import pandas as pd
import numpy as np
import os 

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.utils import processpool_wrapper
from providers.ncbi.uniprot_dbs_dict import UniprotDict

from providers.ncbi.uniprot_api import submit_id_mapping
from providers.ncbi.uniprot_api import check_id_mapping_results_ready
from providers.ncbi.uniprot_api import get_id_mapping_results_link
from providers.ncbi.uniprot_api import get_id_mapping_results_search
from providers.ncbi.dataset import Dataset

import logging
logger = logging.getLogger('iteration')

class SequenceMapping:
    """ 
    Methods and attributes to map target sequences between different id standards
    using the UniProt Rest API. The API supports mapping of different id standards, 
    e.g., UniParc, UniProtKB-AC, RefSeq, EMBL-CDS, etc. 
    However, not all id standards are supported for mapping to each other.
    The standard defines the databases to be used for mapping. A full list of
    available databases can be found at:
    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        target_list (list): The list of target sequences.
        to_type (str): The id standard to map to.
        msg (str): The message to print.
        cores (int): The number of CPU cores to use.
        console (RichConsole): The RichConsole object to print messages.
        uniprot_dict (dict): The dictionary with the identifiers standard information and the target sequences.
        supported (list): The list of supported databases.
        target_types (list): The list of target types.
        mapping_df (pd.DataFrame): The DataFrame with the mapping results.
    """

    def __init__(self, dataset: Dataset, config: Configuration, target_list: list, to_type: str, msg: str = None):
        """
        Initialize the SequenceMapping object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            target_list (list): The list of target sequences.
            to_type (str): The id standard to map to.
            msg (str, optional): The console message to display. Defaults to None.
        """        
        self.target_list = target_list
        self.to_type = to_type
        if msg is not None:
            self.msg = 'Mapping {} to {}'.format(msg, to_type)
        else:
            self.msg = 'Mapping target sequences to {}'.format(to_type)

        # extract needed values from config
        self.out_dir = dataset.metadata_dir
        self.cores = config.arguments['n_cpu']['value']
        self.console = RichConsole()

    def run(self):
        """
        Run the mapping of target sequences to the specified id standard:
            - Create the UniProt dictionary with the target sequences.
            - Create the parallel input for the mapping.
            - Run the mapping.
            - Finalize the mapping.
        Uses parallel processing with processpool_wrapper from utils.py.
        """        
        with self.console.status(self.msg):
            self.create_uniprot_dbs_with_targets()
            parallel_args = self.create_parallel_input()
            mapping_df_list = processpool_wrapper(self.cores, parallel_args, self.do_mapping)
            # concatenate all mapping dfs to one
            self.mapping_df = pd.concat(mapping_df_list, ignore_index=True)
            # add target column
            self.add_target_column()

    def finalize(self): 
        """
        Finalize the mapping of target sequences to the specified database:
            - Add the NCBI column.
            - Log the failed mappings.
            - Save the mapping DataFrame to a CSV file.        
        """        
        # possible existing target column is overwritten
        self.add_target_column()
        self.add_ncbi_column()
        self.log_failed_finalize()

        # write dataframe to csv file
        self.mapping_df.to_csv(os.path.join(self.out_dir,'mapping.csv'), index=False) 
        self.console.print_done('All mapping done. Table saved to mapping.csv')
    
    def get_codes(self, id_type: str = None) -> list:
        """
        Get the list of codes from the mapping DataFrame.

        Args:
            id_type (str, optional): The id standard of the codes to get. Defaults to None using to_type attribute.

        Returns:
            list: The list of codes.
        """        
        if id_type is None:
            id_type = self.to_type        
        return self.mapping_df[id_type].dropna().tolist()
    
    def get_target_to_result_dict(self) -> dict:
        """
        Get a dictionary with the target sequences as keys and the mapped results as values.

        Returns:
            dict: The dictionary with the target sequences and the mapped results.
        """        
        df = self.mapping_df[(self.mapping_df['target'].notna()) &
                              (self.mapping_df[self.to_type].notna())] 
        return df.set_index('target')[self.to_type].to_dict()
    
    def get_targets_and_ncbi_codes(self) -> list[tuple]:
        """
        Get a list of tuples with the target sequences and the NCBI codes.

        Returns:
            list[tuple]: The list of tuples with the target sequences and the NCBI codes.
        """        
        # return list of tuples with target and ncbi_code
        df = self.mapping_df[(self.mapping_df['target'].notna()) &
                              (self.mapping_df['ncbi_code'].notna())]        
        return list(zip(df['target'], df['ncbi_code']))
           
    def create_uniprot_dbs_with_targets(self) -> None:
        """
        Create the UniProt dictionary and assign the target sequences to the different type standards.
        """        
        uniprot_dict = UniprotDict()
        # fill uniprot_dbs with the targets        
        uniprot_dict.assign_target_to_type(self.target_list)
        self.supported = uniprot_dict.get_supported_databases()
        self.target_types = uniprot_dict.get_target_types()
        self.uniprot_dict = uniprot_dict.get_uniprot_dict()

    def create_parallel_input(self, api_limit: int = 1000) -> list[tuple]:
        """
        Create the parallel input for the mapping of target sequences to the specified database.

        Args:
            api_limit (int, optional): The API limit for the number of IDs to map. Defaults to 1000.

        Returns:
            list[tuple]: The list of tuples with the parallel input arguments.
        """        
        parallel_input_args = []

        # 1. Create a parallel workload for each uniprot type
        # This assures that we can do all of the same type with one request
        for key, values in self.uniprot_dict.items():
            if len(values['targets']) != 0:
                from_type = key
                ids = values['targets']
                from_db = values['from_dbs']
                to_db = self.uniprot_dict[self.to_type]['to_dbs']
                
                parallel_input_args.append((ids, from_db, to_db, from_type, self.to_type))
                
        # 2. Split the ids into sublists of more than api_limit
        # This assures that we do not exceed the api limit
        refinded_parallel_input_args = []
        additions = []
        ids_sublists = []
        for parallel_arg in parallel_input_args:
            ids, from_db, to_db, from_type, _ = parallel_arg    

            if len(ids) > api_limit:
                ids_sublists += self.split_recursive(ids, api_limit)
                additions += [(sublist, from_db, to_db, from_type, self.to_type)
                                        for sublist in ids_sublists]
            else:
                refinded_parallel_input_args.append(parallel_arg)

        parallel_input_args = refinded_parallel_input_args + additions
        return parallel_input_args
    
    def split_recursive(self, ids: list, api_limit: int) -> list:
        """
        Split the IDs into sublists of more than the API limit.

        Args:
            ids (list): The list of IDs to split.
            api_limit (int): The API limit for the number of IDs to map.

        Returns:
            list: The list of sublists of IDs.
        """        
        # Base case: if the length of ids is less than or equal to api_limit, 
        # return a list containing ids
        if len(ids) <= api_limit:
            # the trick is to return [ids], a list of a list which we can iterate through
            return [ids]
        
        # Recursive case: split ids into two halves and recursively split each half
        mid = len(ids) // 2
        return self.split_recursive(ids[:mid], api_limit) + self.split_recursive(ids[mid:], api_limit)         

    def do_mapping(self, arg: tuple[list,str,str,str,str]) -> pd.DataFrame:
        """
        Do the mapping of target sequences to the specified database used
        in the parallel processing:
            - Run the mapping API.
            - Extract the mapping from the API results.
            - Create the mapping DataFrame.
            - Make the mapping DataFrame unique.

        Args:
            arg (tuple[list,str,str,str,str]): The tuple with the arguments.
                First element is the list of IDs to map.
                Second element is the database of the id standard to map from.
                Third element is the database of the id standard to map to.
                Fourth element is the database id standard to map from.
                Fifth element is the database id standard to map to.

        Returns:
            pd.DataFrame: The DataFrame with the unique mapping results.
        """ 
        ids, from_db, to_db, from_type, to_type = arg
        
        api_results = self.run_mapping_api(ids, from_db, to_db)
        mapping_extracted = self.extract_mapping_from_api_results(ids, api_results, from_type, to_type)
        mapping_df = self.create_mapping_df(mapping_extracted)
        unique_mapping_df = self.make_unique(mapping_df, [from_type], to_type)
        
        return unique_mapping_df

    def run_mapping_api(self, ids: list, from_db: str, to_db: str) -> list:  
        """
        Run the mapping of target sequences to the specified database using the UniProt
        API in uniprot_api.py.

        Args:
            ids (list): The list of IDs to map.
            from_db (str): The database of the id standard to map from.
            to_db (str): The database of the id standard to map to.

        Returns:
            list: The list of mapping results.
        """        
        #ids, from_db, to_db, from_type, to_type = arg
        # call methods from gcsnap.uniprot_api.py
        job_id = submit_id_mapping(from_db=from_db, to_db=to_db, ids=ids)
        
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)['results']

        if isinstance(results, list):
            return results
        else:           
            return None

    def extract_mapping_from_api_results(self, ids: list,  results: list, from_type: str , to_type: str) -> dict:
        """
        Extract the mapping results from the API results.

        Args:
            ids (list): The list of IDs to map.
            results (list): The list of mapping results.
            from_type (str): The database of the id standard to map from.
            to_type (str): The database of the id standard to map to.

        Returns:
            dict: _description_
        """        
        # the UniProt Rest API returns different formats depending on what db is searched
        # for instance when mapping to EMBL-Genbank, the results is a list with 1 level dict:
            # [{'from': 'A0A0E3MFP2', 'to': 'AKA75089.1'},
            # {'from': 'A0A0E3MFP2', 'to': 'AKA77782.1'}]
        # mapping to UniProtKB returns a list with nested dict:
            # [{'from': 'P05067', 'to': {'entryType': 'UniProtKB reviewed (Swiss-Prot)', 'primaryAccession': 'P05067', 
            
        extracted = {'from_type': from_type,
                    'to_type': to_type,
                    'from_id': [],
                    'to_id': []}
        
        if results is None:
            extracted['from_id'] += ids
            extracted['to_id'] += [pd.NA] * len(ids)
            return extracted
                
        for result in results:
            extracted['from_id'].append(result['from'])
            # get result_id
            if isinstance(results[0]['to'], dict):
                extracted['to_id'].append(result['to']['primaryAccession'])
            else:
                extracted['to_id'].append(result['to'])

        # fill not found ids with nan
        not_found_ids = list(set(ids).difference(set(extracted['from_id'])))
        extracted['from_id'] += not_found_ids
        extracted['to_id'] += [pd.NA] * len(not_found_ids)
            
        return extracted    
        
    def create_mapping_df(self, mapping_dict: dict) -> pd.DataFrame:     
        """
        Create a DataFrame from the mapping dictionary.

        Args:
            mapping_dict (dict): The mapping dictionary.

        Returns:
            pd.DataFrame: The DataFrame with the mapping results.
        """           
        id_types =  self.supported
        mappings = {id_type: [] for id_type in id_types}
        
        from_type = mapping_dict['from_type']
        mappings[from_type] = mapping_dict['from_id']
        to_type = mapping_dict['to_type']
        mappings[to_type] = mapping_dict['to_id']
        
        # fill all cols with nan
        for id_type in list(set(id_types).difference(set([from_type,to_type]))):
            mappings[id_type] += [pd.NA] * len(mapping_dict['from_id'])
                
        return pd.DataFrame.from_dict(mappings)

    def make_unique(self, mapping_df: pd.DataFrame, from_type: list, to_type: str) -> pd.DataFrame:
        """
        Make the mapping DataFrame unique by removing duplicates. If there is a tie of same duplicates,
        select the one with the shorter 'to_type' id.

        Args:
            mapping_df (pd.DataFrame): The DataFrame with the mapping results.
            from_type (list): The list of database id standards to map from.
            to_type (str): The database id standard to map to.

        Returns:
            pd.DataFrame: The DataFrame with the unique mapping results.
        """        
        # make copy before handling data frame
        df = mapping_df.copy()

        # add length colum for sorting
        # Rule of thumb, just take the first one for UniParc.            
        # The longer UniProt-KB are currated automatically. Take the samller one
        # pd.NA is a float, hence len() directly is not possible
        df['Length'] = df[to_type].apply(lambda x: len(x) if pd.notna(x) else 0)
        # sort, default ascending
        df = df.sort_values(by='Length')         

        # do it for every column in from_type
        for col in from_type: 
            # is duplicated detects Nan as duplicated, avoid this
            mask = df[col].notna()
            # Checking for duplicates only on non-NaN values while preserving the original index
            is_duplicated = df.loc[mask, col].duplicated(keep=False).reindex(df.index, fill_value=False)
            # DataFrame of non-unique entries, selecting the first occurrence of each duplicated entry
            # sorting is not changed by group by
            first_non_unique_df = df[is_duplicated].groupby(col, sort=False).head(1)
            # DataFrame of unique entries (i.e., entries that are not duplicated)
            unique_df = df[~is_duplicated]
            # Combine the two DataFrames
            filtered_unique_df = pd.concat([unique_df, first_non_unique_df], ignore_index=True)
        
        # drop sorting column
        filtered_unique_df.drop(columns='Length', inplace=True)   
        return filtered_unique_df
    
    def add_target_column(self) -> None:
        """
        Add the target column to the mapping DataFrame and fill it with the target sequences.
        """        
        # add target column with NA or fill it with NA if present
        self.mapping_df['target'] = pd.NA
        for key in self.target_types:
            targets = self.uniprot_dict[key]['targets']
            # Update the 'target' column with the value from 'key' column where 'key' is in targets
            self.mapping_df.loc[self.mapping_df[key].isin(targets), 'target'] = self.mapping_df[key]

    def add_ncbi_column(self) -> None:
        """
        Add the NCBI column to the mapping.
        """        
        self.mapping_df['ncbi_code'] = pd.NA

        # Update 'ncbi_code' with values from 'RefSeq' if not NA, otherwise from 'EMBL-CDS'
        # np.where(condition, x, y) returns x if condition is True, otherwise y
        self.mapping_df['ncbi_code'] = np.where(self.mapping_df['RefSeq'].notna(), 
                                                self.mapping_df['RefSeq'], self.mapping_df['EMBL-CDS'])

    def log_failed_finalize(self) -> None:
        """
        Log the failed mappings and finalize the mapping of target sequences to the specified database.
        """        
        df = self.mapping_df.copy()

        # RefSeqs where no UniProtKB-AC was found
        ref_no_uniprot = df.loc[(df['RefSeq'].notna()) & (df['UniProtKB-AC'].isna()) ,'target'].to_list()
        if len(ref_no_uniprot) > 0:
            message = '{} RefSeq ids not mapped to UniProtKB-AC but included in NCBI-Code.'.format(len(ref_no_uniprot))
            self.console.print_warning(message)
            for id in ref_no_uniprot:
                logger.warning(f'Target sequence {id}') 

        # those where no UniProtKB-AC was found
        no_uniprot = df.loc[df['UniProtKB-AC'].isna(), 'target'].to_list()
        no_uniprot = list(set(no_uniprot).difference(set(ref_no_uniprot)))
        if len(no_uniprot) > 0:
            message = '{} ids not mapped to UniProtKB-AC.'.format(len(no_uniprot))
            self.console.print_warning(message)
            for id in no_uniprot:
                logger.warning(f'Target sequence {id}')   

        # those where no RefSeq or EMBL-CDS was found:
        no_ncbi = df.loc[df['ncbi_code'].isna(), 'target'].to_list()
        no_ncbi = list(set(no_ncbi).difference(set(no_uniprot)))
        if len(no_ncbi) > 0:
            message = '{} ids not mapped to any NCBI-Code (RefSeq or EMBL-CDS).'.format(len(no_ncbi))
            self.console.print_warning(message)
            for id in no_ncbi:
                logger.warning(f'Target sequence {id}')          

    def log_failed(self) -> None:
        """
        Log the failed mappings of target sequences to the specified id standard.
        """        
        df = self.mapping_df.copy()
        no_hits = df.loc[df[self.to_type].isna(), 'target'].to_list()
        if len(no_hits) > 0:
            message = '{} ids not mapped to {}.'.format(len(no_hits),self.to_type)
            self.console.print_warning(message)
            for id in no_hits:
                logger.warning(f'Target sequence {id}')            

    def merge_mapping_dfs(self, mapping_df: pd.DataFrame, key_column: str = 'UniProtKB-AC',
                          columns_to_merge: list = ['EMBL-CDS']) -> pd.DataFrame:
        """
        Merge the mapping DataFrames on the key column and update the nan columns 
        in the first DataFrame with values from the second DataFrame.

        Args:
            mapping_df (pd.DataFrame): The DataFrame to merge with.
            key_column (str, optional): The key column to merge on. Defaults to 'UniProtKB-AC'.
            columns_to_merge (list, optional): The columns to merge. Defaults to ['EMBL-CDS'].

        Returns:
            pd.DataFrame: The merged DataFrame with all mappings from all id standards.
        """        
        # Merge DataFrames on 'UniProtKB-AC'
        merged_df = pd.merge(self.mapping_df, mapping_df[[key_column] + columns_to_merge], 
                             on=key_column, how='left', suffixes=('', '_y'))
        # Update nan columns in df1 with values from df2
        for col in columns_to_merge:
            merged_df[col] = merged_df[col].fillna(merged_df[f'{col}_y'])
            # Drop the additional '_y' column
            merged_df.drop(columns=[f'{col}_y'], inplace=True)
        self.mapping_df = merged_df.copy()