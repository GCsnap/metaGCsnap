import os
import gzip # to work with .gz
import urllib.request
import time
from datetime import datetime, timedelta

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import WarningToLog

from gcsnap.providers.MGnify.dataset import Dataset

import logging
logger = logging.getLogger('iteration')

class Assemblies:    
    """
    Methods and attributes to download and parse flanking genes given NCBI codes.

    Attributes:
        cores (int): Number of cores to use for parallel processing.
        n_flanking5 (int): Number of flanking genes to extract at the 5' end of target.
        n_flanking3 (int): Number of flanking genes to extract at the 3' end of target.
        exclude_partial (bool): Exclude partial genomic blocks.
        config (Configuration): Configuration object.
        console (RichConsole): Console object.
        assembly_dir (str): Path to store assembly summaries. [I already have all the summaries and assembly metadata]
        targets_and_cds_codes (list): List of tuples with target and ncbi code. [I don't need this kind of mapping, as i already have downloaded everything locally]
        accessions (dict): Dictionary with ncbi codes and assembly accessions.
    """

    def __init__(self, config: Configuration, dataset: Dataset):                     
        """
        Initialize the Assemblies object.

        Args:
            config (Configuration): Configuration object containing the arguments.
            mappings (list[tuple[str,str]]): Contains the target and its ncbi code.
        """        
        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value'] 
        self.n_flanking5 = config.arguments['n_flanking5']['value']  
        self.n_flanking3 = config.arguments['n_flanking3']['value']
        self.exclude_partial = config.arguments['exclude_partial']['value'] 
        self.config = config

        self.console = RichConsole()

        # this will be replaced by assembly --> contig
        self.assembly_dir = dataset.contigs_dir

        try:
            self.viable_cds = dataset.viable_cds
        except AttributeError:
            self.viable_cds = dataset.mgyp_metadata['ERZ_cds_id'].unique().tolist()
            
        self.mgyp_metadata = dataset.mgyp_metadata.fillna('not found')
        self.mgyp_metadata = self.mgyp_metadata[self.mgyp_metadata['ERZ_cds_id'].isin(self.viable_cds)]
        self.assemblies = self.mgyp_metadata['ERZ'].unique().tolist()
        
        self.assembly_file = dataset.contig_file
        self.proteins_file = dataset.contig_proteins_file
        self.gff_file = dataset.contig_gff_file

        # input list with [(target, cds_code)] in GCsnap
        # input list with [(mgyp, coding sequence)] in metaGCsnap
        #self.targets_and_cds_codes = mappings

        # this will be replaced by targets_and_cds_codes --> targets_and_cds_codes
        #self.targets_and_cds_codes = [(row['MGYP'], row['ERZ_cds_id']) for _, row in self.mgyp_metadata.iterrows() if row['ERZ_cds_id'] != 'NA']
        
        self.targets_and_cds_codes = [(row['ERZ_cds_id'], row['ERZ_cds_id']) for _, row in self.mgyp_metadata.iterrows() if row['ERZ_cds_id'] != 'not found']
        self.partial_cds = []
        #### new attributes
        self.contig_file = dataset.contig_gff_file
        
    def get_flanking_genes(self) -> dict:
        """
        Getter for the flanking_genes attribute.

        Returns:
            dict: The flanking gene infrmation.
        """        
        return self.flanking_genes
    
    def run(self) -> None:
        """
        Run the process to download and extract flanking genes for the targets:
            - Find the assembly accessions for the given NCBI codes.
            - Download and extract flanking genes for each target in parallel.
        Uses parallel processing with the processpool_wrapper from utils.py
        """        
        # find the assembly accessions
        #cds_codes = [target[1] for target in self.targets_and_cds_codes]
        self.find_assembly_accessions()
        self.find_mgyc_accessions()
        self.find_mgyp_accessions()
        self.find_contig_accessions()
        self.find_urls()
        
        # download and parse the assembly summary files
        with self.console.status('Extract flanking genes from contigs'):
            dict_list = processpool_wrapper(self.cores, self.targets_and_cds_codes, self.run_each)
            # combine results
            self.flanking_genes = {k: v for d in dict_list for k, v in d.items() 
                                   if v.get('flanking_genes') is not None}
            not_found = {k: v for d in dict_list for k, v in d.items() 
                         if v.get('flanking_genes') is None}
        if not_found:
            #pass
            self.log_not_found(not_found)

        if not self.flanking_genes: 
            msg = 'No flanking genes found for any target sequence. Continuing is not possible.'
            self.console.stop_execution(msg = msg)

    def run_each(self, args: tuple[str,str]) -> dict[str, dict]:
        """
        Run the process to download and extract flanking genes for a single target
        used in the parallel processing.

        Args:
            args (tuple[str,str]): Contains the target and its ncbi code.

        Returns:
            dict[str, dict]: The flanking genes and assembly information.
        """        
        target, cds_code = args
        try:
            # get ERZ corresponding to the contig
            # MGYC
            #accession = self.cds_codes_to_contig[cds_code]
            #print('1',cds_code,accession)
            assembly_accession = self.get_assembly_accession(cds_code)
            assembly_url = self.get_assembly_url(assembly_accession)
            
            mgyc_accession = self.get_mgyc_accession(cds_code)
            mgyp_accession = self.get_mgyp_accession(cds_code)
            contig_accession = self.get_contig_accession(cds_code)
            assembly_file, lines = self.download_and_read_gz_file(mgyc_accession)
            flanking_genes = self.parse_assembly(cds_code, lines)
            
            assembly_metadata = {   'source': 'MGnify', 
                                    'target_cds': cds_code, 
                                    'genomic_region': contig_accession, 
                                    'contig': mgyc_accession, 
                                    'assembly_accession': assembly_accession, 
                                    'assembly_url': assembly_url, 
                                    'target': mgyp_accession, 
                                    'target_source': 'MGnify',
                                    'dna_file': self.assembly_file[mgyc_accession].as_posix(),
                                    'protein_file': self.proteins_file[mgyc_accession].as_posix(),
                                    'gff_file': self.gff_file[mgyc_accession].as_posix() }

            return {target: {'flanking_genes': flanking_genes,
                             'assembly_metadata':  assembly_metadata}}   # here, i want erz-contig, assembly, contig?

        except WarningToLog as e:
            # return None for flanking genes and message, logged later
            return {target: {'flanking_genes': None,
                             'msg': str(e)}}

    def find_assembly_accessions(self) -> None:
        """
        Find the assembly accessions for the given NCBI codes using the EntrezQuery class.

        Args:
            cds_codes (list): The list of NCBI codes.
        Sets:
            dict of {cds:contig id}
        """        
        # get file with assembly links, no logging as its done after run()

        self.assembly_accessions = { row['ERZ_cds_id']: row['ERZ'] for _, row in self.mgyp_metadata.iterrows()}

    def get_assembly_accession(self, cds_code: str) -> str:
        """
        Get the assembly accession for a given NCBI code.

        Args:
            cds_code (str): The NCBI code.

        Raises:
            WarningToLog: If no assembly accession is found.

        Returns:
            str: The assembly accession.
        """        
        accession = self.assembly_accessions.get(cds_code, 'unk')
        if accession == 'unk':
            raise WarningToLog('No assembly accession found for {}'.format(cds_code))
        return accession

    def find_mgyp_accessions(self) -> None:
        """
        Find the contig accessions for the given NCBI codes using the EntrezQuery class.

        Args:
            cds_codes (list): The list of NCBI codes.
        Sets:
            dict of {cds:contig id}
        """        
        # get file with assembly links, no logging as its done after run()

        self.mgyp_accession = { row['ERZ_cds_id']: row['MGYP'] for _, row in self.mgyp_metadata.iterrows()}

    def get_mgyp_accession(self, cds_code: str) -> str:
        """
        Get the contig accession for a given NCBI code.

        Args:
            cds_code (str): The NCBI code.

        Raises:
            WarningToLog: If no assembly accession is found.

        Returns:
            str: The assembly accession.
        """        
        accession = self.mgyp_accession.get(cds_code, 'unk')
        if accession == 'unk':
            raise WarningToLog('No assembly accession found for {}'.format(cds_code))
        return accession

    def find_mgyc_accessions(self) -> None:
        """
        Find the contig accessions for the given NCBI codes using the EntrezQuery class.

        Args:
            cds_codes (list): The list of NCBI codes.
        Sets:
            dict of {cds:contig id}
        """        
        # get file with assembly links, no logging as its done after run()

        self.mgyc_accession = { row['ERZ_cds_id']: row['MGYC'] for _, row in self.mgyp_metadata.iterrows()}

    def get_mgyc_accession(self, cds_code: str) -> str:
        """
        Get the contig accession for a given NCBI code.

        Args:
            cds_code (str): The NCBI code.

        Raises:
            WarningToLog: If no assembly accession is found.

        Returns:
            str: The assembly accession.
        """        
        accession = self.mgyc_accession.get(cds_code, 'unk')
        if accession == 'unk':
            raise WarningToLog('No assembly accession found for {}'.format(cds_code))
        return accession

    def find_contig_accessions(self) -> None:
        """
        Find the contig accessions for the given NCBI codes using the EntrezQuery class.

        Args:
            cds_codes (list): The list of NCBI codes.
        Sets:
            dict of {cds:contig id}
        """        
        # get file with assembly links, no logging as its done after run()

        self.contig_accession = { row['ERZ_cds_id']: row['ERZ_contig'] for _, row in self.mgyp_metadata.iterrows()}

    def get_contig_accession(self, cds_code: str) -> str:
        """
        Get the contig accession for a given NCBI code.

        Args:
            cds_code (str): The NCBI code.

        Raises:
            WarningToLog: If no assembly accession is found.

        Returns:
            str: The assembly accession.
        """        
        accession = self.contig_accession.get(cds_code, 'unk')
        if accession == 'unk':
            raise WarningToLog('No assembly accession found for {}'.format(cds_code))
        return accession

    def find_urls(self) -> None:
        """
        Find the assembly URLs for the given assembly accessions using the AssemblyLinks class.

        Sets:
            dict of {accession:url}
        """        
        # get file with assembly links, no logging as its done after run()

        mgnify_url = 'https://www.ebi.ac.uk/metagenomics/assemblies/{}'
        self.links = {a: mgnify_url.format(a) for a in self.assemblies}

    def get_assembly_url(self, assembly_accession: str) -> str:
        """
        Get the assembly URL for a given assembly accession.

        Args:
            assembly_accession (str): The assembly accession.

        Raises:
            WarningToLog: If no URL is found.

        Returns:
            str: The assembly URL.
        """        
        url = self.links.get(assembly_accession, 'unk')   
        if url == 'unk':
            raise WarningToLog('No url found for accession {}'.format(assembly_accession))
        return url      
    
    # this will be replaced by download_and_read_gz_file --> read_gz_file and url-->mgyc
    def download_and_read_gz_file(self, accession: str, retries: int = 3) -> tuple:
        """
        Wrapper to download and read the assembly file in .gz format.

        Args:
            url (str): The URL of the assembly file.
            retries (int, optional): The number of retries of download. Defaults to 3.

        Raises:
            WarningToLog: If the file was not downloaded properly.

        Returns:
            tuple: The full path of the downloaded file and the content of the file.
        """     

        

        try:
            full_path = self.contig_file[accession]
            # this will be replaced by download_gz_file --> full_path --> self.mgyc_file
            content = self.read_gz_file(full_path)
        
        except (urllib.error.URLError, EOFError, FileNotFoundError, KeyError) as e:        
            raise WarningToLog('Failed to open {} with error {}'.format(full_path, e))

        return full_path, content

    def download_gz_file(self, url: str) -> str:
        """
        Download the .gz file and save it without uncompressing.

        Args:
            url (str): The URL of the assembly file.

        Returns:
            str: The full path of the downloaded file.
        """        
        assembly_label = url.split('/')[-1]  # e.g., GCF_000260135.1_ASM26013v1
        assembly_file_gz = '{}_genomic.gff.gz'.format(assembly_label)
        full_path = os.path.join(self.assembly_dir, assembly_file_gz)

        days = self.age
        if os.path.exists(full_path):
            # time check, download again if older than days
            file_time = datetime.fromtimestamp(os.path.getmtime(full_path))
            if datetime.now() - file_time > timedelta(days=days):
                # Download the .gz file and save it without uncompressing
                with urllib.request.urlopen(url + '/' + assembly_file_gz) as response:
                    with open(full_path, 'wb') as out_file:
                        out_file.write(response.read())
            else:
                return full_path
        else:
            # Download the .gz file and save it without uncompressing
            with urllib.request.urlopen(url + '/' + assembly_file_gz) as response:
                with open(full_path, 'wb') as out_file:
                    out_file.write(response.read())

        return full_path

    def read_gz_file(self, file_path: str) -> list:
        """
        Read the content of a .gz file.

        Args:
            file_path (str): The path of the .gz file.

        Returns:
            list: The content of the file as a list of lines.
        """        
        with gzip.open(file_path, 'rt', encoding='utf-8') as file:
            content = file.read()
        return content.splitlines()                
                
    def delete_assemblies(self) -> None:
        """
        Delete the assembly files in the assembly directory.
        """        
        for filename in os.listdir(self.assembly_dir):
            file_path = os.path.join(self.assembly_dir, filename)
            os.remove(file_path)

    def parse_assembly(self, cds_code: str, lines: list) -> dict:
        """
        Wrapper to extract the genomic context block and to extract the flanking genes.

        Args:
            cds_code (str): The NCBI code.
            lines (list): The content of the assembly file.

        Returns:
            dict: The flanking genes.
        """        
        genomic_context_block = self.extract_genomic_context_block(cds_code, lines)
        return self.parse_genomic_context_block(cds_code, genomic_context_block)

    def extract_genomic_context_block(self, target_cds_code: str, lines: list) -> list:
        """
        Extract first all lines belonging to the scaffold containing the target gene.
        Then extract the genomic context block (n_flanking5 on 5' end and 
        n_flanking3 on 3' end based on the direction of the target) from that scaffold.
        One line of the genomic context block looks like this:
            JQ926483.1	Genbank	CDS	1	668	.	+	0
            ID=cds-AFI40896.1;Parent=gene-VP1;Dbxref=NCBI_GP:AFI40896.1;
            Name=AFI40896.1;end_range=668,.;gbkey=CDS;gene=VP1;partial=true;
            product=RNA-dependent RNA polymerase;protein_id=AFI40896.1;start_range=.,1

        Args:
            target_cds_code (str): The NCBI code of the target gene.
            lines (list): The content of the assembly file.

        Returns:
            list: The lines containing the flanking genes.
        """       
         
         # line numbers where different scaffolds (sequence regions) start
        scaffold_positions = [0] + [index for index, val in enumerate(lines) if 
                            val.startswith('##sequence-region')] + [len(lines)]   
        #print('scaffold_positions',scaffold_positions)
        # line number of target
        target_position = [index for index, val in enumerate(lines) 
                           if 'ID=cds-{}'.format(target_cds_code) in val or
                           'Name={}'.format(target_cds_code) in val or
                           'protein_id={}'.format(target_cds_code) in val]

        if not target_position:
            raise WarningToLog('{} not found in'.format(target_cds_code))

        # select scaffold region with target
        region_start_end = [(scaffold_positions[i], scaffold_positions[i + 1]) for i in 
                range(len(scaffold_positions) - 1) if scaffold_positions[i] 
                <= target_position[0] <= scaffold_positions[i + 1]]
        
        # extract the scaffold and lines containing coding sequence (CDS) information
        scaffold = [line for line in lines[region_start_end[0][0]:region_start_end[0][1]] if 
                    (len(line.split('\t')) >= 3 and line.split('\t')[2] == 'CDS')]
        
        # target position in scaffold
        index_of_target = [index for index, val in enumerate(scaffold) 
                           if 'ID=cds-{}'.format(target_cds_code) in val or
                           'Name={}'.format(target_cds_code) in val or
                           'protein_id={}'.format(target_cds_code) in val][0]
                
        # need to know direction to define what flanking genes to extract
        direction_of_target = scaffold[index_of_target].split('\t')[6]
        
        # define neighbor indizes depending on direction
        if direction_of_target == '+':
            # upper slice index can be out of bounds without error, lower must be non negative
            start = 0 if (index_of_target - self.n_flanking5) < 0 else index_of_target - self.n_flanking5
            end = index_of_target + self.n_flanking3 + 1
            # extract the genomic context 
            genomic_context_block = scaffold[start:end]
        else:
            # upper slice index can be out of bounds without error, lower must be non negative
            start = 0 if (index_of_target - self.n_flanking3) < 0 else index_of_target - self.n_flanking3
            end = index_of_target + self.n_flanking5 + 1
            # extract the genomic context 
            genomic_context_block = scaffold[start:end]
            # reverse if direction of target is '-'
            genomic_context_block = genomic_context_block[::-1]
        
        # exclude partials if desired
        if self.exclude_partial and len(genomic_context_block) < (self.n_flanking5 + self.n_flanking3 + 1):
            #self.partial_cds.append(target_cds_code)
            raise WarningToLog('Partial genomic block for {} excluded! Partial CDS: {}'.format(target_cds_code, len(self.partial_cds)))
        
        return genomic_context_block

    def parse_genomic_context_block(self, target_cds_code: str, genomic_context_block: list) -> dict:
        """
        Parse the genmoic context block to extract the flanking genes information.
        One line of the genomic context block looks like this:
            JQ926483.1	Genbank	CDS	1	668	.	+	0
            ID=cds-AFI40896.1;Parent=gene-VP1;Dbxref=NCBI_GP:AFI40896.1;
            Name=AFI40896.1;end_range=668,.;gbkey=CDS;gene=VP1;partial=true;
            product=RNA-dependent RNA polymerase;protein_id=AFI40896.1;start_range=.,1

        Args:
            target_cds_code (str): The NCBI code of the target gene.
            genomic_context_block (list): The lines containing the flanking genes.

        Returns:
            dict: The extracted flanking genes information.
        """
        # result dictionary
        flanking_genes = GenomicContext.get_empty_flanking_genes()

        # parse the genomic context
        for line in genomic_context_block:   

            line_data = line.split('\t')            
            start = int(line_data[3])
            end = int(line_data[4])
            direction = line_data[6]
            
            if 'cds-' in line:
                cds_code = line_data[8].split('ID=cds-')[1].split(';')[0]
            elif 'Name=' in line:
                cds_code = line_data[8].split('Name=')[1].split(';')[0]
            else:
                cds_code = 'unk'

            if 'pseudo=' not in line and 'product=' in line and 'fragment' not in line:
                prot_name = line_data[8].split('product=')[1].split(';')[0]
            else:
                prot_name = 'pseudogene'

            # it means that this is some kind of fragmented gene (has introns?...) 
            # and so we have to collect the largest interval encompassing it
            if cds_code in flanking_genes['cds_codes'] and flanking_genes['cds_codes'][-1] == cds_code: 
                if start < flanking_genes['starts'][-1]:
                    flanking_genes['starts'][-1] = start
                if end > flanking_genes['ends'][-1]:
                    flanking_genes['ends'][-1] = end
            else:
                if '|' in cds_code:
                    cds_code = cds_code.replace('|','_')

                flanking_genes['cds_codes'].append(cds_code)
                flanking_genes['names'].append(prot_name)
                flanking_genes['starts'].append(start)
                flanking_genes['ends'].append(end)
                flanking_genes['directions'].append(direction)
        
        # index of target in flanking genes
        index_of_target = flanking_genes['cds_codes'].index(target_cds_code)

        # direction of target
        direction_of_target = flanking_genes['directions'][index_of_target]
        
        # add relative starts and ends depending on direction
        if direction_of_target == '+':
            for key in ['starts','ends']:
                lst = [e - flanking_genes['starts'][index_of_target] + 1 for e in flanking_genes[key]]
                flanking_genes['relative_{}'.format(key)] = lst  
        else:
            # order is reversed, old ends determin the starts and vice-versa
            # old end of target is the new base point
            lst = [flanking_genes['ends'][index_of_target] - e + 1 for e in flanking_genes['ends']]
            flanking_genes['relative_starts'] = lst  
            lst = [flanking_genes['ends'][index_of_target] - e + 1 for e in flanking_genes['starts']]
            flanking_genes['relative_ends'] = lst  
            
            # turn + and - directions
            flanking_genes['directions'] = ['+' if d == '-' else '-' for d in flanking_genes['directions']]
            
        return flanking_genes       

    def log_not_found(self, not_found: dict) -> None:
        """
        Write the targets for which no flanking genes were found to the log file.

        Args:
            not_found (dict): The targets for which no flanking genes were found.
        """        
        message = 'No flanking genes found for {} target sequences.'.format(len(not_found))
        self.console.print_warning(message)
        for k,v in not_found.items():
            logger.warning('For target {}: {}'.format(k,v.get('msg'))) 
