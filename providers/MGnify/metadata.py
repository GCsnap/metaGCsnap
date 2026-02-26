import pandas as pd
import bisect
import time
import json
import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# in house imports
import providers.MGnify.parquet_search as ps
import providers.MGnify.helpers as hp

class Metadata():

    ''' this class do all the compute intensive searche son MGnify tabular metadata and updates local ones.'''
    def __init__(self, dataset):

        self.mgyps = dataset.mgyps

        self.MGnify_local = dataset.MGnify_local
        self.output_dir = dataset.output_dir

        self.mgyp_metadata_file = dataset.mgyp_metadata_file
        self.assembly_metadata_file = dataset.assembly_metadata_file
        self.targets_file = dataset.targets_file

        self.metadata_dir = dataset.metadata_dir
    
    def parquetSearch(self, entity, reference, queries):

        # Load delimiters from the specified path.
        reference_path = self.MGnify_local / reference
        delimiters_file = reference_path / 'delimiters.csv'
        delimiters, first_ids = ps.load_delimiters(delimiters_file,entity)

        # Group queries by file and block, storing boundaries.
        queries_by_file, results = ps.group_queries_by_block(queries, delimiters, first_ids, entity=entity)

        start_time = time.time()
        processed_queries = []
        tot_queries = len(queries)
        #tot_queries = len(query_list)
        start_time = time.time()

        for filename, blocks in queries_by_file.items():

            ps.process_file_blocks(reference_path, filename, blocks, results, entity=entity)
            
            processed_queries+=[query for block in blocks.values() for query in block['queries']]
            percentage_done = (len(processed_queries) / tot_queries) * 100
            time_elapsed = time.time() - start_time
            minutes, seconds = divmod(time_elapsed, 60)
            print(f"Processed {len(processed_queries)}/{tot_queries}={percentage_done:.2f}% of {entity} queries. Time elapsed: {int(minutes)}m {seconds:.2f}s", end='\r')

        end_time = time.time()

        print(f"  Time taken: {end_time - start_time:.4f} seconds")

        return results
    
    def get_mgyp_metadata(self):

        def get_p2a_metadata(p2a,protein):

            # for a protein request, get the list of assemblies and metadata
            columns = ['ERZ', 'contig', 'MGYC', 'start', 'end', 'strand', 'type']

            # get data for all the assemblies of the protein
            #metadata = p2a.loc[protein,'contig_id'].split(';')
            metadata = p2a[protein]
            
            if metadata is None:
                return None
            
            metadata = metadata.split(';')

            # Process the metadata list
            processed_data = []
            for item in metadata:
                erz_mgyc, start_end, strand, type_ = item.split(':')
                erz, mgyc = erz_mgyc.split('.')
                start, end = start_end.split('-')
                processed_data.append([erz, item, mgyc, start, end, strand, type_])

            # Create a new DataFrame with the processed data
            metadata = pd.DataFrame(processed_data, columns=columns)
            #metadata=metadata.set_index('ERZ')

            return metadata

        print('Building MGYP metadata')
        # get the protein2assembly mapping
        # we will need to wrap, until final_metadata is defined, 
        # these lines and set up the large tabular data search
        # so far this method does not requires api access
        # the method runs locally
        
        # Load the protein2assembly dictionary from the JSON file
        protein2assembly_file = self.metadata_dir / 'protein2assembly.json'

        if os.path.exists(protein2assembly_file):
            # Load the mgyc2contig dictionary from the JSON file
            with open(protein2assembly_file, 'r') as f:
                protein2assembly = json.load(f)
        else:
            protein2assembly = self.parquetSearch(reference='seq_metadata' ,entity='sequence', queries=self.mgyps)
            with open(protein2assembly_file, 'w') as f:
                json.dump(protein2assembly, f)

        # prepare all the urls to be scraped
        list_target_metadata = []
        missing_mgyps = []
        for p in self.mgyps:
            
            tm = get_p2a_metadata(p2a=protein2assembly, protein=p)

            if tm is not None:
                tm['MGYP'] = p 
                list_target_metadata.append(tm)
            else:
                missing_mgyps.append(p)

        # Concatenate all metadata
        self.mgyp_metadata = pd.concat(list_target_metadata,ignore_index=False)
        self.mgyp_metadata[['start','end']] = self.mgyp_metadata[['start','end']].astype(int)
        self.mgyp_metadata.to_csv( self.mgyp_metadata_file)

        # Save missing MGYPs to a file
        if missing_mgyps:
            missing_mgyp_file = self.metadata_dir / 'missing_mgyp.txt'
            with open(missing_mgyp_file, 'w') as f:
                for mgyp_id in missing_mgyps:
                    f.write(f"{mgyp_id}\n")
            print(f"Missing {len(missing_mgyps)} MGYPs saved to metadata/missing_mgyp.txt")

        # save assemblies
        self.assemblies = self.mgyp_metadata['ERZ'].unique().tolist()
        
    def assign_contigs(self):
        
        print('Assigning ERZ-contigs to MGYCs')
        mgyc2contig_file = self.metadata_dir / 'mgyc2contig.json'

        if os.path.exists(mgyc2contig_file):
            # Load the mgyc2contig dictionary from the JSON file
            with open(mgyc2contig_file, 'r') as f:
                mgyc2contig = json.load(f)
        else:
            mgyc2contig = self.parquetSearch(reference='contig_map' ,entity='contig', queries=self.mgyp_metadata['MGYC'].unique().tolist())
            with open(mgyc2contig_file, 'w') as f:
                json.dump(mgyc2contig, f)

        # Add a new column 'ERZ_contig' by mapping the 'contig' column using the mgyc2contig dictionary
        self.mgyp_metadata['ERZ_contig'] = self.mgyp_metadata['MGYC'].map(mgyc2contig)

        # Save the updated self.mgyp_metadata DataFrame to a CSV file
        self.mgyp_metadata.to_csv(self.mgyp_metadata_file)

        print(f"ERZ-contigs assigned to MGYCs and saved.")

    def assign_cds_from_gff(self,dataset,batch_size=20):
                
        cds_feature = 'ERZ_cds_id'
        self.MGYP2cds_file = self.metadata_dir / 'mgyp2cds.json'

        def inspect_gff(gff_file, table):
            """
            Extracts a sequence from a GFF file based on the provided contig, start, and end positions.
            """

            table[cds_feature]=''

            try:
                col_names = [ "seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
                col_types = { c:int if c in ["start", "end", "phase"] else str for c in col_names }
                gff_df = pd.read_csv(gff_file, sep='\t', header=None, comment='#', names=col_names, dtype=col_types, compression='gzip',)
                
                res = {}
                
                for i, row in table.iterrows():
                    
                    matches = gff_df[ (gff_df['seqid'] == row['ERZ_contig']) & 
                                    (gff_df['start'] == row['start']) & 
                                    (gff_df['end'] == row['end']) ]
                    
                    if not matches.empty:
                        attributes = matches.iloc[0]['attributes']
                        attr_dict = dict(item.split('=') for item in attributes.split(';') if '=' in item)
                        table.loc[i,cds_feature] = attr_dict.get('ID')
                    else:
                        table.loc[i,cds_feature] = ''
                
                return table
            except (EOFError, FileNotFoundError):
                print(f"Error reading or not found GFF file: {gff_file}",end='\r')
                return table

        def process_batch(batch, inspect_func):
            """Process a small batch of tasks."""
            results = []
            for gff_file, table in batch:
                result = inspect_func(gff_file, table)
                results.append(result)
            return results

        def assign_cds_from_json():

            with open(self.MGYP2cds_file) as f:
                raw_dict = json.load(f)

            # reconstruct tuple keys
            lookup_dict = { tuple(k.split(hp.TSEP)): v for k, v in raw_dict.items() }

            def get_key(row):
                TUPLE_COLS = ["MGYP", "ERZ_contig", "start", "end"]
                return tuple(str(row[col]) for col in TUPLE_COLS)

            self.mgyp_metadata[cds_feature] = self.mgyp_metadata.apply(lambda row: lookup_dict.get(get_key(row), None), axis=1)
            self.mgyp_metadata[cds_feature] = self.mgyp_metadata[cds_feature].fillna('NA')

        print('Checking if CDS ids have already been assigned')
        if os.path.exists(self.MGYP2cds_file):
            print('Assigning CDS ids from json')
            assign_cds_from_json()
            self.mgyp_metadata.to_csv( self.mgyp_metadata_file, index=False)
            return
        
        print('Assigning CDS ids from gffs')
        # 1. Prepare tasks (one per assembly)
        tasks = []

        for contig, mgyc_file in dataset.contig_gff_file.items():

            mgyc_table = self.mgyp_metadata[self.mgyp_metadata['MGYC'] == contig ].copy()

            tasks.append((mgyc_file, mgyc_table))
            #inspect_gff(erz_gff_file, erz_table), for testing

        # 2. Parallel batch execution
        total_assemblies = len(tasks)
        processed_assemblies = 0
        start_time = time.time()

        mgyp2cds = []

        with ThreadPoolExecutor() as executor:
            futures = []
            for task_batch in hp.chunked_iterable(tasks, batch_size):
                futures.append(executor.submit(process_batch, task_batch, inspect_gff))
            
            for future in as_completed(futures):
                batch_result = future.result()
                for partial_result in batch_result:
                    mgyp2cds.append(partial_result)
                                        
                processed_assemblies += len(batch_result)
                elapsed_time = time.time() - start_time
                hours, remainder = divmod(elapsed_time, 3600)
                minutes, seconds = divmod(remainder, 60)
                print(f"Processed {processed_assemblies}/{total_assemblies} assemblies. Time elapsed: {int(hours)}h {int(minutes)}m {int(seconds)}s", end='\r')
        
        # 3. Update metadata
        mgyp2cds = pd.concat(mgyp2cds, ignore_index=True)
        
        tuple_mgyp2cds = { (row["MGYP"], row["ERZ_contig"], row["start"], row["end"]): row["ERZ_cds_id"].replace('cds-','') for _, row in mgyp2cds.iterrows()}
        string_mgyp2cds = {hp.TSEP.join(map(str, k)): v for k, v in tuple_mgyp2cds.items()}

        with open(self.MGYP2cds_file, "w") as f:
            json.dump(string_mgyp2cds, f, indent=2)

        assign_cds_from_json()
        self.mgyp_metadata.to_csv( self.mgyp_metadata_file, index=False )