# Connection to MGnify API
# Allows us to query specific values in given fields (e.g.: 'taxon-lineage').

# Standard library imports
import os
import time
import random
import json
import gzip
import zlib
import argparse
from pathlib import Path
from itertools import islice
from threading import Lock
from urllib.request import urlretrieve

# Third-party libraries
import requests
import pandas as pd
import tqdm
from Bio import SeqIO

from jsonapi_client import Session as APISession, Modifier

import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
# In-house libraries
from gcsnap.configuration import Configuration
from gcsnap.providers.MGnify.dataset import Dataset
import gcsnap.providers.MGnify.helpers as hp
from gcsnap.utils import handle_compressed_fasta

class Scraper:

    def __init__(self, dataset: Dataset, config: Configuration):
        
        self.MGnify_API = hp.MGnify_API
        self.output_dir = dataset.output_dir

        ############################################################## read MGPY metadata

        self.targets_file = dataset.targets_file
        self.assembly_metadata_file = dataset.assembly_metadata_file

        self.mgyp_metadata_file = dataset.mgyp_metadata_file 
        self.mgyp_present = self.mgyp_metadata_file.exists()
        if self.mgyp_present:
            self.mgyp_metadata = pd.read_csv(self.mgyp_metadata_file, index_col='ERZ')
        else: 
            raise ValueError(f"MGnify metadata file not found: {self.mgyp_metadata_file}")

        self.assemblies_dir  = dataset.assemblies_dir
        self.proteins_dir = dataset.proteins_dir
        self.assemblies_metadata_dir = dataset.assemblies_metadata_dir
        self.assembly_targets_dir = dataset.assembly_targets_dir
        self.assembly_allfiles_dir = dataset.assembly_allfiles_dir

        ############################################################## target list assemblies
        self.failed_assemblies_file = dataset.failed_assemblies_file
        self.failed_assemblies = []

        if self.failed_assemblies_file.exists():
            with open(self.failed_assemblies_file, 'r') as f:
                self.failed_assemblies = [line.strip() for line in f if line.strip()]

        self.assemblies = dataset.assemblies

        self.failed_assemblies_set = set(self.failed_assemblies)
        
        self.assemblies = [
            assembly for assembly in self.assemblies
            if assembly not in self.failed_assemblies_set
        ]
        
        self.proxies = config.arguments['MGnify_proxies']['value']

    def update_metadata(self):
        
        ''' To be runned after aggregate_assemblies_metadata(). so that the api_metadata and assembly_metadata can be read.'''
        self.targets_file_present = self.targets_file.exists()
        if self.targets_file_present:
            self.targets = pd.read_csv(self.targets_file)
        else:
            raise ValueError("api_metadata.csv file not found.")

        self.assembly_present = self.assembly_metadata_file.exists()
        if self.assembly_present:
            self.assembly_metadata = pd.read_csv(self.assembly_metadata_file)
        else:
            raise ValueError("assembly_metadata.csv file not found.")

    def _robust_api_get(self,url, session, max_retries=5, wait_seconds=10):
        """Robustly GET data from an API endpoint, with retries and wait."""
        for attempt in range(max_retries):
            try:
                response = session.get(url)
                response.raise_for_status()
                return response.json()
            except (requests.exceptions.RequestException, ValueError) as e:
                print(f"Attempt {attempt+1}/{max_retries} failed for {url}: {e}")
                print(f"Waiting {wait_seconds} seconds before retrying...")
                time.sleep(wait_seconds)
        print(f"Giving up on {url} after {max_retries} attempts.")
        return {}

    def _robust_paginated_get(self,start_url, session, max_retries=5, wait_seconds=10):
        """Robustly GET paginated data from an API endpoint."""
        data = []
        next_url = start_url
        while next_url:
            json_data = self._robust_api_get(next_url, session, max_retries, wait_seconds)
            if not json_data or 'data' not in json_data:
                break
            data.extend(json_data['data'])
            next_url = json_data.get('links', {}).get('next')
            time.sleep(0.5)
        return data

    def get_assembly_info(self,assembly):
        """
        Fetches and processes metadata for a given assembly from the MGnify API,
        using robust error handling and retries.
        """
        assembly_file = self.assemblies_metadata_dir / f'{assembly}.csv'
        if assembly_file.exists():
            assembly_info = pd.read_csv(assembly_file, index_col='assembly')
            return assembly_info

        req_session = requests.Session()
        req_session.proxies.update(self.proxies)

        # Analyses (paginated)
        analysis_endpoint = f'{self.MGnify_API}/assemblies/{assembly}/analyses'
        analysis_data = self._robust_paginated_get(analysis_endpoint, req_session)
        if not analysis_data:
            print(f"Failed to retrieve data for assembly {assembly}")
            return pd.DataFrame()

        assembly_df = pd.json_normalize(analysis_data)
        assembly_df = assembly_df.rename(columns={'id': 'analysis.id'})
        assembly_df = assembly_df.sort_values(by='attributes.pipeline-version', ascending=False).head(1)

        # Study
        study_id = assembly_df.loc[assembly_df.index[0], 'relationships.study.data.id']
        study_url = f'{self.MGnify_API}/studies/{study_id}'
        study_res = self._robust_api_get(study_url, req_session)
        if not study_res or 'data' not in study_res:
            print(f"Failed to retrieve study {study_id}")
            study_df = pd.DataFrame()
        else:
            study_df = pd.json_normalize(study_res['data']).rename(columns={'id': 'study.id'})

        # Sample
        sample_id = assembly_df.loc[assembly_df.index[0], 'relationships.sample.data.id']
        sample_url = f'{self.MGnify_API}/samples/{sample_id}'
        sample_res = self._robust_api_get(sample_url, req_session)
        if not sample_res or 'data' not in sample_res:
            print(f"Failed to retrieve sample {sample_id}")
            sample_df = pd.DataFrame()
        else:
            sample_metadata = sample_res['data']['attributes'].get('sample-metadata', [])
            if sample_metadata:
                sample_df = pd.DataFrame([{m['key']: m['value'] for m in sample_metadata}])
            else:
                sample_df = pd.DataFrame()

        assembly_info = pd.concat([assembly_df.reset_index(drop=True), sample_df, study_df.reset_index(drop=True)], axis=1)
        assembly_info['assembly'] = assembly
        assembly_info = assembly_info.set_index('assembly')
        assembly_info.to_csv(assembly_file, index=True)

        return assembly_info

    def get_analysis_results(self,assembly_info, assembly, mgyp):
        """
        Fetches and processes analysis results for a given assembly from the MGnify API,
        using robust error handling and retries.
        """
        analysis_file = self.assembly_targets_dir / f'{assembly}.csv'
        if analysis_file.exists():
            analysisRes = pd.read_csv(analysis_file, index_col=0)
            return analysisRes

        req_session = requests.Session()
        req_session.proxies.update(self.proxies)

        # Downloads (paginated)
        
        analysis = assembly_info['analysis.id'].iloc[0]
        analysis_endpoint = f'{self.MGnify_API}/analyses/{analysis}/downloads'
        analysis_data = self._robust_paginated_get(analysis_endpoint, req_session)
        if not analysis_data:
            print(f"Failed to retrieve downloads for analysis {analysis}")
            return pd.DataFrame()

        analysisRes = pd.json_normalize(analysis_data)
        
        pipeline_feature = 'relationships.pipeline.data.id'
        description_feature = 'attributes.description.label'

        analysisRes.to_csv( self.assembly_allfiles_dir / f'{assembly}.csv' )

        MGnify_pipeline = analysisRes[pipeline_feature].dropna()
        if not MGnify_pipeline.empty:
            version = MGnify_pipeline.loc[MGnify_pipeline.index[0]]
            descriptions = {'5.0': hp.pipelineV5_targets, '4.1': hp.pipelineV4_targets}
            target_descriptions = descriptions.get(version)
            if target_descriptions:
                file_locator = analysisRes[description_feature].isin(target_descriptions)
                analysisRes = analysisRes[file_locator]

        file_endpoint = f'{self.MGnify_API}/analyses/{analysis}/file'
        files_dict = {alias: f'{file_endpoint}/{alias}' for alias in analysisRes['attributes.alias']}
        analysisRes['ftp.url'] = analysisRes['attributes.alias'].map(files_dict)
        analysisRes['MGYP'] = mgyp
        analysisRes['ERZ'] = assembly
        analysisRes.to_csv(analysis_file, index=True)
        return analysisRes

    def aggregate_assemblies_metadata(self):
        """
        Aggregate all files from assembly_info and api_files folders.
        Save the resulting DataFrames and assign them to instance variables.
        """

        self.assemblies = [
            assembly for assembly in self.assemblies
            if assembly not in self.failed_assemblies_set
        ]        

        # Aggregate assembly_info files
        if not os.path.exists(self.assembly_metadata_file):
            self.assembly_metadata_dir = self.output_dir / 'assemblies' / 'metadata'
            assembly_info_files = [self.assembly_metadata_dir / f'{assembly}.csv' for assembly in self.assemblies]
            assembly_info_dfs = [pd.read_csv(file, index_col=0) for file in assembly_info_files]
            self.assembly_metadata = pd.concat(assembly_info_dfs, ignore_index=False)
            self.assembly_metadata.to_csv(self.assembly_metadata_file, index=True)
        else:
            self.assembly_metadata = pd.read_csv(self.assembly_metadata_file, index_col=0)

        print(f"Aggregated assembly metadata files [metadata/assembly_metadata.csv] saved.")

        # Aggregate target files
        if not os.path.exists(self.targets_file):
            self.assembly_targets_dir = self.output_dir / 'assemblies' / 'targets'
            target_files = [self.assembly_targets_dir / f'{assembly}.csv' for assembly in self.assemblies]
            target_files_dfs = [pd.read_csv(file, index_col=0) for file in target_files]
            self.assembly_targets = pd.concat(target_files_dfs, ignore_index=False)
            self.assembly_targets.to_csv(self.targets_file, index=False)
        else:
            self.assembly_targets = pd.read_csv(self.targets_file, index_col=0)

        print(f"Aggregated assembly target files [metadata/target_files.csv] saved.")
        
    def process_assembly(self, assembly):
        
        mgyps = self.mgyp_metadata.loc[assembly, 'MGYP']
        
        if isinstance(mgyps, pd.Series):
            mgyps = mgyps.tolist()
            mgyps = ','.join(mgyps)
        
        assembly_info = self.get_assembly_info(assembly)

        if assembly_info.empty:
            if assembly not in self.failed_assemblies_set:
                self.failed_assemblies.append(assembly)
                self.failed_assemblies_set.add(assembly)
                with open(self.failed_assemblies_file, 'a') as f:
                    f.write(f"{assembly}\n")
            print(f"Failed to retrieve data for assembly {assembly}. Added to failed list.")
            return
        else:
            analysis_files = self.get_analysis_results(assembly_info, assembly, mgyps)
            return
        # Function to process a single assembly using a given proxy
    

        # Helper to create batches from a list
        
    def list_target_files(self):

        def chunks(lst, n):
            for i in range(0, len(lst), n):
                yield lst[i:i + n]

        # Read proxy list from file
        #with open("supplementary/proxies.txt") as f:
        #    raw_proxies = [line.strip() for line in f if line.strip()]
            #proxies = [{'proxies': {'http': proxy, 'https': proxy}} for proxy in raw_proxies]
        proxies = self.proxies
        if not proxies:
            raise ValueError("No proxies found in supplementary/proxies.txt")

        batch_size = len(proxies)
        list_api_metadata=[]
        list_ftp_metadata=[]

        # Sequential version for debugging
        total_assemblies = len(self.assemblies)

        self.assemblies.sort()

        with tqdm.tqdm(total=total_assemblies, desc="Processing assemblies") as pbar:
            for i, assembly in enumerate(self.assemblies):

                start_batch_time = time.time()
                self.process_assembly(assembly)
                elapsed_time = time.time() - start_batch_time

                pbar.update(1)

    def download_quality_check(self, batch_size=32):

        def is_valid_gzip(file_path):
            """Return True if file_path is a valid gzip file."""
            try:
                with gzip.open(file_path, 'rb') as f:
                    while f.read(1024 * 1024):
                        pass
                return True
            except (OSError, EOFError, zlib.error):
                self.targets['corrupted'] = True
                return False

        def check_batch(rows, lock):
            for idx, row in rows:
                file_path = row['local file']
                if is_valid_gzip(file_path):
                    with lock:
                        self.targets.loc[idx, 'checked'] = True
                else:
                    print(f"Corrupted {file_path}")

        print('Checking downloaded files')

        self.assembly_present = self.assembly_metadata_file.exists()
        if self.assembly_present:
            self.assembly_metadata = pd.read_csv(self.assembly_metadata_file)

        if 'downloaded' not in self.targets.columns or 'local file' not in self.targets.columns:
            raise ValueError("Missing 'downloaded' or 'local file' columns in api_metadata. Something went wrong during download.")

        # Add column 'checked' initialized to False
        self.targets['checked'] = False
        self.targets['corrupted'] = False
        lock = Lock()
        
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = []
            total_files = len(self.targets)
            checked_files = 0
            for chunk in hp.chunked_iterable(self.targets.iterrows(), batch_size):
                futures.append(executor.submit(check_batch, chunk, lock))
                concurrent.futures.wait(futures)
                checked_files += len(chunk)
                print(f"Checked {checked_files}/{total_files} files.", end='\r')

        # Save the updated api_metadata DataFrame to its file
                
        self.targets.to_csv(self.targets_file, index=False)

        #n_corrupted_files = self.targets[self.targets['checked'] == False].shape[0]

