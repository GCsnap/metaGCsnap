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
from tqdm import tqdm
import re

# Third-party libraries
import requests
import pandas as pd
from Bio import SeqIO

from jsonapi_client import Session as APISession, Modifier

import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
# In-house libraries
from gcsnap.configuration import Configuration
from providers.MGnify.dataset import Dataset
import providers.MGnify.helpers as hp
from gcsnap.utils import handle_compressed_fasta

def extract_sequences_worker(source_file, target_map):
    """
    Worker function for parallel extraction.
    Returns:
        written_files (list): Paths of files successfully written.
        found_keys (set): The specific keys (substrings) from target_map that were found.
    """
    if not os.path.exists(source_file):
        return [], set()

    # Filter out targets that already exist
    active_targets = {k: v for k, v in target_map.items() if not os.path.exists(v)}
    
    if not active_targets:
        return [], set()

    buffers = {out_path: [] for out_path in active_targets.values()}
    found_keys = set()
    written_files = []

    try:
        with gzip.open(source_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                for substring, out_path in active_targets.items():
                    if substring in record.id:
                        buffers[out_path].append(record)
                        found_keys.add(substring)
                        break
        
        # Write buffered records
        for out_path, records in buffers.items():
            if records:
                try:
                    with gzip.open(out_path, "wt") as out_handle:
                        SeqIO.write(records, out_handle, "fasta")
                    written_files.append(out_path)
                except Exception as e:
                    print(f"Error writing {out_path}: {e}")

    except Exception as e:
        print(f"Error processing source {source_file}: {e}")

    return written_files, found_keys

def convert_to_gff(dna, cds, gff, fann):
    """
    Generates a GFF3 file from a DNA fasta and a Prodigal Protein fasta.
    Includes a subfunction to parse InterPro annotations for protein naming.
    Default name is 'protein of unknown function' if no annotation matches.
    """

    def smart_open(path, mode='rt'):
        if str(path).endswith('.gz'): return gzip.open(path, mode)
        return open(path, mode)

    def get_annotation_map(fann_path):
        """
        Subfunction to read InterPro TSV, calculate coverage, 
        and return the best description for each protein.
        """
        if not os.path.exists(fann_path):
            return {}

        try:
            
            # Read minimal columns needed
            df = pd.read_csv(
                fann_path, 
                sep='\t', 
                usecols=["Protein_accession", "Start_location", "Stop_location", "InterPro_annotations_description"],
                dtype={ "Protein_accession": str, "Start_location": int, "Stop_location": int, "InterPro_annotations_description": str},
                on_bad_lines='skip',
                engine='c' 
            )

            # Drop rows where description is missing
            df = df.dropna(subset=['InterPro_annotations_description'])

            if df.empty:
                return {}

            # Calculate coverage (Length of the match)
            df['coverage'] = df['Stop_location'].astype(int) - df['Start_location'].astype(int)

            # Sort by coverage descending
            df = df.sort_values(by='coverage', ascending=False)

            # Keep only the single best entry (highest coverage) for each protein
            best_hits = df.drop_duplicates(subset='Protein_accession', keep='first')

            # Convert to dictionary { 'ProteinID': 'Description' }

            return dict(zip(best_hits['Protein_accession'], best_hits['InterPro_annotations_description']))

        except Exception as e:
            # Silent fail or warning, return empty dict so defaults are used
            # print(f"Warning: Failed to parse annotations from {fann_path}: {e}")
            return {}

    def parse_fasta_lengths(fasta_path):
        lengths = {}
        try:
            with smart_open(fasta_path, 'rt') as fh:
                seq_id = None
                seq_len = 0
                for line in fh:
                    line = line.rstrip()
                    if not line: continue
                    if line.startswith('>'):
                        if seq_id is not None:
                            lengths[seq_id] = seq_len
                        seq_id = line[1:].split()[0]
                        seq_len = 0
                    else:
                        seq_len += len(line)
                if seq_id is not None:
                    lengths[seq_id] = seq_len
            return lengths
        except (OSError, EOFError):
            return {}

    # --- Main Execution Logic ---

    if os.path.exists(gff):
        return
    
    # 1. Load Annotations
    annotation_map = get_annotation_map(fann)

    # 2. Get contig lengths
    contig_lengths = parse_fasta_lengths(dna)

    cds_entries = []
    
    # 3. Parse Prodigal .faa file and merge with annotations
    try:
        with smart_open(cds, 'rt') as faa_fh:
            for line in faa_fh:
                if not line.startswith('>'): continue
                
                # Header format: >ID # Start # End # Strand # Info
                header_parts = line[1:].strip().split('#')
                header_parts = [h.strip() for h in header_parts]

                fasta_id_full = header_parts[0].split()[0]
                # Infer contig name (everything before the last underscore)
                contig_name = '_'.join(fasta_id_full.split('_')[:-1])
                
                if len(header_parts) < 4: continue

                start, end, strand_flag = header_parts[1:4]
                strand = '+' if strand_flag == '1' else '-'
                
                # Handle Attributes
                raw_attrs = header_parts[-1]
                
                # Standardize ID
                attrs = re.sub(r'ID=[^;]+', f'ID=cds-{fasta_id_full}', raw_attrs)
                if 'ID=' not in attrs:
                    attrs = f'ID=cds-{fasta_id_full};' + attrs

                # Retrieve Best Annotation
                func_desc = annotation_map.get(fasta_id_full)
                
                if func_desc:
                    # Sanitize for GFF3 format
                    safe_desc = func_desc.replace(';', '%3B').replace('=', '%3D').replace('\n', ' ')
                    attrs += f";Name={safe_desc};product={safe_desc}"
                else:
                    # --- UPDATED DEFAULT NAME ---
                    attrs += ";Name=protein of unknown function;product=protein of unknown function"

                phase = 0 
                gff_line = f"{contig_name}\tprodigal\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"
                
                cds_entries.append((contig_name, int(start), gff_line))
                
    except Exception as e:
        print(f"Error processing CDS for GFF {gff}: {e}")
        return

    # 4. Sort and Write GFF
    cds_entries.sort()
    try:
        with gzip.open(gff, 'wt') as out_fh:
            out_fh.write("##gff-version 3\n")
            
            # Write Header (Sequence Regions)
            for contig, length in contig_lengths.items():
                out_fh.write(f"##sequence-region {contig} 1 {length}\n")
            
            # Write Features
            for _, _, gff_line in cds_entries:
                out_fh.write(gff_line)
                
    except Exception as e:
         print(f"Error writing GFF {gff}: {e}")

class Gatherer():

    def __init__(self, dataset: Dataset, config: Configuration):

        self.MGnify_API = hp.MGnify_API
        self.output_dir = dataset.output_dir
        self.proxies = config.arguments['MGnify_proxies']['value']
        self.failed_assemblies = dataset.failed_assemblies

        ##############################################################
        # initialize download paths

        if dataset.targets_file.exists():
            self.download_targets = pd.read_csv( dataset.targets_file )

        self.download_streams = ['assembly', 'assembly_cds', 'assembly_fannot']
        self.assemblies_dir  = dataset.assemblies_dir
        self.proteins_dir = dataset.proteins_dir
        self.annotations_dir = dataset.annotations_dir

        ##############################################################
        # initialize extraction paths
        if dataset.mgyp_metadata_file.exists():
            self.mgyp_metadata = pd.read_csv( dataset.mgyp_metadata_file )
        
        self.mgyp_metadata['ERZ_contig'] = self.mgyp_metadata['ERZ_contig'].str.replace(' NODE','-NODE')
        self.contigs_dir = dataset.contigs_dir
        self.contigs_proteins_dir = dataset.contigs_proteins_dir
        self.contigs_ann_dir = dataset.contigs_ann_dir
        self.contigs_gff_dir = dataset.contigs_gff_dir

        self.failed_file = self.output_dir / 'metadata' / "failed_extractions.txt"
        self.failed_mgycs = set()
        
        if self.failed_file.exists():
            with open(self.failed_file, 'r') as f:
                self.failed_mgycs = set(line.strip() for line in f if line.strip())

    def _schedule_extraction(self):

        self.directories = { 'contigs': self.contigs_dir, 'cds': self.contigs_proteins_dir, 'fannotation': self.contigs_ann_dir, 'gff': self.contigs_gff_dir }

        self.mgycs = self.mgyp_metadata['MGYC'].unique().tolist()
        self.mgyc_to_erz = { row['MGYC']: row['ERZ'] for _,row in self.mgyp_metadata.iterrows() }
        self.mgyc_to_erz_contig = { row['MGYC']: row['ERZ_contig'] for _,row in self.mgyp_metadata.iterrows() }

        # Setup DataFrame
        self.extraction_targets = pd.DataFrame(index=self.mgycs)
        self.extraction_targets['ERZ'] = self.extraction_targets.index.map(self.mgyc_to_erz)
        self.extraction_targets['ERZ_contig'] = self.extraction_targets.index.map(self.mgyc_to_erz_contig)        
        self.extraction_targets = self.extraction_targets[ ~self.extraction_targets['ERZ'].isin( self.failed_assemblies ) ]
        self.extraction_targets.dropna(subset=['ERZ_contig'],inplace=True)
        
        # clean and update variables
        self.mgycs = self.extraction_targets.index.tolist()
        self.erz = self.extraction_targets['ERZ'].unique().tolist()
        self.download_targets = self.download_targets[ self.download_targets['ERZ'].isin( self.erz ) ]

        file_specs = {
            'contig': {'directory': self.contigs_dir, 'format': '.fna.gz', 'column_prefix': 'contig'},
            'cds': {'directory': self.contigs_proteins_dir, 'format': '.faa.gz', 'column_prefix': 'cds'},
            'gff': {'directory': self.contigs_gff_dir, 'format': '.gff.gz', 'column_prefix': 'gff'},
            'fannot':  {'directory': self.contigs_ann_dir, 'format': '_InterPro.tsv.gz', 'column_prefix': 'fannot'}
        }

        for _, specs in file_specs.items():
            
            output_dir = specs['directory']
            file_format = specs['format']
            column_prefix = specs['column_prefix']

            file_map = {mgyc: output_dir / f'{mgyc}{file_format}' for mgyc in self.mgycs}
            self.extraction_targets[f'{column_prefix}_file'] = self.extraction_targets.index.map(file_map)
            self.extraction_targets[f'extracted_{column_prefix}_file'] = self.extraction_targets[f'{column_prefix}_file'].apply(lambda x: os.path.exists(str(x)))

        self.extraction_targets['extraction_failed'] = self.extraction_targets.index.isin(self.failed_mgycs)

    def _schedule_download(self):

        self.download_descriptions = { 'assembly': hp.contigs_description, 'assembly_cds': hp.cds_description, 'assembly_fannot': hp.fannotation_description }
        self.download_directories = { 'assembly': self.assemblies_dir, 'assembly_cds': self.proteins_dir, 'assembly_fannot': self.annotations_dir }

        stream_to_extracted_map = { 'assembly': 'contig', 'assembly_cds': 'cds', 'assembly_fannot': 'fannot' }
        
        for ds in self.download_streams:
            
            idx = self.download_targets[ self.download_targets['attributes.description.label'].isin( self.download_descriptions[ds] )].index
            
            local_file = {e:[] for e in self.erz}
            url = {e:[] for e in self.erz}

            for _, row in self.download_targets.loc[idx].iterrows():
                
                filename = self.download_directories[ds] / row['ftp.url'].split('/')[-1]
                local_file[ row['ERZ'] ].append( str(filename) )
                url[ row['ERZ'] ].append( row['ftp.url'] )

            self.extraction_targets[f'{ds}_local_file'] = self.extraction_targets['ERZ'].map(local_file)
            self.extraction_targets[f'{ds}_url'] = self.extraction_targets['ERZ'].map(url)
            #self.extraction_targets[f'downloaded_{ds}'] = self.extraction_targets[f'{ds}_local_file'].apply( lambda x: all(os.path.exists(str(f)) for f in x) if isinstance(x, list) else (os.path.exists(str(x)) if pd.notna(x) else False))
            extracted_col = f'extracted_{stream_to_extracted_map.get(ds)}_file'
            
            # Function to check download existence
            check_download = lambda x: all(os.path.exists(str(f)) for f in x) if isinstance(x, list) else (os.path.exists(str(x)) if pd.notna(x) else False)
            
            # Set to True if (Extracted File Exists) OR (Downloaded File Exists)
            # We use row-wise apply (axis=1) to access both columns
            self.extraction_targets[f'downloaded_{ds}'] = self.extraction_targets.apply(
                lambda row: True if row.get(extracted_col, False) else check_download(row[f'{ds}_local_file']), 
                axis=1
            )
        self.extraction_targets.sort_values('ERZ',inplace=True)

    def _download_batch_files(self, batch_erzs):
        """
        Internal worker: Downloads raw files. 
        Handles columns that contain LISTS of files (e.g. multiple CDS files per assembly).
        """
        stream_to_status_col = {
            'assembly': 'extracted_contig_file',
            'assembly_cds': 'extracted_cds_file',
            'assembly_fannot': 'extracted_fannot_file'
        }

        # Filter master DF for this batch
        batch_view = self.extraction_targets[self.extraction_targets['ERZ'].isin(batch_erzs)]
        
        tasks = []

        for erz in batch_erzs:
            # Get all MGYC targets belonging to this ERZ
            erz_mgycs = batch_view[batch_view['ERZ'] == erz]
            if erz_mgycs.empty: continue

            for ds in self.download_streams:
                status_col = stream_to_status_col.get(ds)
                
                # Check status: Do we need these files?
                # If we aren't tracking status (e.g. unknown stream) or if work is incomplete:
                if not status_col or not erz_mgycs[status_col].all():
                    
                    # Retrieve the lists
                    # Note: The column contains a list, so .iloc[0] returns that list.
                    url_list = erz_mgycs.iloc[0].get(f'{ds}_url')
                    path_list = erz_mgycs.iloc[0].get(f'{ds}_local_file')
                    
                    # Validation: Ensure they are actually lists and match in length
                    if isinstance(url_list, list) and isinstance(path_list, list):
                        for url, local_path in zip(url_list, path_list):
                            if pd.notna(url) and pd.notna(local_path) and not os.path.exists(local_path):
                                tasks.append((url, local_path))
                                
                    # Fallback for legacy string format (just in case)
                    elif isinstance(url_list, str) and isinstance(path_list, str):
                        if not os.path.exists(path_list):
                            tasks.append((url_list, path_list))

        if not tasks:
            return

        # Define task locally or use self-contained helper
        def download_task(url, path):
            try:
                response = requests.get(url, stream=True, proxies=self.proxies)
                if response.status_code == 200:
                    with open(path, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)
                    return True
            except Exception as e:
                print(f"Failed {url}: {e}")
            return False

        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(download_task, t[0], t[1]) for t in tasks]
            concurrent.futures.wait(futures)

    def _extract_batch_sequences(self, batch_erzs):
        """
        Extracts sequences. 
        Updates 'extracted' status on success.
        Updates 'failed' status on failure and writes to disk.
        """
        
        batch_df = self.extraction_targets[self.extraction_targets['ERZ'].isin(batch_erzs)]
        
        jobs = []
        # We need to track which MGYCs we are attempting in this batch to check for failures later
        # Structure: job_future -> { 'type': 'contig'/'cds', 'erz': erz_name, 'targets': set(mgyc_ids) }
        future_metadata = {} 
        
        # Also track "Batch Level" expected targets per ERZ to calculate failures after all files are processed
        # { 'ERZ_NAME': { 'contig': set(expected_mgycs), 'cds': set(expected_mgycs) } }
        erz_expectations = {}

        with ProcessPoolExecutor(max_workers=8) as executor:

            for erz, group in batch_df.groupby('ERZ'):
                
                erz_expectations[erz] = {'contig': set(), 'cds': set()}

                # --- 1. Contigs (DNA) ---
                # Filter: Not extracted AND Not previously failed
                pending_contigs = group[
                    (~group['extracted_contig_file']) & 
                    (~group['extraction_failed'])
                ]
                
                if not pending_contigs.empty:
                    sources = group.iloc[0]['assembly_local_file']
                    if isinstance(sources, str): sources = [sources]
                    
                    if isinstance(sources, list):
                        target_map = dict(zip(pending_contigs['ERZ_contig'], pending_contigs['contig_file']))
                        # We expect to find these
                        erz_expectations[erz]['contig'] = set(pending_contigs['ERZ_contig'])
                        
                        for src in sources:
                            if pd.notna(src) and os.path.exists(src):
                                future = executor.submit(extract_sequences_worker, src, target_map)
                                jobs.append(future)
                                future_metadata[future] = {'type': 'contig', 'erz': erz}

                # --- 2. CDS (Proteins) ---
                pending_cds = group[
                    (~group['extracted_cds_file']) & 
                    (~group['extraction_failed'])
                ]

                if not pending_cds.empty:
                    sources = group.iloc[0]['assembly_cds_local_file']
                    if isinstance(sources, str): sources = [sources]
                    
                    if isinstance(sources, list):
                        target_map = dict(zip(pending_cds['ERZ_contig'], pending_cds['cds_file']))
                        erz_expectations[erz]['cds'] = set(pending_cds['ERZ_contig'])
                        
                        for src in sources:
                            if pd.notna(src) and os.path.exists(src):
                                future = executor.submit(extract_sequences_worker, src, target_map)
                                jobs.append(future)
                                future_metadata[future] = {'type': 'cds', 'erz': erz}

            # --- Process Results ---
            # We must aggregate found keys per ERZ+Type because one ERZ might have multiple source files
            # { 'ERZ_NAME': { 'contig': set(found_mgycs), 'cds': set(found_mgycs) } }
            found_tracker = {erz: {'contig': set(), 'cds': set()} for erz in erz_expectations}

            for future in as_completed(jobs):
                try:
                    meta = future_metadata[future]
                    created_files, found_keys = future.result()
                    
                    # 1. Update Success State (Immediate)
                    for filepath in created_files:
                        filepath_str = str(filepath)
                        #print(filepath_str)
                        if 'contig' in filepath_str:
                            mask = self.extraction_targets['contig_file'].astype(str) == filepath_str
                            self.extraction_targets.loc[mask, 'extracted_contig_file'] = True
                        if 'proteins' in filepath_str:
                            mask = self.extraction_targets['cds_file'].astype(str) == filepath_str
                            self.extraction_targets.loc[mask, 'extracted_cds_file'] = True
                    
                    # 2. Accumulate found keys for failure checking
                    found_tracker[meta['erz']][meta['type']].update(found_keys)
                    
                except Exception as e:
                    print(f"Extraction worker failed: {e}")
            #print(self.extraction_targets)

            # --- Calculate Failures ---
            new_failures = set()
            
            for erz, expectations in erz_expectations.items():
                # Check Contigs
                expected_contigs = expectations['contig']
                found_contigs = found_tracker[erz]['contig']
                missing_contigs = expected_contigs - found_contigs
                
                # Check CDS
                expected_cds = expectations['cds']
                found_cds = found_tracker[erz]['cds']
                missing_cds = expected_cds - found_cds
                
                # If a target was expected but missing in ALL files for this ERZ, it's a failure
                # We need to map back from 'ERZ_contig' (e.g. substring) to 'MGYC' (index)
                if missing_contigs:
                    # Get MGYCs corresponding to these missing substrings
                    mask = (self.extraction_targets['ERZ'] == erz) & (self.extraction_targets['ERZ_contig'].isin(missing_contigs))
                    new_failures.update(self.extraction_targets[mask].index.tolist())

                if missing_cds:
                    mask = (self.extraction_targets['ERZ'] == erz) & (self.extraction_targets['ERZ_contig'].isin(missing_cds))
                    new_failures.update(self.extraction_targets[mask].index.tolist())

            # --- Update Disk and Memory ---
            if new_failures:
                print(f"  Warning: {len(new_failures)} sequences could not be found. Marking as failed.")
                
                # Update DataFrame
                self.extraction_targets.loc[list(new_failures), 'extraction_failed'] = True
                self.failed_mgycs.update(new_failures)
                
                # Append to file
                with open(self.failed_file, 'a') as f:
                    for mgyc in new_failures:
                        f.write(f"{mgyc}\n")

    def _extract_batch_fannotation(self, batch_erzs):
        """Internal worker: Extracts functional annotations for the batch."""
        return
    
    def _generate_gff_batch(self, batch_erzs, max_workers=8):
        """
        Creates GFF files for the specific batch of ERZs.
        """
        # 1. Filter targets for this batch
        batch_view = self.extraction_targets[self.extraction_targets['ERZ'].isin(batch_erzs)]

        # 2. Identify rows that need GFF creation
        # Condition: GFF missing AND Protein present
        todo_mask = (
            (~batch_view['extracted_gff_file']) & 
            (batch_view['extracted_contig_file']) & 
            (batch_view['extracted_cds_file'])
        )
        
        todo_rows = batch_view[todo_mask]

        if todo_rows.empty:
            print(f'todo_rows is empty for batch {batch_erzs}')
            return

        # 3. Prepare jobs
        jobs = []
        for mgyc, row in todo_rows.iterrows():
            dna = str(row['contig_file'])
            proteins = str(row['cds_file'])
            gff = str(row['gff_file'])
            fannot = str(row['fannot_file'])

            if os.path.exists(dna) and os.path.exists(proteins) and os.path.exists(fannot):
                jobs.append((mgyc, dna, proteins, gff, fannot))

        if not jobs:
            print('no jobs to do')
            return

        # 4. Execute GFF Conversion
        # We rely on the global 'convert_to_gff' function defined outside the class
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_mgyc = {
                executor.submit(convert_to_gff, dna, prot, gff, fannot): mgyc 
                for mgyc, dna, prot, gff, fannot in jobs
            }
            
            for future in as_completed(future_to_mgyc):
                mgyc = future_to_mgyc[future]
                try:
                    future.result()
                    # Update status in the main dataframe
                    self.extraction_targets.loc[mgyc, 'extracted_gff_file'] = True
                except Exception as e:
                    print(f"GFF creation failed for {mgyc}: {e}")

    def _read_interpro(self, interpro_file):
        """
        Helper to safely read InterProScan TSV files. 
        Robustly handles files with or without headers to prevent type errors.
        """
        interpro_columns = [
            "Protein_accession", "Sequence_MD5", "Sequence_length", "Analysis",
            "Signature_accession", "Signature_description", "Start_location",
            "Stop_location", "Score", "Status", "Date",
            "InterPro_annotations_accession", "InterPro_annotations_description",
            "GO_annotations", "Pathways_annotations"
        ]

        try:
            # 1. Read all data as strings first (dtype=str) to prevent immediate int conversion crashes
            df = pd.read_csv(
                interpro_file, 
                sep='\t', 
                header=None, 
                names=interpro_columns, 
                dtype=str,               # Read everything as string initially
                on_bad_lines='skip',
                engine='c'
            )
            
            if df.empty:
                return df

            # 2. Robust Header Removal
            # If the file has a header, the 'Start_location' column will contain the word "Start_location" (or "start_location")
            # We filter for rows where Start_location is actually a digit.
            df = df[df['Start_location'].str.isnumeric()]

            # 3. Safe Type Conversion
            # Now that header rows are gone, we can safely convert columns
            df['Start_location'] = pd.to_numeric(df['Start_location'])
            df['Stop_location'] = pd.to_numeric(df['Stop_location'])
            df['Sequence_length'] = pd.to_numeric(df['Sequence_length'])
            
            # Score is left as object/str because it can be mixed (e.g. floats, scientific notation, or "-")
            
            return df

        except Exception as e:
            print(f"Error reading InterPro file {interpro_file}: {e}")
            return pd.DataFrame(columns=interpro_columns)
    
    def _extract_annotations_batch(self, batch_erzs):
        """
        Extracts InterPro annotations for the batch.
        1. Reads the large Assembly Annotation file (once per ERZ).
        2. Gets Protein IDs from the *already extracted* CDS files.
        3. Filters the large dataframe and saves the subset.
        """
        # Filter for this batch
        batch_view = self.extraction_targets[self.extraction_targets['ERZ'].isin(batch_erzs)]
        # We must group by ERZ to avoid reading the large source file multiple times
        for erz, group in batch_view.groupby('ERZ'):
            
            # Check if we have work to do for this ERZ (any missing annotations?)
            if group['extracted_fannot_file'].all():

                continue

            # Check if source file exists
            source_file = group.iloc[0].get('assembly_fannot_local_file')[0]

            if pd.isna(source_file) or not os.path.exists(str(source_file)):
                print(f"Missing source annotation file ({source_file}) for {erz}, skipping.")
                continue
            
            # --- 1. Load the Big InterPro File ---
            # print(f"  Reading InterPro for {erz}...")
            interpro_df = self._read_interpro(source_file)
            
            if interpro_df.empty:
                continue

            # --- 2. Iterate over MGYCs (Targets) ---
            for mgyc, row in group.iterrows():
                
                # Skip if already done
                if row['extracted_fannot_file']:
                    continue

                extracted_cds_path = row['cds_file']
                output_path = row['fannot_file']

                if not os.path.exists(extracted_cds_path):
                    # We can't filter without protein IDs, so we skip
                    continue

                try:
                    # Get Protein IDs from the extracted CDS file
                    ids = []
                    with gzip.open(extracted_cds_path, "rt") as handle:
                        # Use simple parsing for speed
                        for line in handle:
                            if line.startswith('>'):
                                # >ID description -> ID
                                ids.append(line[1:].split()[0])

                    if not ids:
                        # Create empty file if no proteins found (to mark as done)
                        pd.DataFrame(columns=interpro_df.columns).to_csv(output_path, sep='\t', index=False)
                        self.extraction_targets.loc[mgyc, 'extracted_fannot_file'] = True
                        continue

                    # Filter the big dataframe
                    # subset = interpro_df[interpro_df['Protein_accession'].isin(ids)]
                    # Optimization: Set index for faster lookup if list is large, 
                    # but 'isin' is usually optimized enough.
                    subset = interpro_df.loc[interpro_df['Protein_accession'].isin(ids)]
                    
                    # Save
                    subset.to_csv(output_path, sep='\t', index=False)
                    self.extraction_targets.loc[mgyc, 'extracted_fannot_file'] = True

                except Exception as e:
                    print(f"Error extracting annotations for {mgyc}: {e}")

    def _cleanup_batch_files(self, batch_erzs):
        """
        Deletes the raw downloaded files for the batch to free up disk space.
        Handles LISTS of files (e.g., multiple CDS files per assembly).
        """
        batch_view = self.extraction_targets[self.extraction_targets['ERZ'].isin(batch_erzs)]
        unique_erzs = batch_view.drop_duplicates(subset=['ERZ'])

        for _, row in unique_erzs.iterrows():
            for ds in self.download_streams:
                # Get the value (could be string, list, or NaN)
                local_path_data = row.get(f'{ds}_local_file')
                
                # Normalize to list
                if isinstance(local_path_data, str):
                    file_list = [local_path_data]
                elif isinstance(local_path_data, list):
                    file_list = local_path_data
                else:
                    continue # Skip NaN or other types

                # Iterate and delete
                for local_path in file_list:
                    # Safety check: ensure path is valid string and file exists
                    if pd.notna(local_path) and os.path.exists(local_path):
                        try:
                            os.remove(local_path)
                            #print(f"Deleted {local_path}")
                        except OSError as e:
                            print(f"Warning: Could not delete {local_path}: {e}")

    def run_pipeline(self, batch_size=5):
        """
        Main Pipeline Entry Point with Progress Tracking.
        """

        print("Initializing Pipeline Schedules...")
        self._schedule_extraction()
        self._schedule_download()

        # --- 1. Status Report ---
        total_mgycs = len(self.extraction_targets)
        
        # Count successes
        completed_contigs = self.extraction_targets['extracted_contig_file'].sum()
        completed_cds = self.extraction_targets['extracted_cds_file'].sum()
        completed_gff = self.extraction_targets['extracted_gff_file'].sum()
        completed_fann = self.extraction_targets['extracted_fannot_file'].sum()
        
        # Count failures
        # (Assuming 'extraction_failed' column exists from schedule_extraction)
        failed_count = self.extraction_targets['extraction_failed'].sum()

        # Determine which ERZs actually need work
        # Logic: An ERZ is pending if ANY of its MGYCs are missing ANY file AND are NOT failed.
        # We exclude failed rows entirely from the "To Process" list.
        pending_mask = (
            (
                (~self.extraction_targets['extracted_contig_file']) | 
                (~self.extraction_targets['extracted_cds_file']) |
                (~self.extraction_targets['extracted_gff_file']) |
                (~self.extraction_targets['extracted_fannot_file'])
            ) & 
            (~self.extraction_targets['extraction_failed'])
        )
        
        pending_erzs = self.extraction_targets[pending_mask]['ERZ'].unique().tolist()
        pending_erzs.sort()
        total_pending_erzs = len(pending_erzs)
        total_erzs = len(self.erz)

        print("\n" + "="*45)
        print(f"      PIPELINE STATUS REPORT")
        print("="*45)
        print(f"Total Targets (MGYCs) : {total_mgycs}")
        print(f"  - Failed / Skipped  : {failed_count}")
        print("-" * 45)
        print(f"Completed Contigs     : {completed_contigs}")
        print(f"Completed CDS         : {completed_cds}")
        print(f"Completed GFF         : {completed_gff}")
        print(f"Completed Annotations : {completed_fann}")
        print("-" * 45)
        print(f"Assemblies to Process : {total_pending_erzs}/{total_erzs}")
        print(f"Batch Size            : {batch_size}")
        print("="*45 + "\n")

        if total_pending_erzs == 0:
            print("All targets are already extracted. Pipeline complete.")
            return

        # --- 2. Batch Execution ---
        erz_batches = [
            pending_erzs[i : i + batch_size] 
            for i in range(0, total_pending_erzs, batch_size)
        ]
        
        self.current_chunk = []
        with tqdm(total=len(erz_batches), desc="Pipeline Progress", unit="batch") as pbar:
            
            for chunk in erz_batches:

                # 1. Download (Only if files missing)
                pbar.set_postfix_str("Downloading")
                #print("Downloading")
                self._download_batch_files(chunk)
                
                # 2. Extract (Only if files present)
                pbar.set_postfix_str("Extracting protein sequences")
                #print("Extracting sequences")
                self._extract_batch_sequences(chunk)

                pbar.set_postfix_str("Extracting Annotations")
                #print("Extracting Annotations")
                self._extract_annotations_batch(chunk)
            
                # 3. Cleanup (Always delete raw files to save space)
                pbar.set_postfix_str("Generating genomic features")
                #print("Generating genomic features")
                self._generate_gff_batch(chunk)
    
                # 3. Cleanup (Always delete raw files to save space)
                pbar.set_postfix_str("Cleaning")
                #print("Cleaning")
                self._cleanup_batch_files(chunk)
                
                pbar.update(1)
            #return


        print("\nPipeline Execution Finished.")