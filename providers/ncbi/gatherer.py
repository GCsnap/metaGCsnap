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
import requests
import pandas as pd
from Bio import SeqIO

from jsonapi_client import Session as APISession, Modifier

import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

from gcsnap.configuration import Configuration
from providers.MGnify.dataset import Dataset
import providers.ncbi.helpers as hp
from gcsnap.utils import handle_compressed_fasta


def get_genomic_region(gff_path, protein_id):
    """
    Parses a GFF3 file (compressed or uncompressed) to find the chromosome/contig 
    where a specific protein is encoded.

    Args:
        gff_path (str): Path to the .gff or .gff.gz file.
        protein_id (str): The ID of the protein to search for.

    Returns:
        str: The name of the chromosome/contig (column 1 in GFF), or None if not found.
    """
    
    # Internal helper to handle gzip automatically
    def smart_open(path, mode='rt'):
        if str(path).endswith('.gz'): 
            # 'rt' mode is critical here to read as text (strings) instead of bytes
            return gzip.open(path, mode)
        return open(path, mode)

    # 1. Prepare Regex Patterns
    # We escape the protein_id to ensure special characters (like dots in 'NP_123.1') don't break regex
    safe_id = re.escape(protein_id)
    
    # Matches exact ID=... or Name=... or ID=cds-... (common in Prodigal)
    # The (;|$) ensures we match "ID=Prot1" but NOT "ID=Prot12"
    patterns = [
        re.compile(f"ID={safe_id}(;|$)"),
        re.compile(f"Name={safe_id}(;|$)"),
        re.compile(f"ID=cds-{safe_id}(;|$)"),
        re.compile(f"Parent={safe_id}(;|$)") # Sometimes helpful if looking for CDS linking to Gene
    ]

    found_region = None

    try:
        # 2. Open File (Handles .gz automatically)
        with smart_open(gff_path) as fh:
            
            for line in fh:
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.split('\t')
                
                # Standard GFF3 has 9 columns
                if len(parts) < 9:
                    continue

                # Attributes are in the last column (index 8)
                attributes = parts[8]

                # 3. Check for Match
                # We check all patterns. If any match, we return the SeqID (col 0)
                for pattern in patterns:
                    if pattern.search(attributes):
                        found_region = parts[0]
                        return found_region

    except FileNotFoundError:
        print(f"Error: File not found at {gff_path}")
    except Exception as e:
        print(f"Error reading GFF {gff_path}: {e}")

    return found_region

def extract_cds_from_region(gff_path: str, contig: str) -> tuple[list[str], str]:
    
    """
    Extract CDS protein IDs and sliced GFF content for a given contig.
    
    Args:
        gff_path: Path to .gff or .gff.gz file
        contig: Contig/sequence name (e.g. 'MLJW01000001.1')
    
    Returns:
        - List of protein IDs (e.g. ['OIR19122.1', 'OIR19123.1'])
        - Sliced GFF string (header lines + all features on the contig)
    """
    opener = gzip.open if gff_path.endswith('.gz') else open
    cds_ids = []
    header_lines = []
    feature_lines = []
    
    with opener(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            seq, _, feature, *_, attributes = fields
            if seq != contig:
                continue

            feature_lines.append(line)

            if feature == 'CDS':
                match = re.search(r'protein_id=([^;]+)', attributes)
                if match:
                    cds_ids.append(match.group(1))
    
    sliced_gff = ''.join(header_lines) + ''.join(feature_lines)
    
    return cds_ids, sliced_gff


class Gatherer():

    def __init__(self, dataset: Dataset, config: Configuration):

        #self.MGnify_API = hp.MGnify_API
        self.output_dir = dataset.output_dir
        self.proxies = config.arguments['MGnify_proxies']['value']
        #self.failed_assemblies = dataset.failed_assemblies

        ##############################################################
        # initialize download paths

        if dataset.targets_file.exists():
            self.download_targets = pd.read_csv( dataset.targets_file )

        self.download_streams = ['assembly', 'proteins', 'features']
        self.assemblies_dir  = dataset.assemblies_dir
        self.proteins_dir = dataset.proteins_dir
        self.features_dir = dataset.gff_dir

        self.contigs_dir = dataset.contigs_dir
        self.contigs_proteins_dir = dataset.contigs_proteins_dir
        self.contigs_gff_dir = dataset.contigs_gff_dir

        self.targets = pd.read_csv(dataset.targets_file)
        self.assemblies = self.targets['assembly_accession'].unique()

        self.failed_file = self.output_dir / 'metadata' / "failed_extractions.txt"
        self.failed_assemblies = set()
        
        if self.failed_file.exists():
            with open(self.failed_file, 'r') as f:
                self.failed_assemblies = set(line.strip() for line in f if line.strip())

        self.ncbi_metadata_file = dataset.ncbi_metadata_file
        self.assembly_metadata = dataset.assembly_metadata
        self.mappings_file = dataset.mappings_file
        
    def _schedule_download(self):
        
        # 1. Setup Configuration
        # Map the stream name -> (Target Directory, Description List in input DF)
        stream_config = {
            'assembly':      {'dir': self.assemblies_dir, 'desc': hp.assembly_description},
            'proteins':  {'dir': self.proteins_dir,   'desc': hp.cds_description},
            'features':  {'dir': self.features_dir,   'desc': hp.features_description}
        }

        # 2. Initialize the Assembly-Centric DataFrame
        # Get unique assemblies from the input targets
        
        self.download_targets = pd.DataFrame(index=self.assemblies)

        # 3. Populate Columns for each Stream
        for ds, config in stream_config.items():
            
            # Filter the raw input to get rows for this specific file type
            subset = self.targets[self.targets['file_type'].isin(config['desc'])].copy()
            
            # Define column names
            url_col = f'{ds}_url'
            file_col = f'{ds}_file'
            status_col = f'{ds}_downloaded'

            if subset.empty:
                self.download_targets[url_col] = None
                self.download_targets[file_col] = None
                self.download_targets[status_col] = False
                continue

            # Create mappings from Assembly Accession -> Data
            # 1. FTP URL Mapping
            url_map = subset.set_index('assembly_accession')['ftp_url'].to_dict()
            self.download_targets[url_col] = self.download_targets.index.map(url_map)

            # 2. Local File Path Mapping
            target_dir = config['dir']
            # We use the URL to derive the filename
            path_map = subset.set_index('assembly_accession')['ftp_url'].apply(
                lambda url: str(target_dir / url.split('/')[-1]) if pd.notnull(url) else None
            ).to_dict()
            
            self.download_targets[file_col] = self.download_targets.index.map(path_map)

            # 3. Check Existence (Vectorized)
            self.download_targets[status_col] = self.download_targets[file_col].apply(
                lambda x: os.path.exists(x) if x and pd.notnull(x) else False
            )

        # 4. Global status
        download_cols = [f'{ds}_downloaded' for ds in stream_config.keys()]
        self.download_targets['all_downloaded'] = self.download_targets[download_cols].all(axis=1)
        #self.download_targets['assembly_accession'] = self.download_targets.index
        self.download_targets = self.download_targets.reset_index()
        self.download_targets.rename(columns={'index':'assembly_accession'},inplace=True)

    def _schedule_extraction(self):
        """
        Builds self.extraction_targets: one row per unique genomic_region (contig),
        listing the files that need to be produced:
          - contig_file      : extracted .fna.gz for the contig
          - cds_file         : extracted proteins .faa.gz for the contig
          - gff_file         : sliced .gff.gz for the contig

        The index is genomic_region (== contig id, e.g. MLJW01000001.1).
        Each row also carries the assembly_accession it belongs to, so we know
        which downloaded files to open.
        """

        # ncbi_metadata has one row per ncbi_code (protein), with columns:
        #   ncbi_code, genomic_region, assembly_accession, ...
        meta = self.ncbi_metadata.copy()

        # Drop rows with no genomic region resolved
        meta = meta.dropna(subset=['genomic_region'])

        # Exclude assemblies that previously failed
        meta = meta[~meta['assembly_accession'].isin(self.failed_assemblies)]

        # One row per unique contig — keep assembly_accession alongside
        contigs = (
            meta[['genomic_region', 'assembly_accession']]
            .drop_duplicates(subset='genomic_region')
            .set_index('genomic_region')
        )

        self.extraction_targets = contigs.copy()

        # Map contig -> output file paths
        file_specs = {
            'contig': (self.contigs_dir,          '.fna.gz'),
            'cds':    (self.contigs_proteins_dir,  '.faa.gz'),
            'gff':    (self.contigs_gff_dir,        '.gff.gz'),
        }

        for col_prefix, (out_dir, ext) in file_specs.items():
            file_map = {
                contig: out_dir / f'{contig}{ext}'
                for contig in self.extraction_targets.index
            }
            self.extraction_targets[f'{col_prefix}_file'] = self.extraction_targets.index.map(file_map)
            self.extraction_targets[f'extracted_{col_prefix}_file'] = (
                self.extraction_targets[f'{col_prefix}_file']
                .apply(lambda x: os.path.exists(str(x)))
            )

        self.extraction_targets['extraction_failed'] = (
            self.extraction_targets.index.isin(self.failed_assemblies)
        )

        # Convenience: unique assemblies still needed
        self.unique_contigs   = self.extraction_targets.index.tolist()
        self.unique_assemblies = self.extraction_targets['assembly_accession'].unique().tolist()

    def _download_batch_files(self, batch_asms):
        """
        Internal worker: Downloads raw files for a batch of assemblies.
        Coherent with the assembly-centric download_targets DataFrame.
        """

        # 1. Filter the Master DataFrame for this batch
        # We only want rows corresponding to the assemblies in the current batch
        batch_df = self.download_targets[
            self.download_targets['assembly_accession'].isin(batch_asms)
        ]

        if batch_df.empty:
            return

        tasks = []

        # 2. Collect Tasks for each Stream
        for ds in self.download_streams: # e.g. ['assembly', 'proteins', 'features']
            
            url_col = f'{ds}_url'
            file_col = f'{ds}_file'
            status_col = f'{ds}_downloaded'

            # Safety check: Ensure columns exist (in case a stream was skipped/empty)
            if url_col not in batch_df.columns:
                continue

            # 3. Vectorized Filtering
            # We want rows where:
            # a) The URL is not Null (file exists for this assembly)
            # b) The file is NOT marked as downloaded/extracted yet
            target_rows = batch_df[
                pd.notna(batch_df[url_col]) & 
                (batch_df[status_col] == False)
            ]

            # Extract pairs of (URL, Local Path)
            # zip() allows us to iterate cleanly over the two series
            for url, local_path in zip(target_rows[url_col], target_rows[file_col]):
                if pd.notna(local_path):
                     tasks.append((url, local_path))

        if not tasks:
            return

        # 4. Define the Download Worker
        def download_task(url, path):
            try:
                # Ensure directory exists just in case
                os.makedirs(os.path.dirname(path), exist_ok=True)
                
                # Streaming download
                with requests.get(url, stream=True) as r:
                    r.raise_for_status() # Check for HTTP errors (404, etc.)
                    with open(path, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk: 
                                f.write(chunk)
                return True
            except Exception as e:
                print(f"Failed downloading {url} to {path}: {e}")
                # Clean up partial file if it failed
                if os.path.exists(path):
                    os.remove(path)
                return False

        # 5. Execute in Parallel
        # Adjust max_workers based on your connection/CPU
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(download_task, url, path) for url, path in tasks]
            concurrent.futures.wait(futures)

    def _get_ncbi_metadata(self):

        # 1. Prepare Assembly Metadata (Explode NCBI codes)
        # Extracts the list of proteins per assembly and flattens it
        exploded_df = self.assembly_metadata[['#assembly_accession', 'ncbi_code', 'ftp_path','taxid','organism_name']].explode('ncbi_code')
        exploded_df = exploded_df.reset_index(drop=True)

        # 2. Filter Mappings
        # Only keep mappings that are relevant
        mappings = pd.read_csv(self.mappings_file)
        relevant_mask = mappings['ncbi_code'].isin(exploded_df['ncbi_code'])
        mappings = mappings[relevant_mask]

        # 3. Prepare Feature Files (GFF paths)
        features_files = self.download_targets[['assembly_accession', 'features_file']]

        # 4. Merge Everything
        # Merge mappings with assembly info
        metadata = mappings.merge(exploded_df, on='ncbi_code', how='inner')
        
        # Merge with feature file paths to get the GFF location
        metadata = metadata.merge(
            features_files, 
            left_on='#assembly_accession', 
            right_on='assembly_accession',
            how='left'
        )

        # 5. Extract Genomic Region (Chromosome/Contig)
        # This applies the external 'get_genomic_region' function row-by-row
        print("Extracting genomic regions from GFF files...")
        metadata['genomic_region'] = metadata.apply(
            lambda row: get_genomic_region(row['features_file'], row['ncbi_code']), 
            axis=1
        )

        # 6. Identify Source Database
        # Determine which column (RefSeq, UniProt, etc.) matches the 'target' ID
        source_cols = [
            'UniProtKB-AC', 'UniProtKB-ID', 'GeneID', 
            'RefSeq', 'UniParc', 'EMBL-CDS', 'Ensembl'
        ]
        
        # Create a boolean mask of matches
        matches = metadata[source_cols].eq(metadata['target'], axis=0)
        
        # Assign the column name where the match occurred
        metadata['source'] = matches.idxmax(axis=1)
        
        # Handle cases where no match was found (set to 'Unknown')
        metadata.loc[~matches.any(axis=1), 'source'] = 'Unknown'

        # 7. Cleanup
        # Remove technical columns that are no longer needed
        metadata.drop(columns=['features_file', '#assembly_accession'], inplace=True, errors='ignore')

        self.ncbi_metadata = metadata
        self.ncbi_metadata.to_csv(self.ncbi_metadata_file,index=False)

    def _extract_batch_sequences(self, batch_assemblies):
        """
        For each assembly in the batch:
          1. Find all unique contigs that belong to it (from extraction_targets).
          2. For each contig, extract the contig sequence from the assembly .fna.gz
             and write it to contigs_dir/<contig>.fna.gz.
          3. Run extract_cds_from_region on the assembly .gff.gz to get:
             - sliced GFF  -> saved to contigs_gff_dir/<contig>.gff.gz
             - protein ids -> used to fish sequences from assembly .faa.gz
                           -> saved to contigs_proteins_dir/<contig>.faa.gz
        """

        # Join extraction_targets with download_targets to get local file paths
        dl = self.download_targets.set_index('assembly_accession')

        for asm in batch_assemblies:

            if asm not in dl.index:
                print(f"  [warn] No download record for {asm}, skipping.")
                continue

            dl_row       = dl.loc[asm]
            assembly_fna = dl_row.get('assembly_file')
            assembly_faa = dl_row.get('proteins_file')
            assembly_gff = dl_row.get('features_file')

            # Sanity: all three source files must exist
            missing = [p for p in [assembly_fna, assembly_faa, assembly_gff]
                       if not p or not os.path.exists(str(p))]
            if missing:
                print(f"  [warn] Missing source files for {asm}: {missing}")
                continue

            # All contigs that belong to this assembly and still need work
            asm_contigs = self.extraction_targets[
                (self.extraction_targets['assembly_accession'] == asm) &
                (
                    (~self.extraction_targets['extracted_contig_file']) |
                    (~self.extraction_targets['extracted_cds_file'])    |
                    (~self.extraction_targets['extracted_gff_file'])
                ) &
                (~self.extraction_targets['extraction_failed'])
            ]

            if asm_contigs.empty:
                continue

            # ----------------------------------------------------------------
            # Step 1: Extract contig sequences from assembly .fna.gz
            # ----------------------------------------------------------------
            contigs_needed = asm_contigs[~asm_contigs['extracted_contig_file']].index.tolist()

            if contigs_needed:
                try:
                    opener = gzip.open if str(assembly_fna).endswith('.gz') else open
                    with opener(assembly_fna, 'rt') as fh:
                        for record in SeqIO.parse(fh, 'fasta'):
                            if record.id in contigs_needed:
                                out_path = self.extraction_targets.loc[record.id, 'contig_file']
                                os.makedirs(os.path.dirname(str(out_path)), exist_ok=True)
                                with gzip.open(str(out_path), 'wt') as out_fh:
                                    SeqIO.write(record, out_fh, 'fasta')
                                self.extraction_targets.loc[record.id, 'extracted_contig_file'] = True
                except Exception as e:
                    print(f"  [error] Extracting contigs from {assembly_fna}: {e}")

            # ----------------------------------------------------------------
            # Step 2: For each contig, slice GFF and collect protein IDs
            # ----------------------------------------------------------------
            contig_to_protein_ids = {}

            for contig, row in asm_contigs.iterrows():

                if not row['extracted_gff_file']:
                    try:
                        protein_ids, sliced_gff = extract_cds_from_region(str(assembly_gff), contig)
                        gff_out = row['gff_file']
                        os.makedirs(os.path.dirname(str(gff_out)), exist_ok=True)
                        with gzip.open(str(gff_out), 'wt') as gff_fh:
                            gff_fh.write(sliced_gff)
                        self.extraction_targets.loc[contig, 'extracted_gff_file'] = True
                        contig_to_protein_ids[contig] = set(protein_ids)
                    except Exception as e:
                        print(f"  [error] Slicing GFF for {contig}: {e}")
                else:
                    # GFF already extracted — we still need protein ids for CDS extraction below
                    # Re-derive from already-saved sliced GFF
                    if not row['extracted_cds_file']:
                        try:
                            ids, _ = extract_cds_from_region(str(assembly_gff), contig)
                            contig_to_protein_ids[contig] = set(ids)
                        except Exception as e:
                            print(f"  [error] Re-deriving protein IDs for {contig}: {e}")

            # ----------------------------------------------------------------
            # Step 3: Extract proteins from assembly .faa.gz
            # ----------------------------------------------------------------
            contigs_need_cds = asm_contigs[~asm_contigs['extracted_cds_file']].index.tolist()

            if contigs_need_cds and contig_to_protein_ids:

                # Build a flat map: protein_id -> (contig, output_faa_path)
                # We'll buffer records per contig output file
                buffers = {
                    contig: []
                    for contig in contigs_need_cds
                    if contig in contig_to_protein_ids
                }
                # Invert for fast lookup: protein_id -> contig
                protein_to_contig = {
                    pid: contig
                    for contig, pids in contig_to_protein_ids.items()
                    if contig in buffers
                    for pid in pids
                }

                try:
                    opener = gzip.open if str(assembly_faa).endswith('.gz') else open
                    with opener(assembly_faa, 'rt') as fh:
                        for record in SeqIO.parse(fh, 'fasta'):
                            contig = protein_to_contig.get(record.id)
                            if contig:
                                buffers[contig].append(record)
                except Exception as e:
                    print(f"  [error] Reading proteins from {assembly_faa}: {e}")

                # Write buffered records per contig
                for contig, records in buffers.items():
                    cds_out = self.extraction_targets.loc[contig, 'cds_file']
                    try:
                        os.makedirs(os.path.dirname(str(cds_out)), exist_ok=True)
                        with gzip.open(str(cds_out), 'wt') as out_fh:
                            SeqIO.write(records, out_fh, 'fasta')
                        self.extraction_targets.loc[contig, 'extracted_cds_file'] = True
                    except Exception as e:
                        print(f"  [error] Writing CDS for {contig}: {e}")

    def _cleanup_batch_files(self, batch_assemblies):
        """
        Deletes raw downloaded files (.fna, .faa, .gff) for each assembly in the
        batch, but ONLY once all contig-level outputs that depend on that assembly
        are fully extracted.

        Safety rule: if any contig is still missing an output, the raw source is
        kept intact so a re-run can recover without re-downloading.
        """
        dl = self.download_targets.set_index('assembly_accession')

        for asm in batch_assemblies:

            if asm not in dl.index:
                continue

            # Guard: only clean up when every contig for this assembly is done
            asm_contigs = self.extraction_targets[
                self.extraction_targets['assembly_accession'] == asm
            ]

            all_done = (
                asm_contigs['extracted_contig_file'].all() and
                asm_contigs['extracted_cds_file'].all()    and
                asm_contigs['extracted_gff_file'].all()
            )

            if not all_done:
                print(f"  [skip cleanup] {asm}: not all contigs fully extracted yet.")
                continue

            # Delete each raw stream file
            dl_row = dl.loc[asm]
            for ds in self.download_streams:
                raw_path = dl_row.get(f'{ds}_file')
                if pd.isna(raw_path) or not raw_path:
                    continue
                raw_path = str(raw_path)
                if os.path.exists(raw_path):
                    try:
                        os.remove(raw_path)
                    except OSError as e:
                        print(f"  [warn] Could not delete {raw_path}: {e}")

    def run_pipeline(self, batch_size=5):
        """
        Main Pipeline Entry Point with Progress Tracking.

        Stages:
          1. Schedule downloads  (_schedule_download)
          2. Download raw files  (_download_batch_files)
          3. Build NCBI metadata (_get_ncbi_metadata)  — resolves genomic_region per protein
          4. Schedule extraction (_schedule_extraction) — one row per unique contig
          5. Extract             (_extract_batch_sequences) — contig .fna, sliced .gff, contig .faa
        """

        print("Initializing Download Schedule...")
        self._schedule_download()

        # ----------------------------------------------------------------
        # Phase 1: Download
        # ----------------------------------------------------------------
        pending_dl_mask = (
            (~self.download_targets['assembly_downloaded']) |
            (~self.download_targets['proteins_downloaded']) |
            (~self.download_targets['features_downloaded'])
        )
        pending_dl_assemblies = (
            self.download_targets[pending_dl_mask]['assembly_accession']
            .unique().tolist()
        )
        pending_dl_assemblies.sort()

        total_assemblies    = len(self.assemblies)
        total_dl_pending    = len(pending_dl_assemblies)

        print("\n" + "="*45)
        print(f"      PIPELINE STATUS REPORT  (Download)")
        print("="*45)
        print(f"Total Assemblies      : {total_assemblies}")
        print(f"Completed Assembly    : {self.download_targets['assembly_downloaded'].sum()}")
        print(f"Completed Proteins    : {self.download_targets['proteins_downloaded'].sum()}")
        print(f"Completed Features    : {self.download_targets['features_downloaded'].sum()}")
        print("-" * 45)
        print(f"Assemblies to Download: {total_dl_pending}/{total_assemblies}")
        print(f"Batch Size            : {batch_size}")
        print("="*45 + "\n")

        dl_batches = [
            pending_dl_assemblies[i : i + batch_size]
            for i in range(0, total_dl_pending, batch_size)
        ]

        if total_dl_pending == 0:
            print("All files already downloaded.")
        else:
            with tqdm(total=len(dl_batches), desc="Downloading", unit="batch") as pbar:
                for chunk in dl_batches:
                    self._download_batch_files(chunk)
                    pbar.update(1)

        # ----------------------------------------------------------------
        # Phase 2: Build NCBI metadata (resolves genomic_region per protein)
        # ----------------------------------------------------------------
        print("\nBuilding NCBI metadata (resolving genomic regions)...")
        self._get_ncbi_metadata()

        # ----------------------------------------------------------------
        # Phase 3: Schedule extraction (contig-centric view)
        # ----------------------------------------------------------------
        print("Initializing Extraction Schedule...")
        self._schedule_extraction()

        completed_contigs = self.extraction_targets['extracted_contig_file'].sum()
        completed_cds     = self.extraction_targets['extracted_cds_file'].sum()
        completed_gff     = self.extraction_targets['extracted_gff_file'].sum()
        failed_count      = self.extraction_targets['extraction_failed'].sum()
        total_contigs     = len(self.extraction_targets)

        pending_ex_mask = (
            (
                (~self.extraction_targets['extracted_contig_file']) |
                (~self.extraction_targets['extracted_cds_file'])    |
                (~self.extraction_targets['extracted_gff_file'])
            ) &
            (~self.extraction_targets['extraction_failed'])
        )

        # We batch by assembly (not contig) so we open each large file only once
        pending_ex_assemblies = (
            self.extraction_targets[pending_ex_mask]['assembly_accession']
            .unique().tolist()
        )
        pending_ex_assemblies.sort()
        total_ex_pending = len(pending_ex_assemblies)

        print("\n" + "="*45)
        print(f"      PIPELINE STATUS REPORT  (Extraction)")
        print("="*45)
        print(f"Total Unique Contigs  : {total_contigs}")
        print(f"  - Failed / Skipped  : {failed_count}")
        print("-" * 45)
        print(f"Completed Contigs     : {completed_contigs}")
        print(f"Completed CDS         : {completed_cds}")
        print(f"Completed GFF         : {completed_gff}")
        print("-" * 45)
        print(f"Assemblies to Process : {total_ex_pending}/{len(self.unique_assemblies)}")
        print(f"Batch Size            : {batch_size}")
        print("="*45 + "\n")

        ex_batches = [
            pending_ex_assemblies[i : i + batch_size]
            for i in range(0, total_ex_pending, batch_size)
        ]

        if total_ex_pending == 0:
            print("All targets already extracted. Pipeline complete.")
            return

        with tqdm(total=len(ex_batches), desc="Extracting", unit="batch") as pbar:
            for chunk in ex_batches:
                pbar.set_postfix_str("Extracting sequences + GFF")
                self._extract_batch_sequences(chunk)
                #pbar.set_postfix_str("Cleaning up raw files")
                #self._cleanup_batch_files(chunk)
                pbar.update(1)

        print("\nPipeline Execution Finished.")

        """
        Main Pipeline Entry Point with Progress Tracking.
        """

        print("Initializing Pipeline Schedules...")

        self._schedule_download()

        # --- 1. Status Report ---
        total_mgycs = len(self.download_targets)
        
        # Count successes
        completed_assembly = self.download_targets['assembly_downloaded'].sum()
        completed_proteins = self.download_targets['proteins_downloaded'].sum()
        completed_features = self.download_targets['features_downloaded'].sum()
        
        # Count failures
        # (Assuming 'extraction_failed' column exists from schedule_extraction)
        #failed_count = self.download_targets['extraction_failed'].sum()

        # Determine which ERZs actually need work
        # Logic: An ERZ is pending if ANY of its MGYCs are missing ANY file AND are NOT failed.
        # We exclude failed rows entirely from the "To Process" list.
        pending_mask = (
            (
                (~self.download_targets['assembly_downloaded']) | 
                (~self.download_targets['proteins_downloaded']) |
                (~self.download_targets['features_downloaded']) 
            ) #& (~self.download_targets['extraction_failed'])
        )
        
        pending_assemblies = self.download_targets[pending_mask]['assembly_accession'].unique().tolist()
        pending_assemblies.sort()
        total_pending_assemblies = len(pending_assemblies)
        total_assemblies = len(self.assemblies)

        print("\n" + "="*45)
        print(f"      PIPELINE STATUS REPORT")
        print("="*45)
        print(f"Total Targets (MGYCs) : {total_mgycs}")
        #print(f"  - Failed / Skipped  : {failed_count}")
        print("-" * 45)
        print(f"Completed Contigs     : {completed_assembly}")
        print(f"Completed CDS         : {completed_proteins}")
        print(f"Completed GFF         : {completed_features}")
        print("-" * 45)
        print(f"Assemblies to Process : {total_pending_assemblies}/{total_assemblies}")
        print(f"Batch Size            : {batch_size}")
        print("="*45 + "\n")

        erz_batches = [ pending_assemblies[i : i + batch_size] 
                        for i in range(0, total_pending_assemblies, batch_size) ]
        
        if total_pending_assemblies == 0:
            print("All targets are already downloaded. Pipeline complete.")        
        else:
            with tqdm(total=len(erz_batches), desc="Pipeline Progress", unit="batch") as pbar:
                
                for chunk in erz_batches:

                    # 1. Download (Only if files missing)
                    pbar.set_postfix_str("Downloading")
                    #print("Downloading")
                    self._download_batch_files(chunk)

        # 2. Make metadata
        print("Making NCBI metadata")
        self._get_ncbi_metadata()
                
        with tqdm(total=len(erz_batches), desc="Pipeline Progress", unit="batch") as pbar:

            for chunk in erz_batches:

                # 2. Extract (Only if files present)
                pbar.set_postfix_str("Extracting sequences")
                #print("Extracting sequences")
                #self._extract_batch_sequences(chunk)
    
                # 3. Cleanup (Always delete raw files to save space)
                pbar.set_postfix_str("Cleaning")
                #print("Cleaning")
                self._cleanup_batch_files(chunk)
                
                pbar.update(1)
            #return


        print("\nPipeline Execution Finished.")