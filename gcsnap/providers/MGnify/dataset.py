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
from concurrent.futures import ThreadPoolExecutor, as_completed

# In-house libraries
from gcsnap.configuration import Configuration
from gcsnap.utils import handle_compressed_fasta

import gcsnap.providers.MGnify.helpers as hp

##########################################
### helper for the list_target_files method

class Dataset:
    """
    A class to manipulate and manage local files related to MGnify data.
    """

    def __init__(self, config, basename: str):
        """
        Initialize the MGnifyLocalDB class. Here you just set the files/ variables that are needed at initialization time

        Args:
            output_dir (str or Path): The base directory where local files are stored.
        """

        self.MGnify_local = Path(config.arguments['MGnify_path']['value']) #Path(config.arguments['data_path']['value']) / 'MGnify'
        self.GTDB_local = Path(config.arguments['kraken_path']['value']) #Path(config.arguments['data_path']['value']) / 'kraken' / 'GTDB'
        #self.taxonomy_local = Path(config.arguments['taxonomy_path']['value']) #Path(config.arguments['data_path']['value']) / 'taxonkit' / 'GTDB'
        self.output_dir = Path(config.arguments['MGnify_dir']['value'] ) #[]Path(config.arguments['out_label']['value'] )
        #self.query_fasta = config.targets[0] # this should be modified. orignal GCsnap target is just a list of id, while here we have only one fasta file. ideally, we can give it many fasta files, let's see.
        self.query_basename = basename
        
        # scraper
        self.MGnify_API = hp.MGnify_API

        self.metadata_dir = self.output_dir / 'metadata'
        self.seq_search_dir = self.output_dir / 'mmseq_targets'
        self.tmp_dir = Path( config.arguments['tmp_folder']['value'] )

        for d in [self.metadata_dir, self.seq_search_dir, self.tmp_dir]:
            os.makedirs(d, exist_ok=True)
        
        self.mgyp_metadata_file = self.metadata_dir / hp.mgyp_metadata
        self.targets_file = self.metadata_dir / 'target_files.csv'
        self.assembly_metadata_file = self.metadata_dir / 'assembly_metadata.csv'

        # mmseqs
        self.seq_search = {} # a dict to store all sequence search session instructions

        self.genome_classification = config.arguments['genome_classification']['value']

    def update_metadata(self):
        
        ''' To be runned after metadata.aggregate_assemblies_metadata(). so that the api_metadata and assembly_metadata can be read.'''
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

        self.mgyp_present = self.mgyp_metadata_file.exists()
        if self.mgyp_present:
            self.mgyp_metadata = pd.read_csv(self.mgyp_metadata_file)
        else:
            raise ValueError(f"mgyp_metadata_file.csv file not found.")

    def set_sequences_search(self,query_fasta):
        
        self.reference_db = self.MGnify_local / hp.mgyp_database  # Path to the precomputed reference database

        self.query_fasta = query_fasta  # Path to the query contig FASTA file
        # Path to the query FASTA file
        
        for ext in [".fasta", ".faa", ".fa", ".gz"]:
            if self.query_basename.endswith(ext):
                self.query_basename = self.query_basename[: -len(ext)]
        
        # mmseq directories
        self.query_db = str( self.seq_search_dir / 'query/DB' )  # Path to the query database
        self.result_db = str( self.seq_search_dir / f'result/{self.query_basename}')  # Path to the query database

        for d in [self.query_db,self.result_db]:
            os.makedirs(os.path.dirname(d), exist_ok=True)
        
        self.mmseqs_output_file = self.seq_search_dir / hp.mgyp_search_out  # Path to the output directory
        self.mmseqs_hits_fasta = self.seq_search_dir / hp.mgyp_search_fasta  # Path to the output directory
        self.mmseqs_hits_ids = self.seq_search_dir / hp.mgyp_search_ids  # Path to the output directory

        # Handle compressed or uncompressed FASTA
        self.query_fasta_uncompressed = handle_compressed_fasta(self.query_fasta)

        search_dict = { 'query_fasta_uncompressed': self.query_fasta_uncompressed,
                        'query_basename': self.query_basename,
                        'query_db': self.query_db, 'result_db': self.result_db, 'reference_db': self.reference_db, 
                        'mmseqs_output_file': self.mmseqs_output_file, 'hits_fasta': self.mmseqs_hits_fasta, 'hits_ids': self.mmseqs_hits_ids}

        self.seq_search['focal'] = search_dict

    def update_after_sequences_search(self,ids_file=None):
        
        if ids_file is None:
            mmseqs_output = pd.read_csv(self.mmseqs_output_file, sep='\t', header=None)
            self.mgyps = list(set(mmseqs_output[1].tolist()))
            print("Total hits: {}".format(mmseqs_output.shape[0]))
        else:
            with open(ids_file, 'r') as f:
                self.mgyps = [line.strip() for line in f]

        
        print("Unique hits: {}".format(len(self.mgyps)))

    def set_targets(self,ids):
        
        if isinstance(ids, str) and Path(ids).is_file():
            with open(ids, 'r') as f:
                self.mgyps = [line.strip() for line in f if line.strip()]
        else:
            self.mgyps = list(ids)

        print("Unique hits: {}".format(len(self.mgyps)))

    def set_scraper(self):

        self.mgyp_present = self.mgyp_metadata_file.exists()
        
        if self.mgyp_present:
            self.mgyp_metadata = pd.read_csv(self.mgyp_metadata_file, index_col='ERZ')
        else:
            raise ValueError(f"MGnify metadata file not found: {self.mgyp_metadata_file}")
        
        # retain only the assemblies for which you were able to get the ERZ-contig
        self.assemblies = self.mgyp_metadata.dropna(subset='ERZ_contig').index.unique().tolist()
        self.assemblies.sort()
        
        self.assemblies_dir = self.output_dir / 'assemblies' / 'assemblies' # DNA contigs / assemblies sequences
        self.proteins_dir = self.output_dir / 'assemblies' / 'assemblies_cds' # proteins predicted from assemblies
        self.annotations_dir = self.output_dir / 'assemblies' / 'annotations'  # functional annotations for assembly content
        self.assemblies_metadata_dir = self.output_dir / 'assemblies' / 'metadata' # metadata about the assembly experiment (biome, instrument, etc)
        self.assembly_targets_dir = self.output_dir / 'assemblies' / 'targets' # target files for each assembly, see hp.pipelineV5_targets and hp.pipelineV4_targets
        self.assembly_allfiles_dir = self.output_dir / 'assemblies' / 'all_files' # all files associated to an ERZ entry that can be potentially obtained from MGnify

        for d in [self.assemblies_dir, self.proteins_dir, self.annotations_dir, self.assemblies_metadata_dir, self.assembly_targets_dir, self.assembly_allfiles_dir]:
            os.makedirs(d, exist_ok=True)

        self.targets_file_present = self.targets_file.exists()
        self.assembly_present = self.assembly_metadata_file.exists()

        self.failed_assemblies_file = self.metadata_dir / 'failed_assemblies.txt'
                    
    def update_after_scrape(self):
        
        ''' To be runned after metadata.aggregate_assemblies_metadata(). so that the api_metadata and assembly_metadata can be read.'''
        self.mgyp_present = self.mgyp_metadata_file.exists()
        if self.mgyp_present:
            self.mgyp_metadata = pd.read_csv(self.mgyp_metadata_file)
        else:
            raise ValueError(f"{self.mgyp_metadata_file} file not found.")

        self.targets_file_present = self.targets_file.exists()
        if self.targets_file_present:
            self.targets = pd.read_csv(self.targets_file)

        self.assembly_metadata_file_present = self.assembly_metadata_file.exists()
        if self.assembly_metadata_file_present:
            self.assembly_metadata = pd.read_csv(self.assembly_metadata_file)

        self.failed_assemblies = []
        if self.failed_assemblies_file.exists():
            with open(self.failed_assemblies_file, 'r') as f:
                self.failed_assemblies = [line.strip() for line in f if line.strip()]

    def set_gatherer(self):

        self.contigs_dir = self.output_dir / 'contigs' / 'contigs'
        self.contigs_proteins_dir = self.output_dir / 'contigs' / 'proteins'
        self.contigs_gff_dir = self.output_dir / 'contigs' / 'gff'
        self.contigs_ann_dir = self.output_dir / 'contigs' / 'annotations'
        
        for c in [self.contigs_dir, self.contigs_proteins_dir,self.contigs_gff_dir,self.contigs_ann_dir]:
            os.makedirs(c, exist_ok=True)

    def update_after_gathering(self,gatherer):
        
        #self.contig_file = contigs.contig_file
        #self.contig_proteins_file = contigs.contig_proteins_file
        #self.contig_gff_file = contigs.contig_gff_file

        gat = gatherer.extraction_targets
        condition = (gat['extracted_contig_file'] & gat['extracted_cds_file'] & gat['extracted_gff_file'] & gat['extracted_fannot_file'])
        gat = gat[ condition ]

        self.contig_file = {row.name: row['contig_file'] for _,row in gat.iterrows()}
        self.contig_proteins_file = {row.name: row['cds_file'] for _,row in gat.iterrows()}
        self.contig_gff_file = {row.name: row['gff_file'] for _,row in gat.iterrows()}
        

    def set_contig_binning(self):
        
        #self.query_contig_fasta = self.output_dir / 'contigs.fna.gz'  # Path to the query contig FASTA file
        self.query_contig_fasta = [ str(s) for s in self.contig_file.values() if os.path.exists(s) ]

        self.binning_out_dir = self.output_dir / 'binning'
        self.binning_out_dir.mkdir(parents=True, exist_ok=True)
        self.sourmash_output = self.binning_out_dir / 'contigs.sig'
        
        # Path to the query FASTA file

        for ext in [".fasta", ".faa", ".fa", ".gz"]:
            if self.query_basename.endswith(ext):
                self.query_basename = self.query_basename[: -len(ext)]

    def update_after_contig_binning(self, binning):
        
        self.binning_distance_matrix = binning.distance_matrix_file

        distance_matrix = pd.read_csv(self.binning_distance_matrix, compression='gzip', index_col='contig')
        labels = distance_matrix.index.tolist()
        
        mgyc = self.mgyp_metadata['MGYC'].unique().tolist()
        contig_lengths = {}
        for m in mgyc:
            fasta_path = self.contigs_dir / f"{m}.fna.gz"
            length = 0
            try:
                with gzip.open(fasta_path, "rt") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        length += len(record.seq)
                contig_lengths[m] = length
            except Exception as e:
                contig_lengths[m] = 0

        with open(self.binning_out_dir / 'bins.json', 'r') as f:
            bins = json.load(f)
        
        self.mgyp_metadata['ERZ_contig_assigned'] = self.mgyp_metadata['ERZ_contig'].isin(labels)

        filler_rank = 'metagenomic bin' if self.genome_classification == 'binning' else 'root'
        filler = {'taxon_id':'0', 'rank':filler_rank}
        
        for v,f in filler.items():

            self.mgyp_metadata[v]=f
        
        self.mgyp_metadata['taxon_name'] = 'bin '+self.mgyp_metadata['ERZ_contig'].map(bins).astype(str)
        self.mgyp_metadata['contig_length'] = self.mgyp_metadata['ERZ_contig'].map(contig_lengths)
        self.mgyp_metadata.to_csv(self.mgyp_metadata_file, index=False)

        self.viable_cds = self.mgyp_metadata[self.mgyp_metadata['ERZ_contig_assigned']]['ERZ_cds_id'].unique().tolist()
        self.viable_cds = [ v for v in self.viable_cds if 'cov' in str(v)]


    def set_taxonomic_assignment(self):
        
        #self.query_contig_fasta = self.output_dir / 'contigs.fna.gz'  # Path to the query contig FASTA file
        self.query_contig_fasta = [ str(s) for s in self.contig_file.values()]

        self.taxonomic_out_dir = self.output_dir / 'taxonomy'
        self.taxonomic_out_dir.mkdir(parents=True, exist_ok=True)
        self.kraken2_output = self.taxonomic_out_dir / 'contigs_output.tsv'
        self.kraken2_report = self.taxonomic_out_dir / 'contigs.report'
        
        # Path to the query FASTA file

        for ext in [".fasta", ".faa", ".fa", ".gz"]:
            if self.query_basename.endswith(ext):
                self.query_basename = self.query_basename[: -len(ext)]

        self.kraken2_db = self.GTDB_local 

    def update_after_taxonomic_assignment(self):
        
        self.kraken2_output = self.output_dir / 'taxonomy' / 'contigs_output.tsv'

        taxonomy = pd.read_csv(self.kraken2_output, sep='\t', dtype=str)

        self.mgyp_metadata['ERZ_contig_assigned'] = self.mgyp_metadata['ERZ_contig'].isin(taxonomy['ERZ-contig'].values)

        filler = {'taxon_id':'-1','taxon_name':'exception','length':'-1', 'rank':'exception'}

        for v,f in filler.items():
            r = taxonomy[['ERZ-contig',v]].set_index('ERZ-contig')[v].to_dict()
            
            self.mgyp_metadata[v]=self.mgyp_metadata['ERZ_contig'].map(r)
            self.mgyp_metadata[v]=self.mgyp_metadata[v].fillna(f)
        
        self.mgyp_metadata.rename(columns={'length':'contig_length'}, inplace=True)
        self.mgyp_metadata.to_csv(self.mgyp_metadata_file, index=False)

        self.viable_cds = self.mgyp_metadata[self.mgyp_metadata['ERZ_contig_assigned']]['ERZ_cds_id'].unique().tolist()
