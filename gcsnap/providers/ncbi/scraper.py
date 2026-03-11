import os
import gzip # to work with .gz
import urllib.request
import time
from datetime import datetime, timedelta
import pandas as pd

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.utils import CustomLogger
from gcsnap.utils import processpool_wrapper
from gcsnap.utils import WarningToLog
from gcsnap.genomic_context import GenomicContext

from gcsnap.providers.ncbi.entrez_query import EntrezQuery
from gcsnap.providers.ncbi.dataset import Dataset
from gcsnap.providers.ncbi.assembly_links import AssemblyLinks

import logging
logger = logging.getLogger('iteration')

class Scraper:    
    """
    Methods and attributes to download and parse flanking genes given NCBI codes.

    Attributes:
        cores (int): Number of cores to use for parallel processing.
        config (Configuration): Configuration object.
        console (RichConsole): Console object.
        assembly_dir (str): Path to store assembly summaries.
        targets_and_cds_codes (list): List of tuples with target and ncbi code.
        accessions (dict): Dictionary with ncbi codes and assembly accessions.
    """

    def __init__(self, dataset: Dataset, config: Configuration, mappings: list[tuple[str,str]]):                     
        """
        Initialize the Assemblies object.

        Args:
            config (Configuration): Configuration object containing the arguments.
            mappings (list[tuple[str,str]]): Contains the target and its ncbi code.
        """        
        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value'] 
        self.config = config
        self.console = RichConsole()
        self.age = config.arguments['assemblies_data_update_age']['value'] 

        # input list with [(target, cds_code)]
        self.targets_and_cds_codes = mappings
        self.metadata_dir = dataset.metadata_dir
    
    def run(self) -> None:
        """
        Run the process to download and extract flanking genes for the targets:
            - Find the assembly accessions for the given NCBI codes.
            - Download and extract flanking genes for each target in parallel.
        Uses parallel processing with the processpool_wrapper from utils.py
        """        
        # find the assembly accessions
        
        cds_codes = [target[1] for target in self.targets_and_cds_codes]
        
        # assign the assembly to the ncbi codes
        self.find_accessions(cds_codes)

        # download and parse the assembly summary files
        self.load_summaries()

        #list target files to download
        self.list_targets()

    def load_summaries(self) -> None:
        """
        Load the assembly summaries.
        """        
        # get file with assembly links

        self.assembly_links = AssemblyLinks(out_dir=self.metadata_dir, config=self.config)
        self.assembly_links.run( self.accessions )
        self.assembly_metadata = self.assembly_links.assembly_metadata
   
    def find_accessions(self, cds_codes: list) -> None:
        """
        Find the assembly accessions for the given NCBI codes using the EntrezQuery class.

        Args:
            cds_codes (list): The list of NCBI codes.
        """        
        # get file with assembly links, no logging as its done after run()
        entrez = EntrezQuery(self.config, cds_codes, db='protein', rettype='ipg', 
                             retmode='xml', logging=False)
        
        self.accessions = entrez.run()

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
        accession = self.accessions.get(cds_code, 'unk')
        if accession == 'unk':
            raise WarningToLog('No assembly accession found for {}'.format(cds_code))
        return accession

    def list_targets(self) -> str:
        """
        Get the assembly URL for a given assembly accession.

        Args:
            assembly_accession (str): The assembly accession.

        Raises:
            WarningToLog: If no URL is found.

        Returns:
            str: The assembly URL.
        """        
        
        def expand_row(row):
            ftp_path = row['ftp_path'].rstrip('/')
            assembly_name = ftp_path.split('/')[-1]
            accession = row['#assembly_accession']
            
            # Create three rows
            return [
                {
                    'assembly_accession': accession,
                    'file_type': 'genomic_fna',
                    'ftp_url': f"{ftp_path}/{assembly_name}_genomic.fna.gz"
                },
                {
                    'assembly_accession': accession,
                    'file_type': 'protein_faa',
                    'ftp_url': f"{ftp_path}/{assembly_name}_protein.faa.gz"
                },
                {
                    'assembly_accession': accession,
                    'file_type': 'genomic_gff',
                    'ftp_url': f"{ftp_path}/{assembly_name}_genomic.gff.gz"
                }
            ]

        # Expand all rows
        all_rows = []
        for idx, row in self.assembly_metadata.iterrows():
            
            all_rows.extend(expand_row(row))

        # Create new dataframe
        self.targets = pd.DataFrame(all_rows)

        self.targets.to_csv(self.metadata_dir / 'target_files.csv', index=False)
        
    
