import os
from datetime import datetime, timedelta
import time
# pip install pypdl
from pypdl import Pypdl

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

from gcsnap.providers.ncbi.dataset import Dataset
import pandas as pd
from pathlib import Path
from typing import Union

class AssemblyLinks:    
    """
    Methods and attributes to download and parse the assembly summaries from NCBI
    for RefSeq and Genbank to retreive the links to the assemblies.

    Attributes:
        console (RichConsole): The RichConsole object to print messages.
        cores (int): The number of CPU cores to use.
        assembly_dir (str): The path to store the assembly summaries.
        db_list (list): The list of databases to download the assembly summaries.
        links (dict): The final dictionary with the assembly code and the link to it.
    """

    def __init__(self, out_dir: Union[str, Path], config: Configuration): 
        """
        Initialize the AssemblyLinks object.

        Args:
            out_dir (Union[str, Path]): The output directory for assembly summaries.
            config (Configuration): The Configuration object containing the arguments.
        """           
        self.console = RichConsole()

        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value']
        self.age = config.arguments['assemblies_data_update_age']['value']

        #parent_path = os.path.join(config.arguments['ncbi_dir']['value'],'assemblies')
        
        #if parent_path is None:
            # set path to store assembly summaries
        #    parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        self.metadata_dir = out_dir
        
        # database list
        self.db_list = ['genbank','refseq']
        
        # final dictionaire with {assembly code: link}
        self.links = {}
    
    def run(self,accessions) -> None:    
        """
        Run the download and parsing of the assembly summaries.
        """        

        assembly_metadata = []

        for db in self.db_list:
            
            self.check_and_download(db)  
            assembly_metadata.append(self.parse_summaries(db,accessions))      

        self.assembly_metadata = pd.concat(assembly_metadata)  
        self.assembly_metadata = self.assembly_metadata.reset_index(drop=True)

        self.assembly_metadata.to_pickle(os.path.join(self.metadata_dir,'assembly_metadata.pkl'))

    def get(self) -> dict[str, str]:
        """
        Getter for the links attribute.

        Returns:
            dict[str, str]: The dictionary with the assembly code and the link to it.
        """        
        return self.links              
            
    def check_and_download(self, db) -> None:
        """
        Check if the assembly summary exists and download it if it does not.
        Using time check to download again if older than 14 days.

        Args:
            db (_type_): The database (genbank or refseq) to download the assembly summary.
        """        
        file = os.path.join(self.metadata_dir,'assembly_summary_{0}.txt'.format(db))

        days=self.age
        if os.path.exists(file):
            # time check, download again if older than days
            file_time = datetime.fromtimestamp(os.path.getmtime(file))
            if datetime.now() - file_time > timedelta(days=days):
                # print('Assembly summary is older than {} days, downloading again.'.format(days))
                self.console.print_info('Assembly summary {} is older than {} days, downloading again.'.format(db, days))
                self.download_summary(db)
            else:
                # print('Assembly summary is not older than {} days, not downloading.'.format(days))
                self.console.print_info('Assembly summary {} is not older than {} days, not downloading.'.format(db, days))
        else:
            # print('Assembly summary does not exist, downloading.')
            self.console.print_info('Assembly summary {} does not exist, downloading.'.format(db))
            self.download_summary(db)

    def download_summary(self, db: str) -> None:
        """
        Download the assembly summary from NCBI in parallel with pypdl Downloader.

        Args:
            db (str): The database (genbank or refseq) to download the assembly summary.
        """        
        url = 'https://ftp.ncbi.nlm.nih.gov/genomes/{0}/assembly_summary_{0}.txt'.format(db) 

        # Download multithreaded with pypdl Downloader
        dl = Pypdl()
        dl.start(url=url,
                file_path=str(self.metadata_dir), 
                segments=self.cores, 
                multisegment=True, 
                block=False, 
                display=False)
        
        # retreive total size of the file does not work, the file does not have this information
        # total_size = dl.size
        # print(f"Total size: {total_size}")

        # set it manually (lower than the actual size, but it is enough for the progress bar)
        total_size = 1 * 10**9 if db == self.db_list[0] else 0.1 * 10**9

        # Download multithreaded with pypdl Downloader
        with self.console.progress('Downloading {} assembly summary'.format(db), 
                                   total=total_size) as (progress, task_id):
             
            while not dl.completed:
                current_size = dl.current_size
                progress.update(task_id, completed=current_size)
                time.sleep(1)   

        dl.shutdown() 
        
    def parse_summaries(self, db: str, accessions: list) -> dict[str,str]:
        """
        Parse the assembly summary to extract the assembly code and the link to it.

        Args:
            accessions (list): The list of assembly accessions to parse.
            db (str): The database (genbank or refseq) to parse the assembly summary.

        Returns:
            dict[str,str]: The dictionary with the assembly code and the link to it.
        """        

        file = os.path.join(self.metadata_dir,'assembly_summary_{0}.txt'.format(db))

        summary = pd.read_csv(file,sep='\t',skiprows=1,dtype=str)
        
        assembly_to_proteins = {}
        for protein, assembly in accessions.items():
            assembly_to_proteins.setdefault(assembly, []).append(protein)

        assemblies = list(accessions.values())
        summary = summary[summary['#assembly_accession'].isin(assemblies)]
        summary['ncbi_code'] = summary['#assembly_accession'].map(assembly_to_proteins)

        return summary