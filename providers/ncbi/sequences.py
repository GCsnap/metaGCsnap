import json
from Bio import SeqIO
import gzip
import os

from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.utils import processpool_wrapper

from providers.ncbi.dataset import Dataset


def extract_target_sequences(input_file: str, targets: list) -> dict:
    """
    Read the input FASTA file (gzipped) and return all sequences whose record.id exactly matches any target.
    Returns a dictionary: {id: {'seq': sequence}}
    Handles missing files gracefully by returning an empty dict.
    """
    result = {}
    try:
        with gzip.open(input_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in targets:
                    result[record.id] = {'seq': str(record.seq)}
    except FileNotFoundError:
        return {}
    return result

class Sequences:
    """ 
    Methods and attributes to get the sequences for the flanking genes of the target genes.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cores (int): The number of CPU cores to use.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        sequences (dict): The dictionary with the sequences of the flanking genes.
        console (RichConsole): The RichConsole object to print messages.
    """

    def __init__(self, config: Configuration, gc: GenomicContext, dataset: Dataset) -> None:
        """
        Initialize the Sequences object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
        """        
        self.config = config
        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value'] 

        # set arguments
        self.gc = gc

        self.ncbi_metadata = dataset.ncbi_metadata
        self.proteins_file = dataset.contig_proteins_file

        self.console = RichConsole()

    def get_sequences(self) -> dict:
        """
        Getter for the sequences attribute.

        Returns:
            dict: The dictionary with flanking genes and their sequences.
        """        
        return self.genomic_context        

    def run(self) -> None:
        """
        Run the assignment of sequences to the flanking genes.
            - Find sequences for all flanking genes.
            - Add sequences, tax id and species name to flanking genes.
        Uses parallel processing with the processpool_wrapper from utils.py.
        """        
        # Find sequnces for all cds codes
        #print(self.gc.get_all_cds_codes())
        self.find_sequences(self.gc.get_all_cds_codes())
        self.find_taxonomy()

        # Prepare a list of tuples (target, dict_for_target)
        # here each process gets one target: {} combination
        # henve we have many different processes
        parallel_args = self.gc.get_syntenies_key_value_list()

        with self.console.status('Add sequences, tax id and species name to flanking genes'):
            dict_list = processpool_wrapper(self.cores, parallel_args, self.run_each)
            # combine results
            self.genomic_context = {k: v for d in dict_list for k, v in d.items()}

    def run_each(self, args: tuple[str,dict]) -> dict:
        """
        Run the assignment of sequences to the flanking genes for one target used
        in parallel processing.

        Args:
            args (tuple[str,dict]): The arguments for the sequence assignment.
                First element is the target gene.
                Second element is the dictionary with the flanking genes of the target gene.

        Returns:
            dict: The dictionary with the flanking genes and their sequences.
        """        
        target, content_dict = args
        # update flanking genes with sequence
        sequences = [self.get_sequence(cds_code) for cds_code in content_dict['flanking_genes']['cds_codes']]
        content_dict['flanking_genes']['sequences'] = sequences

        # add species and taxid for target_cds code (first one in the list)
        target_cds = content_dict['assembly_metadata']['target_cds']
        # species in contained twice in the dict

        content_dict['taxonomy'] = self.get_taxonomy(target_cds)
        #content_dict['taxonomy'] = content_dict['flanking_genes']['taxonomy']
        #content_dict['flanking_genes']['taxID'] = self.get_taxid(target_cds)

        return {target: content_dict}

    def find_taxonomy(self):
        
        self.taxonomy = {}

        for _, row in self.ncbi_metadata.iterrows():

            self.taxonomy[row['ncbi_code']]={'taxon_id':row['taxid'],'taxon_name':row['organism_name']}

    def get_taxonomy(self, cds_code: str) -> str:

        try:
            tax = self.taxonomy.get(cds_code, {})   
        except Exception:
            tax = {'taxon_id':'unk','taxon_name':'unk'}

        return tax
    
    def find_sequences(self, cds_codes: list) -> None:
        """
        Find sequences for all flanking genes.
        Assign to each contig its taxonomy. For each genomic context, open the mgyc file and get the sequences. Assign taxonomy like the contig.

        Args:
            cds_codes (list): The list of cds codes to find sequences for.
        """    

        self.sequences = {}

        for v in self.gc.get_syntenies().values():

            genomic_region = v["assembly_metadata"]['genomic_region']
            target_proteins = v["flanking_genes"]["cds_codes"]
            protein_file = self.proteins_file[genomic_region]

            seqs=extract_target_sequences(protein_file,target_proteins)

            self.sequences.update(seqs)

    def get_sequence(self, cds_code: str) -> str:
        """
        Get the sequence for a cds code. 
        If the sequence is not found, a fake sequence is returned.

        Args:
            cds_code (str): The cds code.

        Returns:
            str: The sequence for the cds code.
        """        
        entry = self.sequences.get(cds_code, {})        
        return entry.get('seq', 'FAKESEQUENCEFAKESEQUENCEFAKESEQUENCEFAKESEQUENCE')