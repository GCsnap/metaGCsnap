import json
from Bio import SeqIO
import gzip
import os

from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.utils import processpool_wrapper

from gcsnap.providers.ncbi.dataset import Dataset


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

def extract_dna_sequences(dna_file: str, genomic_region: str, starts: list, ends: list) -> dict:
    """
    Read the contig DNA FASTA file (gzipped) and extract the genomic context and individual
    gene sequences using 1-based genomic coordinates.

    Args:
        dna_file (str): Path to the gzipped FASTA file containing the contig DNA.
        genomic_region (str): Contig record ID to locate within the FASTA file.
        starts (list): 1-based start positions of the flanking genes.
        ends (list): 1-based end positions of the flanking genes.

    Returns:
        dict: {'context': str, 'features': [str, ...]}
              'context' is the full genomic context spanning all flanking genes.
              'features' is a list of individual gene DNA sequences in the same order as starts/ends.
    """
    try:
        with gzip.open(dna_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == genomic_region:
                    contig_seq = record.seq
                    gc_start = min(starts) - 1  # convert 1-based to 0-based
                    gc_end   = max(ends)         # 1-based inclusive end = 0-based exclusive end
                    context  = str(contig_seq[gc_start:gc_end])
                    features = [str(contig_seq[s - 1 : e]) for s, e in zip(starts, ends)]
                    return {'context': context, 'features': features}
    except FileNotFoundError:
        pass
    return {'context': '', 'features': []}

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
        self.find_protein_sequences(self.gc.get_all_cds_codes())
        self.find_genomic_sequences()
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
        dna_sequences = [self.get_sequence(cds_code, 'genomic') for cds_code in content_dict['flanking_genes']['cds_codes']]
        content_dict['flanking_genes']['dna_sequences'] = dna_sequences

        protein_sequences = [self.get_sequence(cds_code, 'proteins') for cds_code in content_dict['flanking_genes']['cds_codes']]
        content_dict['flanking_genes']['protein_sequences'] = protein_sequences

        # add genomic context sequence
        content_dict['context_sequence'] = self.context_sequences.get(target, '')

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
    
    def find_protein_sequences(self, cds_codes: list) -> None:
        """
        Find protein sequences for all flanking genes.

        Args:
            cds_codes (list): The list of cds codes to find sequences for.
        """
        self.protein_sequences = {}

        for v in self.gc.get_syntenies().values():

            genomic_region = v["assembly_metadata"]['genomic_region']
            target_proteins = v["flanking_genes"]["cds_codes"]
            protein_file = self.proteins_file[genomic_region]

            seqs = extract_target_sequences(protein_file, target_proteins)

            self.protein_sequences.update(seqs)

    def find_genomic_sequences(self) -> None:
        """
        Extract DNA sequences for all flanking genes and their genomic contexts.
        For each synteny entry, open the contig DNA file and slice out the context region
        and each individual gene using the stored 1-based genomic coordinates.

        Sets:
            self.genomic_sequences (dict): {cds_code: {'seq': dna_str}} for each flanking gene.
            self.context_sequences (dict): {target: dna_str} for each synteny context.
        """
        self.genomic_sequences = {}
        self.context_sequences = {}

        for k, v in self.gc.get_syntenies().items():
            dna_file       = v["assembly_metadata"]["dna_file"]
            genomic_region = v["assembly_metadata"]["genomic_region"]
            cds_codes      = v["flanking_genes"]["cds_codes"]
            starts         = v["flanking_genes"]["starts"]
            ends           = v["flanking_genes"]["ends"]

            seqs = extract_dna_sequences(dna_file, genomic_region, starts, ends)

            for cds_code, dna_seq in zip(cds_codes, seqs['features']):
                self.genomic_sequences[cds_code] = {'seq': dna_seq}

            self.context_sequences[k] = seqs['context']

    def get_sequence(self, cds_code: str, alphabet: str) -> str:
        """
        Get the sequence for a cds code.
        If the sequence is not found, a fake sequence is returned.

        Args:
            cds_code (str): The cds code.
            alphabet (str): 'genomic' for DNA sequences, 'proteins' for protein sequences.

        Returns:
            str: The sequence for the cds code.
        """
        if alphabet == 'genomic':
            entry = self.genomic_sequences.get(cds_code, {})
        if alphabet == 'proteins':
            entry = self.protein_sequences.get(cds_code, {})
        return entry.get('seq', 'FAKESEQUENCEFAKESEQUENCEFAKESEQUENCEFAKESEQUENCE')