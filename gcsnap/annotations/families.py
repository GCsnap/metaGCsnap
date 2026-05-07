import os

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.utils import processpool_wrapper
from gcsnap.utils import split_dict_chunks

from gcsnap.annotations.mmseqs_cluster import MMseqsCluster

class Families:
    """ 
    Methods and attributes to assign families to flanking genes of the target genes.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cores (int): The number of CPU cores to use.
        out_label (str): The label of the output.
        out_dir (str): The path to store the output.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        syntenies (dict): The dictionary with the syntenies of the target genes.
        families (dict): The dictionary with the families assigned to the flanking genes.
        families_adapted (dict): The dictionary with the adapted families assigned to the flanking genes.
        cluster_list (list): The list of clusters.
        cluster_order (list): The order of the clusters.
        console (RichConsole): The RichConsole object to print messages.
    """

    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str=None, sensitivity: float = 6.0, f_community: float = 0.01):
        """
        Initialize the Families object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_label (str): The label of the output.
            sensitivity (float): MMseqs2 sensitivity parameter.
            f_community (float): MMseqs2 community fraction parameter.
        """
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.sensitivity = sensitivity
        self.f_community = f_community
        # set arguments
        
        if out_label is None:
            self.out_label = config.arguments['out_label']['value']
        else:
            self.out_label = out_label
        
        self.out_dir = os.path.join(out_label, 'genomic_context','sequences')     
        os.makedirs(self.out_dir, exist_ok=True)
        self.gc = gc
        self.syntenies = gc.get_syntenies()
        self.families_adapted = {}

        self.console = RichConsole()

    def get_families(self) -> dict:
        """
        Getter for the families attribute.

        Returns:
            dict: The dictionary with the families assigned to the flanking genes.
        """        
        return self.families_adapted

    def run(self) -> None:
        """
        Run the assignment of families to the flanking genes:
            - Find the clusters with MMseqsCluster.
            - Assign the families to the flanking genes.
            - Adapt the families where its outside possible ranges of the clusters.
        Uses parallel processing with processpool_wrapper from utils.py.
        Only done when more than one target gene is present.
        """        

        if len(self.syntenies.keys()) < 2:
            msg = 'Found assembly for only one target. Continue is not possible'
            self.console.stop_execution(msg)

        # MMseqsCluster creates the directory
        self.find_cluster()

        with self.console.status('Assigning families to flanking genes'):

            # 1. Add the family to the flanking genes
            # do in parallel, however, each needs the self.cluster_list and self.cluster_order
            # is this might be large, we use as many batches as there are cores
            parallel_args = split_dict_chunks(self.syntenies, self.cores)  
            # a list of tuple[dict, list] is returned
            result_list = processpool_wrapper(self.cores, parallel_args, self.assign_families)
            # combine results
            self.families = {k: v for tup in result_list for k, v in tup[0].items()}
            curr_numbers = [num for tup in result_list for num in tup[1]]
            # sort the curr_numbers and remove -1
            curr_numbers = sorted(list(set(curr_numbers)))
            if -1 in curr_numbers:
                curr_numbers.remove(-1)

            # 2. adapt the families where its outside possible ranges
            parallel_args = [(sub_dict, curr_numbers) 
                             for sub_dict in split_dict_chunks(self.families, self.cores)] 

            dict_list = processpool_wrapper(self.cores, parallel_args, self.adapt_families)
            # combine results
            self.families_adapted = {k: v for sub_dict in dict_list for k, v in sub_dict.items()}  

    def find_cluster(self) -> None:
        """Find protein family clusters using MMseqs2."""
        cluster = MMseqsCluster(
            self.config, self.gc, self.out_dir,
            sensitivity=self.sensitivity,
            f_community=self.f_community,
        )
        self.mmseqs_results_file = cluster.mmseqs_results
        self.fasta_file = cluster.fasta_file
        cluster.run()
        self.cluster_list = cluster.get_clusters_list()
        self.cluster_order = cluster.get_cluster_order()

    def assign_families(self, args: tuple[dict]) -> tuple[dict,list]:
        """
        Assign the family numbers to the flanking genes.

        Args:
            args (tuple[dict]): The syntenies of the target genes.

        Returns:
            tuple[dict,list]: The syntenies with the families assigned to the flanking genes and the current numbers.
        """         
        syntenies = args
        # loop over all targets in the chunk
        
        curr_numbers = []
        
        for k in syntenies.keys():
            
            syntenies[k]['flanking_genes']['families'] = []
            
            for i, ncbi_code in enumerate(syntenies[k]['flanking_genes']['cds_codes']):
                
                protein_name = syntenies[k]['flanking_genes']['names'][i]
                try:
                    protein_family = self.cluster_list[self.cluster_order.index(ncbi_code)]
                except:
                    protein_family = -1

                if protein_name == 'pseudogene':
                    # pseudogenes get a negative family number
                    protein_family = -1
                if ncbi_code == k:
                    # the target gene get the highest family number
                    protein_family = max(self.cluster_list)+1

                syntenies[k]['flanking_genes']['families'].append(protein_family)
                curr_numbers.append(protein_family)

                if ncbi_code == syntenies[k]['assembly_metadata']['target_cds']:
                    syntenies[k]['target_family'] = protein_family 

        return (syntenies , curr_numbers) 
    
    def adapt_families(self, args: tuple[dict]) -> dict:
        """
        Adapt the familiy numbers 

        Args:
            args (tuple[dict]): The syntenies of the target genes and assigned families.

        Returns:
            dict: The syntenies with the adapted families assigned to the flanking genes.
        """  
        print('Adapt families')
        syntenies, curr_numbers = args

        return syntenies
    
        for k in syntenies.keys():
            for i, _ in enumerate(syntenies[k]['flanking_genes']['ncbi_codes']):
                
                protein_family = syntenies[k]['flanking_genes']['families'][i]

                if 0 <= protein_family <= max(self.cluster_list):
                    # the pseudogenes have -1, hence we restrict adaption to positiv number
                    # also excluded are the target genes with a higher number than the highest cluster number
                    cluster_range = range(min(self.cluster_list), max(self.cluster_list)+1)
                    # this gives problems with -1. It sorted, hence -1 is at the beginning
                    # and now indexes are + 1, thats why -1 was removed
                    if protein_family != cluster_range[curr_numbers.index(protein_family)]:
                        protein_family = cluster_range[curr_numbers.index(protein_family)]
                        syntenies[k]['flanking_genes']['families'][i] = protein_family

                # add the family of the target gene to the syntenies: 'target_family'
                syntenies[k]['target_family'] = syntenies[k]['flanking_genes']['families'][
                    syntenies[k]['flanking_genes']['ncbi_codes'].index(syntenies[k]['assembly_id']['ERZ-cds'])]

        return syntenies
