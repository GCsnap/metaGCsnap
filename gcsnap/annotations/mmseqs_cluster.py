import os
import subprocess
import numpy as np
# pip install scipy
from scipy.cluster import hierarchy
from scipy.spatial import distance
from infomap import Infomap
import pandas as pd
import json
from collections import Counter

from gcsnap.consts import MMSeqsParams
from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.rich_console import RichConsole

class MMseqsCluster:
    """ 
    Methods and attributes to cluster flanking genes using MMseqs2.
    Needs MMseqs2 installed or the path to the executable set either in the config.yaml 
    or as a CLI argument.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cores (int): The number of CPU cores to use.
        max_evalue (float): The maximum e-value for the search.
        min_coverage (float): The minimum coverage for the search.
        num_iterations (int): The number of iterations for the search.
        mmseqs_executable (str): The path to the MMseqs2 executable.
        default_base (int): The default base for the distance matrix.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        out_dir (str): The path to store the output of MMseqs.
        sensitivity (float): The sensitivity for the search.
        console (RichConsole): The RichConsole object to print messages.
    """

    def __init__(self, config: Configuration, gc: GenomicContext, out_dir: str, engine : str = 'infomap', feature : str = 'pident', sensitivity : float = 6.0, f_community : float = 0.01):
        """
        Initialize the MMseqsCluster object

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_dir (str): The path to store the output.
        """        
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.max_evalue = config.arguments['max_evalue']['value']
        self.min_coverage = config.arguments['min_coverage']['value']
        self.num_iterations = config.arguments['num_iterations']['value']
        self.mmseqs_executable = r'{}'.format(config.arguments['mmseqs_executable_path']['value'])
        self.default_base = 10

        self.mmseqs2_columns = MMSeqsParams.tsv_columns
        # set arguments
        self.gc = gc
        self.n_gcs = len(self.gc.curr_targets)
        # if mmseqs temporary folder is not set, use the output folder
        self.out_dir = out_dir
        print(f'Output directory: {self.out_dir}')
        # check if existing
        if not os.path.isdir(self.out_dir):
            os.mkdir('./'+self.out_dir)            

        self.communities_file = os.path.join(self.out_dir, 'infomap_communities.json')
        self.mmseqs_results = os.path.join(self.out_dir,'flanking_sequences.mmseqs')

        self.sensitivity = sensitivity
        self.feature = feature
        self.engine = engine
        self.f_community = f_community
        
        self.min_size = f_community * self.n_gcs

        if self.min_size < 1:
            self.min_size = 1
        
        self.fasta_file = self.gc.write_to_fasta('flanking_sequences.fasta', 
                                            self.out_dir, exclude_pseudogenes = False)  
        self.console = RichConsole()

    def run(self) -> None:
        """
        Run the clustering of flanking genes using MMseqs2 and Scipy:
            - Prepare data for MMseqs
            - Run MMseqs
            - Extract distance matrix
            - Find clusters with Scipy
            - Mask singleton clusters
        """        
        with self.console.status('Prepare data for MMseqs'):
            self.cluster_order = self.gc.get_fasta_order(exclude_pseudogenes = False) 

        with self.console.status('Running MMseqs'):            
            self.run_mmseqs()           
            
        with self.console.status('Find clusters'):
            if self.engine == 'infomap':
                self.find_communities()
            else:
                self.extract_distance_matrix()
                self.find_clusters()

            #self.mask_singleton_clusters()

    def get_distance_matrix(self) -> np.array:
        """
        Getter for the distance_matrix attribute.

        Returns:
            np.array: The distance matrix.
        """        
        return self.distance_matrix

    def get_clusters_list(self) -> list[int]:
        """
        Getter for the cluster_list attribute.

        Returns:
            list[int]: The list of clusters.
        """        
        return self.cluster_list      

    def get_cluster_order(self) -> list[str]:
        """
        Getter for the cluster_order attribute.

        Returns:
            list[str]: The order of the clusters.
        """        
        return self.cluster_order       

    def run_mmseqs(self) -> None:
        """
        Run MMseqs to cluster flanking genes.

        Raises:
            FileNotFoundError: If MMseqs is not installed or the path to executable is wrongly set.
        """            
        
        if not os.path.isfile(self.mmseqs_results):
            print(f'Running MMseqs with sensitivity {self.sensitivity}, output file: {self.mmseqs_results}')
            try:
                _, stderr = self.mmseqs_command('mmseqs')
                if len(stderr) > 0:
                    raise FileNotFoundError
            except FileNotFoundError:
                try:
                    _, stderr = self.mmseqs_command('mmseqs')
                    if len(stderr) > 0:
                        raise FileNotFoundError
                except FileNotFoundError:
                    try:
                        _, stderr = self.mmseqs_command(self.mmseqs_executable)
                    except:
                        self.console.print_error('No MMseqs installation was found') 
                        self.console.print_hint('Please install MMseqs or add the path to the executable to config.yaml.')
                        self.console.stop_execution()         
        else:
            print(f'MMseqs results already exist: {self.mmseqs_results}')

    def mmseqs_command(self, mmseqs: str) -> tuple:
        """
        Run MMseqs command to execute.

        Args:
            mmseqs (str): Either 'mmseqs' or if not installed, the path to the MMseqs executable.

        Returns:
            tuple: The stdout and stderr of the MMseqs command.
        """        

        format = ','.join(self.mmseqs2_columns)
        # returns stdout,stderr
        self.command = [mmseqs, 
                'easy-search', 
                self.fasta_file, 
                self.fasta_file, 
                self.mmseqs_results, 
                self.out_dir, 
                '-e', str(self.max_evalue), 
                '-s', str(self.sensitivity),
                '-c', str(self.min_coverage),
                '--num-iterations', str(self.num_iterations),
                '--threads', str(self.cores),
                '--format-output', format]

        print(' '.join(self.command))
        result = subprocess.run(self.command, capture_output=True, text=True)        
        return result.stdout, result.stderr       
    
    def extract_distance_matrix(self) -> None:
        """
        Extract the distance matrix from the MMseqs results.
        """        
        # crate base distance matrix
        distance_matrix = [[self.default_base if i!=j else 0 for i in self.cluster_order] 
                        for j in self.cluster_order]
        queries_labels = {query: i for i, query in enumerate(self.cluster_order)}

        # read mmseqs results
        with open(self.mmseqs_results, 'r') as f:
            mmseqs_records = f.readlines()

        for hsp in mmseqs_records:
            hsp = hsp.split()
            if len(hsp) > 0:
                query = hsp[0].split('|')[0]
                query_index = queries_labels[query]
                target = hsp[1].split('|')[0]
                if target != query:
                    target_index = queries_labels[target]
                    distance_matrix[query_index][target_index] = 0
                    distance_matrix[target_index][query_index] = 0

        self.distance_matrix = np.array(distance_matrix)      

    def find_clusters(self, t: int = 0) -> None:
        """
        Find clusters using the distance matrix with Scipy hierarchical clustering.

        Args:
            t (int, optional): The threshold for the clustering. Defaults to 0.
        """        
        distance_matrix = distance.squareform(self.distance_matrix)
        linkage = hierarchy.linkage(distance_matrix, method = 'single')
        clusters = hierarchy.fcluster(linkage, t, criterion = 'distance')
        self.cluster_list = [int(i) for i in clusters]

    def filter_small_communities(self,node_communities):
        """
        Re-assigns nodes from communities smaller than min_size to community 0.
        
        Args:
            node_communities (dict): Original {node_id: community_id} mapping.
            min_size (int): The minimum size for a community to be kept.
            
        Returns:
            dict: A new {node_id: community_id} mapping with small 
                communities re-assigned to 0.
        """ 
        
        # 1. Count the size of each community
        community_counts = Counter(node_communities.values())

        # 2. Identify all communities smaller than min_size
        small_communities = { community  for community, count in community_counts.items() if count <= self.min_size }

        # 3. Create the new dictionary with re-assigned nodes
        new_node_communities = { node: 0 if community in small_communities else community for node, community in node_communities.items() }
        
        return new_node_communities

    def find_communities(self) -> None:
        
        print(f'Finding communities with Infomap')
        if os.path.isfile(self.communities_file):
            print(f'Read communities from file: {self.communities_file}')
            communities = pd.read_json(self.communities_file, typ='series').to_dict()
            ids = list(communities.keys())
        else:
            S = pd.read_csv( self.mmseqs_results , sep='\t', header=None, names=self.mmseqs2_columns)
            S['query'] = S['query'].str.split('|').str[0]
            S['target'] = S['target'].str.split('|').str[0]
            S = S[['query', 'target', self.feature]]
            
            ids = list(set(S['query'].values).union(S['target'].values))

            id_to_node = {p:i for i, p in enumerate(ids)}
            node_to_id = {i:p for p, i in id_to_node.items()}

            im = Infomap("--directed --flow-model directed --seed 42 --silent")

            # Add edges to Infomap
            for _, row in S[['query', 'target', self.feature]].iterrows():
                q,t = id_to_node[row["query"]], id_to_node[row["target"]]
                im.add_link(q, t, row[self.feature])

            # Run the algorithm
            im.run()

            # Extract and print the communities
            communities = {node_to_id[node.node_id]: node.module_id for node in im.nodes}

            C=len(set([node.module_id for node in im.nodes]))
            unassigned_nodes = list(set(self.cluster_order) - set(ids))
            for i,n in enumerate(unassigned_nodes): communities[n] = C+1+i

            communities = self.filter_small_communities(communities)

            #communities = {str(k):str(v) for k,v in communities.items()}

            print(f'Saving communities file: {self.communities_file}')
            # Save communities to file
            with open(self.communities_file, 'w') as f:
                json.dump(communities, f, indent=4)

        self.cluster_list = [communities[node] for node in self.cluster_order]

    def mask_singleton_clusters(self, mask: int = 0) -> None:
        """
        Mask singleton clusters.

        Args:
            mask (int, optional): The value to mask the singleton clusters. Defaults to 0.
        """        
        new_clusters_list = []

        self.n_nodes = len(self.cluster_list)

        for value in self.cluster_list:

            f = list(self.cluster_list).count(value) / self.n_gcs

            if f < self.f_community:
                new_clusters_list.append(mask)
            else:
                new_clusters_list.append(value)

        self.cluster_list = new_clusters_list
  