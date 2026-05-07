import os
import pacmap

import json
import numpy as np
import numpy as np
import pandas as pd

from typing import Union
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture
from scipy.spatial import distance
from scipy.cluster import hierarchy

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext

import logging
logger = logging.getLogger('iteration')

class Operons:
    """ 
    Methods and attributes to cluster operons.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cores (int): The number of CPU cores to use.
        advanced (bool): The advanced operon clustering option.
        max_freq (float): The maximum frequency of a family to be considered.
        min_freq (float): The minimum frequency of a family to be considered.
        families (dict): The dictionary with the families assigned to the flanking genes.
        syntenies (dict): The dictionary with the syntenies of the target genes.
        out_label (str): The label of the output.
        console (RichConsole): The RichConsole object to print messages.
        sorted_targets (list): The list of sorted targets.
        operon_clusters (list): The list of operon clusters.
        clean_coordinates (np.ndarray): The clean coordinates.
        all_coordinates (np.ndarray): The all coordinates.
    """

    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str=None):
        """
        Initialize the Operons object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_label (str): The label of the output.
        """        
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.advanced = config.arguments['operon_cluster_advanced']['value']
        self.max_freq = config.arguments['max_family_freq']['value']
        self.min_freq = config.arguments['min_family_freq']['value']

        # get families and syntenies from the synties object
        self.families = gc.get_families()
        self.syntenies = gc.get_syntenies()
        self.out_label = out_label

        self.console = RichConsole() 

        # data structures created within the class
        # set default values if nothing has to be done in cases only 
        # one target is present in syntenies
        self.sorted_targets = list(self.syntenies.keys())
        self.operon_cluster = [1] * len(self.sorted_targets)      

        if out_label is None:
            self.out_label = config.arguments['out_label']['value']
        else:
            self.out_label = out_label

        self.out_dir = os.path.join(self.out_label, 'operons')
        if not os.path.isdir(self.out_dir):
            os.mkdir(self.out_dir)

        self.distance_matrix_file = os.path.join(self.out_dir, 'distance_matrix.csv.gz')
        self.operons_file = os.path.join(self.out_dir, 'operons.json')

        gc.distance_matrix_file = self.distance_matrix_file
        gc.operons_file = self.operons_file

    def get_operons(self) -> dict:
        """
        Getter for the operon_clusters attribute.

        Returns:
            dict: The syntenies dictionary with the operon clusters.
        """        
        return self.syntenies        

    def run(self,skip=False) -> None:
        """
        Run the clustering of operons. Either standard or advanced clustering.
        """     

        # if syntenies has only one target, the default values are kept
        if len(self.syntenies) > 1:
            if self.advanced:

                with self.console.status('Running advanced operon clustering with PaCMAP'):
                    self.run_advanced()
            else:
                if skip:
                    with self.console.status('Skipping operon clustering'):
                        self.run_skip()
                else:
                    with self.console.status('Running standard operon clustering'):
                        self.run_standard()

        # save results
        with open(self.operons_file, 'w') as f:
            json.dump(self.operon_cluster, f)

        self.distance_matrix.index = self.distance_matrix.index.rename('target')
        self.distance_matrix.to_csv(self.distance_matrix_file,compression='gzip')
        
        # update syntenies with the new operon clusters
        self.update_syntenies()       

    def run_skip(self) -> None:
        """
        Skip operon clustering if you want to display all the genomic contexts in one figure:
        Return a dummy operon clustering.
        """        
        
        trivial_D = np.ones((len(self.sorted_targets), len(self.sorted_targets)))
        np.fill_diagonal(trivial_D, 0)
        self.distance_matrix = pd.DataFrame(trivial_D, index=self.sorted_targets, columns=self.sorted_targets)

        self.operon_cluster  = {target: i+1 for i,target in enumerate(self.sorted_targets)}
        
    def run_advanced(self) -> None:
        """
        Run the advanced clustering of operons with PaCMAP:
            - Find the operon clusters with PaCMAP.
            - Find all PaCMAP coordinates.
            - Update the syntenies with the new operon clusters.
        """
        # get the clusters by excluding the most common families
        res = self.find_operon_clusters_with_PaCMAP(clean = True, coordinates_only = False)
        clean_coords_arr, clusters_arr, self.sorted_targets = res
        # and now all PaCMAP coordinates by using all families.
        # This will be later used for sorting the dendrogram
        res = self.find_operon_clusters_with_PaCMAP(clean = False, coordinates_only = True)
        all_coords_arr, self.sorted_targets = res

        # convert to target-keyed dicts — consistent with standard path
        self.operon_cluster    = {t: int(clusters_arr[i])   for i, t in enumerate(self.sorted_targets)}
        self.clean_coordinates = {t: clean_coords_arr[i]    for i, t in enumerate(self.sorted_targets)}
        self.all_coordinates   = {t: all_coords_arr[i]      for i, t in enumerate(self.sorted_targets)}

        # build distance matrix from pairwise Euclidean distances in 2D PaCMAP space
        # consistent with what run() expects (.index rename, .to_csv)
        n = len(self.sorted_targets)
        D = np.zeros((n, n))
        for i, ti in enumerate(self.sorted_targets):
            for j, tj in enumerate(self.sorted_targets):
                if i < j:
                    d = np.linalg.norm(clean_coords_arr[i] - clean_coords_arr[j])
                    D[i, j] = d
                    D[j, i] = d
        self.distance_matrix = pd.DataFrame(D, index=self.sorted_targets, columns=self.sorted_targets)

    def find_operon_clusters_with_PaCMAP(self, clean: bool = True, 
                                         coordinates_only: bool = True,
                                         max_eps_cost: float = 1.3, 
                                         eps_step: float = 0.02
                                         ) -> Union[tuple[np.ndarray,list],
                                                    tuple[np.ndarray,np.ndarray,list]]:
        """
        Find the operon clusters with PaCMAP.

        Args:
            clean (bool, optional): Flag to clean the families. Defaults to True.
            coordinates_only (bool, optional): Flag to return only the coordinates. Defaults to True.
            max_eps_cost (float, optional): Maximum cost for the epsilon. Defaults to 1.3.
            eps_step (float, optional): Step for the epsilon improvement. Defaults to 0.02.

        Returns:
            Union[tuple[np.ndarray,list], tuple[np.ndarray,np.ndarray,list]]: The coordinates and clusters.
        """        
        
        presence_matrix, sorted_targets, sorted_families = self.get_family_presence_matrix(clean)
        paCMAP_embedding = pacmap.PaCMAP(n_components = 2)
        paCMAP_coordinat = paCMAP_embedding.fit_transform(presence_matrix)          

        if coordinates_only:
            return paCMAP_coordinat, sorted_targets
        else:
            # embed into n-D pacMAP space
            n_dims = len(sorted_families)
            if n_dims < 2:
                n_dims = 2
            
            paCMAP_embedding = pacmap.PaCMAP(n_components = n_dims)
            paCMAP_N_coordinat = paCMAP_embedding.fit_transform(presence_matrix)

		# find clusters in the paCMAP space
		# do this by selecting the best eps based on the number of clusters it creates compared to the number of operons 
		# that are not assigned a clusters (i.e., a given maximum cost)        
        eps = self.calculate_start_eps(paCMAP_N_coordinat)    

        n_clusters = [0]
        n_singletons = [0]
        cost = 0
        logger.info('Refining EPS for DBSCAN clustering based on Gaussian Mixture Model')
        while cost <= max_eps_cost and eps > 0:
            eps = eps - eps_step

            model = DBSCAN(eps = eps)
            model.fit(paCMAP_coordinat)
            clusters = model.fit_predict(paCMAP_coordinat)

            n = len(set(clusters))
            delta_n_clusters = n - n_clusters[-1]
            delta_singletons = list(clusters).count(-1) - n_singletons[-1]
            if delta_n_clusters > 0:
                cost = delta_singletons/delta_n_clusters
            else:
                cost = 0

        logger.info('Calculated EPS: {}'.format(eps))
        logger.info('Cost: {}'.format(cost))  
        logger.info('Number of clusters: {}'.format(n))              

        return paCMAP_coordinat, clusters, sorted_targets                

    def get_family_presence_matrix(self, clean: bool = True) -> tuple[np.ndarray, list, list]:
        """
        Get the family presence matrix.

        Args:
            clean (bool, optional): Flag to clean the families. Defaults to True.

        Returns:
            tuple[np.ndarray, list, list]: The presence matrix, sorted targets, and sorted families.
        """        
        sorted_targets = sorted(list(self.syntenies.keys()))
        sorted_families   = [i for i in sorted(list(self.families.keys())) if (i >= 0)]

        # select only the protein families that are not very frequenct but also not very rare
        if clean and len(sorted_families) > 10:
            families_frequency = [len(self.families[family]['members']) for family in sorted_families]
            families_frequency = [i*100/len(self.syntenies) for i in families_frequency]
            sorted_families	= [family for i, family in enumerate(sorted_families) 
                               if families_frequency[i] <= self.max_freq 
                               and families_frequency[i] >= self.min_freq]

        presence_matrix = [[0 for i in sorted_families] for i in sorted_targets]
        for i, target_i in enumerate(sorted_targets):
            operon_i = self.syntenies[target_i]['flanking_genes']['families']
            for family in operon_i:
                if family in sorted_families:
                    presence_matrix[i][sorted_families.index(family)] += 1

        return np.array(presence_matrix), sorted_targets, sorted_families   

    def calculate_start_eps(self, coordinates: list) -> float:
        """
        Calculate the start epsilon.

        Args:
            coordinates (list): The coordinates.

        Returns:
            float: The start epsilon.
        """        
        distances = []
        for i, vector_i in enumerate(coordinates):
            for j, vector_j in enumerate(coordinates):
                if j>i:
                    dist = np.linalg.norm(np.array(vector_i) - np.array(vector_j))
                    distances.append(dist)

        distances = np.array(distances)
        mixture = GaussianMixture(n_components=2).fit(distances.reshape(-1,1))
        means = mixture.means_.flatten()
        sds = np.sqrt(mixture.covariances_).flatten()

        mean = min(means)
        sd = sds[list(means).index(mean)]

        eps = mean - sd        
        return round(eps, 2)    

    def run_standard(self) -> None:
        """
        Run the standard clustering of operons:
            - Compute the operon distance matrix.
            - Find the clusters in the distance matrix.
        """        
        self.distance_matrix, self.sorted_targets = self.compute_operon_distance_matrix()
        self.operon_cluster = self.find_clusters_in_distance_matrix(self.distance_matrix, t = 0.2)

    def compute_operon_distance_matrix(self) -> tuple[np.ndarray,list]:
        """
        Compute the operon distance matrix.

        Returns:
            tuple[np.ndarray,list]: The distance matrix and sorted targets.
        """  
        sorted_targets  = sorted(list(self.syntenies.keys()))      
        D = np.zeros((len(sorted_targets), len(sorted_targets)))
        
        for i, operon_i in enumerate(sorted_targets):
            for j, operon_j in enumerate(sorted_targets):
                if i < j:
                    vector_i = self.syntenies[operon_i]['flanking_genes']['families']
                    reference_family_i = self.syntenies[operon_i]['target_family']
                    
                    vector_j = self.syntenies[operon_j]['flanking_genes']['families']
                    reference_family_j = self.syntenies[operon_j]['target_family']

                    if len(vector_i) >= len(vector_j):
                        reference_vector = vector_i
                        curr_vector = vector_j
                        index_reference_family_reference_vector = (
                            list(reference_vector).index(reference_family_i))
                        index_reference_family_curr_vector = (
                            list(curr_vector).index(reference_family_j))
                    else:
                        reference_vector = vector_j
                        curr_vector = vector_i
                        index_reference_family_reference_vector = (
                            list(reference_vector).index(reference_family_j))
                        index_reference_family_curr_vector = (
                            list(curr_vector).index(reference_family_i))

                    number_fillings = 0
                    if len(curr_vector) < len(reference_vector):
                        if (index_reference_family_curr_vector < 
                            index_reference_family_reference_vector):
                            for a in range(index_reference_family_reference_vector - 
                                           index_reference_family_curr_vector):
                                curr_vector = np.insert(curr_vector, 0, -1)
                                number_fillings += 1
                        
                        if len(curr_vector) < len(reference_vector):
                            for b in range(len(reference_vector) - len(curr_vector)):
                                curr_vector = np.append(curr_vector, -1)
                                number_fillings += 1

                    families_present = set(np.concatenate([reference_vector,curr_vector]))
                    overlapping_families = set(reference_vector).intersection(set(curr_vector))
                    dist = 1-len(overlapping_families)/len(families_present)

                    D[i,j] = dist
                    D[j,i] = dist

        distance_matrix = pd.DataFrame(D, index=sorted_targets, columns=sorted_targets)
        return distance_matrix, sorted_targets       

    def find_clusters_in_distance_matrix(self, distance_matrix: pd.DataFrame, t: float = 0.0) -> dict: 
        """
        Find the clusters in the distance matrix.

        Args:
            distance_matrix (pd.DataFrame): The distance matrix with targets as index/columns.
            t (float, optional): The threshold for the clustering. Defaults to 0.2.

        Returns:
            dict: Dictionary mapping target IDs to their cluster assignments.
        """        
        # Save the target index before converting to condensed form
        targets = distance_matrix.index
        
        # Convert to condensed distance matrix for scipy
        distance_matrix_condensed = distance.squareform(distance_matrix.values)
        linkage = hierarchy.linkage(distance_matrix_condensed, method = 'single')
        clusters = hierarchy.fcluster(linkage, t, criterion = 'distance')

        # Map targets to their cluster assignments
        return {target: int(cluster) for target, cluster in zip(targets, clusters)}

    def update_syntenies(self) -> None:
        """
        Update the syntenies with the new operon clusters
        """  

        for target,operon_type in self.operon_cluster.items():
            self.syntenies[target]['operon_type'] = operon_type

            if self.advanced:
                self.syntenies[target]['operon_filtered_PaCMAP'] = (
                    list([float(a) for a in self.clean_coordinates[target]]))
                self.syntenies[target]['operon_PaCMAP'] = (
                    list([float(a) for a in self.all_coordinates[target]]))