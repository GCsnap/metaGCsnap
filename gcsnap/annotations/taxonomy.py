import os
import subprocess
import gzip
from pathlib import Path
import pandas as pd
#import pytaxonkit
import csv
import json
from collections import defaultdict
import re
from infomap import Infomap

import gcsnap.providers.MGnify.helpers as hp
from gcsnap.annotations.taxtree import get_dist
from gcsnap.genomic_context import GenomicContext
from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration


class SourMashBinning:

    """
    SourMash binning of metagenomic contigs object.
    """

    def __init__(self, gc: GenomicContext, config: Configuration, similarity_threshold: float = 0.97):
        
        """
        Initialize the SourMashBinning object.
        :param dataset: A dataset object containing paths and parameters for SourMash-distance based binning.
        :param gc: A genomic context object containing paths and parameters for the genomic context.
        :param config: A configuration object containing paths and parameters for the configuration.
        :param similarity_threshold: The similarity threshold for the binning. 0.97 is retained as a species-level threshold.
        """

        self.console = RichConsole()

        self.sourmash_executable = config.arguments['sourmash_executable_path']['value']

        self.gc = gc

        self.contigs_fasta = [ self.gc.syntenies[target]['assembly_metadata']['dna_file'] for target in gc.syntenies.keys() ]
        self.contigs_fasta = [ f for f in self.contigs_fasta if os.path.exists(f)]
        self.gere_cols = { gc.syntenies[target]['assembly_metadata']['dna_file']: gc.syntenies[target]['assembly_metadata']['genomic_region'] for target in gc.syntenies.keys() } 
        
        self.binning_out_dir = self.gc.out_label / 'binning'
        self.binning_out_dir.mkdir(parents=True, exist_ok=True)
        self.signatures = self.binning_out_dir / 'contigs.sig'
        
        # Path to the query FASTA file (e.g., metagenomic contigs)

        # Path to the SourMash output and report file
        self.labels_file = self.binning_out_dir / "labels.txt"
        self.similarity_matrix_file = self.binning_out_dir / "similarity.csv"
        self.containment_matrix_file = self.binning_out_dir / "cont.npy"
        self.bins_distance_matrix_file = self.binning_out_dir / "bins_distance_matrix.csv.gz"
        self.targets_distance_matrix_file = self.binning_out_dir / "distance_matrix.csv.gz"
        
        self.gc.metagenomic_bins_distance_file = self.targets_distance_matrix_file

        self.bins_file =self.binning_out_dir / "bins.json"
        self.tree_file = self.binning_out_dir / "tree.json"
        self.feature = 'sourmash_similarity'
        self.similarity_threshold = similarity_threshold
        
        # Number of threads to use for SourMash
        self.threads = config.arguments['n_cpu']['value']

        region_to_targets = defaultdict(list)
        for target, data in gc.syntenies.items():
            region = data['assembly_metadata']['genomic_region']
            region_to_targets[region].append(target)

        self.region_to_targets = dict(region_to_targets)
        

    def get_distance_matrix(self):
        """
        Get the binning distance matrix.
        """

        try:
            return self.distance_matrix
        except ValueError:
            return pd.read_csv( self.targets_distance_matrix_file, compression='gzip', index_col='contig')

    def get_tree(self):

        return self.tree
        
    def _run_command(self, command: list):
        """
        Helper method to run a shell command and handle errors.
        """

        try:
            subprocess.run(command, capture_output=True, text=True)
        except FileNotFoundError:
            if "sourmash" not in command:
                raise

            adjusted_command = [
                str(self.sourmash_executable) if arg == "sourmash" else arg
                for arg in command
            ]

            try:
                subprocess.run(adjusted_command, capture_output=True, text=True)
            except FileNotFoundError:
                raise RuntimeError(
                    "Sourmash executable not found. Please install Sourmash or update "
                    "the 'sourmash-executable-path' setting in the configuration."
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"Command '{' '.join(adjusted_command)}' failed with error: {e}"
                )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Command '{' '.join(command)}' failed with error: {e}")

    def _create_signatures(self):

        # sourmash sketch dna -p k=31,scaled=1000 contigs/contigs/MGYC*.fna.gz -o sourmash/contigs.sig

        command = [ "sourmash", "sketch", "dna", "-p", "k=31,scaled=1000"] + self.contigs_fasta + ["-o", self.signatures ]
        command = [str(s) for s in command]

        self._run_command(command)

    def _compute_similarity_matrix(self):

        # sourmash compare --containment sourmash/contigs.sig -o sourmash/cont.npy --csv sourmash/similarity.csv  --labels-save sourmash/labels.txt
        command = [ "sourmash", "compare", "--containment", self.signatures, "-o", self.containment_matrix_file, "--csv", self.similarity_matrix_file, "--labels-save", self.labels_file ]
        command = [str(s) for s in command]

        self._run_command(command)

    def _compute_target_distance_matrix(self):

        # 1. Prepare the exact order of old contigs and their new target names
        source_contigs = []
        new_targets = []

        for contig in self.distance_matrix.index:
            # Get the targets from the dict (fallback to the contig name if missing)
            targets = self.region_to_targets.get(contig, [contig])
            for target in targets:
                source_contigs.append(contig)
                new_targets.append(target)

        # 2. Expand the matrix using .loc (this automatically duplicates rows and columns!)
        target_dist_matrix = self.distance_matrix.loc[source_contigs, source_contigs].copy()

        # 3. Replace the old contig names with the new target names
        target_dist_matrix.index = new_targets
        target_dist_matrix.columns = new_targets
        target_dist_matrix.index = target_dist_matrix.index.rename('target')
        target_dist_matrix = 0.5*(target_dist_matrix+target_dist_matrix.T)
        target_dist_matrix.to_csv(self.targets_distance_matrix_file,index=True,compression='gzip')

    def _compute_distance_matrix(self):

        labels = pd.read_csv(self.labels_file, usecols=['sort_order','filename'],index_col='sort_order')
        
        order_to_contig_dict = labels['filename'].to_dict()
        order_to_contig_dict = { k: self.gere_cols[v] for k,v in order_to_contig_dict.items()}
        
        self.sim_matrix = pd.read_csv(self.similarity_matrix_file)
        self.sim_matrix.index += 1

        self.sim_matrix = self.sim_matrix.rename(columns=self.gere_cols, index=order_to_contig_dict)
        self.distance_matrix = 1 - self.sim_matrix
        self.distance_matrix.index = self.distance_matrix.index.rename('contig')
        self.distance_matrix.to_csv(self.bins_distance_matrix_file, compression='gzip')
    
    def _find_species_bins(self):


        if os.path.isfile(self.bins_file):
            print(f'Read metagenomic bins from file.')
            self.bins = pd.read_json(self.bins_file, typ='series').to_dict()
        else:
            sim = 1-pd.read_csv(self.bins_distance_matrix_file, compression='gzip', index_col='contig')
            sim = sim.reset_index().rename(columns={'contig': 'query'})

            S = sim.melt( id_vars=['query'],  var_name='target', value_name=self.feature )
            S=S[S[self.feature]>self.similarity_threshold]

            ids = list(set(S['query'].values).union(S['target'].values))

            only_self_links = (S['query'] == S['target']).all()
            
            if only_self_links:
                self.bins = {c:i+1 for i,c in enumerate(ids)}
            else:
                id_to_node = {p:i for i, p in enumerate(ids)}
                node_to_id = {i:p for p, i in id_to_node.items()}

                im = Infomap("--flow-model undirected --seed 42 --silent")

                # Add edges to Infomap
                for _, row in S[['query', 'target', self.feature]].iterrows():
                    q,t = id_to_node[row["query"]], id_to_node[row["target"]]
                    im.add_link(q, t, row[self.feature])

                # Run the algorithm
                im.run()

                # Extract and print the communities
                self.bins = {node_to_id[node.node_id]: str(node.module_id) for node in im.nodes}
            
                with open(self.bins_file, 'w') as f:
                    json.dump(self.bins, f, indent=2)
    
    def _release_tree(self):
        
        if os.path.exists(self.bins_file):

            with open(self.bins_file, 'r') as f:
                self.bins = json.load(f)

        self.tree = {}

        for contig,bin in self.bins.items():
            self.tree[f'bin {bin}']={'target_members':[],'cds_codes':[]}
            for cds in self.region_to_targets[contig]:
                self.tree[f'bin {bin}']['target_members'].append(cds)
                self.tree[f'bin {bin}']['cds_codes'].append(cds)
        
        self.tree = {'root':self.tree}
        with open(self.tree_file, 'w') as f:
            json.dump(self.tree, f, indent=2)

    def get_sim_matrix(self):
        """
        Get the similarity matrix.
        """
        return self.sim_matrix
    
    def _update_taxonomy(self):

        
        for k,v in self.gc.syntenies.items():
            
            contig_id = v['assembly_metadata']['genomic_region']
            v['bin_type'] = self.bins.get(contig_id,0)
            
            #if v['taxonomy']['taxon_name'] == None
            #    v['taxonomy']['taxon_name']=f'metagenomic bin {self.bins.get(contig_id,"unclassified")}'

    def run(self):
        
        with self.console.status('Launching SourMash binning of genomic regions'):
            pass

        if not os.path.exists(self.signatures):
            with self.console.status('Creating signatures'):
                self._create_signatures()
        else:
            with self.console.status('Signatures file already exists.'):
                pass

        if not os.path.exists(self.similarity_matrix_file):
            with self.console.status('Computing similarity matrix'):
                self._compute_similarity_matrix()
        else:
            with self.console.status('Similarity matrix already exists.'):
                pass

        if not os.path.exists(self.bins_distance_matrix_file):
            with self.console.status('Computing distance matrix'):
                self._compute_distance_matrix()
        else:
            with self.console.status('Distance matrix already exists.'):
                pass

        if not os.path.exists(self.targets_distance_matrix_file):
            with self.console.status('Preparing target level distance matrix'):
                self._compute_target_distance_matrix()
        else:
            with self.console.status('Target level distance matrix already exists.'):
                pass

        with self.console.status('Finding species level bins'):
            self._find_species_bins()
        
        self._update_taxonomy()

        self._release_tree()