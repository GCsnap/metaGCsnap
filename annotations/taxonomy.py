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
import providers.MGnify.helpers as hp
from annotations.taxtree import get_dist
from gcsnap.genomic_context import GenomicContext
from gcsnap.configuration import Configuration
from providers.MGnify.dataset import Dataset

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

        self.sourmash_executable = config.arguments['sourmash_executable_path']['value']

        self.gc = gc

        self.contigs_fasta = [ self.gc.syntenies[target]['assembly_metadata']['dna_file'] for target in gc.syntenies.keys() ]
        self.contigs_fasta = [ f for f in self.contigs_fasta if os.path.exists(f)]
        self.gere_cols = { gc.syntenies[target]['assembly_metadata']['dna_file']: gc.syntenies[target]['assembly_metadata']['genomic_region'] for target in gc.syntenies.keys() } 
        
        self.binning_out_dir = self.gc.out_label / 'binning'
        self.binning_out_dir.mkdir(parents=True, exist_ok=True)
        self.output = self.binning_out_dir / 'contigs.sig'
        
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
        
        # Print initialization details
        print(f'''SourMash binning of metagenomic contigs.''')

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

        command = [ "sourmash", "sketch", "dna", "-p", "k=31,scaled=1000"] + self.contigs_fasta + ["-o", self.output ]
        command = [str(s) for s in command]

        self._run_command(command)

    def _compute_similarity_matrix(self):

        # sourmash compare --containment sourmash/contigs.sig -o sourmash/cont.npy --csv sourmash/similarity.csv  --labels-save sourmash/labels.txt
        command = [ "sourmash", "compare", "--containment", self.output, "-o", self.containment_matrix_file, "--csv", self.similarity_matrix_file, "--labels-save", self.labels_file ]
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

        print(f'Finding species level bins with Infomap:')
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
        
        if not os.path.exists(self.bins_file):

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
        
        print(f'''Launchiung SourMash binning of metagenomic contigs''')

        if not os.path.exists(self.output):
            self._create_signatures()
        else:
            print(f"Output file {self.output} already exists.")

        if not os.path.exists(self.similarity_matrix_file):
            print('Computing similarity matrix...')
            self._compute_similarity_matrix()
        else:
            print(f"Output file {self.similarity_matrix_file} already exists.")

        if not os.path.exists(self.bins_distance_matrix_file):
            print('Computing distance matrix...')
            self._compute_distance_matrix()
        else:
            print(f"Output file {self.bins_distance_matrix_file} already exists.")

        if not os.path.exists(self.targets_distance_matrix_file):
            print('Preparing target level distance matrix...')
            self._compute_target_distance_matrix()
        else:
            print(f"Output file {self.targets_distance_matrix_file} already exists.")

        self._find_species_bins()
        
        self._update_taxonomy()

        self._release_tree()

class Kraken2Taxonomy:

    """
    Kraken2 taxonomic assignment of metagenomic contigs object.
    """

    def __init__(self, dataset: Dataset, gc: GenomicContext, config: Configuration):
        
        """
        Initialize the Kraken2Taxonomy object.
        :param dataset: A dataset object containing paths and parameters for Kraken2.
        """

        # Path to the query FASTA file (e.g., metagenomic contigs)
        self.contigs_fasta = dataset.query_contig_fasta

        # Basename of the query FASTA file
        #self.query_basename = dataset.query_basename

        # Path to the Kraken2 reference database
        self.reference_db = dataset.kraken2_db

        # Path to the Kraken2 output and report file
        self.query_fasta = dataset.taxonomic_out_dir / "contigs.fna.gz"
        self.output = dataset.kraken2_output
        self.report = dataset.kraken2_report
        self.taxonomy_json = dataset.taxonomic_out_dir / "taxonomy.json"
        self.distance_matrix_file = dataset.taxonomic_out_dir / "distance_matrix.csv.gz"

        gc.taxonomic_tree_file = self.taxonomy_json
        gc.taxonomic_distance_file = self.distance_matrix_file
        
        # Number of threads to use for Kraken2
        self.threads = config.arguments['n_cpu']['value']

        # taxonomic database
        self.taxonomy_db = dataset.taxonomy_local

        # Group MGYPs and ERZ_cds by ERZ_contig and create a dictionary
        self.erz_to_mgyp_cds = ( dataset.mgyp_metadata
            .groupby('ERZ_contig')[['ERZ_cds_id', 'MGYP']]
            .apply(lambda df: list(df.itertuples(index=False, name=None)))
            .to_dict()
        )
        
        # Print initialization details
        print(f'''Kraken2 taxonomic assignment of metagenomic contigs with GTDB.''')

    def _prepare_query(self):

        with gzip.open(self.query_fasta, "wt") as outfile:
            for fasta_file in self.contigs_fasta:
                try:
                    with gzip.open(fasta_file, "rt") as infile:
                        for line in infile:
                            outfile.write(line)
                except Exception as e:
                    print(f"Error processing {fasta_file}: {e}")

    def _run_command(self, command: list):
        """
        Helper method to run a shell command and handle errors.
        """
        try:
            subprocess.run(command, check=True, text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Command '{' '.join(command)}' failed with error: {e}")

    def _assign_taxonomy(self):
        """
        Perform Kraken2 search between the query assembly and the reference database.
        """
        command = [ "kraken2", "--db", str(self.reference_db),
                    "--output", str(self.output), "--report", str(self.report),
                    "--threads", f'{self.threads}', "--use-names", str(self.query_fasta) ]

        self._run_command(command)

    def _curate_taxonomy(self):
        
        # Map rank codes to prefixes
        rank_prefix = { 'R': 'r__', 'R1': 'd__', 'P': 'p__', 'C': 'c__','O': 'o__', 'F': 'f__', 'G': 'g__', 'S': 's__'}
        rank_name = { 'R': 'root', 'R1': 'domain', 'P': 'phylum', 'C': 'class','O': 'order', 'F': 'family', 'G': 'genus', 'S': 'species'}

        def reconstruct_lineage(taxon_id):

            if taxon_id not in tax_lookup:
                return None

            lineage_names = []
            lineage_ids = []
            current = taxon_id

            while True:
                node = tax_lookup.get(current)
                if node is None:
                    break
                rank = node['rank']
                prefix = rank_prefix.get(rank, 'x__')  # default unknown ranks
                lineage_names.append(f"{prefix}{node['taxid']}")
                lineage_ids.append(str(current))
                if current == node['parent']:  # reached root
                    break
                current = node['parent']

            lineage_names.reverse()
            lineage_ids.reverse()
            return ';'.join(lineage_names), ';'.join(lineage_ids), rank_name[tax_lookup[taxon_id]['rank']]

        K = pd.read_csv(self.output,sep='\t',header=None, names=['classified', 'genomic_region', 'classification', 'length', 'matches'],dtype=str)
        K[['taxon_name', 'taxon_id']] = K['classification'].str.extract(r'^(.*) \(taxid (\d+)\)$')

        unique_taxon_ids = [taxon_id for taxon_id in K['taxon_id'].dropna().unique() if taxon_id != '0']

        tax_file = self.reference_db / 'ktaxonomy.tsv'
        tax_df = pd.read_csv(tax_file, sep='\t', header=None, usecols=[0,2,4,6,8],names=['node','parent','rank','rankid','taxid'],dtype=str)
        del tax_df['rankid']

        tax_lookup = tax_df.set_index('node').to_dict('index')

        results = []

        for tax_id in unique_taxon_ids:
            lineage = reconstruct_lineage(tax_id)
            if lineage:
                lineage_names, lineage_ids, final_rank = lineage
                results.append((tax_id, lineage_names, lineage_ids, final_rank))

        # If needed, convert to dataframe
        lineage_df = pd.DataFrame(results, columns=['taxon_id', 'lineage_names', 'lineage_ids', 'rank'])
        lineage_df = lineage_df.set_index('taxon_id')

        for l in ['lineage_names', 'lineage_ids', 'rank']:

            K[l]=K['taxon_id'].map(lineage_df[l].to_dict())

        cols = ['classified', 'genomic_region','length','taxon_name', 'taxon_id','rank','lineage_names', 'lineage_ids', 'matches']
        K['length'] = K['length'].astype(int)
        K.sort_values(by='length', ascending=False, inplace=True)
        
        K['taxon_name'].replace('unclassified', 'unclassified contig', inplace=True)
        K['taxon_name'].replace('root', 'unclassified contig', inplace=True)
        K['rank'].fillna('root', inplace=True)
        K['lineage_names'].fillna('r__root', inplace=True)
        K['lineage_ids'].fillna('1', inplace=True)
        K['taxon_id'].replace(0,1, inplace=True)

        K.to_csv(self.output, sep='\t', index=False, columns=cols)

    def _prepare_taxonomy_tree(self):

        tree = {}

        try:
            tax = pd.read_csv(self.output, sep='\t', usecols=['genomic_region', 'lineage_names'], dtype=str)
            tax = tax.dropna(subset=['lineage_names'])
            lineage_to_contigs = tax.groupby('lineage_names')['genomic_region'].apply(list)
        except Exception as e:
            print(f"An error occurred during Kraken2 output reading: {e}. Inspect {self.output}")
            with open(self.taxonomy_json, "w") as json_file:
                json.dump(tree, json_file, indent=2)

        for lineage_str, contig_list in lineage_to_contigs.items():
        
            # --- Process Contig Data ---
            K_values = []
            for c in contig_list:
                if c in self.erz_to_mgyp_cds:
                    K_values.append(self.erz_to_mgyp_cds[c])
                # else:
                    # print(f"Warning: Contig {c} not found in erz_to_mgyp_cds dict.")

            targets = [item[0] for item_list in K_values for item in item_list]

            current_level_dict = tree
            lineage_map = {}
            for part in lineage_str.split(';'):
                if '__' in part:
                    prefix, name = part.split('__', 1)
                    if prefix in hp.tax_ranks_dict:
                        rank_name = hp.tax_ranks_dict[prefix]
                        lineage_map[rank_name] = name

            last_known_level_dict = None

            for rank in hp.tax_ranks:
                if rank in lineage_map:
                    taxon_name = lineage_map[rank]
                    current_level_dict = current_level_dict.setdefault(taxon_name, {})
                    last_known_level_dict = current_level_dict
                else:
                    break
            
            # --- FINAL STEP (MODIFIED) ---
            # Instead of appending the erz_cds dict, merge its lists directly
            # into the 'last_known_level_dict'
            if last_known_level_dict is not None:
                # setdefault will create an empty list [] if the key doesn't exist,
                # then .extend() adds all items from the 'targets' list.
                last_known_level_dict.setdefault('target_members', []).extend(targets)
                last_known_level_dict.setdefault('ncbi_codes', []).extend(targets) # Using 'targets' as per your logic
            else:
                # Handle unclassified
                unclassified_dict = tree.setdefault('unclassified', {})
                unclassified_dict.setdefault('target_members', []).extend(targets)
                unclassified_dict.setdefault('ncbi_codes', []).extend(targets)

        self.tree = tree
        # Save tree as JSON to self.taxonomy_json
        with open(self.taxonomy_json, "w") as json_file:
            json.dump(self.tree, json_file, indent=2)

    def _compute_distance_matrix(self):

        tree = self.get_tree()
        self.distance_matrix = get_dist(tree)
        self.distance_matrix.index = self.distance_matrix.index.rename('target')
        self.distance_matrix.to_csv(self.distance_matrix_file, compression='gzip')

    def get_tree(self):
        """
        Get the taxonomy tree.
        """
        return self.tree

    def get_distance_matrix(self):
        """
        Get the taxonomy distance matrix.
        """
        return self.distance_matrix

    def run(self):
        """
        Full pipeline: create run taxonomic assignment.

        :param output_file: Path to the final output file.
        """
        
        print(f'''Launchiung Kraken2 taxonomic assignment of metagenomic contigs''')

        # Step 1: Perform search
        if not os.path.exists(self.output):
            self._prepare_query()
            self._assign_taxonomy()

        # Step 2: Curate taxonomy
        if os.path.exists(self.output):
            existing_data = pd.read_csv(self.output, sep='\t')
            
            if 'taxon_name' in existing_data.columns:
                print(f"taxon_name column already exists in {self.output}. Skipping taxonomy unfolding.")
            else:
                self._curate_taxonomy()

        # get taxonomy tree
        if not os.path.exists(self.taxonomy_json):
            #self._prepare_taxonomy_json()
            self._prepare_taxonomy_tree()
        else:
            print('Read taxonomy from json file')
            with open(self.taxonomy_json, 'r') as f:
                self.tree = json.load(f)

        # get taxonomy distance matrix
        if not os.path.exists(self.distance_matrix_file):
            self._compute_distance_matrix()