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


import metagenomics.helpers as hp
from metagenomics.taxtree import get_dist
from gcsnap.genomic_context import GenomicContext
from gcsnap.configuration import Configuration


class MMseqsSequenceSearch:
    
    """
    MMseqs2 sequence search object.
    """

    def __init__(self, dataset, session, config: Configuration):
        
        """
        Initialize the MMseqsSequenceSearch object.

        :param query_fasta: Path to the query FASTA file (compressed or uncompressed).
        :param reference_db: Path to the precomputed reference database.
        :param output_dir: Path to the output directory.
        :param tmp_dir: Path to the temporary directory (default: /tmp).
        """

        #self.query_fasta = dataset.query_fasta  # Path to the query FASTA file
        
        # Handle compressed or uncompressed FASTA
        self.query_fasta_uncompressed = dataset.seq_search[session]['query_fasta_uncompressed'] #dataset.query_fasta_uncompressed #self._handle_compressed_fasta()
        self.query_basename = dataset.seq_search[session]['query_basename'] # dataset.query_basename #Path(self.query_fasta).name
        
        # mmseq directories
        self.query_db = dataset.seq_search[session]['query_db'] #dataset.query_db #os.path.join(args.out, 'mmseq','query', 'DB')  # Path to the query database
        self.result_db = dataset.seq_search[session]['result_db'] #dataset.result_db #os.path.join(args.out, 'mmseq','result', f'{self.query_basename}')  # Path to the query database
        self.reference_db = dataset.seq_search[session]['reference_db'] #dataset.reference_db #os.path.join(args.MGnify,hp.mgyp_database)  # Path to the precomputed reference database
        
        self.output_file = dataset.seq_search[session]['mmseqs_output_file'] # dataset.mmseqs_output_file #os.path.join(args.out,hp.mgyc_search_out)  # Path to the output directory
        self.hits_fasta = dataset.seq_search[session]['hits_fasta']
        self.hits_ids = dataset.seq_search[session]['hits_ids']
        self.tmp_dir = dataset.tmp_dir #args.tmp_dir  # Path to the temporary directory
        
        self.threads = config.arguments['n_cpu']['value'] #args.threads  # Number of threads to use for MMseqs2
        self.min_seq_id = config.arguments['min_seq_id']['value']
        self.cov_mode = config.arguments['cov_mode']['value']
        self.max_evalue = config.arguments['max_evalue']['value']
        self.coverage = config.arguments['min_coverage']['value']
        self.sensitivity = config.arguments['sensitivity']['value']

        print(f'''Initializing MMseqs2 sequence search with query: {self.query_fasta_uncompressed}, reference: {self.reference_db}, output: {self.output_file}''')

    def _run_command(self, command: list):
        """
        Helper method to run a shell command and handle errors.
        """
        command=[str(s) for s in command]
        try:
            subprocess.run(command, check=True, text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Command '{' '.join(command)}' failed with error: {e}")

    def _create_query_database(self):
        """
        Create an MMseqs2 database from the query FASTA file.
        """
        command = ["mmseqs", "createdb", self.query_fasta_uncompressed, self.query_db]
        self._run_command(command)

    def _search(self):
        """
        Perform MMseqs2 search between the query database and the reference database.
        """
        command = [ "mmseqs", "search", self.query_db, self.reference_db, self.result_db,
                    self.tmp_dir, "--search-type", "3", "-s", self.sensitivity, '--threads',self.threads, 
                    "--min-seq-id", self.min_seq_id, "--cov-mode",self.cov_mode,
                    "-c", self.coverage, "-e", self.max_evalue, '--split', 4, "-a"]

        command = [str(s) for s in command]
        self._run_command(command)

    def _convert_alignment(self):
        """
        Convert the alignment results to a readable format (e.g., tab-separated).
        """

        output_format = "query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid,qorfstart,qorfend,torfstart,torfend"
        command = ["mmseqs", "convertalis", self.query_db, self.reference_db, self.result_db, self.output_file,"--threads",self.threads,"--format-output", output_format]
        self._run_command(command)

    def _release_hits(self):

        mmseqs_cols = ['query','target','evalue','gapopen','pident','fident','nident','qstart','qend','qlen','tstart','tend','tlen','alnlen','raw','bits','cigar','qseq','tseq','qheader','theader','qaln','taln','qframe','tframe','mismatch','qcov','tcov','qset','qsetid','tset','tsetid','qorfstart','qorfend','torfstart','torfend']
        X=pd.read_csv(self.output_file, sep='\t', dtype=str,names=mmseqs_cols)
        
        S=X[['target','tseq']].drop_duplicates()

        with open(self.hits_fasta, "w") as fasta:
            for _, row in S.iterrows():
                fasta.write(f">{row['target']}\n{row['tseq']}\n")

        # Save just the IDs to a separate file
        with open(self.hits_ids, "w") as idf:
            for target_id in S['target']:
                idf.write(f"{target_id}\n")

    def run_pipeline(self):
        """
        Full pipeline: create query database, run search, and convert results.

        :param output_file: Path to the final output file.
        """
        
        print(f'''Launchiung MMseqs2 sequence search pipeline''')

        # Step 1: Create query database
        self._create_query_database()
    
        # Step 2: Perform search
        self._search()

        # Step 3: Convert alignment results
        self._convert_alignment()

        self._release_hits()

class SourMashBinning:

    """
    SourMash binning of metagenomic contigs object.
    """

    def __init__(self, dataset: Dataset, gc: GenomicContext, config: Configuration, similarity_threshold: float = 0.97):
        
        """
        Initialize the SourMashBinning object.
        :param dataset: A dataset object containing paths and parameters for SourMash-distance based binning.
        :param gc: A genomic context object containing paths and parameters for the genomic context.
        :param config: A configuration object containing paths and parameters for the configuration.
        :param similarity_threshold: The similarity threshold for the binning. 0.97 is retained as a species-level threshold.
        """

        self.sourmash_executable = config.arguments['sourmash_executable_path']['value']
        # Path to the query FASTA file (e.g., metagenomic contigs)
        self.contigs_fasta = dataset.query_contig_fasta

        # Path to the SourMash output and report file
        self.output = dataset.sourmash_output
        self.labels_file = dataset.binning_out_dir / "labels.txt"
        self.similarity_matrix_file = dataset.binning_out_dir / "similarity.csv"
        self.containment_matrix_file = dataset.binning_out_dir / "cont.npy"
        self.distance_matrix_file = dataset.binning_out_dir / "distance_matrix.csv.gz"
        gc.metagenomic_bins_distance_file = self.distance_matrix_file
        self.bins_file =dataset.binning_out_dir / "bins.json"
        self.tree_file = dataset.binning_out_dir / "tree.json"
        self.feature = 'sourmash_similarity'
        self.similarity_threshold = similarity_threshold

        gc.binning_distance_file = self.distance_matrix_file
        
        # Number of threads to use for SourMash
        self.threads = config.arguments['n_cpu']['value']

        # Group MGYPs and ERZ_cds by ERZ_contig and create a dictionary
        self.erz_to_mgyp_cds = ( dataset.mgyp_metadata
            .groupby('ERZ_contig')[['ERZ_cds_id', 'MGYP']]
            .apply(lambda df: list(df.itertuples(index=False, name=None)))
            .to_dict()
        )

        self.mgyc_to_erz_contig = {row['MGYC']:row['ERZ_contig'] for _,row in dataset.mgyp_metadata.iterrows()}
        self.cds_to_erz_contig = (
            dataset.mgyp_metadata
            .groupby('ERZ_contig')['ERZ_cds_id']
            .apply(list)
            .to_dict()
        )
        
        # Print initialization details
        print(f'''SourMash binning of metagenomic contigs.''')

    def get_distance_matrix(self):
        """
        Get the binning distance matrix.
        """
        return self.distance_matrix

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

    def _compute_distance_matrix(self):

        sim_matrix = pd.read_csv(self.similarity_matrix_file)

        mgyc_cols = { path: re.search(r'(MGYC\d+)', path).group(1) for path in self.contigs_fasta }
        sim_matrix = sim_matrix.rename(columns=mgyc_cols)
        self.sim_matrix = sim_matrix.rename(columns=self.mgyc_to_erz_contig)
        self.sim_matrix = self.sim_matrix.rename(index={i:c for i,c in enumerate(self.sim_matrix.columns)})
        self.distance_matrix = 1 - self.sim_matrix
        self.distance_matrix.index = self.distance_matrix.index.rename('contig')
        self.distance_matrix.to_csv(self.distance_matrix_file, compression='gzip')
    
    def _find_species_bins(self):

        print(f'Finding species level bins with Infomap: {self.bins_file}')
        if os.path.isfile(self.bins_file):
            print(f'Read metagenomic bins from file: {self.bins_file}')
            self.bins = pd.read_json(self.bins_file, typ='series').to_dict()
        else:
            sim = 1-pd.read_csv(self.distance_matrix_file, compression='gzip', index_col='contig')
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
                self.bins = {node_to_id[node.node_id]: node.module_id for node in im.nodes}
            
                with open(self.bins_file, 'w') as f:
                    json.dump(self.bins, f)
     
    def _release_tree(self):
        
        if not os.path.exists(self.bins_file):

            with open(self.bins_file, 'r') as f:
                self.bins = json.load(f)

        self.tree = {}

        for contig,bin in self.bins.items():
            self.tree[f'bin {bin}']={'target_members':[],'ncbi_codes':[]}
            for cds in self.cds_to_erz_contig[contig]:
                self.tree[f'bin {bin}']['target_members'].append(cds)
                self.tree[f'bin {bin}']['ncbi_codes'].append(cds)
        
        self.tree = {'root':self.tree}
        with open(self.tree_file, 'w') as f:
            json.dump(self.tree, f)

    def get_sim_matrix(self):
        """
        Get the similarity matrix.
        """
        return self.sim_matrix

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

        if not os.path.exists(self.distance_matrix_file):
            print('Computing distance matrix...')
            self._compute_distance_matrix()
        else:
            print(f"Output file {self.distance_matrix_file} already exists.")

        self._find_species_bins()
        
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

        K = pd.read_csv(self.output,sep='\t',header=None, names=['classified', 'ERZ-contig', 'classification', 'length', 'matches'],dtype=str)
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

        cols = ['classified', 'ERZ-contig','length','taxon_name', 'taxon_id','rank','lineage_names', 'lineage_ids', 'matches']
        K['length'] = K['length'].astype(int)
        K.sort_values(by='length', ascending=False, inplace=True)
        
        K['taxon_name'].replace('unclassified', 'unclassified contig', inplace=True)
        K['taxon_name'].replace('root', 'unclassified contig', inplace=True)
        K['rank'].fillna('root', inplace=True)
        K['lineage_names'].fillna('r__root', inplace=True)
        K['lineage_ids'].fillna('1', inplace=True)
        K['taxon_id'].replace(0,1, inplace=True)

        K.to_csv(self.output, sep='\t', index=False, columns=cols)

    '''
    def _prepare_taxonomy_json(self):

        print('Working on taxonomy tree')
        tree = {}

        if not os.path.exists(self.output):
            print(f"Error: The file '{self.output}' was not found.")
            return tree
 
        # Read the TSV file to build the taxonomic tree

        tax = pd.read_csv(self.output, sep='\t', usecols=['ERZ-contig', 'lineage_names'], dtype=str)
        
        for _, row in tax.iterrows():
            contig_id = row['ERZ-contig']
            lineage_names = row['lineage_names']

            # Split the lineage string by semicolon to get taxonomic levels
            lineage_parts = [part.split('__')[-1] for part in lineage_names.split(';')]
            
            current_level = tree
            
            # Traverse or build the nested dictionary for each taxonomic rank
            for i, part_name in enumerate(lineage_parts):
                if part_name not in current_level:
                    current_level[part_name] = {}
                
                # Check if this is the species level (the last part)
                if i == len(lineage_parts) - 1:
                    # Initialize the final species-level dictionary if it doesn't exist
                    if "ncbi_codes" not in current_level[part_name]:
                        current_level[part_name] = {
                            "ncbi_codes": [],
                            "target_members": []
                        }
                    
                    # Now, gather all associated targets and CDS IDs
                    if contig_id in self.erz_to_mgyp_cds:
                        for target, cds in self.erz_to_mgyp_cds[contig_id]:

                            if not target or str(target).strip() == "" or str(target).lower() == "nan":
                                #target = "not found"
                                continue
                            if not cds or str(cds).strip() == "" or str(cds).lower() == "nan":
                                #cds = "not found"
                                continue

                            # Append the CDS IDs (target members) to the list
                            if cds not in current_level[part_name]["target_members"]:
                                current_level[part_name]["target_members"].append(target)
                                current_level[part_name]["ncbi_codes"].append(target)
                                
                current_level = current_level[part_name]

        tree['root']['Unclassified contig'] = {'target_members': tree['root']['target_members'], 'ncbi_codes': tree['root']['ncbi_codes']}
        
        del tree['root']['target_members']
        del tree['root']['ncbi_codes']
            
        self.tree = tree
        # Save tree as JSON to self.taxonomy_json
        with open(self.taxonomy_json, "w") as json_file:
            json.dump(tree, json_file, indent=2)
    '''

    def _prepare_taxonomy_tree(self):

        tree = {}

        try:
            tax = pd.read_csv(self.output, sep='\t', usecols=['ERZ-contig', 'lineage_names'], dtype=str)
            tax = tax.dropna(subset=['lineage_names'])
            lineage_to_contigs = tax.groupby('lineage_names')['ERZ-contig'].apply(list)
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