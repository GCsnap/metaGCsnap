import os
import time
import pandas as pd

# Import your specific modules here
# Assuming these are available in your python path or relative imports
from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext

from providers.ncbi.dataset import Dataset
from providers.ncbi.scraper import Scraper
from providers.ncbi.sequence_mapping import SequenceMapping
from providers.ncbi.gatherer import Gatherer
from providers.ncbi.assemblies import Assemblies
from providers.ncbi.sequences import Sequences


class Maker:
    
    def __init__(self, config, targets, console):
        """
        Initialize the NCBI provider wrapper.
        
        Args:
            config: The configuration object.
            targets: The targets object containing target lists.
        """
        self.config = config
        self.targets = targets
        self.out_label = config.arguments['out_label']['value']
        self.console = console
        self.dataset = None
        self.gatherer = None
        # Filter targets specific to NCBI (containing 'NCBI')
        self.target_list = [
            t for t in targets.get_targets_dict().get(self.out_label, []) 
            if not 'MGYP' in t
        ]
        
        # Setup ncbi
        self.ncbi_dir = os.path.join(self.out_label, 'providers', 'ncbi')
        os.makedirs(self.ncbi_dir, exist_ok=True)
        
        # Update config with provider directory
        self.config.arguments['ncbi_dir'] = {'value': self.ncbi_dir}

    def _resolve_mappings(self):
        """
        Handles ID mapping (UniProt -> RefSeq -> EMBL).
        Checks for existing mapping file to avoid redundant processing.
        """
        self.mapping_file = os.path.join(self.dataset.metadata_dir, 'mapping.csv')
        
        if not os.path.exists(self.mapping_file):

            print('Creating new mapping file:', self.mapping_file)
            
            # 1. Map to RefSeq
            mapping_a = SequenceMapping(dataset=self.dataset, config=self.config, target_list=self.target_list, to_type='UniProtKB-AC')
            mapping_a.run()
        
            mapping_b = SequenceMapping(dataset=self.dataset, config=self.config, target_list=mapping_a.get_codes(), to_type='RefSeq')
            mapping_b.run()
        
            # Merge RefSeq results
            mapping_a.merge_mapping_dfs(mapping_b.mapping_df, columns_to_merge=['RefSeq'])
            
            # 2. Map to EMBL-CDS
            mapping_c = SequenceMapping(dataset=self.dataset, config=self.config, target_list=mapping_a.get_codes(), to_type='EMBL-CDS')
            mapping_c.run()
        
            # 3. Finalize
            mapping_a.merge_mapping_dfs(mapping_c.mapping_df)
            mapping_a.finalize()
            
            return mapping_a.get_targets_and_ncbi_codes()
        
        else:
            print('Reading existing mapping file:', self.mapping_file)
            self.mappings = pd.read_csv(self.mapping_file)
            
            # Return list of tuples (target, ncbi_code)
            return self.mappings.dropna(subset=['target', 'ncbi_code']) \
                           [['target', 'ncbi_code']] \
                           .apply(tuple, axis=1).tolist()

    def _scrape_metadata(self):

        """Initializes dataset and scrapes metadata."""
        self.dataset = Dataset(self.config, basename='ncbi')
        self.dataset.set_targets(ids=self.target_list)

        self.targets_and_ncbi_codes = self._resolve_mappings()

        self.dataset.set_scraper()
        scraper = Scraper(dataset=self.dataset, config=self.config, mappings= self.targets_and_ncbi_codes)

        # Scrape target files if not already present
        if not (self.dataset.assembly_present and self.dataset.targets_file_present):
            start = time.time()
            scraper.run()
            elapsed = time.time() - start
            minutes = int(elapsed // 60)
            seconds = elapsed % 60
            print(f"list target files took {minutes} minutes {seconds:.2f} seconds")

        self.dataset.update_after_scrape()
        self.console.print_done(' NCBI metadata available')
        
    def _run_gatherer(self):

        """Runs the download pipeline."""
        start = time.time()
        
        self.dataset.set_gatherer()
        self.gatherer = Gatherer(dataset=self.dataset, config=self.config)
        self.gatherer.run_pipeline()
        self.dataset.update_after_gathering(self.gatherer)
        
        self.console.print_done(f"download_files took {time.time() - start:.2f} seconds")

    def _build_context(self):

        """Parses assemblies and sequences to build Genomic Context."""
        self.console.print_working_on('genomic context')
        
        self.gc = GenomicContext(self.config, self.out_label)
        self.gc.curr_targets = self.dataset.ncbip

        # 1. Contigs and Assembly Parsing
        self.assemblies = Assemblies(self.config, self.dataset)
        self.assemblies.run()
        self.gc.update_syntenies(self.assemblies.get_flanking_genes())        

        self.console.print_working_on('taxonomic assignment of contigs')

        # 2. Sequence Information
        self.sequences = Sequences(self.config, self.gc, self.dataset)
        self.sequences.run()
        self.gc.update_syntenies(self.sequences.get_sequences())
        
        # 3. Write to disk
        output_file = os.path.join(self.ncbi_dir, 'genomic_context_information.json')
        self.gc.write_syntenies_to_json(output_file)
        
        return self.gc

    def get_genomic_context(self):
        
        """Main execution method to retrieve ncbi data."""
        
        if not self.target_list:
            print("No ncbi targets found.")
            return None

        output_file = os.path.join(self.ncbi_dir, 'genomic_context_information.json')

        if os.path.exists(output_file):
            self.console.print_done('Genomic context file already present, reading from disk.')
            ncbi_gc = GenomicContext(self.config, self.out_label)
            ncbi_gc.read_syntenies_from_json(output_file)
            return ncbi_gc
        else:
            # 1. Scrape
            
            metadata = self._scrape_metadata()

            # 2. Gather (Download)
            self._run_gatherer()

            # 3. Post-processing metadata
            self.dataset.update_metadata() 

            # 4. Build Context
            ncbi_gc = self._build_context()

        return ncbi_gc