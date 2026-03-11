import os
import time
import pandas as pd

# Import your specific modules here
# Assuming these are available in your python path or relative imports
from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext

from gcsnap.providers.MGnify.dataset import Dataset
from gcsnap.providers.MGnify.metadata import Metadata
from gcsnap.providers.MGnify.scraper import Scraper
from gcsnap.providers.MGnify.gatherer import Gatherer
from gcsnap.providers.MGnify.assemblies import Assemblies as MGnifyAssemblies
from gcsnap.providers.MGnify.sequences import Sequences as MGnifySequences

class Maker:

    def __init__(self, config, targets, console):
        """
        Initialize the MGnify provider wrapper.
        
        Args:
            config: The configuration object.
            targets: The targets object containing target lists.
        """
        self.config = config
        self.targets = targets
        self.out_label = config.arguments['out_label']['value']
        self.console = console
        
        self.gatherer = None
        # Filter targets specific to MGnify (containing 'MGYP')
        self.target_list = [
            t for t in targets.get_targets_dict().get(self.out_label, []) 
            if 'MGYP' in t
        ]
        
        # Setup Directories
        self.mgnify_dir = os.path.join(self.out_label, 'providers', 'MGnify')
        os.makedirs(self.mgnify_dir, exist_ok=True)
        
        # Update config with provider directory
        self.config.arguments['MGnify_dir'] = {'value': self.mgnify_dir}

        # initialize dataset
        self.dataset = Dataset(self.config, basename='MGnify')

    def _scrape_metadata(self):

        """Initializes dataset and scrapes metadata."""
        
        self.dataset.set_targets(ids=self.target_list)

        metadata = Metadata(self.dataset)
        metadata.get_mgyp_metadata()
        metadata.assign_contigs()

        self.dataset.set_scraper()
        self.scraper = Scraper(dataset=self.dataset, config=self.config)

        # Scrape target files if not already present
        if not (self.dataset.assembly_present and self.dataset.targets_file_present):
            start = time.time()
            self.scraper.list_target_files()
            elapsed = time.time() - start
            minutes = int(elapsed // 60)
            seconds = elapsed % 60
            print(f"list target files took {minutes} minutes {seconds:.2f} seconds")

        self.scraper.aggregate_assemblies_metadata()
        self.scraper.update_metadata()
        self.dataset.update_after_scrape()
        self.console.print_done(' MGYPs metadata available')
        
        return metadata

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
        self.gc.curr_targets = self.dataset.mgyps

        # 1. Contigs and Assembly Parsing
        mgnify_assemblies = MGnifyAssemblies(self.config, self.dataset)
        mgnify_assemblies.run()
        self.gc.update_syntenies(mgnify_assemblies.get_flanking_genes())        

        self.console.print_working_on('taxonomic assignment of contigs')

        if self.config.arguments['genome_classification']['value']=='taxonomy':
            
            self.dataset.set_taxonomic_assignment()

            genomes = Kraken2Taxonomy(self.dataset)
            genomes.run()
            
            self.dataset.update_after_taxonomic_assignment(genomes)

        #elif self.config.arguments['genome_classification']['value']=='binning':

        #    self.dataset.set_contig_binning()
        #    genomes = SourMashBinning(self.dataset,self.gc,self.config)
        #    genomes.run()

        #    self.dataset.update_after_contig_binning(genomes)

        # 2. Sequence Information
        self.sequences = MGnifySequences(self.config, self.gc, self.dataset)
        self.sequences.run()
        self.gc.update_syntenies(self.sequences.get_sequences())
        
        # 3. Write to disk
        output_file = os.path.join(self.mgnify_dir, 'genomic_context_information.json')
        self.gc.write_syntenies_to_json(output_file)
        
        return self.gc

    def get_genomic_context(self):
        """Main execution method to retrieve MGnify data."""
        if not self.target_list:
            print("No MGnify targets found.")
            return None

        output_file = os.path.join(self.mgnify_dir, 'genomic_context_information.json')

        if os.path.exists(output_file):
            self.console.print_done('Genomic context file already present, reading from disk.')
            mgnify_gc = GenomicContext(self.config, self.out_label)
            mgnify_gc.read_syntenies_from_json(output_file)
            return mgnify_gc
        else:
            # 1. Scrape
            metadata = self._scrape_metadata()

            # 2. Gather (Download)
            self._run_gatherer()

            # 3. Post-processing metadata
            metadata.assign_cds_from_gff(self.dataset)
            self.dataset.update_metadata() 

            # 4. Build Context
            mgnify_gc = self._build_context()

        return mgnify_gc