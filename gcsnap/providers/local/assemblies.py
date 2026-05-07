import os
import gzip # to work with .gz

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.providers.local.db_handler_assemblies import AssembliesDBHandler
from gcsnap.providers.local.parallel_tools import ParallelTools

from gcsnap.utils import split_list_chunks
from gcsnap.utils import WarningToLog

import logging
logger = logging.getLogger('iteration')

class Assemblies:
    """
    Methods and attributes to query assembly accession and contig info from the
    local database and parse flanking genes from per-assembly GFF files.

    Files are organised as:
        gff-path  (config.yaml / CLI)
        └── <assembly_accession>.gff.gz   (one per assembly, all contigs inside)

        fna-path  (optional, config.yaml / CLI)
        └── <assembly_accession>.fna.gz   (matching nucleotide FASTA)

        db-path   (config.yaml / CLI)
        └── GCsnap/
            ├── assemblies.db
            └── sequences.db

    The mappings table inside ``assemblies.db`` already stores the triplet
    ``(seq_code, assembly_accession, contig)``, so the runtime no longer has to
    scan ``##sequence-region`` headers — we just slice CDS lines whose first
    column equals the known contig.

    Attributes:
        n_flanking5 (int): Number of flanking genes to extract at the 5' end of target.
        n_flanking3 (int): Number of flanking genes to extract at the 3' end of target.
        chunks (int): Number of chunks for parallel execution.
        exclude_partial (bool): Exclude partial genomic blocks.
        database_path (str): Directory containing assemblies.db / sequences.db.
        gff_path (str): Flat folder containing all .gff.gz files.
        fna_path (str | None): Flat folder containing all .fna.gz files (optional).
        config (Configuration): Configuration object.
        console (RichConsole): Console object.
        targets_and_cds_codes (list): Input list of (target, cds_code) tuples.
    """

    def __init__(self, config: Configuration, mappings: list[tuple[str, str]]) -> None:
        """
        Initialize the Assemblies object.

        Args:
            config (Configuration): Configuration object containing the arguments.
            mappings (list[tuple[str,str]]): List of (target, cds_code) tuples.
        """
        # required configuration arguments
        self.n_flanking5 = config.arguments['n_flanking5']['value']
        self.n_flanking3 = config.arguments['n_flanking3']['value']
        # -1 as master rank 0 does not compute
        self.chunks = (config.arguments['n_nodes']['value']
                       * config.arguments['n_ranks_per_node']['value']) - 1
        self.exclude_partial = config.arguments['exclude_partial']['value']
        self.database_path = config.arguments['db_path']['value']
        self.gff_path = config.arguments['gff_path']['value']
        # fna is optional — only needed for DNA extraction
        self.fna_path = config.arguments.get('fna_path', {}).get('value')

        self.config = config
        self.console = RichConsole()

        self.gff_suffix = '.gff.gz'
        self.fna_suffix = '.fna.gz'

        # input list with [(target, cds_code)]
        self.targets_and_cds_codes = mappings

    def get_flanking_genes(self) -> dict:
        """Getter for the flanking_genes attribute."""
        return self.flanking_genes

    def run(self) -> None:
        """
        Split inputs into chunks and run parallel extraction.
        """
        parallel_args = split_list_chunks(self.targets_and_cds_codes, self.chunks)

        with self.console.status('Extracting flanking genes from assemblies'):
            dict_list = ParallelTools.parallel_wrapper(parallel_args, self.run_each)

            self.flanking_genes = {k: v for d in dict_list for k, v in d.items()
                                   if v.get('flanking_genes') is not None}
            not_found = {k: v for d in dict_list for k, v in d.items()
                         if v.get('flanking_genes') is None and v.get('msg', '').startswith('File')}
            partial = {k: v for d in dict_list for k, v in d.items()
                       if v.get('flanking_genes') is None and v.get('msg', '').startswith('Partial')}

        if not_found:
            self.log_not_found(not_found)
        if partial:
            self.log_partial(partial)

        if len(self.targets_and_cds_codes) != len(self.flanking_genes):
            missing_targets = [t[0] for t in self.targets_and_cds_codes
                               if t[0] not in self.flanking_genes
                               and t[0] not in not_found
                               and t[0] not in partial]
            missing = {t: {'msg': 'No assembly info found in database'} for t in missing_targets}
            self.log_not_found(missing)

        if not self.flanking_genes:
            self.console.stop_execution(
                msg='No flanking genes found for any target sequence. Continuing is not possible.')

    def run_each(self, args: list[tuple[str, str]]) -> dict[str, dict]:
        """
        Worker:
            - Resolve cds_code -> (assembly_accession, contig) from mappings table.
            - Resolve assembly_accession -> (url, taxid, species) from assemblies table.
            - Open the per-assembly GFF, slice the contig, extract flanking genes.
        """
        target_tuples = args

        # (target, cds, accession, contig)
        accession_tuples = self.get_assembly_accessions(target_tuples)
        # (target, cds, accession, contig, url, taxid, species)
        info_tuples = self.get_assembly_info(accession_tuples)

        flanking_genes = {}
        for element in info_tuples:
            flanking_genes |= self.read_and_parse_assembly(element)

        return flanking_genes

    def read_and_parse_assembly(self, element: tuple) -> dict:
        """
        Read and parse one assembly file for one target.

        Args:
            element (tuple): (target, cds_code, accession, contig, url, taxid, species)
        """
        target, cds_code, accession, contig, url, taxid, species = element

        try:
            gff_file = self.get_gz_path(accession, 'gff')
            lines = self.read_gz_file(gff_file)
            flanking_genes = self.parse_assembly(cds_code, contig, lines)

            assembly_metadata = {
                'source': 'local',
                'target_cds': cds_code,
                'cds_code': cds_code,
                'genomic_region': contig,
                'contig': contig,
                'assembly_accession': accession,
                'assembly_url': url,
                'target': target,
                'target_source': 'local',
                'gff_file': gff_file,
                'species': species,
                'taxID': taxid,
            }
            # only attach dna_file when an fna folder is configured — sequences.py
            # will skip DNA extraction whenever this key is missing
            if self.fna_path is not None:
                assembly_metadata['dna_file'] = self.get_gz_path(url, 'fna')

            return {target: {'flanking_genes': flanking_genes,
                             'assembly_metadata': assembly_metadata}}

        except WarningToLog as e:
            return {target: {'flanking_genes': None,
                             'msg': str(e)}}

    def get_assembly_accessions(self, target_tuples: list[tuple[str, str]]
                                ) -> list[tuple[str, str, str, str]]:
        """
        Query (assembly_accession, contig) for the cds_codes from the mappings table.

        Returns:
            list of (target, cds_code, accession, contig)
        """
        cds_codes = [element[1] for element in target_tuples]
        assembly_db = AssembliesDBHandler(os.path.join(self.database_path))
        # mappings table: (seq_code, assembly_accession, contig)
        result_tuples = assembly_db.select(cds_codes, request_size=5000)
        # build a lookup so we don't do an O(N*M) cross join
        lookup = {row[0]: (row[1], row[2]) for row in result_tuples}

        out = []
        for target, cds in target_tuples:
            if cds in lookup:
                accession, contig = lookup[cds]
                out.append((target, cds, accession, contig))
        return out

    def get_assembly_info(self, target_tuples: list[tuple[str, str, str, str]]
                          ) -> list[tuple[str, str, str, str, str, str, str]]:
        """
        Query (url, taxid, species) for the assembly accessions and append.

        Returns:
            list of (target, cds, accession, contig, url, taxid, species)
        """
        assembly_accessions = [element[2] for element in target_tuples]

        assembly_db = AssembliesDBHandler(os.path.join(self.database_path))
        # assemblies table: (assembly_accession, url, taxid, species)
        result_tuples = assembly_db.select(assembly_accessions,
                                           table='assemblies',
                                           request_size=20000)
        lookup = {row[0]: (row[1], row[2], row[3]) for row in result_tuples}

        out = []
        for target, cds, accession, contig in target_tuples:
            if accession in lookup:
                url, taxid, species = lookup[accession]
                out.append((target, cds, accession, contig, url, taxid, species))
        return out

    def get_gz_path(self, url: str, source: str) -> str:
        """
        Build the on-disk path of an assembly file based on its URL.

        Args:
            url (str): URL or basename as stored in the assemblies table.
            source (str): One of 'gff' or 'fna'.
        """
        if source == 'fna':
            if self.fna_path is None:
                raise WarningToLog('fna_path not configured but fna file requested')
            gz_dir = self.fna_path
            suffix = self.fna_suffix
        elif source == 'gff':
            gz_dir = self.gff_path
            suffix = self.gff_suffix
        else:
            raise ValueError(f'Unknown source {source}')

        file = os.path.basename(url) + suffix
        return os.path.join(gz_dir, file)

    def read_gz_file(self, file_path: str) -> list:
        """Read a gzipped text file and return its lines."""
        try:
            with gzip.open(file_path, 'rt', encoding='utf-8') as file:
                content = file.read()
            return content.splitlines()
        except FileNotFoundError:
            raise WarningToLog('File {} not found'.format(file_path))

    def parse_assembly(self, cds_code: str, contig: str, lines: list) -> dict:
        """Wrapper around context block extraction + parsing."""
        genomic_context_block = self.extract_genomic_context_block(cds_code, contig, lines)
        return self.parse_genomic_context_block(cds_code, genomic_context_block)

    def extract_genomic_context_block(self, target_cds_code: str,
                                      contig: str, lines: list) -> list:
        """
        Slice the genomic context for the target out of the GFF.

        Because the contig name is now known up-front (from the mappings table),
        we can build the per-contig CDS list in a single pass over the file —
        no need to scan ``##sequence-region`` markers, no scaffold offset
        bookkeeping. This is the main runtime win from carrying contig info in
        the database.
        """
        # one pass: keep CDS lines whose chromosome column matches our contig
        scaffold = []
        for line in lines:
            if not line or line[0] == '#':
                continue
            cols = line.split('\t')
            if len(cols) < 9 or cols[2] != 'CDS':
                continue
            if cols[0] != contig:
                continue
            scaffold.append(line)

        if not scaffold:
            raise WarningToLog('Contig {} not found in GFF for {}'.format(contig, target_cds_code))

        # locate target inside the contig
        target_indices = [i for i, val in enumerate(scaffold)
                          if 'ID=cds-{}'.format(target_cds_code) in val
                          or 'Name={}'.format(target_cds_code) in val
                          or 'protein_id={}'.format(target_cds_code) in val
                          or 'locus_tag={}'.format(target_cds_code) in val]

        if not target_indices:
            raise WarningToLog('{} not found in contig {}'.format(target_cds_code, contig))

        index_of_target = target_indices[0]
        direction_of_target = scaffold[index_of_target].split('\t')[6]

        if direction_of_target == '+':
            start = max(0, index_of_target - self.n_flanking5)
            end = index_of_target + self.n_flanking3 + 1
            genomic_context_block = scaffold[start:end]
        else:
            start = max(0, index_of_target - self.n_flanking3)
            end = index_of_target + self.n_flanking5 + 1
            genomic_context_block = scaffold[start:end]
            # reverse if target points in the - direction
            genomic_context_block = genomic_context_block[::-1]

        if self.exclude_partial and len(genomic_context_block) < (self.n_flanking5 + self.n_flanking3 + 1):
            raise WarningToLog('Partial genomic block for {} excluded!'.format(target_cds_code))

        return genomic_context_block

    def parse_genomic_context_block(self, target_cds_code: str,
                                    genomic_context_block: list) -> dict:
        """Parse the GFF block into the GenomicContext flanking_genes dict."""
        flanking_genes = GenomicContext.get_empty_flanking_genes()

        for line in genomic_context_block:
            line_data = line.split('\t')
            start = int(line_data[3])
            end = int(line_data[4])
            direction = line_data[6]
            attrs = line_data[8]

            if 'cds-' in attrs:
                cds_code = attrs.split('ID=cds-')[1].split(';')[0]
            elif 'Name=' in attrs:
                cds_code = attrs.split('Name=')[1].split(';')[0]
            elif 'locus_tag=' in attrs:
                cds_code = attrs.split('locus_tag=')[1].split(';')[0]
            else:
                cds_code = 'unk'

            if 'pseudo=' not in attrs and 'product=' in attrs and 'fragment' not in attrs:
                prot_name = attrs.split('product=')[1].split(';')[0]
            else:
                prot_name = 'pseudogene'

            # merge fragmented genes (introns / split CDS) into one entry
            if (flanking_genes['cds_codes']
                    and flanking_genes['cds_codes'][-1] == cds_code):
                if start < flanking_genes['starts'][-1]:
                    flanking_genes['starts'][-1] = start
                if end > flanking_genes['ends'][-1]:
                    flanking_genes['ends'][-1] = end
                continue

            if '|' in cds_code:
                cds_code = cds_code.replace('|', '_')

            flanking_genes['cds_codes'].append(cds_code)
            flanking_genes['names'].append(prot_name)
            flanking_genes['starts'].append(start)
            flanking_genes['ends'].append(end)
            flanking_genes['directions'].append(direction)

        # locate target in the parsed list (it may differ from the input cds_code
        # if locus_tag-based parsing collapsed it)
        try:
            index_of_target = flanking_genes['cds_codes'].index(target_cds_code)
        except ValueError:
            raise WarningToLog('Target {} missing from parsed block'.format(target_cds_code))

        direction_of_target = flanking_genes['directions'][index_of_target]

        if direction_of_target == '+':
            for key in ['starts', 'ends']:
                lst = [e - flanking_genes['starts'][index_of_target] + 1
                       for e in flanking_genes[key]]
                flanking_genes['relative_{}'.format(key)] = lst
        else:
            base = flanking_genes['ends'][index_of_target]
            flanking_genes['relative_starts'] = [base - e + 1 for e in flanking_genes['ends']]
            flanking_genes['relative_ends'] = [base - e + 1 for e in flanking_genes['starts']]
            flanking_genes['directions'] = ['+' if d == '-' else '-'
                                            for d in flanking_genes['directions']]

        return flanking_genes

    def log_not_found(self, not_found: dict) -> None:
        message = 'No flanking genes found for {} target sequences.'.format(len(not_found))
        self.console.print_warning(message)
        for k, v in not_found.items():
            logger.warning('For target {}: {}'.format(k, v.get('msg')))

    def log_partial(self, partial: dict) -> None:
        message = ('Partial genomic blocks for {} target sequences excluded. '
                   'Set --exclude-partial False to include.').format(len(partial))
        self.console.print_warning(message)
        for k, v in partial.items():
            logger.warning('For target {}: {}'.format(k, v.get('msg')))
