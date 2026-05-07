import os
import gzip

from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.utils import split_dict_chunks

from gcsnap.providers.local.db_handler_sequences import SequenceDBHandler
from gcsnap.providers.local.parallel_tools import ParallelTools


def extract_dna_sequences(dna_file: str, contig: str,
                          starts: list, ends: list) -> dict:
    """
    Read a per-assembly nucleotide FASTA (gzipped) and extract the genomic
    context plus the individual flanking-gene DNA sequences.

    Mirrors the MGnify provider's helper but works with per-assembly .fna.gz
    files: it scans records until it hits the requested contig, then slices.

    Args:
        dna_file (str): Path to the gzipped FASTA file.
        contig (str): Contig record id to locate inside the FASTA.
        starts (list): 1-based start positions of the flanking genes.
        ends (list): 1-based inclusive end positions of the flanking genes.

    Returns:
        dict: ``{'context': str, 'features': [str, ...]}``. Empty strings on
              failure so the caller doesn't need to handle missing files.
    """
    # local import keeps Bio out of the import path when DNA isn't requested
    try:
        from Bio import SeqIO
    except ImportError:
        return {'context': '', 'features': ['' for _ in starts]}

    try:
        with gzip.open(dna_file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if record.id != contig:
                    continue
                seq = record.seq
                
                gc_start = min(starts) - 1
                gc_end = max(ends)
                context = str(seq[gc_start:gc_end])
                features = [str(seq[s - 1:e]) for s, e in zip(starts, ends)]
                return {'context': context, 'features': features}
    except FileNotFoundError:
        pass
    return {'context': '', 'features': ['' for _ in starts]}


class Sequences:
    """
    Attach sequences and taxonomy to the flanking genes.

    Protein sequences are pulled from the local ``sequences.db``. DNA
    sequences (per-gene + the full genomic context) are extracted on the fly
    from the per-assembly .fna.gz files when an ``fna_path`` is configured —
    matching the MGnify provider's behaviour as closely as possible.

    Attributes:
        database_path (str): Path to the directory containing sequences.db.
        gc (GenomicContext): GenomicContext to enrich.
        chunks (int): Parallel chunk count.
        sequences_dict (dict): {cds_code: protein_seq} from sequences.db.
        sequences (dict): Final enriched syntenies.
    """

    FAKE_PROTEIN = 'FAKESEQUENCEFAKESEQUENCEFAKESEQUENCEFAKESEQUENCE'

    def __init__(self, config: Configuration, gc: GenomicContext):
        self.config = config
        self.database_path = config.arguments['db_path']['value']
        self.chunks = (config.arguments['n_nodes']['value']
                       * config.arguments['n_ranks_per_node']['value']) - 1
        self.gc = gc
        self.console = RichConsole()

    def get_sequences(self) -> dict:
        return self.sequences

    def run(self) -> None:
        """
        Pre-fetch protein sequences from sequences.db, then dispatch DNA /
        taxonomy enrichment in parallel chunks.
        """
        # 1. protein sequences for *all* flanking genes — one DB round-trip
        all_cds_codes = self.gc.get_all_cds_codes()
        self.find_protein_sequences(all_cds_codes)

        # 2. parallel per-target enrichment (DNA + taxonomy lift)
        syntenies = self.gc.get_syntenies()
        parallel_args = split_dict_chunks(syntenies, self.chunks)

        with self.console.status('Add sequences, tax id and species name to flanking genes'):
            dict_list = ParallelTools.parallel_wrapper(parallel_args, self.run_each)
            self.sequences = {k: v for d in dict_list for k, v in d.items()}

    def run_each(self, args: dict) -> dict:
        """
        Per-chunk worker. Pulls protein sequences from the dict cached in
        ``self.sequences_dict``, opens the FNA once per target to extract DNA,
        and lifts taxonomy out of assembly_metadata.
        """
        targets_content = args
        results = {}

        for target, content_dict in targets_content.items():
            cds_codes = content_dict['flanking_genes']['cds_codes']

            # protein sequences (already in cache)
            content_dict['flanking_genes']['protein_sequences'] = [
                self.get_protein_sequence(c) for c in cds_codes
            ]

            # DNA sequences — only if dna_file was attached upstream
            assembly_metadata = content_dict['assembly_metadata']
            dna_file = assembly_metadata.get('dna_file')
            contig = assembly_metadata.get('genomic_region') or assembly_metadata.get('contig')

            if dna_file and contig:
                starts = content_dict['flanking_genes']['starts']
                ends = content_dict['flanking_genes']['ends']
                seqs = extract_dna_sequences(dna_file, contig, starts, ends)
                content_dict['flanking_genes']['dna_sequences'] = seqs['features']
                content_dict['context_sequence'] = seqs['context']
            else:
                # keep the schema stable even when DNA isn't available
                content_dict['flanking_genes']['dna_sequences'] = ['' for _ in cds_codes]
                content_dict['context_sequence'] = ''

            # lift taxonomy out of assembly_metadata into its own top-level key
            content_dict['taxonomy'] = {
                'taxon_id': assembly_metadata.pop('taxID', None),
                'taxon_name': assembly_metadata.pop('species', None),
            }

            results[target] = content_dict

        return results

    def find_protein_sequences(self, cds_codes: list) -> None:
        """
        Bulk-load all needed protein sequences from sequences.db into a dict.
        """
        sequences_db = SequenceDBHandler(os.path.join(self.database_path))
        self.sequences_dict = sequences_db.select_as_dict(cds_codes)

    def get_protein_sequence(self, cds_code: str) -> str:
        """Look up a protein sequence by cds_code, returning a fake on miss."""
        return self.sequences_dict.get(cds_code, self.FAKE_PROTEIN)
