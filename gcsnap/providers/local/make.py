import os

from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext

from gcsnap.providers.local.assemblies import Assemblies
from gcsnap.providers.local.sequences import Sequences


class Maker:
    """
    Local-database provider.

    Retrieves flanking genes and sequences from the pre-built GCsnap SQLite
    databases (``assemblies.db`` + ``sequences.db``) and flat GFF files,
    without any network access.

    Expected config keys (set in config.yaml or via CLI):
        ``db-path``   – directory containing assemblies.db and sequences.db
        ``gff-path``  – flat directory with all *.gff.gz annotation files
        ``n-cpu``     – number of parallel workers (default 4)

    The target IDs supplied via ``--local-targets`` must be protein CDS codes
    that exist in the local database (e.g. ``AFI40896.1``).  Each ID is used
    directly as its own CDS code — no external ID mapping is needed.
    """

    def __init__(self, config: Configuration, targets, console):
        self.config  = config
        self.console = console
        self.out_label = config.arguments['out_label']['value']

        self.target_list = targets.get_provider_targets('local')

        # Provider working directory.
        self.local_dir = os.path.join(self.out_label, 'providers', 'local')
        os.makedirs(self.local_dir, exist_ok=True)

        # Assemblies and Sequences still reference n_nodes / n_ranks_per_node
        # for MPI-style chunking.  Those keys were removed from the main config
        # because MPI is optional; we inject sensible serial fallbacks here so
        # the provider classes work without modification.
        n_cpu = config.arguments.get('n_cpu', {}).get('value', 4)
        config.arguments.setdefault('n_nodes',         {'value': 1})
        config.arguments.setdefault('n_ranks_per_node', {'value': n_cpu + 1})

    # ── public API ────────────────────────────────────────────────────────────

    def get_genomic_context(self) -> GenomicContext | None:
        """
        Return a :class:`~gcsnap.genomic_context.GenomicContext` populated
        with flanking genes and sequences for all local targets.

        If a cached JSON file from a previous run exists it is loaded directly,
        skipping the (potentially slow) DB queries.

        Returns ``None`` when the target list is empty.
        """
        if not self.target_list:
            self.console.print_warning('Local provider: no targets supplied.')
            return None

        cache_file = os.path.join(self.local_dir, 'genomic_context_information.json')

        if os.path.exists(cache_file):
            self.console.print_done(
                'Local provider: cached context found, reading from disk.'
            )
            local_gc = GenomicContext(self.config, self.out_label)
            local_gc.read_syntenies_from_json(cache_file)
            return local_gc

        return self._build_context()

    # ── private helpers ───────────────────────────────────────────────────────

    def _build_context(self) -> GenomicContext:
        """
        Query the local databases and build the GenomicContext.

        Step 1 – Assemblies:  map each target to its CDS code (identity
                               mapping for local IDs) and extract flanking
                               genes from the GFF files via the DB.
        Step 2 – Sequences:   fetch amino-acid sequences for all flanking
                               gene CDS codes from sequences.db.
        Step 3 – Persist:     write syntenies to JSON for future caching.
        """
        self.console.print_working_on(
            f'Local provider: building context for {len(self.target_list)} target(s)'
        )

        self.gc = GenomicContext(self.config, self.out_label)
        self.gc.curr_targets = list(self.target_list)

        # Step 1 – flanking genes
        # Local target IDs are CDS codes directly: mapping is (target, target).
        mappings = [(t, t) for t in self.target_list]

        self.assemblies = Assemblies(self.config, mappings)
        self.assemblies.run()
        self.gc.update_syntenies(self.assemblies.get_flanking_genes())

        # Step 2 – sequences
        self.sequences = Sequences(self.config, self.gc)
        self.sequences.run()
        self.gc.update_syntenies(self.sequences.get_sequences())

        # Step 3 – cache to disk
        cache_file = os.path.join(self.local_dir, 'genomic_context_information.json')
        self.gc.write_syntenies_to_json(cache_file)

        self.console.print_done('Local provider: context ready.')
        return self.gc
