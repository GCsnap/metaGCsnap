"""
GCsnap – main entry point.

Provider flags (at least one required):

  --ncbi-targets   PATH   plain-text file with NCBI / UniProt target IDs
  --mgnify-targets PATH   plain-text file with MGnify (MGYP…) target IDs
  --local-targets  PATH   plain-text file with local-database target IDs

If a provider flag is given but its optional dependencies are not installed,
GCsnap exits immediately with an actionable install message instead of
crashing mid-run.
"""

import sys
import json
from pathlib import Path

from gcsnap.rich_console import RichConsole
from gcsnap.configuration import Configuration
from gcsnap.timing import Timing
from gcsnap.genomic_context import GenomicContext
from gcsnap.utils import CustomLogger

from gcsnap.targets import Target
from gcsnap.annotations.families import Families
from gcsnap.annotations.operons import Operons
from gcsnap.visual.runner import Figures


# ── constants ─────────────────────────────────────────────────────────────────

# Maps CLI provider name → sub-directory used inside <out_label>/providers/
_PROVIDER_DIRS = {
    'ncbi':   'ncbi',
    'mgnify': 'MGnify',
    'local':  'local',
}

# ── helpers ───────────────────────────────────────────────────────────────────

def _read_ids(path: str) -> list | None:
    """
    Read a plain-text file of IDs – one per line, # comments ignored.

    Returns None (instead of raising) when the file does not exist, so the
    caller can print a friendly message and exit cleanly.
    """
    if not Path(path).is_file():
        return None
    with open(path) as fh:
        return [
            line.strip()
            for line in fh
            if line.strip() and not line.startswith('#')
        ]


def _validate_provider_paths(active: dict, config_args: dict) -> None:
    """
    Check that every path required by an active provider is set and exists on
    disk. Exits immediately with an actionable message if a path is missing.
    """
    required = {
        'local':  ['db_path', 'gff_path'],
        'mgnify': ['MGnify_path'],
    }
    config_key_to_flag = {
        'db_path':    'db-path',
        'gff_path':   'gff-path',
        'MGnify_path': 'MGnify-path',
    }
    for provider in active:
        for key in required.get(provider, []):
            value = config_args.get(key, {}).get('value')
            flag  = config_key_to_flag[key]
            if not value:
                print(
                    f'\n[ERROR] Provider "{provider}" requires --{flag} '
                    f'but it is not set in config.yaml.\n'
                    f'        Set {flag} in config.yaml or pass --{flag} on the CLI.\n'
                )
                sys.exit(1)
            if not Path(value).exists():
                print(
                    f'\n[ERROR] Provider "{provider}" requires --{flag} '
                    f'but the path does not exist: {value}\n'
                    f'        Check that the path is correct in config.yaml.\n'
                )
                sys.exit(1)


def _load_provider(name: str):
    """
    Lazy-import a provider Maker class.

    Exits with a clear install hint when optional provider dependencies are
    missing, rather than letting a bare ImportError surface mid-run.
    """
    try:
        if name == 'ncbi':
            from gcsnap.providers.ncbi.make import Maker
        elif name == 'mgnify':
            from gcsnap.providers.MGnify.make import Maker
        elif name == 'local':
            from gcsnap.providers.local.make import Maker
        else:
            raise ValueError(f'Unknown provider: {name}')
        return Maker
    except ImportError as exc:
        print(
            f'\n[ERROR] Provider "{name}" requires additional dependencies '
            f'that are not installed.\n'
            f'        Run:  bash install_providers.sh --{name}\n'
            f'        Details: {exc}\n'
        )
        sys.exit(1)



# ── pipeline ──────────────────────────────────────────────────────────────────

def main():
    """
    GCsnap pipeline.

    A. Parse configuration and arguments.
    B. Dispatch targets to active providers; merge syntenies into a single
       GenomicContext.  If a previous complete run exists, load from disk
       and skip straight to figures.
    C. Find protein families (MMseqs2 clustering).
    D. Find operon / genomic-context types.
    E. Produce figures.
    F. Write outputs.
    """

    # ── initialise logging & console ─────────────────────────────────────
    CustomLogger.configure_loggers()
    console = RichConsole('base')
    console.print_title()

    # ── A. configuration ─────────────────────────────────────────────────
    config = Configuration()
    config.parse_arguments()

    out_label = config.arguments['out_label']['value']

    working_dir = Path(out_label).resolve()
    working_dir.mkdir(parents=True, exist_ok=True)

    CustomLogger.configure_iteration_logger(out_label, str(working_dir))
    RichConsole.set_out_label(out_label)
    config.write_configuration_yaml_log('input_arguments.log', str(working_dir))

    timing = Timing(config.arguments.get('timing', {}).get('value', False))
    t_all  = timing.timer('All steps')

    # ── provider flags ────────────────────────────────────────────────────
    provider_flags = {
        'ncbi':   config.arguments.get('ncbi_targets',   {}).get('value'),
        'mgnify': config.arguments.get('mgnify_targets', {}).get('value'),
        'local':  config.arguments.get('local_targets',  {}).get('value'),
    }
    active = {name: path for name, path in provider_flags.items() if path}

    if not active:
        console.print_warning(
            'No provider targets supplied.\n'
            'Use --ncbi-targets, --mgnify-targets, or --local-targets '
            '(each pointing to a plain-text file of IDs, one per line).'
        )
        sys.exit(1)

    _validate_provider_paths(active, config.arguments)

    # Build the Target object and populate per-provider ID lists.
    targets = Target(config)
    #all_target_ids: list = []
    for provider, path in active.items():
        ids = _read_ids(path)
        if ids is None:
            console.print_warning(
                f'Target file for provider "{provider}" was not found: {path}\n'
                'Please check the path and try again. Quitting.'
            )
            sys.exit(0)
        targets.add_provider_targets(out_label, provider, ids)
        #all_target_ids.extend(ids)

    # ── B. collect genomic context ────────────────────────────────────────
    gc = GenomicContext(config, str(working_dir))

    if gc.all_init_files_present:
        # Fast path: previous complete run – load everything from disk.
        console.print_done('All result files present – loading from disk.')
        gc.load_from_files(targets.all_ids)

    else:
        t_collect = timing.timer('Step 1: Collecting genomic contexts')
        gc.curr_targets = []
        for provider_name in active:
            provider_ids = targets.get_provider_targets(provider_name)
            console.print_working_on(
                f'Provider [{provider_name}] – {len(provider_ids)} target(s)'
            )

            # Fast path: if the provider already wrote its JSON on a previous
            # (possibly partial) run, load it directly without importing the
            # provider or re-running any network/DB queries.
            provider_cache = (
                working_dir / 'providers'
                / _PROVIDER_DIRS[provider_name]
                / 'genomic_context_information.json'
            )
            if provider_cache.is_file():
                console.print_done(
                    f'Provider [{provider_name}] – cache found, loading from disk.'
                )
                provider_gc = GenomicContext(config, str(working_dir))
                try:
                    provider_gc.read_syntenies_from_json(str(provider_cache))
                    provider_gc.curr_targets = provider_ids
                except (json.JSONDecodeError, OSError):
                    console.print_warning(
                        f'Provider [{provider_name}] – cache file is corrupt, re-running provider.'
                    )
                    Maker       = _load_provider(provider_name)
                    maker       = Maker(config, targets, console)
                    provider_gc = maker.get_genomic_context()
            else:
                Maker       = _load_provider(provider_name)
                maker       = Maker(config, targets, console)
                provider_gc = maker.get_genomic_context()

            if provider_gc is not None:
                gc.update_syntenies(provider_gc.syntenies)
                gc.curr_targets.extend(provider_gc.curr_targets)
            else:
                console.print_warning(f'Provider [{provider_name}]: no genomic context returned.')

        gc.write_syntenies_to_json()
        t_collect.stop()

    # ── guard ─────────────────────────────────────────────────────────────
    if len(gc.curr_targets) < 2:
        console.print_warning(
            f'Only {len(gc.curr_targets)} target(s) collected – '
            'need at least 2 to analyse genomic context.'
        )
        sys.exit(0)

    console.print_working_on(
        f'Task {out_label} – {len(gc.curr_targets)} total target(s)'
    )

    if config.arguments['collect_only']['value']:
        console.print_skipped_step(
            'collect_only mode: stopping after context collection.'
        )

    else:

        # ── C. protein families ───────────────────────────────────────────
        t_families = timing.timer('Step 2: Finding protein families')
        families = Families(config, gc, str(working_dir))
        families.run()
        gc.update_syntenies(families.get_families())
        gc.create_and_write_families_summary(families)
        t_families.stop()

        # ── D. operon types ───────────────────────────────────────────────
        t_operons = timing.timer('Step 3: Finding operon types')
        operons = Operons(config, gc, str(working_dir))
        operons.run()
        gc.update_syntenies(operons.get_operons())
        gc.create_and_write_operon_types_summary()
        gc.find_most_populated_operon_types()
        t_operons.stop()

        # ── E. figures ────────────────────────────────────────────────────
        t_figures = timing.timer('Step 4: Figures')
        figures = Figures(config, gc, str(working_dir))
        figures.run()
        t_figures.stop()

        # ── F. outputs ────────────────────────────────────────────────────
        t_output = timing.timer('Step 5: Write outputs')
        gc.write_summary_table(
            f'{working_dir.name}_summary_table.tab',
            str(working_dir),
        )
        gc.write_families_to_json()
        t_output.stop()

    # ── wrap up ───────────────────────────────────────────────────────────
    gc.write_syntenies_to_json()
    CustomLogger.log_to_iteration(
        f'Successfully finished task {out_label} '
        f'with {len(gc.curr_targets)} targets.'
    )
    console.print_done(
        f'Task {out_label} with {len(gc.curr_targets)} targets'
    )
    t_all.stop()
    timing.to_csv(str(working_dir / 'timing.csv'))
    console.print_final()


if __name__ == '__main__':
    main()
