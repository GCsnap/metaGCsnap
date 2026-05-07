import os
import gzip
import glob
import time
import datetime
import sqlite3
import argparse

# import database handlers
from gcsnap.providers.local.db_handler_assemblies import AssembliesDBHandler
from gcsnap.providers.local.db_handler_sequences import SequenceDBHandler

from gcsnap.utils import processpool_wrapper


def parse_gff_for_contigs(gff_path: str) -> dict[str, str]:
    """
    Parse a single GFF file and return {protein_id: contig}.

    Contig = column 0 of each CDS line.
    Protein ID is taken from ``locus_tag=`` first, then ``ID=cds-``
    as fallback (matching what the FAA headers produce).
    """
    result: dict[str, str] = {}
    try:
        opener = gzip.open if gff_path.endswith('.gz') else open
        with opener(gff_path, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                cols = line.split('\t')
                if len(cols) < 9 or cols[2] != 'CDS':
                    continue
                contig = cols[0]
                attrs = cols[8]

                protein_id = None
                #if 'locus_tag=' in attrs:
                #    protein_id = attrs.split('locus_tag=')[1].split(';')[0].strip()
                if 'ID=cds-' in attrs:
                    protein_id = attrs.split('ID=cds-')[1].split(';')[0].strip()

                if protein_id:
                    result[protein_id] = contig
    except FileNotFoundError:
        pass
    return result


def parse_gff_batch(gff_paths: list[str]) -> dict[str, str]:
    """Parse a list of GFF files and merge the results."""
    merged: dict[str, str] = {}
    for p in gff_paths:
        merged.update(parse_gff_for_contigs(p))
    return merged


def build_protein_to_contig(gff_dir: str, assembly_accessions: set[str],
                            n_processes: int) -> dict[str, str]:
    """
    For every assembly in *assembly_accessions*, locate its GFF file in
    *gff_dir* and build a global ``{protein_id: contig}`` lookup.

    GFF filename convention: ``{assembly_accession}.gff.gz``.
    """
    gff_paths = []
    for acc in assembly_accessions:
        candidate = os.path.join(gff_dir, f'{acc}.gff.gz')
        if os.path.exists(candidate):
            gff_paths.append(candidate)

    if not gff_paths:
        return {}

    # split into n_processes roughly-equal chunks
    batches = split_into_parts(gff_paths, max(1, n_processes))
    parallel_args = [b for b in batches if b]
    results = processpool_wrapper(n_processes, parallel_args, parse_gff_batch)

    merged: dict[str, str] = {}
    for d in results:
        merged.update(d)
    return merged


def split_into_batches(lst: list, batch_size: int = 1000):
    for i in range(0, len(lst), batch_size):
        yield lst[i:i + batch_size]


def split_into_parts(lst: list, n: int) -> list[list]:
    q, r = divmod(len(lst), n)
    return [lst[i * q + min(i, r):(i + 1) * q + min(i + 1, r)] for i in range(n)]


def execute_handler(args: tuple):
    handler, batch = args
    sequences, mapping = handler.parse_sequences_from_faa_files(batch)
    return sequences, mapping


def execute_handlers(handler: SequenceDBHandler, batch: list, n_processes: int) -> tuple[list[tuple[str,str]], list[tuple[str,str]]]:
    batches = split_into_parts(batch, n_processes)
    parallel_args = [(handler, subbatch) for subbatch in batches]
    result = processpool_wrapper(n_processes, parallel_args, execute_handler)
    mappings = [item for sublist in result for item in sublist[1]]
    sequences = [item for sublist in result for item in sublist[0]]
    return (sequences, mappings)


def reindex(handler) -> None:
    handler.reindex()


def fetch_sample_rows(db_file: str, table: str, n: int = 5) -> list[tuple]:
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute(f'SELECT * FROM {table} LIMIT {n}')
    rows = cursor.fetchall()
    conn.close()
    return rows


def write_readme(out_dir: str, faa_dir: str, summary_files: list[str],
                 n_assemblies: int, n_sequences: int, n_unique: int,
                 elapsed: str) -> None:
    asse_db = os.path.join(out_dir, 'assemblies.db')
    seq_db  = os.path.join(out_dir, 'sequences.db')

    sample_assemblies = fetch_sample_rows(asse_db, 'assemblies')
    sample_mappings   = fetch_sample_rows(asse_db, 'mappings')
    sample_sequences  = fetch_sample_rows(seq_db,  'sequences')

    def md_table(headers: list[str], rows: list[tuple]) -> str:
        sep = '| ' + ' | '.join(['---'] * len(headers)) + ' |'
        header_row = '| ' + ' | '.join(headers) + ' |'
        data_rows = ['| ' + ' | '.join(str(c)[:40] for c in row) + ' |' for row in rows]
        return '\n'.join([header_row, sep] + data_rows)

    lines = [
        '# GCsnap Database',
        '',
        f'**Created:** {datetime.date.today()}  ',
        f'**Build time:** {elapsed}  ',
        '',
        '## Input resources',
        '',
        f'**FAA directory:** `{faa_dir}`  ',
        '',
        '**Assembly summary files:**',
    ]
    for sf in summary_files:
        lines.append(f'- `{sf}`')

    lines += [
        '',
        '## Statistics',
        '',
        f'- Assemblies processed: {n_assemblies:,}',
        f'- Total protein sequences: {n_sequences:,}',
        f'- Unique protein sequences: {n_unique:,}',
        '',
        '## Database files',
        '',
        '- `assemblies.db` — tables: `assemblies` (accession, url, taxid, species), `mappings` (seq_code → assembly_accession, contig)',
        '- `sequences.db` — table: `sequences` (seq_code, sequence)',
        '',
        '## Sample entries',
        '',
        '### assemblies (first 5 rows)',
        '',
        md_table(['assembly_accession', 'url', 'taxid', 'species'], sample_assemblies),
        '',
        '### mappings (first 5 rows)',
        '',
        md_table(['seq_code', 'assembly_accession', 'contig'], sample_mappings),
        '',
        '### sequences (first 5 rows)',
        '',
        md_table(['seq_code', 'sequence (truncated to 40 chars)'], sample_sequences),
    ]

    readme_path = os.path.join(out_dir, 'README.md')
    with open(readme_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'README written to {readme_path}')


def create_dbs(faa_dir: str, summary_files: list[str], out_dir: str,
               n_processes: int, gff_dir: str | None = None) -> tuple[int, int, int]:

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # open handlers and create tables
    asse_db_handler = AssembliesDBHandler(out_dir, 'assemblies.db')
    asse_db_handler.create_tables()
    seq_db_handler  = SequenceDBHandler(out_dir, 'sequences.db')
    seq_db_handler.create_table()

    print(f'Assembly database created: {os.path.join(out_dir, "assemblies.db")}')
    print(f'Sequence database created: {os.path.join(out_dir, "sequences.db")}')

    batch_size = n_processes * 500
    n_assemblies = 0
    n_sequences  = 0
    n_missing_contig = 0

    # 1. fill assemblies table from each summary file
    for summary_file in summary_files:
        asse_db_handler.insert_assemblies_from_summary(summary_file)
        print(f'Summary {os.path.basename(summary_file)} done')

    # 2. parse all .faa.gz files from the flat faa directory
    file_paths = glob.glob(os.path.join(faa_dir, '*.gz'))

    for batch in split_into_batches(file_paths, batch_size):
        sequence_list, mapping_list = execute_handlers(seq_db_handler, batch, n_processes)
        n_sequences  += len(mapping_list)
        n_assemblies += len(batch)

        # Enrich mappings with contig info from the matching GFFs.
        # Each entry becomes (seq_code, assembly_accession, contig).
        if gff_dir is not None:
            batch_accessions = {acc for _, acc in mapping_list}
            protein_to_contig = build_protein_to_contig(
                gff_dir, batch_accessions, n_processes,
            )

            enriched = []
            for seq_code, acc in mapping_list:
                contig = protein_to_contig.get(seq_code)
                if contig is None:
                    n_missing_contig += 1
                enriched.append((seq_code, acc, contig))
            mapping_list = enriched
        else:
            # No GFF dir provided — keep contig as NULL so the schema still matches.
            mapping_list = [(seq_code, acc, None) for seq_code, acc in mapping_list]

        seq_db_handler.insert_sequences(sequence_list)
        asse_db_handler.insert_mappings(mapping_list)

        print(f'{n_assemblies:,} assemblies and {n_sequences:,} sequences done so far')

    if gff_dir is not None and n_missing_contig:
        print(f'Warning: {n_missing_contig:,} protein IDs had no matching contig in the GFFs.')

    # reindex both databases in parallel
    parallel_args = [asse_db_handler, seq_db_handler]
    processpool_wrapper(2, parallel_args, reindex)

    print('Assembly indexing done')
    print('Sequence indexing done')

    n_unique = seq_db_handler.select_number_of_entries()

    return n_assemblies, n_sequences, n_unique


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create assemblies and sequences databases.')
    parser.add_argument('-p', '--processes', type=int, required=True,
                        help='Number of parallel processes to use.')
    parser.add_argument('--faa-dir', type=str, required=True,
                        help='Flat directory containing all *.gz files.')
    parser.add_argument('--gff-dir', type=str, required=False, default=None,
                        help='Flat directory containing all *.gff.gz files. '
                             'When provided, the mappings table is enriched with '
                             'a contig column so the runtime provider can skip '
                             'scaffold-scanning in the GFF.')
    parser.add_argument('--summary-files', nargs='+', type=str, required=True,
                        help='One or more assembly_summary_*.txt files to load into assemblies.db.')
    parser.add_argument('--out-dir', type=str, required=True,
                        help='Output directory for the GCsnap database (will be created if absent).')
    args = parser.parse_args()

    st = time.time()
    n_assemblies, n_sequences, n_unique = create_dbs(
        faa_dir       = args.faa_dir,
        summary_files = args.summary_files,
        out_dir       = args.out_dir,
        n_processes   = args.processes,
        gff_dir       = args.gff_dir,
    )

    elapsed_time = time.time() - st
    elapsed_str  = str(datetime.timedelta(seconds=round(elapsed_time)))

    print(f'{n_assemblies:,} assemblies with {n_sequences:,} sequences '
          f'({n_unique:,} unique) done in {elapsed_str}')

    write_readme(
        out_dir       = args.out_dir,
        faa_dir       = args.faa_dir,
        summary_files = args.summary_files,
        n_assemblies  = n_assemblies,
        n_sequences   = n_sequences,
        n_unique      = n_unique,
        elapsed       = elapsed_str,
    )
