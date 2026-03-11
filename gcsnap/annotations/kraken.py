#!/usr/bin/env python3
"""
Standalone Kraken2 taxonomic assignment of metagenomic contigs.

Reimplements the Kraken2Taxonomy class from annotations/taxonomy.py
without any dependency on the allSnap internal package structure.

Usage
-----
    python Kraken2Taxonomy.py \\
        --fasta-list  /path/to/fasta_paths.txt \\
        --db          /path/to/kraken2_db \\
        --out-dir     /path/to/output_dir \\
        [--threads    8]

Input files
-----------
--fasta-list  : plain-text file, one gzipped contig FASTA path per line
--db          : directory that is a valid Kraken2 database
                (must contain ktaxonomy.tsv for lineage reconstruction)

Outputs (written to --out-dir)
-------------------------------
contigs.fna.gz          merged query FASTA (all input FASTAs concatenated)
kraken2_output.tsv      raw Kraken2 output, then enriched with lineage columns
kraken2_report.txt      Kraken2 summary report
taxonomy.json           nested taxonomy tree  { "root": { ... } }
distance_matrix.csv.gz  all-vs-all pairwise taxonomic distance matrix
"""

import os
import subprocess
import gzip
import json
import argparse
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Taxonomic rank constants
# (non-binning / GTDB mode — mirrors helpers.tax_ranks('taxonomy'))
# ---------------------------------------------------------------------------

TAX_RANKS = [
    'root', 'domain', 'phylum', 'class',
    'order', 'family', 'genus', 'species', 'strain',
]

# Maps the single-letter prefix used in lineage strings to its rank name.
# E.g.  "d__2157"  →  prefix='d'  →  rank='domain'
TAX_RANKS_DICT = {
    'r': 'root',
    'd': 'domain',
    'p': 'phylum',
    'c': 'class',
    'o': 'order',
    'f': 'family',
    'g': 'genus',
    's': 'species',
}


# ---------------------------------------------------------------------------
# Helper functions (inlined from annotations/taxtree.py)
# ---------------------------------------------------------------------------

def _parse_taxonomy_tree(node, current_path, all_targets_list, target_to_path_map):
    """Recursively traverse the taxonomy JSON tree, collecting targets and paths."""
    if 'target_members' in node:
        for target_id in node['target_members']:
            if target_id not in target_to_path_map:   # keep first occurrence only
                all_targets_list.append(target_id)
                target_to_path_map[target_id] = list(current_path)

    for key, sub_node in node.items():
        if key not in {'target_members', 'ncbi_codes'} and isinstance(sub_node, dict):
            _parse_taxonomy_tree(
                sub_node, current_path + [key],
                all_targets_list, target_to_path_map,
            )


def _taxonomic_distance(path_a, path_b):
    """
    Path distance between two taxonomic paths.

    Finds the Lowest Common Ancestor (LCA) depth, then sums the
    steps from each node up to the LCA.
    """
    lca_depth = 0
    for i in range(min(len(path_a), len(path_b))):
        if path_a[i] == path_b[i]:
            lca_depth = i + 1
        else:
            break
    return (len(path_a) - lca_depth) + (len(path_b) - lca_depth)


def _get_dist(taxonomy_data):
    """
    Build an all-vs-all taxonomic distance matrix from a nested taxonomy tree.

    Args:
        taxonomy_data (dict): must contain a top-level 'root' key
                              (the JSON written by _prepare_taxonomy_tree).

    Returns:
        pd.DataFrame: symmetric distance matrix indexed by target ID,
                      or None if no targets are found.
    """
    all_targets_list  = []
    target_to_path_map = OrderedDict()

    print("Parsing taxonomic tree and mapping targets to paths…")
    _parse_taxonomy_tree(
        taxonomy_data['root'], ['root'],
        all_targets_list, target_to_path_map,
    )

    n = len(all_targets_list)
    if n == 0:
        print("No targets found in 'target_members' keys.")
        return None

    print(f"  Found {n} unique targets.")
    print("Building all-vs-all distance matrix…")

    dist_mat = np.zeros((n, n), dtype=int)
    paths    = [target_to_path_map[t] for t in all_targets_list]

    for i in range(n):
        for j in range(i + 1, n):
            d = _taxonomic_distance(paths[i], paths[j])
            dist_mat[i, j] = d
            dist_mat[j, i] = d

    print("  Matrix complete.")
    df = pd.DataFrame(dist_mat, index=all_targets_list, columns=all_targets_list)
    df.index = df.index.rename('target')
    return df


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class Kraken2Taxonomy:
    """
    Standalone Kraken2 taxonomic assignment pipeline.

    Steps
    -----
    1. Merge all input gzipped FASTA files → contigs.fna.gz
    2. Run Kraken2 to assign taxonomy         → kraken2_output.tsv
    3. Curate output: reconstruct full lineages from ktaxonomy.tsv
    4. Build a nested taxonomy tree           → taxonomy.json
    5. Compute pairwise taxonomic distances   → distance_matrix.csv.gz

    All steps are idempotent: if the output file already exists the
    step is skipped, allowing the pipeline to be resumed after a crash.
    """

    def __init__(
        self,
        fasta_list_file: str,
        kraken2_db: str,
        out_dir: str,
        threads: int = 4,
    ):
        """
        Args:
            fasta_list_file : .txt file with one gzipped FASTA path per line.
            kraken2_db      : Path to the Kraken2 reference database directory.
            out_dir         : Output directory (created if necessary).
            threads         : CPU threads passed to Kraken2.
        """
        with open(fasta_list_file, 'r') as fh:
            self.contigs_fasta = [line.strip() for line in fh if line.strip()]

        self.reference_db = Path(kraken2_db)
        self.out_dir      = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.query_fasta          = self.out_dir / "contigs.fna.gz"
        self.output               = self.out_dir / "kraken2_output.tsv"
        self.report               = self.out_dir / "kraken2_report.txt"
        self.taxonomy_json        = self.out_dir / "taxonomy.json"
        self.distance_matrix_file = self.out_dir / "distance_matrix.csv.gz"

        self.threads = threads

        # set after run()
        self.tree            = None
        self.distance_matrix = None

        print(f"Kraken2Taxonomy initialised.")
        print(f"  FASTA files : {len(self.contigs_fasta)}")
        print(f"  Database    : {self.reference_db}")
        print(f"  Output dir  : {self.out_dir}")
        print(f"  Threads     : {self.threads}")

    # ------------------------------------------------------------------
    # Step 1 — merge FASTA files
    # ------------------------------------------------------------------

    def _prepare_query(self):
        """Concatenate all gzipped contig FASTAs into a single gzipped file."""
        print(f"Merging {len(self.contigs_fasta)} FASTA file(s) → {self.query_fasta}")
        with gzip.open(self.query_fasta, "wt") as outfile:
            for fasta_file in self.contigs_fasta:
                try:
                    with gzip.open(fasta_file, "rt") as infile:
                        for line in infile:
                            outfile.write(line)
                except Exception as exc:
                    print(f"  Warning: skipping {fasta_file} — {exc}")

    # ------------------------------------------------------------------
    # Step 2 — run Kraken2
    # ------------------------------------------------------------------

    def _run_command(self, command: list):
        """Run a subprocess, raise RuntimeError on non-zero exit."""
        try:
            subprocess.run(command, check=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"Command failed: {' '.join(str(c) for c in command)}\n{exc}"
            )

    def _assign_taxonomy(self):
        """Run Kraken2 against the reference database."""
        print("Running Kraken2…")
        command = [
            "kraken2",
            "--db",      str(self.reference_db),
            "--output",  str(self.output),
            "--report",  str(self.report),
            "--threads", str(self.threads),
            "--use-names",
            str(self.query_fasta),
        ]
        self._run_command(command)
        print(f"  Output  : {self.output}")
        print(f"  Report  : {self.report}")

    # ------------------------------------------------------------------
    # Step 3 — curate taxonomy (lineage reconstruction)
    # ------------------------------------------------------------------

    def _curate_taxonomy(self):
        """
        Enrich the raw Kraken2 output with reconstructed lineage columns.

        Reads ktaxonomy.tsv from the Kraken2 database directory and
        traverses the taxonomy tree upward from each assigned taxon to
        build the full lineage string (e.g. "r__1;d__2157;p__1224;…").

        The output TSV is overwritten in-place with the additional columns:
            taxon_name, taxon_id, rank, lineage_names, lineage_ids
        """
        # Rank codes used in Kraken2 output → GTDB-style prefixes / friendly names
        rank_prefix = {
            'R':  'r__', 'R1': 'd__',
            'P':  'p__', 'C':  'c__',
            'O':  'o__', 'F':  'f__',
            'G':  'g__', 'S':  's__',
        }
        rank_name = {
            'R':  'root',   'R1': 'domain',
            'P':  'phylum', 'C':  'class',
            'O':  'order',  'F':  'family',
            'G':  'genus',  'S':  'species',
        }

        def reconstruct_lineage(taxon_id):
            if taxon_id not in tax_lookup:
                return None
            lineage_names, lineage_ids = [], []
            current = taxon_id
            while True:
                node = tax_lookup.get(current)
                if node is None:
                    break
                rk     = node['rank']
                prefix = rank_prefix.get(rk, 'x__')
                lineage_names.append(f"{prefix}{node['taxid']}")
                lineage_ids.append(str(current))
                if current == node['parent']:   # reached root
                    break
                current = node['parent']
            lineage_names.reverse()
            lineage_ids.reverse()
            final_rank = rank_name.get(tax_lookup[taxon_id]['rank'], 'unknown')
            return ';'.join(lineage_names), ';'.join(lineage_ids), final_rank

        print(f"Curating taxonomy from {self.output}…")
        K = pd.read_csv(
            self.output, sep='\t', header=None,
            names=['classified', 'genomic_region', 'classification', 'length', 'matches'],
            dtype=str,
        )
        K[['taxon_name', 'taxon_id']] = K['classification'].str.extract(
            r'^(.*) \(taxid (\d+)\)$'
        )

        unique_taxon_ids = [
            t for t in K['taxon_id'].dropna().unique() if t != '0'
        ]

        # Load the taxonomy lookup table bundled with the Kraken2 DB
        tax_file = self.reference_db / 'ktaxonomy.tsv'
        if not tax_file.exists():
            raise FileNotFoundError(
                f"ktaxonomy.tsv not found in database directory: {self.reference_db}"
            )
        tax_df = pd.read_csv(
            tax_file, sep='\t', header=None,
            usecols=[0, 2, 4, 6, 8],
            names=['node', 'parent', 'rank', 'rankid', 'taxid'],
            dtype=str,
        )
        del tax_df['rankid']
        tax_lookup = tax_df.set_index('node').to_dict('index')

        # Reconstruct lineage for every unique taxon ID
        results = []
        for tax_id in unique_taxon_ids:
            lineage = reconstruct_lineage(tax_id)
            if lineage:
                ln, li, fr = lineage
                results.append((tax_id, ln, li, fr))

        lineage_df = (
            pd.DataFrame(results, columns=['taxon_id', 'lineage_names', 'lineage_ids', 'rank'])
            .set_index('taxon_id')
        )

        for col in ['lineage_names', 'lineage_ids', 'rank']:
            K[col] = K['taxon_id'].map(lineage_df[col].to_dict())

        cols_out = [
            'classified', 'genomic_region', 'length',
            'taxon_name', 'taxon_id', 'rank',
            'lineage_names', 'lineage_ids', 'matches',
        ]
        K['length'] = K['length'].astype(int)
        K.sort_values(by='length', ascending=False, inplace=True)

        K['taxon_name'].replace('unclassified', 'unclassified contig', inplace=True)
        K['taxon_name'].replace('root',         'unclassified contig', inplace=True)
        K['rank'].fillna('root',       inplace=True)
        K['lineage_names'].fillna('r__root', inplace=True)
        K['lineage_ids'].fillna('1',       inplace=True)
        K['taxon_id'].replace(0, 1,        inplace=True)

        K.to_csv(self.output, sep='\t', index=False, columns=cols_out)
        print(f"  Curated output written to {self.output}")

    # ------------------------------------------------------------------
    # Step 4 — build taxonomy tree
    # ------------------------------------------------------------------

    def _prepare_taxonomy_tree(self):
        """
        Build a nested taxonomy tree from the curated Kraken2 output.

        Each leaf node contains 'target_members' (list of contig IDs).
        The tree is wrapped under a top-level 'root' key so that
        _get_dist() can traverse it with a consistent entry point.

        Writes: taxonomy.json
        """
        print("Building taxonomy tree…")
        tree = {}

        try:
            tax = pd.read_csv(
                self.output, sep='\t',
                usecols=['genomic_region', 'lineage_names'],
                dtype=str,
            )
            tax = tax.dropna(subset=['lineage_names'])
            lineage_to_contigs = (
                tax.groupby('lineage_names')['genomic_region'].apply(list)
            )
        except Exception as exc:
            print(f"  Error reading {self.output}: {exc}")
            self.tree = {'root': tree}
            with open(self.taxonomy_json, 'w') as fh:
                json.dump(self.tree, fh, indent=2)
            return

        for lineage_str, contig_list in lineage_to_contigs.items():
            # In standalone mode, targets ARE the contig IDs.
            targets = list(contig_list)

            # Parse the lineage string into a {rank_name: value} map.
            # E.g.  "r__1;d__2157;p__1224"  →  {'root':'1','domain':'2157',...}
            lineage_map = {}
            for part in lineage_str.split(';'):
                if '__' in part:
                    prefix, name = part.split('__', 1)
                    if prefix in TAX_RANKS_DICT:
                        lineage_map[TAX_RANKS_DICT[prefix]] = name

            # Walk down the tree, creating nested dicts as needed.
            current_level    = tree
            last_known_level = None
            for rank in TAX_RANKS:
                if rank in lineage_map:
                    taxon_name       = lineage_map[rank]
                    current_level    = current_level.setdefault(taxon_name, {})
                    last_known_level = current_level
                else:
                    break

            # Attach targets to the deepest resolved node.
            if last_known_level is not None:
                last_known_level.setdefault('target_members', []).extend(targets)
                last_known_level.setdefault('ncbi_codes',     []).extend(targets)
            else:
                unclassified = tree.setdefault('unclassified', {})
                unclassified.setdefault('target_members', []).extend(targets)
                unclassified.setdefault('ncbi_codes',     []).extend(targets)

        # Wrap with 'root' so that _get_dist() can find taxonomy_data['root'].
        self.tree = {'root': tree}
        with open(self.taxonomy_json, 'w') as fh:
            json.dump(self.tree, fh, indent=2)
        print(f"  Taxonomy tree written to {self.taxonomy_json}")

    # ------------------------------------------------------------------
    # Step 5 — distance matrix
    # ------------------------------------------------------------------

    def _compute_distance_matrix(self):
        """Compute and save the all-vs-all taxonomic distance matrix."""
        print("Computing pairwise taxonomic distance matrix…")
        dm = _get_dist(self.tree)
        if dm is not None:
            self.distance_matrix = dm
            self.distance_matrix.to_csv(
                self.distance_matrix_file, compression='gzip'
            )
            print(f"  Distance matrix written to {self.distance_matrix_file}")

    # ------------------------------------------------------------------
    # Public getters
    # ------------------------------------------------------------------

    def get_tree(self):
        """Return the taxonomy tree dict (populated after run())."""
        return self.tree

    def get_distance_matrix(self):
        """Return the distance DataFrame (populated after run())."""
        return self.distance_matrix

    # ------------------------------------------------------------------
    # Pipeline entry point
    # ------------------------------------------------------------------

    def run(self):
        """
        Execute the full pipeline.

        Every step checks whether its output file already exists and
        skips the computation if so (idempotent / resumable).
        """
        print("\n=== Kraken2Taxonomy pipeline ===")

        # --- Step 1 + 2: Run Kraken2 ---
        if not self.output.exists():
            self._prepare_query()
            self._assign_taxonomy()
        else:
            print(f"Kraken2 output already present — skipping: {self.output}")

        # --- Step 3: Curate taxonomy ---
        if self.output.exists():
            header = pd.read_csv(self.output, sep='\t', nrows=0)
            if 'taxon_name' in header.columns:
                print("Lineage columns already present — skipping curation.")
            else:
                self._curate_taxonomy()

        # --- Step 4: Build taxonomy tree ---
        if not self.taxonomy_json.exists():
            self._prepare_taxonomy_tree()
        else:
            print(f"Taxonomy JSON already present — loading: {self.taxonomy_json}")
            with open(self.taxonomy_json, 'r') as fh:
                self.tree = json.load(fh)

        # --- Step 5: Distance matrix ---
        if not self.distance_matrix_file.exists():
            self._compute_distance_matrix()
        else:
            print(f"Distance matrix already present — skipping: {self.distance_matrix_file}")
            self.distance_matrix = pd.read_csv(
                self.distance_matrix_file, compression='gzip', index_col='target'
            )

        print("=== Done ===\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Standalone Kraken2 taxonomic assignment of metagenomic contigs.\n"
            "Reimplements Kraken2Taxonomy from annotations/taxonomy.py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--fasta-list", required=True, metavar="FILE",
        help="Plain-text file listing gzipped contig FASTA paths, one per line.",
    )
    parser.add_argument(
        "--db", required=True, metavar="PATH",
        help="Kraken2 reference database directory (must contain ktaxonomy.tsv).",
    )
    parser.add_argument(
        "--out-dir", required=True, metavar="PATH",
        help="Output directory (created if it does not exist).",
    )
    parser.add_argument(
        "--threads", type=int, default=4, metavar="N",
        help="Number of CPU threads for Kraken2 (default: 4).",
    )
    return parser


def main():
    args = _build_parser().parse_args()

    runner = Kraken2Taxonomy(
        fasta_list_file=args.fasta_list,
        kraken2_db=args.db,
        out_dir=args.out_dir,
        threads=args.threads,
    )
    runner.run()


if __name__ == "__main__":
    main()
