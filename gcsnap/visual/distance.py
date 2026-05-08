"""
Distance-matrix computation for all sort modes.

Each function returns *(matrix, labels)* where *matrix* is a square
``np.ndarray`` of pairwise distances and *labels* is the list of target IDs
(or operon-type IDs) in the same order as the matrix rows/columns.
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd


# ── helpers ──────────────────────────────────────────────────────────────────

def _read_csv_matrix(path: str | os.PathLike, index_col: str) -> pd.DataFrame:
    return pd.read_csv(path, compression='gzip', index_col=index_col)


def _filter_to_targets(dm: pd.DataFrame, targets: list[str]) -> pd.DataFrame:
    """Keep only rows/cols present in *targets* (same pattern as taxonomy matrix)."""
    common = [t for t in targets if t in dm.index]
    return dm.loc[common, common]


# ── per-mode functions ────────────────────────────────────────────────────────

def operon_distance_matrix(
    gc,
    operons: dict,
) -> tuple[np.ndarray, list[str]]:
    """
    Distance matrix for ``sort_mode='operon'``.

    Primary path
        Read the pre-computed ``gc.distance_matrix_file`` (written by the
        ``Operons`` pipeline step as a gzip CSV indexed by ``'target'``), then
        filter to the targets currently in ``gc.syntenies``.

    Fallback
        If the file is absent, compute pairwise Jaccard distances from the
        protein-family structure stored in each operon's
        ``'operon_protein_families_structure'``.
    """
    dm_file = getattr(gc, 'operon_distance_matrix_file', None)
    if dm_file is not None and os.path.exists(dm_file):
        dm = _read_csv_matrix(dm_file, index_col='target')
        dm = _filter_to_targets(dm, list(gc.syntenies.keys()))
        return dm.values, dm.index.tolist()

    # --- fallback: Jaccard over protein-family vectors ---
    print('[visual] distance_matrix_file not found – computing from family overlap.')
    labels: list[str] = []
    vectors: list[list[int]] = []
    for operon_type, odata in operons.items():
        for i, target in enumerate(odata['target_members']):
            labels.append(target)
            vectors.append(odata['operon_protein_families_structure'][i])

    n = len(labels)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = set(vectors[i]), set(vectors[j])
            union = a | b
            dist = 1 - len(a & b) / len(union) if union else 0.0
            matrix[i, j] = matrix[j, i] = dist

    return matrix, labels


def metagenomic_bins_distance_matrix(gc) -> tuple[np.ndarray, list[str]]:
    """
    Distance matrix for ``sort_mode='metagenomic bins'``.

    Reads ``gc.metagenomic_bins_distance_matrix_file``, which is already indexed by
    target (same layout as the operon distance matrix).
    """
    dm = _read_csv_matrix(gc.metagenomic_bins_distance_matrix_file, index_col='target')
    dm = _filter_to_targets(dm, list(gc.syntenies.keys()))
    return dm.values, dm.index.tolist()


def taxonomy_distance_matrix(gc) -> tuple[np.ndarray, list[str]]:
    """
    Distance matrix for ``sort_mode='taxonomy'``.

    Reads ``gc.taxonomic_distance_file`` (target-indexed gzip CSV).
    Falls back to operon distance if the file is not available.
    """
    if not hasattr(gc, 'taxonomic_distance_file') or gc.taxonomic_distance_file is None:
        return operon_distance_matrix(gc, {})
    dm = _read_csv_matrix(gc.taxonomic_distance_file, index_col='target')
    dm = _filter_to_targets(dm, list(gc.syntenies.keys()))
    return dm.values, dm.index.tolist()


def get_distance_matrix(
    gc,
    operons: dict,
    sort_mode: str,
) -> tuple[np.ndarray, list[str]]:
    """
    Dispatch to the correct distance-matrix function based on *sort_mode*.

    Parameters
    ----------
    gc:
        :class:`~gcsnap.genomic_context.GenomicContext` object.
    operons:
        Dict returned by ``gc.get_selected_operons()``.
    sort_mode:
        One of ``'operon'``, ``'metagenomic bins'``, ``'taxonomy'``.

    Returns
    -------
    matrix : np.ndarray
        Square distance matrix.
    labels : list[str]
        Row/column labels (target or contig IDs).
    """
    if sort_mode == 'operon':
        return operon_distance_matrix(gc, operons)
    elif sort_mode == 'metagenomic bins':
        return metagenomic_bins_distance_matrix(gc)
    elif sort_mode == 'taxonomy':
        return taxonomy_distance_matrix(gc)
    else:
        raise ValueError(
            f"Unknown sort_mode '{sort_mode}'. "
            "Expected 'operon', 'metagenomic bins', or 'taxonomy'."
        )
