"""
Dendrogram computation.

One public function: :func:`build_dendrogram`.
"""
from __future__ import annotations

import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial import distance as sp_distance

from gcsnap.visual.types import DendrogramResult


def build_dendrogram(
    dist_matrix: np.ndarray,
    labels: list[str],
    method: str = 'average',
) -> DendrogramResult:
    """
    Compute a hierarchical dendrogram from a square distance matrix.

    Parameters
    ----------
    dist_matrix:
        Square (N × N) distance matrix.
    labels:
        Row/column labels – one per target.
    method:
        Linkage method passed to :func:`scipy.cluster.hierarchy.linkage`.

    Returns
    -------
    :class:`~visual.types.DendrogramResult`

    Notes
    -----
    ``y_positions`` uses scipy's exact leaf placement formula ``(2*i+1)*5`` so
    that scatter dots drawn at those y-values land precisely on the tree arms
    defined by ``icoord`` / ``dcoord``.
    """
    n = len(labels)
    if n < 2:
        # Degenerate case – y_positions / cluster_ids length must match leaf_labels
        return DendrogramResult(
            icoord      = [],
            dcoord      = [],
            leaf_labels = labels,
            y_positions = np.array([(2 * i + 1) * 5 for i in range(n)], dtype=float),
            cluster_ids = np.ones(n, dtype=int),
        )

    condensed = sp_distance.squareform(dist_matrix)
    linkage   = hierarchy.linkage(condensed, method=method)
    dendro    = hierarchy.dendrogram(linkage, labels=labels,
                                     get_leaves=True, no_plot=True)

    leaf_order  = dendro['leaves']           # indices into *labels*
    leaf_labels = [labels[i] for i in leaf_order]

    # y_positions: scipy places the i-th leaf (in display order) at exactly
    # (2*i + 1) * 5.  Using this exact formula keeps the scatter dots perfectly
    # aligned with the tree arms drawn from dendro['icoord'].
    y_positions = np.array([(2 * i + 1) * 5 for i in range(n)], dtype=float)

    # cluster_ids: cut the tree into at most 20 top-level branches so that each
    # branch can be assigned one of the 20 tab20 colours.  fcluster returns
    # labels in *original* order (same as `labels`); reorder to display order.
    n_cut       = min(n, 20)
    raw_clusters = hierarchy.fcluster(linkage, t=n_cut, criterion='maxclust')
    cluster_ids  = np.array([raw_clusters[i] for i in leaf_order], dtype=int)

    return DendrogramResult(
        icoord      = dendro['icoord'],
        dcoord      = dendro['dcoord'],
        leaf_labels = leaf_labels,
        y_positions = y_positions,
        cluster_ids = cluster_ids,
    )
