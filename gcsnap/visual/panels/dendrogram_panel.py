"""
Bokeh dendrogram panel.

Public function: :func:`dendrogram_panel` → ``(figure, syn_den_data dict)``.
The returned ``syn_den_data`` dict carries the leaf y-positions and labels
that the genomic-map panel needs to line up its rows.
"""
from __future__ import annotations

import numpy as np
from bokeh.plotting import figure
from bokeh.models import HoverTool

from gcsnap.visual.types import DendrogramResult, FamilyColorMap


def dendrogram_panel(
    dendro:          DendrogramResult,
    syntenies:       dict,
    operons:         dict,
    sort_mode:       str,
    family_colors:   FamilyColorMap,
    height_factor:   float = 25.0,
    dot_size:        int   = 8,
    target_colors:   dict | None = None,
) -> tuple[figure, dict]:
    """
    Build the Bokeh dendrogram figure that sits to the left of the genomic
    context map.

    Returns
    -------
    fig : bokeh.plotting.figure
    syn_den_data : dict
        ``{'leaf_label': [...], 'y': [...]}`` – consumed by
        :func:`visual.panels.genomic_map.genomic_map_panel`.
    """
    n = len(dendro.leaf_labels)
    height = max(300, int(n * height_factor * 1.2))

    # ── coordinate ranges ────────────────────────────────────────────────
    # Horizontal dendrogram: x = distance (dcoord), y = leaf position (icoord).
    # scipy places leaf i at y = (2*i+1)*5  →  y ∈ [5, (2n-1)*5].
    # We reverse the x-axis so the root (max distance) is on the LEFT and
    # the leaves (distance = 0) are on the RIGHT, flush against the genomic map.
    y_bot = 0.0
    y_top = float((2 * n - 1) * 5 + 5)

    # Compute max merge-height.  Guard against all-identical members (Jaccard=0):
    # scipy produces valid dcoord with all heights = 0, which collapses x_range
    # to (0,0) and makes the figure invisible.  Fall back to 1.0 so tree arms
    # and leaf dots are always rendered with a sensible axis extent.
    if dendro.dcoord:
        max_dist = max(v for row in dendro.dcoord for v in row)
    else:
        max_dist = 0.0
    is_degenerate = (max_dist == 0.0)   # all members identical → flat tree
    if max_dist == 0.0:
        max_dist = 1.0

    x_left  =  max_dist * 1.15   # root side (reversed → high x appears on left)
    x_right = -max_dist * 0.05   # leaf side  (x ≈ 0)

    fig = figure(
        width=300, height=height,
        x_range=(x_left, x_right),   # reversed: root at left, leaves at right
        y_range=(y_bot, y_top),
        toolbar_location='left',
        title='Genomic context dendrogram',
    )

    # ── draw the tree arms (horizontal dendrogram) ────────────────────────
    # Each row in dcoord/icoord describes one merge step as four points:
    #   dcoord → x (distance axis)
    #   icoord → y (leaf-position axis)
    for drow, irow in zip(dendro.dcoord, dendro.icoord):
        fig.line(drow, irow, line_color='black', line_width=1)

    # ── coloured leaf dots ────────────────────────────────────────────────
    # Dots sit at x=0 (distance=0, i.e. the leaf side) and y=exact leaf pos.
    # Priority:
    #   1. target_colors provided (inherited from global interactive figure)
    #      → look up each leaf label; unknown labels fall back to grey.
    #   2. local dendrogram is degenerate (all distances = 0) → all grey.
    #   3. default → tab20 colours from local fcluster branch assignments.
    if target_colors is not None:
        leaf_colors = [
            target_colors.get(lbl, '#aaaaaa') for lbl in dendro.leaf_labels
        ]
    else:
        leaf_colors = _leaf_colors(dendro, degenerate=is_degenerate)

    leaf_source = dict(
        x=[0.0] * n,                    # distance = 0 → right edge of dendrogram
        y=list(dendro.y_positions),      # exact scipy leaf positions (2i+1)*5
        color=leaf_colors,
        leaf_label=dendro.leaf_labels,
    )

    dots = fig.scatter(
        x='x', y='y',
        fill_color='color', line_color='black',
        size=dot_size, source=leaf_source,
    )
    fig.add_tools(HoverTool(tooltips=[('Target', '@leaf_label')], renderers=[dots]))

    # ── axis / grid styling ───────────────────────────────────────────────
    fig.xaxis.major_tick_line_color = None
    fig.xaxis.minor_tick_line_color = None
    fig.xaxis.major_label_text_color = None
    fig.xaxis.axis_line_width = 0
    fig.yaxis.major_tick_line_color = None
    fig.yaxis.minor_tick_line_color = None
    fig.yaxis.major_label_text_color = None
    fig.yaxis.axis_line_width = 0
    fig.grid.visible = False
    fig.outline_line_width = 0

    syn_den_data = {
        'leaf_label': dendro.leaf_labels,
        'y':          list(dendro.y_positions),
    }
    return fig, syn_den_data


# ── internal helpers ──────────────────────────────────────────────────────────

def _leaf_colors(dendro: DendrogramResult, degenerate: bool = False) -> list[str]:
    """
    Return one hex colour string per leaf, derived purely from the dendrogram.

    When *degenerate* is True (all pairwise distances are zero, i.e. every
    member has an identical protein-family structure) all dots are grey –
    coloring by cluster would be meaningless since ``fcluster`` assigns every
    leaf to the same cluster anyway.

    Otherwise ``dendro.cluster_ids`` is a 1-based integer array (in display
    order) produced by cutting the linkage tree into at most 20 top-level
    branches with ``scipy.cluster.hierarchy.fcluster(..., criterion='maxclust')``.
    Leaves in the same branch share a tab20 colour, making cluster membership
    immediately visible.
    """
    if degenerate:
        return ['#aaaaaa'] * len(dendro.cluster_ids)
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('tab20')
    return [_to_hex(cmap(((cid - 1) % 20) / 20)) for cid in dendro.cluster_ids]


def _to_hex(rgba_tuple) -> str:
    import matplotlib.colors as mc
    return mc.to_hex(rgba_tuple)
