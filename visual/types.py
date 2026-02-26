"""
Shared data structures for visual.

Every public function in this package takes and returns these types,
not opaque **kwargs dicts or objects with class-level state.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class FamilyColor:
    """Visual style for one protein family."""
    hex_color:  str   # fill colour  e.g. '#4e9af1'
    line_color: str   # border colour e.g. '#000000'
    line_style: str   # '-' solid | ':' dashed

    # convenience alias used by matplotlib arrow patches
    @property
    def face(self) -> str:
        return self.hex_color

    @property
    def edge(self) -> str:
        return self.line_color


# family_id → FamilyColor
FamilyColorMap = dict[int, FamilyColor]


@dataclass
class DendrogramResult:
    """
    Output of :func:`visual.dendrogram.build_dendrogram`.

    icoord / dcoord are the raw scipy dendrogram coordinates used to draw the
    tree lines.  leaf_labels and y_positions are the ordered list of target IDs
    together with their y-axis positions – both bottom-to-top.

    cluster_ids is a 1-based integer array (same length as leaf_labels, same
    display order) produced by cutting the linkage tree into at most 20
    top-level branches via ``scipy.cluster.hierarchy.fcluster``.  Leaves that
    belong to the same branch share the same cluster_id and therefore the same
    dot colour in the dendrogram panel.
    """
    icoord:      list[list[float]]
    dcoord:      list[list[float]]
    leaf_labels: list[str]    # target IDs in dendrogram leaf order
    y_positions: np.ndarray   # y coordinate for each leaf
    cluster_ids: np.ndarray   # branch/cluster label per leaf (1-based, display order)


@dataclass
class MostCommonContext:
    """
    The consensus genomic context across one set of operons.

    Computed by :func:`visual.panels.most_common.find_most_common_context`.
    """
    families:           list[int]
    avg_starts:         list[float]
    avg_ends:           list[float]
    directions:         list[str]
    family_frequencies: list[float]   # 0-100 %
    avg_sizes:          list[float]
    stdev_sizes:        list[float]
    tm_annotations:     list[str]     # '' when not available
