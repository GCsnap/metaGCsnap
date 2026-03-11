"""
Interactive (Bokeh) genomic-context HTML figure.

Assembles all panels into a single HTML file.

Layout
------
::

    ┌──────────────────────────────────────────────────────┐
    │  [Tab: Co-occurrence] [Tab: Adjacency] [Tab: Most common] │
    ├──────────────────────┬───────────────────────────────┤
    │  Dendrogram          │  Genomic context map  │ Legend│
    └──────────────────────┴───────────────────────┴───────┘

Public function: :func:`draw_interactive`.
"""
from __future__ import annotations
import os

from bokeh.plotting import output_file, save
from bokeh.layouts import gridplot
from bokeh.models import TabPanel, Tabs

from gcsnap.visual.colors import assign_family_colors
from gcsnap.visual.distance import get_distance_matrix
from gcsnap.visual.dendrogram import build_dendrogram
from gcsnap.visual.panels.most_common import find_most_common_context, most_common_panel
from gcsnap.visual.panels.dendrogram_panel import dendrogram_panel
from gcsnap.visual.panels.genomic_map import genomic_map_panel
from gcsnap.visual.panels.network import cooccurrence_panel, adjacency_panel
from gcsnap.visual.panels.legend import legend_panel


def draw_interactive(
    gc,
    family_colors,          # FamilyColorMap (bokeh colours)
    reference_family: int,
    sort_mode: str,
    out_dir: str,
    min_coocc: float = 0.0,
    gc_legend_mode: str = 'assembly',
) -> str:
    """
    Build and save the interactive HTML figure.

    Parameters
    ----------
    gc:
        :class:`~gcsnap.genomic_context.GenomicContext` object.
    family_colors:
        Output of :func:`visual.colors.assign_family_colors` with Bokeh-
        compatible colours.
    reference_family:
        Family id of the query protein.
    sort_mode:
        ``'operon'``, ``'metagenomic bins'``, or ``'taxonomy'``.
    out_dir:
        Directory where the HTML is written.
    min_coocc:
        Minimum co-occurrence weight to keep in the network (0–1).
    gc_legend_mode:
        Y-axis label mode for the genomic map (``'assembly'``, ``'operon'``, …).

    Returns
    -------
    str
        Path to the saved HTML file.
    """
    syntenies        = gc.get_syntenies()
    operons          = gc.get_selected_operons()
    families_summary = gc.get_families()

    # ── 1. Distance matrix + dendrogram ──────────────────────────────────
    dist_matrix, labels = get_distance_matrix(gc, operons, sort_mode)
    dendro = build_dendrogram(dist_matrix, labels)

    # ── 2. Most-common context (top panel) ───────────────────────────────
    ctx = find_most_common_context(operons, syntenies, families_summary)
    mc_fig = most_common_panel(
        ctx, family_colors, families_summary, reference_family,
    )

    # ── 3. Dendrogram panel (left) ────────────────────────────────────────
    den_fig, syn_den_data = dendrogram_panel(
        dendro, syntenies, operons, sort_mode, family_colors,
    )

    # ── 4. Genomic map panel (centre) ────────────────────────────────────
    gc_fig = genomic_map_panel(
        syntenies, operons, family_colors, dendro,
        most_common_fig=mc_fig,
        dendro_fig=den_fig,
        gc_legend_mode=gc_legend_mode,
    )

    # ── 5. Network panels (tabs) ─────────────────────────────────────────
    coocc_fig = cooccurrence_panel(
        operons, families_summary, family_colors, reference_family,
        mc_fig, min_coocc=min_coocc,
    )
    adj_fig = adjacency_panel(
        operons, families_summary, family_colors, reference_family,
        mc_fig, min_coocc=min_coocc,
    )

    # ── 6. Legend panel (right) ───────────────────────────────────────────
    leg_fig = legend_panel(
        family_colors, families_summary, reference_family, gc_fig,
    )

    # ── 7. Assemble layout ────────────────────────────────────────────────
    tab1 = TabPanel(child=coocc_fig, title='Gene co-occurrence network')
    tab2 = TabPanel(child=adj_fig,   title='Gene adjacency network')
    tab3 = TabPanel(child=mc_fig,    title='Most common genomic features')
    tabs = Tabs(tabs=[tab1, tab2, tab3])

    grid = gridplot(
        [[None, tabs, None],
         [den_fig, gc_fig, leg_fig]],
        merge_tools=True,
    )

    # ── 8. Save ───────────────────────────────────────────────────────────
    os.makedirs(out_dir, exist_ok=True)
    out_name = f'{sort_mode.replace(" ", "_")}_interactive.html'
    out_path = os.path.join(out_dir, out_name)
    output_file(out_path, title=f'Genomic context – {sort_mode}')
    save(grid)

    return out_path
