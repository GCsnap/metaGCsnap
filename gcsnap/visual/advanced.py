"""
Advanced interactive HTML figure – one Bokeh tab per group.

Two grouping modes are supported via the *group_by* parameter:

``group_by='operon'`` (default)
    One tab per GC type (from ``gc.get_selected_operons()``).
    Within-group dendrogram uses Jaccard distances over protein-family
    structure vectors (no external distance file required).

``group_by='metagenomic_bins'``
    One tab per SourMash DNA bin (from ``gc.selected_metagenomic_bins``).
    Within-group dendrogram uses the SourMash distance matrix
    (``gc.metagenomic_bins_distance_matrix_file``), falling back to
    Jaccard if the file is absent or a bin has fewer than two members
    present in the matrix.

Each tab layout::

    ┌────────────┬────────────────────────────────────────┐
    │ Dendrogram │  Genomic context map (one row/target)  │
    │  (300 px)  │           (1800 px)                    │
    ├────────────┴────────────────────────────────────────┤  ← total 2100 px
    │  Assembly metadata table   │  Protein families      │
    │       (1300 px)            │     (800 px)           │  ← same height
    └────────────────────────────┴────────────────────────┘

A first "Summary" tab lists all groups with member counts and their
dominant protein families.

Public function: :func:`draw_advanced`.
"""
from __future__ import annotations

import os
from collections import Counter

import numpy as np

from bokeh.plotting import output_file, save
from bokeh.models import (
    TabPanel, Tabs, Div,
    DataTable, TableColumn, ColumnDataSource, HTMLTemplateFormatter,
)
from bokeh.layouts import column, row

from gcsnap.visual.dendrogram import build_dendrogram
from gcsnap.visual.panels.most_common import _model_state
from gcsnap.visual.panels.dendrogram_panel import dendrogram_panel
from gcsnap.visual.panels.genomic_map import genomic_map_panel


# ── public entry point ────────────────────────────────────────────────────────

def draw_advanced(
    gc,
    family_colors,
    reference_family: int,
    sort_mode:        str,
    out_dir:          str,
    gc_legend_mode:   str = 'assembly',
    group_by:         str = 'operon',
    select_groups:    list | None = None,
    out_filename:     str | None = None,
) -> str:
    """
    Build and save a multi-tab interactive HTML figure.

    Parameters
    ----------
    gc:
        :class:`~gcsnap.genomic_context.GenomicContext` object.
    family_colors:
        Output of :func:`visual.colors.assign_family_colors`.
    reference_family:
        Family id of the query protein.
    sort_mode:
        Passed to :func:`dendrogram_panel` for axis labelling only.
    out_dir:
        Directory where the HTML file is written.
    gc_legend_mode:
        Y-axis label mode for the genomic map rows:
        ``'assembly'`` (default), ``'operon'``, or a taxonomy-rank string.
    group_by:
        How to partition targets into tabs.

        * ``'operon'`` (default) – one tab per GC type read from
          ``gc.get_selected_operons()``.  Dendrogram uses Jaccard over
          protein-family structure vectors.
        * ``'metagenomic_bins'`` – one tab per SourMash DNA bin read from
          ``gc.selected_metagenomic_bins``.  Dendrogram uses the SourMash
          distance matrix (``gc.metagenomic_bins_distance_matrix_file``),
          falling back to Jaccard when unavailable.
    select_groups:
        Optional list of group identifiers to include.  Pass integers
        (e.g. ``[1, 2, 5]``) or full key strings (e.g.
        ``['Metagenomic bin 00001']``).  ``None`` (default) includes all
        groups.
    out_filename:
        Output HTML file name (basename only, no directory).  When
        ``None`` (default) the name is ``advanced_{group_by}.html``
        (e.g. ``advanced_operon.html`` or
        ``advanced_metagenomic_bins.html``).

    Returns
    -------
    str
        Absolute path to the saved HTML file.
    """
    syntenies        = gc.get_syntenies()
    families_summary = gc.get_families()

    # ── select & normalise groups ─────────────────────────────────────────
    if group_by == 'metagenomic_bins':
        raw_groups = gc.selected_metagenomic_bins
        mode_label = 'Metagenomic bin'
        page_title = 'Genomic context – per metagenomic bin'
    else:  # 'operon'
        raw_groups = gc.get_selected_operons()
        mode_label = 'GC type'
        page_title = 'Genomic context – per operon type'

    # Default filename: advanced_{group_by}.html (underscores, no spaces)
    fname = out_filename or f'advanced_{group_by.replace(" ", "_")}.html'

    # Ensure every group has 'operon_protein_families_structure'.
    # Bin groups only carry 'target_members'; we reconstruct structures from
    # syntenies so that _family_table and the Jaccard fallback always work.
    groups = _normalize_groups(raw_groups, syntenies)

    # Apply optional group filter (integers or full key strings).
    if select_groups is not None:
        groups = _filter_groups(groups, select_groups)

    # Assign one fixed tab20 colour per group (by sorted key index).
    # Every leaf dot in a tab gets that tab's group colour, making the
    # colour encode group identity consistently across all tabs.
    group_colors = _group_color_map(groups)

    os.makedirs(out_dir, exist_ok=True)

    tabs: list[TabPanel] = [
        _summary_tab(groups, families_summary, family_colors, mode_label=mode_label),
    ]

    for group_type in sorted(groups.keys()):
        tab = _group_tab(
            group_type   = group_type,
            gdata        = groups[group_type],
            syntenies    = syntenies,
            families_summary = families_summary,
            family_colors    = family_colors,
            reference_family = reference_family,
            sort_mode        = sort_mode,
            gc_legend_mode   = gc_legend_mode,
            group_by         = group_by,
            gc               = gc,
            group_color      = group_colors[group_type],
        )
        tabs.append(tab)

    all_tabs = Tabs(tabs=tabs)

    out_path = os.path.join(out_dir, fname)
    output_file(out_path, title=page_title)
    save(all_tabs)
    display_path = os.path.join(os.path.basename(os.path.dirname(out_path)), os.path.basename(out_path))

    return out_path


# ── summary tab ───────────────────────────────────────────────────────────────

def _summary_tab(
    groups:           dict,
    families_summary: dict,
    family_colors,
    mode_label:       str = 'GC type',
) -> TabPanel:
    """First tab: overview table of all groups (operon types or bins)."""
    data: dict[str, list] = {
        'group_type':   [],
        'n_members':    [],
        'top_families': [],
        'color':        [],
    }

    for gt in sorted(groups.keys()):
        gdata = groups[gt]
        n     = len(gdata['target_members'])

        fam_counts: Counter = Counter(
            fam
            for struct in gdata['operon_protein_families_structure']
            for fam in set(struct)
            if fam > 0
        )
        top_fam = fam_counts.most_common(1)[0][0] if fam_counts else 0
        fc      = family_colors.get(top_fam, family_colors.get(0))

        top_names = ', '.join(
            families_summary.get(f, {}).get('name', str(f))
            for f, _ in fam_counts.most_common(3)
        ) or '—'

        data['group_type'].append(gt)
        data['n_members'].append(n)
        data['top_families'].append(top_names)
        data['color'].append(fc.hex_color)

    color_fmt = HTMLTemplateFormatter(
        template=(
            '<span style="background:<%= value %>;display:block;'
            'width:90%;height:18px;border:1px solid #666;"></span>'
        )
    )
    columns = [
        TableColumn(field='group_type',   title=mode_label,              width=160),
        TableColumn(field='n_members',    title='# members',             width=80),
        TableColumn(field='color',        title='Top-family colour',     width=100,
                    formatter=color_fmt),
        TableColumn(field='top_families', title='Top 3 families (by member count)',
                    width=500),
    ]
    table = DataTable(
        source  = ColumnDataSource(data),
        columns = columns,
        width   = 860,
        height  = max(150, min(600, 50 + 28 * len(data['group_type']))),
    )
    div = Div(
        text=(
            f'<h2 style="margin-bottom:4px">'
            f'Genomic context – per {mode_label.lower()}</h2>'
            '<p style="color:#555">Select a tab to explore the detailed genomic '
            f'context of each {mode_label.lower()}.</p>'
        )
    )
    return TabPanel(child=column(div, table), title='Summary')


# Fixed pixel widths for the two panel rows.
# Top row:    _DEN_W (dendrogram) + _GC_W (genomic map) = _TOTAL_W
# Bottom row: _ASM_W (assembly)   + _FAM_W (families)   = _TOTAL_W
_DEN_W   = 300
_GC_W    = 1800
_TOTAL_W = _DEN_W + _GC_W   # 2100 px
_ASM_W   = 1300              # assembly metadata table
_FAM_W   = _TOTAL_W - _ASM_W  # 800 px – protein family table


# ── per-group tab ─────────────────────────────────────────────────────────────

def _group_tab(
    group_type:      str,
    gdata:           dict,
    syntenies:       dict,
    families_summary: dict,
    family_colors,
    reference_family: int,
    sort_mode:        str,
    gc_legend_mode:   str,
    group_by:         str,
    gc,
    group_color:      str | None = None,
) -> TabPanel:
    """
    Build one Bokeh tab for a single group (operon type or metagenomic bin).

    Layout
    ------
    ::

        title_div
        ┌────────────┬────────────────────────────────────────┐
        │ Dendrogram │  Genomic context map                   │
        │ (_DEN_W)   │  (_GC_W)          y_range shared →    │
        ├────────────┴──────────────────┬─────────────────────┤  ← _TOTAL_W px wide
        │  Assembly metadata (_ASM_W)   │  Families (_FAM_W)  │
        │  (one row per target)         │  (one row/family)   │  ← shared height
        └───────────────────────────────┴─────────────────────┘

    For ``group_by='operon'`` the dendrogram uses Jaccard distances over
    protein-family structure vectors.  For ``group_by='metagenomic_bins'``
    it uses the pre-computed SourMash distance matrix (falling back to
    Jaccard when fewer than 2 members share an entry in that matrix).
    """
    members        = gdata['target_members']
    n              = len(members)

    curr_syntenies = {t: syntenies[t] for t in members if t in syntenies}
    curr_operons   = {group_type: gdata}

    # ── within-group distance matrix + dendrogram ─────────────────────────
    dist_matrix, labels = _within_group_distance_matrix(gdata, gc, group_by)
    dendro = build_dendrogram(dist_matrix, labels)

    # ── dendrogram panel (left, _DEN_W px wide) ───────────────────────────
    # All leaves in this tab belong to the same group → uniform group colour.
    target_colors = (
        {t: group_color for t in dendro.leaf_labels}
        if group_color is not None else None
    )
    den_fig, _ = dendrogram_panel(
        dendro, curr_syntenies, curr_operons, sort_mode, family_colors,
        target_colors=target_colors,
    )

    # ── genomic map panel (right, _GC_W px wide, no most-common bar) ──────
    gc_fig = genomic_map_panel(
        curr_syntenies, curr_operons, family_colors, dendro,
        most_common_fig = None,
        dendro_fig      = den_fig,
        gc_legend_mode  = gc_legend_mode,
        width           = _GC_W,
    )

    # ── shared table height ────────────────────────────────────────────────
    n_fam_rows = len({
        fam
        for struct in gdata['operon_protein_families_structure']
        for fam in set(struct)
        if fam > 0
    })
    asm_h        = max(100, min(500, 50 + 28 * n))
    fam_h        = max(100, min(500, 50 + 28 * n_fam_rows))
    table_height = max(asm_h, fam_h)

    # ── assembly metadata table (left, _ASM_W px) ─────────────────────────
    asm_div, asm_table = _assembly_table(
        gdata, curr_syntenies,
        width=_ASM_W, height=table_height,
    )

    # ── protein family table (right, _FAM_W px) ───────────────────────────
    fam_div, fam_table = _family_table(
        gdata, families_summary, family_colors, n,
        width=_FAM_W, height=table_height,
    )

    # ── header div ────────────────────────────────────────────────────────
    title_div = Div(
        text=(
            f'<h3 style="margin:4px 0">{group_type}</h3>'
            f'<p style="color:#555;margin:0">{n} member(s)</p>'
        )
    )

    layout = column(
        title_div,
        row(den_fig, gc_fig),
        row(
            column(asm_div, asm_table),
            column(fam_div, fam_table),
        ),
    )
    return TabPanel(child=layout, title=group_type)


# ── helpers ───────────────────────────────────────────────────────────────────

def _group_color_map(groups: dict) -> dict[str, str]:
    """
    Assign one tab20 colour to each group, indexed by sorted key position.

    ``groups`` is the already-normalised dict (after filtering).  Colours
    are assigned in sorted key order so the same group always gets the same
    colour regardless of the order groups happen to be iterated.

    Returns ``{group_key: hex_color}``.
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mc

    cmap = plt.get_cmap('tab20')
    return {
        key: mc.to_hex(cmap((i % 20) / 20))
        for i, key in enumerate(sorted(groups.keys()))
    }


def _normalize_groups(groups: dict, syntenies: dict) -> dict:
    """
    Ensure every group dict has ``'operon_protein_families_structure'``.

    Operon groups (from ``selected_operons``) already carry this field.
    Bin groups (from ``selected_metagenomic_bins``) only have
    ``'target_members'``; here we reconstruct the family-structure lists
    from *syntenies* so that :func:`_family_table` and the Jaccard fallback
    in :func:`_within_group_distance_matrix` always work correctly.

    Targets absent from *syntenies* are silently dropped from the group.
    """
    normalized: dict = {}
    for label, gdata in groups.items():
        if 'operon_protein_families_structure' in gdata:
            normalized[label] = gdata        # already normalised
        else:
            valid_members = [t for t in gdata['target_members'] if t in syntenies]
            structures    = [
                syntenies[t]['flanking_genes']['families']
                for t in valid_members
            ]
            entry = dict(gdata)
            entry['target_members']                    = valid_members
            entry['operon_protein_families_structure'] = structures
            normalized[label] = entry
    return normalized


def _filter_groups(groups: dict, select_groups: list) -> dict:
    """
    Return the subset of *groups* requested by *select_groups*.

    Each element of *select_groups* can be either:

    * **int** – matched by the zero-padded number that ends the key, e.g.
      ``1`` matches ``'Metagenomic bin 00001'`` or ``'GC Type 00001'``.
    * **str** – matched as an exact key lookup.

    Unknown identifiers are silently skipped.
    """
    result: dict = {}
    for item in select_groups:
        if isinstance(item, int):
            suffix = f'{item:05d}'
            for key in groups:
                if key.endswith(suffix):
                    result[key] = groups[key]
        elif item in groups:
            result[item] = groups[item]
    return result


def _within_group_distance_matrix(
    gdata:    dict,
    gc,
    group_by: str,
) -> tuple[np.ndarray, list[str]]:
    """
    Pairwise distance matrix for the members of a single group.

    ``group_by='operon'``
        Jaccard over protein-family structure vectors (no external file).

    ``group_by='metagenomic_bins'``
        Reads ``gc.metagenomic_bins_distance_matrix_file`` (gzip CSV,
        ``index_col='target'``), filters to this bin's members, and falls
        back to Jaccard when the file is absent or fewer than two members
        are found in the matrix index.
    """
    if group_by == 'metagenomic_bins':
        dm_file = getattr(gc, 'metagenomic_bins_distance_matrix_file', None)
        if dm_file is not None and os.path.exists(str(dm_file)):
            import pandas as pd
            dm      = pd.read_csv(dm_file, compression='gzip', index_col='target')
            members = gdata['target_members']
            common  = [m for m in members if m in dm.index]
            if len(common) >= 2:
                dm_sub = dm.loc[common, common]
                return dm_sub.values, dm_sub.index.tolist()
        # fall through to Jaccard (file absent or < 2 members in matrix)

    return _jaccard_distance_matrix(gdata)


def _jaccard_distance_matrix(gdata: dict) -> tuple[np.ndarray, list[str]]:
    """
    Pairwise Jaccard distances from protein-family structure vectors.

    Distance = 1 - |A ∩ B| / |A ∪ B|.
    Used as the within-group dendrogram distance whenever the external
    distance-matrix file is unavailable.
    Note: many operon types have identical structures → all distances = 0 →
    flat dendrogram.  :func:`~visual.panels.dendrogram_panel.dendrogram_panel`
    handles this gracefully (grey leaf dots, sensible axis range).
    """
    members    = gdata['target_members']
    structures = gdata['operon_protein_families_structure']
    n          = len(members)
    matrix     = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            a     = set(structures[i])
            b     = set(structures[j])
            union = a | b
            d     = 1.0 - len(a & b) / len(union) if union else 0.0
            matrix[i, j] = matrix[j, i] = d

    return matrix, list(members)


def _assembly_table(
    odata:     dict,
    syntenies: dict,
    width:     int = _ASM_W,
    height:    int | None = None,
) -> tuple[Div, DataTable]:
    """
    Bokeh DataTable with assembly / taxonomy metadata for each target.

    Columns: Sequence_ID · source · target_cds · genomic_region ·
    assembly_accession · assembly_url (clickable) · target · target_source ·
    taxon_id · taxon_name · bin

    The ``bin`` column shows the integer SourMash bin assignment when
    present in the syntenies (empty string otherwise), making the column
    non-intrusive for non-metagenomics datasets.

    Parameters
    ----------
    width:
        Total widget width in pixels (default ``_ASM_W = 1300``).
    height:
        Widget height in pixels.  Auto-computed from row count when
        ``None``; pass an explicit value to match a sibling table's height.
    """
    data: dict[str, list] = {
        'sequence_id':        [],
        'source':             [],
        'target_cds':         [],
        'genomic_region':     [],
        'assembly_accession': [],
        'assembly_url':       [],
        'target':             [],
        'target_source':      [],
        'taxon_id':           [],
        'taxon_name':         [],
        'bin':                [],
    }

    for t in odata['target_members']:
        syn  = syntenies.get(t, {})
        meta = syn.get('assembly_metadata', {})
        tax  = syn.get('taxonomy', {})

        data['sequence_id'].append(t)
        data['source'].append(meta.get('source', ''))
        data['target_cds'].append(meta.get('target_cds', ''))
        data['genomic_region'].append(meta.get('genomic_region', ''))
        data['assembly_accession'].append(meta.get('assembly_accession', ''))
        data['assembly_url'].append(
            meta.get('assembly_url', meta.get('assembly_link', ''))
        )
        data['target'].append(meta.get('target', ''))
        data['target_source'].append(meta.get('target_source', ''))
        data['taxon_id'].append(str(tax.get('taxon_id', '')))
        data['taxon_name'].append(tax.get('taxon_name', ''))
        data['bin'].append(str(syn.get('bin', '')))

    url_fmt = HTMLTemplateFormatter(
        template='<a href="<%= value %>" target="_blank" title="<%= value %>">link ↗</a>'
    )
    # Column widths sum to ~1270 px – comfortably within the default _ASM_W=1300.
    columns = [
        TableColumn(field='sequence_id',        title='Sequence ID',        width=130),
        TableColumn(field='source',             title='Source',             width=55),
        TableColumn(field='target_cds',         title='Target CDS',         width=110),
        TableColumn(field='genomic_region',     title='Genomic region',     width=155),
        TableColumn(field='assembly_accession', title='Assembly accession', width=145),
        TableColumn(field='assembly_url',       title='Assembly URL',       width=75,
                    formatter=url_fmt),
        TableColumn(field='target',             title='Target',             width=115),
        TableColumn(field='target_source',      title='Target source',      width=105),
        TableColumn(field='taxon_id',           title='Taxon ID',           width=75),
        TableColumn(field='taxon_name',         title='Taxon name',         width=215),
        TableColumn(field='bin',                title='Bin',                width=45),
    ]

    n_rows       = len(odata['target_members'])
    table_height = height if height is not None else max(100, min(500, 50 + 28 * n_rows))
    table = DataTable(
        source  = ColumnDataSource(data),
        columns = columns,
        width   = width,
        height  = table_height,
    )
    div = Div(text='<b style="font-size:11pt">Assembly metadata:</b>', width=width)
    return div, table


def _family_table(
    odata:            dict,
    families_summary: dict,
    family_colors,
    n_members:        int,
    width:            int = _FAM_W,
    height:           int | None = None,
) -> tuple[Div, DataTable]:
    """
    Bokeh DataTable listing all protein families found in this group,
    sorted by frequency (descending).

    Parameters
    ----------
    width:
        Total widget width in pixels (default ``_FAM_W = 800``).
    height:
        Widget height in pixels.  Auto-computed from row count when
        ``None``; pass an explicit value to match a sibling table's height.
    """
    structures = odata['operon_protein_families_structure']

    fam_counts: Counter = Counter(
        fam
        for struct in structures
        for fam in set(struct)
        if fam > 0
    )

    data: dict[str, list] = {
        'family':    [],
        'color':     [],
        'name':      [],
        'count':     [],
        'frequency': [],
        'structure': [],
    }

    for fam, count in sorted(fam_counts.items(), key=lambda x: -x[1]):
        fc          = family_colors.get(fam, family_colors.get(0))
        fam_summary = families_summary.get(fam, {})

        data['family'].append(fam)
        data['color'].append(fc.hex_color)
        data['name'].append(fam_summary.get('name', 'n.a.'))
        data['count'].append(count)
        data['frequency'].append(round(count * 100.0 / n_members, 1))
        data['structure'].append(_model_state(fam_summary, fam))

    color_fmt = HTMLTemplateFormatter(
        template=(
            '<span style="background:<%= value %>;display:block;'
            'width:90%;height:18px;border:1px solid #aaa;"></span>'
        )
    )
    # Column widths sum to 760 px – fits within the default _FAM_W=800
    # leaving ~40 px for the DataTable's internal padding / scrollbar gutter.
    columns = [
        TableColumn(field='family',    title='Family',    width=60),
        TableColumn(field='color',     title='Colour',    width=60,
                    formatter=color_fmt),
        TableColumn(field='name',      title='Name',      width=250),
        TableColumn(field='count',     title='Count',     width=60),
        TableColumn(field='frequency', title='Freq (%)',  width=80),
        TableColumn(field='structure', title='Structure', width=190),
    ]

    n_rows       = len(data['family'])
    table_height = height if height is not None else max(100, min(500, 50 + 28 * n_rows))
    table = DataTable(
        source  = ColumnDataSource(data),
        columns = columns,
        width   = width,
        height  = table_height,
    )
    div = Div(
        text='<b style="font-size:11pt">Protein families:</b>',
        width=width,
    )
    return div, table
