"""
"Most conserved genomic context" Bokeh panel.

Public functions
----------------
find_most_common_context  – pure computation, returns MostCommonContext
most_common_panel         – builds the Bokeh figure
"""
from __future__ import annotations

import statistics
from collections import Counter

from bokeh.plotting import figure
from bokeh.models import HoverTool, ColumnDataSource, TapTool
from bokeh.models.callbacks import OpenURL

from gcsnap.visual.types import FamilyColorMap, MostCommonContext


# ── computation ───────────────────────────────────────────────────────────────

def find_most_common_context(
    operons:         dict,
    syntenies:       dict,
    families_summary: dict,
) -> MostCommonContext:
    """
    Find the most frequent genomic context arrangement across all operons.

    For each gene position, we count which family appears most often,
    then record average start/end coordinates and directions.
    """
    # Collect all protein-family structures
    all_structures: list[list[int]] = []
    for odata in operons.values():
        all_structures.extend(odata['operon_protein_families_structure'])

    if not all_structures:
        return MostCommonContext([], [], [], [], [], [], [], [])

    # Use the most common structure length as the reference
    len_counts = Counter(len(s) for s in all_structures)
    ref_len = len_counts.most_common(1)[0][0]
    ref_structures = [s for s in all_structures if len(s) == ref_len]

    # Per position: most common family, frequency
    selected_families: list[int] = []
    family_frequencies: list[float] = []
    for pos in range(ref_len):
        families_at_pos = [s[pos] for s in ref_structures]
        most_common_fam, count = Counter(families_at_pos).most_common(1)[0]
        selected_families.append(most_common_fam)
        family_frequencies.append(count * 100.0 / len(ref_structures))

    # Collect coordinate / direction / size data from syntenies
    avg_starts:   list[float] = []
    avg_ends:     list[float] = []
    directions:   list[str]   = []
    avg_sizes:    list[float] = []
    stdev_sizes:  list[float] = []
    tm_annotations: list[str] = []

    # Pick one representative target per operon to read coordinates
    starts_per_pos:  list[list[float]] = [[] for _ in range(ref_len)]
    ends_per_pos:    list[list[float]] = [[] for _ in range(ref_len)]
    dirs_per_pos:    list[list[str]]   = [[] for _ in range(ref_len)]
    sizes_per_pos:   list[list[float]] = [[] for _ in range(ref_len)]
    tm_per_pos:      list[list[str]]   = [[] for _ in range(ref_len)]

    for odata in operons.values():
        for target in odata['target_members']:
            fg = syntenies[target]['flanking_genes']
            if len(fg['cds_codes']) != ref_len:
                continue
            for pos in range(ref_len):
                starts_per_pos[pos].append(fg['relative_starts'][pos])
                ends_per_pos[pos].append(fg['relative_ends'][pos])
                dirs_per_pos[pos].append(fg['directions'][pos])
                size = abs(fg['relative_ends'][pos] - fg['relative_starts'][pos])
                sizes_per_pos[pos].append(size)
                if 'TM_annotations' in fg:
                    tm_per_pos[pos].append(fg['TM_annotations'][pos])

    for pos in range(ref_len):
        s_list = starts_per_pos[pos]
        e_list = ends_per_pos[pos]
        d_list = dirs_per_pos[pos]
        z_list = sizes_per_pos[pos]
        tm_list = tm_per_pos[pos]

        avg_starts.append(statistics.mean(s_list) if s_list else 0.0)
        avg_ends.append(statistics.mean(e_list) if e_list else 0.0)
        directions.append(Counter(d_list).most_common(1)[0][0] if d_list else '+')
        avg_sizes.append(statistics.mean(z_list) if z_list else 0.0)
        stdev_sizes.append(statistics.stdev(z_list) if len(z_list) > 1 else 0.0)
        tm_annotations.append(Counter(tm_list).most_common(1)[0][0] if tm_list else '')

    return MostCommonContext(
        families           = selected_families,
        avg_starts         = avg_starts,
        avg_ends           = avg_ends,
        directions         = directions,
        family_frequencies = family_frequencies,
        avg_sizes          = avg_sizes,
        stdev_sizes        = stdev_sizes,
        tm_annotations     = tm_annotations,
    )


# ── Bokeh figure ──────────────────────────────────────────────────────────────

def most_common_panel(
    ctx:              MostCommonContext,
    family_colors:    FamilyColorMap,
    families_summary: dict,
    reference_family: int,
    width:            int = 2000,
    height:           int = 200,
) -> figure:
    """
    Build the "most conserved gene per position" Bokeh figure.

    The figure is a row of gene-arrow patches, one per position in the
    most-common context.  Hovering shows protein name, TM prediction, etc.
    """
    data = {
        'xs': [], 'ys': [],
        'facecolor': [], 'edgecolor': [],
        'transparency': [],
        'text_x': [], 'text_y': [],
        'protein_family': [],
        'protein_name': [],
        'relative_start': [], 'relative_end': [],
        'protein_size': [],
        'family_frequency': [],
        'found_models': [],
        'model_links': [],
    }

    for i, fam in enumerate(ctx.families):
        dx = ctx.avg_ends[i] - ctx.avg_starts[i] + 1
        direction = ctx.directions[i]

        if direction == '-':
            x_tail = ctx.avg_ends[i]
            dx = -dx
            x_head = x_tail + dx
            x_head_start = x_head + 100
        else:
            x_tail = ctx.avg_starts[i]
            x_head = x_tail + dx
            x_head_start = x_head - 100

        text_x = (x_tail + x_head_start) / 2

        fc = family_colors.get(fam, family_colors.get(0))
        transparency = ctx.family_frequencies[i] / 100.0 if fam != 0 else 0.2

        # Protein info
        fam_summary = families_summary.get(fam, {})
        protein_name   = fam_summary.get('name', 'n.a.') if fam != 0 else 'Non-conserved'
        protein_size   = (f"{round(ctx.avg_sizes[i])} (±{round(ctx.stdev_sizes[i])})"
                          if fam != 0 else 'n.a.')
        freq_str       = f"{round(ctx.family_frequencies[i], 1)}%" if fam != 0 else 'n.a.'
        rel_start      = f"{int(ctx.avg_starts[i]):,}" if fam != 0 else 'n.a.'
        rel_end        = f"{int(ctx.avg_ends[i]):,}" if fam != 0 else 'n.a.'
        model_state    = _model_state(fam_summary, fam)
        model_link     = _swiss_model_url(fam_summary)
        fam_label      = fam if (fam > 0 and fam != reference_family) else ''

        data['xs'].append([x_tail, x_tail, x_head_start, x_head, x_head_start])
        data['ys'].append([0.75, 1.25, 1.25, 1.0, 0.75])
        data['facecolor'].append(fc.hex_color)
        data['edgecolor'].append(fc.line_color)
        data['transparency'].append(transparency)
        data['text_x'].append(text_x)
        data['text_y'].append(1.25)
        data['protein_family'].append(fam_label)
        data['protein_name'].append(protein_name)
        data['relative_start'].append(rel_start)
        data['relative_end'].append(rel_end)
        data['protein_size'].append(protein_size)
        data['family_frequency'].append(freq_str)
        data['found_models'].append(model_state)
        data['model_links'].append(model_link)

    src = ColumnDataSource(data)
    fig = figure(
        width=width, height=height,
        y_range=[0, 4],
        title='Most conserved gene per position',
        toolbar_location='left',
    )

    # Draw solid patches manually (no hover) for fill_alpha support
    for i in range(len(data['xs'])):
        fig.patch(
            data['xs'][i], data['ys'][i],
            fill_color=data['facecolor'][i],
            line_color=data['edgecolor'][i],
            fill_alpha=data['transparency'][i],
            line_alpha=data['transparency'][i],
            line_width=1,
        )

    # Invisible patches for hover
    hover_renderer = fig.patches(
        'xs', 'ys',
        fill_color=None, line_color=None, line_width=0,
        source=src,
        hover_fill_color='white', hover_fill_alpha=0.5,
        hover_line_color='edgecolor',
        selection_fill_color='facecolor',
        selection_line_color='edgecolor',
        nonselection_fill_color='facecolor',
        nonselection_line_color='edgecolor',
        nonselection_fill_alpha=0.2,
    )

    fig.text('text_x', 'text_y', text='protein_family',
             text_baseline='bottom', text_align='center',
             text_font_size={'value': '6pt'}, source=src)

    tooltips = [
        ('Protein name',             '@protein_name'),
        ('Structural model',          '@found_models'),
        ('Frequency at position',     '@family_frequency'),
        ('Median size (stdev)',        '@protein_size'),
        ('Median start',              '@relative_start'),
        ('Median end',                '@relative_end'),
    ]
    fig.add_tools(HoverTool(tooltips=tooltips, renderers=[hover_renderer]))
    fig.add_tools(TapTool(callback=OpenURL(url='@model_links'), renderers=[hover_renderer]))

    # Minimal axis decoration
    fig.yaxis.ticker = [1]
    fig.yaxis.major_tick_line_color = None
    fig.yaxis.minor_tick_line_color = None
    fig.yaxis.major_label_text_color = None
    fig.yaxis.axis_line_width = 0
    fig.grid.visible = False
    fig.outline_line_width = 0

    return fig


# ── private helpers ───────────────────────────────────────────────────────────

def _tm_display(annotation: str) -> tuple[str, str]:
    """Return (tm_text shown on arrow, tm_pred_text shown in tooltip)."""
    if annotation == 'TM':
        return 'TM', 'Yes'
    elif annotation == 'SP':
        return 'SP', 'Contains signal peptide'
    elif annotation:
        return annotation, 'No'
    return '', 'n.a.'


def _model_state(fam_summary: dict, fam: int) -> str:
    if fam <= 0:
        return ''
    raw = fam_summary.get('model_state', '')
    if raw == 'Model exists':
        return 'Yes (click to view in Swiss-Model repository)'
    elif raw == 'Model does not exist':
        return 'No (click to model with Swiss-Model)'
    elif raw:
        return 'Not possible to find'
    return 'n.a.'


def _swiss_model_url(fam_summary: dict) -> str:
    if 'structure' not in fam_summary:
        return 'n.a.'
    structure = fam_summary['structure']
    if not structure:
        uniprot = fam_summary.get('uniprot_code', '')
        return f'https://swissmodel.expasy.org/repository/uniprot/{uniprot}'
    return structure
