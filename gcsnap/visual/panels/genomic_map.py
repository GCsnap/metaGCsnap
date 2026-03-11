"""
Main genomic-context Bokeh panel.

Each row is one target (ordered by the dendrogram).
Columns are the flanking genes drawn as coloured arrows.

Public function: :func:`genomic_map_panel`.
"""
from __future__ import annotations

from bokeh.plotting import figure
from bokeh.models import HoverTool, ColumnDataSource, TapTool
from bokeh.models.callbacks import OpenURL

from gcsnap.visual.types import FamilyColorMap, DendrogramResult


def genomic_map_panel(
    syntenies:         dict,
    operons:           dict,
    family_colors:     FamilyColorMap,
    dendro_result:     DendrogramResult,
    most_common_fig    = None,   # optional – shares x_range when provided
    dendro_fig         = None,   # optional – shares y_range / height when provided
    gc_legend_mode:    str = 'assembly',
    width:             int = 1800,  # fallback width when most_common_fig is None
) -> figure:
    """
    Build the genomic-context map figure.

    Rows are ordered by *dendro_result.leaf_labels* (bottom to top).
    Each gene is drawn as a pentagon arrow patch coloured by protein family.

    Parameters
    ----------
    syntenies:
        ``gc.get_syntenies()`` dict.
    operons:
        ``gc.get_selected_operons()`` dict.
    family_colors:
        Output of :func:`visual.colors.assign_family_colors`.
    dendro_result:
        Output of :func:`visual.dendrogram.build_dendrogram`.
    most_common_fig:
        When provided (non-None), its ``x_range`` and ``width`` are shared so
        panning stays synchronised with the most-common-context bar.
        Pass ``None`` (default) to let Bokeh auto-range the x-axis – used in
        the advanced per-operon view where no most-common panel is shown.
    dendro_fig:
        When provided (non-None), its ``y_range`` and ``height`` are shared so
        rows stay aligned with the dendrogram leaf dots.
    gc_legend_mode:
        What to show on the y-axis: ``'assembly'``, ``'operon'``, or a taxonomy
        rank string.  Defaults to ``'assembly'``.
    width:
        Pixel width used when *most_common_fig* is ``None``.  Ignored otherwise.
    """
    leaf_labels = dendro_result.leaf_labels
    y_positions = dendro_result.y_positions
    n = len(leaf_labels)

    y_step = (y_positions[-1] - y_positions[0]) / max(1, n - 1) if n > 1 else 1.0
    y_half = y_step / 4.0

    # Build the target → operon-type lookup once
    target_to_operon = {
        t: ot
        for ot, odata in operons.items()
        for t in odata['target_members']
    }

    data = {
        'operon':        [],
        'target_id':     [],
        'protein_id':     [],
        'assembly':      [],
        'name':          [],
        'protein_family':[],
        'relative_start':[],
        'relative_end':  [],
        'facecolor':     [],
        'edgecolor':     [],
        'xs':            [],
        'ys':            [],
        'text_x':        [],
        'text_y':        [],
        'tm_text_x':     [],
        'tm_text_y':     [],
        'tm_text':       [],
        'tm_pred_text':  [],
        'assembly_link':   [],
    }

    yticklabels: dict[int, str] = {}

    for i, target in enumerate(leaf_labels):
        curr_y   = float(y_positions[i])
        operon   = target_to_operon.get(target, 'unknown')
        synteny  = syntenies.get(target, {})
        meta     = synteny.get('assembly_metadata', {})
        fg       = synteny.get('flanking_genes', {})
        assembly = meta.get('assembly_accession', '')
        region   = meta.get('genomic_region', target)
        asm_link = meta.get('assembly_link', meta.get('assembly_url', ''))
        taxon    = synteny.get('taxonomy', {}).get('taxon_name', '')

        # y-axis label
        if gc_legend_mode == 'operon':
            yticklabels[int(curr_y)] = operon
        elif gc_legend_mode == 'assembly':
            yticklabels[int(curr_y)] = f'{assembly} | {operon}'
        else:
            yticklabels[int(curr_y)] = taxon or assembly

        cds_codes = fg.get('cds_codes', [])
        for j, gene in enumerate(cds_codes):
            fam       = fg['families'][j]
            gene_name = fg['names'][j]
            dx        = fg['relative_ends'][j] - fg['relative_starts'][j] + 1
            direction = fg['directions'][j]

            if direction == '-':
                x_tail      = fg['relative_ends'][j]
                dx          = -dx
                x_head      = x_tail + dx
                x_head_st   = x_head + 100
            else:
                x_tail      = fg['relative_starts'][j]
                x_head      = x_tail + dx
                x_head_st   = x_head - 100

            text_x = (x_tail + x_head_st) / 2.0

            fc = family_colors.get(fam, family_colors.get(0))

            tm_text    = ''
            tm_pred    = 'n.a.'
            if 'TM_annotations' in fg:
                ann = fg['TM_annotations'][j]
                if ann == 'TM':
                    tm_text, tm_pred = 'TM', 'Yes'
                elif ann == 'SP':
                    tm_text, tm_pred = 'SP', 'Contains signal peptide'
                else:
                    tm_pred = 'No'

            fam_label = fam if (fam > 0) else ''

            data['operon'].append(operon)
            data['target_id'].append(region)
            data['protein_id'].append(gene)
            data['assembly'].append(assembly)
            data['name'].append(gene_name)
            data['protein_family'].append(fam_label)
            data['relative_start'].append(f'{x_tail:,}')
            data['relative_end'].append(f'{x_head:,}')
            data['facecolor'].append(fc.hex_color)
            data['edgecolor'].append(fc.line_color)
            data['xs'].append([x_tail, x_tail, x_head_st, x_head, x_head_st])
            data['ys'].append([curr_y - y_half, curr_y + y_half,
                               curr_y + y_half, curr_y, curr_y - y_half])
            data['text_x'].append(text_x)
            data['text_y'].append(curr_y + y_half)
            data['tm_text_x'].append(text_x)
            data['tm_text_y'].append(curr_y)
            data['tm_text'].append(tm_text)
            data['tm_pred_text'].append(tm_pred)
            data['assembly_link'].append(asm_link)

    src = ColumnDataSource(data)

    # Resolve width / x_range / height / y_range from companion figures
    fig_width   = most_common_fig.width   if most_common_fig is not None else width
    fig_x_range = most_common_fig.x_range if most_common_fig is not None else None
    fig_height  = dendro_fig.height       if dendro_fig       is not None else max(300, n * 30)
    fig_y_range = dendro_fig.y_range      if dendro_fig       is not None else None

    fig_kwargs: dict = dict(
        width            = fig_width,
        height           = fig_height,
        toolbar_location = 'left',
        title            = 'Representative genomic contexts (hover for details)',
    )
    if fig_x_range is not None:
        fig_kwargs['x_range'] = fig_x_range
    if fig_y_range is not None:
        fig_kwargs['y_range'] = fig_y_range

    fig = figure(**fig_kwargs)

    # Draw solid patches (non-interactive) for colour/alpha
    for i in range(len(data['xs'])):
        fig.patch(
            data['xs'][i], data['ys'][i],
            fill_color=data['facecolor'][i],
            line_color=data['edgecolor'][i],
            line_width=1,
        )

    # Invisible patches layer for hover + tap
    fig.patches(
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

    tooltips = [
        ('GC type',            '@operon'),
        ('Input ID',           '@target_id'),
        ('Protein ID',         '@protein_id'),
        ('Assembly',           '@assembly'),
        ('Gene start',         '@relative_start'),
        ('Gene end',           '@relative_end'),
        ('Protein name',       '@name'),
        ('Family code',        '@protein_family'),
        ('Membrane protein',   '@tm_pred_text'),
    ]
    fig.add_tools(HoverTool(tooltips=tooltips))
    fig.add_tools(TapTool(callback=OpenURL(url='@assembly_link')))

    # y-axis labels
    fig.yaxis.ticker = list(yticklabels.keys())
    fig.yaxis.major_label_overrides = {k: v for k, v in yticklabels.items()}
    fig.yaxis.major_tick_line_color = None
    fig.yaxis.axis_line_width = 0
    fig.xaxis.axis_label = 'Position relative to target (bp)'
    fig.grid.visible = False
    fig.outline_line_width = 0

    return fig
