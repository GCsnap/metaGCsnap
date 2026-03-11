"""
Protein-family colour legend Bokeh panel.

Public function: :func:`legend_panel`.
"""
from __future__ import annotations

from bokeh.plotting import figure
from bokeh.models import HoverTool, ColumnDataSource, TapTool
from bokeh.models.callbacks import OpenURL

from gcsnap.visual.types import FamilyColorMap
from gcsnap.visual.panels.most_common import _swiss_model_url, _model_state


def legend_panel(
    family_colors:    FamilyColorMap,
    families_summary: dict,
    reference_family: int,
    genomic_map_fig,                 # shares height
    rescale_height:   bool = False,
) -> figure:
    """
    Build the colour-legend figure that sits to the right of the genomic map.

    Each row is one protein family: a coloured rectangle + family name +
    Swiss-Model link on click.
    """
    families = sorted(family_colors.keys())
    n = len(families)
    height = genomic_map_fig.height if not rescale_height else max(200, n * 20)

    data: dict[str, list] = {
        'xs':          [],
        'ys':          [],
        'facecolor':   [],
        'edgecolor':   [],
        'family_name': [],
        'family_code': [],
        'found_models':[],
        'model_links': [],
        'text_x':      [],
        'text_y':      [],
    }

    for i, fam in enumerate(families):
        fc     = family_colors[fam]
        y      = n - i                       # top-down order
        x_tail, dx = 0.0, 5.0
        x_head = x_tail + dx

        if fam == reference_family:
            label = f'Target: {families_summary.get(fam, {}).get("name", str(fam))}'
        elif fam == 0:
            label = 'Non-conserved gene'
        elif fam == -1:
            label = 'Pseudogene'
        else:
            label = families_summary.get(fam, {}).get('name', str(fam))

        summary = families_summary.get(fam, {})
        data['xs'].append([x_tail, x_tail, x_head, x_head, x_tail])
        data['ys'].append([y - 0.25, y + 0.25, y + 0.25, y - 0.25, y - 0.25])
        data['facecolor'].append(fc.hex_color)
        data['edgecolor'].append(fc.line_color)
        data['family_name'].append(label)
        data['family_code'].append(str(fam) if fam not in (0, -1) else '')
        data['found_models'].append(_model_state(summary, fam))
        data['model_links'].append(_swiss_model_url(summary))
        data['text_x'].append(x_head + 0.5)
        data['text_y'].append(y)

    src = ColumnDataSource(data)

    fig = figure(
        width=400,
        height=height,
        x_range=(-1, 25),
        y_range=(0, n + 1),
        toolbar_location=None,
        title='Protein family legend',
    )

    # Draw patches manually (for linestyle support)
    for i, fam in enumerate(families):
        fc = family_colors[fam]
        kwargs = {}
        if fc.line_style == ':':
            kwargs['line_dash'] = 'dotted'
        fig.patch(
            data['xs'][i], data['ys'][i],
            fill_color=fc.hex_color,
            line_color=fc.line_color,
            **kwargs,
        )

    # Invisible patches for hover / tap
    fig.patches(
        'xs', 'ys',
        fill_color=None, line_color=None, line_width=0,
        source=src,
        hover_fill_color='white', hover_fill_alpha=0.5,
        hover_line_color='edgecolor',
    )

    fig.text('text_x', 'text_y', text='family_name',
             text_align='left', text_baseline='middle',
             text_font_size={'value': '8pt'}, source=src)

    tooltips = [
        ('Family',          '@family_name (@family_code)'),
        ('Structural model','@found_models'),
    ]
    fig.add_tools(HoverTool(tooltips=tooltips))
    fig.add_tools(TapTool(callback=OpenURL(url='@model_links')))

    fig.grid.visible = False
    fig.outline_line_width = 0
    for ax in [fig.xaxis, fig.yaxis]:
        ax.ticker                = []     # remove all tick marks
        ax.major_tick_line_color = None
        ax.minor_tick_line_color = None
        ax.major_label_text_color = None
        ax.axis_line_width       = 0

    return fig
