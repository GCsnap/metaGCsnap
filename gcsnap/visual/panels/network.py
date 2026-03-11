"""
Gene co-occurrence and adjacency network Bokeh panels.

Public functions
----------------
cooccurrence_panel  – gene co-occurrence network
adjacency_panel     – gene adjacency network
empty_network_panel – placeholder when no network can be built
"""
from __future__ import annotations

import numpy as np
import networkx as nx
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from bokeh.plotting import figure
from bokeh.models import HoverTool, TapTool, ColumnDataSource
from bokeh.models.callbacks import OpenURL
from bokeh.models import MultiLine, Scatter
import networkx as nx
from bokeh.plotting import figure, from_networkx, show
from bokeh.colors import RGB

from gcsnap.visual.types import FamilyColorMap
from gcsnap.visual.panels.most_common import _swiss_model_url, _model_state, _tm_display


# ── public API ────────────────────────────────────────────────────────────────

def cooccurrence_panel(
    operons:          dict,
    families_summary: dict,
    family_colors:    FamilyColorMap,
    reference_family: int,
    most_common_fig,
    min_coocc:        float = 0.0,
) -> figure:
    """
    Build the gene co-occurrence network figure.

    Returns an empty panel with a message if co-occurrence is zero.
    """
    try:
        matrix, family_labels = _cooccurrence_matrix(operons, families_summary, min_coocc)
        graph, coords, family_labels = _build_graph(
            matrix, family_labels, family_colors, families_summary
        )
        return _network_figure(
            graph, coords, family_labels, families_summary,
            family_colors, reference_family, most_common_fig,
            title='Gene co-occurrence network',
        )
    except ValueError as exc:
        return empty_network_panel(
            most_common_fig,
            message=f'No co-occurrence detected ({exc})',
        )


def adjacency_panel(
    operons:          dict,
    families_summary: dict,
    family_colors:    FamilyColorMap,
    reference_family: int,
    most_common_fig,
    min_coocc:        float = 0.0,
) -> figure:
    """
    Build the gene adjacency network figure.

    Only edges between genes that appear consecutively in at least one
    genomic context are kept.
    """
    try:
        co_matrix, family_labels = _cooccurrence_matrix(operons, families_summary, min_coocc)
        adj_matrix, _            = _adjacency_matrix(operons, families_summary)
        # Mask co-occurrence edges that are not adjacent
        combined = np.where(adj_matrix > 0, co_matrix, 0)
        graph, coords, family_labels = _build_graph(
            combined, family_labels, family_colors, families_summary
        )
        return _network_figure(
            graph, coords, family_labels, families_summary,
            family_colors, reference_family, most_common_fig,
            title='Gene adjacency network',
        )
    except ValueError as exc:
        return empty_network_panel(
            most_common_fig,
            message=f'No adjacency detected ({exc})',
        )


def empty_network_panel(most_common_fig, message: str = 'No network data') -> figure:
    """A blank placeholder panel with a text message."""
    fig = figure(
        width=most_common_fig.width,
        height=most_common_fig.height,
        title=message,
        toolbar_location='left',
    )
    fig.grid.visible = False
    fig.outline_line_width = 0
    return fig


# ── matrix computation ────────────────────────────────────────────────────────

def _cooccurrence_matrix(
    operons:          dict,
    families_summary: dict,
    min_coocc:        float,
) -> tuple[np.ndarray, list[int]]:
    """
    Count how often each pair of non-zero families appears in the same
    genomic context.  Raises ValueError when the matrix is all-zero.
    """
    family_labels = sorted(f for f in families_summary if f > 0)
    n = len(family_labels)
    idx = {f: i for i, f in enumerate(family_labels)}
    matrix = np.zeros((n, n))

    for odata in operons.values():
        for structure in odata['operon_protein_families_structure']:
            unique = [f for f in set(structure) if f > 0]
            for a in range(len(unique)):
                for b in range(a + 1, len(unique)):
                    i, j = idx[unique[a]], idx[unique[b]]
                    matrix[i, j] += 1
                    matrix[j, i] += 1

    max_val = matrix.max()
    if max_val == 0:
        raise ValueError('co-occurrence matrix is all-zero')

    matrix /= max_val
    matrix = np.where(matrix < min_coocc, 0, matrix)
    matrix = np.where(matrix != 0, (np.exp(matrix * 2) - 1), matrix)

    # Drop rows/cols that are entirely zero
    keep = ~np.all(matrix == 0, axis=1)
    matrix = matrix[keep][:, keep]
    family_labels = [f for f, k in zip(family_labels, keep) if k]

    if len(family_labels) == 0:
        raise ValueError('all families below min_coocc threshold')

    return matrix, family_labels


def _adjacency_matrix(
    operons:          dict,
    families_summary: dict,
) -> tuple[np.ndarray, list[int]]:
    """Count how often families appear consecutively."""
    family_labels = sorted(f for f in families_summary if f > 0)
    n = len(family_labels)
    idx = {f: i for i, f in enumerate(family_labels)}
    matrix = np.zeros((n, n))

    for odata in operons.values():
        for structure in odata['operon_protein_families_structure']:
            for k in range(len(structure) - 1):
                a, b = structure[k], structure[k + 1]
                if a > 0 and b > 0:
                    i, j = idx[a], idx[b]
                    matrix[i, j] += 1
                    matrix[j, i] += 1

    return matrix, family_labels


# ── graph building ────────────────────────────────────────────────────────────

def _build_graph(
    matrix:           np.ndarray,
    family_labels:    list[int],
    family_colors:    FamilyColorMap,
    families_summary: dict,
) -> tuple[nx.Graph, dict, list[int]]:
    """
    Convert a co-occurrence/adjacency matrix into a NetworkX graph with visual
    attributes.

    Isolated nodes (degree == 0) are removed before building the graph so they
    never appear in the figure.  The possibly-shortened *family_labels* list is
    returned as the third element of the tuple.
    """
    # Drop rows/cols for families that have no connections in *this* matrix.
    # (The adjacency panel can produce all-zero rows after masking co-occurrence
    # with adjacency, even if _cooccurrence_matrix already filtered its own zeros.)
    degrees = matrix.sum(axis=1)
    keep_mask = degrees > 1
    if not keep_mask.all():
        matrix        = matrix[np.ix_(keep_mask, keep_mask)]
        family_labels = [f for f, k in zip(family_labels, keep_mask) if k]

    G = nx.from_numpy_array(matrix)
    if len(G.nodes) == 0:
        raise ValueError('graph has no nodes after removing isolated nodes')

    edge_cmap = plt.get_cmap('Greys')
    edge_norm = mcolors.Normalize(vmin=0, vmax=4)

    for u, v, attrs in G.edges(data=True):
        w = attrs['weight']
        rgba = edge_cmap(edge_norm(round(w)))
        G[u][v]['line_width'] = w
        G[u][v]['edge_color'] = RGB(
            int(255 * rgba[0]), int(255 * rgba[1]), int(255 * rgba[2])
        )
        G[u][v]['weight'] = 1  # uniform for layout

    for node in G.nodes:
        fam = family_labels[node]
        fc  = family_colors.get(fam, family_colors.get(0))
        G.nodes[node]['node_label'] = fam
        G.nodes[node]['node_color'] = fc.hex_color

    coords = nx.spring_layout(G)
    return G, coords, family_labels


# ── figure assembly ───────────────────────────────────────────────────────────

def _network_figure(
    graph:            nx.Graph,
    coords:           dict,
    family_labels:    list[int],
    families_summary: dict,
    family_colors:    FamilyColorMap,
    reference_family: int,
    most_common_fig,
    title:            str,
) -> figure:
    xs = [coords[n][0] for n in coords]
    ys = [coords[n][1] for n in coords]
    x_range = (min(xs) - 0.5, max(xs) + 0.5)
    y_range = (min(ys) - 0.5, max(ys) + 0.5)

    y_step = (y_range[1] - y_range[0]) * 0.09

    # Node feature table
    node_data = _node_data(
        graph, coords, family_labels, families_summary,
        family_colors, reference_family, y_step,
    )

    fig = figure(
        width=most_common_fig.width,
        height=most_common_fig.height,
        x_range=x_range,
        y_range=y_range,
        title=title,
        toolbar_location='left',
    )

    renderer = from_networkx(graph, coords, scale=1, center=(0, 0))
    renderer.edge_renderer.glyph = MultiLine(
        line_width='line_width', line_color='edge_color'
    )
    renderer.node_renderer.glyph = Scatter(
        size=22, marker='circle', fill_color='node_color'
    )
    fig.renderers.append(renderer)

    # ── inject tooltip fields into the node renderer's data source ────────
    # HoverTool/TapTool must target this data source; a separate ColumnDataSource
    # is invisible to them and causes "???" in every tooltip field.
    node_ds  = renderer.node_renderer.data_source
    node_ids = list(node_ds.data['index'])   # node indices in data-source order

    node_ds.data['protein_name'] = [
        families_summary.get(family_labels[n], {}).get('name', 'n.a.')
        for n in node_ids
    ]
    node_ds.data['family'] = [
        str(family_labels[n]) if family_labels[n] > 0 else ''
        for n in node_ids
    ]
    _tm_info = [
        _tm_from_function(families_summary.get(family_labels[n], {}), family_labels[n])
        for n in node_ids
    ]
    node_ds.data['tm_type']      = [t[0] for t in _tm_info]
    node_ds.data['tm_pred_text'] = [t[2] for t in _tm_info]
    node_ds.data['found_models'] = [
        _model_state(families_summary.get(family_labels[n], {}), family_labels[n])
        for n in node_ids
    ]
    node_ds.data['model_links'] = [
        _swiss_model_url(families_summary.get(family_labels[n], {}))
        for n in node_ids
    ]

    deg = {t[0]:t[1] for t in graph.degree}
    node_ds.data['degree'] = [deg[n] for n in node_ids]

    # ── text labels (family code + TM) drawn near each node ──────────────
    src = ColumnDataSource(node_data)
    fig.text('text_x', 'text_y', text='family',
             text_baseline='bottom', text_align='center',
             text_font_size={'value': '6pt'}, source=src)
    fig.text('tm_text_x', 'tm_text_y', text='tm_text',
             text_color='white', text_baseline='middle', text_align='center',
             text_font_size={'value': '6pt'}, source=src)

    # ── tools attached to the node renderer so they read the right source ─
    fig.add_tools(HoverTool(
        tooltips=[
            ('Protein name',     '@protein_name'),
            ('Family code',      '@family'),
            ('Degree',           '@degree'),
            ('Membrane protein', '@tm_pred_text @tm_type'),
            ('Structural model', '@found_models'),
        ],
        renderers=[renderer.node_renderer],
    ))
    fig.add_tools(TapTool(
        renderers=[renderer.node_renderer],
        callback=OpenURL(url='@model_links'),
    ))

    for ax in [fig.xaxis, fig.yaxis]:
        ax.major_tick_line_color = None
        ax.minor_tick_line_color = None
        ax.major_label_text_color = None
        ax.axis_line_width = 0
    fig.grid.visible = False
    fig.outline_line_width = 0

    return fig


def _node_data(
    graph:            nx.Graph,
    coords:           dict,
    family_labels:    list[int],
    families_summary: dict,
    family_colors:    FamilyColorMap,
    reference_family: int,
    y_step:           float,
) -> dict:
    data: dict[str, list] = {
        'text_x': [], 'text_y': [],
        'tm_text_x': [], 'tm_text_y': [],
        'family': [], 'tm_text': [], 'tm_type': [], 'tm_pred_text': [],
        'protein_name': [], 'found_models': [], 'model_links': [],
    }
    for node in graph.nodes:
        fam     = family_labels[node]
        c       = coords[node]
        summary = families_summary.get(fam, {})
        func    = summary.get('function', {})

        # TM annotation
        tm_type, tm_text, tm_mode = _tm_from_function(summary, fam)

        data['text_x'].append(c[0])
        data['text_y'].append(c[1] + y_step)
        data['tm_text_x'].append(c[0])
        data['tm_text_y'].append(c[1])
        data['family'].append(fam if fam > 0 and fam != reference_family else '')
        data['tm_text'].append(tm_text)
        data['tm_type'].append(tm_type)
        data['tm_pred_text'].append(tm_mode)
        data['protein_name'].append(summary.get('name', 'n.a.'))
        data['found_models'].append(_model_state(summary, fam))
        data['model_links'].append(_swiss_model_url(summary))

    return data


def _tm_from_function(summary: dict, fam: int) -> tuple[str, str, str]:
    """Return (tm_type, tm_text_on_arrow, tooltip_mode)."""
    if fam > 0 and 'function' in summary:
        fn = summary['function']
        if 'TM_topology' in fn:
            tm_type = fn['TM_topology']
            if tm_type:
                return tm_type, 'TM', 'Yes → type:'
            return '', '', 'No'
        return '', '', ''
    return 'n.a.', '', 'n.a.'
