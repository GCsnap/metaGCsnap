"""
Summary visualization – three cross-linked Bokeh scatter panels.

Panel 1  –  PaCMAP scatter of individual genomic contexts
Panel 2  –  Cluster similarity network (centroids + edge alpha ∝ similarity)
Panel 3  –  Sequence similarity scatter (CLANS map or MMseqs-derived PaCMAP)

The three panels share a single ColumnDataSource, so lasso/tap selections
propagate across all of them automatically.

Public entry point
------------------
draw_summary(gc, out_dir, clans_file=None, out_filename='resume.html') -> str
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pacmap

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import column, gridplot
from bokeh.models import (
    HoverTool, LassoSelectTool, ColumnDataSource,
    DataTable, TableColumn, HTMLTemplateFormatter,
)

from gcsnap.genomic_context import GenomicContext
from gcsnap.consts import MMSeqsParams


# ── public API ────────────────────────────────────────────────────────────────

def draw_summary(
    gc:           GenomicContext,
    out_dir:      str,
    clans_file:   str | None = None,
    out_filename: str = 'resume.html',
) -> str:
    """
    Build the three-panel summary HTML page and write it to *out_dir*.

    Parameters
    ----------
    gc:
        A :class:`~gcsnap.genomic_context.GenomicContext` whose syntenies
        already contain ``operon_filtered_PaCMAP`` (i.e.
        :class:`~gcsnap.annotations.operons.Operons` was run with
        ``operon_cluster_advanced=True``).
    out_dir:
        Directory where the HTML file will be saved.
    clans_file:
        Path to an optional ``.clans`` file.  When ``None`` the function
        looks for a pre-existing MMseqs TSV at
        ``{out_label}/genomic_context/sequences/flanking_sequences.mmseqs``
        and derives sequence-similarity coordinates from it via PaCMAP.
    out_filename:
        Basename of the output HTML file.

    Returns
    -------
    str
        Absolute path to the saved HTML file.

    Raises
    ------
    ValueError
        If ``operon_filtered_PaCMAP`` is not found in the syntenies.
    """
    syntenies = gc.get_syntenies()
    _check_pacmap(syntenies)

    os.makedirs(out_dir, exist_ok=True)

    # Build a fresh operons view directly from syntenies (avoids stale cache
    # when advanced operons were run after the initial Figures setup).
    operons        = _operons_from_syntenies(syntenies)
    cluster_colors = _cluster_colors(operons)
    clans_coords   = _get_clans_coords(gc, syntenies, clans_file)
    source         = _build_source(operons, syntenies, cluster_colors, clans_coords)

    p1 = _pacmap_scatter(source)
    p2 = _cluster_network(source, operons, p1)
    p3 = _clans_scatter(source)

    scatter_row = gridplot([[p1, p2, p3]], merge_tools=True)
    table       = _targets_table(operons, syntenies, cluster_colors)
    layout      = column(scatter_row, table)

    out_path = os.path.join(out_dir, out_filename)
    output_file(out_path)
    save(layout)
    return out_path


# ── data preparation ──────────────────────────────────────────────────────────

def _check_pacmap(syntenies: dict) -> None:
    first = next(iter(syntenies.values()))
    if 'operon_filtered_PaCMAP' not in first:
        raise ValueError(
            'operon_filtered_PaCMAP not found in syntenies. '
            'Run Operons with operon_cluster_advanced=True before calling draw_summary.'
        )


def _operons_from_syntenies(syntenies: dict) -> dict:
    """
    Build an operons dict with PaCMAP data directly from syntenies,
    bypassing any potentially stale gc.selected_operons cache.

    Keys are 'GC Type {:05d}' strings.  Centroids are computed on the fly
    for non-singleton clusters.
    """
    raw: dict = {}
    for target, syn in syntenies.items():
        key = 'GC Type {:05d}'.format(syn['operon_type'])
        if key not in raw:
            raw[key] = {'target_members': [], 'operon_filtered_PaCMAP': []}
        raw[key]['target_members'].append(target)
        raw[key]['operon_filtered_PaCMAP'].append(syn['operon_filtered_PaCMAP'])

    for key, odata in raw.items():
        if '-' not in key:
            odata['operon_centroid_PaCMAP'] = list(
                np.mean(odata['operon_filtered_PaCMAP'], axis=0)
            )

    return raw


def _cluster_colors(operons: dict) -> dict[str, str]:
    """
    Assign a hex colour to every operon type.
    Singletons (``'-'`` in key) → grey.  Others → gist_rainbow spread.
    """
    non_singletons = sorted(k for k in operons if '-' not in k)
    cmap = plt.get_cmap('gist_rainbow')
    norm = mcolors.Normalize(vmin=0, vmax=max(1, len(non_singletons) - 1))

    colors: dict[str, str] = {}
    for i, key in enumerate(non_singletons):
        rgba = cmap(norm(i))
        colors[key] = '#{:02x}{:02x}{:02x}'.format(
            int(255 * rgba[0]), int(255 * rgba[1]), int(255 * rgba[2])
        )
    for key in operons:
        if '-' in key:
            colors[key] = '#aaaaaa'

    return colors


def _build_source(
    operons:        dict,
    syntenies:      dict,
    cluster_colors: dict[str, str],
    clans_coords:   dict[str, tuple[float, float]],
) -> ColumnDataSource:
    """
    Build the shared ColumnDataSource used by all three panels.

    One row per target member; an extra centroid row per non-singleton
    cluster (used by the network panel nodes, invisible in panels 1 and 3).
    """
    data: dict[str, list] = {
        'x':         [],   # PaCMAP x per member
        'y':         [],   # PaCMAP y per member
        'avg_x':     [],   # centroid x (NaN for member rows)
        'avg_y':     [],   # centroid y (NaN for member rows)
        'clans_x':   [],   # sequence-similarity x per member
        'clans_y':   [],   # sequence-similarity y per member
        'facecolor': [],
        'edgecolor': [],
        'size':      [],   # dot radius for member rows
        'node_size': [],   # dot radius for centroid rows (0 for member rows)
        'type':      [],
        'target':    [],
        'species':   [],
    }

    for operon_type, odata in sorted(operons.items()):
        color      = cluster_colors[operon_type]
        is_single  = '-' in operon_type
        dot_size   = 2 if is_single else 5
        edge_color = color if is_single else '#000000'

        for i, target in enumerate(odata['target_members']):
            x, y    = odata['operon_filtered_PaCMAP'][i]
            cx, cy  = clans_coords.get(target, (np.nan, np.nan))
            species = (syntenies[target]
                       .get('taxonomy', {})
                       .get('taxon_name', 'n.a.'))

            data['x'].append(float(x))
            data['y'].append(float(y))
            data['avg_x'].append(np.nan)
            data['avg_y'].append(np.nan)
            data['clans_x'].append(float(cx) if not np.isnan(cx) else np.nan)
            data['clans_y'].append(float(cy) if not np.isnan(cy) else np.nan)
            data['facecolor'].append(color)
            data['edgecolor'].append(edge_color)
            data['size'].append(dot_size)
            data['node_size'].append(0)
            data['type'].append(operon_type)
            data['target'].append(str(target))
            data['species'].append(species)

        # Centroid entry — rendered only by the network panel (avg_x/avg_y columns)
        if not is_single and 'operon_centroid_PaCMAP' in odata:
            cx_c, cy_c = odata['operon_centroid_PaCMAP']
            data['x'].append(np.nan)
            data['y'].append(np.nan)
            data['avg_x'].append(float(cx_c))
            data['avg_y'].append(float(cy_c))
            data['clans_x'].append(np.nan)
            data['clans_y'].append(np.nan)
            data['facecolor'].append(color)
            data['edgecolor'].append('#000000')
            data['size'].append(0)
            data['node_size'].append(10)
            data['type'].append(operon_type)
            data['target'].append('')
            data['species'].append('')

    return ColumnDataSource(data)


# ── CLANS / sequence-similarity coordinates ───────────────────────────────────

def _get_clans_coords(
    gc:         GenomicContext,
    syntenies:  dict,
    clans_file: str | None,
) -> dict[str, tuple[float, float]]:
    if clans_file is not None:
        return _parse_clans_file(clans_file)

    mmseqs_path = _find_mmseqs_file(gc)
    if mmseqs_path is not None:
        return _coords_from_mmseqs(mmseqs_path, syntenies)

    # No sequence similarity data available — Panel 3 will be empty
    return {}


def _find_mmseqs_file(gc: GenomicContext) -> str | None:
    candidate = os.path.join(
        gc.out_label, 'genomic_context', 'sequences', 'flanking_sequences.mmseqs'
    )
    return candidate if os.path.exists(candidate) else None


def _parse_clans_file(path: str) -> dict[str, tuple[float, float]]:
    """Read sequence accessions and 2-D positions from a .clans file."""
    seq_map:   dict[int, str] = {}
    coords:    dict[str, tuple[float, float]] = {}
    seq_count  = 0
    in_seq     = False
    in_pos     = False

    with open(path, 'r') as fh:
        for line in fh:
            if '<seq>' in line:
                in_seq = True
            elif '</seq>' in line:
                in_seq = False
            elif in_seq and line.startswith('>'):
                code = (line[1:].split()[0]
                        .split(':')[0].split('|')[0].split('_#')[0].strip())
                seq_map[seq_count] = code
                seq_count += 1
            elif '<pos>' in line:
                in_pos = True
            elif '</pos>' in line:
                in_pos = False
            elif in_pos and not in_seq:
                parts = line.strip().split()
                idx, x, y = int(parts[0]), float(parts[1]), float(parts[2])
                coords[seq_map[idx]] = (x, -y)   # flip y as in original gcsnap

    return coords


def _coords_from_mmseqs(
    mmseqs_path: str,
    syntenies:   dict,
) -> dict[str, tuple[float, float]]:
    """
    Derive 2-D sequence-similarity coordinates from an existing
    flanking_sequences.mmseqs all-vs-all TSV.

    Only target↔target hits are used; the resulting identity matrix is
    embedded with PaCMAP to produce coordinates comparable to a CLANS map.
    """
    targets        = set(syntenies.keys())
    sorted_targets = sorted(targets)
    idx            = {t: i for i, t in enumerate(sorted_targets)}
    n              = len(sorted_targets)

    df = pd.read_csv(
        mmseqs_path, sep='\t', header=None,
        names=MMSeqsParams.tsv_columns,
        usecols=['query', 'target', 'pident'],
    )
    df = df[df['query'].isin(targets) & df['target'].isin(targets)]

    sim = np.zeros((n, n))
    np.fill_diagonal(sim, 1.0)
    for _, row in df.iterrows():
        i = idx.get(row['query'])
        j = idx.get(row['target'])
        if i is not None and j is not None:
            val = float(row['pident']) / 100.0   # pident is 0–100, normalise to 0–1
            sim[i, j] = max(sim[i, j], val)
            sim[j, i] = max(sim[j, i], val)

    dist       = 1.0 - sim
    embedding  = pacmap.PaCMAP(n_components=2)
    coords_arr = embedding.fit_transform(dist)

    return {
        sorted_targets[i]: (float(coords_arr[i, 0]), float(coords_arr[i, 1]))
        for i in range(n)
    }


# ── edge data for cluster network ─────────────────────────────────────────────

def _edge_source(operons: dict) -> ColumnDataSource:
    """
    Compute pairwise minimum distances between cluster members in PaCMAP space,
    normalise to [0, 1] with a power transform, and return a ColumnDataSource
    of line segments (one per cluster pair) for the network panel.
    """
    non_singletons = sorted(k for k in operons if '-' not in k
                            and 'operon_centroid_PaCMAP' in operons[k])

    xs:     list = []
    ys:     list = []
    alphas: list = []

    for i, oi in enumerate(non_singletons):
        for j, oj in enumerate(non_singletons):
            if i >= j:
                continue
            coords_i = np.array(operons[oi]['operon_filtered_PaCMAP'])
            coords_j = np.array(operons[oj]['operon_filtered_PaCMAP'])
            min_dist = min(
                np.linalg.norm(pi - pj)
                for pi in coords_i for pj in coords_j
            )
            xi, yi = operons[oi]['operon_centroid_PaCMAP']
            xj, yj = operons[oj]['operon_centroid_PaCMAP']
            xs.append([float(xi), float(xj)])
            ys.append([float(yi), float(yj)])
            alphas.append(min_dist)

    # Normalise: nearby clusters → opaque edges, distant → transparent
    if alphas:
        min_d = min(alphas)
        max_d = max(alphas)
        span  = max_d - min_d if max_d > min_d else 1.0
        alphas = [
            round((1.0 - (d - min_d) / span) ** 30, 3)
            for d in alphas
        ]

    return ColumnDataSource({'x': xs, 'y': ys, 'alpha': alphas})


# ── panel builders ────────────────────────────────────────────────────────────

_TOOLTIPS = [
    ('GC type', '@type'),
    ('Species', '@species'),
]


def _pacmap_scatter(source: ColumnDataSource) -> figure:
    p = figure(
        title='Genomic context types/clusters',
        width=500, height=500,
    )
    p.circle(
        'x', 'y', size='size',
        line_color='edgecolor', fill_color='facecolor',
        alpha=1, source=source,
    )
    _clean_axes(p)
    p.add_tools(HoverTool(tooltips=_TOOLTIPS))
    p.add_tools(LassoSelectTool())
    p.background_fill_color = 'lightgrey'
    p.background_fill_alpha = 0.2
    return p


def _cluster_network(
    source:  ColumnDataSource,
    operons: dict,
    ref_fig: figure,
) -> figure:
    non_singletons = [k for k in operons
                      if '-' not in k and 'operon_centroid_PaCMAP' in operons[k]]

    if len(non_singletons) < 2:
        return _empty_panel(
            ref_fig,
            'Genomic context types/clusters similarity network'
            ' (need ≥ 2 clusters)',
        )

    edge_src = _edge_source(operons)

    p = figure(
        title='Genomic context types/clusters similarity network',
        width=ref_fig.width, height=ref_fig.height,
        x_range=ref_fig.x_range, y_range=ref_fig.y_range,
    )
    p.multi_line('x', 'y', color='black', alpha='alpha', source=edge_src)
    p.circle(
        'avg_x', 'avg_y', size='node_size',
        line_color='edgecolor', fill_color='facecolor',
        alpha=1, source=source, name='nodes',
    )
    _clean_axes(p)
    p.add_tools(HoverTool(tooltips=_TOOLTIPS, name='nodes'))
    p.add_tools(LassoSelectTool())
    p.background_fill_color = 'lightgrey'
    p.background_fill_alpha = 0.2
    return p


def _clans_scatter(source: ColumnDataSource) -> figure:
    p = figure(
        title='Sequence similarity cluster (CLANS) map',
        width=500, height=500,
    )
    p.circle(
        'clans_x', 'clans_y', size='size',
        line_color='edgecolor', fill_color='facecolor',
        alpha=1, source=source,
    )
    _clean_axes(p)
    p.add_tools(HoverTool(tooltips=_TOOLTIPS))
    p.add_tools(LassoSelectTool())
    p.background_fill_color = 'lightgrey'
    p.background_fill_alpha = 0.2
    return p


# ── table ─────────────────────────────────────────────────────────────────────

def _targets_table(
    operons:        dict,
    syntenies:      dict,
    cluster_colors: dict[str, str],
) -> DataTable:
    """
    Scrollable DataTable listing every target with its metadata.

    Columns: Target, Source, Target source, Taxon name, Taxon ID,
             GC type, Color (coloured square), Bin type.
    """
    data: dict[str, list] = {
        'target':        [],
        'source':        [],
        'target_source': [],
        'taxon_name':    [],
        'taxon_id':      [],
        'gc_type':       [],
        'color':         [],
        'bin_type':      [],
    }

    for operon_type, odata in sorted(operons.items()):
        color = cluster_colors[operon_type]
        for target in odata['target_members']:
            syn = syntenies[target]
            am  = syn.get('assembly_metadata', {})
            tax = syn.get('taxonomy', {})

            data['target'].append(am.get('target', str(target)))
            data['source'].append(am.get('source', 'n.a.'))
            data['target_source'].append(am.get('target_source', 'n.a.'))
            data['taxon_name'].append(tax.get('taxon_name', 'n.a.'))
            data['taxon_id'].append(str(tax.get('taxon_id', 'n.a.')))
            data['gc_type'].append(operon_type)
            data['color'].append(color)
            data['bin_type'].append(str(syn.get('bin_type', 'n.a.')))

    color_formatter = HTMLTemplateFormatter(
        template=(
            '<span style="color:<%= value %>; font-size:18pt; '
            'text-shadow: 1px 1px 2px #000000;">&#9632;</span>'
        )
    )

    columns = [
        TableColumn(field='target',        title='Target'),
        TableColumn(field='source',        title='Source'),
        TableColumn(field='target_source', title='Target source'),
        TableColumn(field='taxon_name',    title='Taxon name'),
        TableColumn(field='taxon_id',      title='Taxon ID'),
        TableColumn(field='gc_type',       title='GC type'),
        TableColumn(field='color',         title='Color', formatter=color_formatter),
        TableColumn(field='bin_type',      title='Bin type'),
    ]

    return DataTable(
        source      = ColumnDataSource(data),
        columns     = columns,
        width       = 1500,
        height      = 300,
        index_position = 0,
        sortable    = True,
        reorderable = True,
    )


# ── helpers ───────────────────────────────────────────────────────────────────

def _clean_axes(p: figure) -> None:
    for ax in [p.xaxis, p.yaxis]:
        ax.major_tick_line_color  = None
        ax.minor_tick_line_color  = None
        ax.major_label_text_color = None
        ax.axis_line_width        = 0
    p.grid.visible = False


def _empty_panel(ref_fig: figure, message: str) -> figure:
    p = figure(
        width=ref_fig.width, height=ref_fig.height,
        title=message,
    )
    _clean_axes(p)
    p.outline_line_width = 0
    return p
