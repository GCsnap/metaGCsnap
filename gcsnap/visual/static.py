"""
Static (matplotlib) genomic-context figure.

Saves a PNG/PDF with gene arrows coloured by protein family, one row per
operon type, plus a separate legend file.

Public function: :func:`draw_static`.
"""
from __future__ import annotations
import os

import numpy as np
import matplotlib.pyplot as plt

from gcsnap.visual.types import FamilyColorMap


def draw_static(
    syntenies:        dict,
    operons:          dict,
    families_summary: dict,
    family_colors:    FamilyColorMap,
    reference_family: int,
    out_dir:          str,
    out_format:       str = 'png',
) -> None:
    """
    Save the static genomic-context figure and its legend.

    Parameters
    ----------
    syntenies:
        ``gc.get_syntenies()``.
    operons:
        ``gc.get_selected_operons()``.
    families_summary:
        ``gc.get_families()``.
    family_colors:
        Output of :func:`visual.colors.assign_family_colors`.
    reference_family:
        Family id of the query protein.
    out_dir:
        Directory where the files are written.
    out_format:
        ``'png'`` or ``'pdf'``.
    """
    os.makedirs(out_dir, exist_ok=True)
    fig_path    = os.path.join(out_dir, f'genomic_context.{out_format}')
    legend_path = os.path.join(out_dir, f'genomic_context_legend.{out_format}')

    _draw_genomic_context(syntenies, operons, family_colors,
                          reference_family, fig_path, out_format)
    _draw_legend(family_colors, families_summary, reference_family,
                 legend_path, out_format)


# ── genomic context figure ────────────────────────────────────────────────────

def _draw_genomic_context(
    syntenies:        dict,
    operons:          dict,
    family_colors:    FamilyColorMap,
    reference_family: int,
    out_path:         str,
    out_format:       str,
) -> None:
    n_operons = len(operons)
    fig_height = max(2, int(n_operons / 1.5))
    fig, axes = plt.subplots(
        1, 2, figsize=(20, fig_height),
        gridspec_kw={'width_ratios': [4, 1]},
    )
    ax_gc, ax_hist = axes

    curr_y = n_operons
    all_xs:      list[float] = []
    all_pops:    list[float] = []
    yticklabels: list[str]   = []

    for operon_type in sorted(operons.keys()):
        odata   = operons[operon_type]
        target  = odata['target_members'][0]
        fg      = syntenies[target]['flanking_genes']
        meta    = syntenies[target]['assembly_metadata']
        assembly = meta.get('assembly_accession', '')
        taxon   = syntenies[target].get('taxonomy', {}).get('taxon_name', '')

        n_members = len(odata['target_members'])
        pop_pct   = n_members * 100.0 / len(syntenies)

        # Gene arrows
        for i, gene in enumerate(fg.get('cds_codes', [])):
            fam       = fg['families'][i]
            dx        = fg['relative_ends'][i] - fg['relative_starts'][i] + 1
            direction = fg['directions'][i]

            if direction == '-':
                x_tail = fg['relative_ends'][i]
                dx     = -dx
            else:
                x_tail = fg['relative_starts'][i]

            fc = family_colors.get(fam, family_colors.get(0))
            zorder = 1 if fam == 0 else len(fg['cds_codes']) - i + 1

            ax_gc.arrow(
                x_tail, curr_y, dx, 0,
                width=0.5, head_width=0.5,
                length_includes_head=True, head_length=100,
                zorder=zorder,
                facecolor=fc.hex_color,
                edgecolor=fc.line_color,
                linestyle=fc.line_style,
            )

            # Family-code label for non-reference families
            if fam < 0 and fam != reference_family:
                text_x = x_tail + dx / 2
                ax_gc.text(text_x, curr_y + 0.3, str(fam),
                           horizontalalignment='center', fontsize=7)

            all_xs += [fg['relative_starts'][i], fg['relative_ends'][i]]

        # Frequency bar
        ax_hist.arrow(
            0, curr_y, pop_pct, 0,
            width=0.5, head_width=0.5,
            length_includes_head=True, head_length=0,
            facecolor='black', edgecolor='black',
        )
        ax_hist.text(
            pop_pct + 2.5, curr_y,
            f'{pop_pct:.1f}%, ({n_members})',
            verticalalignment='center', fontsize=8,
        )

        all_pops.append(pop_pct)
        yticklabels.append(f'{operon_type} | {assembly}\n({taxon})')
        curr_y -= 1

    # Axis formatting – genomic context side
    if all_xs:
        ax_gc.set_xlim(min(all_xs) - 100, max(all_xs) + 100)
    ax_gc.set_ylim(0, n_operons + 1)
    ax_gc.set_yticks(np.arange(1, n_operons + 1))
    ax_gc.set_yticklabels(yticklabels[::-1], fontsize=10, ha='right')
    ax_gc.spines['right'].set_visible(False)
    ax_gc.spines['left'].set_visible(False)
    ax_gc.tick_params(axis='x', bottom=False, top=False, labelbottom=False)
    ax_gc.tick_params(axis='y', length=0, pad=4)   # no tick marks, labels close

    # Axis formatting – frequency histogram side
    ax_hist.set_xlim(0, max(all_pops, default=100) + 10)
    ax_hist.set_ylim(0, n_operons + 1)
    ax_hist.spines['right'].set_visible(False)
    ax_hist.tick_params(axis='y', left=False, right=False,
                        labelleft=False, labelright=False)
    ax_hist.set_xlabel('Operon frequency (%)')
    ax_hist.set_title(f'{sum(all_pops):.1f}% of all targets represented')

    try:
        plt.tight_layout()
    except Exception:
        pass

    plt.savefig(out_path, format=out_format, transparent=True)
    plt.close('all')
    display_path = os.path.join(os.path.basename(os.path.dirname(out_path)), os.path.basename(out_path))
    print(f' static figure saved → {display_path}')


# ── legend figure ─────────────────────────────────────────────────────────────

def _draw_legend(
    family_colors:    FamilyColorMap,
    families_summary: dict,
    reference_family: int,
    out_path:         str,
    out_format:       str,
) -> None:
    families  = sorted(family_colors.keys())
    n         = len(families)
    fig_height = max(2, int(n / 1.5))

    plt.clf()
    fig, ax = plt.subplots(figsize=(10, fig_height))

    curr_y = n
    for fam in families:
        fc = family_colors[fam]
        plt.arrow(
            0, curr_y, 5, 0,
            width=0.5, head_width=0.5,
            length_includes_head=True, head_length=0.5,
            facecolor=fc.hex_color,
            edgecolor=fc.line_color,
            linestyle=fc.line_style,
        )

        if fam == 0:
            label = 'Non-conserved gene'
        elif fam == reference_family:
            label = f'Target: {families_summary.get(fam, {}).get("name", str(fam))}'
        elif fam == -1:
            label = 'Pseudogene'
        else:
            plt.text(2.25, curr_y + 0.3, str(fam),
                     ha='center', fontsize=8)
            label = families_summary.get(fam, {}).get('name', str(fam))

        plt.text(7, curr_y, label, va='center', fontsize=9)
        curr_y -= 1

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.tick_params(axis='both', left=False, right=False,
                    bottom=False, top=False,
                    labelleft=False, labelright=False, labelbottom=False)
    plt.xlim(0, 50)
    plt.ylim(0, n + 1)

    plt.savefig(out_path, format=out_format, transparent=True)
    plt.close('all')
    display_path = os.path.join(os.path.basename(os.path.dirname(out_path)), os.path.basename(out_path))
    print(f'legend saved → {display_path}')
