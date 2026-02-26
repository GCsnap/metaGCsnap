"""
Top-level entry point for visual.

:class:`Figures` is a drop-in replacement for ``visuals.figures.Figures``.
It has the same ``__init__`` signature and a ``run()`` method, but internally
it delegates to the clean, composable functions in the rest of the package.
"""
from __future__ import annotations
import os

from gcsnap.configuration import Configuration
from gcsnap.genomic_context import GenomicContext
from gcsnap.rich_console import RichConsole

from visual.colors import assign_family_colors
from visual.static import draw_static
from visual.interactive import draw_interactive
from visual.advanced import draw_advanced


class Figures:
    """
    Orchestrates static + interactive figure generation.

    Parameters
    ----------
    config:
        The application :class:`~gcsnap.configuration.Configuration`.
    gc:
        The :class:`~gcsnap.genomic_context.GenomicContext` with all data.
    out_label:
        Base output path (same as the rest of the pipeline).
    starting_directory:
        Unused – kept for API compatibility with the old ``visuals.Figures``.
    """

    def __init__(
        self,
        config: Configuration,
        gc:     GenomicContext,
        out_label: str,
        starting_directory: str | None = None,
    ) -> None:
        self.config    = config
        self.gc        = gc
        self.out_label = out_label

        # Derived from gc once, reused throughout
        self.operons          = gc.get_selected_operons()
        self.syntenies        = gc.get_syntenies()
        self.families_summary = gc.get_families()
        self.most_populated_operon = gc.get_most_populated_operon()

        # The reference family is the family of the first member of the most
        # populated operon – same convention as the old code.
        self.reference_family: int = (
            self.syntenies
            [self.operons[self.most_populated_operon]['target_members'][0]]
            ['target_family']
        )

        self.cmap       = config.arguments['genomic_context_cmap']['value']
        self.out_format = config.arguments['out_format']['value']
        self.out_dir    = os.path.join(out_label, 'visuals')
        self.console    = RichConsole()

        # Family colours computed once; all draw_* functions share this object.
        self.family_colors = assign_family_colors(
            families         = list(self.families_summary.keys()),
            reference_family = self.reference_family,
            cmap_name        = self.cmap,
        )

    # ── public API ────────────────────────────────────────────────────────

    def run(self) -> None:
        """Run both the static and (optionally) the interactive figure."""
        self.run_static()
        if self.config.arguments['interactive']['value']:
            self.run_interactive()

    def run_static(self) -> None:
        """Run only the static (matplotlib) figure."""
        try:
            with self.console.status('Creating static genomic context figure'):
                draw_static(
                    syntenies        = self.syntenies,
                    operons          = self.operons,
                    families_summary = self.families_summary,
                    family_colors    = self.family_colors,
                    reference_family = self.reference_family,
                    out_dir          = self.out_dir,
                    out_format       = self.out_format,
                )
            self.console.print_info(
                f'Static genomic context figures created'
            )
        except Exception as exc:
            self.console.print_warning(f'Static figure failed: {exc}')

    def run_interactive(self, sort_mode: str | None = None) -> str:
        """
        Run only the interactive (Bokeh) figure and return the HTML path.

        Parameters
        ----------
        sort_mode:
            Override the sort mode from config.  Pass ``None`` to use the
            value from ``config.arguments['sort_mode']['value']``.
        """
        if sort_mode is None:
            sort_mode = self.config.arguments['sort_mode']['value']

        try:
            with self.console.status(
                f'Creating interactive figure (sort_mode={sort_mode!r})'
            ):
                path = draw_interactive(
                    gc               = self.gc,
                    family_colors    = self.family_colors,
                    reference_family = self.reference_family,
                    sort_mode        = sort_mode,
                    out_dir          = self.out_dir,
                )
            return path
        except Exception as exc:
            self.console.print_warning(f'Interactive figure failed: {exc}')
            return ''

    def run_advanced(
        self,
        sort_mode:      str | None = None,
        gc_legend_mode: str = 'assembly',
        group_by:       str = 'operon',
        select_groups:  list | None = None,
        out_filename:   str | None = None,
    ) -> str:
        """
        Run the advanced interactive figure and return the HTML path.

        Produces a single multi-tab HTML file with one Bokeh tab per group.

        Parameters
        ----------
        sort_mode:
            Override the sort mode from config.  ``None`` uses the config value.
        gc_legend_mode:
            Y-axis label mode for the genomic map rows.
            ``'assembly'`` (default), ``'operon'``, or a taxonomy rank string.
        group_by:
            How to partition targets into tabs.

            * ``'operon'`` (default) – one tab per GC type.
            * ``'metagenomic_bins'`` – one tab per SourMash DNA bin.
        select_groups:
            Optional subset of groups to include.  Pass integers (e.g.
            ``[1, 2, 5]``) or full key strings.  ``None`` includes all.
        out_filename:
            Output HTML basename.  ``None`` → ``advanced_{group_by}.html``.
        """
        if sort_mode is None:
            sort_mode = self.config.arguments['sort_mode']['value']

        try:
            with self.console.status(
                f'Creating advanced interactive figure '
                f'(group_by={group_by!r}, sort_mode={sort_mode!r})'
            ):
                path = draw_advanced(
                    gc               = self.gc,
                    family_colors    = self.family_colors,
                    reference_family = self.reference_family,
                    sort_mode        = sort_mode,
                    out_dir          = self.out_dir,
                    gc_legend_mode   = gc_legend_mode,
                    group_by         = group_by,
                    select_groups    = select_groups,
                    out_filename     = out_filename,
                )
            display_path = os.path.join(os.path.basename(os.path.dirname(path)), os.path.basename(path))
            self.console.print_info(f'Advanced interactive figure saved to {display_path}')
            return path
        except Exception as exc:
            self.console.print_warning(f'Advanced interactive figure failed: {exc}')
            return ''
