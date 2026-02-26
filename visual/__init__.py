"""
visual – clean reimplementation of the metaGCsnap visualization layer.

Drop-in replacement for the ``visuals/`` package.  The old code is untouched.

Quick start
-----------
::

    from visual import Figures

    figs = Figures(config, gc, out_label)
    figs.run()

Or use the lower-level functions directly:

::

    from visual.colors import assign_family_colors
    from visual.static import draw_static
    from visual.interactive import draw_interactive
"""
from visual.runner import Figures
from visual.advanced import draw_advanced

__all__ = ['Figures', 'draw_advanced']
