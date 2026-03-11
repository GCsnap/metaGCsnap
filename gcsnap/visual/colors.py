"""
Color assignment for protein families.

Single entry point: :func:`assign_family_colors`.
"""
from __future__ import annotations
import random
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from gcsnap.visual.types import FamilyColor, FamilyColorMap


def _random_hex(seed: int) -> str:
    """Deterministic random hex colour seeded on the family id."""
    rng = random.Random(seed)
    r, g, b = rng.randint(50, 200), rng.randint(50, 200), rng.randint(50, 200)
    return f'#{r:02x}{g:02x}{b:02x}'


def assign_family_colors(
    families: list[int],
    reference_family: int,
    cmap_name: str = 'tab20',
) -> FamilyColorMap:
    """
    Assign a :class:`FamilyColor` to every family id.

    Rules
    -----
    - family == 0   → light grey fill, light grey border (unknown / non-conserved)
    - family == -1  → white fill, grey dashed border (pseudogene)
    - family == ref → grey fill, black border (the query protein)
    - everything else → colour from *cmap_name*, black border

    Parameters
    ----------
    families:
        All family ids present in the dataset (including 0 and -1).
    reference_family:
        The family id of the query / reference protein.
    cmap_name:
        Matplotlib colourmap name used for regular families.
    """
    cmap = plt.get_cmap(cmap_name)
    norm = mcolors.Normalize(vmin=0, vmax=max(1, len(families)))

    result: FamilyColorMap = {}
    for i, fam in enumerate(sorted(families)):
        if fam == 0:
            result[fam] = FamilyColor(
                hex_color  = mcolors.to_hex('lightgrey'),
                line_color = mcolors.to_hex('lightgrey'),
                line_style = '-',
            )
        elif fam == -1:
            result[fam] = FamilyColor(
                hex_color  = mcolors.to_hex('white'),
                line_color = mcolors.to_hex('grey'),
                line_style = ':',
            )
        elif fam == reference_family:
            result[fam] = FamilyColor(
                hex_color  = mcolors.to_hex('grey'),
                line_color = mcolors.to_hex('black'),
                line_style = '-',
            )
        else:
            result[fam] = FamilyColor(
                hex_color  = _random_hex(fam),
                line_color = mcolors.to_hex('black'),
                line_style = '-',
            )

    return result
