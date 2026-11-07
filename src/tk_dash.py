"""Dashed outlines on a Tk canvas that look the same on Windows and macOS.

Tk on macOS draws each dash with square ends that stick out half the line's width at both ends, so every dash is a line width longer
than asked and every gap a line width shorter: a 2 px wide outline with dash=(6, 3) had 1 px gaps, which hardly show, and one with
dash=(3, 2) or a 3 px one with (3, 2) none at all. The UMAP browser's section view showed long dashes, or one dash per side. dash()
gives the pattern to pass, so that the dashes and gaps come out as asked on both."""
import sys


def dash(pattern, width):
    """The dash option for an outline `width` pixels wide whose dashes and gaps should be `pattern` (pixels: dash, gap, dash, gap...)."""
    if sys.platform != "darwin" or width <= 1:
        return tuple(pattern)
    return tuple(max(1, round(p - width)) if i % 2 == 0 else round(p + width) for i, p in enumerate(pattern))
