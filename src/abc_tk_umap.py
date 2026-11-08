"""A Tk version of the atlas browser's interactive UMAP viewer (first version: the UMAP, a grid of the run's brain sections, colouring
by group / single group(s) / gene / imputed gene, and the hover link between the two).

Why it is quicker than the matplotlib window: every plot here is one image on a Tk canvas that is drawn directly from the cells' positions
(a point is a few pixels written into an array), so a recolour, a zoom or a pan costs the same however many artists the figure would
have had, and each section panel is a plain canvas of its own with its own zoom and pan. The hover ring is a canvas item.

    show_umap_window_tk(adata, abc_cache, imputed_state=None, adata_backed=None, rois=None, run_folder=None, ensure_imputed=None)

blocks until the window is closed, like the matplotlib one. Not yet in this version: more than one gene at once, the UMAP variants,
the group (family) highlight on hover, and the exports.
"""
import colorsys
import gc
import queue
import math
import re
import sys
from pathlib import Path
import threading
import time
import types
import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk, messagebox

import cv2
import numpy as np
import pandas as pd
from PIL import Image, ImageTk
from scipy.spatial import cKDTree

import abc_core as core
import abc_project
from tk_dash import dash
import abc_tk_ccf
from abc_tk_umap_export import ExportsMixin, hexcolor

LEIDEN_KEY = "leiden"
GREY = np.array([205, 205, 205], np.uint8)
MULTI_BG = (26, 26, 26)           # the plots' background in the gene views (a zero-expression cell is black in the multi-gene view, so it must not melt into it)
CCF_FIT_CELLS_PER_SECTION = 30000   # the CCF atlas fit of a section map uses at most this many of its cells (the lookup of millions is seconds)
MAX_GENES = 6                      # the first three are red / green / blue circles; the next three the same colors as small dots on top
CHANNEL_HEX = ("#e53935", "#43a047", "#4a6cff", "#e53935", "#43a047", "#4a6cff")
CHANNEL_NAMES = ("red", "green", "blue", "red", "green", "blue")
CHANNEL_HEX_LIGHT = ("#00a8cc", "#d81b9a", "#d99a00", "#00a8cc", "#d81b9a", "#d99a00")        # the light theme's inks: cyan, magenta, yellow (a deep one: pale on white)
CHANNEL_NAMES_LIGHT = ("cyan", "magenta", "yellow", "cyan", "magenta", "yellow")
INK_ABSORB = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.3, 1.0]])    # what each ink takes out of (red, green, blue): cyan the red, magenta the green, yellow the blue
INK_PAPER = 230.0                 # the light theme's colour for no expression: a light grey, close to the white background
LIGHT_GREY = np.array([222, 222, 222], np.uint8)
RING = "#ff1493"
MAX_COLORED_CATEGORIES = 400
LEVEL_NAMES = {"class": "Class", "subclass": "Subclass", "supertype": "Supertype", "cluster": "Cluster", "leiden": "Leiden", "neurotransmitter": "Neurotransmitter"}
NOT_IN_TAXONOMY = ("leiden", "neurotransmitter")        # the levels that are not a step down from the one above (so a group's table is just its own row)
TILE_PX = 230
BLINK_PAUSE_MS = 1000              # a highlighted cell (its ring or white dot) starts to blink when it has stood this long...
BLINK_ON_MS, BLINK_OFF_MS = 450, 250       # ...shown this long, hidden this long, over and over (until another cell is highlighted or it is cleared)
MIN_TILE_PX = 50
ROI_LINE_PX = {8: 2, 7: 2, 6: 2, 5: 3, 4: 4, 3: 5, 2: 6, 1: 8}      # (more across than that, as a fit can have: 2)        # the ROI outline's line width at each Section size stop (panels across): thicker as the sections grow
TILE_COLUMNS = (0, 8, 7, 6, 5, 4, 3, 2, 1)    # the Section size slider's stops: how many panels fit across the pane (1/8 of its width ... all of it); 0: the
                                              # size that shows every section in the pane at once, as large as that allows (see fit_tiles)
FIT_STOP = TILE_COLUMNS.index(0)
LARGE_PROJECT_CELLS = 200_000             # a project of this many cells or more starts with its UMAP hidden (Show UMAP)
OUTLINE_PX_LARGEST = 2                    # the CCF outlines' width (pixels) at the largest section map size (one pixel at the others)
DEFAULT_TILE_INDEX = TILE_COLUMNS.index(2)    # (two across)
CELL_SIZES_UM = (20, 30, 40, 50, 60, 70, 80, 90, 100)       # the Cell size slider's stops: the width of a cell's dot in the section panels, in microns
OTHER_CELLS_DOT_RATIO = 0.5       # the dots of the cells of the other kind (e.g. the non-neurons of a neuron run), relative to the others'
DEFAULT_DOT_INDEX = 2                     # (40 microns, DOT_MM: the size the dots have always had)
MARKER_SIZES = (0.5, 1.0, 2.0)            # the Marker size slider's stops: times the size of the hover ring and of the red + marks in the section panels
DEFAULT_MARKER_INDEX = 1
GROUP_SORT_DESCENDING_FIRST = {"count"}     # the group list's columns whose first click sorts largest first (the others: smallest / A first)
LABEL_SIZES = (0.7, 1.0, 1.4, 2.0)         # the UMAP's Label size slider's stops: times the group numbers' font size
LABEL_SIZE_NAMES = ("Small", "Medium", "Large", "Very large")
DEFAULT_LABEL_INDEX = 1
TRANSLUCENT = 0.3                                  # how opaque one cell's dot is with the Translucent box ticked (overlapping dots add up)
OPACITIES = (0.15, 0.25, 0.35, 0.5, 0.75, 1.0)      # the Opacity slider's stops, low to high: how opaque one cell's dot is (overlapping dots add up)
DEFAULT_OPACITY_INDEX = 5                          # (fully opaque: the dots on top hide the ones below)
LABEL_ZOOM_CAP = 5.0               # the group numbers grow with the zoom up to this, then stay the same size
LAYOUT_SETTLE_MS = 150             # dragging the divider (or the window's edge) changes the section pane's width many times a second: the panels are laid out and drawn again this long after it stops
HOVER_LABEL_PX = 170               # the cell tip's table (CellTip): its first column wraps here
HOVER_GAP_PX = 8                   # the space between the cell tip's columns
PREVIEW_ZOOM = True                # a wheel step on a plot shows the last picture scaled, and draws the view for real once the wheel is still
ZOOM_DRAW_MS = 120
PAN_DRAW_MS = 45                   # while a plot is dragged it is drawn again this often (in between, what is on the canvas is moved as it is)
ZOOM_FILL_PARALLEL = 4             # while a zoom settles, this many section maps are drawn at once (in threads of their own: see PointCanvas.draw)
PAN_FILL_GAP_MS = 10               # while the zoomed panel is dragged, its maps are drawn one at a time, this far apart (for the mouse's moves between)
PANEL_PAN_SETTLE_MS = 150          # ...and all of them once the drag has been still this long
SHAPE_SETTLE_MS = 250              # while the view is being zoomed or panned the cells are drawn as plain discs; their shapes come back this long after it stops
LABEL_WEIGHT_RANGE = (0.6, 2.0)    # the size factor of the smallest and the largest group's number, when sized by cell count
LABEL_HIGHLIGHT_FACTOR = 2.0       # the number of the group being marked with '+' signs is drawn this much bigger, for as long as it is marked
HOVER_HOLD_MS = 250                # the pointer must rest this long before the plots dim
GROUP_HOLD_MS = 500                # ...and this long before the cell's group is marked
CELL_TIP_MS = 1000                 # ...and this long before its description comes up in a tip too (CellTip), with a Copy button
CELL_TIP_GRACE_PX = 40             # the tip stays while the pointer is this close to where it came (on its way onto the tip)
GROUP_TIP_MS = 400                 # a row of the group list (or of the Single group(s) legend) is described in a tip after this long
DENSE_UMAP_CELLS = 100_000         # a UMAP of more cells than this starts with Small dots (its Cell size slider's default)
DOT_MM = 0.04                      # a cell's dot in the section panels: 40 um across on the tissue
BIG_DOT_PX = 40                    # dots wider than this are drawn one by one as circles (the maximum filter gets slow)
SCALEBAR_UM = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
HOVER_RADIUS_PX = 9
HOVER_REACH_DOTS = 0.75               # at a high zoom the pointer reaches a cell this many dot widths away (HOVER_RADIUS_PX is the least)
SCALE_UNITS = {("standard", "log"): "log2(1 + expr)", ("standard", "linear"): "expr (raw counts)",
               ("imputed", "log"): "log2 (imputed)", ("imputed", "linear"): "linear (imputed)"}


def rank_color(rank):
    """The browser's colour for a category of this rank (0 = commonest): ten hues stepped by the golden angle, in a vivid and a muted set."""
    hue = ((rank % 10) * 0.6180339887498949) % 1.0
    s, v = (0.90, 0.95) if (rank // 10) % 2 == 0 else (0.45, 0.85)
    return np.array([int(round(c * 255)) for c in colorsys.hsv_to_rgb(hue, s, v)], np.uint8)


def plural_name(name):
    """'subclass' -> 'subclasses', 'cluster' -> 'clusters' (a level's name as the status line says it)."""
    low = name if any(c.isupper() for c in name[1:]) else name.lower()            # ('MMC Class' keeps its capitals)
    if any(c.isdigit() for c in low) or len(low) > 28:                              # (a long column name, 'DF classifications 0.25 0.15 3107')
        return low + " groups"
    if low.endswith("ss"):
        return low + "es"                                                           # (classes)
    return low if low.endswith("s") else low + "s"                                  # (clusters stay; subclass -> subclasses is the 'ss' above)


def is_dark(bgcolor):
    return sum(bgcolor) < 300


def ink_for(color, bgcolor):
    """The colour of a group's number on this background. The numbers are made in a darker shade of the group's colour (see color_all), which is
    for a light background; on a dark one the shade is brightened back."""
    if not is_dark(bgcolor):
        return color
    rgb = [int(color[i:i + 2], 16) for i in (1, 3, 5)]
    return "#%02x%02x%02x" % tuple(min(int(round(v / 0.55)), 255) for v in rgb)


def ink_rgb(norms):
    """The colours of cells shown in the light theme with several genes: each gene is an ink (cyan, magenta, yellow, in the order of the genes)
    that takes its colour out of a light grey paper in proportion to the gene's value (`norms`: one array in 0..1 for each of up to three genes),
    so that two genes together give blue, red or green, and three give a dark grey, and no expression is the paper's grey."""
    absorbed = sum(np.asarray(t, dtype=float)[:, None] * INK_ABSORB[i] for i, t in enumerate(norms))
    return np.rint(INK_PAPER * (1.0 - np.clip(absorbed, 0.0, 1.0))).astype(np.uint8)


LIGHT_STOPS = ((0.0, (230, 230, 230)), (0.25, (253, 212, 158)), (0.5, (252, 141, 89)), (0.75, (215, 48, 31)), (1.0, (127, 0, 0)))


def light_rgb(t):
    """The colours of the light theme's gene scale, for values t in 0..1: a light grey for none (close to a white background), then through
    pale orange and orange to red and a dark red, so that the highest values stand out against the low ones."""
    t = np.clip(np.asarray(t, dtype=float), 0.0, 1.0)
    xs = [s[0] for s in LIGHT_STOPS]
    out = np.stack([np.interp(t, xs, [s[1][c] for s in LIGHT_STOPS]) for c in range(3)], axis=-1)
    return np.rint(out).astype(np.uint8)


def viridis_rgb(t):
    from matplotlib import colormaps
    return (colormaps["viridis"](np.clip(t, 0, 1))[:, :3] * 255).astype(np.uint8)


# ---------------------------------------------------------------------------------------------------------------------- the point canvas

# The browser's category markers (CATEGORY_MARKER_SHAPES): a group's rank picks a color from a cycle of 20, and every 20 ranks the shape changes (and
# after 200, the shapes are hollow), so that groups that share a color still differ. A marker code is its shape's place in SHAPE_MARKERS, plus 10 if hollow.
SHAPE_MARKERS = ('o', 's', '^', '<', 'v', '>', 'D', 'P', 'X', '*')


def _plus_points(arm):
    a, b = (1 - arm) / 2, (1 + arm) / 2
    return [(a, 0), (b, 0), (b, a), (1, a), (1, b), (b, b), (b, 1), (a, 1), (a, b), (0, b), (0, a), (a, a)]


def _x_points(arm):
    pts = np.array(_plus_points(arm)) - 0.5
    rot = pts @ np.array([[np.cos(np.pi / 4), -np.sin(np.pi / 4)], [np.sin(np.pi / 4), np.cos(np.pi / 4)]]).T
    return (rot * (0.5 / np.abs(rot).max()) + 0.5).tolist()


def _star_points():
    return [(0.5 + (0.5 if k % 2 == 0 else 0.2) * np.cos(-np.pi / 2 + k * np.pi / 5), 0.5 + (0.5 if k % 2 == 0 else 0.2) * np.sin(-np.pi / 2 + k * np.pi / 5))
            for k in range(10)]


SHAPE_POLYGONS = {      # each marker but the circle, as a polygon in a unit square (y down)
    's': [(0, 0), (1, 0), (1, 1), (0, 1)], '^': [(0.5, 0), (1, 1), (0, 1)], 'v': [(0, 0), (1, 0), (0.5, 1)],
    '<': [(1, 0), (1, 1), (0, 0.5)], '>': [(0, 0), (1, 0.5), (0, 1)], 'D': [(0.5, 0), (1, 0.5), (0.5, 1), (0, 0.5)],
    'P': _plus_points(0.34), 'X': _x_points(0.26), '*': _star_points()}


def shape_mask(code, d):
    """The marker of this code as a d x d mask (1 inside), the right way up. (A maximum filter writes its kernel turned half way round, so the
    kernels made from these are turned the other way first: see PointCanvas.marker_kernel.)"""
    marker, is_open = SHAPE_MARKERS[code % 10], code >= 10
    S = 8 * d                                              # (drawn eight times as big and shrunk, so that small shapes keep their form)
    big = np.zeros((S, S), np.uint8)
    if marker == 'o':
        cv2.circle(big, (S // 2, S // 2), S // 2, 255, -1)
    else:
        cv2.fillPoly(big, [np.round(np.array(SHAPE_POLYGONS[marker]) * (S - 1)).astype(np.int32)], 255)
    mask = (cv2.resize(big, (d, d), interpolation=cv2.INTER_AREA) >= 128).astype(np.uint8)
    if is_open and d > 2:
        t = max(1, int(round(d / 6)))
        mask = mask & (1 - cv2.erode(mask, np.ones((2 * t + 1, 2 * t + 1), np.uint8), borderType=cv2.BORDER_CONSTANT, borderValue=0))   # (the box's edge is not inside)
    return mask


def draw_big_marker(img, X, Y, c, d, code):
    """One marker of diameter d at (X, Y) drawn straight into the image (for dots so big that few are on screen)."""
    marker, is_open = SHAPE_MARKERS[code % 10], code >= 10
    thick = max(1, d // 12)
    if marker == 'o':
        cv2.circle(img, (X, Y), d // 2, c, thick if is_open else -1)
    else:
        pts = np.round((np.array(SHAPE_POLYGONS[marker]) - 0.5) * d + (X, Y)).astype(np.int32)
        if is_open:
            cv2.polylines(img, [pts], True, c, thick)
        else:
            cv2.fillPoly(img, [pts], c)


def shape_outline(code, cx, cy, size):
    """The flat (x, y, ...) corners of the marker of this code, `size` pixels across and centered on (cx, cy), as drawn into the pictures
    (a disc: 32 corners)."""
    marker = SHAPE_MARKERS[code % 10]
    if marker == 'o':
        pts = [(0.5 + 0.5 * np.cos(a), 0.5 + 0.5 * np.sin(a)) for a in np.linspace(0, 2 * np.pi, 32, endpoint=False)]
    else:
        pts = SHAPE_POLYGONS[marker]
    return [c for x, y in pts for c in (cx + (x - 0.5) * size, cy + (y - 0.5) * size)]


def draw_marker(canvas, cx, cy, code, color, size=12):
    """A marker of this code, `size` pixels across, on a Tk canvas (a legend swatch)."""
    marker, is_open = SHAPE_MARKERS[code % 10], code >= 10
    h = size / 2
    style = {"fill": "", "outline": color, "width": 1.5} if is_open else {"fill": color, "outline": ""}
    if marker == 'o':
        canvas.create_oval(cx - h, cy - h, cx + h, cy + h, **style)
    else:
        canvas.create_polygon(*[c for x, y in SHAPE_POLYGONS[marker] for c in (cx + (x - 0.5) * size, cy + (y - 0.5) * size)], **style)


def stamp_plus(img, px, py, rgb, arm, thick=2):
    """A '+' (`thick` pixels thick, two by default, `arm` pixels each way) at each point, in place. `rgb`: one color (3,) or one per point (n, 3)."""
    h, w = img.shape[:2]
    per_point = np.ndim(rgb) == 2
    first = 1 - thick // 2                                       # (the first of the thick rows, from the point: 0 and 1 for two)
    if arm > 12:                                                 # (a big '+', from a close zoom, where few points are in view: a bar or two for each)
        for i, (X, Y) in enumerate(zip(np.asarray(px).tolist(), np.asarray(py).tolist())):
            c = rgb[i] if per_point else rgb
            x0, x1, ya, yb = max(X - arm, 0), min(X + arm + 1, w), max(Y + first, 0), min(Y + first + thick, h)
            if x1 > x0 and yb > ya:
                img[ya:yb, x0:x1] = c
            y0, y1, xa, xb = max(Y - arm, 0), min(Y + arm + 1, h), max(X + first, 0), min(X + first + thick, w)
            if y1 > y0 and xb > xa:
                img[y0:y1, xa:xb] = c
        return
    for k in range(-arm, arm + 1):
        for j in range(first, first + thick):
            for dx, dy in ((k, j), (j, k)):
                xx, yy = px + dx, py + dy
                ok = (xx >= 0) & (xx < w) & (yy >= 0) & (yy < h)
                img[yy[ok], xx[ok]] = rgb[ok] if per_point else rgb


def blend_dots(img, px, py, rgb, d, opacity, kernel_cache, kernel=None):
    """Draw dots of diameter `d` at (px, py) (already shifted by `d`, as in PointCanvas.draw) on `img` with the given opacity, so that dots
    piled on each other show it: where n dots overlap, the pixel is their average color covering the background with 1 - (1 - opacity) ** n."""
    h, w = img.shape[:2]
    W, H = w + 2 * d, h + 2 * d
    keep = (px >= 0) & (px < W) & (py >= 0) & (py < H)       # (a point just past the margin would land on the wrong row, or past the end)
    if not keep.all():
        px, py, rgb = px[keep], py[keep], rgb[keep]
    flat = py * W + px
    size = H * W
    st = np.empty((H, W, 4), np.float32)
    for c in range(3):
        st[..., c] = np.bincount(flat, weights=rgb[:, c], minlength=size).reshape(H, W)
    st[..., 3] = np.bincount(flat, minlength=size).reshape(H, W)
    if d > 1:
        if kernel is not None:                             # (a marker of some other shape than a disc)
            k = kernel.astype(np.float32)
        else:
            k = kernel_cache.get(d)
            if k is None:
                k = kernel_cache[d] = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (d, d)).astype(np.float32)
        st = cv2.filter2D(st, -1, k)
    st = st[d:d + h, d:d + w]
    n = st[..., 3]
    hit = n > 0.5
    if not hit.any():
        return
    n = n[hit]
    cover = (1.0 - (1.0 - opacity) ** n)[:, None]
    mean = st[hit][:, :3] / n[:, None]
    img[hit] = (img[hit] * (1.0 - cover) + mean * cover + 0.5).astype(np.uint8)


class GeneBox(tk.Text):
    """A one-line text box, looking like an entry, in which the gene names that were rejected can be shown in red (an entry can only color
    all of its text, or none)."""
    BAD = "bad"

    def __init__(self, parent):
        super().__init__(parent, height=1, wrap="none", undo=True, font="TkDefaultFont", relief="solid", borderwidth=1, padx=4, pady=3,
                         highlightthickness=1, highlightcolor="#4a90d9", highlightbackground="#a0a0a0")
        self.tag_configure(self.BAD, foreground="#d00000")
        self.bind("<Tab>", lambda e: (self.tk_focusNext().focus_set(), "break")[1])
        self.bind("<Control-a>", lambda e: (self.tag_add("sel", "1.0", "end-1c"), "break")[1])

    def get_text(self):
        return self.get("1.0", "end-1c")

    def set_text(self, text):
        self.delete("1.0", tk.END)
        self.insert("1.0", text)

    def clear_bad(self):
        self.tag_remove(self.BAD, "1.0", tk.END)

    def mark_bad(self, names):
        """Show the gene names `names` (as typed, any case) in red wherever they are in the box."""
        self.clear_bad()
        bad = {n.lower() for n in names}
        pos = 0
        for part in self.get_text().split(","):
            word = part.strip()
            if word and word.lower() in bad:
                start = pos + part.index(word)
                self.tag_add(self.BAD, f"1.{start}", f"1.{start + len(word)}")
            pos += len(part) + 1


class GeneVar:
    """What a StringVar is to an entry, for a GeneBox: get() and set()."""

    def __init__(self, box):
        self.box = box

    def get(self):
        return self.box.get_text()

    def set(self, text):
        self.box.set_text(text)


def label_weights(counts):
    """The size factor of each group's number from its cell count: 0.6 for the smallest group to 2.0 for the largest, in between by the
    logarithm of the count (so the major groups stand out when there are many)."""
    counts = np.maximum(np.asarray(counts, dtype=float), 1.0)
    if len(counts) == 0 or counts.max() <= counts.min():
        return np.ones(len(counts))
    t = (np.log(counts) - np.log(counts.min())) / (np.log(counts.max()) - np.log(counts.min()))
    return LABEL_WEIGHT_RANGE[0] + (LABEL_WEIGHT_RANGE[1] - LABEL_WEIGHT_RANGE[0]) * t


def expressing_percentile(vals, point, info=None, hi=None):
    """Where `point` stands among the cells that express the gene (a finite value above 0): the percentage of them below it (0 if none express it). It is
    read from the same distribution the histogram draws (see histogram_columns; `info` is its data, expression_histogram(vals), `hi` the top of the range):
    where the values are whole, each value's cells are spread over the stretch it stands for, so a point inside that stretch has part of them below it, as the
    histogram shows, and the reading is the usual percentile rank of a value, counting half of the cells that have it. The cells with no expression are left
    out, as they are below any point above 0 and would only say how few cells express the gene."""
    if info is not None and info[1] is not None and hi is not None:
        u, c = info[1]
        return 100.0 * float(np.interp(point, discrete_bounds(u, hi), np.concatenate(([0.0], np.cumsum(c))), left=0.0, right=float(c.sum()))) / float(c.sum())
    expr = np.isfinite(vals) & (vals > 0)
    ne = int(np.count_nonzero(expr))
    return 100.0 * np.count_nonzero(expr & (vals < point)) / ne if ne else 0.0


def value_at_percentile(info, hi, percent):
    """The value that is the `percent`th percentile of the cells that express the gene (the inverse of expressing_percentile), or None if none do."""
    srt, discrete = info
    if not len(srt):
        return None
    if discrete is not None:
        u, c = discrete
        return float(np.interp(percent / 100.0 * c.sum(), np.concatenate(([0.0], np.cumsum(c))), discrete_bounds(u, hi)))
    return float(np.quantile(srt, percent / 100.0))


HIST_HEIGHT = 34                    # the histogram under a gene's levels slider (pixels)
LINEAR_WHITE_PERCENTILE = 98.0      # in the linear scale the white point of a gene starts at this percentile of the cells that express it
HIST_PAD = 7                        # and room under its axis for the dot that shows the value of the cell the pointer is on (with its halo)
HIST_MAX_DISCRETE = 600             # a gene with no more different values than this (counts are whole numbers: a few dozen) has its cells spread over the stretch each value stands for


def expression_histogram(vals):
    """What a histogram of a gene's expression is drawn from: (the values of the cells that express it (finite, above 0), sorted; the different values and
    how many cells have each, if there are few of them (else None)). The cells with no expression are left out: there are far more of them than of any
    other value, and they would be all there is to see."""
    v = vals[np.isfinite(vals) & (vals > 0)]
    srt = np.sort(v)
    discrete = None
    if len(srt):
        u, c = np.unique(srt, return_counts=True)
        if len(u) <= HIST_MAX_DISCRETE:
            discrete = (u, c)
    return srt, discrete


def discrete_bounds(u, hi):
    """The stretches of the values `u` (few and whole: counts, or the log of 1 + counts) as the histogram spreads them: the edges halfway between one value and the
    next (from 0 below the first), len(u) + 1 of them; the highest value's stretch ends where the track does (`hi`), so that no cell is cut off."""
    points = np.concatenate(([0.0], u))
    end = min(u[-1] + (u[-1] - (u[-2] if len(u) > 1 else 0.0)) / 2, hi)
    return np.concatenate(((points[:-1] + points[1:]) / 2, [end]))


def histogram_columns(info, lo, hi, even, n):
    """How many of the cells in `info` (see expression_histogram) fall in each of `n` equal parts of a slider's track, the track running over the values
    from `lo` to `hi` (a position is the fraction of the range, or its square root unless `even`: see RangeSlider). Where the values are few and whole
    (counts; the log of 1 + counts), each value's cells are spread over the stretch halfway to its neighbors, so that the track is filled and not a
    comb of spikes: a column narrower than that holds its share."""
    srt, discrete = info
    if not len(srt) or hi <= lo:
        return np.zeros(n)
    position = np.linspace(0.0, 1.0, n + 1)
    edges = lo + (position if even else position ** 2) * (hi - lo)
    if discrete is not None:
        u, c = discrete
        return np.diff(np.interp(edges, discrete_bounds(u, hi), np.concatenate(([0.0], np.cumsum(c))), left=0.0, right=float(c.sum())))
    below = np.searchsorted(srt, edges, side="left")
    below[-1] = len(srt)
    return np.diff(below).astype(float)


class RangeSlider(tk.Canvas):
    """A track with two handles, for the lowest and the highest end of a range (0 to 1): the black level (left, a black dot) and the white level
    (right, a white dot) of a color scale. It is as wide as the color bar above it (the track runs between the same margins), calls
    `command(low, high)` as a handle is dragged, and a double-click puts both back at the ends. `low` and `high` are the fractions of the range
    the levels are at; the track is drawn by their square root (a position is the square root of its value), so that the low end, where most of a
    gene's expression is, has most of the track, unless `even` is set (a log scale, which already spreads the values): then by the fractions themselves."""
    MIN_GAP = 0.005           # the handles keep at least this much of the track between them

    def __init__(self, parent, command, margin=10, height=26):
        super().__init__(parent, height=height, highlightthickness=0)
        self.low, self.high, self.margin, self.command, self.dragging = 0.0, 1.0, margin, command, None
        self.even = False                                  # (True: the track is even, a position is the fraction itself)
        self.default = (0.0, 1.0)                          # (where a double-click puts the handles: the whole range, or what the window sets)
        self.bind("<Configure>", lambda e: self.draw())
        self.bind("<ButtonPress-1>", self.press)
        self.bind("<B1-Motion>", self.move)
        self.bind("<ButtonRelease-1>", lambda e: setattr(self, "dragging", None))
        self.bind("<Double-Button-1>", lambda e: self.set(*self.default, notify=True))

    def x_of(self, fraction):
        fraction = min(max(fraction, 0.0), 1.0)
        return self.margin + (fraction if self.even else math.sqrt(fraction)) * max(self.winfo_width() - 2 * self.margin, 1)

    def fraction_of(self, x):
        position = min(max((x - self.margin) / max(self.winfo_width() - 2 * self.margin, 1), 0.0), 1.0)
        return position if self.even else position * position

    def set_even(self, on):
        """An even track (for a log scale) or the square-root one (linear): the handles stay at the same fractions, at other places of the track."""
        if bool(on) != self.even:
            self.even = bool(on)
            self.draw()

    def set(self, low, high, notify=False):
        self.low, self.high = low, high
        self.draw()
        if notify:
            self.command(low, high)

    def draw(self):
        self.delete("all")
        y, r = 11, 6
        x0, x1 = self.x_of(0), self.x_of(1)
        a, b = self.x_of(self.low), self.x_of(self.high)
        self.create_line(x0, y, x1, y, fill="#c8c8c8", width=3, capstyle="round")
        self.create_line(a, y, b, y, fill="#707070", width=3)
        self.create_oval(a - r, y - r, a + r, y + r, fill="black", outline="#707070", width=1)
        self.create_oval(b - r, y - r, b + r, y + r, fill="white", outline="black", width=1)

    def press(self, e):
        a, b = self.x_of(self.low), self.x_of(self.high)
        self.dragging = "low" if abs(e.x - a) <= abs(e.x - b) and not (a == b and e.x > a) else "high"
        self.move(e)

    def move(self, e):
        if self.dragging is None:
            return
        f = self.fraction_of(e.x)
        if self.dragging == "low":
            self.low = min(f, self.high - self.MIN_GAP)
        else:
            self.high = max(f, self.low + self.MIN_GAP)
        self.low, self.high = max(self.low, 0.0), min(self.high, 1.0)
        self.draw()
        self.command(self.low, self.high)


class HoverTip:
    """A small tip beside the pointer while it is over `widget`; `text_fn()` gives the text (asked again as the pointer moves)."""

    def __init__(self, widget, text_fn):
        self.widget, self.text_fn, self.tip, self.label = widget, text_fn, None, None
        widget.bind("<Enter>", self.show, add="+")
        widget.bind("<Motion>", self.show, add="+")
        widget.bind("<Leave>", self.hide, add="+")
        widget.bind("<ButtonPress>", self.show, add="+")

    def show(self, e):
        text = self.text_fn()
        if not text:                                          # (nothing to say where the pointer is now: no tip)
            self.hide()
            return
        if self.tip is None:
            self.tip = tk.Toplevel(self.widget)
            self.tip.wm_overrideredirect(True)
            self.label = tk.Label(self.tip, bg="#ffffe0", relief="solid", borderwidth=1, font=("Segoe UI", 9), padx=4, pady=1, justify=tk.LEFT, wraplength=460)
            self.label.pack()
        self.label.config(text=text)
        # beside the pointer, but kept on the screen: to its left if it would run off the right edge, above it if off the bottom
        self.tip.update_idletasks()
        w, h = self.tip.winfo_reqwidth(), self.tip.winfo_reqheight()
        sw, sh = self.widget.winfo_screenwidth(), self.widget.winfo_screenheight()
        x, y = e.x_root + 12, e.y_root + 18
        if 0 <= e.x_root < sw and x + w > sw:
            x = max(e.x_root - w - 12, 0)
        if 0 <= e.y_root < sh and y + h > sh:
            y = max(e.y_root - h - 8, 0)
        self.tip.wm_geometry(f"+{x}+{y}")

    def hide(self, e=None):
        if self.tip is not None:
            self.tip.destroy()
            self.tip = self.label = None


PIN_OUTLINE = [(5, 2), (11, 2), (10, 8), (13, 10), (13, 12), (3, 12), (3, 10), (6, 8)]       # a push-pin's head and body, upright in a 16 x 20 box


def pin_scale(widget):
    """How much bigger than 16 x 20 pixels the pin is drawn: the screen's pixels per inch over 72, which is what a Mac counts in (its pin looked the
    right size at 16 x 20); a Windows screen is 96 pixels per inch or more, and a pin of the same pixels was small."""
    try:
        return max(1.0, widget.winfo_fpixels("1i") / 72.0)
    except tk.TclError:
        return 1.0


def draw_pin(canvas, pinned, scale=None):
    """A push-pin in `canvas` (16 x 20 times `scale`: see pin_scale): upright, red and filled while the tip is pinned; tilted, gray and hollow while it is not."""
    import math
    scale = pin_scale(canvas) if scale is None else scale
    canvas.delete("all")
    angle = 0.0 if pinned else math.radians(40)
    cx, cy = 8, 11

    def turn(x, y):
        dx, dy = x - cx, y - cy
        return (cx + dx * math.cos(angle) - dy * math.sin(angle)) * scale, (cy + dx * math.sin(angle) + dy * math.cos(angle)) * scale

    fill, line = ("#d03030", "#802020") if pinned else ("#ffffe0", "#808080")
    needle = [turn(8, 12), turn(8, 19)]
    canvas.create_line(*needle[0], *needle[1], fill=line, width=max(2, round(2 * scale)))
    canvas.create_polygon([c for pt in PIN_OUTLINE for c in turn(*pt)], fill=fill, outline=line, width=max(1, round(scale)))


class CellTip:
    """The hovered cell's description in a tip of its own, once the pointer has rested on one cell for
    CELL_TIP_MS, with a button that copies it. Unlike a HoverTip it stays while the pointer goes onto it (to reach the button): it goes when
    the pointer is on another cell more than CELL_TIP_GRACE_PX from where the tip came, or leaves the tip itself."""

    def __init__(self, root, delay_ms=None):
        self.root, self.tip, self.job, self.key, self.anchor = root, None, None, None, None
        self.delay_ms = CELL_TIP_MS if delay_ms is None else delay_ms
        self.drag = None               # while a tip is being dragged: {"top", "dx", "dy", "moved"}, so that it is not taken away from under the pointer
        self.pinned = []               # the tips that were pinned (by moving them or by the pin): they stay, with a close box, until closed or unpinned
        self.bars = {}                 # tip -> (its top bar, its pin picture)
        self.hot = []                  # the rows that react to the pointer: (tip, token, the labels to light, the widgets that make up its area)
        self.poll_job = None           # the pointer is looked at while a tip is up (see poll)
        self.entered = False           # the pointer has been on the floating tip since it came up
        self.last_moved = False        # the last press on a tip went on to move it (so that it was not a click on the pin)
        self.suspended = False         # a plot or the pane is being dragged: no tip, and no waiting for one, until it stops (suspend, resume)
        self.row_hover_fn = None       # fn(token) when the pointer rests on a row that has a token (a 4th item), fn(None) when it leaves it
        self.row_active = self.row_labels = self.row_job = None       # the row the pointer is on: its token and its labels (shown highlighted)

    def pointer_on(self, key, rows_fn, action_fn=None, swatch_fn=None, copy_text_fn=None):
        """The pointer is on a new cell (`key`); `rows_fn()` gives its description, as rows (what, number, name), asked when the tip is shown.
        `action_fn()`, if given, is asked then too: None, or (text, command) for a button beside Copy. So is `swatch_fn()`: None, or
        (rgb, marker code, background rgb), the cell as it is drawn, shown at the top (where cells are close together, which one the
        description is of is not always plain from the ring alone). `copy_text_fn()`, if given, is the text the Copy button copies (else
        the rows, a tab between their columns)."""
        if self.suspended:
            return
        if key == self.key:
            return
        self.key = key
        self.cancel()
        if self.tip is not None and not self.pointer_in_tip() and not self.near_anchor():
            self.hide()
        self.job = self.root.after(self.delay_ms, lambda: self.show(rows_fn, action_fn, swatch_fn, copy_text_fn))

    def suspend(self):
        """A drag has begun: the tip that is up goes (not the ones that were moved), the wait for one is dropped, and none starts until resume."""
        self.suspended = True
        self.cancel()
        self.key = None
        if self.tip is not None and not self.pointer_in_tip():
            self.hide()

    def resume(self):
        self.suspended = False

    def pointer_off(self):
        """The pointer is on no cell (or left the plot): the tip goes, unless the pointer went onto it."""
        self.key = None
        self.cancel()
        if self.tip is not None:
            self.root.after(300, self.hide_unless_pointer_in)

    def cancel(self):
        if self.job is not None:
            self.root.after_cancel(self.job)
            self.job = None

    def pointer_in_tip(self):
        if self.tip is None:
            return False
        try:
            x, y = self.root.winfo_pointerxy()
            tx, ty = self.tip.winfo_rootx(), self.tip.winfo_rooty()
            return tx <= x < tx + self.tip.winfo_width() and ty <= y < ty + self.tip.winfo_height()
        except tk.TclError:
            return False

    def near_anchor(self):
        if self.anchor is None:
            return False
        x, y = self.root.winfo_pointerxy()
        return abs(x - self.anchor[0]) <= CELL_TIP_GRACE_PX and abs(y - self.anchor[1]) <= CELL_TIP_GRACE_PX

    def hide_unless_pointer_in(self):
        if self.drag is None and not self.pointer_in_tip():
            self.hide()

    def make_movable(self, widget):
        """The tip can be dragged by any of its text (not its buttons) with the left button. Once moved it is pinned (see pin)."""
        if isinstance(widget, ttk.Button):
            return
        widget.config(cursor="fleur")
        widget.bind("<ButtonPress-1>", self.drag_start, add="+")
        widget.bind("<B1-Motion>", self.drag_to, add="+")
        widget.bind("<ButtonRelease-1>", self.drag_end, add="+")
        for child in widget.winfo_children():
            self.make_movable(child)

    ROW_COLOR = "#ffe08a"

    POLL_MS = 80

    def start_poll(self):
        if self.poll_job is None:
            self.poll_job = self.root.after(self.POLL_MS, self.poll)

    def poll(self):
        """While a tip is up, look at where the pointer is, instead of waiting to be told: on macOS a tip (a window without a frame, never the
        active one) is sent no Enter, Leave or Motion events until it is clicked, so its rows would not react and it would not go when the pointer
        left it. The rows are lit by where the pointer is over them, and the floating tip goes when the pointer has been on it and then leaves."""
        self.poll_job = None
        live = []
        for h in self.hot:
            try:
                if h[0].winfo_exists():
                    live.append(h)
            except tk.TclError:
                pass
        self.hot = live
        if not live and self.tip is None and not self.pinned:
            return
        try:
            x, y = self.root.winfo_pointerxy()
            hit = None
            for tip, token, labels, area in live:
                boxes = [(w.winfo_rootx(), w.winfo_rooty(), w.winfo_rootx() + w.winfo_width(), w.winfo_rooty() + w.winfo_height())
                         for w in area if w.winfo_ismapped()]
                if boxes and min(b[0] for b in boxes) <= x < max(b[2] for b in boxes) and min(b[1] for b in boxes) <= y < max(b[3] for b in boxes):
                    hit = (token, labels)
                    break
            if hit is not None and self.row_hover_fn is not None and self.drag is None:
                self.row_enter(*hit)
            elif self.row_active is not None:
                self.row_leave()
            if self.tip is not None:
                if self.pointer_in_tip():
                    self.entered = True
                elif self.entered and self.drag is None:
                    self.entered = False
                    self.hide()
        except tk.TclError:
            pass
        self.start_poll()

    def row_enter(self, token, labels):
        """The pointer is on a row of the tip that stands for a group of cells: the row is lit, and the cells are shown (row_hover_fn)."""
        if self.row_job is not None:
            self.root.after_cancel(self.row_job)
            self.row_job = None
        if token == self.row_active and labels[0] is self.row_labels[0]:
            return
        self.row_unlight()
        self.row_active, self.row_labels = token, labels
        for w in labels:
            w.config(bg=self.ROW_COLOR)
        self.row_hover_fn(token)

    def row_leave(self):
        """(a little later: going from one of the row's labels to the next leaves and enters again at once)"""
        if self.row_job is None and self.row_active is not None:
            self.row_job = self.root.after(60, self.row_left)

    def row_left(self):
        self.row_job = None
        self.row_unlight()

    def row_unlight(self):
        if self.row_active is None:
            return
        labels, self.row_active, self.row_labels = self.row_labels, None, None
        for w in labels:
            try:
                w.config(bg="#ffffe0")
            except tk.TclError:
                pass
        if self.row_hover_fn is not None:
            self.row_hover_fn(None)

    def drag_start(self, e):
        top = e.widget.winfo_toplevel()
        self.drag = {"top": top, "dx": e.x_root - top.winfo_x(), "dy": e.y_root - top.winfo_y(), "x": e.x_root, "y": e.y_root, "moved": False}

    def drag_to(self, e):
        d = self.drag
        if d is None:
            return
        if abs(e.x_root - d["x"]) > 3 or abs(e.y_root - d["y"]) > 3:
            d["moved"] = True
        if d["moved"]:
            d["top"].wm_geometry(f"+{e.x_root - d['dx']}+{e.y_root - d['dy']}")

    def drag_end(self, e):
        d, self.drag = self.drag, None
        self.last_moved = bool(d is not None and d["moved"])
        if self.last_moved and d["top"] is self.tip:
            self.pin()

    def pin_clicked(self, tip):
        """A click on the pin (not a drag that began on it): a tip that is up becomes pinned; a pinned one is unpinned, and so goes."""
        if self.last_moved:
            self.last_moved = False
            return
        if tip is self.tip:
            self.pin()
        elif tip in self.pinned:
            self.close_pinned(tip)

    def pin(self):
        """The tip was moved or its pin clicked: it stays where it is until its close box is clicked or it is unpinned, whatever the pointer
        does, and a new tip can come up for the next cell beside it. It can still be moved."""
        tip = self.tip
        self.tip = self.anchor = None
        self.key = None
        tip.unbind("<Leave>")
        bar, pin = self.bars[tip]
        draw_pin(pin, True)
        close = tk.Label(bar, text="\u2715", bg="#ffffe0", fg="gray35", font=("Segoe UI", 9, "bold"), padx=5, cursor="hand2")
        close.pack(side=tk.RIGHT)
        close.bind("<Enter>", lambda e: close.config(bg="#e04040", fg="white"))
        close.bind("<Leave>", lambda e: close.config(bg="#ffffe0", fg="gray35"))
        close.bind("<ButtonRelease-1>", lambda e: self.close_pinned(tip))
        self.pinned.append(tip)

    def close_pinned(self, tip):
        self.row_unlight()
        self.bars.pop(tip, None)
        if tip in self.pinned:
            self.pinned.remove(tip)
        try:
            tip.destroy()
        except tk.TclError:
            pass

    def show(self, rows_fn, action_fn=None, swatch_fn=None, copy_text_fn=None):
        self.job = None
        if self.pointer_in_tip():                       # (on its way to the Copy button: the tip it is going to stays as it is)
            return
        rows = rows_fn()
        if not rows:
            return
        self.hide()
        self.tip = tk.Toplevel(self.root)
        self.tip.wm_overrideredirect(True)
        self.tip.attributes("-topmost", True)
        bg = "#ffffe0"
        box = tk.Frame(self.tip, bg=bg, relief="solid", borderwidth=1)
        box.pack()
        bar = tk.Frame(box, bg=bg)
        bar.pack(side=tk.TOP, fill=tk.X, padx=3, pady=(2, 0))
        scale = pin_scale(bar)
        pin = tk.Canvas(bar, width=round(16 * scale), height=round(20 * scale), bg=bg, highlightthickness=0, cursor="hand2")
        pin.pack(side=tk.LEFT)
        draw_pin(pin, False)
        self.bars[self.tip] = (bar, pin)
        swatch = swatch_fn() if swatch_fn is not None else None
        if swatch is not None:                          # the cell as it is drawn: its color and its group's marker, on the plot's background
            rgb, code, back = swatch[:3]
            token = swatch[3] if len(swatch) > 3 else None         # (what the row_hover_fn is told while the pointer is on this: the cell itself)
            hexc = "#%02x%02x%02x" % tuple(int(c) for c in rgb[:3])
            top = tk.Frame(box, bg=bg)
            top.pack(anchor="e", padx=6, pady=(5, 0))               # (at the right: the tip opens with its top left corner at the pointer, so the symbol was under it)
            caption = tk.Label(top, text="this cell, as drawn", bg=bg, fg="gray35", font=("Segoe UI", 9))
            caption.pack(side=tk.LEFT, padx=(0, 6))
            mark = tk.Canvas(top, width=22, height=22, bg="#%02x%02x%02x" % tuple(int(c) for c in back[:3]),
                             highlightthickness=1, highlightbackground="gray60")
            mark.pack(side=tk.LEFT)
            draw_marker(mark, 12, 12, int(code), hexc, size=14)
            bring = swatch[4] if len(swatch) > 4 else None
            if bring is not None:                       # (to the right of the symbol: the cell comes into view in the plots where it is out of it)
                button = ttk.Button(top, text="Bring into view", command=bring)
                button.pack(side=tk.LEFT, padx=(8, 0))
                if len(swatch) > 5 and not swatch[5]:    # (greyed out while the UMAP is not shown)
                    button.state(["disabled"])
                    HoverTip(button, lambda: "Show the UMAP (the Show UMAP box above the section maps) to use this.")
            if token is not None and self.row_hover_fn is not None:
                for w in (top, caption, mark):
                    w.bind("<Enter>", lambda e, t=token, c=caption: self.row_enter(t, [c]), add="+")
                    w.bind("<Leave>", lambda e: self.row_leave(), add="+")
                self.hot.append((self.tip, token, [caption], [top]))
        # in columns: what (gray), the number (bold, to the right), the name
        grid = tk.Frame(box, bg=bg)
        grid.pack(anchor="w", padx=6, pady=(4, 2))
        for r, row in enumerate(rows):
            what, num, name = row[:3]
            token = row[3] if len(row) > 3 else None
            tk.Label(grid, text=what, bg=bg, fg="gray35", font=("Segoe UI", 9, "bold"), anchor="w", justify=tk.LEFT,
                     wraplength=HOVER_LABEL_PX).grid(row=r, column=0, sticky="nw", padx=(0, HOVER_GAP_PX))
            tk.Label(grid, text=num, bg=bg, font=("Segoe UI", 9, "bold"), anchor="e").grid(row=r, column=1, sticky="ne", padx=(0, HOVER_GAP_PX))
            tk.Label(grid, text=name, bg=bg, font=("Segoe UI", 9), anchor="w", justify=tk.LEFT, wraplength=300).grid(row=r, column=2, sticky="nw")
            if token is not None and self.row_hover_fn is not None:
                labels = [w for w in grid.grid_slaves(row=r)]
                for w in labels:
                    w.bind("<Enter>", lambda e, t=token, ls=labels: self.row_enter(t, ls), add="+")
                    w.bind("<Leave>", lambda e: self.row_leave(), add="+")
                self.hot.append((self.tip, token, labels, labels))
        # copied with a tab between the columns: pasted into a spreadsheet, they stay columns
        text = copy_text_fn() if copy_text_fn is not None else "\n".join("\t".join(row[:3]) for row in rows)
        buttons = tk.Frame(box, bg=bg)
        buttons.pack(anchor="e", padx=6, pady=(0, 6))
        action = action_fn() if action_fn is not None else None
        if action is not None:
            label, command = action
            ttk.Button(buttons, text=label, command=lambda: (self.hide(), command())).pack(side=tk.LEFT, padx=(0, 6))
        copy = ttk.Button(buttons, text="Copy to clipboard")
        copy.config(command=lambda: self.copy(text, copy))
        copy.pack(side=tk.LEFT)
        self.tip.bind("<Leave>", lambda e: self.root.after(100, self.hide_unless_pointer_in))
        self.make_movable(box)
        self.entered = False
        self.start_poll()
        pin.bind("<ButtonRelease-1>", lambda e, t=self.tip: self.pin_clicked(t), add="+")        # (after the drag handlers: they say whether this press moved the tip)
        pin.config(cursor="hand2")
        # just beside the pointer (close, so that it is quickly reached), but kept on the screen
        x0, y0 = self.root.winfo_pointerxy()
        self.anchor = (x0, y0)
        self.tip.update_idletasks()
        w, h = self.tip.winfo_reqwidth(), self.tip.winfo_reqheight()
        sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        x, y = x0 + 6, y0 + 6
        if x + w > sw:
            x = max(x0 - w - 6, 0)
        if y + h > sh:
            y = max(y0 - h - 6, 0)
        self.tip.wm_geometry(f"+{x}+{y}")

    def copy(self, text, button):
        self.root.clipboard_clear()
        self.root.clipboard_append(text)
        button.config(text="Copied")
        self.root.after(1200, lambda: button.winfo_exists() and button.config(text="Copy to clipboard"))

    def hide(self):
        self.row_unlight()
        if self.tip is not None:
            self.bars.pop(self.tip, None)
            try:
                self.tip.destroy()
            except tk.TclError:
                pass
            self.tip = self.anchor = None


class LazyLevels(dict):
    """The group labels (as text) of every cell of the atlas at each level, made the first time a level is asked for: making all of them takes
    half a second, and a run is shown at one level at a time. Used like a dict of arrays: `levels.get(level)` (None if the table has no such level)."""

    def __init__(self, table, names):
        super().__init__()
        self.table, self.names = table, list(names)

    def get(self, key, default=None):
        if key not in self:
            dict.__setitem__(self, key, self.table[key].astype(str).to_numpy() if key in self.names and key in self.table.columns else None)
        value = dict.__getitem__(self, key)
        return default if value is None else value

    def __getitem__(self, key):
        value = self.get(key)
        if value is None:
            raise KeyError(key)
        return value


WINDOW_MARGIN = 0.15                                   # a zoomed section map is made this far (times the pane's size) beyond what is on screen, so that a small move needs no new picture


def grown(d, old, new):
    """An offset `d` from the corner of something of size `old` that is now of size `new`: inside it, it grows with it; before it or beyond it, it stays."""
    return d if d < 0 else d * new / old if d <= old else new + (d - old)


def roi_dash(wd):
    """The dash option for an ROI outline wd pixels wide. Windows draws (3 * wd, 2 * wd) as short dashes, 5 to 10 along a side; macOS
    draws it as asked, one or two dashes a side, so there the dashes are given about as long as the line is wide."""
    if wd <= 2:
        return dash((6, 3), wd)
    if sys.platform == "darwin":
        return dash((max(4, wd + 1), max(3, wd * 2 // 3)), wd)
    return (3 * wd, 2 * wd)


class PointCanvas:
    """Points drawn straight onto a Tk canvas: wheel zoom about the pointer, left- or right-drag pan, never out past the whole (`home`)
    view and never off its edges. `layers` are drawn in order; each is {'x','y','rgb','size'}. Markers (the hover ring), rectangles
    and text labels are canvas items given in data coordinates and repositioned after each draw."""

    MAX_ZOOM = 400.0

    def __init__(self, parent, home, flip_y, width=300, height=300):
        self.canvas = tk.Canvas(parent, width=width, height=height, bg="white", highlightthickness=1, highlightbackground="#bbbbbb")
        self.home = home                            # (xmin, xmax, ymin, ymax)
        self.flip_y = flip_y                        # True: data y grows upward on screen (the UMAP); False: downward (a brain section)
        self.layers = []
        self.cx = (home[0] + home[1]) / 2
        self.cy = (home[2] + home[3]) / 2
        self.scale = 1.0
        self.at_fit = True
        self.photo = None
        self.pending = None
        self.item = self.canvas.create_image(0, 0, anchor="nw")
        self.markers = {}                           # name -> [item, x, y] (x None: hidden)
        self.shaped_key = self.shaped_found = None  # the last cell looked for by shaped_cell, and what was found
        self.rects = []                             # [item, x0, y0, x1, y1]
        self.rect_line = None                       # the width the ROI outlines have been given (see rect_width)
        self.rect_px = 2                            # the width they should have
        self.labels = []                            # [item, x, y]
        self.label_text = {}                        # item -> its text (the group's number)
        self.label_weight = {}                      # item -> its size factor, by the group's cell count
        self.label_color = {}                       # item -> the colour it was given (for a light background: see ink_for)
        self.scale_labels = True                    # the labels' size follows their weight
        self.shape_kernels = {}                     # (marker code, diameter) -> the kernel that grows a point into that marker
        self.shapes_while_moving = False            # True: layers with shapes keep them while the view moves (for small canvases: drawing twice is the slower way)
        self.moving = False                         # the view is being zoomed or panned: layers with shapes are drawn as plain discs
        self.settle_job = None
        self.boosted = None                         # the text of the label drawn bigger now (the highlighted group's), or None
        self.label_factor = 1.0                     # the Label size slider
        self.label_size = 8                         # the group labels' font size: grows with the zoom up to LABEL_ZOOM_CAP, then stops
        self.scalebar = None                        # (line item, text item) once enabled: data units are mm
        self.on_motion = None                       # fn(cx, cy) ; fn(None, None) on leave
        self.on_view = None                         # fn() after a redraw
        self.on_change = None                       # fn(zoom, dx, dy) after the user zoomed or panned (the view, relative to the whole view)
        self.pan_from = None
        self.drag_hook = None                       # fn(True) when a drag of this view begins, fn(False) when it ends
        self.panel_hook = None                      # when set (the whole-panel zoom of the section maps): its wheel, drag and double-click act on the panel, not on this view
        self.panel_zoom = 1.0                       # how much the whole panel is enlarged: the dots grow with it
        self.opacity = 1.0                          # of one cell's dot, below 1 overlapping dots add up (the Opacity slider)
        self.kernels = {}                           # dot diameter -> the disc used to spread dots when they are blended
        self.marker_factor = 1.0                    # times the size of the hover ring and of the + marks
        self.dot_factor = 1.0                       # times the size of the dots that have a size on the tissue
        self.visible_fn = None                      # fn() -> whether this canvas is on screen; one that is not is only drawn once it is
        self.stale = False                          # needs drawing (it was not on screen when asked)
        self.base = None                            # the last picture of the points (H x W x 3)
        self.shown_view = None                      # (cx, cy) the items on the canvas were placed for, and the scale: a drag moves them from there (preview_pan)
        self.shown_scale = None
        self.previewed = (0, 0)                     # how far a drag has moved them since
        self.pan_job = None
        self.zoom_job = None
        self.view_cache = []                        # [(picture, cx, cy, scale)]: pictures drawn before at other zooms; they fill in what zooming out uncovers (see preview_zoom)
        self.cache_look = None                      # what the pictures in view_cache were drawn with: they are dropped when it changes
        self.labels_hidden = False
        self.drawn = None                           # the part of the canvas the picture `base` covers (x0, y0, x1, y1), in the canvas's own (virtual) pixels
        self.vsize = None                           # a virtual canvas (the whole-panel zoom): its size (w, h). The Tk widget is only as big as the part of it that is on screen
        self.vview = None                           # ... and that part: (x0, y0, w, h) of the virtual canvas, which is what the widget shows
        self.off = (0, 0)                           # where the widget's top left is in the virtual canvas
        self.border = None                          # a virtual canvas draws its own border (the widget's would be round the part shown)
        self.dim = False                            # a veil over the points (while a cell is highlighted)
        self.bgcolor = (255, 255, 255)
        self.plus_arm = 3                           # the '+' marks' arm length (pixels) at the whole view
        self.pluses = None                          # (x, y) of the red '+' marks of the highlighted group
        self.compose_pending = None
        self.outlines = None                        # closed outlines (arrays of (x, y) in data units) drawn over the points: the CCF's regions
        self.outline_px = 1                         # how wide the outlines are (pixels): twice as wide at the largest section map size
        self.cell_alpha = 1.0                       # how opaque the cells are: lowered by the CCF atlas slider, so the outlines (always opaque) stand out
        c = self.canvas
        c.bind("<Configure>", lambda e: self.on_resize())
        c.bind("<MouseWheel>", lambda e: self.zoom_at(e.x, e.y, 1.25 if e.delta > 0 else 0.8))
        c.bind("<Button-4>", lambda e: self.zoom_at(e.x, e.y, 1.25))
        c.bind("<Button-5>", lambda e: self.zoom_at(e.x, e.y, 0.8))
        for press, motion, release in (("<ButtonPress-1>", "<B1-Motion>", "<ButtonRelease-1>"), ("<ButtonPress-3>", "<B3-Motion>", "<ButtonRelease-3>"),
                                       ("<ButtonPress-2>", "<B2-Motion>", "<ButtonRelease-2>")):
            c.bind(press, self.pan_start)
            c.bind(motion, self.pan_move)
            c.bind(release, self.pan_end)
        c.bind("<Motion>", lambda e: self.on_motion and self.on_motion(e.x, e.y), add="+")
        c.bind("<Leave>", lambda e: self.on_motion and self.on_motion(None, None))
        c.bind("<Double-Button-1>", lambda e: self.panel_hook.reset() if self.panel_hook else (self.fit(), self.announce()))

    # ---- geometry

    def size(self):
        """The size of the canvas (a virtual canvas: its whole size, which the widget shows a part of)."""
        if self.vsize is not None:
            return self.vsize
        return max(self.canvas.winfo_width(), 2), max(self.canvas.winfo_height(), 2)

    def view_size(self):
        """The size of the part shown (the widget)."""
        if self.vview is not None:
            return self.vview[2], self.vview[3]
        return self.size()

    def set_virtual(self, size):
        """Make this a canvas of `size` (w, h) pixels that is too big for a widget (zoomed 100 times, say: tens of thousands of pixels, which Tk makes
        slowly and may not make at all): the widget is only the part on screen (set_view), and everything is in the big canvas's pixels. None: an
        ordinary canvas again (the caller sets the widget's size)."""
        if size is None:
            if self.vsize is None:
                return
            self.vsize = self.vview = None
            self.off = (0, 0)
            self.canvas.config(highlightthickness=1)
            if self.border is not None:
                self.canvas.delete(self.border)
                self.border = None
            return
        size = (max(int(size[0]), 2), max(int(size[1]), 2))
        if self.vsize is None:
            self.canvas.config(highlightthickness=0)
            self.border = self.canvas.create_rectangle(0, 0, 0, 0, outline="#bbbbbb")
        if size != self.vsize:
            self.vsize = size
            self.vview = None
            if self.at_fit:
                self.fit()
            else:
                self.note_motion()
                self.clamp()
                self.request_draw()

    def set_view(self, x0, y0, w, h, background=None):
        """The part of a virtual canvas the widget shows, and draws (`background`: see draw)."""
        new = (int(x0), int(y0), max(int(w), 1), max(int(h), 1))
        if new == self.vview:
            return
        old_off = self.off
        self.vview = new
        self.off = (new[0], new[1])
        # What is on the canvas stays where it is in the big canvas until it is drawn again (the picture, the hover ring, the labels): the widget is
        # in a new place, so they move the other way. Left alone they would jump (backwards, when dragging) for the moment the new picture takes.
        self.canvas.move("all", old_off[0] - new[0], old_off[1] - new[1])
        self.canvas.config(width=new[2], height=new[3])
        lw, lh = self.vsize
        self.canvas.coords(self.border, -new[0], -new[1], lw - 1 - new[0], lh - 1 - new[1])
        # Drawn now, not when the program is next idle: the widget's new place and its new picture then reach the screen together (put off, the
        # screen could show the widget in its new place with the old picture in it for a moment: a strip of it wrong).
        if self.pending is not None:
            self.canvas.after_cancel(self.pending)
            self.pending = None
        self.draw(background)

    def fit_scale(self):
        w, h = self.size()
        x0, x1, y0, y1 = self.home
        return min(w / max(x1 - x0, 1e-9), h / max(y1 - y0, 1e-9))

    def announce(self):
        if self.on_change:
            self.on_change(self.scale / self.fit_scale(), self.cx - (self.home[0] + self.home[1]) / 2, self.cy - (self.home[2] + self.home[3]) / 2)

    def note_motion(self):
        """The view is changing (zoom, pan, fit, size): until it has been still for SHAPE_SETTLE_MS the cells with shapes are drawn as plain discs,
        which is quicker (one pass instead of one for each shape), and translucent cells are drawn solid (blending them is much slower); then they are
        drawn again as they are."""
        shapes = not self.shapes_while_moving and any(lay.get("marks") is not None for lay in self.layers)
        if not shapes and self.opacity >= 1.0:
            return
        self.moving = True
        if self.settle_job is not None:
            self.canvas.after_cancel(self.settle_job)
        self.settle_job = self.canvas.after(SHAPE_SETTLE_MS, self.settle)

    def settle(self):
        self.settle_job = None
        self.moving = False
        self.request_draw()

    def follow(self, zoom, dx, dy):
        """Take the same view (zoom, and offset from the whole view's center) as the panel that was moved."""
        self.note_motion()
        self.scale = self.fit_scale() * zoom
        self.cx = (self.home[0] + self.home[1]) / 2 + dx
        self.cy = (self.home[2] + self.home[3]) / 2 + dy
        self.at_fit = zoom <= 1.0000001
        self.clamp()
        self.request_draw()

    def fit(self):
        self.note_motion()
        self.at_fit = True
        self.cx, self.cy = (self.home[0] + self.home[1]) / 2, (self.home[2] + self.home[3]) / 2
        self.scale = self.fit_scale()
        self.request_draw()

    def clamp(self):
        w, h = self.size()
        fs = self.fit_scale()
        self.scale = min(max(self.scale, fs), fs * self.MAX_ZOOM)
        x0, x1, y0, y1 = self.home
        hw, hh = w / (2 * self.scale), h / (2 * self.scale)
        self.cx = (x0 + x1) / 2 if 2 * hw >= x1 - x0 else min(max(self.cx, x0 + hw), x1 - hw)
        self.cy = (y0 + y1) / 2 if 2 * hh >= y1 - y0 else min(max(self.cy, y0 + hh), y1 - hh)

    def on_resize(self):
        if self.vsize is not None:                  # (a virtual canvas has the size it was given; the widget is only the part shown)
            return
        if self.at_fit:
            self.fit()
        else:
            self.note_motion()
            self.clamp()
            self.request_draw()

    def to_canvas(self, x, y):
        w, h = self.size()
        sy = -(y - self.cy) if self.flip_y else (y - self.cy)
        return w / 2 + (x - self.cx) * self.scale - self.off[0], h / 2 + sy * self.scale - self.off[1]

    def to_data(self, px, py):
        w, h = self.size()
        px, py = px + self.off[0], py + self.off[1]
        x = self.cx + (px - w / 2) / self.scale
        dy = (py - h / 2) / self.scale
        return x, self.cy + (-dy if self.flip_y else dy)

    def zoom_at(self, px, py, factor):
        if self.panel_hook is not None:
            self.panel_hook.zoom(self.canvas, px, py, factor)
            return
        x, y = self.to_data(px, py)
        before = self.scale
        fs = self.fit_scale()
        self.scale = min(max(self.scale * factor, fs), fs * self.MAX_ZOOM)
        if self.scale == before:
            return
        by_bitmap = self.can_preview_zoom()
        if not by_bitmap:
            self.note_motion()
        self.at_fit = abs(self.scale - fs) < 1e-9 * fs
        w, h = self.size()
        self.cx = x - (px - w / 2) / self.scale
        dy = (py - h / 2) / self.scale
        self.cy = y + dy if self.flip_y else y - dy
        self.clamp()
        if by_bitmap:
            self.preview_zoom()
        else:
            self.request_draw()
        self.announce()

    def bring_into_view(self, x, y, margin=6):
        """Put data point (x, y) in view: nothing changes if it is (`margin` pixels from the edges), else the view is centered on it, as far as the view
        may go (it never leaves the whole view). The zoom stays. Returns whether the view moved."""
        w, h = self.view_size()
        px, py = self.to_canvas(x, y)
        if margin <= px <= w - margin and margin <= py <= h - margin:
            return False
        self.note_motion()
        self.cx, self.cy = x, y
        self.at_fit = False
        self.clamp()
        self.request_draw()
        self.announce()
        return True

    def position_in_widget(self, x, y):
        """Where data point (x, y) is in the widget (to_canvas), also when the widget has not been laid out yet (a map in a pane, scrolled out of
        view, is 1 x 1 to Tk until it is shown): then, from the size it is to have and its whole view."""
        if self.vsize is None and (self.canvas.winfo_width() <= 2 or self.canvas.winfo_height() <= 2):
            w, h = self.canvas.winfo_reqwidth(), self.canvas.winfo_reqheight()
            x0, x1, y0, y1 = self.home
            scale = min(w / max(x1 - x0, 1e-9), h / max(y1 - y0, 1e-9))
            sy = -(y - (y0 + y1) / 2) if self.flip_y else (y - (y0 + y1) / 2)
            return w / 2 + (x - (x0 + x1) / 2) * scale, h / 2 + sy * scale
        return self.to_canvas(x, y)

    def remember_view(self):
        """Keep the picture just drawn, for filling in the edges when a later zoom out reaches beyond what a picture covers: the widest one drawn (the
        first, usually the whole view) and the last three."""
        if self.vsize is not None or self.moving or self.base is None:
            return
        look = (id(self.layers), self.dot_factor, self.opacity, tuple(self.bgcolor), id(self.outlines), self.cell_alpha, self.panel_zoom, self.size())
        if look != self.cache_look:
            self.cache_look, self.view_cache = look, []
        self.view_cache = [e for e in self.view_cache if abs(e[3] / self.scale - 1.0) > 0.02]
        self.view_cache.append((self.base, self.cx, self.cy, self.scale))
        if len(self.view_cache) > 4:
            widest = min(self.view_cache, key=lambda e: e[3])
            self.view_cache = [widest] + [e for e in self.view_cache if e is not widest][-3:]

    def view_coefficients(self, cx, cy, scale):
        """(a, bx, by): the picture that was drawn for the view (cx, cy, scale) seen from the current view: its pixel = a * this pixel + (bx, by)."""
        w, h = self.size()
        a = scale / self.scale
        sign = -1.0 if self.flip_y else 1.0
        return a, w / 2 + (self.cx - cx) * scale - a * w / 2, h / 2 + sign * (self.cy - cy) * scale - a * h / 2

    def can_preview_zoom(self):
        return (PREVIEW_ZOOM and self.vsize is None and self.base is not None and self.shown_view is not None and self.drawn is not None
                and self.base.shape[:2] == (self.size()[1], self.size()[0]))

    def preview_zoom(self):
        """A wheel step shows the last picture enlarged (or reduced) about the pointer, at once; the view is really drawn again when the wheel has been still
        for ZOOM_DRAW_MS. Drawing every cell at every step is what makes zooming a big plot sluggish. (Like the whole-panel zoom of the section maps.)"""
        w, h = self.size()
        a, bx, by = self.view_coefficients(self.shown_view[0], self.shown_view[1], self.shown_scale)      # (old picture pixels per new one)
        bg = tuple(int(v) for v in self.bgcolor)
        # (cv2 and a photo reused rather than made again: 4 ms where PIL's transform and a new photo took 20)
        img = cv2.warpAffine(np.ascontiguousarray(self.base), np.array([[a, 0, bx], [0, a, by]], np.float32), (w, h),
                             flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP, borderMode=cv2.BORDER_CONSTANT, borderValue=bg)
        if bx < 0 or by < 0 or a * w + bx > w or a * h + by > h:           # (zoomed out, or moved: the picture does not reach the edges)
            under = self.earlier_picture(w, h, bg)
            if under is not None:
                mask = cv2.warpAffine(np.full((self.base.shape[0], self.base.shape[1]), 255, np.uint8), np.array([[a, 0, bx], [0, a, by]], np.float32),
                                      (w, h), flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP, borderMode=cv2.BORDER_CONSTANT, borderValue=0) > 127
                under[mask] = img[mask]
                img = under
        picture = Image.fromarray(img)
        if self.photo is not None and (self.photo.width(), self.photo.height()) == (w, h):
            self.photo.paste(picture)
        else:
            self.photo = ImageTk.PhotoImage(picture, master=self.canvas)
        self.canvas.itemconfig(self.item, image=self.photo)
        self.canvas.coords(self.item, 0, 0)
        self.reposition(light=True)                                        # (the ring, the rectangles: at their places in the new view; the labels come back with the real drawing)
        if self.zoom_job is not None:
            self.canvas.after_cancel(self.zoom_job)
        self.zoom_job = self.canvas.after(ZOOM_DRAW_MS, self.request_draw)

    def earlier_picture(self, w, h, bg):
        """An earlier picture, of a wider (or the widest) view, shown at the current view, for the parts the picture being scaled does not cover; or None."""
        if not self.view_cache:
            return None
        x0, x1 = self.cx - w / (2 * self.scale), self.cx + w / (2 * self.scale)
        y0, y1 = self.cy - h / (2 * self.scale), self.cy + h / (2 * self.scale)
        best = None
        for e in self.view_cache:
            if abs(e[3] / self.shown_scale - 1.0) <= 0.02:
                continue                                                 # (the picture already being shown)
            ew, eh = e[0].shape[1], e[0].shape[0]
            covers = (e[1] - ew / (2 * e[3]) <= x0 and e[1] + ew / (2 * e[3]) >= x1 and e[2] - eh / (2 * e[3]) <= y0 and e[2] + eh / (2 * e[3]) >= y1)
            # the sharpest picture that covers the whole view; failing that, the widest there is
            rank = (1, e[3]) if covers else (0, -e[3])
            if best is None or rank > best[0]:
                best = (rank, e)
        if best is None:
            return None
        picture, ecx, ecy, escale = best[1]
        a, bx, by = self.view_coefficients(ecx, ecy, escale)
        return cv2.warpAffine(np.ascontiguousarray(picture), np.array([[a, 0, bx], [0, a, by]], np.float32), (w, h),
                              flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP, borderMode=cv2.BORDER_CONSTANT, borderValue=bg)

    def pan_start(self, e):
        if self.drag_hook is not None:
            self.drag_hook(True)
        if self.panel_hook is not None:
            self.panel_hook.pan_start(e)
            return
        self.pan_from = (e.x, e.y, self.cx, self.cy)

    def pan_move(self, e):
        if self.panel_hook is not None:
            self.panel_hook.pan_move(e)
            return
        if self.pan_from is None:
            return
        x0, y0, cx, cy = self.pan_from
        self.note_motion()
        self.cx = cx - (e.x - x0) / self.scale
        dy = (e.y - y0) / self.scale
        self.cy = cy + dy if self.flip_y else cy - dy
        self.at_fit = False
        self.clamp()
        self.preview_pan()
        self.announce()

    def preview_pan(self):
        """A drag moves what is on the canvas (the picture, the ring, the labels) as it is, at once, and draws the view again at most every PAN_DRAW_MS: drawing
        it at every mouse move, every cell in view, is what makes dragging a big plot slow."""
        if self.shown_view is None or self.scale != self.shown_scale:
            self.request_draw()
            return
        dx = (self.shown_view[0] - self.cx) * self.scale
        dy = ((self.cy - self.shown_view[1]) if self.flip_y else (self.shown_view[1] - self.cy)) * self.scale
        self.canvas.move("all", dx - self.previewed[0], dy - self.previewed[1])
        self.previewed = (dx, dy)
        if self.pan_job is None and self.pending is None:
            self.pan_job = self.canvas.after(PAN_DRAW_MS, self.pan_draw)

    def pan_draw(self):
        self.pan_job = None
        self.request_draw()

    def pan_end(self, e):
        self.pan_from = None
        if self.drag_hook is not None:
            self.drag_hook(False)

    # ---- drawing

    def set_layers(self, layers):
        self.layers = layers
        self.request_draw()

    def request_draw(self):
        if self.pending is None:
            self.pending = self.canvas.after_idle(self.draw)

    def draw(self, background=None):
        """Draw the view. `background`: a function; then the picture is made in a thread (numpy and cv2 let go of Python while they work, so the
        program goes on seeing to the mouse and the wheel meanwhile), and the function is called when it has been put on the canvas, or has been
        dropped because something else was drawn since. It is called at once if there is nothing to draw."""
        self.draw_inner(background)

    def draw_inner(self, background=None):
        self.pending = None
        if self.pan_job is not None:
            self.canvas.after_cancel(self.pan_job)
            self.pan_job = None
        if self.zoom_job is not None:
            self.canvas.after_cancel(self.zoom_job)
            self.zoom_job = None
        w, h = self.size()
        if w <= 2 or h <= 2:
            return background and background()
        if self.vsize is not None and self.vview is None:
            self.stale = True                                                 # (a virtual canvas no part of which is shown yet: place_tiles gives it its view, and that draws)
            return background and background()
        if self.visible_fn is not None and not self.visible_fn():
            self.stale = True
            return background and background()
        self.stale = False
        self.draw_token = getattr(self, "draw_token", 0) + 1
        if self.at_fit:
            self.scale = self.fit_scale()
            self.cx, self.cy = (self.home[0] + self.home[1]) / 2, (self.home[2] + self.home[3]) / 2
        zoom = self.scale / self.fit_scale()
        if self.vsize is not None:                                            # (a virtual canvas: only the part the widget shows is drawn)
            v = self.vview or (0, 0, 3, 3)
            ox, oy, ex, ey = v[0], v[1], v[0] + v[2], v[1] + v[3]
        else:
            ox, oy, ex, ey = 0, 0, w, h
        ox, oy, ex, ey = max(int(ox), 0), max(int(oy), 0), min(int(ex), w), min(int(ey), h)
        if ex - ox < 2 or ey - oy < 2:
            ox, oy, ex, ey = 0, 0, min(w, 3), min(h, 3)
        W, H = w, h                                                           # (the whole canvas; w, h below are the picture's)
        w, h = ex - ox, ey - oy
        cx, cy, scale = self.cx, self.cy, self.scale                          # (what the picture is made from, kept: the view may move while a thread makes it)
        # the data in view: the picture's pixels (ox..ex, oy..ey of the canvas) in data units
        x_lo, x_hi = cx + (ox - W / 2) / scale, cx + (ex - W / 2) / scale
        d_lo, d_hi = (oy - H / 2) / scale, (ey - H / 2) / scale
        y_lo, y_hi = (cy - d_hi, cy - d_lo) if self.flip_y else (cy + d_lo, cy + d_hi)
        opacity = 1.0 if self.moving else self.opacity                        # (solid while the view moves: blending translucent dots is much slower; see note_motion)
        diameters = []
        for lay in self.layers:
            d = lay["drawn_d"] = self.dot_diameter(lay, zoom)       # (kept, with the scale: the highlighted cell's outline is made as big as the dots really were drawn)
            lay["drawn_scale"] = scale
            diameters.append(d)
        layers, flip_y, bgcolor, cell_alpha, outlines = list(self.layers), self.flip_y, self.bgcolor, self.cell_alpha, self.outlines

        def render():
            img = np.empty((h, w, 3), np.uint8)
            img[:] = bgcolor
            for lay, d in zip(layers, diameters):
                x, y = lay["x"], lay["y"]
                pm = d / scale
                m = (x >= x_lo - pm) & (x <= x_hi + pm) & (y >= y_lo - pm) & (y <= y_hi + pm)
                idx = np.flatnonzero(m)
                if not len(idx):
                    continue
                px = np.floor((x[idx] - cx) * scale + W / 2 - ox).astype(np.int32)      # (floor: the same pixel as on the whole canvas, also left of or above the picture)
                dy = (y[idx] - cy) * scale
                py = np.floor((-dy if flip_y else dy) + H / 2 - oy).astype(np.int32)
                rgb = lay["rgb"][idx]
                if lay.get("shape") == "plus":
                    stamp_plus(img, px, py, rgb, d)
                    continue
                if lay.get("marks") is not None and not self.moving:         # (a marker code for each cell: the group's shape; plain discs while the view moves)
                    self.draw_marked(img, px, py, rgb, lay["marks"][idx], d)
                    continue
                if opacity < 1.0 and lay.get("shape") != "half":
                    blend_dots(img, px + d, py + d, rgb, d, opacity, self.kernels)
                    continue
                if d > BIG_DOT_PX:                                             # (very big dots, so few on screen: each drawn as a circle, later ones on top)
                    r = d // 2
                    for X, Y, c in zip(px.tolist(), py.tolist(), rgb.tolist()):
                        cv2.circle(img, (X, Y), r, c, -1)
                    continue
                px, py = px + d, py + d
                # each point writes its number (later = on top) at its pixel; a bigger dot is that number grown by a maximum filter
                lab = np.zeros((h + 2 * d, w + 2 * d), np.float32)
                lab[py, px] = np.arange(1, len(idx) + 1, dtype=np.float32)
                if lay.get("shape") == "half" and d > 1:                     # (the right half of a disc, flat side on the cell)
                    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (d, d))
                    k[:, :d // 2] = 0
                    lab = cv2.dilate(lab, k[:, d // 2:], anchor=(0, d // 2))
                elif d > 1:
                    lab = cv2.dilate(lab, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (d, d)))
                lab = lab[d:d + h, d:d + w]
                hit = lab > 0
                img[hit] = rgb[lab[hit].astype(np.int64) - 1]
            if cell_alpha < 0.999:
                img = (img.astype(np.float32) * cell_alpha + np.asarray(bgcolor, np.float32) * (1.0 - cell_alpha)).astype(np.uint8)
            if outlines:
                self.draw_outlines(img, W, H, ox, oy, (x_lo, x_hi, y_lo, y_hi), (cx, cy, scale), outlines)
            return img

        def put(img):
            self.drawn = (ox, oy, ex, ey)
            self.base = img
            self.compose()
            self.reposition()
            self.shown_view, self.shown_scale, self.previewed = (cx, cy), scale, (0, 0)          # (what is on the canvas is for this view: see preview_pan)
            self.remember_view()
            if self.on_view:
                self.on_view()

        if background is None:
            self.drawn = (ox, oy, ex, ey)           # (set first: the notes below say how big the picture was)
            put(render())
            return
        token, box = self.draw_token, {}

        def work():
            try:
                box["img"] = render()
            except Exception as e:                    # (e.g. the arrays changed under it: the picture is made again on the main thread when it is wanted)
                box["error"] = e

        thread = threading.Thread(target=work, daemon=True)
        thread.start()

        def poll():
            if thread.is_alive():
                self.canvas.after(4, poll)
                return
            if token == getattr(self, "draw_token", 0):          # (nothing has been drawn since: else this picture is for a view that is gone)
                if "img" in box:
                    put(box["img"])
                else:
                    self.stale = True
                    self.request_draw()
            background()
        self.canvas.after(4, poll)

    def set_outlines(self, polygons):
        self.outlines = polygons or None
        self.request_draw()

    def draw_outlines(self, img, w, h, ox, oy, view, at=None, outlines=None):
        """The outlines in view, one pixel wide, dark on a light background and light on a dark one. (w, h: the whole canvas; ox, oy: where the
        picture starts in it; view: the data x0, x1, y0, y1 it shows.)"""
        x0, x1, y0, y1 = view
        cx, cy, scale = at or (self.cx, self.cy, self.scale)
        shown = []
        for p in (self.outlines if outlines is None else outlines):
            lo, hi = p.min(axis=0), p.max(axis=0)
            if hi[0] < x0 or lo[0] > x1 or hi[1] < y0 or lo[1] > y1:
                continue
            px = (p[:, 0] - cx) * scale + w / 2 - ox
            dy = (p[:, 1] - cy) * scale
            py = (-dy if self.flip_y else dy) + h / 2 - oy
            shown.append(np.round(np.column_stack([px, py])).astype(np.int32))
        if shown:
            ink = (215, 215, 215) if sum(self.bgcolor) < 300 else (90, 90, 90)
            cv2.polylines(img, shown, True, ink, self.outline_px, cv2.LINE_AA)

    def marker_kernel(self, code, d):
        """The kernel with which a maximum filter makes a marker of this code and diameter d out of a point (cached)."""
        key = (code, d)
        k = self.shape_kernels.get(key)
        if k is None:
            k = self.shape_kernels[key] = (cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (d, d)) if code == 0
                                           else np.ascontiguousarray(shape_mask(code, d)[::-1, ::-1]))
        return k

    def draw_marked(self, img, px, py, rgb, marks, d):
        """Dots of different shapes: `marks` is a marker code for each point. Later points are on top of earlier ones, whatever their shape (each
        point writes its number at its pixel; each shape grows those numbers with its own kernel; the biggest number wins where they meet)."""
        h, w = img.shape[:2]
        keep = (px >= -d) & (px < w + d) & (py >= -d) & (py < h + d)       # (the caller's margin is d, inclusive, and rounding can put a point
        if not keep.all():                                               # exactly on its far edge: outside the margin of d the arrays below have)
            px, py, rgb, marks = px[keep], py[keep], rgb[keep], marks[keep]
        if d > BIG_DOT_PX:                                               # (very big dots, so few on screen: one by one)
            for X, Y, c, code in zip(px.tolist(), py.tolist(), rgb.tolist(), marks.tolist()):
                draw_big_marker(img, X, Y, c, d, int(code))
            return
        codes = np.unique(marks)
        opacity = 1.0 if self.moving else self.opacity                        # (solid while the view moves: see note_motion)
        if opacity < 1.0:
            for code in codes:
                sel = marks == code
                blend_dots(img, px[sel] + d, py[sel] + d, rgb[sel], d, opacity, self.kernels, kernel=self.marker_kernel(int(code), d))
            return
        if d >= 6 and len(px) <= h * w * d / 6000:                         # (few, big markers: each stamped where it is, in order, which comes out the same as the
            stamps = {int(c): self.marker_kernel(int(c), d)[::-1, ::-1] > 0 for c in codes}       # filter's result but is not a filter over the whole picture)
            off = d - 1 - d // 2
            for X, Y, c, code in zip(px.tolist(), py.tolist(), rgb.tolist(), marks.tolist()):
                x0, y0 = X - off, Y - off
                sx0, sy0, sx1, sy1 = max(-x0, 0), max(-y0, 0), min(d, w - x0), min(d, h - y0)
                if sx1 > sx0 and sy1 > sy0:
                    img[y0 + sy0:y0 + sy1, x0 + sx0:x0 + sx1][stamps[code][sy0:sy1, sx0:sx1]] = c
            return
        px, py = px + d, py + d
        numbers = np.arange(1, len(px) + 1, dtype=np.float32)
        best = None
        for code in codes:
            sel = marks == code
            lab = np.zeros((h + 2 * d, w + 2 * d), np.float32)
            lab[py[sel], px[sel]] = numbers[sel]
            if d > 1:
                lab = cv2.dilate(lab, self.marker_kernel(int(code), d))
            best = lab if best is None else np.maximum(best, lab)
        lab = best[d:d + h, d:d + w]
        hit = lab > 0
        img[hit] = rgb[lab[hit].astype(np.int64) - 1]

    def compose(self):
        """The picture shown: the points, veiled and with the highlighted group's '+' marks when those are on."""
        self.compose_pending = None
        if self.base is None:
            return
        if self.visible_fn is not None and not self.visible_fn():
            self.stale = True
            return
        img = self.base
        if self.dim:
            img = (img >> 1) + (np.array(self.bgcolor, np.uint8) >> 1)
        elif self.pluses is not None:
            img = img.copy()
        if self.pluses is not None and len(self.pluses[0]):
            h, w = img.shape[:2]
            x, y = self.pluses
            ox, oy = self.drawn[0], self.drawn[1]
            W, H = self.size()
            px = np.floor((x - self.cx) * self.scale + W / 2 - ox).astype(np.int32)
            dy = (y - self.cy) * self.scale
            py = np.floor((-dy if self.flip_y else dy) + H / 2 - oy).astype(np.int32)
            # The '+' grows with the zoom as the cells' dots do (the whole view's and the whole panel's zoom together: see dot_diameter), and gets thicker
            # with it; one at the whole view is as it always was.
            total = self.scale / self.fit_scale() * self.panel_zoom
            arm = max(int(round((self.plus_arm + 1) * self.marker_factor * total ** 0.7)), 1)
            keep = (px >= -arm) & (px < w + arm) & (py >= -arm) & (py < h + arm)              # (only those that can show)
            stamp_plus(img, px[keep], py[keep], np.array([255, 0, 0], np.uint8), arm, thick=max(2, round(arm / 6)))
        self.photo = ImageTk.PhotoImage(Image.fromarray(np.ascontiguousarray(img)), master=self.canvas)
        self.canvas.itemconfig(self.item, image=self.photo)
        self.canvas.coords(self.item, self.drawn[0] - self.off[0], self.drawn[1] - self.off[1])

    def request_compose(self):
        if self.pending is None and self.compose_pending is None:
            self.compose_pending = self.canvas.after_idle(self.compose)

    def set_background(self, rgb):
        if tuple(rgb) != tuple(self.bgcolor):
            self.bgcolor = tuple(rgb)
            self.canvas.config(bg="#%02x%02x%02x" % self.bgcolor)
            if self.scalebar:                                   # (white on a dark background)
                ink = "white" if is_dark(self.bgcolor) else "black"
                self.canvas.itemconfig(self.scalebar[0], fill=ink)
                self.canvas.itemconfig(self.scalebar[1], fill=ink)
            for item, color in self.label_color.items():
                self.canvas.itemconfig(item, fill=ink_for(color, self.bgcolor))

    def set_dim(self, on):
        if on != self.dim:
            self.dim = on
            self.request_compose()

    def set_pluses(self, xy):
        if xy is None and self.pluses is None:
            return
        self.pluses = xy
        self.request_compose()

    def refresh(self):
        """Called when this canvas may have come on screen: draw it if it was waiting."""
        if self.stale and self.pending is None:
            self.request_draw()

    # ---- items in data coordinates

    def dot_diameter(self, lay, zoom):
        """The diameter in pixels of a layer's dots (`zoom`: the view's, scale / fit_scale)."""
        if "size_mm" in lay:                                              # (a size on the tissue: the same at any window size)
            return int(max(round(lay["size_mm"] * self.dot_factor * self.scale), 1))
        return int(max(round(lay["size"] * (zoom * self.panel_zoom) ** 0.7), 1))          # (the dots grow as the view, or the whole panel, is zoomed in)

    def shaped_cell(self, x, y):
        """(layer, index in it) of the cell drawn at data position (x, y) with a group's shape (the topmost layer that has one there: the dots of
        the other kind of cell are a layer of their own, smaller), or None: no layer has shapes (cells not drawn by group), or none has a cell there."""
        key = (id(self.layers), x, y)
        if self.shaped_key != key:
            found = None
            tol = 1e-6 * max(self.home[1] - self.home[0], self.home[3] - self.home[2], 1e-12)
            for lay in reversed(self.layers):
                if lay.get("marks") is None or not len(lay["x"]):
                    continue
                dist = np.hypot(lay["x"] - x, lay["y"] - y)
                k = int(np.argmin(dist))
                if dist[k] <= tol:
                    found = (lay, k)
                    break
            self.shaped_key, self.shaped_found = key, found
        return self.shaped_found

    def add_marker(self, name, radius=8):
        item = self.canvas.create_oval(0, 0, 0, 0, outline=RING, width=2, state="hidden")
        # [ring, x, y, radius, the cell redrawn over itself with a darker edge (when cells are drawn by group)]
        shape = self.canvas.create_polygon(0, 0, 0, 0, 0, 0, fill="", outline="black", state="hidden")
        self.markers[name] = [item, None, None, radius, [shape], {"job": None, "off": False}]       # (the last: its blinking, see set_marker)

    MARKER_STYLES = {"ring": ("", "#e00000"), "mirror": ("white", "black")}      # (fill, outline): the cell under the pointer; the same cell in the other plot

    def set_marker(self, name, x, y, style="ring"):
        """When cells are drawn by group the cell is also drawn again over itself, with its shape, size and color and a darker edge (a ring
        does not make a triangle's place clear: its center is not the ring's); see place_marker."""
        m = self.markers[name]
        self.stop_blink(m)
        m[1], m[2] = x, y
        fill, outline = self.MARKER_STYLES[style]
        self.canvas.itemconfig(m[0], fill=fill, outline=outline)
        self.place_marker(m)
        if x is not None:                         # (after a pause it blinks: a single cell is easy to lose among so many)
            m[5]["job"] = self.canvas.after(BLINK_PAUSE_MS, lambda: self.blink_step(m, True))

    def stop_blink(self, m):
        job, m[5]["job"], m[5]["off"] = m[5]["job"], None, False
        if job is not None:
            try:
                self.canvas.after_cancel(job)
            except tk.TclError:
                pass

    def blink_step(self, m, off):
        """Hide the marked cell (`off`) or show it again, and come back for the other half; ends when the marker is cleared or moved."""
        m[5]["job"] = None
        if m[1] is None:
            return
        try:
            m[5]["off"] = off
            self.place_marker(m)
            m[5]["job"] = self.canvas.after(BLINK_OFF_MS if off else BLINK_ON_MS, lambda: self.blink_step(m, not off))
        except tk.TclError:
            m[5]["off"] = False                   # (the window is gone)

    def place_marker(self, m):
        if m[1] is None:
            for item in [m[0]] + m[4]:
                self.canvas.itemconfig(item, state="hidden")
            return
        px, py = self.to_canvas(m[1], m[2])
        r = m[3] * self.marker_factor
        self.canvas.coords(m[0], px - r, py - r, px + r, py + r)
        self.canvas.itemconfig(m[0], state="normal")
        self.canvas.tag_raise(m[0])
        got = self.shaped_cell(m[1], m[2])
        for item in m[4]:                         # (over the ring: the other plot's cell has a white disc for a ring)
            if got is None:
                self.canvas.itemconfig(item, state="hidden")
                continue
            lay, k = got
            code, rgb = int(lay["marks"][k]), lay["rgb"][k]
            # The size it was drawn at. Until the next drawing (after a zoom step the picture is stretched at once, and drawn again a moment later)
            # the picture on screen is that one stretched by the change in scale, so this is too.
            d = lay["drawn_d"] * self.scale / lay["drawn_scale"] if lay.get("drawn_d") else self.dot_diameter(lay, self.scale / self.fit_scale())
            # The cell again, as it was drawn (its shape, size and color) with a darker edge of its color, both in proportion to its size, so they
            # keep their relation at any zoom. A canvas line is centered on its path.
            hollow = code >= 10
            if hollow:
                # The cell's own line is d / 6 wide, inside its edge (up to BIG_DOT_PX; the maximum filter's), or d / 12, centered on its edge
                # (beyond: cv2's, see draw_big_marker). The new line follows the same path and is a third wider, so it covers it and shows as
                # one darker line of the same color family.
                line = max(1, round(d / 6)) if d <= BIG_DOT_PX else max(1, d // 12) + 2       # (cv2 draws a line about 2 px wider than asked)
                wd = max(2.0, 1.35 * line)
                size = d - line if d <= BIG_DOT_PX else d
            else:
                wd = max(2.0, d / 10)                      # (inside the cell's edge: half the line's width in, the whole stays the cell's size)
                size = d - wd
            fill = "#%02x%02x%02x" % tuple(int(c) for c in rgb[:3])
            edge = "#%02x%02x%02x" % tuple(int(c * 0.55) for c in rgb[:3])
            self.canvas.coords(item, *shape_outline(code, px, py, max(size, 2)))
            self.canvas.itemconfig(item, state="normal", fill="" if hollow else fill, outline=edge, width=wd,
                                   joinstyle=tk.MITER)
            self.canvas.tag_raise(item)
        if m[5]["off"]:                           # (the dark half of a blink: placed all the same, so that it is where it should be on coming back)
            for item in [m[0]] + m[4]:
                self.canvas.itemconfig(item, state="hidden")

    def add_rect(self, x0, y0, x1, y1):
        item = self.canvas.create_rectangle(0, 0, 0, 0, outline="black", width=2, dash=dash((6, 3), 2))
        self.rects.append([item, x0, y0, x1, y1])

    def rect_width(self):
        """The ROI outline's line width: set by the Section size slider (see ROI_LINE_PX), thicker as the sections are bigger."""
        return self.rect_px

    def enable_scalebar(self):
        self.scalebar = (self.canvas.create_line(0, 0, 0, 0, fill="black", width=3), self.canvas.create_text(0, 0, anchor="s", font=("Segoe UI", 9, "bold"), fill="black"))

    def place_scalebar(self):
        w, h = self.size()
        view_mm = w / self.scale
        cands = [u for u in SCALEBAR_UM if u <= view_mm * 1000 * 0.22]
        um = cands[-1] if cands else SCALEBAR_UM[0]
        vw, vh = self.view_size()
        x0, y = vw * 0.06, vh * 0.94
        x1 = x0 + um / 1000 * self.scale
        self.canvas.coords(self.scalebar[0], x0, y, x1, y)
        self.canvas.coords(self.scalebar[1], (x0 + x1) / 2, y - 4)
        self.canvas.itemconfig(self.scalebar[1], text=f"{um:g} µm" if um < 1000 else f"{um / 1000:g} mm")
        self.canvas.tag_raise(self.scalebar[0])
        self.canvas.tag_raise(self.scalebar[1])

    def clear_labels(self):
        for it in self.labels:
            self.canvas.delete(it[0])
        self.labels = []
        self.label_text = {}
        self.label_weight = {}
        self.label_color = {}
        self.boosted = None

    def add_label(self, text, x, y, color="black", weight=1.0):
        """`weight`: times the font size (the labels are sized by their group's cell count)."""
        item = self.canvas.create_text(0, 0, text=text, font=("Segoe UI", self.label_size, "bold"), fill=ink_for(color, self.bgcolor) if color.startswith("#") else color)
        self.label_color[item] = color if color.startswith("#") else "#000000"
        self.labels.append([item, x, y])
        self.label_text[item] = text
        self.label_weight[item] = weight
        # its own size now (the view's size has not changed, so update_label_size would not do it)
        self.canvas.itemconfig(item, font=("Segoe UI", self.label_font_size(item, self.label_size or 8), "bold"))

    def boost_label(self, text):
        """Draw the label with this text bigger (the group being marked with '+' signs), and the others as they were; None puts it back."""
        if text == self.boosted:
            return
        self.boosted = text
        self.update_label_size(force=True)
        self.reposition()

    def set_label_factor(self, factor):
        self.label_factor = factor
        self.label_size = None                      # (so the next update applies the new size)
        self.update_label_size()

    def label_font_size(self, item, size):
        if self.scale_labels:
            size = max(int(round(size * self.label_weight.get(item, 1.0))), 1)
        big = self.boosted is not None and self.label_text.get(item) == self.boosted
        return max(int(round(size * LABEL_HIGHLIGHT_FACTOR)), size + 1) if big else size

    def update_label_size(self, force=False):
        zoom = min(self.scale / self.fit_scale(), LABEL_ZOOM_CAP)
        size = max(int(round(8 * self.label_factor * (1 + 0.75 * (zoom - 1)))), 1)
        if size != self.label_size or force:
            self.label_size = size
            for it in self.labels:
                self.canvas.itemconfig(it[0], font=("Segoe UI", self.label_font_size(it[0], size), "bold"))

    def reposition(self, light=False):
        """Put the ring, the rectangles, the scale bar and the group labels where they belong in the view. `light`: all but the labels, which are hidden
        (hundreds of them: moving them and sizing their text at every wheel step is most of what a step cost)."""
        if light:
            if not self.labels_hidden:
                for it, _x, _y in self.labels:
                    self.canvas.itemconfig(it, state="hidden")
                self.labels_hidden = True
        elif self.labels:
            self.update_label_size()
            self.labels_hidden = False
        if self.scalebar:
            self.place_scalebar()
        for m in self.markers.values():
            self.place_marker(m)
        wd = self.rect_width() if self.rects else 2
        for it, x0, y0, x1, y1 in self.rects:
            a, b = self.to_canvas(x0, y0)
            c, d = self.to_canvas(x1, y1)
            self.canvas.coords(it, a, b, c, d)
            if self.rect_line != wd:
                self.canvas.itemconfig(it, width=wd, dash=roi_dash(wd))
            self.canvas.tag_raise(it)
        self.rect_line = wd
        w, h = self.view_size()
        for it, x, y in ([] if light else self.labels):
            px, py = self.to_canvas(x, y)
            if 0 <= px <= w and 0 <= py <= h:
                self.canvas.coords(it, px, py)
                self.canvas.itemconfig(it, state="normal")
                self.canvas.tag_raise(it)
            else:
                self.canvas.itemconfig(it, state="hidden")
        if self.boosted is not None:                  # (the bigger number on top of the others)
            for it, _x, _y in self.labels:
                if self.label_text.get(it) == self.boosted:
                    self.canvas.tag_raise(it)

    def nearest(self, tree, px, py):
        """Index (into the tree's points) of the point within reach of the canvas position, or None. The reach is HOVER_RADIUS_PX, or more when the
        dots have grown with the zoom (cells are then far apart on the screen): most of a dot's width, so the pointer need not be on the center."""
        x, y = self.to_data(px, py)
        biggest = max((lay.get("drawn_d", 0) for lay in self.layers), default=0)
        reach = max(HOVER_RADIUS_PX, HOVER_REACH_DOTS * biggest)
        dist, k = tree.query([x, y], distance_upper_bound=reach / self.scale)
        return None if not np.isfinite(dist) else int(k)


# ---------------------------------------------------------------------------------------------------------------------- the window

class UmapWindow(ExportsMixin):
    def __init__(self, adata, abc_cache, imputed_state, adata_backed, rois, run_folder, ensure_imputed, title, parent=None, host=None, params=None,
                 own=None, panel_zoom=False):
        self.start_panel_zoom = panel_zoom                  # start with Zoom whole panel ticked (and the section maps at their smallest)
        self.adata, self.abc_cache = adata, abc_cache
        self.own = own                      # a dataset of the user's own (abc_tk_own.OwnDataset): no atlas, no sections, no imputed genes
        self.no_umap = bool(adata.uns.get("no_umap"))        # a project without a UMAP (abc_project.py): the section maps and the genes, read from the atlas
        self.atlas_counts = None            # ... whose genes are read from the atlas's own file (abc_project.AtlasCounts)
        self.level_names = dict(own.level_names) if own is not None else dict(LEVEL_NAMES)          # the group types: column -> the words on the screen
        self.level_has_ids = dict(own.level_has_ids) if own is not None else {"neurotransmitter": False}      # (whether a level's groups have numbers to type)
        self.host = host or {}             # the browser's own functions this window uses (settings dialog, UMAP variants, the exports)
        self.run_folder = run_folder
        self.rois_input = list(rois or [])                   # (as given: a subset's window draws the same ROIs)
        self.subset_windows = []                             # the windows Analyze subset opened
        self.params = dict(params) if params else None
        self.viewer_settings = self.host["load_viewer_settings"]() if "load_viewer_settings" in self.host else {"hover_enabled": True, "gene456_shape": "small circle", "umap_circle_size": "Default", "persist_settings": False}
        self.base_title = title
        self.imputed_state = imputed_state if imputed_state is not None else {"adata": None, "load_thread": None, "load_error": None}
        self.ensure_imputed = ensure_imputed
        self.own_root = parent is None
        self.root = tk.Tk() if parent is None else tk.Toplevel(parent)
        self.root.title(title)
        self.root.geometry("1750x950")
        try:
            self.root.state("zoomed")                          # (its final size from the start: nothing is drawn at one size and then again at another)
        except tk.TclError:
            pass
        loading =ttk.Label(self.root, text="Preparing the project window...", font=("Segoe UI", 14), padding=40)
        loading.place(relx=0.5, rely=0.5, anchor="center")     # (in the middle of the window, which already fills the screen)
        self.root.update()
        self.msgs = queue.Queue()
        self.busy = False
        self.gc_was_enabled = gc.isenabled()
        gc.disable()
        self.mode = "all"                  # "all" | "single" | "gene" | "imputed"
        if own is not None:
            names = list(own.level_names)
            self.level = next((k for k in names if k.lower().endswith("subclass")), next((k for k in names if k.lower().endswith("class")), names[0]))
        else:
            self.level = "subclass" if "subclass" in adata.obs.columns else "class"
        self.shown_gene = None
        self.number_arrays = {}            # a dataset of your own: the numeric columns' values, made when first hovered
        self.shown_dataset = "standard"
        self.values = None                 # the shown gene's value for each cell (as coloured), or None
        self.tree = None
        self.hover_job = self.group_job = None
        self.gene_cache = None             # ((genes, imputed, scale), values read) of the genes shown, so a legend click re-uses them
        self.hidden = set()                # the genes switched off by a legend click
        self.multi_rows = None
        self.last_group_click = None       # (group, whether it was switched to gray) of the last click on a legend row: where a shift-click starts
        self.hover_group = None            # the legend row the pointer is on
        self.single_names = {}             # Single group(s): number -> full name of each group
        self.off_groups = set()            # the groups switched to gray by a click on their legend swatch (All groups: their names; Single group(s): their numbers)
        self.last_ids = None               # the groups typed in Single group(s): level, targets, colors, cells
        self.export_cat = None             # what the categorical colours mean (for a saved legend)
        self.label_info = []
        self.view_kind, self.view_token, self.view_title = "all", "all", ""
        self.hover_key = None
        self.ccf_started = False           # start_ccf_lookup has been called (by the first cell tip)
        self.ccf_run = self.ccf_bg = None  # the cells' CCF coordinates (AP, DV, ML), for the ML and DV of the cell tip: see start_ccf_lookup
        self.ccf_state = "none"
        self.variants = {}                 # UMAP parameter key -> (coordinates, Leiden labels), the ones already shown
        try:
            self.load_data(adata_backed, rois)
            if self.host and self.params is not None:
                self.variants[self.host["umap_params_key"](self.params)] = (self.coords.copy(), self.info.get("leiden"))
            loading.destroy()
            self.build()
            self.recolor()
        except BaseException:
            try:
                self.root.destroy()                       # (not left on the screen, blank, when the window could not be made)
            except tk.TclError:
                pass
            raise
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.bind("<Escape>", lambda e: self.close())
        self.poll()
        # (No background loads here: they held Python's one running thread for seconds at a time, just when the section maps were being zoomed. The imputed
        # dataset is loaded when Imputed Gene is chosen, the CCF coordinates when the first cell tip needs them: see coord_rows.)

    def preload_imputed(self):
        """Start loading the imputed gene dataset in the background, if it is already on this computer: the choice "Imputed Gene" is common."""
        start = self.host.get("preload_imputed")
        if start is not None and self.own is None:
            try:
                start(self.imputed_state, self.abc_cache)
            except Exception as e:
                print(f"Could not start loading the imputed genes in the background: {e}", flush=True)

    # ---- data

    def load_data(self, adata_backed, rois):
        a = self.adata
        if self.no_umap:
            self.coords = np.zeros((a.n_obs, 2))                # (no embedding: every cell at one point of a UMAP that is not shown)
        else:
            self.coords = np.asarray(a.obsm["X_umap"], dtype=float)[:, :2]
        self.n = len(self.coords)
        if self.no_umap and adata_backed is not None:
            self.atlas_counts = abc_project.AtlasCounts(adata_backed, core.atlas_positions(a.obs.index, adata_backed.obs.index))
        self.has_counts = "counts" in a.layers or self.own is not None or self.atlas_counts is not None
        self.has_umap_data = not self.no_umap                                   # (a project made without a UMAP has none to show)
        sec_col = core.SECTION_COL
        bsp = None
        if adata_backed is not None and sec_col in a.obs.columns:
            bsp = core.load_section_spatial_coords(adata_backed, self.abc_cache)      # (the whole atlas: the run's cells are taken from it, so the table is read once)
        if sec_col not in a.obs.columns:
            spatial = None
        elif bsp is not None:
            pos = core.atlas_positions(a.obs.index, adata_backed.obs.index)          # (the run's cells in the atlas's table: by position, no lookup of the cell IDs)
            spatial = bsp.iloc[pos].set_axis(a.obs.index) if (pos >= 0).all() else bsp.reindex(a.obs.index)
        else:
            spatial = core.load_section_spatial_coords(a, self.abc_cache)
        self.section_of = a.obs[sec_col].to_numpy() if sec_col in a.obs.columns else None
        ok = spatial is not None and "x" in spatial.columns and "y" in spatial.columns
        self.cx = pd.to_numeric(spatial["x"], errors="coerce").to_numpy() if ok else None
        self.cy = pd.to_numeric(spatial["y"], errors="coerce").to_numpy() if ok else None
        self.region_of_run = spatial["parcellation_index"].to_numpy() if ok and "parcellation_index" in spatial.columns else None      # (Allen's CCF region of each cell)
        self.region_names = None
        if self.region_of_run is not None:
            try:
                self.region_names = core.load_region_names(self.abc_cache)
            except Exception as e:
                print(f"Could not read the names of the CCF regions: {e}", flush=True)
        self.sections = core.sorted_sections_descending(a.obs[sec_col].dropna().unique().tolist()) if ok else []
        self.info = {k: (a.obs[k].astype(str).to_numpy() if k in a.obs.columns else None) for k in self.level_names}
        if self.own is None and self.info.get("neurotransmitter") is None:           # (runs made before it was a column: read from the atlas's table, as the taxonomy was)
            try:
                nt = core.load_annotations_aligned(a.obs_names, self.abc_cache)
                if "neurotransmitter" in nt.columns:
                    self.info["neurotransmitter"] = nt["neurotransmitter"].astype(str).to_numpy()
            except Exception as e:
                print(f"Could not read the neurotransmitter of the cells: {e}", flush=True)
        self.levels = [k for k in self.level_names if self.info[k] is not None]
        # the section panels' background: every cell of each section (from the unfiltered data), when given
        self.bg = None
        self.adata_backed = adata_backed
        if ok and bsp is not None:
            try:
                bsec = core.get_section_labels(adata_backed, self.abc_cache)
                if "x" in bsp.columns:
                    codes, uniques = pd.factorize(bsec.to_numpy())
                    order = np.argsort(codes, kind="stable")
                    starts = np.searchsorted(codes[order], np.arange(len(uniques) + 1))
                    where = {u: order[starts[i]:starts[i + 1]] for i, u in enumerate(uniques)}
                    self.bg = {"ids": adata_backed.obs_names, "sec_arr": bsec.to_numpy(), "z": (pd.to_numeric(bsp["z"], errors="coerce").to_numpy() if "z" in bsp.columns else None),
                               "x": pd.to_numeric(bsp["x"], errors="coerce").to_numpy(), "y": pd.to_numeric(bsp["y"], errors="coerce").to_numpy(),
                               "region": bsp["parcellation_index"].to_numpy() if "parcellation_index" in bsp.columns else None,
                               "where": where, "levels": LazyLevels(bsp, self.level_names)}
            except Exception as e:
                print(f"Could not load the all-cell background for the section panels: {e}", flush=True)
        self.run_idx = {}
        if ok:
            codes, uniques = pd.factorize(self.section_of)
            for i, u in enumerate(uniques):
                if pd.notna(u):
                    self.run_idx[u] = np.flatnonzero(codes == i)
        self.id_of_cell = {}
        self.rois = []
        for r in rois or []:
            try:
                vals = [float(r[k]) for k in ("x_min", "x_max", "y_min", "y_max")]
            except (TypeError, ValueError, KeyError):
                continue
            if all(np.isfinite(vals)):
                self.rois.append((str(r["section"]), *vals))
        self.cache = {}
        self.run_index = pd.Index(a.obs_names)
        self.bg_under = None
        if self.info.get("class") is not None:
            cid = core.extract_leading_numeric_id(pd.Series(self.info["class"])).to_numpy()
            frac = float(np.isin(cid, list(core.NON_NEURON_CLASS_IDS)).mean())
            self.bg_under = "nonneuron" if frac < 0.5 else "neuron"           # (the minority kind is drawn underneath)
        # Show UMAP: off to begin with for a large project, whose section maps then take the room. Not when there are no section maps to take it (a
        # dataset of the user's own, a project without sections): then the UMAP is all there is.
        self.umap_shown = self.has_umap_data and (self.n < LARGE_PROJECT_CELLS or self.own is not None or not self.sections)
        # ...unless the run says which cell type it was made for: the launch screen's choice is in the run folder's name
        name = Path(self.run_folder).name.lower() if self.run_folder is not None else ""
        for prefix, under in (("umap_neurons_", "nonneuron"), ("umap_nonneurons_", "neuron"), ("umap_all-cells_", None), ("umap_all_", None)):         # (umap_all_: the older name of an all-cells run)
            if name.startswith(prefix):
                self.bg_under = under

    def category_codes(self, level):
        """(codes of each run cell, uniques) of a level; and for the all-cell background (codes, uniques) over the same uniques."""
        key = ("codes", level)
        if key not in self.cache:
            codes, uniques = pd.factorize(pd.Series(self.info[level]).replace("nan", np.nan))      # (a cell with no value, such as a non-neuron's neurotransmitter: code -1)
            self.cache[key] = (codes, uniques)
        return self.cache[key]

    def bg_codes(self, level):
        key = ("bgcodes", level)
        if key not in self.cache:
            arr = self.bg["levels"].get(level) if self.bg else None
            if arr is None:
                self.cache[key] = None
            else:
                codes, uniques = self.category_codes(level)
                self.cache[key] = pd.Index(uniques).get_indexer(arr)       # (-1: a group none of the run's cells is in)
        return self.cache[key]

    # ---- the window

    def build(self):
        self.tiles = []                    # (made by build_tiles; empty until then, for the scroll callbacks that can come first)
        self.fit_cache, self.fit_across = {}, None      # (see fit_tiles)
        self.started_in_background = 0
        self.cell_tip = CellTip(self.root)
        self.cell_tip.row_hover_fn = self.tip_row_hover
        self.pw = pw = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True)
        # the section panels
        left = ttk.Frame(pw)
        bars = ttk.Frame(left, padding=(6, 4))
        bars.pack(side=tk.TOP, fill=tk.X)
        # Zoom whole panel from the start, in every project that has section maps: they start at the size it chooses (panel_zoom_tile_index), so that they
        # are made once, not made and then again
        self.zoom_from_start = bool(self.sections)
        self.var_tile = tk.IntVar(bars, self.panel_zoom_tile_index() if self.zoom_from_start else DEFAULT_TILE_INDEX)
        self.tile_stop_before_zoom = DEFAULT_TILE_INDEX if self.zoom_from_start else None
        self.var_dot = tk.IntVar(bars, DEFAULT_DOT_INDEX)
        self.var_marker = tk.IntVar(bars, DEFAULT_MARKER_INDEX)
        self.dot_factors = [um / (DOT_MM * 1000) for um in CELL_SIZES_UM]          # (times DOT_MM)
        marker_names = ("Small", "Medium", "Large")
        self.sliders = [
            self.add_slider(bars, 0, "Section size:", self.var_tile, len(TILE_COLUMNS), self.layout_tiles,
                            lambda i: "fit all" if TILE_COLUMNS[i] == 0 else "full width" if TILE_COLUMNS[i] == 1 else f"1/{TILE_COLUMNS[i]} width",
                            lambda i: (f"Every section in the pane at once, as large as that allows: {self.fit_across} across." if self.fit_across else
                                       "Every section in the pane at once, as large as that allows.") if TILE_COLUMNS[i] == 0 else f"{TILE_COLUMNS[i]} across the pane."),
            self.add_slider(bars, 1, "Cell size:", self.var_dot, len(CELL_SIZES_UM), self.apply_dot_size,
                            lambda i: f"{CELL_SIZES_UM[i]} \u00b5m", None,
                            label_tip="The size of every cell's dot in the section maps: all the neurons or glia shown, not only the highlighted ones."),
            self.add_slider(bars, 2, "Marker size:", self.var_marker, len(MARKER_SIZES), self.apply_marker_size, lambda i: marker_names[i], None,
                            label_tip="The size of the marks on the highlighted cells only (the ring under the pointer and the red + marks), "
                                      "not of the other cells' dots.")]
        # The buttons on three lines (narrow, as tall as the sliders). On macOS an ordinary button's bezel is drawn at the system's one height,
        # a line of text, whatever the button's own: the three lines were cut to one. A Toolbutton's bevel there takes the button's height.
        aqua = self.root.tk.call("tk", "windowingsystem") == "aqua"
        self.stacked_style = "Stacked.Toolbutton" if aqua else "Stacked.TButton"
        ttk.Style(self.root).configure(self.stacked_style, justify="center", padding=(10, 4), anchor="center")       # (room round the text)
        ttk.Button(bars, text="Reset\nslider\ndefault", width=0, style=self.stacked_style,
                   command=self.reset_sliders).grid(row=0, column=0, rowspan=3, sticky="nsw")
        # the checkboxes, in two columns: Hover highlight, ROI only, Zoom whole panel; then CCF atlas overlay, with its slider under it
        self.var_roi_only = tk.BooleanVar(bars, False)
        roi_only_box = ttk.Checkbutton(bars, text="ROI only", variable=self.var_roi_only, command=self.apply_roi_only)
        roi_only_box.grid(row=1, column=4, sticky="w", padx=(6, 6))
        if self.rois:
            HoverTip(roi_only_box, lambda: "The cells outside the ROIs go gray, on a white background: only the ROIs are colored. Saved section maps "
                                           "have the same look (and, as an SVG, the gray cells are pictures, so the file stays small).")
        else:
            roi_only_box.state(["disabled"])
            HoverTip(roi_only_box, lambda: "This project has no ROIs.")
        self.var_show_umap = tk.BooleanVar(bars, self.umap_shown)
        show_umap_box = ttk.Checkbutton(bars, text="Show UMAP", variable=self.var_show_umap, command=self.toggle_umap)
        show_umap_box.grid(row=0, column=5, sticky="w", padx=(6, 6))
        HoverTip(show_umap_box, lambda: ("Show the UMAP beside the section maps. Off, the section maps use its room. A project of "
                                         f"{LARGE_PROJECT_CELLS:,} cells or more starts with it off." if self.has_umap_data
                                         else "This project was made without a UMAP."))
        if not self.has_umap_data:
            show_umap_box.state(["disabled"])                  # (a project made without a UMAP has none to show: greyed out, with the tip saying why)
        self.var_ccf = tk.BooleanVar(bars, False)
        ccf_box = ttk.Checkbutton(bars, text="CCF atlas overlay", variable=self.var_ccf, command=self.apply_ccf)
        ccf_box.grid(row=1, column=5, sticky="w", padx=(6, 6))
        HoverTip(ccf_box, lambda: "Draw the outlines of the regions of the Allen Mouse Brain Common Coordinate Framework (CCF, 2020 version) "
                                  "over the section maps. Placed by the cells' reconstructed coordinates, they are the atlas's own; otherwise "
                                  "fitted to each section from the CCF position of its cells (the first time reads the CCF coordinates of the "
                                  "cells: a minute or so).")
        ccf_row = ttk.Frame(bars)                                                              # (its slider, under it)
        ccf_row.grid(row=2, column=5, sticky="w", padx=(22, 6))
        self.var_ccf_opacity = tk.IntVar(bars, DEFAULT_OPACITY_INDEX)                          # (full: the cells as they are)
        ccf_shown = ttk.Label(ccf_row, width=5)
        self.ccf_scale = self.stop_slider(ccf_row, ccf_shown, self.var_ccf_opacity, len(OPACITIES), self.apply_ccf_opacity, lambda i: f"{OPACITIES[i] * 100:g}%",
                                          length=90)            # (shorter than the others, so that the values of the sliders on the left have room)
        self.ccf_scale.pack(side=tk.LEFT)
        ccf_shown.pack(side=tk.LEFT, padx=(4, 0))
        self.ccf_scale.state(["disabled"])                                                      # (until the atlas is shown)
        HoverTip(self.ccf_scale, lambda: "How opaque the cells are while the CCF atlas is shown: lower it to see the outlines better. "
                                         "The outlines stay fully opaque.")
        self.ccf_outlines = None
        self.panel_zoom = 1.0                                                                  # (the whole-panel zoom: see set_panel_zoom)
        self.zoom_burst = None                                                                 # (a run of wheel steps being shown on a snapshot: see zoom_step)
        self.panel_hook = types.SimpleNamespace(zoom=self.panel_zoom_at, pan_start=self.panel_pan_start, pan_move=self.panel_pan_move, reset=self.panel_zoom_reset)
        self.var_panel_zoom = tk.BooleanVar(bars, self.zoom_from_start)
        panel_zoom_box = ttk.Checkbutton(bars, text="Zoom whole panel", variable=self.var_panel_zoom, command=self.apply_panel_zoom_mode)
        panel_zoom_box.grid(row=2, column=4, sticky="w", padx=(6, 6))
        HoverTip(panel_zoom_box, lambda: "Off: the wheel zooms into each section map, which stays in its panel (clipped by its border). "
                                         "On: the wheel zooms the whole panel, the borders and the sections together, and a drag moves the pane; "
                                         "sections scrolled out of view are not drawn. Double-click puts the panel back.")
        # (the same setting as in the Settings window; but off to begin with for "Browse All Neurons", whose many small maps make it distracting,
        # and there the box is for this window alone: see apply_hover_enabled)
        self.var_hover = tk.BooleanVar(bars, False if self.start_panel_zoom else bool(self.viewer_settings.get("hover_enabled", True)))
        hover_box = ttk.Checkbutton(bars, text="Hover highlight", variable=self.var_hover, command=self.apply_hover_enabled)
        hover_box.grid(row=0, column=4, sticky="w", padx=(6, 6))
        self.hover_box = hover_box
        HoverTip(hover_box, lambda: ("When the pointer rests on a cell, dim the plots and mark every cell of its group with a red +. Off: no "
                                     "marking, but the cell is still ringed and described in its tip (the same setting as in Settings)."
                                     if self.umap_shown else "Hover highlight marks the cell's group in the UMAP as well as in the section maps: "
                                                             "it is available while the UMAP is shown."))
        if not self.umap_shown:
            hover_box.state(["disabled"])
        self.update_marker_slider()
        bars.columnconfigure(6, weight=1)          # (the spare width goes to an empty column on the right: the controls stay where they are)
        self.tile_scroll = tk.Canvas(left, highlightthickness=0, bg="#f4f4f4")
        self.watch_scroll(self.tile_scroll)
        sb = ttk.Scrollbar(left, orient=tk.VERTICAL, command=self.tile_scroll.yview)
        self.hsb = ttk.Scrollbar(left, orient=tk.HORIZONTAL, command=self.tile_scroll.xview)        # (shown only while the panel is zoomed in)
        self.tile_scroll.configure(yscrollcommand=lambda a, b: (sb.set(a, b), self.note_view("y", a), self.root.after_idle(self.refresh_visible)),
                                   xscrollcommand=lambda a, b: (self.hsb.set(a, b), self.note_view("x", a), self.root.after_idle(self.refresh_visible)))
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        sb.bind("<MouseWheel>", lambda e: self.tile_scroll.yview_scroll(-1 if e.delta > 0 else 1, "units"))
        self.tile_scroll.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # (Each section panel is its own item in the scrolled canvas, not a part of one big frame: a Tk window cannot be taller than about 32,000
        # pixels, so with 60 panels as tall as the pane is wide one frame holding them all showed only the first 30 or so, the last cut off.)
        self.layout_job, self.layout_seen = None, False
        self.tile_scroll.bind("<Configure>", lambda e: self.layout_tiles_when_settled())
        self.bind_pane_gestures(self.tile_scroll)
        if self.own is None:                                        # (a dataset of the user's own has no sections)
            pw.add(left, weight=2)
        # the UMAP
        mid = ttk.Frame(pw)
        umap_bar = ttk.Frame(mid, padding=(6, 4))
        umap_bar.pack(side=tk.TOP, fill=tk.X)
        # Cell size: the UMAP circle size of Settings (the same setting, so either changes the other). A dense UMAP starts at Small, unless
        # the size is remembered from an earlier session (Settings' "Remember these settings"): then that.
        self.umap_size_names = list(self.host.get("umap_circle_sizes", {}) or {"Default": 1.0})
        self.umap_default_size = "Small" if self.n > DENSE_UMAP_CELLS and "Small" in self.umap_size_names else "Default"
        if not self.viewer_settings.get("persist_settings", False):
            self.viewer_settings["umap_circle_size"] = self.umap_default_size
        self.var_umap_dot = tk.IntVar(umap_bar, self.umap_size_index(self.viewer_settings.get("umap_circle_size", "Default")))
        self.var_label = tk.IntVar(umap_bar, DEFAULT_LABEL_INDEX)
        self.umap_sliders = [
            self.add_slider(umap_bar, 0, "Cell size:", self.var_umap_dot, len(self.umap_size_names), self.apply_umap_dot_size,
                            lambda i: self.umap_size_names[i], None,
                            label_tip="The size of the cells' dots in the UMAP (the same setting as UMAP circle size in Settings)."),
            self.add_slider(umap_bar, 1, "Label size:", self.var_label, len(LABEL_SIZES), self.apply_label_size,
                            lambda i: LABEL_SIZE_NAMES[i], None)]
        ttk.Button(umap_bar, text="Reset\nslider\ndefault", width=0, style=self.stacked_style,
                   command=self.reset_umap_sliders).grid(row=0, column=0, rowspan=2, sticky="nsw")
        umap_bar.columnconfigure(4, weight=1)      # (the spare width goes to an empty column on the right: the controls stay where they are)
        c = self.coords
        padx = (c[:, 0].max() - c[:, 0].min()) * 0.05 or 1.0
        pady = (c[:, 1].max() - c[:, 1].min()) * 0.05 or 1.0
        self.umap = PointCanvas(mid, (c[:, 0].min() - padx, c[:, 0].max() + padx, c[:, 1].min() - pady, c[:, 1].max() + pady), flip_y=True)
        self.umap.shapes_while_moving = True                             # (no plain discs while zooming and then the shapes again: that is two drawings)
        self.umap.canvas.pack(fill=tk.BOTH, expand=True)
        self.umap.add_marker("hover", 9)
        self.umap.on_motion = self.umap_motion
        self.umap.drag_hook = self.drag_state
        self.umap.on_view = None
        self.umap_tree = None if self.no_umap else cKDTree(c)
        self.mid = mid
        if self.umap_shown:
            pw.add(mid, weight=3)
        side = ttk.Frame(pw, padding=10, width=420)
        side.pack_propagate(False)
        self.build_side(side)
        pw.add(side, weight=0)
        self.side_pane = side
        self.root.update_idletasks()
        self.fit_left_pane(False)                                        # (the panes at their final widths before the section maps are made: they are laid out once)
        self.root.update_idletasks()                                     # (so that the section pane has its new width when the maps are sized from it: until the divider's move is carried out, it still reports the old one, the maps came out 1/5 of the old pane and then 1/6 of the new)
        self.build_tiles()
        self.root.after(400, self.fit_left_pane)                         # (and looked at again once the window has settled: a no-op if all is as it should be)

    def toggle_umap(self):
        """The Show UMAP box."""
        want = bool(self.var_show_umap.get())
        if want and not self.has_umap_data:
            self.var_show_umap.set(False)
            messagebox.showinfo("Show UMAP", "A UMAP is not available for this project: it was made without one.", parent=self.root)
            return
        self.set_umap_shown(want)

    def set_umap_shown(self, on):
        """Show or hide the UMAP pane. Hidden, the section maps take its room (they are laid out again for the wider pane), and Hover highlight, which
        marks cells in the UMAP too, is greyed out."""
        if on == self.umap_shown:
            return
        self.umap_shown = on
        if on:
            self.pw.insert(1, self.mid, weight=3)
        else:
            self.pw.forget(self.mid)
        self.hover_box.state(["!disabled"] if on else ["disabled"])
        if not on:
            self.drop_highlights()
        self.update_marker_slider()
        self.root.update_idletasks()
        self.fit_left_pane(False)

    def hover_on(self):
        """Hover highlight is on: its box is ticked and the UMAP is shown (it marks the group there too)."""
        return bool(self.var_hover.get()) and self.umap_shown

    def fit_left_pane(self, again=True):
        """The section pane gets about a third of the window (the panels fill whatever width it has). Done before the panels are made, so that they are made
        once; looked at again after the window has settled, which moves nothing when all is as it should be (a move is a second drawing of everything)."""
        if not self.root.winfo_viewable() or self.own is not None:
            return
        total = self.pw.winfo_width()
        side_w = self.side_pane.winfo_width() if self.side_pane.winfo_width() > 50 else int(self.side_pane.cget("width")) or 420       # (as wide as it is: the user may have dragged it wider)
        edge = total - side_w                                           # (the side panel keeps its own width, at the right)
        if not self.umap_shown:                                         # (no UMAP shown: the section maps fill what the side panel leaves)
            if abs(self.pw.sashpos(0) - int(edge)) > 2:
                self.pw.sashpos(0, int(edge))
            if again:
                self.root.after(600, lambda: self.fit_left_pane(False))
            return
        want = min(max(total * 0.33, 520), total * 0.5, edge - 300)     # (and the UMAP keeps at least 300)
        if abs(self.pw.sashpos(0) - int(want)) > 2:
            self.pw.sashpos(0, int(want))
        if abs(self.pw.sashpos(1) - int(edge)) > 2:
            self.pw.sashpos(1, int(edge))
        if again:                                                       # (once more after the window has settled)
            self.root.after(600, lambda: self.fit_left_pane(False))

    def build_side(self, f):
        self.top_row = ttk.Frame(f)                                       # Color by (left) and Level (right), each a column of radio buttons
        self.top_row.pack(fill=tk.X)
        self.top_row.columnconfigure(0, weight=1)
        self.top_row.columnconfigure(1, weight=1)
        color_col = ttk.LabelFrame(self.top_row, text="Color by", padding=6)             # (boxed like Level)
        color_col.grid(row=0, column=0, sticky="nsew")
        self.color_label = color_col
        self.var_mode = tk.StringVar(f, "all")
        modes = ttk.Frame(color_col)
        modes.pack(fill=tk.X)
        shown = [(value, text) for value, text in (("all", "All groups"), ("single", "Specific group(s)"), ("gene", "Gene"), ("imputed", "Imputed Gene"))
                 if not (value == "imputed" and self.own is not None)]     # (there is no imputed dataset behind a dataset of the user's own)
        for k, (value, text) in enumerate(shown):
            depends = value in ("all", "single")             # (these two work on the groups of the level chosen in Group by, to the right: an arrow points to it)
            rb = ttk.Radiobutton(modes, text=text + ("  →" if depends else ""), value=value, variable=self.var_mode, command=self.on_mode)
            rb.grid(row=k, column=0, sticky="w")
            if depends:
                HoverTip(rb, lambda: "What this shows depends on the level chosen in Group by, to the right.")
        self.var_translucent = tk.BooleanVar(f, False)
        translucent_box = ttk.Checkbutton(color_col, text="Translucent", variable=self.var_translucent, command=self.apply_opacity)
        translucent_box.pack(anchor="w", pady=(4, 0))                      # (under the Color by buttons, in the room beside the Level list)
        HoverTip(translucent_box, lambda: f"Draw the dots {TRANSLUCENT * 100:g}% opaque instead of fully opaque,\n"
                                          "so that where they overlap you see how many there are.\n\n"
                                          "Zooming and moving the plots\n"
                                          "become much slower while this is on.")
        self.level_frame = ttk.LabelFrame(self.top_row, text="Group by", padding=6)
        self.level_frame.grid(row=0, column=1, sticky="nsew", padx=(8, 0))
        HoverTip(self.level_frame, self.level_tip)
        self.var_level = tk.StringVar(f, self.level)
        self.build_level_radios()
        self.box = ttk.Frame(f)
        self.box.pack(fill=tk.X, pady=(8, 0))
        ttk.Frame(self.box, height=1).pack(fill=tk.X)             # (always something in it: a frame with nothing left in it keeps the height its last contents gave it, as the gene
                                                                  # box's, and the group list under it was pushed down by it)
        # a group number box (Specific group(s))
        self.id_frame = ttk.LabelFrame(self.box, text="Group number(s)", padding=6)
        self.id_hint = ttk.Label(self.id_frame, text="One or more, separated by commas (Enter):", foreground="gray30", wraplength=300)
        self.id_hint.pack(anchor="w")
        self.var_ids = tk.StringVar(f, "")
        self.ids_by_level = {}                              # level -> the group numbers last typed for it (see on_level)
        e = ttk.Entry(self.id_frame, textvariable=self.var_ids)
        e.pack(fill=tk.X, pady=(2, 0))
        e.bind("<Return>", lambda ev: self.recolor())
        self.id_entry = e
        # the gene box
        self.gene_frame = ttk.LabelFrame(self.box, text="Gene(s)", padding=6)
        ttk.Label(self.gene_frame, text="Up to six, separated by commas: the first three red / green / blue circles, the next three the same as small dots:", foreground="gray30", wraplength=300).pack(anchor="w")
        self.gene_entry = GeneBox(self.gene_frame)
        self.var_gene = GeneVar(self.gene_entry)
        self.gene_entry.pack(fill=tk.X)
        self.gene_entry.bind("<KeyRelease>", self.on_gene_key)
        self.gene_entry.bind("<Return>", lambda ev: (self.recolor(), "break")[1])
        self.gene_entry.bind("<Down>", lambda ev: self.focus_suggestions())
        self.gene_warning = ttk.Label(self.gene_frame, text="", foreground="#c00000", wraplength=300, justify=tk.LEFT)       # (shown under the box when a gene is rejected)
        self.sugg = tk.Listbox(self.gene_frame, height=6, exportselection=False, activestyle="none")
        self.sugg.bind("<ButtonRelease-1>", lambda ev: self.pick_suggestion())
        self.sugg.bind("<Return>", lambda ev: self.pick_suggestion())
        sc = ttk.Frame(self.gene_frame)
        sc.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(sc, text="Scale:").pack(side=tk.LEFT)
        self.var_scale = tk.StringVar(f, "log")
        for value, text in (("linear", "Linear"), ("log", "Log")):
            ttk.Radiobutton(sc, text=text, value=value, variable=self.var_scale, command=self.on_scale).pack(side=tk.LEFT, padx=(8, 0))
        th = ttk.Frame(self.gene_frame)
        th.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(th, text="Theme:").pack(side=tk.LEFT)
        self.var_theme = tk.StringVar(f, "light")
        self.theme_choice = None                # (None until Dark or Light is clicked: then what was chosen, kept. Until then several genes are shown dark, one gene light)
        self.theme_buttons = []
        for value, text in (("dark", "Dark"), ("light", "Light")):
            b = ttk.Radiobutton(th, text=text, value=value, variable=self.var_theme, command=self.on_theme)
            b.pack(side=tk.LEFT, padx=(8, 0))
            self.theme_buttons.append(b)
            HoverTip(b, lambda: "Dark: a black background, as for the screen. Light: a white background, for a paper. One gene: a scale from light grey "
                                "(none) to dark red. Several genes: cyan, magenta and yellow inks on light grey, that mix to blue, red and green.")
        self.cbar = tk.Canvas(self.gene_frame, width=300, height=70, highlightthickness=0)
        self.cbar.pack(fill=tk.X, pady=(6, 0))
        self.cbar_photo = None
        self.cbar.bind("<Button-1>", self.on_legend_click)
        self.gene_levels = (0.0, 1.0)          # the black and white levels of a gene, as fractions of its range (the slider under the color bar)
        self.multi_levels = {}                 # the same for each of several genes (gene number -> (low, high)): the sliders in the legend
        self.multi_sliders = []
        self.multi_windows = []
        self.levels_job = None
        self.cbar_slider = RangeSlider(self.gene_frame, self.on_gene_levels)
        self.hist = tk.Canvas(self.gene_frame, height=HIST_HEIGHT + HIST_PAD, highlightthickness=0)
        self.cell_marks = None                                 # {gene number: the value of the cell under the pointer}: a dot on each histogram (see mark_cell_on_histograms)
        self.hist_args = None                                  # (what it shows now, to draw again when its width changes: see redraw_hist)
        self.hist_cache = {}                                   # id of a values array -> (the array, its histogram data): the sliders redraw on every drag
        self.hist.hist_regions, self.hist.hist_photos = [], []
        self.hist.bind("<Configure>", lambda e: self.redraw_hist())
        self.legend_redraw, self.cbar_width, self.cbar_job = None, 0, None      # (how to draw the color bar or the several genes' legend again when the panel changes width)
        self.cbar.bind("<Configure>", lambda e: self.cbar_resized(), add="+")
        self.hist_pointer = [None]
        for canvas in (self.hist, self.cbar):
            canvas.bind("<Motion>", lambda e: self.hist_pointer.__setitem__(0, e), add="+")          # (before the tips' own: they ask for the text as the pointer moves)
            canvas.hist_regions, canvas.hist_photos = canvas.__dict__.get("hist_regions", []), canvas.__dict__.get("hist_photos", [])
            HoverTip(canvas, lambda c=canvas: self.hist_pointer[0] and self.histogram_readout(c, self.hist_pointer[0]))
        HoverTip(self.cbar_slider, lambda: "Drag the black and white handles to set the lowest and the highest value of the color scale. "
                                        + ("The track is even." if self.cbar_slider.even else "The track is by the square root of the value: most of it is for the low values. ")
                                        + " Double-click to reset the handles to where they started.")
        # the legend for Single group(s)
        self.legend = tk.Canvas(f, width=300, height=10, highlightthickness=0)
        self.select_bar = ttk.Frame(f)                      # All / None, above the group list and the Single group(s) legend
        all_button = ttk.Button(self.select_bar, text="Show all", width=10, style="Side.TButton", command=self.select_all_groups)
        all_button.pack(side=tk.LEFT)
        none_button = ttk.Button(self.select_bar, text="Show none", width=10, style="Side.TButton", command=self.select_no_groups)
        none_button.pack(side=tk.LEFT, padx=(6, 0))
        HoverTip(all_button, lambda: "Show all groups")
        HoverTip(none_button, lambda: "Show no groups (every group gray)")
        self.group_list = ttk.Frame(f)                      # the groups' numbers and colors (All groups), scrollable
        self.group_header = tk.Canvas(self.group_list, height=20, highlightthickness=0, cursor="hand2")      # the column headings: click to sort
        self.group_header.pack(side=tk.TOP, fill=tk.X, pady=(0, 2))
        self.group_header.bind("<Button-1>", self.on_group_header_click)
        self.group_sort = ("id", False)                    # (column, descending): kept when the level changes
        self.group_canvas = tk.Canvas(self.group_list, width=300, height=40, highlightthickness=0, bg="white")
        gsb = ttk.Scrollbar(self.group_list, orient=tk.VERTICAL, command=self.group_canvas.yview)
        self.group_canvas.configure(yscrollcommand=gsb.set)
        gsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.group_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.group_canvas.bind("<Button-1>", self.on_group_list_click)
        self.group_canvas.bind("<Motion>", self.on_group_list_motion)
        self.group_canvas.bind("<Leave>", lambda e: self.show_legend_group(None))
        self.group_canvas.config(cursor="hand2")
        self.legend.bind("<Button-1>", self.on_single_legend_click)
        self.legend.bind("<Motion>", self.on_single_legend_motion)
        self.legend.bind("<Leave>", lambda e: self.show_legend_group(None))
        self.legend.config(cursor="hand2")
        self.group_canvas.bind("<MouseWheel>", lambda e: self.group_canvas.yview_scroll(-1 if e.delta > 0 else 1, "units"))
        self.build_buttons(f)
        # (No "Under the pointer" box: the hovered cell is described in a tip (CellTip), and so is a row of the group list, which then has the room.)
        self.hover_job = None
        self.group_tip = CellTip(self.root, delay_ms=GROUP_TIP_MS)
        self.var_status = tk.StringVar(f, "")
        ttk.Label(f, textvariable=self.var_status, foreground="gray30", wraplength=300, justify=tk.LEFT).pack(anchor="w", side=tk.BOTTOM)
        ttk.Label(f, text="wheel: zoom   drag: pan   double-click: whole view", foreground="gray40", wraplength=300).pack(anchor="w", side=tk.BOTTOM, pady=(0, 6))
        self.update_frames()

    def build_level_radios(self):
        for w in self.level_frame.winfo_children():
            w.destroy()
        shown = [k for k in self.level_names if k in self.levels or (self.no_umap and k == LEIDEN_KEY)]
        for k, name in enumerate(shown):
            rb = ttk.Radiobutton(self.level_frame, text=self.level_names[name], value=name, variable=self.var_level, command=self.on_level)
            rb.grid(row=k, column=0, sticky="w")                       # (one column, coarse to fine, top to bottom)
            if name not in self.levels:
                rb.state(["disabled"])                                 # (Leiden clusters come from the UMAP's computation: a project without a UMAP has none)
                HoverTip(rb, lambda: "Leiden clusters are found along with the UMAP; this project was made without one.")
            else:
                HoverTip(rb, self.level_tip)

    def level_tip(self):
        if self.mode in ("gene", "imputed"):
            return ("The group that is marked with red + signs when the pointer rests on a cell: every cell of that cell's group at this level.\n"
                    "It does not change how the cells are colored in this mode.")
        return "The level the cells are colored by, and the group that is marked with red + signs when the pointer rests on a cell."

    def build_buttons(self, f):
        ttk.Style(self.root).configure("Side.TButton", font=("Segoe UI", 9))
        bar = ttk.Frame(f)
        bar.pack(fill=tk.X, side=tk.TOP, pady=(0, 10), before=self.top_row)
        bar.columnconfigure(0, weight=1)
        bar.columnconfigure(1, weight=1)
        self.group_buttons = []
        rows = ((("Settings...", self.open_settings, False), ("Open data folder", self.open_data_folder, False)),
                (("Save section maps", self.save_section_maps, False), ("Save UMAP", self.save_umap, False)),
                (("Export DEGs", self.export_degs, True), ("Export expression", self.export_expression, True)),
                (("Export coordinates", self.export_coordinates, True), ("3D in browser", self.show_3d_browser, True)),
                (("Analyze...", self.open_analyze, False), ("Create subproject", self.analyze_subset, False)),
                (("Close", self.close, False),))
        for r, pair in enumerate(rows):
            for c, (text, cmd, grouped) in enumerate(pair):
                b = ttk.Button(bar, text=text, command=cmd, style="Side.TButton")
                b.grid(row=r, column=c, columnspan=2 if len(pair) == 1 else 1, sticky="ew", padx=2, pady=2)
                if grouped:
                    self.group_buttons.append(b)
                if self.own is not None and text == "Save section maps":
                    b.state(["disabled"])                                 # (no section maps)
                if self.no_umap and text == "Save UMAP":
                    b.state(["disabled"])
                    HoverTip(b, lambda: "This project was made without a UMAP.")
                if text == "Create subproject":
                    HoverTip(b, lambda: "Analyze only the groups switched on in the legend: a project of its own, in a subfolder of this one, with its own "
                                        "PCA and UMAP, opened in a second window.")
        self.update_group_buttons()

    def analyze_subset(self):
        """The Analyze subset button: the groups switched on in the legend, analyzed again on their own (abc_tk_subset.py)."""
        import abc_tk_subset
        abc_tk_subset.analyze_subset(self)

    def update_group_buttons(self):
        """The group exports and 3D views are for the groups typed in Single group(s)."""
        for b in self.group_buttons:
            b.state(["!disabled"] if self.mode == "single" and self.own is None else ["disabled"])

    def update_frames(self):
        self.level_frame.config(text="Group by" if self.mode in ("all", "single") else "Group by (hover)")
        if self.level_has_ids.get(self.level, True):                 # (what the box takes: the groups' numbers, or their names where they have no numbers)
            self.id_frame.config(text="Group number(s)")
            self.id_hint.config(text="One or more, separated by commas (Enter):")
        else:
            self.id_frame.config(text="Group name(s)")
            self.id_hint.config(text="One or more names (a start will do: GABA, Glut), separated by commas (Enter):")
        for w in (self.id_frame, self.gene_frame, self.select_bar, self.legend, self.group_list):
            w.pack_forget()
        if self.mode == "single":
            self.id_frame.pack(fill=tk.X)
            self.id_entry.focus_set()
        elif self.mode in ("gene", "imputed"):
            self.gene_frame.pack(fill=tk.X)
            self.gene_entry.focus_set()

    # ---- the section panels

    def build_tiles(self):
        self.tiles = []
        if not self.sections:
            return
        boxes = {}
        for s in self.sections:
            idx = self.run_idx[s]
            xs, ys = self.cx[idx], self.cy[idx]
            if self.bg is not None and s in self.bg["where"]:          # (the whole section, not only the run's cells in it)
                w = self.bg["where"][s]
                xs, ys = self.bg["x"][w], self.bg["y"][w]
            ok = np.isfinite(xs) & np.isfinite(ys)
            if ok.any():
                boxes[s] = (xs[ok].min(), xs[ok].max(), ys[ok].min(), ys[ok].max())
        span_x = max(b[1] - b[0] for b in boxes.values()) * 1.04 or 1.0          # (one size for every panel, so they share a scale)
        span_y = max(b[3] - b[2] for b in boxes.values()) * 1.04 or 1.0
        self.tile_aspect = span_x / span_y
        for s in self.sections:
            if s not in boxes:
                continue
            b = boxes[s]
            mx, my = (b[0] + b[1]) / 2, (b[2] + b[3]) / 2
            fr = tk.Frame(self.tile_scroll, bg="#f4f4f4")
            item = self.tile_scroll.create_window(0, 0, window=fr, anchor="nw")
            token = core.sanitize_section_token(s)
            title = tk.StringVar(fr, token)
            caption = tk.Label(fr, textvariable=title, bg="#f4f4f4", font=("Segoe UI", 9, "bold"))
            caption.pack(anchor="w")
            self.bind_pane_gestures(fr)
            self.bind_pane_gestures(caption)
            pc = PointCanvas(fr, (mx - span_x / 2, mx + span_x / 2, my - span_y / 2, my + span_y / 2), flip_y=False, width=TILE_PX, height=int(TILE_PX / self.tile_aspect))
            pc.canvas.pack()
            pc.canvas.bind("<Shift-MouseWheel>", lambda e: self.tile_scroll.yview_scroll(-1 if e.delta > 0 else 1, "units"))     # (the plain wheel zooms)
            pc.add_marker("hover", 7)
            pc.plus_arm = 2
            pc.dot_factor = self.dot_factors[self.var_dot.get()]
            pc.opacity = TRANSLUCENT if self.var_translucent.get() else 1.0
            pc.marker_factor = MARKER_SIZES[self.var_marker.get()]
            pc.visible_fn = lambda fr=fr: self.frame_visible(fr)
            pc.shapes_while_moving = True
            pc.panel_zoom = self.panel_zoom
            pc.panel_hook = self.panel_hook if self.var_panel_zoom.get() else None
            if not self.tiles:
                pc.enable_scalebar()                       # (the first section's panel carries the scale bar)
            tile = {"item": item, "label": s, "frame": fr, "pc": pc, "title": title, "token": token, "tree": None, "caption": caption}
            pc.on_motion = lambda px, py, t=tile: self.tile_motion(t, px, py)
            pc.drag_hook = self.drag_state
            pc.on_change = lambda zoom, dx, dy, t=tile: self.link_tiles(t, zoom, dx, dy)
            for r in self.rois:
                if r[0] == str(s) or core.sanitize_section_token(r[0]) == token:
                    pc.add_rect(r[1], r[3], r[2], r[4])
            self.tiles.append(tile)
        self.tile_cols = self.tile_px = 0
        self.layout_tiles()

    def apply_ccf_opacity(self):
        """The CCF atlas slider: the cells of the section maps are drawn again, more or less see-through (the outlines are not)."""
        if not self.var_ccf.get():
            return
        a = OPACITIES[self.var_ccf_opacity.get()]
        for t in self.tiles:
            t["pc"].cell_alpha = a
            if self.frame_visible(t["frame"]):
                t["pc"].request_draw()
            else:
                t["pc"].stale = True

    def apply_ccf(self):
        """The CCF atlas checkbox: show or hide the regions' outlines on the section maps (made the first time it is ticked)."""
        self.ccf_scale.state(["!disabled"] if self.var_ccf.get() else ["disabled"])
        if not self.var_ccf.get():
            for t in self.tiles:
                t["pc"].cell_alpha = 1.0                                   # (the cells as they are again)
                t["pc"].set_outlines(None)
            return
        if self.ccf_outlines is None and not self.load_ccf():
            self.var_ccf.set(False)
            return
        a = OPACITIES[self.var_ccf_opacity.get()]
        for t in self.tiles:
            t["pc"].cell_alpha = a
            t["pc"].set_outlines(self.ccf_outlines.get(t["label"]))

    def load_ccf(self):
        """Make the outlines of every section map: the volume, the cells' CCF positions, each section's fit. True if done."""
        h = self.host
        if not self.tiles:
            return False
        if core.reconstructed_in_use() and self.bg is not None and self.bg.get("z") is not None:
            return self.load_ccf_from_slices()
        if not h["ensure_ccf_downloaded"](self.abc_cache, parent=self.root):
            self.var_status.set("The CCF coordinates were not downloaded, so there is no atlas to draw.")
            return False
        cells, homes = {}, {}
        for t in self.tiles:
            s = t["label"]
            if self.bg is not None and s in self.bg["where"]:          # (every cell of the section: a better fit than the run's own)
                w = self.bg["where"][s]
                ids, x, y = self.bg["ids"][w], self.bg["x"][w], self.bg["y"][w]          # (indexed first, then made an array: converting every id of the atlas for each section is slow)
            else:
                idx = self.run_idx[s]
                ids, x, y = self.adata.obs_names[idx], self.cx[idx], self.cy[idx]
            if len(ids) > CCF_FIT_CELLS_PER_SECTION:                    # (a section's fit needs tens of thousands of cells, not millions)
                keep = np.sort(np.random.default_rng(0).choice(len(ids), CCF_FIT_CELLS_PER_SECTION, replace=False))
                ids, x, y = ids[keep], x[keep], y[keep]
            cells[s] = (np.asarray(ids, dtype=object), x, y)
            homes[s] = t["pc"].home


        def work():
            path = self.abc_cache.get_data_path(directory=h["CCF_VOLUME_DIRECTORY"], file_name=h["CCF_VOLUME_FILE_NAME"])
            volume = h["read_nifti_volume"](path)
            all_ids = np.concatenate([c[0] for c in cells.values()])
            aligned, _shell = h["load_ccf_for_cells"](self.abc_cache, all_ids.tolist())
            xyz = aligned[["x_ccf", "y_ccf", "z_ccf"]].to_numpy(dtype=float)
            at, fitted = 0, {}
            for s, (ids, x, y) in cells.items():
                fitted[s] = (x, y, xyz[at:at + len(ids)])
                at += len(ids)
            return abc_tk_ccf.section_outlines(volume, fitted, homes)

        try:
            self.ccf_outlines = h["run_with_loading_window"]("Drawing the CCF atlas on the section maps (the first read of the CCF coordinates takes about a minute).", work, show_after_s=0)
        except Exception as e:
            self.var_status.set(f"Could not draw the CCF atlas: {e}")
            return False
        empty = [core.sanitize_section_token(s) for s, p in self.ccf_outlines.items() if not p]
        if empty:
            self.var_status.set(f"No CCF position for section {', '.join(empty)}: no outlines there.")
        return True

    def load_ccf_from_slices(self):
        """The outlines when the cells are placed by Allen's reconstructed coordinates: the CCF annotation resampled into that same space, one slice per
        section, is read off as it is (the cells' own region labels are exactly its values at their positions), so nothing is fitted. True if done."""
        h = self.host
        homes = {t["label"]: t["pc"].home for t in self.tiles}
        section_z = {}
        for t in self.tiles:
            w = self.bg["where"].get(t["label"])
            section_z[t["label"]] = float(np.nanmedian(self.bg["z"][w])) if w is not None and len(w) else float("nan")

        def work():
            path = self.abc_cache.get_data_path(directory=core.CCF_DIRECTORY, file_name=core.RECON_ANNOTATION_FILE_NAME)
            return abc_tk_ccf.slice_outlines(h["read_nifti_volume"](path), section_z, homes)

        try:
            self.ccf_outlines = h["run_with_loading_window"]("Drawing the CCF atlas on the section maps.", work, show_after_s=0)
        except Exception as e:
            self.var_status.set(f"Could not draw the CCF atlas: {e}")
            return False
        empty = [core.sanitize_section_token(s) for s, p in self.ccf_outlines.items() if not p]
        if empty:
            self.var_status.set(f"No CCF atlas for section {', '.join(empty)} (Allen did not register it): no outlines there.")
        return True

    def stop_slider(self, parent, shown, var, n, command, value_text, length=130):
        """A slider of `parent` with `n` stops (0 .. n-1), not yet placed: it jumps from stop to stop, shows `value_text(stop)` in the label `shown`, and
        calls `command()` when the stop changes. The same width wherever it is (but for `length`, in pixels)."""
        scale = ttk.Scale(parent, from_=0, to=n - 1, orient=tk.HORIZONTAL, length=length)
        pending = []

        def run():
            pending.clear()
            command()

        def moved(v):
            i = int(round(float(v)))
            if abs(float(v) - i) > 1e-6:
                scale.set(i)                               # (snaps to the stop, which calls this again)
                return
            shown.config(text=value_text(i))
            if var.get() != i:
                var.set(i)
                # The command (laying out and drawing all the section maps again, say) is slow: run inside the drag it held up the mouse events,
                # which the slider then went on working through after the button was let go. So the thumb moves at once, and the work follows.
                if pending:
                    scale.after_cancel(pending[0])
                pending[:] = [scale.after(30, run)]

        scale.config(command=moved)
        scale.set(var.get())
        shown.config(text=value_text(var.get()))
        return scale

    def add_slider(self, parent, row, text, var, n, command, value_text, tip_text, label_tip=None):
        """A slider with `n` stops (0 .. n-1): it jumps from stop to stop, shows `value_text(stop)` beside it, and calls `command()` when the
        stop changes. `tip_text(stop)`, if given, is shown in a tip while the pointer is over the slider, and `label_tip` (a string) over
        its label. Returns the slider."""
        label = ttk.Label(parent, text=text)
        label.grid(row=row, column=1, sticky="w", padx=(12, 6))
        if label_tip is not None:
            HoverTip(label, lambda: label_tip)
        shown = ttk.Label(parent)
        shown.grid(row=row, column=3, sticky="w", padx=(6, 0))
        # The column of values as wide as the longest value its sliders can show, and no wider (a fixed width in characters left empty space
        # beside the shorter ones, room the checkboxes to the right lacked on a small screen); never narrower, so nothing moves as they change.
        font = tkfont.nametofont(str(shown.cget("font")) or "TkDefaultFont")
        widest = max(font.measure(value_text(i)) for i in range(n)) + 12      # (+ the label's padx and its own margin)
        parent.columnconfigure(3, minsize=max(widest, int(parent.grid_columnconfigure(3)["minsize"])))
        scale = self.stop_slider(parent, shown, var, n, command, value_text)
        scale.grid(row=row, column=2, sticky="w")
        scale.companions = (label, shown)                # (greyed out with it: see set_slider_enabled)
        if tip_text is not None:
            HoverTip(scale, lambda: tip_text(var.get()))
        return scale

    @staticmethod
    def set_slider_enabled(scale, on):
        """Grey out (or not) a slider made by add_slider, with its label and the value shown beside it."""
        for w in (scale, *getattr(scale, "companions", ())):
            w.state(["!disabled"] if on else ["disabled"])

    def update_marker_slider(self):
        """The Marker size slider is for the hover marks: greyed out while Hover highlight is off."""
        self.set_slider_enabled(self.sliders[2], self.hover_on())

    def apply_label_size(self):
        """The UMAP's Label size slider: the group numbers drawn on the UMAP, bigger or smaller."""
        self.umap.set_label_factor(LABEL_SIZES[self.var_label.get()])

    def reset_umap_sliders(self):
        for scale, default in zip(self.umap_sliders, (self.umap_size_index(self.umap_default_size), DEFAULT_LABEL_INDEX)):
            scale.set(default)

    def umap_size_index(self, name):
        """The Cell size slider's stop for a UMAP circle size of Settings (the middle one, Default, if it is not one of them)."""
        names = self.umap_size_names
        return names.index(name) if name in names else (names.index("Default") if "Default" in names else 0)

    def apply_umap_dot_size(self):
        """The UMAP's Cell size slider: the UMAP circle size of Settings, saved, and the UMAP drawn again with it."""
        name = self.umap_size_names[self.var_umap_dot.get()]
        if name == self.viewer_settings.get("umap_circle_size", "Default"):
            return
        self.viewer_settings["umap_circle_size"] = name
        if "save_viewer_settings" in self.host:
            self.host["save_viewer_settings"](self.viewer_settings)
        if not self.busy:
            self.recolor()

    def view_background(self):
        """The background of the UMAP and the section maps: white, dark in the dark theme of the gene views."""
        return MULTI_BG if getattr(self, "tiles_dark", False) else (255, 255, 255)

    def reset_sliders(self):
        for scale, default in zip(self.sliders, (self.panel_zoom_tile_index() if self.var_panel_zoom.get() else DEFAULT_TILE_INDEX, DEFAULT_DOT_INDEX, DEFAULT_MARKER_INDEX)):
            scale.set(default)

    def layout_tiles_when_settled(self):
        """The section pane changed size (the divider or the window is being dragged): lay the panels out once it has stopped changing, not at
        every step, which redrew all of them each time."""
        if not self.layout_seen:                                    # (the first time, when the window opens, at once)
            self.layout_seen = True
            self.layout_tiles()
            return
        self.trail(f"pane resized to {self.tile_scroll.winfo_width()}x{self.tile_scroll.winfo_height()}")
        if self.layout_job is not None:
            self.root.after_cancel(self.layout_job)
        self.layout_job = self.root.after(LAYOUT_SETTLE_MS, self.layout_after_settling)

    def layout_after_settling(self):
        self.layout_job = None
        self.layout_tiles()

    def layout_tiles(self):
        """The panels are a third, a half or all of the pane's width (the Section size slider), as many to a row as that gives; the
        pane scrolls to the rest. While the whole panel is zoomed they are bigger than the pane, by up to 100 times: each is then a virtual
        canvas (PointCanvas.set_virtual), of which the widget is only the part in view (place_tiles)."""
        if not getattr(self, "tiles", None):
            return
        cols, base = self.tile_params()
        px = int(base * self.panel_zoom)
        same = cols == self.tile_cols and px == self.tile_px
        self.trail(f"layout_tiles cols={cols} px={px} zoom={self.panel_zoom:.3g} " + ("(unchanged: nothing done)" if same else "(laid out again)"))
        if same:
            return
        self.tile_cols, self.tile_px = cols, px
        big = self.panel_zoom > 1.0
        ph = int(px / self.tile_aspect)
        for t in self.tiles:
            pc = t["pc"]
            pc.outline_px = OUTLINE_PX_LARGEST if cols == 1 else 1            # (the CCF outlines: wider at the largest size)
            pc.rect_px = ROI_LINE_PX.get(cols, 2)
            if pc.rects and pc.rect_line != pc.rect_px:
                pc.reposition()                                                         # (the ROI outlines take the new line width)
            if big:
                pc.set_virtual((px, ph))
            else:
                pc.set_virtual(None)
                pc.canvas.config(width=px, height=ph)
                t["frame"].slot_wh = None
                t["shown"] = True
                if not t["caption"].winfo_manager():                # (not winfo_ismapped: that is only true once Tk has had its idle turn after a pack)
                    t["caption"].pack(anchor="w", before=pc.canvas)
                self.tile_scroll.itemconfig(t["item"], state="normal")
        self.root.update_idletasks()                                  # (so that the panels know their size)
        if big:
            ch = max(t["caption"].winfo_reqheight() for t in self.tiles)
            fw, fh = px, ch + ph
        else:
            fw = max(t["frame"].winfo_reqwidth() for t in self.tiles)
            fh = max(t["frame"].winfo_reqheight() for t in self.tiles)
        for k, t in enumerate(self.tiles):                            # (8 across and 6 down between panels, as the grid had)
            x, y = 4 + (k % cols) * (fw + 8), 3 + (k // cols) * (fh + 6)
            t["frame"].tile_x, t["frame"].tile_y = x, y               # (where it is in the scrolled canvas: see frame_visible)
            t["frame"].slot_wh = (fw, fh) if big else None            # (the room it takes there, when its widget is only the part in view)
            self.tile_scroll.coords(t["item"], x, y)
        rows = -(-len(self.tiles) // cols)
        pad_x, pad_y = self.pane_slack()
        self.tile_scroll.configure(scrollregion=(0, 0, 4 + cols * (fw + 8) + pad_x, 3 + rows * (fh + 6) + pad_y))
        self.root.after_idle(self.refresh_visible)

    def place_tiles(self, max_draws=None, background=None, started=None):
        """While the whole panel is zoomed: each section map's widget is put where the part of the map that is on screen is (and made that big, a
        little more all round), and drawn there; a map not on screen has no widget shown. A map whose widget still covers what is on screen is left.
        `max_draws`: at most this many maps are drawn (those with the most on screen first); while the panel is being dragged, none unless it says.
        `background`: a function; the map is then drawn in a thread and the function called when that is done (see PointCanvas.draw), at most one
        map to a call. `started()` is called as each map's draw is about to begin (before it can end: a map that has nothing to draw ends at once). Returns how many
        more are waiting to be drawn."""
        if self.panel_zoom <= 1.0 or not getattr(self, "tiles", None):
            return 0
        if max_draws is None and getattr(self, "panning", False):
            max_draws = 0                            # (a drag goes on at once, the new parts blank: pan_fill draws them, one map at a time)
        waiting = []
        vw, vh = self.tile_scroll.winfo_width(), self.tile_scroll.winfo_height()
        vl, vt = self.tile_scroll.canvasx(0), self.tile_scroll.canvasy(0)
        mx, my = int(WINDOW_MARGIN * vw), int(WINDOW_MARGIN * vh)
        for t in self.tiles:
            fr, pc = t["frame"], t["pc"]
            if getattr(fr, "slot_wh", None) is None or pc.vsize is None:
                continue
            lw, lh = pc.vsize
            ch = fr.slot_wh[1] - lh                                   # (the caption above the map)
            cx0, cy0 = fr.tile_x, fr.tile_y + ch                      # (where the map is in the scrolled pane)
            nx0, nx1 = max(vl - cx0, 0), min(vl + vw - cx0, lw)       # (what is on screen, in the map's pixels)
            ny0, ny1 = max(vt - cy0, 0), min(vt + vh - cy0, lh)
            if nx1 - nx0 < 1 or ny1 - ny0 < 1:
                if t.get("shown", True):
                    self.tile_scroll.itemconfig(t["item"], state="hidden")
                    t["shown"] = False
                continue
            v = pc.vview
            if t.get("shown", True) and v is not None and v[0] <= nx0 and v[1] <= ny0 and v[0] + v[2] >= nx1 and v[1] + v[3] >= ny1:
                continue
            waiting.append(((nx1 - nx0) * (ny1 - ny0), t, ch, cx0, cy0, lw, lh))
        waiting.sort(key=lambda w: -w[0])
        if max_draws is None:
            max_draws = len(waiting)
        for _area, t, ch, cx0, cy0, lw, lh in waiting[:max_draws]:
            pc = t["pc"]
            x0, x1 = max(vl - mx - cx0, 0), min(vl + vw + mx - cx0, lw)
            y0, y1 = max(vt - my - cy0, 0), min(vt + vh + my - cy0, lh)
            with_caption = y0 <= 0
            # (winfo_manager, not winfo_ismapped: a caption just packed is not "mapped" until Tk's idle turn, and one taken for absent here was
            # then never hidden, and appeared later: the map's picture 31 pixels low in its frame, the vertical hiccup at the end of a zoom)
            if with_caption and not t["caption"].winfo_manager():
                t["caption"].pack(anchor="w", before=pc.canvas)
            elif not with_caption and t["caption"].winfo_manager():
                t["caption"].pack_forget()
            self.tile_scroll.coords(t["item"], cx0 + x0, cy0 + y0 - (ch if with_caption else 0))
            if not t.get("shown", True):
                self.tile_scroll.itemconfig(t["item"], state="normal")
                t["shown"] = True
            if background is not None:
                self.started_in_background += 1
                if started is not None:
                    started()
            pc.set_view(int(x0), int(y0), int(x1) - int(x0), int(y1) - int(y0), background)
        if waiting[:max_draws]:
            # A map's widget moved, resized and given its new picture: Tk does those in passes of its own (the pane moves the widget, the frame packs
            # its caption and canvas again, the canvas draws), and a screen refresh between them showed a moment of the old picture in the new
            # place, a blink. All of them now, together, so that only the finished map reaches the screen.
            self.root.update_idletasks()
        return max(len(waiting) - max_draws, 0)

    # ---- zooming the whole panel (the borders and the sections together)

    MAX_PANEL_ZOOM = 200.0

    def max_panel_zoom(self):
        return self.MAX_PANEL_ZOOM

    def panel_zoom_tile_index(self):
        """The Section size stop the whole-panel zoom starts at: the one that fits every section into the pane (whatever shape the pane is), as large as
        that allows, so that a few sections are not left small in a big pane and many are not left below its edge."""
        return FIT_STOP

    def tile_params(self):
        """(panels across, the width of a panel's map at panel zoom 1) for the Section size stop and the pane as it is now."""
        stop = TILE_COLUMNS[self.var_tile.get()]
        w = max(self.tile_scroll.winfo_width(), 100)
        if stop == 0:
            return self.fit_tiles(w)
        return stop, max((w - 4) // stop - 8 - 2, MIN_TILE_PX)               # (8: the gaps either side of a panel; 2: its canvas border)

    def pane_slack(self):
        """(width, height) of the room at the right and the bottom of the pane that the panels leave empty at panel zoom 1 (a fit that the height limits does not
        fill the width: 0 if they fill it). It stays part of what the pane scrolls over at every zoom: or a zoom close to the unzoomed one, with the content
        hardly wider than the pane, could not scroll as far as it takes to keep the point under the pointer in place, and the content slid instead."""
        cols, base = self.tile_params()
        n = max(len(self.tiles), 1)
        ch = max((t["caption"].winfo_reqheight() for t in self.tiles), default=20)
        fw, fh = base + 2, ch + int(base / self.tile_aspect) + 2                # (a panel's slot at zoom 1: its map with the border it has then)
        w = self.tile_scroll.winfo_width()
        h = self.tile_scroll.winfo_height() + (self.hsb.winfo_reqheight() if self.hsb.winfo_manager() else 0)
        return max(0, w - (4 + cols * (fw + 8))), max(0, h - (3 + -(-n // cols) * (fh + 6)))

    def fit_tiles(self, w):
        """(across, map width) for the largest panels that show every section in the pane (width `w`, its height now) at once: for each number of panels
        across, the panels are as wide as the width allows or, if the rows then run below the bottom, as the height allows (and then do not fill the
        width), and the number across that gives the biggest panels is taken. The same for a pane wide or tall, and for few sections or many. If even the
        smallest panels (MIN_TILE_PX) do not fit, they are that small, as many across as the width holds, and the pane scrolls."""
        # (the height as it is without the scroll bar under the pane, which comes while the panel is zoomed: or the size would change with the zoom)
        h = self.tile_scroll.winfo_height() + (self.hsb.winfo_reqheight() if self.hsb.winfo_manager() else 0)
        n = max(len(self.tiles), 1)
        ch = max((t["caption"].winfo_reqheight() for t in self.tiles), default=20)
        key = (w, h, n, ch, self.tile_aspect)
        hit = self.fit_cache.get(key)
        if hit is None:
            best = None
            for across in range(1, n + 1):
                px = (w - 4) // across - 8 - 2
                if px < MIN_TILE_PX:
                    break
                rows = -(-n // across)
                px = min(px, int(((h - 3) / rows - ch - 8) * self.tile_aspect) + 2)        # (a row is the caption, the map with its border, and a gap of 6; the loop below makes it exact)
                while px >= MIN_TILE_PX and 3 + rows * (ch + int(px / self.tile_aspect) + 2 + 6) > h:
                    px -= 1
                if px >= MIN_TILE_PX and (best is None or px > best[1]):
                    best = (across, px)
            if best is None:
                best = (max(1, (w - 4) // (MIN_TILE_PX + 10)), MIN_TILE_PX)
            self.fit_cache.clear()
            hit = self.fit_cache[key] = best
        self.fit_across = hit[0]
        return hit

    def apply_panel_zoom_mode(self):
        """The Zoom whole panel checkbox. On: every section map is put back to its whole view and the wheel then zooms the panel. Off: the panel is
        put back to its size, and the wheel zooms into each section again."""
        on = self.var_panel_zoom.get()
        for t in self.tiles:
            t["pc"].panel_hook = self.panel_hook if on else None
            if on:
                t["pc"].fit()
        self.cancel_zoom_burst()
        if on:
            self.tile_stop_before_zoom = self.var_tile.get()          # (the section maps start at a size that fits them all: see panel_zoom_tile_index)
            self.sliders[0].set(self.panel_zoom_tile_index())
        else:
            self.set_panel_zoom(1.0)
            if getattr(self, "tile_stop_before_zoom", None) is not None:
                self.sliders[0].set(self.tile_stop_before_zoom)       # (and back to what they were)
                self.tile_stop_before_zoom = None

    def tile_geometry(self, t):
        """(x, y of the panel in the scrolled pane, the caption's height, the map's width and height)."""
        fr, pc = t["frame"], t["pc"]
        ch = t["caption"].winfo_reqheight()
        if pc.vsize is not None:
            lw, lh = pc.vsize
        else:
            lw, lh = pc.canvas.winfo_width(), pc.canvas.winfo_height()
        return getattr(fr, "tile_x", 0), getattr(fr, "tile_y", 0), ch, max(lw, 1), max(lh, 1)

    def spot_in_tiles(self, ax, ay):
        """Where a point of the scrolled pane is in the panels, so that it can be found again after they are laid out bigger or smaller: the nearest panel
        and the point's offset from its map's corner. (The captions and the gaps do not grow with the zoom, the maps do: scaling the pane's
        coordinates all alike puts the point off by the captions and gaps of every row above it, which showed as a hop when the zoom settled.)"""
        best = None
        for k, t in enumerate(self.tiles):
            x, y, ch, lw, lh = self.tile_geometry(t)
            dx = max(x - ax, 0, ax - (x + lw))
            dy = max(y - ay, 0, ay - (y + ch + lh))
            if best is None or dx * dx + dy * dy < best[0]:
                best = (dx * dx + dy * dy, k, ax - x, ay - (y + ch), lw, lh)
        return best[1:] if best else None

    def spot_now(self, spot):
        """Where that spot is now (x, y in the scrolled pane): offsets within the map grow with the map, those outside it (caption, gaps) stay as they were."""
        k, dx, dy, old_w, old_h = spot
        x, y, ch, lw, lh = self.tile_geometry(self.tiles[k])

        return x + grown(dx, old_w, lw), y + ch + grown(dy, old_h, lh)

    def set_panel_zoom(self, zoom, anchor=None, spot=None, draw=True):
        """Enlarge the whole panel to `zoom` times its size (1 to MAX_PANEL_ZOOM): the panels are laid out again, bigger, so that borders, gaps and
        sections grow together, and only those in view are drawn. `anchor` = (x, y in the scrolled pane, x, y in the window): the point that stays
        under the pointer. `spot`: that point already found in the panels (see spot_in_tiles; for a layout that is not the one the anchor's x, y
        are of). `draw` False: only laid out, the maps drawn by the caller (zoom_fill)."""
        if getattr(self, "panning", False):
            self.end_pan()                           # (a drag just before: its drawing put off no longer, as the maps are laid out and drawn here)
        zoom = min(max(zoom, 1.0), self.max_panel_zoom())
        if abs(zoom - self.panel_zoom) < 1e-9 or not self.tiles:
            return
        if spot is None and anchor is not None:
            spot = self.spot_in_tiles(anchor[0], anchor[1])
        self.trail(f"set_panel_zoom {self.panel_zoom:.3g} -> {zoom:.3g}")
        self.panel_zoom = zoom
        for t in self.tiles:
            t["pc"].panel_zoom = zoom
        self.layout_tiles()
        if zoom > 1.0:
            self.hsb.pack(side=tk.BOTTOM, fill=tk.X, before=self.tile_scroll)
        else:
            self.hsb.pack_forget()
        new = [float(v) for v in self.tile_scroll.cget("scrollregion").split()]
        if anchor is not None and (spot is None or not (new[2] > 0 and new[3] > 0)):
            self.zoom_note("no spot to keep under the pointer", spot=spot, scrollregion=new, anchor=anchor)
        if spot is not None and new[2] > 0 and new[3] > 0:
            nx, ny = self.spot_now(spot)
            wx, wy = anchor[2], anchor[3]
            self.tile_scroll.xview_moveto(min(max((nx - wx) / new[2], 0.0), 1.0))
            self.tile_scroll.yview_moveto(min(max((ny - wy) / new[3], 0.0), 1.0))
            self.root.update_idletasks()
            # (the view should now be where the spot is under the pointer, unless that is past an edge: if not, say so in zoom_debug.log, and put it right once)
            want_x = min(max(nx - wx, 0.0), max(new[2] - self.tile_scroll.winfo_width(), 0.0))
            want_y = min(max(ny - wy, 0.0), max(new[3] - self.tile_scroll.winfo_height(), 0.0))
            got_x, got_y = self.tile_scroll.canvasx(0), self.tile_scroll.canvasy(0)
            # (moveto lands a few pixels off, by the canvas's border and its rounding: put right to the pixel, or the real layout is seen to jump
            # against the snapshot it replaces)
            for _ in range(2):
                if abs(got_x - want_x) < 0.5 and abs(got_y - want_y) < 0.5:
                    break
                fx, fy = self.tile_scroll.xview()[0], self.tile_scroll.yview()[0]               # (the fractions it is at now, moved by the error)
                self.tile_scroll.xview_moveto(min(max(fx + (want_x - got_x) / new[2], 0.0), 1.0))
                self.tile_scroll.yview_moveto(min(max(fy + (want_y - got_y) / new[3], 0.0), 1.0))
                got_x, got_y = self.tile_scroll.canvasx(0), self.tile_scroll.canvasy(0)
            if abs(got_x - want_x) > 1.5 or abs(got_y - want_y) > 1.5:
                self.trail(f"view {got_x - want_x:+.1f}, {got_y - want_y:+.1f} px off after zoom {zoom:.3g}")
            if abs(got_x - want_x) > 8 or abs(got_y - want_y) > 8:               # (a few pixels is the canvas's own border: moveto is off by it. More is not.)
                self.zoom_note("the view is not where the zoom put it", want=(want_x, want_y), got=(got_x, got_y), spot=spot, spot_now=(nx, ny), anchor=anchor,
                               scrollregion=new, pane=(self.tile_scroll.winfo_width(), self.tile_scroll.winfo_height()), zoom=(self.panel_zoom, zoom))
                self.tile_scroll.xview_moveto(want_x / new[2])
                self.tile_scroll.yview_moveto(want_y / new[3])
            self.zoom_view = (self.tile_scroll.canvasx(0), self.tile_scroll.canvasy(0))
            self.root.after(400, self.check_zoom_view)
        if draw:
            self.place_tiles()
            self.root.after_idle(self.refresh_visible)

    def check_zoom_view(self):
        """Some time after a zoom: has the view moved since, with no drag or wheel? (The sudden jump to the top left: see zoom_debug.log.)"""
        if self.zoom_burst is not None or getattr(self, "panning", False) or getattr(self, "zoom_view", None) is None:
            return
        now = (self.tile_scroll.canvasx(0), self.tile_scroll.canvasy(0))
        if abs(now[0] - self.zoom_view[0]) > 3 or abs(now[1] - self.zoom_view[1]) > 3:
            self.zoom_note("the view moved by itself after a zoom", was=self.zoom_view, now=now, zoom=self.panel_zoom,
                           scrollregion=self.tile_scroll.cget("scrollregion"), pane=(self.tile_scroll.winfo_width(), self.tile_scroll.winfo_height()))
        self.zoom_view = None

    def time_note(self, what, since, limit=0.04):
        """A step of the zoom that took longer than `limit` seconds goes in the trail (and so in zoom_debug.log with the next note)."""
        took = time.perf_counter() - since
        if took > limit:
            self.trail(f"SLOW {what}: {took * 1000:.0f} ms")

    def trail(self, what):
        """The last few things that moved or re-laid the section pane, with the time, kept for zoom_note to write out with a jump."""
        trail = self.__dict__.setdefault("zoom_trail", [])
        trail.append(f"{time.time() % 1000:8.3f} {what}")
        del trail[:-40]

    def note_view(self, axis, first):
        """The pane's scroll position changed (`first`: the fraction at its top or left): a jump to the very start is put in the trail."""
        last = self.__dict__.setdefault("last_view", {"x": 0.0, "y": 0.0})
        now = float(first)
        if now == 0.0 and last[axis] > 0.0:
            self.trail(f"view {axis} went to 0 (was {last[axis]:.3f})")
            if last[axis] > 0.05:                                # (not a step of a drag or a click on the scroll bar: all at once)
                self.zoom_note(f"the view jumped to the start ({axis}, was {last[axis]:.3f})", panning=bool(getattr(self, "panning", False)),
                               zoom=self.panel_zoom, scrollregion=self.tile_scroll.cget("scrollregion"),
                               pane=(self.tile_scroll.winfo_width(), self.tile_scroll.winfo_height()))
        last[axis] = now

    def watch_scroll(self, canvas):
        """(Diagnostic, for the view that jumps to the top left by itself: see zoom_debug.log.) Every call that scrolls the section pane goes into the trail with
        the function that made it, so that a jump to the start says what did it."""
        def wrap(name, scrolls):
            original = getattr(canvas, name)

            def traced(*args, **kw):
                if scrolls(args):
                    caller = sys._getframe(1)
                    self.trail(f"{name}{tuple(args)[:3]} by {caller.f_code.co_name}:{caller.f_lineno}")
                return original(*args, **kw)
            setattr(canvas, name, traced)
        for name in ("xview_moveto", "yview_moveto", "xview_scroll", "yview_scroll"):
            wrap(name, lambda args: True)
        for name in ("xview", "yview"):                         # (also what the scroll bars call; with no arguments they only ask)
            wrap(name, lambda args: bool(args))

    def zoom_note(self, what, **state):
        """A line in zoom_debug.log (next to the cache) for a zoom that did not go as it should: for finding what makes the view jump."""
        state["trail"] = list(self.__dict__.get("zoom_trail", []))
        try:
            with open(core.CACHE_DIR / "zoom_debug.log", "a", encoding="utf-8") as f:
                f.write(time.strftime("%Y-%m-%d %H:%M:%S") + " " + what + " " + repr(state) + chr(10))
        except Exception:
            pass

    ZOOM_SETTLE_MS = 180

    def panel_zoom_at(self, canvas, px, py, factor):
        """The wheel over a section map (px, py in that map) while the whole panel is zoomed: about the pointer."""
        wx = canvas.winfo_rootx() + px - self.tile_scroll.winfo_rootx()
        wy = canvas.winfo_rooty() + py - self.tile_scroll.winfo_rooty()
        self.zoom_step(wx, wy, factor)

    def zoom_step(self, wx, wy, factor):
        """One wheel step. Laying out and drawing the panels again at every step is jerky, so the zoom works on a picture of the pane: the first
        step takes a snapshot of it, each step shows that snapshot, scaled, at once, and when the wheel has been still for ZOOM_SETTLE_MS the panels
        really are laid out and drawn at the new size (set_panel_zoom). Without a snapshot (no screen to take it from) each step does it all."""
        b = self.zoom_burst
        t_step = time.perf_counter()
        if b is None:
            image = self.snapshot_pane()
            self.time_note("snapshot_pane", t_step)
            if image is None:
                self.set_panel_zoom(self.panel_zoom * factor, (self.tile_scroll.canvasx(wx), self.tile_scroll.canvasy(wy), wx, wy))
                return
            b = self.zoom_burst = {"image": image, "start": self.panel_zoom, "target": self.panel_zoom, "at": (wx, wy),
                                   "content": (self.tile_scroll.canvasx(wx), self.tile_scroll.canvasy(wy)), "label": None, "job": None, "photo": None}
            b["spot"] = self.spot_in_tiles(*b["content"])
            self.remember_pane(b, image)
        elif abs(wx - b["at"][0]) > 2 or abs(wy - b["at"][1]) > 2:
            self.reanchor(b, wx, wy)                                 # (the pointer has moved since the zoom began: the zoom goes on about where it is now)
        b["target"] = min(max(b["target"] * factor, 1.0), self.max_panel_zoom())
        if b.get("fill") is not None:                                # (settling was drawing the maps: it stops, to start again when the wheel does)
            self.root.after_cancel(b["fill"])
            b["fill"] = None
        if b.get("show") is None:                                    # (steps that came in while the last was being shown are shown as one: no backlog)
            b["show"] = self.root.after_idle(self.show_zoom_step)
        if b["job"] is not None:
            self.root.after_cancel(b["job"])
        b["job"] = self.root.after(self.ZOOM_SETTLE_MS, self.settle_zoom)

    def spot_in_layout(self, lay, ax, ay):
        """spot_in_tiles for the layout `lay` (see layout_at) instead of the panels as they are: where the point (ax, ay) of the scrolled pane is in them."""
        best = None
        for k in range(len(self.tiles)):
            x, y = 4 + (k % lay["cols"]) * (lay["fw"] + 8), 3 + (k // lay["cols"]) * (lay["fh"] + 6)
            dx = max(x - ax, 0, ax - (x + lay["mw"]))
            dy = max(y - ay, 0, ay - (y + lay["ch"] + lay["mh"]))
            if best is None or dx * dx + dy * dy < best[0]:
                best = (dx * dx + dy * dy, k, ax - x, ay - (y + lay["ch"]), lay["mw"], lay["mh"])
        return best[1:] if best else None

    def reanchor(self, b, wx, wy):
        """The pointer moved while the zoom was going on (the snapshot is up, and the wheel is turned again): the zoom is about the point of the maps that
        is under the pointer now, in the picture as it is shown. The picture does not jump: the view that keeps that point under the pointer is the one
        it is in."""
        lay, view = b.get("shown_layout"), b.get("shown_view")
        if lay is not None and view is not None:
            b["spot"] = self.spot_in_layout(lay, view[0] + wx, view[1] + wy)
        else:                                                        # (nothing shown yet: the panels are as they were when the zoom began)
            b["spot"] = self.spot_in_tiles(self.tile_scroll.canvasx(wx), self.tile_scroll.canvasy(wy))
        b["at"] = (wx, wy)

    def show_zoom_step(self):
        b = self.zoom_burst
        if b is None:
            return
        b["show"] = None
        t_show = time.perf_counter()
        self.show_zoom_snapshot(b)
        self.time_note("show_zoom_snapshot", t_show)

    def pane_look(self):
        """What the pictures of the pane were made with: they are dropped when it changes (a recolor, a size, ...)."""
        return (tuple(id(t["pc"].layers) for t in self.tiles), self.var_dot.get(), self.var_marker.get(), self.var_translucent.get(), self.var_ccf.get(),
                self.var_roi_only.get(), self.var_tile.get(), self.tile_scroll.winfo_width(), getattr(self, "tiles_dark", False))        # (not its height: the scroll bar that comes with the zoom changes that)

    PANE_CACHE_BYTES = 160 * 1024 * 1024                         # the pictures of the pane kept (see remember_pane)
    PANE_CACHE_MAX = 12

    def remember_pane(self, b, image):
        """Keep the picture of the pane at the zoom a wheel run starts from: the pictures kept make a pyramid, from which show_zoom_snapshot takes,
        for each part of the pane, the one that fits the zoom best. A picture replaces an older one only if it is of about the same zoom AND the
        same part of the panel (a picture of somewhere else is not the same picture). The widest (lowest zoom) one is always kept; of the others the
        newest, to PANE_CACHE_MAX pictures and PANE_CACHE_BYTES."""
        look = self.pane_look()
        if look != getattr(self, "pane_cache_look", None):
            self.pane_cache_look, self.pane_cache = look, []
        b["left"], b["top"] = self.tile_scroll.canvasx(0), self.tile_scroll.canvasy(0)
        if getattr(self, "panning", False):                           # (a drag leaves parts of the maps blank until it stops: not a picture to keep)
            return
        entry = {"z": self.panel_zoom, "left": b["left"], "top": b["top"], "image": image}
        pw, ph = image.size
        same = lambda e: (abs(e["z"] / self.panel_zoom - 1.0) <= 0.02 and abs(e["left"] - entry["left"]) < 0.25 * pw
                          and abs(e["top"] - entry["top"]) < 0.25 * ph)
        self.pane_cache = [e for e in self.pane_cache if not same(e)] + [entry]
        size = lambda e: e["image"].size[0] * e["image"].size[1] * 3
        while len(self.pane_cache) > 1 and (len(self.pane_cache) > self.PANE_CACHE_MAX or sum(size(e) for e in self.pane_cache) > self.PANE_CACHE_BYTES):
            widest = min(self.pane_cache, key=lambda e: e["z"])
            # (Many pictures of one zoom, from going in and out there, must not push out those of the other zooms: the oldest of the
            # most crowded half-octave of zoom goes.)
            band = lambda e: round(2 * math.log2(max(e["z"], 1e-6)))
            crowd = {}
            for e in self.pane_cache:
                crowd[band(e)] = crowd.get(band(e), 0) + 1
            lost = [e for e in self.pane_cache if e is not widest and e is not entry]
            if not lost:
                break
            top = max(crowd[band(e)] for e in lost)
            self.pane_cache.remove(next(e for e in lost if crowd[band(e)] == top))

    def blit_window_region(self, x, y, w, h):
        """Just that part of the window, copied from the window's own device context (about 20 ms: PrintWindow draws the whole window, 170 ms). The same
        picture on any monitor. None if it cannot be had, or is blank."""
        try:
            import ctypes
            from ctypes import wintypes as wt
            u, g = ctypes.WinDLL("user32"), ctypes.WinDLL("gdi32")            # (their own copies of the functions: the argument types set here are not the program's)
            vp, ci = ctypes.c_void_p, ctypes.c_int
            u.GetDC.restype, u.GetDC.argtypes = vp, [vp]
            u.ReleaseDC.argtypes = [vp, vp]
            g.CreateCompatibleDC.restype, g.CreateCompatibleDC.argtypes = vp, [vp]
            g.CreateCompatibleBitmap.restype, g.CreateCompatibleBitmap.argtypes = vp, [vp, ci, ci]
            g.SelectObject.argtypes = [vp, vp]
            g.BitBlt.argtypes = [vp, ci, ci, ci, ci, vp, ci, ci, ctypes.c_uint]
            g.GetDIBits.argtypes = [vp, vp, ctypes.c_uint, ctypes.c_uint, vp, vp, ctypes.c_uint]
            g.DeleteObject.argtypes = [vp]
            g.DeleteDC.argtypes = [vp]

            class Header(ctypes.Structure):
                _fields_ = [("size", wt.DWORD), ("w", wt.LONG), ("h", wt.LONG), ("planes", wt.WORD), ("bits", wt.WORD), ("compression", wt.DWORD),
                            ("image", wt.DWORD), ("xppm", wt.LONG), ("yppm", wt.LONG), ("used", wt.DWORD), ("important", wt.DWORD)]
            hwnd = self.root.winfo_id()
            dc = u.GetDC(hwnd)
            if not dc:
                return None
            mem = g.CreateCompatibleDC(dc)
            bmp = g.CreateCompatibleBitmap(dc, w, h)
            try:
                g.SelectObject(mem, bmp)
                if not g.BitBlt(mem, 0, 0, w, h, dc, x, y, 0x00CC0020):
                    return None
                buf = ctypes.create_string_buffer(w * h * 4)
                g.GetDIBits(mem, bmp, 0, h, buf, ctypes.byref(Header(40, w, -h, 1, 32, 0, 0, 0, 0, 0, 0)), 0)
            finally:
                g.DeleteObject(bmp)
                g.DeleteDC(mem)
                u.ReleaseDC(hwnd, dc)
            part = Image.frombuffer("RGBA", (w, h), buf, "raw", "BGRA", 0, 1).convert("RGB")
            return None if part.getextrema() == ((0, 0), (0, 0), (0, 0)) else part
        except Exception:
            return None

    def print_window_region(self, x, y, w, h):
        """A part of this window's own picture of itself (Windows PrintWindow), x, y from the window's top left: right on any monitor, whatever its scaling.
        (A screen grab is not: it comes out black, or the wrong place, for a window on a monitor other than the main one.) None if it cannot be had."""
        if sys.platform != "win32":
            return None
        quick = self.blit_window_region(x, y, w, h)
        if quick is not None:
            return quick
        try:
            import ctypes
            from ctypes import wintypes as wt
            u, g = ctypes.windll.user32, ctypes.windll.gdi32
            u.GetDC.restype = ctypes.c_void_p
            g.CreateCompatibleDC.restype = ctypes.c_void_p
            g.CreateCompatibleBitmap.restype = ctypes.c_void_p
            hwnd = wt.HWND(self.root.winfo_id())
            rect = wt.RECT()
            u.GetClientRect(hwnd, ctypes.byref(rect))
            W, H = rect.right, rect.bottom
            if W <= 0 or H <= 0:
                return None
            screen = u.GetDC(None)
            mem = g.CreateCompatibleDC(ctypes.c_void_p(screen))
            bmp = g.CreateCompatibleBitmap(ctypes.c_void_p(screen), W, H)
            g.SelectObject(ctypes.c_void_p(mem), ctypes.c_void_p(bmp))
            try:
                if not u.PrintWindow(hwnd, ctypes.c_void_p(mem), 2):                  # (2: the full content, the windows inside it too)
                    return None

                class Header(ctypes.Structure):
                    _fields_ = [("size", wt.DWORD), ("w", wt.LONG), ("h", wt.LONG), ("planes", wt.WORD), ("bits", wt.WORD), ("compression", wt.DWORD),
                                ("image", wt.DWORD), ("xppm", wt.LONG), ("yppm", wt.LONG), ("used", wt.DWORD), ("important", wt.DWORD)]
                header = Header(40, W, -H, 1, 32, 0, 0, 0, 0, 0, 0)
                buf = ctypes.create_string_buffer(W * H * 4)
                g.GetDIBits(ctypes.c_void_p(mem), ctypes.c_void_p(bmp), 0, H, buf, ctypes.byref(header), 0)
            finally:
                g.DeleteObject(ctypes.c_void_p(bmp))
                g.DeleteDC(ctypes.c_void_p(mem))
                u.ReleaseDC(None, ctypes.c_void_p(screen))
            whole = Image.frombuffer("RGBA", (W, H), buf, "raw", "BGRA", 0, 1).convert("RGB")
            part = whole.crop((x, y, x + w, y + h))
            if part.size != (w, h) or part.getextrema() == ((0, 0), (0, 0), (0, 0)):         # (nothing in it: not to be trusted)
                return None
            return part
        except Exception:
            return None

    def snapshot_pane(self):
        """A picture of the section pane as it is on screen now (PIL Image), or None if it cannot be taken."""
        try:
            self.root.update_idletasks()
            c = self.tile_scroll
            x, y, w, h = c.winfo_rootx(), c.winfo_rooty(), c.winfo_width(), c.winfo_height()
            image = self.print_window_region(x - self.root.winfo_rootx(), y - self.root.winfo_rooty(), w, h)
            if image is not None:
                return image
            from PIL import ImageGrab
            try:
                image = ImageGrab.grab(bbox=(x, y, x + w, y + h), all_screens=True)         # (all screens: without it only the primary monitor is captured, and a window on another is black)
            except TypeError:
                image = ImageGrab.grab(bbox=(x, y, x + w, y + h))
            image = image.convert("RGB")
            return image if image.size == (w, h) else image.resize((w, h))
        except Exception:
            return None

    def layout_at(self, z):
        """The layout the panels have (or would have) at panel zoom z, as layout_tiles makes it: the number across and down, the slot of one panel (width, height),
        the caption's height, and the size of its map."""
        cols, base = self.tile_params()
        ch = max(t["caption"].winfo_reqheight() for t in self.tiles)
        px = int(base * z)
        ph = int(px / self.tile_aspect)
        border = 2 if z <= 1.0 else 0                                    # (a map of the unzoomed panel has its widget's border)
        mw, mh = px + border, ph + border
        pad_x, pad_y = self.pane_slack()
        return {"cols": cols, "rows": -(-len(self.tiles) // cols), "fw": mw, "fh": ch + mh, "ch": ch, "mw": mw, "mh": mh, "pad_x": pad_x, "pad_y": pad_y}

    def view_at(self, z, spot, at, lay):
        """Where the pane is scrolled to (left, top) at zoom z when the spot under the pointer is to stay under it."""
        k, dx, dy, old_w, old_h = spot
        x = 4 + (k % lay["cols"]) * (lay["fw"] + 8) + grown(dx, old_w, lay["mw"])
        y = 3 + (k // lay["cols"]) * (lay["fh"] + 6) + lay["ch"] + grown(dy, old_h, lay["mh"])
        total_w = 4 + lay["cols"] * (lay["fw"] + 8) + lay["pad_x"]
        total_h = 3 + lay["rows"] * (lay["fh"] + 6) + lay["pad_y"]
        wx, wy = at
        return (min(max(x - wx, 0.0), max(total_w - self.tile_scroll.winfo_width(), 0.0)),
                min(max(y - wy, 0.0), max(total_h - self.tile_scroll.winfo_height(), 0.0)))

    def pane_maps(self, lay_n, left_n, top_n, lay_c, left_c, top_c, w, h):
        """For each column and row of the pane at the layout lay_n (scrolled to left_n, top_n): where in a picture of the pane at the layout lay_c (scrolled
        to left_c, top_c) the same place is. The maps grow and the captions and gaps do not, so it is not one scale for all."""
        xn = left_n + np.arange(w, dtype=np.float64) - 4
        k = np.clip(np.floor(xn / (lay_n["fw"] + 8)), 0, lay_n["cols"] - 1)
        o = xn - k * (lay_n["fw"] + 8)
        oc = np.where(o < 0, o, np.where(o <= lay_n["mw"], o * (lay_c["mw"] / lay_n["mw"]), lay_c["mw"] + (o - lay_n["mw"])))
        mx = 4 + k * (lay_c["fw"] + 8) + oc - left_c
        yn = top_n + np.arange(h, dtype=np.float64) - 3
        r = np.clip(np.floor(yn / (lay_n["fh"] + 6)), 0, lay_n["rows"] - 1)
        o = yn - r * (lay_n["fh"] + 6)
        ch = lay_n["ch"]
        oc = np.where(o < ch, o, np.where(o <= ch + lay_n["mh"], ch + (o - ch) * (lay_c["mh"] / lay_n["mh"]), ch + lay_c["mh"] + (o - ch - lay_n["mh"])))
        my = 3 + r * (lay_c["fh"] + 6) + oc - top_c
        return mx, my

    def pane_picture(self, image, maps, w, h):
        """`image` (a snapshot of the pane) at the places the maps say; also which pixels of the result it covers. (RGB array, boolean array)"""
        mx, my = maps
        src = np.asarray(image)
        sh, sw = src.shape[:2]
        inside = ((my >= 0) & (my <= sh - 1))[:, None] & ((mx >= 0) & (mx <= sw - 1))[None, :]
        out = cv2.remap(src, np.ascontiguousarray(np.broadcast_to(mx.astype(np.float32)[None, :], (h, w))),
                        np.ascontiguousarray(np.broadcast_to(my.astype(np.float32)[:, None], (h, w))), cv2.INTER_LINEAR,
                        borderMode=cv2.BORDER_CONSTANT, borderValue=self.margin_rgb())
        return out, inside

    def pyramid_order(self, target):
        """The kept pictures in the order they are wanted for a pane at panel zoom `target`: those of at least that zoom first (they are made smaller, so
        sharp), the nearest zoom first; then those of a lower zoom (made bigger, so soft), the nearest first."""
        up = sorted((e for e in getattr(self, "pane_cache", []) if e["z"] >= target * 0.98), key=lambda e: e["z"])
        down = sorted((e for e in getattr(self, "pane_cache", []) if e["z"] < target * 0.98), key=lambda e: -e["z"])
        return up + down

    def show_zoom_snapshot(self, b):
        """The snapshot of the pane as it will be at the zoom so far, over the pane: each part of it moved and grown as the real layout will have it (the maps
        grow, the captions and gaps do not), so that what is drawn for real at the end lands where this was; and the edges it does not reach from an earlier
        picture of the pane."""
        t0 = time.perf_counter()
        w, h = b["image"].size
        z = b["target"]
        lay_n = self.layout_at(z)
        left_n, top_n = self.view_at(z, b["spot"], b["at"], lay_n) if b.get("spot") else (b["left"] * z / b["start"], b["top"] * z / b["start"])
        b["shown_layout"], b["shown_view"] = lay_n, (left_n, top_n)          # (what the picture shown is, for reanchor)
        out = np.empty((h, w, 3), np.uint8)
        out[:] = self.margin_rgb()
        todo = np.ones((h, w), bool)                                  # (the parts of the pane no picture has been put on yet)
        # The picture the zoom started from is always first: it is the pane as it is. The older ones only fill in what it does not reach
        # (they may be out of date, or have blank parts from a drag: the pane is not drawn all over while one goes on).
        first = {"z": b["start"], "left": b["left"], "top": b["top"], "image": b["image"]}
        others = [e for e in self.pyramid_order(z) if e["image"] is not b["image"]]
        widest = min(others, key=lambda e: e["z"]) if others else None
        used = 0
        for e in [first] + [e for e in others if e is not widest] + ([widest] if widest is not None else []):
            if not todo.any():
                break
            if used >= 4 and e is not widest:                        # (the cost of a step: the widest picture, which reaches furthest, is not left out)
                continue
            mx, my = self.pane_maps(lay_n, left_n, top_n, self.layout_at(e["z"]), e["left"], e["top"], w, h)
            sh, sw = e["image"].size[1], e["image"].size[0]
            ix, iy = (mx >= 0) & (mx <= sw - 1), (my >= 0) & (my <= sh - 1)
            if not ix.any() or not iy.any():
                continue
            reach = iy[:, None] & ix[None, :] & todo
            if not reach.any():
                continue                                             # (nothing new: not counted)
            picture, _ = self.pane_picture(e["image"], (mx, my), w, h)
            np.copyto(out, picture, where=reach[:, :, None])
            todo &= ~reach
            used += 1
        t1 = time.perf_counter()
        fresh = b.get("photo") is None or (b["photo"].width(), b["photo"].height()) != (w, h)
        if fresh:
            b["photo"] = ImageTk.PhotoImage(Image.fromarray(out), master=self.tile_scroll)
        else:
            b["photo"].paste(Image.fromarray(out))                   # (the label shows the same photo, now with the new picture in it: nothing else to do)
        t2 = time.perf_counter()
        if b["label"] is None:
            b["label"] = tk.Label(self.tile_scroll, bd=0, highlightthickness=0, cursor="fleur")
            pointer = lambda e: (e.x_root - self.tile_scroll.winfo_rootx(), e.y_root - self.tile_scroll.winfo_rooty())       # (where it is now, not where the zoom began)
            b["label"].bind("<MouseWheel>", lambda e: self.zoom_step(*pointer(e), 1.25 if e.delta > 0 else 0.8))
            b["label"].bind("<Button-4>", lambda e: self.zoom_step(*pointer(e), 1.25))
            b["label"].bind("<Button-5>", lambda e: self.zoom_step(*pointer(e), 0.8))
            b["label"].place(x=0, y=0, width=w, height=h)
        if fresh:
            b["label"].configure(image=b["photo"])
            b["label"].lift()
        t3 = time.perf_counter()
        b["label"].update_idletasks()
        if time.perf_counter() - t0 > 0.04:
            self.time_note(f"snapshot step: maps+remap {1000 * (t1 - t0):.0f}, photo {1000 * (t2 - t1):.0f}, configure {1000 * (t3 - t2):.0f}, "
                           f"update_idletasks {1000 * (time.perf_counter() - t3):.0f} ms; total", t0, 0)

    def settle_zoom(self):
        """The wheel has been still: the panels are laid out and drawn for real at the zoom reached, and the snapshot taken away."""
        b = self.zoom_burst
        if b is None:
            return
        b["job"] = None
        wx, wy = b["at"]
        if b.get("spot") is None:
            self.set_panel_zoom(b["target"], (b["content"][0], b["content"][1], wx, wy))
            self.root.update_idletasks()
            self.refresh_visible()                                  # (what is in view is drawn before the snapshot goes, so that there is no blank between)
            self.root.update_idletasks()
            self.cancel_zoom_burst()
            return
        # The panels are laid out now (quick), and drawn one map at a time with the snapshot still over them, so that the wheel
        # starting again is seen to between the maps: it goes on from the snapshot (which does not depend on the layout under it), and the
        # maps drawn so far are laid out again at the new zoom when it stops. (zoom_fill)
        t_lay = time.perf_counter()
        b["laid"] = True                                            # (before the layout: its idle work must not draw everything, see refresh_visible)
        self.set_panel_zoom(b["target"], (None, None, wx, wy), spot=b["spot"], draw=False)
        self.time_note("settle: set_panel_zoom (layout)", t_lay, 0.02)
        b["settled_at"] = time.perf_counter()
        self.trail("settle: laid out, filling")
        self.zoom_fill()

    def zoom_fill(self):
        """Settling a zoom: one map drawn where it now is (those with the most on screen first), then the next after the events that came in between;
        when all are drawn the snapshot goes."""
        b = self.zoom_burst
        if b is None:
            return
        b["fill"] = None
        if b["job"] is not None:                                    # (the wheel has started again: settle_zoom comes when it stops)
            return
        busy = b.get("busy", 0)                                     # (maps being drawn in their threads: the end of each comes back here)
        if busy >= ZOOM_FILL_PARALLEL:
            return
        t_fill = time.perf_counter()
        self.started_in_background = 0

        def started():
            b["busy"] = b.get("busy", 0) + 1                        # (counted before the draw can end, and so call map_drawn)
            b["maps"] = b.get("maps", 0) + 1
        more = self.place_tiles(max_draws=ZOOM_FILL_PARALLEL - busy, background=lambda: self.map_drawn(b), started=started)
        self.time_note(f"zoom_fill maps (more={more})", t_fill)
        if self.started_in_background or b.get("busy", 0):          # (some are being drawn: the next are started as they end)
            return
        if more:
            b["fill"] = self.root.after(PAN_FILL_GAP_MS, self.zoom_fill)
            return
        t0 = time.perf_counter()
        for t in self.tiles:                                        # (maps off screen when recolored, and maps that are waiting for a draw)
            if t.get("dirty") and self.frame_visible(t["frame"]):
                self.build_tile(t)
            t["pc"].refresh()
            if time.perf_counter() - t0 > 0.03 and t is not self.tiles[-1]:
                b["fill"] = self.root.after(PAN_FILL_GAP_MS, self.zoom_fill)
                return
        t_idle = time.perf_counter()
        self.root.update_idletasks()
        self.time_note("settle: final update_idletasks", t_idle, 0.02)
        if any(t["pc"].pending is not None for t in self.tiles):    # (those draws come at idle: the snapshot goes after them)
            b["fill"] = self.root.after(PAN_FILL_GAP_MS, self.zoom_fill)
            return
        took = time.perf_counter() - b.get("settled_at", time.perf_counter())
        if took > 2.0:                                              # (a slow settle: said in zoom_debug.log, with what it was spent on)
            try:
                with open(core.CACHE_DIR / "zoom_debug.log", "a", encoding="utf-8") as f:
                    f.write(time.strftime("%H:%M:%S") + f" SLOW SETTLE {took:.1f}s: zoom {b['start']:.3g} -> {b['target']:.3g}, {b.get('maps', 0)} maps drawn, "
                            f"{len(self.tiles)} sections, pane {self.tile_scroll.winfo_width()}x{self.tile_scroll.winfo_height()}, mode {self.mode}" + chr(10))
            except Exception:
                pass
        self.time_note("settle: whole (laid out to snapshot gone)", b.get("settled_at", time.perf_counter()), 0.0)
        self.cancel_zoom_burst()
        self.trail("settle: done, snapshot gone")

    def map_drawn(self, b):
        """A map drawn in the background by zoom_fill is on the screen (or dropped): on to the next, unless the wheel has started again, in which case
        settle_zoom goes on from here when it stops."""
        b["busy"] = max(b.get("busy", 0) - 1, 0)
        if self.zoom_burst is b and b["job"] is None and b.get("fill") is None:
            b["fill"] = self.root.after(PAN_FILL_GAP_MS, self.zoom_fill)

    def cancel_zoom_burst(self):
        b, self.zoom_burst = self.zoom_burst, None
        if b is not None:
            for name in ("job", "fill", "show"):
                if b.get(name) is not None:
                    self.root.after_cancel(b[name])
            if b["label"] is not None:
                b["label"].destroy()

    def panel_zoom_reset(self):
        self.trail("panel_zoom_reset: a double-click on a map put the panel back")
        self.cancel_zoom_burst()
        self.set_panel_zoom(1.0)

    def bind_pane_gestures(self, widget):
        """The wheel, a drag and a double-click over `widget` (the pane itself, a panel's frame or its title: not a section map, which has its
        own) act on the whole panel while Zoom whole panel is on. The wheel over the pane itself scrolls it as it always did."""
        widget.bind("<MouseWheel>", self.on_pane_wheel, add="+")
        for press, motion in (("<ButtonPress-1>", "<B1-Motion>"), ("<ButtonPress-2>", "<B2-Motion>"), ("<ButtonPress-3>", "<B3-Motion>")):
            widget.bind(press.replace("Press", "Release"), lambda e: self.drag_state(False), add="+")
            widget.bind(press, lambda e: self.var_panel_zoom.get() and self.panel_pan_start(e), add="+")
            widget.bind(motion, lambda e: self.var_panel_zoom.get() and self.panel_pan_move(e), add="+")
        # (No double-click here: Tk gives a quick second press to a Double-Button-1 binding instead of to the press binding above, so a drag that
        # started within a double-click's time of an earlier click would not have moved the pane. The section maps, which have their own, reset it.)

    def on_pane_wheel(self, e):
        if self.var_panel_zoom.get():
            self.zoom_step(e.x_root - self.tile_scroll.winfo_rootx(), e.y_root - self.tile_scroll.winfo_rooty(), 1.25 if e.delta > 0 else 0.8)
        elif e.widget is self.tile_scroll:
            self.tile_scroll.yview_scroll(-1 if e.delta > 0 else 1, "units")

    def panel_pan_start(self, e):
        self.drag_state(True)
        self.cancel_zoom_burst()
        self.trail("pan: button pressed on the pane")
        self.tile_scroll.scan_mark(e.x_root - self.tile_scroll.winfo_rootx(), e.y_root - self.tile_scroll.winfo_rooty())

    def panel_pan_move(self, e):
        """A drag of the zoomed panel moves it at once: the maps whose pictures no longer cover the screen are not all drawn again there and
        then (which held the drag up as it went past the margin they are drawn with), but one at a time between the mouse's moves (pan_fill),
        the new parts blank meanwhile; and everything once the drag has been still for PANEL_PAN_SETTLE_MS (pan_settle)."""
        if not getattr(self, "panning", False):
            self.trail("pan: drag begins")
        self.panning = True
        self.tile_scroll.scan_dragto(e.x_root - self.tile_scroll.winfo_rootx(), e.y_root - self.tile_scroll.winfo_rooty(), gain=1)
        self.root.after_idle(self.refresh_visible)
        if getattr(self, "pan_settle_job", None) is not None:
            self.root.after_cancel(self.pan_settle_job)
        self.pan_settle_job = self.root.after(PANEL_PAN_SETTLE_MS, self.pan_settle)
        if getattr(self, "pan_fill_job", None) is None:
            self.pan_fill_job = self.root.after(PAN_DRAW_MS, self.pan_fill)

    def pan_fill(self):
        """During a drag of the zoomed panel: one map drawn where it now is (the one with the most on screen), then the mouse's moves are seen
        to before the next."""
        self.pan_fill_job = None
        if not getattr(self, "panning", False):
            return
        if self.place_tiles(max_draws=1):
            self.pan_fill_job = self.root.after(PAN_FILL_GAP_MS, self.pan_fill)

    def pan_settle(self):
        """The drag has been still: everything in view drawn, as without one."""
        self.end_pan()
        self.refresh_visible()

    def end_pan(self):
        self.panning = False
        for name in ("pan_settle_job", "pan_fill_job"):
            job = getattr(self, name, None)
            if job is not None:
                self.root.after_cancel(job)
                setattr(self, name, None)

    def apply_dot_size(self):
        """The Cell size slider: redraw the section panels with bigger or smaller dots (the ones off screen when they come on screen)."""
        f = self.dot_factors[self.var_dot.get()]
        for t in self.tiles:
            t["pc"].dot_factor = f
            if self.frame_visible(t["frame"]):
                t["pc"].request_draw()
            else:
                t["pc"].stale = True

    def apply_opacity(self):
        """The Translucent box: the UMAP and the section panels are drawn again with dots that are see-through (or not)."""
        f = TRANSLUCENT if self.var_translucent.get() else 1.0
        self.umap.opacity = f
        self.umap.request_draw()
        for t in self.tiles:
            t["pc"].opacity = f
            if self.frame_visible(t["frame"]):
                t["pc"].request_draw()
            else:
                t["pc"].stale = True

    def apply_marker_size(self):
        """The Marker size slider: the hover ring and the + marks of the section panels, bigger or smaller."""
        f = MARKER_SIZES[self.var_marker.get()]
        for t in self.tiles:
            pc = t["pc"]
            pc.marker_factor = f
            for m in pc.markers.values():
                pc.place_marker(m)
            if pc.pluses is not None:
                pc.request_compose()

    def link_tiles(self, source, zoom, dx, dy):
        """Zooming or panning one section does the same to all the others."""
        for t in self.tiles:
            if t is not source:
                t["pc"].follow(zoom, dx, dy)

    def section_tile(self, label):
        for t in self.tiles:
            if t["label"] == label:
                return t
        return None

    # ---- colouring

    def on_mode(self):
        if self.busy:
            self.var_mode.set(self.mode)
            return
        mode = self.var_mode.get()
        if mode == "imputed" and self.imputed_state["adata"] is None:
            ok = bool(self.ensure_imputed(self.imputed_state, self.abc_cache)) if self.ensure_imputed else False
            if not ok:
                self.var_mode.set(self.mode)
                return
        self.mode = mode
        self.update_frames()
        self.update_group_buttons()
        self.recolor()

    def on_level(self):
        if self.busy:
            self.var_level.set(self.level)
            return
        old = self.level
        self.level = self.var_level.get()
        if self.level != old:                                          # (each level keeps its own group numbers: those last typed for it, else none)
            self.ids_by_level[old] = self.var_ids.get()
            self.var_ids.set(self.ids_by_level.get(self.level, ""))
            self.update_frames()
        if self.mode in ("all", "single"):
            self.recolor()

    def gene_rgb(self, t):
        """The colour of a gene's value (t in 0..1) in the theme chosen."""
        return light_rgb(t) if self.var_theme.get() == "light" else viridis_rgb(t)

    def on_theme(self):
        """Dark or Light: the gene shown is drawn again in the other theme (nothing is read again)."""
        if self.busy:
            self.var_theme.set("light" if self.var_theme.get() == "dark" else "dark")
            return
        self.theme_choice = self.var_theme.get()          # (chosen: from now on it is kept, for one gene and for several)
        if self.mode in ("gene", "imputed") and self.gene_cache is not None:
            key = self.gene_cache[0]
            self.show_genes(list(key[0]), key[1], key[2], self.gene_cache[1])

    def on_scale(self):
        if self.mode in ("gene", "imputed") and self.var_gene.get().strip():
            self.recolor()

    def recolor(self, keep_off=False):
        if not keep_off:
            self.last_group_click = None
            self.off_groups = set()                      # (a new level, mode or list of groups starts with every group colored)
        self.var_status.set("")
        self.legend_redraw = None                        # (the color bar or gene legend is drawn again only while there is one: see cbar_resized)
        self.cell_marks = None
        self.select_bar.pack_forget()
        self.legend.pack_forget()
        self.group_list.pack_forget()
        if self.mode == "all":
            self.color_all()
        elif self.mode == "single":
            self.color_single()
        else:
            self.color_gene()
            self.group_list.pack_forget()                # (a gene's colors say nothing of the groups: no list, even before a gene is entered and the cells keep their group colors)
            self.legend.pack_forget()
        self.sync_select_bar()

    def dot_scale(self):
        return self.host.get("umap_circle_sizes", {}).get(self.viewer_settings.get("umap_circle_size", "Default"), 1.0)

    def apply(self, rgb_run, order_run, bg_rgb_fn, labels=None, extra_umap_layers=(), dark=False):
        """Show the colouring. `rgb_run`: (n, 3) colour of each run cell; `order_run`: the order to draw them in (lowest first);
        `bg_rgb_fn(section label)` -> (x, y, rgb) of the section panel's background cells, or None."""
        self.rgb_run = rgb_run
        o = order_run
        self.run_cell_rank = self.bg_cell_rank = None   # (the mode that has groups says them: see set_cell_ranks)
        self.tiles_dark = dark
        self.umap_order = o                           # (the order of the UMAP's cells in its first layer: see set_umap_marks)
        self.umap.set_background(self.view_background())
        self.set_tile_background()
        self.umap.set_layers([{"x": self.coords[o, 0], "y": self.coords[o, 1], "rgb": rgb_run[o], "size": 4 * self.dot_scale()}] + list(extra_umap_layers))
        self.umap.clear_labels()
        self.label_info = list(labels or [])
        for text, x, y, color, weight in labels or []:
            self.umap.add_label(text, x, y, color, weight)
        self.umap.reposition()
        self.tile_queue = list(self.tiles)
        self.tile_fn = bg_rgb_fn
        if not getattr(self, "tile_job", None):
            self.tile_job = self.root.after(1, self.tile_step)

    def shape_codes(self, ranks):
        """The marker code of each rank (the browser's rule: the color cycle is 20 ranks, a new shape every cycle, hollow shapes after 200); a
        negative rank (a cell in no listed group) is a plain disc. Made from the browser's own category_rank_shape, so that the screen and the
        saved UMAP agree."""
        table = getattr(self, "shape_table", None)
        if table is None:
            fn = self.host.get("category_rank_shape")
            table = self.shape_table = np.array(
                [(SHAPE_MARKERS.index(m) + (10 if is_open else 0)) for m, is_open in (fn(r) for r in range(2 * 10 * 20 * 2))] if fn else [0], np.int64)
        ranks = np.asarray(ranks)
        return np.where(ranks >= 0, table[np.maximum(ranks, 0) % len(table)], 0)

    def set_umap_marks(self, rank):
        """Give the UMAP's cells the shape of their group (`rank`: each cell's group rank, or -1; in the cells' own order) and draw it again."""
        lay = self.umap.layers[0] if self.umap.layers else None
        order = getattr(self, "umap_order", None)
        if lay is None or order is None or len(rank) != len(order):
            return
        lay["marks"] = self.shape_codes(rank)[order].astype(np.int64)
        self.umap.request_draw()

    def set_cell_ranks(self, run_rank, bg_rank):
        """The group rank of each cell of the run (`run_rank`, -1 for none) and of each cell of the sections' backgrounds (`bg_rank`, or None): the
        UMAP's cells take their group's shape now, and the section panels' when they are made (just after: see apply)."""
        self.run_cell_rank, self.bg_cell_rank = run_rank, bg_rank
        self.set_umap_marks(run_rank)

    def margin_rgb(self):
        """The color of the blank margins a zoom out shows where no picture of the pane reaches (see show_zoom_snapshot): the dark of the plots in the dark
        theme, light gray (the pane's own color) otherwise: dark plots are not shown in a sea of white. The pane itself, between the panels, stays light."""
        return MULTI_BG if getattr(self, "tiles_dark", False) else (244, 244, 244)

    def set_tile_background(self):
        """The section panels' background: the UMAP's (see view_background), and white in ROI only (whose gray cells are for it)."""
        bgc = (255, 255, 255) if self.var_roi_only.get() else self.view_background()
        for t in getattr(self, "tiles", []):
            t["pc"].set_background(bgc)

    def apply_roi_only(self):
        """The ROI only box was ticked or cleared: the panels are made again (those in view now, the others as they come into view)."""
        self.set_tile_background()
        for t in getattr(self, "tiles", []):
            t["dirty"] = True
        self.refresh_visible()

    def frame_visible(self, fr):
        """Whether a section panel is (at least partly) in the scrolled pane."""
        slot = getattr(fr, "slot_wh", None)                          # (while the panel is zoomed: the room the map takes, of which its widget is only the part in view)
        top, h = getattr(fr, "tile_y", fr.winfo_y()), slot[1] if slot else fr.winfo_reqheight()
        vt, vh = self.tile_scroll.canvasy(0), self.tile_scroll.winfo_height()
        if not (top + h > vt and top < vt + vh):
            return False
        left = getattr(fr, "tile_x", None)                          # (sideways too, while the whole panel is zoomed in)
        if left is None:
            return True
        vl, vw = self.tile_scroll.canvasx(0), self.tile_scroll.winfo_width()
        return left + (slot[0] if slot else fr.winfo_reqwidth()) > vl and left < vl + vw

    def refresh_visible(self):
        """Panels that were off screen when they were recolored or moved are made now that they are on screen."""
        b = self.zoom_burst
        if b is not None and b.get("laid"):                    # (a zoom settling: its maps are drawn one at a time, see zoom_fill)
            if b.get("fill") is None and b["job"] is None:
                b["fill"] = self.root.after(PAN_FILL_GAP_MS, self.zoom_fill)
            return
        self.place_tiles()
        if getattr(self, "panning", False):          # (while the panel is dragged: drawn bit by bit, by pan_fill, and all of it once the drag stops)
            return
        for t in self.tiles:
            if t.get("dirty") and self.frame_visible(t["frame"]):
                self.build_tile(t)
            t["pc"].refresh()

    def tile_step(self):
        self.tile_job = None
        t0 = time.perf_counter()
        while self.tile_queue and time.perf_counter() - t0 < 0.03:
            t = self.tile_queue.pop(0)
            t["dirty"] = True
            if self.frame_visible(t["frame"]):
                self.build_tile(t)
        if self.tile_queue:
            self.tile_job = self.root.after(1, self.tile_step)

    def build_tile(self, t):
        """The section's dots, back to front: the cells of the other kind (half size), then the gray cells of the run's kind (those of the
        whole section, then the run's own), and last all the colored ones, so a colored cell is never hidden by a gray one."""
        t["dirty"] = False
        idx = self.run_idx[t["label"]]
        layers = []
        bgl = self.tile_fn(t["label"])
        roi_only = self.var_roi_only.get()
        shapes = True                    # (the cells of a section map have the shape of their group, as in the UMAP; plain discs while the view moves)
        sources = []                     # (x, y, rgb, marks) of the run's kind, the section's cells first (marks: a marker code for each cell, or None)
        if bgl is not None:
            x, y, rgb, w = bgl
            bg_rank = getattr(self, "bg_cell_rank", None)
            mk = self.shape_codes(bg_rank[w]).astype(np.int64) if shapes and bg_rank is not None and not roi_only else None
            if roi_only:                 # (the cells of the whole section go gray: the run's own cells, drawn after them, are the ROIs)
                rgb = np.tile(LIGHT_GREY, (len(x), 1))
            if self.bg_under is not None and self.bg["levels"].get("class") is not None:
                # the cells of the other kind (non-neurons in a neuron run) go underneath, so they never cover the ones being looked at
                ids = core.extract_leading_numeric_id(pd.Series(self.bg["levels"]["class"][w])).to_numpy()
                non_neuron = np.isin(ids, list(core.NON_NEURON_CLASS_IDS))
                other = non_neuron if self.bg_under == "nonneuron" else ~non_neuron
                ox, oy, orgb = x[other], y[other], rgb[other]
                omk = mk[other] if mk is not None else None
                ogray = (orgb == GREY).all(axis=1) | (orgb == LIGHT_GREY).all(axis=1)
                for pick, is_gray in ((ogray, True), (~ogray, False)):
                    layer = {"x": ox[pick], "y": oy[pick], "rgb": orgb[pick], "size_mm": DOT_MM * OTHER_CELLS_DOT_RATIO, "gray": is_gray}
                    if omk is not None:
                        layer["marks"] = omk[pick]
                    layers.append(layer)
                x, y, rgb = x[~other], y[~other], rgb[~other]
                mk = mk[~other] if mk is not None else None
            sources.append((x, y, rgb, mk))
        # the run's own cells of this section, in drawing order (lowest first)
        sel = idx[np.argsort(self.rank_run[idx], kind="stable")]
        run_rank = getattr(self, "run_cell_rank", None)
        sources.append((self.cx[sel], self.cy[sel], self.rgb_run[sel],
                        self.shape_codes(run_rank[sel]).astype(np.int64) if shapes and run_rank is not None else None))
        grays = []
        for x, y, rgb, mk in sources:
            gray = (rgb == GREY).all(axis=1) | (rgb == LIGHT_GREY).all(axis=1)
            grays.append((x, y, rgb, gray, mk))
        for x, y, rgb, gray, mk in grays:
            layer = {"x": x[gray], "y": y[gray], "rgb": rgb[gray], "size_mm": DOT_MM, "gray": True}      # ("gray": saved as pictures in an SVG)
            if mk is not None:
                layer["marks"] = mk[gray]
            layers.append(layer)
        for x, y, rgb, gray, mk in grays:
            layer = {"x": x[~gray], "y": y[~gray], "rgb": rgb[~gray], "size_mm": DOT_MM}
            if mk is not None:
                layer["marks"] = mk[~gray]
            layers.append(layer)
        t["pc"].set_layers(layers)

    def color_all(self):
        level = self.level
        if self.info.get(level) is None:
            self.var_status.set(f"No '{level}' column on this dataset.")
            return
        codes, uniques = self.category_codes(level)
        if len(uniques) == 0:                                  # (no cell has a value at this level, such as the neurotransmitter of a run of non-neurons)
            self.var_status.set(f"No cell of this dataset has a {self.level_names[level].lower()}.")
            return
        counts = np.bincount(codes[codes >= 0], minlength=len(uniques))
        ranked = np.argsort(-counts, kind="stable")[:MAX_COLORED_CATEGORIES]
        lut = np.tile(GREY, (len(uniques), 1))
        for r, u in enumerate(ranked):
            lut[u] = GREY if str(uniques[u]) in self.off_groups else rank_color(r)
        rgb = np.where((codes >= 0)[:, None], lut[np.maximum(codes, 0)], LIGHT_GREY).astype(np.uint8)         # (code -1: a cell with no annotation at this level: gray)
        self.rank_run = np.zeros(self.n, np.int32)
        bgc = self.bg_codes(level) if self.bg else None

        def bg_fn(label):
            if not self.bg:
                return None
            w = self.bg["where"].get(label)
            if w is None:
                return None
            c = bgc[w] if bgc is not None else np.full(len(w), -1)
            col = np.where((c >= 0)[:, None], lut[np.maximum(c, 0)], LIGHT_GREY)
            return self.bg["x"][w], self.bg["y"][w], col.astype(np.uint8), w

        # each coloured group's number at its centre
        cent = pd.DataFrame({"g": codes, "x": self.coords[:, 0], "y": self.coords[:, 1]})
        on = [u for u in ranked if str(uniques[u]) not in self.off_groups]
        cent = cent[cent["g"].isin(on)].groupby("g").median()
        labels = []
        weights = label_weights(counts[cent.index.to_numpy()])
        for (g, row), weight in zip(cent.iterrows(), weights):
            m = re.match(r"^\s*(\d+)", str(uniques[g]))
            dark = (lut[g].astype(float) * 0.55).astype(int)               # (the group's own color, darker)
            labels.append((m.group(1) if m else str(uniques[g]), row["x"], row["y"], "#%02x%02x%02x" % tuple(dark), float(weight)))
        self.apply(rgb, np.arange(self.n), bg_fn, labels)
        self.view_kind, self.view_token = "all", f"all_{level}"
        self.view_title = f"All {plural_name(self.level_names[level])}"
        rank_of = np.full(len(uniques), -1)
        rank_of[ranked] = np.arange(len(ranked))
        self.export_cat = {"rank": np.where(codes >= 0, rank_of[np.maximum(codes, 0)], -1), "title": self.level_names[level],
                           "entries": [(hexcolor(GREY if str(uniques[u]) in self.off_groups else rank_color(r)), str(uniques[u]), int(counts[u])) for r, u in enumerate(ranked)]}
        self.set_cell_ranks(self.export_cat["rank"], np.where(bgc >= 0, rank_of[np.maximum(bgc, 0)], -1) if bgc is not None else None)
        self.shown_values = None
        n_cat = int((counts > 0).sum())
        self.var_status.set(f"Showing all {n_cat} {plural_name(self.level_names[level])}" + (f"; the {len(ranked)} largest colored, the rest gray." if n_cat > len(ranked) else ".") + " Click group in list to select/deselect.")
        self.values = None
        self.draw_group_list(ranked, uniques, counts)

    def draw_legend_rows(self, canvas, rows, x0):
        """The rows of a group legend, in columns: swatch, number, cell count (right-aligned) and the group's name (all names start at the same x).
        `rows`: (swatch color, number, count text, name, grayed out). Returns the width used."""
        font = tkfont.Font(root=canvas, family="Segoe UI", size=9)
        canvas.legend_fonts = (font, tkfont.Font(root=canvas, family="Segoe UI", size=9, weight="bold"))      # (kept: the bold one is for the row of the hovered cell)
        self.bold_row = None
        w_num = max((font.measure(r[1]) for r in rows), default=0)
        w_cnt = max((font.measure(r[2]) for r in rows), default=0)
        x_num, x_cnt, x_name = x0 + 22, x0 + 22 + w_num + 8 + w_cnt, x0 + 22 + w_num + 8 + w_cnt + 10
        self.legend_cols = (x_num, x_cnt - w_cnt, x_cnt, x_name)          # (where the number, count and name columns start; the count's right edge)
        for k, (col, num, cnt, name, off, *more) in enumerate(rows):
            cy = 4 + 22 * k + 11
            fill = "gray55" if off else "black"
            draw_marker(canvas, x0 + 6, cy, more[0] if more else 0, col)                               # (the shape the group's cells have in the UMAP)
            canvas.create_text(x_num, cy, anchor="w", text=num, font=font, fill=fill, tags=f"row{k}")
            canvas.create_text(x_cnt, cy, anchor="e", text=cnt, font=font, fill=fill, tags=f"row{k}")
            canvas.create_text(x_name, cy, anchor="w", text=name, font=font, fill=fill, tags=f"row{k}")
        return x_name + max((font.measure(r[3]) for r in rows), default=0) + 8

    def draw_group_list(self, ranked, uniques, counts):
        """The colored groups, each with its swatch, number, cell count and name (the ones past the colored cap are gray and not listed),
        in the order of the Sort by buttons."""
        self.group_list_args = (ranked, uniques, counts)
        c = self.group_canvas
        c.delete("all")
        rows = []
        for r, u in enumerate(ranked):
            m = re.match(r"^\s*(\d+)", str(uniques[u]))
            off = str(uniques[u]) in self.off_groups
            name = re.sub(r"^\s*\d+\s*", "", str(uniques[u])) if m else ""
            rows.append({"id": int(m.group(1)) if m else 10 ** 9, "color": hexcolor(GREY if off else rank_color(r)),
                         "num": m.group(1) if m else str(uniques[u]), "count": int(counts[u]), "name": name, "off": off, "key": str(uniques[u]),
                         "code": int(self.shape_codes(np.array([r]))[0])})
        column, descending = self.group_sort
        rows.sort(key=lambda r: (r["id"], r["key"]))                      # (ties stay in ID order)
        if column == "count":
            rows.sort(key=lambda r: r["count"], reverse=descending)
        elif column == "name":
            rows.sort(key=lambda r: (r["name"] or r["num"]).lower(), reverse=descending)
        elif descending:
            rows.reverse()
        self.group_rows = [r["key"] for r in rows]
        width = self.draw_legend_rows(c, [(r["color"], r["num"], f"({r['count']})", r["name"], r["off"], r["code"]) for r in rows], 6)
        c.configure(scrollregion=(0, 0, max(300, width), 8 + 22 * len(rows)))
        self.draw_group_header()
        self.group_list.pack(fill=tk.BOTH, expand=True, pady=(6, 0), after=self.box)
        self.sync_select_bar()                              # (packed again, the list is above the buttons: they go back above it)

    def draw_group_header(self):
        """The headings of the group list's columns, over the columns, each with an arrow: dark on the column it is sorted by (pointing the way it
        is sorted), light gray on the others (pointing the way a click on them would sort)."""
        h = self.group_header
        h.delete("all")
        x_num, x_cnt_left, x_cnt, x_name = self.legend_cols
        column, descending = self.group_sort
        font = tkfont.Font(root=h, family="Segoe UI", size=9, weight="bold")
        gap = font.measure(" ")
        for key, text, x, right in (("id", "#", x_num, False), ("count", "Cells", x_cnt, True), ("name", "Name", x_name, False)):
            active = key == column
            down = descending if active else key in GROUP_SORT_DESCENDING_FIRST
            arrow = "\u25bc" if down else "\u25b2"
            color = "gray20" if active else "gray75"
            w_arrow = font.measure(arrow)
            if right:                                     # (the text ends at the column's right edge, with the arrow there)
                h.create_text(x, 10, anchor="e", text=arrow, font=font, fill=color)
                h.create_text(x - w_arrow - gap, 10, anchor="e", text=text, font=font, fill="gray30")
                self.count_head_left = min(x_cnt_left, x - w_arrow - gap - font.measure(text)) - 2      # (where a click counts for Cells, not #)
            else:
                h.create_text(x, 10, anchor="w", text=text, font=font, fill="gray30")
                h.create_text(x + font.measure(text) + gap, 10, anchor="w", text=arrow, font=font, fill=color)

    def on_group_header_click(self, e):
        """Click a heading to sort by that column; click it again to reverse the order."""
        if self.mode != "all" or getattr(self, "group_list_args", None) is None:
            return
        x_num, x_cnt_left, x_cnt, x_name = self.legend_cols
        key = "id" if e.x < getattr(self, "count_head_left", x_cnt_left) else "count" if e.x < x_name else "name"
        column, descending = self.group_sort
        self.group_sort = (key, not descending) if key == column else (key, key in GROUP_SORT_DESCENDING_FIRST)
        self.draw_group_list(*self.group_list_args)
        self.group_canvas.yview_moveto(0)

    def group_keys(self, level, values):
        """What each of `values` (a level's group labels) is called when typed in Specific group(s): its number, or, at a level whose groups have names and no
        numbers (neurotransmitter), the name itself."""
        if self.level_has_ids.get(level, True):
            return core.extract_leading_numeric_id(pd.Series(values)).to_numpy()
        return np.asarray([str(v) for v in values], dtype=object)

    def typed_groups(self, level, uniques):
        """(the groups typed in the box, the entries that name none): numbers at a level that has them, else names (see parse_names)."""
        return self.parse_ids() if self.level_has_ids.get(level, True) else self.parse_names(uniques)

    def parse_names(self, uniques):
        """The group names typed in the box, separated by commas, as the level's own names: whatever the case, and the start of a name will do if it is of
        one group alone ('Glut' is Glut, though Glut-GABA also starts with it: the whole name wins). (names, the entries that name no group or several)"""
        names = [str(u) for u in uniques]
        folded = [n.casefold() for n in names]
        found, bad = [], []
        for tok in (t.strip() for t in self.var_ids.get().split(",")):
            if not tok:
                continue
            key = tok.casefold()
            hits = [i for i, f in enumerate(folded) if f == key] or [i for i, f in enumerate(folded) if f.startswith(key)]
            if len(hits) == 1:
                if names[hits[0]] not in found:
                    found.append(names[hits[0]])
            else:
                bad.append(tok)
        return found, bad

    def parse_ids(self):
        ids, bad = [], []
        for tok in re.split(r"[,\s]+", self.var_ids.get().strip()):
            if not tok:
                continue
            try:
                v = int(tok)
            except ValueError:
                bad.append(tok)
                continue
            if v not in ids:
                ids.append(v)
        return ids, bad

    def color_single(self):
        level = self.level
        if self.info.get(level) is None:
            self.var_status.set(f"No '{level}' column on this dataset.")
            return
        codes, uniques = self.category_codes(level)
        ids, bad = self.typed_groups(level, uniques)
        uid = self.group_keys(level, uniques)
        colors = {t: (GREY if t in self.off_groups else np.array([255, 0, 0], np.uint8) if len(ids) == 1 else rank_color(i)) for i, t in enumerate(ids)}
        lut = np.tile(GREY, (len(uniques), 1))
        for t, col in colors.items():
            lut[uid == t] = col
        rgb = np.where((codes >= 0)[:, None], lut[np.maximum(codes, 0)], LIGHT_GREY).astype(np.uint8)         # (code -1: a cell with no annotation at this level: gray)
        hit = np.zeros(len(uniques), np.int32)
        for i, t in enumerate(ids):
            hit[uid == t] = i + 1
        self.rank_run = np.where(codes >= 0, hit[np.maximum(codes, 0)], 0).astype(np.int32)
        order = np.argsort(self.rank_run, kind="stable")
        bgc = self.bg_codes(level) if self.bg else None

        def bg_fn(label):
            if not self.bg:
                return None
            w = self.bg["where"].get(label)
            if w is None:
                return None
            c = bgc[w] if bgc is not None else np.full(len(w), -1)
            h = np.where(c >= 0, hit[np.maximum(c, 0)], 0)
            col = np.where((h > 0)[:, None], lut[np.maximum(c, 0)], LIGHT_GREY).astype(np.uint8)
            o = np.argsort(h, kind="stable")
            return self.bg["x"][w][o], self.bg["y"][w][o], col[o], w[o]

        self.apply(rgb, order, bg_fn)
        counts = {t: int((uid[np.maximum(codes, 0)] == t)[codes >= 0].sum()) for t in ids}
        full = {t: (str(uniques[np.flatnonzero(uid == t)[0]]) if (uid == t).any() else f"{self.level_names[level]} {t}") for t in ids}
        self.last_ids = {"level": level, "targets": ids, "colors": colors, "counts": counts}
        self.view_kind, self.view_token = "single", f"{level}_{'-'.join(str(t) for t in ids)}"
        self.view_title = f"{self.level_names[level]} {', '.join(str(t) for t in ids)}"
        order_of = {t: i for i, t in enumerate(ids)}
        self.export_cat = {"rank": np.where((codes >= 0) & (hit[np.maximum(codes, 0)] > 0), hit[np.maximum(codes, 0)] - 1, -1), "title": self.level_names[level],
                           "entries": [(hexcolor(colors[t]), full[t], counts[t]) for t in ids]}
        self.set_cell_ranks(self.export_cat["rank"], np.where((bgc >= 0) & (hit[np.maximum(bgc, 0)] > 0), hit[np.maximum(bgc, 0)] - 1, -1)
                            if bgc is not None else None)
        self.shown_values = None
        parts = []
        if ids:
            parts.append("Highlighting " + ", ".join(f"{t} ({counts[t]} cells)" for t in ids if counts[t]))
            if any(not counts[t] for t in ids):
                parts.append("not found: " + ", ".join(str(t) for t in ids if not counts[t]))
        else:
            parts.append(f"Enter one or more {self.level_names[level].lower()} {'numbers' if self.level_has_ids.get(level, True) else 'names'}")        # (the joined status ends with its own period)
        if bad:
            parts.append("ignored: " + ", ".join(bad))
        self.var_status.set("; ".join(parts) + ".")
        self.single_names = full
        self.draw_legend(ids, colors, counts, full)
        self.values = None

    def draw_legend(self, ids, colors, counts, names=None):
        if len(ids) < 2:
            return
        self.legend.pack(fill=tk.X, pady=(6, 0), after=self.box)
        self.sync_select_bar()
        self.legend.delete("all")
        self.legend.config(height=22 * len(ids) + 8)
        numbered = self.level_has_ids.get(self.level, True)
        rows = [("#%02x%02x%02x" % tuple(int(v) for v in colors[t]), str(t) if numbered else "", f"({counts[t]})",
                 re.sub(r"^\s*\d+\s*", "", (names or {}).get(t, "")) if numbered else str(t), t in self.off_groups,
                 int(self.shape_codes(np.array([i]))[0]))
                for i, t in enumerate(ids)]
        self.draw_legend_rows(self.legend, rows, 4)

    # ---- genes

    def symbols(self):
        if self.own is not None:                                       # (its genes, and its numeric columns: either can color the cells)
            return list(self.own.gene_names) + list(self.own.numbers)
        return core.gene_symbol_list(self.adata, self.imputed_state, self.mode == "imputed")

    def show_gene_warning(self, text):
        """The explanation under the gene box of what was wrong with the genes typed (nothing: the line is hidden)."""
        self.gene_warning.config(text=text)
        if text:
            self.gene_warning.pack(fill=tk.X, pady=(4, 0), after=self.gene_entry)
        else:
            self.gene_warning.pack_forget()

    def explain_missing_genes(self, missing, imputed):
        """Why each of the genes `missing` from the shown dataset could not be shown: it may be in the other dataset (switch to it), or in neither."""
        here, there = ("imputed", "standard") if imputed else ("standard", "imputed")
        there_choice = "Gene" if imputed else "Imputed Gene"          # (the choice above that shows the other dataset)
        other_adata = self.adata if imputed else self.imputed_state["adata"]
        other = None
        if other_adata is not None:
            other = {str(g).lower() for g in other_adata.var["gene_symbol"]}
        elsewhere = [g for g in missing if other is not None and g.lower() in other]
        nowhere = [g for g in missing if g not in elsewhere]
        lines = []
        if elsewhere:
            lines.append(f"{', '.join(elsewhere)}: not in the {here} dataset, but it is in the {there} dataset: choose \u201c{there_choice}\u201d above.")
        if nowhere:
            hint = " Check the spelling." if other is not None or self.own is not None else \
                " Check the spelling (the imputed genes can be checked once that dataset is loaded: choose \u201cImputed Gene\u201d)."
            lines.append(f"{', '.join(nowhere)}: not found in the {here} dataset.{hint}")
        return "\n".join(lines)

    def on_gene_key(self, e):
        text = self.var_gene.get()
        if text == getattr(self, "rejected_text", None):          # (the key that asked for the genes, or one that moves the cursor: the red and the warning stay)
            return
        self.rejected_text = None
        self.gene_entry.clear_bad()
        self.show_gene_warning("")
        if "\n" in text:                                       # (a pasted line break: names are separated by commas)
            self.var_gene.set(text.replace("\n", " "))
        if e.keysym in ("Return", "Up", "Down", "Escape", "Tab"):
            return
        token = self.var_gene.get().rsplit(",", 1)[-1].strip()
        matches = core.rank_gene_suggestions(self.symbols(), token, 6) if token else []
        self.sugg.delete(0, tk.END)
        for m in matches:
            self.sugg.insert(tk.END, m)
        if matches:
            self.sugg.config(height=len(matches))
            self.sugg.pack(fill=tk.X, after=self.gene_entry)
        else:
            self.sugg.pack_forget()

    def focus_suggestions(self):
        if self.sugg.winfo_ismapped() and self.sugg.size():
            self.sugg.focus_set()
            self.sugg.selection_clear(0, tk.END)
            self.sugg.selection_set(0)
        return "break"

    def pick_suggestion(self):
        sel = self.sugg.curselection()
        if not sel:
            return
        text = self.var_gene.get()
        prefix = text.rsplit(",", 1)[0].strip() + ", " if "," in text else ""
        self.var_gene.set(prefix + self.sugg.get(sel[0]))
        self.sugg.pack_forget()
        self.gene_entry.focus_set()
        self.recolor()

    def color_gene(self):
        self.sugg.pack_forget()
        text = self.var_gene.get().strip()
        names = core.parse_gene_names(text)
        if not names:
            self.cbar.delete("all")
            self.show_levels_slider(False)
            self.var_status.set("Enter a gene name.")
            self.color_all_quiet()
            return
        if len(names) > MAX_GENES:
            self.var_status.set(f"Enter at most {MAX_GENES} gene names, separated by commas.")
            return
        imputed = self.mode == "imputed"
        if not imputed and not self.has_counts:
            self.var_status.set("No raw counts on this dataset (layers['counts'] missing).")
            return
        source = self.imputed_state["adata"] if imputed else self.adata
        _cols, resolved, missing = self.own.resolve(names) if self.own is not None else core.resolve_gene_names(source, names)
        if missing:
            self.gene_entry.mark_bad(missing)
            self.rejected_text = self.var_gene.get()
            self.show_gene_warning(self.explain_missing_genes(missing, imputed))
            self.var_status.set(f"Not found in the {'imputed' if imputed else 'standard'} dataset: {', '.join(missing)}.")
            return
        self.rejected_text = None
        self.gene_entry.clear_bad()
        self.show_gene_warning("")
        self.var_gene.set(", ".join(resolved))
        scale = self.var_scale.get()
        key = (tuple(resolved), imputed, scale)
        if self.gene_cache is not None and self.gene_cache[0] == key:                   # (a legend click, or the same genes again)
            self.show_genes(resolved, imputed, scale, self.gene_cache[1])
            return
        backed = self.adata_backed if self.bg else None
        bg_ids = self.bg["ids"] if self.bg else None

        def read_values(gene):
            return self.read_run_gene(gene, imputed, scale)

        def read_background(gene):
            """The gene's value in every cell of the sections (the whole unfiltered dataset), or None if it can't be had."""
            if backed is None:
                return None
            try:
                if imputed:
                    full = core._get_full_gene_column(source, core.find_gene_index(source, gene))
                    pos = source.obs_names.get_indexer(bg_ids)
                    out = np.full(len(bg_ids), np.nan)
                    out[pos >= 0] = full[pos[pos >= 0]]
                    return np.exp2(out) - 1.0 if scale == "linear" else out
                bcol = core.find_gene_index(backed, gene)
                if bcol is None:
                    return None
                out = core._get_full_gene_column(backed, bcol)
                return core.log2p1(out) if scale == "log" else out
            except Exception as e:
                print(f"Could not read {gene} for the cells outside the selection: {e}", flush=True)
                return None

        def work():
            return [(read_values(g), read_background(g)) for g in resolved]

        def done(res):
            if isinstance(res, BaseException):
                self.var_status.set(f"Could not read {', '.join(resolved)}: {res}")
                return
            self.keep_gene_settings(self.gene_cache[0] if self.gene_cache is not None else None, key, res)
            self.gene_cache = (key, res)
            self.show_genes(resolved, imputed, scale, res)

        self.job(f"Reading {', '.join(resolved)}...", work, done)

    def zero_color(self):
        """The colour of a cell with no expression (the lowest value) in the gene view on the screen, or None in the other views. Those cells are
        most of a run and all look alike, so a saved SVG holds them as a picture, not as a dot apiece."""
        if self.mode not in ("gene", "imputed"):
            return None
        if self.view_kind == "multi":
            return (230, 230, 230) if self.var_theme.get() == "light" else (0, 0, 0)
        if self.view_kind == "gene":
            return tuple(int(v) for v in self.gene_rgb(np.array([0.0]))[0])
        return None

    def units_text(self, dataset, scale):
        """What the color bar says the values are: for a numeric column of a dataset of your own, just 'value'."""
        shown = self.gene_cache[0][0] if self.gene_cache is not None else ()
        if self.own is not None and len(shown) == 1 and self.own.is_number(shown[0]):
            return "value"
        return self.scale_units()[(dataset, scale)]

    def scale_units(self):
        """What the color bar says the values are."""
        if self.own is not None:
            linear = "expr (normalized)" if self.own.normalization == "cp10k" or self.own.kind != "counts" else "expr (counts)"
            return {("standard", "log"): "log2(1 + expr)", ("standard", "linear"): linear}
        return SCALE_UNITS

    def default_levels_for(self, vals, scale):
        """The black and white levels a gene starts at: the whole range in the log scale; in the linear scale, where the few cells that express a gene most
        strongly make the rest all but black, the white point at the 98th percentile of the cells that express it (LINEAR_WHITE_PERCENTILE)."""
        if scale != "linear":
            return (0.0, 1.0)
        fin = np.isfinite(vals)
        if not fin.any():
            return (0.0, 1.0)
        lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
        v = value_at_percentile(self.histogram_of(vals), hi, LINEAR_WHITE_PERCENTILE)
        return (0.0, min(max((v - lo) / (hi - lo), 0.0), 1.0)) if v is not None and hi > lo else (0.0, 1.0)

    def keep_gene_settings(self, old_key, new_key, per_gene=None):
        """A different set of genes is going to be shown: a gene that was shown before keeps what was set for it (the black and white levels of its
        slider, and whether it was hidden); a new gene starts at the default levels (default_levels_for), shown. Nothing is kept when the dataset or the
        scale changed (the levels are fractions of a range that is not the same any more). `per_gene`: the genes' values, [(run, background)]."""
        names = list(new_key[0])
        self.default_levels = {i: (self.default_levels_for(per_gene[i][0], new_key[2]) if per_gene is not None else (0.0, 1.0)) for i in range(len(names))}
        levels, hidden = {}, set()
        if old_key is not None and old_key[1:] == new_key[1:]:
            old_names = list(old_key[0])
            old_levels = {old_names[0]: self.gene_levels} if len(old_names) == 1 else {g: self.multi_levels.get(i, (0.0, 1.0)) for i, g in enumerate(old_names)}
            hidden_names = {old_names[i] for i in self.hidden if i < len(old_names)} if len(old_names) > 1 else set()
            levels = {i: old_levels[g] for i, g in enumerate(names) if g in old_levels}
            hidden = {i for i, g in enumerate(names) if g in hidden_names}
        self.gene_levels = levels.get(0, self.default_levels[0]) if len(names) == 1 else (0.0, 1.0)
        self.multi_levels = {i: levels.get(i, self.default_levels[i]) for i in range(len(names))} if len(names) > 1 else {}
        self.hidden = hidden if len(names) > 1 else set()

    def read_run_gene(self, gene, imputed, scale):
        """A gene's value in each cell of the run (imputed or standard, on the linear or the log scale). Blocking: on a worker thread."""
        if self.own is not None:
            return self.own.column(gene, scale)
        adata = self.adata
        source = self.imputed_state["adata"] if imputed else adata
        col = core.find_gene_index(source, gene)
        if imputed:
            full = core._get_full_gene_column(source, col)
            pos = source.obs_names.get_indexer(adata.obs_names)
            vals = np.full(adata.n_obs, np.nan)
            found = pos >= 0
            vals[found] = full[pos[found]]
            if scale == "linear":
                vals = np.exp2(vals) - 1.0
        elif self.atlas_counts is not None:                      # (no UMAP: the gene's column of the atlas's own file, for the project's cells)
            vals = self.atlas_counts.column(col).astype(float)
            if scale == "log":
                vals = core.log2p1(vals)
        else:
            expr = adata.layers["counts"][:, col]
            expr = expr.toarray() if hasattr(expr, "toarray") else expr
            vals = np.asarray(expr).ravel().astype(float)
            if scale == "log":
                vals = core.log2p1(vals)
        return vals

    def counts_of_rows(self, rows):
        """The raw counts of the run's cells `rows` (indices into the run) as a cells x genes matrix (scipy sparse or a numpy array): from the run's own
        counts layer, or, for a project without a UMAP, read from the atlas."""
        if self.atlas_counts is not None:
            return self.atlas_counts.rows(rows)
        return self.adata.layers["counts"][rows]

    def counts_adata(self):
        """A copy of the run's AnnData with the raw counts as its X, for the exports (a project without a UMAP: the counts of all its cells read from the
        atlas, so a big project needs the memory for them)."""
        if self.atlas_counts is None:
            m = self.adata.copy()
            m.X = m.layers["counts"].copy()
            return m
        import anndata
        m = anndata.AnnData(X=self.atlas_counts.rows(np.arange(self.adata.n_obs)), obs=self.adata.obs.copy(), var=self.adata.var.copy())
        return m

    def open_analyze(self):
        """The Analyze button: a window to plot genes against each other, per cell or per group."""
        win = getattr(self, "analyze_window", None)
        try:
            if win is not None and win.top.winfo_exists():
                win.top.lift()
                return
        except tk.TclError:
            pass
        from abc_tk_analyze import AnalyzeWindow
        self.analyze_window = AnalyzeWindow(self)

    def color_all_quiet(self):
        mode = self.mode
        self.mode = "all"
        self.color_all()
        self.mode = mode

    def show_genes(self, genes, imputed, scale, per_gene):
        dataset = "imputed" if imputed else "standard"
        if self.mode not in ("gene", "imputed") or (self.mode == "imputed") != imputed:
            return
        if self.theme_choice is None:                         # (not chosen yet: several genes are mixed on a dark background, one gene is on a light one)
            self.var_theme.set("dark" if len(genes) > 1 else "light")
        self.values = [(g, pv[0]) for g, pv in zip(genes, per_gene)]
        self.view_token = "_".join(genes) + ("_imputed" if imputed else "")
        self.export_cat = None
        if len(genes) == 1:
            self.view_kind = "gene"
            self.view_title = f"{genes[0]}{' [imputed]' if imputed else ''}"
            self.shown_values = (genes[0], dataset, scale, per_gene[0][0])
            self.show_single(genes[0], dataset, scale, per_gene[0][0], per_gene[0][1])
        else:
            self.view_kind, self.multi_genes = "multi", list(genes)
            self.view_title = " / ".join(f"{g} ({self.channel_names()[i]})" for i, g in enumerate(genes)) + (" [imputed]" if imputed else "")
            self.show_multi(genes, dataset, scale, per_gene)

    def show_single(self, gene, dataset, scale, vals, bg_vals):
        self.cell_marks = None                          # (the dots of the cell under the pointer are for the genes shown: they come back with the next move)
        fin = np.isfinite(vals)
        vmin, vmax = (float(np.nanmin(vals)), float(np.nanmax(vals))) if fin.any() else (0.0, 1.0)
        # The black and white levels can be moved (the slider under the color bar), in either scale: they cut the range the colors cover.
        self.clear_multi_sliders()
        low, high = self.gene_levels
        whole = (vmin, vmax)                            # (the histogram runs over all of it, whatever the levels)
        vmin, vmax = vmin + low * (vmax - vmin), vmin + high * (vmax - vmin)
        self.shown_range = (vmin, vmax)                 # (the colors' range as shown, with the black and white levels: a saved picture's color bar has it too)
        self.show_levels_slider(True)
        self.cbar_slider.set_even(scale == "log")                                       # (a log scale is already spread out: an even track)
        self.cbar_slider.default = self.default_levels.get(0, (0.0, 1.0))
        self.cbar_slider.set(low, high)
        self.hist_args = (self.histogram_of(vals), whole[0], whole[1], scale == "log", low, high, self.gene_rgb, self.units_text(dataset, scale))
        self.redraw_hist()
        t = np.clip((vals - vmin) / (vmax - vmin), 0, 1) if vmax > vmin else np.zeros_like(vals)
        rgb = self.gene_rgb(np.where(fin, t, 0))
        rgb[~fin] = GREY
        order = np.argsort(np.where(fin, vals, -np.inf), kind="stable")
        self.rank_run = np.argsort(order, kind="stable").astype(np.int64)       # the rank of each cell: highest expression on top

        def bg_fn(label):
            if not self.bg:
                return None
            w = self.bg["where"].get(label)
            if w is None:
                return None
            if bg_vals is None:
                return self.bg["x"][w], self.bg["y"][w], np.tile(LIGHT_GREY, (len(w), 1)), w
            v = bg_vals[w]
            ok = np.isfinite(v)
            tt = np.clip((v - vmin) / (vmax - vmin), 0, 1) if vmax > vmin else np.zeros_like(v)
            col = self.gene_rgb(np.where(ok, tt, 0))
            col[~ok] = LIGHT_GREY
            o = np.argsort(np.where(ok, v, -np.inf), kind="stable")           # (highest on top)
            return self.bg["x"][w][o], self.bg["y"][w][o], col[o], w[o]

        self.apply(rgb, order, bg_fn, dark=self.var_theme.get() == "dark")             # (the dark theme: a black background, as in the multi-gene view)
        info = self.histogram_of(vals)
        readout = (expressing_percentile(vals, vmin, info, whole[1]), expressing_percentile(vals, vmax, info, whole[1]))
        self.legend_redraw = lambda: self.draw_colorbar(vmin, vmax, scale, dataset, readout=readout)
        self.legend_redraw()
        self.var_status.set(f"Showing {gene}{' [imputed]' if dataset == 'imputed' else ''} (every cell of the sections, colored on the selection's scale)."
                            if bg_vals is not None else f"Showing {gene}." if self.own is not None else "(cells outside this run's selection are gray).")

    # ---- several genes: red / green / blue circles for the first three (cyan / magenta / yellow in the light theme), the same three colors as small dots on top for the next three

    def channel_names(self):
        return CHANNEL_NAMES_LIGHT if self.var_theme.get() == "light" else CHANNEL_NAMES

    def channel_hex(self):
        return CHANNEL_HEX_LIGHT if self.var_theme.get() == "light" else CHANNEL_HEX

    def show_multi(self, genes, dataset, scale, per_gene):
        self.cell_marks = None
        self.show_levels_slider(False)                  # (the single gene's slider: each gene has its own in the legend)
        light = self.var_theme.get() == "light"
        mix = (lambda chans: ink_rgb(chans)) if light else (lambda chans: (core.multi_gene_rgb(chans, (0.0, 0.0, 0.0)) * 255).astype(np.uint8))
        n = len(genes)
        norm, wholes = [], []
        for i, (vals, _bg) in enumerate(per_gene):
            fin = np.isfinite(vals)
            lo, hi = (float(np.nanmin(vals)), float(np.nanmax(vals))) if fin.any() else (0.0, 1.0)
            wholes.append((lo, hi))                                   # (the histogram runs over all of it, whatever the levels)
            low, high = self.multi_levels.get(i, (0.0, 1.0))          # (the black and white levels set with the legend's sliders cut the range)
            lo, hi = lo + low * (hi - lo), lo + high * (hi - lo)
            norm.append((np.clip((vals - lo) / (hi - lo), 0, 1) if hi > lo else np.zeros_like(vals), lo, hi))
        self.multi_ranges = [(v[1], v[2]) for v in norm]           # (each gene's range as shown: for the legend of a saved picture)
        self.multi_scale = scale
        hidden = self.hidden
        zero = np.zeros(self.n)
        runn = [zero if i in hidden else norm[i][0] for i in range(n)]
        intensity = lambda arrs: np.maximum.reduce([np.nan_to_num(a, nan=0.0) for a in arrs])
        circ = runn[:3]
        rgb = mix(circ)
        order = np.argsort(intensity(circ), kind="stable")                      # (the most expressing on top)
        self.rank_run = np.argsort(order, kind="stable").astype(np.int64)
        extra = []
        if n > 3:
            dots = runn[3:6]
            inten = intensity(dots)
            vis = np.flatnonzero(inten > 0)
            vis = vis[np.argsort(inten[vis], kind="stable")]
            drgb = mix([d[vis] for d in dots])
            half = self.viewer_settings.get("gene456_shape") == "half circle"       # (the Settings choice: a half circle as big as the base one, or a small circle)
            extra = [{"x": self.coords[vis, 0], "y": self.coords[vis, 1], "rgb": drgb, "size": (4 if half else 2) * self.dot_scale(),
                      "shape": "half" if half else "circle"}]

        def bg_fn(label):
            if not self.bg:
                return None
            w = self.bg["where"].get(label)
            if w is None:
                return None
            chans = []
            for i in range(min(n, 3)):
                bgv = per_gene[i][1]
                if i in hidden or bgv is None:
                    chans.append(np.zeros(len(w)))
                    continue
                lo, hi = norm[i][1], norm[i][2]
                v = bgv[w]
                chans.append(np.clip((np.nan_to_num(v, nan=lo) - lo) / (hi - lo), 0, 1) if hi > lo else np.zeros(len(w)))
            col = mix(chans)
            o = np.argsort(intensity(chans), kind="stable")
            return self.bg["x"][w][o], self.bg["y"][w][o], col[o], w[o]

        self.apply(rgb, order, bg_fn, extra_umap_layers=extra, dark=not light)
        readouts = [(expressing_percentile(per_gene[i][0], norm[i][1], self.histogram_of(per_gene[i][0]), wholes[i][1]),
                     expressing_percentile(per_gene[i][0], norm[i][2], self.histogram_of(per_gene[i][0]), wholes[i][1])) for i in range(n)]
        channel = lambda i: (lambda t, c=i % 3: mix([t if k == c else np.zeros_like(t) for k in range(3)]))          # (a gene's color at a place in its range: its own ink alone)
        legend_ranges = [(v[1], v[2]) for v in norm]
        legend_hists = [(self.histogram_of(per_gene[i][0]),) + wholes[i] + (channel(i), self.units_text(dataset, scale)) for i in range(n)]
        self.legend_redraw = lambda: self.draw_multi_legend(genes, legend_ranges, scale, readouts, legend_hists)
        self.legend_redraw()
        parts = " and ".join(f"'{g}' ({self.channel_names()[i]})" for i, g in enumerate(genes))
        self.var_status.set(f"Showing {parts}{' [imputed]' if dataset == 'imputed' else ''}."
                            + (f" Hidden: {', '.join(genes[i] for i in sorted(hidden))}." if hidden else ""))

    def draw_multi_legend(self, genes, ranges, scale, readouts=None, hists=None):
        c = self.cbar
        self.cbar_photo = None
        self.cbar_width = c.winfo_width()
        c.hist_regions, c.hist_photos = [], []                      # (the histograms are drawn again below)
        h8, h10 = self.line_height(8), self.line_height(10, bold=True)
        name_y = h10 // 2 + 3                                       # (the name and its circle; then the range; where the two levels stand among the cells that express the gene, and their slider)
        range_y = h10 + 4
        cells_y = range_y + h8 + 2
        slider_y = cells_y + h8 + 4
        hist_y = slider_y + 26 + 2                                  # (the histogram under the slider)
        row = hist_y + HIST_HEIGHT + HIST_PAD + 10
        if len(self.multi_sliders) != len(genes):
            c.delete("all")
            self.clear_multi_sliders()                              # (the sliders are kept while they are dragged: only made again for other genes)
            self.multi_sliders = [RangeSlider(c, lambda lo, hi, i=i: self.on_multi_levels(i, lo, hi)) for i in range(len(genes))]
            self.multi_windows = [c.create_window(0, 0, anchor="nw", window=sl, height=26) for sl in self.multi_sliders]
        else:
            c.delete("legend")                                      # (only the texts and circles: the sliders stay where they are, so that a drag in progress goes on)
            c.delete("hist")
        c.config(height=8 + row * len(genes) + 2 * h8 + 8)
        for i, (g, (lo, hi)) in enumerate(zip(genes, ranges)):
            y = 8 + row * i
            off = i in self.hidden
            color = "#999999" if off else self.channel_hex()[i]
            half = i >= 3 and self.viewer_settings.get("gene456_shape") == "half circle"      # (genes 4-6: the Settings choice, a half circle or a small circle)
            r = 6 if i < 3 or half else 3                          # (genes 4-6 are the small dots)
            if half:                                               # (the right half of a disc, as the cells of that gene are drawn)
                c.create_arc(18 - r, y + name_y - r, 18 + r, y + name_y + r, start=-90, extent=180, style=tk.CHORD, tags="legend",
                             **({"outline": "#999999"} if off else {"fill": color, "outline": ""}))
            elif off:
                c.create_oval(18 - r, y + name_y - r, 18 + r, y + name_y + r, outline="#999999", tags="legend")
            else:
                c.create_oval(18 - r, y + name_y - r, 18 + r, y + name_y + r, fill=color, outline="", tags="legend")
            c.create_text(34, y + name_y, anchor="w", text=g, fill=color, font=("Segoe UI", 10, "bold"), tags=("legend", f"row{i}"))      # (centered on the circle)
            c.create_text(34, y + range_y, anchor="nw", text="hidden" if off else f"{'' if scale == 'linear' else 'log2: '}{lo:.1f}\u2013{hi:.1f}",
                          fill="#666666", font=("Segoe UI", 8), tags="legend")
            if readouts is not None and not off:
                for e, x, anchor in ((readouts[i][0], 10, "nw"), (readouts[i][1], max(c.winfo_width(), 200) - 10, "ne")):
                    c.create_text(x, y + cells_y, anchor=anchor, text=f"{e:.1f} percentile", font=("Segoe UI", 8), fill="gray30", tags="legend")
            self.multi_sliders[i].set_even(scale == "log")
            self.multi_sliders[i].default = self.default_levels.get(i, (0.0, 1.0))
            self.multi_sliders[i].set(*self.multi_levels.get(i, (0.0, 1.0)))
            c.coords(self.multi_windows[i], 0, y + slider_y)
            c.itemconfig(self.multi_windows[i], width=max(c.winfo_width(), 200))
            if hists is not None and not off:
                info, whole_lo, whole_hi, color_fn, units = hists[i]
                low, high = self.multi_levels.get(i, (0.0, 1.0))
                self.draw_histogram(c, y + hist_y, max(c.winfo_width(), 200), info, whole_lo, whole_hi, scale == "log", low, high, color_fn, units, key=i)
        width = max(c.winfo_width(), 200)
        hint = c.create_text(width / 2, 8 + row * len(genes) + 2, anchor="n", width=width - 20, justify=tk.CENTER,          # (wrapped: it is wider than a narrow panel)
                             text="click a gene to hide / show it; drag the handles for its black and white levels",
                             fill="#888888", font=("Segoe UI", 8, "italic"), tags="legend")
        c.config(height=c.bbox(hint)[3] + 6)                        # (as tall as the hint, however many lines it took)
        self.multi_rows = (8, row, len(genes))
        self.apply_cell_marks()

    def clear_multi_sliders(self):
        for slider in self.multi_sliders:
            slider.destroy()
        self.multi_sliders = []
        self.multi_windows = []

    def on_multi_levels(self, gene_index, low, high):
        """A legend slider moved: the colors are redrawn (as for the single gene's slider: see on_gene_levels)."""
        self.multi_levels[gene_index] = (low, high)
        if self.levels_job is None:
            self.levels_job = self.root.after(60, self.apply_gene_levels)

    def toggle_group(self, key, shift=False, order=()):
        """Switch a group between its color and gray (a click on its legend swatch); the cells, the swatch and the labels follow. With the shift
        key, every group in the legend from the one last clicked to this one is put in the state the last-clicked one was left in."""
        if self.busy:
            return
        order = list(order)
        last = self.last_group_click
        if shift and last is not None and last[0] in order and key in order:
            a, b = sorted((order.index(last[0]), order.index(key)))
            if last[1]:
                self.off_groups.update(order[a:b + 1])
            else:
                self.off_groups.difference_update(order[a:b + 1])
            self.last_group_click = (key, last[1])
        else:
            self.off_groups.symmetric_difference_update({key})
            self.last_group_click = (key, key in self.off_groups)
        self.recolor(keep_off=True)

    def sync_select_bar(self):
        """All / None, right under the box, so above the group list or the legend. Only where a click on the entries switches them (All groups and
        Single group(s); with a gene the list is only what the cells were colored by before)."""
        self.select_bar.pack_forget()
        if self.mode in ("all", "single") and (self.group_list.winfo_manager() or self.legend.winfo_manager()):
            self.select_bar.pack(fill=tk.X, pady=(6, 0), after=self.box)

    def legend_keys(self):
        """The entries of the legend on the screen, as toggle_group knows them: the group list's (All), or the numbers of Single group(s)."""
        if self.mode == "all":
            return list(getattr(self, "group_rows", []))
        if self.mode == "single":
            return list((self.last_ids or {}).get("targets", []))
        return []

    def select_all_groups(self):
        """The All button: every entry of the legend is shown in its color."""
        if self.busy or not self.legend_keys():
            return
        self.last_group_click = None
        self.off_groups = set()
        self.recolor(keep_off=True)

    def select_no_groups(self):
        """The None button: every entry of the legend is switched off (gray)."""
        keys = self.legend_keys()
        if self.busy or not keys:
            return
        self.last_group_click = None
        self.off_groups = set(keys)
        self.recolor(keep_off=True)

    def on_group_list_click(self, e):
        rows = getattr(self, "group_rows", [])
        i = int((self.group_canvas.canvasy(e.y) - 4) // 22)
        if self.mode == "all" and 0 <= i < len(rows):
            self.toggle_group(rows[i], bool(e.state & 0x1), rows)

    def on_single_legend_click(self, e):
        ids = (self.last_ids or {}).get("targets", [])
        i = int((e.y - 4) // 22)
        if self.mode == "single" and len(ids) > 1 and 0 <= i < len(ids):
            self.toggle_group(ids[i], bool(e.state & 0x1), ids)

    # ---- the legend rows under the pointer

    def group_info_rows(self, level, label):
        """A group as table rows like a cell's: its place in the taxonomy (the levels above it, then its own), how many cells it has in the run and
        in how many sections. None if the run has no such group."""
        values = self.info.get(level)
        if values is None:
            return None
        mask = values == label
        if not mask.any():
            return None
        i0 = int(np.flatnonzero(mask)[0])
        levels = list(self.level_names)
        chain = [level] if level in NOT_IN_TAXONOMY or self.own is not None else [k for k in levels[:levels.index(level) + 1] if k not in NOT_IN_TAXONOMY]
        rows = []
        for k in chain:
            v = self.info[k][i0] if self.info.get(k) is not None else None
            if v is None or str(v) == "nan":
                continue
            m = re.match(r"^\s*(\d+)\s*(.*)$", str(v))
            rows.append((self.level_names[k], m.group(1) if m else "", m.group(2) if m else str(v)))
        rows.append(("Cells", f"{int(mask.sum())}", "in this run"))
        if self.section_of is not None:
            rows.append(("Sections", f"{len(np.unique(self.section_of[mask]))}", ""))
        return rows

    def show_legend_group(self, key):
        """Describe the legend row's group in a tip (key None: the pointer left the rows)."""
        if key == self.hover_group:
            return
        self.hover_group = key
        self.cancel_hover_jobs()
        self.hover_key = None
        self.clear_highlights()
        self.show_info(None)                         # (also puts the bold row back to normal)
        if key is None:
            self.group_tip.pointer_off()
            return
        label = key if self.mode == "all" else self.single_names.get(key)
        rows = self.group_info_rows(self.level, label) if label is not None else None
        if rows is not None:
            if self.mode == "all":
                canvas, names = self.group_canvas, getattr(self, "group_rows", [])
            else:
                canvas, names = self.legend, (self.last_ids or {}).get("targets", [])
            self.bold_legend_row(canvas, names.index(key) if key in names else None)
            self.group_tip.pointer_on(("group", self.mode, self.level, key), lambda: rows)
        else:
            self.group_tip.pointer_off()
        if rows is not None:
            self.set_dim(True)                       # (as when the pointer rests on one of the group's cells)
            self.mark_group(label)

    def on_group_list_motion(self, e):
        rows = getattr(self, "group_rows", [])
        i = int((self.group_canvas.canvasy(e.y) - 4) // 22)
        self.show_legend_group(rows[i] if self.mode == "all" and 0 <= i < len(rows) else None)

    def on_single_legend_motion(self, e):
        ids = (self.last_ids or {}).get("targets", [])
        i = int((e.y - 4) // 22)
        self.show_legend_group(ids[i] if self.mode == "single" and len(ids) > 1 and 0 <= i < len(ids) else None)

    def on_legend_click(self, e):
        if self.mode not in ("gene", "imputed") or self.multi_rows is None or self.gene_cache is None:
            return
        top, row, n = self.multi_rows
        i = int((e.y - top) // row)
        if not 0 <= i < n or len(self.gene_cache[0][0]) < 2 or self.busy:
            return
        self.hidden.symmetric_difference_update({i})
        key = self.gene_cache[0]
        self.show_genes(list(key[0]), key[1], key[2], self.gene_cache[1])

    def on_gene_levels(self, low, high):
        """The levels slider moved: the colors are redrawn (at most about 15 times a second while it is dragged)."""
        self.gene_levels = (low, high)
        if self.levels_job is None:
            self.levels_job = self.root.after(60, self.apply_gene_levels)

    def apply_gene_levels(self):
        self.levels_job = None
        if self.gene_cache is None or self.mode not in ("gene", "imputed"):
            return
        if self.busy:                                           # (something else is being read: try again shortly)
            self.levels_job = self.root.after(200, self.apply_gene_levels)
            return
        key = self.gene_cache[0]
        self.show_genes(list(key[0]), key[1], key[2], self.gene_cache[1])

    def show_levels_slider(self, on):
        if on:
            self.cbar_slider.pack(fill=tk.X, after=self.cbar)
            self.hist.pack(fill=tk.X, after=self.cbar_slider)
        else:
            self.cbar_slider.pack_forget()
            self.hist.pack_forget()

    def histogram_of(self, vals):
        """The data of a gene's histogram (expression_histogram), made once for the values on the screen: the sliders redraw as they are dragged."""
        hit = self.hist_cache.get(id(vals))
        if hit is None or hit[0] is not vals:
            if len(self.hist_cache) >= 8:
                self.hist_cache.pop(next(iter(self.hist_cache)))
            hit = self.hist_cache[id(vals)] = (vals, expression_histogram(vals))
        return hit[1]

    def draw_histogram(self, canvas, y_top, width, info, lo, hi, even, low, high, color_fn, label, key=0):
        """The histogram in `canvas` from y_top down HIST_HEIGHT pixels: along the track of a slider of this width (RangeSlider: the same margins and the same
        mapping), the number of cells in each pixel's stretch of it, as the logarithm of 1 + the number, over the whole range `lo`..`hi`. Each bar has the color
        the cells of its value have (`color_fn(t)`, t the place in the colors' range, between the slider's handles `low` and `high`, fractions of the whole
        range: below the black handle the cells are the first color, above the white one the last); the bars outside the handles are paler. The pointer over
        it says the value and the height (see histogram_readout)."""
        margin = 10
        track = max(int(width) - 2 * margin, 1)
        counts = histogram_columns(info, lo, hi, even, track)
        canvas.hist_regions.append({"x0": margin, "y0": y_top, "n": track, "counts": counts, "lo": lo, "hi": hi, "even": even, "label": label, "key": key,
                                    "color_fn": color_fn, "black": lo + low * (hi - lo), "white": lo + high * (hi - lo)})
        heights = np.log10(1.0 + counts)
        top = float(heights.max())
        if top <= 0:
            return
        rows = HIST_HEIGHT
        px = np.where(counts > 0, np.maximum(np.round(heights / top * (rows - 3)).astype(int), 1), 0)            # (a column with any cell is at least one pixel high)
        centre = (np.arange(track) + 0.5) / track
        value = lo + (centre if even else centre ** 2) * (hi - lo)
        black, white = lo + low * (hi - lo), lo + high * (hi - lo)
        t = np.clip((value - black) / (white - black), 0.0, 1.0) if white > black else np.zeros(track)
        colors = np.asarray(color_fn(t), dtype=np.uint8)
        bar = (rows - 1 - np.arange(rows))[:, None] < px[None, :]                                                # (rows from the bottom up to each column's height)
        picture = np.zeros((rows, track, 4), np.uint8)
        picture[..., :3] = colors[None, :, :]
        a = int(round((low if even else math.sqrt(low)) * track))
        b = int(round((high if even else math.sqrt(high)) * track))
        inside = (np.arange(track) >= a) & (np.arange(track) < max(b, a))
        picture[..., 3] = np.where(bar, np.where(inside, 255, 110)[None, :], 0)
        left, right, above = np.zeros_like(bar), np.zeros_like(bar), np.zeros_like(bar)       # (a gray edge round the bars: a bar of the palest color shows on the pane)
        left[:, 1:], right[:, :-1], above[1:, :] = bar[:, :-1], bar[:, 1:], bar[:-1, :]
        edge = bar & (~left | ~right | ~above)
        picture[edge] = (118, 118, 118, 255)
        picture[..., 3][edge & ~inside[None, :].repeat(rows, axis=0)] = 140
        photo = ImageTk.PhotoImage(Image.fromarray(picture, "RGBA"), master=canvas)
        canvas.hist_photos.append(photo)                                                                          # (kept: Tk shows the picture only while Python holds it)
        canvas.create_image(margin, y_top, anchor="nw", image=photo, tags="hist")
        canvas.create_line(margin, y_top + rows - 1, margin + track, y_top + rows - 1, fill="#a8a8a8", tags="hist")

    def histogram_readout(self, canvas, event):
        """What the pointer is over in a histogram of `canvas` (as draw_histogram left it): the value there and the height of the bar, the number of cells in
        that pixel's stretch of the track (None if it is not over one)."""
        for r in getattr(canvas, "hist_regions", []):
            x, y = canvas.canvasx(event.x) - r["x0"], canvas.canvasy(event.y) - r["y0"]
            if 0 <= x < r["n"] and 0 <= y < HIST_HEIGHT:
                p = (int(x) + 0.5) / r["n"]
                v = r["lo"] + (p if r["even"] else p * p) * (r["hi"] - r["lo"])
                n = float(r["counts"][int(x)])
                cells = "no cells" if n <= 0 else "fewer than 1 cell" if n < 0.5 else f"{n:,.0f} cell{'s' if round(n) != 1 else ''}"
                return f"{r['label']}: {v:.2f}\n{cells} (the bar's height)"
        return None

    def cbar_resized(self):
        """The side panel was made wider or narrower: the color bar (or the legend of several genes, with its histograms) is drawn again for the new width,
        once the panel has stopped changing. (Only for a change of width: drawing it changes its height, which comes back here.)"""
        if self.legend_redraw is None or self.cbar.winfo_width() == self.cbar_width:
            return
        if self.cbar_job is not None:
            self.root.after_cancel(self.cbar_job)
        self.cbar_job = self.root.after(60, self.run_legend_redraw)

    def run_legend_redraw(self):
        self.cbar_job = None
        if self.legend_redraw is not None and self.mode in ("gene", "imputed"):
            self.legend_redraw()

    def redraw_hist(self):
        """The histogram under the single gene's slider, again (a new gene or levels, or the panel changed width)."""
        self.hist.delete("all")
        self.hist.hist_regions, self.hist.hist_photos = [], []
        if self.hist_args is not None and self.hist.winfo_manager():
            info, lo, hi, even, low, high, color_fn, label = self.hist_args
            self.draw_histogram(self.hist, 0, max(self.hist.winfo_width(), 200), info, lo, hi, even, low, high, color_fn, label, key=0)
        self.apply_cell_marks()

    def apply_cell_marks(self):
        """The dots of the cell under the pointer on the histograms (in the single gene's, and in each gene's of the legend): at the place of its value on the
        histogram's axis, which is the slider's track, in the color the cell has on the maps (the histogram's there), with a black rim."""
        marks = self.cell_marks or {}
        for canvas in (self.hist, self.cbar):
            canvas.delete("cellmark")
            for r in getattr(canvas, "hist_regions", []):
                v = marks.get(r["key"])
                if v is None or not np.isfinite(v):
                    continue
                span = r["hi"] - r["lo"]
                f = min(max((v - r["lo"]) / span, 0.0), 1.0) if span > 0 else 0.0
                x = r["x0"] + (f if r["even"] else math.sqrt(f)) * r["n"]
                y = r["y0"] + HIST_HEIGHT - 1
                t = min(max((v - r["black"]) / (r["white"] - r["black"]), 0.0), 1.0) if r["white"] > r["black"] else 0.0
                color = np.asarray(r["color_fn"](np.array([t])), dtype=np.uint8)[0]
                canvas.create_oval(x - 6, y - 6, x + 6, y + 6, fill="white", outline="", tags="cellmark")                   # (a white halo: the dot shows on the dark bars of the dark theme)
                canvas.create_oval(x - 4, y - 4, x + 4, y + 4, fill="#%02x%02x%02x" % tuple(int(c) for c in color[:3]), outline="black", width=1, tags="cellmark")

    def mark_cell_on_histograms(self, info):
        """The cell `info` is under the pointer (None: it is not on one): while genes are shown, a dot on each gene's histogram where that cell's value
        is. A cell of the run has its value in the genes read for the run; any other cell of the sections in those read for the whole atlas."""
        marks = None
        if info is not None and self.mode in ("gene", "imputed") and self.gene_cache is not None and self.view_kind in ("gene", "multi"):
            marks = {}
            i, g = info["umap_idx"], info["bg_idx"]
            for k, (vals, bg) in enumerate(self.gene_cache[1]):
                v = vals[i] if i is not None else bg[g] if g is not None and bg is not None else None
                if v is not None:
                    marks[k] = float(v)
        if marks != self.cell_marks:
            self.cell_marks = marks
            self.apply_cell_marks()
    def line_height(self, size, bold=False):
        """The height in pixels of a line of Segoe UI text of this size on this screen (the scaling of Windows makes it more than the points say)."""
        fonts = self.__dict__.setdefault("line_fonts", {})
        key = (size, bold)
        if key not in fonts:
            fonts[key] = tkfont.Font(root=self.root, family="Segoe UI", size=size, weight="bold" if bold else "normal")
        return fonts[key].metrics("linespace")

    def draw_colorbar(self, vmin, vmax, scale, dataset, readout=None):
        """`readout`: where the black and the white point stand among the cells that express the gene (see expressing_percentile): (black, white)."""
        c = self.cbar
        self.cbar_width = c.winfo_width()
        c.delete("all")
        c.hist_regions, c.hist_photos = [], []                      # (no histograms in it: they are the legend's of several genes)
        h8, h9 = self.line_height(8), self.line_height(9)
        w = max(c.winfo_width(), 200)
        ramp = self.gene_rgb(np.linspace(0, 1, w - 20))[None, :, :]
        img = np.repeat(ramp, 16, axis=0)
        self.cbar_photo = ImageTk.PhotoImage(Image.fromarray(img), master=c)
        c.create_image(10, 8, anchor="nw", image=self.cbar_photo)
        c.create_text(10, 28, anchor="nw", text=f"{vmin:.2f}", font=("Segoe UI", 9))
        c.create_text(w - 10, 28, anchor="ne", text=f"{vmax:.2f}", font=("Segoe UI", 9))
        y = 28 + h9 + 2
        if readout is not None:
            for e, x, anchor in ((readout[0], 10, "nw"), (readout[1], w - 10, "ne")):
                c.create_text(x, y, anchor=anchor, text=f"{e:.1f} percentile", font=("Segoe UI", 8), fill="gray30")
            y += h8 + 4
        c.create_text(w / 2, y + h9 / 2, text=self.units_text(dataset, scale), font=("Segoe UI", 9), fill="gray30")
        c.config(height=int(y + h9 + 8))

    # ---- background work

    def job(self, text, work, done):
        if self.busy:
            return False
        self.busy = True
        self.var_status.set(text)
        self.root.config(cursor="watch")

        def target():
            try:
                res = work()
            except BaseException as e:
                res = e
            self.msgs.put((done, res))

        threading.Thread(target=target, daemon=True).start()
        return True

    def poll(self):
        try:
            while True:
                done, res = self.msgs.get_nowait()
                self.busy = False
                self.root.config(cursor="")
                done(res)
        except queue.Empty:
            pass
        self.ticks = getattr(self, "ticks", 0) + 1
        if self.ticks % 20 == 0 and not self.busy:          # (garbage is collected here, on the Tk thread: a worker thread freeing Tk images fails)
            gc.collect()
        self.poll_id = self.root.after(300, self.poll)

    # ---- hover and the link between the UMAP and the section panels

    def number_text(self, name, i):
        """A numeric column's value in run cell `i`, as text."""
        values = self.number_arrays.get(name)
        if values is None:
            values = self.number_arrays[name] = self.own.number_values(name)
        return f"{values[i]:.4g}" if np.isfinite(values[i]) else "n/a"

    def level_text(self, info):
        """The lines describing the cell under the pointer."""
        lines = []
        sec = info["sec"]
        if sec is not None:
            lines.append(f"Section {core.sanitize_section_token(sec)}")
        region = self.region_of(info)
        if region is not None:
            lines.append(f"Region: {region[0]} {region[1]}" if region[0] != "unassigned" else "Region: unassigned (outside the CCF atlas)")
        i = info["umap_idx"]
        for k in self.level_names:
            if i is not None:
                v = self.info[k][i] if self.info.get(k) is not None else None
            else:
                arr = self.bg["levels"].get(k) if self.bg else None
                v = arr[info["bg_idx"]] if arr is not None else None
            if v is not None and str(v) != "nan":
                lines.append(f"{self.level_names[k]}: {v}")
        if self.own is not None and i is not None:
            for name in self.own.numbers:
                lines.append(f"{name}: {self.number_text(name, i)}")
        if i is not None and self.values is not None:
            for gene, vals in self.values:
                if self.own is not None and self.own.is_number(gene):
                    continue
                lines.append(f"{gene} = {vals[i]:.3g}" if np.isfinite(vals[i]) else f"{gene} = n/a")
        elif i is None:
            lines.append("(not in the UMAP)")
        return "\n".join(lines)

    def region_index(self, info):
        """The number of the CCF region Allen gives the hovered cell, or None if it has none."""
        i = info["umap_idx"]
        if i is not None:
            value = self.region_of_run[i] if self.region_of_run is not None else np.nan
        else:
            arr = self.bg.get("region") if self.bg else None
            value = arr[info["bg_idx"]] if arr is not None else np.nan
        if value is None or not np.isfinite(value):
            return None
        return int(value)

    def region_of(self, info):
        """(acronym, name, structure, division) of the CCF region Allen gives the hovered cell (the registration's own answer, whatever the outlines
        drawn show), or None if there is no table of regions or the cell has none."""
        if self.region_names is None:
            return None
        value = self.region_index(info)
        return None if value is None else self.region_names.get(value)

    def cell_rows(self, info):
        """The hovered cell as table rows: (what, number, name). A group's label '268 CS-PRNr-DR ...' is split into its number and its name."""
        rows = []
        sec = info["sec"]
        if sec is not None:
            rows.append(("Section", core.sanitize_section_token(sec), "", ("section", sec)))          # (a token: the row lights the section's cells)
        region = self.region_of(info)
        if region is not None:
            acronym, name, structure, division = region
            number = self.region_index(info)               # (the rows are tokens: each lights the cells of its region, its structure or its division)
            if acronym == "unassigned":
                rows.append(("Region", "", "unassigned (outside the CCF atlas)", ("region", ("region", number))))
            else:
                rows.append(("Region", acronym, name, ("region", ("region", number))))
                if structure and structure not in (name, "unassigned"):
                    rows.append(("", "", f"in {structure}", ("region", ("structure", structure))))
                if division and division not in (structure, name, "unassigned"):
                    rows.append(("Division", "", division, ("region", ("division", division))))
        i = info["umap_idx"]
        for k in self.level_names:
            if i is not None:
                v = self.info[k][i] if self.info.get(k) is not None else None
            else:
                arr = self.bg["levels"].get(k) if self.bg else None
                v = arr[info["bg_idx"]] if arr is not None else None
            if v is None or str(v) == "nan":
                continue
            m = re.match(r"^\s*(\d+)\s*(.*)$", str(v))
            rows.append((self.level_names[k], m.group(1) if m else "", m.group(2) if m else str(v), (k, v) if self.info.get(k) is not None else None))
        if self.own is not None and i is not None:                      # (the numeric columns of a dataset of your own, under the groups)
            for name in self.own.numbers:
                rows.append((name, self.number_text(name, i), ""))
        rows += self.coord_rows(info)
        if i is not None and self.values is not None:
            for gene, vals in self.values:
                if self.own is not None and self.own.is_number(gene):          # (already listed with the numeric columns)
                    continue
                rows.append((gene, f"{vals[i]:.3g}" if np.isfinite(vals[i]) else "n/a", ""))
        elif i is None:
            rows.append(("", "", "(not in the UMAP)"))
        return rows

    def coord_rows(self, info):
        """The cell's position: its CCF AP, and approximate AP and DV from bregma and ML from the midline (core.paxinos_coords). From the cell's
        own CCF coordinates once they are read (see start_ccf_lookup); until then AP only, from its section's place in the CCF."""
        sec = info["sec"]
        if sec is None:
            return []
        core.start_section_coordinate_load(self.abc_cache)
        z, section_ap = core.section_coords(sec)
        rows = []
        if not self.ccf_started:
            self.start_ccf_lookup()                                  # (the first cell tip: not before, so that nothing reads it while the maps are being drawn and zoomed)
        ccf = None
        if info["umap_idx"] is not None and self.ccf_run is not None:
            ccf = self.ccf_run[info["umap_idx"]]
        elif info["bg_idx"] is not None and self.ccf_bg is not None:
            ccf = self.ccf_bg[info["bg_idx"]]
        if ccf is not None and np.isfinite(ccf).all():
            ap, ml, dv = core.paxinos_coords(float(ccf[0]), float(ccf[1]), float(ccf[2]))
            rows += [("CCF AP", f"{ccf[0]:.2f}", "mm from the front of the CCF"),
                     ("AP from bregma", f"{ap:+.2f}", "mm (CCF offset, approx.)"), ("ML from midline", f"{ml:+.2f}", "mm (CCF; + is right)"),
                     ("DV from bregma", f"{dv:+.2f}", "mm (CCF offset, approx.)")]
        elif section_ap is not None:
            ap, _ml, _dv = core.paxinos_coords(section_ap)
            rows += [("CCF AP", f"{section_ap:.2f}", "mm from the front of the CCF (section average)"),
                     ("AP from bregma", f"{ap:+.2f}", "mm (CCF offset, approx., section average)")]
            if self.ccf_state == "reading":
                rows.append(("ML, DV", "", "reading the cells' CCF coordinates..."))
        return rows

    def cell_swatch(self, info):
        """(rgb, marker code, background rgb): the cell as it is drawn where the pointer is, for its tip; None if that is not known. A cell of
        the run has its own color and its group's marker; a cell of a section that is not in the run, the color its section map gives it."""
        try:
            on_section = info["source"] == "section"
            back = self.umap.bgcolor
            if on_section:
                tile = next((t for t in self.tiles if t["label"] == info["sec"]), None)
                back = tile["pc"].bgcolor if tile is not None else self.view_background()
            i, g = info["umap_idx"], info["bg_idx"]
            if i is not None:
                rank = getattr(self, "run_cell_rank", None)
                return self.rgb_run[i], int(self.shape_codes([rank[i]])[0]) if rank is not None else 0, back
            if g is None or self.tile_fn is None:
                return None
            got = self.tile_fn(info["sec"])
            if got is None:
                return None
            _x, _y, rgb, w = got
            at = np.flatnonzero(w == g)
            if not len(at):
                return None
            rank = getattr(self, "bg_cell_rank", None)
            code = int(self.shape_codes([rank[g]])[0]) if rank is not None and not self.var_roi_only.get() else 0
            return rgb[at[0]], code, back
        except Exception:
            return None                                  # (no swatch, rather than no tip)

    def ccf_action(self):
        """The cell tip's extra button: fetch the CCF coordinates that its ML and DV need, while they are not on this computer."""
        if self.ccf_state == "missing" and "ensure_ccf_downloaded" in self.host:
            return "Get ML/DV...", self.download_ccf_coords
        return None

    def start_ccf_lookup(self):
        """Read the CCF coordinates of the run's cells and of the section maps' other cells (for their ML and DV), in a thread: slow the
        first time (the table of every cell of the atlas: a minute or so), seconds after. Only if they are on this computer: see ccf_action."""
        self.ccf_run = self.ccf_bg = None
        self.ccf_state = "none"
        self.ccf_started = True
        h = self.host
        if self.own is not None or self.abc_cache is None or "load_ccf_for_cells" not in h:
            return
        if core.ccf_files_missing(self.abc_cache):
            self.ccf_state = "missing"
            return
        self.ccf_state = "reading"
        names, bg = self.adata.obs_names, (self.bg["ids"] if self.bg is not None and self.bg.get("ids") is not None else [])
        # (the cells' coordinates are kept in the run's folder, the first time: the table behind them is the whole atlas. Which cells they are for is
        # told by their number and the first and last of them.)
        signature = f"{len(names)}|{names[0] if len(names) else ''}|{names[-1] if len(names) else ''}|{len(bg)}|{bg[0] if len(bg) else ''}|{bg[-1] if len(bg) else ''}"
        saved = Path(self.run_folder) / "ccf_cells.npz" if self.run_folder is not None else None

        def work():
            try:
                if saved is not None and saved.exists():
                    try:
                        with np.load(saved, allow_pickle=False) as z:
                            if str(z["signature"]) == signature:
                                self.ccf_run, self.ccf_bg = z["run"], (z["bg"] if z["bg"].size else None)
                                self.ccf_state = "ready"
                                return
                    except Exception:
                        pass                                            # (unreadable or not for these cells: read again)
                run_ids = list(self.adata.obs_names)                     # (here and not before the thread: for 500,000 cells and every cell of the sections that took 6 seconds with the window frozen)
                bg_ids = list(self.bg["ids"]) if self.bg is not None and self.bg.get("ids") is not None else []
                table, _shell = h["load_ccf_for_cells"](self.abc_cache, run_ids + bg_ids)
                xyz = table[["x_ccf", "y_ccf", "z_ccf"]].to_numpy(dtype=np.float32)       # (AP, DV, ML)
                self.ccf_run, self.ccf_bg = xyz[:len(run_ids)], (xyz[len(run_ids):] if bg_ids else None)
                self.ccf_state = "ready"
                if saved is not None:
                    try:
                        np.savez(saved, run=self.ccf_run, bg=self.ccf_bg if self.ccf_bg is not None else np.zeros((0, 3), np.float32), signature=np.array(signature))
                    except Exception:
                        pass                                            # (a folder that cannot be written to: it is read again next time)
            except Exception as e:
                self.ccf_state = "failed"
                print(f"Could not read the cells' CCF coordinates (for their ML and DV): {e}", flush=True)

        threading.Thread(target=work, daemon=True).start()

    def download_ccf_coords(self):
        # (said when the reconstructed coordinates are in use: a table of much the same size, downloaded earlier, that this is easily taken for)
        note = ("This is a separate table from the reconstructed coordinates already downloaded: those place the cells on the section "
                "maps, but are not CCF positions.") if core.reconstructed_in_use() else ""
        if self.host["ensure_ccf_downloaded"](self.abc_cache, parent=self.root, what="CCF coordinates of every cell (for their ML and DV)",
                                              note=note):
            self.start_ccf_lookup()

    def cell_info(self, i, source):
        return {"key": ("run", i), "umap_idx": i, "bg_idx": None, "sec": self.section_of[i] if self.section_of is not None else None,
                "x": self.cx[i] if self.cx is not None else np.nan, "y": self.cy[i] if self.cy is not None else np.nan, "source": source}

    def bg_cell_info(self, sec, g):
        """A cell of the section that is not in the run (so not in the UMAP)."""
        return {"key": ("bg", g), "umap_idx": None, "bg_idx": g, "sec": sec, "x": self.bg["x"][g], "y": self.bg["y"][g], "source": "section"}

    def bold_legend_row(self, canvas, k):
        """Show the legend row k of `canvas` in bold (None: no row); the row bold before goes back to normal."""
        old = getattr(self, "bold_row", None)
        if old is not None and old == (canvas, k):
            return
        if old is not None:
            try:
                old[0].itemconfigure(f"row{old[1]}", font=old[0].legend_fonts[0])
            except (tk.TclError, AttributeError):
                pass
        self.bold_row = None
        if k is None or not hasattr(canvas, "legend_fonts"):
            return
        canvas.itemconfigure(f"row{k}", font=canvas.legend_fonts[1])
        self.bold_row = (canvas, k)

    def bold_group_of(self, info):
        """The group list's row (or Single group(s) legend's) of the group the cell `info` belongs to, in bold; none if there isn't one."""
        canvas = k = None
        level = self.level
        if info is not None and self.mode in ("all", "single") and self.info.get(level) is not None:
            i = info["umap_idx"]
            if i is not None:
                label = self.info[level][i]
            else:
                arr = self.bg["levels"].get(level) if self.bg else None
                label = arr[info["bg_idx"]] if arr is not None else None
            if label is not None and str(label) != "nan":
                if self.mode == "all":
                    rows = getattr(self, "group_rows", [])
                    if str(label) in rows:
                        canvas, k = self.group_canvas, rows.index(str(label))
                else:
                    ids = (self.last_ids or {}).get("targets", [])
                    m = re.match(r"^\s*(\d+)", str(label))
                    if len(ids) > 1 and not self.level_has_ids.get(self.level, True) and str(label) in ids:
                        canvas, k = self.legend, ids.index(str(label))
                    elif len(ids) > 1 and m and int(m.group(1)) in ids:
                        canvas, k = self.legend, ids.index(int(m.group(1)))
        self.bold_legend_row(canvas, k)

    def show_info(self, info):
        """Ring the cell where the pointer is, and mirror it in the other plot (its description is in its tip: see hover)."""
        self.bold_group_of(info)
        self.mark_cell_on_histograms(info)
        if info is None:
            self.umap.set_marker("hover", None, None)
            for t in self.tiles:
                t["pc"].set_marker("hover", None, None)
            return
        src = info["source"]
        i = info["umap_idx"]
        if i is None:
            self.umap.set_marker("hover", None, None)
        else:
            self.umap.set_marker("hover", *self.coords[i], style="ring" if src == "umap" else "mirror")
        finite = info["x"] is not None and np.isfinite(info["x"]) and np.isfinite(info["y"])
        for t in self.tiles:
            if t["label"] == info["sec"] and finite:
                t["pc"].set_marker("hover", info["x"], info["y"], style="ring" if src == "section" else "mirror")
            else:
                t["pc"].set_marker("hover", None, None)

    def drag_state(self, on):
        """A plot, a section map or the zoomed pane is being dragged (or the drag has ended): no cell tip comes up while it is, and the wait for
        one starts again when it stops (at the cell then under the pointer)."""
        if on:
            self.cell_tip.suspend()
            return
        self.cell_tip.resume()
        info = getattr(self, "last_hover_info", None)
        if info is not None:
            self.cell_tip.pointer_on(info["key"], lambda: self.cell_rows(info), self.ccf_action, lambda: self.swatch_with_token(info))

    def hover(self, info):
        """The pointer is on a cell (info) or on none (None). The ring and description are immediate; once it rests, the plots dim and, a
        little later, the cell's group is marked with red '+' signs (only those two are turned off by the Hover highlight box). After a longer
        rest the description comes up in a tip as well (CellTip), whether or not the box is ticked."""
        self.last_hover_info = info
        if info is None:
            self.cancel_hover_jobs()
            self.hover_key = None
            self.cell_tip.pointer_off()
            self.show_info(None)
            self.clear_highlights()
            return
        self.show_info(info)
        self.cell_tip.pointer_on(info["key"], lambda: self.cell_rows(info), self.ccf_action, lambda: self.swatch_with_token(info))
        if not self.hover_on():
            self.hover_key = None                    # (so that ticking the box again marks the next cell pointed at)
            return
        if info["key"] == self.hover_key:
            return
        self.cancel_hover_jobs()
        self.clear_highlights(keep_dim=True)
        self.hover_key = info["key"]
        self.hover_job = self.root.after(HOVER_HOLD_MS, lambda: self.hover_settled(info))
        self.group_job = self.root.after(GROUP_HOLD_MS, lambda: self.show_group(info))

    def cancel_hover_jobs(self):
        for name in ("hover_job", "group_job"):
            j = getattr(self, name, None)
            if j is not None:
                self.root.after_cancel(j)
                setattr(self, name, None)

    def hover_settled(self, info):
        self.hover_job = None
        self.set_dim(True)

    def set_dim(self, on):
        self.umap.set_dim(on)
        for t in self.tiles:
            t["pc"].set_dim(on)

    def clear_highlights(self, keep_dim=False):
        if not keep_dim:
            self.set_dim(False)
        self.umap.set_pluses(None)
        self.umap.boost_label(None)
        for t in self.tiles:
            t["pc"].set_pluses(None)

    def show_group(self, info):
        """Red '+' signs on every cell of the hovered cell's group (at the Group level), in the UMAP and in every section."""
        self.group_job = None
        level = self.level
        if self.info.get(level) is None:
            return
        i = info["umap_idx"]
        if i is not None:
            label = self.info[level][i]
        else:
            arr = self.bg["levels"].get(level) if self.bg else None
            if arr is None:
                return
            label = arr[info["bg_idx"]]
        if label is None or str(label) == "nan":
            return
        self.mark_group(label, skip=i)

    def swatch_with_token(self, info):
        """cell_swatch, with the cell itself added: the tip lights the cell on the plots while the pointer is on its picture."""
        got = self.cell_swatch(info)
        return None if got is None else tuple(got) + (("cell", info), lambda: self.bring_cell_into_view(info), self.umap_shown)

    def bring_cell_into_view(self, info):
        """The tip's Bring into view: the cell is made visible in the UMAP and in its section's map, in whichever it is out of view (centered where it
        can be); where it is in view nothing changes. Not while the UMAP is hidden (the button is greyed out then; a tip kept up may still have it)."""
        if not self.umap_shown:
            return
        i = info["umap_idx"]
        if i is not None:
            self.umap.bring_into_view(*self.coords[i])
        finite = info["x"] is not None and np.isfinite(info["x"]) and np.isfinite(info["y"])
        tile = next((t for t in self.tiles if t["label"] == info["sec"]), None)
        if tile is None or not finite:
            return
        pc, fr = tile["pc"], tile["frame"]
        if pc.vsize is None and pc.panel_hook is None:                 # (a map zoomed on its own: the cell is brought into the map's view first)
            pc.bring_into_view(info["x"], info["y"])
        px, py = pc.position_in_widget(info["x"], info["y"])
        if pc.vsize is not None and getattr(fr, "slot_wh", None) is not None:        # (the whole panel zoomed: where the map is in the pane, from the layout)
            ax = fr.tile_x + px + pc.off[0]
            ay = fr.tile_y + (fr.slot_wh[1] - pc.vsize[1]) + py + pc.off[1]
        else:                                                           # (the panels at their size: the map is under its caption in its frame)
            ax = getattr(fr, "tile_x", fr.winfo_x()) + px
            ay = getattr(fr, "tile_y", fr.winfo_y()) + (tile["caption"].winfo_reqheight() if tile["caption"].winfo_manager() else 0) + py
        c = self.tile_scroll
        region = [float(v) for v in c.cget("scrollregion").split()]
        vw, vh = c.winfo_width(), c.winfo_height()
        left, top, margin = c.canvasx(0), c.canvasy(0), 6
        if margin <= ax - left <= vw - margin and margin <= ay - top <= vh - margin:
            return                                                     # (in view in the pane: nothing to do)
        if region[2] > vw:
            c.xview_moveto(min(max(ax - vw / 2, 0.0), region[2] - vw) / region[2])
        if region[3] > vh:
            c.yview_moveto(min(max(ay - vh / 2, 0.0), region[3] - vh) / region[3])

    def show_cell_dot(self, info):
        """The exact cell the tip is about, with the white dot of the other plot's copy of a hovered cell (bigger than a red +), in the UMAP and
        in its section's map."""
        self.mark_cell_on_histograms(info)
        i = info["umap_idx"]
        if i is None:
            self.umap.set_marker("hover", None, None)
        else:
            self.umap.set_marker("hover", *self.coords[i], style="mirror")
        finite = info["x"] is not None and np.isfinite(info["x"]) and np.isfinite(info["y"])
        for t in self.tiles:
            if t["label"] == info["sec"] and finite:
                t["pc"].set_marker("hover", info["x"], info["y"], style="mirror")
            else:
                t["pc"].set_marker("hover", None, None)

    def tip_row_hover(self, token):
        """The pointer is on a group's row in a cell tip (token: (level, label)), or left it (None): that group's cells are marked with red '+'
        signs and the rest dimmed, as for Hover highlight, whether or not that is on."""
        self.cancel_hover_jobs()
        self.hover_key = None                       # (so that the cell under the pointer later is marked again, if Hover highlight is on)
        self.clear_highlights()
        if token is None:
            self.show_info(None)                    # (the dot of a cell lit from the tip goes)
            return
        if token[0] == "cell":
            self.show_cell_dot(token[1])
            return
        if token[0] == "section":
            done = self.mark_section(token[1])
        elif token[0] == "region":
            done = self.mark_region(*token[1])
        else:
            done = self.mark_group(token[1], level=token[0])
        if done:
            self.set_dim(True)

    def region_ids(self, kind, key):
        """The numbers of the CCF regions that make up a region (itself), a structure or a division."""
        if kind == "region":
            return [key]
        column = {"structure": 2, "division": 3}[kind]
        return [number for number, names in (self.region_names or {}).items() if names[column] == key]

    def mark_region(self, kind, key):
        """Red '+' signs on every cell in a CCF region, structure or division (as Allen gives it): in the UMAP and in the section maps in view."""
        if self.region_of_run is None and not (self.bg and self.bg.get("region") is not None):
            return False
        ids = self.region_ids(kind, key)
        run_mask = np.isin(self.region_of_run, ids) if self.region_of_run is not None else np.zeros(len(self.coords), bool)
        if run_mask.any():
            self.umap.set_pluses((self.coords[run_mask, 0], self.coords[run_mask, 1]))
        else:
            self.umap.set_pluses(None)
        self.umap.boost_label(None)
        barr = self.bg.get("region") if self.bg else None
        found = bool(run_mask.any())
        for t in self.tiles:
            if not self.frame_visible(t["frame"]):
                t["pc"].set_pluses(None)
                continue
            if barr is not None:
                w = self.bg["where"][t["label"]]
                m = np.isin(barr[w], ids)
                xy = (self.bg["x"][w][m], self.bg["y"][w][m])
            elif self.region_of_run is not None:
                idx = self.run_idx[t["label"]]
                m = np.isin(self.region_of_run[idx], ids)
                xy = (self.cx[idx][m], self.cy[idx][m])
            else:
                xy = ((), ())
            t["pc"].set_pluses(xy if len(xy[0]) else None)
            found = found or len(xy[0]) > 0
        return found

    def mark_section(self, sec):
        """Red '+' signs on every cell of the section `sec`: in the UMAP (those of the run) and in that section's map (all its cells)."""
        if self.section_of is None:
            return False
        run_mask = self.section_of == sec
        if not run_mask.any():
            return False
        self.umap.set_pluses((self.coords[run_mask, 0], self.coords[run_mask, 1]))
        self.umap.boost_label(None)
        for t in self.tiles:
            if t["label"] != sec or not self.frame_visible(t["frame"]):
                t["pc"].set_pluses(None)
                continue
            if self.bg is not None and sec in self.bg["where"]:
                w = self.bg["where"][sec]
                xy = (self.bg["x"][w], self.bg["y"][w])
            else:
                idx = self.run_idx[sec]
                xy = (self.cx[idx], self.cy[idx])
            t["pc"].set_pluses(xy if len(xy[0]) else None)
        return True

    def mark_group(self, label, skip=None, level=None):
        """Red '+' signs on every cell of the group `label` (at `level`, the Group level if not given), in the UMAP and in every section. `skip`:
        a UMAP cell that keeps its own ring instead."""
        level = level or self.level
        run_mask = self.info[level] == label
        if not run_mask.any():
            return False                             # (a group none of the run's cells is in: nothing in the UMAP to relate to)
        if skip is not None:
            run_mask[skip] = False                   # (the hovered cell keeps its own ring)
        self.umap.set_pluses((self.coords[run_mask, 0], self.coords[run_mask, 1]))
        number = re.match(r"^\s*(\d+)", str(label))               # (the group's number is what the UMAP shows at its center: bigger while it is marked)
        self.umap.boost_label(number.group(1) if number and level == self.level else None)
        barr = self.bg["levels"].get(level) if self.bg else None
        for t in self.tiles:
            if not self.frame_visible(t["frame"]):
                t["pc"].pluses = None
                continue
            if barr is not None:
                w = self.bg["where"][t["label"]]
                m = barr[w] == label
                xy = (self.bg["x"][w][m], self.bg["y"][w][m])
            else:
                idx = self.run_idx[t["label"]]
                m = self.info[level][idx] == label
                xy = (self.cx[idx][m], self.cy[idx][m])
            t["pc"].set_pluses(xy if len(xy[0]) else None)
        return True

    def umap_motion(self, px, py):
        if px is None:
            self.hover(None)
            return
        i = self.umap.nearest(self.umap_tree, px, py)
        self.hover(None if i is None else self.cell_info(i, "umap"))

    def tile_motion(self, tile, px, py):
        if px is None:
            self.hover(None)
            return
        if tile["tree"] is None:
            sec = tile["label"]
            if self.bg is not None and sec in self.bg["where"]:         # (any cell of the section, in the run or not)
                w = self.bg["where"][sec]
                x, y = self.bg["x"][w], self.bg["y"][w]
                good = np.flatnonzero(np.isfinite(x) & np.isfinite(y))
                tile["good"] = w[good]
                tile["umap_of"] = self.run_index.get_indexer(self.bg["ids"][tile["good"]])
                pts = np.column_stack([x[good], y[good]])
            else:
                idx = self.run_idx[sec]
                good = np.flatnonzero(np.isfinite(self.cx[idx]) & np.isfinite(self.cy[idx]))
                tile["good"] = idx[good]
                tile["umap_of"] = tile["good"]
                pts = np.column_stack([self.cx[idx][good], self.cy[idx][good]])
            tile["tree"] = cKDTree(pts) if len(good) else False
        if tile["tree"] is False:
            return
        k = tile["pc"].nearest(tile["tree"], px, py)
        if k is None:
            self.hover(None)
            return
        u = int(tile["umap_of"][k])
        self.hover(self.cell_info(u, "section") if u >= 0 else self.bg_cell_info(tile["label"], int(tile["good"][k])))

    # ---- Settings and the UMAP parameters

    def apply_hover_enabled(self):
        """The Hover highlight box above the section maps: the setting of the same name in Settings, saved as it is when that is. Off: no
        dimming and no red '+' signs; the hovered cell's ring and its description in its tip stay. In "Browse All Neurons" (which
        starts with it off) it is for this window only, so the saved setting, for the other projects, is not changed."""
        if not self.start_panel_zoom:
            self.viewer_settings["hover_enabled"] = bool(self.var_hover.get())
            if "save_viewer_settings" in self.host:
                self.host["save_viewer_settings"](self.viewer_settings)
        self.update_marker_slider()
        if not self.hover_on():
            self.drop_highlights()

    def drop_highlights(self):
        """Hover highlight was turned off: the dimming and '+' signs of the cell under the pointer go (its ring and description stay)."""
        self.cancel_hover_jobs()
        self.hover_key = None
        self.clear_highlights()

    def open_settings(self):
        dialog = self.host.get("settings_dialog")
        if dialog is None or self.busy:
            return
        new = dialog(self.root, self.viewer_settings, umap_params=self.params, n_cells=self.n)
        if new is None:
            return
        new_params = new.pop("umap_params", None)        # (belong to this run, not to the user's preferences)
        shape_changed = new.get("gene456_shape") != self.viewer_settings.get("gene456_shape")
        size_changed = new.get("umap_circle_size") != self.viewer_settings.get("umap_circle_size")
        hover_changed = bool(new.get("hover_enabled", True)) != bool(self.viewer_settings.get("hover_enabled", True))
        self.viewer_settings.update(new)
        self.host["save_viewer_settings"](self.viewer_settings)
        if hover_changed or not self.start_panel_zoom:          # (Browse All Neurons: the box keeps its own state unless the setting itself changed)
            self.var_hover.set(bool(self.viewer_settings.get("hover_enabled", True)))
        self.update_marker_slider()
        if not self.hover_on():
            self.drop_highlights()
        if size_changed:                                            # (the UMAP's Cell size slider shows it: set without its own redraw, done below)
            self.umap_sliders[0].set(self.umap_size_index(self.viewer_settings.get("umap_circle_size", "Default")))
        if (size_changed or (shape_changed and self.mode in ("gene", "imputed"))) and not self.busy:
            self.recolor()
        if new_params is not None and self.params is not None:
            self.switch_variant(new_params)

    def switch_variant(self, new_params):
        """Show the UMAP for these parameters: from memory, else from the run's cache, else by calculating it (progress window)."""
        h = self.host
        key = h["umap_params_key"](new_params)
        if key == h["umap_params_key"](self.params):
            return
        entry = self.variants.get(key)
        message = None
        if entry is None:
            folder = self.run_folder if self.run_folder is not None else h["CACHE_DIR"]
            self.var_status.set("Working... recalculating the UMAP; the viewer is locked until it finishes.")
            self.root.config(cursor="watch")
            self.root.update()
            try:
                csv_path, message = h["ensure_umap_variant_csv"](self.adata, new_params, folder)
                if csv_path is None:
                    self.var_status.set(message or "The UMAP was not recalculated.")
                    return
                try:
                    new_coords, new_leiden = h["load_umap_variant"](csv_path, self.adata.obs_names)
                except Exception as e:
                    self.var_status.set(f"Could not read the recalculated UMAP: {e}")
                    return
            finally:
                self.root.config(cursor="")
            entry = (np.asarray(new_coords, dtype=float)[:, :2], None if new_leiden is None else np.asarray(new_leiden, dtype=object).astype(str))
            self.variants[key] = entry
        self.install_variant(entry[0], entry[1], new_params)
        self.var_status.set(message or f"UMAP shown with {h['umap_params_label'](new_params)}.")

    def install_variant(self, new_coords, new_leiden, new_params):
        """Swap the embedding and its Leiden clustering in place and redraw everything derived from either."""
        h = self.host
        self.hover(None)
        self.coords = np.asarray(new_coords, dtype=float)[:, :2].copy()
        if new_leiden is not None:
            self.adata.obs[LEIDEN_KEY] = pd.Categorical(new_leiden)
            self.info["leiden"] = np.asarray(new_leiden).astype(str)
        else:
            self.adata.obs.drop(columns=[LEIDEN_KEY], inplace=True, errors="ignore")
            self.info["leiden"] = None
        for k in [k for k in self.cache if isinstance(k, tuple) and k[1] == "leiden"]:
            del self.cache[k]
        self.levels = [k for k in self.level_names if self.info[k] is not None]
        if self.level not in self.levels:
            self.level = self.levels[0]
            self.var_level.set(self.level)
        self.build_level_radios()
        c = self.coords
        padx = (c[:, 0].max() - c[:, 0].min()) * 0.05 or 1.0
        pady = (c[:, 1].max() - c[:, 1].min()) * 0.05 or 1.0
        self.umap.home = (c[:, 0].min() - padx, c[:, 0].max() + padx, c[:, 1].min() - pady, c[:, 1].max() + pady)
        self.umap.fit()
        self.umap_tree = cKDTree(c)
        self.params = h["normalize_umap_params"](new_params)
        self.root.title(self.window_title())
        if self.run_folder is not None:
            h["write_last_umap_params"](self.run_folder, self.params)
        self.recolor()

    def window_title(self):
        title = "Interactive UMAP viewer"
        h = self.host
        if self.params is not None and not h["umap_params_are_default"](self.params):
            title = f"{title} ({h['umap_params_label'](self.params)})"
        if self.run_folder is not None:
            title = f"{Path(self.run_folder).name} \u2014 {title}"
        return title

    def open_data_folder(self):
        folder = Path(self.run_folder) if self.run_folder is not None else Path(self.host["CACHE_DIR"])
        folder.mkdir(parents=True, exist_ok=True)
        self.host["open_with_default_viewer"](folder)

    # ---- closing

    def close(self):
        if self.gc_was_enabled:
            gc.enable()
        try:
            self.root.after_cancel(self.poll_id)
        except Exception:
            pass
        if self.own_root:
            try:
                self.root.quit()
            except tk.TclError:
                pass
        for win in list(self.subset_windows):
            try:
                win.close()
            except Exception:
                pass
        if self.layout_job is not None:
            try:
                self.root.after_cancel(self.layout_job)
            except Exception:
                pass
        core.release_tk_variables(self)
        self.root.destroy()

    def run(self):
        try:
            self.root.state("zoomed")
        except tk.TclError:
            pass
        self.root.lift()
        self.root.focus_force()
        if self.own_root:
            self.root.mainloop()
        else:
            self.root.wait_window(self.root)


def show_umap_window_tk(adata, abc_cache, imputed_state=None, adata_backed=None, rois=None, run_folder=None, ensure_imputed=None, title=None,
                        parent=None, host=None, params=None, own=None, panel_zoom=False):
    """Open the viewer and return when it is closed. `panel_zoom`: start with Zoom whole panel ticked. `host`: the browser's own functions (see UmapWindow); `params`: the UMAP parameters
    the embedding in `adata` was made with."""
    if "X_umap" not in adata.obsm and not adata.uns.get("no_umap"):            # (a project made without a UMAP has none, and is shown without it)
        print("No UMAP embedding available; skipping the interactive UMAP viewer.")
        return
    win = UmapWindow(adata, abc_cache, imputed_state, adata_backed, rois, run_folder, ensure_imputed, title or "Interactive UMAP viewer", parent,
                     host=host, params=params, own=own, panel_zoom=panel_zoom)
    win.run()
