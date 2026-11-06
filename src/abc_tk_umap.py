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
from pathlib import Path
import threading
import time
import types
import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk

import cv2
import numpy as np
import pandas as pd
from PIL import Image, ImageTk
from scipy.spatial import cKDTree

import abc_core as core
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
MIN_TILE_PX = 50
TILE_COLUMNS = (6, 5, 4, 3, 2, 1)         # the Section size slider's stops: how many panels fit across the pane (1/6 of its width ... all of it)
OUTLINE_PX_LARGEST = 2                    # the CCF outlines' width (pixels) at the largest section map size (one pixel at the others)
DEFAULT_TILE_INDEX = 4                    # (two across)
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
HOVER_LABEL_PX = 110               # the cell tip's table (CellTip): its first column wraps here
HOVER_GAP_PX = 8                   # the space between the cell tip's columns
PREVIEW_ZOOM = True                # a wheel step on a plot shows the last picture scaled, and draws the view for real once the wheel is still
ZOOM_DRAW_MS = 120
PAN_DRAW_MS = 45                   # while a plot is dragged it is drawn again this often (in between, what is on the canvas is moved as it is)
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


def draw_marker(canvas, cx, cy, code, color, size=12):
    """A marker of this code, `size` pixels across, on a Tk canvas (a legend swatch)."""
    marker, is_open = SHAPE_MARKERS[code % 10], code >= 10
    h = size / 2
    style = {"fill": "", "outline": color, "width": 1.5} if is_open else {"fill": color, "outline": ""}
    if marker == 'o':
        canvas.create_oval(cx - h, cy - h, cx + h, cy + h, **style)
    else:
        canvas.create_polygon(*[c for x, y in SHAPE_POLYGONS[marker] for c in (cx + (x - 0.5) * size, cy + (y - 0.5) * size)], **style)


def stamp_plus(img, px, py, rgb, arm):
    """A '+' (two pixels thick, `arm` pixels each way) at each point, in place. `rgb`: one color (3,) or one per point (n, 3)."""
    h, w = img.shape[:2]
    per_point = np.ndim(rgb) == 2
    for k in range(-arm, arm + 1):
        for dx, dy in ((k, 0), (k, 1), (0, k), (1, k)):
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


def percent_below(vals, point):
    """The percentage of the cells (the finite values) below `point`, and of those that express the gene (a value above 0) that are:
    (all cells, expressing cells)."""
    fin = np.isfinite(vals)
    n = int(np.count_nonzero(fin))
    if n == 0:
        return 0.0, 0.0
    below = fin & (vals < point)
    expr = fin & (vals > 0)
    ne = int(np.count_nonzero(expr))
    return 100.0 * np.count_nonzero(below) / n, (100.0 * np.count_nonzero(below & expr) / ne if ne else 0.0)


class RangeSlider(tk.Canvas):
    """A track with two handles, for the lowest and the highest end of a range (0 to 1): the black level (left, a black dot) and the white level
    (right, a white dot) of a color scale. It is as wide as the color bar above it (the track runs between the same margins), calls
    `command(low, high)` as a handle is dragged, and a double-click puts both back at the ends. `low` and `high` are the fractions of the range
    the levels are at; the track is drawn by their square root (a position is the square root of its value), so that the low end, where most of a
    gene's expression is, has most of the track."""
    MIN_GAP = 0.02            # the handles keep at least this much of the track between them

    def __init__(self, parent, command, margin=10, height=26):
        super().__init__(parent, height=height, highlightthickness=0)
        self.low, self.high, self.margin, self.command, self.dragging = 0.0, 1.0, margin, command, None
        self.bind("<Configure>", lambda e: self.draw())
        self.bind("<ButtonPress-1>", self.press)
        self.bind("<B1-Motion>", self.move)
        self.bind("<ButtonRelease-1>", lambda e: setattr(self, "dragging", None))
        self.bind("<Double-Button-1>", lambda e: self.set(0.0, 1.0, notify=True))

    def x_of(self, fraction):
        return self.margin + math.sqrt(min(max(fraction, 0.0), 1.0)) * max(self.winfo_width() - 2 * self.margin, 1)

    def fraction_of(self, x):
        position = min(max((x - self.margin) / max(self.winfo_width() - 2 * self.margin, 1), 0.0), 1.0)
        return position * position

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
        if self.tip is None:
            self.tip = tk.Toplevel(self.widget)
            self.tip.wm_overrideredirect(True)
            self.label = tk.Label(self.tip, bg="#ffffe0", relief="solid", borderwidth=1, font=("Segoe UI", 9), padx=4, pady=1, justify=tk.LEFT, wraplength=460)
            self.label.pack()
        self.label.config(text=self.text_fn())
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


class CellTip:
    """The hovered cell's description in a tip of its own, once the pointer has rested on one cell for
    CELL_TIP_MS, with a button that copies it. Unlike a HoverTip it stays while the pointer goes onto it (to reach the button): it goes when
    the pointer is on another cell more than CELL_TIP_GRACE_PX from where the tip came, or leaves the tip itself."""

    def __init__(self, root, delay_ms=None):
        self.root, self.tip, self.job, self.key, self.anchor = root, None, None, None, None
        self.delay_ms = CELL_TIP_MS if delay_ms is None else delay_ms

    def pointer_on(self, key, rows_fn, action_fn=None):
        """The pointer is on a new cell (`key`); `rows_fn()` gives its description, as rows (what, number, name), asked when the tip is shown.
        `action_fn()`, if given, is asked then too: None, or (text, command) for a button beside Copy."""
        if key == self.key:
            return
        self.key = key
        self.cancel()
        if self.tip is not None and not self.pointer_in_tip() and not self.near_anchor():
            self.hide()
        self.job = self.root.after(self.delay_ms, lambda: self.show(rows_fn, action_fn))

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
        if not self.pointer_in_tip():
            self.hide()

    def show(self, rows_fn, action_fn=None):
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
        # in columns: what (gray), the number (bold, to the right), the name
        grid = tk.Frame(box, bg=bg)
        grid.pack(anchor="w", padx=6, pady=(4, 2))
        for r, (what, num, name) in enumerate(rows):
            tk.Label(grid, text=what, bg=bg, fg="gray35", font=("Segoe UI", 9, "bold"), anchor="w", justify=tk.LEFT,
                     wraplength=HOVER_LABEL_PX).grid(row=r, column=0, sticky="nw", padx=(0, HOVER_GAP_PX))
            tk.Label(grid, text=num, bg=bg, font=("Segoe UI", 9, "bold"), anchor="e").grid(row=r, column=1, sticky="ne", padx=(0, HOVER_GAP_PX))
            tk.Label(grid, text=name, bg=bg, font=("Segoe UI", 9), anchor="w", justify=tk.LEFT, wraplength=300).grid(row=r, column=2, sticky="nw")
        # copied with a tab between the columns: pasted into a spreadsheet, they stay columns
        text = "\n".join("\t".join(row) for row in rows)
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
        if self.tip is not None:
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
        self.rects = []                             # [item, x0, y0, x1, y1]
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

    def set_view(self, x0, y0, w, h):
        """The part of a virtual canvas the widget shows, and draws."""
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
        self.draw()

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

    # ---- drawing

    def set_layers(self, layers):
        self.layers = layers
        self.request_draw()

    def request_draw(self):
        if self.pending is None:
            self.pending = self.canvas.after_idle(self.draw)

    def draw(self):
        self.pending = None
        if self.pan_job is not None:
            self.canvas.after_cancel(self.pan_job)
            self.pan_job = None
        if self.zoom_job is not None:
            self.canvas.after_cancel(self.zoom_job)
            self.zoom_job = None
        w, h = self.size()
        if w <= 2 or h <= 2:
            return
        if self.visible_fn is not None and not self.visible_fn():
            self.stale = True
            return
        self.stale = False
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
        self.drawn = (ox, oy, ex, ey)
        W, H = w, h                                                           # (the whole canvas; w, h below are the picture's)
        w, h = ex - ox, ey - oy
        img = np.empty((h, w, 3), np.uint8)
        img[:] = self.bgcolor
        # the data in view: the picture's pixels (ox..ex, oy..ey of the canvas) in data units
        x_lo, x_hi = self.cx + (ox - W / 2) / self.scale, self.cx + (ex - W / 2) / self.scale
        d_lo, d_hi = (oy - H / 2) / self.scale, (ey - H / 2) / self.scale
        y_lo, y_hi = (self.cy - d_hi, self.cy - d_lo) if self.flip_y else (self.cy + d_lo, self.cy + d_hi)
        opacity = 1.0 if self.moving else self.opacity                        # (solid while the view moves: blending translucent dots is much slower; see note_motion)
        for lay in self.layers:
            if "size_mm" in lay:                                              # (a size on the tissue: the same at any window size)
                d = int(max(round(lay["size_mm"] * self.dot_factor * self.scale), 1))
            else:
                d = int(max(round(lay["size"] * (zoom * self.panel_zoom) ** 0.7), 1))          # (the dots grow as the view, or the whole panel, is zoomed in)
            x, y = lay["x"], lay["y"]
            pm = d / self.scale
            m = (x >= x_lo - pm) & (x <= x_hi + pm) & (y >= y_lo - pm) & (y <= y_hi + pm)
            idx = np.flatnonzero(m)
            if not len(idx):
                continue
            px = np.floor((x[idx] - self.cx) * self.scale + W / 2 - ox).astype(np.int32)      # (floor: the same pixel as on the whole canvas, also left of or above the picture)
            dy = (y[idx] - self.cy) * self.scale
            py = np.floor((-dy if self.flip_y else dy) + H / 2 - oy).astype(np.int32)
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
        if self.cell_alpha < 0.999:
            img = (img.astype(np.float32) * self.cell_alpha + np.asarray(self.bgcolor, np.float32) * (1.0 - self.cell_alpha)).astype(np.uint8)
        if self.outlines:
            self.draw_outlines(img, W, H, ox, oy, (x_lo, x_hi, y_lo, y_hi))
        self.base = img
        self.compose()
        self.reposition()
        self.shown_view, self.shown_scale, self.previewed = (self.cx, self.cy), self.scale, (0, 0)          # (what is on the canvas is for this view: see preview_pan)
        self.remember_view()
        if self.on_view:
            self.on_view()

    def set_outlines(self, polygons):
        self.outlines = polygons or None
        self.request_draw()

    def draw_outlines(self, img, w, h, ox, oy, view):
        """The outlines in view, one pixel wide, dark on a light background and light on a dark one. (w, h: the whole canvas; ox, oy: where the
        picture starts in it; view: the data x0, x1, y0, y1 it shows.)"""
        x0, x1, y0, y1 = view
        shown = []
        for p in self.outlines:
            lo, hi = p.min(axis=0), p.max(axis=0)
            if hi[0] < x0 or lo[0] > x1 or hi[1] < y0 or lo[1] > y1:
                continue
            px = (p[:, 0] - self.cx) * self.scale + w / 2 - ox
            dy = (p[:, 1] - self.cy) * self.scale
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
            zoom = self.scale / self.fit_scale()
            stamp_plus(img, px, py, np.array([255, 0, 0], np.uint8), max(int(min((self.plus_arm + zoom ** 0.4) * self.marker_factor, 14 * self.marker_factor)), 1))
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

    def add_marker(self, name, radius=8):
        item = self.canvas.create_oval(0, 0, 0, 0, outline=RING, width=2, state="hidden")
        self.markers[name] = [item, None, None, radius]

    MARKER_STYLES = {"ring": ("", "#e00000"), "mirror": ("white", "black")}      # (fill, outline): the cell under the pointer; the same cell in the other plot

    def set_marker(self, name, x, y, style="ring"):
        m = self.markers[name]
        m[1], m[2] = x, y
        fill, outline = self.MARKER_STYLES[style]
        self.canvas.itemconfig(m[0], fill=fill, outline=outline)
        self.place_marker(m)

    def place_marker(self, m):
        if m[1] is None:
            self.canvas.itemconfig(m[0], state="hidden")
            return
        px, py = self.to_canvas(m[1], m[2])
        r = m[3] * self.marker_factor
        self.canvas.coords(m[0], px - r, py - r, px + r, py + r)
        self.canvas.itemconfig(m[0], state="normal")
        self.canvas.tag_raise(m[0])

    def add_rect(self, x0, y0, x1, y1):
        item = self.canvas.create_rectangle(0, 0, 0, 0, outline="#ff8c00", width=2, dash=(6, 3))
        self.rects.append([item, x0, y0, x1, y1])

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
        for it, x0, y0, x1, y1 in self.rects:
            a, b = self.to_canvas(x0, y0)
            c, d = self.to_canvas(x1, y1)
            self.canvas.coords(it, a, b, c, d)
            self.canvas.tag_raise(it)
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
        """Index (into the tree's points) of the point within HOVER_RADIUS_PX of the canvas position, or None."""
        x, y = self.to_data(px, py)
        dist, k = tree.query([x, y], distance_upper_bound=HOVER_RADIUS_PX / self.scale)
        return None if not np.isfinite(dist) else int(k)


# ---------------------------------------------------------------------------------------------------------------------- the window

class UmapWindow(ExportsMixin):
    def __init__(self, adata, abc_cache, imputed_state, adata_backed, rois, run_folder, ensure_imputed, title, parent=None, host=None, params=None,
                 own=None, panel_zoom=False):
        self.start_panel_zoom = panel_zoom                  # start with Zoom whole panel ticked (and the section maps at their smallest)
        self.adata, self.abc_cache = adata, abc_cache
        self.own = own                      # a dataset of the user's own (abc_tk_own.OwnDataset): no atlas, no sections, no imputed genes
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
        loading =ttk.Label(self.root, text="Preparing the UMAP viewer...", font=("Segoe UI", 14), padding=40)
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
        self.ccf_run = self.ccf_bg = None  # the cells' CCF coordinates (AP, DV, ML), for the Paxinos ML and DV of the cell tip: see start_ccf_lookup
        self.ccf_state = "none"
        self.variants = {}                 # UMAP parameter key -> (coordinates, Leiden labels), the ones already shown
        self.load_data(adata_backed, rois)
        if self.host and self.params is not None:
            self.variants[self.host["umap_params_key"](self.params)] = (self.coords.copy(), self.info.get("leiden"))
        loading.destroy()
        self.build()
        self.recolor()
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.bind("<Escape>", lambda e: self.close())
        self.poll()
        self.root.after(1500, self.preload_imputed)         # (once the window has been drawn)
        self.root.after(2500, self.start_ccf_lookup)

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
        self.coords = np.asarray(a.obsm["X_umap"], dtype=float)[:, :2]
        self.n = len(self.coords)
        self.has_counts = "counts" in a.layers or self.own is not None
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
        self.cell_tip = CellTip(self.root)
        self.pw = pw = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True)
        # the section panels
        left = ttk.Frame(pw)
        bars = ttk.Frame(left, padding=(6, 4))
        bars.pack(side=tk.TOP, fill=tk.X)
        # Zoom whole panel from the start ("Browse All Neurons"): the section maps start at their smallest, so that they are made once, not made and then again
        self.var_tile = tk.IntVar(bars, TILE_COLUMNS.index(max(TILE_COLUMNS)) if self.start_panel_zoom else DEFAULT_TILE_INDEX)
        self.tile_stop_before_zoom = DEFAULT_TILE_INDEX if self.start_panel_zoom else None
        self.var_dot = tk.IntVar(bars, DEFAULT_DOT_INDEX)
        self.var_marker = tk.IntVar(bars, DEFAULT_MARKER_INDEX)
        self.dot_factors = [um / (DOT_MM * 1000) for um in CELL_SIZES_UM]          # (times DOT_MM)
        marker_names = ("Small", "Medium", "Large")
        self.sliders = [
            self.add_slider(bars, 0, "Section size:", self.var_tile, len(TILE_COLUMNS), self.layout_tiles,
                            lambda i: "full width" if TILE_COLUMNS[i] == 1 else f"1/{TILE_COLUMNS[i]} width", None),
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
        roi_only_box.grid(row=1, column=4, sticky="w", padx=(12, 6))
        if self.rois:
            HoverTip(roi_only_box, lambda: "The cells outside the ROIs go gray, on a white background: only the ROIs are colored. Saved section maps "
                                           "have the same look (and, as an SVG, the gray cells are pictures, so the file stays small).")
        else:
            roi_only_box.state(["disabled"])
            HoverTip(roi_only_box, lambda: "This project has no ROIs.")
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
        self.var_panel_zoom = tk.BooleanVar(bars, self.start_panel_zoom)
        panel_zoom_box = ttk.Checkbutton(bars, text="Zoom whole panel", variable=self.var_panel_zoom, command=self.apply_panel_zoom_mode)
        panel_zoom_box.grid(row=2, column=4, sticky="w", padx=(12, 6))
        HoverTip(panel_zoom_box, lambda: "Off: the wheel zooms into each section map, which stays in its panel (clipped by its border). "
                                         "On: the wheel zooms the whole panel, the borders and the sections together, and a drag moves the pane; "
                                         "sections scrolled out of view are not drawn. Double-click puts the panel back.")
        # (the same setting as in the Settings window; but off to begin with for "Browse All Neurons", whose many small maps make it distracting,
        # and there the box is for this window alone: see apply_hover_enabled)
        self.var_hover = tk.BooleanVar(bars, False if self.start_panel_zoom else bool(self.viewer_settings.get("hover_enabled", True)))
        hover_box = ttk.Checkbutton(bars, text="Hover highlight", variable=self.var_hover, command=self.apply_hover_enabled)
        hover_box.grid(row=0, column=4, sticky="w", padx=(12, 6))
        HoverTip(hover_box, lambda: "When the pointer rests on a cell, dim the plots and mark every cell of its group with a red +. Off: no "
                                    "marking, but the cell is still ringed and described in its tip (the same setting as in Settings).")
        self.update_marker_slider()
        bars.columnconfigure(6, weight=1)          # (the spare width goes to an empty column on the right: the controls stay where they are)
        self.tile_scroll = tk.Canvas(left, highlightthickness=0, bg="#f4f4f4")
        sb = ttk.Scrollbar(left, orient=tk.VERTICAL, command=self.tile_scroll.yview)
        self.hsb = ttk.Scrollbar(left, orient=tk.HORIZONTAL, command=self.tile_scroll.xview)        # (shown only while the panel is zoomed in)
        self.tile_scroll.configure(yscrollcommand=lambda a, b: (sb.set(a, b), self.root.after_idle(self.refresh_visible)),
                                   xscrollcommand=lambda a, b: (self.hsb.set(a, b), self.root.after_idle(self.refresh_visible)))
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
        self.umap.on_view = None
        self.umap_tree = cKDTree(c)
        pw.add(mid, weight=3)
        side = ttk.Frame(pw, padding=10, width=420)
        side.pack_propagate(False)
        self.build_side(side)
        pw.add(side, weight=0)
        self.side_pane = side
        self.root.update_idletasks()
        self.fit_left_pane(False)                                        # (the panes at their final widths before the section maps are made: they are laid out once)
        self.build_tiles()
        self.root.after(400, self.fit_left_pane)                         # (and looked at again once the window has settled: a no-op if all is as it should be)

    def fit_left_pane(self, again=True):
        """The section pane gets about a third of the window (the panels fill whatever width it has). Done before the panels are made, so that they are made
        once; looked at again after the window has settled, which moves nothing when all is as it should be (a move is a second drawing of everything)."""
        if not self.root.winfo_viewable() or self.own is not None:
            return
        total = self.pw.winfo_width()
        side_w = int(self.side_pane.cget("width")) or 420
        edge = total - side_w                                           # (the side panel keeps its own width, at the right)
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
            ttk.Radiobutton(modes, text=text, value=value, variable=self.var_mode, command=self.on_mode).grid(row=k, column=0, sticky="w")
        self.var_translucent = tk.BooleanVar(f, False)
        translucent_box = ttk.Checkbutton(color_col, text="Translucent", variable=self.var_translucent, command=self.apply_opacity)
        translucent_box.pack(anchor="w", pady=(4, 0))                      # (under the Color by buttons, in the room beside the Level list)
        HoverTip(translucent_box, lambda: f"Draw the dots {TRANSLUCENT * 100:g}% opaque instead of fully opaque,\n"
                                          "so that where they overlap you see how many there are.\n\n"
                                          "Zooming and moving the plots\n"
                                          "become much slower while this is on.")
        self.level_frame = ttk.LabelFrame(self.top_row, text="Level", padding=6)
        self.level_frame.grid(row=0, column=1, sticky="nsew", padx=(8, 0))
        HoverTip(self.level_frame, self.level_tip)
        self.var_level = tk.StringVar(f, self.level)
        self.build_level_radios()
        self.box = ttk.Frame(f)
        self.box.pack(fill=tk.X, pady=(8, 0))
        # a group number box (Specific group(s))
        self.id_frame = ttk.LabelFrame(self.box, text="Group number(s)", padding=6)
        ttk.Label(self.id_frame, text="One or more, separated by commas (Enter):", foreground="gray30", wraplength=300).pack(anchor="w")
        self.var_ids = tk.StringVar(f, "")
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
        self.gene_levels = (0.0, 1.0)          # the black and white levels of a gene shown in linear scale, as fractions of its range (the slider under the color bar)
        self.multi_levels = {}                 # the same for each of several genes (gene number -> (low, high)): the sliders in the legend
        self.multi_sliders = []
        self.multi_windows = []
        self.levels_job = None
        self.cbar_slider = RangeSlider(self.gene_frame, self.on_gene_levels)
        HoverTip(self.cbar_slider, lambda: "Drag the black and white handles to set the lowest and the highest value of the color scale. The track is by the square root of the value: "
                                        "most of it is for the low values. Double-click to reset.")
        # the legend for Single group(s)
        self.legend = tk.Canvas(f, width=300, height=10, highlightthickness=0)
        self.select_bar = ttk.Frame(f)                      # All / None, above the group list and the Single group(s) legend
        all_button = ttk.Button(self.select_bar, text="All", width=6, style="Side.TButton", command=self.select_all_groups)
        all_button.pack(side=tk.LEFT)
        none_button = ttk.Button(self.select_bar, text="None", width=6, style="Side.TButton", command=self.select_no_groups)
        none_button.pack(side=tk.LEFT, padx=(6, 0))
        HoverTip(all_button, lambda: "Select all groups")
        HoverTip(none_button, lambda: "Deselect all groups")
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
        for k, name in enumerate(self.levels):
            rb = ttk.Radiobutton(self.level_frame, text=self.level_names[name], value=name, variable=self.var_level, command=self.on_level)
            rb.grid(row=k, column=0, sticky="w")                       # (one column, coarse to fine, top to bottom)
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
        self.level_frame.config(text="Level" if self.mode in ("all", "single") else "Level (hover)")
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
        shown = ttk.Label(parent, width=11)                # (room for the longest value, "full width" or "Very large", and no more)
        shown.grid(row=row, column=3, sticky="w", padx=(6, 0))
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
        self.set_slider_enabled(self.sliders[2], self.var_hover.get())

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
        for scale, default in zip(self.sliders, (DEFAULT_TILE_INDEX, DEFAULT_DOT_INDEX, DEFAULT_MARKER_INDEX)):
            scale.set(default)

    def layout_tiles_when_settled(self):
        """The section pane changed size (the divider or the window is being dragged): lay the panels out once it has stopped changing, not at
        every step, which redrew all of them each time."""
        if not self.layout_seen:                                    # (the first time, when the window opens, at once)
            self.layout_seen = True
            self.layout_tiles()
            return
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
        cols = TILE_COLUMNS[self.var_tile.get()]
        w = max(self.tile_scroll.winfo_width(), 100)
        px = int(max((w - 4) // cols - 8 - 2, MIN_TILE_PX) * self.panel_zoom)          # (8: the gaps either side of a panel; 2: its canvas border)
        if cols == self.tile_cols and px == self.tile_px:
            return
        self.tile_cols, self.tile_px = cols, px
        big = self.panel_zoom > 1.0
        ph = int(px / self.tile_aspect)
        for t in self.tiles:
            pc = t["pc"]
            pc.outline_px = OUTLINE_PX_LARGEST if cols == min(TILE_COLUMNS) else 1            # (the CCF outlines: wider at the largest size)
            if big:
                pc.set_virtual((px, ph))
            else:
                pc.set_virtual(None)
                pc.canvas.config(width=px, height=ph)
                t["frame"].slot_wh = None
                t["shown"] = True
                if not t["caption"].winfo_ismapped():
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
        self.tile_scroll.configure(scrollregion=(0, 0, 4 + cols * (fw + 8), 3 + rows * (fh + 6)))
        self.root.after_idle(self.refresh_visible)

    def place_tiles(self, max_draws=None):
        """While the whole panel is zoomed: each section map's widget is put where the part of the map that is on screen is (and made that big, a
        little more all round), and drawn there; a map not on screen has no widget shown. A map whose widget still covers what is on screen is left.
        `max_draws`: at most this many maps are drawn (those with the most on screen first); while the panel is being dragged, none unless it says.
        Returns how many more are waiting to be drawn."""
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
            if with_caption and not t["caption"].winfo_ismapped():
                t["caption"].pack(anchor="w", before=pc.canvas)
            elif not with_caption and t["caption"].winfo_ismapped():
                t["caption"].pack_forget()
            self.tile_scroll.coords(t["item"], cx0 + x0, cy0 + y0 - (ch if with_caption else 0))
            if not t.get("shown", True):
                self.tile_scroll.itemconfig(t["item"], state="normal")
                t["shown"] = True
            pc.set_view(int(x0), int(y0), int(x1) - int(x0), int(y1) - int(y0))
        if waiting[:max_draws]:
            # A map's widget moved, resized and given its new picture: Tk does those in passes of its own (the pane moves the widget, the frame packs
            # its caption and canvas again, the canvas draws), and a screen refresh between them showed a moment of the old picture in the new
            # place, a blink. All of them now, together, so that only the finished map reaches the screen.
            self.root.update_idletasks()
        return max(len(waiting) - max_draws, 0)

    # ---- zooming the whole panel (the borders and the sections together)

    MAX_PANEL_ZOOM = 100.0

    def max_panel_zoom(self):
        return self.MAX_PANEL_ZOOM

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
            self.tile_stop_before_zoom = self.var_tile.get()          # (the section maps start at their smallest, 1/6 of the width: to zoom into)
            self.sliders[0].set(TILE_COLUMNS.index(max(TILE_COLUMNS)))
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

    def set_panel_zoom(self, zoom, anchor=None):
        """Enlarge the whole panel to `zoom` times its size (1 to MAX_PANEL_ZOOM): the panels are laid out again, bigger, so that borders, gaps and
        sections grow together, and only those in view are drawn. `anchor` = (x, y in the scrolled pane, x, y in the window): the point that stays
        under the pointer."""
        if getattr(self, "panning", False):
            self.end_pan()                           # (a drag just before: its drawing put off no longer, as the maps are laid out and drawn here)
        zoom = min(max(zoom, 1.0), self.max_panel_zoom())
        if abs(zoom - self.panel_zoom) < 1e-9 or not self.tiles:
            return
        spot = self.spot_in_tiles(anchor[0], anchor[1]) if anchor is not None else None
        self.panel_zoom = zoom
        for t in self.tiles:
            t["pc"].panel_zoom = zoom
        self.layout_tiles()
        if zoom > 1.0:
            self.hsb.pack(side=tk.BOTTOM, fill=tk.X, before=self.tile_scroll)
        else:
            self.hsb.pack_forget()
        new = [float(v) for v in self.tile_scroll.cget("scrollregion").split()]
        if spot is not None and new[2] > 0 and new[3] > 0:
            nx, ny = self.spot_now(spot)
            wx, wy = anchor[2], anchor[3]
            self.tile_scroll.xview_moveto(min(max((nx - wx) / new[2], 0.0), 1.0))
            self.tile_scroll.yview_moveto(min(max((ny - wy) / new[3], 0.0), 1.0))
        self.place_tiles()
        self.root.after_idle(self.refresh_visible)

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
        if b is None:
            image = self.snapshot_pane()
            if image is None:
                self.set_panel_zoom(self.panel_zoom * factor, (self.tile_scroll.canvasx(wx), self.tile_scroll.canvasy(wy), wx, wy))
                return
            b = self.zoom_burst = {"image": image, "start": self.panel_zoom, "target": self.panel_zoom, "at": (wx, wy),
                                   "content": (self.tile_scroll.canvasx(wx), self.tile_scroll.canvasy(wy)), "label": None, "job": None, "photo": None}
            b["spot"] = self.spot_in_tiles(*b["content"])
            self.remember_pane(b, image)
        b["target"] = min(max(b["target"] * factor, 1.0), self.max_panel_zoom())
        self.show_zoom_snapshot(b)
        if b["job"] is not None:
            self.root.after_cancel(b["job"])
        b["job"] = self.root.after(self.ZOOM_SETTLE_MS, self.settle_zoom)

    def pane_look(self):
        """What the pictures of the pane were made with: they are dropped when it changes (a recolor, a size, ...)."""
        return (tuple(id(t["pc"].layers) for t in self.tiles), self.var_dot.get(), self.var_marker.get(), self.var_translucent.get(), self.var_ccf.get(),
                self.var_roi_only.get(), self.var_tile.get(), self.tile_scroll.winfo_width())        # (not its height: the scroll bar that comes with the zoom changes that)

    def remember_pane(self, b, image):
        """Keep the picture of the pane at the zoom a wheel run starts from: when the zoom later goes out again it fills in what the scaled picture does not
        reach. The widest (lowest zoom) one is kept, and the last three."""
        look = self.pane_look()
        if look != getattr(self, "pane_cache_look", None):
            self.pane_cache_look, self.pane_cache = look, []
        b["left"], b["top"] = self.tile_scroll.canvasx(0), self.tile_scroll.canvasy(0)
        entry = {"z": self.panel_zoom, "left": b["left"], "top": b["top"], "image": image}
        self.pane_cache = [e for e in self.pane_cache if abs(e["z"] / self.panel_zoom - 1.0) > 0.02] + [entry]
        if len(self.pane_cache) > 4:
            widest = min(self.pane_cache, key=lambda e: e["z"])
            self.pane_cache = [widest] + [e for e in self.pane_cache if e is not widest][-3:]

    def snapshot_pane(self):
        """A picture of the section pane as it is on screen now (PIL Image), or None if it cannot be taken."""
        try:
            from PIL import ImageGrab
            self.root.update_idletasks()
            c = self.tile_scroll
            x, y, w, h = c.winfo_rootx(), c.winfo_rooty(), c.winfo_width(), c.winfo_height()
            image = ImageGrab.grab(bbox=(x, y, x + w, y + h)).convert("RGB")
            return image if image.size == (w, h) else image.resize((w, h))
        except Exception:
            return None

    def layout_at(self, z):
        """The layout the panels have (or would have) at panel zoom z, as layout_tiles makes it: the number across and down, the slot of one panel (width, height),
        the caption's height, and the size of its map."""
        cols = TILE_COLUMNS[self.var_tile.get()]
        w = max(self.tile_scroll.winfo_width(), 100)
        base = max((w - 4) // cols - 8 - 2, MIN_TILE_PX)
        ch = max(t["caption"].winfo_reqheight() for t in self.tiles)
        px = int(base * z)
        ph = int(px / self.tile_aspect)
        border = 2 if z <= 1.0 else 0                                    # (a map of the unzoomed panel has its widget's border)
        mw, mh = px + border, ph + border
        return {"cols": cols, "rows": -(-len(self.tiles) // cols), "fw": mw, "fh": ch + mh, "ch": ch, "mw": mw, "mh": mh}

    def view_at(self, z, spot, at, lay):
        """Where the pane is scrolled to (left, top) at zoom z when the spot under the pointer is to stay under it."""
        k, dx, dy, old_w, old_h = spot
        x = 4 + (k % lay["cols"]) * (lay["fw"] + 8) + grown(dx, old_w, lay["mw"])
        y = 3 + (k // lay["cols"]) * (lay["fh"] + 6) + lay["ch"] + grown(dy, old_h, lay["mh"])
        total_w = 4 + lay["cols"] * (lay["fw"] + 8)
        total_h = 3 + lay["rows"] * (lay["fh"] + 6)
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
                        borderMode=cv2.BORDER_CONSTANT, borderValue=(244, 244, 244))
        return out, inside

    def earlier_entry(self, b):
        """The earlier picture of the pane to fill the edges with: the sharpest at or below the zoom reached, else the nearest above; None if there is none."""
        target, start = b["target"], b["start"]
        best = None
        for e in getattr(self, "pane_cache", []):
            if abs(e["z"] / start - 1.0) <= 0.02:
                continue                                                 # (the snapshot being scaled)
            rank = (1, e["z"]) if e["z"] <= target * 1.02 else (0, -e["z"])
            if best is None or rank > best[0]:
                best = (rank, e)
        return best[1] if best else None

    def show_zoom_snapshot(self, b):
        """The snapshot of the pane as it will be at the zoom so far, over the pane: each part of it moved and grown as the real layout will have it (the maps
        grow, the captions and gaps do not), so that what is drawn for real at the end lands where this was; and the edges it does not reach from an earlier
        picture of the pane."""
        w, h = b["image"].size
        z = b["target"]
        lay_n = self.layout_at(z)
        left_n, top_n = self.view_at(z, b["spot"], b["at"], lay_n) if b.get("spot") else (b["left"] * z / b["start"], b["top"] * z / b["start"])
        out, inside = self.pane_picture(b["image"], self.pane_maps(lay_n, left_n, top_n, self.layout_at(b["start"]), b["left"], b["top"], w, h), w, h)
        if not inside.all():
            e = self.earlier_entry(b)
            if e is not None:
                under, _ = self.pane_picture(e["image"], self.pane_maps(lay_n, left_n, top_n, self.layout_at(e["z"]), e["left"], e["top"], w, h), w, h)
                under[inside] = out[inside]
                out = under
        b["photo"] = ImageTk.PhotoImage(Image.fromarray(out), master=self.tile_scroll)
        if b["label"] is None:
            b["label"] = tk.Label(self.tile_scroll, bd=0, highlightthickness=0, cursor="fleur")
            b["label"].bind("<MouseWheel>", lambda e: self.zoom_step(*b["at"], 1.25 if e.delta > 0 else 0.8))
            b["label"].bind("<Button-4>", lambda e: self.zoom_step(*b["at"], 1.25))
            b["label"].bind("<Button-5>", lambda e: self.zoom_step(*b["at"], 0.8))
            b["label"].place(x=0, y=0, width=w, height=h)
        b["label"].configure(image=b["photo"])
        b["label"].lift()
        b["label"].update_idletasks()

    def settle_zoom(self):
        """The wheel has been still: the panels are laid out and drawn for real at the zoom reached, and the snapshot taken away."""
        b = self.zoom_burst
        if b is None:
            return
        b["job"] = None
        wx, wy = b["at"]
        self.set_panel_zoom(b["target"], (b["content"][0], b["content"][1], wx, wy))
        self.root.update_idletasks()
        self.refresh_visible()                                      # (what is in view is drawn before the snapshot goes, so that there is no blank between)
        self.root.update_idletasks()
        self.cancel_zoom_burst()

    def cancel_zoom_burst(self):
        b, self.zoom_burst = self.zoom_burst, None
        if b is not None:
            if b["job"] is not None:
                self.root.after_cancel(b["job"])
            if b["label"] is not None:
                b["label"].destroy()

    def panel_zoom_reset(self):
        self.cancel_zoom_burst()
        self.set_panel_zoom(1.0)

    def bind_pane_gestures(self, widget):
        """The wheel, a drag and a double-click over `widget` (the pane itself, a panel's frame or its title: not a section map, which has its
        own) act on the whole panel while Zoom whole panel is on. The wheel over the pane itself scrolls it as it always did."""
        widget.bind("<MouseWheel>", self.on_pane_wheel, add="+")
        for press, motion in (("<ButtonPress-1>", "<B1-Motion>"), ("<ButtonPress-2>", "<B2-Motion>"), ("<ButtonPress-3>", "<B3-Motion>")):
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
        self.cancel_zoom_burst()
        self.tile_scroll.scan_mark(e.x_root - self.tile_scroll.winfo_rootx(), e.y_root - self.tile_scroll.winfo_rooty())

    def panel_pan_move(self, e):
        """A drag of the zoomed panel moves it at once: the maps whose pictures no longer cover the screen are not all drawn again there and
        then (which held the drag up as it went past the margin they are drawn with), but one at a time between the mouse's moves (pan_fill),
        the new parts blank meanwhile; and everything once the drag has been still for PANEL_PAN_SETTLE_MS (pan_settle)."""
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
        self.level = self.var_level.get()
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
        self.select_bar.pack_forget()
        self.legend.pack_forget()
        self.group_list.pack_forget()
        if self.mode == "all":
            self.color_all()
        elif self.mode == "single":
            self.color_single()
        else:
            self.color_gene()
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
        if not self.level_has_ids.get(level, True):
            self.var_status.set(f"{plural_name(self.level_names[level]).capitalize()} have names, not numbers: use All groups to see them (click a group in the list to gray it).")
            return
        ids, bad = self.parse_ids()
        codes, uniques = self.category_codes(level)
        uid = core.extract_leading_numeric_id(pd.Series(uniques)).to_numpy()
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
            parts.append(f"Enter one or more {self.level_names[level].lower()} numbers.")
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
        rows = [("#%02x%02x%02x" % tuple(int(v) for v in colors[t]), str(t), f"({counts[t]})",
                 re.sub(r"^\s*\d+\s*", "", (names or {}).get(t, "")), t in self.off_groups, int(self.shape_codes(np.array([i]))[0]))
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
            self.keep_gene_settings(self.gene_cache[0] if self.gene_cache is not None else None, key)
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

    def keep_gene_settings(self, old_key, new_key):
        """A different set of genes is going to be shown: a gene that was shown before keeps what was set for it (the black and white levels of its
        slider, and whether it was hidden); a new gene starts at the full range, shown. Nothing is kept when the dataset or the scale changed (the
        levels are fractions of a range that is not the same any more)."""
        names = list(new_key[0])
        levels, hidden = {}, set()
        if old_key is not None and old_key[1:] == new_key[1:]:
            old_names = list(old_key[0])
            old_levels = {old_names[0]: self.gene_levels} if len(old_names) == 1 else {g: self.multi_levels.get(i, (0.0, 1.0)) for i, g in enumerate(old_names)}
            hidden_names = {old_names[i] for i in self.hidden if i < len(old_names)} if len(old_names) > 1 else set()
            levels = {i: old_levels[g] for i, g in enumerate(names) if g in old_levels}
            hidden = {i for i, g in enumerate(names) if g in hidden_names}
        self.gene_levels = levels.get(0, (0.0, 1.0)) if len(names) == 1 else (0.0, 1.0)
        self.multi_levels = levels if len(names) > 1 else {}
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
        else:
            expr = adata.layers["counts"][:, col]
            expr = expr.toarray() if hasattr(expr, "toarray") else expr
            vals = np.asarray(expr).ravel().astype(float)
            if scale == "log":
                vals = core.log2p1(vals)
        return vals

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
        fin = np.isfinite(vals)
        vmin, vmax = (float(np.nanmin(vals)), float(np.nanmax(vals))) if fin.any() else (0.0, 1.0)
        # In linear scale the black and white levels can be moved (the slider under the color bar): they cut the range the colors cover.
        linear = scale == "linear"
        self.clear_multi_sliders()
        low, high = self.gene_levels if linear else (0.0, 1.0)
        vmin, vmax = vmin + low * (vmax - vmin), vmin + high * (vmax - vmin)
        self.shown_range = (vmin, vmax)                 # (the colors' range as shown, with the black and white levels: a saved picture's color bar has it too)
        self.show_levels_slider(linear)
        if linear:
            self.cbar_slider.set(low, high)
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
        self.draw_colorbar(vmin, vmax, scale, dataset, readout=(percent_below(vals, vmin), percent_below(vals, vmax)) if linear else None)
        self.var_status.set(f"Showing {gene}{' [imputed]' if dataset == 'imputed' else ''} (every cell of the sections, colored on the selection's scale)."
                            if bg_vals is not None else f"Showing {gene}." if self.own is not None else "(cells outside this run's selection are gray).")

    # ---- several genes: red / green / blue circles for the first three (cyan / magenta / yellow in the light theme), the same three colors as small dots on top for the next three

    def channel_names(self):
        return CHANNEL_NAMES_LIGHT if self.var_theme.get() == "light" else CHANNEL_NAMES

    def channel_hex(self):
        return CHANNEL_HEX_LIGHT if self.var_theme.get() == "light" else CHANNEL_HEX

    def show_multi(self, genes, dataset, scale, per_gene):
        self.show_levels_slider(False)
        light = self.var_theme.get() == "light"
        mix = (lambda chans: ink_rgb(chans)) if light else (lambda chans: (core.multi_gene_rgb(chans, (0.0, 0.0, 0.0)) * 255).astype(np.uint8))
        n = len(genes)
        norm = []
        for i, (vals, _bg) in enumerate(per_gene):
            fin = np.isfinite(vals)
            lo, hi = (float(np.nanmin(vals)), float(np.nanmax(vals))) if fin.any() else (0.0, 1.0)
            if scale == "linear":                  # (the black and white levels set with the legend's sliders cut the range)
                low, high = self.multi_levels.get(i, (0.0, 1.0))
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
        readouts = [(percent_below(per_gene[i][0], norm[i][1]), percent_below(per_gene[i][0], norm[i][2])) for i in range(n)] if scale == "linear" else None
        self.draw_multi_legend(genes, [(v[1], v[2]) for v in norm], scale, readouts)
        parts = " and ".join(f"'{g}' ({self.channel_names()[i]})" for i, g in enumerate(genes))
        self.var_status.set(f"Showing {parts}{' [imputed]' if dataset == 'imputed' else ''}."
                            + (f" Hidden: {', '.join(genes[i] for i in sorted(hidden))}." if hidden else ""))

    def draw_multi_legend(self, genes, ranges, scale, readouts=None):
        c = self.cbar
        self.cbar_photo = None
        linear = scale == "linear"
        h8, h10 = self.line_height(8), self.line_height(10, bold=True)
        name_y = h10 // 2 + 3                                       # (the name and its circle; then the range; in linear scale what is below the two levels, and their slider)
        range_y = h10 + 4
        cells_y = range_y + h8 + 2
        slider_y = cells_y + 2 * h8 + 4
        row = slider_y + 26 + 10 if linear else range_y + h8 + 12
        if len(self.multi_sliders) != (len(genes) if linear else 0):
            c.delete("all")
            self.clear_multi_sliders()                              # (the sliders are kept while they are dragged: only made again for other genes)
            if linear:
                self.multi_sliders = [RangeSlider(c, lambda lo, hi, i=i: self.on_multi_levels(i, lo, hi)) for i in range(len(genes))]
                self.multi_windows = [c.create_window(0, 0, anchor="nw", window=sl, height=26) for sl in self.multi_sliders]
        else:
            c.delete("legend")                                      # (only the texts and circles: the sliders stay where they are, so that a drag in progress goes on)
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
            if linear:
                if readouts is not None and not off:
                    for (a, e), x, anchor in ((readouts[i][0], 10, "nw"), (readouts[i][1], max(c.winfo_width(), 200) - 10, "ne")):
                        c.create_text(x, y + cells_y, anchor=anchor, text=f"{a:.1f}% of cells below", font=("Segoe UI", 8), fill="gray30", tags="legend")
                        c.create_text(x, y + cells_y + h8, anchor=anchor, text=f"{e:.1f}% of expressing", font=("Segoe UI", 8), fill="gray30", tags="legend")
                self.multi_sliders[i].set(*self.multi_levels.get(i, (0.0, 1.0)))
                c.coords(self.multi_windows[i], 0, y + slider_y)
                c.itemconfig(self.multi_windows[i], width=max(c.winfo_width(), 200))
        width = max(c.winfo_width(), 200)
        hint = c.create_text(width / 2, 8 + row * len(genes) + 2, anchor="n", width=width - 20, justify=tk.CENTER,          # (wrapped: it is wider than a narrow panel)
                             text="click a gene to hide / show it" + ("; drag the handles for its black and white levels" if linear else ""),
                             fill="#888888", font=("Segoe UI", 8, "italic"), tags="legend")
        c.config(height=c.bbox(hint)[3] + 6)                        # (as tall as the hint, however many lines it took)
        self.multi_rows = (8, row, len(genes))

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
        else:
            self.cbar_slider.pack_forget()

    def line_height(self, size, bold=False):
        """The height in pixels of a line of Segoe UI text of this size on this screen (the scaling of Windows makes it more than the points say)."""
        fonts = self.__dict__.setdefault("line_fonts", {})
        key = (size, bold)
        if key not in fonts:
            fonts[key] = tkfont.Font(root=self.root, family="Segoe UI", size=size, weight="bold" if bold else "normal")
        return fonts[key].metrics("linespace")

    def draw_colorbar(self, vmin, vmax, scale, dataset, readout=None):
        """`readout` (linear scale): the percentages below the black and the white point, ((all, expressing), (all, expressing))."""
        c = self.cbar
        c.delete("all")
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
            for (a, e), x, anchor in ((readout[0], 10, "nw"), (readout[1], w - 10, "ne")):
                c.create_text(x, y, anchor=anchor, text=f"{a:.1f}% of cells below", font=("Segoe UI", 8), fill="gray30")
                c.create_text(x, y + h8, anchor=anchor, text=f"{e:.1f}% of expressing", font=("Segoe UI", 8), fill="gray30")
            y += 2 * h8 + 4
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
        for t in self.tiles:                  # (the table of the sections' z values is read in the background: it may arrive later)
            if t["title"].get() == t["token"]:
                z = core.section_z_text(t["label"], self.abc_cache, wait=False, short=True)
                if z:
                    t["title"].set(f"{t['token']}   {z}")
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
        z = core.section_z_text(sec, self.abc_cache, wait=False) if sec is not None else ""
        if z:
            lines.append(z)
        if i is not None and self.values is not None:
            for gene, vals in self.values:
                if self.own is not None and self.own.is_number(gene):
                    continue
                lines.append(f"{gene} = {vals[i]:.3g}" if np.isfinite(vals[i]) else f"{gene} = n/a")
        elif i is None:
            lines.append("(not in the UMAP)")
        return "\n".join(lines)

    def region_of(self, info):
        """(acronym, name, structure, division) of the CCF region Allen gives the hovered cell (the registration's own answer, whatever the outlines
        drawn show), or None if there is no table of regions or the cell has none."""
        if self.region_names is None:
            return None
        i = info["umap_idx"]
        if i is not None:
            value = self.region_of_run[i] if self.region_of_run is not None else np.nan
        else:
            arr = self.bg.get("region") if self.bg else None
            value = arr[info["bg_idx"]] if arr is not None else np.nan
        if value is None or not np.isfinite(value):
            return None
        return self.region_names.get(int(value))

    def cell_rows(self, info):
        """The hovered cell as table rows: (what, number, name). A group's label '268 CS-PRNr-DR ...' is split into its number and its name."""
        rows = []
        sec = info["sec"]
        if sec is not None:
            rows.append(("Section", core.sanitize_section_token(sec), ""))
        region = self.region_of(info)
        if region is not None:
            acronym, name, structure, division = region
            if acronym == "unassigned":
                rows.append(("Region", "", "unassigned (outside the CCF atlas)"))
            else:
                rows.append(("Region", acronym, name))
                if structure and structure not in (name, "unassigned"):
                    rows.append(("", "", f"in {structure}"))
                if division and division not in (structure, name, "unassigned"):
                    rows.append(("Division", "", division))
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
            rows.append((self.level_names[k], m.group(1) if m else "", m.group(2) if m else str(v)))
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
        """The cell's position: the atlas's own z of its section, and approximate Paxinos AP, ML and DV (core.paxinos_coords). From the cell's
        own CCF coordinates once they are read (see start_ccf_lookup); until then AP only, from its section's place in the CCF."""
        sec = info["sec"]
        if sec is None:
            return []
        core.start_section_coordinate_load(self.abc_cache)
        z, section_ap = core.section_coords(sec)
        rows = [("Atlas z", f"{z:.2f}", "mm (the section's)")] if z is not None else []
        ccf = None
        if info["umap_idx"] is not None and self.ccf_run is not None:
            ccf = self.ccf_run[info["umap_idx"]]
        elif info["bg_idx"] is not None and self.ccf_bg is not None:
            ccf = self.ccf_bg[info["bg_idx"]]
        if ccf is not None and np.isfinite(ccf).all():
            ap, ml, dv = core.paxinos_coords(float(ccf[0]), float(ccf[1]), float(ccf[2]))
            rows += [("Paxinos AP", f"{ap:+.2f}", "mm from bregma (approx.)"), ("Paxinos ML", f"{ml:.2f}", "mm from the midline"),
                     ("Paxinos DV", f"{dv:+.2f}", "mm from bregma")]
        elif section_ap is not None:
            ap, _ml, _dv = core.paxinos_coords(section_ap)
            rows.append(("Paxinos AP", f"{ap:+.2f}", "mm from bregma (approx., the section's)"))
            if self.ccf_state == "reading":
                rows.append(("Paxinos ML, DV", "", "reading the cells' CCF coordinates..."))
        return rows

    def ccf_action(self):
        """The cell tip's extra button: fetch the CCF coordinates that its Paxinos ML and DV need, while they are not on this computer."""
        if self.ccf_state == "missing" and "ensure_ccf_downloaded" in self.host:
            return "Get ML/DV...", self.download_ccf_coords
        return None

    def start_ccf_lookup(self):
        """Read the CCF coordinates of the run's cells and of the section maps' other cells (for their Paxinos ML and DV), in a thread: slow the
        first time (the table of every cell of the atlas: a minute or so), seconds after. Only if they are on this computer: see ccf_action."""
        self.ccf_run = self.ccf_bg = None
        self.ccf_state = "none"
        h = self.host
        if self.own is not None or self.abc_cache is None or "load_ccf_for_cells" not in h:
            return
        if core.ccf_files_missing(self.abc_cache):
            self.ccf_state = "missing"
            return
        self.ccf_state = "reading"
        run_ids = list(self.adata.obs_names)
        bg_ids = list(self.bg["ids"]) if self.bg is not None and self.bg.get("ids") is not None else []

        def work():
            try:
                table, _shell = h["load_ccf_for_cells"](self.abc_cache, run_ids + bg_ids)
                xyz = table[["x_ccf", "y_ccf", "z_ccf"]].to_numpy(dtype=np.float32)       # (AP, DV, ML)
                self.ccf_run, self.ccf_bg = xyz[:len(run_ids)], (xyz[len(run_ids):] if bg_ids else None)
                self.ccf_state = "ready"
            except Exception as e:
                self.ccf_state = "failed"
                print(f"Could not read the cells' CCF coordinates (for their Paxinos ML and DV): {e}", flush=True)

        threading.Thread(target=work, daemon=True).start()

    def download_ccf_coords(self):
        # (said when the reconstructed coordinates are in use: a table of much the same size, downloaded earlier, that this is easily taken for)
        note = ("This is a separate table from the reconstructed coordinates already downloaded: those place the cells on the section "
                "maps, but are not CCF positions.") if core.reconstructed_in_use() else ""
        if self.host["ensure_ccf_downloaded"](self.abc_cache, parent=self.root, what="CCF coordinates of every cell (for their Paxinos ML and DV)",
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
                    if len(ids) > 1 and m and int(m.group(1)) in ids:
                        canvas, k = self.legend, ids.index(int(m.group(1)))
        self.bold_legend_row(canvas, k)

    def show_info(self, info):
        """Ring the cell where the pointer is, and mirror it in the other plot (its description is in its tip: see hover)."""
        self.bold_group_of(info)
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

    def hover(self, info):
        """The pointer is on a cell (info) or on none (None). The ring and description are immediate; once it rests, the plots dim and, a
        little later, the cell's group is marked with red '+' signs (only those two are turned off by the Hover highlight box). After a longer
        rest the description comes up in a tip as well (CellTip), whether or not the box is ticked."""
        if info is None:
            self.cancel_hover_jobs()
            self.hover_key = None
            self.cell_tip.pointer_off()
            self.show_info(None)
            self.clear_highlights()
            return
        self.show_info(info)
        self.cell_tip.pointer_on(info["key"], lambda: self.cell_rows(info), self.ccf_action)
        if not self.var_hover.get():
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

    def mark_group(self, label, skip=None):
        """Red '+' signs on every cell of the group `label` (at the Group level), in the UMAP and in every section. `skip`: a UMAP cell that
        keeps its own ring instead."""
        level = self.level
        run_mask = self.info[level] == label
        if not run_mask.any():
            return False                             # (a group none of the run's cells is in: nothing in the UMAP to relate to)
        if skip is not None:
            run_mask[skip] = False                   # (the hovered cell keeps its own ring)
        self.umap.set_pluses((self.coords[run_mask, 0], self.coords[run_mask, 1]))
        number = re.match(r"^\s*(\d+)", str(label))               # (the group's number is what the UMAP shows at its center: bigger while it is marked)
        self.umap.boost_label(number.group(1) if number else None)
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
        if not self.var_hover.get():
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
        if not self.var_hover.get():
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
    if "X_umap" not in adata.obsm:
        print("No UMAP embedding available; skipping the interactive UMAP viewer.")
        return
    win = UmapWindow(adata, abc_cache, imputed_state, adata_backed, rois, run_folder, ensure_imputed, title or "Interactive UMAP viewer", parent,
                     host=host, params=params, own=own, panel_zoom=panel_zoom)
    win.run()
