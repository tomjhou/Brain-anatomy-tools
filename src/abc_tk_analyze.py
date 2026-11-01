"""The UMAP viewer's Analyze window, with a tab for each analysis. It starts from what the viewer has (the group numbers, the genes, the level,
the gene source), and then goes its own way.

Expression: plots the second gene against the first, for the cells of some groups, per cell or per group.
Selectivity: one dot for each gene, at its mean in one group (x) and the ratio of that to the highest mean among the largest other groups (y)."""
import functools
from pathlib import Path
import hashlib
import re
import textwrap
import threading
import time
import tkinter as tk
from tkinter import ttk

import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import abc_core as core
from abc_tk_umap import rank_color, plural_name
from abc_tk_umap_export import hexcolor

MAX_LABELED_GROUPS = 2000        # a per-group plot names its points up to this many groups (a supertype plot has about 1,200; the 5,000+ clusters are left unnamed)
LARGEST_GROUP_AREA = 1400.0      # a per-group plot: the marker area (points squared) of the group with the most cells; the others by their cell count
SMALLEST_AREA = 3.0              # (so that a tiny group is still a dot)
LABEL_DARKEN = 0.6               # the group numbers beside the dots: the circle's color times this
LABEL_FONT_SIZE = 8              # the group numbers beside the dots: the same for every group
DENSITY_BINS = 160               # the grid the density shading counts the cells on
HOVER_PX = 9                     # a per-cell plot: how close (pixels) the pointer must be to a cell
DEFAULT_OTHER_GROUPS = 5         # the slider starts on this many of the other groups (while it has not been moved)
MIN_PLOTTED_GROUP_CELLS = 10      # a Per group plot can leave out the groups with fewer cells than this
MIN_GROUP_CELLS = 20             # Selectivity ignores groups with fewer cells (their means are too noisy to be the highest of the others)
SELECTIVITY_FLOOR = 0.01         # added to a gene's mean in the group and in the other groups before the ratio, so that "none anywhere else" is finite
IMPUTED_BLOCK_GENES = 20         # the imputed file is read this many genes at a time
MARK_COLOR = "#d62728"           # the genes typed in the Selectivity tab
LABEL_RATIO = 5.0                # a Selectivity plot names every gene whose ratio (y) is above this
EXPRESSION_DOT_AREA = {False: 4.0, True: 1.5}      # the dots of a per-cell plot (marker area, points squared), for panel genes and for imputed ones
SELECTIVITY_DOT_AREA = {False: 16.0, True: 5.0}    # ... and of a Selectivity plot (500 panel genes, but 8,000 imputed ones)
MAX_LEGEND_GROUPS = 15           # a per-cell plot has a legend up to this many groups
POINT_ALPHA = 0.5


def split_entries(text):
    """The entries of a box: separated by commas. A piece made only of numbers may hold several, separated by spaces ('12 15')."""
    out = []
    for piece in re.split(r"[,;]", text or ""):
        piece = piece.strip()
        if not piece:
            continue
        out.extend(piece.split() if re.fullmatch(r"[\d\s]+", piece) else [piece])
    return out


def match_groups(entries, uniques):
    """(the indexes into `uniques` of the groups the entries name, in the order given; the entries that name none). An entry is a group's number
    ('12') or its name ('Glut', 'IT-ET Glut'), whatever the case."""
    labels = [str(u) for u in uniques]
    names = [re.sub(r"^\s*\d+\s*", "", s).casefold() for s in labels]
    ids = core.extract_leading_numeric_id(pd.Series(labels)).to_numpy()
    chosen, bad = [], []
    for e in entries:
        if e.isdigit():
            hit = np.flatnonzero(ids == float(e))
        else:
            hit = np.array([i for i, (s, n) in enumerate(zip(labels, names)) if e.casefold() in (s.casefold(), n)], dtype=int)
        if len(hit) == 0:
            bad.append(e)
        chosen.extend(int(i) for i in hit if int(i) not in chosen)
    return chosen, bad


def plural(name):
    return plural_name(name)


def darker(color, factor=LABEL_DARKEN):
    """The circle's color, darkened so that the number beside it can be read (a pale circle would give a pale number)."""
    from matplotlib.colors import to_hex, to_rgb
    return to_hex(tuple(v * factor for v in to_rgb(color)))


def jittered(values, imputed, scale, rng):
    """Counts are whole numbers, so cells with the same count are on top of each other in a plot. Each is spread over the stretch of the axis its
    count stands for: a count c over c - 0.5 to c + 0.5 (from 0 for none), on the plotted scale. Imputed values are continuous: left as they are."""
    if imputed or not len(values):
        return values
    c = np.rint(np.exp2(values) - 1.0 if scale == "log" else values)
    low, high = np.maximum(c - 0.5, 0.0), c + 0.5
    if scale == "log":
        low, high = core.log2p1(low), core.log2p1(high)
    return low + (high - low) * rng.random(len(values))


def clamp_limits(lim, home, scale):
    """`lim` (low, high) held inside `home`: no wider than it, and moved back in if it has gone out. None if it already is. On a log axis the
    comparison is of the logarithms."""
    flip = lim[0] > lim[1]
    a, b = (lim[1], lim[0]) if flip else (lim[0], lim[1])
    ha, hb = (home[1], home[0]) if home[0] > home[1] else (home[0], home[1])
    if scale == "log":
        if min(a, ha) <= 0:
            return None
        a, b, ha, hb = np.log10(a), np.log10(b), np.log10(ha), np.log10(hb)
    tolerance = 1e-9 * max(hb - ha, 1e-12)
    span = min(b - a, hb - ha)
    na = min(max(a, ha), hb - span)
    nb = na + span
    if abs(na - a) <= tolerance and abs(nb - b) <= tolerance:
        return None
    if scale == "log":
        na, nb = 10 ** na, 10 ** nb
    return (nb, na) if flip else (na, nb)


def both_scales(value, scale):
    """A plotted value in both units, as text: 'linear 3, log2(1 + x) 2'. `scale` is the one it is plotted in ('log' values are log2(1 + x))."""
    if scale == "log":
        log, linear = float(value), float(np.exp2(value) - 1.0)
    else:
        linear, log = float(value), float(np.log2(1.0 + value))
    text = f"{round(linear):d}" if abs(linear - round(linear)) < 1e-6 and abs(linear) < 1e6 else f"{linear:.3g}"
    return f"linear {text}, log2(1 + x) {log:.3g}"


def local_density(x, y, bins=DENSITY_BINS):
    """For each point, how many points share its cell of a bins x bins grid over the points."""
    h, xe, ye = np.histogram2d(x, y, bins=bins)
    ix = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, bins - 1)
    iy = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, bins - 1)
    return h[ix, iy]


def units(imputed, scale, own=None):
    if imputed:
        return "imputed, log2" if scale == "log" else "imputed, linear"
    if own is not None and not (own.kind == "counts" and own.normalization == "none"):
        return "log2(1 + normalized)" if scale == "log" else "normalized expression"
    return "log2(1 + counts)" if scale == "log" else "counts"


def other_groups_used(slider, n_others):
    """How many of the other groups Selectivity compares with, from the slider (0 to 100): the middle (50) is the largest half of them, the left end
    only the largest, the right end all of them."""
    half = max(int(round(n_others / 2)), 1)
    if slider <= 50:
        k = 1 + (slider / 50.0) * (half - 1)
    else:
        k = half + ((slider - 50.0) / 50.0) * (n_others - half)
    return int(min(max(int(round(k)), 1), n_others))


def slider_for_groups(k, n_others):
    """The slider position (0 to 100) at which other_groups_used gives `k` (or the most there are) of `n_others` groups."""
    k = min(max(int(k), 1), max(n_others, 1))
    half = max(int(round(n_others / 2)), 1)
    if k <= half:
        position = 0.0 if half <= 1 else (k - 1) / (half - 1) * 50.0
    else:
        position = 50.0 + (k - half) / max(n_others - half, 1) * 50.0
    for nudge in (0.0, 0.01, -0.01, 0.1, -0.1, 0.5, -0.5):                     # (rounding: the nearest position that really gives k)
        if other_groups_used(position + nudge, n_others) == k:
            return min(max(position + nudge, 0.0), 100.0)
    return min(max(position, 0.0), 100.0)


def shows_working(method):
    """For a method of the window that reads or draws a plot (Analyze, Log / Linear, a box ticked): the Working... note is on the graph from the
    start until nothing is being read or drawn any more (check_working)."""
    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        self.show_working()
        try:
            return method(self, *args, **kwargs)
        finally:
            self.top.after(50, self.check_working)
    return wrapper


def other_group_sizes(counts, target, min_cells=MIN_GROUP_CELLS):
    """The cell counts of the groups other than `target` that have at least `min_cells` cells, largest first (the order Selectivity takes them in)."""
    return [int(counts[g]) for g in np.argsort(-counts, kind="stable") if g != target and counts[g] >= min_cells]


def selectivity_values(means, counts, target, k, floor=SELECTIVITY_FLOOR, min_cells=MIN_GROUP_CELLS):
    """(x, y, how many other groups were used) for every gene: x is the gene's mean in the `target` group, y the ratio of that to the highest
    mean among the `k` largest other groups with at least `min_cells` cells (`floor` is added to both, so that a gene the others do not express at
    all has a large ratio, not an infinite one). `means`: genes x groups; `counts`: the cells of each group."""
    others = [int(g) for g in np.argsort(-counts, kind="stable") if g != target and counts[g] >= min_cells][:k]
    x = means[:, target].astype(float)
    top = means[:, others].max(axis=1).astype(float) if others else np.zeros(len(x))
    return x, (x + floor) / (top + floor), len(others)


def group_means_panel(counts_matrix, codes, n_groups):
    """The mean raw count of every gene in every group, as a genes x groups array, from the run's cells x genes counts."""
    import scipy.sparse as sp
    ok = np.flatnonzero(codes >= 0)
    indicator = sp.csr_matrix((np.ones(len(ok)), (ok, codes[ok])), shape=(len(codes), n_groups))
    sums = indicator.T @ counts_matrix                        # groups x genes
    sums = sums.toarray() if hasattr(sums, "toarray") else np.asarray(sums)
    n = np.bincount(codes[ok], minlength=n_groups).astype(float)
    return (sums / np.maximum(n, 1.0)[:, None]).T.astype(np.float32)


def group_means_imputed(path, positions, level_codes, progress, cancelled, block_genes=IMPUTED_BLOCK_GENES):
    """The mean linear expression (2**v - 1 of the stored log2 values) of every gene of the imputed file in every group, at every group type
    at once: {level: genes x groups}, or None if cancelled. `positions`: each run cell's row in the file (-1: not in it); `level_codes`:
    {level: (the group of each run cell, how many groups)}.
    The file is stored in blocks of 50,000 cells by 5 genes, and reading is what takes the time, so only the blocks that hold run cells are read
    (the run's cells are a few sections of the file's 87 blocks of cells), a block of genes at a time, once for all the group types."""
    import h5py
    import scipy.sparse as sp
    valid = np.flatnonzero(positions >= 0)
    rows_in_file = positions[valid]
    offsets, total = {}, 0
    for level, (_codes, n_groups) in level_codes.items():
        offsets[level] = total
        total += n_groups
    r_idx, c_idx = [], []
    for level, (codes, _n) in level_codes.items():           # one indicator for all the types, the groups of each under those of the one before
        c = codes[valid]
        ok = np.flatnonzero(c >= 0)
        r_idx.append(c[ok] + offsets[level])
        c_idx.append(ok)
    indicator = sp.csc_matrix((np.ones(sum(len(r) for r in r_idx)), (np.concatenate(r_idx), np.concatenate(c_idx))), shape=(total, len(valid)))
    n_cells = np.asarray(indicator.sum(axis=1)).ravel()
    with h5py.File(path, "r") as f:
        x = f["X"]
        size, n_genes = x.chunks[0], x.shape[1]
        block_of = rows_in_file // size
        blocks = []
        for b in np.unique(block_of):
            pick = np.flatnonzero(block_of == b)
            blocks.append((int(b) * size, min((int(b) + 1) * size, x.shape[0]), rows_in_file[pick] - int(b) * size, indicator[:, pick]))
        sums = np.zeros((total, n_genes), np.float32)
        for g0 in range(0, n_genes, block_genes):
            if cancelled():
                return None
            g1 = min(g0 + block_genes, n_genes)
            acc = np.zeros((total, g1 - g0), np.float32)
            for r0, r1, local, ind in blocks:
                linear = np.exp2(x[r0:r1, g0:g1][local].astype(np.float32)) - 1.0
                acc += ind @ linear
            sums[:, g0:g1] = acc
            progress(g1 / n_genes)
    means = (sums / np.maximum(n_cells, 1.0)[:, None]).T.astype(np.float32)           # genes x all the groups
    return {level: np.ascontiguousarray(means[:, offsets[level]:offsets[level] + n_groups]) for level, (_c, n_groups) in level_codes.items()}


def save_means(path, data):
    """Keep the means of every gene in every group (at every group type) for the next time this run is analyzed."""
    arrays = {"names": np.asarray(data["names"], dtype=str)}
    for level, d in data["levels"].items():
        arrays[f"means__{level}"] = d["means"]
        arrays[f"uniques__{level}"] = np.asarray(d["uniques"], dtype=str)
        arrays[f"counts__{level}"] = np.asarray(d["counts"])
    try:
        np.savez(path, **arrays)
    except Exception as e:
        print(f"Could not save the group means ({e}).", flush=True)


def load_saved_means(path, uniques_of):
    """The means saved by save_means, if there are some and they are for these groups, else None."""
    if not path.exists():
        return None
    try:
        with np.load(path, allow_pickle=False) as z:
            levels = {}
            for level, uniques in uniques_of.items():
                if not np.array_equal(z[f"uniques__{level}"], uniques):
                    return None
                levels[level] = {"means": z[f"means__{level}"], "uniques": uniques, "counts": z[f"counts__{level}"]}
            return {"names": z["names"], "levels": levels}
    except Exception as e:
        print(f"The saved group means could not be read ({e}); they are read again.", flush=True)
        return None


class AnalyzeWindow:
    def __init__(self, viewer):
        self.v = viewer
        self.top = top = tk.Toplevel(viewer.root)
        top.title("Analyze")
        height = min(1020, top.winfo_screenheight() - 90)
        width = min(height - 61 + 316, top.winfo_screenwidth() - 40)      # (the graph's pane, right of the 300 px of controls and under the toolbar, is a square)
        top.geometry(f"{width}x{height}")
        top.minsize(900, 600)
        body = ttk.Frame(top, padding=8)
        body.pack(fill=tk.BOTH, expand=True)
        side = ttk.Frame(body, width=290)
        side.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        side.pack_propagate(False)

        # what both tabs use: the group type, the kind of plot, and where the genes come from
        box = ttk.LabelFrame(side, text="Group type", padding=6)
        box.pack(fill=tk.X)
        self.level_keys = [k for k in viewer.level_names if k in viewer.levels]
        self.var_level = tk.StringVar(top, self.v.level_names[viewer.level] if viewer.level in self.level_keys else self.v.level_names[self.level_keys[0]])
        combo = ttk.Combobox(box, textvariable=self.var_level, values=[self.v.level_names[k] for k in self.level_keys], state="readonly")
        combo.pack(fill=tk.X)
        combo.bind("<<ComboboxSelected>>", lambda ev: self.level_changed())

        box = ttk.LabelFrame(side, text="Plot", padding=6)
        box.pack(fill=tk.X, pady=(6, 0))
        self.var_per = tk.StringVar(top, "cell")
        ttk.Radiobutton(box, text="Per cell", value="cell", variable=self.var_per, command=self.per_changed).pack(anchor="w")
        ttk.Radiobutton(box, text="Per group", value="group", variable=self.var_per, command=self.per_changed).pack(anchor="w")
        self.var_exclude_small = tk.BooleanVar(top, True)
        self.chk_exclude = ttk.Checkbutton(box, text=f"Exclude groups with <{MIN_PLOTTED_GROUP_CELLS} cells", variable=self.var_exclude_small, command=self.redraw)
        self.var_jitter = tk.BooleanVar(top, True)
        self.var_density = tk.BooleanVar(top, False)

        # where each tab's genes come from: Expression starts as the viewer's colouring is, Selectivity on the imputed genes (a dataset of the user's own
        # has none)
        can_impute = viewer.own is None
        self.var_source = {"expression": tk.StringVar(top, "imputed" if can_impute and viewer.mode == "imputed" else "panel"),
                           "selectivity": tk.StringVar(top, "imputed" if can_impute else "panel")}
        self.source_shown = {tab: var.get() for tab, var in self.var_source.items()}

        self.var_scale = tk.StringVar(top, viewer.var_scale.get())         # (the viewer's choice to begin with; changing it here leaves the viewer's alone)

        # the two analyses, a tab each (the tab control is as tall as the taller of them, so nothing below it moves)
        outline = tk.Frame(side, highlightthickness=1, highlightbackground="#7f7f7f", highlightcolor="#7f7f7f", bd=0)      # (a line round the tab control, to see where it is)
        outline.pack(fill=tk.X, pady=(8, 0))
        self.tabs = ttk.Notebook(outline)
        self.tabs.pack(fill=tk.X, padx=1, pady=1)
        page = ttk.Frame(self.tabs, padding=6)
        self.tabs.add(page, text="Expression")
        self.source_radios(page, "expression")
        ttk.Label(page, text="Genes, separated by commas: the first is plotted along x, the second up y.", foreground="gray30", wraplength=250).pack(anchor="w")
        self.var_genes = tk.StringVar(top, viewer.var_gene.get())
        self.genes_entry = ttk.Entry(page, textvariable=self.var_genes)
        self.genes_entry.pack(fill=tk.X, pady=(2, 0))
        self.var_genes_hint = tk.StringVar(top, "")
        self.genes_hint = ttk.Label(page, textvariable=self.var_genes_hint, foreground="#b3261e", wraplength=250, justify=tk.LEFT)     # (under the box, where it is seen; the status is at the foot of the window)
        self.var_cells_of = tk.StringVar(top, "")
        self.cells_of_label = ttk.Label(page, textvariable=self.var_cells_of, foreground="gray30", wraplength=250, justify=tk.LEFT)
        self.cells_of_label.pack(anchor="w", pady=(8, 0))
        # (the scale, and the two boxes about the dots of the Expression plot: Selectivity has no use for them)
        scale_row = ttk.Frame(page)
        scale_row.pack(anchor="w", pady=(10, 0))
        ttk.Label(scale_row, text="Scale:", foreground="gray30").pack(side=tk.LEFT, padx=(0, 8))
        for value, text in (("linear", "Linear"), ("log", "Log")):
            ttk.Radiobutton(scale_row, text=text, value=value, variable=self.var_scale, command=self.scale_changed).pack(side=tk.LEFT, padx=(0, 12))
        self.chk_jitter = ttk.Checkbutton(page, text="Jitter counts", variable=self.var_jitter, command=self.redraw)
        self.chk_jitter.pack(anchor="w", pady=(6, 0))
        from abc_tk_umap import HoverTip                  # (the viewer's own tip; imported here, as that module imports this one)
        HoverTip(self.chk_jitter, lambda: "Cells with the same count would sit on top of each other: each is spread over the stretch of the axis its count stands for.")
        self.chk_density = ttk.Checkbutton(page, text="Shade per-cell dots by density", variable=self.var_density, command=self.redraw)
        self.chk_density.pack(anchor="w")

        page = ttk.Frame(self.tabs, padding=6)
        self.tabs.add(page, text="Selectivity")
        self.source_radios(page, "selectivity")
        ttk.Label(page, text="Group number or name (only the first is used):", foreground="gray30", wraplength=250).pack(anchor="w")
        self.var_groups = tk.StringVar(top, viewer.var_ids.get())
        self.groups_entry = ttk.Entry(page, textvariable=self.var_groups)
        self.groups_entry.pack(fill=tk.X, pady=(2, 0))
        ttk.Label(page, text="Genes to mark in red, separated by commas:", foreground="gray30", wraplength=250, justify=tk.LEFT).pack(anchor="w", pady=(8, 0))
        self.var_sel_genes = tk.StringVar(top, viewer.var_gene.get())
        self.sel_genes_entry = ttk.Entry(page, textvariable=self.var_sel_genes)
        self.sel_genes_entry.pack(fill=tk.X, pady=(2, 0))
        ttk.Label(page, text="Compared with the largest other groups:", foreground="gray30", wraplength=250, justify=tk.LEFT).pack(anchor="w", pady=(8, 0))
        self.var_slider = tk.DoubleVar(top, 50.0)
        self.slider_touched = False             # (until it is moved, the slider keeps to DEFAULT_OTHER_GROUPS whatever the group type and group)
        slider = ttk.Scale(page, from_=0, to=100, variable=self.var_slider, orient=tk.HORIZONTAL, command=lambda v: self.slider_moved())
        slider.pack(fill=tk.X)
        self.var_others = tk.StringVar(top, "")
        ttk.Label(page, textvariable=self.var_others, foreground="gray30", wraplength=250, justify=tk.LEFT).pack(anchor="w")
        for var, entry in ((self.var_groups, self.groups_entry), (self.var_genes, self.genes_entry), (self.var_sel_genes, self.sel_genes_entry)):
            var.trace_add("write", lambda *a, e=entry: (e.configure(foreground=""), self.update_labels()))
            entry.bind("<Return>", lambda ev: self.analyze())
        self.sel = None                         # the Selectivity data on the screen (the means of every gene in every group), so the slider only redraws
        self.means_cache = {}
        self.job = None                         # a read of many genes in the background (see start_means)
        self.slider_job = None

        self.btn_analyze = ttk.Button(side, text="Analyze", command=self.analyze)
        self.btn_analyze.pack(fill=tk.X, pady=(10, 0))
        ttk.Button(side, text="Save graph", command=self.save_graph).pack(fill=tk.X, pady=(6, 0))
        ttk.Button(side, text="Export CSV", command=self.export_csv).pack(fill=tk.X, pady=(6, 0))
        self.var_status = tk.StringVar(top, "")
        ttk.Label(side, textvariable=self.var_status, foreground="gray30", wraplength=270, justify=tk.LEFT).pack(anchor="w", pady=(10, 0))

        right = ttk.Frame(body)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.fig = Figure(figsize=(9.0, 8.4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.text(0.5, 0.5, "Choose groups and genes, then click Analyze.", ha="center", va="center", color="gray", transform=self.ax.transAxes)
        self.ax.set_axis_off()
        self.hover = None                       # what the pointer can be over (see plot_expression, plot_selectivity)
        self.background = None
        self.last_motion = 0.0
        self.hover_group = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.table = None                       # the numbers of the graph on the screen (a DataFrame), for Export CSV
        self.last_plot = None                   # the arguments of the plot on the screen (a checkbox draws it again)
        self.cax = None
        self.legend_room = False                # the plot leaves room at the right for a legend
        self.canvas.mpl_connect("resize_event", lambda ev: self.place_axes())
        self.canvas.mpl_connect("draw_event", self.on_draw)
        self.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.canvas.mpl_connect("figure_leave_event", lambda ev: self.hide_hover())
        self.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("button_release_event", self.on_release)
        self.pan_start = None                   # while the plot is being dragged: the pointer and the limits where it began
        self.home_view = None                   # the limits of the plot as first drawn (a double-click goes back to them)
        self.clamping = False
        bar = ttk.Frame(right)
        bar.pack(side=tk.BOTTOM, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, bar)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.working = tk.Label(right, text="Working...", font=("Segoe UI", 20), foreground="#555555", background="white")      # (covers the graph while an analysis is going)
        self.update_labels()
        self.canvas.draw()

    def export_csv(self):
        """Save the numbers of the graph on the screen as a CSV file: per cell, one row for each cell plotted (its group and the two values, and the
        jittered positions it was drawn at); per group, one row for each group (the cells and the two means); Selectivity, one row for each gene."""
        v = self.v
        if self.table is None:
            self.var_status.set("Click Analyze first: there is no graph to export.")
            return
        proposed = v.host["unique_export_paths"](v.target_dir(), "analysis_graph", (".csv",))[".csv"].name
        path = v.ask_save_path("Export the graph's data", proposed, ".csv", "CSV file", parent=self.top)
        if path is None:
            self.var_status.set("Export CSV cancelled.")
            return
        try:
            self.table.to_csv(path, index=False)
            self.var_status.set(f"Saved {len(self.table):,} rows to {Path(path).name} in {Path(path).parent}.")
        except Exception as e:
            self.var_status.set(f"Could not export the data: {e}")
            print(f"Could not export the data: {e}", flush=True)

    def save_graph(self):
        """Save the graph on the screen as a PNG, JPG or SVG (asked as for the viewer's other saves)."""
        v = self.v
        if self.table is None:
            self.var_status.set("Click Analyze first: there is no graph to save.")
            return
        fmt = v.ask_image_format("Save the graph", "", parent=self.top)
        if fmt is None:
            self.var_status.set("Save graph cancelled.")
            return
        paths = v.ask_export_paths(fmt, "analysis_graph", "Save the graph", parent=self.top)
        if paths is None:
            self.var_status.set("Save graph cancelled.")
            return
        path = paths[fmt]
        self.top.config(cursor="watch")
        self.var_status.set("Saving the graph...")
        self.top.update()
        try:
            from abc_tk_umap_export import JPG_ARGS
            self.hide_hover()
            self.fig.savefig(path, bbox_inches="tight", dpi=300 if fmt != "svg" else None, **({"pil_kwargs": JPG_ARGS} if fmt == "jpg" else {}))
            self.var_status.set(f"Saved {path.name} to {path.parent}.")
        except Exception as e:
            self.var_status.set(f"Could not save the graph: {e}")
            print(f"Could not save the graph: {e}", flush=True)
        self.top.config(cursor="")

    # ---- the boxes and the choices above the tabs

    def source_radios(self, page, tab):
        """Panel gene / Imputed gene, at the top of a tab's page (not for a dataset of the user's own: it has no imputed genes)."""
        if self.v.own is not None:
            return
        ttk.Label(page, text="Genes from:", foreground="gray30").pack(anchor="w")          # (above the two, which would not fit beside it)
        row = ttk.Frame(page)
        row.pack(anchor="w", pady=(0, 6))
        for value, text in (("panel", "Panel gene"), ("imputed", "Imputed gene")):
            ttk.Radiobutton(row, text=text, value=value, variable=self.var_source[tab], command=lambda t=tab: self.source_changed(t)).pack(side=tk.LEFT, padx=(0, 10))

    def imputed(self, tab):
        return self.var_source[tab].get() == "imputed"

    def ensure_source(self, tab):
        """True once the tab's genes can be read: for imputed genes the imputed dataset is loaded first, if it is not (as when the viewer's own
        Imputed Gene is chosen). If that is declined or fails the tab goes back to the genes it had."""
        v = self.v
        if self.imputed(tab) and v.imputed_state["adata"] is None:
            ok = bool(v.ensure_imputed(v.imputed_state, v.abc_cache)) if v.ensure_imputed else False
            if not ok:
                self.var_source[tab].set(self.source_shown[tab] if self.source_shown[tab] != "imputed" else "panel")
                return False
        self.source_shown[tab] = self.var_source[tab].get()
        return True

    def source_changed(self, tab):
        """Panel gene / Imputed gene was clicked on a tab."""
        if not self.ensure_source(tab):
            return
        if tab == "selectivity":
            self.sel = None
        self.var_status.set("")

    def per_changed(self):
        """Per cell / Per group: the two boxes are about the dots of single cells, so they are greyed out for a plot of groups."""
        state = "normal" if self.var_per.get() == "cell" else "disabled"
        self.chk_jitter.configure(state=state)
        self.chk_density.configure(state=state)
        if self.var_per.get() == "group":                  # (the box about small groups, just under the radio button, is for a plot of groups)
            self.chk_exclude.pack(anchor="w", padx=(20, 0))
        else:
            self.chk_exclude.pack_forget()

    def level_changed(self):
        self.sel = None                         # (the means are kept for every group type: the next Analyze plots the new one from them)
        self.update_labels()

    def show_genes_hint(self, text):
        """The note under the genes box (nothing hides it)."""
        self.var_genes_hint.set(text)
        if text:
            self.genes_hint.pack(anchor="w", pady=(4, 0), before=self.cells_of_label)
        else:
            self.genes_hint.pack_forget()

    def scale_changed(self):
        """Log / Linear: the plot on the screen is read and drawn again on the new scale."""
        if self.last_plot is not None:
            self.expression()

    def update_labels(self):
        """The small texts that follow the boxes: which groups the Expression tab plots, and how many other groups the slider compares with."""
        n_genes = len(core.parse_gene_names(self.var_genes.get()))
        self.show_genes_hint("Two genes are needed: add a second after a comma." if n_genes == 1 else "")
        entries = split_entries(self.var_groups.get())
        self.var_cells_of.set("Cells of the groups on the Selectivity tab: " + (", ".join(entries) if entries else "every group") + ".")
        try:
            codes, uniques = self.v.category_codes(self.level())
        except Exception:
            return
        counts = np.bincount(codes[codes >= 0], minlength=len(uniques))
        chosen = match_groups(entries[:1], uniques)[0] if entries else []
        sizes = other_group_sizes(counts, chosen[0] if chosen else None)
        n_others = max(len(sizes), 1)
        if not self.slider_touched:
            self.var_slider.set(slider_for_groups(DEFAULT_OTHER_GROUPS, n_others))
        k = other_groups_used(self.var_slider.get(), n_others)
        smallest = sizes[k - 1] if sizes else 0
        self.var_others.set(f"Largest {k} of {n_others} other groups used, each with at least {smallest:,} cells (groups under {MIN_GROUP_CELLS} cells are ignored).")

    def slider_moved(self):
        self.slider_touched = True
        self.update_labels()
        if self.sel is not None and self.slider_job is None:
            self.slider_job = self.top.after(60, self.slider_redraw)

    def slider_redraw(self):
        self.slider_job = None
        if self.sel is not None and self.tabs.index(self.tabs.select()) == 1:
            self.plot_selectivity()

    def analyze(self):
        """The Analyze button: the analysis of the tab that is showing (while a long read is going, the same button cancels it)."""
        if self.job is not None:
            self.job["cancel"] = True
            self.var_status.set("Cancelling...")
            return
        if self.tabs.index(self.tabs.select()) == 0:
            self.expression()
        else:
            self.selectivity()

    def show_working(self):
        """Cover the graph with the Working... note (painted now: reading and drawing a big plot holds the window up)."""
        self.working.place(in_=self.canvas.get_tk_widget(), relx=0, rely=0, relwidth=1, relheight=1)
        self.top.update()

    def check_working(self):
        """Take the Working... note off the graph once nothing is being read or drawn any more (the reads run in the background and draw when done)."""
        if self.v.busy or self.job is not None:
            self.top.after(150, self.check_working)
        else:
            self.working.place_forget()

    # ---- what the boxes ask for

    def level(self):
        names = {self.v.level_names[k]: k for k in self.level_keys}
        return names[self.var_level.get()]

    def groups(self, level):
        """(codes of each run cell at the level, the groups' indexes in `uniques`, uniques, how they were asked: typed or every group) or None
        (the status says why)."""
        codes, uniques = self.v.category_codes(level)
        entries = split_entries(self.var_groups.get())
        if not entries:
            counts = np.bincount(codes[codes >= 0], minlength=len(uniques))
            return codes, [int(i) for i in np.argsort(-counts, kind="stable") if counts[i] > 0], uniques, False
        chosen, bad = match_groups(entries, uniques)
        if bad:
            self.groups_entry.configure(foreground="red")
            self.var_status.set(f"No {self.v.level_names[level].lower()} named {', '.join(bad)} in this run.")
            return None
        return codes, chosen, uniques, True

    def genes(self):
        """The two genes to plot, as spelled in the dataset, or None (the status says why)."""
        v = self.v
        names = core.parse_gene_names(self.var_genes.get())
        if len(names) < 2:
            self.genes_entry.configure(foreground="red")
            self.show_genes_hint("Two genes are needed: enter them separated by a comma.")
            self.var_status.set("Enter two genes on the Expression tab, separated by a comma.")
            return None
        if not self.ensure_source("expression"):
            self.var_status.set("The imputed genes are not available: Panel gene is used.")
        imputed = self.imputed("expression")
        if not imputed and not v.has_counts:
            self.var_status.set("No raw counts on this dataset (layers['counts'] missing).")
            return None
        source = v.imputed_state["adata"] if imputed else v.adata
        _cols, resolved, missing = core.resolve_gene_names(source, names)
        if missing:
            self.genes_entry.configure(foreground="red")
            self.var_status.set(f"Not found in the {'imputed' if imputed else 'standard'} dataset: {', '.join(missing)}.")
            return None
        return resolved

    # ---- Expression

    @shows_working
    def expression(self):
        v = self.v
        if v.busy:
            self.var_status.set("The viewer is busy reading something: try again in a moment.")
            return
        level = self.level()
        picked = self.groups(level)
        if picked is None:
            return
        genes = self.genes()
        if genes is None:
            return
        codes, chosen, uniques, typed = picked
        imputed, scale, per = self.imputed("expression"), self.var_scale.get(), self.var_per.get()
        use = genes[:2]
        extra = f" (only the first two of {len(genes)} are plotted)" if len(genes) > 2 else ""

        def work():
            return [v.read_run_gene(g, imputed, scale) for g in use]

        def done(res):
            if isinstance(res, BaseException):
                self.var_status.set(f"Could not read {', '.join(use)}: {res}")
                return
            self.last_plot = (level, codes, chosen, uniques, typed, use, res, imputed, scale)
            self.plot_expression(*self.last_plot, per)
            self.var_status.set(self.var_status.get() + extra)

        self.var_status.set(f"Reading {', '.join(use)}...")
        if not v.job(f"Reading {', '.join(use)}...", work, done):
            self.var_status.set("The viewer is busy reading something: try again in a moment.")

    @shows_working
    def redraw(self):
        """A box was ticked or cleared: the same genes and groups drawn again (nothing is read again)."""
        if self.last_plot is not None:
            self.plot_expression(*self.last_plot, self.var_per.get())

    def plot_expression(self, level, codes, chosen, uniques, typed, genes, values, imputed, scale, per):
        gx, gy = genes
        vx, vy = values
        ok = np.isfinite(vx) & np.isfinite(vy)
        position = {g: i for i, g in enumerate(chosen)}
        color_of = {g: hexcolor(rank_color(i)) for g, i in position.items()}
        in_group = np.isin(codes, chosen) & ok
        fig, ax = self.fig, self.ax
        if self.cax is not None:
            self.cax.remove()
            self.cax = None
        ax.clear()
        ax.set_axis_on()
        self.home_view = None                               # (the new plot's own limits are not held to the old plot's)
        self.hover = None
        self.hover_group = None
        label_of = {g: (re.match(r"^\s*(\d+)", str(uniques[g])).group(1) if re.match(r"^\s*(\d+)", str(uniques[g])) else str(uniques[g])) for g in chosen}
        full_of = {g: str(uniques[g]).strip() for g in chosen}
        unit = units(imputed, scale, self.v.own)
        own = self.v.own
        continuous = imputed or (own is not None and not (own.kind == "counts" and own.normalization == "none"))      # (values that are not whole counts are not jittered)
        if per == "cell":
            cell_group = codes[in_group]
            px, py = vx[in_group], vy[in_group]
            note = ""
            if self.var_jitter.get():
                rng = np.random.default_rng(0)
                px, py = jittered(px, continuous, scale, rng), jittered(py, continuous, scale, rng)
                note = " These values are continuous: not jittered." if continuous else " Counts jittered within their step."
            columns = {"cell": np.asarray(self.v.adata.obs_names)[in_group], self.v.level_names[level]: [full_of[g] for g in cell_group],
                       gx: vx[in_group], gy: vy[in_group]}
            if self.var_jitter.get():
                columns[f"{gx} (plotted, jittered)"], columns[f"{gy} (plotted, jittered)"] = px, py
            self.table = pd.DataFrame(columns)
            shade = self.var_density.get() and len(px) > 0
            if shade:
                from matplotlib.colors import LogNorm
                dens = local_density(px, py)
                order = np.argsort(dens, kind="stable")                      # (the dense dots on top)
                sc = ax.scatter(px[order], py[order], c=dens[order], cmap="viridis", norm=LogNorm(vmin=1, vmax=max(float(dens.max()), 2.0)), s=EXPRESSION_DOT_AREA[imputed],
                                linewidths=0, rasterized=True)
                self.cax = ax.inset_axes([1.03, 0.0, 0.03, 1.0])
                fig.colorbar(sc, cax=self.cax, label="cells in the same small square of the graph")
            else:
                colors = np.array([color_of[g] for g in cell_group]) if len(cell_group) else []
                ax.scatter(px, py, c=colors, s=EXPRESSION_DOT_AREA[imputed], alpha=POINT_ALPHA, linewidths=0, rasterized=True)
            self.hover = {"x": px, "y": py, "group": cell_group, "full": full_of, "color": color_of,
                          "ox": vx[in_group], "oy": vy[in_group], "genes": (gx, gy), "scale": scale}          # (the cell's own values: the dots may be jittered)
            self.pluses = ax.scatter([], [], marker="+", c="black", s=40, linewidths=0.9, animated=True, zorder=5)
            self.tip = ax.annotate("", xy=(0, 0), xytext=(14, 14), textcoords="offset points", fontsize=9, animated=True, zorder=6,
                                   bbox={"boxstyle": "round,pad=0.3", "fc": "#ffffe0", "ec": "#888888"}, annotation_clip=False)
            self.tip.set_visible(False)
            if len(chosen) <= MAX_LEGEND_GROUPS and not shade:
                from matplotlib.lines import Line2D
                handles = [Line2D([], [], marker="o", linestyle="", color=color_of[g], label=f"{full_of[g]} ({int((cell_group == g).sum()):,})")
                           for g in chosen if (cell_group == g).any()]
                if handles:
                    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, frameon=True, title=self.v.level_names[level])
            ax.set_xlabel(f"{gx} ({unit})")
            ax.set_ylabel(f"{gy} ({unit})")
            ax.set_title(f"{int(in_group.sum()):,} cells in {len(set(cell_group.tolist()))} {plural(self.v.level_names[level])}")
            self.var_status.set(f"{int(in_group.sum()):,} cells plotted." + note)
        else:
            xs, ys, names, cols, counts, shown_groups = [], [], [], [], [], []
            self.hover = None
            for g in chosen:
                m = in_group & (codes == g)
                if not m.any() or (self.var_exclude_small.get() and m.sum() < MIN_PLOTTED_GROUP_CELLS):
                    continue
                xs.append(float(np.mean(vx[m])))
                ys.append(float(np.mean(vy[m])))
                names.append(label_of[g])
                shown_groups.append(g)
                cols.append(color_of[g])
                counts.append(int(m.sum()))
            self.table = pd.DataFrame({self.v.level_names[level]: [full_of[g] for g in shown_groups], "cells": counts, f"mean {gx}": xs, f"mean {gy}": ys})
            areas = np.maximum(np.array(counts, dtype=float) / max(counts) * LARGEST_GROUP_AREA, SMALLEST_AREA) if counts else []
            ax.scatter(xs, ys, c=cols, s=areas, alpha=0.85, linewidths=0.4, edgecolors="black")
            if len(xs) <= MAX_LABELED_GROUPS:
                for x, y, n, a, c in zip(xs, ys, names, areas, cols):
                    ax.annotate(n, (x, y), xytext=(np.sqrt(a) / 2 + 2, 0), textcoords="offset points", fontsize=LABEL_FONT_SIZE,
                                ha="left", va="center", color=darker(c))
            ax.set_xlabel(f"mean {gx} ({unit})")
            ax.set_ylabel(f"mean {gy} ({unit})")
            ax.set_title(f"{len(xs)} {plural(self.v.level_names[level])}")
            self.var_status.set(f"{len(xs)} groups plotted (the mean of each, over {sum(counts):,} cells; the circle's area follows the cell count).")
        ax.grid(True, alpha=0.25)
        ax.margins(0.08)                                      # (room for the big dots and their numbers at the edges)
        self.legend_room = per == "cell" and len(chosen) <= MAX_LEGEND_GROUPS and not self.var_density.get()      # (a legend, not the colour bar: its room is measured)
        self.place_axes()
        self.canvas.draw()
        self.remember_view()

    # ---- Selectivity: every gene's mean in one group against how much higher it is than in the other groups

    @shows_working
    def selectivity(self):
        v = self.v
        level = self.level()
        codes, uniques = v.category_codes(level)
        entries = split_entries(self.var_groups.get())
        if not entries:
            self.groups_entry.configure(foreground="red")
            self.var_status.set("Enter a group on this tab: Selectivity compares it with the others.")
            return
        chosen, bad = match_groups(entries[:1], uniques)
        if bad or not chosen:
            self.groups_entry.configure(foreground="red")
            self.var_status.set(f"No {self.v.level_names[level].lower()} named {entries[0]} in this run.")
            return
        note_source = ""
        if not self.ensure_source("selectivity"):
            note_source = " The imputed genes are not available: Panel gene is used."
        imputed = self.imputed("selectivity")
        if not imputed and not v.has_counts:
            self.var_status.set("No raw counts on this dataset (layers['counts'] missing).")
            return
        target = chosen[0]
        counts = np.bincount(codes[codes >= 0], minlength=len(uniques))
        if counts[target] < MIN_GROUP_CELLS:
            self.groups_entry.configure(foreground="red")
            self.var_status.set(f"{str(uniques[target]).strip()} has only {int(counts[target])} cells in this run: a group needs at least {MIN_GROUP_CELLS}.")
            return
        meta = dict(level=level, target=target, uniques=uniques, imputed=imputed,
                    note=(f" (only the first, {entries[0]}, is used)" if len(entries) > 1 else "") + note_source)
        data = self.ensure_means(imputed, meta)
        if data is not None:
            self.finish_means(data, meta)

    def ensure_means(self, imputed, meta):
        """The mean of every gene in every group, for every group type at once (so that changing the group type needs no new read): from memory,
        else from the file saved for this run, else read now in the background (None is returned and the plot follows when it is done)."""
        v = self.v
        by_level = {k: v.category_codes(k) for k in self.level_keys}
        digest = hashlib.md5(core._index_fingerprint(v.adata.obs_names).encode())
        for k, (c, _u) in by_level.items():
            digest.update(k.encode() + c.astype(np.int32).tobytes())
        fingerprint = digest.hexdigest()[:16]                  # (the cells, and the group each is in at every level)
        key = (imputed, fingerprint)
        if key in self.means_cache:
            return self.means_cache[key]
        uniques_of = {k: np.asarray([str(u) for u in u_], dtype=str) for k, (_c, u_) in by_level.items()}
        counts_of = {k: np.bincount(c[c >= 0], minlength=len(by_level[k][1])) for k, (c, _u) in by_level.items()}
        if imputed:
            source = v.imputed_state["adata"]
            path = source.filename
            positions = source.obs_names.get_indexer(v.adata.obs_names)
            names = source.var["gene_symbol"].astype(str).to_numpy()
            disk = core.CACHE_DIR / f"analyze_means_{fingerprint}.npz"
            loaded = load_saved_means(disk, uniques_of)
            if loaded is not None:
                self.means_cache[key] = loaded
                return loaded
            level_codes = {k: (c, len(u_)) for k, (c, u_) in by_level.items()}

            def work(progress, cancelled):
                means = group_means_imputed(path, positions, level_codes, progress, cancelled)
                if means is None:
                    return None
                data = {"names": names, "levels": {k: {"means": means[k], "uniques": uniques_of[k], "counts": counts_of[k]} for k in means}}
                save_means(disk, data)
                return data
        else:
            names = v.adata.var["gene_symbol"].astype(str).to_numpy()
            matrix = v.adata.layers["counts"] if v.own is None else None
            means_of = (lambda c, n: group_means_panel(matrix, c, n)) if v.own is None else v.own.group_means      # (a dataset of the user's own: its normalized values)

            def work(progress, cancelled):
                return {"names": names, "levels": {k: {"means": means_of(c, len(u_)), "uniques": uniques_of[k], "counts": counts_of[k]}
                                                   for k, (c, u_) in by_level.items()}}
        self.start_means(work, key, meta)
        return None

    def start_means(self, work, key, meta):
        """Read the means in a background thread (the window keeps going; Analyze is Cancel meanwhile), then draw."""
        state = {"cancel": False, "progress": 0.0, "done": False, "result": None, "error": None}
        self.job = state

        def run():
            try:
                state["result"] = work(lambda f: state.__setitem__("progress", f), lambda: state["cancel"])
            except BaseException as e:
                state["error"] = e
            state["done"] = True

        threading.Thread(target=run, daemon=True).start()
        self.btn_analyze.config(text="Cancel")
        label = "the imputed genes (the first time for this run: a minute or two)" if meta["imputed"] else "the genes"
        started = time.perf_counter()

        def poll():
            if not state["done"]:
                self.var_status.set(f"Reading {label}: {state['progress'] * 100:.0f}%  ({time.perf_counter() - started:.0f} s)")
                self.top.after(300, poll)
                return
            self.job = None
            self.btn_analyze.config(text="Analyze")
            if state["error"] is not None:
                self.var_status.set(f"Could not read the genes: {state['error']}")
            elif state["result"] is None:
                self.var_status.set("Cancelled.")
            else:
                self.means_cache[key] = state["result"]
                self.finish_means(state["result"], meta)

        poll()

    def finish_means(self, data, meta):
        self.sel = dict(meta, data=data)
        self.plot_selectivity()

    def plot_selectivity(self):
        s = self.sel
        level, target, uniques, imputed = s["level"], s["target"], s["uniques"], s["imputed"]
        data = s["data"]
        means, counts = data["levels"][level]["means"], data["levels"][level]["counts"]
        n_others = max(len(other_group_sizes(counts, target)), 1)
        k = other_groups_used(self.var_slider.get(), n_others)
        x, y, used = selectivity_values(means, counts, target, k)
        smallest = other_group_sizes(counts, target)[used - 1] if used else 0
        keep = np.isfinite(x) & np.isfinite(y) & (x > 0)
        all_names = data["names"]
        names, x, y = all_names[keep], x[keep], y[keep]
        fig, ax = self.fig, self.ax
        if self.cax is not None:
            self.cax.remove()
            self.cax = None
        ax.clear()
        ax.set_axis_on()
        self.home_view = None
        self.hover = None
        self.hover_group = None
        label = str(uniques[target]).strip()
        ax.scatter(x, y, s=SELECTIVITY_DOT_AREA[bool(imputed)], c="#3b6ea8", alpha=0.7, linewidths=0, rasterized=True)
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        marked, not_in_data, not_shown = [], [], []
        position = {str(n).casefold(): i for i, n in enumerate(names)}
        in_dataset = {str(n).casefold() for n in all_names}
        for g in core.parse_gene_names(self.var_sel_genes.get()):
            if g.casefold() in position:
                marked.append(position[g.casefold()])
            elif g.casefold() in in_dataset:
                not_shown.append(g)
            else:
                not_in_data.append(g)
        if not_in_data:
            self.sel_genes_entry.configure(foreground="red")
        if marked:
            ax.scatter(x[marked], y[marked], s=46, c=MARK_COLOR, linewidths=0.6, edgecolors="white", zorder=4)
            for i in marked:
                ax.annotate(names[i], (x[i], y[i]), xytext=(6, 4), textcoords="offset points", fontsize=LABEL_FONT_SIZE + 1, fontweight="bold",
                            color=MARK_COLOR, zorder=5)
        named = np.flatnonzero(y > LABEL_RATIO)                                            # (every gene whose ratio is above LABEL_RATIO is named)
        for i in named:
            if i not in marked:
                ax.annotate(names[i], (x[i], y[i]), xytext=(5, 3), textcoords="offset points", fontsize=LABEL_FONT_SIZE, color="#1d3b5c")
        own = self.v.own
        unit = ("mean imputed expression (linear)" if imputed else "mean normalized expression" if own is not None and not (own.kind == "counts" and own.normalization == "none")
                else "mean count per cell")
        ax.set_xlabel(f"{unit} in {label}")
        ax.set_ylabel(textwrap.fill(f"ratio to the highest of the {used} largest other {plural(self.v.level_names[level])}", 42))
        ax.set_title(f"{len(x):,} genes: {label}")
        ax.grid(True, alpha=0.25, which="both")
        self.hover = {"kind": "genes", "x": x, "y": y, "names": names}
        self.table = pd.DataFrame({"gene": names, ax.get_xlabel().replace("\n", " "): x, ax.get_ylabel().replace("\n", " "): y})
        self.pluses = ax.scatter([], [], marker="+", c="black", s=40, linewidths=0.9, animated=True, zorder=5)
        self.tip = ax.annotate("", xy=(0, 0), xytext=(14, 14), textcoords="offset points", fontsize=9, animated=True, zorder=6,
                               bbox={"boxstyle": "round,pad=0.3", "fc": "#ffffe0", "ec": "#888888"}, annotation_clip=False)
        self.tip.set_visible(False)
        ax.margins(0.08)
        self.legend_room = False
        self.place_axes()
        self.canvas.draw()
        self.remember_view()
        dropped = int((~keep).sum())
        extra = ""
        if not_in_data:
            extra += f" Not in the dataset: {', '.join(not_in_data)}."
        if not_shown:
            extra += f" None in {label}, so not on the graph: {', '.join(not_shown)}."
        self.var_status.set(f"{len(x):,} genes plotted" + (f" ({dropped:,} with none in {label} are left out)" if dropped else "")
                            + f", each against the {used} largest other {plural(self.v.level_names[level])} (at least {smallest:,} cells each). "
                            f"{len(named):,} with a ratio above {LABEL_RATIO:g} are named." + s.get("note", "") + extra
                            + " Scroll to zoom, double-click to reset.")

    # ---- zoom: the wheel, about the pointer; a double-click goes back (the toolbar's zoom and pan tools work as well)

    def remember_view(self):
        self.home_view = (self.ax.get_xlim(), self.ax.get_ylim())
        # (clearing the axes drops their callbacks, so this is done after each plot: zooming out stops at this view, and panning stays inside it)
        self.ax.callbacks.connect("xlim_changed", self.keep_in_view)
        self.ax.callbacks.connect("ylim_changed", self.keep_in_view)

    def keep_in_view(self, _ax):
        """The limits changed (the wheel, the toolbar's zoom or pan): the view may not be bigger than the first one drawn, nor leave it."""
        if self.home_view is None or self.clamping:
            return
        self.clamping = True
        try:
            for lim, home, scale, setter in ((self.ax.get_xlim(), self.home_view[0], self.ax.get_xscale(), self.ax.set_xlim),
                                             (self.ax.get_ylim(), self.home_view[1], self.ax.get_yscale(), self.ax.set_ylim)):
                clamped = clamp_limits(lim, home, scale)
                if clamped is not None:
                    setter(*clamped)
        finally:
            self.clamping = False

    def reset_view(self):
        if self.home_view is not None:
            self.ax.set_xlim(self.home_view[0])
            self.ax.set_ylim(self.home_view[1])
            self.place_axes()
            self.canvas.draw_idle()

    def on_press(self, ev):
        """A double-click goes back to the first view; a drag with the left button pans (unless one of the toolbar's own tools is on)."""
        if ev.inaxes is not self.ax or self.home_view is None:
            return
        if ev.dblclick:
            self.reset_view()
        elif ev.button == 1 and not self.toolbar.mode:
            self.pan_start = (ev.x, ev.y, self.ax.get_xlim(), self.ax.get_ylim())
            self.hide_hover()
            self.canvas.get_tk_widget().config(cursor="fleur")

    def on_release(self, ev):
        if self.pan_start is None:
            return
        self.pan_start = None
        self.canvas.get_tk_widget().config(cursor="")
        self.place_axes()                                         # (the tick labels may have changed width)
        self.canvas.draw_idle()

    def pan_to(self, ev):
        """The plot follows the pointer from where it was pressed (whatever the axes' scales; the view cannot leave the first one: keep_in_view)."""
        x0, y0, xlim, ylim = self.pan_start
        box = self.ax.bbox
        for lim, delta, size, scale, setter in ((xlim, ev.x - x0, box.width, self.ax.get_xscale(), self.ax.set_xlim),
                                                (ylim, ev.y - y0, box.height, self.ax.get_yscale(), self.ax.set_ylim)):
            a, b = (np.log10(lim[0]), np.log10(lim[1])) if scale == "log" else lim
            shift = -delta / max(size, 1.0) * (b - a)
            a, b = a + shift, b + shift
            setter(*((10 ** a, 10 ** b) if scale == "log" else (a, b)))
        self.canvas.draw_idle()

    def on_scroll(self, ev):
        if ev.inaxes is not self.ax or ev.xdata is None or self.home_view is None:
            return
        factor = 0.8 if ev.button == "up" else 1.25
        for lim, centre, setter, scale in ((self.ax.get_xlim(), ev.xdata, self.ax.set_xlim, self.ax.get_xscale()),
                                           (self.ax.get_ylim(), ev.ydata, self.ax.set_ylim, self.ax.get_yscale())):
            if scale == "log":
                a, b, c = np.log10(lim[0]), np.log10(lim[1]), np.log10(centre)
                setter(10 ** (c - (c - a) * factor), 10 ** (c + (b - c) * factor))
            else:
                setter(centre - (centre - lim[0]) * factor, centre + (lim[1] - centre) * factor)
        self.place_axes()                                         # (the tick labels change width: the margin follows)
        self.canvas.draw_idle()

    def place_axes(self):
        """The graph is a square, as big as the window allows: its side is the smaller of the room across and the room down (margins for the
        labels, and at the right for a legend when there is one)."""
        W, H = self.fig.get_size_inches() * self.fig.dpi
        grow = self.fig.dpi / 100.0                                       # (the text is bigger on a screen with scaling: so are the margins)
        left, bottom, top = 110.0 * grow, 65.0 * grow, 40.0 * grow
        right = (250.0 if self.legend_room else 30.0) * grow
        for _ in range(3):
            side = max(min(W - left - right, H - bottom - top), 50.0)
            y0 = bottom + (H - bottom - top - side) / 2                   # (centered up and down when the width is what limits it)
            self.ax.set_position([left / W, y0 / H, side / W, side / H])
            try:                                                          # (then what the axes' text really takes is measured, and the margins widened to fit)
                box = self.ax.get_tightbbox(self.canvas.get_renderer())
            except Exception:
                return
            over = (max(0.0, 4.0 - box.x0), max(0.0, 4.0 - box.y0), max(0.0, box.x1 - (W - 4.0)), max(0.0, box.y1 - (H - 4.0)))
            if max(over) < 1.0:
                return
            left, bottom, right, top = left + over[0], bottom + over[1], right + over[2], top + over[3]

    # ---- the pointer over a per-cell plot: the cell's group in a tip, and a plus on every cell of that group

    def on_draw(self, _event):
        """The figure was drawn: keep it as the background the tip and the plusses are blitted onto."""
        self.background = self.canvas.copy_from_bbox(self.fig.bbox)
        self.hover_group = None
        if self.hover is not None:
            self.pluses.set_offsets(np.empty((0, 2)))
            self.tip.set_visible(False)

    def hide_hover(self):
        if self.hover is None or self.hover_group is None or self.background is None:
            return
        self.hover_group = None
        self.canvas.restore_region(self.background)
        self.canvas.blit(self.fig.bbox)

    def on_motion(self, ev):
        if self.pan_start is not None:
            self.pan_to(ev)
            return
        h = self.hover
        if h is None or self.background is None:
            return
        now = time.perf_counter()
        if now - self.last_motion < 0.03:                      # (at most about thirty looks a second)
            return
        self.last_motion = now
        if ev.inaxes is not self.ax or ev.xdata is None:
            self.hide_hover()
            return
        if h.get("kind") == "genes":                              # (Selectivity: the gene under the pointer, whatever the axes' scales)
            pts = self.ax.transData.transform(np.column_stack([h["x"], h["y"]]))
            d2 = (pts[:, 0] - ev.x) ** 2 + (pts[:, 1] - ev.y) ** 2
            i = int(np.argmin(d2)) if len(d2) else -1
            if i < 0 or d2[i] > HOVER_PX ** 2:
                self.hide_hover()
                return
            self.canvas.restore_region(self.background)
            self.hover_group = i
            self.tip.set_text(f"{h['names'][i]}\n{h['x'][i]:.3g}, {h['y'][i]:.3g}")
            self.tip.xy = (h["x"][i], h["y"][i])
            self.tip.set_visible(True)
            self.ax.draw_artist(self.tip)
            self.canvas.blit(self.fig.bbox)
            return
        x0, x1 = self.ax.get_xlim()
        y0, y1 = self.ax.get_ylim()
        box = self.ax.bbox
        sx, sy = box.width / abs(x1 - x0), box.height / abs(y1 - y0)
        d2 = ((h["x"] - ev.xdata) * sx) ** 2 + ((h["y"] - ev.ydata) * sy) ** 2
        i = int(np.argmin(d2)) if len(d2) else -1
        if i < 0 or d2[i] > HOVER_PX ** 2:
            self.hide_hover()
            return
        g = h["group"][i]
        self.canvas.restore_region(self.background)
        if g != self.hover_group:
            self.hover_group = g
            same = h["group"] == g
            self.hover_count = int(same.sum())
            self.pluses.set_offsets(np.column_stack([h["x"][same], h["y"][same]]))
        lines = [f"{h['full'][g]}", f"{self.hover_count:,} cells"]
        for gene, value in zip(h["genes"], (h["ox"][i], h["oy"][i])):
            lines.append(f"{gene}: {both_scales(value, h['scale'])}")
        self.tip.set_text("\n".join(lines))
        self.tip.xy = (h["x"][i], h["y"][i])
        self.tip.set_visible(True)
        self.ax.draw_artist(self.pluses)
        self.ax.draw_artist(self.tip)
        self.canvas.blit(self.fig.bbox)
