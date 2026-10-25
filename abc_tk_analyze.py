"""The UMAP viewer's Analyze window: plot genes against each other, for the cells of some groups, per cell or per group.

It starts from what the viewer has (the group numbers, the genes, the level), and then goes its own way. The groups and the genes are
boxes with several entries separated by commas. "Expression" plots the second gene against the first; "Selectivity" is not made yet."""
import re
import time
import tkinter as tk
from tkinter import ttk

import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import abc_core as core
from abc_tk_umap import LEVEL_NAMES, rank_color
from abc_tk_umap_export import hexcolor

MAX_LABELED_GROUPS = 400         # a per-group plot names its points up to this many groups
LARGEST_GROUP_AREA = 1400.0      # a per-group plot: the marker area (points squared) of the group with the most cells; the others by their cell count
SMALLEST_AREA = 3.0              # (so that a tiny group is still a dot)
LABEL_DARKEN = 0.6               # the group numbers beside the dots: the circle's color times this
LABEL_FONT_SIZE = 8              # the group numbers beside the dots: the same for every group
DENSITY_BINS = 160               # the grid the density shading counts the cells on
HOVER_PX = 9                     # a per-cell plot: how close (pixels) the pointer must be to a cell
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


def plural(level):
    name = LEVEL_NAMES[level].lower()
    return name + ("es" if name.endswith("s") else "s")


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


def local_density(x, y, bins=DENSITY_BINS):
    """For each point, how many points share its cell of a bins x bins grid over the points."""
    h, xe, ye = np.histogram2d(x, y, bins=bins)
    ix = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, bins - 1)
    iy = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, bins - 1)
    return h[ix, iy]


def units(imputed, scale):
    if imputed:
        return "imputed, log2" if scale == "log" else "imputed, linear"
    return "log2(1 + counts)" if scale == "log" else "counts"


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

        box = ttk.LabelFrame(side, text="Groups", padding=6)
        box.pack(fill=tk.X)
        ttk.Label(box, text="Numbers or names, separated by commas. Empty: every group.", foreground="gray30", wraplength=260).pack(anchor="w")
        self.var_groups = tk.StringVar(top, viewer.var_ids.get())
        self.groups_entry = ttk.Entry(box, textvariable=self.var_groups)
        self.groups_entry.pack(fill=tk.X, pady=(2, 0))
        box = ttk.LabelFrame(side, text="Genes", padding=6)
        box.pack(fill=tk.X, pady=(8, 0))
        ttk.Label(box, text="Separated by commas: the first is plotted along x, the second up y.", foreground="gray30", wraplength=260).pack(anchor="w")
        self.var_genes = tk.StringVar(top, viewer.var_gene.get())
        self.genes_entry = ttk.Entry(box, textvariable=self.var_genes)
        self.genes_entry.pack(fill=tk.X, pady=(2, 0))
        for var, entry in ((self.var_groups, self.groups_entry), (self.var_genes, self.genes_entry)):
            var.trace_add("write", lambda *a, e=entry: e.configure(foreground=""))
            entry.bind("<Return>", lambda ev: self.expression())

        box = ttk.LabelFrame(side, text="Group type", padding=6)
        box.pack(fill=tk.X, pady=(8, 0))
        self.level_keys = [k for k in LEVEL_NAMES if k in viewer.levels]
        self.var_level = tk.StringVar(top, LEVEL_NAMES[viewer.level] if viewer.level in self.level_keys else LEVEL_NAMES[self.level_keys[0]])
        ttk.Combobox(box, textvariable=self.var_level, values=[LEVEL_NAMES[k] for k in self.level_keys], state="readonly").pack(fill=tk.X)

        box = ttk.LabelFrame(side, text="Plot", padding=6)
        box.pack(fill=tk.X, pady=(8, 0))
        self.var_per = tk.StringVar(top, "cell")
        ttk.Radiobutton(box, text="Per cell", value="cell", variable=self.var_per).pack(anchor="w")
        ttk.Radiobutton(box, text="Per group", value="group", variable=self.var_per).pack(anchor="w")
        self.var_jitter = tk.BooleanVar(top, True)
        self.var_density = tk.BooleanVar(top, False)
        ttk.Checkbutton(box, text="Jitter counts (cells with the same count are spread out)", variable=self.var_jitter,
                        command=self.redraw).pack(anchor="w", pady=(8, 0))
        ttk.Checkbutton(box, text="Shade per-cell dots by density", variable=self.var_density, command=self.redraw).pack(anchor="w")

        row = ttk.Frame(side)
        row.pack(fill=tk.X, pady=(10, 0))
        row.columnconfigure(0, weight=1)
        row.columnconfigure(1, weight=1)
        ttk.Button(row, text="Expression", command=self.expression).grid(row=0, column=0, sticky="ew", padx=(0, 3))
        sel = ttk.Button(row, text="Selectivity", state="disabled")
        sel.grid(row=0, column=1, sticky="ew", padx=(3, 0))
        self.var_status = tk.StringVar(top, "")
        ttk.Label(side, textvariable=self.var_status, foreground="gray30", wraplength=270, justify=tk.LEFT).pack(anchor="w", pady=(10, 0))

        right = ttk.Frame(body)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.fig = Figure(figsize=(9.0, 8.4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.text(0.5, 0.5, "Choose groups and genes, then click Expression.", ha="center", va="center", color="gray", transform=self.ax.transAxes)
        self.ax.set_axis_off()
        self.hover = None                       # per-cell plot: what the pointer can be over (see plot_expression)
        self.background = None
        self.last_motion = 0.0
        self.hover_group = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.last_plot = None                   # the arguments of the plot on the screen (a checkbox draws it again)
        self.cax = None
        self.legend_room = False                # the plot leaves room at the right for a legend
        self.canvas.mpl_connect("resize_event", lambda ev: self.place_axes())
        self.canvas.mpl_connect("draw_event", self.on_draw)
        self.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.canvas.mpl_connect("figure_leave_event", lambda ev: self.hide_hover())
        bar = ttk.Frame(right)
        bar.pack(side=tk.BOTTOM, fill=tk.X)
        NavigationToolbar2Tk(self.canvas, bar).update()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.draw()

    # ---- what the boxes ask for

    def level(self):
        names = {LEVEL_NAMES[k]: k for k in self.level_keys}
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
            self.var_status.set(f"No {LEVEL_NAMES[level].lower()} named {', '.join(bad)} in this run.")
            return None
        return codes, chosen, uniques, True

    def genes(self):
        """The two genes to plot, as spelled in the dataset, or None (the status says why)."""
        v = self.v
        names = core.parse_gene_names(self.var_genes.get())
        if len(names) < 2:
            self.genes_entry.configure(foreground="red")
            self.var_status.set("Enter two genes, separated by a comma.")
            return None
        imputed = v.mode == "imputed"
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
        imputed, scale, per = v.mode == "imputed", v.var_scale.get(), self.var_per.get()
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
        self.hover = None
        self.hover_group = None
        label_of = {g: (re.match(r"^\s*(\d+)", str(uniques[g])).group(1) if re.match(r"^\s*(\d+)", str(uniques[g])) else str(uniques[g])) for g in chosen}
        full_of = {g: str(uniques[g]).strip() for g in chosen}
        unit = units(imputed, scale)
        if per == "cell":
            cell_group = codes[in_group]
            px, py = vx[in_group], vy[in_group]
            note = ""
            if self.var_jitter.get():
                rng = np.random.default_rng(0)
                px, py = jittered(px, imputed, scale, rng), jittered(py, imputed, scale, rng)
                note = " Imputed values are continuous: not jittered." if imputed else " Counts jittered within their step."
            shade = self.var_density.get() and len(px) > 0
            if shade:
                from matplotlib.colors import LogNorm
                dens = local_density(px, py)
                order = np.argsort(dens, kind="stable")                      # (the dense dots on top)
                sc = ax.scatter(px[order], py[order], c=dens[order], cmap="viridis", norm=LogNorm(vmin=1, vmax=max(float(dens.max()), 2.0)), s=4,
                                linewidths=0, rasterized=True)
                self.cax = ax.inset_axes([1.03, 0.0, 0.03, 1.0])
                fig.colorbar(sc, cax=self.cax, label="cells in the same small square of the graph")
            else:
                colors = np.array([color_of[g] for g in cell_group]) if len(cell_group) else []
                ax.scatter(px, py, c=colors, s=4, alpha=POINT_ALPHA, linewidths=0, rasterized=True)
            self.hover = {"x": px, "y": py, "group": cell_group, "full": full_of, "color": color_of}
            self.pluses = ax.scatter([], [], marker="+", c="black", s=40, linewidths=0.9, animated=True, zorder=5)
            self.tip = ax.annotate("", xy=(0, 0), xytext=(14, 14), textcoords="offset points", fontsize=9, animated=True, zorder=6,
                                   bbox={"boxstyle": "round,pad=0.3", "fc": "#ffffe0", "ec": "#888888"}, annotation_clip=False)
            self.tip.set_visible(False)
            if len(chosen) <= MAX_LEGEND_GROUPS and not shade:
                from matplotlib.lines import Line2D
                handles = [Line2D([], [], marker="o", linestyle="", color=color_of[g], label=f"{full_of[g]} ({int((cell_group == g).sum()):,})")
                           for g in chosen if (cell_group == g).any()]
                if handles:
                    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, frameon=True, title=LEVEL_NAMES[level])
            ax.set_xlabel(f"{gx} ({unit})")
            ax.set_ylabel(f"{gy} ({unit})")
            ax.set_title(f"{int(in_group.sum()):,} cells in {len(set(cell_group.tolist()))} {plural(level)}")
            self.var_status.set(f"{int(in_group.sum()):,} cells plotted." + note)
        else:
            xs, ys, names, cols, counts = [], [], [], [], []
            self.hover = None
            for g in chosen:
                m = in_group & (codes == g)
                if not m.any():
                    continue
                xs.append(float(np.mean(vx[m])))
                ys.append(float(np.mean(vy[m])))
                names.append(label_of[g])
                cols.append(color_of[g])
                counts.append(int(m.sum()))
            areas = np.maximum(np.array(counts, dtype=float) / max(counts) * LARGEST_GROUP_AREA, SMALLEST_AREA) if counts else []
            ax.scatter(xs, ys, c=cols, s=areas, alpha=0.85, linewidths=0.4, edgecolors="black")
            if len(xs) <= MAX_LABELED_GROUPS:
                for x, y, n, a, c in zip(xs, ys, names, areas, cols):
                    ax.annotate(n, (x, y), xytext=(np.sqrt(a) / 2 + 2, 0), textcoords="offset points", fontsize=LABEL_FONT_SIZE,
                                ha="left", va="center", color=darker(c))
            ax.set_xlabel(f"mean {gx} ({unit})")
            ax.set_ylabel(f"mean {gy} ({unit})")
            ax.set_title(f"{len(xs)} {plural(level)}")
            self.var_status.set(f"{len(xs)} groups plotted (the mean of each, over {sum(counts):,} cells; the circle's area follows the cell count).")
        ax.grid(True, alpha=0.25)
        ax.margins(0.08)                                      # (room for the big dots and their numbers at the edges)
        self.legend_room = per == "cell" and (len(chosen) <= MAX_LEGEND_GROUPS or self.var_density.get())
        self.place_axes()
        self.canvas.draw()

    def place_axes(self):
        """The graph is a square, as big as the window allows: its side is the smaller of the room across and the room down (margins for the
        labels, and at the right for a legend when there is one)."""
        W, H = self.fig.get_size_inches() * self.fig.dpi
        left, bottom, top = 85.0, 65.0, 40.0
        right = 250.0 if self.legend_room else 30.0
        side = max(min(W - left - right, H - bottom - top), 50.0)
        y0 = bottom + (H - bottom - top - side) / 2                      # (centered up and down when the width is what limits it)
        self.ax.set_position([left / W, y0 / H, side / W, side / H])

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
            self.pluses.set_offsets(np.column_stack([h["x"][same], h["y"][same]]))
            self.tip.set_text(f"{h['full'][g]}\n{int(same.sum()):,} cells")
        self.tip.xy = (h["x"][i], h["y"][i])
        self.tip.set_visible(True)
        self.ax.draw_artist(self.pluses)
        self.ax.draw_artist(self.tip)
        self.canvas.blit(self.fig.bbox)
