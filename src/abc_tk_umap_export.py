"""The exports and 3D views of the Tk UMAP viewer (abc_tk_umap.py): Save UMAP, Save section maps, Export DEGs, Export expression, Export
coordinates, 3D view and 3D in browser. They are the same outputs as the matplotlib viewer's (the same file contents, and the same file
names as the ones proposed); the statistics, dialogs and 3D viewers are the browser's own functions, which it hands to the window as `host`.

`ExportsMixin` is mixed into UmapWindow: it reads the window's state (what is colored by what, the cells, the view) and writes files. Each
export first shows a Save dialog with the file name it would use (in the run's folder, or the folder of the last export), so that the name
can be changed; cancelling the dialog cancels the export.
"""
import math
import re
import threading
from pathlib import Path
from tkinter import filedialog, messagebox

import numpy as np
import pandas as pd

import abc_core as core

EXPORT_FONTSIZE = 14
EXPORT_AXES_PX = 6400                         # the UMAP's plotting square in a saved file, at the save DPI
EXPORT_MARGIN_IN = 2.5
EXPORT_LEGEND_MAX_ROWS = 35
SECTION_MAPS_FULL_WIDTH_PX = 4000            # a panel at the largest Section map size (one across) is this many pixels wide in a saved PNG; smaller sizes in proportion
PANEL_PAD_IN, PANEL_TITLE_IN = 0.08, 0.25     # the space round a panel in a saved section map figure, and above it for its title
PANELS_ACROSS_CHOICES = ("Auto (near-square)", "As on screen", "1", "2", "3", "4", "5", "6", "8", "10", "12")
JPG_ARGS = {"quality": 85, "subsampling": 2, "optimize": True}      # a saved JPG: quality 85, colors at half resolution: on a section map, 2 to 4 times smaller than the PNG
JPG_MAX_SIDE_PX = 65500                       # the JPEG format cannot hold a picture more than 65,535 pixels on a side
SVG_PICTURE_DPI = 200                         # the pictures inside a saved SVG (the blank and gray cells): a picture costs the same however many cells it holds, so a big one loses to vector dots below about 70,000 cells
UMAP_POINT_SIZE = 10                          # matplotlib's s (an area) of a dot in a saved UMAP
POINT_ALPHA = 0.75


def hexcolor(rgb):
    return "#%02x%02x%02x" % tuple(int(v) for v in rgb)


def parse_panels_across(text):
    """The choice typed or picked in the "Panels across" box: 'auto', 'screen', or a whole number of columns."""
    text = (text or "").strip().lower()
    if text.isdigit() and int(text) >= 1:
        return int(text)
    return "screen" if text.startswith("as on") else "auto"


def section_map_columns(n, panel_w_in, panel_h_in, choice, screen_cols):
    """The number of panels across in a saved section map figure. 'screen': as the pane has them. A whole number: that many. 'auto': the number
    that makes the figure closest to a square, and if it cannot be square, taller than wide is preferred to wider than tall (a small extra cost
    for wider: it decides between layouts about as far from a square, but never picks a tall thin strip over a wide one that is much closer to a
    square); a layout that leaves few empty places is preferred to one that leaves many. Only the panels count in the shape: the legend band under
    them goes as far past the square as it needs to. The panels' size is set by the Section map size slider, so the figure is as big as that makes it."""
    pad, title = PANEL_PAD_IN, PANEL_TITLE_IN
    if choice == "screen":
        return max(1, screen_cols)
    if isinstance(choice, int):
        return min(max(choice, 1), n)
    best = None
    for cols in range(1, n + 1):
        rows = -(-n // cols)
        width = cols * (panel_w_in + pad)
        aspect = width / (rows * (panel_h_in + pad + title))
        cost = abs(math.log(aspect)) + (0.12 if aspect > 1.0 else 0.0) + 0.15 * (cols * rows - n) / (cols * rows)
        if best is None or cost < best[0] - 1e-9:
            best = (cost, cols)
    return best[1]


def sibling(path, tail, ext):
    """`path`, a file name that ends in `ext`, with `tail` put before the ending. The name may have dots inside it (a UMAP's md0.5_pc100), so
    Path.with_suffix, which takes the text after the last dot for an ending, must not be used on it."""
    return path.with_name(path.name[:-len(ext)] + tail + ext)


def safe_token(text, default="view"):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_") or default


class ExportsMixin:
    # ---- shared helpers

    def target_dir(self):
        d = Path(self.run_folder) if self.run_folder is not None else Path(self.host["CACHE_DIR"])
        d.mkdir(parents=True, exist_ok=True)
        return d

    def ask_save_path(self, title, proposed, extension, kind, parent=None):
        """The file the user picks in a Save dialog, or None if the dialog was cancelled. The dialog offers `proposed` (a file name) in the run's
        folder, or in the folder of the last export of this window; a name typed without the extension gets it. The dialog itself asks before
        replacing a file."""
        folder = getattr(self, "export_dir", None)
        if folder is None or not Path(folder).is_dir():
            folder = self.target_dir()
        chosen = filedialog.asksaveasfilename(
            parent=parent or self.root, title=title, initialdir=str(folder), initialfile=proposed, defaultextension=extension,
            filetypes=[(kind, "*" + extension), ("All files", "*.*")])
        if not chosen:
            return None
        path = Path(chosen)
        if path.suffix.lower() != extension:       # (a name typed with some other ending keeps it, and gets this one after it)
            path = path.with_name(path.name + extension)
        self.export_dir = path.parent
        return path

    def ask_image_format(self, title, svg_note, panels_across=False, summary=None, parent=None):
        """"png", "jpg" or "svg", chosen in a small dialog over `parent` (the window by default; None is returned if it was cancelled or closed). The last
        choice made in this window is offered first. `svg_note` is added under the SVG choice."""
        import tkinter as tk
        from tkinter import ttk
        owner = parent or self.root
        dlg = tk.Toplevel(owner)
        dlg.title(title)
        dlg.transient(owner)
        dlg.resizable(False, False)
        var = tk.StringVar(dlg, getattr(self, "export_format", "png") if getattr(self, "export_format", "png") in ("png", "jpg", "svg") else "png")
        box = ttk.Frame(dlg, padding=16)
        box.pack()
        ttk.Label(box, text="Save as:").pack(anchor="w")
        ttk.Radiobutton(box, text="PNG: a picture", value="png", variable=var).pack(anchor="w", pady=(6, 0))
        png_info = tk.StringVar(dlg, "")                 # (what the picture will be, under its choice; the SVG's goes at the end of its description)
        if summary is not None:
            ttk.Label(box, textvariable=png_info, foreground="gray30", wraplength=340, justify=tk.LEFT).pack(anchor="w", padx=(22, 0))
        ttk.Radiobutton(box, text="JPG: a picture in a much smaller file", value="jpg", variable=var).pack(anchor="w", pady=(4, 0))
        jpg_info = tk.StringVar(dlg, "The same picture as the PNG, compressed: typically 2 to 4 times smaller, with slight blurring and color fringes at the edges of dots and text.")
        ttk.Label(box, textvariable=jpg_info, foreground="gray30", wraplength=340, justify=tk.LEFT).pack(anchor="w", padx=(22, 0))
        ttk.Radiobutton(box, text="SVG: vector graphics, for editing in a drawing program", value="svg", variable=var).pack(anchor="w", pady=(4, 0))
        svg_info = tk.StringVar(dlg, svg_note or "")
        if svg_note or summary is not None:
            ttk.Label(box, textvariable=svg_info, foreground="gray30", wraplength=340, justify=tk.LEFT).pack(anchor="w", padx=(22, 0))
        across = tk.StringVar(dlg, getattr(self, "export_panels_across", PANELS_ACROSS_CHOICES[0]))
        if panels_across:                          # (section maps: how the panels are arranged in the saved picture)
            line = ttk.Frame(box)
            line.pack(anchor="w", pady=(12, 0))
            ttk.Label(line, text="Panels across:").pack(side=tk.LEFT)
            ttk.Combobox(line, textvariable=across, values=PANELS_ACROSS_CHOICES, width=18).pack(side=tk.LEFT, padx=(8, 0))
            ttk.Label(box, text="Auto makes the panels close to a square (taller rather than wider); the legend is under them.", foreground="gray30",
                      wraplength=340, justify=tk.LEFT).pack(anchor="w")
            if summary is not None:                                       # (the layout and the size of each file the choice makes, updated as it is changed)
                def refresh(*args):
                    try:
                        png_text, svg_text, jpg_text = summary(across.get())
                    except Exception:
                        png_text, svg_text, jpg_text = "", "", ""
                    png_info.set(png_text)
                    svg_info.set(((svg_note + " ") if svg_note else "") + svg_text)
                    jpg_info.set("The same picture as the PNG, compressed: typically 2 to 4 times smaller, with slight blurring and color fringes at the edges of dots and text." + jpg_text)
                across.trace_add("write", refresh)
                refresh()
        result = []

        def finish(value):
            result.append(value)
            dlg.destroy()
        row = ttk.Frame(box)
        row.pack(anchor="e", pady=(14, 0))
        ttk.Button(row, text="OK", command=lambda: finish(var.get())).pack(side=tk.LEFT)
        ttk.Button(row, text="Cancel", command=lambda: finish(None)).pack(side=tk.LEFT, padx=(6, 0))
        dlg.bind("<Return>", lambda e: finish(var.get()))
        dlg.bind("<Escape>", lambda e: finish(None))
        dlg.protocol("WM_DELETE_WINDOW", lambda: finish(None))
        dlg.update_idletasks()
        try:                                        # (in the middle of the window, not in a corner of the screen)
            x = owner.winfo_rootx() + (owner.winfo_width() - dlg.winfo_reqwidth()) // 2
            y = owner.winfo_rooty() + (owner.winfo_height() - dlg.winfo_reqheight()) // 3
            dlg.geometry(f"+{max(x, 0)}+{max(y, 0)}")
            dlg.grab_set()
        except tk.TclError:
            pass
        owner.wait_window(dlg)
        choice = result[0] if result else None
        if choice:
            self.export_format = choice
            if panels_across:
                self.export_panels_across = across.get()
        return choice

    def ask_export_paths(self, fmt, stem, title, parent=None):
        """The file to save for the format `fmt` ("png", "jpg" or "svg"): {format: path}, or None if the Save dialog was cancelled."""
        ext = "." + fmt
        proposed = self.host["unique_export_paths"](self.target_dir(), stem, (ext,))[ext].name
        chosen = self.ask_save_path(title, proposed, ext, fmt.upper() + " image", parent=parent)
        return None if chosen is None else {fmt: chosen}

    def confirm_replacing(self, paths):
        """True unless some of `paths` (files the Save dialog did not ask about, saved along with the one it did) exist and the user declines
        to replace them."""
        existing = [p for p in paths if Path(p).exists()]
        if not existing:
            return True
        return messagebox.askyesno("Replace files?", "These files exist and would be replaced:\n\n" + "\n".join(str(p) for p in existing)
                                   + "\n\nReplace them?", icon="warning", parent=self.root)

    def begin_wait(self, message):
        """The hourglass at once, wherever the pointer is: the widgets that have a cursor of their own (the legend's hand, ...) get it too, since a
        window's cursor does not reach them. The pending events are handled first, so that the Save dialog that was just closed is gone from the screen
        before the long work begins. `end_wait` puts the cursors back."""
        self.wait_cursors = []
        stack = [self.root]
        while stack:
            w = stack.pop()
            for child in w.winfo_children():
                if child.winfo_class() == "Toplevel":
                    continue
                stack.append(child)
                try:
                    own = str(child.cget("cursor"))
                    if own:
                        self.wait_cursors.append((child, own))
                        child.config(cursor="watch")
                except Exception:
                    pass
        self.root.config(cursor="watch")
        self.var_status.set(message)
        self.root.update()

    def end_wait(self):
        for w, own in getattr(self, "wait_cursors", []):
            try:
                w.config(cursor=own)
            except Exception:
                pass
        self.wait_cursors = []
        self.root.config(cursor="")

    def view_token_safe(self):
        token = safe_token(getattr(self, "view_token", "view"))
        params_token = self.host["umap_params_token"](self.params) if self.params is not None else ""
        return f"{token}_{params_token}" if params_token else token

    def level_ids(self, level):
        return core.extract_leading_numeric_id(pd.Series(self.info[level])).to_numpy()

    def group_choice(self, action):
        """(level, targets) of the groups highlighted in Single group(s) that have cells in this run, or None (saying why)."""
        sel = getattr(self, "last_ids", None)
        if self.mode != "single" or not sel or not sel["targets"]:
            self.var_status.set(f"Enter at least one group number that exists in this run before {action}.")
            return None
        targets = [t for t in sel["targets"] if sel["counts"].get(t, 0) > 0]
        if not targets:
            self.var_status.set(f"Enter at least one group number that exists in this run before {action}.")
            return None
        return sel["level"], targets

    def background(self, text, work, done):
        """Run `work` (no Tk in it) off the main thread; `done(result)` runs back on it."""
        return self.job(text, work, done)

    # ---- Save UMAP

    def view_limits(self):
        c = self.umap
        w, h = c.size()
        return (c.cx - w / (2 * c.scale), c.cx + w / (2 * c.scale)), (c.cy - h / (2 * c.scale), c.cy + h / (2 * c.scale))

    def export_channels(self):
        """(the colours, the names) of the three genes' channels for the legend of a saved UMAP, in the theme on the screen."""
        if self.var_theme.get() == "light":
            return ("#00a8cc", "#d81b9a", "#d99a00") * 2, ("cyan", "magenta", "yellow") * 2
        return ("red", "green", "#4040ff", "red", "green", "#4040ff"), ("red", "green", "blue", "red", "green", "blue")

    def build_umap_figure(self):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.lines import Line2D
        h = self.host
        dpi = h["UMAP_SAVE_DPI"]
        axes_in = EXPORT_AXES_PX / dpi
        total = axes_in + 2 * EXPORT_MARGIN_IN
        fig = Figure(figsize=(total, total))
        FigureCanvasAgg(fig)
        frac, margin = axes_in / total, EXPORT_MARGIN_IN / total
        ax = fig.add_axes([margin, margin, frac, frac])
        ax.set_aspect("equal", adjustable="box")
        xlim, ylim = self.view_limits()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_facecolor(tuple(v / 255 for v in self.umap.bgcolor))                    # (the background on the screen, black to white)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("UMAP1", fontsize=EXPORT_FONTSIZE)
        ax.set_ylabel("UMAP2", fontsize=EXPORT_FONTSIZE)
        ax.set_title(getattr(self, "view_title", ""), fontsize=EXPORT_FONTSIZE)
        size = UMAP_POINT_SIZE * self.dot_scale() ** 2
        coords = self.coords
        kind = getattr(self, "view_kind", "all")
        entries = None
        if kind == "gene":
            gene, dataset, scale, vals = self.shown_values
            fin = np.isfinite(vals)
            vmin, vmax = getattr(self, "shown_range", None) or ((float(np.nanmin(vals)), float(np.nanmax(vals))) if fin.any() else (0.0, 1.0))      # (as shown: with the levels of the slider)
            order = np.argsort(np.where(fin, vals, -np.inf), kind="stable")
            cmap = self.gene_cmap()
            lowest = (fin & (vals <= vmin))[order]                              # (the cells with no expression: most of them, all alike)
            sc = None
            for part, raster in ((lowest, True), (~lowest, False)):
                if part.any():
                    sc = ax.scatter(coords[order[part], 0], coords[order[part], 1], c=vals[order[part]], cmap=cmap, vmin=vmin, vmax=vmax, s=size,
                                    alpha=POINT_ALPHA, linewidths=0, rasterized=raster)
            cb = fig.colorbar(sc, ax=ax)
            cb.set_label(gene, fontsize=EXPORT_FONTSIZE)
            cb.ax.tick_params(labelsize=EXPORT_FONTSIZE)
        elif kind == "multi":
            zero = np.array(self.zero_color() or (-1, -1, -1))
            for lay in self.umap.layers:
                blank = (lay["rgb"] == zero).all(axis=1)                          # (no expression of any gene: a picture in an SVG)
                dot = size * (0.25 if lay.get("shape") != "half" and lay["size"] < 4 * self.dot_scale() else 1.0)
                for part, raster in ((blank, True), (~blank, False)):
                    if part.any():
                        ax.scatter(lay["x"][part], lay["y"][part], c=lay["rgb"][part] / 255.0, s=dot, linewidths=0, alpha=None, rasterized=raster)
            handles = [Line2D([0], [0], marker="o", linestyle="", color=color, label=f"{g} ({name})")
                       for i, (g, color, name) in enumerate(zip(self.multi_genes, *self.export_channels()))
                       if i not in self.hidden]
            fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), bbox_transform=ax.transAxes, fontsize=EXPORT_FONTSIZE, frameon=True)
        else:
            cat = getattr(self, "export_cat", None)
            if cat is None:
                ax.scatter(coords[:, 0], coords[:, 1], c="lightgray", s=size, alpha=POINT_ALPHA, linewidths=0, rasterized=True)
            else:
                rest = cat["rank"] < 0
                ax.scatter(coords[rest, 0], coords[rest, 1], c="lightgray", s=size, alpha=POINT_ALPHA, linewidths=0, rasterized=True)
                for r, (color, label, count) in enumerate(cat["entries"]):
                    m = cat["rank"] == r
                    if not m.any():
                        continue
                    marker, is_open = h["category_rank_shape"](r)
                    if is_open:
                        ax.scatter(coords[m, 0], coords[m, 1], facecolors="none", edgecolors=color, marker=marker, s=size, alpha=POINT_ALPHA, linewidths=1.0)
                    else:
                        ax.scatter(coords[m, 0], coords[m, 1], c=color, marker=marker, s=size, alpha=POINT_ALPHA, linewidths=0)
                entries = cat
            from abc_tk_umap import ink_for
            by_count = self.var_label_scale.get()                 # (the "Label font size scales with cell count" checkbox: the saved UMAP follows it)
            for text, x, y, color, weight in getattr(self, "label_info", []):
                if xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]:
                    ax.text(x, y, text, fontsize=EXPORT_FONTSIZE * 0.6 * (weight if by_count else 1.0), color=ink_for(color, self.umap.bgcolor), fontweight="bold", ha="center", va="center", clip_on=True)
            if entries is not None:
                handles = self.legend_handles(entries, sort="id")
                ncol = max(1, math.ceil(len(handles) / EXPORT_LEGEND_MAX_ROWS))
                fig.legend(handles=handles, title=entries["title"], loc="center left", bbox_to_anchor=(1.02, 0.5), bbox_transform=ax.transAxes,
                           fontsize=EXPORT_FONTSIZE * 0.85, title_fontsize=EXPORT_FONTSIZE, ncol=ncol, frameon=True, handlelength=1.0,
                           handleheight=1.0, labelspacing=0.4, columnspacing=1.5)
        return fig

    def section_map_plan(self, choice):
        """What a saved section map will be for this choice of panels across ('auto', 'screen' or a number): the layout, and the size of the picture.
        A panel's size in the picture follows the Section map size slider and nothing else (not the window or the screen): at full width it is
        SECTION_MAPS_FULL_WIDTH_PX pixels wide, at one half of the width half of that, and so on; its shape is the shape of the panels on the screen."""
        dpi = self.host["SECTION_MAPS_SAVE_DPI"]
        n = len(self.tiles)
        screen_w, screen_h = self.tiles[0]["pc"].size()
        pw_px = round(SECTION_MAPS_FULL_WIDTH_PX / max(1, self.tile_cols))
        ph_px = round(pw_px * screen_h / screen_w)
        pw_in, ph_in = pw_px / dpi, ph_px / dpi
        legend = self.section_legend()                                  # (what the view on the screen says its colors are: a band under the panels)
        cols = section_map_columns(n, pw_in, ph_in, choice, max(1, self.tile_cols))
        rows = -(-n // cols)
        width = cols * (pw_in + PANEL_PAD_IN)
        band_h = self.legend_band_height(legend, width)
        height = rows * (ph_in + PANEL_PAD_IN + PANEL_TITLE_IN) + band_h
        return {"cols": cols, "rows": rows, "pw_in": pw_in, "ph_in": ph_in, "width": width, "height": height, "band_h": band_h, "legend": legend,
                "panel_px": (pw_px, ph_px), "px": (round(width * dpi), round(height * dpi))}

    def section_map_summary(self, choice_text):
        """(the PNG's line, the SVG's line) the save dialog shows for this choice of panels across: the layout and the size of each file."""
        plan = self.section_map_plan(parse_panels_across(choice_text))
        w_px, h_px = plan["px"]
        mp = w_px * h_px / 1e6
        png = (f"{plan['cols']} panel{'s' if plan['cols'] > 1 else ''} across and {plan['rows']} down" + (", the legend under them" if plan["band_h"] else "") + ". "
               f"Each panel is {plan['panel_px'][0]:,} x {plan['panel_px'][1]:,} pixels (set by the Section map size slider, not by the window). "
               f"The picture is {w_px:,} x {h_px:,} pixels ({mp:.0f} megapixels, about {w_px * h_px * 4 / 1e6:,.0f} MB of memory to make).")
        if mp > 250:
            png += " That is very large: use more panels across, or a smaller Section map size."
        svg = f"The page is {plan['width']:.1f} x {plan['height']:.1f} inches."
        jpg = (f" This one is too big for a JPG, which can be at most {JPG_MAX_SIDE_PX:,} pixels on a side: use more panels across, or a smaller Section map size."
               if max(w_px, h_px) > JPG_MAX_SIDE_PX else "")
        return png, svg, jpg

    def gene_cmap(self):
        """The colormap of a single gene in the theme on the screen."""
        from matplotlib import colormaps
        from matplotlib.colors import LinearSegmentedColormap
        from abc_tk_umap import LIGHT_STOPS
        if self.var_theme.get() == "light":
            return LinearSegmentedColormap.from_list("light_gene", [(s[0], tuple(v / 255 for v in s[1])) for s in LIGHT_STOPS])
        return colormaps["viridis"]

    def section_legend(self):
        """What the legend band of a saved section map holds, as the view on the screen has it (None if there is nothing to say):
        {"kind": "groups", "title", "handles", "labels"} | {"kind": "colorbar", ...} | {"kind": "genes", "rows"}."""
        kind = getattr(self, "view_kind", "all")
        cat = getattr(self, "export_cat", None)
        if kind in ("all", "single") and cat is not None:
            handles = self.legend_handles(cat, sort="screen" if kind == "all" else "given")
            return {"kind": "groups", "title": cat["title"], "handles": handles, "labels": [h.get_label() for h in handles]}
        if kind == "gene" and getattr(self, "shown_values", None) is not None:
            gene, dataset, scale, vals = self.shown_values
            vmin, vmax = getattr(self, "shown_range", None) or ((float(np.nanmin(vals)), float(np.nanmax(vals))) if np.isfinite(vals).any() else (0.0, 1.0))
            return {"kind": "colorbar", "gene": gene, "units": self.units_text(dataset, scale), "vmin": vmin, "vmax": vmax, "cmap": self.gene_cmap()}
        if kind == "multi" and getattr(self, "multi_ranges", None):
            colors, names = self.export_channels()
            linear = getattr(self, "multi_scale", "log") == "linear"
            rows = [(g, colors[i], names[i], f"{'' if linear else 'log2: '}{lo:.1f}\u2013{hi:.1f}") for i, (g, (lo, hi)) in enumerate(zip(self.multi_genes, self.multi_ranges))
                    if i not in self.hidden]
            return {"kind": "genes", "rows": rows}
        return None

    def legend_grid(self, spec, width):
        """(columns, rows) of the entries of a group legend in a band `width` inches wide."""
        chars = max((len(s) for s in spec["labels"]), default=10)
        ncol = max(1, int((width - 0.5) // (0.55 + 0.072 * chars)))
        return ncol, -(-len(spec["labels"]) // ncol)

    def legend_band_height(self, spec, width):
        """The height in inches of the band under the panels for this legend (0: none)."""
        if spec is None:
            return 0.0
        if spec["kind"] == "groups":
            return 0.55 + 0.2 * self.legend_grid(spec, width)[1]
        if spec["kind"] == "colorbar":
            return 1.3
        return 0.5 + 0.26 * len(spec["rows"])

    def draw_legend_band(self, fig, spec, fw, fh, band_h):
        """The legend in the band along the bottom of the figure (band_h inches high)."""
        from matplotlib.lines import Line2D
        ax = fig.add_axes([0, 0, 1, band_h / fh])
        ax.axis("off")
        if spec["kind"] == "groups":
            ncol = self.legend_grid(spec, fw)[0]
            ax.legend(handles=spec["handles"], loc="upper left", bbox_to_anchor=(0.01, 0.98), ncol=ncol, frameon=False, fontsize=10, title=spec["title"],
                      title_fontsize=11, handlelength=1.0, handleheight=1.0, labelspacing=0.4, columnspacing=1.5, alignment="left")
        elif spec["kind"] == "colorbar":
            from matplotlib.cm import ScalarMappable
            from matplotlib.colors import Normalize
            cax = fig.add_axes([0.5 / fw, 0.45 / fh, min(0.4 * fw, 8.0) / fw, 0.2 / fh])
            cb = fig.colorbar(ScalarMappable(norm=Normalize(spec["vmin"], spec["vmax"]), cmap=spec["cmap"]), cax=cax, orientation="horizontal")
            cb.set_label(f"{spec['gene']}   ({spec['units']})", fontsize=11)
            cb.ax.tick_params(labelsize=10)
        else:
            ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=color, label=f"{g} ({name})   {rng}") for g, color, name, rng in spec["rows"]],
                      loc="upper left", bbox_to_anchor=(0.01, 0.98), frameon=False, fontsize=10, alignment="left")

    def legend_handles(self, cat, sort="id"):
        from matplotlib.lines import Line2D
        rows = []
        for r, (color, label, count) in enumerate(cat["entries"]):
            marker, is_open = self.host["category_rank_shape"](r)
            rows.append((color, label, count, marker, is_open))
        if sort == "given":                                              # (as the entries are: the groups in the order they were typed)
            pass
        elif sort == "screen":                                           # (the order of the list on the screen)
            place = {k: i for i, k in enumerate(getattr(self, "group_rows", []))}
            rows.sort(key=lambda e: place.get(e[1], math.inf))
        elif sort == "count":
            rows.sort(key=lambda e: e[2], reverse=True)
        else:
            def key(e):
                m = re.match(r"^(\d+)", e[1])
                return int(m.group(1)) if m else math.inf
            rows.sort(key=key)
        return [Line2D([], [], marker=m, linestyle="none", markersize=9, markerfacecolor="none" if o else c, markeredgecolor=c,
                       markeredgewidth=1.2 if o else 0.5, label=f"{lab}  ({n} cells)") for c, lab, n, m, o in rows]

    def save_umap(self):
        if self.busy:
            return
        token, dpi = self.view_token_safe(), self.host["UMAP_SAVE_DPI"]
        fmt = self.ask_image_format("Save the UMAP", "The gray cells, and in a gene view the cells with no expression, are pictures inside it, to keep the file small; "
                                                     "the colored cells stay vector graphics.")
        if fmt is None:
            self.var_status.set("Save UMAP cancelled.")
            return
        cat = getattr(self, "export_cat", None)
        has_legend = getattr(self, "view_kind", "all") in ("all", "single") and cat is not None
        paths = self.ask_export_paths(fmt, f"umap_{token}", "Save the UMAP" + (" (its legend is saved with it, as _legend)" if has_legend else ""))
        if paths is None:
            self.var_status.set("Save UMAP cancelled.")
            return
        legend_paths = {f: sibling(path, "_legend", "." + f) for f, path in paths.items()}
        if has_legend and not self.confirm_replacing(list(legend_paths.values())):
            self.var_status.set("Save UMAP cancelled.")
            return
        self.begin_wait("Saving the UMAP...")
        try:
            fig = self.build_umap_figure()
            written = []
            for f, path in paths.items():
                self.var_status.set(f"Saving {f}...")
                self.root.update()
                fig.savefig(path, bbox_inches="tight", dpi=SVG_PICTURE_DPI if f == "svg" else dpi, **({"pil_kwargs": JPG_ARGS} if f == "jpg" else {}))          # (for an SVG the dpi is that of the pictures in it)
                written.append(path)
                if has_legend:
                    from matplotlib.figure import Figure
                    from matplotlib.backends.backend_agg import FigureCanvasAgg
                    handles = self.legend_handles(cat, sort="count")
                    lf = Figure()
                    FigureCanvasAgg(lf)
                    lf.legend(handles=handles, title=cat["title"], loc="center", fontsize=EXPORT_FONTSIZE * 0.85, title_fontsize=EXPORT_FONTSIZE,
                              ncol=max(1, math.ceil(len(handles) / EXPORT_LEGEND_MAX_ROWS)), frameon=True, handlelength=1.0, handleheight=1.0,
                              labelspacing=0.4, columnspacing=1.5)
                    lf.savefig(legend_paths[f], bbox_inches="tight", dpi=None if f == "svg" else dpi, **({"pil_kwargs": JPG_ARGS} if f == "jpg" else {}))
                    written.append(legend_paths[f])
            folders = list(dict.fromkeys(str(p.parent) for p in written))
            self.var_status.set(f"Saved {', '.join(p.name for p in written)} to {' and '.join(folders)}.")
        except Exception as e:
            self.var_status.set(f"Could not save the UMAP: {e}")
            print(f"Could not save the UMAP: {e}", flush=True)
        self.end_wait()

    # ---- Save section maps

    def save_section_maps(self):
        if self.busy or not self.tiles:
            if not self.tiles:
                self.var_status.set("No section panels to save.")
            return
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.patches import Rectangle
        h = self.host
        fmt = self.ask_image_format("Save the section maps", "The gray cells, and in a gene view the cells with no expression, are pictures inside it, to keep the file small; the colored cells and the "
                                                              "outlines stay vector graphics.", panels_across=True, summary=self.section_map_summary)
        if fmt is None:
            self.var_status.set("Save section maps cancelled.")
            return
        if fmt == "jpg":
            big = max(self.section_map_plan(parse_panels_across(getattr(self, "export_panels_across", "")))["px"])
            if big > JPG_MAX_SIDE_PX:
                self.var_status.set(f"A JPG can be at most {JPG_MAX_SIDE_PX:,} pixels on a side and this picture would be {big:,}: use more panels across, "
                                    "a smaller Section map size, or a PNG.")
                return
        roi_only = "_ROI_only" if self.var_roi_only.get() else ""        # (the ROI only box above the panels: the name says what the maps show)
        ccf = "_CCF-overlay" if any(t["pc"].outlines for t in self.tiles) else ""        # (and the CCF atlas box)
        paths = self.ask_export_paths(fmt, f"{self.view_token_safe()}_section_maps{roi_only}{ccf}", "Save the section maps")
        if paths is None:
            self.var_status.set("Save section maps cancelled.")
            return
        self.begin_wait("Saving the section maps...")
        try:
            written = []
            n = len(self.tiles)
            pw_px, ph_px = self.tiles[0]["pc"].size()
            plan = self.section_map_plan(parse_panels_across(getattr(self, "export_panels_across", "")))
            legend, cols, rows, pw_in, ph_in, band_h = plan["legend"], plan["cols"], plan["rows"], plan["pw_in"], plan["ph_in"], plan["band_h"]
            pad, row_extra = PANEL_PAD_IN, PANEL_TITLE_IN
            for fmt, chosen in paths.items():             # (each format is made on its own: the SVG holds its gray cells as pictures)
                fig = Figure(figsize=(cols * (pw_in + pad), rows * (ph_in + pad + row_extra) + band_h))
                FigureCanvasAgg(fig)
                fw, fh = fig.get_size_inches()
                for k, t in enumerate(self.tiles):
                    if t.get("dirty") or not t["pc"].layers:
                        self.build_tile(t)
                    pc = t["pc"]
                    r, c = divmod(k, cols)
                    x0 = c * (pw_in + pad) / fw
                    y0 = 1 - (r + 1) * (ph_in + pad + row_extra) / fh
                    ax = fig.add_axes([x0, y0, pw_in / fw, ph_in / fh])
                    ax.set_aspect("equal", adjustable="box")
                    w, hh = pc.size()
                    xl = (pc.cx - w / (2 * pc.scale), pc.cx + w / (2 * pc.scale))
                    yl = (pc.cy - hh / (2 * pc.scale), pc.cy + hh / (2 * pc.scale))
                    ax.set_xlim(xl)
                    ax.set_ylim(yl[1], yl[0])                            # (a section's y grows downward)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_facecolor(tuple(v / 255 for v in pc.bgcolor))
                    ax.set_title(t["title"].get(), fontsize=9, loc="left", pad=3)
                    pts_per_mm = 72.0 * pw_in / (xl[1] - xl[0])
                    from abc_tk_umap import DOT_MM, GREY, LIGHT_GREY, SHAPE_MARKERS
                    zero_rgb = np.array(self.zero_color()) if self.zero_color() is not None else None
                    for lay in pc.layers:
                        dot_mm = lay.get("size_mm", DOT_MM) * pc.dot_factor          # (a dot's width on the tissue, as on the screen: the Cell size slider)
                        m = (lay["x"] >= xl[0]) & (lay["x"] <= xl[1]) & (lay["y"] >= yl[0]) & (lay["y"] <= yl[1])
                        x, y, rgb = lay["x"][m], lay["y"][m], lay["rgb"][m]
                        shown = rgb if pc.cell_alpha >= 0.999 else np.rint(rgb * pc.cell_alpha + np.asarray(pc.bgcolor) * (1.0 - pc.cell_alpha))      # (dimmed, as on the screen)
                        marks = lay["marks"][m] if lay.get("marks") is not None else None         # (a marker code for each cell: the shape of its group)
                        # (in an SVG the gray cells, many thousands of them, are pictures: a vector dot apiece would make the file huge. They are told by
                        # their color, so that the gray ones of a layer that is mostly colored, or the other way round, are found too.)
                        if fmt == "svg" and not lay.get("gray"):
                            gray = (rgb == GREY).all(axis=1) | (rgb == LIGHT_GREY).all(axis=1)
                            if zero_rgb is not None:
                                gray |= (rgb == zero_rgb).all(axis=1)                  # (and the cells with no expression, in a gene view)
                            parts = ((gray, True), (~gray, False))
                        else:
                            parts = ((np.ones(len(x), bool), fmt == "svg" and bool(lay.get("gray"))),)
                        for sel, raster in parts:
                            if not sel.any():
                                continue
                            for code in (np.unique(marks[sel]) if marks is not None else [0]):
                                sub = sel & (marks == code) if marks is not None else sel
                                marker, is_open = SHAPE_MARKERS[int(code) % 10], int(code) >= 10
                                if is_open:
                                    ax.scatter(x[sub], y[sub], facecolors="none", edgecolors=shown[sub] / 255.0, marker=marker, s=(dot_mm * pts_per_mm) ** 2,
                                               linewidths=0.5, rasterized=raster)
                                else:
                                    ax.scatter(x[sub], y[sub], c=shown[sub] / 255.0, marker=marker, s=(dot_mm * pts_per_mm) ** 2, linewidths=0,
                                               rasterized=raster)
                    if pc.outlines:                                      # (the CCF atlas box: the regions' outlines, as on the screen)
                        from matplotlib.collections import LineCollection
                        seen = [np.vstack([p, p[:1]]) for p in pc.outlines
                                if p[:, 0].max() >= xl[0] and p[:, 0].min() <= xl[1] and p[:, 1].max() >= yl[0] and p[:, 1].min() <= yl[1]]
                        ink = "#d7d7d7" if sum(pc.bgcolor) < 300 else "#5a5a5a"
                        ax.add_collection(LineCollection(seen, colors=ink, linewidths=0.4 * pc.outline_px))        # (as wide as on the screen)
                    for _it, a, b, c2, d in pc.rects:
                        ax.add_patch(Rectangle((a, b), c2 - a, d - b, fill=False, edgecolor="#ff8c00", linestyle=(0, (4, 2)), linewidth=1.2))
                    if pc.scalebar:
                        from abc_tk_umap import SCALEBAR_UM
                        view_mm = xl[1] - xl[0]
                        cands = [u for u in SCALEBAR_UM if u <= view_mm * 1000 * 0.22]
                        um = cands[-1] if cands else SCALEBAR_UM[0]
                        ink = "white" if sum(pc.bgcolor) < 300 else "black"
                        bx0 = xl[0] + 0.06 * view_mm
                        by = yl[1] - 0.06 * (yl[1] - yl[0])
                        ax.plot([bx0, bx0 + um / 1000], [by, by], color=ink, linewidth=2.5, solid_capstyle="butt")
                        ax.text(bx0 + um / 2000, by - 0.01 * view_mm, f"{um:g} µm" if um < 1000 else f"{um / 1000:g} mm", color=ink, fontsize=8,
                                ha="center", va="bottom")
                if legend is not None:
                    self.draw_legend_band(fig, legend, fw, fh, band_h)
                self.var_status.set(f"Saving {fmt}...")
                self.root.update()
                fig.savefig(chosen, dpi=h["SECTION_MAPS_SAVE_DPI"], **({"pil_kwargs": JPG_ARGS} if fmt == "jpg" else {}))         # (for an SVG the dpi is that of the pictures in it)
                written.append(chosen)
            folders = list(dict.fromkeys(str(p.parent) for p in written))
            self.var_status.set(f"Saved {', '.join(p.name for p in written)} to {' and '.join(folders)} "
                                f"({cols} panels across, {rows} down; {fw:.0f} x {fh:.0f} in).")
        except Exception as e:
            self.var_status.set(f"Could not save section maps: {e}")
            print(f"Could not save section maps: {e}", flush=True)
        self.end_wait()

    # ---- Export DEGs / expression

    def imputed_for_export(self):
        """The imputed dataset (loaded now, with the browser's own prompt, if it is not), or None."""
        if self.imputed_state["adata"] is None and self.ensure_imputed:
            self.ensure_imputed(self.imputed_state, self.abc_cache)
        return self.imputed_state["adata"]

    def ask_csv_paths(self, what, prefix, token, has_counts, has_imputed):
        """{"merfish": path, "imputed": path} for the .csv files an export of `what` will write: a Save dialog for each of the two datasets that
        has data (a file is only written if its dataset has something to put in it), each proposing `<prefix>_<token>_<dataset>.csv`. None if a
        dialog was cancelled (nothing is exported then)."""
        paths = {}
        for key, present in (("merfish", has_counts), ("imputed", has_imputed)):
            if not present:
                continue
            proposed = self.host["unique_export_paths"](self.target_dir(), f"{prefix}_{token}_{key}", (".csv",))[".csv"].name
            path = self.ask_save_path(f"Save the {what} of the {'MERFISH' if key == 'merfish' else 'imputed'} dataset", proposed, ".csv", "CSV file")
            if path is None:
                return None
            paths[key] = path
        return paths

    def export_degs(self):
        if self.busy:
            return
        chosen = self.group_choice("exporting DE genes")
        if chosen is None:
            return
        level, targets = chosen
        h = self.host
        level_name = LEVEL_LONG[level]
        choice = h["prompt_deg_comparison_mode"](self.root, f"DEGs for {level_name} {', '.join(str(t) for t in targets)}", level_name)
        if choice is None:
            self.var_status.set("DEG export cancelled.")
            return
        mode, reference = choice["mode"], choice["reference"]
        eff = [t for t in targets if t != reference] if mode == "reference" else list(targets)
        if not eff:
            self.var_status.set(f"Every entered {level_name.lower()} ID equals the reference ({reference}); nothing to compare.")
            return
        ids_arr = self.level_ids(level)
        adata, sc, has_counts = self.adata, h["sc"], self.has_counts
        imputed_full = self.imputed_for_export()
        token = safe_token(f"{level}_{'-'.join(str(t) for t in eff)}" + (f"_vsRef{reference}" if mode == "reference" else ""), level)
        paths = self.ask_csv_paths("DE genes", "DEG", token, has_counts, imputed_full is not None)
        if paths is None:
            self.var_status.set("DEG export cancelled.")
            return

        def work():
            written, skipped = [], []
            if has_counts:
                m = adata.copy()
                m.X = m.layers["counts"].copy()
                sc.pp.normalize_total(m)
                sc.pp.log1p(m)
                frames = []
                for t in eff:
                    if mode == "reference":
                        pair = (ids_arr == t) | (ids_arr == reference)
                        comp, gm = m[pair], ids_arr[pair] == t
                    else:
                        comp, gm = m, ids_arr == t
                    if not gm.any() or gm.all():
                        continue
                    df = h["compute_de_genes"](comp, gm, log1p_base=None)
                    df.insert(0, "group_id", t)
                    if mode == "reference":
                        df.insert(1, "reference_id", reference)
                    frames.append(df)
                if frames:
                    out = pd.concat(frames, ignore_index=True)
                    path = paths["merfish"]
                    out.to_csv(path, index=False)
                    written.append((path, len(out), len(frames)))
                else:
                    skipped.append("MERFISH (no specified ID had both its own cells and a group to compare with)")
            else:
                skipped.append("MERFISH (no counts layer on this dataset)")
            if imputed_full is not None:
                pos = imputed_full.obs_names.get_indexer(adata.obs_names)
                found = pos >= 0
                if found.any():
                    sub = h["materialize_subset"](imputed_full, pos[found])
                    iid = ids_arr[found]
                    frames = []
                    for t in eff:
                        if mode == "reference":
                            pair = (iid == t) | (iid == reference)
                            comp, gm = sub[pair], iid[pair] == t
                        else:
                            comp, gm = sub, iid == t
                        if not gm.any() or gm.all():
                            continue
                        df = h["compute_de_genes"](comp, gm, log1p_base=2)
                        df.insert(0, "group_id", t)
                        if mode == "reference":
                            df.insert(1, "reference_id", reference)
                        frames.append(df)
                    if frames:
                        out = pd.concat(frames, ignore_index=True)
                        path = paths["imputed"]
                        out.to_csv(path, index=False)
                        written.append((path, len(out), len(frames)))
                    else:
                        skipped.append("imputed (no specified ID had both its own cells and a group to compare with)")
                else:
                    skipped.append("imputed (none of this run's cells found in the imputed dataset)")
            else:
                skipped.append("imputed (dataset not loaded)")
            return written, skipped

        def done(res):
            if isinstance(res, BaseException):
                self.var_status.set(f"Could not export DE genes: {res}")
                return
            written, skipped = res
            if written:
                summary = "; ".join(f"{p.name} ({n} rows, {k} of {len(eff)} IDs)" for p, n, k in written)
                msg = (f"Exported {summary}" + (f" vs. reference {level_name.lower()} {reference}" if mode == "reference" else "")
                       + f" to {', '.join(sorted({str(p.parent) for p, _n, _k in written}))}.")
                if skipped:
                    msg += f" Skipped: {'; '.join(skipped)}."
                self.var_status.set(msg)
            else:
                self.var_status.set(f"Could not export DE genes: {'; '.join(skipped) or 'nothing to export'}.")

        self.background("Computing DE genes...", work, done)

    def export_expression(self):
        if self.busy:
            return
        chosen = self.group_choice("exporting expression")
        if chosen is None:
            return
        level, targets = chosen
        h = self.host
        ids_arr = self.level_ids(level)
        adata, sc, has_counts = self.adata, h["sc"], self.has_counts
        imputed_full = self.imputed_for_export()
        token = safe_token(f"{level}_{'-'.join(str(t) for t in targets)}", level)
        paths = self.ask_csv_paths("expression", "EXPR", token, has_counts, imputed_full is not None)
        if paths is None:
            self.var_status.set("Expression export cancelled.")
            return

        def work():
            written, skipped = [], []
            if has_counts:
                m = adata.copy()
                m.X = m.layers["counts"].copy()
                sc.pp.normalize_total(m)
                sc.pp.log1p(m)
                df = pd.DataFrame({"gene_symbol": m.var["gene_symbol"].values, "ensembl_id": m.var_names})
                included = []
                for t in targets:
                    gm = ids_arr == t
                    if gm.any():
                        df[f"mean_log2_expr_{t}"] = h["mean_log2_expression"](m, gm, log1p_base=None).values
                        included.append(t)
                if included:
                    df = df.sort_values("gene_symbol", kind="stable").reset_index(drop=True)
                    path = paths["merfish"]
                    df.to_csv(path, index=False)
                    written.append((path, len(df), len(included)))
                else:
                    skipped.append("MERFISH (no specified ID had any matching cells)")
            else:
                skipped.append("MERFISH (no counts layer on this dataset)")
            if imputed_full is not None:
                pos = imputed_full.obs_names.get_indexer(adata.obs_names)
                found = pos >= 0
                if found.any():
                    sub = h["materialize_subset"](imputed_full, pos[found])
                    iid = ids_arr[found]
                    df = pd.DataFrame({"gene_symbol": sub.var["gene_symbol"].values, "ensembl_id": sub.var_names})
                    included = []
                    for t in targets:
                        gm = iid == t
                        if gm.any():
                            df[f"mean_log2_expr_{t}"] = h["mean_log2_expression"](sub, gm, log1p_base=2).values
                            included.append(t)
                    if included:
                        df = df.sort_values("gene_symbol", kind="stable").reset_index(drop=True)
                        path = paths["imputed"]
                        df.to_csv(path, index=False)
                        written.append((path, len(df), len(included)))
                    else:
                        skipped.append("imputed (no specified ID had any matching cells)")
                else:
                    skipped.append("imputed (none of this run's cells found in the imputed dataset)")
            else:
                skipped.append("imputed (dataset not loaded)")
            return written, skipped

        def done(res):
            if isinstance(res, BaseException):
                self.var_status.set(f"Could not export expression: {res}")
                return
            written, skipped = res
            if written:
                msg = ("Exported " + "; ".join(f"{p.name} ({n} rows, {k} of {len(targets)} IDs)" for p, n, k in written)
                       + f" to {', '.join(sorted({str(p.parent) for p, _n, _k in written}))}.")
                if skipped:
                    msg += f" Skipped: {'; '.join(skipped)}."
                self.var_status.set(msg)
            else:
                self.var_status.set(f"Could not export expression: {'; '.join(skipped) or 'nothing to export'}.")

        self.background("Computing expression...", work, done)

    # ---- coordinates and the 3D views

    def selected_group_cells(self, level, targets):
        """(table, scope): one row per cell of the groups: the whole brain when the all-cell table is loaded, else only this run's cells."""
        bg = self.bg
        frames = []
        barr = bg["levels"].get(level) if bg else None
        use_bg = bg is not None and barr is not None and bg.get("sec_arr") is not None
        if use_bg:
            ids = core.extract_leading_numeric_id(pd.Series(barr)).to_numpy()
            for t in targets:
                rows = np.flatnonzero(ids == t)
                if not rows.size:
                    continue
                f = pd.DataFrame({"cell_label": bg["ids"][rows], "group_id": t, "brain_section_label": bg["sec_arr"][rows]})
                for tl in ("class", "subclass", "supertype", "cluster"):
                    v = bg["levels"].get(tl)
                    if v is not None:
                        f[tl] = v[rows]
                f["section_x_mm"] = bg["x"][rows]
                f["section_y_mm"] = bg["y"][rows]
                if bg.get("z") is not None:
                    f["section_z_mm"] = bg["z"][rows]
                frames.append(f)
        else:
            ids = self.level_ids(level)
            for t in targets:
                rows = np.flatnonzero(ids == t)
                if not rows.size:
                    continue
                f = pd.DataFrame({"cell_label": self.adata.obs_names[rows], "group_id": t})
                if self.section_of is not None:
                    f["brain_section_label"] = self.section_of[rows]
                for tl in ("class", "subclass", "supertype", "cluster"):
                    if self.info.get(tl) is not None:
                        f[tl] = self.info[tl][rows]
                if self.cx is not None:
                    f["section_x_mm"] = self.cx[rows]
                    f["section_y_mm"] = self.cy[rows]
                frames.append(f)
        table = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["cell_label", "group_id"])
        umap_rows = self.run_index.get_indexer(table["cell_label"])
        in_run = umap_rows >= 0
        table["in_run"] = in_run
        table["umap1"] = np.where(in_run, self.coords[np.maximum(umap_rows, 0), 0], np.nan)
        table["umap2"] = np.where(in_run, self.coords[np.maximum(umap_rows, 0), 1], np.nan)
        return table, ("whole brain" if use_bg else "this run only")

    def attach_ccf(self, table, shell_points=0):
        h = self.host
        if not h["ensure_ccf_downloaded"](self.abc_cache, parent=self.root):
            return None, False
        aligned, shell = h["run_with_loading_window"](
            "Reading the CCF coordinates (the first read of this 4-million-cell table takes about a minute).",
            lambda: h["load_ccf_for_cells"](self.abc_cache, table["cell_label"].tolist(), shell_points))
        for source, name in (("x_ccf", "ccf_x_mm"), ("y_ccf", "ccf_y_mm"), ("z_ccf", "ccf_z_mm")):
            table[name] = aligned[source].to_numpy()
        names = h["load_parcellation_names"](self.abc_cache)
        if names is not None:
            parcel = aligned["parcellation_index"].to_numpy()
            parcel = np.where(np.isfinite(parcel), parcel, -1).astype(int)
            for column in ("division", "structure", "substructure"):
                table[f"ccf_{column}"] = names[column].reindex(parcel).to_numpy()
        return shell, True

    def export_coordinates(self):
        if self.busy:
            return
        chosen = self.group_choice("exporting coordinates")
        if chosen is None:
            return
        level, targets = chosen
        stem = f"coordinates_{level}_{'-'.join(str(t) for t in targets)}"
        proposed = self.host["unique_export_paths"](self.target_dir(), stem, (".csv",))[".csv"].name
        path = self.ask_save_path("Save the coordinates of the cells of those groups", proposed, ".csv", "CSV file")
        if path is None:
            self.var_status.set("Coordinates export cancelled.")
            return
        self.root.config(cursor="watch")
        self.var_status.set("Gathering cells...")
        self.root.update()
        try:
            table, scope = self.selected_group_cells(level, targets)
            if table.empty:
                self.var_status.set("None of those IDs has any cells to export.")
            else:
                _shell, have_ccf = self.attach_ccf(table)
                front = ["cell_label", "group_id", "brain_section_label", "class", "subclass", "supertype", "cluster", "section_x_mm", "section_y_mm",
                         "section_z_mm", "ccf_x_mm", "ccf_y_mm", "ccf_z_mm", "ccf_division", "ccf_structure", "ccf_substructure", "in_run", "umap1", "umap2"]
                table = table[[c for c in front if c in table.columns]]
                table.to_csv(path, index=False)
                self.var_status.set(f"Saved {len(table):,} cells ({scope})" + ("" if have_ccf else " (no CCF columns: that download was declined)") + f" to {path}.")
        except Exception as e:
            self.var_status.set(f"Could not export coordinates: {e}")
            print(f"Could not export coordinates: {e}", flush=True)
        self.root.config(cursor="")

    def ask_3d_jitter(self):
        """Ask whether to jitter the section-axis coordinate of the cells (they sit on the planes of the sections, which looks layered): the
        spread in mm (0 = none), or None if cancelled. The last answer is the default."""
        import tkinter as tk
        from tkinter import ttk
        top = tk.Toplevel(self.root)
        top.title("3D in browser")
        top.transient(self.root)
        top.resizable(False, False)
        on = tk.BooleanVar(top, getattr(self, "jitter_3d_on", True))
        spread = tk.StringVar(top, getattr(self, "jitter_3d_mm", "0.15"))
        body = ttk.Frame(top, padding=12)
        body.pack()
        ttk.Checkbutton(body, text="Randomly jitter the section-axis (z) position of the cells", variable=on).grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Label(body, text="Spread, +/- mm:").grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(body, textvariable=spread, width=6).grid(row=1, column=1, sticky="w", padx=4, pady=(8, 0))
        ttk.Label(body, foreground="#555555", wraplength=300, justify="left",
                  text="The cells lie on the planes of the sections, which gives a layered look. Jitter spreads them out along that axis "
                       "(by default +/- 0.15 mm). Only the picture changes.").grid(row=2, column=0, columnspan=3, sticky="w", pady=(6, 0))
        result = []

        def ok():
            try:
                mm = max(0.0, float(spread.get())) if on.get() else 0.0
            except ValueError:
                spread.set("0.15")
                return
            self.jitter_3d_on, self.jitter_3d_mm = on.get(), spread.get()
            result.append(mm)
            top.destroy()
        buttons = ttk.Frame(body)
        buttons.grid(row=3, column=0, columnspan=3, sticky="e", pady=(10, 0))
        ttk.Button(buttons, text="Cancel", command=top.destroy).pack(side="right")
        ttk.Button(buttons, text="Open", command=ok).pack(side="right", padx=(0, 6))
        top.bind("<Return>", lambda e: ok())
        top.bind("<Escape>", lambda e: top.destroy())
        top.grab_set()
        self.root.wait_window(top)
        return result[0] if result else None

    def prepare_3d(self, action, jitter_mm=0.0):
        chosen = self.group_choice(action)
        if chosen is None:
            return None
        level, targets = chosen
        h = self.host
        table, scope = self.selected_group_cells(level, targets)
        if table.empty:
            self.var_status.set("None of those IDs has any cells to show.")
            return None
        shell, have_ccf = self.attach_ccf(table, shell_points=h["THREED_SHELL_POINTS"])
        if have_ccf:
            columns = ("ccf_z_mm", "ccf_x_mm", "ccf_y_mm")
            axis_labels = ("left - right (CCF z, mm)", "anterior - posterior (CCF x, mm)", "dorsal - ventral (CCF y, mm)")
            shell_xyz = shell[:, [2, 0, 1]]
            outline = h["brain_outline_segments"](shell_xyz)
        else:
            columns = ("section_x_mm", "section_z_mm", "section_y_mm")
            axis_labels = ("section x (mm)", "section position (mm)", "section y (mm)")
            shell_xyz, outline = None, None
        level_name = LEVEL_LONG[level]
        groups = []
        for t in targets:
            rows = table[table["group_id"] == t]
            hover = None
            if have_ccf and "ccf_structure" in rows.columns:
                hover = (rows["ccf_structure"].fillna("?").astype(str) + " | section " + rows["brain_section_label"].astype(str).str.rsplit(".", n=1).str[-1]
                         + " | " + rows["cell_label"].astype(str)).to_numpy()
            color = self.last_ids["colors"].get(t)
            xyz = rows[list(columns)].to_numpy(dtype=float)
            if jitter_mm > 0:                    # (the middle column is the one across the sections)
                xyz[:, 1] += np.random.default_rng().uniform(-jitter_mm, jitter_mm, len(xyz))
            groups.append({"label": f"{level_name} {t} ({len(rows):,})", "color": hexcolor(color) if color is not None else "red",
                           "xyz": xyz, "hover": hover})
        structures = table["ccf_structure"].to_numpy() if have_ccf and "ccf_structure" in table.columns else None
        return {"groups": groups, "shell_xyz": shell_xyz, "outline": outline, "axis_labels": axis_labels, "structures": structures,
                "title": f"3D view: {level_name.lower()} {', '.join(str(t) for t in targets)} ({scope})", "have_ccf": have_ccf, "level": level,
                "targets": targets}

    def show_3d(self):
        if self.busy:
            return
        self.root.config(cursor="watch")
        self.var_status.set("Preparing the 3D view...")
        self.root.update()
        try:
            p = self.prepare_3d("opening the 3D view")
            if p is not None:
                self.view3d = getattr(self, "view3d", [])
                self.view3d.append(self.host["open_3d_group_view"](p["groups"], p["shell_xyz"], p["outline"], p["axis_labels"], p["title"]))
                self.var_status.set("3D view opened: drag to rotate, click legend entries to hide/show a group."
                                    + ("" if p["have_ccf"] else " (CCF unavailable: section coordinates, so the sections look stacked.)"))
        except Exception as e:
            self.var_status.set(f"Could not open the 3D view: {e}")
            print(f"Could not open the 3D view: {e}", flush=True)
        self.root.config(cursor="")

    def show_3d_browser(self):
        if self.busy:
            return
        try:
            import plotly  # noqa: F401
        except ImportError:
            self.var_status.set("The browser 3D view needs plotly: run  python -m pip install plotly  and restart.")
            return
        if self.group_choice("opening the 3D view") is None:
            return
        jitter = self.ask_3d_jitter()
        if jitter is None:
            return
        self.root.config(cursor="watch")
        self.var_status.set("Preparing the 3D view...")
        self.root.update()
        try:
            p = self.prepare_3d("opening the 3D view", jitter)
            if p is not None:
                h = self.host
                target = self.target_dir()
                path = h["unique_export_paths"](target, f"3d_{p['level']}_{'-'.join(str(t) for t in p['targets'])}", (".html",))[".html"]
                edges = h["build_3d_edges"](self.abc_cache, p["structures"], self.root) if p["have_ccf"] else None
                self.var_status.set("Writing the 3D page...")
                self.root.update()
                h["write_3d_group_html"](p["groups"], p["shell_xyz"], p["outline"], p["axis_labels"], p["title"], path, edges=edges)
                h["open_with_default_viewer"](path)
                self.var_status.set(f"Opened {path.name} in your browser (saved in {target}).")
        except Exception as e:
            self.var_status.set(f"Could not open the 3D view in the browser: {e}")
            print(f"Could not open the 3D view in the browser: {e}", flush=True)
        self.root.config(cursor="")


LEVEL_LONG = {"class": "Class", "subclass": "Subclass", "supertype": "Supertype", "cluster": "Cluster", "leiden": "Leiden", "neurotransmitter": "Neurotransmitter"}
