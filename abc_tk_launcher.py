"""A Tk version of the atlas browser's launch screen: the output folder, then one of "Open Brain Browser / ROI Chooser", "Load
Existing ROI Dataset" or "Open 3D Viewer". Which cells to analyze (Neurons, Non-neurons or All cells) is asked later, by
choose_run_options(), once the sections and ROIs are confirmed. It is a function of the same shape as the matplotlib one (prompt_startup_panel in
ABC_atlas_browser.py): the browser passes in the few of its own functions it needs (recent folders, the check of a run folder, the 3D
viewer launcher, the settings window), so this file has no start-up side effects of its own.

    outcome = prompt_startup_panel_tk(out_folder, abc_cache, helpers)

`outcome` is (action, cell_type_selection, load_folder, out_folder) with action 'new_run' or 'load_existing', or None if the window was closed
or Exit was clicked. The cell type is None for a new run (not chosen yet) and unused for an existing one (the run says which it is).
"""
import os
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

CELL_TYPES = (("Neurons", "Neurons"), ("NonNeurons", "Non-neurons"), ("All", "All cells"))      # (value, the words on the screen)
ACCENT = "#1f5fa8"
RECOVERABLE = {"ROI coordinates", "Processed AnnData"}      # (cached files of a run that can be made again: not required to load it)


SUBSAMPLE_CAPS = (("200,000 cells (recommended)", 200_000), ("No limit: every cell", None))
LARGE_PROJECT_CELLS = 200_000            # a selection of more cells than this is a "large project": the dialog asks how to handle it


def choose_run_options(current="Neurons", counts=None, folder_for=None):
    """The one window between confirming the sections and ROIs and the long analysis that follows, so that the analysis then runs without
    asking anything more: which cells to analyze, how to handle a large project (a selection of more than LARGE_PROJECT_CELLS cells: use a
    random subset, or all), the project name that is added to the folder name, and whether to make the optional subclass / supertype plots.
    `counts`: {"Neurons": n, "NonNeurons": n, "All": n}, the cells in the selection, if known (they are shown, and the large-project question
    is asked only when the chosen kind has more than LARGE_PROJECT_CELLS; if not known it is always asked).
    Returns a dict {"cell_type": "Neurons" | "NonNeurons" | "All", "cap": the most cells to use, or None for no limit, "label": the project
    name (may be empty), "optional_plots": bool}, or None if cancelled.
    `folder_for(cell_type, project_name, cap)`, if given, returns the path of the project's folder for those choices (`cap`: the most cells
    to use, None for no limit): it is shown, read-only, under the project name, and follows the choices as they are made."""
    root = tk.Tk()
    root.title("Analysis options")
    root.resizable(False, False)
    cell_type = tk.StringVar(root, current)
    cap_choice = tk.StringVar(root, SUBSAMPLE_CAPS[0][0])
    custom_cap = tk.StringVar(root, "")
    label = tk.StringVar(root, "")
    optional = tk.BooleanVar(root, False)
    answer = {"value": None}
    f = ttk.Frame(root, padding=22)
    f.pack()

    ttk.Label(f, text="Which cells do you want to analyze?", font=("Segoe UI", 12, "bold")).pack(anchor="w")
    ttk.Label(f, text="The cells of this kind in the sections and ROIs you chose are made into the UMAP.", foreground="gray30",
              wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 4))
    for value, text in CELL_TYPES:
        shown = text if not counts or counts.get(value) is None else f"{text}  ({counts[value]:,} cells)"
        ttk.Radiobutton(f, text=shown, value=value, variable=cell_type).pack(anchor="w", pady=(3, 0))

    slot = ttk.Frame(f)                                           # (where the large-project question goes, when it is asked: see show_large)
    slot.pack(anchor="w")
    large = ttk.Frame(slot)
    ttk.Label(large, text="How to handle large projects", font=("Segoe UI", 12, "bold")).pack(anchor="w", pady=(16, 0))
    large_text = tk.StringVar(large, "")
    ttk.Label(large, textvariable=large_text, foreground="gray30", wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 4))
    for text, _cap in SUBSAMPLE_CAPS:
        ttk.Radiobutton(large, text=text, value=text, variable=cap_choice).pack(anchor="w", pady=(3, 0))
    custom_row = ttk.Frame(large)
    custom_row.pack(anchor="w", pady=(3, 0))
    ttk.Radiobutton(custom_row, text="At most", value="custom", variable=cap_choice).pack(side=tk.LEFT)
    custom_entry = ttk.Entry(custom_row, textvariable=custom_cap, width=10)
    custom_entry.pack(side=tk.LEFT, padx=(6, 6))
    ttk.Label(custom_row, text="cells (randomly chosen)").pack(side=tk.LEFT)
    custom_cap.trace_add("write", lambda *args: cap_choice.set("custom") if custom_cap.get() else None)         # (typing a number picks "At most")

    ttk.Label(f, text="Project name (optional)", font=("Segoe UI", 12, "bold")).pack(anchor="w", pady=(16, 0))
    ttk.Label(f, text="Added to the end of the project's folder name, to remember what it is, e.g. RMTg-VTA.", foreground="gray30",
              wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 4))
    ttk.Entry(f, textvariable=label, width=44).pack(anchor="w")
    if folder_for is not None:
        ttk.Label(f, text="The project will be saved in this folder:", foreground="gray30").pack(anchor="w", pady=(10, 2))
        folder_text = tk.StringVar(f, "")
        ttk.Entry(f, textvariable=folder_text, width=70, state="readonly").pack(anchor="w", fill=tk.X)

    ttk.Checkbutton(f, text="Also make optional subclass / supertype UMAP plots and spatial maps (may be very time consuming)",
                    variable=optional).pack(anchor="w", pady=(16, 0))

    row = ttk.Frame(f)
    row.pack(fill=tk.X, pady=(18, 0))

    def is_large():
        n = counts.get(cell_type.get()) if counts else None
        return n is None or n > LARGE_PROJECT_CELLS

    def show_large(*args):
        if is_large():
            n = counts.get(cell_type.get()) if counts else None
            large_text.set((f"This selection has {n:,} cells. " if n is not None else "") +
                           "A UMAP of that many cells takes a very long time. Use a random subset of this size, or all of them:")
            large.pack(anchor="w")
        else:
            large.pack_forget()

    def chosen_cap(complain=False):
        """The most cells the choices allow (None: no limit); LARGE_PROJECT_CELLS for a selection that is not large (it is not limited by it).
        With a custom number that is not a number: None if not complaining, else the complaint and 0."""
        if not is_large():
            return LARGE_PROJECT_CELLS
        cap = dict(SUBSAMPLE_CAPS).get(cap_choice.get(), "custom")
        if cap == "custom":
            try:
                cap = int(custom_cap.get().replace(",", "").strip())
                if cap < 1:
                    raise ValueError
            except ValueError:
                if complain:
                    messagebox.showerror("How to handle large projects", "Enter a whole number of cells, or pick one of the other choices.", parent=root)
                    return 0
                return "unknown"
        return cap

    def show_folder(*args):
        if folder_for is None:
            return
        cap = chosen_cap()
        folder_text.set(str(folder_for(cell_type.get(), label.get().strip(), LARGE_PROJECT_CELLS if cap == "unknown" else cap)))

    def done(ok):
        if not ok:
            answer["value"] = None
            root.destroy()
            return
        cap = chosen_cap(complain=True)
        if cap == 0:
            return
        answer["value"] = {"cell_type": cell_type.get(), "cap": cap, "label": label.get().strip(), "optional_plots": bool(optional.get())}
        root.destroy()

    ok_button = ttk.Button(row, text="Continue", command=lambda: done(True))
    ok_button.pack(side=tk.RIGHT)
    ttk.Button(row, text="Cancel", command=lambda: done(False)).pack(side=tk.RIGHT, padx=(0, 8))
    root.bind("<Return>", lambda e: done(True))
    root.bind("<Escape>", lambda e: done(False))
    root.protocol("WM_DELETE_WINDOW", lambda: done(False))
    cell_type.trace_add("write", show_large)
    for var in (cell_type, label, cap_choice, custom_cap):
        var.trace_add("write", show_folder)
    show_folder()
    root.update_idletasks()
    show_large()
    root.update_idletasks()
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry(f"+{(sw - root.winfo_width()) // 2}+{(sh - root.winfo_height()) // 4}")
    root.attributes("-topmost", True)
    ok_button.focus_force()
    root.mainloop()
    return answer["value"]


RUN_FILE = "umap_coords.csv"        # the one file a run needs (see required_cached_run_files in ABC_atlas_browser.py): a folder with it is a run


def scan_runs(folder):
    """The saved runs in `folder`: [(name, path, saved: the time umap_coords.csv was written)], the most recent first. A subfolder without
    umap_coords.csv is not a run (the folder of old outputs, say) and is left out. Slow on a network drive: call it from a thread."""
    runs = []
    try:
        entries = [e for e in os.scandir(folder) if e.is_dir()]
    except OSError:
        return runs
    for e in entries:
        try:
            runs.append((e.name, Path(e.path), os.stat(os.path.join(e.path, RUN_FILE)).st_mtime))
        except OSError:
            continue
    runs.sort(key=lambda r: r[2], reverse=True)
    return runs


class RunPicker:
    """A list of the saved runs in a folder, by name to start with (a click on the Run or the Saved heading sorts by that, and again the other way
    round): click one to highlight it (the most recent starts highlighted), then OK (or Enter, or a double-click) to take it; Cancel (or Escape)
    to go back. Browse... picks a run folder from anywhere else with the system's dialog. A box above the list
    narrows it to the names that hold what is typed. `result` is the chosen folder (a Path), or None."""

    def __init__(self, parent, folder):
        import queue
        import threading
        self.result = None
        self.folder = folder
        self.scanned = False                  # the folders have been looked into
        self.runs = []
        self.shown = []                       # the runs in the list now (the ones that match the filter)
        self.top = top = tk.Toplevel(parent)
        top.title("Load an existing run")
        top.geometry(f"{min(1300, top.winfo_screenwidth() - 80)}x{min(760, top.winfo_screenheight() - 120)}")
        top.minsize(640, 300)
        self.sort = ("name", False)           # (the heading sorted by, and whether it runs downward)
        body = ttk.Frame(top, padding=12)
        body.pack(fill=tk.BOTH, expand=True)
        ttk.Label(body, text=f"Saved runs in {folder}", font=("Segoe UI", 11)).pack(anchor="w")
        row = ttk.Frame(body)
        row.pack(fill=tk.X, pady=(8, 6))
        ttk.Label(row, text="Filter:").pack(side=tk.LEFT)
        self.var_filter = tk.StringVar(top, "")
        entry = ttk.Entry(row, textvariable=self.var_filter)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(6, 0))
        self.var_filter.trace_add("write", lambda *a: self.fill())
        frame = ttk.Frame(body)
        frame.pack(fill=tk.BOTH, expand=True)
        import tkinter.font as tkfont
        style = ttk.Style(top)
        style.configure("Runs.Treeview", rowheight=tkfont.nametofont("TkDefaultFont").metrics("linespace") + 12)      # (room around the text, at any screen scaling)
        self.tree = ttk.Treeview(frame, columns=("name", "saved"), show="headings", selectmode="browse", style="Runs.Treeview")
        self.tree.heading("name", command=lambda: self.sort_by("name"))
        self.tree.heading("saved", command=lambda: self.sort_by("saved"))
        self.mark_sort()
        self.tree.column("name", width=720, stretch=True)
        self.tree.column("saved", width=150, stretch=False)
        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tree.bind("<Double-Button-1>", lambda e: self.ok() if self.tree.identify_row(e.y) else None)
        self.var_status = tk.StringVar(top, "Looking for runs...")
        ttk.Label(body, textvariable=self.var_status, foreground="gray30").pack(anchor="w", pady=(6, 0))
        buttons = ttk.Frame(body)
        buttons.pack(fill=tk.X, pady=(10, 0))
        ttk.Button(buttons, text="Browse…", command=self.browse).pack(side=tk.LEFT)
        ttk.Button(buttons, text="Cancel", command=self.cancel).pack(side=tk.RIGHT)
        self.btn_ok = ttk.Button(buttons, text="OK", command=self.ok, state="disabled")
        self.btn_ok.pack(side=tk.RIGHT, padx=(0, 8))
        top.bind("<Return>", lambda e: self.ok())
        top.bind("<Escape>", lambda e: self.cancel())
        top.protocol("WM_DELETE_WINDOW", self.cancel)
        self.tree.bind("<<TreeviewSelect>>", lambda e: self.btn_ok.config(state="normal" if self.tree.selection() else "disabled"))
        top.transient(parent)
        entry.focus_set()
        # the folders are looked into on a thread (a network drive can take a while), and the list is filled when that is done
        self.found = queue.Queue()
        threading.Thread(target=lambda: self.found.put(scan_runs(folder)), daemon=True).start()
        self.poll()

    def poll(self):
        try:
            self.runs = self.found.get_nowait()
        except Exception:
            try:
                self.top.after(100, self.poll)
            except tk.TclError:
                pass
            return
        self.scanned = True
        self.fill()

    def mark_sort(self):
        """The headings: an arrow on the one the list is sorted by."""
        column, descending = self.sort
        for key, text in (("name", "Run"), ("saved", "Saved")):
            self.tree.heading(key, text=text + (("  ▼" if descending else "  ▲") if key == column else ""))

    def sort_by(self, column):
        """A click on a heading: sort by it (and the other way round if it already is)."""
        self.sort = (column, not self.sort[1]) if column == self.sort[0] else (column, column == "saved")      # (dates start with the newest)
        self.mark_sort()
        self.fill()

    def fill(self):
        """The runs that match the filter, in the sorted order; the run highlighted before stays so, else the most recent one is."""
        if not self.scanned:
            return
        keep = None
        if self.tree.selection() and self.shown:
            keep = self.shown[int(self.tree.selection()[0])][1]
        text = self.var_filter.get().strip().casefold()
        column, descending = self.sort
        self.shown = sorted([r for r in self.runs if text in r[0].casefold()],
                            key=(lambda r: r[0].casefold()) if column == "name" else (lambda r: r[2]), reverse=descending)
        self.tree.delete(*self.tree.get_children())
        import time
        for i, (name, _path, saved) in enumerate(self.shown):
            self.tree.insert("", tk.END, iid=str(i), values=(name, time.strftime("%Y-%m-%d %H:%M", time.localtime(saved))))
        if self.shown:
            paths = [r[1] for r in self.shown]
            at = paths.index(keep) if keep in paths else max(range(len(self.shown)), key=lambda i: self.shown[i][2])
            self.tree.selection_set(str(at))
            self.tree.focus(str(at))
            self.tree.see(str(at))
        self.btn_ok.config(state="normal" if self.shown else "disabled")
        if not self.runs:
            self.var_status.set("No saved runs in this folder. Use Browse… to pick a run folder from somewhere else.")
        else:
            self.var_status.set(f"{len(self.shown)} of {len(self.runs)} runs" if text else f"{len(self.runs)} runs, the most recent first.")

    def ok(self):
        picked = self.tree.selection()
        if picked and self.shown:
            self.result = self.shown[int(picked[0])][1]
            self.top.destroy()

    def cancel(self):
        self.result = None
        self.top.destroy()

    def browse(self):
        chosen = filedialog.askdirectory(parent=self.top, title="Select an existing run folder to load", initialdir=str(self.folder))
        if chosen:
            self.result = Path(chosen)
            self.top.destroy()

    def show(self):
        """Open it over its parent and wait for it to be closed; returns the chosen folder or None."""
        self.top.update_idletasks()
        parent = self.top.master
        x = parent.winfo_rootx() + max((parent.winfo_width() - self.top.winfo_width()) // 2, 0)
        y = parent.winfo_rooty() + 40
        self.top.geometry(f"+{x}+{y}")
        self.top.grab_set()
        self.top.wait_window()
        return self.result


def prompt_startup_panel_tk(out_folder, abc_cache, helpers):
    """`helpers`: load_recent_output_folders(), remember_recent_output_folder(folder), required_cached_run_files(folder) ->
    {label: Path}, open_3d_viewer(abc_cache), settings_dialog(abc_cache)."""
    root = tk.Tk()
    root.title("AllenABC UMAP Explorer")
    root.resizable(False, False)
    result = {"outcome": None}

    outer = ttk.Frame(root, padding=(28, 20, 28, 22))
    outer.pack()
    top = ttk.Frame(outer)
    top.pack(fill=tk.X)
    ttk.Label(top, text="AllenABC UMAP Explorer", font=("Segoe UI", 20, "bold")).pack(side=tk.LEFT)
    ttk.Button(top, text="Settings…", command=lambda: helpers["settings_dialog"](abc_cache)).pack(side=tk.RIGHT, anchor="n")

    ttk.Label(outer, text="Output folder (where runs are saved and loaded from):", font=("Segoe UI", 11)).pack(anchor="w", pady=(18, 4))
    row = ttk.Frame(outer)
    row.pack(fill=tk.X)
    folder_var = tk.StringVar(root, str(out_folder))

    def recent_items():
        typed = folder_var.get().strip().strip('"')
        items = list(helpers["load_recent_output_folders"]())
        norm = lambda p: os.path.normcase(os.path.abspath(os.path.expanduser(p)))
        if typed and norm(typed) not in {norm(p) for p in items}:
            items = [typed] + items
        return items

    folder_box = ttk.Combobox(row, textvariable=folder_var, width=70, postcommand=lambda: folder_box.config(values=recent_items()))
    folder_box.pack(side=tk.LEFT, fill=tk.X, expand=True)
    folder_box.config(values=recent_items())

    def browse():
        raw = folder_var.get().strip().strip('"')
        initial = raw if raw and Path(os.path.expanduser(raw)).is_dir() else str(out_folder)
        chosen = filedialog.askdirectory(parent=root, title="Choose the output folder", initialdir=initial)
        if chosen:
            folder_var.set(os.path.normpath(chosen))

    ttk.Button(row, text="Browse…", command=browse).pack(side=tk.LEFT, padx=(8, 0))

    def committed_folder():
        """The folder in the box, made if it doesn't exist; None (after saying why) if it can't be."""
        raw = folder_var.get().strip().strip('"')
        if not raw:
            messagebox.showerror("Output folder", "Enter an output folder, or click Browse… to pick one.", parent=root)
            return None
        path = Path(os.path.expanduser(raw))
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            messagebox.showerror("Output folder", f"Can't create or access that folder:\n\n{path}\n\n{e}", parent=root)
            return None
        return path

    def new_run():
        folder = committed_folder()
        if folder is None:
            return
        helpers["remember_recent_output_folder"](folder)
        result["outcome"] = ("new_run", None, None, folder)
        root.destroy()

    def load_existing():
        folder = committed_folder()
        if folder is None:
            return
        chosen_path = RunPicker(root, folder).show()
        if chosen_path is None:
            return
        required = helpers["required_cached_run_files"](chosen_path)
        missing = [f"{label}: {path}" for label, path in required.items() if label not in RECOVERABLE and not path.exists()]
        if missing:
            messagebox.showerror("Missing cached files", "This folder is missing required cached file(s):\n\n" + "\n".join(missing), parent=root)
            return                                          # (stay here, so another folder can be picked)
        helpers["remember_recent_output_folder"](folder)
        result["outcome"] = ("load_existing", "All", chosen_path, folder)
        root.destroy()

    def viewer3d():
        helpers["open_3d_viewer"](abc_cache)                # (returns when that window is closed; this one is still here)

    ttk.Separator(outer).pack(fill=tk.X, pady=(20, 14))
    style = ttk.Style(root)
    style.configure("Big.TButton", font=("Segoe UI", 12), padding=(10, 9))
    for text, command in (("Open Brain Browser/ROI Chooser", new_run), ("Load Existing ROI Dataset…", load_existing),
                          ("Open 3D Viewer…", viewer3d), ("Exit", root.destroy)):
        ttk.Button(outer, text=text, style="Big.TButton", command=command).pack(fill=tk.X, pady=4)

    root.protocol("WM_DELETE_WINDOW", root.destroy)
    root.update_idletasks()
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry(f"+{(sw - root.winfo_width()) // 2}+{(sh - root.winfo_height()) // 3}")
    root.mainloop()
    return result["outcome"]


if __name__ == "__main__":
    print(prompt_startup_panel_tk(
        Path.cwd(), None, {"load_recent_output_folders": lambda: [], "remember_recent_output_folder": lambda f: None,
                           "required_cached_run_files": lambda p: {}, "open_3d_viewer": lambda c: print("3D viewer"),
                           "settings_dialog": lambda c: print("settings")}))
