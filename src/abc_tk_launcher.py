"""A Tk version of the atlas browser's launch screen: the output folder, then one of "Choose ROI(s) for project", "Load
saved ROI project" or "Open 3D Viewer". Which cells to analyze (Neurons, Non-neurons or All cells) is asked later, by
choose_run_options(), once the sections and ROIs are confirmed. It is a function of the same shape as the matplotlib one (prompt_startup_panel in
ABC_atlas_browser.py): the browser passes in the few of its own functions it needs (recent folders, the check of a run folder, the 3D
viewer launcher, the settings window), so this file has no start-up side effects of its own.

    outcome = prompt_startup_panel_tk(out_folder, abc_cache, helpers)

`outcome` is (action, cell_type_selection, load_folder, out_folder) with action 'new_run', 'load_existing', 'own_dataset' or 'all_neurons' (make
the "All neurons" project: every neuron, every section, with no UMAP; an existing one comes back with its folder as load_folder), or None if
the window was closed or Exit was clicked. The cell type is None for a new run (not chosen yet) and unused for an existing one (the run says which it is).
"""
import contextlib
import gc
import os
import queue
import sys
import subprocess
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

CELL_TYPES = (("Neurons", "Neurons"), ("NonNeurons", "Non-neurons"), ("All", "All cells"))      # (value, the words on the screen)
ACCENT = "#1f5fa8"
RECOVERABLE = {"ROI coordinates", "Processed AnnData"}      # (cached files of a run that can be made again: not required to load it)


SUBSAMPLE_CAPS = (("200,000 cells (recommended)", 200_000), ("No limit: every cell", None))
ALL_NEURONS_LABEL = "All neurons"         # "Browse All Neurons": every neuron of every section, without a UMAP (a project of seconds), saved as ALL_NEURONS_FOLDER
ALL_NEURONS_FOLDER = "browse_neurons_allsections_All-neurons"        # (the name a project of Neurons / all sections / this label, without a UMAP, gets)
LARGE_PROJECT_CELLS = 200_000          # a selection of more cells than this is a "large project": the dialog asks how to handle it


def bring_to_front(root):
    """Puts a window that is just opening in front of the others, with the keyboard focus. Windows won't let a program take the foreground
    unless the user's last click or key was in it, and the browser takes a while to start (time enough to click somewhere else), so
    without this its first window can open behind the others, with only its taskbar button flashing. Being "always on top" is allowed
    anyway: on for a moment brings the window in front, then off again (else it would stay above everything, even its own file dialogs)."""
    def front():
        root.attributes("-topmost", True)
        root.lift()
        root.focus_force()
        root.after(400, lambda: root.winfo_exists() and root.attributes("-topmost", False))
    root.after(50, front)                  # (once the window is on screen)


def choose_run_options(current="Neurons", counts=None, folder_for=None, whole_brain=False, count_state=None):
    """The one window between confirming the sections and ROIs and the long analysis that follows, so that the analysis then runs without
    asking anything more: which cells to analyze, whether to make a UMAP at all (not making one, the default for a large project, takes seconds,
    whatever the cells), how to handle a large project that is to have one (a selection of more than LARGE_PROJECT_CELLS cells: use a random subset, or
    all), the project name that is added to the folder name, and whether to make the optional subclass / supertype plots.
    `counts`: {"Neurons": n, "NonNeurons": n, "All": n}, the cells in the selection, if known (they are shown, and the large-project question
    is asked only when the chosen kind has more than LARGE_PROJECT_CELLS; if not known it is always asked).
    Returns a dict {"cell_type": "Neurons" | "NonNeurons" | "All", "make_umap": bool, "cap": the most cells to use in the UMAP, or None for no limit
    (always None without a UMAP), "label": the project name (may be empty), "optional_plots": bool}, or None if cancelled.
    `count_state`: {"done": bool, "counts": {...} or None, "error": exception or None}, filled in by a thread that is counting the cells now: the
    window opens at once, shows "counting..." for the numbers, keeps Continue off, and fills the numbers in (and settles the large-project question)
    when the counting is done. The caller's `folder_for` may read the same dict.
    `whole_brain`: no sections or ROIs were picked, so every section is used: the window says so (it is the question that used to be asked at the console).
    `folder_for(cell_type, project_name, cap, make_umap)`, if given, returns the path of the project's folder for those choices (`cap`: the most
    cells to use, None for no limit): it is shown, read-only, under the project name, and follows the choices as they are made."""
    root = tk.Tk()
    root.title("Analysis options")
    root.resizable(False, False)
    cell_type = tk.StringVar(root, current)
    cap_choice = tk.StringVar(root, SUBSAMPLE_CAPS[0][0])
    custom_cap = tk.StringVar(root, "")
    label = tk.StringVar(root, "")
    optional = tk.BooleanVar(root, False)
    skip_umap = tk.BooleanVar(root, False)                        # (set from the number of cells when it is known: see sync_skip)
    skip_touched = {"yes": False}                                 # (the box was ticked or unticked by hand: it is no longer set from the number of cells)
    answer = {"value": None}
    f = ttk.Frame(root, padding=22)
    f.pack()

    ttk.Label(f, text="Which cells do you want to analyze?", font=("Segoe UI", 12, "bold")).pack(anchor="w")
    ttk.Label(f, text="The cells of this kind in the sections and ROIs you chose are shown in the section maps, and made into a UMAP unless you skip it.",
              foreground="gray30", wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 4))
    if whole_brain:                                               # (nothing was picked: this stands in for a question, which Cancel answers with no)
        ttk.Label(f, text="Whole brain: no sections or ROIs were picked, so every section is used. Cancel returns to the start screen.",
                  foreground="#a05000", font=("Segoe UI", 10, "bold"), wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 6))
    known = {"counts": counts, "pending": count_state is not None and not count_state["done"]}          # (what is known about the selection's cells now)
    kind_buttons = {}

    def kind_text(value, text):
        if known["pending"]:
            return f"{text}  (counting...)"
        c = known["counts"]
        return text if not c or c.get(value) is None else f"{text}  ({c[value]:,} cells)"

    for value, text in CELL_TYPES:
        kind_buttons[value] = ttk.Radiobutton(f, text=kind_text(value, text), value=value, variable=cell_type)
        kind_buttons[value].pack(anchor="w", pady=(3, 0))

    ttk.Label(f, text="UMAP", font=("Segoe UI", 12, "bold")).pack(anchor="w", pady=(16, 0))
    skip_box = ttk.Checkbutton(f, text="Skip the UMAP (the default for a large project)", variable=skip_umap,
                               command=lambda: (skip_touched.update(yes=True), show_large(), show_folder()))
    skip_box.pack(anchor="w", pady=(3, 0))
    ttk.Label(f, text="Without a UMAP the project opens in seconds, however many cells it has, and shows every cell of the selection in the section "
                      "maps. Gene expression is read from the atlas as it is wanted. There is then no UMAP to show, and no Leiden clusters; a UMAP "
                      "cannot be added to the project later.", foreground="gray30", wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 0))

    slot = ttk.Frame(f)                                           # (where the large-project question goes: always there, greyed until it applies: see show_large)
    slot.pack(anchor="w")
    large = ttk.Frame(slot)
    large.pack(anchor="w")
    large_widgets = []
    ttk.Label(large, text="How to handle large projects (with a UMAP)", font=("Segoe UI", 12, "bold")).pack(anchor="w", pady=(16, 0))
    large_text = tk.StringVar(large, "")
    ttk.Label(large, textvariable=large_text, foreground="gray30", wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 4))
    first_text, last_text = SUBSAMPLE_CAPS[0][0], SUBSAMPLE_CAPS[-1][0]            # (the recommended number, then a number of your own, then no limit)
    b = ttk.Radiobutton(large, text=first_text, value=first_text, variable=cap_choice)
    b.pack(anchor="w", pady=(3, 0))
    large_widgets.append(b)
    custom_row = ttk.Frame(large)
    custom_row.pack(anchor="w", pady=(3, 0))
    b = ttk.Radiobutton(custom_row, text="At most", value="custom", variable=cap_choice)
    b.pack(side=tk.LEFT)
    large_widgets.append(b)
    custom_entry = ttk.Entry(custom_row, textvariable=custom_cap, width=10)
    custom_entry.pack(side=tk.LEFT, padx=(6, 6))
    large_widgets.append(custom_entry)
    ttk.Label(custom_row, text="cells (randomly chosen)").pack(side=tk.LEFT)
    custom_cap.trace_add("write", lambda *args: cap_choice.set("custom") if custom_cap.get() else None)         # (typing a number picks "At most")
    b = ttk.Radiobutton(large, text=last_text, value=last_text, variable=cap_choice)
    b.pack(anchor="w", pady=(3, 0))
    large_widgets.append(b)

    ttk.Label(f, text="Project name (optional)", font=("Segoe UI", 12, "bold")).pack(anchor="w", pady=(16, 0))
    ttk.Label(f, text="Added to the end of the project's folder name, to remember what it is, e.g. RMTg-VTA.", foreground="gray30",
              wraplength=400, justify=tk.LEFT).pack(anchor="w", pady=(2, 4))
    ttk.Entry(f, textvariable=label, width=44).pack(anchor="w")
    if folder_for is not None:
        ttk.Label(f, text="The project will be saved in this folder:", foreground="gray30").pack(anchor="w", pady=(10, 2))
        folder_text = tk.StringVar(f, "")
        ttk.Entry(f, textvariable=folder_text, width=70, state="readonly").pack(anchor="w", fill=tk.X)

    optional_box = ttk.Checkbutton(f, text="Also make optional subclass / supertype UMAP plots and spatial maps (may be very time consuming)",
                                   variable=optional)
    optional_box.pack(anchor="w", pady=(16, 0))

    row = ttk.Frame(f)
    row.pack(fill=tk.X, pady=(18, 0))

    def count_of_kind():
        c = known["counts"]
        return c.get(cell_type.get()) if c else None

    def is_large():
        """Whether the large-project question applies: when the number of cells is not known (the count failed), it is asked regardless. While the
        count is still going, it does not apply yet (Continue is off until it is done)."""
        if known["pending"]:
            return False
        n = count_of_kind()
        return n is None or n > LARGE_PROJECT_CELLS

    def sync_skip(*args):
        """Until the box has been ticked or unticked by hand, it follows the selection: ticked for a large one (the default), not for one a UMAP is
        quick for."""
        if not skip_touched["yes"] and not known["pending"]:
            skip_umap.set(is_large())

    def show_large(*args):
        """The large-project block is always there (so that nothing moves when the count arrives): greyed while counting, when the selection is
        small, and while the UMAP is skipped; active when it is large and a UMAP is made."""
        sync_skip()
        optional_box.state(["disabled"] if skip_umap.get() else ["!disabled"])
        if skip_umap.get():
            large_text.set("No limit is needed without a UMAP: every cell of the selection is shown.")
            for w in large_widgets:
                w.state(["disabled"])
            return
        if known["pending"]:
            large_text.set("Counting the cells in your selection...")
        elif is_large():
            n = count_of_kind()
            large_text.set((f"This selection has {n:,} cells. " if n is not None else "") +
                           "A UMAP of that many cells takes a very long time. Use a pseudorandom subset of this size (the same cells each time), or all of them:")
        else:
            large_text.set(f"This selection has {count_of_kind():,} cells, so there is no need to limit it (that is asked for more than "
                           f"{LARGE_PROJECT_CELLS:,}).")
        for w in large_widgets:
            w.state(["!disabled"] if is_large() else ["disabled"])

    def chosen_cap(complain=False):
        """The most cells the choices allow (None: no limit); LARGE_PROJECT_CELLS for a selection that is not large (it is not limited by it).
        With a custom number that is not a number: None if not complaining, else the complaint and 0. Without a UMAP there is no limit: None."""
        if skip_umap.get():
            return None
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
        folder_text.set(str(folder_for(cell_type.get(), label.get().strip(), LARGE_PROJECT_CELLS if cap == "unknown" else cap, not skip_umap.get())))

    def done(ok):
        if not ok:
            answer["value"] = None
            root.destroy()
            return
        if known["pending"]:                                          # (the count is not done: Continue is off, and so is Enter)
            return
        cap = chosen_cap(complain=True)
        if cap == 0:
            return
        answer["value"] = {"cell_type": cell_type.get(), "make_umap": not skip_umap.get(), "cap": cap, "label": label.get().strip(),
                           "optional_plots": bool(optional.get()) and not skip_umap.get()}
        root.destroy()

    ok_button = ttk.Button(row, text="Continue", command=lambda: done(True))
    ok_button.pack(side=tk.RIGHT)
    if known["pending"]:
        ok_button.state(["disabled"])

    def counted():
        """The counting is done (or failed): the numbers, the large-project block and the folder follow, and Continue is on."""
        known["pending"] = False
        known["counts"] = None if count_state["error"] is not None else count_state["counts"]
        for value, text in CELL_TYPES:
            kind_buttons[value].config(text=kind_text(value, text))
        show_large()
        show_folder()
        ok_button.state(["!disabled"])

    def watch_count():
        if not count_state["done"]:
            root.after(100, watch_count)
            return
        try:
            counted()
        except tk.TclError:                                           # (the window was closed meanwhile)
            pass

    ttk.Button(row, text="Cancel", command=lambda: done(False)).pack(side=tk.RIGHT, padx=(0, 8))
    root.bind("<Return>", lambda e: done(True))
    root.bind("<Escape>", lambda e: done(False))
    root.protocol("WM_DELETE_WINDOW", lambda: done(False))
    cell_type.trace_add("write", show_large)
    for var in (cell_type, label, cap_choice, custom_cap, skip_umap):
        var.trace_add("write", show_folder)
    show_folder()
    root.update_idletasks()
    show_large()
    root.update_idletasks()
    if known["pending"]:
        watch_count()
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry(f"+{(sw - root.winfo_width()) // 2}+{(sh - root.winfo_height()) // 4}")
    root.attributes("-topmost", True)
    ok_button.focus_force()
    root.mainloop()
    return answer["value"]


SUBSET_PREFIX = "subset_"          # the subfolders Analyze subset makes in a run (abc_tk_subset.py), each a run of its own
MAX_SUBSET_DEPTH = 8               # how many levels of subsets inside subsets the list follows
ROI_FILE = "roi_coords.csv"          # the ROIs a run was made from: a folder with this but no RUN_FILE is a run that stopped before its UMAP was made
RUN_FILE = "umap_coords.csv"        # the one file a run needs (see required_cached_run_files in ABC_atlas_browser.py): a folder with it is a run


def scan_runs(folder):
    """The saved runs in `folder`: [(name, path, opened: when it was last opened, or else saved)], the most recent first. A run counts as opened
    when it was made (the time umap_coords.csv was written) and each time it is opened after that (see abc_core.record_run_opened). A subfolder
    without umap_coords.csv is not a run (the folder of old outputs, say) and is left out. Slow on a network drive: call it from a thread."""
    import abc_core
    opened = abc_core.run_opened_times()
    runs = []
    try:
        entries = [e for e in os.scandir(folder) if e.is_dir()]
    except OSError:
        return runs
    import abc_project
    for e in entries:
        try:
            saved = os.stat(os.path.join(e.path, RUN_FILE)).st_mtime
        except OSError:
            # A project made without a UMAP: its project file says so, and it is a finished project.
            if not abc_project.has_umap(e.path):
                try:
                    saved = os.stat(os.path.join(e.path, abc_project.PROJECT_FILE)).st_mtime
                    runs.append((e.name, Path(e.path), max(saved, opened.get(abc_core._run_key(e.path), 0.0)), 0, None, True))
                except OSError:
                    pass
                continue
            # No UMAP: a run that stopped before it was made. If it saved its ROIs it is listed, as unfinished, to be finished from them.
            try:
                if e.name.startswith("umap_"):
                    runs.append((e.name, Path(e.path), os.stat(os.path.join(e.path, ROI_FILE)).st_mtime, 0, None, False))
            except OSError:
                pass
            continue
        runs.append((e.name, Path(e.path), max(saved, opened.get(abc_core._run_key(e.path), 0.0)), 0, None, True))
        _add_subsets(e.path, 1, runs, opened, abc_core)
    runs.sort(key=lambda r: r[2], reverse=True)
    return runs


def _add_subsets(path, depth, runs, opened, abc_core):
    """The subsets (Analyze subset in the viewer) kept in subfolders of the run at `path`, and theirs in turn: [(name, folder, opened, depth, the folder
    of the project it was made from)] added to `runs`."""
    if depth > MAX_SUBSET_DEPTH:
        return
    try:
        entries = [e for e in os.scandir(path) if e.is_dir() and e.name.startswith(SUBSET_PREFIX)]
    except OSError:
        return
    for e in entries:
        try:
            saved = os.stat(os.path.join(e.path, RUN_FILE)).st_mtime
        except OSError:
            continue
        runs.append((e.name, Path(e.path), max(saved, opened.get(abc_core._run_key(e.path), 0.0)), depth, Path(path), True))
        _add_subsets(e.path, depth + 1, runs, opened, abc_core)


class RunPicker:
    """A list of the saved runs in a folder, by name to start with (a click on the Run or the Saved heading sorts by that, and again the other way
    round): click one to highlight it (the one opened most recently starts highlighted), then OK (or Enter, or a double-click) to take it; Cancel (or Escape)
    to go back. Browse... picks a run folder from anywhere else with the system's dialog. A box above the list
    narrows it to the names that hold what is typed. `result` is the chosen folder (a Path), or None."""

    def __init__(self, parent, folder, helpers=None):
        import queue
        import threading
        self.result = None
        self.folder = folder
        self.helpers = helpers or {}
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
        self.tree.column("saved", width=170, stretch=False)
        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tree.tag_configure("incomplete", foreground="#a05000")
        self.tree.bind("<Double-Button-1>", self.double_clicked)
        self.menu = tk.Menu(top, tearoff=False)
        self.menu.add_command(label="Delete project folder and files…", command=self.delete_project)
        for press in (("<Button-3>", "<Button-2>", "<Control-Button-1>") if top.tk.call("tk", "windowingsystem") == "aqua" else ("<Button-3>",)):
            self.tree.bind(press, self.right_clicked)
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
        for key, text in (("name", "Run"), ("saved", "Last opened")):
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
        # (the projects in the sorted order, each followed by its subsets, sorted the same way: a project that matches the filter brings its subsets,
        # and a subset that matches brings the projects it is under)
        children, parent_of = {}, {}
        for r in self.runs:
            children.setdefault(r[4], []).append(r)
            parent_of[r[1]] = r[4]
        show = set()
        for r in self.runs:
            if text in r[0].casefold():
                stack = [r]
                while stack:
                    item = stack.pop()
                    show.add(item[1])
                    stack.extend(children.get(item[1], []))
                up = parent_of.get(r[1])
                while up is not None:
                    show.add(up)
                    up = parent_of.get(up)
        key = (lambda r: r[0].casefold()) if column == "name" else (lambda r: r[2])
        self.shown = []

        def emit(parent):
            for r in sorted(children.get(parent, []), key=key, reverse=descending):
                if r[1] in show:
                    self.shown.append(r)
                    emit(r[1])
        emit(None)
        self.tree.delete(*self.tree.get_children())
        import time
        for i, (name, _path, saved, depth, _parent, complete) in enumerate(self.shown):
            shown_name = "      " * depth + "\u21b3 " + name if depth else name           # (a subset: under its project, indented)
            if not complete:
                shown_name += "     (incomplete: no UMAP yet; OK finishes it from its saved ROIs)"
            self.tree.insert("", tk.END, iid=str(i), values=(shown_name, time.strftime("%Y-%m-%d %H:%M", time.localtime(saved))),
                             tags=() if complete else ("incomplete",))
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
            self.var_status.set(f"{len(self.shown)} of {len(self.runs)} projects" if text else f"{len(self.runs)} projects, subsets under the project they were made from.")

    def double_clicked(self, event):
        if not self.tree.identify_row(event.y):
            return None
        self.ok()
        return "break"                  # (the window is gone: the tree's own double-click handler, which runs after this one, would fail on it)

    def right_clicked(self, event):
        """A right click on a project: it is highlighted, and the menu with Delete is offered."""
        row = self.tree.identify_row(event.y)
        if not row:
            return None
        self.tree.selection_set(row)
        self.tree.focus(row)
        try:
            self.menu.tk_popup(event.x_root, event.y_root)
        finally:
            self.menu.grab_release()
        return "break"

    def delete_project(self):
        """Delete the highlighted project's folder (with the subsets kept inside it) and the processed data file saved for it in the cache, after a
        second question that says exactly what goes."""
        picked = self.tree.selection()
        if not picked or not self.shown:
            return
        name, path = self.shown[int(picked[0])][:2]
        path = Path(path)
        inside = [r for r in self.runs if r[1] != path and path in Path(r[1]).parents]
        files = []                                              # (the processed data saved in the cache for the project and for each subset)
        find = self.helpers.get("required_cached_run_files")
        if find is not None:
            for r in [self.shown[int(picked[0])]] + inside:
                try:
                    f = find(r[1]).get("Processed AnnData")
                except Exception:
                    f = None
                if f is not None and Path(f).exists():
                    files.append(Path(f))
        text = f"Permanently delete the project \"{name}\"?\n\nThe folder, with everything in it:\n    {path}\n"
        if inside:
            text += f"\nincluding its {len(inside)} subset{'s' if len(inside) != 1 else ''}:\n" + "\n".join("    " + r[0] for r in inside[:8]) + (
                f"\n    … and {len(inside) - 8} more" if len(inside) > 8 else "") + "\n"
        if files:
            text += "\nand the saved processed data in the cache:\n" + "\n".join("    " + f.name for f in files[:8]) + (
                f"\n    … and {len(files) - 8} more" if len(files) > 8 else "") + "\n"
        text += "\nThis cannot be undone."
        with question():
            yes = messagebox.askyesno("Delete project folder and files", text, icon="warning", default="no", parent=self.top)
        if not yes:
            return
        import shutil
        self.top.config(cursor="watch")
        self.var_status.set(f"Deleting {name}...")
        self.top.update_idletasks()
        problems = []
        try:
            shutil.rmtree(path)
        except OSError as e:
            problems.append(f"{path}: {e}")
        for f in files:
            try:
                f.unlink()
            except OSError as e:
                problems.append(f"{f}: {e}")
        self.top.config(cursor="")
        if problems:
            messagebox.showerror("Could not delete everything", "These could not be deleted:\n\n" + "\n".join(problems), parent=self.top)
        gone = {r[1] for r in self.runs if (r[1] == path or path in Path(r[1]).parents) and not Path(r[1]).exists()}
        self.runs = [r for r in self.runs if r[1] not in gone]
        self.fill()

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


QUESTIONS = {"open": 0}          # the launcher's own questions on the screen now: the update question waits for them (two at once would pile up)


@contextlib.contextmanager
def question():
    """Around a messagebox question: while it is up, the automatic update waits for it (see start_auto_update)."""
    QUESTIONS["open"] += 1
    try:
        yield
    finally:
        QUESTIONS["open"] -= 1


def start_auto_update(root, button, helpers):
    """With "Download updates automatically" on (the Settings window's Common tab), the updates the check at launch found are downloaded without a
    question, and the browser restarts (see app_update.AutoUpdate, which the other apps use too). It waits while a dialog has the input or a question
    is up, does nothing while the button is busy (the user's own check), and does not restart while the atlas is downloading: that would stop the
    download. The button then says "Updated: restart to use it"."""
    import app_update

    def blocked():
        if QUESTIONS["open"]:
            return True
        try:
            top = root.grab_current()
            return top is not None and top is not root and top.winfo_exists()
        except (tk.TclError, KeyError):
            return False

    def status(text):
        button.config(text=text or getattr(button, "idle_text", "Check for updates"), state=tk.DISABLED if text == "Updating..." else tk.NORMAL)
        if text and text != "Updating...":
            button.idle_text = text

    def restart_now():
        root.destroy()
        restart()

    app_update.AutoUpdate(root, {
        "busy": lambda: str(button.cget("state")) == "disabled",
        "blocked": blocked,
        "status": status,
        "keep_open": lambda: ("the atlas is still downloading (" + helpers["atlas_download_running"]() + "), and restarting would stop it")
                      if (helpers or {}).get("atlas_download_running", lambda: None)() else None,
        "restart": restart_now,
    }).start()


def check_for_update(root, button, changes=None, helpers=None):
    """The launch screen's Update… button: looks on GitHub for newer commits (in a thread: the network), says what is new, and on a yes
    brings the apps' folder up to date (git pull for a clone, else GitHub's ZIP: see app_update.py) and offers to restart the browser.
    `changes`: those the check at launch already found: the question is asked at once, without looking again. `helpers`: as for
    prompt_startup_panel_tk (the restart question says if the atlas download is running)."""
    import app_update
    done = queue.Queue()

    def work(job, *args):
        try:
            done.put((True, job(*args)))
        except app_update.UpdateError as e:
            done.put((False, str(e)))
        except Exception as e:                      # (anything unforeseen: still reported, not a button stuck on "Checking…")
            done.put((False, f"{type(e).__name__}: {e}"))

    def run(job, *args, then, label):
        button.config(text=label, state=tk.DISABLED)
        root.config(cursor="watch")
        threading.Thread(target=work, args=(job,) + args, daemon=True).start()

        def poll():
            try:
                ok, value = done.get_nowait()
            except queue.Empty:
                root.after(100, poll)
                return
            then(ok, value)
        poll()

    def finish(text=None, error=False):
        button.config(text=getattr(button, "idle_text", "Check for updates"), state=tk.NORMAL)
        root.config(cursor="")
        if text:
            with question():
                (messagebox.showerror if error else messagebox.showinfo)("Update", text, parent=dialog_parent())

    def dialog_parent():
        """The window a question goes over: the dialog that has the input (e.g. the project list), else the launch screen; a question over
        the launch screen would be hidden behind that dialog, and its answer blocked."""
        try:
            top = root.grab_current()
            if top is not None and top.winfo_exists():
                return top
        except (tk.TclError, KeyError):
            pass
        return root

    def checked(ok, value):
        if not ok:
            return finish(value, True)
        if not value:
            button.idle_text = "Check for updates"
            return finish("You have the newest version.")
        shown = "\n".join("  • " + line for line in value[:12]) + (f"\n  ... and {len(value) - 12} more" if len(value) > 12 else "")
        news = f"{len(value)} new change{'s' if len(value) != 1 else ''} on GitHub:\n\n{shown}\n\n"
        reason = app_update.no_pull_reason(folder)
        if reason:                                  # (a git clone, and no git: what is new, and how to get it)
            return finish(news + reason)
        root.config(cursor="")
        with question():
            yes = messagebox.askyesno("Update", news + "Download them now?", parent=dialog_parent())
        if not yes:
            return finish()
        run(app_update.pull, folder, then=pulled, label="Updating…")

    def pulled(ok, value):
        if not ok:
            return finish(value, True)
        button.idle_text = "Check for updates"
        finish()
        downloading = (helpers or {}).get("atlas_download_running", lambda: None)()
        if downloading:
            # (restarting would pause the download of the atlas until the next start: so no restart is offered while it is going)
            with question():
                messagebox.showinfo("Update", "Updated. The new version runs the next time the browser starts.\n\n"
                                    f"The atlas is still downloading ({downloading}), so the browser is not restarted now: that would stop the "
                                    "download. Close the browser after the download has finished.", parent=dialog_parent())
            return
        with question():
            yes = messagebox.askyesno("Update", "Updated. The new version runs after the browser restarts.\n\nRestart it now?", parent=dialog_parent())
        if yes:
            top = dialog_parent()
            if top is not root:
                top.destroy()                               # (a dialog waiting for its window to close: it ends first)

            def go():
                root.destroy()
                restart()
            root.after(100, go)

    try:
        folder = app_update.repo_root()
    except app_update.UpdateError as e:
        return finish(str(e), True)
    if changes:
        button.config(state=tk.DISABLED)
        checked(True, changes)
        return
    run(app_update.check, folder, then=checked, label="Checking…")


def restart():
    """Starts the browser again in the same console window (app_update.relaunch() would start it with none, hiding its progress output),
    and ends this one."""
    args = [sys.executable, os.path.abspath(sys.argv[0])] + sys.argv[1:]
    os.environ["ABC_RESTARTED"] = "1"               # (the new copy waits a while for this one to end, and so for its hold on instance_lock to be let go)
    sys.stdout.flush()
    sys.stderr.flush()
    if os.name == "nt":
        subprocess.Popen(args, close_fds=True)
        os._exit(0)                                  # (not sys.exit: background loading threads would hold it open)
    os.execv(sys.executable, args)                  # (macOS / Linux: the same process, so the same Terminal window)


def prompt_startup_panel_tk(out_folder, abc_cache, helpers):
    """`helpers`: load_recent_output_folders(), remember_recent_output_folder(folder), required_cached_run_files(folder) ->
    {label: Path}, open_3d_viewer(abc_cache), settings_dialog(abc_cache), and (optional) atlas_download_running() -> how far the atlas
    download has got, as text, while it runs (else None)."""
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
    if not getattr(sys, "frozen", False):                   # (the packaged .exe is not a folder from GitHub: nothing to update)
        update_button = ttk.Button(top, text="Check for updates", command=lambda: check_for_update(root, update_button, helpers=helpers))
        update_button.pack(side=tk.RIGHT, anchor="n", padx=(0, 6))

        def update_found(changes):
            """The check at launch found newer commits: the Check for updates button becomes the Update available one (its dialog lists them and asks).
            With "Download updates automatically" on they are downloaded at once instead, without a question, and the browser restarts; with it off
            nothing else is shown."""
            update_button.idle_text = f"Update available ({len(changes)} change{'s' if len(changes) != 1 else ''})"
            if str(update_button.cget("state")) != "disabled":
                update_button.config(text=update_button.idle_text)
            if app_update.auto_update_enabled():
                start_auto_update(root, update_button, helpers)
        import app_update
        app_update.watch_for_update(root, update_found)

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
        chosen_path = RunPicker(root, folder, helpers).show()
        if chosen_path is None:
            return
        chosen_path = Path(chosen_path)
        import abc_project
        if not (chosen_path / RUN_FILE).exists() and (chosen_path / ROI_FILE).exists() and abc_project.has_umap(chosen_path):
            # a run that stopped before its UMAP was made: finished from the ROIs it saved (the cells, the project name and the size limit are in its name)
            with question():
                finish_it = messagebox.askyesno(
                    "Finish this project?",
                    f"{chosen_path.name}\n\nThis project was never finished: its ROIs were saved, but the UMAP was not made.\n\n"
                    "Make it now, from the saved ROIs? It takes as long as a new project (the cells, the project name and the size limit "
                    "are taken from the folder's name).", parent=root)
            if not finish_it:
                return
            helpers["remember_recent_output_folder"](folder)
            result["outcome"] = ("resume", None, chosen_path, folder)
            root.destroy()
            return
        required = helpers["required_cached_run_files"](chosen_path)
        missing = [f"{label}: {path}" for label, path in required.items() if label not in RECOVERABLE and not path.exists()]
        if missing:
            messagebox.showerror("Missing cached files", "This folder is missing required cached file(s):\n\n" + "\n".join(missing), parent=root)
            return                                          # (stay here, so another folder can be picked)
        helpers["remember_recent_output_folder"](folder)
        result["outcome"] = ("load_existing", "All", chosen_path, folder)
        root.destroy()

    def own_dataset():
        folder = committed_folder()
        if folder is None:
            return
        chosen = filedialog.askopenfilename(parent=root, title="Open a dataset of your own (.h5ad)", filetypes=[("AnnData files", "*.h5ad"), ("All files", "*.*")])
        if not chosen:
            return
        result["outcome"] = ("own_dataset", None, Path(chosen), folder)
        root.destroy()

    def viewer3d():
        helpers["open_3d_viewer"](abc_cache)                # (returns when that window is closed; this one is still here)

    def existing_all_neurons(folder):
        """The "All neurons" project in `folder`, made on an earlier day, or None. (A UMAP of 500,000 of the neurons, which this button used to make, is
        a project like another: it is opened with Load Existing, and this button now makes the project of all of them.)"""
        import abc_project
        return folder / ALL_NEURONS_FOLDER if abc_project.read_project(folder / ALL_NEURONS_FOLDER) is not None else None

    def make_all_neurons(folder):
        helpers["remember_recent_output_folder"](folder)
        result["outcome"] = ("all_neurons", None, None, folder)
        root.destroy()

    def all_neurons():
        folder = committed_folder()
        if folder is None:
            return
        # made on an earlier day: simply opened, like any saved run
        existing = existing_all_neurons(folder)
        if existing is not None:
            required = helpers["required_cached_run_files"](existing)
            missing = [f"{label}: {path}" for label, path in required.items() if label not in RECOVERABLE and not path.exists()]
            if missing:
                messagebox.showerror("Missing cached files", "The saved \"All neurons\" project is missing required cached file(s):\n\n" + "\n".join(missing), parent=root)
                return
            helpers["remember_recent_output_folder"](folder)
            result["outcome"] = ("all_neurons", None, existing, folder)           # (with the folder: opened, not made)
            root.destroy()
            return
        # not made yet: made now, without asking (a few seconds: see create_no_umap_project; the atlas still downloading, the first time, is waited for in
        # a window of its own, which says so)
        make_all_neurons(folder)

    ttk.Separator(outer).pack(fill=tk.X, pady=(20, 14))
    style = ttk.Style(root)
    # On macOS an ordinary button's bezel is drawn at the system's one height, a line of text in its own size, whatever the button's: these
    # bigger buttons were bands of empty space round a bezel too short for them. A Toolbutton's bevel there takes the button's height (as the
    # UMAP viewer's "Reset slider default" buttons do), so the buttons are as tall as their padding makes them, and the gaps are only the pady.
    aqua = root.tk.call("tk", "windowingsystem") == "aqua"
    big = "Big.Toolbutton" if aqua else "Big.TButton"
    style.configure(big, font=("Segoe UI", 12), padding=(10, 9), anchor="center")
    ttk.Button(outer, text="Browse All Neurons", style=big, command=all_neurons).pack(fill=tk.X, pady=4)
    # choosing ROIs and loading a saved ROI project go together: set apart by space above and below the two
    ttk.Frame(outer, height=22).pack(fill=tk.X)
    ttk.Button(outer, text="Choose ROI(s) for project", style=big, command=new_run).pack(fill=tk.X, pady=4)
    ttk.Button(outer, text="Load saved ROI project…", style=big, command=load_existing).pack(fill=tk.X, pady=4)
    ttk.Frame(outer, height=22).pack(fill=tk.X)
    pair = ttk.Frame(outer)                                         # (the two less used ones share a row)
    pair.pack(fill=tk.X, pady=4)
    pair.columnconfigure((0, 1), weight=1, uniform="pair")
    ttk.Button(pair, text="Open .h5ad dataset…", style=big, command=own_dataset).grid(row=0, column=0, sticky="ew", padx=(0, 4))
    ttk.Button(pair, text="Open 3D Viewer…", style=big, command=viewer3d).grid(row=0, column=1, sticky="ew", padx=(4, 0))
    def exit_app():
        """Exit (or the window closed): while the atlas is still downloading, warn first, since closing the program stops the download."""
        running = helpers.get("atlas_download_running", lambda: None)()
        if running:
            with question():
                leave = messagebox.askyesno(
                    "The atlas is still downloading",
                    f"The atlas data is still downloading ({running}).\n\n"
                    "Exiting now pauses the download: it continues from where it stopped the next time the program is started.\n\n"
                    "Exit anyway?", icon="warning", default="no", parent=root)
            if not leave:
                return
        root.destroy()

    ttk.Button(outer, text="Exit", style=big, command=exit_app).pack(fill=tk.X, pady=4)

    root.protocol("WM_DELETE_WINDOW", exit_app)
    root.update_idletasks()
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry(f"+{(sw - root.winfo_width()) // 2}+{(sh - root.winfo_height()) // 3}")
    bring_to_front(root)
    root.after(50, helpers.get("launch_screen_shown", lambda: None))           # (the window is up: the start-up time is logged)
    root.mainloop()
    # Anything of this window still waiting for the cyclic collector is freed now, here, on the main thread: freed later by a worker thread
    # (the analysis's threads call gc.collect()), Tcl ends the program ("Tcl_AsyncDelete: async handler deleted by the wrong thread").
    try:
        root.destroy()
    except tk.TclError:
        pass
    root = update_button = None
    gc.collect()
    return result["outcome"]


if __name__ == "__main__":
    print(prompt_startup_panel_tk(
        Path.cwd(), None, {"load_recent_output_folders": lambda: [], "remember_recent_output_folder": lambda f: None,
                           "required_cached_run_files": lambda p: {}, "open_3d_viewer": lambda c: print("3D viewer"),
                           "settings_dialog": lambda c: print("settings")}))
