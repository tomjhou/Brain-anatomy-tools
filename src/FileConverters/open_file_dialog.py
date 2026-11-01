"""A reusable "open a saved file" dialog for Tk programs: a toggle at the top picks one of two panels, and Open acts on the panel that is showing.

  Recent <things>       a list of the recently used files (with optional extra columns, filled in by a background thread), double-click opens
  Browse for a <thing>  a "Choose file" button, a box for a typed or pasted path, and an area a file can be dragged onto from Explorer

The module knows nothing about any program: it is given the recent paths and a few optional callbacks, and returns what was chosen. It does not
open the file, check that it exists or remember anything (the caller keeps its own list of recent files), so it can be used as it is by any program.
Dragging needs tkinterdnd2 and a window made with it (TkinterDnD.Tk); without that the dialog works the same, minus the dragging.

    from open_file_dialog import ask_open_file
    result = ask_open_file(root, recent_paths, title="Open saved project", noun="project", filetypes=[("Projects", "*.json"), ("All files", "*.*")],
                           columns=[("name", "Project", 280, "w", True), ("sections", "Sections", 110, "e", False), ("folder", "Folder", 400, "w", True)],
                           details=lambda path: {"sections": count(path)}, on_remove=forget_recent, allow_delete=True)
    if result:
        path, from_list = result
"""
import os
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

SELECTED_BG = "#2f6fd4"      # the toggle button of the panel showing
IDLE_BG = "#e8e8e8"
IDLE_HOVER_BG = "#d8d8d8"
DEFAULT_COLUMNS = (("name", "File", 280, "w", True), ("folder", "Folder", 400, "w", True))


def ask_open_file(parent, recent, *, title="Open", noun="file", filetypes=(("All files", "*.*"),), columns=DEFAULT_COLUMNS, details=None,
                  on_remove=None, allow_delete=False, delete_note="", drop=True, initial_dir=None):
    """Show the dialog over `parent` and wait for it. Returns (path, from_list) for the file chosen (from_list: it was picked in the list, not
    typed, browsed for or dropped), or None if the dialog was cancelled.

    recent        the recently used paths, newest first (the list is shown first when it is not empty)
    noun          what a file is called in the dialog's words ("project", "recording", ...)
    filetypes     for the Browse button's file dialog
    columns       the list's columns: (key, heading, width, anchor, stretches). The keys "name" and "folder" are filled in from the path; any
                  other key is filled in from `details`
    details       details(path) -> {key: text}, called for each recent file in a background thread (reading a file can be slow, on a network
                  drive) while the dialog is already up; an OSError makes the row's last detail column read "(missing)"
    on_remove     on_remove(path): the caller forgets that path; when given, the dialog has a "Remove from list" button (and it is called after a
                  file is deleted too)
    allow_delete  show a "Delete file" button: the file is deleted from the disk, after asking
    delete_note   a few words after the Delete button's text, e.g. "raw files remain intact"
    drop          allow a file to be dragged in from Windows Explorer (if tkinterdnd2 and the window allow it)
    initial_dir   where the Browse dialog starts (default: the folder of the newest recent file)"""
    paths = list(recent)
    a_noun = ("an " if noun[:1].lower() in "aeiou" else "a ") + noun          # ("an experiment", "a project")
    dlg = tk.Toplevel(parent)
    r, g, b = dlg.winfo_rgb(dlg.cget("background"))
    light = (0.299 * r + 0.587 * g + 0.114 * b) / 65535 > 0.5         # (the program may use a dark theme: greys that read on it)
    muted, hint = ("#555555", "gray30") if light else ("#a8a8a8", "#b8b8b8")
    dlg.title(title)
    dlg.transient(parent)
    dlg.geometry(f"+{parent.winfo_rootx() + 120}+{parent.winfo_rooty() + 100}")
    mode = tk.StringVar(dlg, "recent" if paths else "browse")
    toggle = ttk.Frame(dlg, padding=(10, 10, 10, 6))
    toggle.pack(fill=tk.X)
    toggle_buttons = {}
    panels = ttk.Frame(dlg, padding=(10, 0))
    panels.pack(fill=tk.BOTH, expand=True)
    recent_panel = ttk.Frame(panels)
    browse_panel = ttk.Frame(panels)

    # ---- the file panel
    var_file = tk.StringVar(dlg)
    top = ttk.Frame(browse_panel)
    top.pack(fill=tk.X, pady=(6, 0))

    def browse():
        start = initial_dir or (os.path.dirname(paths[0]) if paths else "")
        file = filedialog.askopenfilename(parent=dlg, title=title, initialdir=start or None, filetypes=list(filetypes))
        if file:
            var_file.set(os.path.normpath(file))
    ttk.Button(top, text="Choose file", command=browse).pack(side=tk.LEFT)
    entry = ttk.Entry(top, textvariable=var_file, width=70)
    entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(6, 0))
    drop_note = tk.Label(browse_panel, text=f"...or drag {a_noun}{'' if noun == 'file' else ' file'} here from Windows Explorer", relief=tk.GROOVE, bd=2, fg=muted, height=5)
    drop_note.pack(fill=tk.BOTH, expand=True, pady=(10, 4))
    if drop:
        try:
            from tkinterdnd2 import DND_FILES

            def dropped(event):
                names = parent.tk.splitlist(event.data)
                if names:
                    var_file.set(os.path.normpath(names[0]))
                return event.action
            for target in (browse_panel, drop_note, entry):
                target.drop_target_register(DND_FILES)
                target.dnd_bind("<<Drop>>", dropped)
        except Exception as e:                               # (no tkinterdnd2, or the window was not made with it)
            print(f"Dragging a file into the dialog is not available ({e}).", flush=True)

    # ---- the recent-files panel
    keys = [c[0] for c in columns]
    detail_keys = [k for k in keys if k not in ("name", "folder")]
    tree = ttk.Treeview(recent_panel, columns=keys, show="headings", selectmode="browse", height=min(max(len(paths), 4), 12))
    for key, heading, width, anchor, stretch in columns:
        tree.heading(key, text=heading)
        tree.column(key, width=width, anchor=anchor, stretch=stretch)

    def row_values(path, fill):
        return tuple(os.path.basename(path) if k == "name" else os.path.dirname(path) if k == "folder" else fill for k in keys)
    # The dialog comes up at once, from the names alone; the details are filled in, row by row, by a thread.
    for i, p in enumerate(paths):
        tree.insert("", tk.END, iid=str(i), values=row_values(p, "..."))
    found = queue.Queue()

    def read_details():
        for i, p in enumerate(paths):
            try:
                found.put((i, details(p) if details else {}))
            except OSError:
                found.put((i, {k: "" for k in detail_keys} | ({detail_keys[-1]: "(missing)"} if detail_keys else {})))
            except Exception:
                found.put((i, {k: "" for k in detail_keys}))

    if detail_keys:
        threading.Thread(target=read_details, daemon=True).start()

        def show_details():
            try:
                while True:
                    i, values = found.get_nowait()
                    if tree.exists(str(i)):
                        for k in detail_keys:
                            tree.set(str(i), k, values.get(k, ""))
            except queue.Empty:
                pass
            except tk.TclError:
                return                                   # (the dialog was closed)
            if dlg.winfo_exists():
                dlg.after(100, show_details)
        dlg.after(50, show_details)
    tree.pack(fill=tk.BOTH, expand=True)
    ttk.Label(recent_panel, text=f"Double-click {a_noun}, or select it and click Open.", foreground=hint).pack(anchor="w", pady=(4, 0))
    if paths:
        tree.selection_set("0")
        tree.focus("0")
    chosen = []                                             # [(path, from the list)]

    def choose(event=None):
        if mode.get() == "browse":
            typed = var_file.get().strip().strip('"')
            if not typed:
                return
            chosen.append((typed, False))
        else:
            sel = tree.selection()
            if not sel:
                return
            chosen.append((paths[int(sel[0])], True))
        dlg.destroy()
    entry.bind("<Return>", choose)
    tree.bind("<Double-1>", choose)
    tree.bind("<Return>", choose)

    def forget():
        sel = tree.selection()
        if sel:
            on_remove(paths[int(sel[0])])
            tree.delete(sel[0])

    def delete():
        sel = tree.selection()
        if not sel:
            return
        p = paths[int(sel[0])]
        if not messagebox.askyesno(f"Delete {noun}", f"Delete the {noun} file{chr(10)}{chr(10)}{p}{chr(10)}{chr(10)}from the disk? Only this file goes"
                                   + (f" ({delete_note})" if delete_note else "") + ". This can't be undone.", icon="warning", default="no", parent=dlg):
            return
        try:
            if os.path.isfile(p):
                os.remove(p)
        except OSError as e:
            messagebox.showerror("Can't delete", str(e), parent=dlg)
            return
        if on_remove:
            on_remove(p)
        tree.delete(sel[0])

    # ---- the buttons below: Remove / Delete belong to the list, so they show with it only
    bar = ttk.Frame(dlg, padding=10)
    bar.pack(fill=tk.X)
    list_buttons = ttk.Frame(bar)
    list_buttons.pack(side=tk.LEFT)
    if on_remove:
        ttk.Button(list_buttons, text=f"Remove {noun} from list", command=forget).pack(side=tk.LEFT)
    if allow_delete:
        ttk.Button(list_buttons, text=f"Delete {noun} file" + (f" ({delete_note})" if delete_note else ""), command=delete).pack(side=tk.LEFT, padx=8)
    ttk.Button(bar, text="Cancel", command=dlg.destroy).pack(side=tk.RIGHT)       # (packed first: the rightmost)
    ttk.Button(bar, text="Open", command=choose).pack(side=tk.RIGHT, padx=(0, 8))

    # ---- the toggle: the panel showing is the one Open acts on
    def show_mode(which):
        mode.set(which)
        for key, button in toggle_buttons.items():
            on = key == which
            button.configure(bg=SELECTED_BG if on else IDLE_BG, fg="white" if on else "black", relief=tk.SUNKEN if on else tk.RAISED,
                             activebackground=SELECTED_BG if on else IDLE_HOVER_BG, activeforeground="white" if on else "black")
        recent_panel.pack_forget()
        browse_panel.pack_forget()
        if which == "recent":
            recent_panel.pack(fill=tk.BOTH, expand=True)
            list_buttons.pack(side=tk.LEFT)
            tree.focus_set()
        else:
            browse_panel.pack(fill=tk.BOTH, expand=True)
            list_buttons.pack_forget()
            entry.focus_set()
    for key, text in (("recent", f"Recent {noun}s"), ("browse", f"Browse for {a_noun}")):
        b = tk.Button(toggle, text=text, font=("Segoe UI", 10, "bold"), bd=2, padx=16, pady=4, cursor="hand2", command=lambda k=key: show_mode(k))
        b.pack(side=tk.LEFT, padx=(0, 6))
        toggle_buttons[key] = b
    show_mode(mode.get())
    dlg.bind("<Escape>", lambda e: dlg.destroy())
    dlg.grab_set()
    parent.wait_window(dlg)
    return chosen[0] if chosen else None
