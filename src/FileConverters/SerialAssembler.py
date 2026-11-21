"""The start window: two steps, one after the other.

  1. Load / detect images - choose Zeiss .czi slide scans, find the sections on them, put them in order (section_loader).
  2. Align images         - the alignment tool (jpx_aligner); after step 1 it opens with those sections, in memory.

Run with a file name (a .jpx, or a project saved by step 1) to go straight to the alignment tool."""
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import messagebox, ttk

import neuroinfo_stack as ns

ns.ensure_packages({"PIL": "pillow"})

# --- updates: as in the other apps (app_update.py, in src: the setting is one for all the apps on this computer)

AUTO_UPDATE_NOTE = ("Applies to all the Brain-anatomy-tools apps on this computer. The assembler updates only while this start window is\n"
                    "showing, never while a tool is open. Off: the Check for updates button only says that updates are available.")
UPDATES = {"started": False, "done": False, "changes": None, "auto_tried": False}    # (the check at launch: once per run, not per start window)


def _updater():
    """app_update, or None: the packaged .exe is not a folder from GitHub (nothing to update), and app_update is in src, above this folder."""
    if getattr(sys, "frozen", False):
        return None
    try:
        import app_update
        return app_update
    except ImportError:
        return None


def _watch_for_updates(app_update):
    """At the first start window: looks on GitHub for newer commits in a background thread (its answer is kept in UPDATES, for whichever
    start window is showing when it comes). Silent on any failure (offline, ...): the button reports those when it is used."""
    if UPDATES["started"]:
        return
    UPDATES["started"] = True

    def work():
        try:
            UPDATES["changes"] = app_update.check(app_update.repo_root())
        except Exception:
            UPDATES["changes"] = None
        UPDATES["done"] = True
    threading.Thread(target=work, daemon=True).start()


def restart():
    """Starts the assembler again (in the same console window: its progress output stays visible) and ends this run."""
    args = [sys.executable, os.path.abspath(sys.argv[0])] + sys.argv[1:]
    sys.stdout.flush()
    sys.stderr.flush()
    if os.name == "nt":
        subprocess.Popen(args, close_fds=True)
        os._exit(0)                     # (not sys.exit: background threads would hold it open)
    os.execv(sys.executable, args)


def _starts(root, enabled):
    """The Start buttons on or off: off while an update replaces the program files (a tool started then could load half of the old
    version and half of the new)."""
    for b in getattr(root, "start_buttons", []):
        b.config(state=tk.NORMAL if enabled else tk.DISABLED)


def _idle_text(changes):
    return f"Update available ({len(changes)} change{'s' if len(changes) != 1 else ''})" if changes else "Check for updates"


def check_for_update(root, button, app_update):
    """The start window's Check for updates button: looks on GitHub for newer commits (in a thread: the network), says what is new, and on a
    yes brings the apps' folder up to date (git pull for a clone, else GitHub's ZIP: see app_update.py) and offers to restart the assembler.
    Changes the check at launch found are asked about at once, without looking again."""
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
        _starts(root, job is not app_update.pull)
        root.config(cursor="watch")
        threading.Thread(target=work, args=(job,) + args, daemon=True).start()

        def poll():
            if not root.winfo_exists():
                return
            try:
                ok, value = done.get_nowait()
            except queue.Empty:
                root.after(100, poll)
                return
            then(ok, value)
        poll()

    def finish(text=None, error=False):
        button.config(text=_idle_text(UPDATES["changes"]), state=tk.NORMAL)
        _starts(root, True)
        root.config(cursor="")
        if text:
            (messagebox.showerror if error else messagebox.showinfo)("Update", text, parent=root)

    def checked(ok, value):
        if not ok:
            return finish(value, True)
        UPDATES["changes"] = value
        if not value:
            return finish("You have the newest version.")
        shown = "\n".join("  • " + line for line in value[:12]) + (f"\n  ... and {len(value) - 12} more" if len(value) > 12 else "")
        news = f"{len(value)} new change{'s' if len(value) != 1 else ''} on GitHub:\n\n{shown}\n\n"
        reason = app_update.no_pull_reason(folder)
        if reason:                                  # (a git clone, and no git: what is new, and how to get it)
            return finish(news + reason)
        root.config(cursor="")
        if not messagebox.askyesno("Update", news + "Download them now?", parent=root):
            return finish()
        run(app_update.pull, folder, then=pulled, label="Updating…")

    def pulled(ok, value):
        if not ok:
            return finish(value, True)
        UPDATES["changes"] = None
        finish()
        if messagebox.askyesno("Update", "Updated. The new version runs after the assembler restarts.\n\nRestart it now?", parent=root):
            root.destroy()
            restart()

    try:
        folder = app_update.repo_root()
    except app_update.UpdateError as e:
        return finish(str(e), True)
    if UPDATES["changes"]:
        button.config(state=tk.DISABLED)
        checked(True, UPDATES["changes"])
        return
    run(app_update.check, folder, then=checked, label="Checking…")


def settings_dialog(root, app_update):
    """The start window's Settings: for now, only whether updates are downloaded automatically (saved at once, for all the apps)."""
    dlg = tk.Toplevel(root)
    dlg.title("Settings")
    dlg.resizable(False, False)
    dlg.transient(root)
    body = ttk.Frame(dlg, padding=16)
    body.pack(fill=tk.BOTH, expand=True)
    box = ttk.LabelFrame(body, text="Updates", padding=(10, 6))
    box.pack(fill=tk.X)
    auto = tk.BooleanVar(dlg, app_update.auto_update_enabled())
    ttk.Checkbutton(box, text="Download updates automatically at start-up, then restart", variable=auto,
                    command=lambda: app_update.set_auto_update(auto.get())).pack(anchor=tk.W)
    ttk.Label(box, text=AUTO_UPDATE_NOTE, foreground="gray30", justify=tk.LEFT).pack(anchor=tk.W, pady=(2, 0))
    ttk.Button(body, text="Close", command=dlg.destroy).pack(anchor=tk.E, pady=(12, 0))
    dlg.bind("<Escape>", lambda e: dlg.destroy())
    dlg.geometry(f"+{root.winfo_rootx() + 40}+{root.winfo_rooty() + 40}")
    dlg.grab_set()
    return dlg


def _update_buttons(root, header, app_update):
    """Settings and Check for updates, at the top right of the start window; the check at launch shown on the button as it comes in, and
    with "Download updates automatically" on, its updates downloaded (see app_update.AutoUpdate) and the assembler restarted."""
    ttk.Button(header, text="Settings…", command=lambda: settings_dialog(root, app_update)).pack(side=tk.RIGHT, anchor="n")
    button = ttk.Button(header, text=_idle_text(UPDATES["changes"]))
    button.configure(command=lambda: check_for_update(root, button, app_update))
    button.pack(side=tk.RIGHT, anchor="n", padx=(0, 6))
    _watch_for_updates(app_update)

    def blocked():
        try:
            top = root.grab_current()
            return top is not None and top is not root and top.winfo_exists()
        except (tk.TclError, KeyError):
            return False

    def status(text):
        button.config(text=text or _idle_text(UPDATES["changes"]), state=tk.DISABLED if text == "Updating..." else tk.NORMAL)
        _starts(root, text != "Updating...")

    def restart_now():
        root.destroy()
        restart()

    def wait():
        if not root.winfo_exists():
            return
        if not UPDATES["done"]:
            root.after(500, wait)
            return
        if str(button.cget("state")) != "disabled":
            button.config(text=_idle_text(UPDATES["changes"]))
        if UPDATES["changes"] and app_update.auto_update_enabled() and not UPDATES["auto_tried"]:
            UPDATES["auto_tried"] = True
            app_update.AutoUpdate(root, {"busy": lambda: str(button.cget("state")) == "disabled", "blocked": blocked, "status": status,
                                         "keep_open": lambda: None, "restart": restart_now}).start()
    root.after(500, wait)


def choose():
    """The start window. Returns "detect", "align" or None (closed)."""
    root = tk.Tk()
    root.title("Section alignment")
    root.resizable(False, False)
    choice = []
    outer = ttk.Frame(root, padding=28)
    outer.pack()
    header = ttk.Frame(outer)
    header.pack(fill=tk.X)
    ttk.Label(header, text="Brain section assembler", font=("Segoe UI", 18, "bold")).pack(side=tk.LEFT, anchor="w")
    app_update = _updater()
    if app_update is not None:
        _update_buttons(root, header, app_update)
    ttk.Label(outer, text="Steps 1-3 load and order the sections. Step 4 picks up where they leave off.", foreground="gray30").pack(anchor="w", pady=(2, 18))

    def card(number, title, text, key, primary):
        f = tk.Frame(outer, bd=1, relief="solid", padx=16, pady=12, bg="white")
        f.pack(fill=tk.X, pady=6)
        tk.Label(f, text=number, font=("Segoe UI", 14, "bold"), fg="#1f5fa8", bg="white", width=9).pack(side=tk.LEFT)
        mid = tk.Frame(f, bg="white")
        mid.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(8, 16))
        tk.Label(mid, text=title, font=("Segoe UI", 13, "bold"), bg="white", anchor="w").pack(fill=tk.X)
        tk.Label(mid, text=text, bg="white", fg="gray30", anchor="w", justify="left", wraplength=360).pack(fill=tk.X)
        b = ttk.Button(f, text="Start", width=10, command=lambda: (choice.append(key), root.destroy()))
        b.pack(side=tk.RIGHT)
        root.start_buttons = getattr(root, "start_buttons", []) + [b]
        if primary:
            b.focus_set()

    card("Steps 1-3", "Load / detect / order images",
         "Choose slide scans (.czi), find the sections on them, and put the sections in order. Or open a .jpx stack to reorder its sections.", "detect", True)
    card("Step 4", "Align images",
         "Line up the sections and export the stack. Opens the sections from steps 1-3, or a JPX you already have.", "align", False)
    ttk.Separator(outer).pack(fill=tk.X, pady=(14, 8))
    ttk.Style(root).configure("Exit.TButton", font=("Segoe UI", 13))
    ttk.Button(outer, text="Exit", width=24, style="Exit.TButton", command=root.destroy).pack(side=tk.BOTTOM, ipady=6, pady=(0, 4))   # (centered)
    root.mainloop()
    return choice[0] if choice else None


def run_aligner(stack=None, path=None, open_dialog=False, back_to_start=False):
    """The alignment tool, until its window is closed. With back_to_start, closing it returns here (to show the start window again)."""
    import jpx_aligner as J
    app = J.Aligner(path if stack is None else None)
    app.return_to_launcher = back_to_start
    if stack is not None:
        app.show_stack(stack, path)
    elif open_dialog and path is None:
        app.root.after(200, app.open_dialog)
    app.run()
    return app.back_to_ordering


def detect_and_align():
    """Steps 1-3 (the loader), then step 4 (the alignment tool). The alignment tool's "Back to step 3" hides its window and shows the loader
    at the ordering tab with the same project; the loader's Next then brings the same alignment window back, as it was (zoom, brightness,
    the pictures loaded, loading carrying on meanwhile), with the sections in the new order. Only if the project changed in a way that
    needs it (outlines redrawn, another JPX) is a new alignment window made, from the project (its saved alignment is kept)."""
    import jpx_aligner
    import section_loader
    import section_project as sp
    aligner, resume = None, None
    while True:
        loader = section_loader.SectionLoader()
        if resume:
            loader.resume(resume, 2)
        stack = loader.run()
        if stack is None:                        # (the loader was closed: the hidden alignment window goes too)
            if aligner is not None:
                aligner.shutdown()
            return
        project = loader.project
        ordered = None
        if aligner is not None and aligner.reuse_info is not None:
            ordered = sp.reuse_order(aligner.reuse_info, project, aligner.stack, stack)
        if ordered is not None:
            aligner.show_again(ordered, project.z_spacing_um)
        else:
            if aligner is not None:
                aligner.shutdown()
            aligner = jpx_aligner.Aligner(None)
            aligner.return_to_launcher = True
            aligner.show_stack(stack, project.path)
            aligner.set_start_brightness(loader.aligner_brightness())   # (the ordering tab's Auto / Full scale)
            aligner.reuse_info = sp.reuse_info(project, stack)
        aligner.run()
        if not aligner.back_to_ordering:         # (the alignment window was closed: back to the start window)
            return
        resume = project.path


def warn_if_h5py_missing():
    """Up front, before any work: say so when h5py (which reads NeuroInfo case files) could not be installed. Returns
    False if the user chose to stop."""
    if not ns.H5PY_PROBLEM:
        return True
    from tkinter import messagebox
    root = tk.Tk()
    root.withdraw()
    try:
        return messagebox.askyesno(
            "h5py is missing",
            f"{ns.H5PY_PROBLEM}" + chr(10) * 2 + "Without it, NeuroInfo case files (.h5) cannot be read, so the alignment of a "
            "NeuroInfo stack cannot be recovered. Everything else works." + chr(10) * 2 + "Continue anyway?",
            icon="warning", parent=root)
    finally:
        root.destroy()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not warn_if_h5py_missing():
        return
    if argv:                                  # a file named: straight to the alignment tool
        run_aligner(path=argv[0])
        return
    while True:                               # (the start window comes back whenever a tool is closed)
        what = choose()
        if what == "align":
            run_aligner(open_dialog=True, back_to_start=True)
        elif what == "detect":
            detect_and_align()
        else:
            break
    sys.stdout.flush()                        # (leave at once: background threads may be mid-task)
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":
    try:      # (its own icon and taskbar button: app_identity.py is in src, above this folder; not in the packaged .exe)
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))     # (after this folder: its modules come first)
        from app_identity import use_app_identity
        use_app_identity("section_assembler")
        import previous_install             # (installed from a new folder: offer to bring over the earlier copy's settings and files)
        previous_install.offer("section_assembler")
    except ImportError:
        pass
    main()
