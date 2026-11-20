"""After the tools are installed from a new folder while an earlier copy is still on the computer: the first program opened from the new
folder offers to bring over what the earlier copy had, so that nothing is downloaded again and no results seem to be lost.

The installers (INSTALL_APPS_WINDOWS.bat via check_previous_install.ps1, INSTALL_APPS_MAC.command) find the earlier copy, ask whether the
new folder is to be the installed one, and if so leave a note (NOTE_PATH: "previous=<folder>" and "new=<folder>"). offer(), called by
each program as it starts, finds that note from the new folder and offers:
  - the settings that are kept in a copy's own folder (the ABC and CCF browsers' settings, the recent output folders; the others are
    kept per user already, and so are shared by every copy);
  - the ABC atlas download (up to about 60 GB): moved beside the new copy (on the same drive: a rename, at once), or used where it is;
  - the ABC browser's output folder (the user's results), in the same way.
Anything in the Downloads folder is strongly recommended to be moved: Downloads gets cleared out, by hand or automatically.

Also kept here, per user (in the settings file the programs share, SETTINGS_PATH): where the atlas download is, so that every copy
finds it, and a new copy never downloads it again (saved_atlas_location, remember_atlas_location).

Needs only tkinter and the standard library: it runs before anything else is loaded."""
import json
import os
import subprocess
import sys

TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))              # this copy's top folder (src's parent)
SETTINGS_PATH = os.path.join(os.path.expanduser("~"), ".brain_anatomy_tools_settings.json")    # (shared with app_update.py)
NOTE_PATH = os.path.join(os.path.expanduser("~"), ".brain_anatomy_tools_switch.txt")          # (written by the installers)
ATLAS_NAME = "abc_atlas_cache"                       # the atlas download's default folder, beside a copy (DEFAULT_ATLAS_CACHE_DIR)
OUTPUT_NAME = "ABC_atlas_browser_output"             # the ABC browser's default output folder, beside a copy (DEFAULT_OUTPUT_FOLDER)
ATLAS_LOCATION_FILE = os.path.join("src", "cache_local", "abc_atlas_cache_location.json")   # (a copy's own record of where the atlas is)
RECENT_FILE = os.path.join("src", "cache_local", "recent_output_folders.json")
SETTINGS_FILES = (os.path.join("src", "viewer_settings.json"),                 # the settings kept in a copy's own folder
                  os.path.join("src", "cache_local", "ccf_browser_settings.json"),
                  os.path.join("src", "cache_local", "3d_viewer_recent_csv.json"),
                  os.path.join("src", "cache_local", "3d_viewer_last_csv.txt"))


# --- the settings file the programs share

def _read_settings():
    try:
        with open(SETTINGS_PATH, encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _write_setting(key, value):
    data = _read_settings()
    data[key] = value
    try:
        with open(SETTINGS_PATH, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    except OSError:
        pass


def saved_atlas_location(copy_file):
    """Where the atlas download is, as a path text, or None: this copy's own record (`copy_file`, its abc_atlas_cache_location.json)
    first, else the one kept per user, which any copy saved. Of the two, one that is an existing folder is preferred (a record of a folder
    that is gone is still returned when it is the only one, so that the caller can say so)."""
    found = []
    try:
        with open(copy_file, encoding="utf-8") as f:
            found.append(json.load(f).get("path"))
    except (OSError, ValueError, AttributeError):
        pass
    found.append(_read_settings().get("atlas_location"))
    found = [str(p) for p in found if p]
    return next((p for p in found if os.path.isdir(p)), found[0] if found else None)


def remember_atlas_location(folder):
    """Keep `folder` as the atlas download's place for every copy (the caller keeps its own record as well)."""
    _write_setting("atlas_location", str(folder))


# --- folders

def downloads_folders():
    """The Downloads folder(s): where the system says it is (on Windows it can be moved, e.g. to another drive), and the usual one."""
    found = []
    if sys.platform == "win32":
        try:
            import ctypes
            import uuid
            from ctypes import wintypes

            class GUID(ctypes.Structure):
                _fields_ = [("Data1", wintypes.DWORD), ("Data2", wintypes.WORD), ("Data3", wintypes.WORD), ("Data4", ctypes.c_ubyte * 8)]
            g = uuid.UUID("374DE290-123F-4565-9164-39C4925E467B")              # FOLDERID_Downloads
            guid = GUID(g.fields[0], g.fields[1], g.fields[2], (ctypes.c_ubyte * 8)(*g.bytes[8:]))
            path = ctypes.c_wchar_p()
            if ctypes.windll.shell32.SHGetKnownFolderPath(ctypes.byref(guid), 0, None, ctypes.byref(path)) == 0:
                found.append(path.value)
                ctypes.windll.ole32.CoTaskMemFree(path)
        except Exception:
            pass
        found.append(os.path.join(os.environ.get("USERPROFILE") or os.path.expanduser("~"), "Downloads"))
    else:
        found.append(os.path.join(os.path.expanduser("~"), "Downloads"))
    return [d for d in found if d]


def _key(path):
    return os.path.normcase(os.path.realpath(path)).rstrip("\\/") + os.sep


def in_downloads(path):
    """Whether `path` is the Downloads folder or in it (compared by the real path, without case where the system ignores it)."""
    here = _key(path)
    return any(here.startswith(_key(d)) for d in downloads_folders())


def same_folder(a, b):
    """Whether two paths are the same folder on disk: whatever their spelling (case, short names, links, a drive letter or the network
    name of the same share). False if either is not there."""
    try:
        return os.path.samefile(a, b)
    except OSError:
        return False


def same_drive(path, folder):
    """Whether `path` can be moved into `folder` by a rename (the same drive): no copying, at once."""
    try:
        return os.stat(path).st_dev == os.stat(folder).st_dev
    except OSError:
        return False


def folder_size(path):
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.path.getsize(os.path.join(root, name))
            except OSError:
                pass
    return total


def size_text(n):
    return f"{n / 1e9:.1f} GB" if n >= 1e9 else f"{n / 1e6:.0f} MB" if n >= 1e6 else "under 1 MB"


# --- the note the installers leave, and what it offers

def read_note():
    """(previous folder, new folder) from the installers' note, or None."""
    try:
        with open(NOTE_PATH, encoding="utf-8-sig") as f:
            pairs = dict(line.rstrip("\r\n").split("=", 1) for line in f if "=" in line)
    except OSError:
        return None
    previous, new = pairs.get("previous", "").strip(), pairs.get("new", "").strip()
    return (previous, new) if previous and new else None


def clear_note():
    try:
        os.remove(NOTE_PATH)
    except OSError:
        pass


def plan(top=TOP):
    """What there is to bring over to this copy (`top`) from the previous one the installers noted, or None (no note, or the note is not
    about this copy: then nothing is done, and the note is left for the copy it is about). A dict:
      previous, new: the two folders;
      settings: the settings files the previous copy has;
      atlas, output: each None (none found, or already beside the new copy) or a dict: path, size, in_downloads, movable (on the same
        drive, with no folder of that name beside the new copy yet), target (where a move puts it), why (why it can't be moved)."""
    note = read_note()
    if note is None or not same_folder(note[1], top):
        return None
    previous = note[0]
    if not os.path.isdir(previous) or same_folder(previous, top):
        clear_note()                                             # (nothing to bring over: the earlier copy is gone, or is this one)
        return None
    new_parent, old_parent = os.path.dirname(top), os.path.dirname(previous)
    result = {"previous": previous, "new": top, "previous_in_downloads": in_downloads(previous),
              "settings": [f for f in SETTINGS_FILES + (RECENT_FILE,) if os.path.isfile(os.path.join(previous, f))]}

    def item(path, name):
        if not path or not os.path.isdir(path) or same_folder(os.path.dirname(path), new_parent):
            return None                                          # (none, or beside the new copy already: it is found there as it is)
        target = os.path.join(new_parent, name)
        why = None
        if os.path.exists(target):
            why = f"a folder named {name} is beside the new copy already"
        elif not same_drive(path, new_parent):
            why = "it is on another drive than the new copy (moving it would mean copying it all)"
        return {"path": path, "size": folder_size(path), "in_downloads": in_downloads(path), "movable": why is None,
                "target": target, "why": why}

    atlas = saved_atlas_location(os.path.join(previous, ATLAS_LOCATION_FILE))
    if not atlas or not os.path.isdir(atlas):
        atlas = os.path.join(old_parent, ATLAS_NAME)
    result["atlas"] = item(atlas, ATLAS_NAME)
    # The output folder the previous copy was using: the one its launch screen showed (the most recent of its recent output folders,
    # as DEFAULT_OUTPUT_FOLDER's use in ABC_atlas_browser.py), which may be anywhere; the default one beside it only if it has none.
    output = os.path.join(old_parent, OUTPUT_NAME)
    try:
        with open(os.path.join(previous, RECENT_FILE), encoding="utf-8") as f:
            recent = json.load(f)
        if isinstance(recent, list) and recent and os.path.isdir(str(recent[0])):
            output = str(recent[0])
    except (OSError, ValueError):
        pass
    result["output"] = item(output, os.path.basename(os.path.normpath(output)) or OUTPUT_NAME)
    return result


def previous_copy_running(previous):
    """Whether a program of the previous copy is running (from its private Python), which would hold its files open."""
    try:
        if sys.platform == "win32":
            script = os.path.join(TOP, "src", "venv_in_use.ps1")
            done = subprocess.run(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", script, "-Venv",
                                   os.path.join(previous, ".venv")], capture_output=True, text=True, timeout=30,
                                  creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            return bool(done.stdout.strip())
        done = subprocess.run(["pgrep", "-f", os.path.join(previous, ".venv-macos", "bin") + os.sep], capture_output=True, timeout=10)
        return done.returncode == 0
    except Exception:
        return False                                             # (unknown: a move that then fails says so)


def carry_out(p, settings=True, move_atlas=False, move_output=False):
    """Bring over what was chosen. Returns (lines saying what was done, lines saying what could not be). Each part on its own: one that
    fails leaves the others done, and itself as it was."""
    done, failed = [], []
    new, previous = p["new"], p["previous"]
    moved = {}                                                   # old path -> new path, for the records that name them
    for name, item, chosen in (("atlas download", p.get("atlas"), move_atlas), ("output folder", p.get("output"), move_output)):
        if item is None or not chosen or not item["movable"]:
            continue
        try:
            os.rename(item["path"], item["target"])
            moved[item["path"]] = item["target"]
            done.append(f"Moved the {name} to {item['target']}.")
        except OSError as e:
            failed.append(f"The {name} could not be moved ({e.strerror or e}); it is still at {item['path']}.")

    def renamed(path):
        """`path` where it is now, if it was in (or was) a folder that moved."""
        for old, to in moved.items():
            if _key(path).startswith(_key(old)):
                rel = os.path.relpath(os.path.realpath(path), os.path.realpath(old))
                return to if rel == "." else os.path.join(to, rel)
        return path

    os.makedirs(os.path.join(new, "src", "cache_local"), exist_ok=True)
    if settings and p["settings"]:
        copied = 0
        for rel in p["settings"]:
            try:
                with open(os.path.join(previous, rel), "rb") as f:
                    data = f.read()
                if rel == RECENT_FILE:                           # (the folders that moved, at their new place)
                    data = json.dumps([renamed(str(x)) for x in json.loads(data.decode("utf-8"))], indent=2).encode("utf-8")
                with open(os.path.join(new, rel), "wb") as f:
                    f.write(data)
                copied += 1
            except (OSError, ValueError) as e:
                failed.append(f"Could not copy {rel} ({e}).")
        if copied:
            done.append(f"Copied the settings ({copied} file{'s' if copied != 1 else ''}).")
    atlas = p.get("atlas")
    if atlas is not None:                                        # the atlas, moved or not: where every copy (and this one) finds it
        where = moved.get(atlas["path"], atlas["path"])
        remember_atlas_location(where)
        try:
            with open(os.path.join(new, ATLAS_LOCATION_FILE), "w", encoding="utf-8") as f:
                json.dump({"path": where}, f, indent=2)
        except OSError:
            pass
        if atlas["path"] not in moved:
            done.append(f"The atlas download is used where it is: {where}.")
    output = p.get("output")
    if output is not None and output["path"] not in moved:      # (kept where it is: in the recent output folders, so it is found)
        recent = os.path.join(new, RECENT_FILE)
        try:
            with open(recent, encoding="utf-8") as f:
                folders = json.load(f)
        except (OSError, ValueError):
            folders = []
        if not any(same_folder(str(x), output["path"]) for x in folders):
            try:
                with open(recent, "w", encoding="utf-8") as f:
                    json.dump([output["path"]] + list(folders), f, indent=2)
            except OSError:
                pass
        done.append(f"The output folder is used where it is: {output['path']}.")
    return done, failed


# --- asking

DOWNLOADS_WARNING = ("Files in Downloads are often cleared out, by hand or automatically ("
                     + ("Windows' Storage Sense, and clean-up apps, delete files there" if sys.platform == "win32" else
                        "clean-up apps, and optimized storage, remove files there")
                     + "), and the atlas download and your results would go with them. Moving them out is strongly recommended.")


def offer(program="", top=TOP):
    """Called by each program as it starts: if the installers noted a previous copy and this is the new one, ask what to bring over, and
    do it. Never raises (it must not stop the program starting); does nothing in the packaged program, which is not a copy of this folder."""
    if getattr(sys, "frozen", False):
        return
    try:
        p = plan(top)
        if p is None:
            return
        _ask(p)
    except Exception as e:
        print(f"Could not offer to bring over the previous copy's settings and files ({e}).")


def _ask(p):
    import tkinter as tk
    from tkinter import ttk
    root = tk.Tk()
    root.withdraw()
    root.title("Previous copy found")
    root.resizable(False, False)
    try:
        root.attributes("-topmost", True)
    except tk.TclError:
        pass
    frame = ttk.Frame(root, padding=(20, 16, 20, 14))
    frame.pack(fill=tk.BOTH, expand=True)
    wrap = 760
    text = ("The Brain anatomy tools are now installed from:\n    " + p["new"] + "\n\nThe previous copy is still at:\n    " + p["previous"]
            + "\n\nChoose what the new copy should use from it:")
    ttk.Label(frame, text=text, wraplength=wrap, justify=tk.LEFT).pack(anchor="w")
    choices = {}

    def item_rows(key, name):
        it = p.get(key)
        if it is None:
            return
        box = ttk.LabelFrame(frame, text=f"{name} ({size_text(it['size'])}), at {it['path']}", padding=(10, 4))
        box.pack(fill=tk.X, pady=(10, 0))
        if it["in_downloads"]:                                # (the warning, by the folder it is about and the choice it argues for)
            ttk.Label(box, text="This is in your Downloads folder. " + DOWNLOADS_WARNING, wraplength=wrap - 30, justify=tk.LEFT,
                      foreground="firebrick").pack(anchor="w", pady=(0, 4))
        var = choices[key] = tk.StringVar(root, "move" if it["movable"] and it["in_downloads"] else "keep")
        move_text = f"Move it beside the new copy, to {it['target']} (instant: nothing is copied)"
        if it["in_downloads"] and it["movable"]:
            move_text += "  (strongly recommended)"
        mv = ttk.Radiobutton(box, text=move_text, value="move", variable=var)
        mv.pack(anchor="w")
        if not it["movable"]:
            mv.state(["disabled"])
            ttk.Label(box, text=f"(It can't be moved from here: {it['why']}. Move it yourself if you want it elsewhere.)",
                      foreground="gray35", wraplength=wrap - 30, justify=tk.LEFT).pack(anchor="w")
        keep_text = "Leave it where it is; the new copy uses it from there"
        if it["in_downloads"]:
            keep_text += "  (not recommended: it is in Downloads)"
        ttk.Radiobutton(box, text=keep_text, value="keep", variable=var).pack(anchor="w")

    item_rows("atlas", "The ABC atlas download")
    item_rows("output", "Your ABC browser results")

    if p["settings"]:                                         # (last: the folders above are what the question is mostly about)
        choices["settings"] = tk.BooleanVar(root, True)
        ttk.Checkbutton(frame, text="Also copy the settings (the ABC and CCF browsers' settings, the recent output folders)",
                        variable=choices["settings"]).pack(anchor="w", pady=(12, 0))

    answer = tk.StringVar(root, "")
    row = ttk.Frame(frame)
    row.pack(anchor="e", pady=(16, 0))
    # (two buttons, not a button for each way of doing it: the choices above already say move or keep for each folder)
    ttk.Button(row, text="OK", default="active", command=lambda: answer.set("go")).pack(side=tk.LEFT)
    ttk.Button(row, text="Not now", command=lambda: answer.set("later")).pack(side=tk.LEFT, padx=(8, 0))
    ttk.Label(frame, text="Nothing is deleted, whichever you choose. Not now asks again next time.",
              foreground="gray35", wraplength=wrap, justify=tk.LEFT).pack(anchor="w", pady=(8, 0))
    root.protocol("WM_DELETE_WINDOW", lambda: answer.set("later"))
    root.bind("<Return>", lambda e: answer.set("go"))
    root.bind("<Escape>", lambda e: answer.set("later"))
    root.update_idletasks()
    w, h = root.winfo_reqwidth(), root.winfo_reqheight()
    root.geometry(f"+{max((root.winfo_screenwidth() - w) // 2, 0)}+{max((root.winfo_screenheight() - h) // 3, 0)}")
    root.deiconify()
    root.lift()
    root.focus_force()
    root.wait_variable(answer)
    picked = {k: v.get() for k, v in choices.items()}
    root.destroy()
    import centered_dialog
    if answer.get() != "go":
        return
    moving = picked.get("atlas") == "move" or picked.get("output") == "move"
    if moving and previous_copy_running(p["previous"]):
        centered_dialog.info("Close the previous copy first",
                             "A program of the previous copy is open, and has its files open, so nothing was moved.\n\n"
                             "Close the Brain anatomy tools opened from\n    " + p["previous"] + "\nthen open this program again: it asks again.")
        return
    done, failed = carry_out(p, settings=picked.get("settings", False), move_atlas=picked.get("atlas") == "move",
                             move_output=picked.get("output") == "move")
    if not failed:
        clear_note()
    left = []
    for key in ("atlas", "output"):
        it = p.get(key)
        if it is not None and picked.get(key) != "move" and it["in_downloads"]:
            left.append(it["path"])
    message = "\n".join(done) if done else "Nothing needed changing."
    if failed:
        message += "\n\nNot done:\n" + "\n".join(failed) + "\n\nOpen this program again to try those once more."
    if left:
        message += "\n\nStill in Downloads (at risk of being cleared out):\n" + "\n".join("    " + x for x in left)
    message += ("\n\nThe previous copy is still at\n    " + p["previous"] + "\nOnce you have checked that everything is here, you can delete it.")
    centered_dialog.info("Done" if not failed else "Partly done", message)
