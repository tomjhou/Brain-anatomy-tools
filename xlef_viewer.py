"""XLEF viewer: browse the images of a Leica LAS X experiment (.xlef) and view their channels together.

The .xlef lists one .xlif (XML metadata) per image; each .xlif names the per-channel TIFFs. Only the small .xlif files are read to build the
list. An image's TIFFs are opened when it is selected, memory-mapped, and only the rows needed for the current view are read, so a zoomed-out
view of a 680 MB image reads a few tens of MB. OneDrive/SharePoint data: a file that is online-only is downloaded IN FULL as soon as any byte
of it is read, so the list marks such images and the viewer asks before opening one.

Usage:  python xlef_viewer.py [file.xlef]

Metadata used: <Frame File=...> per channel (and z), ChannelScalingInfo (LAS X contrast), WideFieldChannelInfo UserDefName (channel
names), DimensionDescription (size, pixel size). Test headless: build XlefViewer(), root.withdraw(), patch canvas.winfo_width/height,
SETTINGS_PATH and messagebox, pump root.update(), grab ImageTk.getimage(app.photo).
"""
import os
import re
import sys
import json
import queue
import threading
import xml.etree.ElementTree as ET
from urllib.parse import unquote
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import tifffile
from PIL import Image, ImageTk

SETTINGS_PATH = os.path.join(os.path.expanduser("~"), ".xlef_viewer_settings.json")

# Windows file attributes of OneDrive/SharePoint "online-only" placeholders (as in tifs_to_jp2.py).
_CLOUD_ONLY_ATTRS = 0x00400000 | 0x00040000 | 0x00001000  # RECALL_ON_DATA_ACCESS | RECALL_ON_OPEN | OFFLINE

LUT_COLORS = {"red": (1, 0, 0), "green": (0, 1, 0), "blue": (0, 0, 1), "cyan": (0, 1, 1), "magenta": (1, 0, 1),
              "yellow": (1, 1, 0), "gray": (1, 1, 1), "grey": (1, 1, 1), "white": (1, 1, 1)}
OVERVIEW_MAX = 2048      # longest side of the overview kept in memory per channel
MAX_ZOOM = 16.0          # screen pixels per image pixel


def is_online_only(path):
    """True for a cloud placeholder that isn't stored on this computer. Checking the attributes doesn't trigger a download."""
    try:
        return bool(getattr(os.stat(path), "st_file_attributes", 0) & _CLOUD_ONLY_ATTRS)
    except OSError:
        return False


# ---------------------------------------------------------------------------------------------------------------------- metadata

class Channel:
    def __init__(self, name, color, max_value, black, white):
        self.name, self.color, self.max_value = name, color, max_value
        self.black, self.white = black, white        # LAS X's display range, in pixel values


class ImageEntry:
    """One image (.xlif): its size, pixel size, channels and the TIFF of each channel (and z plane)."""

    def __init__(self, name, xlif_path):
        self.name, self.xlif_path = name, xlif_path
        self.width = self.height = 0
        self.nz = 1
        self.pixel_um = None
        self.channels = []
        self.files = {}          # (z, c) -> TIFF path
        self.problem = None      # why it can't be shown, if it can't

    def all_files(self):
        return sorted(set(self.files.values()))

    def online_only_files(self):
        return [f for f in self.all_files() if is_online_only(f)]


class FolderEntry:
    def __init__(self, name):
        self.name, self.children = name, []


def _resolve(ref, base_dir):
    """A Leica reference ("..%5Csection1%20snc--c00.tif") as a path, relative to the folder of the file that holds it."""
    rel = unquote(ref).replace("\\", os.sep).replace("/", os.sep)
    return os.path.normpath(os.path.join(base_dir, rel))


def _children(element, base_dir):
    """The entries an <Element>'s <Children> refer to (other .xlif/.xlef files)."""
    out = []
    kids = element.find("Children")
    if kids is None:
        return out
    for ref in kids.findall("Reference"):
        f = ref.get("File")
        if not f:
            continue
        path = _resolve(f, base_dir)
        if os.path.basename(path).lower().startswith("iomanagerconfigu"):
            continue
        entry = load_entry(path)
        if entry is not None:
            out.append(entry)
    return out


def load_entry(path):
    """An .xlif/.xlef file as an ImageEntry, a FolderEntry, or None if it holds neither."""
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as e:
        entry = ImageEntry(os.path.splitext(os.path.basename(path))[0], path)
        entry.problem = f"Can't read {os.path.basename(path)}: {e}"
        return entry
    element = root.find("Element")
    if element is None:
        return None
    name = element.get("Name") or os.path.splitext(os.path.basename(path))[0]
    base_dir = os.path.dirname(path)
    image = element.find("Data/Image")
    if image is not None:
        return _image_entry(name, path, element, image, base_dir)
    folder = FolderEntry(name)
    folder.children = _children(element, base_dir)
    return folder


def _image_entry(name, path, element, image, base_dir):
    entry = ImageEntry(name, path)
    desc = image.find("ImageDescription")
    dims = {d.get("DimID"): d for d in desc.iter("DimensionDescription")} if desc is not None else {}
    try:
        entry.width = int(dims["1"].get("NumberOfElements"))
        entry.height = int(dims["2"].get("NumberOfElements"))
        if "3" in dims:
            entry.nz = int(dims["3"].get("NumberOfElements"))
        length, unit = float(dims["1"].get("Length")), dims["1"].get("Unit", "m")
        scale = {"m": 1e6, "mm": 1e3, "um": 1.0, "µm": 1.0}.get(unit, 1e6)
        entry.pixel_um = length * scale / entry.width if entry.width else None
    except (KeyError, TypeError, ValueError):
        entry.problem = "The image size is missing from the metadata."
    if any(k not in ("1", "2", "3") for k in dims):
        entry.problem = entry.problem or ("This image has dimensions beyond x, y and z (time, tiles, ...), which this viewer doesn't show: "
                                          + ", ".join(sorted(dims)))

    # Channel names as typed in LAS X (widefield: WideFieldChannelInfo UserDefName; confocal: the dye name), and LAS X's display range.
    names = {}
    for info in element.iter("WideFieldChannelInfo"):
        names.setdefault(info.get("Channel"), info.get("UserDefName") or info.get("FluoCubeName") or "")
    names = list(names.values())
    scaling = list(element.iter("ChannelScalingInfo"))
    descs = list(desc.iter("ChannelDescription")) if desc is not None else []
    for i, cd in enumerate(descs):
        bits = int(cd.get("Resolution") or 16)
        top = 2 ** bits - 1
        nm = names[i] if i < len(names) and names[i] else ""
        if not nm:
            for prop in cd.iter("ChannelProperty"):
                if prop.findtext("Key") == "DyeName" and prop.findtext("Value"):
                    nm = prop.findtext("Value")
        lut = (cd.get("LUTName") or "Gray")
        nm = nm or f"Channel {i + 1} ({lut})"
        color = LUT_COLORS.get(lut.lower(), (1, 1, 1))
        black, white = 0, top
        if i < len(scaling):
            try:
                b, w = float(scaling[i].get("BlackValue")), float(scaling[i].get("WhiteValue"))
                if 0 <= b < w <= 1:
                    black, white = round(b * top), max(round(w * top), round(b * top) + 1)
            except (TypeError, ValueError):
                pass
        entry.channels.append(Channel(nm, color, top, black, white))

    # The TIFFs: one <Frame File=...> per channel (and z plane). File names carry --cNN (and --zNN); otherwise channels vary fastest.
    frames = [f.get("File") for f in element.iter("Frame") if f.get("File")]
    nch = max(len(entry.channels), 1)
    for k, ref in enumerate(frames):
        fpath = _resolve(ref, base_dir)
        low = os.path.basename(fpath).lower()
        mc, mz = re.search(r"--c(\d+)", low), re.search(r"--z(\d+)", low)
        c = int(mc.group(1)) if mc else k % nch
        z = int(mz.group(1)) if mz else (k // nch if len(frames) > nch else 0)
        entry.files[(z, c)] = fpath
    if not entry.problem:
        if not entry.files:
            entry.problem = "The metadata names no image files (perhaps the data is in a .lof file, which this viewer doesn't read)."
        elif any(f.lower().endswith(".lof") for f in entry.files.values()):
            entry.problem = "The data is in a .lof file, which this viewer doesn't read."
    return entry


def load_experiment(path):
    """The .xlef as a FolderEntry (its images and sub-folders)."""
    entry = load_entry(path)
    if isinstance(entry, FolderEntry):
        entry.name = os.path.splitext(os.path.basename(path))[0]
        return entry
    folder = FolderEntry(os.path.splitext(os.path.basename(path))[0])
    if entry is not None:
        folder.children.append(entry)
    return folder


# ---------------------------------------------------------------------------------------------------------------------- pixels

class OpenImage:
    """The selected image's TIFFs, memory-mapped (or read in full when a TIFF can't be mapped, e.g. compressed)."""

    def __init__(self, entry):
        self.entry = entry
        self.planes = {}         # (z, c) -> 2-D array
        self.overviews = {}      # (z, c) -> downsampled array
        self.ov_step = max(1, int(np.ceil(max(entry.width, entry.height) / OVERVIEW_MAX)))

    def plane(self, z, c):
        key = (z, c)
        if key not in self.planes:
            path = self.entry.files.get(key)
            if path is None:
                return None
            try:
                a = tifffile.memmap(path, mode="r")
            except (ValueError, OSError):
                a = tifffile.imread(path)
            while a.ndim > 2:        # a z stack in one file: pick the plane; extra sample axes: take the first
                a = a[z if a.shape[0] == self.entry.nz else 0]
            self.planes[key] = a
        return self.planes[key]

    def overview(self, z, c):
        key = (z, c)
        if key not in self.overviews:
            p = self.plane(z, c)
            if p is None:
                return None
            s = self.ov_step
            self.overviews[key] = np.ascontiguousarray(p[::s, ::s])
        return self.overviews[key]

    def close(self):
        self.planes.clear()
        self.overviews.clear()


def make_luts(channels, settings):
    """One (65536, 3) uint16 lookup table per shown channel: pixel value -> its RGB contribution (0-255)."""
    luts = []
    x = np.arange(65536, dtype=np.float32)
    for ch, st in zip(channels, settings):
        if not st["on"]:
            luts.append(None)
            continue
        lo, hi = st["lo"], max(st["hi"], st["lo"] + 1)
        v = np.clip((x - lo) / (hi - lo), 0, 1) * 255
        luts.append(np.stack([v * k for k in ch.color], axis=1).astype(np.uint16))
    return luts


def composite(arrays, luts):
    """The arrays (one per channel, same shape, or None) mixed into an RGB uint8 image by their lookup tables (additive)."""
    shape = next((a.shape for a in arrays if a is not None), (1, 1))
    acc = np.zeros(shape + (3,), np.uint16)
    for a, lut in zip(arrays, luts):
        if a is None or lut is None:
            continue
        if a.dtype == np.uint8:
            a = a.astype(np.uint16)
        elif a.dtype != np.uint16:
            a = np.clip(a, 0, 65535).astype(np.uint16)
        acc += lut[a]
    return np.minimum(acc, 255).astype(np.uint8)


def nice_length(um):
    """The largest 1/2/5 x 10^n not above um."""
    if um <= 0:
        return 0
    e = 10 ** np.floor(np.log10(um))
    for m in (5, 2, 1):
        if m * e <= um:
            return m * e
    return e


# ---------------------------------------------------------------------------------------------------------------------- the window

class RangeSlider(tk.Canvas):
    """One slider with two handles: black level (left) and white level (right). Drag a handle to move it, drag the bar between
    them to shift both, click elsewhere on the track to bring the nearer handle there. Shift while dragging moves 10x slower.
    command(lo, hi) is called on every change made with the mouse (not on set())."""

    PAD, TRACK_Y, HANDLE_W = 8, 24, 7

    def __init__(self, parent, max_value, lo, hi, color, command=None):
        bg = ttk.Style(parent).lookup("TFrame", "background") or "SystemButtonFace"
        super().__init__(parent, height=36, highlightthickness=0, background=bg, cursor="sb_h_double_arrow")
        self.max_value, self.color, self.command = max_value, color, command
        self.lo, self.hi = lo, hi
        self.drag = None
        self.bind("<Configure>", lambda e: self._draw())
        self.bind("<ButtonPress-1>", self._press)
        self.bind("<B1-Motion>", self._motion)
        self.bind("<ButtonRelease-1>", lambda e: setattr(self, "drag", None))

    def set(self, lo, hi):
        self.lo, self.hi = int(lo), int(hi)
        self._draw()

    def _span(self):
        return max(self.winfo_width() - 2 * self.PAD, 1)

    def _x(self, v):
        return self.PAD + v / self.max_value * self._span()

    def _value(self, x):
        return (x - self.PAD) / self._span() * self.max_value

    def _draw(self):
        self.delete("all")
        y, hw = self.TRACK_Y, self.HANDLE_W
        x0, x1 = self.PAD, self.PAD + self._span()
        xl, xh = self._x(self.lo), self._x(self.hi)
        self.create_rectangle(x0, y - 2, x1, y + 2, fill="#bbb", outline="")
        self.create_rectangle(xl, y - 4, xh, y + 4, fill=self.color, outline="#555")
        for x, fill in ((xl, "black"), (xh, "white")):
            self.create_polygon(x - hw / 2, y + 9, x + hw / 2, y + 9, x + hw / 2, y - 5, x, y - 9, x - hw / 2, y - 5,
                                fill=fill, outline="#333")
        self.create_text(x0, 1, anchor=tk.NW, text=f"black {self.lo}", fill="#333")
        self.create_text(x1, 1, anchor=tk.NE, text=f"white {self.hi}", fill="#333")

    def _press(self, event):
        xl, xh = self._x(self.lo), self._x(self.hi)
        reach = self.HANDLE_W
        if abs(event.x - xh) <= reach and (event.x >= xh or abs(event.x - xh) < abs(event.x - xl)):
            which = "hi"
        elif abs(event.x - xl) <= reach:
            which = "lo"
        elif xl < event.x < xh:
            which = "both"
        else:      # on the track outside the range: bring the nearer handle here
            which = "lo" if event.x < xl else "hi"
            self._apply(which, round(self._value(event.x)) - getattr(self, which))
        self.drag = (which, event.x, self.lo, self.hi)

    def _motion(self, event):
        if not self.drag:
            return
        which, x0, lo0, hi0 = self.drag
        dv = (event.x - x0) / self._span() * self.max_value * (0.1 if event.state & 0x0001 else 1)
        self._apply(which, round(dv), lo0, hi0)

    def _apply(self, which, dv, lo=None, hi=None):
        """Move handle(s) `which` by dv from (lo, hi) (default: where they are now)."""
        lo, hi, top = self.lo if lo is None else lo, self.hi if hi is None else hi, self.max_value
        if which == "both":
            dv = min(max(dv, -lo), top - hi)
            lo, hi = lo + dv, hi + dv
        elif which == "lo":
            lo = min(max(lo + dv, 0), hi - 1)
        else:
            hi = max(min(hi + dv, top), lo + 1)
        if (lo, hi) != (self.lo, self.hi):
            self.lo, self.hi = lo, hi
            self._draw()
            if self.command:
                self.command(lo, hi)


class XlefViewer:
    def __init__(self, path=None):
        self.root = tk.Tk()
        self.root.title("XLEF viewer")
        self.settings = self._load_settings()
        self.root.geometry(self.settings.get("geometry", "1400x900"))
        self.entries = {}            # tree item -> ImageEntry
        self.img = None              # OpenImage
        self.z = 0
        self.ch_settings = []        # per channel: {"on", "lo", "hi"}
        self.scale = 1.0             # screen pixels per image pixel
        self.cx = self.cy = 0.0      # image point at the canvas centre
        self.fit_scale = 1.0
        self.photo = None
        self.drag = None
        self.redraw_job = None
        self.generation = 0
        self.results = queue.Queue()
        self.loading = False
        self.last_view_rgb = None
        self._build()
        self.root.protocol("WM_DELETE_WINDOW", self.quit)
        self.root.after(50, self._poll_results)
        if path:
            self.root.after(10, lambda: self.open_xlef(path))

    # --- layout

    def _build(self):
        menu = tk.Menu(self.root)
        fm = tk.Menu(menu, tearoff=0)
        fm.add_command(label="Open .xlef...", accelerator="Ctrl+O", command=self.ask_open)
        fm.add_command(label="Save view as PNG...", command=self.save_view)
        fm.add_separator()
        fm.add_command(label="Quit", command=self.quit)
        menu.add_cascade(label="File", menu=fm)
        vm = tk.Menu(menu, tearoff=0)
        vm.add_command(label="Fit to window", accelerator="Home", command=self.fit)
        vm.add_command(label="Actual pixels (1:1)", command=lambda: self.zoom_to(1.0))
        menu.add_cascade(label="View", menu=vm)
        self.root.config(menu=menu)

        panes = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        panes.pack(fill=tk.BOTH, expand=True)
        left = ttk.Frame(panes, width=320)
        panes.add(left, weight=0)
        right = ttk.Frame(panes)
        panes.add(right, weight=1)

        top = ttk.Frame(left)
        top.pack(fill=tk.X, padx=4, pady=4)
        ttk.Button(top, text="Open .xlef...", command=self.ask_open).pack(side=tk.LEFT)
        self.file_label = ttk.Label(left, text="No file open", wraplength=300, foreground="#555")
        self.file_label.pack(fill=tk.X, padx=4)

        tree_frame = ttk.Frame(left)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.tree = ttk.Treeview(tree_frame, show="tree", selectmode="browse")
        sb = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.bind("<<TreeviewSelect>>", self.on_select)
        self.tree.tag_configure("cloud", foreground="#2060b0")
        self.tree.tag_configure("bad", foreground="#999")
        ttk.Label(left, text="☁ = online-only in OneDrive (downloaded when opened)", foreground="#2060b0").pack(anchor=tk.W, padx=4)

        self.z_frame = ttk.Frame(left)
        ttk.Label(self.z_frame, text="Z plane").pack(side=tk.LEFT)
        self.z_var = tk.IntVar(value=0)
        self.z_scale = tk.Scale(self.z_frame, from_=0, to=0, orient=tk.HORIZONTAL, variable=self.z_var, command=self.on_z)
        self.z_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

        self.ch_frame = ttk.LabelFrame(left, text="Channels")
        self.ch_frame.pack(fill=tk.X, padx=4, pady=4)
        btns = ttk.Frame(left)
        btns.pack(fill=tk.X, padx=4, pady=(0, 4))
        ttk.Button(btns, text="Auto contrast", command=self.auto_all).pack(side=tk.LEFT)
        ttk.Button(btns, text="LAS X contrast", command=self.lasx_all).pack(side=tk.LEFT, padx=4)

        self.canvas = tk.Canvas(right, background="black", highlightthickness=0, cursor="fleur")
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.status = ttk.Label(right, text="Open an .xlef file (File > Open, or drop it on RUN_XLEF_VIEWER.bat).", anchor=tk.W)
        self.status.pack(fill=tk.X)
        self.cursor_label = ttk.Label(right, text="", anchor=tk.W, foreground="#444")
        self.cursor_label.pack(fill=tk.X)

        c = self.canvas
        c.bind("<Configure>", lambda e: self.on_resize())
        c.bind("<ButtonPress-1>", self.on_press)
        c.bind("<B1-Motion>", self.on_drag)
        c.bind("<ButtonRelease-1>", self.on_release)
        c.bind("<MouseWheel>", self.on_wheel)
        c.bind("<Button-4>", lambda e: self.on_wheel(e, 120))
        c.bind("<Button-5>", lambda e: self.on_wheel(e, -120))
        c.bind("<Motion>", self.on_motion)
        c.bind("<Double-Button-1>", lambda e: self.fit())
        self.root.bind("<Control-o>", lambda e: self.ask_open())
        self.root.bind("<Home>", lambda e: self.fit())
        for k in range(1, 10):
            self.root.bind(str(k), lambda e, k=k: self.toggle_channel(k - 1))

    # --- settings

    def _load_settings(self):
        try:
            with open(SETTINGS_PATH, encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}

    def _save_settings(self):
        try:
            self.settings["geometry"] = self.root.geometry()
            with open(SETTINGS_PATH, "w", encoding="utf-8") as f:
                json.dump(self.settings, f, indent=1)
        except OSError:
            pass

    def quit(self):
        self._save_settings()
        self.generation += 1
        if self.img:
            self.img.close()
        self.root.destroy()

    # --- opening

    def ask_open(self):
        last = self.settings.get("last_file", "")
        path = filedialog.askopenfilename(parent=self.root, title="Open a Leica experiment",
                                          initialdir=os.path.dirname(last) if last else None,
                                          filetypes=[("Leica experiment", "*.xlef *.xlif"), ("All files", "*.*")])
        if path:
            self.open_xlef(path)

    def open_xlef(self, path):
        path = os.path.normpath(path)
        try:      # reads only the .xlef and .xlif files (a few hundred KB each, even when they must be downloaded)
            experiment = load_experiment(path)
        except Exception as e:      # a bad file must not take the window down
            messagebox.showerror("XLEF viewer", f"Can't read {path}:\n{e}", parent=self.root)
            return
        self.settings["last_file"] = path
        self._close_image()
        self.tree.delete(*self.tree.get_children())
        self.entries.clear()
        self._fill_tree("", experiment)
        self.root.title(f"XLEF viewer - {os.path.basename(path)}")
        self.file_label.config(text=path)
        n = len(self.entries)
        n_cloud = sum(1 for e in self.entries.values() if e.online_only_files())
        self.status.config(text=f"{n} image{'s' * (n != 1)}" + (f", {n_cloud} online-only" if n_cloud else "")
                           + ". Click one to view it.")
        first = next(iter(self.entries), None)
        if first and not self.entries[first].online_only_files():
            self.tree.selection_set(first)
            self.tree.see(first)

    def _fill_tree(self, parent, folder):
        for child in sorted(folder.children, key=lambda e: _natural_key(e.name)):
            if isinstance(child, FolderEntry):
                item = self.tree.insert(parent, tk.END, text=child.name, open=True)
                self._fill_tree(item, child)
                continue
            cloud = bool(child.online_only_files())
            text = child.name + ("  ☁" if cloud else "")
            tags = ("bad",) if child.problem else (("cloud",) if cloud else ())
            item = self.tree.insert(parent, tk.END, text=text, tags=tags)
            self.entries[item] = child

    def _close_image(self):
        self.generation += 1
        if self.img:
            self.img.close()
        self.img = None
        self.photo = None
        self.last_view_rgb = None
        self.canvas.delete("all")

    def on_select(self, event=None):
        sel = self.tree.selection()
        entry = self.entries.get(sel[0]) if sel else None
        if entry is None or (self.img and self.img.entry is entry):
            return
        if entry.problem:
            self._close_image()
            self.status.config(text=f"{entry.name}: {entry.problem}")
            return
        cloud = entry.online_only_files()
        if cloud:
            mb = sum(os.path.getsize(f) for f in cloud) / 1e6
            if not messagebox.askyesno("Download from OneDrive?",
                                       f"{len(cloud)} of {entry.name}'s files ({mb:,.0f} MB) are online-only.\n\n"
                                       "Viewing the image makes OneDrive download them in full. Continue?", parent=self.root):
                self.status.config(text=f"{entry.name} not opened (online-only).")
                return
        self.show_entry(entry)

    def show_entry(self, entry):
        same_channels = (len(entry.channels) == len(self.ch_settings) and self.img is not None
                         and [c.name for c in entry.channels] == [c.name for c in self.img.entry.channels])
        self._close_image()
        self.img = OpenImage(entry)
        self.z = min(self.z, entry.nz - 1)
        if not same_channels:      # keep the user's contrast when stepping through images of the same channels
            self.ch_settings = [{"on": True, "lo": c.black, "hi": c.white} for c in entry.channels]
        self._build_channel_controls()
        if entry.nz > 1:
            self.z_scale.config(to=entry.nz - 1)
            self.z_var.set(self.z)
            self.z_frame.pack(fill=tk.X, padx=4, before=self.ch_frame)
        else:
            self.z_frame.pack_forget()
        self.fit(redraw=False)
        self._load_overviews()

    def _load_overviews(self):
        """Read the overview of each channel in the background (the first read of a file is the slow part), showing each as it comes."""
        img, z, gen = self.img, self.z, self.generation
        self.loading = True
        n = len(img.entry.channels)
        mb = sum(os.path.getsize(f) for f in img.entry.all_files() if os.path.exists(f)) / 1e6

        def work():
            for c in range(n):
                if gen != self.generation:
                    return
                self.results.put(("progress", gen, f"Reading {img.entry.name}: channel {c + 1} of {n} ..."))
                try:
                    img.overview(z, c)
                except Exception as e:
                    self.results.put(("error", gen, f"Can't read channel {c + 1} of {img.entry.name}: {e}"))
                    return
                self.results.put(("overview", gen, None))
            self.results.put(("loaded", gen, f"{mb:,.0f} MB on disk"))

        threading.Thread(target=work, daemon=True).start()
        self.status.config(text=f"Reading {img.entry.name} ...")

    # --- channel controls

    def _build_channel_controls(self):
        for w in self.ch_frame.winfo_children():
            w.destroy()
        self.ch_widgets = []
        for i, (ch, st) in enumerate(zip(self.img.entry.channels, self.ch_settings)):
            f = ttk.Frame(self.ch_frame)
            f.pack(fill=tk.X, padx=2, pady=2)
            color = "#%02x%02x%02x" % tuple(int(200 * k + 30) for k in ch.color)
            on = tk.BooleanVar(value=st["on"])
            row = ttk.Frame(f)
            row.pack(fill=tk.X)
            tk.Label(row, text="  ", background=color).pack(side=tk.LEFT, padx=(0, 4))
            ttk.Checkbutton(row, text=f"{i + 1}. {ch.name}", variable=on,
                            command=lambda i=i, v=on: self._set(i, "on", v.get())).pack(side=tk.LEFT)
            ttk.Button(row, text="Auto", width=5, command=lambda i=i: self.auto(i)).pack(side=tk.RIGHT)
            slider = RangeSlider(f, ch.max_value, st["lo"], st["hi"], color, command=lambda lo, hi, i=i: self._set_lo_hi(i, lo, hi))
            slider.pack(fill=tk.X, pady=(2, 0))
            self.ch_widgets.append((on, slider))

    def _set(self, i, key, value):
        if i >= len(self.ch_settings) or self.ch_settings[i][key] == value:
            return
        self.ch_settings[i][key] = value
        self.schedule_redraw()

    def _set_lo_hi(self, i, lo, hi):
        if i < len(self.ch_settings) and (self.ch_settings[i]["lo"], self.ch_settings[i]["hi"]) != (lo, hi):
            self.ch_settings[i]["lo"], self.ch_settings[i]["hi"] = lo, hi
            self.schedule_redraw()

    def toggle_channel(self, i):
        if self.img and i < len(self.ch_settings) and not isinstance(self.root.focus_get(), tk.Entry):
            on = self.ch_widgets[i][0]
            on.set(not on.get())
            self._set(i, "on", on.get())

    def _set_range(self, i, lo, hi):
        self.ch_widgets[i][1].set(lo, hi)
        self.ch_settings[i]["lo"], self.ch_settings[i]["hi"] = lo, hi
        self.schedule_redraw()

    def auto(self, i):
        """Contrast from the overview: 0.1-99.9 percentile of the non-zero pixels (stitched scans have zero outside the tiles)."""
        if not self.img or (self.z, i) not in self.img.overviews:
            return
        a = self.img.overviews[(self.z, i)]
        v = a[a > 0]
        if v.size == 0:
            return
        lo, hi = np.percentile(v, [0.1, 99.9])
        self._set_range(i, int(lo), max(int(hi), int(lo) + 1))

    def auto_all(self):
        for i in range(len(self.ch_settings)):
            self.auto(i)

    def lasx_all(self):
        if self.img:
            for i, ch in enumerate(self.img.entry.channels):
                self._set_range(i, ch.black, ch.white)

    def on_z(self, value):
        z = int(float(value))
        if self.img and z != self.z:
            self.z = z
            self.generation += 1
            self._load_overviews()

    # --- view geometry

    def canvas_size(self):
        return max(self.canvas.winfo_width(), 1), max(self.canvas.winfo_height(), 1)

    def clamp_view(self):
        """No zooming out past the fit view, no panning past the image edges."""
        if not self.img:
            return
        e = self.img.entry
        cw, ch = self.canvas_size()
        self.fit_scale = min(cw / e.width, ch / e.height)
        self.scale = min(max(self.scale, self.fit_scale), MAX_ZOOM)
        for attr, size, view in (("cx", e.width, cw), ("cy", e.height, ch)):
            half = view / self.scale / 2
            v = getattr(self, attr)
            v = size / 2 if half * 2 >= size else min(max(v, half), size - half)
            setattr(self, attr, v)

    def fit(self, redraw=True):
        if not self.img:
            return
        self.scale = 0
        self.cx, self.cy = self.img.entry.width / 2, self.img.entry.height / 2
        self.clamp_view()
        if redraw:
            self.schedule_redraw(0)

    def zoom_to(self, scale, sx=None, sy=None):
        """Zoom to `scale`, keeping the image point under screen point (sx, sy) (default: the centre) in place."""
        if not self.img:
            return
        cw, ch = self.canvas_size()
        sx = cw / 2 if sx is None else sx
        sy = ch / 2 if sy is None else sy
        ix, iy = self.screen_to_image(sx, sy)
        self.scale = scale
        self.clamp_view()
        self.cx = ix - (sx - cw / 2) / self.scale
        self.cy = iy - (sy - ch / 2) / self.scale
        self.clamp_view()
        self.schedule_redraw(0)

    def screen_to_image(self, sx, sy):
        cw, ch = self.canvas_size()
        return self.cx + (sx - cw / 2) / self.scale, self.cy + (sy - ch / 2) / self.scale

    def on_resize(self):
        if self.img:
            at_fit = abs(self.scale - self.fit_scale) < 1e-9
            self.clamp_view()
            if at_fit:
                self.fit(redraw=False)
            self.schedule_redraw()

    def on_wheel(self, event, delta=None):
        d = event.delta if delta is None else delta
        self.zoom_to(self.scale * (1.25 ** (d / 120)), event.x, event.y)

    def on_press(self, event):
        self.drag = (event.x, event.y, self.cx, self.cy)

    def on_drag(self, event):
        if not self.drag or not self.img:
            return
        x0, y0, cx0, cy0 = self.drag
        self.cx = cx0 - (event.x - x0) / self.scale
        self.cy = cy0 - (event.y - y0) / self.scale
        self.clamp_view()
        self.redraw(preview_only=True)

    def on_release(self, event):
        self.drag = None
        self.schedule_redraw(0)

    def on_motion(self, event):
        if not self.img:
            return
        e = self.img.entry
        ix, iy = self.screen_to_image(event.x, event.y)
        if not (0 <= ix < e.width and 0 <= iy < e.height):
            self.cursor_label.config(text="")
            return
        x, y = int(ix), int(iy)
        parts = [f"x {x}, y {y} px"]
        if e.pixel_um:
            parts.append(f"({x * e.pixel_um:,.1f}, {y * e.pixel_um:,.1f} µm)")
        vals = []
        for c, ch in enumerate(e.channels):
            p = self.img.planes.get((self.z, c))     # only planes already open: hovering must not open files
            if p is not None and not self.loading:
                vals.append(f"{ch.name}: {int(p[y, x])}")
        self.cursor_label.config(text="   ".join(parts + vals))

    # --- drawing

    def schedule_redraw(self, delay=30):
        if self.redraw_job:
            self.root.after_cancel(self.redraw_job)
        self.redraw_job = self.root.after(delay, self._redraw_job)

    def _redraw_job(self):
        self.redraw_job = None
        self.redraw()

    def _view_rect(self, step):
        """The visible part of the image, in units of `step` image pixels (i0, j0, i1, j1: columns/rows of a [::step, ::step] array),
        and where it goes on the canvas (left, top, width, height)."""
        e = self.img.entry
        cw, ch = self.canvas_size()
        ix0, iy0 = self.screen_to_image(0, 0)
        ix1, iy1 = self.screen_to_image(cw, ch)
        i0, j0 = max(0, int(ix0 // step)), max(0, int(iy0 // step))
        i1, j1 = min(-(-e.width // step), int(-(-ix1 // step))), min(-(-e.height // step), int(-(-iy1 // step)))
        left = round((i0 * step - self.cx) * self.scale + cw / 2)
        top = round((j0 * step - self.cy) * self.scale + ch / 2)
        w = max(1, round((i1 - i0) * step * self.scale))
        h = max(1, round((j1 - j0) * step * self.scale))
        return (i0, j0, i1, j1), (left, top, w, h)

    def redraw(self, preview_only=False):
        """Draw at once from the overview; if the view needs more detail than the overview has, read it from the files in the background."""
        if not self.img:
            return
        n = len(self.img.entry.channels)
        luts = make_luts(self.img.entry.channels, self.ch_settings)
        s = self.img.ov_step
        ovs = [self.img.overviews.get((self.z, c)) if luts[c] is not None else None for c in range(n)]
        if all(lut is None for lut in luts):        # every channel switched off
            self.canvas.delete("all")
            self._draw_scale_bar()
            return
        if any(o is not None for o in ovs):
            (i0, j0, i1, j1), place = self._view_rect(s)
            if i1 > i0 and j1 > j0:
                crop = [o[j0:j1, i0:i1] if o is not None else None for o in ovs]
                self._show(composite(crop, luts), *place, smooth=self.scale * s < 1)
        step = max(1, int(1 / self.scale))
        if preview_only or self.loading or step >= s:
            return
        self.generation += 1
        gen, img, z = self.generation, self.img, self.z
        on = [c for c in range(n) if luts[c] is not None]
        (i0, j0, i1, j1), place = self._view_rect(step)
        smooth = self.scale * step < 1

        def work():
            try:
                arrays = [None] * n
                for c in on:
                    if gen != self.generation:
                        return
                    arrays[c] = img.plane(z, c)[j0 * step:j1 * step:step, i0 * step:i1 * step:step]
                self.results.put(("detail", gen, (composite(arrays, luts), *place, smooth)))
            except Exception as e:
                self.results.put(("error", gen, f"Can't read {img.entry.name}: {e}"))

        threading.Thread(target=work, daemon=True).start()

    def _show(self, rgb, left, top, w, h, smooth):
        im = Image.fromarray(rgb)
        if im.size != (w, h):
            im = im.resize((w, h), Image.BILINEAR if smooth else Image.NEAREST)
        self.last_view_rgb = im
        self.photo = ImageTk.PhotoImage(im, master=self.root)
        self.canvas.delete("all")
        self.canvas.create_image(left, top, anchor=tk.NW, image=self.photo)
        self._draw_scale_bar()
        self._update_status()

    def _draw_scale_bar(self):
        e = self.img.entry
        if not e.pixel_um:
            return
        cw, ch = self.canvas_size()
        um = nice_length(cw / 5 / self.scale * e.pixel_um)
        px = um / e.pixel_um * self.scale
        x, y = 16, ch - 20
        self.canvas.create_rectangle(x, y, x + px, y + 5, fill="white", outline="black")
        label = f"{um / 1000:g} mm" if um >= 1000 else f"{um:g} µm"
        self.canvas.create_text(x + px / 2, y - 3, text=label, fill="white", anchor=tk.S)

    def _update_status(self):
        e = self.img.entry
        size = f"{e.width} x {e.height} px"
        if e.pixel_um:
            size += f", {e.pixel_um:.3f} µm/px"
        z = f", z {self.z + 1}/{e.nz}" if e.nz > 1 else ""
        if not self.loading:
            self.status.config(text=f"{e.name}   {size}{z}   zoom {self.scale * 100:.1f}%")

    def _poll_results(self):
        try:
            while True:
                kind, gen, data = self.results.get_nowait()
                if kind in ("progress", "error", "loaded", "overview") and self.img is None:
                    continue
                if kind == "progress":
                    self.status.config(text=data)
                elif kind == "overview" and gen == self.generation:
                    self.redraw(preview_only=True)
                elif kind == "loaded" and gen == self.generation:
                    self.loading = False
                    self.redraw()
                elif kind == "detail" and gen == self.generation and not self.drag:
                    self._show(*data)
                elif kind == "error":
                    self.loading = False
                    self.status.config(text=data)
        except queue.Empty:
            pass
        self.root.after(30, self._poll_results)

    def save_view(self):
        if self.last_view_rgb is None:
            return
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".png", filetypes=[("PNG", "*.png")],
                                            initialfile=f"{self.img.entry.name}.png")
        if path:
            self.last_view_rgb.save(path)

    def run(self):
        self.root.mainloop()


def _natural_key(s):
    """Sort "Section2" before "Section10"."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


if __name__ == "__main__":
    XlefViewer(sys.argv[1] if len(sys.argv) > 1 else None).run()
