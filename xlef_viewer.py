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
import time
import collections
import functools
import shutil
import hashlib
import uuid
import queue
import threading
import concurrent.futures
import xml.etree.ElementTree as ET
from urllib.parse import quote, unquote
import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk, filedialog, messagebox

import numpy as np
import tifffile
from PIL import Image, ImageTk

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "FileConverters"))
from open_file_dialog import ask_open_file
import leica_luts           # (LAS X's colour maps, sampled: see lut_table)
from tooltip import Tooltip, place_on_monitor, monitor_area

SETTINGS_PATH = os.path.join(os.path.expanduser("~"), ".xlef_viewer_settings.json")
THUMB_CACHE = os.path.join(os.environ.get("LOCALAPPDATA") or os.path.expanduser("~"), "JhouLab", "xlef_viewer_thumbs")
MAX_RECENT = 20
MOVING_PAUSE = 0.25          # s without a wheel click or drag before the view counts as at rest (see XlefViewer._moving)
HIRES_NOTE_DELAY = 0.5       # s of waiting for a sharper view before "Loading high-resolution images..." is shown over the image
KEPT_IMAGES_PERCENT = 10      # share of the computer's memory (%) for the pictures of recently viewed and preloaded images, which
                              # make switching to one instant (a setting); 500 MB if the memory's size can't be read
THUMB_SIZE = 128         # longest side of a thumbnail, screen pixels
THUMB_MAX_READ = 64e6    # a TIFF that can't be memory-mapped (compressed) is read whole for its thumbnail only up to this size (bytes)

# Windows file attributes of OneDrive/SharePoint "online-only" placeholders (as in tifs_to_jp2.py).
_CLOUD_ONLY_ATTRS = 0x00400000 | 0x00040000 | 0x00001000  # RECALL_ON_DATA_ACCESS | RECALL_ON_OPEN | OFFLINE

LUT_COLORS = {"red": (1, 0, 0), "green": (0, 1, 0), "blue": (0, 0, 1), "cyan": (0, 1, 1), "magenta": (1, 0, 1),
              "yellow": (1, 1, 0), "gray": (1, 1, 1), "grey": (1, 1, 1), "white": (1, 1, 1)}
# LAS X's colour maps (LUTName), in the order of its own list. A custom colour is stored as "Gradient (r,g,b)" (black to that colour).
LUT_ORDER = ["Gray", "Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "R&B", "Spectrum", "SpectrumLog", "GlowDark", "Glow",
             "Glow (O)", "Glow (U)", "Glow (O&U)", "Gray (O&U)", "P.Color 1", "P.Color 2", "P.Color 3", "P.Color 4", "P.Color 5",
             "P.Color 6", "P.Color 7", "Geo (L&S)", "Geo (L)", "Geo (S)", "Cyan Cool", "Rainbow"]
GRADIENT = re.compile(r"Gradient\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)", re.I)


@functools.lru_cache(maxsize=64)
def lut_table(name):
    """LAS X colour map `name` as (table, under, over): table = (256, 3) floats 0-1 from the dark end to the bright one; under/over =
    the colour (0-1) LAS X marks underexposed (0) / overexposed (maximum) pixels with in its (U)/(O) maps, else None. Plain colours and
    "Gradient (r,g,b)" go from black to that colour; the others come from leica_luts (sampled from LAS X); an unknown name is gray."""
    t = np.linspace(0, 1, 256)[:, None]
    m = GRADIENT.fullmatch(name.strip())
    if m:
        return t * (np.array([int(v) for v in m.groups()]) / 255), None, None
    if name.strip().lower() in LUT_COLORS:
        return t * np.array(LUT_COLORS[name.strip().lower()], float), None, None
    hexes = "".join(leica_luts.SAMPLED.get(name.strip(), ""))
    if not hexes:
        return t * np.ones(3), None, None
    samples = np.array([[int(hexes[k + j:k + j + 2], 16) for j in (0, 2, 4)] for k in range(0, len(hexes), 6)], float) / 255
    under = samples[0].copy() if "U)" in name else None            # ("Glow (U)", "Glow (O&U)", "Gray (O&U)")
    over = samples[-1].copy() if "(O" in name else None
    if under is not None:
        samples[:2] = samples[2]                                     # (the mark is applied exactly, not spread over the dark end)
    if over is not None:
        samples[-2:] = samples[-3]
    x = np.linspace(0, 1, len(samples))
    table = np.stack([np.interp(np.linspace(0, 1, 256), x, samples[:, k]) for k in range(3)], axis=1)
    return table, under, over


def lut_ui_color(name):
    """One colour (0-1) standing for a map, for the channel's swatch and slider: a plain colour itself, a map at 85% of its range."""
    table, _, _ = lut_table(name)
    return tuple(float(v) for v in table[int(0.85 * 255)])
OVERVIEW_MAX = 2048      # longest side of the overview kept in memory per channel
MAX_ZOOM = 16.0          # screen pixels per image pixel
PIXEL_GRID_COLOR = "#484848"     # lines along the pixels' edges at MAX_ZOOM, once the wheel tries to zoom in further


# ---------------------------------------------------------------------------------------------------------------------- color themes

# Light is Windows' own look (the native ttk theme, system colors); dark recolors ttk's "clam" theme. Everything this viewer draws itself
# takes its colors from PAL, the palette in use. (Windows draws the menu bar, message boxes and file dialogs itself, in its own colors.)
THEMES = {
    "light": {"bg": "SystemButtonFace", "fg": "SystemButtonText", "field": "SystemWindow", "field_fg": "SystemWindowText",
              "muted": "#555555", "link": "#2060b0", "warn": "#c03030", "bad": "#999999",
              "dim": "#7a7a7a", "link_dim": "#7f9fcb",
              "idle": "#e8e8e8", "idle_hover": "#d8d8d8", "idle_fg": "black", "selected": "#2f6fd4",
              "grid": "#f4f4f4", "tile": "white", "tile_line": "#d0d0d0", "tile_sel": "#cfe0fa", "tile_text": "black",
              "placeholder": "#e6e6e6", "placeholder_fg": "#777777", "download_bar": "#707070",
              "track": "#bbbbbb", "slider_text": "#333333", "handle_line": "#333333",
              "off_bg": "#c8c8c8", "off_fg": "#808080", "off_track": "#d4d4d4", "off_bar": "#b8b8b8", "off_black": "#8a8a8a",
              "off_white": "#e4e4e4", "off_line": "#999999"},
    "dark": {"bg": "#2b2b2b", "fg": "#e6e6e6", "field": "#1e1e1e", "field_fg": "#e6e6e6",
             "muted": "#a8a8a8", "link": "#6ea8ff", "warn": "#ff6b6b", "bad": "#7a7a7a",
             "dim": "#9a9a9a", "link_dim": "#4d76b3",
             "idle": "#3c3c3c", "idle_hover": "#4a4a4a", "idle_fg": "#e6e6e6", "selected": "#2f6fd4",
             "grid": "#232323", "tile": "#363636", "tile_line": "#4a4a4a", "tile_sel": "#28466f", "tile_text": "#e6e6e6",
             "placeholder": "#444444", "placeholder_fg": "#9a9a9a", "download_bar": "#ffffff",
             "track": "#5a5a5a", "slider_text": "#d0d0d0", "handle_line": "#9a9a9a",
             "off_bg": "#4a4a4a", "off_fg": "#8a8a8a", "off_track": "#404040", "off_bar": "#5a5a5a", "off_black": "#6a6a6a",
             "off_white": "#9a9a9a", "off_line": "#666666"},
}
PAL = dict(THEMES["dark"])


def style_theme(root, name, native):
    """Make `name` the theme in use: PAL, the ttk styles (live, for widgets already made), and the option database (for Tk widgets made
    from now on, e.g. dialogs). `native` is the ttk theme Windows started with (the light theme's)."""
    PAL.clear()
    PAL.update(THEMES[name])
    style = ttk.Style(root)
    root.option_clear()
    if name == "light":
        style.theme_use(native)
    else:
        p = PAL
        style.theme_use("clam")
        style.configure(".", background=p["bg"], foreground=p["fg"], fieldbackground=p["field"], bordercolor="#444444",
                        darkcolor=p["bg"], lightcolor="#3a3a3a", troughcolor=p["field"], selectbackground=p["selected"],
                        selectforeground="white", insertcolor=p["fg"], focuscolor=p["selected"])
        style.map(".", background=[("disabled", p["bg"]), ("active", p["idle_hover"])], foreground=[("disabled", "#777777")])
        style.configure("TButton", background=p["idle"], bordercolor="#555555", padding=(6, 2))
        style.map("TButton", background=[("disabled", p["bg"]), ("pressed", "#555555"), ("active", p["idle_hover"])])
        style.configure("Treeview", background=p["field"], fieldbackground=p["field"], foreground=p["field_fg"])
        style.map("Treeview", background=[("selected", p["selected"])], foreground=[("selected", "white")])
        style.configure("TCheckbutton", indicatorbackground=p["field"], indicatorforeground=p["fg"])
        style.map("TCheckbutton", indicatorbackground=[("disabled", p["bg"]), ("active", p["idle_hover"])])
        style.configure("TRadiobutton", indicatorbackground=p["field"], indicatorforeground=p["fg"])
        style.configure("TLabelframe", bordercolor="#555555")
        style.configure("TScrollbar", background=p["idle"], troughcolor=p["field"], arrowcolor=p["fg"], bordercolor="#444444")
        style.configure("TEntry", fieldbackground=p["field"], foreground=p["field_fg"])
        style.configure("TSeparator", background="#555555")
        style.configure("Sash", sashthickness=6, gripcount=0, background=p["bg"])
        for key, value in (("*Background", p["bg"]), ("*Foreground", p["fg"]), ("*activeBackground", p["idle_hover"]),
                           ("*activeForeground", p["fg"]), ("*highlightBackground", p["bg"]), ("*troughColor", p["field"]),
                           ("*Text.Background", p["field"]), ("*Text.Foreground", p["field_fg"]), ("*insertBackground", p["fg"]),
                           ("*selectBackground", p["selected"]), ("*selectForeground", "white"), ("*Entry.Background", p["field"])):
            root.option_add(key, value)
    for name_, key in (("Muted.TLabel", "muted"), ("Link.TLabel", "link"), ("Warn.TLabel", "warn")):
        style.configure(name_, foreground=PAL[key])


def style_toggle(button, on):
    """A toggle button (List / Thumbnails, the annotation tools): the one in use filled in."""
    p = PAL
    button.configure(bg=p["selected"] if on else p["idle"], fg="white" if on else p["idle_fg"], relief=tk.SUNKEN if on else tk.RAISED,
                     activebackground=p["selected"] if on else p["idle_hover"], activeforeground="white" if on else p["idle_fg"],
                     highlightbackground=p["bg"])


def recolor(widget, old, new):
    """Tk (not ttk) widgets made under the old palette: every color option that has an old palette color gets the new one. Widgets with
    their own colors (the image, the slide, the channel swatches, ...) are marked `keep_colors` and left alone, with their children."""
    if getattr(widget, "keep_colors", False):
        return
    swap = {}
    for key, value in old.items():
        if key in new and isinstance(value, str):
            swap.setdefault(str(value).lower(), new[key])
    for opt in ("background", "foreground", "activebackground", "activeforeground", "highlightbackground", "troughcolor",
                "insertbackground", "selectbackground", "disabledforeground"):
        try:
            value = str(widget.cget(opt)).lower()
        except tk.TclError:
            continue
        if value in swap:
            widget.configure(**{opt: swap[value]})
    for child in widget.winfo_children():
        recolor(child, old, new)


@functools.lru_cache(maxsize=None)
def _network_drive(root):
    """The drive (e.g. "M:\\", or a "\\\\server\\share\\" path) is on the network (asked of Windows once per drive)."""
    if root.startswith(("\\\\", "//")):
        return True
    try:
        import ctypes
        return ctypes.windll.kernel32.GetDriveTypeW(root) == 4            # (DRIVE_REMOTE)
    except (AttributeError, OSError):
        return False


def is_online_only(path):
    """True for a cloud placeholder that isn't stored on this computer. Checking the attributes doesn't trigger a download. A file on
    a network drive is never one (OneDrive keeps its files on a local drive), and isn't asked about: over a slow connection, each
    question is a round trip to the server."""
    drive = os.path.splitdrive(os.path.abspath(path))[0]
    if drive and _network_drive(drive.rstrip("\\/") + "\\"):
        return False
    try:
        return bool(getattr(os.stat(path), "st_file_attributes", 0) & _CLOUD_ONLY_ATTRS)
    except OSError:
        return False


def stored_bytes(path):
    """How much of a file is actually stored on this computer (for an online-only OneDrive file: 0, growing as it is downloaded),
    asked of Windows without opening the file (so no download is started). None if it can't be told."""
    try:
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.GetCompressedFileSizeW.argtypes = [wintypes.LPCWSTR, ctypes.POINTER(wintypes.DWORD)]
        k32.GetCompressedFileSizeW.restype = wintypes.DWORD
        high = wintypes.DWORD(0)
        low = k32.GetCompressedFileSizeW(path, ctypes.byref(high))
        if low == 0xFFFFFFFF and ctypes.get_last_error():
            return None
        return (high.value << 32) + low
    except (OSError, AttributeError):
        return None


AUTO_DOWNLOAD_MAX = 50e6     # an online-only file up to this size (bytes) is downloaded when needed without asking (pyramid levels are)


def size_text(n_bytes):
    """A file size for people: "640 MB", or from 1,000 MB on "1.4 GB"."""
    mb = n_bytes / 1e6
    return f"{mb:,.0f} MB" if mb < 999.5 else f"{mb / 1e3:,.1f} GB"


@functools.lru_cache(maxsize=None)
def readable_color(rgb, background, min_contrast=3.0):
    """A channel's color (rgb, 0-1) as text on `background` (rgb, 0-1): itself, or mixed toward white (dark background) or black (light
    background) until its contrast with the background reaches min_contrast (WCAG's ratio; 3 is enough for a short label). "#rrggbb"."""
    def luminance(c):
        lin = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in c]
        return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    lb = luminance(background)
    toward = (1.0, 1.0, 1.0) if lb < 0.5 else (0.0, 0.0, 0.0)
    for k in range(21):
        t = k / 20
        c = tuple(v + (w - v) * t for v, w in zip(rgb, toward))
        lc = luminance(c)
        if (max(lc, lb) + 0.05) / (min(lc, lb) + 0.05) >= min_contrast:
            break
    return "#%02x%02x%02x" % tuple(round(255 * v) for v in c)


REVEAL_LABEL = ("Open folder in Windows Explorer" if sys.platform == "win32" else "Show in Finder" if sys.platform == "darwin"
                else "Open folder in file manager")


def reveal_in_file_manager(path, select=True):
    """Open the folder of `path` in the system's file manager (Windows Explorer, Finder, or what xdg-open opens), with `path` selected
    where the system can do that (select: path is a file to point at; else a folder to open)."""
    import subprocess
    path = os.path.normpath(path)
    if sys.platform == "win32":
        if select:
            subprocess.Popen(f'explorer /select,"{path}"')
        else:
            os.startfile(path)
    elif sys.platform == "darwin":
        subprocess.Popen(["open", "-R", path] if select else ["open", path])
    else:
        subprocess.Popen(["xdg-open", os.path.dirname(path) if select else path])


WINDOW_FRAME = (16, 40)      # (pixels the title bar and borders add to a window's size, beyond what Tk's geometry counts)


def on_screen_geometry(widget, geometry):
    """A saved window geometry ("WxH+X+Y") moved, and if need be shrunk, so that the whole window is on the monitor it is mostly on
    (or the nearest one: a monitor since unplugged), inside its work area (not under the taskbar). Just "WxH": as it is."""
    m = re.fullmatch(r"(\d+)x(\d+)([+-]-?\d+)([+-]-?\d+)", geometry.strip())
    if not m:
        return geometry
    w, h = int(m.group(1)), int(m.group(2))
    x, y = int(m.group(3).replace("+-", "-").replace("+", "")), int(m.group(4).replace("+-", "-").replace("+", ""))
    left, top, right, bottom = monitor_area(widget, x + w // 2, y + h // 2)
    w = max(200, min(w, right - left - WINDOW_FRAME[0]))
    h = max(150, min(h, bottom - top - WINDOW_FRAME[1]))
    x = max(left, min(x, right - w - WINDOW_FRAME[0]))
    y = max(top, min(y, bottom - h - WINDOW_FRAME[1]))
    return f"{w}x{h}+{x}+{y}"


def total_memory():
    """The computer's physical memory (bytes), or None if it can't be read (read once)."""
    if sys.platform != "win32":
        try:
            return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
        except (ValueError, OSError, AttributeError):
            return None
    import ctypes

    class MemoryStatus(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong), ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong), ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong), ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong), ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    status = MemoryStatus()
    status.dwLength = ctypes.sizeof(status)
    return status.ullTotalPhys if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)) else None


def needs_asking(path):
    """True for a file that may only be read after asking: online-only AND bigger than AUTO_DOWNLOAD_MAX (reading any of it makes
    OneDrive download all of it). os.stat only: no download."""
    try:
        return is_online_only(path) and os.path.getsize(path) > AUTO_DOWNLOAD_MAX
    except OSError:
        return False


# ---------------------------------------------------------------------------------------------------------------------- metadata

class Channel:
    def __init__(self, name, lut, max_value, black, white, named=True):
        self.name, self.max_value = name, max_value
        self.lut = self.saved_lut = lut      # LAS X's colour map (LUTName): may be changed in the viewer; saved_lut is the file's
        self.named = named       # the name is from LAS X (False: made up here, e.g. "Channel 1 (Blue)")
        self.black, self.white = black, white        # display range, in pixel values: LAS X's, or as changed in the viewer
        self.saved_black, self.saved_white = black, white    # (as in the file)

    @property
    def color(self):
        """The channel's colour for its swatch and slider (0-1 RGB; for a multi-colour map, one colour standing for it)."""
        return lut_ui_color(self.lut)


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
        self.pyramid = None      # .xlif of LAS X's first reduced copy (_pmd_0), if it made a pyramid; each level names the next
        self.full_bytes = None   # its channel TIFFs' size on disk (bytes), once measured (see XlefViewer._measure_sizes)
        self.pages = {}          # (z, c) -> page of its file, when the channels are the pages of one OME-TIFF (see _ome_pages)
        self.annotations = []    # what was drawn on the image in LAS X (Annotation objects, hidden ones too, as in the file)
        self.annotations_changed = False     # edited in the viewer and not saved yet
        self.stage_um = None     # (x, y) of the image's CENTRE on the microscope stage, um (its <Tile PosX/PosY>), to place it on the slide
        self.saved_name = name   # its name in the file (self.name may be a rename not saved yet)
        self.parent_file = None  # the .xlef/.xlcf whose <Children> refer to this image's .xlif (None: unknown)
        self.flipped = self.saved_flipped = False    # "Flip this section": as changed in the viewer / as in the file (see write_flip)

    @property
    def path(self):
        return self.xlif_path

    @property
    def luts_changed(self):
        return any(c.lut != c.saved_lut for c in self.channels)

    @property
    def contrast_changed(self):
        return any((c.black, c.white) != (c.saved_black, c.saved_white) for c in self.channels)

    @property
    def flip_changed(self):
        return self.flipped != self.saved_flipped

    @property
    def unsaved(self):
        return (self.annotations_changed or self.name != self.saved_name or self.luts_changed or self.contrast_changed
                or self.flip_changed)

    def all_files(self):
        return sorted(set(self.files.values()))

    def online_only_files(self):
        return [f for f in self.all_files() if is_online_only(f)]


class FolderEntry:
    """A LAS X folder (.xlcf) or a nested experiment (.xlef). `path` is its file (None for a folder listed from the disk because its
    own file couldn't be read: it can't be renamed)."""

    def __init__(self, name, path=None):
        self.name, self.children, self.path = name, [], path
        self.saved_name = name
        self.parent_file = None
        self.annotations_changed = False
        self.order_changed = False   # what it lists, or their order, changed in the viewer (dragged, deleted) and not saved yet

    @property
    def unsaved(self):
        return self.name != self.saved_name or self.order_changed


def _resolve(ref, base_dir):
    """A Leica reference ("..%5Csection1%20snc--c00.tif") as a path, relative to the folder of the file that holds it."""
    rel = unquote(ref).replace("\\", os.sep).replace("/", os.sep)
    return os.path.normpath(os.path.join(base_dir, rel))


def _children(element, base_dir, container=None):
    """The entries an <Element>'s <Children> refer to (other .xlif/.xlef files); each remembers `container` (the file referring to it)."""
    out = []
    kids = element.find("Children")
    if kids is None:
        return out
    refs = []
    for ref in kids.findall("Reference"):
        f = ref.get("File")
        if not f:
            continue
        path = _resolve(f, base_dir)
        if not os.path.basename(path).lower().startswith("iomanagerconfigu"):
            refs.append((ref, path))
    # (read together, not one after another: on a network drive each file costs a round trip to the server, and over wifi those
    # added up to seconds; the order stays the file's)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        entries = list(pool.map(load_entry, [path for _, path in refs]))
    for (ref, _), entry in zip(refs, entries):
        if entry is not None:
            entry.parent_file = container
            entry.ref_uuid = ref.get("UUID") or ""
            out.append(entry)
    return out


def _folder_from_disk(xlcf_path):
    """A LAS X folder whose folder file (.xlcf) can't be read (one was found empty: 0 bytes): its contents as LAS X lays them out on
    disk instead, i.e. the images' .xlif files next to the .xlcf (in the folder's Metadata folder) and each subfolder's own .xlcf
    (<folder>/<sub>/Metadata/<sub>.xlcf). None if nothing is found."""
    meta = os.path.dirname(xlcf_path)
    here = os.path.dirname(meta)
    name = os.path.splitext(os.path.basename(xlcf_path))[0]
    try:        # (LAS X's references are lower case: the folder's name as it is on disk)
        name = next((d for d in os.listdir(os.path.dirname(here)) if d.lower() == os.path.basename(here).lower()), name)
    except OSError:
        pass
    folder = FolderEntry(f"{name} (folder file unreadable: images found on disk)")
    try:
        for f in sorted(os.listdir(meta)):
            low = f.lower()
            if low.endswith(".xlif") and "_pmd_" not in low and not low.endswith("_histo.xlif") and not low.startswith("iomanagerconfigu"):
                entry = load_entry(os.path.join(meta, f))
                if entry is not None:
                    folder.children.append(entry)
        for sub in sorted(os.listdir(here)):
            sub_meta = os.path.join(here, sub, "Metadata")
            if sub.lower() == "metadata" or not os.path.isdir(sub_meta):
                continue
            for f in sorted(os.listdir(sub_meta)):
                if f.lower().endswith(".xlcf"):
                    entry = load_entry(os.path.join(sub_meta, f))
                    if entry is not None:
                        folder.children.append(entry)
    except OSError:
        return None
    return folder if folder.children else None


def load_entry(path):
    """An .xlif/.xlef file as an ImageEntry, a FolderEntry, or None if it holds neither."""
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as e:
        if path.lower().endswith(".xlcf"):
            folder = _folder_from_disk(path)
            if folder is not None:
                return folder
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
    folder = FolderEntry(name, path)
    folder.children = _children(element, base_dir, path)
    return folder


def _ome_pages(entry, element, dims, descs, base_dir):
    """The other way LAS X stores an image: all its channels (and z planes) as the pages of one OME-TIFF, named by <OMETiffBlock
    File=... Offset=... Size=...> (or a few such blocks, each a part of the data). Which page holds which channel comes from the
    metadata alone (each channel's BytesInc, and the z dimension's, over one page's bytes): the TIFF is not opened, so an online-only
    one is not downloaded. Sets entry.files and entry.pages."""
    blocks = []
    for b in element.iter("OMETiffBlock"):
        try:
            blocks.append((int(b.get("Offset") or 0), int(b.get("Size") or 0), _resolve(b.get("File"), base_dir)))
        except (TypeError, ValueError):
            continue
    if not blocks or not entry.width or not entry.height:
        return
    try:
        page_bytes = int(dims["2"].get("BytesInc")) * entry.height          # (one row's bytes, times the rows)
        z_inc = int(dims["3"].get("BytesInc")) if "3" in dims else 0
    except (KeyError, TypeError, ValueError):
        return
    if page_bytes <= 0:
        return
    for c, cd in enumerate(descs):
        for z in range(entry.nz):
            try:
                at = int(cd.get("BytesInc") or 0) + z * z_inc
            except ValueError:
                return
            for offset, size, fpath in blocks:
                if offset <= at < offset + max(size, page_bytes) and (at - offset) % page_bytes == 0:
                    entry.files[(z, c)] = fpath
                    entry.pages[(z, c)] = (at - offset) // page_bytes
                    break


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
        names.setdefault(info.get("Channel"), info.get("UserDefName") or "")     # (FluoCubeName is the filter cube, e.g. "MICA_WF")
    # In the order of their numbers (2000, 2001, ...), not as listed: LAS X lists the active channel's info first, again.
    names = [names[k] for k in sorted(names, key=lambda k: int(k) if (k or "").isdigit() else 1 << 30)]
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
        named = bool(nm)
        nm = nm or f"Channel {i + 1} ({lut})"
        black, white = 0, top
        if i < len(scaling):
            try:
                b, w = float(scaling[i].get("BlackValue")), float(scaling[i].get("WhiteValue"))
                if 0 <= b <= w <= 1:
                    black, white = round(b * top), round(w * top)
                    if white == black and b < w:
                        white = black + 1
            except (TypeError, ValueError):
                pass
        entry.channels.append(Channel(nm, lut, top, black, white, named))

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
    if not frames:
        _ome_pages(entry, element, dims, descs, base_dir)
    kids = element.find("Children")
    for ref in (kids if kids is not None else []):
        if "_pmd_" in unquote(ref.get("File", "")).lower():
            entry.pyramid = _resolve(ref.get("File"), base_dir)
            break
    entry.annotations = [a for a in (Annotation.from_element(an) for an in element.iter("Annotation")) if a is not None]
    mine = next((at for at in element.iter("Attachment") if at.get("Name") == VIEWER_ATTACHMENT), None)
    entry.flipped = entry.saved_flipped = mine is not None and mine.get("FlipX") == "1"
    tile = next(element.iter("Tile"), None)
    try:
        if tile is not None:
            entry.stage_um = (float(tile.get("PosX")) * 1e6, float(tile.get("PosY")) * 1e6)
    except (TypeError, ValueError):
        pass
    if not entry.problem:
        if not entry.files:
            entry.problem = "The metadata names no image files (perhaps the data is in a .lof file, which this viewer doesn't read)."
        elif any(f.lower().endswith(".lof") for f in entry.files.values()):
            entry.problem = "The data is in a .lof file, which this viewer doesn't read."
    return entry


def count_images(path, _seen=None):
    """How many images an .xlef holds, counting those inside its folders (.xlcf, or a nested .xlef). Reads only those container files,
    not the images' own .xlif files."""
    seen = set() if _seen is None else _seen
    key = os.path.normcase(os.path.abspath(path))
    if key in seen:
        return 0
    seen.add(key)
    kids = ET.parse(path).getroot().find("Element/Children")
    n = 0
    for ref in (kids if kids is not None else []):
        if not ref.get("File"):
            continue
        child = _resolve(ref.get("File"), os.path.dirname(path))
        low = os.path.basename(child).lower()
        if low.startswith("iomanagerconfigu"):
            continue
        if low.endswith((".xlcf", ".xlef")):
            try:
                n += count_images(child, seen)
            except (OSError, ET.ParseError):
                if low.endswith(".xlcf"):           # (an unreadable folder file: count what is on disk, see _folder_from_disk)
                    folder = _folder_from_disk(child)
                    n += _count_folder(folder) if folder is not None else 0
        else:
            n += 1
    return n


def _count_folder(folder):
    return sum(_count_folder(c) if isinstance(c, FolderEntry) else 1 for c in folder.children)


# LAS X annotation types (<Annotation type=...> in an image's "Annotations" attachment), worked out from LAS X 1.15's own drawing of them.
ANNOTATION_KINDS = {128: "measure", 2: "line", 256: "arrow", 16: "rectangle", 32: "ellipse", 512: "text"}
TWO_POINT = ("measure", "line", "arrow")     # (drawn from one end to the other; an arrow's head is at its second point)
ANNOTATION_TYPES = {kind: t for t, kind in ANNOTATION_KINDS.items()}
ELLIPSE_POINTS = 100     # LAS X writes an ellipse as 100 points, at angles 0, 2pi/99, ..., 2pi (the last repeats the first)


class Annotation:
    """One LAS X annotation, kept in LAS X's own terms so that it can be edited and written back: its XML attributes (all of them, so
    that those this viewer doesn't use are written back unchanged) and its vertices.
    Geometry: metres from the image's CENTRE, y down. A vertex is scaled (transScalingX/Y), rotated by transRotation (degrees, clockwise
    on screen) and moved by transTransX/Y. A text's move is its top-left corner. Unknown types are shown as an open line through their
    vertices ("other") and written back as they were."""

    def __init__(self, attrs, verts):
        self.attrs = dict(attrs)
        self.verts = [tuple(v) for v in verts]

    @classmethod
    def from_element(cls, el):
        try:
            return cls(el.attrib, [(float(v.get("x")), float(v.get("y"))) for v in el.iter("Vertex")])
        except (TypeError, ValueError):
            return None

    @classmethod
    def new(cls, kind, color, center_m, verts, rotation=0.0, text=""):
        """A new annotation, with the attributes LAS X gives its own."""
        attrs = {"color": str(argb_of(color)), "fontName": "Arial", "fontSize": "10", "fontStyle": "Regular",
                 "type": str(ANNOTATION_TYPES[kind]), "name": "", "visibleLabel": "True", "visible": "True", "transRotation": "0",
                 "transTransX": "0", "transTransY": "0", "transScalingX": "1", "transScalingY": "1"}
        if kind == "measure":
            attrs["showLineEnds"] = "True"
        if kind == "arrow":
            attrs["linewidth"] = "6"     # (as LAS X writes its arrows)
        if kind == "text":
            attrs["text"] = text
            attrs["showBorder"] = "False"
            verts = [(0.0, 0.0)] * 4
        a = cls(attrs, verts)
        a.set_pose(center_m, rotation)
        return a

    def _f(self, key, default=0.0):
        try:
            return float(self.attrs.get(key) or default)
        except ValueError:
            return default

    @property
    def kind(self):
        try:
            return ANNOTATION_KINDS.get(int(self.attrs.get("type") or 0), "other")
        except ValueError:
            return "other"

    @property
    def visible(self):
        return self.attrs.get("visible", "True") != "False"

    @property
    def color(self):
        return "#%06x" % ((int(self.attrs.get("color") or -1) & 0xFFFFFFFF) & 0xFFFFFF)

    @color.setter
    def color(self, hex_color):
        self.attrs["color"] = str(argb_of(hex_color))

    @property
    def text(self):
        return (self.attrs.get("text") or "").replace("\r\n", "\n")

    @text.setter
    def text(self, value):
        self.attrs["text"] = value.replace("\r\n", "\n").replace("\n", "\r\n")

    @property
    def font(self):
        style = self.attrs.get("fontStyle") or ""
        return (self.attrs.get("fontName") or "Arial", int(self._f("fontSize", 10)), "bold" if "Bold" in style else "normal")

    @property
    def center(self):
        return self._f("transTransX"), self._f("transTransY")

    @property
    def rotation(self):
        return self._f("transRotation")

    def set_pose(self, center_m, rotation=None):
        self.attrs["transTransX"], self.attrs["transTransY"] = _num(center_m[0]), _num(center_m[1])
        if rotation is not None:
            self.attrs["transRotation"] = _num(rotation % 360)

    def local(self):
        """The vertices with the scaling applied (metres, before rotation and move)."""
        sx, sy = self._f("transScalingX", 1), self._f("transScalingY", 1)
        return [(x * sx, y * sy) for x, y in self.verts]

    def set_local(self, verts):
        """New vertices (metres, before rotation and move); the scaling is folded into them."""
        self.verts = [tuple(v) for v in verts]
        self.attrs["transScalingX"] = self.attrs["transScalingY"] = "1"

    def to_world(self, x, y):
        """A local point (metres) -> metres from the image centre."""
        r = np.radians(self.rotation)
        c, s = np.cos(r), np.sin(r)
        tx, ty = self.center
        return tx + x * c - y * s, ty + x * s + y * c

    def to_local(self, x, y):
        r = np.radians(self.rotation)
        c, s = np.cos(r), np.sin(r)
        dx, dy = x - self.center[0], y - self.center[1]
        return dx * c + dy * s, -dx * s + dy * c

    def points(self):
        """Where it is drawn: metres from the image centre (text: its top-left corner only)."""
        if self.kind == "text":
            return [self.center]
        return [self.to_world(x, y) for x, y in self.local()]

    def half_size(self):
        """Half width and half height of a rectangle or ellipse, in its own (unrotated) frame, metres."""
        loc = self.local()
        return max((abs(x) for x, _ in loc), default=0.0), max((abs(y) for _, y in loc), default=0.0)

    def set_half_size(self, a, b):
        a, b = max(a, 1e-7), max(b, 1e-7)
        if self.kind == "ellipse":
            t = np.linspace(0, 2 * np.pi, ELLIPSE_POINTS)
            self.set_local(list(zip(a * np.cos(t), b * np.sin(t))))
        else:
            self.set_local([(-a, -b), (a, -b), (a, b), (-a, b)])

    def length_label(self):
        """A line's length, as LAS X shows it, or ""."""
        if self.kind != "measure" or len(self.verts) != 2 or self.attrs.get("visibleLabel") == "False":
            return ""
        (x0, y0), (x1, y1) = self.local()
        return f"{np.hypot(x1 - x0, y1 - y0) * 1e6:.2f} µm"

    def to_xml(self, indent, nl):
        """The <Annotation> element, as LAS X writes it."""
        attrs = " ".join(f'{k}="{_xml_attr(v)}"' for k, v in self.attrs.items())
        lines = [f"{indent}<Annotation {attrs}>"]
        lines += [f'{indent}  <Vertex x="{_num(x)}" y="{_num(y)}" />' for x, y in self.verts]
        lines.append(f"{indent}</Annotation>")
        return nl.join(lines)


ANNOTATIONS_BLOCK = re.compile(r'([ \t]*)(<Attachment\b[^>]*\bName="Annotations"[^>]*?)(/>|>.*?</Attachment>)', re.S)
BACKUP_SUFFIX = ".before_xlef_viewer.bak"
VIEWER_ATTACHMENT = "JhouXlefViewer"     # this viewer's own <Attachment> in an image's .xlif: what LAS X has no place for (see write_flip)


def write_annotations(xlif_path, annotations):
    """Write `annotations` into the image's .xlif: ONLY its "Annotations" attachment is replaced (or added, at the end of the <Image>
    element, where LAS X puts it); every other byte of the file stays as it was (encoding, BOM, line endings). A copy of the file as it
    was before the viewer first changed it is kept next to it (BACKUP_SUFFIX). The new file is checked (it must read back, with these
    annotations) and written to a temporary file first, then moved over the old one. Returns the backup's path."""
    with open(xlif_path, "rb") as f:
        raw = f.read()
    bom = raw.startswith(b"\xef\xbb\xbf")
    text = raw.decode("utf-8-sig")
    nl = "\r\n" if "\r\n" in text else "\n"
    m = ANNOTATIONS_BLOCK.search(text)
    if m:
        indent, opening = m.group(1), m.group(2) + ">"
        start, end = m.start(), m.end()
    else:
        close = re.search(r"([ \t]*)</Image>", text)
        if close is None:
            raise ValueError("no <Image> element in the file")
        indent = close.group(1) + "  "
        opening = '<Attachment Name="Annotations" Application="LAS AF" Version="1.0">'
        start = end = close.start()
    body = [indent + opening] + [a.to_xml(indent + "  ", nl) for a in annotations] + [indent + "</Attachment>"]
    block = nl.join(body) + ("" if m else nl)
    new_text = text[:start] + block + text[end:]
    # (check before writing anything: the file must read back, with these annotations)
    root = ET.fromstring(new_text.encode("utf-8"))
    back = [Annotation.from_element(e) for e in root.iter("Annotation")]
    if len(back) != len(annotations) or any(b is None or b.attrs != a.attrs or len(b.verts) != len(a.verts)
                                            for a, b in zip(annotations, back)):
        raise ValueError("the annotations did not read back as written")
    return _write_text(xlif_path, new_text, bom)


def _real_case(path):
    """The path with its last part spelt as on disk (Windows ignores case, LAS X's references are lower case: "section1 snc.xlif" for
    "Section1 SNC.xlif"). Unchanged if it isn't found."""
    folder, name = os.path.split(path)
    try:
        return os.path.join(folder, next((n for n in os.listdir(folder) if n.lower() == name.lower()), name))
    except OSError:
        return path


def _read_text(path):
    """A LAS X XML file's text, and whether it starts with a byte-order mark (kept when it is written back)."""
    with open(path, "rb") as f:
        raw = f.read()
    return raw.decode("utf-8-sig"), raw.startswith(b"\xef\xbb\xbf")


def _write_text(path, new_text, bom):
    """Write a LAS X file the careful way: a copy of it as it was before the viewer first changed it is kept next to it (BACKUP_SUFFIX),
    and the new text goes to a temporary file first, which then replaces the old one. Returns the backup's path."""
    path = _real_case(path)          # (so that the file keeps its name exactly as it is: LAS X's references are lower case)
    backup = path + BACKUP_SUFFIX
    if not os.path.exists(backup):
        shutil.copy2(path, backup)
    tmp = path + ".xlef_viewer_tmp"
    with open(tmp, "wb") as f:
        f.write((b"\xef\xbb\xbf" if bom else b"") + new_text.encode("utf-8"))
    os.replace(tmp, path)
    return backup


ELEMENT_NAME = re.compile(r'(<Element\b[^>]*?\bName=")([^"]*)(")')
PREVIEW_MEMORY_NOTE = ("For the pictures of the images viewed recently, and of the ones loaded in the background so that they\n"
                       "open at once. Memory is only taken as images are loaded.")
RENAME_NOTE = ("Off: a rename changes only the name LAS X shows (in its .xlif / .xlcf).\n"
               "On: the image's .xlif and channel TIFFs, or the folder, get the new name too.")
FLIP_NOTE = ("On: Save writes it into the section's .xlif, in a field of this viewer's own (JhouXlefViewer), so it goes\n"
             "with the files to other computers. LAS X keeps such a field when it saves the file again (tried with\n"
             "LAS X 6.2), but this isn't guaranteed: a future version of LAS X could drop it, or object to it.\n"
             "Off: it is kept on this computer only (for each section's .xlif path), and LAS X's files are left as they are.\n"
             "Either way, a flip already saved in a section's .xlif is shown.")
BAD_NAME_CHARS = '\\/:*?"<>|'


def with_element_name(text, new_name):
    """The file's text with its <Element Name="..."> (the name LAS X shows) changed; checked to read back."""
    m = ELEMENT_NAME.search(text)
    if m is None:
        raise ValueError("the file has no element name")
    new_text = text[:m.start(2)] + _xml_attr(new_name) + text[m.end(2):]
    if ET.fromstring(new_text.encode("utf-8")).find("Element").get("Name") != new_name:
        raise ValueError("the new name did not read back")
    return new_text


def write_luts(path, luts):
    """Write the channels' colour maps into an image's .xlif, as LAS X does: only the LUTName of each <ChannelDescription> changes
    (in order; e.g. "Cyan", or "Gradient (255,119,104)" for a custom colour). Checked to read back; backed up (see _write_text)."""
    text, bom = _read_text(path)
    tags = list(re.finditer(r"<ChannelDescription\b[^>]*>", text))
    if len(tags) != len(luts):
        raise ValueError(f"the file has {len(tags)} channels, not {len(luts)}")
    out, last = [], 0
    for m, lut in zip(tags, luts):
        tag = m.group(0)
        new_tag, n = re.subn(r'\bLUTName="[^"]*"', f'LUTName="{_xml_attr(lut)}"', tag)
        if n == 0:
            new_tag = tag[:-1].rstrip("/").rstrip() + f' LUTName="{_xml_attr(lut)}"' + ("/>" if tag.endswith("/>") else ">")
        out += [text[last:m.start()], new_tag]
        last = m.end()
    new_text = "".join(out) + text[last:]
    back = [cd.get("LUTName") for cd in ET.fromstring(new_text.encode("utf-8")).iter("ChannelDescription")]
    if back != list(luts):
        raise ValueError("the colour maps did not read back as written")
    return _write_text(path, new_text, bom)


def write_contrast(path, ranges, tops):
    """Write the channels' display ranges (black, white in pixel values; tops: each channel's largest value) into an image's .xlif, as
    LAS X keeps them: BlackValue/WhiteValue (0-1) of each <ChannelScalingInfo> (in order, one per channel), Automatic="0". Nothing
    else changes. Checked to read back; backed up (see _write_text)."""
    text, bom = _read_text(path)
    tags = list(re.finditer(r"<ChannelScalingInfo\b[^>]*>", text))
    if len(tags) != len(ranges):
        raise ValueError(f"the file has contrast for {len(tags)} channels, not {len(ranges)}")
    out, last = [], 0
    for m, (lo, hi), top in zip(tags, ranges, tops):
        tag = m.group(0)
        for key, v in (("BlackValue", f"{lo / top:.12f}"), ("WhiteValue", f"{hi / top:.12f}"), ("Automatic", "0")):
            tag, n = re.subn(rf'\b{key}="[^"]*"', f'{key}="{v}"', tag)
            if n == 0:
                tag = tag[:-1].rstrip("/").rstrip() + f' {key}="{v}"' + (" />" if tag.endswith("/>") else ">")
        out += [text[last:m.start()], tag]
        last = m.end()
    new_text = "".join(out) + text[last:]
    back = [(round(float(e.get("BlackValue")) * top), round(float(e.get("WhiteValue")) * top))
            for e, top in zip(ET.fromstring(new_text.encode("utf-8")).iter("ChannelScalingInfo"), tops)]
    if back != [tuple(r) for r in ranges]:
        raise ValueError("the contrast did not read back as written")
    return _write_text(path, new_text, bom)


def write_flip(path, flipped):
    """Write "Flip this section" into an image's .xlif, in this viewer's own attachment (LAS X has no flip of its own; it keeps an
    attachment it doesn't know when it saves the file again):  <Attachment Name="JhouXlefViewer" Application="JhouXlefViewer"
    Version="1.0" FlipX="1" />  (FlipX="0" when unflipped). Added before </Image> if the file has none. Checked to read back; backed up
    (see _write_text)."""
    text, bom = _read_text(path)
    value = "1" if flipped else "0"
    m = re.search(rf'<Attachment\b[^>]*\bName="{VIEWER_ATTACHMENT}"[^>]*>', text)
    if m:
        tag, n = re.subn(r'\bFlipX="[^"]*"', f'FlipX="{value}"', m.group(0))
        if n == 0:
            tag = tag[:-1].rstrip("/").rstrip() + f' FlipX="{value}"' + (" />" if tag.endswith("/>") else ">")
        new_text = text[:m.start()] + tag + text[m.end():]
    else:
        close = re.search(r"([ \t]*)</Image>", text)
        if close is None:
            raise ValueError("no <Image> element in the file")
        nl = "\r\n" if "\r\n" in text else "\n"
        new_text = (text[:close.start()] + close.group(1) + "  " + f'<Attachment Name="{VIEWER_ATTACHMENT}" Application='
                    f'"{VIEWER_ATTACHMENT}" Version="1.0" FlipX="{value}" />' + nl + text[close.start():])
    back = [at.get("FlipX") for at in ET.fromstring(new_text.encode("utf-8")).iter("Attachment") if at.get("Name") == VIEWER_ATTACHMENT]
    if back != [value]:
        raise ValueError("the flip did not read back as written")
    return _write_text(path, new_text, bom)


def write_element_name(path, new_name):
    """Rename an image or folder the way LAS X shows it: only the Name in its own .xlif/.xlcf changes; file names stay."""
    text, bom = _read_text(path)
    return _write_text(path, with_element_name(text, new_name), bom)


def _leica_ref(target, base_dir):
    """A reference to `target` from a file in base_dir, written as LAS X writes them: relative, lower case, backslashes, URL-encoded
    (".%5Cmetadata%5Csection%2017.xlif")."""
    rel = os.path.relpath(target, base_dir).replace("/", "\\")
    if not rel.startswith(".."):
        rel = ".\\" + rel
    return quote(rel.lower(), safe="")


def _retarget(text, base_dir, moved):
    """The text with every File="..." reference to a moved file (moved: normcase(old path) -> new path) pointing at its new path."""
    def fix(m):
        target = moved.get(os.path.normcase(_resolve(m.group(2), base_dir)))
        return m.group(0) if target is None else m.group(1) + _leica_ref(target, base_dir) + m.group(3)
    return re.sub(r'(\bFile=")([^"]*)(")', fix, text)


def plan_disk_rename(item, new_name):
    """Everything renaming an image or folder ON DISK involves, worked out and checked before anything is done:
    [("text", path, new text), ("move", old path, new path)] (texts first, then moves, in this order). Raises ValueError if it can't
    be done safely (a name already taken, an unexpected layout, the referring file unknown).
      image: its .xlif and the channel TIFFs named after it ("<name>--C00.tif") get the new name; its .xlif's references to the TIFFs
             and its <Element Name> are changed, and so is the reference to it in the folder (or experiment) file that lists it.
      folder: its directory and its .xlcf ("<name>/Metadata/<name>.xlcf") get the new name, its <Element Name> too, and the reference
             to it in the file that lists it. What is inside keeps its names (references inside are relative, so they stay right).
    Pyramid files keep their names: they are found through references, not names."""
    if any(c in new_name for c in BAD_NAME_CHARS) or new_name != new_name.strip() or new_name.endswith("."):
        raise ValueError(f'a file name can\'t contain {BAD_NAME_CHARS} or end in a space or a dot')
    if not item.parent_file or not os.path.isfile(item.parent_file):
        raise ValueError("the file that lists it is unknown or unreadable")

    def free(old, new):
        if os.path.normcase(old) != os.path.normcase(new) and os.path.exists(new):
            raise ValueError(f"{os.path.basename(new)} already exists")
    moves = []
    if isinstance(item, ImageEntry):
        own = _real_case(item.xlif_path)
        old = os.path.splitext(os.path.basename(own))[0]
        new_own = os.path.join(os.path.dirname(own), new_name + ".xlif")
        free(own, new_own)
        for f in item.all_files():
            f = _real_case(f)
            b = os.path.basename(f)
            if b.lower().startswith(old.lower() + "--"):
                nf = os.path.join(os.path.dirname(f), new_name + b[len(old):])
                free(f, nf)
                moves.append((f, nf))
        text, _ = _read_text(own)
        own_text = with_element_name(_retarget(text, os.path.dirname(own), {os.path.normcase(a): b for a, b in moves}), new_name)
        moves.append((own, new_own))
        listed = {os.path.normcase(own): new_own}
    else:
        own = _real_case(item.path) if item.path else None
        meta = os.path.dirname(own or "")
        here = _real_case(os.path.dirname(meta))
        old = os.path.basename(here)
        if (not own or os.path.basename(meta).lower() != "metadata"
                or os.path.splitext(os.path.basename(own))[0].lower() != old.lower()):
            raise ValueError("its files are not laid out as LAS X does (<name>/Metadata/<name>.xlcf)")
        new_here = os.path.join(os.path.dirname(here), new_name)
        free(here, new_here)
        text, _ = _read_text(own)
        own_text = with_element_name(text, new_name)
        moves += [(own, os.path.join(meta, new_name + ".xlcf")), (here, new_here)]
        listed = {os.path.normcase(own): os.path.join(new_here, "Metadata", new_name + ".xlcf")}
    parent_text, _ = _read_text(item.parent_file)
    new_parent_text = _retarget(parent_text, os.path.dirname(item.parent_file), listed)
    if new_parent_text == parent_text:
        raise ValueError("the file that lists it doesn't refer to it")
    return [("text", own, own_text), ("text", item.parent_file, new_parent_text)] + [("move", a, b) for a, b in moves]


def carry_out(plan):
    """Do a plan_disk_rename() plan; if any step fails, undo the steps done (texts restored, moves reversed) and raise."""
    done = []
    try:
        for kind, a, b in plan:
            if kind == "text":
                old_bytes = open(a, "rb").read()
                _write_text(a, b, old_bytes.startswith(b"\xef\xbb\xbf"))
                done.append(("text", a, old_bytes))
            elif kind == "mkdir":
                if not os.path.isdir(a):
                    os.makedirs(a)
                    done.append(("mkdir", a, None))
            else:
                shutil.move(a, b)             # (a rename on the same drive; between drives a copy, then the original deleted)
                done.append(("move", a, b))
    except Exception:
        for kind, a, b in reversed(done):
            try:
                if kind == "text":
                    with open(a, "wb") as f:
                        f.write(b)
                elif kind == "mkdir":
                    os.rmdir(a)
                else:
                    shutil.move(b, a)
            except OSError:
                pass
        raise


def experiment_home(container):
    """The folder a LAS X experiment or folder keeps its things in: the .xlef's own folder, or for a folder's .xlcf (<name>/Metadata/
    <name>.xlcf) the <name> folder."""
    d = os.path.dirname(container)
    return os.path.dirname(d) if container.lower().endswith(".xlcf") and os.path.basename(d).lower() == "metadata" else d


def owned_files(xlif_path):
    """Every file an image's .xlif brings with it: itself, its channel TIFFs (or .lof), and what it refers to in turn (LAS X's reduced
    copies, _pmd_N, each an .xlif with TIFFs of its own; histograms). Read from the references, not guessed from names."""
    out, todo = [], [os.path.normpath(xlif_path)]
    while todo:
        path = todo.pop()
        if any(os.path.normcase(path) == os.path.normcase(p) for p in out):
            continue
        out.append(path)
        if not path.lower().endswith(".xlif"):
            continue
        text, _ = _read_text(path)
        for ref in re.findall(r'\bFile="([^"]*)"', text):
            target = _resolve(ref, os.path.dirname(path))
            if not os.path.basename(target).lower().startswith("iomanagerconfigu"):
                todo.append(target)
    return out


def _moved_text(text, old_dir, new_dir, moved):
    """An .xlif's text for its new place: every File="..." reference (read from old_dir) written relative to new_dir, pointing at the
    file's new path if it moves too (moved: normcase(old path) -> new path)."""
    def fix(m):
        old = _resolve(m.group(2), old_dir)
        return m.group(1) + _leica_ref(moved.get(os.path.normcase(old), old), new_dir) + m.group(3)
    return re.sub(r'(\bFile=")([^"]*)(")', fix, text)


def plan_move_image(entry, src_root, dst_container, dst_root):
    """Moving an image into another experiment (or folder), files and all, as a carry_out() plan: its .xlif goes to the destination's
    Metadata folder, the files beside it (channel TIFFs) to the destination's folder, LAS X's reduced copies (in the experiment's Pyramid
    folder) to the destination experiment's Pyramid folder; the references inside the moved .xlif files are rewritten for their new
    places. src_root/dst_root: the two experiments' folders. Checked first: raises ValueError (nothing done) if a file is missing,
    lies outside its experiment, or its new name is taken. Returns (plan, {normcase(old path): new path})."""
    meta = os.path.dirname(entry.xlif_path)          # (an image's .xlif: <folder>/Metadata/<name>.xlif)
    src_home = os.path.dirname(meta) if os.path.basename(meta).lower() == "metadata" else meta
    dst_home = experiment_home(dst_container)
    files = owned_files(entry.xlif_path)
    moved = {}
    for f in files:
        f = _real_case(f)
        if not os.path.exists(f):
            raise ValueError(f"{os.path.basename(f)} is missing")
        pyramid = os.path.join(src_root, "Pyramid")
        if os.path.normcase(f).startswith(os.path.normcase(pyramid + os.sep)):
            new = os.path.join(dst_root, "Pyramid", os.path.relpath(f, pyramid))
        elif os.path.normcase(f).startswith(os.path.normcase(src_home + os.sep)):
            new = os.path.join(dst_home, os.path.relpath(f, src_home))
        elif os.path.normcase(f).startswith(os.path.normcase(src_root + os.sep)):
            new = os.path.join(dst_root, os.path.relpath(f, src_root))
        else:
            raise ValueError(f"{os.path.basename(f)} lies outside its experiment's folder")
        if os.path.exists(new):
            raise ValueError(f"{new} already exists")
        moved[os.path.normcase(f)] = new
    plan = [("mkdir", d, None) for d in sorted({os.path.dirname(n) for n in moved.values()})]
    for f in files:
        if f.lower().endswith(".xlif"):
            new = moved[os.path.normcase(_real_case(f))]
            text, _ = _read_text(f)
            plan.append(("text", _real_case(f), _moved_text(text, os.path.dirname(f), os.path.dirname(new), moved)))
    plan += [("move", _real_case(f), moved[os.path.normcase(_real_case(f))]) for f in files]
    return plan, moved


def children_text(container, listed, children):
    """The container's (.xlef/.xlcf) text with the images and folders it shows set to `children`, in that order: [(path, uuid)].
    listed: the paths of those it showed when it was read. Its other references (LAS X's settings file, ...) keep their places; the
    shown ones take the places of those shown before (all together, where the first was). Checked to read back."""
    text, _ = _read_text(container)
    base = os.path.dirname(container)
    if len(re.findall(r"<Children\b", text)) != 1:
        raise ValueError("its list of contents is not laid out as expected")
    nl = "\r\n" if "\r\n" in text else "\n"
    m = re.search(r"([ \t]*)<Children\s*/>", text)
    if m:                                    # (an empty list: <Children />)
        text = text[:m.start()] + f"{m.group(1)}<Children>{nl}{m.group(1)}</Children>" + text[m.end():]
    head_indent = re.search(r"([ \t]*)<Children>", text).group(1)
    m = re.search(r"<Children>(.*?)[ \t]*</Children>", text, re.S)
    refs = re.findall(r"[ \t]*<Reference\b[^>]*/>", m.group(1))
    indent = re.match(r"[ \t]*", refs[0]).group(0) if refs else head_indent + "  "
    shown = {os.path.normcase(p) for p in listed} | {os.path.normcase(p) for p, _ in children}

    def target(ref):
        f = re.search(r'\bFile="([^"]*)"', ref)
        return os.path.normcase(_resolve(f.group(1), base)) if f else ""
    new_lines = [f'{indent}<Reference File="{_leica_ref(p, base)}" UUID="{_xml_attr(u)}" />' for p, u in children]
    out, last = [], 0                # (the shown ones fill the places of those shown before, in the new order; any more go after)
    for ref in refs:
        if target(ref) in shown:
            if new_lines:
                out.append(new_lines.pop(0))
                last = len(out)
        else:
            out.append(ref)
    out[last:last] = new_lines
    new_text = text[:m.start(1)] + nl + "".join(line + nl for line in out) + text[m.end(1):]
    kids = ET.fromstring(new_text.encode("utf-8")).find("Element/Children")
    back = [os.path.normcase(_resolve(r.get("File", ""), base)) for r in (kids if kids is not None else [])]
    if [b for b in back if b in shown] != [os.path.normcase(p) for p, _ in children]:
        raise ValueError("its new list of contents did not read back as written")
    return new_text


def send_to_recycle_bin(paths):
    """Move files to the Windows Recycle Bin (they can be restored from there). Raises OSError if Windows refuses."""
    import ctypes
    from ctypes import wintypes

    class SHFILEOPSTRUCTW(ctypes.Structure):
        _fields_ = [("hwnd", wintypes.HWND), ("wFunc", ctypes.c_uint), ("pFrom", wintypes.LPCWSTR), ("pTo", wintypes.LPCWSTR),
                    ("fFlags", ctypes.c_ushort), ("fAnyOperationsAborted", wintypes.BOOL), ("hNameMappings", ctypes.c_void_p),
                    ("lpszProgressTitle", wintypes.LPCWSTR)]
    FO_DELETE, FOF_SILENT, FOF_NOCONFIRMATION, FOF_ALLOWUNDO, FOF_NOERRORUI = 3, 0x4, 0x10, 0x40, 0x400
    op = SHFILEOPSTRUCTW(None, FO_DELETE, "\0".join(os.path.abspath(p) for p in paths) + "\0\0", None,
                         FOF_SILENT | FOF_NOCONFIRMATION | FOF_ALLOWUNDO | FOF_NOERRORUI, False, None, None)
    err = ctypes.windll.shell32.SHFileOperationW(ctypes.byref(op))
    if err or op.fAnyOperationsAborted:
        raise OSError(f"Windows could not move them to the Recycle Bin (error {err:#x})")


def has_recycle_bin(path):
    """Does the drive of `path` have a Recycle Bin? Local drives do; network drives don't (Windows would delete outright there)."""
    import ctypes
    drive = os.path.splitdrive(os.path.abspath(path))[0] + "\\"
    return ctypes.windll.kernel32.GetDriveTypeW(drive) in (2, 3)        # (DRIVE_REMOVABLE, DRIVE_FIXED)


def argb_of(hex_color):
    """'#rrggbb' -> LAS X's color number (opaque ARGB as a signed 32-bit integer, e.g. red = -65536)."""
    v = 0xFF000000 | int(hex_color.lstrip("#"), 16)
    return v - (1 << 32) if v >= 1 << 31 else v


def _num(v):
    """A number as LAS X writes them (e.g. 0.00090446, 1.121371E-06)."""
    return f"{float(v):.7G}"


def _xml_attr(value):
    return (str(value).replace("&", "&amp;").replace('"', "&quot;").replace("<", "&lt;").replace(">", "&gt;")
            .replace("\r", "&#xD;").replace("\n", "&#xA;"))


def _filetime(hex_value):
    """A LAS X time stamp (Windows FILETIME, hexadecimal: 100 ns steps since 1601) as a local datetime."""
    import datetime
    t = datetime.datetime(1601, 1, 1, tzinfo=datetime.timezone.utc) + datetime.timedelta(microseconds=int(hex_value, 16) // 10)
    return t.astimezone()


def _clock(t):
    """'1:49:40.365 PM' (LAS X's way)."""
    return f"{t.hour % 12 or 12}:{t:%M:%S}.{t.microsecond // 1000:03d} {'AM' if t.hour < 12 else 'PM'}"


def _dec(v, digits=1):
    return f"{float(v):.{digits}f}"


def image_properties(entry, experiment_path):
    """What LAS X's Properties window shows for an image, read from its .xlif: {"header": [(label, value)], "sections": [(title,
    [lines])], "stamps": [(frame, relative s, clock, date)]}. Parts LAS X stores only as codes this viewer can't name are left out."""
    root = ET.parse(entry.xlif_path).getroot()
    el = root.find("Element")
    n = len(entry.channels)
    header = [("Image", entry.name), ("File Location", experiment_path or entry.xlif_path)]
    mem = el.find("Memory")
    if mem is not None and mem.get("Size"):
        header.append(("Size", f"{int(mem.get('Size')) / 1e6:,.2f} MB"))
    stamps_el = el.find(".//TimeStampList")
    times = [_filetime(h) for h in (stamps_el.text or "").split()] if stamps_el is not None else []
    if times:
        header += [("Start Time", f"{times[0].month}/{times[0].day}/{times[0].year} {_clock(times[0])}"),
                   ("End Time", f"{times[-1].month}/{times[-1].day}/{times[-1].year} {_clock(times[-1])}")]
        frames = max(1, len(times) // max(n, 1))
        header.append(("Total Exposures", f"{len(times)} ({n} channels, {frames} frame{'s' * (frames != 1)})"))
    sections = []
    # Sample
    pattern = el.find(".//SamplePattern")
    dyes = []                                    # (name, fluorochrome), from the dye list (Name / Fluochrome / LutName)
    for d in el.iter():
        if d.find("Fluochrome") is not None and d.find("Name") is not None and d.find("LutName") is not None:
            dyes.append((d.findtext("Name"), d.findtext("Fluochrome")))
    if pattern is not None:
        lines = [", ".join(x for x in (pattern.findtext("ProjectName"), pattern.findtext("SampleName")) if x)]
        carrier = pattern.find("CarrierInfo")
        if carrier is not None:
            try:
                count = int(carrier.findtext("MultiCarrierCount") or 1)
                w, h = float(carrier.findtext("Width")) * 1000, float(carrier.findtext("Height")) * 1000
                cw, ch = float(carrier.findtext("CoverGlassWidth")) * 1000, float(carrier.findtext("CoverGlassHeight")) * 1000
                mask = int(carrier.findtext("MultiCarrierMask") or 0)
                positions = [str(k + 1) for k in range(count) if mask >> k & 1]
                lines.append(f"{count} x Slide {w:.0f}x{h:.0f}, Coverslip {cw:.0f}mm x {ch:.0f}mm"
                             + (f", Position{'s' * (len(positions) != 1)} {', '.join(positions)}" if positions else ""))
            except (TypeError, ValueError):
                pass
        medium, index = pattern.findtext("MountingMedium"), pattern.findtext("RefractiveIndex")
        if medium:
            lines.append(f"in {medium}" + (f", n={float(index):g}" if index else ""))
        if dyes:
            lines.append("Staining: " + ", ".join(nm if nm == fl or not fl else f"{nm} | {fl}" for nm, fl in dyes))
        sections.append(("Sample", [l for l in lines if l]))
    # Objective
    cam = next((c for c in el.iter("ATLCameraSettingDefinition") if c.get("ObjectiveName")), None)
    if cam is not None:
        obj = re.sub(r"(\d+)\s*x\b", r"\1x", re.sub(r"\s+", " ", cam.get("ObjectiveName", "")).strip())
        na, immersion = cam.get("NumericalAperture"), cam.get("Immersion")
        if na and "/" not in obj:               # (the name usually has the aperture and immersion already: "... 10x/0.32 DRY")
            obj += f"/{float(na):.2f}" + (f" {immersion}" if immersion else "")
        sections.append(("Objective", [obj]))
    # Fluorescence settings: one per channel (the first WideFieldChannelInfo of each channel number)
    infos = {}
    for info in el.iter("WideFieldChannelInfo"):
        infos.setdefault(info.get("Channel"), info)
    for k, info in enumerate(list(infos.values())[:n]):
        nm, fl = dyes[k] if k < len(dyes) else (info.get("UserDefName") or entry.channels[k].name, "")
        lines = ["Dyes:", "    " + (nm if nm == fl or not fl else f"{nm} | {fl}")]
        leds = [f"    LED {info.get(f'ILLEDWavelength{j}')}: {_dec(info.get(f'ILLEDIntensity{j}', 0))}%"
                for j in range(int(info.get("ILLEDMaxWavelengths") or info.get("ILLEDMaxUsedChannels") or 0))
                if info.get(f"ILLEDActiveState{j}") == "1"]
        if leds:
            lines += ["Illumination:"] + leds
        cams = list(info.iter("IndividualCameraInfo"))
        if cams:
            lines.append(f"Exposure {float(cams[0].get('ExposureTime', 0)) * 1000:.1f}ms")
            gains = [_dec(c.get("Gain", 0)) for c in cams]
            lines.append("Gain " + (gains[0] if len(set(gains)) == 1 else " / ".join(gains)))
        sections.append((f"Fluorescence Setting {k + 1} (WF)", lines))
    # Image detector
    if entry.width and entry.pixel_um:
        sections.append(("Image Detector", [f"Image: {entry.width}px x {entry.height}px",
                                            f"Field of view: {entry.width * entry.pixel_um / 1000:.3f}mm x {entry.height * entry.pixel_um / 1000:.3f}mm",
                                            f"Pixelsize: {entry.pixel_um * 1000:.0f}nm"]))
    if entry.annotations:
        sections.append(("Annotations", [f"{len(entry.annotations)} (lines, shapes and text drawn in LAS X or here)"]))
    stamps = [(k + 1, f"{(t - times[0]).total_seconds():.3f}", _clock(t), f"{t.month}/{t.day}/{t.year}") for k, t in enumerate(times)]
    return {"header": header, "sections": sections, "stamps": stamps}


def load_experiment(path):
    """The .xlef as a FolderEntry (its images and sub-folders)."""
    entry = load_entry(path)
    if isinstance(entry, FolderEntry):
        entry.name = entry.saved_name = os.path.splitext(os.path.basename(path))[0]     # (the project: not renamed in the viewer)
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
        self.overviews = {}      # (z, c) -> downsampled array, in the image's own pixel values
        self.ov_step = max(1, int(np.ceil(max(entry.width, entry.height) / OVERVIEW_MAX)))
        self.pyramid = self.usable_level(OVERVIEW_MAX)
        if self.pyramid is not None:
            self.ov_step = self.pyramid[1]

    def quick_level(self, size):
        """(level, step) of a smaller pyramid level, at least `size` across (the fit view's size on screen), to show first while the
        overview's own level is read; None if there is no smaller one."""
        quick = self.usable_level(size) if self.pyramid is not None else None
        return quick if quick is not None and quick[1] > self.pyramid[1] else None

    def usable_level(self, size):
        """(level, step) when a view can come from LAS X's pyramid: the smallest level at least `size` across, stored on this computer,
        8-bit, exactly 1/step of the image's size. None otherwise (no pyramid, a z stack, a level file that needs asking, ...)."""
        e = self.entry
        if e.nz != 1 or not e.pyramid:
            return None
        try:
            level = pyramid_level(e, size)
            if level is None or not level.width:
                return None
            step = round(e.width / level.width)
            if step < 2 or step & (step - 1) or e.width // step != level.width or e.height // step != level.height:
                return None
            files = [level.files.get((0, c)) for c in range(len(e.channels))]
            if not all(files) or any(needs_asking(f) for f in files):
                return None
            return level, step
        except (OSError, ET.ParseError):
            return None

    def plane(self, z, c):
        key = (z, c)
        if key not in self.planes:
            path = self.entry.files.get(key)
            if path is None:
                return None
            page = self.entry.pages.get(key)
            try:
                a = tifffile.memmap(path, mode="r", page=page)
            except (ValueError, OSError):
                a = tifffile.imread(path, key=page if page is not None else 0) if page is not None else tifffile.imread(path)
            while a.ndim > 2:        # a z stack in one file: pick the plane; extra sample axes: take the first
                a = a[z if a.shape[0] == self.entry.nz else 0]
            self.planes[key] = a
        return self.planes[key]

    def overview(self, z, c, progress=None):
        """Channel c's overview for plane z (read once). progress(fraction), if given, is called as it is read from the image itself."""
        key = (z, c)
        if key not in self.overviews:
            small = self._from_pyramid(c) if self.pyramid is not None else None
            if small is None:
                p = self.plane(z, c)
                if p is None:
                    return None
                s = self.ov_step
                small = np.empty((-(-p.shape[0] // s), -(-p.shape[1] // s)), p.dtype)
                for i in range(0, small.shape[0], 64):          # (in chunks of rows, so that the progress can be shown)
                    small[i:i + 64] = p[i * s:(i + 64) * s:s, ::s]
                    if progress is not None:
                        progress(min(1.0, (i + 64) / small.shape[0]))
            self.overviews[key] = small
        return self.overviews[key]

    def _from_pyramid(self, c, level_step=None):
        """Channel c from a pyramid level (default: the overview's; a file of a few MB), converted back to the image's own values with
        LAS X's fixed gamma (see PYRAMID_GAMMA): the image's own files are not touched. None if the level can't be used (the overview is
        then read from the image itself)."""
        level, s = level_step or self.pyramid
        try:
            small = tifffile.imread(level.files[(0, c)])
        except (OSError, ValueError, KeyError):
            return None
        while small.ndim > 2:
            small = small[0]
        if small.dtype != np.uint8:
            return None
        return pyramid_to_original(small, self.entry.channels[c].max_value)

    def has_all(self, z):
        """All channels' overviews for plane z are in (at whatever step)."""
        return all((z, c) in self.overviews for c in range(len(self.entry.channels)))

    def complete(self, z):
        """All channels' overviews for plane z are in, at the overview's own (sharp) step."""
        return (all((z, c) in self.overviews for c in range(len(self.entry.channels)))
                and (self.pyramid is None or self.ov_step == self.pyramid[1]))

    def nbytes(self):
        return sum(a.nbytes for a in self.overviews.values())

    def close(self):
        self.planes.clear()
        self.overviews.clear()


def read_rows(path, step, wait, progress, chunk=64, page=None):
    """Every step-th row and column of a TIFF (step 1: all of it), read in chunks of rows so that a background read can pause between
    them: wait() is called before each chunk and blocks while it should pause; anything but "go" from it ends the read and is returned.
    progress(fraction) is called after each chunk. Returns the array, or wait()'s answer, or "skip" for a compressed image that is too
    big to read whole in the background."""
    try:
        src = tifffile.memmap(path, mode="r", page=page)          # (page: one of an OME-TIFF's channels; None: the file's image)
    except (ValueError, OSError):
        if step != 1:
            return "skip"
        a = tifffile.imread(path, key=page or 0)            # (a small compressed pyramid level: read whole)
        while a.ndim > 2:
            a = a[0]
        progress(1.0)
        return a
    while src.ndim > 2:
        src = src[0]
    h, w = src.shape
    rows = -(-h // step)
    out = np.empty((rows, -(-w // step)), src.dtype)
    for i in range(0, rows, chunk):
        state = wait()
        if state != "go":
            return state
        out[i:i + chunk] = src[i * step:(i + chunk) * step:step, ::step]
        progress(min(1.0, (i + chunk) / rows))
    return out


def make_luts(channels, settings):
    """Per shown channel, its lookup tables: [(color component 0-2, (65536,) uint16 table: pixel value -> that component, 0-255)], only
    for the components its color uses (red uses one, yellow two), since each table lookup is the slow part of drawing. None if off."""
    luts = []
    x = np.arange(65536, dtype=np.float32)
    for ch, st in zip(channels, settings):
        if not st["on"]:
            luts.append(None)
            continue
        lo, hi = st["lo"], max(st["hi"], st["lo"] + 1)
        table, under, over = lut_table(getattr(ch, "lut", "Gray"))
        idx = np.rint(np.clip((x - lo) / (hi - lo), 0, 1) * 255).astype(np.intp)
        rgb = table[idx] * 255                          # (65536, 3): the channel's colour map through its black/white levels
        if under is not None:
            rgb[0] = under * 255                        # (LAS X's (U) maps: zero marked)
        if over is not None:
            rgb[int(ch.max_value):] = over * 255        # ((O) maps: the maximum value marked)
        luts.append([(i, rgb[:, i].astype(np.uint16)) for i in range(3) if rgb[:, i].any()])
    return luts


def composite(arrays, luts):
    """The arrays (one per channel, same shape, or None) mixed into an RGB uint8 image by their lookup tables (additive)."""
    shape = next((a.shape for a in arrays if a is not None), (1, 1))
    planes = [None, None, None]
    for a, lut in zip(arrays, luts):
        if a is None or lut is None:
            continue
        if a.dtype not in (np.uint8, np.uint16):
            a = np.clip(a, 0, 65535).astype(np.uint16)
        for i, table in lut:
            if planes[i] is None:
                planes[i] = table[a]
            else:
                planes[i] += table[a]
    out = np.zeros(shape + (3,), np.uint8)
    for i, p in enumerate(planes):
        if p is not None:
            out[..., i] = np.minimum(p, 255, out=p)
    return out


# How pixels are drawn when the view is zoomed in past them (one image pixel bigger than one screen pixel). Compared on section 17 at
# 16x: bilinear turns a bright pixel into a plus; Lanczos rings; B-spline and wide Gaussians dim small bright spots. Mitchell-Netravali
# is round, with almost no ringing.
ENLARGE_METHODS = {"mitchell": "Mitchell-Netravali (smooth, no halos; default)", "bicubic": "Bicubic (sharper, faint halos)",
                   "gaussian": "Gaussian (softer)", "bilinear": "Bilinear", "nearest": "Hard pixels (the true pixels as squares)"}


def _cubic(B, C):
    def k(x):
        x = np.abs(x)
        return np.where(x < 1, ((12 - 9 * B - 6 * C) * x ** 3 + (-18 + 12 * B + 6 * C) * x ** 2 + (6 - 2 * B)) / 6,
                        np.where(x < 2, ((-B - 6 * C) * x ** 3 + (6 * B + 30 * C) * x ** 2 + (-12 * B - 48 * C) * x
                                         + (8 * B + 24 * C)) / 6, 0))
    return k


_KERNELS = {"mitchell": _cubic(1 / 3, 1 / 3), "bicubic": _cubic(0, 0.5), "gaussian": lambda x: np.exp(-0.5 * (x / 0.5) ** 2)}


def _taps(n_in, n_out, kernel):
    """For each of n_out output pixels along one axis: the 4 input pixels it is made of and their weights (summing to 1)."""
    x = (np.arange(n_out) + 0.5) * (n_in / n_out) - 0.5            # (where each output pixel's centre falls in the input)
    idx = np.floor(x).astype(int)[:, None] + np.arange(-1, 3)[None, :]
    w = kernel(x[:, None] - idx).astype(np.float32)
    w /= w.sum(axis=1, keepdims=True)
    return np.clip(idx, 0, n_in - 1), w


@functools.lru_cache(maxsize=16)
def _weights(n_in, n_out, method):
    """The (n_out, n_in) matrix that resamples one axis (mostly zeros: 4 weights per row). Kept: redraws at one zoom (a slider, a
    channel) reuse it."""
    idx, w = _taps(n_in, n_out, _KERNELS[method])
    m = np.zeros((n_out, n_in), np.float32)
    np.add.at(m, (np.repeat(np.arange(n_out), 4), idx.ravel()), w.ravel())       # (add: at the edges taps repeat a pixel)
    return m


def enlarge(rgb, w, h, method="mitchell"):
    """An RGB uint8 array enlarged to w x h (a PIL image), the way ENLARGE_METHODS names: each output pixel a weighted sum of the 4 x 4
    nearest input pixels (cubic or Gaussian weights), done as two matrix products per color (rows, then columns), which is fast
    (~10-25 ms for a 1200 x 800 view). Bilinear and hard pixels come from Pillow."""
    if method not in _KERNELS:
        return Image.fromarray(rgb).resize((w, h), Image.NEAREST if method == "nearest" else Image.BILINEAR)
    rows, cols = _weights(rgb.shape[0], h, method), _weights(rgb.shape[1], w, method)
    a = rgb.astype(np.float32)
    out = np.empty((h, w, 3), np.uint8)
    for c in range(3):
        out[..., c] = np.clip(rows @ a[..., c] @ cols.T + 0.5, 0, 255)
    return Image.fromarray(out)


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
    """One slider with two handles: black level and white level (left/right, or bottom/top when vertical). Drag a handle to move it,
    drag the bar between them to shift both, click elsewhere on the track to bring the nearer handle there. Shift while dragging moves
    10x slower. command(lo, hi) is called on every change made with the mouse (not on set())."""

    HANDLE_W = 7
    MIN_GAP = 12             # (px: the handles are drawn at least this far apart; their values may still meet)

    def __init__(self, parent, max_value, lo, hi, color, command=None, vertical=False):
        bg = PAL["bg"]
        size = {"width": 36, "height": 200} if vertical else {"height": 36}
        super().__init__(parent, highlightthickness=0, background=bg, cursor="sb_v_double_arrow" if vertical else "sb_h_double_arrow",
                         **size)
        self.max_value, self.color, self.command, self.vertical = max_value, color, command, vertical
        self.pad = 16 if vertical else 8        # (vertical: room for the values above and below the track)
        self.active = True                      # (False: drawn grey, for a channel that is switched off; it can still be moved)
        self.lo, self.hi = lo, hi
        self.drag = None
        self.bind("<Configure>", lambda e: self._draw())
        self.bind("<ButtonPress-1>", self._press)
        self.bind("<B1-Motion>", self._motion)
        self.bind("<ButtonRelease-1>", lambda e: setattr(self, "drag", None))

    def set(self, lo, hi):
        self.lo, self.hi = int(lo), int(hi)
        self._draw()

    def set_active(self, active):
        self.active = active
        self._draw()

    def set_color(self, color):
        self.color = color
        self._draw()

    def _dimmed(self):
        """The bar's color for a channel switched off: its own color, 60% of the way to the background (so it is still told apart)."""
        (r, g, b), (R, G, B) = self.winfo_rgb(self.color), self.winfo_rgb(PAL["bg"])
        mix = lambda a, z: int((0.4 * a + 0.6 * z) / 257)
        return "#%02x%02x%02x" % (mix(r, R), mix(g, G), mix(b, B))

    def _length(self):
        return self.winfo_height() if self.vertical else self.winfo_width()

    def _span(self):
        return max(self._length() - 2 * self.pad, 1)

    def _x(self, v):
        """Where value v is along the track (0 = the black end)."""
        return self.pad + v / self.max_value * self._span()

    def _value(self, x):
        return (x - self.pad) / self._span() * self.max_value

    def _shown(self):
        """Where the two handles are drawn along the track: at their values, but never closer than MIN_GAP (about a handle's width), as
        LAS X does, so that both can be seen and grabbed even when the values are equal; pushed apart about their middle, and kept on
        the track."""
        al, ah = self._x(self.lo), self._x(self.hi)
        a0, a1 = self.pad, self.pad + self._span()
        gap = min(self.MIN_GAP, a1 - a0)
        if ah - al < gap:
            mid = min(max((al + ah) / 2, a0 + gap / 2), a1 - gap / 2)
            al, ah = mid - gap / 2, mid + gap / 2
        return al, ah

    def _along(self, event):
        """The mouse position along the track, measured from the black end."""
        return self._length() - event.y if self.vertical else event.x

    def _pt(self, a, t):
        """Canvas point at `a` along the track and `t` across it (from the track's middle line)."""
        if self.vertical:
            return self.winfo_width() / 2 + t, self._length() - a
        return a, 24 + t

    def _draw(self):
        self.delete("all")
        hw = self.HANDLE_W
        a0, a1 = self.pad, self.pad + self._span()
        al, ah = self._shown()
        on = self.active
        p = PAL
        self.configure(background=p["bg"])
        text_color = p["slider_text"] if on else p["off_fg"]
        self.create_rectangle(*self._pt(a0, -2), *self._pt(a1, 2), fill=p["track"] if on else p["off_track"], outline="")
        self.create_rectangle(*self._pt(al, -4), *self._pt(ah, 4), fill=self.color if on else self._dimmed(),
                              outline=p["handle_line"] if on else p["off_line"])
        for a, fill in ((al, "black" if on else p["off_black"]), (ah, "white" if on else p["off_white"])):
            shape = [(-hw / 2, 9), (hw / 2, 9), (hw / 2, -5), (0, -9), (-hw / 2, -5)]       # (a pointer, across the track)
            self.create_polygon(*[c for da, t in shape for c in self._pt(a + da, t)], fill=fill,
                                outline=p["handle_line"] if on else p["off_line"])
        if self.vertical:      # (just the numbers: the white handle's above, the black handle's below)
            w, h = self.winfo_width(), self._length()
            self.create_text(w / 2, 1, anchor=tk.N, text=str(self.hi), fill=text_color, font=("Segoe UI", 8))
            self.create_text(w / 2, h - 1, anchor=tk.S, text=str(self.lo), fill=text_color, font=("Segoe UI", 8))
        else:
            self.create_text(a0, 1, anchor=tk.NW, text=f"black {self.lo}", fill=text_color)
            self.create_text(a1, 1, anchor=tk.NE, text=f"white {self.hi}", fill=text_color)

    def _press(self, event):
        x = self._along(event)
        self.before_press = (self.lo, self.hi)
        self.presses = getattr(self, "presses", 0) + 1         # (the viewer makes each press one step of Undo)
        xl, xh = self._shown()          # (the handles where they are drawn, which is where they are clicked)
        reach = self.HANDLE_W
        if min(abs(x - xl), abs(x - xh)) <= reach:
            which = "hi" if abs(x - xh) < abs(x - xl) else "lo"
        elif xl < x < xh:
            which = "both"
        else:      # on the track outside the range: bring the nearer handle here
            which = "lo" if x < xl else "hi"
            self._apply(which, round(self._value(x)) - getattr(self, which))
        self.drag = (which, x, self.lo, self.hi)

    def _motion(self, event):
        if not self.drag:
            return
        which, x0, lo0, hi0 = self.drag
        dv = (self._along(event) - x0) / self._span() * self.max_value * (0.1 if event.state & 0x0001 else 1)
        self._apply(which, round(dv), lo0, hi0)

    def undo_press(self):
        """Put back the handles as they were before the last click (whose first press, of a double click, may have moved one)."""
        lo, hi = getattr(self, "before_press", (self.lo, self.hi))
        self.drag = None
        self._apply("lo", lo - self.lo)
        self._apply("hi", hi - self.hi)

    def _apply(self, which, dv, lo=None, hi=None):
        """Move handle(s) `which` by dv from (lo, hi) (default: where they are now)."""
        lo, hi, top = self.lo if lo is None else lo, self.hi if hi is None else hi, self.max_value
        if which == "both":
            dv = min(max(dv, -lo), top - hi)
            lo, hi = lo + dv, hi + dv
        elif which == "lo":
            lo = min(max(lo + dv, 0), hi)       # (black and white may meet: the image is then two-level)
        else:
            hi = max(min(hi + dv, top), lo)
        if (lo, hi) != (self.lo, self.hi):
            self.lo, self.hi = lo, hi
            self._draw()
            if self.command:
                self.command(lo, hi)


def _thumb_cache_path(entry, size):
    """Where the thumbnail of `entry` is cached: named by its files' paths, sizes and dates (os.stat only: no download)."""
    h = hashlib.sha1(f"{size}|lasx-contrast|{'|'.join(f'{c.lut}:{c.black}-{c.white}' for c in entry.channels)}".encode())     # (the second part changes whenever thumbnails are drawn differently)
    for f in entry.all_files():
        st = os.stat(f)
        h.update(f"{os.path.normcase(f)}|{st.st_size}|{st.st_mtime_ns}".encode())
    return os.path.join(THUMB_CACHE, h.hexdigest() + ".png")


PYRAMID_GAMMA = 1.8      # LAS X's 8-bit pyramid levels: level = 255 * (original / max) ** (1 / 1.8) (measured 2026-10-29: fit 1.804,
                         # median error 4.6 of 4095 over 6 sections x 4 channels, LAS X 1.15)


def pyramid_to_original(small, max_value):
    """A pyramid level's 8-bit values as the image's own values (0..max_value), undoing LAS X's gamma."""
    lut = np.round(max_value * (np.arange(256) / 255) ** PYRAMID_GAMMA).astype(np.uint16)
    return lut[small]


_PYRAMID_LEVELS = {}      # normcase(.xlif path) -> its ImageEntry (None: no such file)


def _pyramid_entry(path):
    """A pyramid level's .xlif, read once and then remembered: opening an image, its quick picture, its thumbnail and the background
    loader each walk the levels, and on a network drive each read is a round trip to the server. (LAS X's levels don't change.)"""
    key = os.path.normcase(path)
    if key not in _PYRAMID_LEVELS:
        _PYRAMID_LEVELS[key] = load_entry(path) if os.path.isfile(path) else None
    return _PYRAMID_LEVELS[key]


def pyramid_level(entry, size):
    """The smallest level of LAS X's pyramid (reduced 8-bit copies, _pmd_0 = 1/2 size, _pmd_1 = 1/4, ...) that is still at least `size`
    on its longest side, as an ImageEntry, or None if there is no pyramid. Reads only the levels' small .xlif files."""
    best, path = None, entry.pyramid
    for _ in range(16):
        if not path:
            break
        level = _pyramid_entry(path)
        if not isinstance(level, ImageEntry) or level.problem or len(level.channels) != len(entry.channels):
            break
        if best is not None and max(level.width, level.height) < size:
            break
        best, path = level, level.pyramid
    return best


def _thumb_source(entry, size):
    """What to make the thumbnail from: the pyramid level for `size` when its files can be read without asking (on this computer, or
    online-only but small: see needs_asking), else the image itself if its files can; (None, "cloud") when neither can."""
    for source in (pyramid_level(entry, size), entry):
        if source is None or not source.files:
            continue
        files = [source.files.get((0, c)) for c in range(len(source.channels))]
        if not any(f and needs_asking(f) for f in files):
            return source, files
    return None, "cloud"


def make_thumbnail(entry, size=THUMB_SIZE):
    """(PIL image, "ok"), or (None, "cloud") when the files are online-only (they are NOT read), or (None, "none").
    Made from the smallest fitting level of LAS X's pyramid if there is one (a few KB per channel), else from every n-th row and column
    of the image itself, and shown with the image's contrast (LAS X's, or as changed in the viewer), as the main view first shows it."""
    if entry.problem or not entry.files or not entry.width:
        return None, "none"
    try:
        cache = _thumb_cache_path(entry, size)
        if os.path.isfile(cache):
            return Image.open(cache).convert("RGB"), "ok"
        source, files = _thumb_source(entry, size)
    except (OSError, ET.ParseError):
        return None, "none"
    if source is None:
        return None, files
    step = max(1, max(source.width, source.height) // size)       # (rounded down: the resize below makes it exactly `size`)
    arrays = []
    for c, f in enumerate(files):
        if f is None:
            arrays.append(None)
            continue
        page = source.pages.get((0, c))                         # (one of an OME-TIFF's pages, or None)
        try:
            a = tifffile.memmap(f, mode="r", page=page)
        except (ValueError, OSError):
            if os.path.getsize(f) > THUMB_MAX_READ:
                return None, "none"
            a = tifffile.imread(f, key=page or 0)
        while a.ndim > 2:
            a = a[0]
        small = np.ascontiguousarray(a[::step, ::step])
        del a
        if source is not entry and small.dtype == np.uint8:      # a pyramid level: back to the image's own values
            small = pyramid_to_original(small, entry.channels[len(arrays)].max_value)
        arrays.append(small)
    luts = make_luts(entry.channels, [{"on": True, "lo": c.black, "hi": c.white} for c in entry.channels])
    im = Image.fromarray(composite(arrays, luts))
    im.thumbnail((size, size), Image.BILINEAR)
    try:
        os.makedirs(THUMB_CACHE, exist_ok=True)
        im.save(cache)
    except OSError:
        pass
    return im, "ok"


SLIDE_HEIGHT = 260       # the slide overview's height above the list, screen pixels


class SlideOverview(ttk.Frame):
    """The slide (LAS X's "Sample Overview" image) with the outline of each image of the experiment on it, the open one highlighted;
    clicking an outline calls on_click(item), and resting the pointer on one shows name_of(item) beside it. Images are placed by their
    stage positions: an image's position and the overview's are both of their CENTRES, in stage um (worked out on Tyler's BF567 slide:
    the outlines then sit on the tissue)."""

    def __init__(self, parent, on_click, name_of=str):
        super().__init__(parent)
        self.on_click, self.name_of = on_click, name_of
        self.tip = None           # (the name's note while the pointer is on an outline: [window, label, item])
        self.canvas = tk.Canvas(self, height=SLIDE_HEIGHT, background="#202020", highlightthickness=0)
        self.canvas.keep_colors = True        # (dark in both themes, like LAS X's)
        self.canvas.pack(fill=tk.X)
        self.extent = None        # (left, top, width, height) of the overview on the stage, um
        self.boxes = []           # (item, (x0, y0, x1, y1) stage um)
        self.picture = None       # PIL image of the slide
        self.photo = None
        self.current = None
        self.canvas.bind("<Configure>", lambda e: self.draw())
        self.canvas.bind("<ButtonPress-1>", self._click)
        self.canvas.bind("<Motion>", self._hover)
        self.canvas.bind("<Leave>", lambda e: self._hide_tip())

    def set_slide(self, extent, boxes):
        self.extent, self.boxes, self.picture, self.current = extent, boxes, None, None
        self._hide_tip()
        self.draw()

    def set_picture(self, image):
        self.picture = image
        self.draw()

    def set_current(self, item):
        self.current = item
        self.draw()

    def _transform(self):
        """(scale, offset x, offset y): canvas = offset + scale * (stage um - overview corner)."""
        cw, ch = max(self.canvas.winfo_width(), 1), max(self.canvas.winfo_height(), 1)
        _, _, w, h = self.extent
        k = min((cw - 8) / w, (ch - 8) / h)
        return k, (cw - k * w) / 2, (ch - k * h) / 2

    def draw(self):
        c = self.canvas
        c.delete("all")
        if not self.extent:
            return
        k, ox, oy = self._transform()
        left, top, w, h = self.extent
        if self.picture is not None:
            size = (max(1, round(k * w)), max(1, round(k * h)))
            self.photo = ImageTk.PhotoImage(self.picture.resize(size, Image.BILINEAR), master=c)
            c.create_image(ox, oy, anchor=tk.NW, image=self.photo)
        else:
            c.create_rectangle(ox, oy, ox + k * w, oy + k * h, outline="#666")
            c.create_text(ox + k * w / 2, oy + k * h / 2, text="slide\n(loading or\nonline-only)", fill="#888", justify=tk.CENTER)
        for item, (x0, y0, x1, y1) in self.boxes:          # (the open one last, on top)
            if item != self.current:
                c.create_rectangle(ox + k * (x0 - left), oy + k * (y0 - top), ox + k * (x1 - left), oy + k * (y1 - top),
                                   outline="#e0b030", width=1, tags=("box", f"i:{item}"))
        for item, (x0, y0, x1, y1) in self.boxes:
            if item == self.current:
                c.create_rectangle(ox + k * (x0 - left), oy + k * (y0 - top), ox + k * (x1 - left), oy + k * (y1 - top),
                                   outline="#ff3030", width=3, tags=("box", f"i:{item}"))

    def _item_at(self, event):
        """The image whose outline the pointer is in (overlapping outlines: the one whose centre is nearest), or None."""
        if not self.extent:
            return None
        k, ox, oy = self._transform()
        left, top, _, _ = self.extent
        x, y = left + (event.x - ox) / k, top + (event.y - oy) / k
        hits = [(abs((x0 + x1) / 2 - x) + abs((y0 + y1) / 2 - y), item) for item, (x0, y0, x1, y1) in self.boxes
                if x0 <= x <= x1 and y0 <= y <= y1]
        return min(hits)[1] if hits else None

    def _click(self, event):
        item = self._item_at(event)
        if item is not None:
            self._hide_tip()
            self.on_click(item)

    def _hover(self, event):
        """The pointer on an outline: the section's name beside the pointer (and a hand: a click opens it)."""
        item = self._item_at(event)
        self.canvas.configure(cursor="hand2" if item is not None else "")
        if item is None:
            self._hide_tip()
            return
        if self.tip is None:
            self.tip = PointerTip(self.canvas)
        self.tip.show(self.name_of(item), event.x_root, event.y_root)

    def _hide_tip(self):
        if self.tip is not None:
            self.tip.hide()
            self.tip = None


class PointerTip:
    """A small note beside the pointer that follows it (for things that aren't widgets of their own: outlines on the slide, a cloud
    after a name in the list); show() moves it and sets its text, hide() removes it."""

    def __init__(self, widget):
        self.win = tk.Toplevel(widget)
        self.win.wm_overrideredirect(True)
        self.label = tk.Label(self.win, background="#ffffe0", foreground="black", relief=tk.SOLID, borderwidth=1, padx=4, pady=2,
                              justify=tk.LEFT, wraplength=360)
        self.label.pack()
        self.text = None

    def show(self, text, x_root, y_root):
        if text != self.text:
            self.label.configure(text=text)
            self.text = text
        place_on_monitor(self.win, x_root + 14, y_root + 16, avoid_y=y_root)        # (kept on the pointer's monitor)

    def hide(self):
        self.win.destroy()


def slide_layout(entries):
    """For the slide overview: (overview entry, its extent (left, top, width, height) in stage um, [(item, box in stage um)]), or None
    when the experiment has no overview image with a stage position. `entries` is {item: ImageEntry}."""
    def extent(e):
        w, h = e.width * e.pixel_um, e.height * e.pixel_um
        return e.stage_um[0] - w / 2, e.stage_um[1] - h / 2, w, h
    placed = {i: e for i, e in entries.items() if e.stage_um and e.width and e.pixel_um}
    overviews = [e for e in placed.values() if "overview" in e.name.lower()]
    if not overviews:
        return None
    ov = max(overviews, key=lambda e: e.width * e.height * e.pixel_um ** 2)
    boxes = []
    for item, e in placed.items():
        if e is not ov:
            x0, y0, w, h = extent(e)
            boxes.append((item, (x0, y0, x0 + w, y0 + h)))
    return ov, extent(ov), boxes


class ThumbGrid(ttk.Frame):
    """The images as a grid of thumbnails (folders as headings), in the order of the list. Clicking one calls on_click(item); the arrow keys
    move through them. Thumbnails arrive later, with set_thumb()."""

    GAP, NAME_H = 8, 30

    def __init__(self, parent, on_click, size=THUMB_SIZE, on_menu=None, on_drop=None, name_tip=None):
        super().__init__(parent)
        self.on_click, self.size, self.on_menu, self.on_drop = on_click, size, on_menu, on_drop
        self.name_tip, self.tip = name_tip, None     # (name_tip(item): the note to show while the pointer is on its name, or None)
        self.press = None         # (a thumbnail pressed: (item, x, y, dragging yet); dropped on another: on_drop(item, target, after))
        self.canvas = tk.Canvas(self, highlightthickness=0, background=PAL["grid"], takefocus=True)
        sb = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=sb.set)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.rows = []            # ("project", name, its list item), ("folder", name) or ("image", item, name, tag)
        self.thumbs = {}          # item -> (PhotoImage or None, state)
        self.dimmed = set()       # items whose names are dimmed (nothing of the image in memory yet: see set_dimmed)
        self.selected = None
        self.cols = 1
        self.tile_y = {}          # item -> y of its tile
        self.canvas.bind("<Configure>", lambda e: self.layout())
        self.canvas.bind("<ButtonPress-1>", self._click)
        self.canvas.bind("<B1-Motion>", self._drag)
        self.canvas.bind("<ButtonRelease-1>", self._release)
        self.canvas.bind("<Button-3>", self._menu)
        self.canvas.bind("<Motion>", self._motion)
        self.canvas.bind("<Leave>", lambda e: self._hide_tip())
        self.canvas.bind("<MouseWheel>", lambda e: self.canvas.yview_scroll(int(-e.delta / 120), "units"))
        for key, d in (("Left", -1), ("Right", 1), ("Up", "-row"), ("Down", "row")):
            self.canvas.bind(f"<{key}>", lambda e, d=d: self._step(d))

    def _name_color(self, item, tag):
        dim = item in self.dimmed
        return (PAL["bad"] if tag == "bad" else PAL["link_dim" if dim else "link"] if tag == "cloud"
                else PAL["dim"] if dim else PAL["tile_text"])

    def set_dimmed(self, item, dim):
        """Dim an image's name (as in the list, while nothing of it is in memory), or not."""
        if (item in self.dimmed) == dim:
            return
        (self.dimmed.add if dim else self.dimmed.discard)(item)
        tag = next((r[3] for r in self.rows if r[0] == "image" and r[1] == item), None)
        if tag is not None:
            self.canvas.itemconfigure(f"name:{item}", fill=self._name_color(item, tag))

    def _motion(self, event):
        """The pointer on an image's name: its note (see name_tip) beside the pointer, if it has one."""
        c = self.canvas
        x, y = c.canvasx(event.x), c.canvasy(event.y)
        text = None
        if self.name_tip is not None and self.press is None:
            for obj in c.find_overlapping(x - 1, y - 1, x + 1, y + 1):
                name = next((t[5:] for t in c.gettags(obj) if t.startswith("name:")), None)
                if name is not None:
                    text = self.name_tip(name)
                    break
        if text is None:
            self._hide_tip()
            return
        if self.tip is None:
            self.tip = PointerTip(c)
        self.tip.show(text, event.x_root, event.y_root)

    def _hide_tip(self):
        if self.tip is not None:
            self.tip.hide()
            self.tip = None

    def set_items(self, rows, message="", keep=False):
        """The rows to show; `message` is shown instead when there are none (e.g. "Loading..."). keep: the thumbnails (and selection)
        of images still among them stay."""
        items = {r[1] for r in rows if r[0] == "image"}
        thumbs = {i: t for i, t in self.thumbs.items() if i in items} if keep else {}
        selected = self.selected if keep and self.selected in items else None
        self.rows, self.thumbs, self.selected, self.message, self.progress = rows, thumbs, selected, message, {}
        if not keep:
            self.canvas.yview_moveto(0)
        self.layout()

    def set_progress(self, item, fraction):
        """A thin blue bar along the bottom of the item's tile, while it loads in the background (fraction None: no bar)."""
        if fraction is None:
            self.progress.pop(item, None)
        else:
            self.progress[item] = fraction
        self._draw_progress(item)

    def _draw_progress(self, item):
        c = self.canvas
        c.delete(f"prog:{item}")
        frac = getattr(self, "progress", {}).get(item)
        box = c.coords(f"bg:{item}")
        if frac is None or not box:
            return
        x0, y0, x1, y1 = box
        c.create_rectangle(x0 + 4, y1 - 6, x1 - 4, y1 - 3, fill=PAL["placeholder"], outline="", tags=(f"prog:{item}",))
        c.create_rectangle(x0 + 4, y1 - 6, x0 + 4 + (x1 - x0 - 8) * frac, y1 - 3, fill="#3a8fff", outline="", tags=(f"prog:{item}",))

    def set_thumb(self, item, image, state):
        photo = ImageTk.PhotoImage(image, master=self.canvas) if image is not None else None
        self.thumbs[item] = (photo, state)
        self._draw_picture(item)

    def set_selected(self, item):
        if item == self.selected:
            return
        old, self.selected = self.selected, item
        for it in (old, item):
            if it is not None:
                self.canvas.itemconfigure(f"bg:{it}", fill=self._bg(it), outline=self._outline(it))
        self.see(item)

    def see(self, item):
        y = self.tile_y.get(item)
        total = float(self.canvas.cget("scrollregion").split()[3]) if self.canvas.cget("scrollregion") else 0
        if y is None or total <= 0:
            return
        top, bottom = (f * total for f in self.canvas.yview())
        tile_h = self.size + self.NAME_H + self.GAP
        if y < top or y + tile_h > bottom:
            self.canvas.yview_moveto(max(0, y - self.GAP) / total)

    def _bg(self, item):
        return PAL["tile_sel"] if item == self.selected else PAL["tile"]

    def _outline(self, item):
        return PAL["selected"] if item == self.selected else PAL["tile_line"]

    def layout(self):
        c = self.canvas
        c.delete("all")
        c.configure(background=PAL["grid"])
        self.tile_y.clear()
        if not self.rows and getattr(self, "message", ""):
            c.create_text(max(c.winfo_width(), 1) / 2, 40, text=self.message, fill=PAL["placeholder_fg"], font=("Segoe UI", 11))
            c.configure(scrollregion=(0, 0, 1, 1))
            return
        width = max(c.winfo_width(), 1)
        tile_w = self.size + self.GAP
        self.cols = max(1, (width - self.GAP) // (tile_w + self.GAP))
        x0 = self.GAP + max(0, (width - self.GAP - self.cols * (tile_w + self.GAP)) // 2)
        y, col = self.GAP, 0
        tile_h = self.size + self.NAME_H
        for row in self.rows:
            if row[0] == "project":                 # (a project: a line across, then its name, larger)
                if col:
                    y, col = y + tile_h + self.GAP, 0
                if y > self.GAP:
                    y += 2 * self.GAP
                c.create_line(self.GAP, y, width - self.GAP, y, fill=PAL["selected"], width=3)
                t = c.create_text(self.GAP, y + 6, anchor=tk.NW, text=row[1], font=("Segoe UI", 12, "bold"), fill=PAL["tile_text"],
                                  width=max(width - 2 * self.GAP, 40), tags=(f"proj:{row[2]}",))
                y = c.bbox(t)[3] + 8
                continue
            if row[0] == "folder":
                if col:
                    y, col = y + tile_h + self.GAP, 0
                t = c.create_text(self.GAP, y, anchor=tk.NW, text=row[1], font=("Segoe UI", 10, "bold"), fill=PAL["tile_text"],
                                  width=max(width - 2 * self.GAP, 40))            # (a long path wraps)
                y = c.bbox(t)[3] + 6
                continue
            _, item, name, tag = row
            x = x0 + col * (tile_w + self.GAP)
            self.tile_y[item] = y
            tags = ("tile", f"i:{item}")
            c.create_rectangle(x, y, x + tile_w, y + tile_h, fill=self._bg(item), outline=self._outline(item), width=2,
                               tags=tags + (f"bg:{item}",))
            short = (name if len(name) <= 34 else name[:32] + "…") + ("  ☁" if tag == "cloud" else "")
            c.create_text(x + tile_w / 2, y + self.size + 4, anchor=tk.N, text=short, width=tile_w - 4, justify=tk.CENTER,
                          fill=self._name_color(item, tag), tags=tags + (f"name:{item}",))
            self._draw_picture(item, x, y)
            self._draw_progress(item)
            col += 1
            if col == self.cols:
                y, col = y + tile_h + self.GAP, 0
        if col:
            y += tile_h + self.GAP
        c.configure(scrollregion=(0, 0, width, y))

    def _draw_picture(self, item, x=None, y=None):
        c = self.canvas
        if x is None:
            box = c.coords(f"bg:{item}")
            if not box:
                return
            x, y = box[0], box[1]
        c.delete(f"pic:{item}")
        tags = ("tile", f"i:{item}", f"pic:{item}")
        cx, cy = x + (self.size + self.GAP) / 2, y + 4 + self.size / 2
        photo, state = self.thumbs.get(item, (None, None))
        if photo is not None:
            c.create_image(cx, cy, image=photo, tags=tags)
        else:
            text = {"cloud": "☁\nonline-only", "none": "no preview"}.get(state, "…")
            c.create_rectangle(cx - self.size / 2 + 4, cy - self.size / 2 + 4, cx + self.size / 2 - 4, cy + self.size / 2 - 4,
                               fill=PAL["placeholder"], outline="", tags=tags)
            c.create_text(cx, cy, text=text, justify=tk.CENTER, fill=PAL["link"] if state == "cloud" else PAL["placeholder_fg"], tags=tags)

    def _item_at(self, event):
        for obj in self.canvas.find_overlapping(self.canvas.canvasx(event.x), self.canvas.canvasy(event.y),
                                                self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)):
            for t in self.canvas.gettags(obj):
                if t.startswith("i:"):
                    return t[2:]
        return None

    def _click(self, event):
        self._hide_tip()
        self.canvas.focus_set()
        item = self._item_at(event)
        self.press = (item, event.x, event.y, False) if item is not None else None
        if item is not None:
            self.on_click(item)

    def _drop_at(self, event, dragged):
        """The thumbnail under the pointer while dragging, and whether the pointer is on its right half: (target, after, its box)."""
        target = self._item_at(event)
        box = self.canvas.coords(f"bg:{target}") if target is not None and target != dragged else None
        if not box:
            return None
        return target, self.canvas.canvasx(event.x) > (box[0] + box[2]) / 2, box

    def _drag(self, event):
        """Dragging a thumbnail (after a few pixels): a bar before or after the thumbnail it would be dropped next to."""
        if self.press is None or self.on_drop is None:
            return
        item, x0, y0, started = self.press
        if not started and abs(event.x - x0) + abs(event.y - y0) < 8:
            return
        self.press = (item, x0, y0, True)
        self.canvas.configure(cursor="fleur")
        self.canvas.delete("dropline")
        if event.y < 16:
            self.canvas.yview_scroll(-1, "units")
        elif event.y > self.canvas.winfo_height() - 16:
            self.canvas.yview_scroll(1, "units")
        hit = self._drop_at(event, item)
        if hit:
            target, after, box = hit
            x = box[2] + self.GAP / 2 if after else box[0] - self.GAP / 2
            self.canvas.create_line(x, box[1], x, box[3], fill=PAL["selected"], width=3, tags="dropline")

    def _release(self, event):
        press, self.press = self.press, None
        if press is None or not press[3]:
            return
        self.canvas.delete("dropline")
        self.canvas.configure(cursor="")
        hit = self._drop_at(event, press[0])
        if hit:
            self.on_drop(press[0], hit[0], hit[1])

    def _menu(self, event):
        """Right-click on a thumbnail, or on a project's name: on_menu(item, x_root, y_root) (the viewer's menu for that section or
        project, as in the list)."""
        self._hide_tip()
        x, y = self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)
        for obj in self.canvas.find_overlapping(x - 1, y - 1, x + 1, y + 1):
            pid = next((t[5:] for t in self.canvas.gettags(obj) if t.startswith("proj:")), None)
            if pid is not None and self.on_menu is not None:
                self.on_menu(pid, event.x_root, event.y_root)
                return
        item = self._item_at(event)
        if item is not None and self.on_menu is not None:
            box = self.canvas.coords(f"bg:{item}")            # (the menu's top at the middle of the thumbnail's picture)
            y = event.y_root
            if box:
                picture_mid = box[1] + 4 + self.size / 2
                y = int(self.canvas.winfo_rooty() + picture_mid - self.canvas.canvasy(0))
            self.on_menu(item, event.x_root, y)

    def _step(self, d):
        items = [r[1] for r in self.rows if r[0] == "image"]
        if not items:
            return
        d = -self.cols if d == "-row" else self.cols if d == "row" else d
        i = items.index(self.selected) + d if self.selected in items else 0
        self.on_click(items[min(max(i, 0), len(items) - 1)])


class XlefViewer:
    def __init__(self, path=None):
        try:      # with tkinterdnd2, .xlef files can be dragged onto the window and into the Open dialog
            from tkinterdnd2 import DND_FILES, TkinterDnD
            self.root = TkinterDnD.Tk()
            self.dnd_files = DND_FILES
        except Exception:
            self.root = tk.Tk()
            self.dnd_files = None
        self.root.title("XLEF viewer")
        self.settings = self._load_settings()
        if "recent_files" not in self.settings:       # (settings written before there was a recent list)
            self.settings["recent_files"] = [self.settings["last_file"]] if self.settings.get("last_file") else []
        self.settings.pop("last_file", None)
        self.native_theme = ttk.Style(self.root).theme_use()        # (Windows' own: the light theme's)
        self.theme = self.settings.get("theme", "dark")
        if self.theme not in THEMES:
            self.theme = "dark"
        style_theme(self.root, self.theme, self.native_theme)
        self.root.configure(background=PAL["bg"])
        self.root.geometry(on_screen_geometry(self.root, self.settings.get("geometry", "1400x900")))
        self.normal_geometry = None  # (the window's size and position when last neither maximized nor minimized: what is saved)
        self.root.bind("<Configure>", self._window_moved, add="+")
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
        self.shrunk = {}             # (view..., channel) -> that channel's overview shrunk to the screen, for the view on screen
        self.projects = {}           # tree item of each open experiment (a top-level branch of the list) -> {"path", "name"}
        self.slide_project = None    # (the project whose slide overview shows)
        self.download_gen = 0        # (bumped when a project is closed: its downloads stop)
        self.kept_images = collections.OrderedDict()     # ImageEntry -> OpenImage, for images seen recently (see _keep_image)
        self.kept_bytes = 0          # (their size, kept up to date here, for the background loader's thread to read)
        self.full_res_ok = None      # (see show_entry)
        self.contrast_step, self.contrast_job = None, None     # (the slider drag being made: one Undo step; see _set_lo_hi)
        self.history = []            # changes that Undo goes back through: ("ann", entry), ("lut", [(entry, channel, old map)]),
                                     # ("range", [(entry, channel, old black, old white)]),
                                     # ("move", (item, old parent, old index)) or ("delete", (item, entry, old parent, old index))
        self.home = {}               # list item -> the project it was read from (an image dragged to another moves there on Save)
        self.deleted = {}            # list item -> ImageEntry, for images deleted and not saved yet (out of the list; Save recycles them)
        self.drag_item = None        # (a name pressed in the list: (item, y, dragging yet); see _tree_drag)
        self.drop_line = None
        self.enlarge_method = self.settings.get("enlarge_method", "mitchell")     # (zoomed-in pixels: see ENLARGE_METHODS)
        if self.enlarge_method not in ENLARGE_METHODS:
            self.enlarge_method = "mitchell"
        self.asking_download = False         # (the download question is open: see _full_res_allowed)
        self.moving_until, self.rest_job = 0.0, None     # (see _moving)
        self.folders = {}            # list item -> FolderEntry (the list's folders; images are in self.entries)
        self.item_cloud = set()      # list items of images with online-only files (a cloud after the name)
        self.rename_job = self.rename_box = None
        self.preload_progress = {}   # list item -> percent, while the background loader reads it
        self.open_progress = None    # (list item, percent) while the preview of the image opened is read
        self.preload_full = []       # images the background loader left out because the memory for previews was full
        self.thumbs_busy = False     # (the thumbnail thread is running: the background loader waits for it while they are shown)
        self.thumbs_shown = False    # (the thumbnail grid is showing, not the list: the thumbnail thread pauses otherwise)
        self.detail_cache = None     # the sharp view's pixels for the view on screen: {"key", "arrays"} (see redraw)
        self.download_queue = []     # list items whose full resolution is to be downloaded from OneDrive (see download_full_res)
        self.download_progress = {}  # list item -> "waiting" or percent, while it is queued or downloading
        self.download_thread = self.download_stop = None
        self.canvas_message = ""
        self.thumbs_go = threading.Event()      # cleared while a chosen image loads: the thumbnail thread waits, so the image comes first
        self.thumbs_go.set()
        self._build()
        self.root.protocol("WM_DELETE_WINDOW", self.quit)
        self.root.after(50, self._poll_results)
        if path:
            self.root.after(10, lambda: self.open_xlef(path))
        elif self.settings["recent_files"]:
            self.root.after(200, self.ask_open)

    # --- layout

    def _build(self):
        menu = tk.Menu(self.root)
        fm = tk.Menu(menu, tearoff=0)
        fm.add_command(label="Open (recent or browse)...", accelerator="Ctrl+O", command=self.ask_open)
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

        # Left, top to bottom: the file's name, the slide overview (when the experiment has one), a bar (Open / Save on the left,
        # List / Thumbnails on the right), the list.
        self.file_label = ttk.Label(left, text="No file open", wraplength=300, style="Muted.TLabel")
        self.file_label.pack(fill=tk.X, padx=4, pady=(4, 0))
        self.slide = SlideOverview(left, on_click=self.select_item,       # (packed when the experiment has a slide overview)
                                   name_of=lambda item: self.entries[item].name if item in self.entries else "")
        self.slide_gen = 0
        self.list_bar = ttk.Frame(left)
        self.list_bar.pack(fill=tk.X, padx=4, pady=(4, 0))
        ttk.Button(self.list_bar, text="Open...", command=self.ask_open).pack(side=tk.LEFT)
        self.save_button = ttk.Button(self.list_bar, text="Save", command=self.save_annotations, state=tk.DISABLED)
        self.save_button.pack(side=tk.LEFT, padx=(4, 0))
        Tooltip(self.save_button, "Write all unsaved changes (annotations, color maps, contrast, renames, the list's order, deleted images: "
                                  "marked with a red * in the list) into LAS X's files, "
                                  "keeping a copy of each file as it was. Nothing is saved until you click this.")
        # The List / Thumbnails toggle: real buttons, the one showing filled in (as the Open dialog's toggle), since ttk's flat Toolbutton
        # style looks like a label on Windows.
        self.view_mode = tk.StringVar(value=self.settings.get("view_mode", "list"))
        self.view_buttons = {}
        for value, text, tip in (("thumbs", "Thumbnails", "Show the images as a grid of thumbnails"),
                                 ("list", "List", "Show the images as a list of names")):      # (packed from the right)
            b = tk.Button(self.list_bar, text=text, bd=2, padx=10, pady=1, cursor="hand2",
                          command=lambda v=value: (self.view_mode.set(v), self.show_view_mode()))
            b.pack(side=tk.RIGHT, padx=(3, 0))
            Tooltip(b, tip)
            self.view_buttons[value] = b

        tree_frame = ttk.Frame(left)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.list_frame = ttk.Frame(tree_frame)
        self.tree = ttk.Treeview(self.list_frame, show="tree", selectmode="browse", columns=("size",))
        self.tree.column("#0", stretch=True)
        self.tree.column("size", width=72, minwidth=60, stretch=False, anchor=tk.E)   # (full resolution's size: see _show_sizes)
        sb = ttk.Scrollbar(self.list_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=lambda *a: (sb.set(*a), self._place_list_bars()))
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.list_bars = {}          # list item -> the thin progress bar laid over the bottom of its row (see _place_list_bar)
        self.unsaved_marks = {}      # list item -> the red * laid over its row after its name (see _place_unsaved_mark)
        self.tree.bind("<Configure>", lambda e: self._place_list_bars(), add="+")
        for ev in ("<<TreeviewOpen>>", "<<TreeviewClose>>"):
            self.tree.bind(ev, lambda e: self.root.after_idle(self._place_list_bars), add="+")
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.bind("<<TreeviewSelect>>", self.on_select)
        self.tree.bind("<<TreeviewSelect>>", lambda e: self.root.after_idle(self._place_list_bars), add="+")   # (marks' colors)
        self.tree.bind("<ButtonPress-1>", self._tree_press, add="+")
        self.tree.bind("<B1-Motion>", self._tree_drag)
        self.tree.bind("<ButtonRelease-1>", self._tree_drop, add="+")
        self.tree.bind("<Motion>", self._tree_motion, add="+")
        self.tree.bind("<Leave>", lambda e: self._hide_tree_tip(), add="+")
        self.tree.bind("<Button-3>", self._tree_menu, add="+")
        self.tree_tip = None
        self.tree.bind("<Double-Button-1>", lambda e: self._cancel_rename_job(), add="+")
        self.thumb_grid = ThumbGrid(tree_frame, on_click=self.select_item, on_menu=lambda item, x, y: self._item_menu(item, x, y),
                                    on_drop=self._grid_drop,
                                    name_tip=lambda item: self._cloud_tip_text(self.entries[item]) if item in self.item_cloud else None)
        self.thumb_gen = 0           # bumped for each experiment opened: thumbnails made for an earlier one are dropped
        self.thumbs_started = False
        self._tree_tag_colors()
        ttk.Label(left, text="☁ = online-only in OneDrive (asks before downloading files over 50 MB)", style="Link.TLabel").pack(anchor=tk.W, padx=4)
        self.show_view_mode()

        # Above the image: the annotation tools (as LAS X's). Nothing is written to LAS X's files until Save is clicked.
        bar = ttk.Frame(right)
        bar.pack(fill=tk.X, padx=4, pady=(4, 0))
        self.show_annotations = tk.BooleanVar(value=self.settings.get("show_annotations", True))
        self.annotations_check = ttk.Checkbutton(bar, text="Show annotations", variable=self.show_annotations,
                                                 command=self.on_annotations_toggle)
        self.annotations_check.pack(side=tk.LEFT, padx=(0, 8))
        Tooltip(self.annotations_check, "Show the lines, shapes and text drawn on this image in LAS X (and here). The tools after it "
                                        "draw new ones.")
        self.tool_buttons = {}
        for tool, text, tip in (("select", "Select", "Select an annotation to move it (drag it), resize it (its corners), rotate it (the "
                                 "round handle) or delete it. Double-click a text to change its words or color. Dragging elsewhere moves the view."),
                                ("rectangle", "Rectangle", "Draw a rectangle: drag from one corner to the other"),
                                ("ellipse", "Ellipse", "Draw an ellipse or circle: drag across it"),
                                ("measure", "Measure", "Draw a measuring line (its length is shown): drag from one end to the other"),
                                ("line", "Line", "Draw a line: drag from one end to the other"),
                                ("arrow", "Arrow", "Draw an arrow: drag from its tail to its head"),
                                ("text", "Text", "Add text: click where its top-left corner goes")):
            b = tk.Button(bar, text=text, bd=2, padx=8, pady=1, cursor="hand2", command=lambda t=tool: self.set_tool(t))
            b.pack(side=tk.LEFT, padx=(0, 2))
            Tooltip(b, tip)
            self.tool_buttons[tool] = b
        ttk.Separator(bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=6)
        # New annotations' colors, for this run only (red and white again when the viewer starts): shapes and lines, and text.
        self.ann_color, self.text_color = "#ff0000", "#ffffff"
        self.settings.pop("annotation_color", None)        # (remembered by earlier versions)
        self.color_button = tk.Button(bar, text="   ", bg=self.ann_color, activebackground=self.ann_color, bd=2, cursor="hand2",
                                      command=self.choose_color)
        self.color_button.keep_colors = True
        self.color_button.pack(side=tk.LEFT, padx=(0, 2))
        Tooltip(self.color_button, "Color of new shapes and lines (red at first), and of the selected one. With Text chosen, or a text "
                                   "selected: the color of new text (white at first), and of the selected text.")
        for text, cmd, tip in (("Delete", self.delete_annotation, "Delete the selected annotation (or press Delete)"),
                               ("Undo", self.undo_annotation, "Undo the last change: annotations or color maps, on any image (Ctrl+Z)")):
            b = ttk.Button(bar, text=text, command=cmd)
            b.pack(side=tk.LEFT, padx=(4, 0))
            Tooltip(b, tip)
        self.unsaved_label = ttk.Label(bar, text="", style="Warn.TLabel")          # ("unsaved changes": Save is left of the list)
        settings_button = ttk.Button(bar, text="⚙ Settings", command=self.open_settings)
        settings_button.pack(side=tk.RIGHT)
        Tooltip(settings_button, "Settings: color theme, zoomed-in pixels, renaming")
        self.flip_x = tk.BooleanVar(value=self.settings.get("flip_left_right", False))
        flip_check = ttk.Checkbutton(bar, text="Flip all sections", variable=self.flip_x, command=self.on_flip_toggle)
        flip_check.pack(side=tk.RIGHT, padx=(0, 8))
        Tooltip(flip_check, "Show every section mirrored, left for right (with its annotations); a section with \"Flip this section\" "
                            "ticked is flipped back. Only the view is flipped: the files, and annotations drawn while flipped, keep the "
                            "image's own left and right.")
        self.unsaved_label.pack(side=tk.RIGHT, padx=6)
        self.tool, self.selected, self.ann_drag, self.save_warned = "select", None, None, False

        # Right of the image: the channels side by side (checkbox, name, vertical black/white slider, Auto), then the buttons for all.
        view_row = ttk.Frame(right)
        view_row.pack(fill=tk.BOTH, expand=True)
        side = ttk.Frame(view_row)
        side.pack(side=tk.RIGHT, fill=tk.Y, padx=4, pady=4)        # (packed before the image, so a narrow window shrinks the image)
        self.z_frame = ttk.Frame(side)
        ttk.Label(self.z_frame, text="Z plane").pack(side=tk.LEFT)
        self.z_var = tk.IntVar(value=0)
        self.z_scale = tk.Scale(self.z_frame, from_=0, to=0, orient=tk.HORIZONTAL, variable=self.z_var, command=self.on_z)
        self.z_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

        self.split_button = ttk.Button(side, text="Merge channels" if self.settings.get("separate_channels", True) else "Separate channels",
                                       command=self.toggle_split)
        self.split_button.pack(fill=tk.X, pady=(0, 4))
        Tooltip(self.split_button, "Show each channel switched on in a panel of its own, side by side. Zooming or panning one moves "
                                   "them all.")
        self.flip_one = tk.BooleanVar(value=False)
        self.flip_one_check = ttk.Checkbutton(side, text="Flip this section", variable=self.flip_one, command=self.on_flip_one_toggle,
                                              state=tk.DISABLED)
        self.flip_one_check.pack(anchor=tk.W, pady=(0, 4))
        Tooltip(self.flip_one_check, "Show this section mirrored, left for right, e.g. one mounted the other way up from the rest "
                                     "(Save writes it into the section's .xlif, so it goes with the files; or, see Settings, it is kept on this "
                                     "computer). With \"Flip all "
                                     "sections\" ticked too, it is shown unflipped. Only the view is flipped, not the image.")
        self.ch_frame = ttk.LabelFrame(side, text="Channels")
        self.ch_frame.pack(fill=tk.BOTH, expand=True)
        btns = ttk.Frame(side)
        btns.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(btns, text="Auto contrast", command=self.auto_all).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(btns, text="LAS X contrast", command=self.lasx_all).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=(4, 0))

        self.canvas = self.main_canvas = tk.Canvas(view_row, background="black", highlightthickness=0, cursor="fleur")
        self.canvas.keep_colors = True        # (black in both themes)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # Channels separated (see _layout_panels): the image's canvases are a grid of panels in panel_frame, one per channel shown,
        # all with the one view (scale, cx, cy: in a panel's own coordinates); self.canvas is the first, and the one being drawn
        # while _each_canvas draws on each in turn.
        self.split = self.settings.get("separate_channels", True)      # (on by default; remembered)
        self.canvases, self.panel_channels, self.panel_cols, self.drawing_panel = [self.canvas], [None], 1, None
        self.panel_frame = tk.Frame(view_row, background="#3a3a3a")
        self.panel_frame.keep_colors = True
        self.panel_frame.bind("<Configure>", lambda e: self._layout_panels())
        self.pixel_grid = False      # (at the largest zoom, after trying to zoom in further: a line along every pixel's edge)
        self.status = ttk.Label(right, text="Open an .xlef file (File > Open, or drag it onto this window or onto RUN_XLEF_VIEWER.bat).",
                                anchor=tk.W)
        self.status.pack(fill=tk.X)
        self.cursor_row = ttk.Frame(right)          # (under the pointer: its position, then each channel's value, in its color)
        self.cursor_row.pack(fill=tk.X)
        self.cursor_label = ttk.Label(self.cursor_row, text="", anchor=tk.W, style="Muted.TLabel")
        self.cursor_label.pack(side=tk.LEFT)
        self.cursor_values = []                      # ([tk.Label per channel], the channels' colors they were made for)

        self._bind_view_canvas(self.canvas)
        self.root.bind("<Delete>", lambda e: None if isinstance(e.widget, (tk.Entry, tk.Text)) else self.delete_annotation())
        self.root.bind("<Control-z>", lambda e: None if isinstance(e.widget, tk.Text) else self.undo_annotation())
        self.root.bind("<Escape>", lambda e: (setattr(self, "selected", None), self.set_tool("select"), self._draw_annotations()))
        self.set_tool("select")
        self.root.bind("<Control-o>", lambda e: self.ask_open())
        self.root.bind("<Home>", lambda e: self.fit())
        for k in range(1, 10):
            self.root.bind(str(k), lambda e, k=k: self.toggle_channel(k - 1))
        if self.dnd_files:
            self.tree.drop_target_register(self.dnd_files)
            self.tree.dnd_bind("<<Drop>>", self.on_drop)

    def _bind_view_canvas(self, c):
        """The mouse on an image canvas (the one, or each panel when the channels are separated: the same view in each)."""
        c.bind("<Configure>", lambda e: self.on_resize())
        c.bind("<ButtonPress-1>", self.on_press)
        c.bind("<B1-Motion>", self.on_drag)
        c.bind("<ButtonRelease-1>", self.on_release)
        c.bind("<MouseWheel>", self.on_wheel)
        c.bind("<Button-4>", lambda e: self.on_wheel(e, 120))
        c.bind("<Button-5>", lambda e: self.on_wheel(e, -120))
        c.bind("<Motion>", self.on_motion)
        c.bind("<Leave>", lambda e: self._hover_handles(None, None))
        c.bind("<Double-Button-1>", self.ann_double)
        for b in (2, 3):          # the middle and right buttons always move the view, whatever the tool
            c.bind(f"<ButtonPress-{b}>", self.pan_start)
            c.bind(f"<B{b}-Motion>", self.pan_move)
            c.bind(f"<ButtonRelease-{b}>", self.pan_end)
        if self.dnd_files:
            c.drop_target_register(self.dnd_files)
            c.dnd_bind("<<Drop>>", self.on_drop)

    # --- separate channels

    def toggle_split(self):
        """The "Separate channels" button: each channel switched on in a panel of its own (or, again, all mixed in one)."""
        self.split = self.settings["separate_channels"] = not self.split
        self._save_settings()
        self.split_button.config(text="Merge channels" if self.split else "Separate channels")
        self._layout_panels(force=True)

    def _layout_panels(self, force=False):
        """The image's canvases: the one, or (channels separated) a grid of panels, one per channel switched on, in as many columns as
        make the image largest in them. Rebuilt only when that changes (the channels shown, or the best number of columns)."""
        if self.split and self.img is None:
            return                          # (between two images: kept as they are)
        want = [c for c, st in enumerate(self.ch_settings) if st["on"]] if self.split else []
        was_fit = self.img is not None and abs(self.scale - self.fit_scale) < 1e-9
        if not want:
            if self.canvases != [self.main_canvas]:
                for c in self.canvases:
                    c.destroy()
                self.panel_frame.pack_forget()
                self.main_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
                self.canvases, self.canvas, self.panel_channels, self.panel_cols = [self.main_canvas], self.main_canvas, [None], 1
                self.photo = None
                self.root.after_idle(self._relaid_out, was_fit)
            return
        area = self.panel_frame if self.panel_frame.winfo_ismapped() else self.main_canvas
        w, h = max(area.winfo_width(), 1), max(area.winfo_height(), 1)
        e, n = self.img.entry, len(want)
        cols = max(range(1, n + 1), key=lambda k: min(w / k / e.width, h / -(-n // k) / e.height))
        if not force and want == self.panel_channels and cols == self.panel_cols:
            return
        rows = -(-n // cols)
        for c in self.canvases:
            if c is not self.main_canvas:
                c.destroy()
        for k in range(max(n, 9)):
            self.panel_frame.grid_columnconfigure(k, weight=1 if k < cols else 0, uniform="panel" if k < cols else "")
            self.panel_frame.grid_rowconfigure(k, weight=1 if k < rows else 0, uniform="panel" if k < rows else "")
        panels = []
        for k in range(n):
            c = tk.Canvas(self.panel_frame, background="black", highlightthickness=0, cursor=self.main_canvas.cget("cursor"))
            c.keep_colors = True
            c.grid(row=k // cols, column=k % cols, sticky="nsew", padx=(0 if k % cols == 0 else 2, 0), pady=(0 if k < cols else 2, 0))
            self._bind_view_canvas(c)
            panels.append(c)
        if self.main_canvas.winfo_ismapped() or not self.panel_frame.winfo_ismapped():
            self.main_canvas.pack_forget()
            self.panel_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.canvases, self.canvas, self.panel_channels, self.panel_cols = panels, panels[0], want, cols
        self.photo = None
        self.root.after_idle(self._relaid_out, was_fit)

    def _relaid_out(self, was_fit):
        """After the canvases changed (separated, merged, a panel more or less): a view that filled the old canvas fills the new one;
        a zoomed-in view keeps its zoom. (The one canvas, shown again, keeps its old size: no <Configure> comes to do this.)"""
        self.root.update_idletasks()                    # (the canvases' new sizes)
        if self.img is None:
            return
        if was_fit:
            self.fit()
        else:
            self.clamp_view()
            self.schedule_redraw(0)

    def _each_canvas(self, draw):
        """draw() once per canvas of the image, with self.canvas set to each in turn (self.drawing_panel: its index)."""
        if self.drawing_panel is not None:
            draw()
            return
        keep = self.canvas
        try:
            for k, c in enumerate(self.canvases):
                self.canvas, self.drawing_panel = c, k
                draw()
        finally:
            self.canvas, self.drawing_panel = keep, None

    def _mix(self, arrays, luts):
        """The channels' arrays mixed into one RGB picture per canvas: all of them together, or (separated) each on its own."""
        if self.panel_channels == [None]:
            return [composite(arrays, luts)]
        return [composite([a if i == c else None for i, a in enumerate(arrays)], [t if i == c else None for i, t in enumerate(luts)])
                for c in self.panel_channels]

    def _draw_panel_label(self):
        """Separated: the channel's number and name at the top left of its panel, in its color."""
        c = self.panel_channels[self.drawing_panel or 0]
        if c is None or self.img is None or c >= len(self.img.entry.channels):
            return
        ch = self.img.entry.channels[c]
        t = self.canvas.create_text(8, 6, anchor=tk.NW, text=f"{c + 1}  {ch.name}", fill=self._swatch_color(c),
                                    font=("Segoe UI", 10, "bold"))
        x0, y0, x1, y1 = self.canvas.bbox(t)
        self.canvas.tag_lower(self.canvas.create_rectangle(x0 - 4, y0 - 2, x1 + 4, y1 + 2, fill="#101010", outline=""), t)

    def _draw_pixel_grid(self):
        """At the largest zoom, once the wheel has tried to zoom in further: a line along the edges of every pixel (over the image only)."""
        if not self.pixel_grid or self.img is None or self.scale < MAX_ZOOM - 1e-9:
            return
        e = self.img.entry
        cw, ch = self.canvas_size()
        ix0, iy0, ix1, iy1 = self._view_corners()
        sx = lambda x: round(self._fx((x - self.cx) * self.scale + cw / 2))
        sy = lambda y: round((y - self.cy) * self.scale + ch / 2)
        top, bottom = max(0, sy(0)), min(ch, sy(e.height))
        left, right = max(0, min(sx(0), sx(e.width))), min(cw, max(sx(0), sx(e.width)))
        for x in range(max(0, int(np.ceil(ix0))), min(e.width, int(ix1)) + 1):
            self.canvas.create_line(sx(x), top, sx(x), bottom, fill=PIXEL_GRID_COLOR, tags="pixelgrid")
        for y in range(max(0, int(np.ceil(iy0))), min(e.height, int(iy1)) + 1):
            self.canvas.create_line(left, sy(y), right, sy(y), fill=PIXEL_GRID_COLOR, tags="pixelgrid")

    # --- settings window, color theme

    def open_settings(self):
        """The Settings window: the color theme (applied at once), zoomed-in pixels, whether renaming also renames files on disk,
        whether "Flip this section" is saved into LAS X's files, and the memory for previews. Remembered."""
        if getattr(self, "settings_window", None) is not None and self.settings_window.winfo_exists():
            self.settings_window.lift()
            return
        dlg = self.settings_window = tk.Toplevel(self.root)
        dlg.title("Settings")
        dlg.transient(self.root)
        dlg.resizable(False, False)
        dlg.geometry(f"+{self.root.winfo_rootx() + 200}+{self.root.winfo_rooty() + 120}")
        body = ttk.Frame(dlg, padding=12)
        body.pack(fill=tk.BOTH, expand=True)
        box = ttk.LabelFrame(body, text="Color theme", padding=(10, 6))
        box.pack(fill=tk.X)
        choice = tk.StringVar(dlg, self.theme)
        for value, text in (("dark", "Dark"), ("light", "Light")):
            ttk.Radiobutton(box, text=text, value=value, variable=choice, command=lambda: self.apply_theme(choice.get())).pack(anchor=tk.W)
        box = ttk.LabelFrame(body, text="Zoomed-in pixels", padding=(10, 6))
        box.pack(fill=tk.X, pady=(10, 0))
        method = tk.StringVar(dlg, self.enlarge_method)

        def set_method():
            self.enlarge_method = self.settings["enlarge_method"] = method.get()
            self._save_settings()
            self.redraw()                   # (drawn again at once: from the pixels already here)
        for value, text in ENLARGE_METHODS.items():
            ttk.Radiobutton(box, text=text, value=value, variable=method, command=set_method).pack(anchor=tk.W)
        box = ttk.LabelFrame(body, text="Renaming", padding=(10, 6))
        box.pack(fill=tk.X, pady=(10, 0))
        on_disk = tk.BooleanVar(dlg, self.settings.get("rename_files", False))

        def set_on_disk():
            self.settings["rename_files"] = on_disk.get()
            self._save_settings()
        ttk.Checkbutton(box, text="Also rename the files and folders on disk", variable=on_disk, command=set_on_disk).pack(anchor=tk.W)
        ttk.Label(box, text=RENAME_NOTE, style="Muted.TLabel", justify=tk.LEFT).pack(anchor=tk.W, pady=(2, 0))
        box = ttk.LabelFrame(body, text="Flip this section", padding=(10, 6))
        box.pack(fill=tk.X, pady=(10, 0))
        in_files = tk.BooleanVar(dlg, self._flip_in_files())
        ttk.Checkbutton(box, text="Save it in a custom field of LAS X's files", variable=in_files,
                        command=lambda: self.set_flip_in_files(in_files.get())).pack(anchor=tk.W)
        ttk.Label(box, text=FLIP_NOTE, style="Muted.TLabel", justify=tk.LEFT).pack(anchor=tk.W, pady=(2, 0))
        box = ttk.LabelFrame(body, text="Memory for previews", padding=(10, 6))
        box.pack(fill=tk.X, pady=(10, 0))
        row = ttk.Frame(box)
        row.pack(anchor=tk.W)
        percent = tk.IntVar(dlg, self.settings.get("preview_memory_percent", KEPT_IMAGES_PERCENT))
        amount = ttk.Label(row)

        def set_percent(*_):
            try:
                value = min(80, max(1, int(percent.get())))
            except (tk.TclError, ValueError):
                return                      # (while typing)
            more = value > self.settings.get("preview_memory_percent", KEPT_IMAGES_PERCENT)
            self.settings["preview_memory_percent"] = value
            self._save_settings()
            amount.config(text=f"% of this computer's memory ({size_text(self.kept_limit())})")
            if more and self.preload_full:
                self._projects_changed()    # (the background loader goes on with the sections it left out)
            elif not more:
                self._keep_within_limit()
        ttk.Spinbox(row, from_=1, to=80, increment=5, width=4, textvariable=percent).pack(side=tk.LEFT)
        amount.pack(side=tk.LEFT, padx=(4, 0))
        percent.trace_add("write", set_percent)
        set_percent()
        ttk.Label(box, text=PREVIEW_MEMORY_NOTE, style="Muted.TLabel", justify=tk.LEFT).pack(anchor=tk.W, pady=(2, 0))
        ttk.Button(body, text="Close", command=dlg.destroy).pack(anchor=tk.E, pady=(10, 0))
        dlg.bind("<Escape>", lambda e: dlg.destroy())

    def _tree_tag_colors(self):
        self.tree.tag_configure("cloud", foreground=PAL["link"])
        self.tree.tag_configure("bad", foreground=PAL["bad"])
        self.tree.tag_configure("dim", foreground=PAL["dim"])             # (an image with nothing of it in memory yet)
        self.tree.tag_configure("cloud_dim", foreground=PAL["link_dim"])
        self.tree.tag_configure("project", font=("Segoe UI", 10, "bold"))
        for bar in getattr(self, "list_bars", {}).values():
            bar.configure(background=PAL["placeholder"])
        if getattr(self, "unsaved_marks", None) or getattr(self, "list_bars", None):
            self.root.after_idle(self._place_list_bars)     # (the marks and bars in the theme's colors)

    def apply_theme(self, name):
        """Switch the color theme, for everything already on screen too, and remember it."""
        if name == self.theme or name not in THEMES:
            return
        old = dict(PAL)
        style_theme(self.root, name, self.native_theme)
        self.theme = self.settings["theme"] = name
        self.root.configure(background=PAL["bg"])
        recolor(self.root, old, PAL)
        for w in self.root.winfo_children():               # (other windows, e.g. Settings)
            if isinstance(w, tk.Toplevel):
                w.configure(background=PAL["bg"])
        self._tree_tag_colors()
        self.show_view_mode()
        self.set_tool(self.tool)
        if self.img:
            for i in range(len(self.ch_widgets)):
                self._show_channel_state(i)
        self.thumb_grid.layout()
        self._save_settings()

    def on_drop(self, event):
        names = self.root.tk.splitlist(event.data)
        if names:
            self.root.after(10, lambda: self.open_xlef(names[0]))
        return event.action

    # --- settings

    def _load_settings(self):
        try:
            with open(SETTINGS_PATH, encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}

    def _save_settings(self):
        try:
            if self.root.state() == "normal":
                self.normal_geometry = self.root.geometry()
            if self.normal_geometry:
                self.settings["geometry"] = self.normal_geometry     # (not a maximized window's: see _window_moved)
            with open(SETTINGS_PATH, "w", encoding="utf-8") as f:
                json.dump(self.settings, f, indent=1)
        except OSError:
            pass

    def _window_moved(self, event):
        """The window resized or moved: remembered when it is a normal window (maximized, Tk gives the full screen's size with the
        position from before: restored as a normal window that would hang off the screen)."""
        if event.widget is self.root and self.root.state() == "normal":
            self.normal_geometry = self.root.geometry()

    def quit(self):
        if not self.unsaved_ok("Closing the viewer", "Close anyway?"):
            return
        self._save_settings()
        self.generation += 1
        if self.img:
            self.img.close()
        self.root.destroy()

    # --- opening

    def add_recent(self, path):
        recent = [p for p in self.settings["recent_files"] if os.path.normcase(p) != os.path.normcase(path)]
        self.settings["recent_files"] = [path] + recent[:MAX_RECENT - 1]
        self._save_settings()

    def remove_recent(self, path):
        self.settings["recent_files"] = [p for p in self.settings["recent_files"] if os.path.normcase(p) != os.path.normcase(path)]
        self._save_settings()

    def ask_open(self):
        """Open an experiment: one from the list of recently opened ones, or an .xlef browsed for, typed or dragged in (open_file_dialog's
        dialog)."""
        def details(path):
            when = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(path)))
            return {"images": str(count_images(path)), "when": when}
        result = ask_open_file(
            self.root, self.settings["recent_files"], title="Open a Leica experiment", noun="experiment",
            filetypes=[("Leica experiment", "*.xlef *.xlif"), ("All files", "*.*")],
            columns=[("name", "Experiment", 280, "w", True), ("images", "Images", 70, "e", False), ("folder", "Folder", 420, "w", True),
                     ("when", "Modified", 130, "w", False)],
            details=details, on_remove=self.remove_recent, drop=self.dnd_files is not None)
        if result:
            path, from_list = result
            if from_list and not os.path.isfile(path):
                messagebox.showerror("Can't open", f"This experiment no longer exists:\n{path}", parent=self.root)
                self.remove_recent(path)
                return
            self.open_xlef(path)

    def open_xlef(self, path):
        path = os.path.normpath(path)
        if not os.path.isfile(path):
            messagebox.showerror("XLEF viewer", f"File not found:\n{path}", parent=self.root)
            return
        for pid, project in self.projects.items():         # (open already: just go to it)
            if os.path.normcase(project["path"]) == os.path.normcase(path):
                self.tree.selection_set(pid)
                self.tree.see(pid)
                return
        # Several projects can be open: each is a branch of the list (and a group of thumbnails). The new one's metadata (the .xlef and
        # .xlif files, a few hundred KB each) is read in the background, so that the window stays alive on a slow drive;
        # _experiment_read() takes it from there.
        name = os.path.splitext(os.path.basename(path))[0]
        pid = self.tree.insert("", tk.END, text=f"{name}   (loading...)", open=True, tags=("project",))
        self.projects[pid] = {"path": path, "name": name}
        self.tree.see(pid)
        if self.img is None:
            self.show_message("Loading...")
        if not self.entries:
            self.thumb_grid.set_items([], message="Loading...")
        self._update_title()
        self.status.config(text=f"Reading {os.path.basename(path)} ...")

        def work():
            try:
                self.results.put(("experiment", pid, (path, load_experiment(path), None)))
            except Exception as e:      # a bad file must not take the window down
                self.results.put(("experiment", pid, (path, None, e)))

        threading.Thread(target=work, daemon=True).start()

    def _experiment_read(self, pid, path, experiment, error):
        if pid not in self.projects:            # (closed while it was read)
            return
        if error is not None:
            self.tree.delete(pid)
            self.projects.pop(pid)
            self._update_title()
            if self.img is None:
                self.show_message("")
            self.status.config(text=f"Can't read {os.path.basename(path)}")
            messagebox.showerror("XLEF viewer", f"Can't read {path}:\n{error}", parent=self.root)
            return
        self.add_recent(path)
        self.tree.item(pid, text=self.projects[pid]["name"])
        self.projects[pid]["root"] = experiment
        experiment.path = path
        self._fill_tree(pid, experiment)
        self._measure_sizes(pid)
        mine = self._entries_of(pid)
        n_cloud = sum(1 for i in mine if i in self.item_cloud)        # (as _fill_tree found them: not asked of each file again)
        self.status.config(text=f"{self.projects[pid]['name']}: {len(mine)} image{'s' * (len(mine) != 1)}"
                                + (f", {n_cloud} online-only" if n_cloud else "") + ". Click one to view it.")
        self._projects_changed()
        first = next(iter(mine), None)
        if first and not self._files_to_ask(mine[first]):
            self.tree.selection_set(first)          # (it opens, and the thumbnails wait for it: see thumbs_go)
            self.tree.see(first)
            self.on_select()
        elif self.img is None:
            self.show_message("Select an image" if mine else "No images in this experiment")
            self._start_slide(pid)

    def _measure_sizes(self, pid):
        """The size of each of the project's images at full resolution (its channel TIFFs), read in the background from the file
        system only (os.stat: an online-only file is not downloaded), then shown in the list (see _show_sizes)."""
        entries = list(self._entries_of(pid).values())

        def work():
            for e in entries:
                total = 0
                for f in e.all_files():
                    try:
                        total += os.path.getsize(f)
                    except OSError:
                        pass
                e.full_bytes = total
            self.results.put(("sizes", pid, None))

        threading.Thread(target=work, daemon=True).start()

    def _show_sizes(self):
        """The size column of the list: each image's full resolution, and the total of what each folder and project lists now
        (images dragged in or out, or deleted, counted as they are shown)."""
        def total(node):
            entry = self.entries.get(node)
            if entry is not None:
                b = entry.full_bytes
            else:
                sizes = [total(child) for child in self.tree.get_children(node)]
                b = sum(x for x in sizes if x is not None) if any(x is not None for x in sizes) else None
            if node:
                self.tree.set(node, "size", size_text(b) if b is not None else "")
            return b
        total("")

    def _projects_changed(self):
        """After a project was opened or closed: the thumbnail grid, the title, and the background work (thumbnails, first pictures)
        start again over all open projects' images (what was done already comes from the cache or memory)."""
        self.thumb_gen += 1
        self.thumbs_started = False
        self.preload_full = []
        self.thumb_grid.set_items(self._grid_rows_all(), message="No images" if self.projects else "Open an experiment", keep=True)
        if self.view_mode.get() == "thumbs":
            self._start_thumbnails()
        self._start_preload()
        self._update_title()

    def _update_title(self):
        names = [p["name"] for p in self.projects.values()]
        self.root.title("XLEF viewer" + (f" - {names[0]}" if len(names) == 1 else f" - {len(names)} experiments" if names else ""))
        self.file_label.config(text=next(iter(self.projects.values()))["path"] if len(names) == 1
                               else f"{len(names)} experiments open" if names else "No file open")

    def _project_of(self, item):
        while item and self.tree.parent(item):
            item = self.tree.parent(item)
        return item if item in self.projects else None

    def _entries_of(self, pid):
        return {i: e for i, e in self.entries.items() if self._project_of(i) == pid}

    def _project_items(self, pid):
        """Every item under project `pid` in the list, at any depth."""
        items = []

        def walk(node):
            for child in self.tree.get_children(node):
                items.append(child)
                walk(child)
        walk(pid)
        return items

    def _unsaved_names(self, pid):
        """The names of project `pid`'s images and folders (and the project itself) with unsaved changes, and of its deleted images."""
        root = self.projects[pid].get("root")
        objects = [self.entries.get(i) or self.folders.get(i) for i in self._project_items(pid)]
        changed = [o.name for o in objects if o is not None and o.unsaved]
        changed += [self.projects[pid]["name"]] if root is not None and root.unsaved else []
        changed += [e.name + " (deleted)" for i, e in self.deleted.items() if self.home.get(i) == pid]
        return changed

    def close_project(self, pid, ask=True):
        """Close one open project (asking first if any of its images or folders have unsaved changes)."""
        if pid not in self.projects:
            return
        items = self._project_items(pid)
        crossing = [i for i in self.home if self.tree.exists(i) and (self.home[i] == pid) != (self._project_of(i) == pid)
                    and i not in self.deleted]
        if ask and crossing:
            messagebox.showinfo("Close project", f"Images were dragged between {self.projects[pid]['name']} and another project. Save "
                                                 "those moves, or undo them, before closing it.", parent=self.root)
            return
        changed = self._unsaved_names(pid)
        if ask and changed and not messagebox.askyesno(
                "Unsaved changes", f"{', '.join(changed[:5])}{' and more' if len(changed) > 5 else ''} "
                                   f"{'has' if len(changed) == 1 else 'have'} unsaved changes.\n\nClosing "
                                   f"{self.projects[pid]['name']} discards them. Close anyway?",
                icon="warning", default="no", parent=self.root):
            return
        closing = {self.entries[i] for i in items if i in self.entries}
        if self.img is not None and self.img.entry in closing:
            self._close_image()
            self.show_message("Select an image")
        for e in closing:
            kept = self.kept_images.pop(e, None)
            if kept is not None:
                kept.close()
        self._kept_changed()
        gone = [i for i in self.deleted if self.home.get(i) == pid]
        for i in gone:
            if self.tree.exists(i):
                self.tree.delete(i)                 # (deleted images are kept out of the list, detached, until Save)
        items += gone
        for i in items:
            self.entries.pop(i, None)
            self.folders.pop(i, None)
            self.deleted.pop(i, None)
            self.home.pop(i, None)
            self.item_cloud.discard(i)
            self.preload_progress.pop(i, None)
            self.download_progress.pop(i, None)
            if i in self.download_queue:
                self.download_queue.remove(i)
        self.tree.delete(pid)
        self.projects.pop(pid)
        if self.slide_project == pid:
            self.slide_project = None
            self.slide_gen += 1
            self.slide.pack_forget()
        self._projects_changed()

    def _label_item(self, item):
        """An image's or folder's name in the list: its name (a cloud after an image with online-only files), a red * after it (see _place_unsaved_mark) while it has
        unsaved changes (a rename, or edited annotations), and how far the reading of its preview has got ("45% preview": in the background, or as it is opened) while that runs. An
        image with nothing of it in memory yet (not open, not kept or preloaded: see _loaded_entries) is dimmed."""
        obj = self.entries.get(item) or self.folders.get(item)
        if obj is None or not self.tree.exists(item):
            return
        if item in self.entries:
            cloud = item in self.item_cloud
            if obj.problem:
                tags = ("bad",)
            elif obj in self._loaded_entries():
                tags = ("cloud",) if cloud else ()
            else:
                tags = ("cloud_dim",) if cloud else ("dim",)
            self.tree.item(item, tags=tags)
            self.thumb_grid.set_dimmed(item, tags in (("cloud_dim",), ("dim",)))
        text = obj.name + ("  ☁" if item in self.item_cloud else "")
        if item in self.preload_progress:
            text += f"  {self.preload_progress[item]}% preview"
        elif self.open_progress is not None and self.open_progress[0] == item:
            text += f"  {self.open_progress[1]}% preview"
        if item in self.download_progress:
            p = self.download_progress[item]
            text += "  ⇣ waiting" if p == "waiting" else f"  ⇣ {p}%"
        self.tree.item(item, text=text)
        self._place_list_bar(item)
        self._place_unsaved_mark(item)

    def _list_bar_fraction(self, item):
        """How far the bar under an image's name is (0-1) and its color: its preview being read (blue), else its full resolution
        being downloaded from OneDrive (white; grey in the light theme); or (None, None)."""
        if item in self.preload_progress:
            return self.preload_progress[item] / 100, "#3a8fff"
        if self.open_progress is not None and self.open_progress[0] == item:
            return self.open_progress[1] / 100, "#3a8fff"
        if item in self.download_progress:
            p = self.download_progress[item]
            return (0 if p == "waiting" else p / 100), PAL["download_bar"]
        return None, None

    def _place_list_bar(self, item):
        """The thin bar under an image's name in the list while its preview is read (blue, as along a thumbnail's bottom) or its full
        resolution downloaded (white): a small canvas laid over the list, moved whenever the rows move (scrolling, folders opened or
        closed, resizing)."""
        frac, color = self._list_bar_fraction(item) if item in self.entries else (None, None)
        bar = self.list_bars.get(item)
        box = self.tree.bbox(item, "#0") if frac is not None and self.tree.exists(item) else ""      # (the name's column)
        if not box:
            if bar is not None:
                if frac is None or not self.tree.exists(item):
                    bar.destroy()
                    del self.list_bars[item]
                else:
                    bar.place_forget()           # (scrolled out of view, or in a closed folder: back when it shows again)
            return
        x, y, w, h = box
        start = x + 20 * (self._tree_depth(item) + 1)            # (where the name starts, as in _tree_motion)
        width = max(40, min(x + w, self.tree.winfo_width()) - start - 8)     # (to the right edge of the list, whatever the name's length)
        if bar is None:
            bar = self.list_bars[item] = tk.Canvas(self.tree, height=3, highlightthickness=0, borderwidth=0,
                                                   background=PAL["placeholder"])
            bar.create_rectangle(0, 0, 0, 3, fill=color, outline="", tags=("fill",))
            self._pass_mouse_to_list(bar)
        bar.coords("fill", 0, 0, width * frac, 3)
        bar.itemconfigure("fill", fill=color)
        bar.configure(width=width)
        bar.place(x=start, y=y + h - 4, width=width, height=3)

    def _place_list_bars(self):
        for item in list(self.list_bars):
            self._place_list_bar(item)
        for item in list(self.unsaved_marks):
            self._place_unsaved_mark(item)

    def _pass_mouse_to_list(self, widget):
        """The mouse on something laid over the list (a progress bar, an unsaved mark) acts on the list under it: selecting,
        dragging, the menu, the wheel, hover notes."""
        for ev in ("<ButtonPress-1>", "<B1-Motion>", "<ButtonRelease-1>", "<Double-Button-1>", "<Button-3>", "<Motion>"):
            widget.bind(ev, lambda e, ev=ev, w=widget: self.tree.event_generate(
                ev, x=e.x + w.winfo_x(), y=e.y + w.winfo_y(), rootx=e.x_root, rooty=e.y_root))
        widget.bind("<MouseWheel>", lambda e, w=widget: self.tree.event_generate(
            "<MouseWheel>", x=e.x + w.winfo_x(), y=e.y + w.winfo_y(), rootx=e.x_root, rooty=e.y_root, delta=e.delta))

    def _is_unsaved(self, item):
        if item in self.projects:
            root = self.projects[item].get("root")
            return root is not None and root.unsaved
        obj = self.entries.get(item) or self.folders.get(item)
        return obj is not None and obj.unsaved

    def _unsaved_what(self, item):
        """What isn't saved yet of this image, folder or project, for the note on its red *: e.g. ["contrast", "color maps"]."""
        obj = self.projects[item].get("root") if item in self.projects else (self.entries.get(item) or self.folders.get(item))
        if obj is None:
            return []
        what = []
        if obj.name != obj.saved_name:
            what.append(f"renamed (was \"{obj.saved_name}\")")
        if getattr(obj, "annotations_changed", False):
            what.append("annotations")
        if isinstance(obj, ImageEntry):
            what += ["color maps"] * obj.luts_changed + ["contrast"] * obj.contrast_changed + ["flip"] * obj.flip_changed
        if getattr(obj, "order_changed", False):
            what.append("the order of what it holds")
        return what

    def _unsaved_tip(self, item, event):
        """The pointer on a row's red *: a note saying what isn't saved yet. True if it is on one."""
        mark = self.unsaved_marks.get(item)
        if mark is None or not mark.winfo_manager():
            return False
        left = int(float(mark.place_info().get("x", 0)))
        if not left - 3 <= event.x <= left + mark.winfo_reqwidth() + 3:
            return False
        what = self._unsaved_what(item)
        if self.tree_tip is None:
            self.tree_tip = PointerTip(self.tree)
        self.tree_tip.show("Unsaved changes" + (": " + ", ".join(what) if what else "") +
                           ".\nSave writes them into LAS X's files.", event.x_root, event.y_root)
        return True

    def _place_unsaved_mark(self, item):
        """A red * at the end of the row's name while an image, folder or project has unsaved changes (the list can't color part of a
        row, so it is a label laid over the list, moved with the row as the progress bars are). On a name too long for the list, it
        sits at the list's right edge, over the name's cut-off end. The pointer on it: a note saying what isn't saved (_unsaved_tip)."""
        mark = self.unsaved_marks.get(item)
        if not self.tree.exists(item) or not self._is_unsaved(item):
            if mark is not None:
                mark.destroy()
                del self.unsaved_marks[item]
            return
        box = self.tree.bbox(item, "#0")
        if mark is None:
            mark = self.unsaved_marks[item] = tk.Label(self.tree, text="*", bd=0, padx=1, pady=0, font=("Segoe UI", 12, "bold"))
            mark.keep_colors = True
            self._pass_mouse_to_list(mark)
        if not box:
            mark.place_forget()                 # (scrolled out of view, or in a closed folder)
            return
        x, y, w, h = box
        font = tkfont.Font(font=("Segoe UI", 10, "bold")) if item in self.projects else tkfont.nametofont("TkDefaultFont")
        start = x + 20 * (self._tree_depth(item) + 1)            # (where the name starts, as in _tree_motion)
        left = min(start + font.measure(self.tree.item(item, "text")) + 3, x + w - mark.winfo_reqwidth() - 2)
        style = ttk.Style()
        selected = item in self.tree.selection()
        bg = (style.lookup("Treeview", "background", ("selected",)) or PAL["selected"]) if selected else             (style.lookup("Treeview", "background") or PAL["field"])
        mark.configure(fg="#ffb4b4" if selected else PAL["warn"], bg=bg)
        if left < start:
            mark.place_forget()                 # (no room at all in a very narrow list)
            return
        mark.place(x=left, y=y + (h - mark.winfo_reqheight()) // 2)

    def _fill_tree(self, parent, folder):
        """The folder's images and sub-folders, in LAS X's order (the order its file lists them in, which dragging changes)."""
        for child in folder.children:
            if isinstance(child, FolderEntry):
                item = self.tree.insert(parent, tk.END, text=child.name, open=True)
                self.folders[item] = child
                self.home[item] = self._project_of(parent)
                self._fill_tree(item, child)
                continue
            cloud = bool(child.online_only_files())
            tags = ("bad",) if child.problem else (("cloud",) if cloud else ())
            item = self.tree.insert(parent, tk.END, text=child.name, tags=tags)
            self.entries[item] = child
            self.home[item] = self._project_of(parent)
            if cloud:
                self.item_cloud.add(item)
            self._label_item(item)

    # --- list / thumbnail view

    def show_view_mode(self):
        thumbs = self.thumbs_shown = self.view_mode.get() == "thumbs"
        for value, b in self.view_buttons.items():
            on = value == self.view_mode.get()
            style_toggle(b, on)
        (self.list_frame.pack_forget if thumbs else self.thumb_grid.pack_forget)()
        (self.thumb_grid if thumbs else self.list_frame).pack(fill=tk.BOTH, expand=True)
        self.settings["view_mode"] = self.view_mode.get()
        sel = self.tree.selection()
        if thumbs:
            self._start_thumbnails()
            if sel:
                self.root.after_idle(lambda: self.thumb_grid.see(sel[0]))
        elif sel:
            self.tree.see(sel[0])

    def _grid_rows_all(self):
        """The thumbnail grid's rows: each open project under a heading of its own (a line across, its name), then its folders and
        images."""
        rows = []
        for pid in self.tree.get_children(""):
            if pid in self.projects:
                rows.append(("project", self.projects[pid]["name"], pid))
                rows.extend(self._grid_rows(pid))
        return rows

    def _grid_rows(self, parent, path=()):
        """The list's rows for the thumbnail grid: folders as headings, then their images. A heading is the folder's whole path in the
        experiment ("vgf62 ctb in ipac / s1"): the grid is flat, and folders of the same name ("s1", "s2") are common."""
        rows = []
        for item in self.tree.get_children(parent):
            if item in self.entries:
                tag = "bad" if self.entries[item].problem else "cloud" if item in self.item_cloud else ""
                rows.append(("image", item, self.entries[item].name, tag))
            else:
                folder = path + (self.tree.item(item, "text"),)
                rows.append(("folder", " / ".join(folder)))
                rows.extend(self._grid_rows(item, folder))
        return rows

    def _start_thumbnails(self):
        """Make the thumbnails in a background thread, in list order. Only images stored on this computer are read (and only every n-th
        row of them); online-only ones get a placeholder. Made once per experiment, and cached on disk."""
        if self.thumbs_started or not self.entries:
            return
        self.thumbs_started = True
        done = self.thumb_grid.thumbs                # (thumbnails already made, kept when another project was opened or closed)
        gen, todo = self.thumb_gen, [(i, e) for i, e in self.entries.items() if done.get(i, (None, None))[1] != "ok"]
        self.thumbs_busy = True

        def work():
            try:
                make_all()
            finally:
                if gen == self.thumb_gen:
                    self.thumbs_busy = False

        def make_all():
            for item, entry in todo:
                # An image being opened comes first: wait (reading nothing) until it is shown, sharp level included. While the list
                # is showing instead of the thumbnails, wait too: the background loader goes first then.
                while not self.thumbs_go.wait(0.2) or not self.thumbs_shown:
                    if gen != self.thumb_gen:
                        return
                    if not self.thumbs_shown:
                        time.sleep(0.2)
                if gen != self.thumb_gen:
                    return
                try:
                    image, state = make_thumbnail(entry)
                except Exception:
                    image, state = None, "none"
                self.results.put(("thumb", gen, (item, image, state)))

        threading.Thread(target=work, daemon=True).start()

    def _start_preload(self):
        """Load each image's first picture (the pyramid level that fills the window, as when it is opened) in the background, in list
        order, after the thumbnails (while they are showing; before them while the list is), and keep it with the recently viewed images: picking it then shows it at once, and only its sharper
        level is still read. It only fills memory that is free (the share set in Settings: see kept_limit), and stops when that is
        used up, saying so in the status bar. An image being opened comes first: the loader pauses, mid-file if need be, until that image is shown, then goes on where it stopped. Images
        whose files would need asking (large and online-only: see needs_asking) are skipped."""
        gen, todo = self.thumb_gen, list(self.entries.items())
        cw, ch = self.canvas_size()

        def wait(entry):
            """Wait while thumbnails are made (when they are showing) or an image is opened; "stop" (another experiment), "skip"
            (this image needs no preloading any more: it is open, or kept) or "go"."""
            while (self.thumbs_busy and self.thumbs_shown) or not self.thumbs_go.is_set():
                if gen != self.thumb_gen:
                    return "stop"
                time.sleep(0.05)
            if gen != self.thumb_gen:
                return "stop"
            if entry in self.kept_images or (self.img is not None and self.img.entry is entry):
                return "skip"
            return "go"

        def work():
            for item, entry in todo:
                state = wait(entry)
                if state == "stop":
                    return
                if state == "skip" or entry.problem or not entry.width or entry.nz != 1:
                    continue
                try:
                    img = OpenImage(entry)
                    level = img.quick_level(max(entry.width, entry.height) * min(cw / entry.width, ch / entry.height)) or img.pyramid
                    n = len(entry.channels)
                    if level is not None:
                        files, step, pages = [level[0].files[(0, c)] for c in range(n)], level[1], level[0].pages
                        size = level[0].width * level[0].height * 2 * n
                    else:
                        files, step, pages = [entry.files[(0, c)] for c in range(n)], img.ov_step, entry.pages
                        size = -(-entry.width // step) * -(-entry.height // step) * 2 * n
                    if any(needs_asking(f) for f in files):
                        continue                        # (large online-only files: only after asking, when opened)
                    if self.kept_bytes + size > self.kept_limit():
                        self.results.put(("preload_full", gen, [e for i, e in todo[todo.index((item, entry)):] if e.nz == 1]))
                        return                          # (memory for kept images is full)
                    arrays, last = [], -1
                    for c, f in enumerate(files):
                        def progress(frac, c=c):
                            nonlocal last
                            pct = int(100 * (c + frac) / n)
                            if pct != last:
                                last = pct
                                self.results.put(("preload_progress", gen, (item, pct)))
                        a = read_rows(f, 1 if level is not None else step, lambda: wait(entry), progress, page=pages.get((0, c)))
                        if not isinstance(a, np.ndarray):
                            break                       # ("stop" or "skip")
                        arrays.append(pyramid_to_original(a, entry.channels[c].max_value) if level is not None and a.dtype == np.uint8
                                      else a)
                    if len(arrays) == n:
                        self.results.put(("preloaded", gen, (item, entry, img, arrays, step)))
                    else:
                        self.results.put(("preload_progress", gen, (item, None)))
                        if a == "stop":
                            return
                except Exception:
                    self.results.put(("preload_progress", gen, (item, None)))

        threading.Thread(target=work, daemon=True).start()

    def _preloaded(self, item, entry, img, arrays, step):
        """A picture loaded in the background: kept with the recently viewed images (as the first to go), unless it is no longer
        needed or there is no room."""
        self.preload_progress.pop(item, None)
        self._label_item(item)
        self.thumb_grid.set_progress(item, None)
        if entry in self.kept_images or (self.img is not None and self.img.entry is entry):
            return
        img.overviews = {(0, c): a for c, a in enumerate(arrays)}
        img.ov_step = step
        if self.kept_bytes + img.nbytes() > self.kept_limit():
            return
        self.kept_images[entry] = img
        self.kept_images.move_to_end(entry, last=False)         # (dropped before any image that was actually looked at)
        self._kept_changed()

    def _kept_changed(self):
        self.kept_bytes = sum(k.nbytes() for k in self.kept_images.values())
        self._show_loaded()

    def _loaded_entries(self):
        """The images with a picture in memory: kept (recently viewed, or preloaded) or open."""
        return set(self.kept_images) | ({self.img.entry} if self.img is not None else set())

    def _show_loaded(self):
        """Dim or undim the names of images whose picture came into memory or left it since the last call (see _label_item)."""
        now = self._loaded_entries()
        for entry in now ^ getattr(self, "shown_loaded", set()):
            item = self._item_of(entry)
            if item is not None:
                self._label_item(item)
        self.shown_loaded = now

    def _start_slide(self, pid):
        """Show the slide overview of project `pid` above the list (the project of the image selected), when it has one: outlines at
        once, the picture of the slide made in the background (like a thumbnail: from LAS X's pyramid if possible, never from
        online-only files; cached)."""
        self.slide_gen += 1
        self.slide_project = pid
        layout = slide_layout(self._entries_of(pid)) if pid else None
        if layout is None:
            self.slide.pack_forget()
            return
        ov, extent, boxes = layout
        self.slide_entry = ov
        self.slide.set_slide(extent, boxes)
        self.slide.pack(fill=tk.X, padx=4, pady=(4, 0), before=self.list_bar)
        sel = self.tree.selection()
        if sel:
            self.slide.set_current(sel[0])
        self._make_slide_picture()

    def _make_slide_picture(self):
        """The picture of the slide (self.slide_entry), made in the background, with its contrast and color maps as they are now."""
        ov, gen = self.slide_entry, self.slide_gen

        def work():
            try:
                image, _ = make_thumbnail(ov, size=2 * SLIDE_HEIGHT)
            except Exception:
                image = None
            self.results.put(("slide", gen, image))

        threading.Thread(target=work, daemon=True).start()

    def select_item(self, item):
        """A click (or arrow key) in the thumbnail grid: select that image in the list, which opens it."""
        if self.tree.exists(item):
            self.tree.selection_set(item)
            self.tree.see(item)
            self.thumb_grid.set_selected(item)

    def _set_open_progress(self, progress):
        """(item, percent) of the image opened, or None once it is read (or left): shown after its name in the list."""
        old, self.open_progress = self.open_progress, progress
        for item in {p[0] for p in (old, progress) if p is not None}:
            self._label_item(item)
            if item not in self.preload_progress:
                self.thumb_grid.set_progress(item, progress[1] / 100 if progress is not None and progress[0] == item else None)

    def _close_image(self):
        self.generation += 1
        self._set_open_progress(None)
        if self.img:
            if not self.loading and self.img.complete(self.z):
                self._keep_image(self.img)
            else:
                self.img.close()
        self.img = None
        self._show_loaded()
        self.photo = None
        self.flip_one.set(False)
        self.flip_one_check.configure(state=tk.DISABLED)        # (no section shown)
        self.last_view_rgb = None
        self.shrunk.clear()
        self.selected = self.ann_drag = self.hover_ann = None
        self._each_canvas(lambda: self.canvas.delete("all"))
        self.cursor_label.config(text="")
        self._show_cursor_values([])
        self.canvas_message = ""
        self.detail_pending = False
        self.detail_cache = None
        self.thumbs_go.set()            # (no image is loading now)
        self._update_annotations_check()

    def _keep_image(self, img):
        """Keep an image that is being left, fully loaded, so that coming back to it needs no reading: its pictures (overviews) stay in
        memory, newest last, up to kept_limit in all (the longest unseen go first); its full-resolution files are closed (they
        are opened again if it is zoomed into)."""
        img.planes.clear()
        self.kept_images[img.entry] = img
        self._keep_within_limit()

    def _keep_within_limit(self):
        """Drop kept images, the longest unseen first, until they fit in the memory for previews."""
        while sum(k.nbytes() for k in self.kept_images.values()) > self.kept_limit() and self.kept_images:
            self.kept_images.popitem(last=False)[1].close()
        self._kept_changed()

    def forget_kept_images(self):
        for img in self.kept_images.values():
            img.close()
        self.kept_images.clear()
        self._kept_changed()

    def on_select(self, event=None):
        sel = self.tree.selection()
        if sel:
            self.thumb_grid.set_selected(sel[0])
            project = self._project_of(sel[0])
            if project != self.slide_project:
                self._start_slide(project)
            self.slide.set_current(sel[0])
        entry = self.entries.get(sel[0]) if sel else None
        if entry is None or (self.img and self.img.entry is entry):
            return
        if entry.problem:
            self._close_image()
            self.show_message(f"{entry.name} can't be shown.\n{entry.problem}")      # (in the image area too: not a blank)
            self.status.config(text=f"{entry.name}: {entry.problem}")
            return
        ask = [] if entry in self.kept_images else self._files_to_ask(entry)
        if ask:
            size = size_text(sum(os.path.getsize(f) for f in ask))
            if not messagebox.askyesno("Download from OneDrive?",
                                       f"{len(ask)} of {entry.name}'s files ({size}) are online-only.\n\n"
                                       "Viewing the image makes OneDrive download them in full. Continue?", parent=self.root):
                self.status.config(text=f"{entry.name} not opened (online-only).")
                return
        self.show_entry(entry, full_res_ok=True if ask else None)

    def _files_to_ask(self, entry):
        """The large online-only files (see needs_asking) that showing `entry` would read: none when its pyramid can be used (its
        full-resolution files are then read only when zooming in, and asked about then), else its own."""
        if entry.nz == 1 and entry.pyramid and OpenImage(entry).pyramid is not None:
            return []
        return [f for f in entry.all_files() if needs_asking(f)]

    def show_entry(self, entry, full_res_ok=None):
        """Show `entry`. full_res_ok: may its large online-only full-resolution files be downloaded when zooming in needs them (True:
        agreed to already; None: ask then; False: no)."""
        self.full_res_ok = full_res_ok
        same_channels = (len(entry.channels) == len(self.ch_settings) and self.img is not None
                         and [c.name for c in entry.channels] == [c.name for c in self.img.entry.channels])
        shown = [st["on"] for st in self.ch_settings] if same_channels else [True] * len(entry.channels)
        kept = self.kept_images.pop(entry, None)        # (seen recently, or loaded in the background: its picture is in memory)
        self._kept_changed()
        self._close_image()
        self.z = min(self.z, entry.nz - 1)
        instant = kept is not None and kept.has_all(self.z)
        if instant:
            self.img = kept
        else:
            if kept is not None:
                kept.close()
            n = len(entry.channels)
            self.show_message(f"Loading channel 1 of {n}..." if n > 1 else "Loading...")
            self.thumbs_go.clear()          # (the thumbnails wait until this image is shown: see _start_thumbnails)
            self.img = OpenImage(entry)
        self._show_loaded()
        # Each image opens with its own LAS X contrast; which channels are switched on carries over to images of the same channels.
        self.ch_settings = [{"on": on, "lo": c.black, "hi": c.white} for c, on in zip(entry.channels, shown)]
        self.flip_one.set(self._section_flipped(entry))
        self.flip_one_check.configure(state=tk.NORMAL)
        self._build_channel_controls()
        self._update_annotations_check()
        if entry.nz > 1:
            self.z_scale.config(to=entry.nz - 1)
            self.z_var.set(self.z)
            self.z_frame.pack(fill=tk.X, padx=4, before=self.ch_frame)
        else:
            self.z_frame.pack_forget()
        self.fit(redraw=False)
        if instant and self.img.complete(self.z):
            self.loading = False
            self.redraw()
            self._update_cursor()
        elif instant:                   # (its first picture was loaded in the background: shown now, the sharper level read next)
            self.loading = True
            self.redraw(preview_only=True)
            self.thumbs_go.clear()
            self._load_overviews()
        else:
            self._load_overviews()

    def _load_overviews(self):
        """Read the overview of each channel in the background (the first read of a file is the slow part), showing each as it comes."""
        img, z, gen = self.img, self.z, self.generation
        self.loading = True
        self._update_cursor()
        n = len(img.entry.channels)

        def on_disk():
            """For the status bar once it is read (asked of the files there, in the background: a round trip each on a network drive)."""
            mb = (img.entry.full_bytes if img.entry.full_bytes is not None else
                  sum(os.path.getsize(f) for f in img.entry.all_files() if os.path.exists(f))) / 1e6
            return f"{mb:,.0f} MB on disk"
        # With a pyramid, a level about the size of the fit view comes first (a quarter of the data, or less); the overview's own,
        # sharper level is read next and replaces it whole, once all its channels are in.
        quick = img.quick_level(max(img.entry.width, img.entry.height) * self.fit_scale)
        have_first = img.pyramid is not None and img.has_all(z)     # (first picture already loaded in the background: only sharpen)
        item = self._item_of(img.entry)
        last = [-1]

        def report(fraction):
            """How far the preview has got, after the name in the list ("45% preview"), as the background loader shows it."""
            pct = int(100 * fraction)
            if pct != last[0]:
                last[0] = pct
                self.results.put(("open_progress", gen, (item, pct)))

        def work():
            try:
                read(quick)
            except Exception as e:       # (always report: the thumbnails wait for this image)
                self.results.put(("error", gen, f"Can't read {img.entry.name}: {e}"))

        def read(quick):
            if quick is not None or have_first:
                # (the list's "% preview" is for the first picture only; the sharper level after it is read silently)
                first = have_first
                if not have_first:
                    report(0)
                    for c in range(n if quick is not None else 0):
                        if gen != self.generation:
                            return
                        self.results.put(("progress", gen, f"Reading {img.entry.name}: channel {c + 1} of {n} ..."))
                        a = img._from_pyramid(c, quick)
                        if a is None:
                            break
                        self.results.put(("quick", gen, (z, c, a, quick[1])))
                        report((c + 1) / n)
                    else:
                        first = quick is not None
                if first:
                    self.results.put(("open_progress", gen, None))
                sharp = []
                for c in range(n):
                    if gen != self.generation:
                        return
                    self.results.put(("progress", gen, f"Reading {img.entry.name}: sharper view, channel {c + 1} of {n} ..."))
                    sharp.append(img._from_pyramid(c))
                    if not first:
                        report((c + 1) / n)
                if all(a is not None for a in sharp):
                    self.results.put(("sharp", gen, (z, sharp, img.pyramid[1])))
                self.results.put(("loaded", gen, on_disk()))
                return
            for c in range(n):
                if gen != self.generation:
                    return
                self.results.put(("progress", gen, f"Reading {img.entry.name}: channel {c + 1} of {n} ..."))
                try:
                    img.overview(z, c, lambda f, c=c: report((c + f) / n))
                except Exception as e:
                    self.results.put(("error", gen, f"Can't read channel {c + 1} of {img.entry.name}: {e}"))
                    return
                self.results.put(("overview", gen, None))
            self.results.put(("loaded", gen, on_disk()))

        threading.Thread(target=work, daemon=True).start()
        self.status.config(text=f"Reading {img.entry.name} ...")

    # --- channel controls

    def _build_channel_controls(self):
        for w in self.ch_frame.winfo_children():
            w.destroy()
        self.ch_widgets = []
        for i, (ch, st) in enumerate(zip(self.img.entry.channels, self.ch_settings)):
            # One column per channel: its swatch (its color and number; click to switch the channel on/off, grey when off), then its
            # black/white slider standing up, with the white level above and the black level below.
            f = ttk.Frame(self.ch_frame)
            f.pack(side=tk.LEFT, fill=tk.Y, padx=2, pady=2)
            color = self._swatch_color(i)
            on = tk.BooleanVar(value=st["on"])
            swatch = tk.Label(f, text=str(i + 1), width=3, font=("Segoe UI", 9, "bold"), cursor="hand2", relief=tk.RAISED, bd=1)
            swatch.pack(pady=(0, 2))
            swatch.bind("<Button-1>", lambda e, i=i: self.toggle_channel(i))
            swatch.bind("<Button-3>", lambda e, i=i: self.pick_channel_color(i))
            tip = ((ch.name + "\n" if ch.named else "") + f"Click to show or hide this channel (or press {i + 1}).\n"
                   "Right-click to change its color.")
            Tooltip(swatch, tip)
            slider = RangeSlider(f, ch.max_value, st["lo"], st["hi"], color, command=lambda lo, hi, i=i: self._set_lo_hi(i, lo, hi),
                                 vertical=True)
            slider.pack(fill=tk.Y, expand=True)
            slider.bind("<Button-3>", lambda e, i=i: self.pick_channel_color(i))
            slider.bind("<Double-Button-1>", lambda e, i=i, s=slider: (s.undo_press(), self.pick_channel_color(i)))
            swatch.color = color
            swatch.keep_colors = True       # (styled by _show_channel_state)
            self.ch_widgets.append((on, slider, swatch))
            self._show_channel_state(i)
        if self.split:
            self._layout_panels()           # (this image's channels, each in a panel)

    def _show_channel_state(self, i):
        """A channel's swatch in its color (number in black or white, whichever reads), or grey when it is switched off; its slider too."""
        on, slider, swatch = self.ch_widgets[i]
        rgb = self.img.entry.channels[i].color
        light = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] > 0.5
        if on.get():
            swatch.configure(bg=swatch.color, fg="black" if light else "white", relief=tk.RAISED)
        else:
            swatch.configure(bg=PAL["off_bg"], fg=PAL["off_fg"], relief=tk.SUNKEN)
        slider.set_active(on.get())

    def _swatch_color(self, i):
        """The channel's color as the swatch and slider show it (a little muted, as LAS X's colors always were here)."""
        rgb = self.img.entry.channels[i].color
        return "#%02x%02x%02x" % tuple(int(200 * k + 30) for k in rgb)

    def pick_channel_color(self, i):
        """Right-click on a channel's swatch or slider: choose its LAS X colour map (see LutChooser)."""
        if self.img and i < len(self.img.entry.channels):
            LutChooser(self, i)

    def _set(self, i, key, value):
        if i >= len(self.ch_settings) or self.ch_settings[i][key] == value:
            return
        self.ch_settings[i][key] = value
        self.schedule_redraw()

    def _set_lo_hi(self, i, lo, hi):
        """A slider moved: the image's contrast changed (an unsaved change: Save writes it into its .xlif; Undo takes back the whole
        drag at once)."""
        if i < len(self.ch_settings) and (self.ch_settings[i]["lo"], self.ch_settings[i]["hi"]) != (lo, hi):
            step = (id(self.img.entry), i, getattr(self.ch_widgets[i][1], "presses", 0))
            self._contrast_change([(i, lo, hi)], new_step=step != self.contrast_step)
            self.contrast_step = step

    def _contrast_change(self, changes, new_step=True):
        """Set the shown image's contrast: changes = [(channel, black, white)], kept in its channels (so it stays when coming back to
        the image) and remembered for Undo as one step (new_step False: part of the step before, as a slider drag)."""
        e = self.img.entry
        changes = [(i, lo, hi) for i, lo, hi in changes if (e.channels[i].black, e.channels[i].white) != (lo, hi)]
        if not changes:
            return
        if new_step or not self.history or self.history[-1][0] != "range":
            self.history = self.history[-199:] + [("range", [(e, i, e.channels[i].black, e.channels[i].white) for i, _, _ in changes])]
        else:                           # (the same drag: the step keeps the values from before it)
            have = {i for _, i, _, _ in self.history[-1][1]}
            self.history[-1][1].extend((e, i, e.channels[i].black, e.channels[i].white) for i, _, _ in changes if i not in have)
        for i, lo, hi in changes:
            e.channels[i].black, e.channels[i].white = lo, hi
            self.ch_widgets[i][1].set(lo, hi)
            self.ch_settings[i]["lo"], self.ch_settings[i]["hi"] = lo, hi
        self.schedule_redraw()
        if self.contrast_job is not None:
            self.root.after_cancel(self.contrast_job)
        self.contrast_job = self.root.after(400, lambda: self._contrast_settled({e}))    # (not on every step of a drag)

    def _contrast_settled(self, changed):
        """After contrast of these images changed (and the slider rests): the list's "unsaved" marks, Save, their thumbnails."""
        self.contrast_job = None
        for e in changed:
            self._label_item(self._item_of(e))
        self._update_annotations_check()
        self._refresh_thumbnails([self._item_of(e) for e in changed])

    def toggle_channel(self, i):
        if self.img and i < len(self.ch_settings) and not isinstance(self.root.focus_get(), tk.Entry):
            on = self.ch_widgets[i][0]
            on.set(not on.get())
            self._show_channel_state(i)
            self._set(i, "on", on.get())
            if self.split:
                self._layout_panels()       # (a panel more, or one less)

    def _auto_range(self, i):
        """Contrast from the overview: 0.1-99.9 percentile of the non-zero pixels (stitched scans have zero outside the tiles); or None."""
        if not self.img or (self.z, i) not in self.img.overviews:
            return None
        a = self.img.overviews[(self.z, i)]
        v = a[a > 0]
        if v.size == 0:
            return None
        lo, hi = np.percentile(v, [0.1, 99.9])
        return i, int(lo), max(int(hi), int(lo) + 1)

    def auto(self, i):
        r = self._auto_range(i)
        if r:
            self._contrast_change([r])

    def auto_all(self):
        if self.img:
            self._contrast_change([r for r in map(self._auto_range, range(len(self.ch_settings))) if r])

    def lasx_all(self):
        """Back to the contrast in the image's LAS X file (as last saved)."""
        if self.img:
            self._contrast_change([(i, ch.saved_black, ch.saved_white) for i, ch in enumerate(self.img.entry.channels)])

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
        self.pixel_grid = False
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
        sx = self._fx(sx)
        before = (self.scale, self.cx, self.cy)
        self.scale = scale
        self.clamp_view()
        if self.scale < MAX_ZOOM - 1e-9:
            self.pixel_grid = False
        self.cx = ix - (sx - cw / 2) / self.scale
        self.cy = iy - (sy - ch / 2) / self.scale
        self.clamp_view()
        if np.allclose((self.scale, self.cx, self.cy), before, rtol=0, atol=1e-6):
            return          # (already at the zoom limit: a redraw would only blink, blurry then sharp again)
        self.schedule_redraw(0)

    def screen_to_image(self, sx, sy):
        cw, ch = self.canvas_size()
        sx = self._fx(sx)
        return self.cx + (sx - cw / 2) / self.scale, self.cy + (sy - ch / 2) / self.scale

    def _fx(self, x):
        """Flip left-right: canvas x <-> where it would be unflipped (the view is worked out unflipped, and mirrored when drawn)."""
        return self.canvas_size()[0] - x if self._flipped() else x

    def _view_corners(self):
        """The image points at the canvas's top-left and bottom-right corners, as if unflipped (ix0 <= ix1 either way)."""
        cw, ch = self.canvas_size()
        return (self.cx - cw / 2 / self.scale, self.cy - ch / 2 / self.scale,
                self.cx + cw / 2 / self.scale, self.cy + ch / 2 / self.scale)

    def on_flip_toggle(self):
        self.settings["flip_left_right"] = self.flip_x.get()
        if self.img:
            self.redraw()

    @staticmethod
    def _flip_key(entry):
        return os.path.normcase(os.path.normpath(entry.xlif_path))

    def _flip_in_files(self):
        """Settings: "Flip this section" is saved into the section's .xlif (see write_flip), not kept on this computer."""
        return self.settings.get("flip_in_files", True)

    def _flip_note(self, entry):
        """This computer's note of the section's flip, by its .xlif's path (in the settings: "flipped_sections", "unflipped_sections"):
        True, False, or None (none). Kept while "Flip this section" isn't saved into the files, or by an earlier version."""
        key = self._flip_key(entry)
        return (True if key in self.settings.get("flipped_sections", []) else
                False if key in self.settings.get("unflipped_sections", []) else None)

    def _set_flip_note(self, entry, flipped):
        """Note the section's flip on this computer (None: forget the note)."""
        key = self._flip_key(entry)
        for name in ("flipped_sections", "unflipped_sections"):
            if key in self.settings.get(name, []):
                self.settings[name] = [k for k in self.settings[name] if k != key]
        if flipped is not None:
            self.settings.setdefault("flipped_sections" if flipped else "unflipped_sections", []).append(key)

    def _section_flipped(self, entry):
        """Whether the section is shown flipped. Saved into the files (the setting): as in its .xlif, or as changed here and not saved
        yet; a note on this computer (from an earlier version, or from when the setting was off) counts as such a change, until Save
        writes it into the file. Else: this computer's note, or, without one, the .xlif's."""
        note = self._flip_note(entry)
        if not self._flip_in_files():
            return entry.flipped if note is None else note
        if note is not None and note == entry.saved_flipped:
            self._set_flip_note(entry, None)        # (the file says so already)
        elif note is not None and note != entry.flipped:
            entry.flipped = note
            self._label_item(self._item_of(entry))
            self._update_annotations_check()
        return entry.flipped

    def on_flip_one_toggle(self):
        """"Flip this section": written into its .xlif by Save (until then, marked unsaved), or, with that setting off, noted on this
        computer at once."""
        if not self.img:
            return
        e = self.img.entry
        if self._flip_in_files():
            e.flipped = self.flip_one.get()
            self._set_flip_note(e, None)
            self._label_item(self._item_of(e))
            self._update_annotations_check()
        else:
            self._set_flip_note(e, None if self.flip_one.get() == e.flipped else self.flip_one.get())
        self.redraw()

    def set_flip_in_files(self, on):
        """The setting: "Flip this section" saved into the files, or kept on this computer. Turned off, flips not saved yet are noted
        on this computer instead (and are no longer unsaved); turned on, the notes become unsaved changes (see _section_flipped)."""
        self.settings["flip_in_files"] = on
        if not on:
            for item, e in self.entries.items():
                if e.flip_changed:
                    self._set_flip_note(e, e.flipped)
                    e.flipped = e.saved_flipped
                    self._label_item(item)
        self._save_settings()
        if self.img:
            self.flip_one.set(self._section_flipped(self.img.entry))
            self.redraw()
        self._update_annotations_check()

    def _flipped(self):
        """The view is mirrored: "Flip all sections" or "Flip this section", not both (both flip it back)."""
        return self.flip_x.get() != (self.img is not None and self.flip_one.get())

    def show_message(self, text):
        """Words in the middle of the image area ("Loading...", ...) until a picture replaces them (any drawing clears the canvas)."""
        self.canvas_message = text
        self._each_canvas(lambda: self.canvas.delete("all"))
        if text:
            cw, ch = self.canvas_size()
            self.canvas.create_text(cw / 2, ch / 2, text=text, fill="#b0b0b0", font=("Segoe UI", 16), tags="message",
                                    width=max(cw - 60, 200), justify=tk.CENTER)

    def on_resize(self):
        if self.photo is None and getattr(self, "canvas_message", ""):
            self.show_message(self.canvas_message)          # (kept in the middle)
        if self.img:
            at_fit = abs(self.scale - self.fit_scale) < 1e-9
            self.clamp_view()
            if at_fit:
                self.fit(redraw=False)
            self.schedule_redraw()

    def on_wheel(self, event, delta=None):
        d = event.delta if delta is None else delta
        if d > 0 and self.img is not None and self.scale >= MAX_ZOOM - 1e-9:
            if not self.pixel_grid:             # (at the largest zoom already: show the pixels' edges instead)
                self.pixel_grid = True
                self.schedule_redraw(0)
            return
        if d < 0 and self.pixel_grid:           # (back from the grid: the same zoom, smooth again)
            self.pixel_grid = False
            self.schedule_redraw(0)
            return
        self._moving()
        self.zoom_to(self.scale * (1.25 ** (d / 120)), event.x, event.y)

    def on_press(self, event):
        self.canvas.focus_set()
        if self.ann_press(event):
            return
        self.pan_start(event)

    def on_drag(self, event):
        if self.ann_drag:
            self.ann_motion(event)
        else:
            self.pan_move(event)

    def on_release(self, event):
        if self.ann_drag:
            self.ann_release(event)
        else:
            self.pan_end(event)

    def _hover_handles(self, sx, sy):
        """The mouse at (sx, sy), or off the image (None): the hand pointer over a handle, the move pointer over an outline, and the
        handles of the shape or line under the mouse shown (hollow) when it isn't the selected one."""
        on = sx is not None and self.img is not None and bool(self.img.entry.pixel_um) and self.show_annotations.get()
        grab = self._handle_at(sx, sy) if on else None
        hover = grab[1] if grab else self._hover_target(sx, sy) if on else None
        if hover is self.selected:
            hover = None
        over = None if grab is None else "move" if grab[0] == "move" else "handle"
        if over != getattr(self, "over_handle", None):
            self.over_handle = over
            self._update_cursor()
        if hover is not getattr(self, "hover_ann", None):
            self.hover_ann = hover
            self._draw_annotations()

    def pan_start(self, event):
        self.drag = (event.x, event.y, self.cx, self.cy)

    def pan_move(self, event):
        if not self.drag or not self.img:
            return
        x0, y0, cx0, cy0 = self.drag
        self._moving()
        self.cx = cx0 - (event.x - x0) / self.scale * (-1 if self._flipped() else 1)
        self.cy = cy0 - (event.y - y0) / self.scale
        self.clamp_view()
        self.redraw(preview_only=True)

    def pan_end(self, event):
        if self.drag:
            self.drag = None
            self.schedule_redraw(0)

    def on_motion(self, event):
        if not self.img:
            return
        self._hover_handles(event.x, event.y)
        e = self.img.entry
        ix, iy = self.screen_to_image(event.x, event.y)
        if not (0 <= ix < e.width and 0 <= iy < e.height):
            self.cursor_label.config(text="")
            self._show_cursor_values([])
            return
        x, y = int(ix), int(iy)
        parts = [f"x {x}, y {y} px"]
        if e.pixel_um:
            parts.append(f"({x * e.pixel_um:,.1f}, {y * e.pixel_um:,.1f} µm)")
        vals = []
        for c, ch in enumerate(e.channels):
            p = self.img.planes.get((self.z, c))     # only planes already open: hovering must not open files
            if p is not None and not self.loading:
                vals.append((c, f"{ch.name}: {int(p[y, x])}"))
        self.cursor_label.config(text="   ".join(parts))
        self._show_cursor_values(vals)

    def _show_cursor_values(self, vals):
        """After the pointer's position: each channel's value under it, "DAPI: 1234", in the channel's color (made lighter or darker
        when it would hardly show on the theme's background: see readable_color)."""
        labels = self.cursor_values
        while len(labels) < len(vals):
            lb = tk.Label(self.cursor_row, text="", anchor=tk.W, padx=6, bd=0, font=tkfont.nametofont("TkDefaultFont"))
            lb.keep_colors = True
            labels.append(lb)
        bg = self.root.winfo_rgb(PAL["bg"])
        bg = tuple(v / 65535 for v in bg)
        for k, lb in enumerate(labels):
            if k < len(vals):
                c, text = vals[k]
                lb.configure(text=text, fg=readable_color(self.img.entry.channels[c].color, bg), bg=PAL["bg"])
                if not lb.winfo_manager():
                    lb.pack(side=tk.LEFT)
            elif lb.winfo_manager():
                lb.pack_forget()

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
        ix0, iy0, ix1, iy1 = self._view_corners()
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
            self._each_canvas(lambda: self.canvas.delete("all"))
            self._draw_annotations()
            self._draw_scale_bar()
            self._draw_load_note()
            return
        step = max(1, int(1 / self.scale))
        need = [c for c in range(n) if luts[c] is not None]           # (the channels to read for a sharp view, if one is needed)
        detail = not preview_only and not self.loading and step < s
        if detail:
            # The sharp view's pixels for this very view are still here (only the contrast, or a channel, changed): recolor them
            # at once instead of drawing the blurrier overview first and reading the files again (that blinked).
            (i0, j0, i1, j1), place = self._view_rect(step)
            key = (id(self.img), self.z, i0, j0, i1, j1, step, place)
            cache = self.detail_cache
            if cache is not None and cache["key"] != key and self._kept_covers(step, need):
                # The new view lies within the sharp pixels already here, at the detail it needs (zooming in, or panning within
                # them): draw it from those; reading the files again would bring the same pixels.
                overlay = self._detail_overlay(luts)
                if overlay is not None:
                    self._put(*overlay)
                    self.detail_pending = False
                    self._update_cursor()
                    return
            if cache is not None and cache["key"] == key:
                arrays = [a if luts[c] is not None else None for c, a in enumerate(cache["arrays"])]
                if any(a is not None for a in arrays):
                    self._show(self._mix(arrays, luts), *place, smooth=self.scale * step < 1)
                need = [c for c in need if cache["arrays"][c] is None]
                if not need:
                    self.detail_pending = False
                    self._update_cursor()
                    return
                ovs = [None] * n                  # (skip the overview: the sharp view stays while the missing channel is read)
        if any(o is not None for o in ovs):
            (i0, j0, i1, j1), place = self._view_rect(s)
            if i1 > i0 and j1 > j0:
                crop = [o[j0:j1, i0:i1] if o is not None else None for o in ovs]
                rows, cols = next(a.shape for a in crop if a is not None)
                left, top, w, h = place
                if (cols, rows) != (i1 - i0, j1 - j0):      # (a pyramid level drops the last partial block at the right and bottom edges)
                    w, h = max(1, round(cols * s * self.scale)), max(1, round(rows * s * self.scale))
                if cols > w and rows > h:
                    # Zoomed out: shrink each channel to the screen size first (averaging, Pillow handles 16-bit), so that the slow part,
                    # mixing the channels, works on only as many pixels as are shown. The shrunk channels are kept for this view, so a
                    # redraw for a slider, or for the next channel arriving, shrinks nothing again.
                    view = (id(self.img), self.z, i0, j0, i1, j1, w, h)
                    for c, a in enumerate(crop):
                        if a is not None:
                            key = view + (c,)
                            if key not in self.shrunk:
                                self.shrunk[key] = np.asarray(Image.fromarray(np.ascontiguousarray(a)).resize((w, h), Image.BOX))
                            crop[c] = self.shrunk[key]
                    for key in [k for k in self.shrunk if k[:-1] != view]:
                        del self.shrunk[key]
                self._show(self._mix(crop, luts), left, top, w, h, smooth=self.scale * s < 1, overlay=self._detail_overlay(luts))
        if not detail:
            if not preview_only:
                self.detail_pending = False
                self._update_cursor()
            return
        if not self._full_res_allowed():
            return
        self.generation += 1
        self.detail_pending = True           # (an hourglass until the sharper view is in)
        self._update_cursor()
        gen, img, z = self.generation, self.img, self.z
        (i0, j0, i1, j1), place = self._view_rect(step)
        key = (id(img), z, i0, j0, i1, j1, step, place)
        smooth = self.scale * step < 1

        def work():
            try:
                arrays = [None] * n
                for c in need:
                    if gen != self.generation:
                        return
                    arrays[c] = np.ascontiguousarray(img.plane(z, c)[j0 * step:j1 * step:step, i0 * step:i1 * step:step])
                self.results.put(("detail", gen, (key, arrays, place, smooth)))
            except Exception as e:
                self.results.put(("error", gen, f"Can't read {img.entry.name}: {e}"))

        threading.Thread(target=work, daemon=True).start()

    def _moving(self):
        """The view is being zoomed or dragged: until it rests for MOVING_PAUSE, enlarged pixels are drawn the fast way (bilinear),
        then the view is drawn once more, the way chosen in Settings (Mitchell-Netravali takes 40-90 ms for a large view: too slow
        for every wheel click)."""
        self.moving_until = time.monotonic() + MOVING_PAUSE
        if self.rest_job is not None:
            self.root.after_cancel(self.rest_job)
        self.rest_job = self.root.after(int(MOVING_PAUSE * 1000) + 20, self._at_rest)

    def _at_rest(self):
        self.rest_job = None
        if self.img is not None and self.enlarge_method in _KERNELS and not self.drag:
            self.redraw()                   # (from the pixels already here: nothing is read again)

    def _draw_method(self):
        """How zoomed-in pixels are enlarged: as chosen in Settings (bilinear while moving, for speed); with the pixel grid shown, as
        squares of one color each, filling the grid's cells."""
        if self.pixel_grid:
            return "nearest"
        return "bilinear" if time.monotonic() < self.moving_until and self.enlarge_method in _KERNELS else self.enlarge_method

    def _mix_sig(self):
        """What the mixed colours depend on: each channel on/off, its levels and its colour map, and the panels (see _mix)."""
        return (tuple((s["on"], s["lo"], s["hi"]) for s in self.ch_settings), tuple(c.lut for c in self.img.entry.channels),
                tuple(self.panel_channels))

    def _detail_overlay(self, luts):
        """While a new view waits for its sharp pixels: the sharp pixels of the last sharp view, for the part of it still on screen,
        scaled to the new view: (PIL image, canvas x, canvas y), to lay over the zoomed-out picture; or None. Zooming in, they cover the
        whole new view; zooming out or panning, part of it. (Not when a channel shown now wasn't read for it: it would lack that color.)"""
        cache = self.detail_cache
        if cache is None or cache["key"][0] != id(self.img) or cache["key"][1] != self.z:
            return None
        _, _, i0, j0, i1, j1, st, _ = cache["key"]
        arrays = cache["arrays"]
        on = [c for c, lut in enumerate(luts) if lut is not None]
        if not on or any(arrays[c] is None for c in on):
            return None
        sig = self._mix_sig()
        if cache.get("sig") != sig:            # (mixed once per slider setting, not on every wheel click)
            cache["rgb"] = self._mix([a if luts[c] is not None else None for c, a in enumerate(arrays)], luts)
            cache["sig"] = sig
        rgbs = cache["rgb"]
        rows, cols = rgbs[0].shape[:2]
        x0, y0 = i0 * st, j0 * st                           # (image pixel of rgb[0, 0])
        cw, ch = self.canvas_size()
        vx0, vy0, vx1, vy1 = self._view_corners()
        a0, b0 = max(0, int((vx0 - x0) // st)), max(0, int((vy0 - y0) // st))
        a1, b1 = min(cols, int(-(-(vx1 - x0) // st))), min(rows, int(-(-(vy1 - y0) // st)))
        if a1 <= a0 or b1 <= b0:
            return None
        w, h = max(1, round((a1 - a0) * st * self.scale)), max(1, round((b1 - b0) * st * self.scale))
        ims = []
        for rgb in rgbs:                                    # (one per canvas)
            part = np.ascontiguousarray(rgb[b0:b1, a0:a1])
            ims.append(enlarge(part, w, h, self._draw_method()) if w > a1 - a0 or h > b1 - b0       # (zoomed in past the pixels)
                       else Image.fromarray(part).resize((w, h), Image.BILINEAR))
        return ims, round((x0 + a0 * st - self.cx) * self.scale + cw / 2), round((y0 + b0 * st - self.cy) * self.scale + ch / 2)

    def _kept_covers(self, step, channels):
        """The kept sharp pixels (detail_cache) cover all of the image that is in view, at `step` or finer, for these channels."""
        cache = self.detail_cache
        if cache is None or cache["key"][0] != id(self.img) or cache["key"][1] != self.z or cache["key"][6] > step:
            return False
        if any(cache["arrays"][c] is None for c in channels):
            return False
        _, _, i0, j0, i1, j1, st, _ = cache["key"]
        rows, cols = next(a.shape for a in cache["arrays"] if a is not None)
        e = self.img.entry
        vx0, vy0, vx1, vy1 = self._view_corners()
        vx0, vy0, vx1, vy1 = max(vx0, 0), max(vy0, 0), min(vx1, e.width), min(vy1, e.height)     # (the part of the image in view)
        x0, y0 = i0 * st, j0 * st
        return x0 <= vx0 and y0 <= vy0 and x0 + cols * st >= min(vx1, e.width - st) and y0 + rows * st >= min(vy1, e.height - st)

    def _show(self, rgbs, left, top, w, h, smooth, overlay=None):
        """rgbs: one RGB array per canvas (see _mix), all of this view."""
        ims = []
        for k, rgb in enumerate(rgbs):
            im = Image.fromarray(rgb)
            if im.size != (w, h):           # (smooth: shrinking; else zoomed in past the pixels: enlarged as chosen in Settings)
                im = im.resize((w, h), Image.BILINEAR) if smooth else enlarge(rgb, w, h, self._draw_method())
            if overlay is not None and k < len(overlay[0]):     # (sharper pixels on top, where there are some for this view)
                sharp, ox, oy = overlay[0][k], overlay[1], overlay[2]
                im.paste(sharp, (ox - left, oy - top))
            ims.append(im)
        self._put(ims, left, top)

    def _put(self, ims, left, top):
        """The finished pictures of the view (PIL images, one per canvas), at (left, top) on each, with the pixel grid, the channel's
        name (separated), the annotations, scale bar and notes. Flip left-right: mirrored here, and placed mirrored."""
        if self._flipped():
            ims = [im.transpose(Image.FLIP_LEFT_RIGHT) for im in ims]
            left = self.canvas_size()[0] - left - ims[0].width
        self.last_view = (ims, left, top, self.canvas_size(), self.panel_cols)
        self.last_view_rgb = ims[0]
        self.photos = [ImageTk.PhotoImage(im, master=self.root) for im in ims]
        self.photo = self.photos[0]

        def draw():
            c = self.canvas
            c.delete("all")
            if self.drawing_panel < len(self.photos):
                c.create_image(left, top, anchor=tk.NW, image=self.photos[self.drawing_panel])
            self._draw_pixel_grid()
            self._draw_panel_label()
            self._draw_annotations()
            self._draw_scale_bar()
            self._draw_load_note()
            self._draw_hires_note_here()      # (still waiting for the sharper view: the note stays over the new picture)
        self._each_canvas(draw)
        self._update_status()

    def _load_note(self):
        """While an image's channels are still arriving: "Loading channel 2 of 4..." (the picture is not complete yet); else ""."""
        img = self.img
        if img is None or not self.loading:
            return ""
        n = len(img.entry.channels)
        have = sum((self.z, c) in img.overviews for c in range(n))
        if have >= n:
            return ""
        return f"Loading channel {have + 1} of {n}..." if n > 1 else "Loading..."

    def _draw_load_note(self):
        self._each_canvas(self._draw_load_note_here)

    def _draw_load_note_here(self):
        """The note in the middle of the canvas; separated, on the panel of the channel being read (on the first panel if that channel
        is switched off), and "Pending..." on the panels of channels still to come (not empty: not read yet)."""
        c = self.canvas
        c.delete("loadnote")
        note = self._load_note()
        if not note:
            return
        if self.panel_channels != [None]:
            reading = next(k for k in range(len(self.img.entry.channels)) if (self.z, k) not in self.img.overviews)
            panel = self.panel_channels.index(reading) if reading in self.panel_channels else 0
            here = self.drawing_panel or 0
            if here != panel:
                if (self.z, self.panel_channels[here]) in self.img.overviews:
                    return
                note = "Pending..."
        cw, ch = self.canvas_size()
        t = c.create_text(cw / 2, ch / 2, text=note, fill="#e0e0e0", font=("Segoe UI", 16), tags="loadnote")
        x0, y0, x1, y1 = c.bbox(t)
        c.tag_lower(c.create_rectangle(x0 - 14, y0 - 8, x1 + 14, y1 + 8, fill="#202020", outline="#606060", tags="loadnote"), t)

    # --- annotations: showing and editing (saved into the image's .xlif only by the Save button)

    def on_annotations_toggle(self):
        self.settings["show_annotations"] = self.show_annotations.get()
        if not self.show_annotations.get():
            self.selected = None
            self.set_tool("select")
        self._draw_annotations()

    def _update_annotations_check(self):
        n = sum(a.visible for a in self.img.entry.annotations) if self.img else 0
        self.annotations_check.configure(text=f"Show annotations ({n or 'none'})")
        changed = any(o.unsaved for o in self._all_objects())
        for item in set(getattr(self, "unsaved_marks", {})) | set(self.projects):
            self._place_unsaved_mark(item)          # (gone after Save or Undo; a project's own, after its rename)
        self.save_button.configure(state=tk.NORMAL if changed else tk.DISABLED)
        self.unsaved_label.configure(text="unsaved changes" if changed else "")

    def m_to_screen(self, p):
        """Metres from the image centre (annotation coordinates) -> canvas pixels."""
        e = self.img.entry
        cw, ch = self.canvas_size()
        px, py = e.width / 2 + p[0] * 1e6 / e.pixel_um, e.height / 2 + p[1] * 1e6 / e.pixel_um
        return self._fx((px - self.cx) * self.scale + cw / 2), (py - self.cy) * self.scale + ch / 2

    def screen_to_m(self, sx, sy):
        e = self.img.entry
        ix, iy = self.screen_to_image(sx, sy)
        return (ix - e.width / 2) * e.pixel_um * 1e-6, (iy - e.height / 2) * e.pixel_um * 1e-6

    def _handles(self, a):
        """The selected annotation's handles: [(name, screen point)]. Rectangle/ellipse: 4 corners and "rotate" (above the top edge);
        line: its 2 ends; others: none (they can be moved)."""
        if a.kind in ("rectangle", "ellipse"):
            ha, hb = a.half_size()
            out = [((sx, sy), self.m_to_screen(a.to_world(sx * ha, sy * hb))) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
            top = self.m_to_screen(a.to_world(0, -hb))
            centre = self.m_to_screen(a.center)
            d = np.hypot(top[0] - centre[0], top[1] - centre[1]) or 1
            out.append(("rotate", (top[0] + (top[0] - centre[0]) / d * 25, top[1] + (top[1] - centre[1]) / d * 25)))
            return out
        if a.kind in TWO_POINT and len(a.verts) == 2:
            return [(("end", i), self.m_to_screen(p)) for i, p in enumerate(a.points())]
        return []

    def _draw_annotations(self):
        self._each_canvas(self._draw_annotations_here)
        self._show_color()                  # (a text selected, or not: the color button shows text's color, or shapes')

    def _draw_annotations_here(self):
        """The image's annotations on top of it, in image coordinates (they move and zoom with it); text and line widths stay the same
        size on screen. The selected one gets its handles."""
        c = self.canvas
        c.delete("ann")
        if not self.img or not self.show_annotations.get():
            return
        for k, a in enumerate(self.img.entry.annotations):
            if not a.visible:
                continue
            tags = ("ann", f"a{k}")
            pts = [self.m_to_screen(p) for p in a.points()]
            col = a.color
            if a.kind == "text":
                t = c.create_text(*pts[0], text=a.text, fill=col, font=a.font, anchor=tk.NE if self._flipped() else tk.NW,
                                  tags=tags)
                if a.attrs.get("showBorder") == "True":
                    c.create_rectangle(*c.bbox(t), outline=col, tags=tags)
            elif a.kind in ("rectangle", "ellipse") and len(pts) >= 3:
                c.create_polygon(*[v for p in pts for v in p], outline=col, fill="", width=2, tags=tags)
            elif len(pts) >= 2:
                c.create_line(*[v for p in pts for v in p], fill=col, width=2, tags=tags,
                              **({"arrow": tk.LAST, "arrowshape": (14, 16, 5)} if a.kind == "arrow" else {}))
                if a.kind == "measure" and len(pts) == 2:
                    (x0, y0), (x1, y1) = pts
                    length = max(np.hypot(x1 - x0, y1 - y0), 1e-9)
                    nx, ny = -(y1 - y0) / length * 6, (x1 - x0) / length * 6       # (6 px either side, across the line)
                    if a.attrs.get("showLineEnds") == "True":
                        for x, y in pts:
                            c.create_line(x - nx, y - ny, x + nx, y + ny, fill=col, width=2, tags=tags)
                    label = a.length_label()
                    if label:
                        c.create_text((x0 + x1) / 2 + 10, (y0 + y1) / 2, text=label, fill=col, font=a.font, anchor=tk.W, tags=tags)
        a = self.selected
        if a is not None and a in self.img.entry.annotations:
            k = self.img.entry.annotations.index(a)
            box = c.bbox(f"a{k}")
            if box and a.kind in ("text", "other"):
                c.create_rectangle(box[0] - 3, box[1] - 3, box[2] + 3, box[3] + 3, outline="white", dash=(3, 3), tags="ann")
            for name, (x, y) in self._handles(a):
                if name == "rotate":
                    top = self.m_to_screen(a.to_world(0, -a.half_size()[1]))
                    c.create_line(*top, x, y, fill="white", tags="ann")
                    c.create_oval(x - 5, y - 5, x + 5, y + 5, fill="white", outline="black", tags="ann")
                else:
                    c.create_rectangle(x - 4, y - 4, x + 4, y + 4, fill="white", outline="black", tags="ann")
        h = getattr(self, "hover_ann", None)
        if h is not None and h is not a and h.visible and h in self.img.entry.annotations:
            b = Annotation(h.attrs, h.verts)
            self._normalize_rotation(b)              # (where its handles will be once it is grabbed: see _handle_at)
            for name, (x, y) in self._handles(b):
                if name != "rotate":                 # (hollow: not selected yet)
                    c.create_rectangle(x - 4, y - 4, x + 4, y + 4, fill="", outline="white", width=2, tags="ann")

    def _handle_at(self, sx, sy):
        """What a press at canvas point (sx, sy) grabs, whatever the tool, as (name, annotation), or None: a handle within 7 px (the
        selected annotation's, rotate included, else a corner or line end of another visible one, the one drawn last, as its handles
        would be once it is selected); else ("move", annotation) for an outline within 6 px (a text: its box), the selected one's first."""
        grab = self._handle_point_at(sx, sy)
        if grab is not None or not self.img or not self.show_annotations.get():
            return grab
        anns = self.img.entry.annotations
        for a in ([self.selected] if self.selected in anns else []) + anns[::-1]:
            if a.visible and self._on_outline(a, sx, sy):
                return "move", a
        return None

    def _on_outline(self, a, sx, sy):
        """Whether canvas point (sx, sy) is within 6 px of annotation a's outline (a rectangle's edges, an ellipse's circumference, a
        line), or on a text's box."""
        if a.kind == "text":
            box = self.canvas.bbox(f"a{self.img.entry.annotations.index(a)}")
            return bool(box) and box[0] - 3 <= sx <= box[2] + 3 and box[1] - 3 <= sy <= box[3] + 3
        tol = 6 / self.scale * self.img.entry.pixel_um * 1e-6          # (6 screen px, in metres)
        p = self.screen_to_m(sx, sy)
        if a.kind in ("rectangle", "ellipse"):
            ha, hb = a.half_size()
            x, y = a.to_local(*p)
            if a.kind == "rectangle":
                ox, oy = abs(x) - ha, abs(y) - hb
                d = np.hypot(max(ox, 0), max(oy, 0)) if ox > 0 or oy > 0 else -max(ox, oy)
            else:                                    # (|g| / |grad g| for g = x²/a² + y²/b² - 1: close to the true distance near it)
                g = (x / (ha or 1e-12)) ** 2 + (y / (hb or 1e-12)) ** 2 - 1
                d = abs(g) / (2 * np.hypot(x / (ha or 1e-12) ** 2, y / (hb or 1e-12) ** 2) or 1e-30)
            return d <= tol
        pts = a.points()
        for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
            dx, dy = x1 - x0, y1 - y0
            t = min(max(((p[0] - x0) * dx + (p[1] - y0) * dy) / (dx * dx + dy * dy or 1), 0), 1)
            if np.hypot(p[0] - x0 - t * dx, p[1] - y0 - t * dy) <= tol:
                return True
        return False

    def _handle_point_at(self, sx, sy):
        """The handle within 7 px of canvas point (sx, sy), as (name, annotation), or None (see _handle_at)."""
        if not self.img or not self.show_annotations.get():
            return None
        anns = self.img.entry.annotations
        a = self.selected
        if a is not None and a in anns:
            for name, (x, y) in self._handles(a):
                if np.hypot(sx - x, sy - y) <= 7:
                    return name, a
        for a in reversed(anns):
            if a is self.selected or not a.visible:
                continue
            b = Annotation(a.attrs, a.verts)
            self._normalize_rotation(b)              # (a copy: selecting it normalizes it, which renames its corners)
            for name, (x, y) in self._handles(b):
                if name != "rotate" and np.hypot(sx - x, sy - y) <= 7:
                    return name, a
        return None

    def _hover_target(self, sx, sy):
        """The unselected shape or line whose handles are shown while the mouse is over it (so that they can be found and grabbed, with
        any tool), or None: a shape whose box (an ellipse's too, as its corners are its handles) is within 10 px of (sx, sy), the
        smallest of several; else a line within 10 px."""
        tol = 10 / self.scale * self.img.entry.pixel_um * 1e-6
        p = self.screen_to_m(sx, sy)
        best = None
        for k, a in enumerate(self.img.entry.annotations):
            if a is self.selected or not a.visible or not self._handles(a):
                continue
            if a.kind in ("rectangle", "ellipse"):
                ha, hb = a.half_size()
                x, y = a.to_local(*p)
                if abs(x) <= ha + tol and abs(y) <= hb + tol:
                    best = min(best or (np.inf, 0), (ha * hb, -k))
                continue
            (x0, y0), (x1, y1) = a.points()
            dx, dy = x1 - x0, y1 - y0
            t = min(max(((p[0] - x0) * dx + (p[1] - y0) * dy) / (dx * dx + dy * dy or 1), 0), 1)
            if np.hypot(p[0] - x0 - t * dx, p[1] - y0 - t * dy) <= tol:
                best = min(best or (np.inf, 0), (0.0, -k))
        return None if best is None else self.img.entry.annotations[-best[1]]

    def _hit(self, sx, sy):
        """The visible annotation at canvas point (sx, sy), or None: the smallest one containing it (shapes count by their inside),
        a line within 6 px, a text by its box; of several, the one drawn last (on top)."""
        tol = 6 / self.scale * self.img.entry.pixel_um * 1e-6          # (6 screen px, in metres)
        p = self.screen_to_m(sx, sy)
        best = None
        for k, a in enumerate(self.img.entry.annotations):
            if not a.visible:
                continue
            if a.kind == "text":
                box = self.canvas.bbox(f"a{k}")
                if box and box[0] - 3 <= sx <= box[2] + 3 and box[1] - 3 <= sy <= box[3] + 3:
                    best = min(best or (np.inf, None), (0.0, -k))
                continue
            if a.kind in ("rectangle", "ellipse"):
                ha, hb = a.half_size()
                x, y = a.to_local(*p)
                inside = ((abs(x) <= ha + tol and abs(y) <= hb + tol) if a.kind == "rectangle"
                          else (x / (ha + tol)) ** 2 + (y / (hb + tol)) ** 2 <= 1)
                if inside:
                    best = min(best or (np.inf, None), (ha * hb, -k))
                continue
            pts = a.points()
            for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
                dx, dy = x1 - x0, y1 - y0
                t = min(max(((p[0] - x0) * dx + (p[1] - y0) * dy) / (dx * dx + dy * dy or 1), 0), 1)
                if np.hypot(p[0] - x0 - t * dx, p[1] - y0 - t * dy) <= tol:
                    best = min(best or (np.inf, None), (0.0, -k))
                    break
        return None if best is None else self.img.entry.annotations[-best[1]]

    def _full_res_allowed(self):
        """Zoomed in past the pyramid: the full-resolution files are needed. If some are large and online-only (see needs_asking):
          - being downloaded already (or waiting to be): it goes to the front of the queue, and the sharp view is drawn when it is in;
          - else ask, once per image while the viewer runs (one question at a time: redraws while it is open don't ask again).
            Download (this section, or all) queues the download; Cancel is remembered: the view stays at the pyramid's sharpness
            (right-click on the section still downloads it)."""
        if self.full_res_ok is True:
            return True
        e = self.img.entry
        ask = [f for f in {e.files.get((self.z, c)) for c in range(len(e.channels))} if f and needs_asking(f)]
        if not ask:
            return True
        item = next((i for i, x in self.entries.items() if x is e), None)

        def not_now(note):
            self.detail_pending = False
            self._update_cursor()
            self.status.config(text=note)
            return False
        if item in self.download_progress:
            if item in self.download_queue and self.download_queue[0] != item:
                self.download_queue.remove(item)
                self.download_queue.insert(0, item)             # (the one being looked at comes next)
            return not_now(f"{e.name}: downloading its full resolution from OneDrive; the sharp view follows when it is in")
        if getattr(e, "full_res_declined", False):
            return not_now(f"{e.name}: full resolution not downloaded (online-only; right-click it to download)")
        if self.asking_download:
            return False
        self.asking_download = True
        try:
            size = size_text(sum(os.path.getsize(f) for f in ask))
            cloud = [i for i in self._entries_of(self._project_of(item)) if i in self.item_cloud and i not in self.download_progress]
            choice = ask_download(self.root, "Download from OneDrive?",
                                  f"Zooming further needs high-resolution files ({size} for this section) that are online-only.",
                                  all_label=f"Download for all sections ({size_text(self._cloud_bytes(cloud))})")
        finally:
            self.asking_download = False
        if choice is None:
            e.full_res_declined = True
            return not_now(f"{e.name}: full resolution not downloaded (online-only; right-click it to download)")
        self.download_full_res([item] if choice == "this" else [item] + [i for i in cloud if i != item])     # (this one first)
        return not_now(f"{e.name}: downloading its full resolution from OneDrive; the sharp view follows when it is in")

    def _update_cursor(self):
        """The pointer over the image: an hourglass while an image loads or a sharper view of it is read, a hand over an annotation's
        handle and the move arrows over its outline (they can be grabbed, whatever the tool), else the tool's."""
        busy = self.loading or getattr(self, "detail_pending", False)
        for c in self.canvases:
            over = getattr(self, "over_handle", None)
            c.configure(cursor="watch" if busy else "hand2" if over == "handle" else "fleur" if over == "move" or self.tool == "select"
                        else "crosshair")
        self._hires_note_timer()

    def _hires_note_timer(self):
        """While a sharper view is read (after the image itself has loaded): "Loading high-resolution images..." over the image once it
        has taken HIRES_NOTE_DELAY (a short wait needs no words: the hourglass says enough); gone when it is in."""
        waiting = getattr(self, "detail_pending", False) and not self.loading
        job = getattr(self, "hires_job", None)
        if waiting and job is None and not getattr(self, "hires_note", False):
            self.hires_job = self.root.after(int(HIRES_NOTE_DELAY * 1000), self._show_hires_note)
        elif not waiting:
            if job is not None:
                self.root.after_cancel(job)
                self.hires_job = None
            if getattr(self, "hires_note", False):
                self.hires_note = False
                self._each_canvas(lambda: self.canvas.delete("hiresnote"))

    def _show_hires_note(self):
        self.hires_job = None
        if getattr(self, "detail_pending", False) and not self.loading and self.img is not None:
            self.hires_note = True
            self._each_canvas(self._draw_hires_note_here)

    def _draw_hires_note_here(self):
        c = self.canvas
        c.delete("hiresnote")
        if not getattr(self, "hires_note", False):
            return
        cw, ch = self.canvas_size()
        t = c.create_text(cw / 2, ch / 2, text="Loading high-resolution images...", fill="#e0e0e0", font=("Segoe UI", 16),
                          tags="hiresnote")
        x0, y0, x1, y1 = c.bbox(t)
        c.tag_lower(c.create_rectangle(x0 - 14, y0 - 8, x1 + 14, y1 + 8, fill="#202020", outline="#606060", tags="hiresnote"), t)

    def set_tool(self, tool):
        """Choose a tool. It stays chosen for one annotation after another (Select, or Escape, to go back to selecting)."""
        self.tool = tool
        for name, b in self.tool_buttons.items():
            on = name == tool
            style_toggle(b, on)
        self._update_cursor()
        self._show_color()
        if tool != "select" and not self.show_annotations.get():
            self.show_annotations.set(True)
            self.on_annotations_toggle()

    def _remember(self):
        """Keep the image's annotations as they are now, for Undo (which goes back through annotation and colour-map changes alike,
        the newest first: see undo_annotation)."""
        e = self.img.entry
        e.undo = getattr(e, "undo", [])[-49:] + [[Annotation(a.attrs, a.verts) for a in e.annotations]]
        self.history = self.history[-199:] + [("ann", e)]

    def _forget_remember(self, e):
        """A change that _remember() prepared for didn't happen after all."""
        e.undo.pop()
        if self.history and self.history[-1] == ("ann", e):
            self.history.pop()

    def _changed(self):
        e = self.img.entry
        if not e.annotations_changed:
            e.annotations_changed = True
            self._label_item(next((i for i, x in self.entries.items() if x is e), None))
        self._update_annotations_check()
        self._draw_annotations()

    def _normalize_rotation(self, a):
        """LAS X stores upright shapes with a rotation of 90/180/270 and their sizes swapped. Before such a shape is resized, its rotation
        is folded into its vertices, so that its handles sit where they look like they should."""
        r = a.rotation % 360
        if a.kind in ("rectangle", "ellipse") and r and abs(r - round(r / 90) * 90) < 1e-6:
            world = [(x, y) for x, y in (a.to_world(*v) for v in a.local())]
            cx, cy = a.center
            a.set_pose(a.center, 0)
            a.set_local([(x - cx, y - cy) for x, y in world])
            ha, hb = a.half_size()
            a.set_half_size(ha, hb)

    def ann_press(self, event):
        """A press on the image, when annotations are shown: grabs a handle or outline under the mouse and starts moving it (with any
        tool: see _handle_at), else
        starts drawing (with a drawing tool), or selects an annotation and starts moving it (Select). True if it was taken (else the
        press pans the view)."""
        if not self.img or not self.img.entry.pixel_um or not self.show_annotations.get():
            return False
        p = self.screen_to_m(event.x, event.y)
        grab = self._handle_at(event.x, event.y)
        if grab is not None:                         # (whatever the tool: a handle is grabbed, not drawn from)
            name, a = grab
            if a is not self.selected:
                self._normalize_rotation(a)
                self.selected = a
            self._remember()
            self.ann_drag = (name, p, Annotation(a.attrs, a.verts), a)
            self._draw_annotations()
            return True
        if self.tool == "text":
            answer = ask_text(self.root, "New text", "", color=self.text_color)
            text, color = answer or ("", None)
            if text:
                self.text_color = color                 # (the next text's color too)
                self._remember()
                a = Annotation.new("text", color, p, [], text=text)
                self.img.entry.annotations.append(a)
                self.selected = a
                self._changed()
            return True
        if self.tool in ("rectangle", "ellipse") + TWO_POINT:
            self._remember()
            verts = [(0.0, 0.0), (0.0, 0.0)] if self.tool in TWO_POINT else []
            a = Annotation.new(self.tool, self.ann_color, p, verts)
            if self.tool not in TWO_POINT:
                a.set_half_size(0, 0)
            self.img.entry.annotations.append(a)
            self.selected = a
            self.ann_drag = ("draw", p, a)
            return True
        hit = self._hit(event.x, event.y)
        if hit is not None and hit is not self.selected:
            self._normalize_rotation(hit)        # (same shape on screen; its handles then sit on its own corners and top)
        self.selected = hit
        self._draw_annotations()
        if hit is None:
            return False
        self._remember()
        self.ann_drag = ("move", p, Annotation(hit.attrs, hit.verts), hit)
        return True

    def ann_motion(self, event):
        mode, p0 = self.ann_drag[0], self.ann_drag[1]
        p = self.screen_to_m(event.x, event.y)
        if mode == "draw":
            a = self.ann_drag[2]
            if a.kind in TWO_POINT:
                a.set_pose(((p0[0] + p[0]) / 2, (p0[1] + p[1]) / 2), 0)
                a.set_local([((p0[0] - p[0]) / 2, (p0[1] - p[1]) / 2), ((p[0] - p0[0]) / 2, (p[1] - p0[1]) / 2)])
            else:
                a.set_pose(((p0[0] + p[0]) / 2, (p0[1] + p[1]) / 2), 0)
                a.set_half_size(abs(p[0] - p0[0]) / 2, abs(p[1] - p0[1]) / 2)
        else:
            orig, a = self.ann_drag[2], self.ann_drag[3]
            if mode == "move":
                a.set_pose((orig.center[0] + p[0] - p0[0], orig.center[1] + p[1] - p0[1]))
            elif mode == "rotate":                  # the handle points the shape's top edge at the mouse
                dx, dy = p[0] - orig.center[0], p[1] - orig.center[1]
                a.set_pose(orig.center, np.degrees(np.arctan2(dx, -dy)))
            elif mode[0] == "end":                  # a line's end: the other end stays
                i = mode[1]
                ends = orig.points()
                ends[i] = p
                (x0, y0), (x1, y1) = ends
                a.set_pose(((x0 + x1) / 2, (y0 + y1) / 2), 0)
                a.set_local([((x0 - x1) / 2, (y0 - y1) / 2), ((x1 - x0) / 2, (y1 - y0) / 2)])
            else:                                   # a corner of a rectangle/ellipse: the opposite corner stays
                sx, sy = mode
                ha, hb = orig.half_size()
                ox, oy = -sx * ha, -sy * hb
                lx, ly = orig.to_local(*p)
                a.set_pose(orig.to_world((lx + ox) / 2, (ly + oy) / 2))
                a.set_half_size(abs(lx - ox) / 2, abs(ly - oy) / 2)
        self._draw_annotations()

    def ann_release(self, event):
        mode = self.ann_drag[0]
        a = self.ann_drag[2] if mode == "draw" else self.ann_drag[3]
        self.ann_drag = None
        e = self.img.entry
        if mode == "draw":
            tiny = (np.hypot(*np.subtract(*a.points())) if a.kind in TWO_POINT else 2 * min(a.half_size())) * 1e6 / e.pixel_um * self.scale < 4
            if tiny:                                # a click without a drag: no shape
                e.annotations.remove(a)
                self._forget_remember(e)
                self.selected = None
                self._draw_annotations()
            else:
                self._changed()
            return
        before = e.undo[-1][e.annotations.index(a)] if a in e.annotations else None
        if before is not None and before.attrs == a.attrs and before.verts == a.verts:
            self._forget_remember(e)                # (nothing moved)
        else:
            self._changed()

    def ann_double(self, event):
        """Double-click: edit a text's words; elsewhere, fit the image to the window."""
        a = self._hit(event.x, event.y) if self.img and self.show_annotations.get() else None
        if a is not None and a.kind == "text":
            answer = ask_text(self.root, "Edit text", a.text, color=a.color)
            if answer is not None and answer != (a.text, a.color):
                if answer[1] != a.color:
                    self.text_color = answer[1]         # (the color picked last is the next text's)
                self._remember()
                a.text, a.color = answer
                self._changed()
            return
        self.fit()

    def delete_annotation(self, event=None):
        if self.img and self.selected in self.img.entry.annotations:
            self._remember()
            self.img.entry.annotations.remove(self.selected)
            self.selected = None
            self._changed()

    def undo_annotation(self, event=None):
        """Undo (Ctrl+Z): the last change, whichever image it was on: annotations, contrast, or colour maps (one section's, or all
        sections' at once). What is undone is still unsaved (Save writes the result)."""
        while self.history:
            kind, data = self.history.pop()
            if kind == "ann":
                e = data
                if not getattr(e, "undo", None):
                    continue
                e.annotations = e.undo.pop()
                e.annotations_changed = True
                if self.img is not None and self.img.entry is e:
                    self.selected = None
                    self._changed()
                else:
                    self._label_item(self._item_of(e))
                    self._update_annotations_check()
                self.status.config(text=f"Undone: an annotation change on {e.name}")
                return
            if kind == "move":
                item, parent, index = data
                if not (self.tree.exists(item) and self.tree.exists(parent)) or item in self.deleted:
                    continue
                crossed = self._project_of(item) != self._project_of(parent)
                self._place(item, parent, index)
                self._structure_changed(members=crossed)
                self.tree.see(item)
                self.status.config(text=f"Undone: moving {(self.entries.get(item) or self.folders.get(item)).name}")
                return
            if kind == "delete":
                item, e, parent, index = data
                if item not in self.deleted or not self.tree.exists(parent):
                    continue
                del self.deleted[item]
                self.entries[item] = e
                self._place(item, parent, index)
                self._structure_changed(members=True)
                self.tree.selection_set(item)
                self.tree.see(item)
                self.status.config(text=f"Undone: deleting {e.name}")
                return
            if kind == "range":
                for e, i, lo, hi in reversed(data):
                    e.channels[i].black, e.channels[i].white = lo, hi
                    if self.img is not None and self.img.entry is e and i < len(self.ch_widgets):
                        self.ch_widgets[i][1].set(lo, hi)
                        self.ch_settings[i]["lo"], self.ch_settings[i]["hi"] = lo, hi
                self.contrast_step = None
                self.schedule_redraw()
                self._contrast_settled({e for e, _, _, _ in data})
                self.status.config(text=f"Undone: contrast change on {data[0][0].name}")
                return
            if kind == "lut":
                for e, i, old in data:
                    e.channels[i].lut = old
                self._luts_changed({e for e, _, _ in data})
                names = sorted({e.name for e, _, _ in data})
                self.status.config(text=f"Undone: color map change on {names[0]}" + (f" and {len(names) - 1} more" if len(names) > 1 else ""))
                return
        self.status.config(text="Nothing to undo")

    def _item_of(self, entry):
        return next((i for i, x in self.entries.items() if x is entry), None)

    def _luts_changed(self, changed):
        """After colour maps of these images changed (picked, or undone): the image on screen and its channel controls, the list's
        "unsaved" marks, Save, and their thumbnails."""
        if self.img is not None and self.img.entry in changed:
            for i in range(len(self.ch_widgets)):
                on, slider, swatch = self.ch_widgets[i]
                swatch.color = self._swatch_color(i)
                slider.set_color(swatch.color)
                self._show_channel_state(i)
            self.redraw()
        items = [self._item_of(e) for e in changed]
        for item in items:
            self._label_item(item)
        self._update_annotations_check()
        self._refresh_thumbnails(items)

    def _refresh_thumbnails(self, items):
        """Make these images' thumbnails again (their colour maps or contrast changed), in the background; the slide's picture above
        the list too, if it is one of them."""
        ov = getattr(self, "slide_entry", None)
        if ov is not None and self.slide.winfo_manager() and any(self.entries.get(i) is ov for i in items):
            self.slide_gen += 1                 # (a picture still being made with the old contrast is dropped)
            self._make_slide_picture()
        if self.thumbs_started:
            gen = self.thumb_gen

            def work():
                for item in items:
                    if item is not None and gen == self.thumb_gen:
                        try:
                            image, state = make_thumbnail(self.entries[item])
                        except Exception:
                            image, state = None, "none"
                        self.results.put(("thumb", gen, (item, image, state)))
            threading.Thread(target=work, daemon=True).start()

    def _color_is_text(self):
        """The color button stands for text's color (with Text chosen, or a text selected), else for shapes' and lines'."""
        a = getattr(self, "selected", None)
        return self.tool == "text" or (self.tool == "select" and a is not None and a.kind == "text")

    def _show_color(self):
        color = self.text_color if self._color_is_text() else self.ann_color
        if self.color_button.cget("bg") != color:
            self.color_button.configure(bg=color, activebackground=color)

    def choose_color(self):
        from tkinter import colorchooser
        text = self._color_is_text()
        _, color = colorchooser.askcolor(color=self.text_color if text else self.ann_color, parent=self.root,
                                         title="Text color" if text else "Annotation color")
        if not color:
            return
        if text:
            self.text_color = color
        else:
            self.ann_color = color
        self._show_color()
        a = self.selected
        if self.img and a in self.img.entry.annotations and (a.kind == "text") == text:
            self._remember()
            a.color = color
            self._changed()

    def save_annotations(self):
        """The Save button: write every change not saved yet into LAS X's files (only ever when it is clicked): edited annotations, color
        maps, contrast and flips (into each image's .xlif) and renames (the Name in the .xlif / .xlcf; with "Also rename the files and folders on disk" in Settings,
        the files and folders too, see plan_disk_rename), the list's order (dragged, deleted: see _save_structure). Each changed file is
        backed up first (BACKUP_SUFFIX)."""
        items = [(i, o) for i, o in list(self.entries.items()) + list(self.folders.items())
                 + [(pid, p.get("root")) for pid, p in self.projects.items()] if o is not None and o.unsaved]
        if not items and not self.deleted:
            return
        notes = [(i, o) for i, o in items if o.annotations_changed]
        maps = [(i, o) for i, o in items if isinstance(o, ImageEntry) and o.luts_changed]
        ranges = [(i, o) for i, o in items if isinstance(o, ImageEntry) and o.contrast_changed]
        flips = [(i, o) for i, o in items if isinstance(o, ImageEntry) and o.flip_changed]
        renames = [(i, o) for i, o in items if o.name != o.saved_name]
        on_disk = bool(renames) and self.settings.get("rename_files", False)
        away = [o for i, o in self.entries.items() if self.home.get(i) not in (None, self._project_of(i))]
        reordered = any(getattr(o, "order_changed", False) for _, o in items)
        structure = reordered or bool(self.deleted)
        if not self.save_warned or self.deleted or away:
            n_del, n_away = len(self.deleted), len(away)
            what = ([f"the annotations of {len(notes)} image{'s' * (len(notes) != 1)}"] if notes else []) + \
                   ([f"the color maps of {len(maps)} image{'s' * (len(maps) != 1)}"] if maps else []) + \
                   ([f"the contrast of {len(ranges)} image{'s' * (len(ranges) != 1)}"] if ranges else []) + \
                   ([f"the flip of {len(flips)} image{'s' * (len(flips) != 1)}"] if flips else []) + \
                   ([f"{len(renames)} rename{'s' * (len(renames) != 1)}" + (" (files and folders on disk too)" if on_disk else "")]
                    if renames else []) + (["the new order of the list"] if reordered else [])
            extra = ((f"\n\n{n_away} image{'s' * (n_away != 1)} dragged into another project move{'s' * (n_away == 1)} there, files "
                      "and all." if away else "")
                     + (f"\n\nDeleted: {', '.join(e.name for e in list(self.deleted.values())[:5])}"
                        f"{' and more' if n_del > 5 else ''}. Their files go to the Recycle Bin (on a network drive, which has "
                        f"none: into a \"Deleted by XLEF viewer\" folder in the experiment's folder)." if n_del else ""))
            if not messagebox.askokcancel(
                    "Save changes into LAS X's files",
                    f"This writes {' and '.join(what)} into LAS X's files. A copy of each file as it is now is kept next to it "
                    f"(ending in {BACKUP_SUFFIX}).{extra}\n\nClose this experiment in LAS X first: LAS X may otherwise overwrite "
                    f"these changes with its own copy.", parent=self.root):
                return
            self.save_warned = True
        failed = []
        for item, e in notes:
            try:
                write_annotations(e.xlif_path, e.annotations)
                e.annotations_changed = False
            except Exception as err:
                failed.append(f"annotations of {e.name}: {err}")
        for item, e in maps:
            try:
                write_luts(e.xlif_path, [c.lut for c in e.channels])
                for c in e.channels:
                    c.saved_lut = c.lut
            except Exception as err:
                failed.append(f"color maps of {e.name}: {err}")
        for item, e in ranges:
            try:
                write_contrast(e.xlif_path, [(c.black, c.white) for c in e.channels], [c.max_value for c in e.channels])
                for c in e.channels:
                    c.saved_black, c.saved_white = c.black, c.white
            except Exception as err:
                failed.append(f"contrast of {e.name}: {err}")
        for item, e in flips:
            try:
                write_flip(e.xlif_path, e.flipped)
                e.saved_flipped = e.flipped
                self._set_flip_note(e, None)
            except Exception as err:
                failed.append(f"flip of {e.name}: {err}")
        touched, binned, aside = set(), 0, 0
        if structure or (renames and on_disk):
            # (nothing may hold the files open while they are moved or renamed: the image on screen, kept images, background readers)
            self._close_image()
            self.forget_kept_images()
            self.thumb_gen += 1
            self.slide_gen += 1
            time.sleep(0.3)
        if structure:
            touched, binned, aside = self._save_structure(failed)
        if renames and on_disk:
            ordered = ([r for r in renames if isinstance(r[1], ImageEntry)] +          # (images first, then folders, deepest first)
                       sorted([r for r in renames if isinstance(r[1], FolderEntry)], key=lambda r: -self._tree_depth(r[0])))
            for item, o in ordered:
                try:
                    carry_out(plan_disk_rename(o, o.name))
                    o.saved_name = o.name
                except Exception as err:
                    failed.append(f"{o.saved_name} -> {o.name}: {err}")
        else:
            for item, o in renames:
                try:
                    if o.path is None:
                        raise ValueError("its folder file is unreadable")
                    write_element_name(o.path, o.name)
                    o.saved_name = o.name
                except Exception as err:
                    failed.append(f"{o.saved_name} -> {o.name}: {err}")
        for item, _ in items:
            self._label_item(item)
        self._update_annotations_check()
        done = max(len(items) + binned + aside - len(failed), 0)
        self.status.config(text=f"Saved {done} change{'s' * (done != 1)} into LAS X's files" + (f"; {len(failed)} failed" if failed else "")
                                + (f". {binned} deleted image{'s' * (binned != 1)} moved to the Recycle Bin" if binned else "")
                                + (f". {aside} deleted image{'s' * (aside != 1)} moved to the \"Deleted by XLEF viewer\" folder "
                                   "(this drive has no Recycle Bin)" if aside else ""))
        if failed:
            messagebox.showerror("Some changes not saved", "These could not be saved (their files were left as they were):\n\n"
                                 + "\n".join(failed), parent=self.root)
        reread = set(touched)
        if renames and on_disk:
            for _, o in renames:          # (anything that failed is dropped: those projects are read again from the files)
                o.name = o.saved_name
            for _, e in notes:
                e.annotations_changed = False
            reread |= {self._project_of(i) for i, _ in renames} - {None}
        if reread:
            self.history = [h for h in self.history if h[0] in ("ann", "lut", "range")]
        for pid in [p for p in self.projects if p in reread]:      # (paths or contents changed: read those experiments again)
            path = self.projects[pid]["path"]
            self.close_project(pid, ask=False)
            self.open_xlef(path)

    # --- the cloud after a name: what it means, and downloading full resolution

    def _cloud_tip_text(self, entry):
        return f"Full resolution image of {entry.name} is in cloud. Zoom in or right-click to download"

    def _tree_motion(self, event):
        """The pointer on the name of an image with a cloud, or on the cloud: a note saying what the cloud means. On a project's name:
        its full path, and whether it has unsaved changes."""
        item = self.tree.identify_row(event.y)
        if self._unsaved_tip(item, event):
            return
        if item in self.projects and self.tree.bbox(item):
            x, y, w, h = self.tree.bbox(item)
            start = x + 20                                       # (after the open mark)
            if not start - 8 <= event.x <= start + tkfont.nametofont("TkDefaultFont").measure(self.tree.item(item, "text")) + 8:
                self._hide_tree_tip()
                return
            changed = self._unsaved_names(item)
            note = (f"Unsaved changes: {', '.join(changed[:5])}{' and more' if len(changed) > 5 else ''}" if changed
                    else "No unsaved changes")
            if self.tree_tip is None:
                self.tree_tip = PointerTip(self.tree)
            self.tree_tip.show(f"{os.path.normpath(self.projects[item]['path'])}\n{note}", event.x_root, event.y_root)
            return
        if item not in self.item_cloud or not self.tree.bbox(item):
            self._hide_tree_tip()
            return
        x, y, w, h = self.tree.bbox(item)
        font = tkfont.nametofont("TkDefaultFont")
        start = x + 20 * (self._tree_depth(item) + 1)            # (where the name starts: the tree's indent per level, the open mark)
        cloud = start + font.measure(self.entries[item].name + "  ")
        if not start - 8 <= event.x <= cloud + font.measure("☁") + 8:       # (anywhere on the name or its cloud)
            self._hide_tree_tip()
            return
        if self.tree_tip is None:
            self.tree_tip = PointerTip(self.tree)
        self.tree_tip.show(self._cloud_tip_text(self.entries[item]), event.x_root, event.y_root)

    def _hide_tree_tip(self):
        if self.tree_tip is not None:
            self.tree_tip.hide()
            self.tree_tip = None

    def _cloud_bytes(self, items):
        return sum(os.path.getsize(f) for i in items for f in self.entries[i].online_only_files())

    def download_full_res(self, items):
        """Download these sections' full-resolution files from OneDrive in the background, one section after another (progress after
        the name in the list; the cloud goes when a section is done). OneDrive downloads a whole file as soon as any byte of it is
        read, so the first byte of each file is read."""
        items = [i for i in items if i in self.item_cloud and i not in self.download_queue and i not in self.download_progress]
        for i in items:
            self.download_queue.append(i)
            self.download_progress[i] = "waiting"
            self._label_item(i)
        self._start_download_thread()
        if items:
            self.status.config(text=f"Downloading full resolution of {len(items)} section{'s' * (len(items) != 1)} from OneDrive ...")

    def _start_download_thread(self):
        if self.download_queue and self.download_thread is None:
            gen, stop = self.download_gen, threading.Event()
            self.download_stop = stop
            self.download_thread = threading.Thread(target=self._download_work, args=(gen, stop), daemon=True)
            self.download_thread.start()

    def _download_work(self, gen, stop):
        try:
            while self.download_queue and not stop.is_set() and gen == self.download_gen:
                item = self.download_queue.pop(0)
                entry = self.entries.get(item)
                files = entry.online_only_files() if entry is not None else []
                sizes = [os.path.getsize(f) for f in files]
                done = 0
                total, last = max(sum(sizes), 1), -1
                for f, size in zip(files, sizes):
                    if stop.is_set() or gen != self.download_gen:
                        break
                    # The read (in a thread of its own) returns only when OneDrive has downloaded the whole file; meanwhile, how much
                    # of it is on this computer already is asked 4 times a second, for the percentage.
                    failed = []

                    def fetch(f=f):
                        try:
                            with open(f, "rb") as fh:
                                fh.read(1)
                        except OSError as e:
                            failed.append(e)
                    reader = threading.Thread(target=fetch, daemon=True)
                    reader.start()
                    while True:
                        here = min(stored_bytes(f) or 0, size)
                        pct = int(100 * (done + here) / total)
                        if pct != last:
                            last = pct
                            self.results.put(("download", gen, (item, pct)))
                        reader.join(0.25)
                        if not reader.is_alive():
                            break
                    if failed:
                        self.results.put(("download_error", gen, f"Couldn't download {os.path.basename(f)}: {failed[0]}"))
                        break
                    done += size
                self.results.put(("downloaded", gen, item))
        finally:
            self.results.put(("download_idle", gen, None))

    def stop_downloads(self):
        if self.download_stop is not None:
            self.download_stop.set()
        for i in self.download_queue:
            self.download_progress.pop(i, None)
            self._label_item(i)
        self.download_queue.clear()
        self.status.config(text="Downloads stopped (the file being downloaded finishes first)")

    def _downloaded(self, item):
        self.download_progress.pop(item, None)
        entry = self.entries.get(item)
        if entry is not None and not entry.online_only_files():
            self.item_cloud.discard(item)
            self.thumb_grid.rows = self._grid_rows_all()          # (its thumbnail's name: no longer blue, no cloud)
            self.thumb_grid.layout()
            if self.img is not None and self.img.entry is entry:
                self.schedule_redraw(0)          # (zoomed in, waiting for it: the sharp view can be read now)
        self._label_item(item)
        if not self.download_queue:
            stopped = self.download_stop is not None and self.download_stop.is_set()
            self.status.config(text="Downloads stopped" if stopped else "Full-resolution downloads finished")

    def _tree_menu(self, event):
        """Right-click in the list: the menu for that row, its top at the middle of the row (so it is plain which row it is for)."""
        item = self.tree.identify_row(event.y)
        if item:
            box = self.tree.bbox(item)
            y = self.tree.winfo_rooty() + box[1] + box[3] // 2 if box else event.y_root
            self._item_menu(item, event.x_root, y)

    def _project_menu(self, menu, pid):
        if pid in self.projects:
            menu.add_command(label=f"Close project {self.projects[pid]['name']}", command=lambda: self.close_project(pid))

    def _item_menu(self, item, x_root, y_root):
        """Right-click on a section (list or thumbnail): download its full resolution, or every section's."""
        self._hide_tree_tip()
        menu = tk.Menu(self.root, tearoff=0)
        if item in self.entries:
            if item in self.item_cloud and item not in self.download_progress:
                menu.add_command(label=f"Download full resolution this section ({size_text(self._cloud_bytes([item]))})",
                                 command=lambda: self.download_full_res([item]))
            else:
                menu.add_command(label="Download full resolution this section" + (" (downloading)" if item in self.download_progress
                                                                                  else " (on this computer)"), state=tk.DISABLED)
        pid = self._project_of(item)
        if item in self.projects:                   # (a project's name: closing it, or downloading its sections)
            self._project_menu(menu, pid)
            menu.add_separator()
        cloud = [i for i in self._entries_of(pid) if i in self.item_cloud and i not in self.download_progress]
        if cloud:
            menu.add_command(label=f"Download full resolution all sections ({len(cloud)} in the cloud, "
                                   f"{size_text(self._cloud_bytes(cloud))})", command=lambda: self.download_full_res(cloud))
        else:
            queued = any(i in self.download_progress for i in self.entries)
            menu.add_command(label="Download full resolution all sections (" + ("all already queued" if queued else "none in the cloud")
                                   + ")", state=tk.DISABLED)
        if self.download_queue or self.download_thread is not None:
            menu.add_separator()
            menu.add_command(label="Stop downloading", command=self.stop_downloads)
        if item not in self.projects and pid is not None:
            menu.add_separator()
            self._project_menu(menu, pid)
        menu.add_separator()
        if item in self.entries:
            menu.add_command(label="Delete section", command=lambda: self.delete_section(item))
        target = self._disk_location(item)
        if target is not None:
            menu.add_command(label=REVEAL_LABEL, command=lambda: reveal_in_file_manager(*target))
        menu.add_command(label="Properties", command=lambda: self.show_properties(item))
        try:
            menu.tk_popup(int(x_root), int(y_root))
        finally:
            menu.grab_release()

    def _disk_location(self, item):
        """Where a list row is on disk, for "Open folder in Windows Explorer": (path, select it). An image: its first TIFF (in the
        folder with its channels); a folder: the folder itself (its .xlcf is in its Metadata folder); a project: its .xlef."""
        if item in self.entries:
            e = self.entries[item]
            files = e.all_files()
            path = files[0] if files else e.xlif_path
            return (path, True) if path and os.path.exists(path) else None
        if item in self.projects:
            path = self.projects[item]["path"]
            return (path, True) if os.path.exists(path) else None
        folder = self.folders.get(item)
        path = getattr(folder, "path", None)
        if not path:
            return None
        here = os.path.dirname(path)
        if os.path.basename(here).lower() == "metadata":
            here = os.path.dirname(here)
        return (here, False) if os.path.isdir(here) else None

    def show_properties(self, item):
        """Right-click > Properties: what LAS X's Properties window shows for an image (see image_properties), or for a folder or a
        project its file and how many images it holds. Read-only text that can be selected and copied."""
        obj = self.entries.get(item) or self.folders.get(item)
        pid = self._project_of(item)
        experiment = self.projects[pid]["path"] if pid else None
        try:
            if isinstance(obj, ImageEntry) and not obj.problem:
                props = image_properties(obj, experiment)
            else:
                name = self.projects[item]["name"] if item in self.projects else obj.name if obj is not None else "?"
                count = []

                def walk(node):
                    for child in self.tree.get_children(node):
                        count.append(child in self.entries)
                        walk(child)
                walk(item)
                path = experiment if item in self.projects else getattr(obj, "path", None) or getattr(obj, "xlif_path", None)
                props = {"header": [("Name", name), ("File", path or "(none)"), ("Images", str(sum(count))),
                                    ("Folders", str(len(count) - sum(count)))] + ([("Problem", obj.problem)] if getattr(obj, "problem", None) else []),
                         "sections": [], "stamps": []}
        except Exception as e:
            messagebox.showerror("Properties", f"The metadata could not be read:\n{e}", parent=self.root)
            return
        win = tk.Toplevel(self.root)
        win.title(f"Properties - {props['header'][0][1]}")
        win.geometry(f"640x760+{self.root.winfo_rootx() + 120}+{self.root.winfo_rooty() + 60}")
        frame = ttk.Frame(win)
        frame.pack(fill=tk.BOTH, expand=True)
        text = tk.Text(frame, wrap=tk.WORD, font=("Segoe UI", 9), padx=12, pady=10, relief=tk.FLAT, borderwidth=0)
        bar = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=text.yview)
        text.configure(yscrollcommand=bar.set)
        bar.pack(side=tk.RIGHT, fill=tk.Y)
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        text.tag_configure("label", font=("Segoe UI", 9))
        text.tag_configure("value", font=("Segoe UI", 9, "bold"))
        text.tag_configure("title", font=("Segoe UI", 9, "bold underline"), spacing1=10)
        text.tag_configure("rule", font=("Segoe UI", 2), background=PAL["muted"])
        text.tag_configure("table", font=("Consolas", 9))
        text.configure(tabs=("150p",))
        for label, value in props["header"]:
            text.insert(tk.END, f"{label}:\t", "label")
            text.insert(tk.END, f"{value}\n", "value")
        if props["sections"] or props["stamps"]:
            text.insert(tk.END, " \n", "rule")
        for title, lines in props["sections"]:
            text.insert(tk.END, f"{title}\n", "title")
            text.insert(tk.END, "\n".join(lines) + "\n")
        if props["stamps"]:
            text.insert(tk.END, "Time Stamps\n", "title")
            rows = [("Frame", "Relative Time (s)", "Absolute Time (h:m:s.ms)", "Date")] + props["stamps"]
            for r in rows:
                text.insert(tk.END, f"{str(r[0]):<7}{r[1]:<19}{r[2]:<26}{r[3]}\n", "table")
        text.configure(state=tk.DISABLED)
        ttk.Button(win, text="Close", command=win.destroy).pack(anchor=tk.E, padx=10, pady=8)
        win.bind("<Escape>", lambda e: win.destroy())
        return text

    def _tree_depth(self, item):
        depth = 0
        while self.tree.parent(item):
            item = self.tree.parent(item)
            depth += 1
        return depth

    # --- renaming in the list: click a name that is already selected (as in Windows Explorer)

    def _cancel_rename_job(self):
        if self.rename_job is not None:
            self.root.after_cancel(self.rename_job)
            self.rename_job = None

    def _tree_press(self, event):
        """A click in the list: on the name of an image or folder that is already selected, start renaming it (after a moment, so that
        a double-click, which opens or closes a folder, doesn't)."""
        self._cancel_rename_job()
        item = self.tree.identify_row(event.y)
        self.drag_item = (item, event.y, False) if item in self.entries or item in self.folders else None
        if (item and item in self.tree.selection() and (item in self.entries or item in self.folders)
                and "text" in self.tree.identify_element(event.x, event.y)):
            self.rename_job = self.root.after(600, lambda: self.start_rename(item))

    def start_rename(self, item):
        """An edit box over the item's name: Enter or clicking elsewhere keeps the new name (shown at once, written by Save), Escape
        keeps the old one."""
        self.rename_job = None
        obj = self.entries.get(item) or self.folders.get(item)
        box = self.tree.bbox(item)
        if obj is None or not box:
            return
        if self.rename_box is not None:
            self.rename_box.destroy()
        x, y, w, h = box
        indent = 20 * (self._tree_depth(item) + 1)          # (the tree's indent per level, and the open/close mark)
        e = self.rename_box = ttk.Entry(self.tree)
        e.insert(0, obj.name)
        e.select_range(0, tk.END)
        e.place(x=x + indent - 4, y=y, width=max(w - indent + 4, 80), height=h)
        e.focus_set()
        done = []

        def finish(keep):
            if done:
                return
            done.append(1)
            new = e.get().strip()
            e.destroy()
            self.rename_box = None
            self.tree.focus_set()
            if keep and new and new != obj.name:
                self.rename(item, new)
        e.bind("<Return>", lambda ev: finish(True))
        e.bind("<KP_Enter>", lambda ev: finish(True))
        e.bind("<Escape>", lambda ev: (finish(False), "break")[1])
        e.bind("<FocusOut>", lambda ev: finish(True))

    def rename(self, item, new):
        """Give an image or folder a new name: shown at once (marked "unsaved"), written into LAS X's files by Save."""
        obj = self.entries.get(item) or self.folders.get(item)
        bad = [c for c in BAD_NAME_CHARS if c in new]
        if bad:
            messagebox.showerror("Can't rename", f"A name can't contain {' '.join(bad)}.", parent=self.root)
            return
        if isinstance(obj, FolderEntry) and obj.path is None:
            messagebox.showerror("Can't rename", "This folder's own file is unreadable, so it can't be renamed here.", parent=self.root)
            return
        obj.name = new
        self._label_item(item)
        self.thumb_grid.rows = self._grid_rows_all()       # (the thumbnails stay; only their names and headings change)
        self.thumb_grid.layout()
        self._update_annotations_check()
        if self.img is not None and self.img.entry is obj:
            self._update_status()
        self.status.config(text=f"Renamed to {new}: not saved yet (Save writes it into LAS X's files)")

    # --- the list's order: dragging names (or thumbnails), deleting images. Shown at once, written by Save (see _save_structure).

    def _all_objects(self):
        """Everything that can have unsaved changes: the images, the folders, and the projects (their order of contents)."""
        return (list(self.entries.values()) + list(self.folders.values())
                + [p["root"] for p in self.projects.values() if p.get("root") is not None])

    def _container(self, item):
        """The FolderEntry a list item that holds others (a project or a folder) stands for, or None."""
        return self.projects[item].get("root") if item in self.projects else self.folders.get(item)

    def _shown_children(self, item):
        return [self.entries.get(k) or self.folders.get(k) for k in self.tree.get_children(item)]

    def _place(self, item, parent, index):
        """Put `item` in the list under `parent`, at position `index` among parent's other children."""
        others = [k for k in self.tree.get_children(parent) if k != item]
        self.tree.move(item, parent, self.tree.index(others[index]) if index < len(others) else tk.END)

    def _structure_changed(self, members=False):
        """After images were dragged or deleted (or that was undone): the folders and projects that now list something other than
        their files do ("unsaved" after the name; Save), the thumbnails' order, and (members: images came or went) the slide."""
        for item, obj in list(self.folders.items()) + [(pid, p.get("root")) for pid, p in self.projects.items()]:
            if obj is None:
                continue
            obj.order_changed = self._shown_children(item) != obj.children
            if item in self.projects:
                self.tree.item(item, text=self.projects[item]["name"])
                self._place_unsaved_mark(item)
            else:
                self._label_item(item)
        self._update_annotations_check()
        self.thumb_grid.rows = self._grid_rows_all()
        self.thumb_grid.layout()
        self.root.after_idle(self._place_list_bars)         # (rows moved: their progress bars with them)
        self._show_sizes()                                  # (the folders' and projects' totals)
        if members and self.slide_project in self.projects:
            self._start_slide(self.slide_project)

    def delete_section(self, item):
        """Right-click > Delete section: out of the list at once; Save takes it out of LAS X's files and moves its files (its .xlif,
        channel TIFFs and reduced copies) to the Recycle Bin. Undo (Ctrl+Z) brings it back until then."""
        entry = self.entries.get(item)
        parent = self.tree.parent(item)
        container = self._container(parent)
        if entry is None or container is None or container.path is None:
            messagebox.showerror("Can't delete", "The file that lists this image is unreadable, so it can't be deleted here.",
                                 parent=self.root)
            return
        index = self.tree.index(item)
        if self.img is not None and self.img.entry is entry:
            self._close_image()
            self.show_message("Select an image")
        kept = self.kept_images.pop(entry, None)
        if kept is not None:
            kept.close()
            self._kept_changed()
        if item in self.download_queue:
            self.download_queue.remove(item)
        self.tree.detach(item)
        self.entries.pop(item)
        self.deleted[item] = entry
        self.history = self.history[-199:] + [("delete", (item, entry, parent, index))]
        self._structure_changed(members=True)
        self.status.config(text=f"Deleted {entry.name}: not saved yet (Save moves its files to the Recycle Bin; Ctrl+Z undoes it)")

    def _drop_ok(self, item, parent):
        """Can `item` go under `parent`? Not a folder into itself (or into what it holds), nor into another project (only images move
        between projects), nor anything to or from a folder whose own file is unreadable."""
        p = parent
        while p:
            if p == item:
                return False
            p = self.tree.parent(p)
        if item in self.folders and self._project_of(parent) != self._project_of(item):
            return False
        target, source = self._container(parent), self._container(self.tree.parent(item))
        return not (target is None or target.path is None or source is None or source.path is None)

    def _drop_target(self, item, y):
        """Where dropping `item` at height y of the list puts it: (parent, index among the parent's other children), or None. On the
        upper or lower half of a name: before or after it; on the middle of a folder's name, or on a project's: into it, at the end."""
        row = self.tree.identify_row(y)
        box = self.tree.bbox(row) if row and row != item else None
        if not box:
            return None
        rel = (y - box[1]) / max(box[3], 1)
        if row in self.projects or (row in self.folders and 0.25 <= rel <= 0.75):
            parent, index = row, len([k for k in self.tree.get_children(row) if k != item])
        else:
            parent = self.tree.parent(row)
            others = [k for k in self.tree.get_children(parent) if k != item]
            index = others.index(row) + (rel >= 0.5)
        return (parent, index) if parent and self._drop_ok(item, parent) else None

    def _show_drop(self, target):
        """A line in the list where a dragged name would go (None: no line)."""
        if target is None:
            if self.drop_line is not None:
                self.drop_line.place_forget()
            return
        parent, index = target
        others = [k for k in self.tree.get_children(parent) if k != self.drag_item[0]]
        if index < len(others):
            box, edge = self.tree.bbox(others[index]), 0
        else:
            last = others[-1] if others else parent
            while others and self.tree.item(last, "open") and self.tree.get_children(last):
                last = self.tree.get_children(last)[-1]
            box, edge = self.tree.bbox(last), 1
        if not box:
            self._show_drop(None)
            return
        if self.drop_line is None:
            self.drop_line = tk.Frame(self.tree, height=3, background=PAL["selected"])
        x = 20 * (self._tree_depth(parent) + 1) if parent else 0
        self.drop_line.place(x=x, y=box[1] + edge * box[3] - 1, width=max(self.tree.winfo_width() - x, 20))

    def _tree_drag(self, event):
        """Dragging a name in the list (after a few pixels: a short wobble is still a click): a line shows where it would go."""
        if self.drag_item is None:
            return
        item, y0, started = self.drag_item
        if not started:
            if abs(event.y - y0) < 6:
                return
            self._cancel_rename_job()
            self._hide_tree_tip()
            self.drag_item = (item, y0, True)
            self.tree.configure(cursor="fleur")
        if event.y < 12:
            self.tree.yview_scroll(-1, "units")
        elif event.y > self.tree.winfo_height() - 12:
            self.tree.yview_scroll(1, "units")
        self._show_drop(self._drop_target(item, event.y))
        return "break"

    def _tree_drop(self, event):
        drag, self.drag_item = self.drag_item, None
        if drag is None or not drag[2]:
            return
        self._show_drop(None)
        self.tree.configure(cursor="")
        target = self._drop_target(drag[0], event.y)
        if target is not None:
            self.move_item(drag[0], *target)
        return "break"

    def _grid_drop(self, item, target, after):
        """A thumbnail dropped before (or after) another: it goes there in the list too."""
        parent = self.tree.parent(target)
        others = [k for k in self.tree.get_children(parent) if k != item]
        if target in others and self._drop_ok(item, parent):
            self.move_item(item, parent, others.index(target) + after)

    def move_item(self, item, parent, index):
        """Move an image (or a folder, within its project) to `parent`, at `index` among the others there: shown at once, written by
        Save (into another project: its files move there too). Undoable."""
        old = (self.tree.parent(item), self.tree.index(item))
        if old == (parent, index):
            return
        crossed = self._project_of(item) != self._project_of(parent)
        self._place(item, parent, index)
        self.history = self.history[-199:] + [("move", (item,) + old)]
        self._structure_changed(members=crossed)
        self.tree.see(item)
        obj = self.entries.get(item) or self.folders.get(item)
        away = self._project_of(item) != self.home.get(item)
        self.status.config(text=f"Moved {obj.name}: not saved yet (Save writes the new order"
                                + (", and moves its files into that project" if away else "") + "; Ctrl+Z undoes it)")

    def _save_structure(self, failed):
        """Save, for the list's order: images dragged into another project move there, files and all (see plan_move_image); each
        folder and project whose contents or their order changed gets its list of references rewritten (see children_text); deleted
        images go to the Recycle Bin (on a drive that has none, a network drive, into a "Deleted by XLEF viewer" folder of their
        experiment instead). Returns the projects to read again, and how many images were deleted where (bin, folder)."""
        touched, binned, kept_aside = set(), 0, 0
        old_path = {}           # ImageEntry -> its .xlif's path before it moved (how the list it left still refers to it)
        # 1. images into another project: their files first (a container then refers only to what is really there)
        for item, entry in list(self.entries.items()):
            pid, home = self._project_of(item), self.home.get(item)
            if pid is None or home is None or pid == home:
                continue
            touched |= {pid, home}
            try:
                plan, moved = plan_move_image(entry, os.path.dirname(self.projects[home]["path"]), self._container(
                    self.tree.parent(item)).path, os.path.dirname(self.projects[pid]["path"]))
                carry_out(plan)
            except Exception as err:
                failed.append(f"moving {entry.name}: {err}")
                folder = next((k for k, o in [(home, self._container(home))] + list(self.folders.items())
                               if o is not None and entry in o.children), home)
                self._place(item, folder, self._container(folder).children.index(entry)
                            if entry in self._container(folder).children else 0)
                continue
            new = lambda p: moved.get(os.path.normcase(p), p)       # noqa: E731
            old_path[entry] = entry.xlif_path
            entry.xlif_path = new(entry.xlif_path)
            entry.files = {k: new(f) for k, f in entry.files.items()}
            entry.pyramid = new(entry.pyramid) if entry.pyramid else None
        # 2. the lists of references
        written = set()
        for item, obj in list(self.folders.items()) + [(pid, p.get("root")) for pid, p in self.projects.items()]:
            if obj is None:
                continue
            shown = self._shown_children(item)
            if shown == obj.children:
                continue
            touched.add(self._project_of(item) if item not in self.projects else item)
            try:
                text = children_text(obj.path, [old_path.get(c, c.path) for c in obj.children],
                                     [(c.path, getattr(c, "ref_uuid", "") or str(uuid.uuid1())) for c in shown])
                _write_text(obj.path, text, _read_text(obj.path)[1])
                written.add(item)
                obj.children = shown
                obj.order_changed = False
                for c in shown:
                    c.parent_file = obj.path
            except Exception as err:
                failed.append(f"the list of {obj.name}: {err}")
        # 3. deleted images, once nothing refers to them
        for item, entry in list(self.deleted.items()):
            parent = next((k for k, o in [(p, self._container(p)) for p in self.projects] + list(self.folders.items())
                           if o is not None and entry.xlif_path in (getattr(c, "xlif_path", None) for c in o.children)), None)
            if parent is not None:          # (its list was not rewritten: it stays)
                failed.append(f"deleting {entry.name}: the file that lists it was not changed")
                continue
            touched.add(self.home.get(item))
            try:
                files = [_real_case(f) for f in owned_files(entry.xlif_path) if os.path.exists(f)]
                root = os.path.dirname(self.projects[self.home[item]]["path"])
                if has_recycle_bin(entry.xlif_path):
                    send_to_recycle_bin(files)
                    binned += 1
                else:
                    aside = os.path.join(root, "Deleted by XLEF viewer")
                    carry_out([("mkdir", os.path.dirname(os.path.join(aside, os.path.relpath(f, root))), None) for f in files]
                              + [("move", f, os.path.join(aside, os.path.relpath(f, root))) for f in files])
                    kept_aside += 1
                del self.deleted[item]
            except Exception as err:
                failed.append(f"deleting {entry.name}: {err}")
        return touched - {None}, binned, kept_aside

    def unsaved_ok(self, doing, question):
        """True if there are no unsaved changes (renames, annotations), or the user agrees to lose them (nothing is ever saved without
        Save)."""
        changed = [o.name for o in self._all_objects() if o.unsaved] + [e.name + " (deleted)" for e in self.deleted.values()]
        if not changed:
            return True
        names = ", ".join(changed[:5]) + (f" and {len(changed) - 5} more" if len(changed) > 5 else "")
        return messagebox.askyesno("Unsaved changes",
                                   f"{names} {'has' if len(changed) == 1 else 'have'} unsaved changes.\n\n{doing} discards them. {question}", icon="warning", default="no", parent=self.root)

    def _draw_scale_bar(self):
        self._each_canvas(self._draw_scale_bar_here)

    def _draw_scale_bar_here(self):
        """At the bottom left of the canvas; separated, on the first panel only (top left: the panels share the view's scale)."""
        e = self.img.entry
        if not e.pixel_um or (self.drawing_panel or 0) != 0:
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
            self.status.config(text=f"{e.name}   {size}{z}   zoom {self.scale * 100:.1f}%" + ("      " + self._preload_full_text()
                                                                                          if self.preload_full else ""))

    def _preload_full_text(self):
        n = sum(1 for e in self.preload_full if e not in self.kept_images and not (self.img is not None and self.img.entry is e))
        return (f"⚠ The memory for previews ({size_text(self.kept_limit())}) is full: {n} more section"
                f"{'s were' if n != 1 else ' was'} not preloaded. It can be increased in ⚙ Settings.")

    def _preload_full_note(self):
        """The background loader stopped with the memory for previews full: say so (and again after each image's own line)."""
        if self.img is None or not self.loading:
            self.status.config(text=self._preload_full_text() if self.img is None else self.status.cget("text").split("      ⚠")[0]
                               + "      " + self._preload_full_text())

    def kept_limit(self):
        """Memory (bytes) for the pictures of recently viewed and preloaded images: the share of the computer's memory set in Settings."""
        total = total_memory()
        return self.settings.get("preview_memory_percent", KEPT_IMAGES_PERCENT) / 100 * total if total else 500e6

    def _poll_results(self):
        try:
            while True:
                kind, gen, data = self.results.get_nowait()
                if kind == "thumb":
                    if gen == self.thumb_gen:
                        self.thumb_grid.set_thumb(*data)
                    continue
                if kind == "preload_full":
                    if gen == self.thumb_gen:
                        self.preload_full = data
                        self._preload_full_note()
                    continue
                if kind in ("preload_progress", "preloaded"):
                    if gen == self.thumb_gen:
                        if kind == "preloaded":
                            self._preloaded(*data)
                        else:
                            item, pct = data
                            if pct is None:
                                self.preload_progress.pop(item, None)
                            else:
                                self.preload_progress[item] = pct
                            self._label_item(item)
                            self.thumb_grid.set_progress(item, None if pct is None else pct / 100)
                    continue
                if kind in ("download", "downloaded", "download_error", "download_idle"):
                    if kind == "download_idle":
                        self.download_thread = None
                        if gen == self.download_gen and not (self.download_stop and self.download_stop.is_set()):
                            self._start_download_thread()             # (sections queued while the thread was finishing)
                    elif gen != self.download_gen:
                        pass
                    elif kind == "download":
                        self.download_progress[data[0]] = data[1]
                        self._label_item(data[0])
                    elif kind == "downloaded":
                        self._downloaded(data)
                    else:
                        self.status.config(text=data)
                    continue
                if kind == "experiment":
                    self._experiment_read(gen, *data)          # (gen: the project's tree item)
                    continue
                if kind == "sizes":
                    self._show_sizes()
                    continue
                if kind in ("loaded", "error") and gen == self.generation:
                    self.thumbs_go.set()                # (the chosen image is in: the thumbnails go on)
                if kind == "slide":
                    if gen == self.slide_gen and data is not None:
                        self.slide.set_picture(data)
                    continue
                if kind in ("progress", "error", "loaded", "overview", "quick", "sharp") and self.img is None:
                    continue
                if kind == "progress":
                    self.status.config(text=data)
                elif kind == "overview" and gen == self.generation:
                    self.redraw(preview_only=True)
                elif kind in ("quick", "sharp") and gen == self.generation:
                    # (the overview's arrays and its step change together, so that every channel shown is at one step)
                    if kind == "quick":
                        z, c, a, step = data
                        if self.img.ov_step != step:
                            self.img.overviews, self.img.ov_step = {}, step
                        self.img.overviews[(z, c)] = a
                    else:
                        z, arrays, step = data
                        self.img.overviews = {(z, c): a for c, a in enumerate(arrays)}
                        self.img.ov_step = step
                    self.shrunk.clear()
                    self.redraw(preview_only=True)
                elif kind == "open_progress" and gen == self.generation:
                    self._set_open_progress(data)
                elif kind == "loaded" and gen == self.generation:
                    self._set_open_progress(None)
                    self.loading = False
                    self.redraw()
                    self._update_cursor()
                elif kind == "detail" and gen == self.generation and not self.drag:
                    # (kept for this view, joined to what was kept of it before, and mixed with the sliders as they are NOW)
                    key, arrays, place, smooth = data
                    if self.detail_cache is not None and self.detail_cache["key"] == key:
                        arrays = [a if a is not None else old for a, old in zip(arrays, self.detail_cache["arrays"])]
                    luts = make_luts(self.img.entry.channels, self.ch_settings)
                    rgb = self._mix([a if lut is not None else None for a, lut in zip(arrays, luts)], luts)
                    self.detail_cache = {"key": key, "arrays": arrays, "rgb": rgb,          # (mixed once: reused while zooming)
                                         "sig": self._mix_sig()}
                    self._show(rgb, *place, smooth)
                    self.detail_pending = False
                    self._update_cursor()
                elif kind == "error":
                    self._set_open_progress(None)
                    self.loading = False
                    self.detail_pending = False
                    self._update_cursor()
                    self.status.config(text=data)
        except queue.Empty:
            pass
        self.root.after(30, self._poll_results)

    def save_view(self):
        if self.last_view_rgb is None:
            return
        ims, left, top, (cw, ch), cols = self.last_view
        if len(ims) > 1:                    # (channels separated: the panels as they are on screen, side by side)
            rows = -(-len(ims) // cols)
            tiled = Image.new("RGB", (cols * cw + (cols - 1) * 2, rows * ch + (rows - 1) * 2), (58, 58, 58))
            for k, im in enumerate(ims):
                panel = Image.new("RGB", (cw, ch))
                panel.paste(im, (left, top))
                tiled.paste(panel, ((k % cols) * (cw + 2), (k // cols) * (ch + 2)))
            self.last_view_rgb = tiled
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".png", filetypes=[("PNG", "*.png")],
                                            initialfile=f"{self.img.entry.name}.png")
        if path:
            self.last_view_rgb.save(path)

    def run(self):
        self.root.mainloop()


class LutChooser:
    """Right-click on a channel's swatch or slider: choose its colour map, as in LAS X's dialog. Left, LAS X's maps (LUT_ORDER) with a
    strip showing each; right, a colour picker (hue x saturation, and brightness) for any colour, stored as LAS X does, as
    "Gradient (r,g,b)". The choice shows on the image at once; This section / All sections (the same channel in every section of
    the project) applies it (undoable, written by Save); Cancel, or closing the window, puts the channel back as it was."""

    ROW, STRIP_W, SQUARE = 20, 100, 180

    def __init__(self, viewer, channel):
        self.v, self.i = viewer, channel
        self.entry = viewer.img.entry
        self.ch = self.entry.channels[channel]
        self.original = self.ch.lut
        win = self.win = tk.Toplevel(viewer.root)
        win.title(f"Color map - channel {channel + 1} ({self.ch.name}), {self.entry.name}")
        win.transient(viewer.root)
        win.resizable(False, False)
        win.geometry(f"+{viewer.root.winfo_rootx() + 220}+{viewer.root.winfo_rooty() + 120}")
        body = ttk.Frame(win, padding=10)
        body.pack(fill=tk.BOTH, expand=True)
        left = ttk.Frame(body)
        left.pack(side=tk.LEFT, fill=tk.Y)
        self.list = tk.Canvas(left, width=self.STRIP_W + 130, height=self.ROW * 16, highlightthickness=0, background=PAL["field"])
        bar = ttk.Scrollbar(left, orient=tk.VERTICAL, command=self.list.yview)
        self.list.configure(yscrollcommand=bar.set, scrollregion=(0, 0, 1, self.ROW * len(LUT_ORDER)))
        self.list.pack(side=tk.LEFT, fill=tk.Y)
        bar.pack(side=tk.LEFT, fill=tk.Y)
        self.strips = []
        for k, name in enumerate(LUT_ORDER):
            table, under, over = lut_table(name)
            rgb = (table[np.linspace(0, 255, self.STRIP_W).astype(int)] * 255).astype(np.uint8)
            if under is not None:
                rgb[:2] = under * 255
            if over is not None:
                rgb[-2:] = over * 255
            photo = ImageTk.PhotoImage(Image.fromarray(np.repeat(rgb[None], self.ROW - 6, axis=0)), master=self.list)
            self.strips.append(photo)
            y = k * self.ROW
            self.list.create_image(4, y + 3, anchor=tk.NW, image=photo, tags=(f"row{k}",))
            self.list.create_text(self.STRIP_W + 12, y + self.ROW / 2, anchor=tk.W, text=name, fill=PAL["field_fg"], tags=(f"row{k}",))
        self.list.bind("<Button-1>", lambda e: self._pick_row(self.list.canvasy(e.y)))
        self.list.bind("<MouseWheel>", lambda e: self.list.yview_scroll(int(-e.delta / 120), "units"))
        right = ttk.Frame(body, padding=(14, 0, 0, 0))
        right.pack(side=tk.LEFT, fill=tk.Y)
        ttk.Label(right, text="Any color:").pack(anchor=tk.W)
        self.square = tk.Canvas(right, width=self.SQUARE, height=self.SQUARE, highlightthickness=0, cursor="crosshair")
        self.square.pack(pady=(2, 6))
        h, s = np.meshgrid(np.linspace(0, 255, self.SQUARE), np.linspace(255, 0, self.SQUARE))
        hsv = np.stack([h, s, np.full_like(h, 255)], axis=-1).astype(np.uint8)
        self.square_photo = ImageTk.PhotoImage(Image.fromarray(hsv, "HSV").convert("RGB"), master=self.square)
        self.square.create_image(0, 0, anchor=tk.NW, image=self.square_photo)
        for ev in ("<Button-1>", "<B1-Motion>"):
            self.square.bind(ev, self._pick_square)
        ttk.Label(right, text="Brightness").pack(anchor=tk.W)
        self.value = tk.DoubleVar(win, 100)
        ttk.Scale(right, from_=10, to=100, variable=self.value, command=lambda v: self._custom()).pack(fill=tk.X)
        self.preview = tk.Canvas(right, width=self.SQUARE, height=14, highlightthickness=0)
        self.preview.pack(pady=(8, 2))
        self.name_label = ttk.Label(right, text="")
        self.name_label.pack(anchor=tk.W)
        self.hs = None                              # (the custom colour's hue and saturation, 0-1)
        m = GRADIENT.fullmatch(self.original.strip())
        if m:
            import colorsys
            h_, s_, v_ = colorsys.rgb_to_hsv(*(int(c) / 255 for c in m.groups()))
            self.hs = (h_, s_)
            self.value.set(max(10, v_ * 100))
        buttons = ttk.Frame(win, padding=(10, 0, 10, 10))
        buttons.pack(fill=tk.X)
        ttk.Button(buttons, text="Cancel", command=self.cancel).pack(side=tk.RIGHT)
        ttk.Button(buttons, text="All sections", command=lambda: self.apply(True)).pack(side=tk.RIGHT, padx=(0, 6))
        ttk.Button(buttons, text="This section", command=lambda: self.apply(False)).pack(side=tk.RIGHT, padx=(0, 6))
        ttk.Label(buttons, text="Use for:").pack(side=tk.RIGHT, padx=(0, 8))
        win.protocol("WM_DELETE_WINDOW", self.cancel)
        win.bind("<Escape>", lambda e: self.cancel())
        self._show(self.original)
        win.grab_set()

    def _show(self, name):
        """Show `name` as the choice: on the image, the channel's swatch and slider, the list's highlight and the preview bar."""
        self.ch.lut = name
        v, i = self.v, self.i
        if v.img is not None and v.img.entry is self.entry and i < len(v.ch_widgets):
            on, slider, swatch = v.ch_widgets[i]
            swatch.color = v._swatch_color(i)
            slider.set_color(swatch.color)
            v._show_channel_state(i)
            v.redraw()
        self.list.delete("sel")
        if name in LUT_ORDER:
            k = LUT_ORDER.index(name)
            self.list.create_rectangle(1, k * self.ROW + 1, self.STRIP_W + 128, (k + 1) * self.ROW - 1, outline=PAL["selected"], width=2,
                                       tags="sel")
        table, _, _ = lut_table(name)
        rgb = (table[np.linspace(0, 255, self.SQUARE).astype(int)] * 255).astype(np.uint8)
        self.bar_photo = ImageTk.PhotoImage(Image.fromarray(np.repeat(rgb[None], 14, axis=0)), master=self.preview)
        self.preview.delete("all")
        self.preview.create_image(0, 0, anchor=tk.NW, image=self.bar_photo)
        self.name_label.configure(text=name)

    def _pick_row(self, y):
        k = int(y // self.ROW)
        if 0 <= k < len(LUT_ORDER):
            self._show(LUT_ORDER[k])

    def _pick_square(self, event):
        x, y = min(max(event.x, 0), self.SQUARE - 1), min(max(event.y, 0), self.SQUARE - 1)
        self.hs = (x / (self.SQUARE - 1), 1 - y / (self.SQUARE - 1))
        self._custom()

    def _custom(self):
        """The colour picked in the square, at the brightness set: "Gradient (r,g,b)"."""
        if self.hs is None:
            return
        import colorsys
        r, g, b = (round(c * 255) for c in colorsys.hsv_to_rgb(self.hs[0], self.hs[1], float(self.value.get()) / 100))
        self._show(f"Gradient ({r},{g},{b})")

    def apply(self, everywhere):
        """Use the choice for this section, or for the same channel of every section of this project (an image whose channel has the
        map already is left as it is). Undoable (Ctrl+Z), and written into each image's .xlif by Save."""
        name = self.ch.lut
        changes = [(self.entry, self.i, self.original)] if name != self.original else []
        if everywhere:
            pid = self.v._project_of(self.v._item_of(self.entry))
            for e in self.v._entries_of(pid).values() if pid else []:
                if e is not self.entry and len(e.channels) > self.i and e.channels[self.i].lut != name:
                    changes.append((e, self.i, e.channels[self.i].lut))
                    e.channels[self.i].lut = name
        self.win.destroy()
        if changes:
            self.v.history = self.v.history[-199:] + [("lut", changes)]
            self.v._luts_changed({e for e, _, _ in changes})
            n = len({e for e, _, _ in changes})
            self.v.status.config(text=f"Channel {self.i + 1} shown as {name} in {n} section{'s' * (n != 1)}: not saved yet "
                                      f"(Save writes it; Ctrl+Z undoes it)")

    def cancel(self):
        self._show(self.original)
        self.win.destroy()


def ask_download(parent, title, message, all_label="Download for all sections"):
    """Ask about downloading full resolution: "this" (this section), "all" (all sections), or None (Cancel, Escape, closed)."""
    dlg = tk.Toplevel(parent)
    dlg.title(title)
    dlg.transient(parent)
    dlg.resizable(False, False)
    body = ttk.Frame(dlg, padding=14)
    body.pack(fill=tk.BOTH, expand=True)
    ttk.Label(body, text=message, wraplength=440, justify=tk.LEFT).pack(anchor=tk.W)
    bar = ttk.Frame(body)
    bar.pack(fill=tk.X, pady=(14, 0))
    result = []

    def pick(value):
        result.append(value)
        dlg.destroy()
    ttk.Button(bar, text="Cancel", command=lambda: pick(None)).pack(side=tk.RIGHT)
    ttk.Button(bar, text=all_label, command=lambda: pick("all")).pack(side=tk.RIGHT, padx=(0, 6))
    this = ttk.Button(bar, text="Download for this section", command=lambda: pick("this"))
    this.pack(side=tk.RIGHT, padx=(0, 6))
    this.focus_set()
    dlg.bind("<Return>", lambda e: pick("this"))
    dlg.bind("<Escape>", lambda e: pick(None))
    dlg.geometry(f"+{parent.winfo_rootx() + 160}+{parent.winfo_rooty() + 140}")
    dlg.grab_set()
    parent.wait_window(dlg)
    return result[0] if result else None


def ask_text(parent, title, initial, color=None):
    """A small dialog for an annotation's text (several lines allowed). The text, or None if cancelled. Given a `color`, the dialog
    also has a Color button (the color picker), and gives (text, color), or None if cancelled."""
    dlg = tk.Toplevel(parent)
    dlg.title(title)
    dlg.transient(parent)
    dlg.geometry(f"+{parent.winfo_rootx() + 160}+{parent.winfo_rooty() + 140}")    # (over the viewer, on its monitor)
    box = tk.Text(dlg, width=40, height=5, font=("Arial", 11), undo=True)
    box.insert("1.0", initial)
    box.pack(fill=tk.BOTH, expand=True, padx=8, pady=(8, 4))
    ttk.Label(dlg, text="Enter starts a new line; Ctrl+Enter is OK.", style="Muted.TLabel").pack(anchor=tk.W, padx=8)
    result = []
    chosen = [color]

    def ok(event=None):
        result.append(box.get("1.0", "end-1c"))
        dlg.destroy()
        return "break"
    bar = ttk.Frame(dlg)
    bar.pack(fill=tk.X, padx=8, pady=8)
    if color is not None:
        def pick():
            from tkinter import colorchooser
            _, c = colorchooser.askcolor(color=chosen[0], parent=dlg, title="Text color")
            if c:
                chosen[0] = c
                swatch.configure(bg=c, activebackground=c)
        ttk.Label(bar, text="Color").pack(side=tk.LEFT)
        swatch = tk.Button(bar, text="   ", bg=color, activebackground=color, bd=2, cursor="hand2", command=pick)
        swatch.keep_colors = True
        swatch.pack(side=tk.LEFT, padx=(4, 0))
        dlg.color_button = swatch           # (for tests)
    ttk.Button(bar, text="Cancel", command=dlg.destroy).pack(side=tk.RIGHT)
    ttk.Button(bar, text="OK", command=ok).pack(side=tk.RIGHT, padx=(0, 6))
    box.bind("<Control-Return>", ok)
    dlg.bind("<Escape>", lambda e: dlg.destroy())
    box.focus_set()
    dlg.grab_set()
    parent.wait_window(dlg)
    if not result:
        return None
    return result[0] if color is None else (result[0], chosen[0])


if __name__ == "__main__":
    XlefViewer(sys.argv[1] if len(sys.argv) > 1 else None).run()
