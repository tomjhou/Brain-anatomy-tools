"""
Generate a UMAP embedding from the Allen Institute MERFISH (C57BL6J-638850) dataset.

Pipeline:
1. Load the raw MERFISH AnnData in backed mode (as in the user's snippet).
2. Prompt for a cell type subset (Neurons / NonNeurons / All) and a brain
   section subset. Section selection first tries a clickable GUI (a grid of
   per-section spatial thumbnails, with non-neurons drawn gray and neurons
   colored by class on top), falling back to a console prompt if no
   interactive display is available. In the GUI, double-clicking a thumbnail
   opens a single-section ROI picker (after a confirmation, since it takes a
   few seconds) to hand-draw one or more rectangular regions of interest;
   ROIs can be collected across multiple sections in one session and, if
   any are drawn, take priority over the whole-section selection.
3. Bring counts into memory (backed='r' data must be loaded before most scanpy ops).
4. Merge in cell metadata (class/subclass/cluster/section labels) from the ABC cache.
5. Filter to the selected cell type, section(s), and ROIs (if any).
6. Subsample if needed, then standard scanpy preprocessing: filter, normalize, log-transform.
7. PCA -> neighbors -> UMAP.
8. Save the UMAP coordinates and a plot to disk (plus, if ROIs were drawn, a
   CSV of their coordinates and a PNG map showing them on their sections).
9. Highlight subclass 268 on the full UMAP.

Every output file for a run is saved into its own subfolder of `out_folder`,
named 'umap_{run_suffix}' (the cell type/section/ROI selection encoded as a
filename-safe string), with short fixed names within it (umap_coords.csv,
class_plot.png, subclass268.png, and, if ROIs were drawn, roi_coords.csv and
roi_map.png; if the optional subclass/supertype plots are generated, also,
for each of 'subclass'/'supertype', a 'umap_by_{subclass,supertype}'
subfolder of UMAP group*.png files and a 'spatial_maps_{subclass,supertype}'
subfolder of spatial group*.png files) — so distinct runs don't clutter
`out_folder` with many similarly-prefixed files.
Every plot PNG is saved alongside a same-named .svg twin.

Adjust `abc_cache` setup at the top to match however you're already initializing it
(e.g. via `abc_atlas_access.abc_atlas_cache.AbcProjectCache`).
"""

import os
import re
import sys
import gc
import json
import math
import time
import socket
import hashlib
import shutil
import warnings
import platform
import threading
import subprocess
from pathlib import Path
import importlib.util
from app_identity import use_app_identity

use_app_identity('abc_atlas_browser')   # its own taskbar icon and button, not Python's (Windows)

_STARTUP_T0 = time.perf_counter()


def startup_note(message):
    """One line on the console while the app is starting up, with the seconds since it started: the first window can take a while to appear
    (the libraries are large), and a silent console looks like a hang."""
    print(f"[{time.perf_counter() - _STARTUP_T0:5.1f}s] {message}", flush=True)


def _make_windows_dpi_aware():
    """On Windows, declare this process DPI-aware before any Tk window gets
    created. Without this, under Windows display scaling (125%/150%/...), Tk
    reports screen dimensions (winfo_screenheight() etc.) in a different
    pixel space than what actually ends up on screen, so sizing a window
    from those numbers (see compute_figsize_for_screen_height) produces a
    window visibly smaller than intended — the two spaces disagree by
    roughly the scaling factor. Must run before the first Tk root/figure is
    created, so this is called at import time, right here."""
    if sys.platform != 'win32':
        return
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)  # PROCESS_SYSTEM_DPI_AWARE
    except Exception:
        try:
            import ctypes
            ctypes.windll.user32.SetProcessDPIAware()
        except Exception:
            pass


# Before check_libraries_or_exit(), which can open a dialog: if any Tk window
# exists before this call, Tk keeps reporting the scaled screen size for the
# rest of the process (e.g. 3413x1440 instead of 5120x2160 at 150%), and every
# window opens that much too large.
_make_windows_dpi_aware()


# Checked before the imports below, at every launch. Import name -> pip
# package spec (they differ for some). Specs are plain argv strings, not
# shell-quoted; _display_spec adds quotes where a user would need them.
REQUIRED_LIBRARIES = {
    'numpy': 'numpy',
    'pandas': 'pandas',
    'matplotlib': 'matplotlib',
    'PIL': 'pillow',
    'cv2': 'opencv-python',
    'scipy': 'scipy',
    'anndata': 'anndata',
    'scanpy': 'scanpy',
    # (installed from GitHub's zip of the repository's main branch, so that git is not needed: pip reads the zip like a package file)
    'abc_atlas_access': 'abc_atlas_access[notebooks] @ https://github.com/alleninstitute/abc_atlas_access/archive/refs/heads/main.zip',
}
# The app runs without these; only the named feature is unavailable.
# Import name -> (pip spec, what it's needed for).
OPTIONAL_LIBRARIES = {
    'igraph': ('igraph', 'Leiden clustering of new UMAP runs'),
    'pyarrow': ('pyarrow', 'Parquet file support'),
    'plotly': ('plotly', 'the 3D view in the browser'),
    'skimage': ('scikit-image', 'the translucent-surfaces background of the browser 3D view'),
    'tkinterdnd2': ('tkinterdnd2', 'dragging a CSV file into the 3D viewer window'),
}
# Install anything missing automatically, without asking (into the interpreter
# running this script), instead of just listing the commands and quitting. Set False to get
# the old list-the-commands-and-quit behavior.
AUTO_INSTALL_LIBRARIES = True
# The packaged program (packaging/atlas_entry.py) has every library inside and no pip to install more with.
RUNNING_PACKAGED = getattr(sys, 'frozen', False)


def _display_spec(spec):
    """A pip spec as it should be typed in a shell: quoted if it has spaces."""
    return f'"{spec}"' if ' ' in spec else spec


def _is_missing_library(name):
    try:
        return importlib.util.find_spec(name) is None
    except (ImportError, ValueError):
        return True


def _install_libraries(specs):
    """pip-install each spec, one pip call apiece (so one failure doesn't
    block the rest), into the interpreter running this script. Output goes to
    the console and, when Tk is available, a small progress window. pip runs
    in a subprocess whose output a reader thread drains into a queue, so the
    window keeps repainting even during pip's long quiet stretches (building
    wheels, resolving dependencies). Returns [(spec, error_tail), ...] for
    the ones that failed."""
    import queue
    ui = None
    try:
        import tkinter as tk
        root = tk.Tk()
        root.title("Installing Python libraries")
        root.attributes('-topmost', True)
        root.protocol("WM_DELETE_WINDOW", lambda: None)  # closing mid-pip would only hide progress
        label = tk.Label(root, text="", justify='left', anchor='w', width=80, padx=12, pady=12)
        label.pack()
        root.update()
        ui = (root, label)
    except Exception:
        ui = None

    def show(text):
        nonlocal ui
        if ui is None:
            return
        try:
            ui[1].config(text=text)
            ui[0].update()
        except Exception:
            ui = None  # window gone; carry on with console output only

    failures = []
    creationflags = subprocess.CREATE_NO_WINDOW if sys.platform == 'win32' else 0
    for i, spec in enumerate(specs, 1):
        header = f"Installing {_display_spec(spec)} ({i} of {len(specs)})..."
        print(header)
        show(header + "\n\nThis can take a few minutes; the console shows pip's full output.")
        try:
            proc = subprocess.Popen(
                [sys.executable, '-m', 'pip', 'install', spec],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
                creationflags=creationflags,
            )
        except OSError as e:
            failures.append((spec, str(e)))
            continue
        lines = queue.Queue()

        def drain(stream=proc.stdout, out=lines):
            for line in stream:
                out.put(line)
            out.put(None)  # end-of-output sentinel

        threading.Thread(target=drain, daemon=True).start()
        tail = []
        while True:
            try:
                line = lines.get(timeout=0.1)
            except queue.Empty:
                line = ''
            if line is None:
                break
            if line:
                print(line, end='')
                tail.append(line.rstrip())
            # Repainted on every pass, including empty ones, so the window
            # stays responsive during pip's quiet stretches.
            show(header + ("\n\n" + tail[-1][:100] if tail else ""))
        if proc.wait() != 0:
            failures.append((spec, "\n".join(tail[-8:])))
    if ui is not None:
        try:
            ui[0].destroy()
        except Exception:
            pass
    importlib.invalidate_caches()  # so find_spec/import see what pip just wrote
    return failures


def check_libraries_or_exit():
    """Look for every library this script needs, without importing them
    (importlib.util.find_spec only locates them, so this costs milliseconds).
    Anything missing, required or optional, is pip-installed at once
    (AUTO_INSTALL_LIBRARIES), with a progress window. If installing fails: missing required libraries show
    the manual install commands and quit, since the imports below would fail
    anyway; missing optional ones only ask whether to continue. Manual
    commands are plain pip; we can't tell whether the user runs conda or a
    base Python, so the message says to run them in whichever environment
    launches this script, and names that interpreter. Everything is also
    printed to the console, where it can be copied."""
    missing_required = [name for name in REQUIRED_LIBRARIES if _is_missing_library(name)]
    missing_optional = [name for name in OPTIONAL_LIBRARIES if _is_missing_library(name)]
    if not missing_required and not missing_optional:
        return

    failures = []
    if AUTO_INSTALL_LIBRARIES:
        # No question first: whoever starts this from a .bat file wants it to just work.
        print("These Python libraries are not installed, so they will be installed now:")
        for name in missing_required:
            print(f"    {name}")
        for name in missing_optional:
            print(f"    {name} (optional: {OPTIONAL_LIBRARIES[name][1]})")
        specs = ([REQUIRED_LIBRARIES[name] for name in missing_required]
                 + [OPTIONAL_LIBRARIES[name][0] for name in missing_optional])
        print("Installing missing Python libraries into:", sys.executable)
        failures = _install_libraries(specs)
        # Re-check rather than trusting pip's exit codes alone: a package can
        # install fine yet still not be importable under the name we expect.
        missing_required = [name for name in REQUIRED_LIBRARIES if _is_missing_library(name)]
        missing_optional = [name for name in OPTIONAL_LIBRARIES if _is_missing_library(name)]
        if not missing_required and not missing_optional:
            return

    lines = []
    if AUTO_INSTALL_LIBRARIES:
        lines += ["Automatic installation did not succeed for everything.", ""]
        for spec, error_tail in failures:
            lines += [f"pip install {_display_spec(spec)} failed:"]
            lines += [f"    {row}" for row in error_tail.splitlines()[-4:]]
            lines += [""]
    if missing_required:
        lines += ["These required Python libraries are not installed:", ""]
        lines += [f"    {name}" for name in missing_required]
        lines += [""]
    if missing_optional:
        lines += ["These optional Python libraries are not installed:", ""]
        lines += [f"    {name} (needed for {OPTIONAL_LIBRARIES[name][1]})" for name in missing_optional]
        lines += [""]
    lines += ["Install them with:", ""]
    lines += [f"    python -m pip install {_display_spec(REQUIRED_LIBRARIES[name])}" for name in missing_required]
    lines += [f"    python -m pip install {OPTIONAL_LIBRARIES[name][0]}" for name in missing_optional]
    lines += ["", "Run these in the same Python environment you use to launch this app",
              "(e.g. after 'conda activate <env>' if you use conda). This app is running from:",
              f"    {sys.executable}"]
    if missing_required:
        lines += ["", "The app will now quit."]
    else:
        lines += ["", "Continue without them?"]
    message = "\n".join(lines)
    print("\n" + message + "\n")

    title = "Missing Python libraries"
    try:
        import centered_dialog                                  # (in the middle of the screen: see centered_dialog.py)
        if missing_required:
            centered_dialog.info(title, message)
            proceed = False
        else:
            proceed = centered_dialog.ask_yes_no(title, message)
    except Exception:
        # No GUI available: fall back to the console.
        if missing_required:
            proceed = False
        else:
            try:
                proceed = input("Continue? [y/N] ").strip().lower() in ('y', 'yes')
            except EOFError:
                proceed = False
    if not proceed:
        sys.exit(1)


if not RUNNING_PACKAGED:
    check_libraries_or_exit()

startup_note("Loading the libraries (numpy, pandas, scanpy, matplotlib, ...). The first start after a reboot is the slowest.")
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from PIL import Image
import matplotlib
# Every window in this app is built on Tk: its sizing (compute_figsize_for_
# screen_height, normalize_tk_scaling), blitting, cursors and dialogs all
# assume the TkAgg backend. Without this, matplotlib picks whichever backend
# it finds first, e.g. Qt when PyQt is installed, and Qt's own display
# scaling then makes every window open too large. Must come before pyplot.
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.widgets import Button, TextBox
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
# Aliased — pathlib.Path is already imported as plain `Path` and used
# throughout this file for filesystem paths; this is matplotlib's own,
# unrelated Path (a vector shape, used below for the multi-gene overlay's
# selectable marker shapes).
from matplotlib.path import Path as MplPath

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
from abc_atlas_access.abc_atlas_cache.cloud_cache import OutdatedManifestWarning

# The data and rendering layer lives in abc_core.py (importable without starting this app); everything moved there is imported back
# under its own name, so the rest of this file is unchanged. Must come after the library check and the matplotlib setup above.
import abc_project
from abc_core import (
    BREGMA_CCF_ML_MM,
    CACHE_DIR, CCF_COORDS_FILE_NAME, CCF_DIRECTORY, CCF_FILES, CCF_PARCELLATION_DIRECTORY, CCF_PARCELLATION_FILE_NAME, COORD_SETTINGS,
    DEFAULT_BREGMA_CCF_AP_MM, DEFAULT_BREGMA_CCF_DV_MM, DEFAULT_GREY_RGBA, FRAME_STATE, RECON_FILES, coordinate_frame_name,
    file_complete, reconstructed_in_use, record_run_opened as core_record_run_opened, IMPUTED_DATASET_DIRECTORY,
    IMPUTED_DATASET_FILE_NAME, NON_NEURON_CLASS_IDS, DEFAULTS_DIR, SECTION_COL, SETTINGS_FILE, SPATIAL_CSV_CHUNK_SIZE,
    build_neuron_class_colors, ccf_files_missing, class_id_to_color, ensure_spatial_load_started, extract_leading_numeric_id,
    find_gene_index, generate_and_cache_section_image, generate_section_grid_image, get_section_labels, load_annotations_aligned,
    load_coordinate_settings, load_imputed_adata, load_rois_csv, load_section_spatial_coords, metadata_file_status, multi_gene_rgb,
    render_gene_expression_array, render_section_level_array, sanitize_section_token, save_coordinate_settings,
    scatter_gray_then_colored, sorted_sections_descending, start_section_coordinate_load, write_run_statistics,
)



# Every interactive window in this file has its own hand-rolled pan/zoom/
# hover controls (right-mouse-drag panning, scroll-to-zoom, ROI dragging,
# hover-for-cell-info, ...) — matplotlib's default NavigationToolbar2Tk
# (home/pan/zoom/save icons + coordinate readout, packed as a separate Tk
# widget below the canvas) is both redundant and, once a window is sized
# large relative to its original layout (see maximize_figure_window), prone
# to visually crowding/overlapping the custom button row at the bottom of
# the figure — those buttons are positioned in the *figure's* own fractional
# coordinate space, which has no notion of the toolbar's separate,
# fixed-pixel-height Tk strip. Disabling it removes that whole class of
# layout conflict. Must be set before any figure is created.
matplotlib.rcParams['toolbar'] = 'None'


def _has_internet_connection(host='allen-brain-cell-atlas.s3.amazonaws.com', port=443, timeout=2.0):
    """Best-effort, fast (few-second) connectivity probe — a plain TCP
    connect attempt, not an HTTP/S3 request. AbcProjectCache.from_cache_dir
    (below) picks an S3-backed cache over a local-only one any time the
    cache directory is merely *writable* — true for basically every normal
    setup, even when every dataset this session will ever touch is already
    fully downloaded — so it doesn't matter whether the data is cached: the
    S3 client still gets constructed and still calls out to list/verify
    the manifest. Without internet, that call doesn't fail fast: boto3's
    own default connect/read timeouts (60s each) plus its default retry
    count meant the app sat there for minutes, apparently frozen, before
    any window ever appeared — this whole check exists to detect that
    upfront, in seconds, and skip straight to the local-only path instead.
    The probe is of the atlas's own bucket, not of s3.amazonaws.com: the
    general S3 address can be reachable while the bucket's own name does
    not resolve (a flaky DNS or VPN), which is the case that failed."""
    try:
        socket.create_connection((host, port), timeout=timeout).close()
        return True
    except OSError:
        return False



CACHE_DIR.mkdir(parents=True, exist_ok=True)

# A crash that ends the program outright (a native abort, as when Tk is used from the wrong thread) leaves its Python traceback here, since the console
# window closes with the program.
try:
    import faulthandler
    _crash_log = open(CACHE_DIR / 'crash_log.txt', 'a', buffering=1)
    faulthandler.enable(_crash_log, all_threads=True)
except Exception:
    pass


# Defaults for the interactive UMAP viewer's own Settings dialog. Gene456_
# shape is validated against GENE456_SHAPE_REGISTRY (defined further down,
# alongside the rest of the multi-gene overlay's own tunables) inside
# load_viewer_settings itself, not here — referencing it only from within
# a function body means definition order in the file doesn't matter, since
# by the time load_viewer_settings is actually called (opening the
# interactive viewer), the whole module has already finished loading.
DEFAULT_VIEWER_SETTINGS = {
    'hover_enabled': True,
    'gene456_shape': 'small circle',
    'umap_circle_size': 'Default',
    'colorize_other_cells': True,   # in a neuron project the non-neurons (and in a non-neuron project the neurons) of the section maps are colored too; off: always gray
    'group_colors': 'cycle',        # 'cycle': every group colored, the look repeating past 400; 'lineage': a group's color comes from its parent's
    # Itself a setting, with a deliberately asymmetric persistence rule —
    # see load_viewer_settings/save_viewer_settings's own docstrings.
    'persist_settings': False,
}


def load_viewer_settings():
    """The interactive UMAP viewer's own settings (hover-highlighting on/
    off, the multi-gene overlay's genes-4-6 marker shape, the overall UMAP
    dot size), starting from DEFAULT_VIEWER_SETTINGS and layering in
    whatever's on disk.

    'persist_settings' itself is the one exception to "nothing persists by
    default": it's always read from the file if present, regardless of its
    own value, so a user who opts into persistence once stays opted in on
    every later launch without the choice to persist needing to itself be
    persisted through some separate mechanism. The other settings are only
    loaded from disk when the file's own 'persist_settings' says to —
    otherwise they're left at their hardcoded defaults for this session,
    same as if the file didn't exist.

    Never raises — any read/parse problem just falls back to defaults, the
    same as no file existing at all."""
    settings = dict(DEFAULT_VIEWER_SETTINGS)
    try:
        with open(SETTINGS_FILE, 'r', encoding='utf-8') as f:
            on_disk = json.load(f)
    except Exception:
        return settings
    if not isinstance(on_disk, dict):
        return settings
    settings['persist_settings'] = bool(on_disk.get('persist_settings', False))
    if settings['persist_settings']:
        if 'hover_enabled' in on_disk:
            settings['hover_enabled'] = bool(on_disk['hover_enabled'])
        if on_disk.get('gene456_shape') in GENE456_SHAPE_REGISTRY:
            settings['gene456_shape'] = on_disk['gene456_shape']
        if on_disk.get('umap_circle_size') in UMAP_CIRCLE_SIZE_REGISTRY:
            settings['umap_circle_size'] = on_disk['umap_circle_size']
        if 'colorize_other_cells' in on_disk:
            settings['colorize_other_cells'] = bool(on_disk['colorize_other_cells'])
        if on_disk.get('group_colors') in ('cycle', 'lineage'):
            settings['group_colors'] = on_disk['group_colors']
    return settings


def save_viewer_settings(settings):
    """Writes `settings` to SETTINGS_FILE — mirroring load_viewer_
    settings' own asymmetric rule: 'persist_settings' is always written,
    but the other settings are only included when settings[
    'persist_settings'] is True. Without that second half, turning
    persistence off would still leave whatever was last saved sitting in
    the file, silently reappearing (looking like persistence was somehow
    still partially on) the next time persistence got turned back on,
    instead of the user's *current* choices at the time they re-enabled
    it.

    Best-effort: a write failure (read-only install location, etc.) is
    swallowed rather than surfaced — failing to persist a preference isn't
    worth interrupting the session over."""
    to_write = {'persist_settings': bool(settings.get('persist_settings', False)),
                'bregma_ap_mm': COORD_SETTINGS['bregma_ap_mm'], 'bregma_dv_mm': COORD_SETTINGS['bregma_dv_mm'],
                'reconstructed_coords': bool(COORD_SETTINGS['reconstructed_coords'])}   # (always kept: see COORD_SETTINGS)
    if to_write['persist_settings']:
        to_write['hover_enabled'] = bool(settings.get('hover_enabled', True))
        to_write['gene456_shape'] = settings.get('gene456_shape', GENE456_SHAPE_DEFAULT)
        to_write['umap_circle_size'] = settings.get('umap_circle_size', UMAP_CIRCLE_SIZE_DEFAULT)
        to_write['group_colors'] = settings.get('group_colors', 'cycle')
        to_write['colorize_other_cells'] = bool(settings.get('colorize_other_cells', True))
    try:
        SETTINGS_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(SETTINGS_FILE, 'w', encoding='utf-8') as f:
            json.dump(to_write, f, indent=2)
    except Exception as e:
        print(f"Warning: could not save viewer settings ({e})")







load_coordinate_settings()











SECTION_Z_LABEL_FONTSIZE = 8    # points: the z label at the top-right of each thumbnail in the section grid (was 6)




def build_common_settings_tab(tab):
    """The "Common" tab of the Settings window (what the launcher and the project browser share): where bregma is in the CCF (for the AP, DV and ML from
    bregma), and whether cells are placed by Allen's reconstructed coordinates. Fills the frame `tab`; returns read(), which gives the values as a dict
    {'bregma_ap_mm', 'bregma_dv_mm', 'reconstructed_coords'}, or raises ValueError with the message to show when an entry is not acceptable."""
    import tkinter as tk
    from tkinter import ttk
    bregma_var = tk.StringVar(tab, f"{COORD_SETTINGS['bregma_ap_mm']:g}")
    bregma_dv_var = tk.StringVar(tab, f"{COORD_SETTINGS['bregma_dv_mm']:g}")
    frame = ttk.Frame(tab, padding=14)
    frame.pack(fill=tk.BOTH)
    row = ttk.Frame(frame)
    row.pack(anchor="w", pady=(0, 0))
    ttk.Label(row, text="Bregma in the CCF (mm):   AP").pack(side=tk.LEFT)
    ttk.Entry(row, textvariable=bregma_var, width=7).pack(side=tk.LEFT, padx=6)
    ttk.Label(row, text="DV").pack(side=tk.LEFT)
    ttk.Entry(row, textvariable=bregma_dv_var, width=7).pack(side=tk.LEFT, padx=6)
    ttk.Label(row, text="ML").pack(side=tk.LEFT)
    ml_entry = ttk.Entry(row, width=7)                          # (the midline: fixed, shown for completeness)
    ml_entry.insert(0, f"{BREGMA_CCF_ML_MM:g}")
    ml_entry.state(["readonly"])
    ml_entry.pack(side=tk.LEFT, padx=6)

    def reset_bregma():
        bregma_var.set(f"{DEFAULT_BREGMA_CCF_AP_MM:g}")
        bregma_dv_var.set(f"{DEFAULT_BREGMA_CCF_DV_MM:g}")

    reset = ttk.Button(row, text="Reset all", command=reset_bregma)
    reset.pack(side=tk.LEFT, padx=(8, 0))
    reset_tip = {'win': None}

    def show_reset_tip(e):
        hide_reset_tip()
        tip = tk.Toplevel(reset)
        tip.wm_overrideredirect(True)
        tip.attributes('-topmost', True)
        tk.Label(tip, text=f"Back to the default: AP {DEFAULT_BREGMA_CCF_AP_MM:g} mm, DV {DEFAULT_BREGMA_CCF_DV_MM:g} mm\n"
                           "(the International Brain Laboratory's estimate of bregma in the CCF)", bg="#ffffe0", relief="solid",
                 borderwidth=1, justify=tk.LEFT, padx=4, pady=1).pack()
        tip.wm_geometry(f"+{e.x_root + 12}+{e.y_root + 18}")
        reset_tip['win'] = tip

    def hide_reset_tip(e=None):
        if reset_tip['win'] is not None:
            reset_tip['win'].destroy()
            reset_tip['win'] = None

    reset.bind("<Enter>", show_reset_tip)
    reset.bind("<Leave>", hide_reset_tip)
    reset.bind("<ButtonPress>", hide_reset_tip, add="+")
    ttk.Label(frame, foreground="gray35", justify=tk.LEFT, wraplength=460, text=(
        "Offsets used to convert common coordinate framework (CCF) brain coordinates to standard bregma-centered (e.g. Paxinos) values. There is no exact mapping between Allen CCF "
        "and Paxinos coordinates: the CCF has no skull, and is not leveled with bregma and lambda at one height (it is pitched down about 5 degrees relative to a flat skull). "
        "Hence, the conversion is approximate. The defaults are the International Brain Laboratory's (IBL's) estimate of bregma in the CCF, the average of five experimenters' picks of "
        "the bregma voxel. To use the bregma position of another registration instead (for example Chon et al. 2019), enter it here.")).pack(
        anchor="w", pady=(6, 0))

    def add_citation(text, links):
        """A line of citation text followed by clickable links [(text, url)]."""
        line = ttk.Frame(frame)
        line.pack(anchor="w", pady=(2, 0))
        ttk.Label(line, foreground="gray35", text=text).pack(side=tk.LEFT)
        for name, url in links:
            link = tk.Label(line, text=name, foreground="#0645ad", cursor="hand2", font=("Segoe UI", 9, "underline"))
            link.pack(side=tk.LEFT, padx=(6, 0))
            link.bind("<Button-1>", lambda e, u=url: __import__("webbrowser").open(u))

    add_citation("Defaults: International Brain Laboratory, iblatlas (Python package):", [("PyPI", "https://pypi.org/project/iblatlas")])
    add_citation("Chon et al. 2019, Nat Commun 10:5067:", [("paper", "https://doi.org/10.1038/s41467-019-13057-w"),
                                                           ("PubMed", "https://pubmed.ncbi.nlm.nih.gov/?term=10.1038%2Fs41467-019-13057-w")])
    ttk.Separator(frame).pack(fill=tk.X, pady=(14, 8))
    recon_var = tk.BooleanVar(tab, bool(COORD_SETTINGS['reconstructed_coords']))
    ttk.Checkbutton(frame, text="Place cells by Allen's reconstructed coordinates", variable=recon_var).pack(anchor="w")
    ttk.Label(frame, foreground="gray35", justify=tk.LEFT, wraplength=460, text=(
        "May be more biologically accurate: it corrects possible tissue distortions in the raw coordinates, and the CCF atlas then lines "
        "up with the cells. One-time download (about 260 MB); applies to the next run you open.")).pack(anchor="w", pady=(2, 0))

    def read():
        try:
            bregma = float(bregma_var.get())
            bregma_dv = float(bregma_dv_var.get())
            if not (0.0 < bregma < 13.5 and -5.0 < bregma_dv < 8.0):
                raise ValueError
        except ValueError:
            raise ValueError("Enter bregma's position as numbers of mm: AP between 0 and 13.5, DV between -5 and 8.")
        return {'bregma_ap_mm': bregma, 'bregma_dv_mm': bregma_dv, 'reconstructed_coords': bool(recon_var.get())}

    return read


def show_coordinate_settings_dialog(abc_cache):
    """The launcher's Settings button: the Settings window (see prompt_viewer_settings_dialog), opened on its Common tab."""
    new = prompt_viewer_settings_dialog(None, load_viewer_settings(), abc_cache=abc_cache, tab="Common")
    if new is not None:
        save_viewer_settings(new)


# Where the ABC atlas download lives. By default it's next to this repository,
# in the parent folder. If it isn't there, the user picks a location once (see
# choose_atlas_cache_location) and the choice is remembered in this file.
DEFAULT_ATLAS_CACHE_DIR = DEFAULTS_DIR / 'abc_atlas_cache'
ATLAS_CACHE_LOCATION_FILE = CACHE_DIR / 'abc_atlas_cache_location.json'
# Free space a drive should have to hold the atlas download. Everything this app downloads, the imputed gene dataset included, is about
# 60 GB; the rest is room for a newer atlas release (a revised file is downloaded beside the old one, up to about 9 GB) and to spare.
ATLAS_CACHE_SPACE_NEEDED_GB = 80
ATLAS_CORE_GB = 9                  # ...of which the core atlas (the cell table and the raw expression data) every session needs
ATLAS_IMPUTED_GB = 50              # ...and the imputed gene dataset, only when imputed genes are asked for


def format_gb(n_bytes):
    return f"{n_bytes / 1e9:,.0f} GB"


def list_local_drives():
    """[(root, free_bytes, total_bytes)] for every local drive, sorted by root.
    On Windows, fixed and removable drives (network, CD and RAM drives are
    left out); elsewhere, just the filesystem root. Drives that can't be read,
    like an empty card reader, are skipped."""
    roots = []
    if sys.platform == 'win32':
        import ctypes
        import string
        DRIVE_REMOVABLE, DRIVE_FIXED = 2, 3
        bitmask = ctypes.windll.kernel32.GetLogicalDrives()
        for i, letter in enumerate(string.ascii_uppercase):
            if bitmask & (1 << i):
                root = f'{letter}:\\'
                if ctypes.windll.kernel32.GetDriveTypeW(root) in (DRIVE_FIXED, DRIVE_REMOVABLE):
                    roots.append(root)
    else:
        roots.append('/')
    drives = []
    for root in sorted(roots):
        try:
            usage = shutil.disk_usage(root)
        except OSError:
            continue
        drives.append((root, usage.free, usage.total))
    return drives


def free_bytes_at(path):
    """Free space on the drive holding `path` (which may not exist yet), or None."""
    path = Path(path)
    for candidate in [path, *path.parents]:
        if candidate.exists():
            try:
                return shutil.disk_usage(candidate).free
            except OSError:
                return None
    return None


def suggest_atlas_cache_dir(drives):
    """The default location if its drive has room; otherwise abc_atlas_cache at the root of the roomiest local drive that does have room; and
    if none has, the default location after all. (Moving to another drive that is short of room too gains little, and on macOS and Linux the
    only drive listed is the filesystem root, /, where a folder cannot be made: the suggestion was one that failed.)"""
    needed = ATLAS_CACHE_SPACE_NEEDED_GB * 1e9
    free_default = free_bytes_at(DEFAULT_ATLAS_CACHE_DIR)
    if free_default is not None and free_default >= needed:
        return DEFAULT_ATLAS_CACHE_DIR
    roomy = [d for d in drives if d[1] >= needed]
    if roomy:
        return Path(max(roomy, key=lambda d: d[1])[0]) / 'abc_atlas_cache'
    return DEFAULT_ATLAS_CACHE_DIR


def looks_like_atlas_cache(folder):
    folder = Path(folder)
    return (folder / '_downloaded_data.json').exists() or (folder / 'metadata').is_dir()


def normalize_chosen_atlas_dir(folder):
    """A folder picked with Browse is used as-is if it already holds an atlas
    download or is named abc_atlas_cache; otherwise abc_atlas_cache is created
    inside it, so picking e.g. E:\\ gives E:\\abc_atlas_cache."""
    folder = Path(folder)
    if folder.name.lower() == 'abc_atlas_cache' or looks_like_atlas_cache(folder):
        return folder
    return folder / 'abc_atlas_cache'


def usable_atlas_dir(folder):
    """`folder` if it can hold the atlas download as it is (it holds one or is empty), else the one subfolder of it that
    holds a download (the cache was moved into a subfolder), else None: it is a folder with other things in it, and the
    download must not be dumped among them."""
    folder = Path(folder)
    try:
        if looks_like_atlas_cache(folder) or not any(folder.iterdir()):
            return folder
        inside = [c for c in folder.iterdir() if c.is_dir() and looks_like_atlas_cache(c)]
    except OSError:
        return None
    return inside[0] if len(inside) == 1 else None


def prompt_atlas_cache_location(suggested, drives):
    """Small Tk dialog: an editable path pre-filled with `suggested`, a Browse
    button, and the free space on the chosen path's drive. Not tkinter's own
    askdirectory, which can't pre-select a folder that doesn't exist yet.
    Returns the chosen Path, or None if cancelled."""
    import tkinter as tk
    from tkinter import filedialog

    needed = ATLAS_CACHE_SPACE_NEEDED_GB * 1e9
    result = {'path': None}
    root = tk.Tk()
    root.title("Choose ABC atlas download location")
    root.attributes('-topmost', True)
    root.resizable(False, False)
    frame = tk.Frame(root, padx=16, pady=12)
    frame.pack(fill='both', expand=True)

    # (no drive with the room for all of it: say what does fit, so that taking the default is still sound advice)
    short = not any(free >= needed for _root, free, _total in drives)
    shortfall = (f"No drive on this computer has {ATLAS_CACHE_SPACE_NEEDED_GB} GB free. That is fine for most use: the core atlas "
                 f"needs about {ATLAS_CORE_GB} GB. The imputed gene dataset (about {ATLAS_IMPUTED_GB} GB) may not fit.\n\n") if short else ""
    tk.Label(frame, justify='left', wraplength=620, text=(
        "No ABC atlas download was found at the default location:\n"
        f"    {DEFAULT_ATLAS_CACHE_DIR}\n\n"
        "Choose where it should go. If this is a new location, the atlas files are "
        f"downloaded there as needed, which can take about {ATLAS_CACHE_SPACE_NEEDED_GB} GB "
        f"({ATLAS_CORE_GB} GB for the core atlas, {ATLAS_IMPUTED_GB} GB for the imputed gene dataset), "
        "so pick a drive with plenty of free space. If you already have a download "
        "elsewhere, choose that folder instead. This is asked only once. The local "
        "drives and their free space are also listed in the console. (The program's own "
        "working files, a few GB, are kept next to the program, not here.)\n\n"
        f"{shortfall}"
        "If you are unsure, accepting the default is likely the best option."
    )).pack(anchor='w')

    path_var = tk.StringVar(value=str(suggested))
    row = tk.Frame(frame, pady=10)
    row.pack(fill='x')
    entry = tk.Entry(row, textvariable=path_var, width=60)
    entry.pack(side='left', fill='x', expand=True)

    def browse():
        current = Path(path_var.get().strip() or str(suggested))
        start = next((str(p) for p in [current, *current.parents] if p.exists()), None)
        chosen = filedialog.askdirectory(parent=root, initialdir=start,
                                         title="Choose a folder for the ABC atlas download")
        if chosen:
            path_var.set(str(normalize_chosen_atlas_dir(chosen)))

    tk.Button(row, text="Browse...", command=browse).pack(side='left', padx=(8, 0))

    space_label = tk.Label(frame, justify='left', anchor='w')
    space_label.pack(fill='x')

    def update_space(*_args):
        text = path_var.get().strip()
        free = free_bytes_at(text) if text else None
        if free is None:
            space_label.config(text="Free space: unknown (drive not found)", fg='firebrick')
        elif looks_like_atlas_cache(text):
            space_label.config(text=f"Existing atlas download found here. Free space: {format_gb(free)}", fg='darkgreen')
        elif free < needed:
            space_label.config(text=f"Free space: {format_gb(free)} (less than the recommended "
                                    f"{ATLAS_CACHE_SPACE_NEEDED_GB} GB)", fg='firebrick')
        else:
            space_label.config(text=f"Free space: {format_gb(free)}", fg='darkgreen')

    path_var.trace_add('write', update_space)
    update_space()

    buttons = tk.Frame(frame, pady=(8))
    buttons.pack(fill='x')

    def ok(_event=None):
        text = path_var.get().strip()
        if not text:
            return
        if free_bytes_at(text) is None:
            # Keep the dialog open rather than failing to create the folder.
            from tkinter import messagebox
            messagebox.showerror("Drive not found", f"Can't find the drive for:\n{text}", parent=root)
            return
        result['path'] = normalize_chosen_atlas_dir(text)       # (a typed folder is not dumped into either)
        root.destroy()

    def cancel(_event=None):
        root.destroy()

    tk.Button(buttons, text="Cancel", width=10, command=cancel).pack(side='right')
    tk.Button(buttons, text="OK", width=10, command=ok).pack(side='right', padx=(0, 8))
    root.bind('<Return>', ok)
    root.bind('<Escape>', cancel)
    root.protocol('WM_DELETE_WINDOW', cancel)
    root.update_idletasks()
    x = (root.winfo_screenwidth() - root.winfo_width()) // 2
    y = (root.winfo_screenheight() - root.winfo_height()) // 3
    root.geometry(f'+{max(x, 0)}+{max(y, 0)}')
    entry.focus_set()
    entry.icursor('end')
    root.mainloop()
    return result['path']


def choose_atlas_cache_location():
    """Where the ABC atlas download lives, in this order: the location saved by
    an earlier choice (if that folder still exists), the default location next
    to this repository (if it exists), or else ask the user once. Asking lists
    every local drive and its free space on the console, then opens
    prompt_atlas_cache_location pre-filled with suggest_atlas_cache_dir. The
    answer is saved to ATLAS_CACHE_LOCATION_FILE. Exits if cancelled."""
    try:
        saved = json.loads(ATLAS_CACHE_LOCATION_FILE.read_text(encoding='utf-8')).get('path')
    except (OSError, ValueError, AttributeError):
        saved = None
    if saved:
        if Path(saved).is_dir():
            usable = usable_atlas_dir(saved)
            if usable is not None:
                if usable != Path(saved):
                    print(f"The ABC atlas download was found inside the saved location, in {usable}.")
                    ATLAS_CACHE_LOCATION_FILE.write_text(json.dumps({'path': str(usable)}, indent=2), encoding='utf-8')
                return usable
            print(f"The saved ABC atlas location {saved} holds other files and no atlas download, so it will not be used.")
        else:
            print(f"The saved ABC atlas location {saved} no longer exists.")
    if DEFAULT_ATLAS_CACHE_DIR.is_dir() and usable_atlas_dir(DEFAULT_ATLAS_CACHE_DIR) is not None:
        return usable_atlas_dir(DEFAULT_ATLAS_CACHE_DIR)

    drives = list_local_drives()
    print(f"No ABC atlas download found at {DEFAULT_ATLAS_CACHE_DIR}.")
    print(f"It can need about {ATLAS_CACHE_SPACE_NEEDED_GB} GB. Local drives:")
    for drive_root, free, total in drives:
        print(f"    {drive_root:6s} {format_gb(free):>10s} free of {format_gb(total)}")
    suggested = suggest_atlas_cache_dir(drives)
    print(f"Suggested location: {suggested}")

    try:
        chosen = prompt_atlas_cache_location(suggested, drives)
    except Exception as e:
        print(f"Could not show the location dialog ({e}).")
        try:
            typed = input(f"ABC atlas location [{suggested}]: ").strip()
        except EOFError:
            typed = ''
        chosen = normalize_chosen_atlas_dir(typed) if typed else suggested
    if chosen is None:
        print("No ABC atlas location chosen; exiting.")
        sys.exit(0)

    try:
        chosen.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"Could not create {chosen} ({e}); exiting.")
        sys.exit(1)
    ATLAS_CACHE_LOCATION_FILE.write_text(json.dumps({'path': str(chosen)}, indent=2), encoding='utf-8')
    print(f"Using {chosen} for the ABC atlas download (saved to {ATLAS_CACHE_LOCATION_FILE}).")
    return chosen


def offer_atlas_release_update(cache):
    """If a newer ABC atlas release exists online than the one this cache last
    used, ask (dialog, or console if no GUI) whether to switch to it. Switching
    is remembered by the library itself (_manifest_last_used.txt in the cache),
    so later runs load the new release without asking. Declining isn't
    remembered: the question comes back every startup until they switch."""
    current = cache.current_manifest
    latest = cache.latest_manifest_file
    if not current or current == latest:
        return

    message = "\n".join([
        "A newer ABC atlas release is available online:",
        f"    in use:  {current}",
        f"    newest:  {latest}",
        "",
        "Switching may re-download any data files that changed between releases (already-downloaded files are kept), "
        "and results saved from the older release may no longer match new runs.",
        "",
        "Switch to the newest release?",
    ])
    print("\n" + message + "\n")
    try:
        import centered_dialog
        switch = centered_dialog.ask_yes_no("Newer ABC atlas release", message)
    except Exception:
        # No GUI available: fall back to the console.
        try:
            switch = input("Switch? [y/N] ").strip().lower() in ('y', 'yes')
        except EOFError:
            switch = False

    if switch:
        print(f"Switching ABC atlas release to {latest}.")
        with warnings.catch_warnings():
            # load_latest_manifest re-warns that files may need re-downloading,
            # which the dialog above already said.
            warnings.simplefilter('ignore', OutdatedManifestWarning)
            cache.load_latest_manifest()
    else:
        print(f"Keeping ABC atlas release {current}.")


import botocore.exceptions
startup_note("Libraries loaded. Looking for the ABC atlas download...")
download_base = choose_atlas_cache_location()
startup_note(f"Atlas download: {download_base}. Checking the internet connection...")
abc_cache = None
offline_reason = "No internet connection was detected."
if _has_internet_connection():
    startup_note("Online: reading the atlas manifest (this also checks for a newer release)...")
    try:
        with warnings.catch_warnings():
            # The library's own "a newer release exists" warning; offer_atlas_release_update asks instead.
            warnings.simplefilter('ignore', OutdatedManifestWarning)
            abc_cache = AbcProjectCache.from_cache_dir(download_base)
    except (OSError, botocore.exceptions.BotoCoreError) as e:      # (the connection failed after all: the local cache is used, as when offline)
        offline_reason = f"The Allen atlas server could not be reached ({type(e).__name__})."
atlas_is_online = abc_cache is not None
if abc_cache is not None:
    offer_atlas_release_update(abc_cache)
else:
    print(f"{offline_reason} Using the local ABC atlas cache only "
          "(skipping the S3 manifest check that would otherwise stall startup for minutes).")
    # NotImplementedError (not a hang) if this session ends up needing a
    # file that was never downloaded — see LocalCache's own _download_
    # file/_download_manifest, both of which just raise instead of trying
    # to reach the network — so a genuinely offline machine with no cache
    # at all still fails fast here, just with a real Python traceback
    # instead of a silent freeze.
    abc_cache = AbcProjectCache.from_local_cache(download_base)
    try:
        import centered_dialog
        centered_dialog.info("Offline", offline_reason + "\n\nThe browser is working offline, with the atlas data already downloaded to this "
                             "computer. Anything not downloaded yet cannot be loaded, and newer atlas releases are not checked for.")
    except Exception:
        pass                                              # (no GUI: the console line above is the notice)


def serialize_downloads(cache):
    """Make a download of any one file wait for another of the same file that is already running. The library writes a download straight
    into its final file, so two at once (e.g. the launch screen's background read of the cell metadata and the browser's own load of it)
    would write over each other, fail the checksum and start again. The second caller now waits, and then finds the file there."""
    locks, guard = {}, threading.Lock()

    def wrap(name):
        original = getattr(cache, name)

        def one_at_a_time(*args, **kwargs):
            names = list(args) + [kwargs[k] for k in ('directory', 'file_name') if k in kwargs]       # (called by position or by keyword)
            with guard:
                lock = locks.setdefault(tuple(names), threading.Lock())
            with lock:
                return original(*args, **kwargs)
        setattr(cache, name, one_at_a_time)
    for method in ('get_data_path', 'get_metadata_path'):
        wrap(method)


serialize_downloads(abc_cache)

# The big files every session needs. The first time on a computer they are downloaded in the background from the moment the program
# starts (the raw expression data alone is ~8 GB: ten minutes or more), so that the wait is over sooner, and anything that needs them
# waits for the download rather than starting another (see serialize_downloads).
ATLAS_CORE_FILES = (("the cell metadata table", 'MERFISH-C57BL6J-638850', 'cell_metadata_with_cluster_annotation', 'metadata'),
                    ("the raw MERFISH expression data", 'MERFISH-C57BL6J-638850', 'C57BL6J-638850/raw', 'data'))
# ...and, first, the reconstructed coordinates and their atlas, when the setting asks for them (it does unless changed): a quarter of a GB, that
# the browser needs as soon as it opens. Asked for separately (prepare_coordinate_frame) they were a second download beside the big one,
# and the two progress bars in the console wrote over each other.
ATLAS_RECON_LABELS = ("the reconstructed cell coordinates", "the atlas regions that go with them")
# ...and the CCF coordinates of the cells (a third of a GB) and the CCF region names: the AP from bregma of the sections, shown in the section grid and the
# single-section readout, is read from them, and so are the CCF atlas overlay and the cell tips' ML and DV. Fetched with the rest, not when first asked for.
ATLAS_CCF_LABELS = ("the CCF cell coordinates", "the CCF region names")
atlas_prefetch = {'thread': None, 'total_gb': 0.0, 'files': [], 'keys': set()}


def atlas_download_running():
    """While the download started by start_atlas_prefetch is still going: how far it has got, as text ("3.1 of 9.1 GB"); else None. (For the
    launch screen's Exit: closing the program stops the download, and a file cut short is downloaded again from the start.)"""
    thread = atlas_prefetch['thread']
    if thread is None or not thread.is_alive():
        return None
    got = 0
    for path, _size in atlas_prefetch['files']:
        try:
            got += path.stat().st_size if path is not None else 0
        except OSError:
            pass
    total = atlas_prefetch['total_gb']
    return f"{got / 1e9:.1f} of {total:.1f} GB" if total else "under way"


def start_atlas_prefetch(cache):
    """Start downloading the core atlas files that are missing, on a background thread (nothing if all are there, or offline).
    Fills atlas_prefetch with the thread and the size, and tells the user once what is going on."""
    core_files = []
    for label, directory, file_name, kind in ATLAS_CORE_FILES:
        path, size, present = metadata_file_status(cache, directory, file_name)
        if present is False:
            core_files.append((label, directory, file_name, kind, size, path))
    ccf = []
    for label, (directory, file_name) in zip(ATLAS_CCF_LABELS, CCF_FILES):
        path, size, present = metadata_file_status(cache, directory, file_name)
        if present is False:
            ccf.append((label, directory, file_name, 'metadata', size, path))
    if not core_files and not ccf:
        return                     # (the reconstructed coordinates alone, a setting changed later, are asked for when needed, as before)
    recon = []
    if core_files and COORD_SETTINGS['reconstructed_coords']:
        for label, (directory, file_name, kind) in zip(ATLAS_RECON_LABELS, RECON_FILES):
            path, size, present = metadata_file_status(cache, directory, file_name)
            if present is False:
                recon.append((label, directory, file_name, kind, size, path))
    missing = recon + ccf + core_files                  # (the small ones first, needed as soon as the browser opens; the big files last)
    for _label, directory, file_name, _kind, size, path in missing:
        atlas_prefetch['files'].append((path, size))
        atlas_prefetch['keys'].add((directory, file_name))
    total_gb = sum(m[4] or 0 for m in missing) / 1e9
    atlas_prefetch['total_gb'] = total_gb

    def work():
        for _label, directory, file_name, kind, _size, _path in missing:          # (one after another: one progress bar in the console)
            try:
                (cache.get_data_path if kind == 'data' else cache.get_metadata_path)(directory=directory, file_name=file_name)
            except Exception as e:
                print(f"Could not download {file_name} ({type(e).__name__}: {e}); it will be tried again when it is needed.")
    atlas_prefetch['thread'] = threading.Thread(target=work, daemon=True)
    atlas_prefetch['thread'].start()

    def size_text(n):
        return f"{n / 1e9:.1f} GB" if n >= 1e8 else f"{max(n / 1e6, 1):.0f} MB"
    parts = "\n".join(f"    {m[0]}: {size_text(m[4])}" if m[4] else f"    {m[0]}" for m in missing)
    if total_gb >= 1:
        message = (f"The atlas data is not on this computer yet, so it is being downloaded now (one time only):\n\n{parts}\n\n"
                   f"That is about {total_gb:.0f} GB and can take 10 minutes or much longer, depending on the connection. Progress is shown in "
                   "the console window. You can carry on to the launch screen; opening the browser will wait for the download to finish. "
                   "Keep this program open until then: an interrupted file starts again from the beginning.")
    else:                                                       # (only the small files, such as the CCF coordinates, were missing)
        message = (f"A little atlas data is not on this computer yet, so it is being downloaded now (one time only):\n\n{parts}\n\n"
                   f"That is about {total_gb * 1000:.0f} MB. Progress is shown in the console window; you can carry on.")
    print("\n" + message + "\n")
    try:
        import centered_dialog                                  # (in the middle of the screen: a native box over a hidden window was partly off it on macOS)
        centered_dialog.info("Downloading the atlas", message)
    except Exception:
        pass                                              # (no GUI: the console text above is the notice)


if atlas_is_online:
    try:
        start_atlas_prefetch(abc_cache)
    except Exception as e:
        print(f"Could not check whether the atlas data needs downloading ({e}).")
startup_note("Atlas manifest ready.")

CELL_TYPE_OPTIONS = {'1': 'Neurons', '2': 'NonNeurons', '3': 'All'}
CELL_TYPE_SUFFIXES = {'Neurons': 'neurons', 'NonNeurons': 'nonneurons', 'All': 'all-cells'}      # (the words in a run folder's name: umap_all-cells_sections-...)


# obs column holding scanpy's own Leiden clustering (Step 6), computed from
# this run's neighbor graph. Deliberately *not* called 'cluster': that name
# is already taken by the Allen ABC taxonomy's finest annotation level (see
# metadata_cols), which is a published label, not something derived here.
LEIDEN_KEY = 'leiden'
# Passed straight through to sc.tl.leiden's own `resolution` argument —
# scanpy's own default (1.0) if left at 1.0 here. Roughly, cluster *count*
# scales with resolution: higher values favor more, smaller clusters; lower
# values favor fewer, larger ones. There's no formula for "resolution X
# gives Y clusters" (it depends on the actual neighbor graph), so tuning
# this is trial and error — rerun Step 6 (or offer_leiden_backfill for a
# cached run) after changing it and check adata.obs[LEIDEN_KEY].nunique().
# Typical exploratory range is roughly 0.3 (coarse) to 2.0 (fine); scanpy's
# docs suggest 0.4-1.2 for a "usual" range, but that's dataset-dependent.
LEIDEN_RESOLUTION = 1.2

# ===========================================================================
# TUNABLE CONSTANTS — INTERACTIVE UMAP VIEWER
# ===========================================================================
# Everything here belongs to the "Interactive UMAP viewer" window opened at
# Step 9 (show_interactive_umap_window) — the one with the brain-section
# panels on the left, the UMAP scatter in the middle, and the sidebar
# controls on the right. Edit freely; all are read live, so a change takes
# effect on the next launch with no other edits needed.
#
# NOT collected here, deliberately:
#   * The section/ROI picker windows' own constants, which live in
#     prompt_section_selection_gui / prompt_subregion_selection. Two names
#     appear in both places with *different* values, so the viewer's copies
#     are prefixed VIEWER_ below to keep them distinct.
#   * This window's layout fractions (GRID_LEFT, GAP, CBAR_WIDTH,
#     QUERY_TOP_Y, the STATUS_PANEL_* family, ...) — they only make sense
#     next to the layout code that consumes them.
#   * Sizes derived from the ones here (UMAP_HIGHLIGHT_BASE_SIZE and
#     friends) and anything derived from the screen-dependent font size;
#     both stay put, and pick up changes made here automatically.
#   * AREA_BOTTOM / GRID_RIGHT / UMAP_LEFT / UMAP_RIGHT / SIDEBAR_LEFT /
#     SIDEBAR_WIDTH / HANDLE_BOTTOM: despite the naming these are mutable
#     layout *state*, reassigned when the resize handles are dragged.

# --- UMAP scatter appearance ---
UMAP_POINT_SIZE = 10       # matplotlib `s` (an AREA) for each cell, at full zoom-out
UMAP_LABEL_FONTSIZE = 10   # centroid ID labels, at full zoom-out; grows with the dots
UMAP_GROUP_DIAMETER_RATIO = 1.25  # family-highlight '+' size vs the cell's own CURRENT diameter

# When True, categorical modes ('All'/'Single Subclass' etc.) use
# category_rank_shape's marker shapes on the *live* UMAP too, not just
# the saved PNG/SVG (which always uses shapes — see render_export_figure).
# Costs one scatter() call per (shape, filled/open) combination actually
# in use instead of one call total, but only at redraw time (a mode/level/
# query change) — the far more frequent zoom/pan ticks already redraw
# from a cached bitmap regardless of how many artists made up the last
# real one, so the expected performance cost of this is minor.
UMAP_USE_SHAPES_ON_SCREEN = True

# On-screen ID legend (draw_id_legend, in cbar_ax) for 'All <level>s' mode,
# at whichever level is currently selected (class/subclass/supertype/
# cluster/leiden). How many rows fit is computed dynamically from the
# window's own current height (see legend_max_rows_for_current_size,
# inside show_interactive_umap_window) rather than a fixed count here —
# a 4k monitor can comfortably show far more rows than a laptop screen.
# Rows beyond that count are dropped smallest-cell-count first (the
# ranking compute_ranked_category_colors already produces), not shrunk to
# fit — past a certain row count the text stops being legible regardless
# of how small it's squeezed; the export-only legend (build_export_
# legend_entries) is what covers the rest for a level with far more
# categories than could ever fit here (subclass/supertype/cluster
# routinely do) — see redraw_all_subclasses' own comment.

# Per-cell opacity. Below 1.0, overlapping cells accumulate into visibly
# darker regions, so a dense cluster reads as dense rather than as one flat
# patch of color — the trade is that an isolated cell is fainter. 1.0
# disables the effect. Applies to the cells only: the hover ring and the
# family-highlight '+' markers stay fully opaque, so they still read clearly
# against whatever is underneath them.
UMAP_POINT_ALPHA = 0.75

# Baseline (every gene at its low end) color for the multi-gene (2 or 3
# genes, red/green/blue) overlay — see multi_gene_rgb — shared by the UMAP
# and the section panels, so a cell's color means the same thing (0 =
# black, rising toward red/green/blue as expression increases) regardless
# of which plot it's read off. Both plots' own *axes* background is set
# separately, to MULTI_GENE_UMAP_FACECOLOR (a dark grey, not this pure
# black) — a zero-expression cell (this color) needs to stay visibly a
# cell, not melt into the space around it the way it would against an
# identically-black background. Applied in redraw_multi_genes directly for
# the UMAP's own `ax`, and via set_section_panel_facecolor for every
# section panel (whose normal, every-other-mode background is pure black —
# SECTION_PANEL_FACECOLOR — restored by clear_colorbar on leaving
# multi-gene mode).
MULTI_GENE_LOW_EXPRESSION_COLOR = (0.0, 0.0, 0.0)
MULTI_GENE_UMAP_FACECOLOR = (0.1, 0.1, 0.1)
# Genes 4-6 (UMAP only, so far): a second red/green/blue overlay, drawn as
# small filled circles centered on top of the circle layer above, instead of
# its own separate plot — see redraw_multi_genes. Opaque (alpha=1). This is
# specifically the 'circle' shape's own diameter multiplier (see GENE456_
# SHAPE_REGISTRY below, where each shape carries its own) — kept smaller
# than the base circle so the base layer's own color still shows around
# its edges. Diameter (not area — matplotlib `s` is area, so this gets
# squared before use) multiplier on the circle layer's own current size.
MULTI_GENE_PLUS_SIZE_DIAMETER_MULTIPLIER = 0.5
# Cells with exactly zero expression across all of genes 4-6 simply don't
# get an overlay dot drawn at all (rather than an opaque black one) — see
# redraw_multi_genes' own visibility mask.
MULTI_GENE_PLUS_ZORDER = 1.5  # above the circle layer's default zorder (1), below hover highlights (5.2+)
# Selectable marker shapes for the overlay layer — user-facing names to
# a (matplotlib marker, own diameter multiplier) pair. The marker half is
# whatever scatter(marker=...) actually wants: 'circle' is a plain built-
# in marker code; 'half-circle' is a custom matplotlib Path (scatter's
# marker= argument accepts an arbitrary Path, not just single-character
# codes) — Path.wedge(-90, 90) is a wedge spanning that range, which is
# exactly a semicircle: flat (diameter) edge on the left, curved edge on
# the right. The size half is per-shape rather than one shared constant:
# 'circle' stays smaller than the base layer's own circle (so that base
# color still shows around its edges — the original design intent), while
# 'half-circle' is drawn at the *same* radius as the base circle, so its
# curved edge lines up exactly with the base circle's own — covering its
# right half outright rather than sitting as a smaller shape somewhere
# inside it. A new shape just needs a new entry here; nothing else needs
# to know how it's actually built or sized. Picked via the Settings dialog
# (prompt_viewer_settings_dialog) and stored in viewer_settings/
# SETTINGS_FILE as 'gene456_shape' — see load_viewer_settings's own
# docstring for that setting's persistence rule.
#
# Path.wedge(-90, 90)'s own raw vertices span radius 1.0 from the origin —
# *not* the same scale matplotlib's built-in markers use. 'o' (any single-
# character marker code) gets normalized by MarkerStyle to a path whose
# own bounding box is [-0.5, 0.5] on each axis, radius 0.5 — verified
# directly (MarkerStyle('o').get_path().transformed(...).get_extents())
# rather than assumed, since scatter's `s` scales a marker by its own raw
# path coordinates, not by some shape-independent notion of "radius".
# Passing the wedge unscaled at the same size_multiplier as 'circle' drew
# it at roughly double the circle's own radius instead of matching it.
# Pre-scaling the vertices by 0.5 here (once, at import time — not
# per-draw) brings it onto the same radius-0.5 scale 'o' already uses, so
# a size_multiplier of 1.0 really does mean "same radius as the base
# circle" rather than needing its own, separately-tuned fudge factor.
_HALF_CIRCLE_RAW_PATH = MplPath.wedge(-90, 90)
GENE456_SHAPE_REGISTRY = {
    'small circle': {'marker': 'o', 'size_multiplier': MULTI_GENE_PLUS_SIZE_DIAMETER_MULTIPLIER},
    'half circle': {
        'marker': MplPath(_HALF_CIRCLE_RAW_PATH.vertices * 0.5, _HALF_CIRCLE_RAW_PATH.codes),
        'size_multiplier': 1.0,
    },
}
GENE456_SHAPE_DEFAULT = 'small circle'

# Selectable overall UMAP dot diameter (a preference, independent of zoom —
# see update_umap_dot_size, which multiplies this together with the
# current zoom multiplier for the main scatter, the hover ring, and the
# group highlight, but deliberately *not* the centroid ID labels' own
# font size, which stays tied to zoom alone). Handy for a sparse dataset
# (bump it up so individual cells are easier to see/click) or a very dense
# one (shrink it so points overlap less).
UMAP_CIRCLE_SIZE_REGISTRY = {
    'Extra small': 0.5,
    'Small': 0.7,
    'Default': 1.0,
    'Large': 1.4,
    'Extra large': 2.0,
}
UMAP_CIRCLE_SIZE_DEFAULT = 'Default'
# The UMAP axes' own background for every *other* mode (categorical
# taxonomy levels, single gene) — restored explicitly by clear_colorbar
# whenever leaving multi-gene mode. Needed because Axes.clear() does *not*
# reset facecolor on its own (confirmed directly — it's one of the few
# Artist properties clear() leaves alone), so without this, a dark
# MULTI_GENE_UMAP_FACECOLOR set while looking at two or three genes stayed
# on screen indefinitely after switching to 'All <level>s'/'Specified
# <level>(s)', where there's no gradient-of-brightness cell color to need
# a dark background for in the first place.
UMAP_AXES_FACECOLOR = 'white'

# How fast dots (and the centroid labels) grow as you zoom in.
#   0.0 = constant on-screen size at every zoom level
#   1.0 = grows in step with the zoom, so a cluster looks identical at any zoom
# Applied to the DIAMETER: at zoom ratio Z, diameter is 1 + RATE*(Z-1).
# matplotlib's `s` is an area, so the dots square this; the label font, being
# a linear dimension, takes it as-is — which is what keeps the two in
# constant proportion.
UMAP_ZOOM_DOT_GROWTH_RATE = 0.75

# How many of a level's categories — ranked by cell count, largest first —
# get a real, distinct color (and a centroid ID label) in 'All <level>s'
# mode. Below this cut, categories are gray with no label on the UMAP
# itself. class/subclass rarely exceed this; supertype/cluster/leiden
# routinely have hundreds to thousands of categories, and matplotlib has no
# way to make that many colors mutually distinguishable regardless of what
# palette generates them — this keeps the ones with the most cells (and
# therefore the most visual weight) identifiable, rather than making every
# category equally hard to tell apart.
#
# The section panels use the *same* ranking but don't gray out the excess —
# they cycle back through the same N colors instead (see
# set_section_colors_categorical), since a flat gray region there reads as
# "nothing here" rather than "many small categories here". A recycled color
# no longer uniquely identifies one category, but hovering any cell still
# reports its exact one regardless of color.
UMAP_MAX_COLORED_CATEGORIES = 400  # category_rank_color/category_rank_shape's own combined cycle length

# Experimental: when True, each centroid ID label's font size scales with
# the 1/4 POWER of that category's own cell count relative to the median
# among labeled categories (see cell_weighted_median_count), instead of
# every label sharing one flat UMAP_LABEL_FONTSIZE — a category at the
# median gets scale 1 (unchanged), one with 10x the median gets scale
# 10**0.25 ≈ 1.78, one with 100x gets ≈3.16 (clamped down to MAX_SCALE
# below), one with 1/10th gets ≈0.56, one with 1/100th gets ≈0.32 (clamped
# up to MIN_SCALE below). Gentler than a straight log10 (which would put
# 10x at scale 2 and 1/10th at scale 0 exactly) — cell counts across
# categories routinely span several orders of magnitude, and the 1/4
# power keeps the largest few categories' labels from dwarfing everything
# else while still keeping tiny categories legible without hitting the
# floor immediately.
# Clamped by the two _SCALE constants below so a very rare or very dominant
# category can't shrink to unreadable or balloon over its neighbors, and
# further adjusted by zoom level — see umap_label_fontsize, the only place
# that clamping actually happens; MIN/MAX_SCALE are just its inputs. Toggle
# UMAP_LABEL_SIZE_BY_CELL_COUNT off if the whole effect reads as more
# distracting than informative.
UMAP_LABEL_SIZE_BY_CELL_COUNT = True
UMAP_LABEL_SIZE_BY_CELL_COUNT_MIN_SCALE = 0.5
UMAP_LABEL_SIZE_BY_CELL_COUNT_MAX_SCALE = 2.5

# --- Brain-section panels (left) ---
SECTION_BACKGROUND_BASE_SIZE = 3     # `s` for each background cell, at full zoom-out
SECTION_POINT_ALPHA = 0.75           # as UMAP_POINT_ALPHA, for the panels' cells
SECTION_ZOOM_DOT_GROWTH_RATE = 0.5   # as UMAP_ZOOM_DOT_GROWTH_RATE, for the panels' dots
SECTION_UNKNOWN_COLOR = 'dimgray'    # cells with no value at the current level
# Panel axes background for every mode except multi-gene (which swaps it to
# MULTI_GENE_UMAP_FACECOLOR — see that constant's own comment; a genuinely
# black background there would make a zero-expression cell, itself black,
# disappear into it). set at panel-build time and restored by clear_
# colorbar whenever leaving multi-gene mode — see set_section_panel_facecolor.
SECTION_PANEL_FACECOLOR = 'black'
VIEWER_SPAN_PERCENTILE = 90          # percentile of section extents setting the shared panel scale
PANEL_PROGRESS_INTERVAL = 10         # log a progress line every N panels while building the grid

# --- Section-panel home-view disk cache (interactive viewer only) ---
# One rendered PNG per (run, section, level) — the "All <level>s" mode's
# categorical background coloring at each panel's fully-zoomed-out (home)
# extent, persisted into the run's own folder (see show_interactive_umap_
# window's own run_folder) so re-showing it — even in a brand new session,
# on a different machine sharing that folder — is instant instead of a
# real ~80-panel vector re-render (which is what was taking ~30s).
SECTION_HOME_CACHE_DPI = 150
# Rendered at whatever size *this* panel's own current on-screen box
# actually needs (see show_section_home_cache_or_scatter), not a fixed
# constant — a run with few sections (each panel then gets a much bigger
# box than one of 60+ sections sharing the same grid area) needs a much
# larger cache than a many-sections run to look sharp, and a fixed size
# tuned for one case was either blurry (too small) or needlessly slow to
# render/load/zoom (too big) for the other.
#
# The size safety margin used when computing that target (see show_
# section_home_cache_or_scatter) is deliberately ZOOM_BITMAP_ONLY_MAX_
# MULTIPLIER itself (defined below), not an independently-tuned constant:
# the cached bitmap is what the zoom-preview system shows, magnified, for
# every ratio between 1.0x (home) and that threshold before switching to a
# real per-cell redraw — so it already has to survive being shown at up to
# that multiplier's own on-screen size just to do its existing job,
# regardless of any cross-machine concern. Reusing that same constant
# (rather than a separate one that happens to start at the same value)
# keeps that relationship explicit and keeps them in sync automatically if
# the zoom threshold ever changes.
# Not screen-size independent by construction any more (the cache used to
# be rendered at one fixed size regardless of what asked for it) — instead
# self-healing: show_section_home_cache_or_scatter compares an existing
# cache's own pixel dimensions against what the *current* panel/window
# actually needs on every load (a cheap PNG header read, not a real
# decode) and regenerates it if it falls short. A smaller window/screen
# than whatever generated the cache gets a fast hit (the cache already
# exceeds its need); one needing more resolution pays the render cost
# once, and every later load — on any machine sharing the run folder —
# benefits from the now-larger cache.
# Bump this whenever render_section_home_view_png's own output would
# change (colors, dot size/style, point selection, ...) in a way that
# makes an already-cached PNG show something subtly wrong — there's no
# way to detect that automatically, so a stale cache would otherwise just
# keep being served as if still valid. Bumped 2 -> 3 for the dot_size_
# scale fix (cell dots baked into the cache were rendering undersized by
# 1/ZOOM_BITMAP_ONLY_MAX_MULTIPLIER^2, making the live, unaffected hover
# family-highlight '+' markers look disproportionately bigger next to
# them) — the resolution-based hit-check alone wouldn't have caught this,
# since it only compares pixel dimensions, not what's drawn inside them.
SECTION_HOME_CACHE_VERSION = 3

# --- Section-panel scale bar (first panel only; see build_section_scalebar) ---
SCALEBAR_TARGET_FRACTION = 0.22      # of the panel's current view width
SCALEBAR_MARGIN_FRACTION = 0.06      # inset from the panel's own edges, as a fraction of its view
SCALEBAR_TEXT_GAP_FRACTION = 0.02    # extra gap above the bar, for its length label
# A 1-2-5 sequence in micrometers, the ABC atlas's own native unit divided by
# 1000 (see HOVER_CELL_RADIUS_DATA_UNITS's comment on why the atlas's 'x'/'y'
# columns are actually millimeters). Covers a single MERFISH cell's own
# width (~10 um) up to several whole sections (~50 mm) side by side.
SCALEBAR_NICE_LENGTHS_UM = (
    1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000,
)
SECTION_MAPS_SAVE_DPI = 200           # matches save_roi_map/save_group_spatial_maps's own PNGs

# --- Interaction timing (milliseconds unless noted) ---
# How long a scroll burst must go quiet before the cheap stand-in bitmaps are
# swapped back for a real, full-quality redraw. Longer = fewer expensive
# redraws while scrolling, at the cost of the sharp view arriving later.
ZOOM_PREVIEW_SETTLE_MS_UMAP = 800
ZOOM_PREVIEW_SETTLE_MS_PANEL = 800
# Below this zoom ratio (home_span / current_span — 1.0 at fully zoomed
# out, 3.0 means the current view spans a third of the home extent), the
# UMAP's own zoom/pan settle (end_zoom_previews) skips its real, full-
# resolution fig.canvas.draw() entirely and just stays on the cheap
# zoom-preview bitmap already on screen — a magnified crop, capped at
# whatever resolution it was captured at, rather than a fresh re-render.
# At or above it, the settle still does a real draw, but first filters the
# UMAP scatter down to only the points inside the current view (see
# filter_main_scatter_to_viewport) — every real draw used to reprocess the
# *entire* ~200k-cell point set through matplotlib's transform/clip
# pipeline regardless of zoom level, which is almost certainly why zooming
# stayed slow even with very little actually visible on screen. Also
# governs the section panels' own zoom/pan settle the same way, against
# section_zoom_ratio (the one shared zoom level all ~80 panels zoom
# together at) instead — see filter_all_section_scatters_to_viewport.
ZOOM_BITMAP_ONLY_MAX_MULTIPLIER = 2.0
VIEWER_HOVER_HOLD_MS = 250       # cursor must settle this long before the hovered cell is looked up
GROUP_HOVER_HOLD_MS = 500        # ...and this long before its whole family is highlighted
LAYOUT_REDRAW_MIN_INTERVAL = 0.03  # SECONDS; caps redraws while dragging the resize handles

# --- Saving ("Save UMAP" button) ---
# Resolution of the saved PNG. The SVG is vector and ignores this — but note
# an SVG of a large UMAP stores every cell as its own element, so for a
# several-hundred-thousand-cell run that file is big and slow to open; the
# PNG is the practical one at those sizes.
UMAP_SAVE_DPI = 450

# --- Rendering detail ---
# How far inside its own boundary a zoom/pan preview snapshot is cropped, to
# keep the axes' antialiased spine out of the bitmap — without this the
# baked-in border draws as a hard rectangle over the layer beneath it.
SPINE_INSET_PX = 2

# --- Status-bar text ---
HOVER_DEFAULT_MESSAGE = 'Hover over any cell to see its class/subclass/supertype/cluster.'
# Between the fields of the hover status line. A semicolon (not a comma)
# because several of these values contain commas of their own.
HOVER_FIELD_SEP = ';   '
# ===========================================================================



# Persisted, most-recent-first list of output folders the user has run
# against, surfaced as a dropdown on the startup panel (prompt_startup_panel)
# so switching between projects is one click rather than a re-typed path.
RECENT_OUTPUT_FOLDERS_PATH = CACHE_DIR / 'recent_output_folders.json'
MAX_RECENT_OUTPUT_FOLDERS = 10


def _normalized_folder_key(folder):
    """Case- and separator-normalized absolute form of `folder`, so
    'M:/proj', 'M:\\proj\\' and 'm:\\PROJ' all de-duplicate to one entry in
    the recent-folders list on Windows. Falls back to the raw string if the
    path can't be normalized (e.g. an unmounted drive)."""
    try:
        return os.path.normcase(os.path.abspath(os.path.expanduser(str(folder))))
    except Exception:
        return str(folder)


def load_recent_output_folders():
    """The recently-used output folders, most-recent-first, at most
    MAX_RECENT_OUTPUT_FOLDERS and de-duplicated. Empty list on first run or
    any read/parse failure — this list is a convenience, never load-bearing,
    so a broken file just means "no history yet"."""
    try:
        with open(RECENT_OUTPUT_FOLDERS_PATH, encoding='utf-8') as f:
            data = json.load(f)
    except Exception:
        return []
    if not isinstance(data, list):
        return []
    out, seen = [], set()
    for entry in data:
        s = str(entry)
        key = _normalized_folder_key(s)
        if key not in seen:
            seen.add(key)
            out.append(s)
    return out[:MAX_RECENT_OUTPUT_FOLDERS]


def remember_recent_output_folder(folder):
    """Move `folder` to the front of the persisted recent-folders list
    (de-duplicated via _normalized_folder_key, capped at
    MAX_RECENT_OUTPUT_FOLDERS) and write it back. Best-effort: a write
    failure is swallowed. Returns the updated list."""
    folder = str(folder)
    key = _normalized_folder_key(folder)
    updated = [folder] + [p for p in load_recent_output_folders()
                          if _normalized_folder_key(p) != key]
    updated = updated[:MAX_RECENT_OUTPUT_FOLDERS]
    try:
        with open(RECENT_OUTPUT_FOLDERS_PATH, 'w', encoding='utf-8') as f:
            json.dump(updated, f, indent=2)
    except Exception:
        pass
    return updated


def prompt_cell_type_selection():
    print("Select cell type subset to compute UMAP for:")
    for key, label in CELL_TYPE_OPTIONS.items():
        print(f"  {key}. {label}")
    while True:
        choice = input("Enter 1, 2, or 3: ").strip()
        if choice in CELL_TYPE_OPTIONS:
            return CELL_TYPE_OPTIONS[choice]
        print("Invalid selection; please enter 1, 2, or 3.")


def materialize_subset(adata, mask):
    """Subset adata by a boolean mask and bring the result into memory.
    Backed AnnData objects can't be .copy()'d directly (anndata raises,
    asking for a filename); .to_memory() on the subset view instead reads
    only the kept rows from disk, which is what keeps memory use down when
    filtering runs right after a backed load."""
    subset = adata[mask]
    return subset.to_memory() if adata.isbacked else subset.copy()


def filter_by_cell_type(adata, selection, class_col='class'):
    """Subset adata to the selected cell type, based on the leading numeric
    ID in the 'class' column (non-neuron classes are 30-34)."""
    if selection == 'All':
        return adata
    if class_col not in adata.obs.columns:
        print(f"Warning: '{class_col}' column not found; cannot filter by cell type. "
              "Proceeding with all cells.")
        return adata

    class_ids = adata.obs[class_col].astype(str).str.strip().str.extract(r'^(\d+)')[0]
    class_ids = pd.to_numeric(class_ids, errors='coerce')
    valid = class_ids.notna()
    is_non_neuron = valid & class_ids.isin(NON_NEURON_CLASS_IDS)

    if selection == 'NonNeurons':
        keep = is_non_neuron
    else:  # 'Neurons'
        keep = valid & ~is_non_neuron

    filtered = materialize_subset(adata, keep.to_numpy())
    print(f"Filtered to {selection}: {filtered.n_obs} of {adata.n_obs} cells kept.")
    return filtered










def parse_index_selection(raw):
    """Parse a comma-separated list of indices and/or dash-ranges, e.g.
    '28-33, 35' -> [28, 29, 30, 31, 32, 33, 35]. Raises ValueError on any
    malformed token. Preserves first-seen order and drops duplicates."""
    indices = []
    for token in raw.split(','):
        token = token.strip()
        if not token:
            continue
        if '-' in token:
            start_str, sep, end_str = token.partition('-')
            start, end = int(start_str.strip()), int(end_str.strip())
            if start > end:
                start, end = end, start
            indices.extend(range(start, end + 1))
        else:
            indices.append(int(token))
    return list(dict.fromkeys(indices))


def prompt_section_selection(section_series):
    """Console-based section picker. Returns (selected_labels, sections_suffix,
    rois); selected_labels is None when all sections are selected. rois is
    always [] here — ROI drawing needs the GUI, so it's unavailable in the
    console fallback.

    A future GUI can bypass this prompt entirely and call filter_by_sections()
    directly with a list of section labels chosen by clicking.
    """
    if section_series is None:
        print("No section metadata available; proceeding with all sections.")
        return None, 'allsections', []

    counts = section_series.value_counts()
    unique_sections = sorted_sections_descending(counts.index.tolist())
    print(f"Available brain sections ({len(unique_sections)}):")
    for i, s in enumerate(unique_sections, start=1):
        print(f"  {i}. {s} ({counts[s]} cells)")
    print("Enter section number(s) to analyze — comma-separated, dash-ranges allowed "
          "(e.g. 28-33, 35), or 'all':")
    while True:
        raw = input("> ").strip()
        if raw.lower() == 'all':
            return None, 'allsections', []
        try:
            idxs = parse_index_selection(raw)
        except ValueError:
            idxs = []
        if idxs and all(1 <= i <= len(unique_sections) for i in idxs):
            selected = [unique_sections[i - 1] for i in idxs]
            tokens = '-'.join(sanitize_section_token(s) for s in selected)
            return selected, f'sections-{tokens}', []
        print(f"Invalid input; enter comma-separated numbers and/or dash-ranges between "
              f"1 and {len(unique_sections)}, or 'all'.")










CCF_CACHE_VERSION = 1
THREED_MAX_POINTS = 60000     # points drawn across all groups in the in-app window; more are randomly thinned
THREED_WEB_MAX_POINTS = 300000  # the browser (WebGL) view copes with far more
THREED_SHELL_POINTS = 60000   # whole-brain sample: its silhouette is the outline, a subset is the faint cloud
THREED_CLOUD_POINTS = 25000
THREED_BACKGROUND = '#12141a'
# The browser view draws each cell as a small lit sphere (a real mesh, so overlapping cells shade and occlude one
# another and read as 3D) rather than as a flat marker sprite, which merges into a flat patch where cells overlap.
# Meshes cost file size and GPU time per cell, so the sphere detail steps down with the number of cells drawn and
# gives way to outlined flat markers past THREED_WEB_SPHERE_MAX_POINTS.
THREED_WEB_SPHERE_RADIUS_FRACTION = 0.003   # sphere radius as a fraction of the plotted data's longest span
THREED_WEB_SMOOTH_SPHERE_MAX_POINTS = 6000  # up to this many cells: 42-vertex spheres (round when zoomed in)
THREED_WEB_SPHERE_MAX_POINTS = 40000        # up to this many: 12-vertex spheres; beyond, outlined flat markers
# Edges for the browser view: the Allen CCF annotation volume (every 25 um voxel labelled with its brain region, the same
# parcellation_index the cells carry) gives the true brain surface, the ventricles and any region's boundary.
CCF_VOLUME_DIRECTORY = 'Allen-CCF-2020'
CCF_VOLUME_FILE_NAME = 'annotation_25'      # ~3.6 MB download; 528 x 320 x 456 voxels, axes (AP, DV, ML) as x/y/z_ccf
CCF_VOLUME_VOXEL_MM = 0.025
CCF_VENTRICLE_DIVISIONS = ('VL', 'V3', 'V4', 'AQ', 'c')   # lateral, third, fourth ventricle, aqueduct, central canal
CCF_EDGES_CACHE_VERSION = 1
# White-matter background objects, each a set of CCF fibre-tract labels picked by the parcellation table's `structure`
# or `substructure` acronyms: (column, acronyms, colour). The corpus callosum is every callosal part (forceps, genu,
# body, splenium, external capsule); the commissures are the anterior (olfactory and temporal limbs), posterior,
# hippocampal, cerebellar, superior and inferior collicular commissures and the trapezoid body; the decussations are the
# supraoptic, dorsal and ventral tegmental and pyramidal ones.
CCF_WHITE_MATTER = {
    'corpus callosum': ('structure', ('cc',), '#ffffff'),
    'commissures': ('substructure', ('act', 'aco', 'pc', 'dhc', 'vhc', 'cbc', 'csc', 'cic', 'tb'), '#ffd9a0'),
    'decussations': ('substructure', ('sup', 'dtd', 'vtd', 'pyd'), '#b8f2e6'),
}
THREED_EDGE_SLICE_MM = 0.5            # contour spacing (every direction) for the brain and ventricle line edges
THREED_EDGE_REGION_SLICE_MM = 0.25    # ... and for a cells' region, which is usually much smaller
THREED_EDGE_REGION_MAX = 4            # at most this many regions (the ones holding most of the shown cells) get edges
THREED_EDGE_REGION_MIN_SHARE = 0.05   # ... and only those holding at least this share of them
THREED_EDGE_BRAIN_COLOR = '#c9d1e0'
THREED_EDGE_VENTRICLE_COLOR = '#5aa9ff'
THREED_EDGE_REGION_COLORS = ('#ffd166', '#06d6a0', '#ef8354', '#bb86fc')

# Plotly paints translucent geometry in data order, not by depth: with the opacity slider below 100% a far sphere can
# land on top of a near one (and a sphere's far side shows through its near side), depending on the viewing direction.
# All the cells therefore share ONE mesh (so one draw order covers every group), and this page script keeps that
# mesh's triangles sorted back-to-front along the camera's view direction (painter's algorithm) whenever the camera or
# the opacity changes. It also owns group visibility: a legend click hides a group by dropping its triangles, since
# the groups no longer have a trace each (the legend entries are empty stand-in traces). The sort is by view
# direction, not true perspective distance, which is plenty for ordering overlapping cells.
# Several scripts of the page want to answer clicks on the legend (the cells' groups, the envelopes). Plotly takes only ONE
# decision from all the legend-click handlers (the last one's return value, in practice), so they cannot each register
# their own: this router, run first, is the only handler, and asks each script's route in turn; a route returns false if
# it dealt with the click (plotly's own toggle is then skipped), true if the click is not its business. A double click
# (which would isolate one entry, bypassing the scripts' own bookkeeping) is turned off whenever any route exists.
THREED_WEB_LEGEND_ROUTER_SCRIPT = """
(function () {
  var gd = document.getElementById('{plot_id}');
  // The busy ("hourglass") cursor: gd.__busyStart() is called the moment a control is used, so the person sees at once that
  // the click registered, and the cursor goes back to normal when the figure has been quiet for a moment (plotly's events
  // after each redraw keep it going while the work continues). The first wait covers plotly's own short delay before it acts
  // on a legend click.
  var style = document.createElement('style');
  style.textContent = 'body.three-d-busy, body.three-d-busy * { cursor: wait !important; }';
  document.head.appendChild(style);
  var quiet = 0, longest = 0;
  function busyEnd() { clearTimeout(quiet); clearTimeout(longest); document.body.classList.remove('three-d-busy'); }
  function keepBusy(delay) {
    if (!document.body.classList.contains('three-d-busy')) return;
    clearTimeout(quiet);
    quiet = setTimeout(busyEnd, delay);
  }
  gd.__busyStart = function () {
    document.body.classList.add('three-d-busy');
    keepBusy(1000);
    clearTimeout(longest);
    longest = setTimeout(busyEnd, 15000);   // never stuck on, whatever happens
  };
  ['plotly_afterplot', 'plotly_restyle', 'plotly_relayout'].forEach(function (name) {
    gd.on(name, function () { keepBusy(450); });
  });
  gd.__legendRoutes = [];
  gd.on('plotly_legendclick', function (e) {
    var carryOn = true;
    gd.__legendRoutes.forEach(function (route) { if (route(e) === false) carryOn = false; });
    return carryOn;
  });
  gd.on('plotly_legenddoubleclick', function () { return gd.__legendRoutes.length ? false : true; });
})();
"""


THREED_WEB_MESH_SCRIPT = """
(function () {
  var gd = document.getElementById('{plot_id}');
  var meshIndex = __MESH_INDEX__, legendIndices = __LEGEND_INDICES__, vertexStarts = __VERTEX_STARTS__;
  // _fullData, not data: plotly keeps the packed (base64) input in data and the decoded typed arrays in _fullData.
  var full = gd._fullData[meshIndex];
  var o = {x: full.x, y: full.y, z: full.z, i: Int32Array.from(full.i), j: Int32Array.from(full.j),
           k: Int32Array.from(full.k)};
  var vertexCount = o.x.length, faceCount = o.i.length;
  var vertexGroup = new Uint8Array(vertexCount);
  vertexStarts.forEach(function (start, group) { for (var v = start; v < vertexCount; v++) vertexGroup[v] = group; });
  var hidden = legendIndices.map(function () { return false; });
  var busy = false, timer = null, forced = false, lastDirection = null;
  // A restyle makes plotly re-plot the 3D scene from the camera stored in the layout, which lags behind the live one while a
  // scroll-zoom or drag is still settling, so the view would snap back a little. Handing the live camera to the layout
  // first avoids that.
  function syncCamera() {
    try {
      var scene = gd._fullLayout.scene, live = scene._scene.getCamera();
      scene.camera = live;
      if (gd.layout.scene) gd.layout.scene.camera = live;
    } catch (err) { /* the layout keeps its own camera */ }
  }
  // Per-face data, computed once: each face's centre (what the Slice cut tests) and its group (what the legend hides).
  var faceCx = new Float32Array(faceCount), faceCy = new Float32Array(faceCount), faceCz = new Float32Array(faceCount);
  var faceGroup = new Uint8Array(faceCount);
  for (var face = 0; face < faceCount; face++) {
    var fa = o.i[face], fb = o.j[face], fc = o.k[face];
    faceCx[face] = (o.x[fa] + o.x[fb] + o.x[fc]) / 3;
    faceCy[face] = (o.y[fa] + o.y[fb] + o.y[fc]) / 3;
    faceCz[face] = (o.z[fa] + o.z[fb] + o.z[fc]) / 3;
    faceGroup[face] = vertexGroup[fa];
  }
  // The back-to-front order of ALL the faces for one viewing direction. Hiding a group or moving the Slice cut only
  // removes faces from this order, never reorders the rest, so it is computed once per direction and filtered after that
  // (a slider drag then costs a pass over the faces, not a sort). The sort itself packs each face's depth and index into
  // one number so the browser's native numeric sort can do it, far faster than a comparator function.
  var sortedAll = null, sortedDirection = null;
  function sortAll(ex, ey, ez) {
    var depth = new Float64Array(vertexCount);
    for (var v = 0; v < vertexCount; v++) depth[v] = o.x[v] * ex + o.y[v] * ey + o.z[v] * ez;
    var key = new Float64Array(faceCount), low = Infinity, high = -Infinity;
    for (var f = 0; f < faceCount; f++) {
      var d = depth[o.i[f]] + depth[o.j[f]] + depth[o.k[f]];
      key[f] = d;
      if (d < low) low = d;
      if (d > high) high = d;
    }
    var bits = Math.ceil(Math.log2(faceCount + 1)), shift = Math.pow(2, bits), top = Math.pow(2, 53 - bits) - 1;
    var scale = high > low ? top / (high - low) : 0;
    var packed = new Float64Array(faceCount);
    for (var g = 0; g < faceCount; g++) packed[g] = Math.floor((key[g] - low) * scale) * shift + g;
    packed.sort();   // far (small) first, near drawn last
    var order = new Uint32Array(faceCount);
    for (var n = 0; n < faceCount; n++) order[n] = packed[n] - Math.floor(packed[n] / shift) * shift;
    return order;
  }
  function rebuild(force) {
    var scene = gd._fullLayout.scene || {};
    var camera = scene.camera || {eye: {x: 1.25, y: 1.25, z: 1.25}, center: {x: 0, y: 0, z: 0}};
    var flipZ = gd.layout.scene && gd.layout.scene.zaxis && gd.layout.scene.zaxis.autorange === 'reversed' ? -1 : 1;
    var ex = camera.eye.x - camera.center.x, ey = camera.eye.y - camera.center.y,
        ez = (camera.eye.z - camera.center.z) * flipZ;
    // The back-to-front order depends on the viewing direction only, so a zoom (which leaves it unchanged) needs no new
    // sort; skipping it keeps the page responsive right after a scroll. Anything else (cut, group, opacity) forces one.
    var length = Math.sqrt(ex * ex + ey * ey + ez * ez) || 1, direction = [ex / length, ey / length, ez / length];
    if (!force && lastDirection && direction[0] * lastDirection[0] + direction[1] * lastDirection[1] +
        direction[2] * lastDirection[2] > 0.99999) return;
    lastDirection = direction;
    var translucent = gd._fullData[meshIndex].opacity < 1;
    if (!force && !translucent) return;   // opaque cells need no re-sorting: a camera move changes nothing here
    var base = null;
    if (translucent) {
      if (!sortedAll || sortedDirection[0] * direction[0] + sortedDirection[1] * direction[1] +
          sortedDirection[2] * direction[2] < 0.99999) {
        sortedAll = sortAll(ex, ey, ez);
        sortedDirection = direction;
      }
      base = sortedAll;
    }
    var cut = gd.__sliceCut;   // {side, x, y, z} while the Slice checkbox is ticked, else absent
    var keep = new Uint32Array(faceCount), count = 0;
    for (var n = 0; n < faceCount; n++) {
      var f = base ? base[n] : n;
      if (hidden[faceGroup[f]]) continue;
      if (cut && (cut.side * (faceCx[f] - cut.x) < 0 || faceCy[f] < cut.y || faceCz[f] < cut.z)) continue;
      keep[count++] = f;
    }
    var ni = new Int32Array(count), nj = new Int32Array(count), nk = new Int32Array(count);
    for (var m = 0; m < count; m++) { var src = keep[m]; ni[m] = o.i[src]; nj[m] = o.j[src]; nk[m] = o.k[src]; }
    busy = true;
    syncCamera();
    Plotly.restyle(gd, {i: [ni], j: [nj], k: [nk]}, [meshIndex]).then(function () { busy = false; }, function () { busy = false; });
  }
  function schedule(force) {
    forced = forced || force !== false;   // only the camera handler passes false
    clearTimeout(timer);
    timer = setTimeout(function () { var f = forced; forced = false; rebuild(f); }, 60);
  }
  gd.__cellsTranslucent = function () { return gd._fullData[meshIndex].opacity < 1; };   // (the Spin button asks)
  gd.__rebuildCells = function () { rebuild(true); };  // the Slice checkbox's script calls this when its cut moves
  gd.on('plotly_relayout', function (e) {
    if (busy || !e) return;
    if (Object.keys(e).some(function (name) { return name.indexOf('scene.camera') === 0; })) schedule(false);
  });
  gd.on('plotly_restyle', function (e) {
    if (busy || !e || !e[0]) return;
    if ('opacity' in e[0] && (!e[1] || e[1].indexOf(meshIndex) >= 0)) schedule();
  });
  gd.__legendRoutes.push(function (e) {   // (see THREED_WEB_LEGEND_ROUTER_SCRIPT)
    var group = legendIndices.indexOf(e.curveNumber);
    if (group < 0) return true;  // outline / cloud entries keep the normal toggle
    hidden[group] = !hidden[group];
    Plotly.restyle(gd, {visible: hidden[group] ? 'legendonly' : true}, [e.curveNumber]);
    schedule();
    return false;
  });
})();
"""


# The "Background" radio group over the plot: Cloud (the faint whole-brain cloud plus the three wall outlines), Lines
# (slice-contour edges of the brain, the ventricles and the cells' region) or Surfaces (translucent shells of the same).
THREED_WEB_EDGE_SCRIPT = """
(function () {
  var gd = document.getElementById('{plot_id}');
  var roles = __ROLES__;
  var switchedOff = {};   // traces the user has switched off with a legend click (some start that way): they stay off
  __OFF__.forEach(function (index) { switchedOff[index] = true; });
  var modes = [['cloud', 'Cloud'], ['lines', 'Lines']];
  if (roles.surfaces.length) modes.push(['surfaces', 'Surfaces']);
  var box = document.createElement('div');
  box.className = 'background-box';   // THREED_WEB_SPIN_SCRIPT adds its axes checkbox to this box
  box.style.cssText = 'position:absolute;left:14px;top:46px;z-index:20;padding:6px 10px;' +
    'background:rgba(29,32,41,0.88);color:#e6e8ee;border:1px solid #3a3f4d;border-radius:6px;font:12px sans-serif;';
  box.appendChild(document.createTextNode('Background: '));
  function apply(mode) {
    var show = {cloud: mode === 'cloud', walls: mode === 'cloud', lines: mode === 'lines', surfaces: mode === 'surfaces'};
    Object.keys(roles).forEach(function (role) {
      roles[role].forEach(function (index) {   // note which are off now (a hidden mode's traces keep their old note)
        var state = gd.data[index].visible;
        if (state === 'legendonly') switchedOff[index] = true; else if (state === true || state === undefined) delete switchedOff[index];
      });
      if (!roles[role].length) return;
      var visible = roles[role].map(function (index) { return show[role] ? (switchedOff[index] ? 'legendonly' : true) : false; });
      Plotly.restyle(gd, {visible: visible}, roles[role]);
    });
  }
  modes.forEach(function (mode, n) {
    var label = document.createElement('label');
    label.style.cssText = 'margin-right:10px;cursor:pointer;';
    var radio = document.createElement('input');
    radio.type = 'radio'; radio.name = 'background-mode'; radio.checked = mode[0] === '__DEFAULT__';
    radio.addEventListener('change', function () {   // the button is painted first, then the (slow) redraw
      if (gd.__busyStart) gd.__busyStart();
      setTimeout(function () { apply(mode[0]); }, 30);
    });
    label.appendChild(radio);
    label.appendChild(document.createTextNode(' ' + mode[1]));
    box.appendChild(label);
  });
  gd.style.position = 'relative';
  gd.appendChild(box);
})();
"""


# The "slice" checkbox under the Background radio group: while it is ticked, three sliders place a coronal, a horizontal
# and a sagittal plane, and only what is on
# the posterior side of the coronal plane, the ventral side of the horizontal plane and the far side of the sagittal
# plane (the side away from the camera, which follows the camera when it is rotated) is shown: a cutaway whose corner you move.
# Opening positions: the coronal plane at the anterior edge and the horizontal plane at the dorsal edge of the plotted
# range (so they hide nothing) and the sagittal plane at the brain's midline (5.7 mm), which hides the near half.
# Plotly has no clipping planes, so the page re-filters the geometry itself: triangles whose centre is cut away are
# dropped from the shell meshes (and from the cells' mesh, through gd.__sliceCut, see THREED_WEB_MESH_SCRIPT), and
# points that are cut away are blanked in the line and marker traces.
THREED_WEB_SLICE_SCRIPT = """
(function () {
  var gd = document.getElementById('{plot_id}');
  var cellMesh = __MESH_INDEX__;
  var traces = [];   // what is re-filtered: {index, kind: 'mesh' | 'points', original arrays}
  gd._fullData.forEach(function (t, index) {
    if (index === cellMesh || !t.x || t.x.length < 2) return;
    if (t.type === 'mesh3d') {
      traces.push({index: index, kind: 'mesh', x: t.x, y: t.y, z: t.z, i: Int32Array.from(t.i), j: Int32Array.from(t.j),
                   k: Int32Array.from(t.k), lastKept: -1, modified: false});
    } else if (t.type === 'scatter3d') {
      // A marker-only trace can simply lose its cut points (fewer points to send is faster), together with their hover
      // text; a line trace needs its gaps, so there the cut points are blanked instead.
      traces.push({index: index, kind: 'points', x: Array.from(t.x), y: Array.from(t.y), z: Array.from(t.z), lastKept: -1,
                   modified: false, compact: t.mode === 'markers', text: t.text ? Array.from(t.text) : null});
    }
  });

  var MIDLINE_MM = 5.7;
  // A restyle makes plotly re-plot the 3D scene from the camera stored in the layout, which lags behind the live one while a
  // scroll-zoom or drag is still settling, so the view would snap back a little. Handing the live camera to the layout
  // first avoids that.
  function syncCamera() {
    try {
      var scene = gd._fullLayout.scene, live = scene._scene.getCamera();
      scene.camera = live;
      if (gd.layout.scene) gd.layout.scene.camera = live;
    } catch (err) { /* the layout keeps its own camera */ }
  }
  var axes = ['xaxis', 'yaxis', 'zaxis'].map(function (name) { return gd._fullLayout.scene[name].range; });
  var low = axes.map(function (r) { return Math.min(r[0], r[1]); }), high = axes.map(function (r) { return Math.max(r[0], r[1]); });
  // The page's scene axes are x = left-right (CCF z), y = anterior-posterior (CCF x), z = dorsal-ventral (CCF y), so
  // position[0] is the sagittal plane, [1] the coronal and [2] the horizontal one.
  var position = [Math.min(Math.max(MIDLINE_MM, low[0]), high[0]), low[1], low[2]];
  // The sagittal plane keeps the side away from the viewer, so which side that is follows the camera: the eye's CCF z is
  // found by carrying the eye (the origin of view space) back through plotly's view and model matrices (the data is
  // scaled by scene.dataScale before the model matrix acts on it).
  var scene = gd._fullLayout.scene._scene, farSide = 1;   // +1: keep x above the sagittal plane, -1: keep x below it
  function matMul(a, b) {
    var out = new Float64Array(16);
    for (var c = 0; c < 4; c++) for (var r = 0; r < 4; r++) {
      var sum = 0;
      for (var k = 0; k < 4; k++) sum += a[k * 4 + r] * b[c * 4 + k];
      out[c * 4 + r] = sum;
    }
    return out;
  }
  function eyeX() {   // Gauss-Jordan on [M | I] for M = view x model, then M^-1 applied to the origin
    var cp = scene.glplot.cameraParams, m = matMul(cp.view, cp.model), n = 4, a = [];
    for (var r = 0; r < n; r++) {
      a.push([]);
      for (var c = 0; c < n; c++) a[r].push(m[c * 4 + r]);
      for (var d = 0; d < n; d++) a[r].push(r === d ? 1 : 0);
    }
    for (var col = 0; col < n; col++) {
      var pivot = col;
      for (var r2 = col + 1; r2 < n; r2++) if (Math.abs(a[r2][col]) > Math.abs(a[pivot][col])) pivot = r2;
      var tmp = a[col]; a[col] = a[pivot]; a[pivot] = tmp;
      var div = a[col][col];
      for (var c2 = 0; c2 < 2 * n; c2++) a[col][c2] /= div;
      for (var r3 = 0; r3 < n; r3++) {
        if (r3 === col) continue;
        var f = a[r3][col];
        for (var c3 = 0; c3 < 2 * n; c3++) a[r3][c3] -= f * a[col][c3];
      }
    }
    // inverse column 3 (the image of the origin) is the eye in scaled data coordinates; its w is entry [3][7]
    return a[0][7] / a[3][7] / scene.dataScale[0];
  }
  function updateSide() { var side = eyeX() < position[0] ? 1 : -1; var changed = side !== farSide; farSide = side; return changed; }
  function cutAway(x, y, z) { return farSide * (x - position[0]) < 0 || y < position[1] || z < position[2]; }

  // --- the controls, under the Background radio group
  var box = document.createElement('div');
  box.style.cssText = 'position:absolute;left:14px;top:86px;z-index:20;padding:6px 10px;' +
    'background:rgba(29,32,41,0.88);color:#e6e8ee;border:1px solid #3a3f4d;border-radius:6px;font:12px sans-serif;';
  box.className = 'slice-box';   // THREED_WEB_SPIN_SCRIPT moves it down under the View box
  var checkboxLabel = document.createElement('label');
  checkboxLabel.style.cssText = 'cursor:pointer;';
  var checkbox = document.createElement('input');
  checkbox.type = 'checkbox';
  checkboxLabel.appendChild(checkbox);
  checkboxLabel.appendChild(document.createTextNode(' slice'));
  var sliderRows = document.createElement('div');
  sliderRows.style.cssText = 'display:none;';
  var names = ['Sagittal (left to right)', 'Coronal (anterior to posterior)', 'Horizontal (dorsal to ventral)'];
  [1, 2, 0].forEach(function (axis) {   // shown as coronal, horizontal, sagittal
    var name = names[axis];
    var row = document.createElement('div');
    row.style.cssText = 'margin-top:6px;';
    var text = document.createElement('div');
    var input = document.createElement('input');
    input.type = 'range'; input.min = low[axis]; input.max = high[axis]; input.step = 0.05; input.value = position[axis];
    input.style.cssText = 'width:170px;margin:3px 0 0 0;display:block;';
    function label() { text.textContent = name + ': ' + position[axis].toFixed(2) + ' mm'; }
    label();
    input.addEventListener('input', function () { position[axis] = parseFloat(input.value); label(); sliderMoved(); });
    row.appendChild(text); row.appendChild(input); sliderRows.appendChild(row);
  });
  box.appendChild(checkboxLabel); box.appendChild(sliderRows);

  // --- re-filtering the geometry (one pass at a time; a slider moved meanwhile triggers one more pass afterwards).
  // What makes a drag quick: (1) each shell's face centres are computed once; (2) an object that is not drawn (switched off in
  // the legend) is skipped, and brought up to date when it is switched on; (3) an object whose kept set did not change (the
  // plane is not crossing it) is not sent to plotly at all, since pushing a mesh to the graphics card is the slow part;
  // (4) the shells are brought up to date on every pass along with the cells (the cells' own mesh is filtered by
  // THREED_WEB_MESH_SCRIPT through gd.__sliceCut); while a pass is under way, slider moves are collected into one more pass.
  var active = false, busy = false, pending = false, timer = 0, heavyReady = true, heavyTimer = 0;
  var cutState = {side: 1, x: 0, y: 0, z: 0};
  function centres(t) {
    if (t.cx) return;
    var n = t.i.length;
    t.cx = new Float32Array(n); t.cy = new Float32Array(n); t.cz = new Float32Array(n);
    for (var f = 0; f < n; f++) {
      var a = t.i[f], b = t.j[f], c = t.k[f];
      t.cx[f] = (t.x[a] + t.x[b] + t.x[c]) / 3; t.cy[f] = (t.y[a] + t.y[b] + t.y[c]) / 3; t.cz[f] = (t.z[a] + t.z[b] + t.z[c]) / 3;
    }
  }
  function apply(restore, withMeshes) {
    cutState.side = farSide; cutState.x = position[0]; cutState.y = position[1]; cutState.z = position[2];
    gd.__sliceCut = restore ? null : cutState;
    var meshes = {i: [], j: [], k: []}, meshIds = [], points = {x: [], y: [], z: [], text: []}, pointIds = [];
    traces.forEach(function (t) {
      if (restore) {   // put back whatever was filtered
        if (!t.modified) return;
        t.modified = false; t.lastKept = -1;
        if (t.kind === 'mesh') { meshes.i.push(t.i); meshes.j.push(t.j); meshes.k.push(t.k); meshIds.push(t.index); }
        else { points.x.push(t.x); points.y.push(t.y); points.z.push(t.z); points.text.push(t.text); pointIds.push(t.index); }
        return;
      }
      if (gd._fullData[t.index].visible !== true) { t.lastKept = -1; return; }   // not drawn: updated when switched on
      if (t.kind === 'mesh') {
        if (!withMeshes) return;
        centres(t);
        var n = t.i.length, keep = new Uint32Array(n), count = 0;
        for (var f = 0; f < n; f++) {
          if (!(cutState.side * (t.cx[f] - cutState.x) < 0 || t.cy[f] < cutState.y || t.cz[f] < cutState.z)) keep[count++] = f;
        }
        if (count === t.lastKept) return;   // the same faces as already shown
        var ni = new Int32Array(count), nj = new Int32Array(count), nk = new Int32Array(count);
        for (var m = 0; m < count; m++) { ni[m] = t.i[keep[m]]; nj[m] = t.j[keep[m]]; nk[m] = t.k[keep[m]]; }
        t.lastKept = count; t.modified = true;
        meshes.i.push(ni); meshes.j.push(nj); meshes.k.push(nk); meshIds.push(t.index);
      } else {
        var total = t.x.length, nx = [], ny = [], nz = [], nt = t.text ? [] : null, blanked = 0;
        for (var q = 0; q < total; q++) {
          var x = t.x[q], y = t.y[q], z = t.z[q];
          var blank = x === null || y === null || z === null || cutState.side * (x - cutState.x) < 0 ||
                      y < cutState.y || z < cutState.z;
          if (blank) {
            blanked++;
            if (t.compact) continue;
            nx.push(null); ny.push(null); nz.push(null);
            if (nt) nt.push(t.text[q]);
          } else {
            nx.push(x); ny.push(y); nz.push(z);
            if (nt) nt.push(t.text[q]);
          }
        }
        if (blanked === t.lastKept) return;
        t.lastKept = blanked; t.modified = true;
        points.x.push(nx); points.y.push(ny); points.z.push(nz); points.text.push(nt); pointIds.push(t.index);
      }
    });
    busy = true;
    syncCamera();
    var jobs = [];
    if (meshIds.length) jobs.push(Plotly.restyle(gd, meshes, meshIds));
    if (pointIds.length) jobs.push(Plotly.restyle(gd, points, pointIds));
    if (gd.__rebuildCells) gd.__rebuildCells();
    function done() { busy = false; if (pending) { pending = false; request(); } }
    Promise.all(jobs).then(done, done);
  }
  function run() { timer = 0; if (!active) return; if (busy) { pending = true; return; } apply(false, heavyReady); }
  function request() { if (active && !timer) timer = setTimeout(run, 40); }
  function sliderMoved() {   // the shells follow the slider along with the cells (one pass at a time, newest position wins)
    heavyReady = true;
    request();
  }
  gd.__sliceReload = function (indices) {   // some meshes were given new geometry: take it as their original
    traces.forEach(function (t) {
      if (indices.indexOf(t.index) < 0 || t.kind !== 'mesh') return;
      var full = gd._fullData[t.index];
      t.x = full.x; t.y = full.y; t.z = full.z;
      t.i = Int32Array.from(full.i); t.j = Int32Array.from(full.j); t.k = Int32Array.from(full.k);
      t.cx = null; t.lastKept = -1; t.modified = false;
    });
    if (active) { heavyReady = true; request(); }
  };
  gd.on('plotly_restyle', function (e) {   // an object was switched on or off: bring the newly drawn ones up to date
    if (active && e && e[0] && 'visible' in e[0]) { heavyReady = true; request(); }
  });

  gd.on('plotly_relayout', function (e) {   // the camera moved: the far side of the sagittal plane may have changed
    // With opaque cells (100% opacity) the cut keeps the side chosen when Slice was ticked, so a camera move costs no
    // redraw at all; the far side follows the camera only while the cells are translucent.
    var translucent = !!(gd.__cellsTranslucent && gd.__cellsTranslucent());
    if (active && translucent && e && Object.keys(e).some(function (key) { return key.indexOf('scene.camera') === 0; })) {
      [120, 600].forEach(function (delay) {   // (plotly updates its matrices just after the event, when it next paints)
        setTimeout(function () { if (active && updateSide()) request(); }, delay);
      });
    }
  });
  checkbox.addEventListener('change', function () {
    active = checkbox.checked;
    gd.__sliceOn = active;   // (the Spin button asks)
    if (active) updateSide();
    sliderRows.style.display = active ? 'block' : 'none';
    if (gd.__busyStart) gd.__busyStart();
    if (active) { heavyReady = true; request(); } else { clearTimeout(timer); timer = 0; setTimeout(function () { apply(true, true); }, 30); }
  });
  gd.style.position = 'relative';
  gd.appendChild(box);
})();
"""


# The "Spin" button, to the right of the Background radio buttons: while it is on, the camera turns smoothly around the
# vertical axis (the azimuth advances, the angle above the horizontal and the distance stay as they are, so the user's own
# tilt and zoom are kept, and a drag or scroll in between simply carries on from there) until the button is clicked again.
# The frames set the live camera directly (plotly's scene.setViewport) rather than through Plotly.relayout, which would redraw the
# whole figure and run the page's other camera handlers on every frame; the layout's camera is brought up to date with
# a cheap camera event now and then, only while the cells are translucent (so the depth order and the Slice view's far side follow), and a proper relayout when the spin stops.
THREED_WEB_SPIN_SCRIPT = """
(function () {
  var gd = document.getElementById('{plot_id}');
  var SECONDS_PER_TURN = 24;
  var scene = gd._fullLayout.scene._scene;
  var box = gd.querySelector('.background-box');
  if (!box) {   // a page without background shells has no radio box: make one for the axes checkbox alone
    box = document.createElement('div');
    box.style.cssText = 'position:absolute;left:14px;top:46px;z-index:20;padding:6px 10px;' +
      'background:rgba(29,32,41,0.88);color:#e6e8ee;border:1px solid #3a3f4d;border-radius:6px;font:12px sans-serif;';
    gd.style.position = 'relative';
    box.appendChild(document.createTextNode('Background: '));
    gd.appendChild(box);
  }
  // The axes checkbox, in the Background box: the grid, tick labels and axis titles of the scene on and off.
  var axesLabel = document.createElement('label');
  axesLabel.style.cssText = 'margin-left:4px;cursor:pointer;';
  var axesBox = document.createElement('input');
  axesBox.type = 'checkbox'; axesBox.checked = true;
  axesBox.addEventListener('change', function () {
    var visible = axesBox.checked;
    Plotly.relayout(gd, {'scene.xaxis.visible': visible, 'scene.yaxis.visible': visible, 'scene.zaxis.visible': visible});
  });
  axesLabel.appendChild(axesBox);
  axesLabel.appendChild(document.createTextNode(' Axes'));
  box.appendChild(axesLabel);
  var button = document.createElement('button');
  button.type = 'button';
  button.textContent = 'Spin';
  var IDLE = 'margin-left:6px;padding:1px 10px;font:12px sans-serif;cursor:pointer;border-radius:4px;' +
    'background:#2b2f3b;color:#e6e8ee;border:1px solid #5a6072;';
  var ACTIVE = IDLE.replace('background:#2b2f3b', 'background:#3b6fd4').replace('#5a6072', '#8fb0ff');
  button.style.cssText = IDLE;

  var spinning = false, last = 0, lastSync = 0;
  // Toggling an object (a legend click, a radio button) makes plotly re-plot the scene from the camera stored in the
  // layout, which is stale while the spin runs, so the view would jump back to where the spin began. The live camera is
  // handed to the layout just before any such change.
  function adoptLiveCamera() {
    try {
      var live = scene.getCamera();
      gd._fullLayout.scene.camera = live;
      if (gd.layout.scene) gd.layout.scene.camera = live;
    } catch (err) { /* the layout keeps its own camera */ }
  }
  // (A press on the legend is caught as a DOM event rather than with gd.on('plotly_legendclick'): plotly takes the LAST
  // legend-click handler's return value to decide whether to go on with its own toggle, and a handler returning nothing
  // here overrode the cells' handler (THREED_WEB_MESH_SCRIPT) that cancels it, so a group switched off came straight back.)
  gd.addEventListener('mousedown', function (event) {
    if (event.target && event.target.closest && event.target.closest('.legend')) adoptLiveCamera();
  }, true);
  ['restyle', 'relayout', 'update'].forEach(function (name) {
    var original = Plotly[name];
    Plotly[name] = function () { adoptLiveCamera(); return original.apply(this, arguments); };
  });
  function syncLayout(full) {   // hand the live camera to the layout and let the page's camera handlers run
    var live = scene.getCamera();
    if (full) { Plotly.relayout(gd, {'scene.camera': live}); return; }
    gd._fullLayout.scene.camera = live;   // cheap: no redraw, just the camera event the handlers listen to
    if (gd.layout.scene) gd.layout.scene.camera = live;
    gd.emit('plotly_relayout', {'scene.camera': live});
  }
  // The cells' back-to-front sort (and, with Slice on, the cut's far side) matters only while the cells are translucent, and
  // the envelopes' own sort only while one is switched on, so
  // that is the only reason to catch up while spinning (like a mouse drag, which catches up when it ends). With opaque
  // cells nothing is done per frame beyond moving the camera, which is what keeps the spin smooth.
  function needsCatchUp() {
    return !!(gd.__cellsTranslucent && gd.__cellsTranslucent()) || !!(gd.__envelopesNeedSort && gd.__envelopesNeedSort());
  }
  function frame(now) {
    if (!spinning) return;
    var dt = Math.min(0.1, (now - last) / 1000);
    last = now;
    var camera = scene.getCamera(), c = camera.center, e = camera.eye;
    var dx = e.x - c.x, dy = e.y - c.y, angle = 2 * Math.PI * dt / SECONDS_PER_TURN;
    var cos = Math.cos(angle), sin = Math.sin(angle);
    scene.setViewport({camera: {eye: {x: c.x + dx * cos - dy * sin, y: c.y + dx * sin + dy * cos, z: e.z},
                                center: c, up: camera.up, projection: camera.projection},
                       aspectratio: gd._fullLayout.scene.aspectratio});
    if (now - lastSync > 700 && needsCatchUp()) { lastSync = now; syncLayout(false); }
    requestAnimationFrame(frame);
  }
  function setSpinning(on) {
    spinning = on;
    button.style.cssText = spinning ? ACTIVE : IDLE;
    button.textContent = spinning ? 'Stop spin' : 'Spin';
    if (spinning) {
      last = performance.now(); lastSync = last;
      requestAnimationFrame(frame);
    } else {
      syncLayout(true);
    }
  }
  button.addEventListener('click', function () { setSpinning(!spinning); });

  // The "View" box, under the Background box: Spin, and Top, Sagittal, Anterior, Posterior and Diagonal put the camera on that side of the brain, at
  // the distance it is now (so the zoom is kept), with anterior up in the top view and left in the sagittal one. The axes
  // are x = left - right, y = anterior - posterior (growing backwards), z = dorsal - ventral (shown reversed, so up on the screen).
  var views = document.createElement('div');
  var viewsTop = 46 + box.offsetHeight + 6;   // (under the Background box)
  views.style.cssText = 'position:absolute;left:14px;top:' + viewsTop + 'px;z-index:20;padding:6px 10px;' +
    'background:rgba(29,32,41,0.88);color:#e6e8ee;border:1px solid #3a3f4d;border-radius:6px;font:12px sans-serif;';
  views.appendChild(document.createTextNode('View: '));
  [['Top', [0, 0, 1], [0, -1, 0]], ['Sagittal', [1, 0, 0], [0, 0, 1]], ['Anterior', [0, -1, 0], [0, 0, 1]], ['Posterior', [0, 1, 0], [0, 0, 1]],
   ['Diagonal', [1, -1, 1], [0, 0, 1]]].forEach(function (view, n) {
    var b = document.createElement('button');
    b.type = 'button';
    b.textContent = view[0];
    b.style.cssText = IDLE.replace('margin-left:6px', 'margin-left:' + (n ? 4 : 0) + 'px');
    b.addEventListener('click', function () {
      if (spinning) setSpinning(false);
      var camera = scene.getCamera(), c = camera.center, e = camera.eye;
      var d = Math.sqrt(Math.pow(e.x - c.x, 2) + Math.pow(e.y - c.y, 2) + Math.pow(e.z - c.z, 2)) || 2;
      var v = view[1], norm = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      Plotly.relayout(gd, {'scene.camera': {eye: {x: c.x + d * v[0] / norm, y: c.y + d * v[1] / norm, z: c.z + d * v[2] / norm},
                                            center: c, up: {x: view[2][0], y: view[2][1], z: view[2][2]},
                                            projection: camera.projection}});
    });
    views.appendChild(b);
  });
  views.appendChild(button);   // Spin belongs with the views
  gd.appendChild(views);
  var sliceBox = gd.querySelector('.slice-box');   // the slice checkbox and its sliders sit under the View box, not beneath it
  if (sliceBox) sliceBox.style.top = (viewsTop + (views.offsetHeight || 32) + 6) + 'px';
})();
"""


# A cell group switched off in the legend has its dot drawn as an open circle in its own colour (so "off" does not depend on
# how bright the colour was to begin with). The legend's swatches for the background objects (the surfaces, the line edges
# in Lines mode, the envelopes) are faint, since they take the objects' own low opacity, so switching one on or off barely
# changed them: each is now drawn as a square in the object's own colour, filled while it is on and open while it is off.
# Drawn over plotly's own legend after every plot (plotly redraws the legend itself, which would drop anything added before).
THREED_WEB_LEGEND_MARK_SCRIPT = """
(function () {
  var gd = document.getElementById('{plot_id}');
  var SVG = 'http://www.w3.org/2000/svg';
  function mark(parent, tag, attributes) {
    var element = document.createElementNS(SVG, tag);
    Object.keys(attributes).forEach(function (key) { element.setAttribute(key, attributes[key]); });
    element.setAttribute('class', 'legend-mark');
    element.style.pointerEvents = 'none';
    parent.appendChild(element);
  }
  // A click on a legend entry is acted on by plotly only after its short wait for a possible double click, and then after the
  // redraw: so the entry is switched in the picture AT ONCE, ahead of the real state (`pending`: trace number -> state),
  // and the real state takes over when the redraw is done (or after a moment, if the click turned out not to count).
  var pending = {};
  function shownState(trace) { return pending[trace.index] !== undefined ? pending[trace.index] : trace.visible; }
  function decorate() {
    Array.prototype.forEach.call(gd.querySelectorAll('.legend g.traces'), function (group) {
      var data = group.__data__, trace = data && data[0] && data[0].trace;
      if (trace && /^[\\u200b\\s]*$/.test(trace.name || '')) {   // a blank-named stand-in trace: just a gap between legend groups
        group.style.setProperty('visibility', 'hidden');
        return;
      }
      if (trace && pending[trace.index] !== undefined) {
        if (pending[trace.index] === trace.visible) delete pending[trace.index];   // the real state has caught up
        else group.style.opacity = pending[trace.index] === 'legendonly' ? 0.5 : 1;
      }
      Array.prototype.forEach.call(group.querySelectorAll('.legend-mark'), function (old) { old.remove(); });
      if (trace && !/ (surface|edges|envelope)$/.test(trace.name || '') && trace.marker &&
          typeof trace.marker.color === 'string') {
        // A cell group: its dot becomes an open circle in its own colour while the group is switched off.
        Array.prototype.forEach.call(group.querySelectorAll('.legendpoints path'), function (dot) {
          var saved = dot.getAttribute('data-original-style');
          if (shownState(trace) === 'legendonly') {
            if (saved === null) dot.setAttribute('data-original-style', dot.getAttribute('style') || '');
            dot.style.setProperty('fill', 'none', 'important');
            dot.style.setProperty('stroke', trace.marker.color, 'important');
            dot.style.setProperty('stroke-width', '2px', 'important');
          } else if (saved !== null) {   // switched back on: put back exactly what plotly had drawn
            dot.setAttribute('style', saved);
            dot.removeAttribute('data-original-style');
          }
        });
        return;
      }
      if (!trace || !/ (surface|edges|envelope)$/.test(trace.name || '')) return;
      // A background object (surface, edges, envelope): plotly's own swatch is replaced by a square in the object's own
      // colour, filled while it is on and open (outline only) while it is off, like the cell groups' dot and open circle.
      var color = (typeof trace.color === 'string' && trace.color) || (trace.line && trace.line.color) ||
                  (trace.marker && typeof trace.marker.color === 'string' && trace.marker.color) || '#c9d1e0';
      var layers = group.querySelector('.layers');
      if (layers) layers.style.setProperty('display', 'none');
      var swatch = group.querySelector('.legendpoints'), box = swatch && swatch.getBBox ? swatch.getBBox() : null;
      var size = 12, x = box && box.width ? box.x : 14, y = box && box.height ? box.y : -6;
      var off = shownState(trace) === 'legendonly';
      mark(group, 'rect', {x: x + 1, y: y + 1, width: size - 2, height: size - 2, rx: 2, ry: 2,
                           fill: off ? 'none' : color, stroke: color, 'stroke-width': 1.6});
    });
  }
  gd.on('plotly_afterplot', decorate);
  gd.on('plotly_restyle', function () { setTimeout(decorate, 0); });
  // (On the press, not the 'click' event: plotly puts a cover over the page while the mouse is down and takes the release
  // itself, so the browser never sends the toggle a 'click' of its own.)
  gd.addEventListener('mousedown', function (event) {   // (capture phase: before plotly's own handling of the press)
    if (event.button !== 0) return;
    var toggle = event.target && event.target.closest && event.target.closest('.legendtoggle');
    var group = toggle && toggle.closest('g.traces'), data = group && group.__data__;
    var trace = data && data[0] && data[0].trace;
    if (!trace || trace.index === undefined) return;
    pending[trace.index] = shownState(trace) === 'legendonly' ? true : 'legendonly';
    decorate();
    if (gd.__busyStart) gd.__busyStart();
    setTimeout(function () { delete pending[trace.index]; decorate(); }, 2500);   // (a click that did not count: back to the truth)
  }, true);
  decorate();
})();
"""


# Opening settings of the browser view's two background controls. The Background opacity slider scales every background
# trace (cloud, wall outlines, edges, shells) from its own 100% look; the Brain surface opacity slider below sets the
# brain shell's own base, and the shell's opacity is that base x the Background slider's factor (so the defaults give
# the shell 0.5 x 0.2 = 10%).
THREED_WEB_BACKGROUND_DEFAULT = 0.2
THREED_WEB_BRAIN_SURFACE_DEFAULT = 0.5
THREED_WEB_SLIDER_ROW_PX = 36   # height of one slider row at the bottom of the page

# The opacity controls at the bottom-left, top to bottom: Brain surface, Background and Cell opacity, all the same kind of
# control (a plain 0-100% slider that updates the figure as it is dragged, with 100% at the right end). Background scales
# every background trace (cloud, wall outlines, edges, shells) from its own 100% look; Brain surface sets the brain shell's
# own base, and the shell's opacity is that base x the Background slider's factor (so the defaults give the shell
# 0.5 x 0.2 = 10%); Cell opacity sets the cells' traces and also scales the cell groups' envelope surfaces, which have
# their own Envelope opacity slider (up to 100%) and two sliders, cells inside and smoothness, that swap in another
# precomputed surface. Changes are batched (at most one redraw every few dozen ms).
THREED_WEB_OPACITY_SCRIPT = """
(function () {
  var gd = document.getElementById('{plot_id}');
  var config = __CONFIG__;   // {cells: [trace indices], envelopes: {mesh, standIns, colors, base, levels, ...} | null, background, brain}
  var env = config.envelopes;
  var state = {cells: 1.0, background: config.background ? config.background.value : 1.0,
               brain: config.brain ? config.brain.value : 1.0, envelope: env ? env.base : 1.0,
               coverage: env ? env.defaultCoverage : 0, smooth: env ? env.defaultSmooth : 0};
  var dirty = {}, timer = 0;
  var stack = document.createElement('div');
  stack.style.cssText = 'position:absolute;right:14px;bottom:10px;z-index:20;color:#e6e8ee;font:11px sans-serif;';
  function addSlider(title, key, options) {
    options = options || {};
    var block = document.createElement('div');
    block.style.cssText = 'height:__ROW__px;';
    var text = document.createElement('div');
    var input = document.createElement('input');
    input.type = 'range'; input.step = 1;
    input.min = options.min !== undefined ? options.min : 0;
    input.max = options.max !== undefined ? options.max : 100;
    input.value = options.toSlider ? options.toSlider(state[key]) : Math.round(state[key] * 100);
    input.style.cssText = 'width:150px;margin:3px 0 0 0;display:block;';
    function label() { text.textContent = title + ': ' + (options.show ? options.show(state[key]) : Math.round(state[key] * 100) + '%'); }
    input.addEventListener('input', function () {
      state[key] = options.fromSlider ? options.fromSlider(parseInt(input.value, 10)) : input.value / 100;
      label(); dirty[key] = true;
      if (!timer) timer = setTimeout(flush, 30);
    });
    label();
    block.appendChild(text); block.appendChild(input); stack.appendChild(block);
  }
  function flush() {
    timer = 0;
    var changed = dirty; dirty = {};
    if (changed.background) {
      var b = config.background;
      var opacities = b.indices.map(function (index, n) {
        return config.brain && index === config.brain.index ? state.brain * state.background
                                                            : Math.round(state.background * b.bases[n] * 1000) / 1000;
      });
      Plotly.restyle(gd, {opacity: opacities}, b.indices);
    } else if (changed.brain) {
      Plotly.restyle(gd, {opacity: state.brain * (config.background ? state.background : 1)}, [config.brain.index]);
    }
    if (changed.cells && config.cells.length) Plotly.restyle(gd, {opacity: state.cells}, config.cells);
    if ((changed.cells || changed.envelope) && env) {   // the envelopes' opacity is their own slider's value times the cells'
      Plotly.restyle(gd, {opacity: Math.round(state.envelope * state.cells * 1000) / 1000}, [env.mesh]);
      relistEnvelopes(false);   // (see-through envelopes list only their near-facing triangles)
    }
    if ((changed.coverage || changed.smooth) && env) {   // swap every envelope's surface for the chosen pair's
      loadShapes();
      relistEnvelopes(true);
    }
  }
  // --- the envelopes: one mesh holding every envelope (a colour per vertex), its triangles listed far to near.
  // Plotly draws a translucent mesh's triangles in the order they are listed, so with the listing as the surface-maker left
  // it a nearer part can be painted before a farther one and end up hidden by it; and separate traces would be painted
  // one after another whatever their depth. So the triangles of every switched-on envelope are re-listed together, far to
  // near for the current viewing direction, whenever the camera has stopped moving, an envelope is switched on or off,
  // or a surface is swapped. The Slice cut, if on, takes the re-listed mesh as its own starting point (gd.__sliceReload).
  var shapes = [], offsets = [], shown = env ? env.colors.map(function () { return false; }) : [];
  var merged = null, sortTimer = 0;
  function loadShapes() {   // every envelope's surface at the chosen coverage and smoothness, merged into one set of vertices
    var x = [], y = [], z = [], color = [];
    offsets = [];
    shapes = env.levels.map(function (rows, g) {
      var row = rows[Math.min(state.smooth, rows.length - 1)], mesh = row[Math.min(state.coverage, row.length - 1)];
      offsets.push(x.length);
      for (var v = 0; v < mesh.x.length; v++) { x.push(mesh.x[v]); y.push(mesh.y[v]); z.push(mesh.z[v]); color.push(env.colors[g]); }
      return {x: mesh.x, y: mesh.y, z: mesh.z, i: Int32Array.from(mesh.i), j: Int32Array.from(mesh.j), k: Int32Array.from(mesh.k)};
    });
    merged = {x: x, y: y, z: z, color: color};
  }
  function viewDirection() {
    var camera = gd._fullLayout.scene._scene.getCamera();
    var flipZ = gd.layout.scene && gd.layout.scene.zaxis && gd.layout.scene.zaxis.autorange === 'reversed' ? -1 : 1;
    return [camera.eye.x - camera.center.x, camera.eye.y - camera.center.y, (camera.eye.z - camera.center.z) * flipZ];
  }
  // While an envelope is see-through its far side only confuses the picture (it is drawn behind, yet plotly's blending lets it
  // show through and, where the two sides overlap, appear to cover the near one), so then only the triangles facing the
  // viewer are listed. Which way is "facing" depends on how the surface-maker wound its triangles and on the axis flips, so
  // it is found from the surface itself: the nearest few percent of its triangles must be facing the viewer.
  function visibleFaces(g, direction, cull, out) {   // appends envelope g's listed triangles (merged vertex numbers) and depths to out
    var shape = shapes[g], vertexCount = shape.x.length, faceCount = shape.i.length, depth = new Float64Array(vertexCount);
    for (var v = 0; v < vertexCount; v++) depth[v] = shape.x[v] * direction[0] + shape.y[v] * direction[1] + shape.z[v] * direction[2];
    var key = new Float64Array(faceCount), facing = new Float64Array(faceCount), order = new Uint32Array(faceCount);
    for (var f = 0; f < faceCount; f++) {
      var a = shape.i[f], b = shape.j[f], c = shape.k[f];
      key[f] = depth[a] + depth[b] + depth[c]; order[f] = f;
      var ux = shape.x[b] - shape.x[a], uy = shape.y[b] - shape.y[a], uz = shape.z[b] - shape.z[a];
      var wx = shape.x[c] - shape.x[a], wy = shape.y[c] - shape.y[a], wz = shape.z[c] - shape.z[a];
      facing[f] = (uy * wz - uz * wy) * direction[0] + (uz * wx - ux * wz) * direction[1] + (ux * wy - uy * wx) * direction[2];
    }
    var sign = 1;
    if (cull && faceCount) {
      order.sort(function (p, q) { return key[p] - key[q]; });
      var nearest = Math.max(1, Math.floor(faceCount * 0.03)), vote = 0;
      for (var m = faceCount - nearest; m < faceCount; m++) vote += facing[order[m]] > 0 ? 1 : -1;
      sign = vote >= 0 ? 1 : -1;
    }
    var base = offsets[g];
    for (var n = 0; n < faceCount; n++) {
      if (cull && facing[n] * sign <= 0) continue;
      out.i.push(shape.i[n] + base); out.j.push(shape.j[n] + base); out.k.push(shape.k[n] + base); out.key.push(key[n]);
    }
  }
  function relistEnvelopes(withGeometry) {
    if (!env) return Promise.resolve();
    var direction = viewDirection(), cull = state.envelope * state.cells < 0.995;
    var all = {i: [], j: [], k: [], key: []};
    shown.forEach(function (on, g) { if (on) visibleFaces(g, direction, cull, all); });
    var count = all.key.length, order = new Uint32Array(count);
    for (var n = 0; n < count; n++) order[n] = n;
    order.sort(function (p, q) { return all.key[p] - all.key[q]; });   // far (small) first
    var ni = new Int32Array(count), nj = new Int32Array(count), nk = new Int32Array(count);
    for (var m = 0; m < count; m++) { var src = order[m]; ni[m] = all.i[src]; nj[m] = all.j[src]; nk[m] = all.k[src]; }
    var update = {i: [ni], j: [nj], k: [nk]};
    if (withGeometry) { update.x = [merged.x]; update.y = [merged.y]; update.z = [merged.z]; update.vertexcolor = [merged.color]; }
    return Plotly.restyle(gd, update, [env.mesh]).then(function () { if (gd.__sliceReload) gd.__sliceReload([env.mesh]); });
  }
  function scheduleRelist() { clearTimeout(sortTimer); sortTimer = setTimeout(function () { relistEnvelopes(false); }, 150); }
  gd.__envelopesNeedSort = function () {   // (the Spin button asks: it catches up only while something needs it)
    return shown.some(function (on) { return on; });
  };
  if (env) {
    loadShapes();
    relistEnvelopes(true);   // (nothing is switched on yet: this puts the mesh in step with the page's own state)
    gd.on('plotly_relayout', function (e) {
      if (e && Object.keys(e).some(function (name) { return name.indexOf('scene.camera') === 0; })) scheduleRelist();
    });
    gd.__legendRoutes.push(function (e) {   // a click on an envelope's legend entry (see THREED_WEB_LEGEND_ROUTER_SCRIPT)
      var g = env.standIns.indexOf(e.curveNumber);
      if (g < 0) return true;
      shown[g] = !shown[g];
      Plotly.restyle(gd, {visible: shown[g] ? true : 'legendonly'}, [e.curveNumber]);
      relistEnvelopes(false);
      return false;
    });
  }
  if (config.brain) addSlider('Brain surface opacity', 'brain');
  if (config.background) addSlider('Background opacity', 'background');
  if (config.cells.length) addSlider('Cell opacity', 'cells');
  if (env) {
    addSlider('Envelope opacity', 'envelope');
    if (env.coverages.length > 1) {
      addSlider('Envelope: cells inside', 'coverage', {
        min: 1, max: env.coverages.length, toSlider: function (level) { return level + 1; },
        fromSlider: function (value) { return value - 1; },
        show: function (level) { return env.coverages[level] + '%'; }});
    }
    if (env.smoothSteps > 1) {
      addSlider('Envelope: smoothness', 'smooth', {
        min: 1, max: env.smoothSteps, toSlider: function (level) { return level + 1; },
        fromSlider: function (value) { return value - 1; },
        show: function (level) { return (level + 1) + ' of ' + env.smoothSteps; }});
    }
  }
  gd.style.position = 'relative';
  gd.appendChild(stack);
})();
"""






def prepare_coordinate_frame(abc_cache, parent=None):
    """Decides which coordinates this session places the cells by (see FRAME_STATE in abc_core.py): Allen's reconstructed ones when the setting
    asks for them (or the run being opened was made with them), once their files are in the cache (asked about, and downloaded, if not); the
    raw ones otherwise, including when the download is declined."""
    want = (FRAME_STATE['override'] == 'reconstructed') if FRAME_STATE['override'] is not None else bool(COORD_SETTINGS['reconstructed_coords'])
    FRAME_STATE['ready'] = False
    if not want:
        return
    if ensure_ccf_downloaded(abc_cache, parent=parent, files=RECON_FILES, what="reconstructed coordinates and the matching atlas"):
        FRAME_STATE['ready'] = True
    else:
        print("The reconstructed coordinates are not available, so the cells are placed by the raw coordinates this time.")


def write_run_frame(run_folder):
    """Records in a run's folder which coordinates its ROIs were drawn on (see read_run_frame)."""
    try:
        (Path(run_folder) / 'coordinate_frame.txt').write_text(coordinate_frame_name() + '\n', encoding='utf-8')
    except Exception as e:
        print(f"Could not record the coordinates of this run ({e}).")


def read_run_frame(run_folder):
    """'raw' or 'reconstructed': the coordinates a saved run's ROIs were drawn on. A run saved before this was recorded used the raw ones."""
    try:
        text = (Path(run_folder) / 'coordinate_frame.txt').read_text(encoding='utf-8').strip()
        return 'reconstructed' if text == 'reconstructed' else 'raw'
    except Exception:
        return 'raw'


def ensure_ccf_downloaded(abc_cache, parent=None, files=CCF_FILES, what="CCF coordinates", note=""):
    """True once every CCF file (or each of `files`, as in RECON_FILES) is in the local cache, downloading
    them first if need be: asks (naming the size, and saying `note` if given) before downloading, then shows a progress window fed
    by the bytes landing on disk. False if the user declined or cancelled, or the download failed."""
    missing = ccf_files_missing(abc_cache, files)
    if not missing:
        return True
    kinds = {(entry[0], entry[1]): (entry[2] if len(entry) > 2 else 'metadata') for entry in files}
    total_bytes = sum(size or 0 for _d, _f, size, _p in missing)
    size_text = f"{total_bytes / 1e6:.0f} MB" if total_bytes else "a large file"
    # Already on their way in the atlas download at start-up (start_atlas_prefetch): not asked about, and not downloaded a second time beside it
    # (the per-file lock of serialize_downloads makes the download below wait for that one, and then find the files there).
    queued = atlas_download_running() is not None and all((d, f) in atlas_prefetch['keys'] for d, f, _s, _p in missing)
    if not queued:
        question = (f"This needs the Allen {what}, a one-time download of "
                    f"{size_text} into your atlas cache folder.\n\n" + (f"{note}\n\n" if note else "") + "Download it now?")
        try:
            if parent is None:                          # (no window to put it over: in the middle of the screen, see centered_dialog.py)
                import centered_dialog
                answer = centered_dialog.ask_yes_no("Download CCF coordinates?", question)
            else:
                from tkinter import messagebox
                answer = messagebox.askyesno("Download CCF coordinates?", question, parent=parent)
            if not answer:
                return False
        except Exception as e:
            print(f"Could not ask about the CCF download ({e}); not downloading.")
            return False

    outcome = {'error': None, 'done': False}

    def worker():
        try:
            for directory, file_name, _size, _path in missing:
                if kinds.get((directory, file_name)) == 'data':
                    abc_cache.get_data_path(directory=directory, file_name=file_name)
                else:
                    abc_cache.get_metadata_path(directory=directory, file_name=file_name)
        except Exception as e:
            outcome['error'] = e
        finally:
            outcome['done'] = True

    threading.Thread(target=worker, daemon=True).start()
    restore_cursors = set_wait_cursor_on_open_figures()
    proc_fig, cancel_flag, set_progress = show_processing_dialog(
        f"Downloading the {what} ({size_text}); a one-time download.")
    started = time.perf_counter()
    try:
        set_progress(0.0)
        while not outcome['done'] and not cancel_flag['cancelled']:
            got = 0
            for _d, _f, _size, local_path in missing:
                try:
                    got += local_path.stat().st_size if local_path is not None else 0
                except OSError:
                    pass
            if total_bytes:
                set_progress(min(1.0, got / total_bytes))
                proc_fig._set_status(f"{got / 1e6:.0f} of {total_bytes / 1e6:.0f} MB  "
                                     f"({time.perf_counter() - started:.0f} s)")
            proc_fig.canvas.flush_events()
            time.sleep(0.05)
    finally:
        close_processing_dialog(proc_fig)
        restore_cursors()
    if cancel_flag['cancelled']:
        print(f"The {what} download was cancelled; it continues in the background and will be there next time.")
        return False
    if outcome['error'] is not None:
        print(f"The {what} download failed: {outcome['error']}")
        return False
    return not ccf_files_missing(abc_cache, files)        # (the files asked for: not the default CCF ones, which a download of others never fetches)


def load_ccf_for_cells(abc_cache, cell_labels, shell_points=0):
    """(table, shell): `table` is a DataFrame aligned to `cell_labels` with float32 x_ccf/y_ccf/z_ccf and
    parcellation_index (NaN for a cell the CCF file lacks); `shell` is a (shell_points, 3) float32 array
    of x/y/z for a random sample of cells across the whole brain (empty when shell_points is 0). Reads
    the big CSV once, then a pickle cache next to the spatial cache. Blocking and slow the first time:
    call it from a worker thread."""
    csv_path = abc_cache.get_metadata_path(directory=CCF_DIRECTORY, file_name=CCF_COORDS_FILE_NAME)
    status = os.stat(csv_path)
    key = f"{csv_path}|{status.st_mtime_ns}|{status.st_size}"
    cache_path = CACHE_DIR / f"ccf_coords_v{CCF_CACHE_VERSION}_{hashlib.md5(key.encode('utf-8')).hexdigest()[:16]}.pkl"
    table = None
    if cache_path.exists():
        try:
            table = pd.read_pickle(cache_path)
        except Exception as e:
            print(f"CCF cache unreadable ({e}); re-reading the CSV.")
    if table is None:
        chunks = []
        for chunk in pd.read_csv(csv_path, index_col=0, dtype={'cell_label': str}, chunksize=SPATIAL_CSV_CHUNK_SIZE):
            chunks.append(chunk)
        table = pd.concat(chunks, copy=False)
        table = table.rename(columns={'x': 'x_ccf', 'y': 'y_ccf', 'z': 'z_ccf'})
        for column in ('x_ccf', 'y_ccf', 'z_ccf'):
            table[column] = table[column].astype(np.float32)
        table['parcellation_index'] = table['parcellation_index'].astype(np.float32)
        try:
            table.to_pickle(cache_path)
        except Exception as e:
            print(f"Could not save the CCF cache ({e}).")
    aligned = table.reindex(pd.Index(cell_labels))
    shell = np.empty((0, 3), dtype=np.float32)
    if shell_points:
        coords_all = table[['x_ccf', 'y_ccf', 'z_ccf']].to_numpy()
        keep = np.flatnonzero(np.isfinite(coords_all).all(axis=1))
        pick = np.random.default_rng(0).choice(keep, size=min(shell_points, len(keep)), replace=False)
        shell = coords_all[pick]
    return aligned, shell


def load_parcellation_names(abc_cache):
    """DataFrame indexed by parcellation_index with the CCF region names (division / structure /
    substructure acronyms), or None if the lookup can't be read."""
    try:
        path = abc_cache.get_metadata_path(directory=CCF_PARCELLATION_DIRECTORY,
                                           file_name=CCF_PARCELLATION_FILE_NAME)
        return pd.read_csv(path).set_index('parcellation_index')[['division', 'structure', 'substructure']]
    except Exception as e:
        print(f"Could not read the CCF region names ({e}).")
        return None


def brain_outline_segments(shell_xyz, bins=80):
    """The whole brain's silhouette seen along each axis, from a (n, 3) sample of cells: {'z': [(m, 2) arrays in
    (x, y)], 'x': [... in (y, z)], 'y': [... in (x, z)]}, each array one closed outline. Built from an occupancy
    image of the projected points (closed, filled and smoothed), then contoured."""
    from scipy import ndimage
    from matplotlib.figure import Figure
    shell = np.asarray(shell_xyz, dtype=float)
    outlines = {}
    for normal, (a, b) in {'z': (0, 1), 'x': (1, 2), 'y': (0, 2)}.items():
        histogram, edges_a, edges_b = np.histogram2d(shell[:, a], shell[:, b], bins=bins)
        occupied = ndimage.binary_fill_holes(ndimage.binary_closing(histogram > 0, iterations=2))
        smooth = ndimage.gaussian_filter(occupied.astype(float), 1.0)
        centers_a = (edges_a[:-1] + edges_a[1:]) / 2
        centers_b = (edges_b[:-1] + edges_b[1:]) / 2
        contour_axes = Figure().subplots()  # a throwaway Agg figure: only its contour paths are wanted
        contours = contour_axes.contour(centers_a, centers_b, smooth.T, levels=[0.5])
        outlines[normal] = [np.asarray(segment) for segment in contours.allsegs[0] if len(segment) > 3]
    return outlines


def read_nifti_volume(path):
    """The voxel array of a (gzipped) NIfTI-1 file as [i, j, k] (the file's own axis order), without needing
    nibabel: this only has to cope with the few plain integer/float volumes the atlas ships."""
    import gzip
    import struct
    with gzip.open(path, 'rb') as handle:
        raw = handle.read()
    header = raw[:348]
    dims = struct.unpack('<8h', header[40:56])
    shape = tuple(dims[1:1 + dims[0]])
    code, = struct.unpack('<h', header[70:72])
    dtypes = {2: '<u1', 4: '<i2', 8: '<i4', 16: '<f4', 64: '<f8', 256: '<i1', 512: '<u2', 768: '<u4'}
    if code not in dtypes:
        raise ValueError(f"unsupported NIfTI data type {code}")
    offset = int(struct.unpack('<f', header[108:112])[0])
    count = int(np.prod(shape))
    return np.frombuffer(raw, dtype=dtypes[code], offset=offset, count=count).reshape(shape[::-1]).T


def coarse_smoothed_mask(mask, step, sigma):
    """`mask` (bool, [i, j, k]) cropped to its own bounding box (so a small region costs little), reduced by `step`
    per axis with any-voxel-set kept (so a thin ventricle survives the reduction), padded by one voxel so surfaces
    close, then lightly blurred: (float32 array or None if the mask is empty, voxel_mm, origin_mm, pad), where
    origin_mm is the (3,) position in mm of the first unpadded voxel's centre and voxel `n` of the array sits at
    origin_mm + (n - pad) * voxel_mm on each axis."""
    from scipy import ndimage
    voxel_mm = step * CCF_VOLUME_VOXEL_MM
    pad = 1
    bounds = []
    for axis in range(3):
        others = tuple(a for a in range(3) if a != axis)
        present = np.flatnonzero(mask.any(axis=others))
        if not len(present):
            return None, voxel_mm, np.zeros(3), pad
        low = max(0, int(present[0]) - 2 * step) // step * step
        high = min(mask.shape[axis], int(present[-1]) + 1 + 2 * step)
        bounds.append((low, low + (high - low) // step * step))
    cropped = mask[bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1], bounds[2][0]:bounds[2][1]].astype(np.uint8)
    if step > 1:
        n0, n1, n2 = (n // step for n in cropped.shape)
        cropped = cropped.reshape(n0, step, n1, step, n2, step).max(axis=(1, 3, 5))
    padded = np.pad(cropped, pad).astype(np.float32)
    blurred = ndimage.gaussian_filter(padded, sigma) if sigma else padded
    origin_mm = np.array([(low + (step - 1) / 2) * CCF_VOLUME_VOXEL_MM for low, _high in bounds])
    return blurred, voxel_mm, origin_mm, pad


def slice_contour_lines(mask, spacing_mm, step=2, sigma=0.7):
    """Outline of a 3D label mask drawn as contour lines on slices every `spacing_mm` along all three axes: (x, y, z)
    float32 arrays with NaN between separate lines, in the 3D view's own axis order (left-right = CCF z, anterior-
    posterior = CCF x, dorsal-ventral = CCF y). Empty arrays if the mask is empty."""
    import contourpy
    empty = np.empty(0, dtype=np.float32)
    blurred, voxel_mm, origin_mm, pad = coarse_smoothed_mask(mask, step, sigma)
    if blurred is None:
        return empty, empty, empty
    coords = [np.arange(n) * voxel_mm + origin_mm[axis] - pad * voxel_mm for axis, n in enumerate(blurred.shape)]
    stride = max(1, int(round(spacing_mm / voxel_mm)))
    pieces = []
    for axis in range(3):
        other = [a for a in range(3) if a != axis]  # slice[b, c] with b < c
        first = int(np.ceil(origin_mm[axis] / spacing_mm - 1e-9))  # slices sit on multiples of spacing_mm in CCF mm
        for position in np.arange(first, 1e9, 1.0):
            mm = position * spacing_mm
            index = int(round((mm - coords[axis][0]) / voxel_mm))
            if index >= blurred.shape[axis]:
                break
            if index < 0:
                continue
            plane = np.take(blurred, index, axis=axis)
            if plane.max() < 0.5:
                continue
            generator = contourpy.contour_generator(x=coords[other[1]], y=coords[other[0]], z=plane)
            for line in generator.lines(0.5):
                points = np.empty((len(line), 3), dtype=np.float32)
                points[:, axis] = coords[axis][index]
                points[:, other[0]] = line[:, 1]
                points[:, other[1]] = line[:, 0]
                pieces.append(points)
                pieces.append(np.full((1, 3), np.nan, dtype=np.float32))
    if not pieces:
        return empty, empty, empty
    stacked = np.vstack(pieces)
    return stacked[:, 2], stacked[:, 0], stacked[:, 1]  # (ML, AP, DV) = (z_ccf, x_ccf, y_ccf)


THREED_ENVELOPE_MAX_GROUPS = 8         # envelopes are built for at most this many groups (the ones with most cells)
THREED_ENVELOPE_MIN_CELLS = 30         # a group with fewer cells has no meaningful envelope
# An envelope has two independent settings, both precomputed into the page and chosen with its own slider: how many of the
# group's cells lie inside the surface (the share), and how smooth the surface is (the width, in grid steps, of the blur
# applied to the cells' density before the surface is cut; a small width follows every bump, a large one rounds them off).
THREED_ENVELOPE_COVERAGES = (0.80, 0.88, 0.92, 0.95, 0.97, 0.99)   # shares of the cells inside, ascending (tightest first)
THREED_ENVELOPE_SMOOTHNESS = (0.4, 0.5, 0.6, 0.7, 1.0, 1.5)                 # blur widths, roughest first
THREED_ENVELOPE_DEFAULT_COVERAGE = 3   # (95%) indexes of the starting positions of the page's two sliders
THREED_ENVELOPE_DEFAULT_SMOOTHNESS = 1
THREED_ENVELOPE_OPACITY = 1.0          # where the page's envelope opacity slider starts (it goes up to 100%)


def group_envelope_meshes(points):
    """Surfaces enveloping a group of cells, for every smoothness and coverage of the settings above: a list over
    smoothness of lists over coverage of (vertices (n, 3) float32 in the points' own axis order, faces (m, 3) int32), or
    None when there is not enough to wrap. Each is an iso-surface of the cells' smoothed density, drawn at the level that
    keeps that share of the cells inside (so it hugs the cloud and ignores stragglers, unlike a convex hull, which a few
    outliers would stretch). A surface that cannot be made repeats a neighbour; the convex hull (a single surface)
    is the fallback where scikit-image is missing."""
    from scipy import ndimage
    points = np.asarray(points, dtype=float)
    points = points[np.isfinite(points).all(axis=1)]
    if len(points) < THREED_ENVELOPE_MIN_CELLS:
        return None
    if len(points) > 300000:   # more than this adds nothing to a smoothed shape
        points = points[np.random.default_rng(0).choice(len(points), 300000, replace=False)]
    low, high = points.min(axis=0), points.max(axis=0)
    cell = max(0.1, float(np.ptp(points, axis=0).max()) / 60)   # grid step, in the points' own units (mm)
    pad = int(np.ceil(3 * max(THREED_ENVELOPE_SMOOTHNESS))) + 1
    shape = np.ceil((high - low) / cell).astype(int) + 1 + 2 * pad
    index = np.floor((points - low) / cell).astype(int) + pad
    counts = np.zeros(shape, dtype=np.float32)
    np.add.at(counts, tuple(index.T), 1.0)
    grid = []
    try:
        from skimage import measure
        for sigma in THREED_ENVELOPE_SMOOTHNESS:
            density = ndimage.gaussian_filter(counts, sigma)
            row = []
            for share in THREED_ENVELOPE_COVERAGES:
                try:
                    level = float(np.percentile(density[tuple(index.T)], 100 * (1 - share)))
                    if level <= 0:
                        raise ValueError("empty density")
                    vertices, faces, _n, _v = measure.marching_cubes(density, level=level, spacing=(cell,) * 3)
                    if not len(faces):
                        raise ValueError("no surface")
                    row.append(((vertices + low - pad * cell).astype(np.float32), faces.astype(np.int32)))
                except Exception:
                    row.append(row[-1] if row else None)
            grid.append(row)
        known = [m for row in grid for m in row if m is not None]
        if known:   # (a surface that failed before any other in its row takes the first one that worked)
            return [[m if m is not None else known[0] for m in row] for row in grid]
    except Exception:
        pass
    try:   # fallback: the convex hull, a single surface
        from scipy.spatial import ConvexHull
        hull = ConvexHull(points)
        return [[(points[hull.vertices].astype(np.float32), _hull_faces(hull))]]
    except Exception:
        return None


def _hull_faces(hull):
    """Faces of a scipy ConvexHull re-indexed to its own vertex list (hull.vertices order), as an int32 (m, 3) array."""
    remap = {int(v): i for i, v in enumerate(hull.vertices)}
    return np.array([[remap[int(v)] for v in simplex] for simplex in hull.simplices], dtype=np.int32)


def surface_mesh(mask, step, sigma=1.0):
    """The surface of a 3D label mask as (vertices (n, 3) float32 in the 3D view's axis order, faces (m, 3) int32),
    by marching cubes at `step` times the volume's own voxel size; None if the mask is empty or scikit-image is
    not installed."""
    try:
        from skimage import measure
    except ImportError:
        return None
    blurred, voxel_mm, origin_mm, pad = coarse_smoothed_mask(mask, step, sigma)
    if blurred is None or blurred.max() < 0.5:
        return None
    vertices, faces, _normals, _values = measure.marching_cubes(blurred, level=0.5, spacing=(voxel_mm,) * 3)
    vertices = vertices - pad * voxel_mm + origin_mm
    return vertices[:, [2, 0, 1]].astype(np.float32), faces.astype(np.int32)  # (ML, AP, DV); a cyclic swap keeps winding


def region_edge_geometry(volume, labels, slice_mm, line_step, surface_step):
    """{'lines': (x, y, z), 'surface': (vertices, faces) or None} for the voxels of `volume` whose label is in `labels`."""
    labels = np.asarray(list(labels), dtype=np.int64)
    lookup = np.zeros(int(volume.max()) + 1, dtype=bool)  # a table lookup, far faster than np.isin over ~77 M voxels
    lookup[labels[(labels >= 0) & (labels < len(lookup))]] = True
    mask = lookup[volume]
    return {'lines': slice_contour_lines(mask, slice_mm, step=line_step),
            'surface': surface_mesh(mask, surface_step)}


def choose_cell_regions(structures):
    """The CCF structure names holding most of the shown cells, biggest first: at most THREED_EDGE_REGION_MAX of them,
    each with at least THREED_EDGE_REGION_MIN_SHARE of the cells (but always the biggest one)."""
    names = pd.Series(structures).dropna().astype(str)
    names = names[~names.isin(['unassigned', 'brain-unassigned', 'fiber tracts-unassigned', '?'])]
    if names.empty:
        return []
    counts = names.value_counts()
    keep = [name for name, count in counts.items() if count / len(names) >= THREED_EDGE_REGION_MIN_SHARE]
    return (keep or [counts.index[0]])[:THREED_EDGE_REGION_MAX]


def ccf_edges_missing(region_names, want_surfaces=True):
    """(store, missing_keys): the edges cache's contents and which of 'brain', 'ventricles' and ('region', name) it
    still lacks (or has only without surfaces, back when scikit-image was not installed). Nothing is computed."""
    cache_path = CACHE_DIR / f"ccf_edges_v{CCF_EDGES_CACHE_VERSION}.pkl"
    store = {}
    if cache_path.exists():
        try:
            store = pd.read_pickle(cache_path)
        except Exception as e:
            print(f"CCF edges cache unreadable ({e}); recomputing.")
    try:
        import skimage  # noqa: F401
        have_skimage = True
    except ImportError:
        have_skimage = False
    wanted = (['brain', 'ventricles'] + [('white', name) for name in CCF_WHITE_MATTER]
              + [('region', name) for name in region_names])
    missing = [key for key in wanted
               if store.get(key) is None or (want_surfaces and have_skimage and not store[key]['has_skimage'])]
    return store, missing, have_skimage


def load_ccf_edges(abc_cache, parcel_names, region_names, want_surfaces=True):
    """Edge geometry for the browser 3D view, from the cache file where already computed (the volume itself is only
    read, and downloaded if need be, when something is missing). Returns {'brain': g, 'ventricles': g, 'regions':
    {name: g}} with each g as in region_edge_geometry; a region name with no voxels is left out. Blocking and slow
    the first time (about half a minute): call it from a worker thread."""
    store, missing, have_skimage = ccf_edges_missing(region_names, want_surfaces)
    if missing:
        path = abc_cache.get_data_path(directory=CCF_VOLUME_DIRECTORY, file_name=CCF_VOLUME_FILE_NAME)
        volume = read_nifti_volume(path)
        for key in missing:
            if key == 'brain':
                geometry = region_edge_geometry(volume, np.arange(1, int(volume.max()) + 1), THREED_EDGE_SLICE_MM, 2, 4)
            elif key == 'ventricles':
                labels = parcel_names.index[parcel_names['division'].isin(CCF_VENTRICLE_DIVISIONS)]
                geometry = region_edge_geometry(volume, labels, THREED_EDGE_SLICE_MM, 1, 2)
            elif key[0] == 'white':
                column, acronyms, _color = CCF_WHITE_MATTER[key[1]]
                labels = parcel_names.index[parcel_names[column].isin(acronyms)]
                geometry = region_edge_geometry(volume, labels, THREED_EDGE_SLICE_MM, 1, 2)
            else:
                labels = parcel_names.index[parcel_names['structure'] == key[1]]
                geometry = region_edge_geometry(volume, labels, THREED_EDGE_REGION_SLICE_MM, 1, 2)
            geometry['has_skimage'] = have_skimage
            store[key] = geometry
        try:
            pd.to_pickle(store, CACHE_DIR / f"ccf_edges_v{CCF_EDGES_CACHE_VERSION}.pkl")
        except Exception as e:
            print(f"Could not save the CCF edges cache ({e}).")
    return {'brain': store['brain'], 'ventricles': store['ventricles'],
            'white_matter': {name: store[('white', name)] for name in CCF_WHITE_MATTER
                             if len(store[('white', name)]['lines'][0])},
            'regions': {name: store[('region', name)] for name in region_names if len(store[('region', name)]['lines'][0])}}


def fair_keep_counts(sizes, budget):
    """How many points each group may keep so the total stays within `budget`: small groups keep everything, and
    only the groups bigger than their share are thinned (a 500-cell group is not cut because another has 200k)."""
    keep, remaining_budget = {}, budget
    order = sorted(range(len(sizes)), key=lambda i: sizes[i])
    for position, index in enumerate(order):
        keep[index] = min(sizes[index], int(remaining_budget / (len(order) - position)))
        remaining_budget -= keep[index]
    return keep


def unit_sphere_mesh(subdivisions):
    """(vertices, faces) of a unit sphere built by subdividing an icosahedron: 12 vertices / 20 faces for 0
    subdivisions, 42 / 80 for 1. Vertices are float32 (n, 3); faces int32 (m, 3), wound outward."""
    t = (1 + 5 ** 0.5) / 2
    vertices = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
                (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2),
             (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5),
             (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    vertices = [np.asarray(v, dtype=float) / np.linalg.norm(v) for v in vertices]
    for _ in range(subdivisions):
        midpoints, refined = {}, []

        def midpoint(a, b):
            key = (min(a, b), max(a, b))
            if key not in midpoints:
                m = vertices[a] + vertices[b]
                vertices.append(m / np.linalg.norm(m))
                midpoints[key] = len(vertices) - 1
            return midpoints[key]

        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = refined
    return np.asarray(vertices, dtype=np.float32), np.asarray(faces, dtype=np.int32)


def sphere_mesh_at(centers, radius, unit_vertices, unit_faces):
    """One mesh holding a sphere of `radius` at each row of `centers` (n, 3): (x, y, z, i, j, k) for plotly's Mesh3d."""
    centers = np.asarray(centers, dtype=np.float32)
    per_sphere = len(unit_vertices)
    vertices = (centers[:, None, :] + unit_vertices[None, :, :] * np.float32(radius)).reshape(-1, 3)
    offsets = (np.arange(len(centers), dtype=np.int32) * per_sphere)[:, None, None]
    faces = (unit_faces[None, :, :] + offsets).reshape(-1, 3)
    return vertices[:, 0], vertices[:, 1], vertices[:, 2], faces[:, 0], faces[:, 1], faces[:, 2]


def open_3d_group_view(groups, shell_xyz, outline, axis_labels, title, invert_z=True, max_points=THREED_MAX_POINTS):
    """A rotatable 3D scatter window on a dark panel (drag with the left mouse button to rotate, right button or
    wheel to zoom). `groups`: list of {'label', 'color', 'xyz'} with xyz an (n, 3) array in plotting order
    (x, y, z); `shell_xyz`: optional whole-brain sample drawn as a faint cloud; `outline`: optional
    brain_outline_segments() result drawn on the three back walls. Click a legend entry to hide or show that
    item ('wall projections' is off until clicked and drops each group's shadow onto the walls). Buttons along the
    bottom jump to standard views and change the dot size. Groups are thinned fairly to `max_points` in total.
    Returns the figure, already shown without blocking."""
    import mpl_toolkits.mplot3d  # noqa: F401  (registers the '3d' projection)
    from matplotlib.lines import Line2D
    from matplotlib.widgets import Button
    background = THREED_BACKGROUND
    rng = np.random.default_rng(0)
    cleaned = []
    for group in groups:
        points = np.asarray(group['xyz'], dtype=float)
        cleaned.append(points[np.isfinite(points).all(axis=1)])
    total = sum(len(points) for points in cleaned)
    keep = fair_keep_counts([len(points) for points in cleaned], max_points)

    fig = plt.figure(figsize=(10, 9), facecolor=background)
    ax3 = fig.add_subplot(111, projection='3d', facecolor=background)
    items = []  # [(label, color, [artists])] in legend order
    size_scale = {'value': 1.0}
    group_artists = []  # (artist, base size)

    cloud = None
    if shell_xyz is not None and len(shell_xyz):
        cloud_points = np.asarray(shell_xyz)[:THREED_CLOUD_POINTS]
        cloud = ax3.scatter(cloud_points[:, 0], cloud_points[:, 1], cloud_points[:, 2], s=2, c='#aeb4c0',
                            alpha=0.16, linewidths=0, depthshade=False)
    # Limits from everything drawn, so the outline walls sit exactly at the data's own extent.
    drawn = [np.asarray(shell_xyz)] if shell_xyz is not None and len(shell_xyz) else []
    drawn += [points for points in cleaned if len(points)]
    everything = np.vstack(drawn) if drawn else np.zeros((1, 3))
    low, high = everything.min(axis=0), everything.max(axis=0)
    spans = np.maximum(high - low, 1e-6)
    wall_x = low[0]
    wall_y = high[1]
    wall_z = high[2] if invert_z else low[2]  # the floor is the ventral side when z is inverted

    outline_artists = []
    if outline:
        for normal, segments in outline.items():
            for segment in segments:
                if normal == 'z':
                    line = ax3.plot(segment[:, 0], segment[:, 1], zs=wall_z, zdir='z')
                elif normal == 'x':
                    line = ax3.plot(segment[:, 0], segment[:, 1], zs=wall_x, zdir='x')
                else:
                    line = ax3.plot(segment[:, 0], segment[:, 1], zs=wall_y, zdir='y')
                line[0].set_color('#8a93a6')
                line[0].set_linewidth(1.0)
                outline_artists.append(line[0])
    if outline_artists:
        items.append(('brain outline', '#8a93a6', outline_artists))
    if cloud is not None:
        items.append(('whole brain (cloud)', '#aeb4c0', [cloud]))

    shadow_artists = []
    shown = 0
    for index, (group, points) in enumerate(zip(groups, cleaned)):
        if keep[index] < len(points):
            points = points[rng.choice(len(points), size=max(1, keep[index]), replace=False)]
        shown += len(points)
        artist = ax3.scatter(points[:, 0], points[:, 1], points[:, 2], s=7, color=group['color'], alpha=0.9,
                             linewidths=0, depthshade=True)
        group_artists.append((artist, 7))
        items.append((group['label'], group['color'], [artist]))
        if len(points):
            for xs, ys, zs, zdir in ((points[:, 0], points[:, 1], wall_z, 'z'), (points[:, 0], points[:, 2], wall_y, 'y'),
                                     (points[:, 1], points[:, 2], wall_x, 'x')):
                shadow = ax3.scatter(xs, ys, zs=zs, zdir=zdir, s=2, color=group['color'], alpha=0.25, linewidths=0,
                                     depthshade=False)
                shadow.set_visible(False)
                shadow_artists.append(shadow)
    if shadow_artists:
        items.append(('wall projections', '0.85', shadow_artists))

    light = '0.85'
    for axis, label in zip((ax3.xaxis, ax3.yaxis, ax3.zaxis), axis_labels):
        axis.set_pane_color((0.07, 0.08, 0.10, 1.0))
        axis._axinfo['grid']['color'] = (1.0, 1.0, 1.0, 0.10)
        axis.set_tick_params(colors='0.65')
        axis.label.set_color(light)
        axis.set_label_text(label)
    if invert_z:
        ax3.invert_zaxis()
    ax3.set_box_aspect(tuple(spans))  # equal scale on all three axes: no distortion of the brain
    ax3.set_title(title if shown >= total else f"{title}  (showing {shown:,} of {total:,} cells)", color='white')

    handles = [Line2D([0], [0], marker='o', linestyle='', color=color, markersize=7) for _l, color, _a in items]
    legend = ax3.legend(handles, [label for label, _c, _a in items], loc='upper left', title='click to hide/show',
                        facecolor='#1d2029', edgecolor='#3a3f4d', labelcolor=light)
    legend.get_title().set_color(light)
    pick_targets = {}
    for index, (text, handle) in enumerate(zip(legend.get_texts(), legend.legend_handles)):
        for part in (text, handle):
            part.set_picker(6)
            pick_targets[part] = (index, text, handle)
    for index, (label, _color, artists) in enumerate(items):
        if label == 'wall projections':  # starts off, so its entry starts dimmed
            legend.get_texts()[index].set_alpha(0.35)
            legend.legend_handles[index].set_alpha(0.35)

    def on_pick(event):
        entry = pick_targets.get(event.artist)
        if entry is None:
            return
        index, text, handle = entry
        artists = items[index][2]
        visible = not artists[0].get_visible()
        for artist in artists:
            artist.set_visible(visible)
        text.set_alpha(1.0 if visible else 0.35)
        handle.set_alpha(1.0 if visible else 0.35)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('pick_event', on_pick)

    # Standard viewpoints and dot size; the axes are x = left/right, y = front/back, z = up/down.
    def resize_dots(factor):
        size_scale['value'] *= factor
        for artist, base in group_artists:
            artist.set_sizes([base * size_scale['value']])
        fig.canvas.draw_idle()

    actions = (('Front', lambda: ax3.view_init(elev=0, azim=-90)), ('Side', lambda: ax3.view_init(elev=0, azim=0)),
               ('Top', lambda: ax3.view_init(elev=90, azim=-90)), ('Default', lambda: ax3.view_init(elev=30, azim=-60)),
               ('Dots -', lambda: resize_dots(0.8)), ('Dots +', lambda: resize_dots(1.25)))
    fig._view_buttons = []
    for i, (name, action) in enumerate(actions):
        button_ax = fig.add_axes([0.03 + i * 0.105, 0.015, 0.095, 0.04])
        button = Button(button_ax, name, color='#2a2e3a', hovercolor='#454b5c')
        button.label.set_color('white')
        button.on_clicked(lambda event, a=action: (a(), fig.canvas.draw_idle()))
        fig._view_buttons.append(button)  # kept alive: a Button with no reference stops responding
    try:
        fig.canvas.manager.set_window_title(title)
    except Exception:
        pass
    fig.show()
    return fig


def write_3d_group_html(groups, shell_xyz, outline, axis_labels, title, path, invert_z=True,
                        max_points=THREED_WEB_MAX_POINTS, edges=None):
    """The same view as open_3d_group_view as one self-contained interactive HTML file (plotly, WebGL), written to
    `path`: true depth handling, smooth with hundreds of thousands of points, hover text per cell where a group
    gives 'hover' strings, and the legend toggles each group. Raises ImportError if plotly is not installed."""
    import plotly.graph_objects as go
    from matplotlib.colors import to_hex
    rng = np.random.default_rng(0)
    cleaned = []
    for group in groups:
        points = np.asarray(group['xyz'], dtype=float)
        finite = np.isfinite(points).all(axis=1)
        hover = group.get('hover')
        cleaned.append((points[finite], None if hover is None else np.asarray(hover, dtype=object)[finite]))
    total = sum(len(points) for points, _h in cleaned)
    keep = fair_keep_counts([len(points) for points, _h in cleaned], max_points)

    figure = go.Figure()
    # The page opens on Surfaces when there are shells (scikit-image was available), else Lines, else the cloud.
    default_mode = 'cloud'
    if edges:
        geometries = [edges['brain'], edges['ventricles'], *edges['white_matter'].values(), *edges['regions'].values()]
        if any(g['surface'] is not None for g in geometries):
            default_mode = 'surfaces'
        elif any(len(g['lines'][0]) for g in geometries):
            default_mode = 'lines'
    roles = {'cloud': [], 'walls': [], 'lines': [], 'surfaces': []}  # trace indices the Background radio group switches
    background_base = {}  # background trace index -> its own opacity at 100% on the Background opacity slider
    legend_off = []       # background traces that start switched off (shown by a click on their legend entry)
    if outline:
        xs, ys, zs = [], [], []
        low_high = np.asarray(shell_xyz) if shell_xyz is not None and len(shell_xyz) else np.vstack([p for p, _h in cleaned])
        wall_x, wall_y = low_high[:, 0].min(), low_high[:, 1].max()
        wall_z = low_high[:, 2].max() if invert_z else low_high[:, 2].min()
        for normal, segments in outline.items():
            for segment in segments:
                n = len(segment)
                if normal == 'z':
                    cx, cy, cz = segment[:, 0], segment[:, 1], np.full(n, wall_z)
                elif normal == 'x':
                    cx, cy, cz = np.full(n, wall_x), segment[:, 0], segment[:, 1]
                else:
                    cx, cy, cz = segment[:, 0], np.full(n, wall_y), segment[:, 1]
                xs += list(cx) + [None]
                ys += list(cy) + [None]
                zs += list(cz) + [None]
        roles['walls'].append(len(figure.data))
        background_base[len(figure.data)] = 1.0
        figure.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines', line=dict(color='#8a93a6', width=2),
                                      name='brain outline', hoverinfo='skip', visible=default_mode == 'cloud'))
    if shell_xyz is not None and len(shell_xyz):
        cloud = np.asarray(shell_xyz)[:THREED_CLOUD_POINTS * 2]
        roles['cloud'].append(len(figure.data))
        background_base[len(figure.data)] = 1.0
        figure.add_trace(go.Scatter3d(
            x=np.round(cloud[:, 0], 3), y=np.round(cloud[:, 1], 3), z=np.round(cloud[:, 2], 3), mode='markers',
            marker=dict(size=1.2, color='#aeb4c0', opacity=0.12), name='whole brain (cloud)', hoverinfo='skip',
            visible=default_mode == 'cloud'))
    if edges:
        # Added before the cells (so a translucent shell is drawn first), hidden until the radio group picks them.
        def break_at_nan(values):  # plotly breaks a line at None, not NaN
            return [None if np.isnan(v) else float(v) for v in np.round(np.asarray(values, dtype=float), 3)]

        edge_sets = [('brain', THREED_EDGE_BRAIN_COLOR, THREED_WEB_BRAIN_SURFACE_DEFAULT, edges['brain']),
                     ('ventricles', THREED_EDGE_VENTRICLE_COLOR, 0.35, edges['ventricles'])]
        for name, geometry in edges['white_matter'].items():
            edge_sets.append((name, CCF_WHITE_MATTER[name][2], 0.35, geometry))
        for region_index, (name, geometry) in enumerate(edges['regions'].items()):
            edge_sets.append((name, THREED_EDGE_REGION_COLORS[region_index % len(THREED_EDGE_REGION_COLORS)], 0.30,
                              geometry))
        # Only the brain and the ventricles start switched on; every other object (white-matter tracts, the cells'
        # regions) starts off, and a click on its legend entry shows it.
        off_names = {name for name, _c, _o, _g in edge_sets if name not in ('brain', 'ventricles')}
        for name, color, shell_opacity, geometry in edge_sets:
            xs, ys, zs = geometry['lines']
            if len(xs):
                line_opacity = 0.5 if name == 'brain' else 1.0  # the dense brain mesh is drawn faded
                roles['lines'].append(len(figure.data))
                background_base[len(figure.data)] = line_opacity
                figure.add_trace(go.Scatter3d(
                    x=break_at_nan(xs), y=break_at_nan(ys), z=break_at_nan(zs), mode='lines',
                    line=dict(color=color, width=2), opacity=line_opacity,
                    name=f"{name} edges", hoverinfo='skip',
                    visible=('legendonly' if name in off_names else True) if default_mode == 'lines' else False))
                if name in off_names:
                    legend_off.append(len(figure.data) - 1)
        for name, color, shell_opacity, geometry in edge_sets:
            if geometry['surface'] is not None:
                vertices, faces = geometry['surface']
                roles['surfaces'].append(len(figure.data))
                background_base[len(figure.data)] = shell_opacity
                figure.add_trace(go.Mesh3d(
                    x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2], i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                    color=color, opacity=shell_opacity, flatshading=False,
                    lighting=dict(ambient=0.6, diffuse=0.7, specular=0.1, roughness=0.8),
                    name=f"{name} surface", showlegend=True, hoverinfo='skip',
                    visible=('legendonly' if name in off_names else True) if default_mode == 'surfaces' else False))
                if name in off_names:
                    legend_off.append(len(figure.data) - 1)
    # An envelope surface for each of the biggest cell groups, in the legend like the other surfaces but off until its
    # entry is clicked (not part of the Background radio group, so switching that never shows or hides them).
    import re
    by_size = sorted(range(len(cleaned)), key=lambda i: -len(cleaned[i][0]))[:THREED_ENVELOPE_MAX_GROUPS]
    def add_legend_spacer(count):
        """An empty legend row (an invisible stand-in trace with a blank name; the page hides it): one line of space between
        the legend's groups, which plotly has no way to ask for. `count` makes each spacer's name different."""
        figure.add_trace(go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', name='​' * count,
                                      marker=dict(size=1, color='rgba(0,0,0,0)'), hoverinfo='skip', showlegend=True))

    envelope_levels, envelope_colors, envelope_names = [], [], []   # per envelope: its [smoothness][coverage] surfaces etc.
    envelope_pieces = []   # the default surfaces, merged below into ONE mesh
    for group_index in sorted(by_size):
        levels = group_envelope_meshes(cleaned[group_index][0])
        if levels is None:
            continue
        smooth_row = levels[min(THREED_ENVELOPE_DEFAULT_SMOOTHNESS, len(levels) - 1)]
        envelope_pieces.append(smooth_row[min(THREED_ENVELOPE_DEFAULT_COVERAGE, len(smooth_row) - 1)])
        try:
            envelope_color = to_hex(groups[group_index]['color'])
        except Exception:
            envelope_color = '#ff0000'
        envelope_colors.append([int(envelope_color[n:n + 2], 16) for n in (1, 3, 5)])
        envelope_names.append(re.sub(r'\s*\([\d,]+\)\s*$', '', str(groups[group_index]['label'])) + ' envelope')
        envelope_levels.append([[{'x': np.round(v[:, 0], 2).tolist(), 'y': np.round(v[:, 1], 2).tolist(),
                                  'z': np.round(v[:, 2], 2).tolist(), 'i': f[:, 0].tolist(), 'j': f[:, 1].tolist(),
                                  'k': f[:, 2].tolist()} for v, f in row] for row in levels])
    envelope_mesh = None       # trace index of the merged mesh; its triangles are listed by THREED_WEB_OPACITY_SCRIPT
    envelope_stand_ins = []    # trace indices of the legend entries (empty stand-in traces, the page script drives them)
    if envelope_pieces:
        # All the envelopes are ONE mesh (a colour per vertex): plotly draws separate traces one after another whatever
        # their depth, so one envelope would always paint over another; in a single mesh the page script can list
        # every triangle of every switched-on envelope far to near. The mesh starts with no triangles (every envelope
        # starts off); the legend entries are empty traces that the script reads clicks for.
        offsets = np.cumsum([0] + [len(v) for v, _f in envelope_pieces])[:-1]
        envelope_mesh = len(figure.data)
        figure.add_trace(go.Mesh3d(
            x=np.round(np.concatenate([v[:, 0] for v, _f in envelope_pieces]), 2),
            y=np.round(np.concatenate([v[:, 1] for v, _f in envelope_pieces]), 2),
            z=np.round(np.concatenate([v[:, 2] for v, _f in envelope_pieces]), 2),
            i=np.zeros(0, dtype=np.int32), j=np.zeros(0, dtype=np.int32), k=np.zeros(0, dtype=np.int32),
            vertexcolor=np.concatenate([np.tile(color, (len(v), 1)) for (v, _f), color in
                                        zip(envelope_pieces, envelope_colors)]).astype(int),
            opacity=THREED_ENVELOPE_OPACITY, flatshading=False,
            lighting=dict(ambient=0.6, diffuse=0.7, specular=0.1, roughness=0.8),
            name='envelopes', showlegend=False, hoverinfo='skip'))
        if any(shell_traces_added for shell_traces_added in (roles['surfaces'], roles['lines'])):
            add_legend_spacer(1)   # between the background objects and the envelopes
        for name, color in zip(envelope_names, envelope_colors):
            envelope_stand_ins.append(len(figure.data))
            figure.add_trace(go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', name=name,
                                          marker=dict(size=8, color='#%02x%02x%02x' % tuple(color)), hoverinfo='skip',
                                          showlegend=True, visible='legendonly'))
    if roles['surfaces'] or roles['lines'] or envelope_mesh is not None:
        add_legend_spacer(2)   # between the objects above and the cell groups
    planned = sum(keep.values())
    # Hover text on a mesh is per vertex, so it is repeated for every vertex of every sphere: with it, smoother
    # (more-vertex) spheres would multiply the file size, so they are only used when no group carries hover text.
    has_hover = any(hover is not None for _p, hover in cleaned)
    if planned <= THREED_WEB_SMOOTH_SPHERE_MAX_POINTS and not has_hover:
        sphere_subdivisions = 1
    elif planned <= THREED_WEB_SPHERE_MAX_POINTS:
        sphere_subdivisions = 0
    else:
        sphere_subdivisions = None  # too many for a mesh each: outlined flat markers instead
    if sphere_subdivisions is not None:
        unit_vertices, unit_faces = unit_sphere_mesh(sphere_subdivisions)
        everything = np.vstack(([np.asarray(shell_xyz)] if shell_xyz is not None and len(shell_xyz) else [])
                               + [points for points, _h in cleaned if len(points)])
        sphere_radius = float(np.ptp(everything, axis=0).max()) * THREED_WEB_SPHERE_RADIUS_FRACTION
    shown = 0
    group_traces = []  # indices of the cells' own traces: what the opacity slider adjusts (not outline or cloud)
    mesh_groups = []   # (label, color, points, hover text) per group, for the single combined sphere mesh
    for index, (group, (points, hover)) in enumerate(zip(groups, cleaned)):
        if keep[index] < len(points):
            chosen = rng.choice(len(points), size=max(1, keep[index]), replace=False)
            points = points[chosen]
            hover = None if hover is None else hover[chosen]
        shown += len(points)
        text = hover if hover is not None and len(points) <= 50000 else None
        try:
            color = to_hex(group['color'])
        except Exception:
            color = 'red'
        if sphere_subdivisions is not None:
            mesh_groups.append((group['label'], color, points, text))
            continue
        group_traces.append(len(figure.data))
        figure.add_trace(go.Scatter3d(
            x=np.round(points[:, 0], 3), y=np.round(points[:, 1], 3), z=np.round(points[:, 2], 3), mode='markers',
            marker=dict(size=2.6, color=color, opacity=1.0, line=dict(color=THREED_BACKGROUND, width=1)),
            name=group['label'],
            text=None if text is None else list(text), hoverinfo='name+text' if text is not None else 'name'))
    post_script = None
    if mesh_groups:
        # One mesh for every group (see THREED_WEB_MESH_SCRIPT for why); a group's colour rides on a per-vertex
        # intensity through a stepped colour scale, so the spheres stay flat-coloured and lit.
        pieces, vertex_starts, vertex_total = [], [], 0
        for group_index, (label, color, points, text) in enumerate(mesh_groups):
            xs, ys, zs, face_i, face_j, face_k = sphere_mesh_at(points, sphere_radius, unit_vertices, unit_faces)
            vertex_starts.append(vertex_total)
            if text is not None:
                hover_text = np.char.add(label + '<br>', np.asarray(text, dtype=str)).astype(object)
            else:
                hover_text = np.full(len(points), label, dtype=object)
            pieces.append((xs, ys, zs, face_i + vertex_total, face_j + vertex_total, face_k + vertex_total,
                           np.full(len(xs), group_index, dtype=np.int8),
                           np.repeat(hover_text, len(unit_vertices))))
            vertex_total += len(xs)
        n_groups = len(mesh_groups)
        colorscale = []
        for group_index, (_label, color, _points, _text) in enumerate(mesh_groups):
            colorscale += [[group_index / n_groups, color], [(group_index + 1) / n_groups, color]]
        mesh_index = len(figure.data)
        figure.add_trace(go.Mesh3d(
            x=np.concatenate([p[0] for p in pieces]), y=np.concatenate([p[1] for p in pieces]),
            z=np.concatenate([p[2] for p in pieces]), i=np.concatenate([p[3] for p in pieces]),
            j=np.concatenate([p[4] for p in pieces]), k=np.concatenate([p[5] for p in pieces]),
            intensity=np.concatenate([p[6] for p in pieces]), intensitymode='vertex', colorscale=colorscale,
            cmin=-0.5, cmax=n_groups - 0.5, showscale=False, opacity=1.0, flatshading=False,
            lighting=dict(ambient=0.45, diffuse=0.9, specular=0.35, roughness=0.5, fresnel=0.15),
            name='cells', showlegend=False, text=list(np.concatenate([p[7] for p in pieces])), hoverinfo='text'))
        group_traces.append(mesh_index)
        legend_indices = []
        for label, color, _points, _text in mesh_groups:  # empty stand-in traces: the legend entries the script drives
            legend_indices.append(len(figure.data))
            figure.add_trace(go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', name=label,
                                          marker=dict(size=8, color=color), hoverinfo='skip', showlegend=True))
        post_script = (THREED_WEB_MESH_SCRIPT.replace('__MESH_INDEX__', str(mesh_index))
                       .replace('__LEGEND_INDICES__', json.dumps(legend_indices))
                       .replace('__VERTEX_STARTS__', json.dumps(vertex_starts)))
    scene = dict(aspectmode='data', xaxis_title=axis_labels[0], yaxis_title=axis_labels[1], zaxis_title=axis_labels[2])
    if invert_z:
        scene['zaxis'] = dict(autorange='reversed', title=axis_labels[2])
    brain_surface_index = (roles['surfaces'][0] if edges and edges['brain']['surface'] is not None and roles['surfaces']
                           else None)  # the brain's own surface is always the first one added
    # The opacity sliders (THREED_WEB_OPACITY_SCRIPT) float over the lower-right corner of the plot, which is empty, so
    # they take nothing from the height available to the 3D scene.
    row = THREED_WEB_SLIDER_ROW_PX
    if background_base:
        # Background traces open at the Background slider's default fraction of their own opacity.
        for index in background_base:
            figure.data[index].opacity = round(THREED_WEB_BACKGROUND_DEFAULT * background_base[index], 3)
    slider_rows = (1 if group_traces else 0) + (1 if background_base else 0) + (1 if brain_surface_index is not None else 0)
    bottom_margin = 10
    figure.update_layout(
        template='plotly_dark', scene=scene, margin=dict(l=0, r=0, t=50, b=bottom_margin),
        title=title if shown >= total else f"{title}  (showing {shown:,} of {total:,} cells)",
        legend=dict(itemsizing='constant'))
    scripts = [THREED_WEB_LEGEND_ROUTER_SCRIPT] + ([post_script] if post_script else [])
    if roles['lines'] or roles['surfaces']:
        scripts.append(THREED_WEB_EDGE_SCRIPT.replace('__ROLES__', json.dumps(roles)).replace('__OFF__', json.dumps(legend_off)).replace('__DEFAULT__', default_mode))
    if slider_rows or envelope_mesh is not None:
        opacity_config = {
            'cells': list(group_traces),
            'envelopes': ({'mesh': envelope_mesh, 'standIns': envelope_stand_ins, 'colors': envelope_colors,
                           'base': THREED_ENVELOPE_OPACITY, 'levels': envelope_levels,
                           'coverages': [round(100 * c) for c in THREED_ENVELOPE_COVERAGES][:len(envelope_levels[0][0])],
                           'smoothSteps': len(envelope_levels[0]),
                           'defaultCoverage': min(THREED_ENVELOPE_DEFAULT_COVERAGE, len(envelope_levels[0][0]) - 1),
                           'defaultSmooth': min(THREED_ENVELOPE_DEFAULT_SMOOTHNESS, len(envelope_levels[0]) - 1)}
                          if envelope_mesh is not None else None),
            'background': ({'indices': list(background_base), 'bases': [background_base[i] for i in background_base],
                            'value': THREED_WEB_BACKGROUND_DEFAULT} if background_base else None),
            'brain': ({'index': brain_surface_index, 'value': background_base[brain_surface_index]}
                      if brain_surface_index is not None else None)}
        scripts.append(THREED_WEB_OPACITY_SCRIPT.replace('__CONFIG__', json.dumps(opacity_config))
                       .replace('__ROW__', str(row)))
    scripts.append(THREED_WEB_SLICE_SCRIPT.replace('__MESH_INDEX__', str(mesh_index if mesh_groups else -1)))
    scripts.append(THREED_WEB_SPIN_SCRIPT)
    scripts.append(THREED_WEB_LEGEND_MARK_SCRIPT)
    figure.write_html(str(path), include_plotlyjs=True, full_html=True, config={'displaylogo': False},
                      post_script='\n'.join(scripts) or None)
    return Path(path)




def is_data_file_cached(abc_cache, directory, file_name):
    """Best-effort check for whether `directory`/`file_name` (as would be
    passed to abc_cache.get_data_path) is already present in the local
    cache — unlike get_data_path itself, this never triggers a download for
    a file that isn't. Returns True/False, or None if the check itself
    failed (e.g. abc_atlas_access's internal cache/manifest shape changed
    in some future version) — callers should treat None as "unknown,
    assume not cached" rather than erroring, since there's no public
    abc_atlas_access API for this as of this writing; reaching into
    AbcProjectCache's underlying S3CloudCache/LocalCache and its manifest
    is the only way to answer this without paying for the download it's
    trying to avoid."""
    try:
        file_attributes = abc_cache.cache._manifest.get_file_attributes(
            directory=directory, file_name=file_name,
        )
        return file_complete(abc_cache, file_attributes)          # (a download cut short is there but not whole: see file_complete)
    except Exception:
        return None






def set_wait_cursor_on_open_figures():
    """Show the busy cursor over every open pyplot window. Returns a function
    that puts each window's previous cursor back, which may be a resize or
    other custom cursor rather than the arrow. Works on the Tk widgets
    directly, so a window whose own code tracks its cursor state isn't
    left out of sync. Windows that aren't Tk are skipped."""
    saved = []
    # Gcf rather than plt.figure(num), which would also change pyplot's
    # current figure.
    from matplotlib._pylab_helpers import Gcf
    for manager in Gcf.get_all_fig_managers():
        try:
            widget = manager.canvas.get_tk_widget()
            saved.append((widget, widget.cget('cursor')))
            widget.configure(cursor='watch')
            widget.update_idletasks()
        except Exception:
            pass

    def restore():
        for widget, cursor in saved:
            try:
                widget.configure(cursor=cursor)
            except Exception:
                pass  # window closed during the load
    return restore


LOAD_TIMINGS_FILE = CACHE_DIR / 'load_timings.json'


def load_recorded_duration(key):
    """Seconds the last successful run of the named load took on this
    machine, or None. Lets a load whose progress can't be measured (opening
    a huge h5ad reads all of obs with no hook to report on) show an
    estimated bar instead of an empty one."""
    try:
        return float(json.loads(LOAD_TIMINGS_FILE.read_text(encoding='utf-8')).get(key))
    except Exception:
        return None


def record_load_duration(key, seconds):
    try:
        try:
            timings = json.loads(LOAD_TIMINGS_FILE.read_text(encoding='utf-8'))
        except Exception:
            timings = {}
        timings[key] = round(seconds, 2)
        LOAD_TIMINGS_FILE.write_text(json.dumps(timings, indent=2), encoding='utf-8')
    except Exception:
        pass  # a missing estimate next time isn't worth failing a load over


def imputed_file_attributes(abc_cache):
    """(expected_size_bytes, local_path) for the imputed dataset, from the
    cache manifest (same private-API reach as is_data_file_cached, for the
    same reason), or (None, None) if it can't be determined."""
    try:
        attrs = abc_cache.cache._manifest.get_file_attributes(
            directory=IMPUTED_DATASET_DIRECTORY, file_name=IMPUTED_DATASET_FILE_NAME)
        return int(attrs.file_size), Path(attrs.local_path)
    except Exception:
        return None, None


def start_imputed_load_thread(imputed_state, abc_cache):
    """Starts loading the imputed gene dataset into imputed_state['adata'] on a thread of its own, and keeps the thread in
    imputed_state['load_thread'] (whoever wants the dataset meanwhile waits for it: see ensure_imputed_gene_dataset_loaded)."""
    def worker():
        try:
            imputed_state['adata'] = load_imputed_adata(abc_cache)
        except Exception as e:
            imputed_state['load_error'] = e
            # Cleared (not left pointing at this now-finished thread) so a
            # later attempt retries fresh instead of treating this as
            # permanent — same reasoning as ensure_spatial_load_started.
            imputed_state['load_thread'] = None

    thread = threading.Thread(target=worker, daemon=True)
    imputed_state['load_thread'] = thread
    thread.start()


def preload_imputed_gene_dataset(imputed_state, abc_cache):
    """Starts loading the imputed gene dataset in the background, so that it is ready by the time someone chooses Imputed Gene (the UMAP
    viewer calls this once its window is drawn). Only when the file is already in the local cache: this never starts a download of its own
    (a download is what ensure_imputed_gene_dataset_loaded asks first about). Returns True if a load was started."""
    if imputed_state['adata'] is not None or imputed_state['load_thread'] is not None:
        return False
    if is_data_file_cached(abc_cache, IMPUTED_DATASET_DIRECTORY, IMPUTED_DATASET_FILE_NAME) is not True:
        return False
    print("Loading the imputed gene dataset in the background...", flush=True)
    start_imputed_load_thread(imputed_state, abc_cache)
    return True


def ensure_imputed_gene_dataset_loaded(imputed_state, abc_cache):
    """Loads the imputed gene dataset into imputed_state['adata'] if it
    isn't already (or already in flight). Warns the user first that this
    can take a long time — but only when it actually might: skipped when
    is_data_file_cached() can confirm the file is already sitting in the
    local cache, so re-selecting 'Imputed Gene' (once it's been downloaded
    once — this session, or an earlier run) doesn't nag about a download
    that isn't going to happen. When the cached-or-not check itself can't
    be answered (is_data_file_cached() returns None — see its own
    docstring), this still warns, since an unnecessary prompt is a much
    smaller cost than an unannounced multi-minute download.

    `imputed_state` must be a dict with 'adata'/'load_thread'/'load_error'
    keys — shared by reference across every caller that wants to reuse the
    same load rather than each paying for their own (see
    prompt_section_selection_gui's own imputed_state, which callers should
    pass through/reuse rather than create a fresh one, same reasoning as
    ensure_spatial_load_started's spatial_cache).

    Returns True once imputed_state['adata'] is ready to use, False if the
    user declined the warning, cancelled, or the load failed — callers
    should leave whatever was selected before unchanged in that case."""
    if imputed_state['adata'] is not None:
        return True
    if imputed_state['load_thread'] is None:
        already_cached = is_data_file_cached(
            abc_cache, IMPUTED_DATASET_DIRECTORY, IMPUTED_DATASET_FILE_NAME,
        )
        if already_cached is not True:
            try:
                import centered_dialog
                confirmed = centered_dialog.ask_yes_no(
                    "Load imputed gene dataset?",
                    "The imputed gene dataset covers far more genes than the standard "
                    "500-gene MERFISH panel, but is much larger. The first time it's "
                    "loaded this session (and cached to disk, if not already), this can "
                    "take a long time. Continue?",
                )
            except Exception as e:
                print(f"Could not show confirmation dialog ({e}); not loading the imputed dataset.")
                return False
            if not confirmed:
                return False

        start_imputed_load_thread(imputed_state, abc_cache)

    thread = imputed_state['load_thread']
    # Hourglass over the app's other windows (e.g. the section window that
    # asked for the dataset) for the whole load. Opened before the dialog, so
    # the dialog keeps its normal arrow and Cancel still looks clickable.
    restore_cursors = set_wait_cursor_on_open_figures()
    proc_fig, cancel_flag, set_progress = show_processing_dialog(
        "Loading imputed gene dataset (first use this session; this can take a while)..."
    )
    try:
        # What can be measured, and what can't:
        #  - Download (first use only): the file grows on disk as it
        #    arrives (abc_atlas_access writes it in place), so bytes on disk
        #    over the manifest's size is a real fraction.
        #  - Everything after that (hash check of a fresh download, then
        #    opening the file, which loads all of obs): no hook exists. If an
        #    earlier open of an already-complete file was timed on this
        #    machine, show an estimate from that, parked below 100% and
        #    switching to the sweeping animation if it runs over; otherwise
        #    just the animation. The animation moves without claiming any
        #    amount of progress, which an empty bar did not.
        expected_bytes, local_path = imputed_file_attributes(abc_cache)
        loop_started = time.perf_counter()
        open_phase_started = None
        opened_from_complete_file = None  # was the file already whole when we started waiting?
        pulse = 0.0
        set_progress(0.0)
        while thread.is_alive() and not cancel_flag['cancelled']:
            elapsed = time.perf_counter() - loop_started
            try:
                size = local_path.stat().st_size if local_path is not None else None
            except OSError:
                size = 0
            if expected_bytes and size is not None and size < expected_bytes:
                if opened_from_complete_file is None:
                    opened_from_complete_file = False
                proc_fig._set_status(
                    f"Downloading: {size / 1e9:.1f} of {expected_bytes / 1e9:.1f} GB  ({elapsed:.0f} s)")
                set_progress(size / expected_bytes)
            else:
                if open_phase_started is None:
                    open_phase_started = time.perf_counter()
                    if opened_from_complete_file is None:
                        opened_from_complete_file = True
                phase_elapsed = time.perf_counter() - open_phase_started
                estimate = load_recorded_duration('imputed_open_s') if opened_from_complete_file else None
                what = "Opening dataset" if opened_from_complete_file else "Verifying download and opening"
                if estimate and phase_elapsed < estimate:
                    proc_fig._set_status(f"{what}…  ({phase_elapsed:.0f} s, about {estimate:.0f} s expected)")
                    set_progress(min(0.95, phase_elapsed / estimate))
                else:
                    proc_fig._set_status(f"{what}…  ({phase_elapsed:.0f} s)")
                    pulse += 0.04
                    proc_fig._set_pulse(pulse)
            proc_fig.canvas.flush_events()
            time.sleep(0.05)
        # Timed only when the file was already complete on arrival: after a
        # fresh download this phase also includes hashing the whole file,
        # which would make every later estimate wildly pessimistic.
        if (open_phase_started is not None and opened_from_complete_file
                and imputed_state['adata'] is not None and not cancel_flag['cancelled']):
            record_load_duration('imputed_open_s', time.perf_counter() - open_phase_started)
    finally:
        close_processing_dialog(proc_fig)
        restore_cursors()

    if cancel_flag['cancelled']:
        # As with spatial_cache's own cancel handling, the thread keeps
        # running in the background rather than being killed outright — a
        # later selection just finds imputed_state['adata'] already there
        # (or retries, if it had failed by then) instead of paying for a
        # second full load.
        print("Cancelled — the imputed dataset load continues in the background; try again shortly.")
        return False
    if imputed_state['adata'] is None:
        print(f"Could not load imputed gene dataset: {imputed_state['load_error']}")
        return False
    print(f"Loaded imputed gene dataset with {imputed_state['adata'].n_obs} cells x "
          f"{imputed_state['adata'].n_vars} genes.")
    return True





# 10 mutually distinguishable colors (Kelly/Boynton-style max-contrast set,
# avoiding low-saturation/near-white/near-black entries) for
# save_group_spatial_maps(): assigned by a category's *position* within
# its group of up to 10, so the same 10 colors — in the same order — are
# reused across every group's figure, unlike class_id_to_color()'s
# golden-angle hue stepping (great for spreading dozens of IDs apart on
# average, but with only 10 in play at once, some pairs land close enough
# in hue to look similar).
SPATIAL_MAP_COLORS = [
    '#e6194b',  # red
    '#3cb44b',  # green
    '#4363d8',  # blue
    '#f58231',  # orange
    '#911eb4',  # purple
    '#42d4f4',  # cyan
    '#f032e6',  # magenta
    '#9a6324',  # brown
    '#000075',  # navy
    '#808000',  # olive
]




GOLDEN_ANGLE = 0.6180339887498949  # shared with class_id_to_color's own local copy above

# Circle first, then nine more of matplotlib's filled markers — paired
# with CATEGORY_COLORS_PER_CYCLE (20) and the filled/open doubling in
# category_rank_shape below, 10 shapes x 20 colors x 2 fill styles = 400
# combinations before the whole cycle repeats.
CATEGORY_MARKER_SHAPES = ('o', 's', '^', '<', 'v', '>', 'D', 'P', 'X', '*')
CATEGORY_COLORS_PER_CYCLE = 20  # 10 hues x 2 saturation/value families — see category_rank_color


def category_rank_color(rank):
    """Deterministic color for a category's rank (0 = most common, by
    whatever ordering assigned that rank — see compute_ranked_category_
    colors/redraw_subclass, its own two callers) within a 20-color cycle:
    10 hues stepped by the golden angle (same mechanism as class_id_to_
    color), generated twice — once vivid, once muted.

    A single hue-only ring (class_id_to_color's own approach, used here
    until this was added) assigns well over 100 *numerically* distinct
    hues at fixed saturation/value, but human hue discrimination at fixed
    saturation/value tops out around a dozen-ish steps — past that,
    hues that are numerically well-separated still read as "the same
    color". Adding a second saturation/value family roughly doubles how
    many colors stay individually recognizable before the cycle repeats.
    Paired with category_rank_shape (below): once the 20-color cycle
    itself repeats, the shape changes too, so two categories that land on
    the same color again still read as different at a glance."""
    hue_index = rank % (CATEGORY_COLORS_PER_CYCLE // 2)
    family = (rank // (CATEGORY_COLORS_PER_CYCLE // 2)) % 2  # 0 = vivid, 1 = muted
    hue = (hue_index * GOLDEN_ANGLE) % 1.0
    saturation, value = (0.90, 0.95) if family == 0 else (0.45, 0.85)
    r, g, b = mcolors.hsv_to_rgb((hue, saturation, value))
    return mcolors.to_hex((r, g, b))


def category_rank_shape(rank):
    """(marker, is_open) for a category's rank — paired with category_
    rank_color (above). Every full color cycle (CATEGORY_COLORS_PER_CYCLE
    = 20 ranks) advances to the next marker in CATEGORY_MARKER_SHAPES (10
    of them); once every shape has been used with every color (20 * 10 =
    200 combinations), the same 200 repeat again as *open* (hollow)
    markers before the whole sequence truly repeats at rank 400
    (20 * 10 * 2) — see CATEGORY_MARKER_SHAPES's own comment. Derived from
    len(CATEGORY_MARKER_SHAPES) rather than a hardcoded count, so this
    stays correct if that tuple's length ever changes again."""
    combos_per_fill = CATEGORY_COLORS_PER_CYCLE * len(CATEGORY_MARKER_SHAPES)  # 200
    wrapped = rank % (combos_per_fill * 2)  # 400
    shape_index = (wrapped // CATEGORY_COLORS_PER_CYCLE) % len(CATEGORY_MARKER_SHAPES)
    is_open = wrapped >= combos_per_fill
    return CATEGORY_MARKER_SHAPES[shape_index], is_open






def raise_figure_window(fig):
    """Best-effort attempt to bring a figure's window to the front and give
    it input focus. Matters most for figures opened from inside another
    window's event handler (nested plt.show(), e.g. a confirmation dialog
    popped up mid double-click) — some backends don't automatically focus
    the new window, so clicks can silently land on the window underneath
    instead, making the new window's buttons look unresponsive."""
    try:
        window = fig.canvas.manager.window
    except AttributeError:
        return
    try:
        # Force the window to actually be created/mapped now, before the
        # focus-stealing calls below — on some backends the window object
        # exists but isn't visible yet until plt.show() runs its own show(),
        # which would be too late for lift()/raise_() to have any effect.
        fig.canvas.manager.show()
    except Exception:
        pass
    try:  # Tk (TkAgg)
        window.attributes('-topmost', 1)
        window.attributes('-topmost', 0)
        window.lift()
        window.focus_force()
    except Exception:
        pass
    try:  # Qt (QtAgg/Qt5Agg)
        window.raise_()
        window.activateWindow()
    except Exception:
        pass




# Rough allowance for OS window chrome (title bar + borders) that sits
# outside the matplotlib canvas but still counts toward the window's actual
# on-screen footprint — without this, a window whose *content* is exactly
# height_frac of the screen ends up slightly taller than that once its
# title bar is added, and can't fully fit; and centering math based only on
# the content size places it slightly off from the window's true center.
WINDOW_CHROME_HEIGHT_PX = 80
WINDOW_CHROME_WIDTH_PX = 20
# Bare desktop kept on every side of a newly-opened window (see
# fit_figure_window_to_work_area/center_figure_window) — enough that the
# window border itself stays grabbable for resizing rather than sitting
# flush against, or past, the screen edge.
WINDOW_EDGE_MARGIN_PX = 40

def compute_figsize_for_screen_height(aspect_wh, height_frac=0.8, default=(10, 7), width_frac=0.95):
    """Compute a (width, height) figsize in inches so the resulting
    *window* (content plus estimated chrome) is `height_frac` of the screen
    height tall at matplotlib's default figure DPI, keeping the given
    width/height aspect ratio. Sizing the figure *before* creating it this
    way (rather than resizing the window after plt.subplots()) avoids a
    real bug that approach had: resizing the Tk window post-creation
    doesn't reliably force matplotlib to relay out its axes at the new size
    before the window is shown, which clipped the right/bottom of the plot
    and pushed the button axes below the visible area. Falls back to
    `default` if screen dimensions can't be determined.

    The height-driven width is also capped at `width_frac` of the screen's
    own width (height reduced to match, preserving aspect_wh) — a wide
    aspect_wh (e.g. the UMAP window's 2.4:1) combined with a tall/narrow
    screen could otherwise compute a width wider than the screen itself,
    since height_frac alone says nothing about how wide the screen is."""
    try:
        import tkinter as tk
        root = tk.Tk()
        root.withdraw()
        screen_w_px, screen_h_px = root.winfo_screenwidth(), root.winfo_screenheight()
        root.destroy()
    except Exception:
        return default
    dpi = matplotlib.rcParams.get('figure.dpi', 100)
    content_h_px = max(screen_h_px * height_frac - WINDOW_CHROME_HEIGHT_PX, 1)
    height_in = content_h_px / dpi
    width_in = height_in * aspect_wh
    max_width_px = max(screen_w_px * width_frac - WINDOW_CHROME_WIDTH_PX, 1)
    max_width_in = max_width_px / dpi
    if width_in > max_width_in:
        width_in = max_width_in
        height_in = width_in / aspect_wh
    return width_in, height_in


def compute_ui_fontsize(button_height_frac=0.05):
    """Font size (points) shared by every button/text across the section
    grid picker, the single-section ROI picker, and the 'Processing...'
    dialog. Computed once, here, from a representative on-screen button
    height, rather than separately per window from that window's own
    figure size and button-height fraction — those two things used to
    differ enough between windows (the grid and single-section pickers use
    different button-height fractions of their own figure; the processing
    dialog uses a fixed, non-screen-relative figsize) to produce visibly
    different sizes, even though the two pickers' figures actually end up
    the same height in inches (compute_figsize_for_screen_height's height
    only depends on screen height and height_frac, not aspect ratio)."""
    ref_height_in = compute_figsize_for_screen_height(1.0, default=(10, 7))[1]
    button_height_pt = button_height_frac * ref_height_in * 72
    return max(8, min(16, button_height_pt * 0.45))


UI_BUTTON_FONTSIZE = compute_ui_fontsize()


def measured_text_width_pt(text, fontsize):
    """Width, in points, of `text` as matplotlib would actually render it at
    `fontsize` — measured via a throwaway off-screen Agg figure rather than
    estimated from character count, so it's exact regardless of font
    metrics. Used by compute_horizontal_radio_layout() to lay out hand-rolled
    radio widgets without needing to force a real draw of whatever
    (possibly large/expensive) figure they'll actually end up on —
    FigureCanvasAgg.get_renderer() lazily creates a renderer on its own,
    unlike the interactive TkAgg canvas's _get_renderer(), which returns None
    until that figure has actually been drawn once."""
    probe_fig = Figure()
    FigureCanvasAgg(probe_fig)
    probe_text = probe_fig.text(0, 0, text, fontsize=fontsize)
    bbox = probe_text.get_window_extent(renderer=probe_fig.canvas.get_renderer())
    return bbox.width * 72 / probe_fig.dpi


_TEXT_WIDTH_CACHE = {}




def compute_horizontal_radio_layout(options, fontsize, left_margin_pt=4.0, dot_label_gap_pt=4.0,
                                     option_gap_pt=14.0, right_margin_pt=4.0):
    """Lay out a hand-rolled horizontal radio widget (dot, label, dot, label,
    ... left to right) using each label's actually-measured rendered width,
    rather than fixed axes-fraction positions for the dots/labels (the
    previous approach in each of this file's 3 radio widgets — grid picker's
    Classes/Gene, single-section picker's Groups/Gene and Linear/Log).

    A fixed axes-fraction layout only ever fits the one screen size/font it
    happened to be tuned against: the dots are a fixed size in *points* (see
    each call site's own dot_size comment) and the labels' rendered width is
    also fixed in points (font metrics don't scale with window size), while
    a fixed *fraction* of the axes shrinks right along with the axes on a
    smaller window or narrower figure — so the same fractional gap that
    looked fine on one screen can end up smaller than the dot's own radius,
    or smaller than a label's own rendered width, on another. That's what
    let a label overlap its own dot, or (on a narrow enough window) overlap
    the *next* dot: the second dot's fixed 0.62-fraction position had no
    relationship to how wide the first label actually rendered.

    Returns (dot_x_pt, label_x_pt, total_width_pt): the first two are tuples
    of point offsets from the axes' left edge — not yet axes-fraction, since
    that depends on how wide the caller ends up making the axes (it should
    be at least `total_width_pt`, converted via the figure's width in
    points, to guarantee no overlap; the caller may make it wider still, in
    which case these positions are simply left-packed with extra room to
    their right)."""
    dot_radius_pt = fontsize / 2
    dot_x_pt = []
    label_x_pt = []
    cursor_pt = left_margin_pt
    for i, opt in enumerate(options):
        dot_center_pt = cursor_pt + dot_radius_pt
        dot_x_pt.append(dot_center_pt)
        label_start_pt = dot_center_pt + dot_radius_pt + dot_label_gap_pt
        label_x_pt.append(label_start_pt)
        cursor_pt = label_start_pt + measured_text_width_pt(opt, fontsize)
        if i < len(options) - 1:
            cursor_pt += option_gap_pt
    total_width_pt = cursor_pt + right_margin_pt
    return tuple(dot_x_pt), tuple(label_x_pt), total_width_pt






def normalize_tk_scaling(fig):
    """Force Tk's own internal pixel-per-point scaling to match matplotlib's
    dpi assumption (dpi pixels/inch = dpi/72 pixels/point). Without this,
    even with the process declared DPI-aware (see _make_windows_dpi_aware),
    Tk independently applies its own DPI-derived scaling (observed ~2x on a
    150%-scaled 4K display) when it actually realizes the canvas widget's
    size — and TkAgg then auto-adjusts the figure to match whatever size Tk
    actually rendered, silently inflating the window well past the size
    compute_figsize_for_screen_height intended. That's what was pushing the
    bottom button row off-screen: not a missing height allowance, but the
    whole window ending up bigger than requested. Must run right after the
    figure/window is created, before anything else is laid out."""
    try:
        window = fig.canvas.manager.window
        dpi = fig.dpi
        window.tk.call('tk', 'scaling', dpi / 72.0)
    except Exception:
        pass


def _windows_monitor_work_area(window):
    """The work area (screen bounds minus the taskbar, whichever edge it's
    docked to) of whichever monitor `window` is actually on, as (left, top,
    width, height) in physical pixels — or None on failure/non-Windows.

    Windows' own 'zoomed' window state (ShowWindow(SW_MAXIMIZE) — what
    window.state('zoomed') triggers) is *supposed* to already respect this
    automatically, but doing it explicitly avoids relying on that: a window
    that's SYSTEM (not per-monitor) DPI-aware, as this process declares
    itself (see _make_windows_dpi_aware), can end up with 'zoomed' computed
    against stale or wrong-monitor metrics, letting the bottom of the window
    land under the taskbar. Querying the monitor's rcWork rectangle directly
    via MonitorFromWindow/GetMonitorInfoW and setting geometry to exactly
    that is unambiguous regardless of DPI-awareness quirks, and — using
    MonitorFromWindow rather than the primary screen's metrics — correct on
    whichever monitor the window actually ends up on, not just the primary
    one."""
    if sys.platform != 'win32':
        return None
    try:
        import ctypes
        from ctypes import wintypes

        class MONITORINFO(ctypes.Structure):
            _fields_ = [
                ('cbSize', wintypes.DWORD),
                ('rcMonitor', wintypes.RECT),
                ('rcWork', wintypes.RECT),
                ('dwFlags', wintypes.DWORD),
            ]

        user32 = ctypes.windll.user32
        # HMONITOR is a pointer-sized handle — on 64-bit Windows that's 8
        # bytes. ctypes defaults an undeclared restype to c_int (4 bytes,
        # signed), which silently truncates/corrupts a real 64-bit handle
        # value rather than raising — GetMonitorInfoW then either fails
        # outright on the corrupted handle, or (worse) succeeds against
        # whatever monitor that garbage handle happens to coincidentally
        # resolve to, producing a plausible-looking but wrong rectangle.
        # Declaring these explicitly is what makes the handle round-trip
        # correctly.
        user32.MonitorFromWindow.restype = ctypes.c_void_p
        user32.MonitorFromWindow.argtypes = [wintypes.HWND, wintypes.DWORD]
        user32.GetMonitorInfoW.restype = wintypes.BOOL
        user32.GetMonitorInfoW.argtypes = [ctypes.c_void_p, ctypes.POINTER(MONITORINFO)]

        hwnd = window.winfo_id()
        MONITOR_DEFAULTTONEAREST = 2
        monitor = user32.MonitorFromWindow(hwnd, MONITOR_DEFAULTTONEAREST)
        if not monitor:
            return None

        info = MONITORINFO()
        info.cbSize = ctypes.sizeof(MONITORINFO)
        if not user32.GetMonitorInfoW(monitor, ctypes.byref(info)):
            return None
        r = info.rcWork
        left, top, right, bottom = r.left, r.top, r.right, r.bottom
        width, height = right - left, bottom - top
        if width <= 0 or height <= 0:
            return None
        return left, top, width, height
    except Exception:
        return None


def maximize_figure_window(fig):
    """Best-effort: maximize a figure's window to fill the current monitor's
    work area (screen bounds minus the taskbar) — explicitly, via
    _windows_monitor_work_area, rather than Windows' own 'zoomed' window
    state, which is supposed to already exclude the taskbar but doesn't
    reliably here (see that function's docstring). Falls back to 'zoomed'
    (or Tk's '-zoomed' attribute, for X11/some window managers) if the work
    area can't be determined — e.g. non-Windows, or the ctypes calls fail.

    Maximizing programmatically like this doesn't reliably deliver a
    correctly-sized resize_event to matplotlib the way an interactive drag-
    resize does — see on_configure in prompt_subregion_selection, which
    this relies on the caller having already wired up (via the toplevel's
    own <Configure> binding, not just matplotlib's resize_event) to actually
    pick up the resulting size once Tk finishes applying it."""
    try:
        window = fig.canvas.manager.window
    except AttributeError:
        return
    work_area = _windows_monitor_work_area(window)
    if work_area is not None:
        left, top, width, height = work_area
        try:
            # Pass 1: content-sized to the *entire* work area — this alone
            # overshoots the window's true *outer* (decorated) bounds past
            # the work area's bottom edge, under the taskbar, by however
            # tall the title bar is: Tk's own docs are explicit that wm
            # geometry's WxH excludes window decorations (border/title bar)
            # entirely, unlike most other toolkits' "maximize", which sizes
            # the outer window to fit the work area. This is exactly what
            # let the button row end up under the taskbar even though the
            # work area itself (confirmed via _windows_monitor_work_area)
            # was already correct — real OS maximize doesn't have this
            # problem since it operates on the outer window directly.
            window.geometry(f"{width}x{height}+{left}+{top}")
            # Needs the window actually realized before its decoration size
            # can be measured (pass 2, below) — deiconify+update_idletasks
            # mirrors center_figure_window's own two-pass approach, just
            # done eagerly here instead of deferred via after_idle, since
            # this function is expected to complete synchronously.
            window.deiconify()
            window.update_idletasks()
            # Pass 2: same trick as center_figure_window's own two-pass
            # centering — winfo_rootx()/winfo_rooty() give the *content*
            # area's screen position, winfo_x()/winfo_y() give the *outer*
            # (decorated) window's, so the gap between them is the real
            # title bar/border chrome. Shrinking the requested content size
            # by that amount is what keeps the window's true outer bounds
            # within the work area instead of extending past it.
            chrome_w = max(0, window.winfo_rootx() - window.winfo_x()) * 2
            chrome_h = max(0, window.winfo_rooty() - window.winfo_y())
            if chrome_w > 0 or chrome_h > 0:
                adj_width = max(1, width - chrome_w)
                adj_height = max(1, height - chrome_h)
                window.geometry(f"{adj_width}x{adj_height}+{left}+{top}")
            return
        except Exception:
            pass
    try:
        window.state('zoomed')
    except Exception:
        try:
            window.attributes('-zoomed', True)
        except Exception:
            pass


def assert_figure_content_size(fig, fig_w_in, fig_h_in):
    """Best-effort: explicitly (re)assert the *canvas widget's* own pixel
    size, computed from (fig_w_in, fig_h_in) at the figure's own dpi — same
    reasoning as maximize_figure_window (whose own target is instead the
    monitor's work area, not a specific figsize): normalize_tk_scaling's fix
    for Tk's own DPI-derived auto-scaling only reliably takes effect once
    the window is actually realized, which otherwise happens lazily
    whenever show_figure_blocking() finally maps it — deiconifying/
    realizing it eagerly here instead (same early timing normalize_tk_
    scaling itself already requires) and reasserting the intended size right
    after catches and corrects any residual mismatch, instead of leaving the
    window to visibly snap to a different (usually smaller) size right
    after it first appears, with no later maximize_figure_window call
    around to paper over it the way the other pickers' own windows get.

    Deliberately resizes the *canvas* (via its own Tk widget's width/height
    options — the same thing FigureCanvasTkAgg sets once at construction)
    rather than calling `window.geometry("WxH")` directly on the toplevel:
    an explicit `wm geometry WxH` call turns off Tk's automatic geometry
    propagation for that window, so if the canvas subsequently resized
    itself again on its own (exactly the behavior this is working around),
    the *window* would no longer follow it — leaving the window pinned at
    this call's size while the canvas shrinks inside it, i.e. a visible gap
    below/right of the actual content. Resizing the canvas widget directly
    leaves propagation on, so the (packed, fill+expand) toplevel keeps
    tracking the canvas's true size the same way it already does for every
    later resize, including whatever internal one this is meant to catch."""
    try:
        canvas_widget = fig.canvas.get_tk_widget()
        dpi = fig.dpi
        w = max(1, int(round(fig_w_in * dpi)))
        h = max(1, int(round(fig_h_in * dpi)))
        canvas_widget.configure(width=w, height=h)
        window = fig.canvas.manager.window
        window.deiconify()
        window.update_idletasks()
    except Exception:
        pass


def _centered_position(area_left, area_top, area_w, area_h, win_w, win_h,
                        margin_px=WINDOW_EDGE_MARGIN_PX):
    """Top-left corner that centers a win_w x win_h window in the given work
    area, clamped so it never extends past either edge.

    The clamp is the point: a plain centered `area_left + (area_w - win_w)
    // 2` goes *negative* for a window wider than the work area, and the
    previous max(0, ...) guard only stopped it from running off the left —
    the overflow simply moved to the right edge instead, which is what put
    the right side of the UMAP window off-screen. Here an oversized window
    is pinned to the margin on the top/left, so at least its own controls
    and the resize border on that side stay reachable."""
    x = area_left + (area_w - win_w) // 2
    y = area_top + (area_h - win_h) // 2
    x = max(area_left + margin_px, min(x, area_left + area_w - win_w - margin_px))
    y = max(area_top + margin_px, min(y, area_top + area_h - win_h - margin_px))
    return int(x), int(y)


def fit_figure_window_to_work_area(fig, margin_px=WINDOW_EDGE_MARGIN_PX):
    """Shrink `fig` if its window wouldn't fit on the monitor it actually
    opened on, leaving at least `margin_px` of bare desktop on every side so
    the window border stays grabbable for resizing.

    compute_figsize_for_screen_height (which picked this figure's size
    before it existed) can only consult a throwaway Tk root's
    winfo_screenwidth()/height() — the *primary* monitor's full bounds,
    taskbar included. That's the wrong number in two common cases: the
    window opens on a different, smaller monitor, or the taskbar/DPI
    scaling eats enough that the nominally-fitting size doesn't. Either way
    the result was a window whose right edge ran off the screen. Once the
    window really exists, MonitorFromWindow gives the true work area of the
    monitor it's on, so the size can be corrected before anything is laid
    out against it.

    Aspect ratio is deliberately *not* preserved: each dimension is clamped
    independently, since the goal is to fit the available desktop, and this
    window's own layout is fully responsive to whatever size it ends up
    with (see its on_figure_resize/apply_region_layout)."""
    try:
        window = fig.canvas.manager.window
        work_area = _windows_monitor_work_area(window)
        if work_area is not None:
            _area_left, _area_top, area_w, area_h = work_area
        else:
            area_w, area_h = window.winfo_screenwidth(), window.winfo_screenheight()
        dpi = fig.dpi
        # Budget is for the *outer* (decorated) window, so the chrome has to
        # come out before comparing against the figure's own content size.
        max_content_w = area_w - 2 * margin_px - WINDOW_CHROME_WIDTH_PX
        max_content_h = area_h - 2 * margin_px - WINDOW_CHROME_HEIGHT_PX
        cur_w_in, cur_h_in = fig.get_size_inches()
        new_w_in = min(cur_w_in, max_content_w / dpi)
        new_h_in = min(cur_h_in, max_content_h / dpi)
        if new_w_in <= 0 or new_h_in <= 0:
            return  # implausible work-area reading — leave the figure alone
        if new_w_in < cur_w_in or new_h_in < cur_h_in:
            fig.set_size_inches(new_w_in, new_h_in)
            assert_figure_content_size(fig, new_w_in, new_h_in)
    except Exception:
        pass


def center_figure_window(fig):
    """Best-effort: center a figure's window on screen (position only,
    doesn't touch its size). Positions in two passes.

    Pass 1 (immediate): uses the figure's own size in inches (known
    precisely — we just computed it via compute_figsize_for_screen_height)
    converted to pixels via its DPI, plus the flat WINDOW_CHROME_* guess for
    title bar/border overhead, rather than querying the Tk window's
    winfo_width()/winfo_height(): those can still reflect a stale,
    not-yet-realized size right after creation, which previously centered
    large windows as if they were tiny — anchoring near the screen's center
    point and then letting the window's real, much larger size overflow
    past it toward the bottom right. This pass just gets the window roughly
    in place immediately, instead of at Tk's default spawn point, while it
    is being realized.

    Pass 2 (deferred, via after_idle): WINDOW_CHROME_* is a single flat
    estimate and doesn't scale with OS DPI scaling, so on a scaled display
    the real chrome can be taller than assumed, leaving the window visibly
    off-center (low). Once Tk has actually mapped the window, its true
    decorated size can be measured instead of guessed: winfo_rootx()/
    winfo_rooty() give the screen position of the window's *content* area,
    while winfo_x()/winfo_y() give the screen position of the *outer*
    window (as tracked by Tk) — the gap between them is the real title
    bar/border chrome, a standard Tkinter trick for recovering WM
    decoration size. This has to be deferred (scheduled instead of read
    immediately) because right after creation those values are still stale
    placeholders, same reason pass 1 doesn't use them directly — after_idle
    guarantees Tk has processed the map/configure events first. Both
    show_figure_blocking()'s and show_processing_dialog()'s event pumping
    (flush_events()/update()) give this callback a chance to run.

    Both passes center within the monitor's *work area* (via
    _windows_monitor_work_area — same query maximize_figure_window uses),
    not its full screen bounds: centering in the full screen height put the
    window's midpoint below the work area's actual visual center by half
    the taskbar's height, since the taskbar eats into the bottom (or
    whichever edge it's docked to) without the centering math accounting
    for it. Falls back to winfo_screenwidth()/height() (full screen, taskbar
    included) if the work-area query fails — e.g. non-Windows."""
    try:
        window = fig.canvas.manager.window
        fig_w_in, fig_h_in = fig.get_size_inches()
        dpi = fig.dpi
        w = int(fig_w_in * dpi) + WINDOW_CHROME_WIDTH_PX
        h = int(fig_h_in * dpi) + WINDOW_CHROME_HEIGHT_PX
        work_area = _windows_monitor_work_area(window)
        if work_area is not None:
            area_left, area_top, screen_w, screen_h = work_area
        else:
            area_left, area_top = 0, 0
            screen_w, screen_h = window.winfo_screenwidth(), window.winfo_screenheight()
        x, y = _centered_position(area_left, area_top, screen_w, screen_h, w, h)
        window.geometry(f"+{x}+{y}")

        def _recenter_with_real_geometry():
            try:
                window.update_idletasks()
                content_w, content_h = window.winfo_width(), window.winfo_height()
                if content_w <= 1 or content_h <= 1:
                    return  # window still not realized; keep the pass-1 estimate
                chrome_w = max(0, window.winfo_rootx() - window.winfo_x()) * 2
                chrome_h = max(0, window.winfo_rooty() - window.winfo_y())
                real_w = content_w + chrome_w
                real_h = content_h + chrome_h
                # Shrink-to-fit, from *measured* geometry rather than the
                # predicted size fit_figure_window_to_work_area works from.
                # That prediction multiplies figsize by fig.dpi and adds the
                # flat WINDOW_CHROME_* estimate — both of which can be wrong
                # (Tk's own DPI-derived scaling, real WM decoration size), and
                # when they are, a window that "fits" on paper still opens
                # with its right edge past the screen. Here the true decorated
                # size is already known, so any overflow is exact and can be
                # taken straight back off the canvas.
                over_w = max(0, real_w - (screen_w - 2 * WINDOW_EDGE_MARGIN_PX))
                over_h = max(0, real_h - (screen_h - 2 * WINDOW_EDGE_MARGIN_PX))
                if over_w > 0 or over_h > 0:
                    new_content_w = max(1, content_w - over_w)
                    new_content_h = max(1, content_h - over_h)
                    # Canvas widget (not window.geometry) for the same reason
                    # assert_figure_content_size documents: an explicit
                    # `wm geometry WxH` turns off Tk's geometry propagation,
                    # leaving the window pinned while the canvas resizes
                    # inside it. set_size_inches keeps matplotlib's own idea
                    # of the size in step, so any resize_event handler
                    # relayouts against the corrected size.
                    fig.set_size_inches(new_content_w / fig.dpi, new_content_h / fig.dpi)
                    fig.canvas.get_tk_widget().configure(width=new_content_w, height=new_content_h)
                    window.update_idletasks()
                    real_w = new_content_w + chrome_w
                    real_h = new_content_h + chrome_h
                x2, y2 = _centered_position(area_left, area_top, screen_w, screen_h, real_w, real_h)
                window.geometry(f"+{x2}+{y2}")
            except Exception:
                pass

        window.after_idle(_recenter_with_real_geometry)
    except Exception:
        pass


def show_figure_blocking(fig):
    """Show `fig` and block until *that specific figure* is closed — not
    plt.show()'s usual behavior of blocking until every open figure is
    closed. This matters whenever a figure is opened while another one
    (e.g. the section grid picker) is still open: with plain plt.show(),
    closing the nested window wouldn't return control, since the call is
    still waiting on the outer figure too — the nested window's buttons
    would appear to do nothing (they do work; the call just never returns)."""
    raise_figure_window(fig)
    fig.canvas.manager.show()
    fig.canvas.draw_idle()
    reraise_timer_id = None
    tk_widget = None
    try:
        # A window manager's own re-stacking in response to another window
        # closing right as this one opens (e.g. a "Processing..." dialog
        # being dismissed just before this call) can lag behind that close
        # by a moment — landing *after* the raise_figure_window() call
        # above and silently undoing it, which read as this window
        # intermittently dropping behind another one shortly after
        # appearing. Raising again a beat later catches that.
        tk_widget = fig.canvas.get_tk_widget()
        reraise_timer_id = tk_widget.after(200, lambda: raise_figure_window(fig))
    except Exception:
        pass
    while plt.fignum_exists(fig.number):
        fig.canvas.flush_events()
        time.sleep(0.02)
    if reraise_timer_id is not None:
        # Without this, closing the window within that 200ms window left
        # this timer to fire after Tk had already destroyed it — not
        # something raise_figure_window's own try/except can catch, since
        # the failure happens in Tcl's own "after" dispatch, before any of
        # our Python code (including that except) gets a chance to run:
        # 'invalid command name "...<lambda>" ("after" script)'.
        try:
            tk_widget.after_cancel(reraise_timer_id)
        except Exception:
            pass


def show_processing_dialog(message):
    """Non-blocking 'Processing...' popup (a Tk window, like the rest of the app) with a progress bar and Cancel button, for covering genuinely slow
    work (not just asking permission up front, which just makes the user wait through two delays instead of one). Does NOT start its own blocking
    event loop — the caller must periodically call win.canvas.flush_events() (e.g. while polling a worker thread) to keep the window responsive
    and let a Cancel click actually register, and is responsible for calling close_processing_dialog(win) when the work finishes or is cancelled.
    (The `canvas.flush_events` name is from when this was a matplotlib window: the callers still use it.)

    Also attaches win._set_pulse(position) (a "still working" animation for unmeasurable work) and win._set_status(text) (a line under the bar) for
    callers that want them.

    Returns (win, cancel_flag, set_progress): cancel_flag['cancelled'] becomes True once Cancel is clicked (or the window is closed);
    set_progress(fraction) fills the bar to `fraction` (0.0-1.0, clamped); set_progress(None) resets the bar to empty, for work whose progress
    isn't (yet) known."""
    import tkinter as tk
    from tkinter import ttk
    own_root = None
    parent = getattr(tk, '_default_root', None)
    try:
        if parent is not None:
            parent.winfo_exists()
    except tk.TclError:
        parent = None
    if parent is None:
        own_root = parent = tk.Tk()
        own_root.withdraw()
    top = tk.Toplevel(parent)
    top.title("Processing")
    top.attributes('-topmost', True)
    top.resizable(False, False)
    cancel_flag = {'cancelled': False}

    def on_cancel():
        cancel_flag['cancelled'] = True

    top.protocol("WM_DELETE_WINDOW", on_cancel)
    frame = ttk.Frame(top, padding=(20, 16, 20, 14))
    frame.pack()
    ttk.Label(frame, text=message, justify='center', wraplength=460).pack(pady=(0, 12))
    bar = ttk.Progressbar(frame, mode='determinate', length=440, maximum=1000)
    bar.pack()
    status = tk.StringVar(top, "")
    ttk.Label(frame, textvariable=status, foreground='gray35').pack(pady=(6, 10))
    ttk.Button(frame, text="Cancel", command=on_cancel).pack()
    top.bind('<Escape>', lambda e: on_cancel())

    animating = {'on': False}

    def set_progress(fraction):
        if animating['on']:
            bar.stop()
            bar.configure(mode='determinate')
            animating['on'] = False
        bar['value'] = 0 if fraction is None else 1000 * max(0.0, min(1.0, fraction))

    def set_pulse(position):
        """Indeterminate "still working" animation, for work whose progress can't be measured: a block sweeping along the bar. It moves by itself while the
        caller keeps pumping events (the `position` is no longer needed)."""
        if not animating['on']:
            bar.configure(mode='indeterminate')
            bar.start(15)
            animating['on'] = True

    def set_status(text):
        status.set(text)

    class Window:
        pass

    win = Window()
    win._set_pulse = set_pulse
    win._set_status = set_status
    win._top, win._own_root = top, own_root

    class Canvas:
        @staticmethod
        def flush_events():
            try:
                top.update()
            except tk.TclError:
                pass                                             # (closed while polling)

    win.canvas = Canvas()
    top.update_idletasks()
    w, h = top.winfo_reqwidth(), top.winfo_reqheight()
    top.geometry(f"+{max(0, (top.winfo_screenwidth() - w) // 2)}+{max(0, (top.winfo_screenheight() - h) // 3)}")
    top.lift()
    top.update()
    return win, cancel_flag, set_progress


def close_processing_dialog(win):
    """Close the window made by show_processing_dialog (and the hidden root it had to make, if there was no window to hang it on)."""
    try:
        win._top.destroy()
    except Exception:
        pass
    if getattr(win, '_own_root', None) is not None:
        try:
            win._own_root.destroy()
        except Exception:
            pass


def run_with_loading_window(message, work, show_after_s=0.4):
    """Run `work()` on a background thread and return its result (re-raising
    anything it raised), showing a small "loading" window if it takes longer
    than `show_after_s`. Exists because a few steps between a click and the
    next window appearing are long blocking reads (opening the 7+ GB raw h5ad
    loads all ~4.3 M rows of obs even in backed mode; cold off disk that is
    many seconds), during which the previous window is already gone and
    nothing else is on screen — so the app looks hung when it isn't.

    The window is only created after `show_after_s`, so a fast (already
    cached) step doesn't flash a window open and shut. With 0 it is made
    and painted before `work` starts, so a worker that holds the GIL for
    seconds (big pandas lookups) cannot keep it from appearing. It shows `message`, an
    indeterminate bar, and elapsed seconds, repainted by this thread's own
    polling loop (Tk stays on the main thread; only `work` runs elsewhere —
    work must not touch any GUI object). No Cancel button: a blocking file
    read can't be interrupted, and a button that did nothing would be worse
    than none. If Tk isn't available it just waits, printing the message."""
    outcome = {}

    def runner():
        try:
            outcome['value'] = work()
        except BaseException as e:  # re-raised on the calling thread below
            outcome['error'] = e

    print(message)
    thread = threading.Thread(target=runner, daemon=True)
    started = time.perf_counter()
    ui = None
    if show_after_s <= 0:
        thread_started = False
    else:
        thread.start()
        thread_started = True

    while not thread_started or thread.is_alive():
        elapsed = time.perf_counter() - started
        if ui is None and elapsed >= show_after_s:
            try:
                import tkinter as tk
                from tkinter import ttk
                root = tk.Tk()
                root.title("Loading")
                root.attributes('-topmost', True)
                root.protocol("WM_DELETE_WINDOW", lambda: None)  # nothing to cancel; ignore close
                tk.Label(root, text=message, justify='left', padx=18, pady=(14), wraplength=460).pack()
                bar = ttk.Progressbar(root, mode='indeterminate', length=420)
                bar.pack(padx=18)
                bar.start(12)
                elapsed_label = tk.Label(root, text="", fg='gray35', pady=10)
                elapsed_label.pack()
                root.update_idletasks()
                w, h = root.winfo_reqwidth(), root.winfo_reqheight()
                root.geometry(f"+{(root.winfo_screenwidth() - w) // 2}+{(root.winfo_screenheight() - h) // 3}")
                ui = (root, elapsed_label, bar)
            except Exception:
                ui = False  # no usable Tk; just wait
        if ui:
            try:
                ui[1].config(text=f"{int(elapsed)} s")
                ui[0].update()
            except Exception:
                ui = False
        if not thread_started:
            thread.start()                       # (after the first paint of the window)
            thread_started = True
        time.sleep(0.03)
    thread.join()
    if ui:
        try:
            ui[2].stop()                       # (its pending step would otherwise run on the destroyed window and print an error)
            # Hidden, and the screen brought up to date, before it is destroyed: on macOS a window destroyed with no event handled after it
            # stays on the screen, and the step that follows (minutes of PCA and UMAP) handles none, so "Waiting for the atlas" stayed up.
            ui[0].withdraw()
            ui[0].update()
            ui[0].destroy()
        except Exception:
            pass
    if 'error' in outcome:
        raise outcome['error']
    return outcome.get('value')


class StepWindow:
    """One "Computing" window for the whole of a new run's computation (Steps 2-8: annotations, filtering, preprocessing, PCA, UMAP,
    clustering, plots, saving), showing the step under way and its time. The long calls run on a background thread through run(), this
    thread repainting the window meanwhile (unrepainted, it would freeze, "Not Responding" on Windows); log_status shows each step's
    message in it too. Opened by the first run(), closed by close() (before the viewer opens, or when the launch screen comes back).
    Not "always on top": the questions a run can ask on the way (how many cells to use, a failed step) must not be hidden behind it."""

    def __init__(self):
        self.ui = None
        self.started = None

    def _ensure(self):
        if self.ui is not None:
            return self.ui
        try:
            import tkinter as tk
            from tkinter import ttk
            root = tk.Tk()
            root.title("Computing")
            root.protocol("WM_DELETE_WINDOW", lambda: None)            # (nothing to cancel)
            tk.Label(root, text="Computing the new run. This can take a while; the console shows the details.",
                     justify='left', padx=18, pady=(14), wraplength=460).pack(anchor='w')
            step = tk.Label(root, text="", justify='left', padx=18, wraplength=460, font=('TkDefaultFont', 10, 'bold'))
            step.pack(anchor='w')
            bar = ttk.Progressbar(root, mode='indeterminate', length=420)
            bar.pack(padx=18, pady=(8, 0))
            bar.start(12)
            elapsed = tk.Label(root, text="", fg='gray35', pady=10)
            elapsed.pack()
            root.update_idletasks()
            w, h = root.winfo_reqwidth(), root.winfo_reqheight()
            root.geometry(f"+{(root.winfo_screenwidth() - w) // 2}+{(root.winfo_screenheight() - h) // 3}")
            root.lift()
            self.ui = (root, step, bar, elapsed)
        except Exception:
            self.ui = False                                                 # (no usable Tk: the console alone)
        return self.ui

    def set(self, message):
        """Show `message` as the step under way (and repaint). Only from the main thread: Tk is not to be touched from another."""
        if threading.current_thread() is not threading.main_thread() or not self.ui:
            return
        try:
            self.ui[1].config(text=message)
            self.started = time.perf_counter()
            self.ui[3].config(text="")
            self.ui[0].update()
        except Exception:
            self.ui = False

    def run(self, message, work):
        """`work()` on a background thread, the window showing `message` and the time it has taken; returns its result (re-raising what it
        raised). `work` must not touch any GUI object."""
        ui = self._ensure()
        if ui:
            self.set(message)
        outcome = {}

        def runner():
            try:
                outcome['value'] = work()
            except BaseException as e:
                outcome['error'] = e
        # While the work runs on its thread: matplotlib draws with the Agg backend (TkAgg made a window for each figure, on the wrong thread: the program
        # died), and garbage collection waits (a collection run on this thread that frees a Tk object ends the program the same way).
        previous_backend = matplotlib.get_backend()
        gc_was_on = gc.isenabled()
        try:
            plt.switch_backend('agg')
        except Exception:
            previous_backend = None
        # The leftover Tk objects of the windows closed so far (the section grid, the options window) are freed here, on the Tk thread: pynndescent
        # (the neighbor graph of Step 6) calls gc.collect() itself, which would otherwise free them on the worker thread, and Tcl ends the program
        # ("Tcl_AsyncDelete: async handler deleted by the wrong thread"). Garbage collection being off does not stop that call.
        gc.collect()
        gc.disable()
        try:
            thread = threading.Thread(target=runner, daemon=True)
            thread.start()
            while thread.is_alive():
                if self.ui:
                    try:
                        self.ui[3].config(text=f"{int(time.perf_counter() - self.started)} s")
                        self.ui[0].update()
                    except Exception:
                        self.ui = False
                time.sleep(0.03)
            thread.join()
        finally:
            if gc_was_on:
                gc.enable()
            if previous_backend and previous_backend.lower() != 'agg':
                try:
                    plt.switch_backend(previous_backend)
                except Exception:
                    pass
        if 'error' in outcome:
            raise outcome['error']
        return outcome.get('value')

    def close(self):
        if self.ui:
            try:
                self.ui[2].stop()
                self.ui[0].withdraw()
                self.ui[0].update()                    # (gone from the screen now, also on macOS: see run_with_loading_window)
                self.ui[0].destroy()
            except Exception:
                pass
        self.ui = None


compute_window = StepWindow()












class UserCancelledSelection(Exception):
    """Raised when the user explicitly cancels a top-level GUI picker
    (closes the window, clicks Cancel/Exit, declines a whole-brain
    confirmation) — distinct from the window failing to open/run at all
    (a real error), so the main script's session loop can return to the
    startup panel for another session instead of falling back to a
    console prompt, which is what any *other* exception from these same
    call sites means."""




def format_scalebar_length(length_um):
    """'500 \u00b5m' below 1 mm, otherwise 'N mm' (e.g. '2 mm', '1.5 mm')."""
    if length_um < 1000:
        return f'{length_um:g} \u00b5m'
    return f'{length_um / 1000:g} mm'


def nice_scalebar_length_um(view_width_mm, target_fraction=SCALEBAR_TARGET_FRACTION,
                             nice_lengths_um=SCALEBAR_NICE_LENGTHS_UM):
    """The largest length (in micrometers, from `nice_lengths_um`) that's no
    more than `target_fraction` of a view this wide (in millimeters — the
    ABC atlas's own native data units; see HOVER_CELL_RADIUS_DATA_UNITS's
    comment elsewhere in this file). Used for a scale bar that reads a
    clean, round physical length at any zoom level: it shrinks as the view
    narrows and grows as it widens, always snapping to one of these fixed
    lengths, so it can never overflow the view or shrink away to nothing."""
    target_um = view_width_mm * 1000 * target_fraction
    candidates = [um for um in nice_lengths_um if um <= target_um]
    return candidates[-1] if candidates else nice_lengths_um[0]


def compute_scalebar_geometry(ax, margin_fraction=SCALEBAR_MARGIN_FRACTION,
                               text_gap_fraction=SCALEBAR_TEXT_GAP_FRACTION):
    """Where a scale bar belongs in `ax`'s lower-left corner right now: a
    dict with the bar's two endpoints, the label's position, and its text —
    all in `ax`'s own data coordinates, so the result is a real vector
    object under savefig too, not a fixed-pixel overlay drawn after the
    fact. Positioned via transAxes -> transData (screen-relative margins,
    not a data-space offset) rather than assuming which data direction is
    'right'/'up' on screen, so this works whether or not the axes' y-axis
    is inverted — e.g. every section panel's tissue-orientation flip (see
    invert_yaxis() elsewhere in this file). Returns None if the axes has no
    usable width yet (e.g. mid-setup, before its limits are set)."""
    x0, x1 = ax.get_xlim()
    view_width_mm = abs(x1 - x0)
    if not np.isfinite(view_width_mm) or view_width_mm <= 0:
        return None
    length_um = nice_scalebar_length_um(view_width_mm)
    length_mm = length_um / 1000
    inv = ax.transData.inverted()
    corner_x, corner_y = inv.transform(ax.transAxes.transform((margin_fraction, margin_fraction)))
    # A second point just to the right of the corner in axes-fraction terms
    # (screen-right, always) — comparing its data-x against corner_x gives
    # the sign that means 'rightward' in *this* axes' own data coordinates,
    # without assuming x is never inverted.
    probe_x, _probe_y = inv.transform(ax.transAxes.transform((margin_fraction + 0.01, margin_fraction)))
    x_sign = 1 if probe_x >= corner_x else -1
    end_x = corner_x + x_sign * length_mm
    text_y = inv.transform(ax.transAxes.transform((margin_fraction, margin_fraction + text_gap_fraction)))[1]
    return {
        'bar_x': (corner_x, end_x), 'bar_y': (corner_y, corner_y),
        'text_xy': (corner_x + x_sign * length_mm / 2, text_y),
        'label': format_scalebar_length(length_um),
    }


def apply_scalebar_geometry(line, text, geometry):
    """Push `geometry` (see compute_scalebar_geometry; None hides the bar)
    onto an existing scale bar's Line2D/Text — the update half of
    build_section_scalebar, for a bar that has to keep tracking a live,
    still-changing view rather than being drawn once for a static export."""
    if geometry is None:
        line.set_visible(False)
        text.set_visible(False)
        return
    line.set_visible(True)
    text.set_visible(True)
    line.set_data(geometry['bar_x'], geometry['bar_y'])
    text.set_position(geometry['text_xy'])
    text.set_text(geometry['label'])


def build_section_scalebar(ax, color='white', fontsize=8, geometry=None, zorder=6, animated=False):
    """Add a scale bar (a Line2D + Text, both real vector artists in `ax`'s
    own data coordinates — not a fixed-pixel overlay) to `ax`, sized and
    positioned from `geometry` (see compute_scalebar_geometry), or freshly
    computed from `ax`'s current view if omitted. Returns (line, text): a
    caller that will keep tracking a live, changing view hangs onto them
    and re-applies new geometry via apply_scalebar_geometry; a one-off
    static export can just discard the return value.

    animated=True (the interactive viewer's own live scale bar; left False
    — the default — for a one-off static export, which has no zoom-preview
    bitmaps or cached images to sit on top of) marks both artists animated
    and gives them a very high zorder, so draw_animated_overlays can
    re-stamp them on top of *any* current panel content — a live scatter,
    a zoom-preview bitmap, or a cached home-view PNG — instead of the bar
    only ever being correct when the panel happens to be showing its real,
    freshly-drawn scatter. Without this, whatever the bar looked like at
    the moment a bitmap snapshot was taken got baked into that bitmap's
    own pixels, stretching/shifting as the bitmap was zoomed, and
    sometimes doubling up with the live bar (see snapshot_axes_region's
    own note on why it hides this pair before capturing)."""
    if geometry is None:
        geometry = compute_scalebar_geometry(ax)
    line = Line2D([0, 0], [0, 0], color=color, linewidth=2.5, solid_capstyle='butt', zorder=zorder)
    ax.add_line(line)
    text = ax.text(0, 0, '', color=color, ha='center', va='bottom', fontsize=fontsize, zorder=zorder)
    if animated:
        line.set_animated(True)
        text.set_animated(True)
    apply_scalebar_geometry(line, text, geometry)
    return line, text


def visible_points_only(offsets, facecolors, sizes, xlim, ylim):
    """`(offsets, facecolors, sizes)` filtered down to just the points inside
    `xlim`/`ylim` (each may be given in either order — a decreasing ylim,
    e.g. every section panel's own inverted y-axis, is handled the same as
    an increasing one). `facecolors`/`sizes` are matplotlib scatter's own
    get_facecolors()/get_sizes() — each is either one row that broadcasts to
    every point (left untouched) or one row per point (filtered the same as
    `offsets`).

    Exists for render_section_maps_export_figure: a background_artist
    always holds a *whole section's* points — many more than are inside the
    current, possibly zoomed-in view — and relying on the new axes' own
    clip-path to hide the rest at render/save time turned out not to be
    enough. It works in every renderer this was actually drawn in (Agg for
    the PNG, a browser for the SVG), but at least one common SVG viewer the
    result also gets opened in doesn't apply a matplotlib clip-path to a
    large scatter's individual points, which showed as this panel's own
    off-screen cells bleeding into whatever sits next to it in the grid.
    Filtering the *data* here removes the dependence on any renderer
    honoring the clip at all, and, incidentally, makes for a much smaller
    file than shipping a whole section's worth of points."""
    x0, x1 = sorted(xlim)
    y0, y1 = sorted(ylim)
    visible = ((offsets[:, 0] >= x0) & (offsets[:, 0] <= x1)
               & (offsets[:, 1] >= y0) & (offsets[:, 1] <= y1))
    facecolors = facecolors if len(facecolors) == 1 else facecolors[visible]
    sizes = sizes if len(sizes) == 1 else sizes[visible]
    return offsets[visible], facecolors, sizes


def unique_export_paths(target_dir, stem, extensions):
    """{ext: path} for a new export named `stem` (no leading '.' on the
    extensions), without overwriting an earlier one: the first save uses
    `stem` bare; every save after that appends the lowest unused '_N'
    suffix, checked across every extension in `extensions` together so a
    PNG and SVG saved in the same call always share one suffix. A suffix
    the user appended themselves after the number (e.g. renaming a copy to
    '..._2_final.png') doesn't confuse the numbering — only the digits
    immediately after '_N' are read; anything past that is ignored."""
    base_exists = any((target_dir / f'{stem}{ext}').exists() for ext in extensions)
    if not base_exists:
        return {ext: target_dir / f'{stem}{ext}' for ext in extensions}
    pattern = re.compile(rf'^{re.escape(stem)}_(\d+)(?:_.*)?$')
    highest = 0
    for ext in extensions:
        for existing in target_dir.glob(f'{stem}_*{ext}'):
            m = pattern.match(existing.stem)
            if m:
                highest = max(highest, int(m.group(1)))
    n = highest + 1
    return {ext: target_dir / f'{stem}_{n}{ext}' for ext in extensions}






def render_section_home_view_png(xs, ys, colors, is_gray, home_xlim, home_ylim, long_edge_in, dot_size_scale=1.0):
    """Off-screen render of one section panel's *complete* (unfiltered)
    background point set at its home extent, colored per `colors` (gray
    cells drawn first/underneath, same convention as apply_panel_colors_
    with_gray_behind — `is_gray` decides that draw order the same way).

    `dot_size_scale` (an AREA multiplier, applied directly to `s`) exists
    because `long_edge_in` is deliberately bigger than this panel's actual
    on-screen size (see show_section_home_cache_or_scatter's own comment,
    ZOOM_BITMAP_ONLY_MAX_MULTIPLIER safety margin) — but matplotlib's `s`
    is a fixed physical size in points^2, independent of the figure's own
    size, so a cell dot rendered at the *same* `s` into a bigger-than-
    needed figure ends up covering a *smaller fraction* of it, and shows
    up smaller still once that oversized image is shrunk back down via
    imshow to the panel's real on-screen size. Caught directly (not just
    suspected) via the physical-unit math: the shrink factor works out to
    exactly 1/(long_edge_in / this panel's own unscaled on-screen size) —
    i.e. 1/ZOOM_BITMAP_ONLY_MAX_MULTIPLIER with the current formula — the
    same relationship the caller uses to compute `dot_size_scale` as that
    ratio squared, canceling it back out. Left at 1.0 (no correction) by
    any other caller that renders at a panel's own true size 1:1.

    Deliberately does *not* invert the y-axis the way the live section
    panels do (see sec_ax.invert_yaxis()'s own comment on the atlas's y-
    increases-downward convention) — this produces a plain, naturally-
    oriented image (like the raw xs/ys, not pre-flipped), the same way a
    real scatter's own raw data isn't pre-flipped either; the *live*
    axes' own inversion is what correctly orients it once displayed via
    imshow (see show_section_home_cache_or_scatter). Inverting here too
    would flip it twice — this image would come out correctly oriented as
    its own standalone picture, but upside down once shown on the already-
    inverted live axes, which is exactly what happened before this
    comment existed.

    Uses a fresh, throwaway Figure+FigureCanvasAgg at a fixed DPI
    (SECTION_HOME_CACHE_DPI) and a `long_edge_in` the *caller* computes
    from this panel's own current on-screen size (see show_section_home_
    cache_or_scatter) — entirely independent of the live interactive
    window's own rendering, so the result looks identical regardless of
    which machine's screen happened to generate it (see section_home_
    cache_path's own docstring), while still being sized appropriately for
    however large this panel actually needs to be. Returns an RGBA uint8
    array, ready for Image.fromarray(...).save(...) or direct imshow use."""
    x0, x1 = home_xlim
    y0, y1 = home_ylim
    width_data, height_data = abs(x1 - x0), abs(y1 - y0)
    if width_data <= 0 or height_data <= 0:
        figsize = (long_edge_in, long_edge_in)
    elif width_data >= height_data:
        figsize = (long_edge_in, long_edge_in * height_data / width_data)
    else:
        figsize = (long_edge_in * width_data / height_data, long_edge_in)
    fig = Figure(figsize=figsize, dpi=SECTION_HOME_CACHE_DPI)
    FigureCanvasAgg(fig)
    fig.patch.set_facecolor(SECTION_PANEL_FACECOLOR)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(SECTION_PANEL_FACECOLOR)
    order = np.argsort(~np.asarray(is_gray), kind='stable')
    ax.scatter(xs[order], ys[order], c=np.asarray(colors, dtype=object)[order],
               s=SECTION_BACKGROUND_BASE_SIZE * dot_size_scale, alpha=SECTION_POINT_ALPHA, linewidths=0)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba()).copy()


# Taxonomy levels a single-section window's Groups view can be colored by.
# 'class' is the image generate_and_cache_section_image() already produces;
# the others are rendered on demand by render_section_level_array().
SECTION_COLOR_LEVELS = ('class', 'subclass', 'supertype', 'cluster')










def roi_selection_unchanged(path, rois, whole_sections):
    """True if the roi_coords.csv at `path` already records exactly this
    session's `rois` and `whole_sections` (order-independent) — lets the
    caller skip re-prompting/rewriting roi_coords.csv and roi_map.png when
    nothing about the ROI/whole-section selection actually changed since
    they were last saved (e.g. reusing a cached run's exact selection).
    False if `path` doesn't exist or can't be parsed as an ROI CSV."""
    path = Path(path)
    if not path.exists():
        return False
    try:
        existing_rois, existing_whole_sections = load_rois_csv(path)
    except Exception:
        return False

    def roi_key(roi):
        return (
            str(roi['section']),
            round(float(roi['x_min']), 6), round(float(roi['x_max']), 6),
            round(float(roi['y_min']), 6), round(float(roi['y_max']), 6),
        )

    try:
        rois_match = {roi_key(r) for r in rois} == {roi_key(r) for r in existing_rois}
    except (TypeError, ValueError):
        return False
    whole_sections_match = set(map(str, whole_sections)) == set(map(str, existing_whole_sections))
    return rois_match and whole_sections_match




def required_cached_run_files(run_folder):
    """The 3 files a previously computed run can draw on for the "load
    existing" startup path to skip straight to the interactive UMAP viewer:
    the run's own umap_coords.csv (embedding + metadata) and roi_coords.csv
    (ROI/whole-section rectangles), plus the full processed AnnData Step 8
    wrote to CACHE_DIR (named after the run folder itself, not stored inside
    it). Returns a dict of label -> Path; callers check .exists() on each.

    Only 'UMAP coordinates' is actually required. 'ROI coordinates' is never
    written at all for a whole-brain run, and 'Processed AnnData' can be
    rebuilt from the other two (see rebuild_processed_adata) — so the
    startup panel treats both as recoverable rather than blocking on them."""
    run_folder = Path(run_folder)
    if not abc_project.has_umap(run_folder):
        return {'Project file': run_folder / abc_project.PROJECT_FILE}           # (a project without a UMAP is only its project file, and its ROIs if any)
    return {
        'UMAP coordinates': run_folder / 'umap_coords.csv',
        'ROI coordinates': run_folder / 'roi_coords.csv',
        'Processed AnnData': CACHE_DIR / f'{run_folder.name}_processed.h5ad',
    }


def leiden_labels_from_csv(series):
    """Leiden labels read back out of umap_coords.csv, normalized to plain
    integer-like strings ('0', '1', ...).

    A CSV round trip doesn't preserve these as written: pandas re-infers the
    column, so '0'/'1' come back as int64 — or, if any row is blank (a cell
    present in the CSV but not in the subset being loaded), as float64,
    whose str() is '0.0'. Left alone that quietly produces a second set of
    category names that don't match the '0'/'1' a freshly computed run
    yields, so the same cluster would look like two different ones
    depending on whether the run came from cache. Blanks stay None
    (unknown) rather than becoming the string 'nan'."""
    numeric = pd.to_numeric(series, errors='coerce')
    out = series.astype(object).where(series.notna(), None)
    is_num = numeric.notna()
    out[is_num] = numeric[is_num].astype('int64').astype(str)
    return out.to_numpy()


def rebuild_processed_adata(run_folder, adata_backed, processed_h5ad_path):
    """Rebuilds (and re-caches) a run's processed AnnData from its saved
    umap_coords.csv plus the raw backed h5ad, for when the .h5ad Step 8
    wrote has been deleted but the run folder itself is intact.

    Everything needed is already in the CSV — the exact cell IDs the run
    settled on, their ABC taxonomy columns, the Leiden clustering, and the
    UMAP embedding itself — so none of the expensive work (subsample,
    preprocessing, PCA, neighbors, UMAP, Leiden) has to happen again, and
    the section/ROI picker doesn't need to be reopened just to re-derive a
    selection the CSV already records. Mirrors what the main pipeline's own
    found_existing_umap branch does with a cached CSV; kept separate rather
    than shared with it because that branch is threaded through the full
    Step 0-8 flow (subsample caps, cell-type filters, section selection)
    that this path deliberately skips.

    Returns the rebuilt AnnData, having also written it back to
    `processed_h5ad_path` so the next load is a plain read again."""
    csv_path = Path(run_folder) / 'umap_coords.csv'
    print(f"Rebuilding processed AnnData from {csv_path}...")
    # converters={0: str} keeps the ~19-digit cell-ID index as text — as
    # float64 pandas silently rounds off its trailing digits.
    umap_coords = pd.read_csv(csv_path, index_col=0, converters={0: str})
    keep_mask = adata_backed.obs.index.isin(umap_coords.index)
    n_found = int(keep_mask.sum())
    if n_found == 0:
        raise RuntimeError(
            f"None of the {len(umap_coords)} cell IDs in {csv_path} were found in the raw dataset."
        )
    if n_found < len(umap_coords):
        print(f"Warning: only {n_found} of {len(umap_coords)} cells from {csv_path} "
              "were found in the raw dataset; rebuilding with those.")
    # Reads only the matched cells off disk, rather than materializing the
    # full multi-million-cell dataset just to discard most of it.
    adata = materialize_subset(adata_backed, keep_mask)
    for meta_col in ('class', 'subclass', 'supertype', 'cluster', SECTION_COL):
        if meta_col in umap_coords.columns and meta_col not in adata.obs.columns:
            adata.obs[meta_col] = umap_coords.loc[adata.obs.index, meta_col].values
    if LEIDEN_KEY in umap_coords.columns:
        adata.obs[LEIDEN_KEY] = pd.Categorical(
            leiden_labels_from_csv(umap_coords.loc[adata.obs.index, LEIDEN_KEY])
        )
    adata.obsm['X_umap'] = umap_coords.loc[adata.obs.index, ['UMAP1', 'UMAP2']].to_numpy()
    # Same as the main pipeline's cache-hit branch: adata.X here is still
    # the untouched raw counts as read from disk (materialize_subset never
    # normalizes or transforms), which is exactly what the viewer's Gene
    # mode wants from layers['counts'].
    adata.layers['counts'] = adata.X.copy()
    processed_h5ad_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing rebuilt AnnData to {processed_h5ad_path}...")
    adata.write_h5ad(processed_h5ad_path)
    print(f"Rebuilt processed AnnData for {adata.n_obs} cells.")
    return adata


def show_error_dialog(title, message):
    """Best-effort modal error dialog, falling back to the console if Tk
    isn't usable. Same withdrawn-root/topmost pattern as this file's other
    tkinter prompts (offer_leiden_backfill's own askyesno, the ROI file
    pickers), with one addition: an already-running Tk root is reused rather
    than a second one created alongside it, since these can fire while a
    matplotlib window (itself a Tk app) is open."""
    print(f"{title}: {message}")
    try:
        import tkinter as tk
        from tkinter import messagebox
        existing_root = getattr(tk, '_default_root', None)
        root = existing_root
        if root is None:
            root = tk.Tk()
            root.withdraw()
        root.attributes('-topmost', True)
        messagebox.showerror(title, message)
        if existing_root is None:
            root.destroy()
    except Exception as e:
        print(f"(Could not show the error dialog: {e})")


def leiden_failure_message(exc):
    """User-facing explanation for a failed Leiden run, with an install hint
    when the cause is a missing optional package.

    scanpy's Leiden support needs igraph (and, for the non-'igraph' flavors,
    leidenalg) — neither is a hard dependency of scanpy itself, so a missing
    one is by far the most likely reason this fails, and it's fixable in one
    command. ModuleNotFoundError carries the module in .name; a plain
    ImportError raised from inside scanpy may not, so the message text is
    checked as a fallback."""
    lines = [f"Leiden clustering could not run:", "", f"{type(exc).__name__}: {exc}", ""]
    if isinstance(exc, ImportError):
        missing = getattr(exc, 'name', None)
        if not missing:
            text = str(exc).lower()
            missing = next((pkg for pkg in ('leidenalg', 'igraph') if pkg in text), None)
        if missing:
            lines += [f"This needs the optional '{missing}' package. Install it with:",
                       "", f"    pip install {missing}", ""]
        else:
            lines += ["This needs scanpy's optional clustering packages. Install them with:",
                       "", "    pip install igraph leidenalg", ""]
    lines.append("Everything else is unaffected: the UMAP and all ABC taxonomy levels "
                  "(class/subclass/supertype/cluster) still work. Only the Leiden level "
                  "will be unavailable for this run.")
    return "\n".join(lines)


def prompt_deg_comparison_mode(parent_window, title, level_name):
    """Modal dialog for the interactive UMAP viewer's 'Export DEGs' button,
    offering two comparison modes: each specified ID versus every other
    cell (export_de_genes' original, still-default behavior — 'others'),
    or each specified ID versus one particular reference ID's own cells
    only ('reference').

    Returns {'mode': 'others' or 'reference', 'reference': int or None}
    once OK is clicked with valid input, or None if cancelled (Cancel, the
    window's own close button, or Escape).

    A real Toplevel, not tkinter.messagebox — this needs radio buttons and
    a text entry, which messagebox doesn't offer — parented to and modal
    over `parent_window` (the viewer's own Tk window, already running its
    own mainloop via show_figure_blocking) via transient()/grab_set()/
    wait_window(), the standard pattern for a modal child dialog. This
    deliberately does *not* create a second, competing Tk() root the way
    this file's earlier tkinter.messagebox-based prompts do (those run
    before any window/mainloop exists yet; this one has to coexist with
    one that's already running)."""
    import tkinter as tk
    from tkinter import messagebox

    result = {'value': None}
    dialog = tk.Toplevel(parent_window)
    dialog.title(title)
    dialog.transient(parent_window)
    dialog.resizable(False, False)

    mode_var = tk.StringVar(value='others')
    level_word = level_name.lower()
    # 'es' for names ending in a sibilant ('class'/'subclass' ->
    # classes/subclasses, the standard English rule for words ending in
    # s/x/z/ch/sh), plain 's' otherwise (supertype, cluster, leiden) — same
    # rule mode_display_label's own 'Single Subclass' branch uses, so this
    # reads correctly for every level without a per-level special case.
    level_word_plural = level_word + ('es' if level_word.endswith('s') else 's')
    tk.Radiobutton(dialog, text=f"...selected {level_word_plural} vs others", variable=mode_var,
                   value='others', anchor='w').pack(fill='x', padx=12, pady=(12, 2))
    tk.Radiobutton(dialog, text=f"...selected {level_word_plural} versus reference", variable=mode_var,
                   value='reference', anchor='w').pack(fill='x', padx=12, pady=(2, 2))

    ref_frame = tk.Frame(dialog)
    tk.Label(ref_frame, text=f"Reference {level_word} ID:").pack(side='left')
    ref_entry = tk.Entry(ref_frame, width=10)
    ref_entry.pack(side='left', padx=(6, 0))
    ref_frame.pack(fill='x', padx=32, pady=(0, 12))

    def update_ref_entry_state(*_args):
        ref_entry.configure(state='normal' if mode_var.get() == 'reference' else 'disabled')

    mode_var.trace_add('write', update_ref_entry_state)
    update_ref_entry_state()  # starts disabled — 'others' is the default selection above

    button_frame = tk.Frame(dialog)
    button_frame.pack(pady=(0, 12))

    def on_ok():
        mode = mode_var.get()
        if mode == 'reference':
            raw = ref_entry.get().strip()
            try:
                reference = int(raw)
            except ValueError:
                # Left open rather than destroyed — same "fix it and try
                # again" flow a real form validation error implies, instead
                # of silently discarding what was already typed elsewhere
                # in the dialog.
                messagebox.showerror("Invalid reference ID", f"'{raw}' is not a valid integer ID.", parent=dialog)
                return
            result['value'] = {'mode': 'reference', 'reference': reference}
        else:
            result['value'] = {'mode': 'others', 'reference': None}
        dialog.destroy()

    def on_cancel():
        result['value'] = None
        dialog.destroy()

    tk.Button(button_frame, text="OK", command=on_ok, width=10, default='active').pack(side='left', padx=6)
    tk.Button(button_frame, text="Cancel", command=on_cancel, width=10).pack(side='left', padx=6)

    dialog.protocol("WM_DELETE_WINDOW", on_cancel)
    dialog.bind('<Return>', lambda event: on_ok())
    dialog.bind('<Escape>', lambda event: on_cancel())

    # Centered over the parent window rather than wherever Tk happens to
    # place a fresh Toplevel (typically the screen's own top-left corner).
    # update_idletasks() first — the dialog's own winfo_width()/height()
    # only reflect its real, laid-out size once the window manager has
    # actually sized the widgets above, not immediately after creating them.
    dialog.update_idletasks()
    px, py = parent_window.winfo_rootx(), parent_window.winfo_rooty()
    pw, ph = parent_window.winfo_width(), parent_window.winfo_height()
    dw, dh = dialog.winfo_width(), dialog.winfo_height()
    dialog.geometry(f"+{px + max(0, (pw - dw) // 2)}+{py + max(0, (ph - dh) // 2)}")

    dialog.grab_set()
    dialog.focus_set()
    dialog.wait_window()  # blocks here (this window's own local event loop) until destroy() above
    return result['value']


# ---------------------------------------------------------------------------
# UMAP parameter variants
# ---------------------------------------------------------------------------
# The embedding a run was first computed with (Step 6) uses scanpy's usual
# parameters. The viewer's Settings dialog can recompute it with different
# ones to compare how they change the picture; each (n_neighbors, min_dist,
# n_pcs) combination is cached in the run folder under a name that carries
# the parameters, so trying a combination a second time is instant.
UMAP_PARAM_DEFAULTS = {'n_neighbors': 15, 'min_dist': 0.5, 'n_pcs': 50}
UMAP_PARAM_MAX_PCS = 100         # Step 6 stores 50; a variant asking for more recomputes PCA (see principal_components_for_variant)
UMAP_PARAMS_FILE_NAME = 'umap_params.json'
UMAP_VARIANT_ESTIMATE_KEY = 'umap_variant_s'
UMAP_VARIANT_METADATA_COLS = ('class', 'subclass', 'supertype', 'cluster', SECTION_COL)


def normalize_umap_params(params):
    """`params` as {'n_neighbors': int, 'min_dist': float, 'n_pcs': int},
    with any missing key taken from UMAP_PARAM_DEFAULTS."""
    merged = dict(UMAP_PARAM_DEFAULTS)
    merged.update({k: v for k, v in dict(params).items() if k in merged})
    return {'n_neighbors': int(merged['n_neighbors']), 'min_dist': float(merged['min_dist']),
            'n_pcs': int(merged['n_pcs'])}


def umap_params_key(params):
    p = normalize_umap_params(params)
    return (p['n_neighbors'], round(p['min_dist'], 6), p['n_pcs'])


def umap_params_are_default(params):
    return umap_params_key(params) == umap_params_key(UMAP_PARAM_DEFAULTS)


def umap_params_token(params):
    """Filename-safe tag for a parameter set ('nn15_md0.5_pc50'); empty for
    the defaults, so everything made with the default embedding keeps the
    names it always had."""
    if umap_params_are_default(params):
        return ''
    p = normalize_umap_params(params)
    return f"nn{p['n_neighbors']}_md{p['min_dist']:g}_pc{p['n_pcs']}"


def umap_params_label(params):
    p = normalize_umap_params(params)
    return f"n_neighbors={p['n_neighbors']}, min_dist={p['min_dist']:g}, n_pcs={p['n_pcs']}"


def umap_variant_csv_path(run_folder, params):
    """Where the embedding for `params` lives: the run's own umap_coords.csv
    for the defaults, umap_coords_<token>.csv otherwise."""
    token = umap_params_token(params)
    name = 'umap_coords.csv' if not token else f'umap_coords_{token}.csv'
    return Path(run_folder) / name


def parse_umap_params(n_neighbors, min_dist, n_pcs, n_cells=None):
    """Validates the Settings dialog's three text fields. Returns
    (params, None) or (None, message for the user)."""
    try:
        nn = int(str(n_neighbors).strip())
    except ValueError:
        return None, "Neighbors (n_neighbors) must be a whole number."
    try:
        md = float(str(min_dist).strip())
    except ValueError:
        return None, "min_dist must be a number between 0 and 1."
    try:
        pcs = int(str(n_pcs).strip())
    except ValueError:
        return None, "n_pcs must be a whole number."
    nn_max = max(2, (n_cells - 1)) if n_cells else None
    if nn < 2 or (nn_max is not None and nn > nn_max):
        return None, f"n_neighbors must be between 2 and {nn_max if nn_max else 'the number of cells'}."
    if not (0.0 <= md <= 1.0):
        return None, "min_dist must be between 0 and 1."
    if not (2 <= pcs <= UMAP_PARAM_MAX_PCS):
        return None, f"n_pcs must be between 2 and {UMAP_PARAM_MAX_PCS}."
    return {'n_neighbors': nn, 'min_dist': md, 'n_pcs': pcs}, None


def read_last_umap_params(run_folder):
    """The parameters the run was last viewed with (UMAP_PARAMS_FILE_NAME in
    the run folder); the defaults when absent or unreadable."""
    try:
        data = json.loads((Path(run_folder) / UMAP_PARAMS_FILE_NAME).read_text(encoding='utf-8'))
        return normalize_umap_params(data)
    except Exception:
        return dict(UMAP_PARAM_DEFAULTS)


def write_last_umap_params(run_folder, params):
    try:
        (Path(run_folder) / UMAP_PARAMS_FILE_NAME).write_text(
            json.dumps(normalize_umap_params(params), indent=2), encoding='utf-8')
    except Exception as e:
        print(f"Warning: could not save the UMAP parameters to {run_folder} ({e}).")


def load_umap_variant(csv_path, obs_names):
    """(embedding as an (n, 2) float array, Leiden labels or None) from a
    variant CSV, aligned to `obs_names`."""
    table = pd.read_csv(csv_path, index_col=0, converters={0: str})
    table = table.reindex(list(obs_names))
    if table['UMAP1'].isna().any() or table['UMAP2'].isna().any():
        raise RuntimeError(f"{csv_path} does not cover every cell in this run.")
    leiden = leiden_labels_from_csv(table[LEIDEN_KEY]) if LEIDEN_KEY in table.columns else None
    return table[['UMAP1', 'UMAP2']].to_numpy(dtype=float), leiden


def principal_components_for_variant(adata, run_folder, n_pcs, report=None):
    """The run's PCA coordinates (cells x >= n_pcs). Normally just the
    X_pca Step 6 stored in the processed AnnData. A run rebuilt from its CSV
    has none, so it is recomputed by replaying Step 5/6 on layers['counts']
    (the same recipe, so the result matches what the run would have had) and
    saved as pca.npy in the run folder for next time."""
    import anndata
    import scanpy as sc
    if 'X_pca' in adata.obsm and adata.obsm['X_pca'].shape[1] >= n_pcs:
        return np.asarray(adata.obsm['X_pca'])
    pca_path = Path(run_folder) / 'pca.npy' if run_folder is not None else None
    if pca_path is not None and pca_path.exists():
        saved = np.load(pca_path)
        if saved.shape[0] == adata.n_obs and saved.shape[1] >= n_pcs:
            return saved
    if 'counts' not in adata.layers:
        raise RuntimeError("This run has no PCA coordinates or raw counts to recompute them from.")
    if report:
        report("Computing principal components")
    counts = adata.layers['counts']
    work = anndata.AnnData(X=counts.copy(), obs=pd.DataFrame(index=adata.obs_names),
                           var=pd.DataFrame(index=adata.var_names))
    from scipy.sparse import issparse
    if issparse(work.X):
        work.X.data = np.nan_to_num(work.X.data, nan=0.0)
    else:
        work.X = np.nan_to_num(work.X, nan=0.0)
    sc.pp.filter_genes(work, min_cells=5)
    sc.pp.normalize_total(work)
    sc.pp.log1p(work)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='zero-centering a sparse array/matrix densifies it',
                                category=UserWarning)
        sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(UMAP_PARAM_MAX_PCS, work.n_vars - 1, work.n_obs - 1), svd_solver='arpack')
    pcs = np.asarray(work.obsm['X_pca'], dtype=np.float32)
    if pca_path is not None:
        try:
            np.save(pca_path, pcs)
        except Exception as e:
            print(f"Warning: could not save {pca_path} ({e}).")
    return pcs


def compute_umap_variant_csv(adata, params, run_folder, report=None):
    """Computes the embedding (and its Leiden clustering, which depends on
    the neighbor graph) for `params` and writes the variant CSV. Returns a
    message about Leiden if it could not be computed, else None. Safe to
    run on a worker thread: touches no GUI state and never modifies
    `adata`."""
    import anndata
    import scanpy as sc
    report = report or (lambda text: None)
    params = normalize_umap_params(params)
    csv_path = umap_variant_csv_path(run_folder, params)
    pcs = principal_components_for_variant(adata, run_folder, params['n_pcs'], report)
    n_pcs = min(params['n_pcs'], pcs.shape[1])
    # Only the PCs go into the working AnnData: the neighbor graph is built
    # from X_pca alone, and this keeps the expression matrix out of it.
    work = anndata.AnnData(np.zeros((adata.n_obs, 1), dtype=np.float32))
    work.obs_names = adata.obs_names
    work.obsm['X_pca'] = np.ascontiguousarray(pcs[:, :n_pcs])
    report("Step 1 of 3: building the neighbor graph")
    sc.pp.neighbors(work, n_neighbors=params['n_neighbors'], use_rep='X_pca')
    report("Step 2 of 3: computing the UMAP embedding")
    sc.tl.umap(work, min_dist=params['min_dist'])
    report("Step 3 of 3: clustering (Leiden)")
    leiden_problem = None
    try:
        sc.tl.leiden(work, key_added=LEIDEN_KEY, flavor='igraph', n_iterations=2,
                     resolution=LEIDEN_RESOLUTION)
    except Exception as e:
        leiden_problem = str(e)

    table = pd.DataFrame(work.obsm['X_umap'], index=adata.obs.index, columns=['UMAP1', 'UMAP2'])
    for col in UMAP_VARIANT_METADATA_COLS:
        if col in adata.obs.columns:
            table[col] = adata.obs[col].values
    if LEIDEN_KEY in work.obs.columns:
        table[LEIDEN_KEY] = work.obs[LEIDEN_KEY].astype(str).values
    tmp_path = csv_path.with_name(csv_path.name + '.tmp')
    table.to_csv(tmp_path)
    os.replace(tmp_path, csv_path)  # a half-written file is never mistaken for a cached variant
    return leiden_problem


_UMAP_VARIANT_JOBS = {}


def ensure_umap_variant_csv(adata, params, run_folder):
    """Makes sure the variant CSV for `params` exists, computing it behind a
    progress window if it doesn't. Returns (csv_path or None, message): None
    when the user cancelled or the computation failed. A cancelled
    computation keeps running in the background and still lands in the
    cache, so asking again shortly afterwards just finds the result."""
    csv_path = umap_variant_csv_path(run_folder, params)
    if csv_path.exists():
        return csv_path, None
    job = _UMAP_VARIANT_JOBS.get(csv_path)
    if job is None or job['done']:
        job = {'status': 'Starting', 'error': None, 'leiden_problem': None, 'done': False}

        def worker():
            try:
                job['leiden_problem'] = compute_umap_variant_csv(
                    adata, params, run_folder, lambda text: job.update(status=text))
            except Exception as e:
                job['error'] = e
            finally:
                job['done'] = True

        _UMAP_VARIANT_JOBS[csv_path] = job
        # pynndescent calls gc.collect() itself mid-computation. That would
        # run on the worker thread, which must not be the one to free the
        # leftover Tk objects of already-closed windows (Tkinter raises "main
        # thread is not in main loop"), and the viewer keeps automatic gc off
        # so such garbage is still waiting. Collected here, on the Tk thread,
        # first.
        gc_started = time.perf_counter()
        gc.collect()
        print(f"[umap-variant] pre-collected garbage on the main thread in "
              f"{time.perf_counter() - gc_started:.2f} s.", flush=True)
        threading.Thread(target=worker, daemon=True).start()

    estimate_key = f"{UMAP_VARIANT_ESTIMATE_KEY}_{adata.n_obs}"
    estimate = load_recorded_duration(estimate_key)
    started = time.perf_counter()
    restore_cursors = set_wait_cursor_on_open_figures()
    proc_fig, cancel_flag, set_progress = show_processing_dialog(
        f"Computing the UMAP with {umap_params_label(params)}.\n"
        "The first calculation in a session takes longest; results are cached for next time.")
    pulse = 0.0
    try:
        set_progress(0.0)
        while not job['done'] and not cancel_flag['cancelled']:
            elapsed = time.perf_counter() - started
            if estimate and elapsed < estimate:
                proc_fig._set_status(f"{job['status']}…  ({elapsed:.0f} s, about {estimate:.0f} s expected)")
                set_progress(min(0.95, elapsed / estimate))
            else:
                proc_fig._set_status(f"{job['status']}…  ({elapsed:.0f} s)")
                pulse += 0.04
                proc_fig._set_pulse(pulse)
            proc_fig.canvas.flush_events()
            time.sleep(0.05)
    finally:
        close_processing_dialog(proc_fig)
        restore_cursors()
    if not job['done']:
        return None, "UMAP recalculation cancelled (it will finish in the background and be cached)."
    if job['error'] is not None:
        return None, f"UMAP recalculation failed: {job['error']}"
    record_load_duration(estimate_key, time.perf_counter() - started)
    return csv_path, (f"Leiden clustering unavailable: {job['leiden_problem']}" if job['leiden_problem'] else None)


_umap_warmup_started = False




def build_viewer_settings_tab(tab, current_settings, umap_params=None, n_cells=None):
    """The "Project browser" tab of the Settings window: hover-highlighting on/off, the multi-gene overlay's genes-4-6 marker shape (GENE456_SHAPE_REGISTRY),
    the overall UMAP dot size (UMAP_CIRCLE_SIZE_REGISTRY), how groups are colored, whether these choices persist across sessions (see load_viewer_settings/
    save_viewer_settings's own docstrings for that setting's own asymmetric persistence rule) and, when `umap_params` (the run's current n_neighbors/min_dist/
    n_pcs) is given, those three, validated against `n_cells`. Fills the frame `tab`; returns read(), which gives the chosen settings as a dict (with
    'umap_params' when they were offered) or raises ValueError with the message to show."""
    import tkinter as tk

    hover_var = tk.BooleanVar(tab, value=current_settings['hover_enabled'])
    tk.Checkbutton(tab, text="Hover highlight", variable=hover_var, anchor='w').pack(
        fill='x', padx=12, pady=(12, 6))

    size_frame = tk.Frame(tab)
    tk.Label(size_frame, text="UMAP circle size:").pack(side='left')
    size_var = tk.StringVar(tab, value=current_settings.get('umap_circle_size', UMAP_CIRCLE_SIZE_DEFAULT))
    tk.OptionMenu(size_frame, size_var, *UMAP_CIRCLE_SIZE_REGISTRY.keys()).pack(side='left', padx=(6, 0))
    size_frame.pack(fill='x', padx=12, pady=(0, 6))

    shape_frame = tk.Frame(tab)
    tk.Label(shape_frame, text="Genes 4-6 shape:").pack(side='left')
    shape_var = tk.StringVar(tab, value=current_settings.get('gene456_shape', GENE456_SHAPE_DEFAULT))
    tk.OptionMenu(shape_frame, shape_var, *GENE456_SHAPE_REGISTRY.keys()).pack(side='left', padx=(6, 0))
    shape_frame.pack(fill='x', padx=12, pady=(0, 6))

    colors_var = tk.StringVar(tab, value=current_settings.get('group_colors', 'cycle'))
    colors_frame = tk.LabelFrame(tab, text="Group colors (when there are more groups than looks)")
    for value, text in (('cycle', "Cycle: every group colored; past the 400th largest the colors and shapes repeat (hover tells which)"),
                        ('lineage', "Lineage: a group's color comes from its parent group, so related groups look alike (experimental)")):
        tk.Radiobutton(colors_frame, text=text, value=value, variable=colors_var, anchor='w', justify='left').pack(fill='x', padx=6)
    colors_frame.pack(fill='x', padx=12, pady=(0, 6))

    other_var = tk.BooleanVar(tab, value=bool(current_settings.get('colorize_other_cells', True)))
    tk.Checkbutton(tab, text="Colorize non-neurons in neuron project and neurons in non-neuron project", variable=other_var, anchor='w', justify='left',
                   wraplength=420).pack(fill='x', padx=12, pady=(0, 6))

    persist_var = tk.BooleanVar(tab, value=current_settings['persist_settings'])
    tk.Checkbutton(
        tab, text="Remember these settings across sessions", variable=persist_var, anchor='w',
    ).pack(fill='x', padx=12, pady=(0, 12))

    umap_vars = {}
    if umap_params is not None:
        umap_frame = tk.LabelFrame(tab, text="UMAP calculation (kept with this run)")
        field_rows = (
            ('n_neighbors', "Neighbors (n_neighbors):", "2 or more; larger favors broad structure"),
            ('min_dist', "Minimum distance (min_dist):", "0 to 1; smaller packs cells tighter"),
            ('n_pcs', "Principal components (n_pcs):", f"2 to {UMAP_PARAM_MAX_PCS}"),
        )
        for row, (key, label, hint) in enumerate(field_rows):
            tk.Label(umap_frame, text=label, anchor='w').grid(row=row, column=0, sticky='w', padx=(8, 4), pady=2)
            umap_vars[key] = tk.StringVar(tab, value=f"{umap_params[key]:g}")
            tk.Entry(umap_frame, textvariable=umap_vars[key], width=7).grid(row=row, column=1, pady=2)
            tk.Label(umap_frame, text=f"default {UMAP_PARAM_DEFAULTS[key]:g} - {hint}", fg='gray35',
                     anchor='w').grid(row=row, column=2, sticky='w', padx=(6, 8), pady=2)
        tk.Label(umap_frame, justify='left', fg='gray35', anchor='w',
                 text="Changing these recalculates the UMAP (about 30 s the first time; each\n"
                      "combination is cached with the run, so revisiting one is instant)."
                 ).grid(row=3, column=0, columnspan=3, sticky='w', padx=8, pady=(4, 2))

        def reset_umap_fields():
            for key, var in umap_vars.items():
                var.set(f"{UMAP_PARAM_DEFAULTS[key]:g}")

        tk.Button(umap_frame, text="Reset to defaults", command=reset_umap_fields).grid(
            row=4, column=0, columnspan=3, sticky='w', padx=8, pady=(2, 8))
        umap_frame.pack(fill='x', padx=12, pady=(0, 12))

    def read():
        chosen_umap_params = None
        if umap_vars:
            chosen_umap_params, problem = parse_umap_params(
                umap_vars['n_neighbors'].get(), umap_vars['min_dist'].get(), umap_vars['n_pcs'].get(), n_cells)
            if problem:
                raise ValueError(problem)
        value = {
            'hover_enabled': hover_var.get(),
            'gene456_shape': shape_var.get(),
            'umap_circle_size': size_var.get(),
            'group_colors': colors_var.get(),
            'colorize_other_cells': other_var.get(),
            'persist_settings': persist_var.get(),
        }
        if chosen_umap_params is not None:
            value['umap_params'] = chosen_umap_params
        return value

    return read


def prompt_viewer_settings_dialog(parent_window, current_settings, umap_params=None, n_cells=None, abc_cache=None, tab="Project browser"):
    """The Settings window of the launcher and of the project browser (the interactive viewer): one dialog with two tabs, "Common" (build_common_settings_tab:
    bregma and the coordinates, applied and saved here) and "Project browser" (build_viewer_settings_tab), opened on `tab`. Same modal-child pattern as
    prompt_deg_comparison_mode — see its own docstring for why a real Toplevel rather than tkinter.messagebox, and why transient()/grab_set()/wait_window()
    rather than a second, competing Tk() root; `parent_window` None (the launcher has no window to hang it on) makes a hidden root of its own.

    `umap_params`/`n_cells`: see build_viewer_settings_tab (they belong to the run, so the caller must keep 'umap_params' out of the persisted settings).
    `abc_cache`: when given, a change of the common settings starts the read of the sections' coordinates again.
    `current_settings` is read but never mutated — the caller decides what to do with the result. Returns the project browser's settings dict once OK is
    clicked (the common ones are already in COORD_SETTINGS and saved), or None if cancelled (Cancel, the window's own close button, or Escape)."""
    import tkinter as tk
    from tkinter import ttk, messagebox

    own_root = None
    if parent_window is None:
        own_root = parent_window = tk.Tk()
        own_root.withdraw()
        own_root.attributes('-topmost', True)
    result = {'value': None}
    dialog = tk.Toplevel(parent_window)
    dialog.title("Settings")
    if own_root is not None:
        dialog.attributes('-topmost', True)
    else:
        dialog.transient(parent_window)
    dialog.resizable(False, False)

    notebook = ttk.Notebook(dialog)
    notebook.pack(fill='both', padx=8, pady=(8, 0))
    common_tab, viewer_tab = tk.Frame(notebook), tk.Frame(notebook)
    notebook.add(common_tab, text="Common")
    notebook.add(viewer_tab, text="Project browser")
    read_common = build_common_settings_tab(common_tab)
    read_viewer = build_viewer_settings_tab(viewer_tab, current_settings, umap_params, n_cells)
    notebook.select(common_tab if tab == "Common" else viewer_tab)

    button_frame = tk.Frame(dialog)
    button_frame.pack(pady=12)

    def on_ok():
        for tab_frame, read in ((common_tab, read_common), (viewer_tab, read_viewer)):
            try:
                result[tab_frame] = read()
            except ValueError as problem:
                notebook.select(tab_frame)
                messagebox.showerror("Settings", str(problem), parent=dialog)
                return
        result['value'] = result.pop(viewer_tab)
        COORD_SETTINGS.update(result.pop(common_tab))
        dialog.destroy()

    def on_cancel():
        result['value'] = None
        dialog.destroy()

    tk.Button(button_frame, text="OK", command=on_ok, width=10, default='active').pack(side='left', padx=6)
    tk.Button(button_frame, text="Cancel", command=on_cancel, width=10).pack(side='left', padx=6)

    dialog.protocol("WM_DELETE_WINDOW", on_cancel)
    dialog.bind('<Return>', lambda event: on_ok())
    dialog.bind('<Escape>', lambda event: on_cancel())

    # Centered over the parent window — see prompt_deg_comparison_mode's
    # own identical block for why update_idletasks() has to run first.
    dialog.update_idletasks()
    if own_root is None:
        px, py = parent_window.winfo_rootx(), parent_window.winfo_rooty()
        pw, ph = parent_window.winfo_width(), parent_window.winfo_height()
        dw, dh = dialog.winfo_width(), dialog.winfo_height()
        dialog.geometry(f"+{px + max(0, (pw - dw) // 2)}+{py + max(0, (ph - dh) // 2)}")

    dialog.grab_set()
    dialog.focus_set()
    try:
        dialog.wait_window()
        if result['value'] is not None:
            save_coordinate_settings()
            if abc_cache is not None:
                start_section_coordinate_load(abc_cache)
        return result['value']
    finally:
        if own_root is not None:
            # The window's Tk variables are freed here, on this thread, while Tk is still up: left to the garbage collector they can be freed later
            # on a background thread, which prints "main thread is not in main loop" once for each of them.
            read_common = read_viewer = notebook = common_tab = viewer_tab = dialog = None
            result.clear()
            import gc
            gc.collect()
            own_root.destroy()


TAXONOMY_COLUMNS = ('class', 'subclass', 'supertype', 'cluster')


def backfill_taxonomy_columns(adata, abc_cache, csv_path, processed_h5ad_path=None):
    """For a saved run whose cells carry no class / subclass / supertype / cluster annotations (its pipeline's
    metadata merge did not happen or failed, so neither umap_coords.csv nor the processed AnnData has them, and the
    viewer reports "No 'subclass' column"): fills them in from the atlas's cell metadata, which the spatial cache
    serves in about a second, and writes them back to the CSV and the processed .h5ad so it is a one-time repair.
    Returns True if anything was added. Touches no GUI object, so it can run behind run_with_loading_window."""
    missing = [c for c in TAXONOMY_COLUMNS if c not in adata.obs.columns]
    if not missing:
        return False
    print(f"This saved run has no {', '.join(missing)} annotations; adding them from the atlas cell metadata.")
    meta = load_section_spatial_coords(adata, abc_cache)  # aligned to adata.obs.index
    if meta is None:
        print("Could not read the atlas cell metadata; the viewer will lack those levels.")
        return False
    added = [c for c in missing if c in meta.columns]
    if not added:
        return False
    for column in added:
        adata.obs[column] = meta[column].to_numpy()
    matched = int(meta[added[0]].notna().sum())
    print(f"Added {', '.join(added)} for {matched} of {adata.n_obs} cells.")
    try:
        table = pd.read_csv(csv_path, index_col=0, converters={0: str})
        for column in added:
            table[column] = adata.obs[column].reindex(table.index).to_numpy()
        front = [c for c in ('UMAP1', 'UMAP2') if c in table.columns]
        table = table[front + [c for c in TAXONOMY_COLUMNS if c in table.columns]
                      + [c for c in table.columns if c not in front and c not in TAXONOMY_COLUMNS]]
        table.to_csv(csv_path)
        print(f"Saved the annotations to {csv_path}.")
    except Exception as e:
        print(f"Warning: could not update {csv_path} with the annotations ({e}).")
    try:
        if processed_h5ad_path is not None and Path(processed_h5ad_path).exists():
            adata.write_h5ad(processed_h5ad_path)
            print(f"Saved the annotations to {processed_h5ad_path}.")
    except Exception as e:
        print(f"Warning: could not update {processed_h5ad_path} with the annotations ({e}).")
    return True


def offer_leiden_backfill(adata, csv_path, processed_h5ad_path=None):
    """For a cached run whose umap_coords.csv predates Leiden clustering:
    offer to compute it now and save it back, so the viewer's Leiden level
    works without recomputing the whole run.

    Only the clustering is recomputed — the saved UMAP embedding is left
    exactly as it was, so the layout the user already knows doesn't shift
    underneath them. That still means redoing normalize/log1p/PCA/neighbors
    (Leiden needs a neighbor graph, and cached runs don't store one), which
    is why this asks rather than just doing it; it's a fraction of a full
    run's cost, but not free.

    A no-op returning False if the clustering is already present, if the
    user declines, or if anything fails — every caller carries on with an
    adata that simply has no Leiden column, which the level selector
    already reports as unavailable."""
    if LEIDEN_KEY in adata.obs.columns:
        return False
    message = (
        f"This run's saved data has no Leiden clustering — it was computed before "
        f"Leiden was added (or Leiden could not run at the time).\n\n"
        f"Compute it now for these {adata.n_obs:,} cells and save it back?\n\n"
        f"The existing UMAP layout is kept unchanged; only the clustering is "
        f"computed, which still requires redoing PCA and the neighbor graph and "
        f"may take a few minutes.\n\n"
        f"Choosing No just leaves the Leiden level unavailable for this run."
    )
    try:
        import centered_dialog
        confirmed = centered_dialog.ask_yes_no("Compute Leiden clustering?", message)
    except Exception as e:
        print(f"Could not show the Leiden prompt ({e}); asking on the console instead.")
        confirmed = prompt_yes_no("Compute Leiden clustering for this cached run? (y/n): ")
    if not confirmed:
        print("Skipping Leiden clustering; the Leiden level will be unavailable for this run.")
        return False

    try:
        print("Computing Leiden clustering for this cached run...")
        leiden_start = time.perf_counter()
        # A copy, so nothing here touches the caller's own `adata` — only
        # the resulting label column is copied back at the end.
        work = adata.copy()
        # `adata.X` is NOT reliably raw counts at this point — this is
        # exactly what was producing the reported "Input X contains NaN"
        # failure. For a run loaded from an *existing* processed .h5ad
        # (as opposed to one rebuilt fresh by rebuild_processed_adata,
        # which leaves .X alone), .X is Step 5's own already-normalized/
        # log1p'd/scaled matrix — the main pipeline scales adata.X in
        # place and only *afterward* is it ever written to disk, so a
        # freshly-loaded run's .X is that scaled (z-scored, hence
        # negative-valued) matrix, not raw counts. Re-running
        # normalize_total/log1p (below) on already-scaled data is where
        # the NaN actually came from: log1p(x) for x <= -1 is NaN
        # outright, and scaled data routinely dips well below -1. The
        # *actual* raw counts survive separately in layers['counts'] —
        # both this loaded-h5ad case (see the main pipeline's own
        # `adata.layers['counts'] = adata.X.copy()`, saved before scaling
        # overwrites .X) and the rebuild_processed_adata case (which
        # populates the same layer, redundantly but harmlessly, since
        # its own .X already *is* raw) — so it's always present here.
        if 'counts' in work.layers:
            work.X = work.layers['counts'].copy()
        # Still defensive from here, unlike the main pipeline's own
        # preprocessing (Step 5), which never reaches PCA with either of
        # these problems because it runs sc.pp.filter_cells(min_counts=20)/
        # filter_genes(min_cells=5) first. This path skips that (the whole
        # point is reusing the already-saved UMAP embedding/cell selection
        # unchanged — see the docstring), so `work.X` (now the raw counts
        # from layers['counts'], per above) is unfiltered, with no gene
        # filtering applied at all. Two ways that reaches PCA broken:
        #  - An actual NaN entry in the raw counts (e.g. a corrupted/missing
        #    MERFISH measurement) — sklearn's PCA rejects any NaN outright
        #    ("Input X contains NaN"), which is the exact failure this was
        #    written for.
        #  - A cell with zero total counts — turns normalize_total's
        #    per-cell division into 0/0 = NaN, which is the same failure
        #    one step removed.
        # Both are rare (a handful of cells at most), so zeroing stray NaNs
        # and dropping any now-all-zero cells costs this recompute nothing
        # meaningful, instead of failing Leiden for every cell over a
        # problem in a few.
        from scipy.sparse import issparse
        if issparse(work.X):
            if work.X.nnz and np.isnan(work.X.data).any():
                n_nan = int(np.isnan(work.X.data).sum())
                print(f"Warning: {n_nan} NaN count value(s) found; treating as zero.")
                work.X.data = np.nan_to_num(work.X.data, nan=0.0)
        else:
            n_nan = int(np.isnan(work.X).sum())
            if n_nan:
                print(f"Warning: {n_nan} NaN count value(s) found; treating as zero.")
                work.X = np.nan_to_num(work.X, nan=0.0)
        n_before = work.n_obs
        sc.pp.filter_cells(work, min_counts=1)
        if work.n_obs < n_before:
            print(f"Warning: dropping {n_before - work.n_obs} cell(s) with zero total counts "
                  "(would otherwise turn into NaN during normalization).")
        if work.n_obs == 0:
            raise RuntimeError("every cell has zero total counts; nothing left to cluster")
        sc.pp.normalize_total(work)
        sc.pp.log1p(work)
        n_comps = min(50, work.n_vars - 1, work.n_obs - 1)
        sc.tl.pca(work, n_comps=n_comps, svd_solver='arpack')
        sc.pp.neighbors(work, n_neighbors=15, n_pcs=n_comps)
        sc.tl.leiden(work, key_added=LEIDEN_KEY, flavor='igraph', n_iterations=2, resolution=LEIDEN_RESOLUTION)
        # Assigned by index, not by position (work.obs[...].values) — work
        # can now be a strict subset of adata (see the zero-count filtering
        # above), so a positional assignment would silently misalign every
        # row after the first dropped cell. Cells not in work (dropped
        # above) are left as NaN, same as any other cell Leiden couldn't
        # place.
        adata.obs[LEIDEN_KEY] = pd.Series(np.nan, index=adata.obs.index, dtype=object)
        adata.obs.loc[work.obs.index, LEIDEN_KEY] = work.obs[LEIDEN_KEY].astype(str).values
        del work
        gc.collect()
        print(f"Leiden found {adata.obs[LEIDEN_KEY].nunique()} clusters in "
              f"{format_duration(time.perf_counter() - leiden_start)}.")
    except Exception as e:
        show_error_dialog("Leiden clustering failed", leiden_failure_message(e))
        adata.obs.drop(columns=[LEIDEN_KEY], inplace=True, errors='ignore')
        return False

    # Written back so this is a one-time cost per run, not per session.
    # Best-effort: a failure to save doesn't invalidate the clustering
    # that's already on `adata` and usable in this session.
    try:
        csv_path = Path(csv_path)
        if csv_path.exists():
            umap_coords = pd.read_csv(csv_path, index_col=0, converters={0: str})
            umap_coords[LEIDEN_KEY] = adata.obs[LEIDEN_KEY].astype(str).reindex(umap_coords.index).values
            umap_coords.to_csv(csv_path)
            print(f"Saved Leiden clustering to {csv_path}.")
    except Exception as e:
        print(f"Warning: could not update {csv_path} with the Leiden column ({e}).")
    try:
        if processed_h5ad_path is not None and Path(processed_h5ad_path).exists():
            adata.write_h5ad(processed_h5ad_path)
            print(f"Saved Leiden clustering to {processed_h5ad_path}.")
    except Exception as e:
        print(f"Warning: could not update {processed_h5ad_path} with the Leiden column ({e}).")
    return True


def get_cached_raw_backed_adata(raw_backed_cache, abc_cache):
    """The raw, unfiltered MERFISH h5ad, opened once in backed mode and
    reused across every session in this process — both the 'new_run' path
    (Step 0) and the 'load_existing' path need this exact same read for
    adata_backed (each section's full spatial background in the
    interactive UMAP viewer), and re-opening it fresh every session the
    user starts (see the startup panel's own session loop) was pure
    redundant work: backed mode only reads the header, not X, but it's
    still a real file-system round trip and h5py setup cost each time.
    `raw_backed_cache` is a plain {} the caller keeps alive across
    sessions; this populates 'adata_backed' in it on first use."""
    if 'adata_backed' not in raw_backed_cache:
        print("Step 0: Locating raw MERFISH h5ad file...")
        h5ad_path = abc_cache.get_data_path(
            directory='MERFISH-C57BL6J-638850',
            file_name='C57BL6J-638850/raw'
        )
        print(f"Step 0: Opening {h5ad_path} in backed mode...")
        raw_backed_cache['adata_backed'] = anndata.read_h5ad(h5ad_path, backed='r')
    else:
        print("Step 0: Reusing already-opened raw MERFISH h5ad from an earlier session this run.")
    return raw_backed_cache['adata_backed']


def build_3d_edges(abc_cache, structures, parent=None):
    """The browser view's edge geometry (brain surface, ventricles, white-matter tracts, and the region(s) holding most of the shown cells'
    `structures`, none if there are none) from the CCF annotation volume, or None (the view then just has its cloud).
    Asks before the one-time download of that volume; the first build takes about half a minute, and is cached after
    that. Runs behind a loading window, so call it from the thread that owns the Tk windows."""
    names = load_parcellation_names(abc_cache)
    if names is None:
        return None
    regions = choose_cell_regions(structures) if structures is not None and len(structures) else []
    _store, missing, have_skimage = ccf_edges_missing(regions)
    if missing:
        _path, size, present = metadata_file_status(abc_cache, CCF_VOLUME_DIRECTORY, CCF_VOLUME_FILE_NAME)
        if present is not True:
            try:
                from tkinter import messagebox
                size_text = f"{size / 1e6:.1f} MB" if size else "a small file"
                if not messagebox.askyesno(
                        "Download the CCF annotation volume?",
                        f"The brain edges need the Allen CCF annotation volume, a one-time download of {size_text} "
                        f"into your atlas cache folder. Download it now?", parent=parent):
                    return None
            except Exception as e:
                print(f"Could not ask about the CCF volume download ({e}); no edges.")
                return None
    try:
        return run_with_loading_window(
            "Building the brain edges (the first time takes about half a minute; they are cached after that).",
            lambda: load_ccf_edges(abc_cache, names, regions, want_surfaces=have_skimage))
    except Exception as e:
        print(f"Could not build the brain edges ({e}); the view has its cloud only.")
        return None


# ---- the standalone 3D viewer (opening screen -> "Open 3D Viewer…") ------------------------------------------------
THREED_VIEWER_DIR = CACHE_DIR / '3d_viewer'   # where its pages are written
COORDINATE_LEVEL_NAMES = {'class': 'Class', 'subclass': 'Subclass', 'supertype': 'Supertype', 'cluster': 'Cluster'}


def read_coordinates_for_3d(path):
    """A coordinates_*.csv (the viewer's "Export coordinates" file) arranged for the 3D page: {'groups', 'structures',
    'axis_labels', 'have_ccf', 'title', 'n_cells'}. One group per group_id (named after the level in the file name, e.g.
    coordinates_subclass_268-270.csv, else the taxonomy column the ids come from); a file without group_id is one
    group. Uses the CCF positions when the file has them, else the section positions (which look like stacked
    slices). Raises ValueError for a file that has neither."""
    path = Path(path)
    table = pd.read_csv(path, dtype={'cell_label': str})
    ccf_columns, section_columns = ('ccf_z_mm', 'ccf_x_mm', 'ccf_y_mm'), ('section_x_mm', 'section_z_mm', 'section_y_mm')
    if all(c in table.columns for c in ccf_columns):
        columns, have_ccf = ccf_columns, True   # left-right, front-back, up-down
        axis_labels = ('left - right (CCF z, mm)', 'anterior - posterior (CCF x, mm)', 'dorsal - ventral (CCF y, mm)')
    elif all(c in table.columns for c in section_columns):
        columns, have_ccf = section_columns, False
        axis_labels = ('section x (mm)', 'section position (mm)', 'section y (mm)')
    else:
        raise ValueError("this file has no ccf_x_mm / ccf_y_mm / ccf_z_mm columns (nor section_x_mm / section_y_mm / "
                         "section_z_mm): is it a coordinates_*.csv from the viewer's Export coordinates button?")
    if table.empty:
        raise ValueError("this file has no cells")
    level = None
    named = re.match(r'coordinates_(class|subclass|supertype|cluster)_', path.name)
    if named:
        level = named.group(1)
    elif 'group_id' in table.columns:
        for candidate in COORDINATE_LEVEL_NAMES:
            if candidate in table.columns and set(
                    extract_leading_numeric_id(table[candidate].dropna().drop_duplicates()).dropna().astype(int)
            ) >= set(pd.to_numeric(table['group_id'], errors='coerce').dropna().astype(int)):
                level = candidate
                break
    level_name = COORDINATE_LEVEL_NAMES.get(level, 'Cells')
    hover = None
    if have_ccf and 'ccf_structure' in table.columns:
        section = (table['brain_section_label'].astype(str).str.rsplit('.', n=1).str[-1]
                   if 'brain_section_label' in table.columns else '')
        hover = (table['ccf_structure'].fillna('?').astype(str) + ' | section ' + section
                 + ' | ' + table['cell_label'].astype(str)).to_numpy() if 'cell_label' in table.columns else None
    if 'group_id' in table.columns:
        ids = list(pd.unique(table['group_id'].dropna()))
    else:
        ids = [None]
    groups = []
    for rank, group_id in enumerate(ids):
        rows = table.index[table['group_id'] == group_id] if group_id is not None else table.index
        label = (f"{level_name} {group_id:g} ({len(rows):,})" if isinstance(group_id, (int, float, np.integer, np.floating))
                 else f"{level_name} {group_id} ({len(rows):,})") if group_id is not None else f"{path.stem} ({len(rows):,})"
        groups.append({'label': label, 'color': 'red' if len(ids) == 1 else category_rank_color(rank),
                       'xyz': table.loc[rows, list(columns)].to_numpy(dtype=float),
                       'hover': None if hover is None else hover[np.asarray(rows)]})
    structures = table['ccf_structure'].to_numpy() if have_ccf and 'ccf_structure' in table.columns else None
    return {'groups': groups, 'structures': structures, 'axis_labels': axis_labels, 'have_ccf': have_ccf,
            'title': f"3D viewer: {path.name}", 'n_cells': len(table)}


def build_3d_viewer_page(abc_cache, csv_path=None, parent=None):
    """Writes the interactive 3D page for the cells in `csv_path` (none when it is None or empty: the brain alone) into
    THREED_VIEWER_DIR and returns its path, or None if the CCF download it needs was declined. Raises on a bad file or
    a missing plotly. Blocking and, the first time, slow: run it from the thread that owns the Tk windows."""
    import plotly  # noqa: F401  (fail early, before any of the slow steps)
    loaded = read_coordinates_for_3d(csv_path) if csv_path else None
    have_ccf = loaded['have_ccf'] if loaded else True
    groups = loaded['groups'] if loaded else []
    axis_labels = loaded['axis_labels'] if loaded else (
        'left - right (CCF z, mm)', 'anterior - posterior (CCF x, mm)', 'dorsal - ventral (CCF y, mm)')
    shell_xyz, outline, edges = None, None, None
    if have_ccf:
        if not ensure_ccf_downloaded(abc_cache, parent=parent):
            return None
        _aligned, shell = run_with_loading_window(
            "Reading the CCF coordinates (the first read of this 4-million-cell table takes about a minute).",
            lambda: load_ccf_for_cells(abc_cache, [], THREED_SHELL_POINTS))
        shell_xyz = shell[:, [2, 0, 1]]
        outline = brain_outline_segments(shell_xyz)
        edges = build_3d_edges(abc_cache, loaded['structures'] if loaded else None, parent)
    THREED_VIEWER_DIR.mkdir(parents=True, exist_ok=True)
    stem = f"3d_{Path(csv_path).stem}" if csv_path else "3d_brain"
    path = unique_export_paths(THREED_VIEWER_DIR, stem, ('.html',))['.html']
    title = loaded['title'] if loaded else "3D viewer: brain (no cells loaded)"
    write_3d_group_html(groups, shell_xyz, outline, axis_labels, title, path, edges=edges)
    return path


THREED_LAST_CSV_FILE = CACHE_DIR / '3d_viewer_last_csv.txt'   # the launcher's last cells file, for the next run


THREED_RECENT_CSV_FILE = CACHE_DIR / '3d_viewer_recent_csv.json'   # the launcher's recently opened cells files
THREED_RECENT_CSV_MAX = 12


def recent_3d_csvs():
    """The recently opened coordinates CSVs (newest first) that still exist, as a list of path strings."""
    try:
        listed = json.loads(THREED_RECENT_CSV_FILE.read_text(encoding='utf-8'))
        return [item for item in listed if isinstance(item, str) and Path(item).is_file()]
    except (OSError, ValueError):
        return []


def add_recent_3d_csv(path_text):
    """Put `path_text` first in the recent list (once), keeping at most THREED_RECENT_CSV_MAX; a failure to save is ignored."""
    if not path_text:
        return
    try:
        try:
            listed = json.loads(THREED_RECENT_CSV_FILE.read_text(encoding='utf-8'))
        except (OSError, ValueError):
            listed = []
        key = os.path.normcase(os.path.normpath(path_text))
        listed = [path_text] + [item for item in listed if isinstance(item, str)
                                and os.path.normcase(os.path.normpath(item)) != key]
        THREED_RECENT_CSV_FILE.parent.mkdir(parents=True, exist_ok=True)
        THREED_RECENT_CSV_FILE.write_text(json.dumps(listed[:THREED_RECENT_CSV_MAX], indent=1), encoding='utf-8')
    except OSError:
        pass


def saved_3d_csv():
    """The path last saved for the 3D viewer launcher, whether or not that file can be found now ('' if none)."""
    try:
        return THREED_LAST_CSV_FILE.read_text(encoding='utf-8').strip()
    except OSError:
        return ''


def remembered_3d_csv():
    """The coordinates CSV the 3D viewer launcher last had in its box, if it still exists, else ''."""
    text = saved_3d_csv()
    try:
        return text if text and Path(text).is_file() else ''
    except OSError:
        return ''


def remember_3d_csv(path_text):
    """Remember `path_text` (or forget, if empty) as the launcher's cells file; a failure to save is not worth an error."""
    try:
        THREED_LAST_CSV_FILE.parent.mkdir(parents=True, exist_ok=True)
        THREED_LAST_CSV_FILE.write_text(path_text or '', encoding='utf-8')
    except OSError:
        pass


def open_3d_viewer_launcher(abc_cache):
    """The window behind the opening screen's "Open 3D Viewer…" button: a box at the top for a coordinates_*.csv (drag
    it in from Windows Explorer, paste its path, or Browse): dropping or browsing a file writes the interactive page
    and opens it in the browser at once, and "Open brain only / reopen" does the same for whatever the box holds (the
    brain alone when it is empty). The window closes itself once the page has opened (an error leaves it open, with the
    message). Returns when the window is closed."""
    import tkinter as tk
    from tkinter import filedialog, ttk
    dnd_files = None
    try:
        from tkinterdnd2 import DND_FILES, TkinterDnD
        root = TkinterDnD.Tk()
        dnd_files = DND_FILES
    except Exception:
        root = tk.Tk()   # no tkinterdnd2: choose the file with Browse or paste its path
    # (Every variable below is created with master=root: the opening screen has its own Tk window, and a variable made
    # without a master would live in THAT window's interpreter, so this window's widgets would never see it.)
    root.title("3D viewer")
    root.geometry("+120+120")
    body = ttk.Frame(root, padding=12)
    body.pack(fill='both', expand=True)
    ttk.Label(body, text="Cells: a coordinates_*.csv from the viewer's Export coordinates (the arrow lists recent files)").pack(
        anchor='w')
    row = ttk.Frame(body)
    row.pack(fill='x', pady=(4, 0))
    path_var = tk.StringVar(master=root, value=remembered_3d_csv())   # last run's file, there to reopen with one click
    # Whatever lands in the box is remembered at once, whether it was dropped, chosen with Browse, pasted or typed (and
    # whether or not the page then opens); only a file that still exists is offered again next time.
    path_var.trace_add('write', lambda *_args: remember_3d_csv(path_var.get().strip().strip('"')))
    # An editable drop-down: type or paste a path, or pick one of the recently opened files from the list.
    entry = ttk.Combobox(row, textvariable=path_var, values=recent_3d_csvs(), width=70)
    entry.pack(side='left', fill='x', expand=True)
    missing_saved = saved_3d_csv() if not path_var.get() else ''
    status_var = tk.StringVar(
        master=root, value=("Last run's file is in the box: press Open to reopen it, or drop another." if path_var.get()
               else (f"The file used last time cannot be found now (is its drive connected?): {missing_saved}"
                     if missing_saved else "No cells loaded: the page will show the brain alone.")))

    def set_status(text):
        status_var.set(text)
        root.update_idletasks()

    busy = {'on': False}
    open_button = {}   # the "Open brain only / reopen" button (made below), greyed out while a page is being prepared

    def open_view():
        if busy['on']:
            return
        wanted = path_var.get().strip().strip('"')
        if wanted and not Path(wanted).is_file():
            set_status(f"Not a file: {wanted}")
            return
        busy['on'] = True
        if 'widget' in open_button:   # greyed out, so it is clear the page is already being made (a drop starts it too)
            open_button['widget'].state(['disabled'])
        set_status("Preparing the 3D page… (the first time can take a minute or two)")
        try:
            page = build_3d_viewer_page(abc_cache, wanted or None, parent=root)
            if page is None:
                set_status("Cancelled: the CCF coordinates are needed for the brain background.")
            else:
                open_with_default_viewer(page)
                remember_3d_csv(wanted)
                add_recent_3d_csv(wanted)
                print(f"Opened {page.name} in your browser (saved in {page.parent}).")
                root.after(200, root.destroy)   # the page is open: this window has done its job
        except ImportError:
            set_status("The 3D viewer needs plotly: run  python -m pip install plotly  and try again.")
        except Exception as e:
            set_status(f"Could not open the 3D view: {e}")
            print(f"Could not open the 3D view: {e}")
        finally:
            busy['on'] = False
            try:
                if 'widget' in open_button:
                    open_button['widget'].state(['!disabled'])
            except tk.TclError:
                pass   # (the window has already closed)

    def browse():
        current = path_var.get().strip().strip('"')
        start = Path(current).parent if current and Path(current).parent.is_dir() else None
        chosen = filedialog.askopenfilename(
            title="Choose a coordinates CSV", parent=root, initialdir=str(start) if start else None,
            filetypes=[("Coordinates CSV", "coordinates*.csv"), ("CSV files", "*.csv"), ("All files", "*.*")])
        if chosen:
            path_var.set(os.path.normpath(chosen))
            open_view()

    def on_drop(event):
        paths = root.tk.splitlist(event.data)
        if paths:
            path_var.set(os.path.normpath(paths[0]))
            set_status("Dropped " + Path(paths[0]).name + (" (the first of several: one file at a time)"
                                                           if len(paths) > 1 else ""))
            root.after(50, open_view)
        return event.action

    ttk.Button(row, text="Browse…", command=browse).pack(side='left', padx=(6, 0))
    ttk.Button(row, text="Clear", command=lambda: (path_var.set(''), remember_3d_csv(''), set_status(
        "No cells loaded: the page will show the brain alone."))).pack(side='left', padx=(4, 0))
    entry.bind('<Return>', lambda e: open_view())
    if dnd_files is not None:
        entry.drop_target_register(dnd_files)
        entry.dnd_bind('<<Drop>>', on_drop)
        root.drop_target_register(dnd_files)
        root.dnd_bind('<<Drop>>', on_drop)
        hint = "Drag a coordinates_*.csv from Explorer into the box above (or anywhere on this window)."
    else:
        hint = "Choose the file with Browse… or paste its path (drag and drop needs the tkinterdnd2 package)."
    ttk.Label(body, text=hint, foreground='#555555').pack(anchor='w', pady=(6, 0))
    ttk.Label(body, textvariable=status_var, foreground='#7a3b00', wraplength=640, justify='left').pack(
        anchor='w', pady=(10, 0))
    buttons = ttk.Frame(body)
    buttons.pack(fill='x', pady=(12, 0))
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side='right')
    open_button['widget'] = ttk.Button(buttons, text="Open brain only / reopen", command=open_view)
    open_button['widget'].pack(side='right', padx=(0, 6))
    entry.focus_set()
    root.lift()
    try:
        root.mainloop()
    except Exception as e:
        print(f"The 3D viewer window stopped: {e}")
    finally:
        try:
            root.destroy()
        except Exception:
            pass


def prompt_startup_panel(out_folder):
    """First screen shown at launch (abc_tk_launcher.prompt_startup_panel_tk): the output folder, then either open the section/ROI picker
    (the full pipeline, Step 0 onward), load a previously computed run's cached files directly into the interactive UMAP viewer, or open
    the 3D viewer. (The matplotlib version of this screen was removed: the cell type is no longer asked here at all.)

    Returns (action, cell_type_selection, load_folder, chosen_out_folder): action is 'new_run', 'load_existing', 'own_dataset' or
    'all_neurons' (make the "All neurons" project: see abc_tk_launcher.ALL_NEURONS_FOLDER); cell_type_selection is
    None for a new run (asked after the sections and ROIs are confirmed) and unused for 'load_existing'; load_folder is the chosen run
    folder (a Path) when action is 'load_existing', else None; chosen_out_folder is the output folder the user settled on (a Path, already
    created) — the caller should adopt it as the working directory for this session.

    Exits the process if the window is closed without choosing an action. Raises if the Tk launch screen can't be loaded (the caller
    falls back to a console prompt)."""
    start_section_coordinate_load(abc_cache)   # (the table of sections' z values, read in the background for the windows to come)
    compute_window.close()                     # (a run that stopped partway, with an error, may have left its "Computing" window up)
    from abc_tk_launcher import prompt_startup_panel_tk
    outcome = prompt_startup_panel_tk(out_folder, abc_cache, {
        'load_recent_output_folders': load_recent_output_folders,
        'remember_recent_output_folder': remember_recent_output_folder,
        'required_cached_run_files': required_cached_run_files,
        'open_3d_viewer': open_3d_viewer_launcher,
        'settings_dialog': show_coordinate_settings_dialog,
        'atlas_download_running': atlas_download_running,
    })
    if outcome is None:
        print("No selection made; exiting.")
        sys.exit(0)
    return outcome


def prompt_section_selection_gui(adata, abc_cache, section_series, out_folder=None, cell_type_selection='All',
                                 imputed_state=None):
    """Clickable, zoomable/pannable grid of per-section spatial thumbnails
    (a single cached bitmap, not one Axes per section — see
    generate_section_grid_image). Click a thumbnail to toggle it
    (multi-select), scroll to zoom/pan the whole grid, then press 'Confirm
    Selection'. Double-click a thumbnail to open a single-section ROI picker
    for that section (after a confirmation, since it takes a few seconds to
    load) — this lets you collect ROIs across multiple sections in one
    session. 'Load ROIs' opens a file picker (defaulting to `out_folder`,
    where ROI CSVs get saved) to merge in a previously saved set. Returns
    (selected_labels, sections_suffix, rois, imputed_adata): the first two
    have the same contract as prompt_section_selection(); rois is a list of
    dicts ({'section','x_min','x_max','y_min','y_max'}), empty if none were
    drawn; imputed_adata is the imputed-gene-expression AnnData if 'Imputed
    Gene' was ever selected this session (see ensure_imputed_gene_dataset_
    loaded) — already loaded and ready to reuse, e.g. by a later interactive
    UMAP viewer — or None if it never was.
    Raises RuntimeError if a GUI can't be shown (caller should fall back to
    the console prompt in that case).

    `cell_type_selection` is the same 'Neurons'/'NonNeurons'/'All' choice
    prompt_cell_type_selection() collects before this is even called — note
    that `adata` here is called with, at this point in the pipeline, still
    unfiltered by it (filter_by_cell_type() doesn't run until after this
    whole picker session ends), so every cell type is still actually
    present in every section shown here regardless of this choice; it's
    only threaded through to each single-section window's gene view, where
    it controls which type (if any) gets greyed out rather than colored by
    expression — see render_gene_expression_array's own docstring."""
    # the Tk picker and single-section window (abc_tk_grid.py, abc_tk_picker.py)
    from abc_tk_grid import prompt_section_selection_tk
    if imputed_state is None:
        imputed_state = {'adata': None, 'load_thread': None, 'load_error': None}
    outcome = prompt_section_selection_tk(adata, abc_cache, section_series, out_folder, cell_type_selection, imputed_state,
                                          ensure_imputed=ensure_imputed_gene_dataset_loaded,
                                          ccf_host={'ensure_ccf_downloaded': ensure_ccf_downloaded, 'load_ccf_for_cells': load_ccf_for_cells})
    if outcome[0] == 'cancelled':
        raise UserCancelledSelection("Section picker cancelled.")
    if outcome[0] != 'confirmed':
        raise RuntimeError("section picker window was closed without confirming a selection")
    _status, chosen, chosen_rois, n_sections = outcome
    # (Nothing picked means the whole brain: the options window that follows says so, instead of a question here.)
    if not chosen or len(chosen) == n_sections:
        return None, 'allsections', chosen_rois, imputed_state['adata']
    tokens = '-'.join(sanitize_section_token(s) for s in chosen)
    return chosen, f'sections-{tokens}', chosen_rois, imputed_state['adata']


def filter_by_sections(adata, section_series, selected_sections, section_col=SECTION_COL):
    """Subset adata to cells whose section label is in `selected_sections`
    (a list of labels), also stashing the section label onto adata.obs for
    downstream use. `selected_sections=None` means "all sections" (no-op)."""
    if section_series is not None and section_col not in adata.obs.columns:
        adata.obs[section_col] = section_series.reindex(adata.obs.index)
    if selected_sections is None:
        return adata
    if section_col not in adata.obs.columns:
        print(f"Warning: '{section_col}' column not available; cannot filter by section. "
              "Proceeding with all sections.")
        return adata
    keep = adata.obs[section_col].isin(selected_sections)
    filtered = materialize_subset(adata, keep.to_numpy())
    print(f"Filtered to {len(selected_sections)} section(s): {filtered.n_obs} of {adata.n_obs} cells kept.")
    return filtered








# ==============================================================================
# Gene entry, shared by the single-section ROI picker and the interactive UMAP
# viewer. Both have a text box taking comma-separated gene names (1-3 here,
# via MAX_GENE_NAMES below), with autocomplete for the name being typed;
# everything about *entering* genes lives here so the two behave identically.
# Rendering stays per-window — the interactive UMAP viewer's own redraw_gene
# uses a higher, viewer-local cap (INTERACTIVE_MULTI_GENE_MAX_NAMES, up to 6)
# instead of MAX_GENE_NAMES, since it's the only renderer that does anything
# with genes 4-6 (drawn as a '+' overlay — see redraw_multi_genes).
# ==============================================================================

MAX_GENE_NAMES = 3
TOO_MANY_GENES_MESSAGE = f"Enter at most {MAX_GENE_NAMES} gene names, separated by commas."
# Color of a gene name that didn't resolve, shown in the box itself until the
# next edit.
GENE_NAME_ERROR_COLOR = 'red'
# Faded text/outline color for controls disabled while a window is busy. Only
# text and outlines fade, not backgrounds: a paler fill reads as highlighted
# rather than disabled.
DISABLED_CONTROL_COLOR = '0.65'






def _gene_name_segments(text, bad_names):
    """Split `text` into drawable (segment, is_bad) pieces for per-name
    coloring. A segment never ends in whitespace: text extents are measured
    from ink, so a trailing space has no width and chaining the next piece
    from that edge would drop it. Whitespace is carried onto the start of
    the following piece instead, where it measures correctly."""
    bad = {b.casefold() for b in bad_names}
    pieces = []
    for piece in re.split(r'(,)', text):
        if not piece:
            continue
        if piece == ',':
            pieces.append((',', False))
            continue
        body = piece.rstrip()
        if body:
            pieces.append((body, body.strip().casefold() in bad))
        trailing = piece[len(body):]
        if trailing:
            pieces.append((trailing, False))
    segments, carry = [], ''
    for seg, is_bad in pieces:
        if not seg.strip():
            carry += seg
            continue
        segments.append((carry + seg, is_bad))
        carry = ''
    return segments


def mark_invalid_gene_names(textbox, bad_names):
    """Show the names in `bad_names` in GENE_NAME_ERROR_COLOR inside the box,
    leaving the rest in the box's normal text color. Clear with
    clear_gene_name_marks().

    A TextBox draws its text as one Text artist, in one color, so this makes
    that artist transparent (it still positions the typing cursor) and draws
    the text as a chain of pieces, each anchored at the right edge of the
    one before and colored per name. Measured against the normal single-color
    rendering, every glyph lands within 1 px at 200 dpi. The pieces are
    children of the box's axes, so they're drawn wherever the box is, blitted
    or not. Call again after changing the box's font size to rebuild them."""
    _remove_gene_name_mark_artists(textbox)
    bad_names = list(bad_names)
    textbox._gene_name_marks = {'bad': bad_names, 'artists': [], 'override_color': None}
    if not bad_names:
        textbox.text_disp.set_alpha(None)
        return
    text_disp = textbox.text_disp
    base_color = text_disp.get_color()
    text_disp.set_alpha(0.0)
    previous = None
    for segment, is_bad in _gene_name_segments(text_disp.get_text(), bad_names):
        style = dict(fontsize=text_disp.get_fontsize(), fontfamily=text_disp.get_fontfamily(),
                     color=GENE_NAME_ERROR_COLOR if is_bad else base_color, annotation_clip=False)
        if previous is None:
            artist = textbox.ax.annotate(segment, xy=text_disp.get_position(), xycoords=text_disp.get_transform(),
                                         va=text_disp.get_va(), ha='left', **style)
        else:
            artist = textbox.ax.annotate(segment, xy=(1, 0), xycoords=previous, va='bottom', ha='left', **style)
        artist._gene_name_is_bad = is_bad
        textbox._gene_name_marks['artists'].append(artist)
        previous = artist


def _remove_gene_name_mark_artists(textbox):
    marks = getattr(textbox, '_gene_name_marks', None)
    if marks:
        for artist in marks['artists']:
            artist.remove()
        marks['artists'] = []


def clear_gene_name_marks(textbox):
    """Undo mark_invalid_gene_names(): remove the colored pieces and show the
    box's own text again. Safe to call when nothing is marked."""
    _remove_gene_name_mark_artists(textbox)
    textbox._gene_name_marks = None
    textbox.text_disp.set_alpha(None)










def region_naming(sections_suffix, rois):
    """The sections and ROIs parts of a run folder's name: (sections_suffix, region_suffix, roi_sections). With ROIs they name the sections
    the ROIs are on and carry a short hash of the exact ROI bounds, which keeps the name unique per distinct ROI configuration without making
    it unreasonably long; without, the sections_suffix is the picker's (e.g. 'allsections') and the region_suffix is empty."""
    if not rois:
        return sections_suffix, '', None
    roi_sections = sorted_sections_descending(list({roi['section'] for roi in rois}))
    roi_key = json.dumps(
        sorted(
            [{k: (round(v, 3) if isinstance(v, float) else v) for k, v in roi.items()} for roi in rois],
            key=lambda r: (r['section'], r['x_min'], r['y_min']),
        )
    )
    if reconstructed_in_use():
        roi_key = 'reconstructed|' + roi_key                    # (the same bounds, drawn on Allen's reconstructed coordinates, select other cells)
    roi_hash = hashlib.md5(roi_key.encode()).hexdigest()[:8]
    return 'sections-' + '-'.join(sanitize_section_token(s) for s in roi_sections), f'_rois{len(rois)}-{roi_hash}', roi_sections


def count_selected_cells(adata, abc_cache, section_series, selected_sections, rois):
    """How many cells a selection holds, without filtering anything: {'Neurons': n, 'NonNeurons': n, 'All': n} for the sections / ROIs chosen
    (as filter_by_sections / filter_by_rois would select them: ROIs, if any, replace the whole-section choice), from the cached annotation and
    coordinate tables, so it takes a few seconds."""
    annotations = load_annotations_aligned(adata.obs.index, abc_cache)
    sections = section_series.to_numpy() if section_series is not None else annotations[SECTION_COL].to_numpy()
    if rois:
        spatial = load_section_spatial_coords(adata, abc_cache)
        xs = pd.to_numeric(spatial['x'], errors='coerce').to_numpy()
        ys = pd.to_numeric(spatial['y'], errors='coerce').to_numpy()
        keep = np.zeros(len(xs), dtype=bool)
        for roi in rois:
            keep |= (sections == roi['section']) & (xs >= roi['x_min']) & (xs <= roi['x_max']) & (ys >= roi['y_min']) & (ys <= roi['y_max'])
    elif selected_sections is not None:
        keep = pd.Series(sections).isin(selected_sections).to_numpy()
    else:
        keep = np.ones(len(sections), dtype=bool)
    # The class number is read from each distinct class name once (about 34 of them), and only for the cells in the selection: matching the
    # pattern against every one of the 4 million cells took three seconds.
    codes, names = pd.factorize(annotations['class'].to_numpy()[keep])
    name_ids = pd.to_numeric(pd.Series(names).astype(str).str.strip().str.extract(r'^(\d+)')[0], errors='coerce').to_numpy(dtype=float)
    class_ids = np.where(codes >= 0, name_ids[np.maximum(codes, 0)], np.nan) if len(name_ids) else np.full(len(codes), np.nan)
    valid = ~np.isnan(class_ids)
    is_non_neuron = valid & np.isin(class_ids, list(NON_NEURON_CLASS_IDS))
    return {'All': int(keep.sum()), 'Neurons': int((valid & ~is_non_neuron).sum()), 'NonNeurons': int(is_non_neuron.sum())}


def filter_by_rois(adata, abc_cache, rois, section_col=SECTION_COL):
    """Subset adata to cells falling within any of `rois` (each a dict with
    'section', 'x_min', 'x_max', 'y_min', 'y_max'). A cell must match both
    the ROI's section and its spatial bounds. No-op if rois is empty/None.
    Assumes adata.obs[section_col] is already populated (e.g. by a prior
    filter_by_sections() call)."""
    if not rois:
        return adata
    if section_col not in adata.obs.columns:
        print(f"Warning: '{section_col}' column not available; cannot filter by ROIs.")
        return adata
    spatial = load_section_spatial_coords(adata, abc_cache)
    if spatial is None or 'x' not in spatial.columns or 'y' not in spatial.columns:
        print("Warning: spatial x/y coordinates not available; cannot filter by ROIs.")
        return adata

    xs = pd.to_numeric(spatial['x'], errors='coerce').to_numpy()
    ys = pd.to_numeric(spatial['y'], errors='coerce').to_numpy()
    sections = adata.obs[section_col].to_numpy()

    keep = np.zeros(len(xs), dtype=bool)
    for roi in rois:
        section_mask = sections == roi['section']
        bounds_mask = (
            (xs >= roi['x_min']) & (xs <= roi['x_max'])
            & (ys >= roi['y_min']) & (ys <= roi['y_max'])
        )
        keep |= section_mask & bounds_mask

    filtered = materialize_subset(adata, keep)
    n_sections = len({roi['section'] for roi in rois})
    print(f"Filtered to {len(rois)} ROI(s) across {n_sections} section(s): "
          f"{filtered.n_obs} of {adata.n_obs} cells kept.")
    return filtered


def save_roi_map(adata, abc_cache, section_series, rois, output_path):
    """Save a grid of per-section scatter plots — one panel per section that
    has at least one ROI, each showing all of that section's cells (gray
    non-neurons, colored neurons) plus its ROI rectangle(s) outlined in red —
    as a visual reference for where the ROIs were drawn."""
    spatial = load_section_spatial_coords(adata, abc_cache)
    if spatial is None or 'x' not in spatial.columns or 'y' not in spatial.columns:
        print("Warning: spatial coordinates not available; cannot save ROI map.")
        return

    print(f"Generating ROI map ...")
    rois_by_section = {}
    for roi in rois:
        rois_by_section.setdefault(roi['section'], []).append(roi)
    sections = sorted_sections_descending(list(rois_by_section.keys()))
    n = len(sections)

    ncols = max(1, math.ceil(math.sqrt(n)))
    nrows = math.ceil(n / ncols)
    # (a plain Agg figure, not pyplot's: this runs on the "Computing" window's worker thread, where a pyplot figure made a Tk window and the program
    # died with "main thread is not in main loop" while saving the map)
    fig = Figure(figsize=(ncols * 4, nrows * 4))
    FigureCanvasAgg(fig)
    axes = fig.subplots(nrows, ncols, squeeze=False)
    axes_flat = axes.ravel()

    section_values = section_series.to_numpy()
    xs_all = pd.to_numeric(spatial['x'], errors='coerce').to_numpy()
    ys_all = pd.to_numeric(spatial['y'], errors='coerce').to_numpy()
    class_ids_all = None
    if 'class' in spatial.columns:
        class_ids_all = extract_leading_numeric_id(spatial['class']).to_numpy()
    point_colors_all, neuron_mask_all = build_neuron_class_colors(class_ids_all)

    for i, section in enumerate(sections):
        ax = axes_flat[i]
        mask = section_values == section
        xs, ys = xs_all[mask], ys_all[mask]
        valid = ~(np.isnan(xs) | np.isnan(ys))
        xs, ys = xs[valid], ys[valid]
        if neuron_mask_all is not None:
            is_neuron = neuron_mask_all[mask][valid]
            colors = point_colors_all[mask][valid]
            scatter_gray_then_colored(ax, xs, ys, is_neuron, colors, gray_size=1, colored_size=2)
        else:
            ax.scatter(xs, ys, s=1, c='steelblue', linewidths=0)

        section_rois = rois_by_section[section]
        for roi in section_rois:
            rect = Rectangle(
                (roi['x_min'], roi['y_min']),
                roi['x_max'] - roi['x_min'], roi['y_max'] - roi['y_min'],
                # Explicit zorder so the outline always draws over the
                # scatter dots (their default zorder ties with a patch's,
                # which would otherwise leave the ordering to insertion
                # order and risk the line looking broken where dots overlap it).
                linewidth=2, edgecolor='red', facecolor='none', zorder=5,
            )
            ax.add_patch(rect)

        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        roi_word = 'ROI' if len(section_rois) == 1 else 'ROIs'
        ax.set_title(f"{sanitize_section_token(section)} ({len(section_rois)} {roi_word})", fontsize=9)

    for j in range(n, len(axes_flat)):
        axes_flat[j].axis('off')

    print(f"Saving ROI map to {output_path} ...", end='', flush=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    fig.savefig(output_path.with_suffix('.svg'))
    print(f" Done saving")


def save_group_spatial_maps(adata, abc_cache, section_series, rois, whole_sections, output_dir,
                             group_col='subclass', cell_type_selection='All'):
    """For each section with at least one ROI or a whole-section pick,
    find which `group_col` categories (e.g. 'subclass' or 'supertype')
    actually occur there among cells of `cell_type_selection` ('Neurons',
    'NonNeurons', or 'All') — within the ROI bounds specifically for a
    section with ROI(s) (a category elsewhere in that section but outside
    every ROI doesn't count), or anywhere in the section for a
    whole-section pick (no ROI to restrict to). Take the union of
    qualifying categories across all such sections, split it into groups of
    up to 10, and save one spatial-scatter grid figure per group into
    `output_dir` (one panel per relevant section, same grid layout as
    save_roi_map). Every cell not in the current group's *qualifying*
    categories — including any cell of the excluded cell type — is drawn
    gray so the section's outline stays visible; each section's ROI
    rectangle(s), if any, are outlined in red on top. Each group's up to 10
    categories get one of SPATIAL_MAP_COLORS each, assigned by position
    within the group, so the same 10 colors are reused (in the same order)
    across every group's figure rather than each group inventing its own
    palette. No-op (prints a message) if there are no relevant sections or
    no qualifying categories."""
    roi_sections = sorted_sections_descending(list({roi['section'] for roi in rois})) if rois else []
    relevant_sections = sorted_sections_descending(list(set(roi_sections) | set(whole_sections or [])))
    if not relevant_sections:
        relevant_sections = sorted_sections_descending(list(section_series.dropna().unique()))
        print(f"No ROI/whole-section picks; defaulting to all {len(relevant_sections)} section(s) for {group_col} spatial maps.")

    spatial = load_section_spatial_coords(adata, abc_cache)
    if spatial is None or 'x' not in spatial.columns or 'y' not in spatial.columns or group_col not in spatial.columns:
        print(f"Warning: spatial x/y/{group_col} columns not available; cannot save {group_col} spatial maps.")
        return

    rois_by_section = {}
    for roi in rois:
        rois_by_section.setdefault(roi['section'], []).append(roi)

    section_values = section_series.reindex(spatial.index).to_numpy()
    xs_all = pd.to_numeric(spatial['x'], errors='coerce').to_numpy()
    ys_all = pd.to_numeric(spatial['y'], errors='coerce').to_numpy()
    group_ids_all = extract_leading_numeric_id(spatial[group_col]).to_numpy()
    group_labels_all = spatial[group_col].astype(str).to_numpy()

    # Restrict everything below to the selected cell type, same rule as
    # filter_by_cell_type(): non-neuron classes are 30-34, everything else
    # valid counts as neuron. Cells excluded here are treated exactly like
    # cells with no `group_col` category at all — always gray, never
    # counted toward a section's qualifying categories.
    if cell_type_selection != 'All' and 'class' in spatial.columns:
        class_ids_all = extract_leading_numeric_id(spatial['class']).to_numpy()
        valid_class = ~np.isnan(class_ids_all)
        is_non_neuron_all = valid_class & np.isin(class_ids_all, list(NON_NEURON_CLASS_IDS))
        if cell_type_selection == 'NonNeurons':
            cell_type_mask_all = valid_class & is_non_neuron_all
        else:  # 'Neurons'
            cell_type_mask_all = valid_class & ~is_non_neuron_all
    else:
        cell_type_mask_all = np.ones(len(xs_all), dtype=bool)

    # For each relevant section, the set of category IDs that qualify there
    # (within-ROI only, if the section has ROIs) — plus a global id->label
    # lookup and the union of all qualifying IDs across every section.
    section_qualifying = {}
    qualifying_ids = set()
    id_to_label = {}
    for section in relevant_sections:
        sec_mask = section_values == section
        if section in rois_by_section:
            target_mask = np.zeros(len(xs_all), dtype=bool)
            for roi in rois_by_section[section]:
                target_mask |= (
                    sec_mask
                    & (xs_all >= roi['x_min']) & (xs_all <= roi['x_max'])
                    & (ys_all >= roi['y_min']) & (ys_all <= roi['y_max'])
                )
        else:
            target_mask = sec_mask
        target_mask = target_mask & cell_type_mask_all
        ids_here = group_ids_all[target_mask]
        valid_ids = set(ids_here[~np.isnan(ids_here)].tolist())
        section_qualifying[section] = valid_ids
        qualifying_ids.update(valid_ids)
        for cid, label in zip(group_ids_all[sec_mask & cell_type_mask_all], group_labels_all[sec_mask & cell_type_mask_all]):
            if not np.isnan(cid) and cid not in id_to_label:
                id_to_label[cid] = label

    if not qualifying_ids:
        print(f"No {group_col} categories found within the selected ROIs/whole sections; skipping spatial maps.")
        return

    sorted_ids = sorted(qualifying_ids)
    groups = [sorted_ids[i:i + 10] for i in range(0, len(sorted_ids), 10)]

    ncols = max(1, math.ceil(math.sqrt(len(relevant_sections))))
    nrows = math.ceil(len(relevant_sections) / ncols)

    # A shared axis span (data units) — separate for x and y — and a
    # per-section centroid, so every panel is drawn at the same physical
    # scale instead of each autoscaling to fill its cell with just its own
    # data — otherwise an anatomically small section (e.g. the olfactory
    # bulb) gets zoomed to look just as large on the page as a full coronal
    # section. Width and height are sized independently (rather than one
    # isotropic span reused for both) because sections are typically wider
    # than tall — a single shared span would end up sized to the width (the
    # larger of the two for most sections), leaving a lot of unnecessary
    # blank margin above and below every panel even though its horizontal
    # fit was already tight. Each is sized to the SPAN_PERCENTILE-th
    # percentile of section extents (not the strict max) plus 5% padding —
    # using the single largest section would leave most panels surrounded
    # by a lot of blank margin; a percentile trades that off deliberately,
    # clipping the few sections above it in this figure only (the full data
    # is still used for actual processing) for noticeably less wasted space
    # in the common case. Every panel still shares the exact same (width,
    # height) pair, so set_aspect('equal') still shrinks every panel by the
    # same amount within its grid cell — which is what keeps row gaps even.
    SPAN_PERCENTILE = 100  # temporarily testing without percentile clipping — was 90
    section_centroids = {}
    half_widths, half_heights = [], []
    for section in relevant_sections:
        sec_mask = section_values == section
        xs_sec, ys_sec = xs_all[sec_mask], ys_all[sec_mask]
        valid_sec = ~(np.isnan(xs_sec) | np.isnan(ys_sec))
        xs_sec, ys_sec = xs_sec[valid_sec], ys_sec[valid_sec]
        if len(xs_sec) == 0:
            continue
        x_min, x_max = xs_sec.min(), xs_sec.max()
        y_min, y_max = ys_sec.min(), ys_sec.max()
        section_centroids[section] = ((x_min + x_max) / 2, (y_min + y_max) / 2)
        half_widths.append((x_max - x_min) / 2)
        half_heights.append((y_max - y_min) / 2)
    shared_half_width = np.percentile(half_widths, SPAN_PERCENTILE) * 1.05 if half_widths else 1.0
    shared_half_height = np.percentile(half_heights, SPAN_PERCENTILE) * 1.05 if half_heights else 1.0

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Clear out any group*.png/.svg files from a previous run in this folder
    # first — otherwise a run that produces fewer groups than a prior one
    # (e.g. fewer qualifying subclasses this time) would leave a stale
    # leftover file that looks like part of the current output.
    for pattern in ('group*.png', 'group*.svg'):
        for stale in output_dir.glob(pattern):
            stale.unlink()

    for group_idx, group_ids in enumerate(groups, start=1):
        color_by_cid = {cid: SPATIAL_MAP_COLORS[i % len(SPATIAL_MAP_COLORS)] for i, cid in enumerate(group_ids)}

        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 4), squeeze=False)
        axes_flat = axes.ravel()

        for i, section in enumerate(relevant_sections):
            ax = axes_flat[i]
            sec_mask = section_values == section
            xs, ys = xs_all[sec_mask], ys_all[sec_mask]
            ids_here = group_ids_all[sec_mask]
            in_cell_type = cell_type_mask_all[sec_mask]
            valid = ~(np.isnan(xs) | np.isnan(ys))
            xs, ys, ids_here, in_cell_type = xs[valid], ys[valid], ids_here[valid], in_cell_type[valid]

            # Only this section's qualifying categories (within its own ROI,
            # or anywhere if it's a whole-section pick), and only cells of
            # the selected cell type, get colored here — everything else
            # (a group category that qualified via a *different* section,
            # the wrong cell type, or no category at all) falls back to gray.
            active_group = [cid for cid in group_ids if cid in section_qualifying.get(section, set())]
            colored_mask = in_cell_type & np.isin(ids_here, active_group)

            ax.scatter(xs[~colored_mask], ys[~colored_mask], s=1, color=DEFAULT_GREY_RGBA, linewidths=0)
            for cid in active_group:
                cid_mask = colored_mask & (ids_here == cid)
                if not cid_mask.any():
                    continue
                ax.scatter(xs[cid_mask], ys[cid_mask], s=2, color=color_by_cid[cid], linewidths=0, alpha=0.5)

            for roi in rois_by_section.get(section, []):
                rect = Rectangle(
                    (roi['x_min'], roi['y_min']),
                    roi['x_max'] - roi['x_min'], roi['y_max'] - roi['y_min'],
                    linewidth=2, edgecolor='red', facecolor='none', zorder=5,
                )
                ax.add_patch(rect)

            cx, cy = section_centroids.get(section, (0.0, 0.0))
            ax.set_xlim(cx - shared_half_width, cx + shared_half_width)
            ax.set_ylim(cy - shared_half_height, cy + shared_half_height)
            ax.set_aspect('equal')
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(section, fontsize=9)

        for j in range(len(relevant_sections), len(axes_flat)):
            axes_flat[j].axis('off')

        # One shared legend for the whole figure — every panel uses the same
        # cid -> color mapping, so a per-panel legend would just repeat itself.
        legend_handles = [
            Line2D([0], [0], marker='o', linestyle='none', markersize=6,
                   color=color_by_cid[cid], label=id_to_label.get(cid, str(cid)))
            for cid in group_ids
        ]
        fig.legend(handles=legend_handles, loc='lower center', ncol=min(5, len(legend_handles)),
                   fontsize=7, bbox_to_anchor=(0.5, -0.02))

        group_path = output_dir / f'group{group_idx:02d}.png'
        fig.tight_layout(rect=(0, 0.06, 1, 1))  # leave room at the bottom for the legend
        fig.savefig(group_path, dpi=200, bbox_inches='tight')
        fig.savefig(group_path.with_suffix('.svg'), bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {group_col} spatial map group {group_idx}/{len(groups)} "
              f"({len(group_ids)} {group_col} categories) to {group_path}.")


def format_duration(seconds):
    """Format a duration in seconds as 'Hh MMm SS.Ss', dropping leading
    all-zero units (e.g. 45.3s stays as '45.3s', 125s becomes '2m 5.0s')."""
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours >= 1:
        return f"{int(hours)}h {int(minutes)}m {secs:.1f}s"
    if minutes >= 1:
        return f"{int(minutes)}m {secs:.1f}s"
    return f"{secs:.1f}s"


def resolve_out_folder(default_folder):
    """Ensure `default_folder` exists (creating it if needed) and return it
    for use as the output directory. If it can't be created or accessed —
    e.g. a network drive that isn't currently mounted on this machine —
    falls back to a native folder-picker dialog so the user can choose an
    alternate location instead of the script just crashing with an OSError.
    Exits the process if the user cancels that dialog."""
    try:
        Path(default_folder).mkdir(parents=True, exist_ok=True)
        return default_folder
    except OSError as e:
        print(f"Could not create/access the default output folder '{default_folder}' ({e}).")

    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        chosen = filedialog.askdirectory(
            title="Default output folder unavailable — choose an alternate output folder",
        )
        root.destroy()
    except Exception as e:
        print(f"Could not show folder picker ({e}).")
        chosen = None

    if not chosen:
        print("No output folder chosen; exiting.")
        sys.exit(0)
    Path(chosen).mkdir(parents=True, exist_ok=True)
    return chosen


def prompt_yes_no(question, default=None):
    """Loop until the user answers y/n to `question` (which should include
    its own trailing '(y/n): ' cue). Returns True for yes, False for no.
    If `default` is given (True/False) and the user just hits enter with no
    input, returns `default` instead of re-prompting."""
    while True:
        choice = input(question).strip().lower()
        if choice == '' and default is not None:
            return default
        if choice in ('y', 'yes'):
            return True
        if choice in ('n', 'no'):
            return False
        print("Please enter 'y' or 'n'.")


def confirm_whole_brain():
    """Ask for confirmation before proceeding with no sections or ROIs
    selected, which means processing every cell in the whole brain — a
    much slower run than a typical subset. Returns True if the caller
    should proceed."""
    return prompt_yes_no(
        "No sections or ROIs selected: this will process the entire brain, "
        "which is time-consuming. Continue? (y/n): "
    )


def prompt_subsample_choice(n_obs, default_cap):
    """Console prompt shown when the filtered cell count `n_obs` exceeds
    `default_cap`: let the user choose to use the full filtered dataset (no
    subsampling), restrict to `default_cap`, or subsample to a custom
    value. Returns the chosen subsample target (an int less than n_obs), or
    None to mean 'use every filtered cell, don't subsample'."""
    print(f"Filtered dataset has {n_obs} cells, which exceeds the default "
          f"subsample cap of {default_cap}.")
    while True:
        choice = input(
            f"Use the (f)ull dataset, (r)estrict to {default_cap} (default), "
            "or enter a (c)ustom number to subsample to? [f/r/c]: "
        ).strip().lower()
        if choice in ('f', 'full'):
            return None
        if choice in ('r', 'restrict', ''):
            return default_cap
        if choice in ('c', 'custom'):
            while True:
                raw = input(f"Enter number of cells to subsample to (1-{n_obs - 1}): ").strip()
                try:
                    n = int(raw)
                except ValueError:
                    print("Please enter a valid integer.")
                    continue
                if 1 <= n < n_obs:
                    return n
                print(f"Enter a value between 1 and {n_obs - 1}, "
                      "or choose 'f' to use the full dataset instead.")
        print("Please enter 'f', 'r', or 'c'.")


def prompt_roi_run_label():
    """Optional short, freeform label the user can add to an ROI run's
    folder name — purely a memory aid (e.g. 'hippocampus') alongside the
    auto-generated ROI hash, which is unique but not remotely readable.
    Blank input (just pressing Enter) skips it entirely, leaving the folder
    named exactly as it always was.

    Appended to region_suffix *before* run_suffix/session_prefix/run_folder
    are built from it, so the label becomes an ordinary part of this run's
    identity, the same way subsample_target already is — no changes needed
    anywhere else in the pipeline (early_cache_hit's glob, the later exact
    csv_path check, 'Load existing run', etc. all just see a longer folder
    name and work unmodified). The one consequence, also already true of
    subsample_target: re-running the identical ROI selection later only
    hits this same cached folder if the *same* label (or none) is given
    again — a different or missing label is a different folder name, hence
    a fresh run. Reusing a labeled run on purpose is what the startup
    panel's 'Load existing run' picker is for."""
    return sanitize_run_label(input("Optional label for this ROI run's folder name (blank to skip): "))


def parse_unfinished_run(folder):
    """What a run folder's name says about how it was made, for finishing a run that stopped before its UMAP was made: {'cell_type', 'label', 'cap'}
    (the cells' kind, the project name or '', and the limit of cells or None), or None if the name is not one a run makes."""
    name = Path(folder).name
    m = re.match(r'^umap_(neurons|nonneurons|all-cells)_(.*)$', name)
    if m is None:
        return None
    cell_type = {v: k for k, v in CELL_TYPE_SUFFIXES.items()}[m.group(1)]
    rest = m.group(2)
    cap = None
    capped = re.search(r'_sub(\d+)$', rest)
    if capped:
        cap = int(capped.group(1))
        rest = rest[:capped.start()]
    # (the part after the sections, and after the ROIs' count and hash if there are ROIs: the project name)
    label = re.match(r'^(?:allsections|sections-[0-9-]+)(?:_rois\d+-[0-9a-f]{8})?(?:_(.+))?$', rest)
    if label is None:
        return None
    return {'cell_type': cell_type, 'label': label.group(1) or '', 'cap': cap}


def sanitize_run_label(raw):
    """A label typed for an ROI run's folder name, made filename-safe (see prompt_roi_run_label); '' for none."""
    raw = (raw or '').strip()
    if not raw:
        return ''
    # Alphanumerics kept as-is; every run of anything else (spaces,
    # punctuation) collapsed to a single hyphen, and any leading/trailing
    # hyphen trimmed — readable word separators instead of sanitize_
    # section_token's own all-punctuation-stripped style, which is fine for
    # a bare trailing section number but would squash a multi-word label
    # like "left hippocampus" into one illegible run of letters.
    label = re.sub(r'[^A-Za-z0-9]+', '-', raw).strip('-')
    return label


def confirm_overwrite(path):
    """If `path` already exists, ask the user whether to overwrite it.
    Returns True if the caller should proceed with writing the file."""
    path = Path(path)
    if not path.exists():
        return True
    return prompt_yes_no(f"File '{path}' already exists. Overwrite? (y/n): ")


PIPELINE_STEPS = 9          # a run's steps, as log_status numbers them: 1 the sections and ROIs, 2-8 the computing (see StepWindow), 9 the viewer


def log_status(message):
    """print() prefixed with the elapsed time, in seconds (3 decimals, i.e.
    millisecond resolution), since `selection_confirmed_at` (a
    time.perf_counter() value set once the user confirms their section/ROI
    selection in the main script below) — so the sequence of pipeline steps
    that follow can be timed against each other at a glance. Only meant for
    status updates printed after that point; earlier messages (data
    loading, the picker itself) use plain print()."""
    elapsed_s = time.perf_counter() - selection_confirmed_at
    # "Step 6:" -> "Step 6/9:" (and "Steps 4-7:" -> "Steps 4-7/9:"), so that how far along a run is reads at a glance
    message = re.sub(r'^(Steps? \d+(?:-\d+)?):', rf'\1/{PIPELINE_STEPS}:', message)
    print(f"[+{elapsed_s:9.3f} s] {message}")
    if message.startswith("Step 9"):
        compute_window.close()                         # (the computing is done: the viewer opens next)
    else:
        compute_window.set(message)


def open_with_default_viewer(path):
    """Open `path` with the OS's default application for its file type
    (e.g. Photos on Windows for a .png) — a separate process entirely from
    this script, so the image stays open (and interactive: zoom, close
    whenever the user wants) regardless of what this script does afterward
    or whether it's still running. Best-effort; prints a warning instead of
    raising if the platform isn't recognized or the launch fails."""
    try:
        system = platform.system()
        if system == 'Windows':
            os.startfile(path)
        elif system == 'Darwin':
            subprocess.Popen(['open', str(path)])
        else:
            subprocess.Popen(['xdg-open', str(path)])
    except Exception as e:
        print(f"Could not open {path} in the default viewer ({e}).")


def darken_color(color, factor=0.55):
    """Return an RGB tuple that is a darker version of `color` (any matplotlib
    color spec), by scaling its RGB channels toward black."""
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)


def add_category_id_labels(adata, ax, color_key, fontsize=6, darken_factor=0.55):
    """Label each category of adata.obs[color_key] at its UMAP centroid with
    its leading numeric ID (e.g. '30 Astro-Epen' -> '30'), in a font colored
    as a darker version of that category's plotted color. Each label is
    clipped to `ax`'s own box (clip_on=True) — a no-op for this function's
    original static, full-view callers (every centroid is always within a
    full-view axes' own bounds), but essential for a zoomed-in interactive
    view (see show_interactive_umap_window's 'All Subclasses' mode): without
    it, a label whose centroid falls outside the current xlim/ylim still
    renders at its mapped screen position regardless, which — for a window
    where another panel sits immediately outside that axes' box — can bleed
    text into whatever's next door instead of just vanishing off-screen.

    Returns the list of Text artists created (empty if there was nothing to
    label), so an interactive caller can keep hold of them — the interactive
    viewer rescales their font size on zoom, in step with its dots."""
    if color_key not in adata.obs.columns or 'X_umap' not in adata.obsm:
        return []
    obs_col = adata.obs[color_key]
    if not isinstance(obs_col.dtype, pd.CategoricalDtype):
        obs_col = obs_col.astype('category')
    categories = obs_col.cat.categories
    colors = adata.uns.get(f'{color_key}_colors')
    if colors is None or len(colors) != len(categories):
        print(f"No stored colors found for '{color_key}'; skipping ID labels.")
        return []

    coords = adata.obsm['X_umap']
    obs_values = obs_col.to_numpy()
    labels = []
    for category, color in zip(categories, colors):
        mask = obs_values == category
        if not mask.any():
            continue
        centroid = coords[mask].mean(axis=0)
        match = re.match(r'^\s*(\d+)', str(category))
        label = match.group(1) if match else str(category)
        labels.append(ax.text(
            centroid[0], centroid[1], label,
            fontsize=fontsize, color=darken_color(color, darken_factor),
            fontweight='bold', ha='center', va='center', clip_on=True,
        ))
    return labels


def save_umap_by_group(adata, color_key, save_dpi, output_dir, group_size=20):
    """Split adata.obs[color_key]'s categories into groups of up to
    `group_size`, saving one UMAP plot per group (group01.png/.svg,
    group02.png/.svg, ... into `output_dir`), each with only that group's
    cells colored and every other cell left at scanpy's na_color gray.

    This exists because scanpy silently colors *every* cell uniform gray —
    no per-category color at all — once a categorical column has more
    categories than it's willing to assign distinct colors to (observed
    around 103; well past the ~300+ subclasses or ~1000+ supertypes this
    atlas has). Restricting each individual plot's color column to just one
    group of categories (everything else set to NaN via the Categorical
    constructor's `categories=` argument, which scanpy then renders as
    na_color) keeps each plot's category count comfortably under that
    limit — it's also a lot less cluttered than trying to tell 300+ colors
    apart on a single plot regardless of whether scanpy can render them.

    No-op (prints a message) if `color_key` isn't in adata.obs. `save_dpi`
    should match the resolution the rest of Step 7's plots are saved at.

    When more than one group is needed (i.e. more than `group_size`
    categories), and every expected group*.png already exists in
    `output_dir`, asks whether to skip regenerating them entirely — with a
    single group this is fast enough not to bother asking."""
    if color_key not in adata.obs.columns:
        print(f"'{color_key}' not available on adata.obs; skipping grouped UMAP plots.")
        return

    obs_col = adata.obs[color_key]
    if not isinstance(obs_col.dtype, pd.CategoricalDtype):
        obs_col = obs_col.astype('category')

    # Sort categories by their leading numeric ID (same convention used by
    # save_group_spatial_maps) so groups are stable/meaningful rather than
    # whatever order pandas happened to assign.
    ids = extract_leading_numeric_id(pd.Series(obs_col.cat.categories))
    sorted_categories = [obs_col.cat.categories[i] for i in ids.sort_values(kind='stable').index]
    groups = [sorted_categories[i:i + group_size] for i in range(0, len(sorted_categories), group_size)]

    output_dir = Path(output_dir)
    # Regenerating is cheap with a single group (<= group_size categories),
    # but with many groups it means one sc.pl.umap() call per group, which
    # adds up — so when there's more than one, offer to skip entirely if
    # every expected output file is already sitting there from a previous
    # run, rather than silently redoing all of it every time.
    if len(groups) > 1:
        expected_files = [output_dir / f'group{i:02d}.png' for i in range(1, len(groups) + 1)]
        if all(f.exists() for f in expected_files):
            if prompt_yes_no(
                f"Found {len(expected_files)} existing {color_key} UMAP group plot(s) in "
                f"{output_dir}; skip regenerating them? (y/n): "
            ):
                print(f"Skipping {color_key} UMAP group plots; using the existing files in {output_dir}.")
                return

    output_dir.mkdir(parents=True, exist_ok=True)
    # Same reasoning as save_group_spatial_maps: clear stale files from a
    # previous run before writing, in case this run produces fewer groups.
    for pattern in ('group*.png', 'group*.svg'):
        for stale in output_dir.glob(pattern):
            stale.unlink()

    temp_col = f'__{color_key}_group'
    for group_idx, group_categories in enumerate(groups, start=1):
        # Categorical(..., categories=group_categories) maps any value not
        # in that list to NaN — exactly the "everything else gray" split
        # this function is for.
        adata.obs[temp_col] = pd.Categorical(obs_col, categories=group_categories)
        # temp_col is the same *name* on every iteration, but a different
        # set of categories each time (a different-sized group, generally)
        # — scanpy's sc.pl.umap doesn't always notice that and recompute:
        # if adata.uns[f'{temp_col}_colors'] is still sitting there from
        # the *previous* group's call, it can get reused as-is rather than
        # regenerated for this group's own (often different) category
        # count. Harmless when consecutive groups both have the full
        # group_size categories, but the last group of a level is usually
        # smaller — e.g. this file's subclass level ends in a 2-category
        # final group right after four 20-category ones — and reusing a
        # stale 20-color array there is exactly what produced add_
        # category_id_labels' "No stored colors found ... skipping ID
        # labels" warning: colors existed, just for the wrong-sized
        # category list. Clearing it first forces a fresh, correctly-sized
        # array every time.
        adata.uns.pop(f'{temp_col}_colors', None)
        ax = sc.pl.umap(adata, color=temp_col, size=2, show=False)
        add_category_id_labels(adata, ax, temp_col)
        group_path = output_dir / f'group{group_idx:02d}.png'
        ax.figure.savefig(group_path, dpi=save_dpi, bbox_inches='tight')
        ax.figure.savefig(group_path.with_suffix('.svg'), bbox_inches='tight')
        plt.close(ax.figure)
        print(f"Saved {color_key} UMAP group {group_idx}/{len(groups)} "
              f"({len(group_categories)} categories) to {group_path}.")
    del adata.obs[temp_col]


# Thresholds shared by every "Export DEGs" call — a gene only makes the
# exported list if its nominal (not multiple-testing-corrected — see
# compute_de_genes' own docstring for why) p-value is below DE_P_VALUE_MAX
# *and* its linear fold change clears one of the two DE_FOLD_CHANGE_*
# bounds (up or down).
DE_P_VALUE_MAX = 0.05
DE_FOLD_CHANGE_MIN_UP = 1.1
DE_FOLD_CHANGE_MAX_DOWN = 0.9






def mean_log2_expression(work_adata, mask, log1p_base=None):
    """Mean log2 expression — the mean *of the already-log-transformed
    values* (not the log of the mean) for the cells where `mask` is True,
    one value per gene, as a Series indexed by work_adata.var_names.

    `work_adata.X` must already be log-transformed (natural-log log1p if
    `log1p_base` is None, matching scanpy's own sc.pp.log1p default, or
    already-log2 if `log1p_base` is 2) — this only reads it, it doesn't
    normalize or log anything itself, since callers need different
    treatment before this point (raw counts need normalizing first; the
    imputed dataset is already log2 on disk and must NOT be logged again).
    Converted to log2 units here regardless of the source base, so the
    result always means what its name says no matter which caller it came
    from. Shared by compute_de_genes (its own mean_log2_group/mean_log2_
    rest columns) and export_expression's plain per-ID average, which
    needs the same number with no statistical test alongside it."""
    X = work_adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.asarray(X)
    to_log2 = np.log2(log1p_base if log1p_base is not None else np.e)
    mask = np.asarray(mask)
    return pd.Series(X[mask].mean(axis=0) * to_log2, index=work_adata.var_names)


def compute_de_genes(work_adata, group_mask, log1p_base=None):
    """Differentially expressed genes for the cells where `group_mask` is
    True versus the rest of `work_adata`'s own cells, via scanpy's own
    sc.tl.rank_genes_groups(method='wilcoxon') — the Wilcoxon rank-sum
    test, matching what Seurat's FindMarkers uses by default, since this
    lab's earlier DE work was done there (comparable results, not a switch
    to a different statistical convention just because the tool changed).

    `work_adata.X` must already be log-transformed — this only runs the
    test, it doesn't normalize or log anything itself, since the two
    callers need different treatment before this point (raw counts need
    normalizing first; the imputed dataset is already log2 on disk and
    must NOT be normalized/logged again). `log1p_base` tells scanpy what
    base that log used: None for natural-log log1p (scanpy's own
    sc.pp.log1p default), or 2 for already-log2 data. This matters — it's
    what scanpy's own fold-change math inverts the log with internally
    (adata.uns['log1p']['base']) — get it wrong and every reported fold
    change is quietly computed against the wrong base, still self-
    consistent-looking but not actually the linear-scale ratio it claims
    to be.

    Returns a DataFrame with columns gene_symbol/ensembl_id/mean_log2_group/
    mean_log2_rest/log2_fold_change/fold_change/p_value/p_value_adj,
    already filtered to p_value < DE_P_VALUE_MAX and fold_change outside
    [DE_FOLD_CHANGE_MAX_DOWN, DE_FOLD_CHANGE_MIN_UP], sorted by fold_change
    descending — most upregulated first, most downregulated last, not
    grouped by effect size regardless of direction. p_value_adj
    (Benjamini-Hochberg) is reported for reference alongside the nominal
    p_value the filter actually uses, since with hundreds to thousands of
    genes tested at once, a flat p<0.05 on the nominal value alone will
    let some fraction through by chance — worth checking that column
    before treating a borderline hit as real.

    Filtering on the *nominal* p-value (not p_value_adj) is a deliberate
    choice to match the caller's exact spec, not an oversight — the
    multiple-testing-corrected value is included specifically so it's not
    lost, not to quietly override that choice."""
    work_adata = work_adata.copy()  # rank_genes_groups writes into .uns; never mutate the caller's own object
    if work_adata.X.dtype == np.float16:
        # The imputed dataset (and possibly a counts layer) can be stored
        # as float16 on disk to save space — scipy.sparse won't even
        # construct a float16 matrix at all, and scanpy's own rank_genes_
        # groups rejects a *dense* float16 array too (raises
        # NotImplementedError('float16'), which is genuinely all the info
        # that exception carries — no context, just the dtype's own str()).
        # float32 is plenty of precision for a log-expression value; this
        # only ever needs to be a working copy for this one DE run, not a
        # permanent upcast of anything on disk.
        work_adata.X = work_adata.X.astype(np.float32)
    if log1p_base is not None:
        work_adata.uns['log1p'] = {'base': log1p_base}
    else:
        work_adata.uns.pop('log1p', None)  # ensures scanpy assumes natural log, not a stale base carried over from elsewhere
    group_mask = np.asarray(group_mask)
    work_adata.obs['_de_group'] = pd.Categorical(np.where(group_mask, 'specified', 'rest'))

    sc.tl.rank_genes_groups(work_adata, groupby='_de_group', groups=['specified'], reference='rest',
                             method='wilcoxon')
    result = sc.get.rank_genes_groups_df(work_adata, group='specified')

    # (result['logfoldchanges'] needs no such conversion — scanpy's own
    # rank_genes_groups already reports that in log2 unconditionally, using
    # log1p_base internally only to invert the log correctly before
    # recomputing the ratio; see this function's own docstring.)
    mean_group = mean_log2_expression(work_adata, group_mask, log1p_base)
    mean_rest = mean_log2_expression(work_adata, ~group_mask, log1p_base)

    # result['names'] comes from work_adata.var_names — in this atlas's own
    # convention that's the Ensembl gene ID (ENSMUSG...), not the common
    # name (Th, Drd2, ...), which instead lives in the separate
    # var['gene_symbol'] column (same column find_gene_index already
    # searches elsewhere in this file). Both are exported, clearly
    # distinguished by name, rather than the exported column being called
    # 'gene_symbol' while actually holding the Ensembl ID.
    gene_symbol = work_adata.var['gene_symbol'].reindex(result['names']).values

    out = pd.DataFrame({
        'gene_symbol': gene_symbol,
        'ensembl_id': result['names'].values,
        'mean_log2_group': mean_group.loc[result['names']].values,
        'mean_log2_rest': mean_rest.loc[result['names']].values,
        'log2_fold_change': result['logfoldchanges'].values,
        'fold_change': 2.0 ** result['logfoldchanges'].values,
        'p_value': result['pvals'].values,
        'p_value_adj': result['pvals_adj'].values,
    })
    out = out[
        (out['p_value'] < DE_P_VALUE_MAX)
        & ((out['fold_change'] > DE_FOLD_CHANGE_MIN_UP) | (out['fold_change'] < DE_FOLD_CHANGE_MAX_DOWN))
    ]
    return out.sort_values('fold_change', ascending=False).reset_index(drop=True)


def show_own_dataset_window(own, abc_cache, project_folder=None):
    """The Tk UMAP viewer on a dataset of the user's own (see abc_tk_own.py); returns when it is closed."""
    from abc_tk_umap import show_umap_window_tk
    host_names = (
        'CACHE_DIR', 'UMAP_SAVE_DPI', 'SECTION_MAPS_SAVE_DPI', 'THREED_SHELL_POINTS', 'load_viewer_settings', 'save_viewer_settings',
        'prompt_viewer_settings_dialog', 'umap_params_key', 'umap_params_label', 'umap_params_are_default', 'umap_params_token',
        'normalize_umap_params', 'write_last_umap_params', 'ensure_umap_variant_csv', 'load_umap_variant', 'open_with_default_viewer',
        'unique_export_paths', 'prompt_deg_comparison_mode', 'materialize_subset', 'compute_de_genes', 'mean_log2_expression',
        'ensure_ccf_downloaded', 'load_ccf_for_cells', 'load_parcellation_names', 'run_with_loading_window', 'brain_outline_segments',
        'open_3d_group_view', 'build_3d_edges', 'write_3d_group_html', 'category_rank_shape', 'visible_points_only', 'sc',
        'read_nifti_volume', 'CCF_VOLUME_DIRECTORY', 'CCF_VOLUME_FILE_NAME')
    host = {name: globals()[name] for name in host_names}
    host['settings_dialog'] = prompt_viewer_settings_dialog
    host['umap_circle_sizes'] = UMAP_CIRCLE_SIZE_REGISTRY
    title = f"{project_folder.name} \u2014 {own.title}" if project_folder is not None else own.title
    show_umap_window_tk(own.viewer_adata(), abc_cache, run_folder=project_folder, title=title, host=host, own=own)


def show_interactive_umap_window(adata, abc_cache, imputed_state=None, adata_backed=None, rois=None,
                                  run_folder=None, panel_zoom=False):
    """Interactive UMAP viewer, opened once the embedding (adata.obsm[
    'X_umap']) is computed. A vertical hand-rolled radio on the right picks
    what colors the scatter:
      - All Subclasses: every cell colored by its 'subclass' category (via
        scanpy's own categorical coloring — same convention, including its
        color-count cap, as the old optional group plots) with each
        category's ID number labeled at its centroid (add_category_id_
        labels) instead of a legend, which wouldn't fit hundreds of
        categories anyway. No query needed.
      - Single Subclass: one subclass ID, typed below, highlighted red
        against gray — the old hardcoded-268 static plot, now with any ID.
      - Gene: raw counts from adata.layers['counts'], log1p'd — not
        adata.X itself, which by this point has been normalized/scaled/
        clipped for PCA and no longer holds interpretable expression
        values.
      - Imputed Gene: the imputed dataset's log2 values, matched to these
        cells by ID, same convention as the single-section picker's gene
        view.
    A text box below the radio (hidden for All Subclasses, which has
    nothing to type) takes the corresponding query and Enter (or the Show
    button) re-colors the plot; the gene box autocompletes against
    whichever gene panel the current mode implies. Gene/Imputed Gene's
    colorbar draws into its own permanently-reserved, fixed-size axes
    (cax=, not ax=) — so the main plot's size never changes switching
    between color-by modes, colorbar or not.

    Left of the UMAP scatter is a brain-section grid: one real Axes per
    section actually present in `adata` (post filter/subsample — i.e.
    exactly the sections that ended up in this run's ROI/whole-section
    selection), each showing that section's cells in gray, all sharing one
    physical scale (so an anatomically small section isn't zoomed to look
    as large as a full coronal one — see the SPAN_PERCENTILE block below)
    but with an independent pan/zoom per panel. Hovering a UMAP point
    (after a short hold, like prompt_subregion_selection's own hover-to-
    identify-cell — deferred, not looked up on every raw mouse-move event,
    so a ~200k-cell nearest-point search never runs more often than the
    mouse actually settles) highlights that cell in both the UMAP itself
    and, if its section is currently in the grid, that one section panel
    (every other panel's highlight is hidden) — and fills in a status line
    below with its class/subclass/supertype/cluster. Both the UMAP scatter
    and every section panel support mouse-wheel zoom and left-drag pan,
    each independently (see pannable_axes/on_scroll_zoom/on_press_pan
    below) — zooming out is clamped per-axes to that axes' own "home"
    (fully-zoomed-out) extent, computed once at setup.

    `adata_backed`, if given, is the original unfiltered backed AnnData
    from Step 1 (the same one save_roi_map/save_group_spatial_maps already
    use, for the same reason) — each section panel's background shows
    *every* cell of that section (via get_section_labels/load_section_
    spatial_coords on `adata_backed`, not `adata`), not just the ones that
    survived this run's own filtering/subsampling. Omit for a panel
    background made of only this run's own cells, same as before this was
    added. `rois`, if given, is the same list of {'section','x_min',
    'x_max','y_min','y_max'} dicts prompt_section_selection_gui collects —
    drawn as dashed orange rectangles on the matching section panel(s),
    same convention as that picker's own ROI borders; entries for whole-
    section picks (blank x_min/x_max/y_min/y_max) are skipped, since
    there's no meaningful rectangle to draw for "the whole section".

    `imputed_state`, if given, is the same {'adata':, 'load_thread':,
    'load_error':} dict prompt_section_selection_gui uses (see
    ensure_imputed_gene_dataset_loaded) — reused here so 'Imputed Gene'
    doesn't pay for a second load if the picker session already loaded it.
    Pass None (a fresh dict is created) if it wasn't, or the picker's GUI
    wasn't used at all (console fallback).

    Blocks until the window is closed, same as every other window in this
    file (see show_figure_blocking)."""
    if imputed_state is None:
        imputed_state = {'adata': None, 'load_thread': None, 'load_error': None}

    no_umap = bool(adata.uns.get('no_umap'))                       # (a project made without a UMAP: the section maps and the genes)
    if 'X_umap' not in adata.obsm and not no_umap:
        print("No UMAP embedding available; skipping the interactive UMAP viewer.")
        return

    # The Tk viewer (abc_tk_umap.py). It opens on the UMAP parameters the run was last viewed with, like the matplotlib one below.
    tk_params = dict(UMAP_PARAM_DEFAULTS)
    if run_folder is not None and not no_umap:
        last_params = read_last_umap_params(run_folder)
        if not umap_params_are_default(last_params):
            try:
                variant_coords, variant_leiden = load_umap_variant(umap_variant_csv_path(run_folder, last_params), adata.obs_names)
                adata.obsm['X_umap'] = variant_coords
                if variant_leiden is not None:
                    adata.obs[LEIDEN_KEY] = pd.Categorical(variant_leiden)
                else:
                    adata.obs.drop(columns=[LEIDEN_KEY], inplace=True, errors='ignore')
                tk_params = last_params
            except Exception as e:
                log_status(f"Step 9: Could not reopen the saved UMAP variant ({e}); using the default embedding.")
    tk_title = "Section map viewer" if no_umap else "Interactive UMAP viewer"
    if not no_umap and not umap_params_are_default(tk_params):
        tk_title = f"{tk_title} ({umap_params_label(tk_params)})"
    if run_folder is not None:
        tk_title = f"{Path(run_folder).name} — {tk_title}"
    from abc_tk_umap import show_umap_window_tk
    # The functions of this file the Tk viewer uses (it cannot import this file: it is the program).
    host_names = (
        'CACHE_DIR', 'UMAP_SAVE_DPI', 'SECTION_MAPS_SAVE_DPI', 'THREED_SHELL_POINTS', 'load_viewer_settings', 'save_viewer_settings',
        'prompt_viewer_settings_dialog', 'umap_params_key', 'umap_params_label', 'umap_params_are_default', 'umap_params_token',
        'normalize_umap_params', 'write_last_umap_params', 'ensure_umap_variant_csv', 'load_umap_variant', 'open_with_default_viewer',
        'unique_export_paths', 'prompt_deg_comparison_mode', 'materialize_subset', 'compute_de_genes', 'mean_log2_expression',
        'ensure_ccf_downloaded', 'load_ccf_for_cells', 'load_parcellation_names', 'run_with_loading_window', 'brain_outline_segments',
        'open_3d_group_view', 'build_3d_edges', 'write_3d_group_html', 'category_rank_shape', 'visible_points_only', 'sc',
        'read_nifti_volume', 'CCF_VOLUME_DIRECTORY', 'CCF_VOLUME_FILE_NAME', 'UMAP_PARAM_DEFAULTS')
    host = {name: globals()[name] for name in host_names}
    host['settings_dialog'] = prompt_viewer_settings_dialog
    host['preload_imputed'] = preload_imputed_gene_dataset
    host['umap_circle_sizes'] = UMAP_CIRCLE_SIZE_REGISTRY
    show_umap_window_tk(adata, abc_cache, imputed_state=imputed_state, adata_backed=adata_backed, rois=rois, run_folder=run_folder,
                        ensure_imputed=ensure_imputed_gene_dataset_loaded, title=tk_title, host=host, params=None if no_umap else tk_params,
                        panel_zoom=panel_zoom)           # (a project without a UMAP has no UMAP parameters: Settings then does not offer them)
    return


def open_no_umap_project(abc_cache, run_folder, backed, imputed_state, section_series=None, mask=None, panel_zoom=False):
    """Open a project without a UMAP (its folder holds project.json saying so, and roi_coords.csv if it has ROIs). Its cells are picked out of the
    atlas's tables from what project.json records (cell type, sections, ROIs), and their genes are read from the atlas's own file when wanted (see
    abc_project.py): there is no processed file. `mask`: the cells, if the caller has just found them."""
    global selection_confirmed_at
    run_folder = Path(run_folder)
    info = abc_project.read_project(run_folder)
    cell_type, sections, rois = abc_project.selection_of(info)
    if rois:
        run_frame = read_run_frame(run_folder)
        if run_frame != coordinate_frame_name():
            print(f"This project's ROIs were drawn on the {run_frame} coordinates, so its cells are placed by those (the setting is for new projects).")
            FRAME_STATE['override'] = run_frame
            prepare_coordinate_frame(abc_cache)
    if abc_project.release_differs(info, abc_cache):
        print(f"Note: this project was made from the atlas release {(info.get('atlas') or {}).get('release')}; the cache is now on "
              f"{abc_project.atlas_release(abc_cache)}. Its cells are picked out of the current tables.")

    def build():
        picked = mask if mask is not None else abc_project.select_cells(backed, abc_cache, cell_type, sections, rois, section_series)
        return abc_project.build_adata(backed, abc_cache, picked)
    adata = run_with_loading_window(f"Opening {run_folder.name}.\nPicking out its cells from the atlas.", build)
    selection_confirmed_at = time.perf_counter()                 # (log_status counts the time from here)
    log_status(f"Opening {run_folder.name} ({adata.n_obs:,} cells, no UMAP).")
    run_with_loading_window("Writing this project's statistics.txt (it had none).",
                            lambda: write_run_statistics(run_folder, adata, abc_cache, rois=rois, overwrite=False,
                                                         note="made when this project without a UMAP was opened and had none"))
    log_status("Opening the viewer...")
    try:
        show_interactive_umap_window(adata, abc_cache, imputed_state=imputed_state, adata_backed=backed, rois=rois, run_folder=run_folder,
                                     panel_zoom=panel_zoom)
    except Exception as e:
        log_status(f"Could not open the viewer: {e}")
    log_status("Viewer closed.")


def create_no_umap_project(abc_cache, out_folder, backed, section_series, imputed_state, cell_type_selection, sections_suffix, region_suffix,
                           whole_sections_selected, rois, run_label, panel_zoom=False):
    """Make a project without a UMAP (the options window's "skip the UMAP"): a folder with project.json and, for ROIs, roi_coords.csv and the ROI
    map; then open it. No subsampling, no preprocessing, no PCA or UMAP is done, so it takes seconds, whatever the number of cells."""
    global selection_confirmed_at
    name = f'browse_{CELL_TYPE_SUFFIXES[cell_type_selection]}_{sections_suffix}{region_suffix}'
    run_folder = Path(out_folder) / name
    log_status(f"Making the project {name} (no UMAP).")
    mask = abc_project.select_cells(backed, abc_cache, cell_type_selection, whole_sections_selected, rois, section_series)
    n_cells = int(mask.sum())
    if n_cells == 0:
        show_error_dialog("Nothing to show", "No cells match this selection (cell type, sections and ROIs), so no project was made.")
        return
    run_folder.mkdir(parents=True, exist_ok=True)
    core_record_run_opened(run_folder)
    if rois:
        write_run_frame(run_folder)
    roi_section_set = {roi['section'] for roi in rois}
    whole_only_sections = [s for s in (whole_sections_selected or []) if s not in roi_section_set] if rois else []
    if rois or whole_only_sections:
        whole_section_rows = [{'section': s, 'x_min': np.nan, 'x_max': np.nan, 'y_min': np.nan, 'y_max': np.nan} for s in whole_only_sections]
        pd.DataFrame(list(rois) + whole_section_rows)[['section', 'x_min', 'x_max', 'y_min', 'y_max']].to_csv(run_folder / 'roi_coords.csv', index=False)
    abc_project.write_project(run_folder, cell_type_selection, whole_sections_selected, rois, label=run_label, has_umap=False, n_cells=n_cells,
                              abc_cache=abc_cache, frame=coordinate_frame_name())
    log_status(f"{n_cells:,} cells in the project.")
    if rois:
        try:
            save_roi_map(backed, abc_cache, section_series, rois, run_folder / 'roi_map.png')
        except Exception as e:
            log_status(f"Could not save the ROI map ({e}); the project is made without it.")
    open_no_umap_project(abc_cache, run_folder, backed, imputed_state, section_series=section_series, mask=mask, panel_zoom=panel_zoom)


def run_session(abc_cache, session_state, raw_backed_cache, imputed_state=None):
    """One full user session: show the startup panel, then either open
    the section/ROI picker plus the full UMAP-computation pipeline, or
    load a previously computed run directly — ending, either way, with
    the interactive UMAP viewer. Returns normally once that viewer
    window closes, or the user cancels partway through (a
    UserCancelledSelection from the section/ROI picker) — never calls
    sys.exit() itself; only prompt_startup_panel's own "closed without
    choosing" case does that, which is deliberately the one place that
    ends the whole program (see the session loop below, which just
    calls this repeatedly). raw_backed_cache is a plain {} the caller
    keeps alive across calls, so the raw MERFISH h5ad (see
    get_cached_raw_backed_adata) is only ever opened once per process,
    not once per session. imputed_state is likewise kept alive by the
    caller, so the imputed gene dataset (see
    ensure_imputed_gene_dataset_loaded) is loaded at most once per process;
    None gives this session its own, as before.

    session_state is a mutable dict carrying {'out_folder': Path} across
    calls: the startup panel lets the user change the output folder, and
    writing it back here (rather than returning it) means every one of this
    function's several early-return paths picks up the new folder for the
    next session without each having to thread it through."""
    # log_status (defined at module level, outside this function) reads
    # selection_confirmed_at as a plain global — it used to be a real
    # module-level variable, set once by this same top-level script code;
    # now that the assignments below run inside this function, they'd
    # otherwise create a *local* invisible to log_status instead of
    # updating the global it actually reads.
    global selection_confirmed_at
    if imputed_state is None:
        imputed_state = {'adata': None, 'load_thread': None, 'load_error': None}
    out_folder = session_state['out_folder']
    try:
        startup_action, cell_type_selection, load_existing_folder, chosen_out_folder = prompt_startup_panel(out_folder)
        if chosen_out_folder is not None:
            out_folder = session_state['out_folder'] = chosen_out_folder
    except Exception as e:
        print(f"GUI startup panel unavailable ({e}); falling back to console prompt.")
        cell_type_selection = prompt_cell_type_selection()
        startup_action, load_existing_folder = 'new_run', None

    # "Browse All Neurons" opens the viewer with Zoom whole panel ticked; an All neurons dataset already made is opened like any saved run.
    browse_all_neurons = startup_action == 'all_neurons'
    if browse_all_neurons and load_existing_folder is not None:
        startup_action = 'load_existing'

    # The Tk launch screen doesn't ask for the cell type: it is asked once the sections and ROIs are confirmed (below), and the picker
    # shows every cell until then.
    if cell_type_selection is not None:
        subset_suffix = CELL_TYPE_SUFFIXES[cell_type_selection]
        print(f"Selected cell type subset: {cell_type_selection}")

    if startup_action == 'own_dataset':
        # A dataset of the user's own (an .h5ad of single-cell or single-nucleus RNA-seq): the UMAP viewer without the atlas (abc_tk_own.py).
        from abc_tk_own import prepare_own_dataset
        try:
            prepared = prepare_own_dataset(Path(load_existing_folder), run_with_loading_window, CACHE_DIR, out_folder)
        except Exception as e:
            show_error_dialog("Could not open the dataset", f"{type(e).__name__}: {e}")
            return
        if prepared is not None:
            own, project_folder = prepared
            selection_confirmed_at = time.perf_counter()                 # (log_status counts the time from here)
            log_status(f"Opening {own.path.name} in the UMAP viewer (project folder {project_folder})...")
            try:
                show_own_dataset_window(own, abc_cache, project_folder)
            except Exception as e:
                log_status(f"Could not open the UMAP viewer: {e}")
            log_status("UMAP viewer closed.")
        return

    from abc_tk_own import read_project, reopen_own_project
    if startup_action == 'load_existing':
        core_record_run_opened(load_existing_folder)                 # (the Load Existing list starts on the run opened last)
    if startup_action == 'load_existing' and read_project(load_existing_folder) is not None:
        # The saved run is the project folder of a dataset of the user's own: it opens from the dataset it names, with no questions.
        try:
            prepared = reopen_own_project(load_existing_folder, run_with_loading_window, CACHE_DIR)
        except Exception as e:
            show_error_dialog("Could not open the dataset", f"{type(e).__name__}: {e}")
            return
        if prepared is not None:
            own, project_folder = prepared
            selection_confirmed_at = time.perf_counter()                 # (log_status counts the time from here)
            log_status(f"Opening {own.path.name} in the UMAP viewer (project folder {project_folder})...")
            try:
                show_own_dataset_window(own, abc_cache, project_folder)
            except Exception as e:
                log_status(f"Could not open the UMAP viewer: {e}")
            log_status("UMAP viewer closed.")
        return

    resume_info = resume_selection = None
    if startup_action == 'resume':
        # A run that stopped before its UMAP was made: its ROIs are in its folder, and its name says the rest. Made again as a new run would be.
        resume_info = parse_unfinished_run(load_existing_folder)
        try:
            resume_rois, resume_whole = load_rois_csv(Path(load_existing_folder) / 'roi_coords.csv')
        except Exception as e:
            show_error_dialog("Could not finish the project", f"Its ROIs could not be read: {e}")
            return
        if resume_info is None:
            show_error_dialog("Could not finish the project", f"The folder's name ({Path(load_existing_folder).name}) is not one that a run makes, so how it was made is not known.")
            return
        if resume_whole:
            resume_selection = (resume_whole, 'sections-' + '-'.join(sanitize_section_token(s) for s in resume_whole), resume_rois, None)
        else:
            resume_selection = (None, 'allsections', resume_rois, None)
    FRAME_STATE['override'] = None
    prepare_coordinate_frame(abc_cache)
    if startup_action == 'resume' and resume_selection[2] and read_run_frame(load_existing_folder) != coordinate_frame_name():
        # (the ROIs were drawn on the other coordinates: they select the cells by those)
        FRAME_STATE['override'] = read_run_frame(load_existing_folder)
        prepare_coordinate_frame(abc_cache)
    if startup_action == 'load_existing':
        # Skips section/ROI selection and UMAP computation (Steps 1-8) entirely
        # — the chosen folder's own cached files (verified by prompt_startup_
        # panel before it let this choice through) already have everything
        # Step 9 needs. Still needs the raw backed h5ad (same read Step 0 below
        # does for the 'new_run' path) since that's what the viewer's
        # adata_backed argument uses for each section's full spatial
        # background, independent of any specific run's own selection —
        # get_cached_raw_backed_adata reuses it across sessions instead of
        # re-opening the file every time.
        print(f"Loading cached run from {load_existing_folder}...")
        if not abc_project.has_umap(load_existing_folder):
            # A project made without a UMAP: nothing to load but the atlas's own tables (see open_no_umap_project).
            backed, section_series = run_with_loading_window(
                f"Loading the MERFISH brain atlas (4+ million cells).\nOpening {Path(load_existing_folder).name}.",
                lambda: (lambda b: (b, get_section_labels(b, abc_cache)))(get_cached_raw_backed_adata(raw_backed_cache, abc_cache)))
            open_no_umap_project(abc_cache, load_existing_folder, backed, imputed_state, section_series=section_series, panel_zoom=browse_all_neurons)
            return
        required_files = required_cached_run_files(load_existing_folder)
        processed_h5ad_path = required_files['Processed AnnData']
        roi_csv_path = required_files['ROI coordinates']

        # Same reasoning as the 'new_run' path below: both reads can be slow
        # and nothing else is on screen while they run.
        def load_raw_and_processed():
            backed = get_cached_raw_backed_adata(raw_backed_cache, abc_cache)
            processed = None
            if processed_h5ad_path.exists():
                print(f"Loading processed AnnData from {processed_h5ad_path}...")
                processed = anndata.read_h5ad(processed_h5ad_path)
            return backed, processed

        adata_backed, adata = run_with_loading_window(
            f"Loading saved run {Path(load_existing_folder).name}.", load_raw_and_processed)
        if adata is None:
            # Deleted (or never kept) — rebuildable from this run's own
            # umap_coords.csv without reopening the section/ROI picker, since
            # that CSV already records exactly which cells the run covered.
            print(f"{processed_h5ad_path} not found; rebuilding it from this run's saved UMAP coordinates.")
            try:
                adata = rebuild_processed_adata(load_existing_folder, adata_backed, processed_h5ad_path)
            except Exception as e:
                print(f"Could not rebuild the processed AnnData ({e}); returning to the startup panel.")
                return
        if roi_csv_path.exists():
            rois, _whole_sections = load_rois_csv(roi_csv_path)
            print(f"Loaded {len(rois)} ROI(s) from {roi_csv_path}.")
            run_frame = read_run_frame(load_existing_folder)
            if rois and run_frame != coordinate_frame_name():
                print(f"This run's ROIs were drawn on the {run_frame} coordinates, so its cells are placed by those (the setting is for new runs).")
                FRAME_STATE['override'] = run_frame
                prepare_coordinate_frame(abc_cache)
        else:
            # Not written at all for a whole-brain run (no ROIs, no explicit
            # whole-section picks) — see the main script's own roi_csv_path-
            # writing block further down, guarded by "if rois or whole_only_
            # sections". Not an error; just means this run covers everything.
            rois = []
            print(f"Warning: {roi_csv_path} not found; assuming this run covers the whole brain (no ROIs).")

        # A run saved without its cell-type annotations (class/subclass/...) would show "No 'subclass' column" for
        # every level: restore them from the atlas metadata, once, and save them back.
        run_with_loading_window(
            "Checking this saved run's cell-type annotations (adds any that are missing).",
            lambda: backfill_taxonomy_columns(adata, abc_cache, Path(load_existing_folder) / 'umap_coords.csv',
                                              processed_h5ad_path))

        # Runs computed before Leiden was added have no clustering saved —
        # offer to backfill it rather than silently showing the Leiden level
        # as unavailable. Declining is fine; nothing below depends on it.
        offer_leiden_backfill(adata, Path(load_existing_folder) / 'umap_coords.csv', processed_h5ad_path)

        # Reference point for log_status()'s elapsed-time prefix, same as the
        # normal flow sets right after its own section/ROI selection step.
        selection_confirmed_at = time.perf_counter()
        # A run saved before statistics.txt existed gets one now (a file that is there is left as it is).
        run_with_loading_window("Writing this run's statistics.txt (it had none).",
                                lambda: write_run_statistics(Path(load_existing_folder), adata, abc_cache, rois=rois, overwrite=False,
                                                             note="made when this saved run was opened and had none"))
        log_status("Step 9: Opening the interactive UMAP viewer...")
        try:
            show_interactive_umap_window(
                adata, abc_cache,
                imputed_state=imputed_state,
                adata_backed=adata_backed, rois=rois, run_folder=Path(load_existing_folder), panel_zoom=browse_all_neurons,
            )
        except Exception as e:
            log_status(f"Step 9: Could not open the interactive UMAP viewer: {e}")
        log_status("Step 9: Interactive UMAP viewer closed.")
        return

    # ---------------------------------------------------------------------------
    # 0. Load data (as provided)
    # ---------------------------------------------------------------------------
    # abc_cache = ...  # assumed already created/configured elsewhere in your workflow

    # Both reads are slow on a cold start and the startup panel is already
    # gone, so run them behind a loading window (see run_with_loading_window)
    # rather than leaving the screen empty. get_section_labels is nearly free
    # when brain_section_label is already in adata.obs, but reads the whole
    # metadata CSV when it isn't.
    def load_raw_and_section_labels():
        loaded = get_cached_raw_backed_adata(raw_backed_cache, abc_cache)
        return loaded, get_section_labels(loaded, abc_cache)

    # The first time on a computer the atlas is still downloading (gigabytes, started at program start): far longer than reading it.
    downloading = atlas_prefetch['thread'] is not None and atlas_prefetch['thread'].is_alive() and atlas_prefetch['total_gb'] >= 1
    if downloading:
        loading_message = (f"Waiting for the MERFISH brain atlas to finish downloading ({atlas_prefetch['total_gb']:.0f} GB, one time only).\n"
                           "This can take 10 minutes or more; see the console window for progress.")
    else:
        loading_message = "Loading the MERFISH brain atlas (4+ million cells).\nThe first load after startup can take a minute."
    adata, section_series = run_with_loading_window(loading_message, load_raw_and_section_labels)
    print(f"Step 0: Loaded backed AnnData with {adata.n_obs} cells x {adata.n_vars} genes.")

    # ---------------------------------------------------------------------------
    # Section subset selection
    # ---------------------------------------------------------------------------
    all_neurons_run = startup_action == 'all_neurons'
    resume_run = startup_action == 'resume'
    if all_neurons_run:
        # "Browse All Neurons" (abc_tk_launcher): every section, no ROIs, so there is nothing to pick.
        selected_sections, sections_suffix, rois, imputed_adata_from_picker = None, 'allsections', [], None
    elif resume_run:
        selected_sections, sections_suffix, rois, imputed_adata_from_picker = resume_selection
        print(f"Finishing {Path(load_existing_folder).name} from its {len(rois)} saved ROI(s).")
    else:
        try:
            print("Attempting to open the clickable section-thumbnail picker...")
            selected_sections, sections_suffix, rois, imputed_adata_from_picker = prompt_section_selection_gui(
                adata, abc_cache, section_series, out_folder=out_folder, cell_type_selection=cell_type_selection or 'All',
                imputed_state=imputed_state,
            )
        except UserCancelledSelection:
            print("Section/ROI selection cancelled; returning to the startup panel.")
            return
        except Exception as e:
            print(f"GUI section picker unavailable ({e}); falling back to console prompt.")
            selected_sections, sections_suffix = prompt_section_selection(section_series)
            rois = []
            imputed_adata_from_picker = None
    # One window now settles everything the analysis would otherwise stop to ask (cell type, the most cells in the UMAP, the ROI run's
    # label, the optional plots), so that it can then run unattended. Without a window (no display), the questions come one by one at the
    # console as they always did (run_options stays None).
    run_options = None
    whole_brain = selected_sections is None and not rois           # (no sections and no ROIs picked: every section is used)
    # The cells of the selection are counted on a thread while the options window is open (it shows "counting..." until they are known).
    count_state = {'done': False, 'counts': None, 'error': None}

    def count_in_background():
        try:
            count_state['counts'] = count_selected_cells(adata, abc_cache, section_series, selected_sections, rois)
            print(f"Cells in the selection: {count_state['counts']}")
        except Exception as e:
            count_state['error'] = e
            print(f"Could not count the cells in the selection ({e}); the large-project question will be asked regardless.")
        finally:
            count_state['done'] = True

    if all_neurons_run:
        # Every neuron of every section, named "All neurons", with no UMAP (a project of seconds: create_no_umap_project): nothing to ask, nothing to count.
        from abc_tk_launcher import ALL_NEURONS_LABEL
        run_options = {'cell_type': 'Neurons', 'make_umap': False, 'cap': None, 'label': ALL_NEURONS_LABEL, 'optional_plots': False}
        count_state['done'] = True
    elif resume_run:
        # The cells, the project name and the limit are in the folder's name: nothing to ask.
        run_options = {'cell_type': resume_info['cell_type'], 'cap': resume_info['cap'], 'label': resume_info['label'], 'optional_plots': False}
        count_state['done'] = True
    else:
        threading.Thread(target=count_in_background, daemon=True).start()
    if not (all_neurons_run or resume_run):
        try:
            from abc_tk_launcher import choose_run_options
            def folder_for(cell_type, project_name, cap, make_umap=True):
                """The folder a project with these choices is saved in (as the analysis below will name it), for the options window to show."""
                sections_part, region_part, _roi_sections = region_naming(sections_suffix, rois)
                name = f'{"umap" if make_umap else "browse"}_{CELL_TYPE_SUFFIXES[cell_type]}_{sections_part}{region_part}'
                if sanitize_run_label(project_name):
                    name += f'_{sanitize_run_label(project_name)}'
                counted = count_state['counts']
                n = counted.get(cell_type) if counted else None
                if make_umap and cap is not None and n is not None and n > cap:
                    name += f'_sub{cap}'
                return Path(out_folder) / name

            run_options = choose_run_options(cell_type_selection or 'Neurons', folder_for=folder_for, whole_brain=whole_brain, count_state=count_state)
            if run_options is None:
                print("Analysis options cancelled; returning to the startup panel.")
                return
        except Exception as e:
            print(f"Analysis options window unavailable ({e}); the questions will be asked here instead.")
            if whole_brain and not confirm_whole_brain():
                print("Whole-brain run not confirmed; returning to the startup panel.")
                return
            if cell_type_selection is None:
                cell_type_selection = prompt_cell_type_selection()
    if run_options is not None:
        cell_type_selection = run_options['cell_type']
        print(f"Analysis options: {run_options}")
    subset_suffix = CELL_TYPE_SUFFIXES[cell_type_selection]
    print(f"Selected cell type subset: {cell_type_selection}")
    # Reference point for log_status()'s elapsed-time prefix, used from here on
    # for every status update in the processing pipeline below.
    selection_confirmed_at = time.perf_counter()
    if selected_sections is None:
        log_status("Selected sections: all")
    else:
        log_status(f"Selected sections: {selected_sections}")
    # Saved off before the ROI override below can replace selected_sections, so
    # the ROI CSV (further down) can still record whole-section picks even when
    # ROIs were also drawn on other sections in the same picker session.
    whole_sections_selected = selected_sections

    # ---------------------------------------------------------------------------
    # If any ROIs were drawn (double-click on a thumbnail in the grid picker),
    # they take priority over the whole-section selection above: the run
    # processes exactly the ROI cells, across however many sections they span.
    # ---------------------------------------------------------------------------
    sections_suffix, region_suffix, roi_sections = region_naming(sections_suffix, rois)
    if rois:
        selected_sections = roi_sections
        log_status(f"Using {len(rois)} ROI(s) across {len(roi_sections)} section(s) "
                   "(overrides whole-section selection above).")
    # The project name (typed in the options window; at the console, asked for ROI runs only, as it always was) ends the folder name.
    run_label = sanitize_run_label(run_options['label']) if run_options is not None else (prompt_roi_run_label() if rois else '')
    if run_label:
        region_suffix += f'_{run_label}'
        log_status(f"Adding the project name '{run_label}' to this run's folder name.")

    if run_options is not None and not run_options.get('make_umap', True):
        # The options window said to skip the UMAP: no subsampling, preprocessing or embedding (and so no size limit): a project that is only
        # its cells, picked out of the atlas's tables (see create_no_umap_project).
        create_no_umap_project(abc_cache, out_folder, adata, section_series, imputed_state, cell_type_selection, sections_suffix, region_suffix,
                               whole_sections_selected, rois, run_label, panel_zoom=all_neurons_run)
        return

    # ---------------------------------------------------------------------------
    # Early cache-hit check: look for an existing UMAP run matching this exact
    # section/cell-type/ROI selection *before* paying for Step 2's metadata
    # merge (a full multi-million-row CSV read+join) or Step 3's filtering
    # (which reads the matched cells off disk) — both become pure waste if
    # we're about to just reload a previously computed embedding anyway.
    # ---------------------------------------------------------------------------
    # subsample_target isn't known yet at this point — it's normally decided by
    # Step 4, from the post-filter cell count this is specifically trying to
    # avoid computing — so this can't look for one exact folder name the way
    # csv_path does further down (once subsample_target *is* known). Instead it
    # globs for any run folder for this selection, with any subsample cap or
    # none, and only acts when that's unambiguous (exactly one match); with
    # zero or multiple matches this just falls through to the unchanged Step
    # 2/3/4 pipeline, which finds the *exact* match (if any) the normal way,
    # once subsample_target is actually known.
    base_run_suffix = f'{subset_suffix}_{sections_suffix}{region_suffix}'
    run_dir_pattern = re.compile(rf'^umap_{re.escape(base_run_suffix)}(?:_sub(\d+))?$')
    early_cache_candidates = []
    for candidate in Path(out_folder).glob(f'umap_{base_run_suffix}*'):
        match = run_dir_pattern.match(candidate.name)
        if match is not None and candidate.is_dir() and (candidate / 'umap_coords.csv').exists():
            early_cache_candidates.append((candidate, int(match.group(1)) if match.group(1) else None))

    if run_options is not None:
        # The options window has already settled the cap, so only a run made with that cap is a match: a single existing run with another cap (the
        # 200,000 of an earlier time, say) must not be reloaded in place of the 800,000 just asked for. The cells in the selection are known (the
        # window waited for the count), so whether this run would be capped at all is too: it is when there are more cells than the cap.
        cap = run_options['cap']
        counted = count_state['counts'].get(run_options['cell_type']) if count_state.get('counts') else None
        if counted is not None:
            expected = cap if (cap is not None and counted > cap) else None
            early_cache_candidates = [(c, s) for c, s in early_cache_candidates if s == expected]
        else:
            early_cache_candidates = [(c, s) for c, s in early_cache_candidates if s == cap]

    early_cache_hit = len(early_cache_candidates) == 1
    subsample_target = None
    if early_cache_hit:
        early_cache_dir, subsample_target = early_cache_candidates[0]
        log_status(f"Found a single existing UMAP run matching this selection at {early_cache_dir} "
                   "before filtering — skipping Step 2's metadata merge and Step 3's filtering "
                   "entirely (see Step 2/3/4 below).")

    # ---------------------------------------------------------------------------
    # 1. Keep the backed AnnData for now — defer the full disk read
    # ---------------------------------------------------------------------------
    # adata stays backed (X is not read from disk yet) through metadata merging
    # and filtering below. Both only touch adata.obs, which is a real in-memory
    # DataFrame even in backed mode, so they don't need X materialized. This
    # means the actual disk read (in materialize_subset(), used by the filters
    # below and by subsampling further down) only ever pulls in the cells that
    # survive filtering/subsampling, instead of loading the entire multi-million
    # -cell dataset into memory just to filter or subsample most of it away
    # immediately after.
    adata_backed = adata  # keep the original, unfiltered backed AnnData for the ROI map below

    # ---------------------------------------------------------------------------
    # 2. Attach cell type metadata (class/subclass/cluster/section), if available
    # ---------------------------------------------------------------------------
    # The ABC atlas ships separate metadata tables (cell_metadata.csv) with
    # cluster/class/subclass/section annotations keyed by cell_label. If you
    # already have this loaded elsewhere, merge it in here.
    if early_cache_hit:
        log_status("Step 2: Skipped — reusing the cached run found above; its own saved "
                   "umap_coords.csv already carries whichever metadata columns it needs "
                   "(see Step 2/3/4's reuse branch further down).")
        metadata_cols = []
        color_key = None
    else:
        log_status("Step 2: Attempting to load and merge cell type metadata...")
        try:
            # The table is read from the CSV only the first time (a gigabyte; the better part of a minute): after that from a cache that
            # is already in the order of adata.obs, so the columns are simply added (see load_annotations_aligned).
            annotations = compute_window.run(f"Step 2/{PIPELINE_STEPS}: Reading the cell types of the atlas...",
                                             lambda: load_annotations_aligned(adata.obs.index, abc_cache))
            metadata_cols = [c for c in annotations.columns if c not in adata.obs.columns]
            # Adding columns to adata.obs (rather than subsetting adata to the
            # intersection first) keeps adata backed — subsetting a backed AnnData
            # with .copy() isn't allowed (see materialize_subset's note), and would
            # force a full disk read here regardless, before any filtering has had
            # a chance to shrink what needs to be read.
            n_matched = int(annotations[SECTION_COL if SECTION_COL in annotations.columns else annotations.columns[0]].notna().sum())
            new_obs = adata.obs.copy()
            for c in metadata_cols:
                new_obs[c] = annotations[c].to_numpy()
            adata.obs = new_obs
            color_key = 'class' if 'class' in adata.obs.columns else None
            log_status(f"Step 2: Merged metadata for {n_matched} of {adata.n_obs} cells. Coloring by '{color_key}'.")
        except Exception as e:
            log_status(f"Step 2: Could not load/merge cell type metadata ({e}); proceeding without it.")
            metadata_cols = []
            color_key = None
            # Surfaced as a dialog too: the run otherwise carries on and is saved with no class/subclass/supertype/
            # cluster columns, and the viewer later says "No 'subclass' column" with the cause long gone from view.
            show_error_dialog("Cell-type annotations could not be attached",
                              f"Merging the atlas cell-type metadata failed:\n\n{e}\n\nThe run continues without "
                              "class/subclass/supertype/cluster. They are added automatically the next time this saved "
                              "run is opened.")

    # ---------------------------------------------------------------------------
    # 3. Filter to the selected cell type, section(s), and ROIs (if any)
    # ---------------------------------------------------------------------------
    cells_after_filter = None
    if early_cache_hit:
        log_status("Step 3: Skipped — the cached run found above already implies this exact "
                   "cell-type/section/ROI filter; materializing straight from its saved "
                   "umap_coords.csv further down picks out the same cells directly.")
    else:
        section_desc = 'all sections' if selected_sections is None else f'{len(selected_sections)} selected section(s)'
        log_status(f"Step 3: Filtering to cell type subset '{cell_type_selection}' and {section_desc}...")

        def filter_cells(adata=adata):
            adata = filter_by_cell_type(adata, cell_type_selection)
            adata = filter_by_sections(adata, section_series, selected_sections)
            if rois:
                adata = filter_by_rois(adata, abc_cache, rois)
            return adata
        adata = compute_window.run(f"Step 3/{PIPELINE_STEPS}: Picking out the cells ({cell_type_selection}, {section_desc})...", filter_cells)
        log_status(f"Step 3: {adata.n_obs} cells remain after filtering.")
        cells_after_filter = int(adata.n_obs)           # (before a random subset, if one is taken: it goes into statistics.txt)

    # ---------------------------------------------------------------------------
    # 4. Decide whether/how much to subsample
    # ---------------------------------------------------------------------------
    # The chosen cap becomes part of the output folder name below, so it has to
    # be decided before that folder is created — a run capped at a different
    # number of cells is a different run, not something to merge into a folder
    # an earlier run with a different cap already claimed.
    N_SUBSAMPLE_DEFAULT = 200_000
    if early_cache_hit:
        # subsample_target was already set (from the matched folder's own
        # '_sub{N}' suffix, or left None if it had none) by the early cache-hit
        # check above — recomputing it here from adata.n_obs would be wrong
        # anyway, since adata is still the full, unfiltered dataset at this
        # point (Step 3 was skipped) rather than the post-filter count this
        # decision is normally based on.
        log_status(f"Step 4: Skipped — reusing subsample_target={subsample_target} from the matched run.")
    else:
        subsample_target = None  # None means "use every filtered cell, no subsampling"
        if run_options is not None:
            cap = run_options['cap']
            if cap is not None and adata.n_obs > cap:
                subsample_target = cap
                log_status(f"Step 4: {adata.n_obs} cells exceed the limit of {cap}: using a pseudorandom subset of {cap} (fixed seed: the same cells each time).")
        elif adata.n_obs > N_SUBSAMPLE_DEFAULT:
            subsample_target = prompt_subsample_choice(adata.n_obs, N_SUBSAMPLE_DEFAULT)

    run_suffix = f'{subset_suffix}_{sections_suffix}{region_suffix}'
    if subsample_target is not None:
        run_suffix += f'_sub{subsample_target}'
    # Every output file for this run is grouped into its own subfolder, named
    # after what used to be their shared filename prefix, so a directory
    # listing of out_folder shows one entry per run instead of a handful of
    # similarly-prefixed files per run.
    session_prefix = f'umap_{run_suffix}'
    run_folder = Path(out_folder) / session_prefix
    run_folder.mkdir(parents=True, exist_ok=True)
    csv_path = run_folder / 'umap_coords.csv'
    try:
        abc_project.write_project(run_folder, cell_type_selection, whole_sections_selected, rois, label=run_label, has_umap=True,
                                 cap=subsample_target, abc_cache=abc_cache, frame=coordinate_frame_name())
    except Exception as e:
        print(f"Could not write {abc_project.PROJECT_FILE} ({e}); the project works without it.")

    # Whole-section picks (clicked to select the entire section, no ROI drawn
    # on it) are only worth recording here if the run also has actual ROIs —
    # otherwise every section in the run is already a whole-section pick, and
    # that's already fully captured by the run's folder name.
    roi_section_set = {roi['section'] for roi in rois}
    whole_only_sections = [s for s in (whole_sections_selected or []) if s not in roi_section_set] if rois else []

    core_record_run_opened(run_folder)
    roi_csv_path = run_folder / 'roi_coords.csv'
    if rois:
        write_run_frame(run_folder)
    # Skip the confirm-overwrite prompt (and the write itself) entirely when
    # roi_coords.csv already records exactly this session's ROI/whole-section
    # selection — most relevantly when reusing a cached run's selection
    # unchanged, where re-prompting to overwrite an identical file is just
    # noise.
    roi_selection_is_unchanged = roi_selection_unchanged(roi_csv_path, rois, whole_only_sections)

    if rois or whole_only_sections:
        if roi_selection_is_unchanged:
            log_status(f"ROI/whole-section selection unchanged since last save; skipping {roi_csv_path}.")
        elif run_options is not None or confirm_overwrite(roi_csv_path):       # (with the options window, files of this run are simply replaced)
            # Whole-section rows get blank x_min/x_max/y_min/y_max, marking
            # "every cell in this section" as distinct from an actual
            # sub-region — load_rois_csv() reads that blankness back out to
            # route them to the right place (selected sections vs. rois) when
            # this file is reloaded via the picker's 'Load ROIs' button.
            whole_section_rows = [
                {'section': s, 'x_min': np.nan, 'x_max': np.nan, 'y_min': np.nan, 'y_max': np.nan}
                for s in whole_only_sections
            ]
            roi_csv_rows = pd.DataFrame(list(rois) + whole_section_rows)
            roi_csv_rows[['section', 'x_min', 'x_max', 'y_min', 'y_max']].to_csv(roi_csv_path, index=False)
            log_status(f"Saved {len(rois)} ROI(s) and {len(whole_only_sections)} whole section(s) to {roi_csv_path}.")
        else:
            log_status(f"Skipped overwriting {roi_csv_path}.")

    if rois:
        roi_map_path = run_folder / 'roi_map.png'
        if roi_selection_is_unchanged and roi_map_path.exists() and roi_map_path.with_suffix('.svg').exists():
            log_status(f"ROI selection unchanged since last save; skipping {roi_map_path}.")
        elif run_options is not None or confirm_overwrite(roi_map_path):
            # Uses adata_backed, the reference to the original, unfiltered
            # AnnData saved off in Step 1 before `adata` gets reassigned by
            # filtering/subsampling below, so the map shows each ROI's full
            # section for context, not just the cells that end up surviving
            # downstream.
            save_roi_map(adata_backed, abc_cache, section_series, rois, roi_map_path)
        else:
            log_status(f"Skipped overwriting {roi_map_path}.")

    # Check for a previously computed UMAP coordinates CSV for this exact cell
    # type + section + subsample-cap selection. If found, skip the heavy
    # subsample/preprocessing/PCA/neighbors/UMAP steps and just reuse it. Also
    # true already if early_cache_hit (Steps 2/3/4 above were skipped on the
    # strength of this same check, run early) — csv_path.exists() is
    # necessarily also true then (same file, found via a different route), but
    # checking early_cache_hit first avoids a redundant stat() call.
    log_status("Looking for existing UMAP for current selection")
    found_existing_umap = early_cache_hit or csv_path.exists()

    if found_existing_umap:
        # -----------------------------------------------------------------
        # 4-7 (skipped). Reuse a previously saved UMAP coordinates CSV.
        # -----------------------------------------------------------------
        log_status(f"Steps 4-7: Found existing UMAP coordinates at {csv_path}; "
                   "skipping subsample/preprocessing/PCA/neighbors/UMAP and reusing it.")
        # converters={0: str} keeps the cell-ID index column as text; otherwise
        # pandas infers it as float64 and silently rounds off the trailing
        # digits of these ~19-digit cell IDs (float64 only has ~15-17 sig figs).
        umap_coords = pd.read_csv(csv_path, index_col=0, converters={0: str})
        keep_mask = adata.obs.index.isin(umap_coords.index)
        # materialize_subset() reads only the matched cells from disk if adata
        # is still backed (e.g. the whole-brain/all-cell-types case, where none
        # of the filters in Step 3 triggered a load — or early_cache_hit, where
        # Step 3 never ran at all) instead of materializing the full dataset
        # just to throw most of it away right here.
        adata = materialize_subset(adata, keep_mask)
        # Only actually adds anything when Step 2's merge was skipped
        # (early_cache_hit) — umap_coords.csv already carries whichever of
        # these Step 8 saved, so there's no need to re-read/re-join the full
        # multi-million-row metadata CSV just to get them back. A no-op
        # (columns already present) in the normal found-late case, where Step
        # 2's own merge already put them there.
        for meta_col in ('class', 'subclass', 'supertype', 'cluster', SECTION_COL):
            if meta_col in umap_coords.columns and meta_col not in adata.obs.columns:
                adata.obs[meta_col] = umap_coords.loc[adata.obs.index, meta_col].values
        # Restored from the CSV rather than recomputed: this branch skips
        # sc.pp.neighbors entirely, so the graph Leiden needs doesn't exist
        # here. Simply absent for runs whose CSV predates this column (or
        # whose Step 6 couldn't run it) — nothing downstream requires it.
        if LEIDEN_KEY in umap_coords.columns and LEIDEN_KEY not in adata.obs.columns:
            adata.obs[LEIDEN_KEY] = pd.Categorical(
                umap_coords.loc[adata.obs.index, LEIDEN_KEY].astype(str).values
            )
        metadata_cols = [c for c in ('class', 'subclass', 'supertype', 'cluster', SECTION_COL)
                          if c in adata.obs.columns]
        color_key = 'class' if 'class' in adata.obs.columns else None
        adata.obsm['X_umap'] = umap_coords.loc[adata.obs.index, ['UMAP1', 'UMAP2']].to_numpy()
        # The interactive UMAP viewer's Gene mode (Step 9) needs raw counts —
        # normally saved off in Step 5, which every cache-hit path (early or
        # late) skips entirely; adata.X is still the untouched raw counts read
        # from disk here (materialize_subset above never normalizes/transforms
        # it), so this just claims that as layers['counts'] rather than
        # re-deriving it.
        adata.layers['counts'] = adata.X.copy()
        # Same backfill offer as the 'load_existing' startup path — this
        # branch reuses an equally old CSV, just reached from the new-run
        # flow instead. Done after layers['counts'] is claimed above, so the
        # copy it works on already carries everything a fresh run would.
        offer_leiden_backfill(adata, csv_path, CACHE_DIR / f'{session_prefix}_processed.h5ad')
        log_status(f"Aligned {adata.n_obs} cells to the saved UMAP embedding. "
                   f"Coloring by '{color_key}'.")
    else:
        if subsample_target is not None:
            log_status(f"Step 4: Subsampling from {adata.n_obs} to {subsample_target} cells...")
            rng = np.random.RandomState(0)
            keep_mask = np.zeros(adata.n_obs, dtype=bool)
            keep_mask[rng.choice(adata.n_obs, size=subsample_target, replace=False)] = True
            # Same reasoning as the cache-hit branch above: if nothing in Step 3
            # triggered a load (e.g. the whole-brain/all-cell-types case), adata
            # is still backed here, so this is where the actual disk read
            # happens — and it only reads the subsampled cells, not the full
            # unfiltered dataset first.
            adata = compute_window.run(f"Step 4/{PIPELINE_STEPS}: Reading {subsample_target:,} cells' expression from the atlas...",
                                       lambda: materialize_subset(adata, keep_mask))
            log_status(f"Step 4: Subsampled to {adata.n_obs} cells.")
        elif adata.isbacked:
            log_status(f"Step 4: Using all {adata.n_obs} filtered cells (no subsampling); loading into memory...")
            adata = compute_window.run(f"Step 4/{PIPELINE_STEPS}: Reading {adata.n_obs:,} cells' expression from the atlas...", adata.to_memory)
        else:
            log_status(f"Step 4: Using all {adata.n_obs} filtered cells (no subsampling).")

        # ---------------------------------------------------------------------------
        # 5. Preprocessing
        # ---------------------------------------------------------------------------
        # Basic QC filtering
        log_status("Step 5: Filtering cells (min_counts=20) and genes (min_cells=5)...")

        def qc_filter():
            sc.pp.filter_cells(adata, min_counts=20)
            sc.pp.filter_genes(adata, min_cells=5)
        compute_window.run(f"Step 5/{PIPELINE_STEPS}: Filtering cells and genes...", qc_filter)
        log_status(f"Step 5: After filtering, {adata.n_obs} cells x {adata.n_vars} genes remain.")

        # Normalize + log-transform (MERFISH panels are small, ~500 genes,
        # so we skip highly-variable-gene selection and use all genes for PCA)
        log_status("Step 5: Normalizing, log-transforming, and scaling...")

        def normalize():
            adata.layers['counts'] = adata.X.copy()          # (the raw counts, kept)
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            # sc.pp.scale zero-centers the data, which densifies the sparse matrix.
            # Harmless here since the MERFISH panel is small (~500 genes vs. tens of
            # thousands for full transcriptome data), so just silence the warning.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    'ignore', message='zero-centering a sparse array/matrix densifies it',
                    category=UserWarning,
                )
                sc.pp.scale(adata, max_value=10)
        compute_window.run(f"Step 5/{PIPELINE_STEPS}: Normalizing, log-transforming and scaling...", normalize)
        log_status("Step 5: Preprocessing complete.")

        # ---------------------------------------------------------------------------
        # 6. PCA -> neighbors -> UMAP
        # ---------------------------------------------------------------------------
        log_status("Step 6: Running PCA (n_comps=50)...")
        pca_start = time.perf_counter()
        compute_window.run(f"Step 6/{PIPELINE_STEPS}: Running PCA...", lambda: sc.tl.pca(adata, n_comps=50, svd_solver='arpack'))
        log_status(f"Step 6: PCA complete in {format_duration(time.perf_counter() - pca_start)}.")
        log_status("Step 6: Computing neighbor graph (n_neighbors=15, n_pcs=50)...")
        compute_window.run(f"Step 6/{PIPELINE_STEPS}: Computing the neighbor graph...", lambda: sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50))
        log_status("Step 6: Running UMAP (this is usually the longest step)...")
        umap_start = time.perf_counter()
        compute_window.run(f"Step 6/{PIPELINE_STEPS}: Running UMAP (this is usually the longest step)...", lambda: sc.tl.umap(adata))
        log_status(f"Step 6: UMAP complete in {format_duration(time.perf_counter() - umap_start)}.")

        # Leiden clustering — scanpy's own de-novo grouping of these cells,
        # entirely independent of the Allen ABC taxonomy columns (class/
        # subclass/supertype/cluster) merged in Step 2. Stored under
        # LEIDEN_KEY rather than 'leiden' straight into obs so it can never
        # be confused with ABC's own 'cluster' column, which is an
        # annotation, not a computed clustering.
        #
        # Reuses the neighbor graph built just above, so this is cheap
        # relative to everything preceding it. Non-fatal if it can't run:
        # leidenalg/igraph are optional extras, and a missing one shouldn't
        # cost the user the whole (already-computed) UMAP.
        log_status("Step 6: Running Leiden clustering...")
        leiden_start = time.perf_counter()
        try:
            # flavor/n_iterations pinned explicitly: leaving flavor unset
            # emits a FutureWarning in scanpy 1.12 about the default
            # changing, and 'igraph' with n_iterations=2 is what that
            # warning steers callers toward.
            compute_window.run(f"Step 6/{PIPELINE_STEPS}: Finding clusters (Leiden)...",
                               lambda: sc.tl.leiden(adata, key_added=LEIDEN_KEY, flavor='igraph', n_iterations=2,
                                                    resolution=LEIDEN_RESOLUTION))
            log_status(f"Step 6: Leiden found {adata.obs[LEIDEN_KEY].nunique()} clusters in "
                       f"{format_duration(time.perf_counter() - leiden_start)}.")
        except Exception as e:
            log_status(f"Step 6: Leiden clustering unavailable ({e}); continuing without it.")
            # Surfaced as a dialog too, not just in the console: the run
            # otherwise carries on to a perfectly normal-looking viewer that
            # silently has no Leiden level, with the only clue several
            # screens back in the terminal.
            show_error_dialog("Leiden clustering failed", leiden_failure_message(e))

    # ---------------------------------------------------------------------------
    # 7. Plot
    # ---------------------------------------------------------------------------
    log_status(f"Step 7: Saving UMAP plot to {run_folder}...")
    sc.settings.figdir = run_folder
    DEFAULT_DPI_SAVE = 150  # scanpy's built-in default
    SAVE_DPI = 3 * DEFAULT_DPI_SAVE  # triple saved-PNG resolution
    # sc.set_figure_params (scanpy=True, its own default) doesn't just apply
    # dpi_save — it also overwrites several *on-screen* rcParams (figure.dpi
    # defaults to 80, among others scanpy styles) for the rest of the process,
    # not just for the plots this step is about to save. Step 9's own window
    # below targets a fixed *pixel* size on screen regardless of dpi, but its
    # text is sized in points (rendered as points * dpi/72 pixels) — so a
    # lingering rcParams['figure.dpi'] of 80 instead of matplotlib's true
    # default of 100 made that same point size render visibly smaller relative
    # to the window there, only in a run that actually went through this step.
    # Snapshotting rcParams now and restoring them once this step's own
    # plotting is done (right before Step 8) keeps that styling scoped to just
    # these saved files, so Step 9 always sees the same rcParams regardless of
    # whether this step ran at all — matching "load existing", which skips it
    # entirely.
    rc_params_before_step7 = dict(matplotlib.rcParams)
    sc.set_figure_params(dpi_save=SAVE_DPI)

    main_plot_path = run_folder / 'class_plot.png'
    if color_key is not None:
        ax = sc.pl.umap(adata, color=color_key, size=2, show=False)
        log_status(f"Step 7: Labeling each '{color_key}' group with its ID number...")
        add_category_id_labels(adata, ax, color_key)
    else:
        ax = sc.pl.umap(adata, size=2, show=False)

    ax.figure.savefig(main_plot_path, dpi=SAVE_DPI, bbox_inches='tight')
    ax.figure.savefig(main_plot_path.with_suffix('.svg'), bbox_inches='tight')
    log_status(f"Step 7: Plot saved to {main_plot_path}.")
    plt.close(ax.figure)
    # Not auto-opened in a viewer window (see Step 9's interactive viewer for
    # the on-screen equivalent) — this just exists as a file on disk, same as
    # the optional subclass/supertype group plots below, until/unless someone
    # opens it themselves.

    # Optional: a separate set of plots of every subclass/supertype (not just
    # whichever single subclass gets highlighted interactively in Step 9 below)
    # colored on the UMAP, plus their spatial location.
    if 'subclass' in adata.obs.columns:
        if run_options['optional_plots'] if run_options is not None else prompt_yes_no(
            "Step 7: Generate optional plots of UMAP subclasses, and spatial plot of subclass/supertype? (y/N): ",
            default=False,
        ):
            # Still generated and saved to disk for the record even though
            # Step 9's interactive viewer now covers the same coloring on
            # demand — never auto-opened in a window (same as main_plot_path
            # above, as of this turn); just a file until someone opens it.
            # One UMAP plot per group isn't optional here — scanpy renders a
            # categorical column entirely gray, with no per-category colors at
            # all, once it has more categories than colors it's willing to
            # assign (subclass/supertype both have far more than that); see
            # save_umap_by_group's docstring. Groups of 20 also just reads
            # better than trying to tell hundreds of colors apart at once.
            for group_col in ('subclass', 'supertype'):
                if group_col not in adata.obs.columns:
                    log_status(f"Step 7: No '{group_col}' column available; skipping {group_col} UMAP plots.")
                    continue
                umap_group_dir = run_folder / f'umap_by_{group_col}'
                log_status(f"Step 7: Plotting {group_col} UMAP groups to {umap_group_dir}...")
                save_umap_by_group(adata, group_col, SAVE_DPI, umap_group_dir)

            # Also show where these subclasses — and, separately, supertypes —
            # actually sit in space, restricted to the sections that were part
            # of this run's ROI/whole-section selection, not the whole,
            # unfiltered dataset. Uses adata_backed (the original, unfiltered
            # backed AnnData from Step 1) rather than the processed/filtered/
            # subsampled `adata`, so each section's full cell population is
            # available for spatial context, same as save_roi_map above.
            roi_section_set = {roi['section'] for roi in rois} if rois else set()
            whole_sections_for_maps = [s for s in (whole_sections_selected or []) if s not in roi_section_set]
            for group_col in ('subclass', 'supertype'):
                if group_col not in adata.obs.columns:
                    log_status(f"Step 7: No '{group_col}' column available; skipping {group_col} spatial maps.")
                    continue
                group_maps_dir = run_folder / f'spatial_maps_{group_col}'
                log_status(f"Step 7: Generating {group_col} spatial maps in {group_maps_dir}...")
                save_group_spatial_maps(
                    adata_backed, abc_cache, section_series, rois, whole_sections_for_maps, group_maps_dir,
                    group_col=group_col, cell_type_selection=cell_type_selection,
                )
        else:
            log_status("Step 7: Skipped the optional subclass/supertype plots.")

    # Restore whatever sc.set_figure_params changed above — see its own comment.
    matplotlib.rcParams.update(rc_params_before_step7)

    # ---------------------------------------------------------------------------
    # 8. Save results
    # ---------------------------------------------------------------------------
    log_status("Step 8: Assembling UMAP coordinates table...")
    umap_coords = pd.DataFrame(
        adata.obsm['X_umap'],
        index=adata.obs.index,
        columns=['UMAP1', 'UMAP2']
    )
    for c in metadata_cols:
        if c in adata.obs.columns:
            umap_coords[c] = adata.obs[c].values
    # Saved alongside the ABC taxonomy columns above, so a cached run can be
    # reloaded with its clustering intact. Absent whenever Step 6 skipped or
    # failed (leidenalg/igraph missing), and absent from every CSV written
    # before this column existed — hence a membership check here and on
    # every read, rather than assuming it's present.
    if LEIDEN_KEY in adata.obs.columns:
        umap_coords[LEIDEN_KEY] = adata.obs[LEIDEN_KEY].values

    if found_existing_umap:
        log_status(f"Step 8: Reused existing UMAP coordinates from {csv_path}; not rewriting it.")
    else:
        log_status(f"Step 8: Writing UMAP coordinates to {csv_path}...")
        umap_coords.to_csv(csv_path)

    # Optionally save the processed AnnData for reuse (kept in the local cache dir, not out_folder)
    h5ad_out_path = CACHE_DIR / f'{session_prefix}_processed.h5ad'
    log_status(f"Step 8: Writing processed AnnData to {h5ad_out_path}...")
    compute_window.run(f"Step 8/{PIPELINE_STEPS}: Saving the run...", lambda: adata.write_h5ad(h5ad_out_path))

    log_status(f"Done. UMAP computed for {adata.n_obs} cells ({cell_type_selection}, {sections_suffix}).")
    log_status(f"Saved: {main_plot_path}, "
               f"{csv_path}, "
               f"{h5ad_out_path}")

    # ---------------------------------------------------------------------------
    # 9. Interactive UMAP viewer (subclass highlight / gene / imputed gene)
    # ---------------------------------------------------------------------------
    compute_window.run(f"Step 8/{PIPELINE_STEPS}: Writing this run's statistics.txt...",          # (in the same window, not one of its own)
                       lambda: write_run_statistics(run_folder, adata, abc_cache, rois=rois, whole_sections=whole_only_sections if rois else None,
                                                    cells_selected=cells_after_filter, overwrite=True, note="made with this run"))
    log_status("Step 9: Opening the interactive UMAP viewer...")
    try:
        show_interactive_umap_window(
            adata, abc_cache,
            # The shared state already holds whatever the picker loaded
            # (imputed_adata_from_picker is that same object, when set).
            imputed_state=imputed_state,
            adata_backed=adata_backed, rois=rois, run_folder=run_folder, panel_zoom=browse_all_neurons,
        )
    except Exception as e:
        log_status(f"Step 9: Could not open the interactive UMAP viewer: {e}")
    log_status("Step 9: Interactive UMAP viewer closed.")


raw_backed_cache = {}
# The imputed gene dataset, once loaded, kept for the life of the process —
# same reasoning as raw_backed_cache. Passed down to every window that can
# load it (section picker, ROI picker, UMAP viewer), so a later session or a
# reloaded previous run reuses it instead of loading it again. Safe to share:
# it's opened read-only (backed='r'), and consumers only read it or copy
# subsets out of it.
session_imputed_state = {'adata': None, 'load_thread': None, 'load_error': None}
# Seeded with the last folder used (persisted across launches — see
# remember_recent_output_folder), or the historical default on a first run,
# then overridden per-session by whatever the startup panel's output-folder
# box is set to (kept in this dict so it persists across sessions — see
# run_session).
_recent_folders = load_recent_output_folders()
# First launch: a folder next to this repository (so, next to the default atlas
# location), which needs no drive letter or user name and isn't cloud-synced
# the way Documents often is. After that, the most recently used folder.
DEFAULT_OUTPUT_FOLDER = (Path.home() / 'Documents' if RUNNING_PACKAGED else DEFAULTS_DIR) / 'ABC_atlas_browser_output'
_initial_out_folder = resolve_out_folder(_recent_folders[0] if _recent_folders else str(DEFAULT_OUTPUT_FOLDER))
session_state = {'out_folder': _initial_out_folder}
startup_note("Opening the launch screen...")
while True:
    run_session(abc_cache, session_state, raw_backed_cache, session_imputed_state)
