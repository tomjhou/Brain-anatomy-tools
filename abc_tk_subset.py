"""Analyze subset: the groups switched on in the viewer's legend become a project of their own, in a subfolder of the viewer's run, with their own
PCA / UMAP / Leiden, opened in a second viewer window.

    analyze_subset(viewer)

The cells are those of the groups still colored in the All groups list, or in the Single group(s) legend (at the level shown). The raw counts of
the run are normalized, scaled and embedded again, from scratch, so the UMAP shows the structure within the subset instead of the one the whole run
had. The subfolder holds umap_coords.csv, roi_coords.csv (the run's) and subset.json (which groups), like a run, and the processed AnnData is kept
in the cache folder under the subfolder's name, as for a run. Selecting the same groups again opens what is saved."""
import hashlib
import json
import re
import time
import warnings
from pathlib import Path
from tkinter import messagebox

import numpy as np
import pandas as pd

import abc_core as core

MIN_SUBSET_CELLS = 100             # fewer cells than this have too little to embed
PCA_COMPONENTS = 50                # as for a run: what Analyze subset offers to begin with
MAX_PCA_COMPONENTS = 150           # and the most it will take
NEIGHBORS = 15
METADATA_COLUMNS = ("class", "subclass", "supertype", "cluster", core.SECTION_COL, "neurotransmitter")


def legend_selection(v):
    """(level, mask, names, message): the cells of the groups switched on in the legend on the screen. `message` says why there is no selection
    (nothing to choose from, none switched on, or all of them): then mask is None."""
    level = v.level
    if v.mode == "all":
        keys = list(getattr(v, "group_rows", []))
        on = [k for k in keys if k not in v.off_groups]
    elif v.mode == "single":
        keys = list((v.last_ids or {}).get("targets", []))
        on = [k for k in keys if k not in v.off_groups]
    else:
        return None, None, [], "Choose All groups or Single group(s) first: the subset is made of the groups switched on in their legend."
    if not keys:
        return None, None, [], "There are no groups in the legend to take a subset from."
    if not on:
        return None, None, [], "No group is selected. Switch on the groups to analyze in the legend (click an entry, or use All), then try again."
    if len(on) == len(keys):
        return None, None, [], ("Every group is selected, so the subset would be the whole run. Switch off the groups to leave out (click an entry, "
                                "or use None and then switch some on), then try again.")
    codes, uniques = v.category_codes(level)
    if v.mode == "all":
        chosen = np.isin(np.asarray([str(u) for u in uniques]), np.asarray([str(k) for k in on]))
    else:
        numbers = core.extract_leading_numeric_id(pd.Series(uniques)).to_numpy()
        chosen = np.isin(numbers, np.asarray(on))
    mask = (codes >= 0) & chosen[np.maximum(codes, 0)]
    names = [str(k) for k in on]
    return level, mask, names, ""


def subset_folder(run_folder, level, names, n_components=PCA_COMPONENTS):
    """The subfolder of the run for these groups: named after the level and the groups (their numbers, when there are few), with a short hash of the
    exact set so that two different selections never share one. The number of PCA components is in the name too (pc50): the same groups analyzed
    with another number are another subset."""
    numbers = [re.match(r"^\s*(\d+)", n) for n in names]
    if all(numbers) and len(names) <= 4:
        token = "-".join(sorted((m.group(1) for m in numbers), key=int))
    else:
        token = f"{len(names)}groups"
    extra = [f"pc{n_components}"]
    digest = hashlib.md5(("|".join([Path(run_folder).name, level] + sorted(names) + extra)).encode("utf-8")).hexdigest()[:8]
    return Path(run_folder) / "_".join(["subset", level, token] + extra + [digest])


def ask_pca_components(parent, n_cells, n_groups):
    """The window that asks how many principal components to use (PCA_COMPONENTS to begin with, at most MAX_PCA_COMPONENTS). Returns the number, or
    None if it was cancelled."""
    import gc
    import tkinter as tk
    from tkinter import ttk
    dlg = tk.Toplevel(parent)
    dlg.title("Analyze subset")
    dlg.transient(parent)
    dlg.resizable(False, False)
    var = tk.StringVar(dlg, str(PCA_COMPONENTS))
    result = {"value": None}
    body = ttk.Frame(dlg, padding=16)
    body.pack()
    ttk.Label(body, text=f"{n_cells:,} cells in {n_groups} group(s) will be analyzed again on their own.", wraplength=380, justify=tk.LEFT).pack(anchor="w")
    row = ttk.Frame(body)
    row.pack(anchor="w", pady=(12, 0))
    ttk.Label(row, text="Dimensions (principal components) of the PCA:").pack(side=tk.LEFT)
    spin = ttk.Spinbox(row, from_=2, to=MAX_PCA_COMPONENTS, textvariable=var, width=6)
    spin.pack(side=tk.LEFT, padx=(8, 0))
    ttk.Label(body, foreground="gray30", wraplength=380, justify=tk.LEFT, text=(
        f"From 2 to {MAX_PCA_COMPONENTS}; {PCA_COMPONENTS} is what the whole run uses. More keep more detail, and more noise. The PCA cannot use more "
        f"dimensions than the cells or genes allow.")).pack(anchor="w", pady=(6, 0))

    def ok():
        try:
            n = int(var.get())
            if not 2 <= n <= MAX_PCA_COMPONENTS:
                raise ValueError
        except ValueError:
            messagebox.showerror("Analyze subset", f"Enter a whole number from 2 to {MAX_PCA_COMPONENTS}.", parent=dlg)
            return
        result["value"] = n
        dlg.destroy()

    buttons = ttk.Frame(body)
    buttons.pack(fill=tk.X, pady=(14, 0))
    ttk.Button(buttons, text="Cancel", command=dlg.destroy).pack(side=tk.RIGHT)
    ttk.Button(buttons, text="Analyze", command=ok).pack(side=tk.RIGHT, padx=(0, 8))
    dlg.bind("<Return>", lambda e: ok())
    dlg.bind("<Escape>", lambda e: dlg.destroy())
    dlg.update_idletasks()
    dlg.geometry(f"+{parent.winfo_rootx() + 120}+{parent.winfo_rooty() + 120}")
    spin.focus_set()
    spin.selection_range(0, tk.END)
    dlg.grab_set()
    parent.wait_window(dlg)
    value = result["value"]
    var = spin = row = body = buttons = dlg = None                 # (the variable is freed here, on this thread, while Tk is up)
    gc.collect()
    return value


def embed_subset(parent, mask, n_components=PCA_COMPONENTS):
    """The subset's AnnData: the run's cells in the groups, from their raw counts, normalized, scaled and embedded (PCA, neighbors, UMAP, Leiden).
    Runs on a worker thread: it touches no window."""
    import anndata
    import scanpy as sc
    source = parent.adata
    rows = np.flatnonzero(mask)
    # (a new AnnData from the raw counts and the cells' annotations, not a copy of the run's: whatever else it holds, the embedding, the graph, its scaled
    # matrix, would only be in the way)
    obs = source.obs.iloc[rows].copy()
    obs = obs.drop(columns=[c for c in ("n_counts", core_leiden_key()) if c in obs.columns])
    sub = anndata.AnnData(X=source.layers["counts"][rows].copy(), obs=obs, var=source.var.copy())
    sc.pp.filter_genes(sub, min_cells=5)
    sub.layers["counts"] = sub.X.copy()
    sc.pp.normalize_total(sub)
    sc.pp.log1p(sub)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="zero-centering a sparse array/matrix densifies it", category=UserWarning)
        sc.pp.scale(sub, max_value=10)
    n_comps = int(min(n_components, sub.n_obs - 1, sub.n_vars - 1))
    sc.tl.pca(sub, n_comps=n_comps, svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=min(NEIGHBORS, sub.n_obs - 1), n_pcs=n_comps)
    sc.tl.umap(sub)
    try:
        sc.tl.leiden(sub, key_added=core_leiden_key(), flavor="igraph", n_iterations=2, resolution=1.0)
    except Exception as e:
        print(f"Leiden clustering is not available for the subset ({e}); continuing without it.", flush=True)
    return sub


def core_leiden_key():
    from abc_tk_umap import LEIDEN_KEY
    return LEIDEN_KEY


def save_subset(sub, folder, parent_folder, level, names, n_components=PCA_COMPONENTS):
    """umap_coords.csv, the run's roi_coords.csv, subset.json, and the processed AnnData (in the cache folder, as for a run)."""
    folder.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(sub.obsm["X_umap"], index=sub.obs.index, columns=["UMAP1", "UMAP2"])
    for column in METADATA_COLUMNS + (core_leiden_key(),):
        if column in sub.obs.columns:
            table[column] = sub.obs[column].values
    table.to_csv(folder / "umap_coords.csv")
    for name in ("roi_coords.csv", "coordinate_frame.txt"):          # (the ROIs, and the coordinates they were drawn on: the subset is opened with the same)
        source = Path(parent_folder) / name
        if source.exists():
            (folder / name).write_bytes(source.read_bytes())
    with open(folder / "subset.json", "w", encoding="utf-8") as f:
        json.dump({"parent": str(parent_folder), "level": level, "groups": names, "cells": int(sub.n_obs), "pca_components": int(n_components), "made": time.strftime("%Y-%m-%d %H:%M")}, f, indent=2)
    for frame in (sub.obs, sub.var):                      # (pandas' own string types are not written by every anndata: plain text is)
        frame.index = frame.index.astype(object)
        for column in frame.columns:
            if isinstance(frame[column].dtype, pd.StringDtype):
                frame[column] = frame[column].astype(object)
    sub.write_h5ad(core.CACHE_DIR / f"{folder.name}_processed.h5ad")


def analyze_subset(v):
    """The Analyze subset button of the viewer `v`."""
    if v.busy:
        return
    level, mask, names, message = legend_selection(v)
    if mask is None:
        messagebox.showwarning("Analyze subset", message, parent=v.root)
        return
    if int(mask.sum()) < MIN_SUBSET_CELLS:
        messagebox.showwarning("Analyze subset", f"The selected groups have only {int(mask.sum())} cells in this run: a subset needs at least "
                                                 f"{MIN_SUBSET_CELLS}. Switch on more groups and try again.", parent=v.root)
        return
    if v.run_folder is None:
        messagebox.showwarning("Analyze subset", "This viewer has no project folder to keep the subset in.", parent=v.root)
        return
    n_components = ask_pca_components(v.root, int(mask.sum()), len(names))
    if n_components is None:
        return
    folder = subset_folder(v.run_folder, level, names, n_components)
    if v.own is not None:
        return analyze_own_subset(v, level, mask, names, folder, n_components)
    processed = core.CACHE_DIR / f"{folder.name}_processed.h5ad"
    host = v.host
    import anndata
    if (folder / "umap_coords.csv").exists() and processed.exists():
        sub = host["run_with_loading_window"](f"Opening the saved subset {folder.name}.", lambda: anndata.read_h5ad(processed), show_after_s=0)
    else:
        n = int(mask.sum())
        sub = host["run_with_loading_window"](
            f"Analyzing {n:,} cells of {len(names)} {v.level_names[level].lower()} group(s): normalizing, PCA ({n_components} dimensions), neighbors, UMAP "
            f"and clustering (a minute or two). Saved in {folder}.",
            lambda: _embed_and_save(v, mask, folder, level, names, n_components), show_after_s=0)
    core.record_run_opened(folder)                       # (the Load Existing list starts on the project opened last)
    from abc_tk_umap import UmapWindow
    params = dict(host.get("UMAP_PARAM_DEFAULTS") or v.params or {})
    title = f"{Path(v.run_folder).name} / {folder.name} — Interactive UMAP viewer (subset of {sub.n_obs:,} cells)"
    win = UmapWindow(sub, v.abc_cache, v.imputed_state, v.adata_backed, getattr(v, "rois_input", []), folder, v.ensure_imputed, title,
                     parent=v.root, host=host, params=params or None)
    v.subset_windows.append(win)
    try:
        win.root.state("zoomed")
    except Exception:
        pass
    win.root.lift()


def analyze_own_subset(v, level, mask, names, folder, n_components=PCA_COMPONENTS):
    """Analyze subset on a dataset of the user's own: the cells' log2 expression is kept as it is (nothing is normalized again), a UMAP is made from
    the subset's most variable genes (abc_tk_own.make_umap), and the subset is saved in a project folder of its own (abc_tk_own.write_subset_project)."""
    import abc_tk_own as own_module
    if v.own.kind == "scaled":
        messagebox.showwarning("Analyze subset", "The values of this dataset were scaled (some are negative), so a UMAP of a subset cannot be made from them.",
                               parent=v.root)
        return
    host = v.host
    saved = (own_module.read_project(folder) is not None and (folder / own_module.SUBSET_EXPRESSION_FILE).exists()
             and (folder / "umap_coords.csv").exists())
    if saved:
        sub = host["run_with_loading_window"](f"Opening the saved subset {folder.name}.", lambda: own_module.load_own_subset(folder), show_after_s=0)
    else:
        rows = np.flatnonzero(mask)

        def work():
            part = v.own.subset(rows)
            part.umap = own_module.make_umap(part, n_components)
            own_module.write_subset_project(folder, part, v.run_folder, level, names, n_components)
            return part

        sub = host["run_with_loading_window"](
            f"Analyzing {len(rows):,} cells of {len(names)} {v.level_names[level].lower()} group(s): the most variable genes, PCA ({n_components} dimensions), neighbors and UMAP "
            f"(a minute or two). Saved in {folder}.", work, show_after_s=0)
    core.record_run_opened(folder)
    from abc_tk_umap import UmapWindow
    title = f"{Path(v.run_folder).name} / {folder.name} \u2014 {sub.title} (subset)"
    win = UmapWindow(sub.viewer_adata(), v.abc_cache, None, None, [], folder, None, title, parent=v.root, host=host, params=None, own=sub)
    v.subset_windows.append(win)
    try:
        win.root.state("zoomed")
    except Exception:
        pass
    win.root.lift()


def _embed_and_save(v, mask, folder, level, names, n_components):
    sub = embed_subset(v, mask, n_components)
    save_subset(sub, folder, v.run_folder, level, names, n_components)
    return sub
