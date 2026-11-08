"""The project file of an Allen atlas project (project.json), and what a project without a UMAP needs.

Every project the browser makes has a project.json in its folder: what its cells are (cell type, sections, ROIs), the atlas release it
was made from and whether it has a UMAP. A project without a UMAP ("no-UMAP") is only that file (and the ROI file, if it has ROIs): its
cells are picked out of the atlas's tables each time it is opened, and the gene expression of its cells is read from the atlas's own
file when a gene is wanted (AtlasCounts), so nothing of the size of a UMAP run's processed file is kept.

Projects made before this file existed have none; they are recognized by umap_coords.csv, as they always were (read_project is None).
This module needs abc_core only, so that it can be used (and tested) without starting the browser."""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import abc_core as core

PROJECT_FILE = "project.json"
KIND = "abc_atlas"                 # (a dataset of the user's own has own_dataset.json: see abc_tk_own)
SCHEMA = 1
DATASET = "MERFISH-C57BL6J-638850"


# ---------------------------------------------------------------------------------------------------------------------- the file

def read_project(folder):
    """The project file of `folder` as a dict, or None if it has none (a project made before there was one, or not a project)."""
    try:
        with open(Path(folder) / PROJECT_FILE, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) and data.get("kind") == KIND else None


def has_umap(folder):
    """Whether the project in `folder` has a UMAP: False only if its project file says so (all earlier projects have one)."""
    info = read_project(folder)
    return True if info is None else bool(info.get("has_umap", True))


def atlas_release(abc_cache):
    """A short text naming the atlas release the cache is on, or None."""
    try:
        return str(abc_cache.current_manifest)
    except Exception:
        return None


def _clean_rois(rois):
    out = []
    for r in rois or []:
        out.append({"section": str(r["section"]), **{k: float(r[k]) for k in ("x_min", "x_max", "y_min", "y_max")}})
    return out


def write_project(folder, cell_type, sections, rois, label="", has_umap=True, cap=None, n_cells=None, abc_cache=None, frame=None):
    """Write `folder`/project.json and return what it holds. `sections`: the whole sections picked (None: every section); `rois`: the list of
    {'section', 'x_min', 'x_max', 'y_min', 'y_max'} (ROIs, when there are any, are the cells; the sections are then only for show); `cap`: the
    most cells the UMAP was made of (a UMAP project only); `frame`: the coordinate frame the ROIs were drawn on. The time the project was first
    made stays as it was."""
    folder = Path(folder)
    previous = read_project(folder) or {}
    info = {
        "kind": KIND,
        "schema": SCHEMA,
        "created": previous.get("created") or time.strftime("%Y-%m-%d %H:%M:%S"),
        "has_umap": bool(has_umap),
        "label": label or "",
        "atlas": {"dataset": DATASET, "release": atlas_release(abc_cache) if abc_cache is not None else previous.get("atlas", {}).get("release")},
        "selection": {
            "cell_type": cell_type,
            "sections": None if sections is None else [str(s) for s in sections],
            "rois": _clean_rois(rois),
            "frame": frame if frame is not None else previous.get("selection", {}).get("frame"),
        },
        "umap": {"cap": None if cap is None else int(cap), "seed": 0} if has_umap else None,
        "n_cells": None if n_cells is None else int(n_cells),
    }
    folder.mkdir(parents=True, exist_ok=True)
    tmp = folder / (PROJECT_FILE + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2)
    tmp.replace(folder / PROJECT_FILE)
    return info


def release_differs(info, abc_cache):
    """True if the atlas release now in the cache is not the one the project was made from (both known)."""
    saved = (info.get("atlas") or {}).get("release")
    now = atlas_release(abc_cache)
    return bool(saved and now and saved != now)


# ---------------------------------------------------------------------------------------------------------------------- the cells

def select_cells(backed, abc_cache, cell_type, sections, rois, section_series=None):
    """Which cells of the atlas (`backed`) a project's selection holds: a boolean array over the atlas's cells. The same cells that
    filter_by_cell_type, filter_by_sections and filter_by_rois pick out, but from the annotation and coordinate tables alone, so no expression is
    read. `cell_type`: 'Neurons', 'NonNeurons' or 'All'; ROIs, if there are any, replace the whole-section choice."""
    annotations = core.load_annotations_aligned(backed.obs.index, abc_cache)
    if section_series is not None:
        section = section_series.to_numpy()
    else:
        section = annotations[core.SECTION_COL].to_numpy()
    if rois:
        spatial = core.load_section_spatial_coords(backed, abc_cache)
        xs = pd.to_numeric(spatial["x"], errors="coerce").to_numpy()
        ys = pd.to_numeric(spatial["y"], errors="coerce").to_numpy()
        keep = np.zeros(len(xs), dtype=bool)
        for roi in rois:
            keep |= (section == roi["section"]) & (xs >= roi["x_min"]) & (xs <= roi["x_max"]) & (ys >= roi["y_min"]) & (ys <= roi["y_max"])
    elif sections is not None:
        keep = pd.Series(section).isin(sections).to_numpy()
    else:
        keep = np.ones(len(section), dtype=bool)
    if cell_type == "All":
        return keep
    # (the class number is read from each distinct class name once, about 34 of them: matching all 4 million cells took seconds)
    codes, names = pd.factorize(annotations["class"].to_numpy())
    ids = pd.to_numeric(pd.Series(names).astype(str).str.strip().str.extract(r"^(\d+)")[0], errors="coerce").to_numpy(dtype=float)
    class_id = np.where(codes >= 0, ids[np.maximum(codes, 0)], np.nan) if len(ids) else np.full(len(codes), np.nan)
    valid = ~np.isnan(class_id)
    non_neuron = valid & np.isin(class_id, list(core.NON_NEURON_CLASS_IDS))
    kind = non_neuron if cell_type == "NonNeurons" else (valid & ~non_neuron)
    return keep & kind


def selection_of(info):
    """(cell_type, sections, rois) of a project file's dict."""
    s = info.get("selection") or {}
    return s.get("cell_type") or "Neurons", s.get("sections"), s.get("rois") or []


def build_adata(backed, abc_cache, mask):
    """The cells of `mask` as an AnnData with their annotations (class, subclass, supertype, cluster, neurotransmitter, section) and the atlas's
    genes, but no expression and no embedding (adata.uns['no_umap'] is True): what the viewer opens a project without a UMAP on. The gene
    values come from AtlasCounts."""
    import anndata
    annotations = core.load_annotations_aligned(backed.obs.index, abc_cache)
    obs = annotations.loc[np.asarray(mask)].copy()
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = pd.Categorical(obs[c])
    adata = anndata.AnnData(obs=obs, var=backed.var.copy())
    adata.uns["no_umap"] = True
    return adata


# ---------------------------------------------------------------------------------------------------------------------- the genes

class AtlasCounts:
    """The raw counts of a project's cells, read from the atlas's file (`backed`: its backed AnnData, a 7.6 GB sparse matrix of 4 million cells) when
    wanted: nothing of them is kept. `positions`: the row of each project cell in that file."""

    def __init__(self, backed, positions):
        self.backed = backed
        self.positions = np.asarray(positions, dtype=np.int64)
        self.n_obs = len(self.positions)

    def column(self, col):
        """Gene number `col`'s count in every project cell. The gene's whole column is read from the file once (a pass over it; seconds to a minute) and
        kept in memory, so every later look at the gene is instant."""
        return core._get_full_gene_column(self.backed, col)[self.positions]

    def rows(self, rows):
        """The counts of project cells `rows` (indices into the project, in the order given) as a scipy CSR matrix of cells x genes."""
        import scipy.sparse as sp
        rows = np.asarray(rows)
        pos = self.positions[rows]
        order = np.argsort(pos, kind="stable")
        sub = self.backed[pos[order]].to_memory()
        x = sub.X if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
        inverse = np.empty(len(order), dtype=np.int64)
        inverse[order] = np.arange(len(order))
        return x[inverse].tocsr()

    def group_means(self, codes, n_groups, chunk_rows=250_000, progress=None, cancelled=None):
        """The mean raw count of every gene in every group, as a genes x groups float32 array (as abc_tk_analyze.group_means_panel does it for a run's own
        counts). `codes`: each project cell's group (-1: none). The file is read in blocks of rows that hold project cells, once. `progress(fraction)`
        is called as it goes; returns None if `cancelled()` says to stop."""
        import scipy.sparse as sp
        codes = np.asarray(codes)
        order = np.argsort(self.positions, kind="stable")
        pos = self.positions[order]
        grp = codes[order]
        n_genes = self.backed.n_vars
        sums = np.zeros((n_groups, n_genes), dtype=np.float64)
        counts = np.bincount(grp[grp >= 0], minlength=n_groups).astype(np.float64)
        total = self.backed.n_obs
        start = 0
        while start < total:
            stop = min(start + chunk_rows, total)
            lo, hi = np.searchsorted(pos, [start, stop])
            if hi > lo:
                if cancelled is not None and cancelled():
                    return None
                block = self.backed.X[start:stop]
                block = block if sp.issparse(block) else sp.csr_matrix(block)
                g = grp[lo:hi]
                ok = np.flatnonzero(g >= 0)
                if len(ok):
                    picked = block[pos[lo:hi][ok] - start]
                    indicator = sp.csr_matrix((np.ones(len(ok)), (g[ok], np.arange(len(ok)))), shape=(n_groups, len(ok)))
                    part = indicator @ picked
                    sums += part.toarray() if hasattr(part, "toarray") else np.asarray(part)
            start = stop
            if progress is not None:
                progress(min(start / total, 1.0))
        return (sums / np.maximum(counts, 1.0)[:, None]).T.astype(np.float32)
