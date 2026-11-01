"""Open a dataset of your own (an .h5ad file: single-cell or single-nucleus RNA-seq from Seurat or scanpy) in the UMAP viewer.

There is no atlas behind such a dataset and no position of the cells in a section: the viewer shows the UMAP only (no section maps, no imputed
genes, none of the Allen-specific exports). What the viewer needs from the file:
  * a UMAP (obsm 'X_umap'): read from the file, or made here (scanpy) when the file has none, and kept in the cache folder;
  * group columns (obs): every categorical column with more than one value, and the integer columns that are clusters; each is a "level";
  * expression, in one of these forms (see choose_source): raw counts, which are normalized here (the choice is asked), or values that are
    already log-normalized (taken as natural logarithms, as scanpy writes them).
Whatever the form, the expression is kept in memory as one sparse matrix of log2(1 + x) values, so that every gene is read in a moment."""
import hashlib
import os
import re
import tkinter as tk
from pathlib import Path
from tkinter import ttk

import numpy as np
import pandas as pd
import scipy.sparse as sp

LN2 = float(np.log(2.0))
PER_CELL = 1e4                      # "per-cell scaling" brings every cell to this many counts
NORMALIZATIONS = (("cp10k", "Scale each cell to 10,000 counts, then log2(1 + x)  (what Seurat's LogNormalize does)"),
                  ("none", "No scaling: log2(1 + counts)"))
KIND_TEXT = {"counts": "raw counts", "lognorm": "log-normalized values (taken as natural logarithms)",
             "linear": "normalized values (not logarithms)", "scaled": "scaled values (some are negative)"}
MAX_GROUPS = 3000                   # a column with more distinct values than this is not a group column
CLUSTER_WORDS = re.compile(r"cluster|leiden|louvain|ident|type|class|group|label|cell", re.IGNORECASE)
NOT_GROUPS = re.compile(r"ncount|nfeature|^df[._]|pann|percent", re.IGNORECASE)
HIERARCHY = ("class", "subclass", "supertype")          # shown first, in this order, when a column is named so (MMC_Class, ...)


# ---- looking into the file

def classify(sample):
    """What the values in `sample` (a dense piece of the matrix) are: 'counts' (whole numbers), 'lognorm' (small decimals), 'linear' (large
    decimals: normalized, not logarithms) or 'scaled' (some negative: z-scores, residuals)."""
    if sample.size == 0 or float(sample.min()) < 0:
        return "scaled"
    if np.allclose(sample, np.rint(sample)):
        return "counts"
    return "lognorm" if float(sample.max()) <= 25 else "linear"


def sample_rows(source, n=2000):
    part = source[:n]
    return part.toarray() if hasattr(part, "toarray") else np.asarray(part)


def choose_source(a):
    """Where the expression to show is in the AnnData `a` (opened in backed mode): (where, matrix, gene names, kind). Raw counts are preferred (they
    can be normalized the way the user wants), then log-normalized values, then others; among equals, a layer named counts, then raw, then X."""
    found = []
    if "counts" in a.layers:
        found.append(("layers['counts']", a.layers["counts"], a.var_names))
    if a.raw is not None:
        found.append(("raw", a.raw.X, a.raw.var_names))
    found.append(("X", a.X, a.var_names))
    ranked = []
    for where, matrix, names in found:
        ranked.append((where, matrix, names, classify(sample_rows(matrix))))
    for kind in ("counts", "lognorm", "linear", "scaled"):
        for r in ranked:
            if r[3] == kind:
                return r
    return ranked[0]


def find_umap_key(a):
    for key in a.obsm.keys():
        if "umap" in key.lower():
            return key
    return None


def find_group_columns(obs):
    """The columns of `obs` to offer as group types, in the order to show them: the class / subclass / supertype ones first, then the others
    in the file's own order."""
    cols = []
    for c in obs.columns:
        s = obs[c]
        try:
            n = int(s.nunique(dropna=True))
        except TypeError:
            continue
        if n < 2:
            continue
        categorical = str(s.dtype) in ("category", "object", "bool", "string") or str(s.dtype).startswith("str")
        if categorical and n <= MAX_GROUPS:
            cols.append(c)
        elif pd.api.types.is_integer_dtype(s) and n <= 1000 and CLUSTER_WORDS.search(str(c)) and not NOT_GROUPS.search(str(c)):
            cols.append(c)

    def rank(c):
        low = str(c).lower()
        for i, word in enumerate(HIERARCHY):
            if low.endswith(word):
                return (0, i)
        return (1, cols.index(c))

    return sorted(cols, key=rank)


def find_number_columns(obs, groups):
    """The numeric columns of `obs` that are not group types (counts per cell, percentages, scores such as a classification's confidence): each can
    color the cells like a gene, and is in the hover. A column whose values are all the same is left out."""
    cols = []
    for c in obs.columns:
        if c in groups:
            continue
        s = obs[c]
        if str(s.dtype) == "bool" or not pd.api.types.is_numeric_dtype(s):
            continue
        values = pd.to_numeric(s, errors="coerce").to_numpy(float)
        values = values[np.isfinite(values)]
        if len(values) and float(np.ptp(values)) > 0:
            cols.append(str(c))
    return cols


def level_display_name(column):
    text = re.sub(r"_+|(?<!\d)\.|\.(?!\d)", " ", str(column)).strip()          # (underscores, and dots that are not decimal points, become spaces)
    return text[:1].upper() + text[1:]


def describe(path):
    """A look into the file, without reading its expression matrix: what the open dialog says about it."""
    import anndata
    a = anndata.read_h5ad(path, backed="r")
    where, matrix, names, kind = choose_source(a)
    key = find_umap_key(a)
    groups = find_group_columns(a.obs)
    info = {"path": str(path), "n_cells": int(a.n_obs), "n_genes": len(names), "where": where, "kind": kind, "umap_key": key,
            "groups": groups, "numbers": find_number_columns(a.obs, groups)}
    try:
        a.file.close()
    except Exception:
        pass
    return info


# ---- the dataset in memory

class OwnDataset:
    """The cells of the file (their group columns and UMAP) and its expression as one sparse cells x genes matrix of log2(1 + x) values."""

    def __init__(self, path, obs, umap, gene_names, log2, kind, normalization, groups, numbers=()):
        self.path = Path(path)
        self.obs = obs
        self.umap = umap                              # (n, 2) or None while it has to be made
        self.gene_names = np.asarray([str(g) for g in gene_names])
        self.log2 = log2                              # CSC float32, cells x genes
        self.kind = kind
        self.normalization = normalization
        self.groups = list(groups)
        self.numbers = [c for c in numbers if c in obs.columns]            # the numeric columns of the cells (see find_number_columns)
        self._number_names = {c.casefold(): c for c in self.numbers}
        self.level_names = {c: level_display_name(c) for c in self.groups}
        self.level_has_ids = {c: bool(obs[c].astype(str).str.match(r"^\s*\d").all()) for c in self.groups}
        self._index = None

    @property
    def title(self):
        return f"{self.path.name} — {len(self.obs):,} cells"

    def gene_index(self, name):
        if self._index is None:
            self._index = {}
            for i, g in enumerate(self.gene_names):
                self._index.setdefault(g.casefold(), i)
        return self._index.get(str(name).strip().casefold())

    def number_name(self, name):
        """The numeric column called `name` (any case), as the dataset spells it, or None."""
        return self._number_names.get(str(name).strip().casefold())

    def is_number(self, name):
        """Whether `name` is a numeric column (and not a gene)."""
        return self.gene_index(name) is None and self.number_name(name) is not None

    def number_values(self, name):
        """A numeric column's value in every cell."""
        return pd.to_numeric(self.obs[self.number_name(name)], errors="coerce").to_numpy(float)

    def resolve(self, names):
        """Each of `names` as the dataset spells it, a gene or a numeric column: (cols, canonical names, the ones that are neither)."""
        canonical, missing = [], []
        for name in names:
            j = self.gene_index(name)
            if j is not None:
                canonical.append(str(self.gene_names[j]))
            elif self.number_name(name) is not None:
                canonical.append(self.number_name(name))
            else:
                canonical.append(name)
                missing.append(name)
        return [None] * len(names), canonical, missing

    def column(self, gene, scale):
        """A gene's value in every cell: log2(1 + x) for the Log scale, x for Linear (normalized expression, not logs). A numeric column's own
        values, on either scale."""
        j = self.gene_index(gene)
        if j is None and self.number_name(gene) is not None:
            return self.number_values(gene)
        if j is None:
            raise KeyError(gene)
        col = self.log2[:, j].toarray().ravel().astype(float)
        return np.exp2(col) - 1.0 if scale == "linear" else col

    def group_means(self, codes, n_groups, step=400):
        """The mean linear expression of every gene in every group (genes x groups), from the group of each cell (`codes`, -1: none)."""
        ok = np.flatnonzero(codes >= 0)
        indicator = sp.csr_matrix((np.ones(len(ok)), (codes[ok], ok)), shape=(n_groups, len(codes)))
        n = np.maximum(np.bincount(codes[ok], minlength=n_groups).astype(float), 1.0)
        out = np.zeros((self.log2.shape[1], n_groups), np.float32)
        for g0 in range(0, self.log2.shape[1], step):
            block = np.exp2(self.log2[:, g0:g0 + step].toarray()) - 1.0
            out[g0:g0 + step] = ((indicator @ block) / n[:, None]).T
        return out

    def subset(self, rows):
        """The dataset of the cells at `rows` (their positions): the same genes, no UMAP yet. The groups that still tell cells apart are kept."""
        obs = self.obs.iloc[rows].copy()
        keep = [c for c in self.groups if obs[c].nunique() > 1] or list(self.groups)
        return OwnDataset(self.path, obs, None, self.gene_names, self.log2.tocsr()[rows].tocsc(), self.kind, self.normalization, keep, self.numbers)

    def viewer_adata(self):
        """An AnnData for the viewer: the cells with their group columns and UMAP, and the genes' names (the expression is read through this object)."""
        import anndata
        obs = pd.DataFrame({c: pd.Categorical(self.obs[c].astype(str)) for c in self.groups}, index=self.obs.index)
        var = pd.DataFrame({"gene_symbol": self.gene_names}, index=pd.Index([str(i) for i in range(len(self.gene_names))]))
        a = anndata.AnnData(obs=obs, var=var)
        a.obsm["X_umap"] = np.asarray(self.umap, dtype=float)
        return a


def to_log2(matrix, kind, normalization):
    """The expression matrix as a CSC matrix of log2(1 + x) values (float32)."""
    m = matrix.tocsr().astype(np.float32) if sp.issparse(matrix) else sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
    if kind == "counts":
        if normalization == "cp10k":
            total = np.asarray(m.sum(axis=1)).ravel()
            total[total == 0] = 1.0
            m.data *= np.repeat((PER_CELL / total).astype(np.float32), np.diff(m.indptr))
        m.data = np.log2(1.0 + m.data).astype(np.float32)
    elif kind == "lognorm":
        m.data = (m.data / LN2).astype(np.float32)
    elif kind == "linear":
        m.data = np.log2(1.0 + m.data).astype(np.float32)
    return m.tocsc()


def load_own(path, info, normalization):
    """Read the file: the cells' columns, the UMAP if it has one, and the expression. Blocking (it takes about a minute for 25,000 cells): call
    it from a worker thread or behind a loading window."""
    import anndata
    a = anndata.read_h5ad(path, backed="r")
    where, matrix, names, kind = choose_source(a)
    umap = np.asarray(a.obsm[info["umap_key"]], dtype=float)[:, :2] if info["umap_key"] else None
    obs = a.obs.copy()
    log2 = to_log2(matrix[:], kind, normalization)
    gene_names = list(names)
    try:
        a.file.close()
    except Exception:
        pass
    return OwnDataset(path, obs, umap, gene_names, log2, kind, normalization, info["groups"], info.get("numbers", []))


# ---- a UMAP for a file that has none

def umap_cache_path(own, cache_dir):
    st = os.stat(own.path)
    key = f"{own.path}|{st.st_size}|{st.st_mtime_ns}|{own.kind}|{own.normalization}"
    return Path(cache_dir) / f"own_umap_{hashlib.md5(key.encode('utf-8')).hexdigest()[:16]}.npy"


def make_umap(own, n_components=50):
    """A UMAP of the cells (scanpy): the 2,000 most variable genes, scaled, `n_components` principal components (50 unless asked: no more than the
    cells and genes allow), 15 neighbors. A few minutes for 25,000 cells (the first time scanpy compiles its UMAP code, which adds a minute)."""
    import anndata
    import scanpy as sc
    x = own.log2.tocsr().astype(np.float32)
    x.data *= np.float32(LN2)                                   # (scanpy's own convention for log values: natural logarithms)
    a = anndata.AnnData(X=x, obs=pd.DataFrame(index=pd.Index([str(i) for i in range(x.shape[0])])))
    sc.pp.highly_variable_genes(a, n_top_genes=min(2000, x.shape[1] - 1), flavor="seurat")
    a = a[:, a.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=min(int(n_components), a.n_vars - 1, a.n_obs - 1))
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=a.obsm["X_pca"].shape[1])
    sc.tl.umap(a)
    return np.asarray(a.obsm["X_umap"], dtype=float)


def ensure_umap(own, cache_dir):
    """The dataset's UMAP: the file's own, else one made earlier and kept, else one made now (and kept)."""
    if own.umap is not None:
        return own.umap
    path = umap_cache_path(own, cache_dir)
    if path.exists():
        saved = np.load(path)
        if saved.shape == (len(own.obs), 2):
            own.umap = saved
            return saved
    own.umap = make_umap(own)
    try:
        Path(cache_dir).mkdir(parents=True, exist_ok=True)
        np.save(path, own.umap)
    except OSError as e:
        print(f"Could not save the UMAP ({e}); it will be made again next time.", flush=True)
    return own.umap


# ---- the dialog, and the whole of opening a file

def ask_options(info):
    """The open dialog: what was found in the file, and how to normalize counts. Returns the normalization ('cp10k' or 'none'), or None."""
    root = tk.Tk()
    root.title("Open your own dataset")
    root.resizable(False, False)
    result = {"value": None}
    body = ttk.Frame(root, padding=(22, 18, 22, 18))
    body.pack()
    ttk.Label(body, text=Path(info["path"]).name, font=("Segoe UI", 13, "bold")).pack(anchor="w")
    ttk.Label(body, text=f"{info['n_cells']:,} cells, {info['n_genes']:,} genes", foreground="gray30").pack(anchor="w", pady=(2, 10))
    lines = [f"Expression: {KIND_TEXT[info['kind']]}, from {info['where']}.",
             "UMAP: found in the file." if info["umap_key"] else "UMAP: none in the file. It will be made here (a few minutes the first time, then kept).",
             "Groups: " + (", ".join(level_display_name(c) for c in info["groups"]) if info["groups"] else "none found (the UMAP will be one color).")]
    for text in lines:
        ttk.Label(body, text=text, wraplength=560, justify=tk.LEFT).pack(anchor="w", pady=2)
    var = tk.StringVar(root, NORMALIZATIONS[0][0])
    if info["kind"] == "counts":
        box = ttk.LabelFrame(body, text="Normalization of the counts", padding=8)
        box.pack(fill=tk.X, pady=(12, 0))
        for value, text in NORMALIZATIONS:
            ttk.Radiobutton(box, text=text, value=value, variable=var).pack(anchor="w")
    elif info["kind"] == "scaled":
        ttk.Label(body, text="These values were scaled (some are negative), so they are shown as they are: log scale only.", foreground="#a05000",
                  wraplength=560).pack(anchor="w", pady=(10, 0))
    row = ttk.Frame(body)
    row.pack(fill=tk.X, pady=(16, 0))

    def done(ok):
        result["value"] = var.get() if ok else None
        root.destroy()

    ttk.Button(row, text="Cancel", command=lambda: done(False)).pack(side=tk.RIGHT)
    ok_button = ttk.Button(row, text="OK", command=lambda: done(True))
    ok_button.pack(side=tk.RIGHT, padx=(0, 8))
    root.bind("<Return>", lambda e: done(True))
    root.bind("<Escape>", lambda e: done(False))
    root.protocol("WM_DELETE_WINDOW", lambda: done(False))
    root.update_idletasks()
    root.geometry(f"+{(root.winfo_screenwidth() - root.winfo_width()) // 2}+{(root.winfo_screenheight() - root.winfo_height()) // 3}")
    root.attributes("-topmost", True)
    ok_button.focus_force()
    root.mainloop()
    return result["value"]


# ---- the project folder: where a dataset of your own keeps what is made from it

PROJECT_FILE = "own_dataset.json"            # in a project folder: says where the dataset is and how it was normalized
STATISTICS_FILE = "statistics.txt"


def same_path(a, b):
    norm = lambda p: os.path.normcase(os.path.abspath(str(p)))
    return norm(a) == norm(b)


def read_project(folder):
    """The project file of `folder`, or None if it is not the project folder of a dataset of your own."""
    import json
    try:
        with open(Path(folder) / PROJECT_FILE, encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) and data.get("source") else None
    except (OSError, ValueError):
        return None


def project_folder_for(out_folder, path):
    """The project folder for the dataset at `path`: one made for it before, else a new name in `out_folder` from the file's (with a few letters of
    the path's hash if another dataset of the same name is already there). The folder is not made here."""
    out = Path(out_folder)
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", Path(path).stem).strip("._") or "dataset"
    base = out / f"own_{stem}"
    for candidate in [base] + [out / f"own_{stem}_{hashlib.md5(str(path).encode('utf-8')).hexdigest()[:6]}"]:
        project = read_project(candidate)
        if project is None and not candidate.exists():
            return candidate
        if project is not None and same_path(project["source"], path):
            return candidate
    return out / f"own_{stem}_{hashlib.md5((str(path) + str(os.stat(path).st_mtime_ns)).encode('utf-8')).hexdigest()[:6]}"


def write_project(folder, own, info):
    """Make the project folder (or bring it up to date): the project file, the UMAP with the cells' group columns, and a statistics file."""
    import json
    import time
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / PROJECT_FILE, "w", encoding="utf-8") as f:
        json.dump({"source": str(own.path), "normalization": own.normalization, "kind": own.kind, "expression_from": info["where"],
                   "n_cells": int(len(own.obs)), "n_genes": int(len(own.gene_names)), "umap_from": "the file" if info["umap_key"] else "made by the browser",
                   "groups": own.groups, "numbers": own.numbers, "written": time.strftime("%Y-%m-%d %H:%M")}, f, indent=2)
    table = pd.DataFrame({"cell_label": own.obs.index.astype(str), "UMAP1": own.umap[:, 0], "UMAP2": own.umap[:, 1]})
    for c in own.groups:
        table[c] = own.obs[c].astype(str).to_numpy()
    for c in own.numbers:
        table[c] = own.number_values(c)
    table.to_csv(folder / "umap_coords.csv", index=False)
    norm = dict(NORMALIZATIONS).get(own.normalization, own.normalization) if own.kind == "counts" else "not applied: " + KIND_TEXT[own.kind]
    lines = ["Statistics of this project", "==========================", f"Project folder:  {folder.name}",
             f"Written:         {time.strftime('%Y-%m-%d %H:%M')}", "", "Dataset", "-------", f"File:            {own.path}",
             f"Cells:           {len(own.obs):,}", f"Genes:           {len(own.gene_names):,}",
             f"Expression:      {KIND_TEXT[own.kind]}, from {info['where']}", f"Normalization:   {norm}",
             f"UMAP:            {'from the file' if info['umap_key'] else 'made by the browser (most variable genes, PCA, neighbors, UMAP)'}", "", "Groups", "------"]
    for c in own.groups:
        lines.append(f"{own.level_names[c]:<28} {own.obs[c].nunique():,} values")
    (folder / STATISTICS_FILE).write_text("\n".join(lines) + "\n", encoding="utf-8")


SUBSET_EXPRESSION_FILE = "expression_log2.npz"      # in the folder of a subset: its cells' expression, so that it opens without the dataset's file
SUBSET_GENES_FILE = "genes.txt"


def write_subset_project(folder, sub, parent_folder, level, names, n_components=50):
    """Make the project folder of a subset (Analyze subset in the viewer): the project file (it names the dataset's file as the parent's does, and says
    which groups the subset is), the cells' UMAP and groups, their expression (log2 values, a sparse file) and the genes' names, and statistics."""
    import json
    import time
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    parent = read_project(parent_folder) or {}
    sp.save_npz(folder / SUBSET_EXPRESSION_FILE, sub.log2.tocsc())
    (folder / SUBSET_GENES_FILE).write_text("\n".join(sub.gene_names) + "\n", encoding="utf-8")
    table = pd.DataFrame({"cell_label": sub.obs.index.astype(str), "UMAP1": sub.umap[:, 0], "UMAP2": sub.umap[:, 1]})
    for c in sub.groups:
        table[c] = sub.obs[c].astype(str).to_numpy()
    for c in sub.numbers:
        table[c] = sub.number_values(c)
    table.to_csv(folder / "umap_coords.csv", index=False)
    when = time.strftime("%Y-%m-%d %H:%M")
    with open(folder / PROJECT_FILE, "w", encoding="utf-8") as f:
        json.dump({"source": str(sub.path), "normalization": sub.normalization, "kind": sub.kind, "expression_from": parent.get("expression_from", ""),
                   "n_cells": int(len(sub.obs)), "n_genes": int(len(sub.gene_names)), "umap_from": "made by the browser (this subset)",
                   "groups": sub.groups, "numbers": sub.numbers, "written": when,
                   "subset": {"parent": str(parent_folder), "level": level, "groups": list(names), "cells": int(len(sub.obs)),
                              "pca_components": int(n_components)}}, f, indent=2)
    lines = ["Statistics of this subset", "=========================", f"Project folder:  {folder.name}", f"Written:         {when}", "",
             f"Made from:       {Path(parent_folder).name}", f"Groups ({level}): " + ", ".join(names), f"Cells:           {len(sub.obs):,}",
             f"Genes:           {len(sub.gene_names):,}", f"File:            {sub.path}", f"UMAP:            the most variable genes, PCA ({n_components} components), neighbors, UMAP, on these cells alone"]
    (folder / STATISTICS_FILE).write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_own_subset(folder):
    """The dataset of a subset folder (write_subset_project), from its own files."""
    saved = read_project(folder)
    folder = Path(folder)
    table = pd.read_csv(folder / "umap_coords.csv", dtype={"cell_label": str})
    groups = [c for c in saved.get("groups", []) if c in table.columns]
    numbers = [c for c in saved.get("numbers", []) if c in table.columns]
    obs = table.set_index("cell_label")[groups + numbers]
    obs.index.name = None
    genes = (folder / SUBSET_GENES_FILE).read_text(encoding="utf-8").split("\n")
    genes = genes[:-1] if genes and genes[-1] == "" else genes
    log2 = sp.load_npz(folder / SUBSET_EXPRESSION_FILE).tocsc().astype(np.float32)
    return OwnDataset(saved["source"], obs, table[["UMAP1", "UMAP2"]].to_numpy(dtype=float), genes, log2, saved.get("kind", "counts"),
                      saved.get("normalization", NORMALIZATIONS[0][0]), groups, numbers)


def saved_umap(folder, own):
    """The UMAP kept in the project folder, if its cells are those of the dataset, else None."""
    try:
        table = pd.read_csv(Path(folder) / "umap_coords.csv", usecols=["cell_label", "UMAP1", "UMAP2"], dtype={"cell_label": str})
    except (OSError, ValueError):
        return None
    if len(table) != len(own.obs) or not (table["cell_label"].to_numpy() == own.obs.index.astype(str).to_numpy()).all():
        return None
    return table[["UMAP1", "UMAP2"]].to_numpy(dtype=float)


def locate_dialog(old_path):
    """The dataset of a project has moved: ask where it is now. Returns the new path, or None."""
    from tkinter import filedialog, messagebox
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    messagebox.showinfo("Dataset not found", f"The dataset of this project is not where it was:\n\n{old_path}\n\nChoose where it is now.", parent=root)
    chosen = filedialog.askopenfilename(parent=root, title=f"Where is {Path(old_path).name}?", initialdir=str(Path(old_path).parent) if Path(old_path).parent.exists() else None,
                                        filetypes=[("AnnData files", "*.h5ad"), ("All files", "*.*")])
    root.destroy()
    return Path(chosen) if chosen else None


def prepare_own_dataset(path, run_blocking, cache_dir, out_folder):
    """Everything between choosing the file and opening the viewer: the project folder for the file (one made earlier is reused, with the
    normalization it was made with; else a look into the file, the dialog and a new folder), the load, and the UMAP if there is none.
    `run_blocking(message, work)` runs `work()` behind a loading window and returns its result. Returns (the OwnDataset, the project folder), or None."""
    path = Path(path)
    folder = project_folder_for(out_folder, path)
    saved = read_project(folder)
    info = run_blocking(f"Looking into {path.name}...", lambda: describe(path))
    if saved is not None and saved.get("normalization") in dict(NORMALIZATIONS):
        normalization = saved["normalization"]
        print(f"This dataset has a project folder already ({folder.name}): using it, with the normalization it was made with.", flush=True)
    else:
        normalization = ask_options(info)
        if normalization is None:
            return None
    return open_with(path, info, normalization, folder, run_blocking, cache_dir)


def open_with(path, info, normalization, folder, run_blocking, cache_dir):
    own = run_blocking(f"Loading {Path(path).name} ({info['n_cells']:,} cells, {info['n_genes']:,} genes).\nAbout a minute from a network drive.",
                       lambda: load_own(path, info, normalization))
    if own.umap is None:
        own.umap = saved_umap(folder, own)
    if own.umap is None:
        if own.kind == "scaled":
            raise ValueError("this file has no UMAP, and its values are scaled, so a UMAP cannot be made from them")
        run_blocking("This file has no UMAP: making one (the most variable genes, PCA, neighbors, UMAP).\nA few minutes the first time; it is kept in the project folder.",
                     lambda: ensure_umap(own, cache_dir))
    write_project(folder, own, info)
    return own, folder


def reopen_own_project(folder, run_blocking, cache_dir):
    """Open a project folder made by prepare_own_dataset: no questions, but where the dataset is if it has moved. Returns (OwnDataset, folder) or None."""
    folder = Path(folder)
    saved = read_project(folder)
    if saved.get("subset") and (folder / SUBSET_EXPRESSION_FILE).exists():
        return run_blocking(f"Opening {folder.name}.", lambda: load_own_subset(folder)), folder          # (a subset keeps its own expression: no file to find)
    path = Path(saved["source"])
    if not path.exists():
        path = locate_dialog(path)
        if path is None:
            return None
    info = run_blocking(f"Looking into {path.name}...", lambda: describe(path))
    return open_with(path, info, saved.get("normalization", NORMALIZATIONS[0][0]), folder, run_blocking, cache_dir)
