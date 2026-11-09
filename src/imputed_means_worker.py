"""The read of the imputed gene file that Analyze's Selectivity needs: the mean linear expression of every gene in every group.

This module imports only numpy, scipy and h5py, so that it can be run as a helper process (a script of its own: the main program cannot be started again
by a process it spawns). Reading the 50 GB file is decompression, which does not run in parallel in threads (HDF5 holds one lock), but does in separate
processes: abc_tk_analyze starts several of these, each doing every Nth block of genes, and they write their sums into one shared file.

Run as:  python imputed_means_worker.py JOB.npz SUMS.npy WORKER_INDEX N_WORKERS PROGRESS_FILE
"""
import sys
import traceback

import numpy as np


def prepare(positions, codes_list, n_groups_list):
    """What the read needs from the run's cells. `positions`: each run cell's row in the file (-1: not in it); `codes_list`, `n_groups_list`: for each group type,
    the group of each run cell (-1 for none) and how many groups there are. Returns (the cells in the file, their rows, the indicator of groups by cells, where each
    type's groups start in it, the total of groups)."""
    import scipy.sparse as sp
    valid = np.flatnonzero(positions >= 0)
    rows_in_file = np.asarray(positions)[valid]
    offsets, total = [], 0
    for n_groups in n_groups_list:
        offsets.append(total)
        total += int(n_groups)
    r_idx, c_idx = [], []
    for codes, offset in zip(codes_list, offsets):           # one indicator for all the types, the groups of each under those of the one before
        c = np.asarray(codes)[valid]
        ok = np.flatnonzero(c >= 0)
        r_idx.append(c[ok].astype(np.int64) + offset)
        c_idx.append(ok)
    indicator = sp.csc_matrix((np.ones(sum(len(r) for r in r_idx), np.float32), (np.concatenate(r_idx), np.concatenate(c_idx))), shape=(total, len(valid)))
    return valid, rows_in_file, indicator, offsets, total


def read_blocks(path, prepared, gene_blocks, sink, progress, cancelled):
    """Sum the linear expression (2**v - 1 of the stored log2 values) of the genes of each `gene_blocks` entry (g0, g1) over each group: sink(g0, g1, sums), sums
    being groups x genes. Only the blocks of cells that hold run cells are read. progress(done, of); False if cancelled."""
    import h5py
    valid, rows_in_file, indicator, offsets, total = prepared
    with h5py.File(path, "r") as f:
        x = f["X"]
        size = x.chunks[0]
        block_of = rows_in_file // size
        blocks = []
        for b in np.unique(block_of):
            pick = np.flatnonzero(block_of == b)
            blocks.append((int(b) * size, min((int(b) + 1) * size, x.shape[0]), rows_in_file[pick] - int(b) * size, indicator[:, pick]))
        for n, (g0, g1) in enumerate(gene_blocks):
            if cancelled():
                return False
            acc = np.zeros((total, g1 - g0), np.float32)
            for r0, r1, local, ind in blocks:
                linear = np.exp2(x[r0:r1, g0:g1][local].astype(np.float32)) - 1.0
                acc += ind @ linear
            sink(g0, g1, acc)
            progress(n + 1, len(gene_blocks))
    return True


def gene_block_list(n_genes, block_genes):
    return [(g0, min(g0 + block_genes, n_genes)) for g0 in range(0, n_genes, block_genes)]


def main(argv):
    job_path, sums_path, index, n_workers, progress_path = argv[0], argv[1], int(argv[2]), int(argv[3]), argv[4]
    with np.load(job_path, allow_pickle=False) as z:
        path = str(z["path"])
        n_groups_list = [int(n) for n in z["n_groups"]]
        positions = z["positions"]
        codes_list = [z[f"codes_{i}"] for i in range(len(n_groups_list))]
        block_genes = int(z["block_genes"])
    sums = np.load(sums_path, mmap_mode="r+")
    mine = gene_block_list(sums.shape[1], block_genes)[index::n_workers]
    prepared = prepare(positions, codes_list, n_groups_list)

    def sink(g0, g1, acc):
        sums[:, g0:g1] = acc

    def progress(done, of):
        with open(progress_path, "w") as f:
            f.write(str(done))

    ok = read_blocks(path, prepared, mine, sink, progress, lambda: False)
    sums.flush()
    return 0 if ok else 1


if __name__ == "__main__":
    try:
        code = main(sys.argv[1:])
    except Exception:
        traceback.print_exc()
        code = 2
    sys.exit(code)
