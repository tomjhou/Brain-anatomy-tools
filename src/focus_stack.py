"""Extended focus (EFI) from the planes of a z-stack, with numpy only: each pixel is taken from the plane where it is sharpest, the plane
with the most fine detail around it (the mean of |Laplacian| over a small square). Channels are done one at a time (each its own
choice of plane), as the caller asks for them.

    focus_stack(planes)                  the planes (2-D arrays of one shape, e.g. memory maps): the extended-focus picture
    BlockFocusReader(readers, w, h)      rows of it at a step, as the viewer's row readers give them (see xlef_viewer._RowReader)
"""
import collections
import threading

import numpy as np

WINDOW = 9           # the side of the square the detail is averaged over, in pixels of the picture being made
CHUNK = 1024         # rows done at a time by focus_stack (with WINDOW rows above and below), so that memory stays small


def sharpness(a, window=WINDOW):
    """The local amount of fine detail of `a`: |Laplacian| averaged over a window x window square, same shape (edges repeated)."""
    a = np.asarray(a, np.float32)
    p = np.pad(a, 1, mode="edge")
    lap = np.abs(4 * a - p[:-2, 1:-1] - p[2:, 1:-1] - p[1:-1, :-2] - p[1:-1, 2:])
    h = window // 2
    q = np.pad(lap, h, mode="edge")
    s = np.zeros((q.shape[0] + 1, q.shape[1] + 1), np.float64)      # (an integral image: each window's sum from 4 corners)
    np.cumsum(np.cumsum(q, 0), 1, out=s[1:, 1:])
    rows, cols = a.shape
    k = 2 * h + 1
    return (s[k:k + rows, k:k + cols] - s[:rows, k:k + cols] - s[k:k + rows, :cols] + s[:rows, :cols]).astype(np.float32)


def focus_stack(planes, window=WINDOW, chunk=CHUNK):
    """The extended-focus picture of `planes` (a list of 2-D arrays of one shape and type): each pixel from the plane that is sharpest
    there (see sharpness). Done CHUNK rows at a time, each with `window` rows of the planes above and below it (so the choice at a
    chunk's edge is the same as if it were done whole)."""
    planes = list(planes)
    if len(planes) == 1:
        return np.array(planes[0])
    rows = planes[0].shape[0]
    out = np.empty(planes[0].shape, planes[0].dtype)
    for r0 in range(0, rows, chunk):
        r1 = min(rows, r0 + chunk)
        a0, a1 = max(0, r0 - window), min(rows, r1 + window)
        best = best_s = None
        for p in planes:
            s = sharpness(p[a0:a1], window)[r0 - a0:r1 - a0]
            v = np.asarray(p[r0:r1])
            if best is None:
                best, best_s = np.array(v), s
            else:
                better = s > best_s
                best[better] = v[better]
                best_s[better] = s[better]
        out[r0:r1] = best
    return out


class BlockFocusReader:
    """Rows of the extended-focus picture of several planes at a step, as xlef_viewer._RowReader gives them: row(r, x0, n, step) is
    every step-th pixel of row r, n of them from column x0 (r and x0 multiples of step). The sharpness needs the rows around a pixel,
    so it is worked out a block of BLOCK rows at a time (with WINDOW rows and columns more around it, where the image has them), from
    `readers` (one per plane, each with that row()), and the block kept for its other rows; several threads may ask at once (a block
    asked for while it is being made waits for it)."""

    BLOCK = 64
    KEEP = 48

    def __init__(self, readers, width, height):
        self.readers, self.width, self.height = readers, width, height
        self.out_dtype = readers[0].out_dtype
        self.blocks = collections.OrderedDict()
        self.lock = threading.Lock()

    def _block(self, step, b, x0, n):
        key = (step, b, x0, n)
        with self.lock:
            got = self.blocks.get(key)
            mine = got is None
            if mine:
                got = self.blocks[key] = [threading.Event(), None, None]
            else:
                self.blocks.move_to_end(key)
        if mine:
            try:
                got[1] = self._make(step, b, x0, n)
            except Exception as e:
                got[2] = e
                with self.lock:
                    self.blocks.pop(key, None)
            finally:
                got[0].set()
                with self.lock:
                    while len(self.blocks) > self.KEEP:
                        old_key = next(iter(self.blocks))
                        if old_key == key:
                            break
                        self.blocks.pop(old_key)
        else:
            got[0].wait()
        if got[2] is not None:
            raise got[2]
        return got[1]

    def _make(self, step, b, x0, n):
        """Block b (rows b * BLOCK ... of the picture at this step), columns x0 // step ... + n: [rows, n]."""
        last_row = -(-self.height // step) - 1
        last_col = -(-self.width // step)
        R0 = b * self.BLOCK
        R1 = min(last_row + 1, R0 + self.BLOCK)
        A0, A1 = max(0, R0 - WINDOW), min(last_row + 1, R1 + WINDOW)
        X = x0 // step
        C0, C1 = max(0, X - WINDOW), min(last_col, X + n + WINDOW)
        planes = []
        for reader in self.readers:
            a = np.empty((A1 - A0, C1 - C0), self.out_dtype)
            for k in range(A0, A1):
                a[k - A0] = reader.row(k * step, C0 * step, C1 - C0, step)
            planes.append(a)
        whole = focus_stack(planes)
        return whole[R0 - A0:R1 - A0, X - C0:X - C0 + n]

    def row(self, r, x0, n, step):
        R = r // step
        b = R // self.BLOCK
        return self._block(step, b, x0, n)[R - b * self.BLOCK]

    def cost(self, n, step):
        return sum(reader.cost(n + 2 * WINDOW, step) for reader in self.readers)

    def close(self):
        for reader in self.readers:
            reader.close()
        with self.lock:
            self.blocks.clear()
