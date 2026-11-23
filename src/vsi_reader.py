"""Olympus cellSens / VS200 .vsi images for the XLEF viewer, read with slideio (pip install slideio).

A .vsi is a small file of metadata; its pixels are in .ets files in the folder "_<name>_" next to it (stack1/, stack10000/, ...), as
JPEG 2000 tiles of 512 x 512 with a pyramid of smaller copies. slideio gives each scanned region of the slide as a scene (several
sections in one .vsi: a scene each), with its channels, z planes and pixel size; this module adds what the viewer needs on top:
  - describe(path): each scene's size, channels (their names, and the display limits set in cellSens, from the .vsi's metadata), z
    planes, pixel size, and the size of its .ets file;
  - overview(...): a scene's channel at every step-th pixel (from the pyramid level that has the detail, as slideio picks it);
  - BlockReader: rows of a channel at a step, as the viewer's row readers give them (see xlef_viewer._RowReader), read in blocks of
    BLOCK x BLOCK pixels at that step, each once.
Steps should be powers of 2 (see pow2_below / pow2_above): slideio then reads them from a pyramid level as they are, in place to a hundredth
of a pixel. Any other step is scaled from the next coarser level, and comes out shifted by up to almost a pixel at that step (measured
2026-11-22: 0.9 px at step 15, i.e. 13 full-resolution pixels; 0.2 px at step 3), which showed as the picture jumping when zoomed. A part
of the image is read as a whole number of steps across, for the same reason (scaled to fit, a width that isn't was stretched).
One slideio Slide per file, shared by all the threads that read from it (checked 2026-11-22: 8 threads reading one scene at once got
the same pixels as reading one block at a time, and in a third of the time)."""
import json
import math
import os
import re
import threading

import numpy as np

try:
    import slideio
except ImportError:              # (the viewer still runs: a .vsi then says what to install)
    slideio = None

MISSING = ("Reading Olympus .vsi files needs the Python package slideio, which isn't installed.\n\n"
           "Starting the viewer from its desktop or Start menu shortcut installs it (if the computer is online); "
           "or install it yourself with:  python -m pip install slideio\nThen restart the viewer.")

BLOCK = 512                      # the side of the blocks BlockReader reads, in pixels at its step (the .ets tiles are 512 x 512)

# The color map for a channel, from its name (cellSens names them by filter or dye): the first pattern that matches.
_LUTS = ((r"dapi|hoechst|\b405\b|blue", "Blue"),
         (r"fitc|gfp|\b488\b|cy2|alexa ?488|green", "Green"),
         (r"cy3\b|tritc|\b5[3-9]\d\b|rfp|mcherry|tdtomato|red", "Red"),
         (r"cy5\b|\b6[3-9]\d\b|far.?red", "Magenta"),
         (r"cy7\b|\b7[3-9]\d\b|\b8\d\d\b", "Cyan"))

_SLIDES = {}
_SLIDES_LOCK = threading.Lock()


def available():
    return slideio is not None


def _slide(path):
    """The file's slideio Slide, opened once (and shared: see the module's note)."""
    if slideio is None:
        raise RuntimeError(MISSING)
    key = os.path.normcase(os.path.abspath(path))
    with _SLIDES_LOCK:
        s = _SLIDES.get(key)
        if s is None:
            s = _SLIDES[key] = slideio.open_slide(path, "VSI")
        return s


def _scene(path, index):
    return _slide(path).get_scene(index)


def lut_for(name):
    low = (name or "").lower()
    return next((lut for pattern, lut in _LUTS if re.search(pattern, low)), "Gray")


# --------------------------------------------------------------------------------------------------- metadata (the .vsi's tag tree)

def _volumes(raw):
    """The 'Multidimensional image volume' nodes of the .vsi's metadata (slideio's JSON of its tag tree), in the file's order: one per
    stack of the slide (the label, the overview, focus maps, and each scanned region)."""
    try:
        root = json.loads(raw)
    except (TypeError, ValueError):
        return []
    out = []

    def walk(node):
        if node.get("name") == "Multidimensional image volume":
            out.append(node)
        v = node.get("value")
        if isinstance(v, list):
            for k in v:
                walk(k)
    walk(root)
    return out


def _leaves(node):
    v = node.get("value")
    if isinstance(v, list):
        for k in v:
            yield from _leaves(k)
    else:
        yield node.get("name"), v


def _volume_info(volume):
    """(stack name, [(black, white) per channel, in order]) of a volume: cellSens's display limits, "(0.000000, 16812.000000)"."""
    name, limits = None, []
    for key, value in _leaves(volume):
        if key == "Stack name" and name is None:
            name = value
        elif key == "Display limits":
            m = re.match(r"\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*\)", str(value))
            if m:
                limits.append((float(m.group(1)), float(m.group(2))))
    return name, limits


def _ets_file(path, volume_index, n_volumes, n_image_stacks):
    """The .ets file of the stack that is volume `volume_index` of the metadata, or None. cellSens numbers the stacks' folders
    stack1 (the label: the first volume), then stack10000, 10001, ... for the volumes after it, in order; a .vsi of one image (a
    projection, say) has its one stack in stack1."""
    folder = os.path.join(os.path.dirname(path), "_" + os.path.splitext(os.path.basename(path))[0] + "_")
    names = [f"stack{10000 + volume_index - 1}"] if volume_index > 0 else []
    if n_image_stacks == 1:
        names.append("stack1")
    for name in names:
        d = os.path.join(folder, name)
        try:
            ets = [os.path.join(d, f) for f in os.listdir(d) if f.lower().endswith(".ets")]
        except OSError:
            continue
        if ets:
            return max(ets, key=os.path.getsize)
    return None


def describe(path):
    """Each scene (scanned region) of the .vsi: a dict of name, index, width, height, nz, pixel_um, z_um, dtype, channels [(name,
    max value, black, white)], ets (its .ets file, or None) and bytes (that file's size, or None). Reads only the metadata."""
    s = _slide(path)
    volumes = _volumes(s.raw_metadata)
    infos = [_volume_info(v) for v in volumes]
    scenes = []
    for i in range(s.num_scenes):
        sc = s.get_scene(i)
        _, _, w, h = sc.rect
        n = sc.num_channels
        bits = [sc.get_channel_significant_bits(c) or 16 for c in range(n)]
        names = [sc.get_channel_name(c) or f"Channel {c + 1}" for c in range(n)]
        scenes.append({"name": sc.name or f"Scene {i + 1}", "index": i, "width": w, "height": h, "nz": max(1, sc.num_z_slices),
                       "pixel_um": (sc.resolution[0] * 1e6) if sc.resolution and sc.resolution[0] else None,
                       "z_um": (sc.z_resolution * 1e6) if sc.z_resolution else None,
                       "dtype": np.dtype(str(sc.get_channel_data_type(0))), "names": names, "bits": bits, "volume": None})
    # Each scene's volume in the metadata: the next one of the same stack name that has a display limit for each of its channels.
    used = set()
    for sc in scenes:
        for k, (name, limits) in enumerate(infos):
            if k not in used and name == sc["name"] and len(limits) >= len(sc["names"]):
                used.add(k)
                sc["volume"] = k
                break
    image_stacks = sum(1 for name, limits in infos if limits)
    for sc in scenes:
        k = sc.pop("volume")
        limits = infos[k][1] if k is not None else []
        sc["channels"] = []
        for c, (name, b) in enumerate(zip(sc.pop("names"), sc.pop("bits"))):
            top = 2 ** b - 1
            black, white = limits[c] if c < len(limits) else (0, top)
            sc["channels"].append((name, top, int(max(0, black)), int(min(top, max(black + 1, white)))))
        sc["ets"] = _ets_file(path, k, len(infos), image_stacks if len(scenes) == 1 else 0) if k is not None else None
        try:
            sc["bytes"] = os.path.getsize(sc["ets"]) if sc["ets"] else None
        except OSError:
            sc["bytes"] = None
    return scenes


# ------------------------------------------------------------------------------------------------------------------------ pixels

def read(path, index, c, z, rect, size=None):
    """Channel c of plane z of scene `index`, the part rect = (x, y, w, h) of it in full-resolution pixels, at `size` (w, h) (None:
    full resolution): slideio reads it from the pyramid level that has that much detail, scaled to `size`."""
    sc = _scene(path, index)
    a = sc.read_block(tuple(int(v) for v in rect), size=(0, 0) if size is None else tuple(int(v) for v in size),
                      channel_indices=[c], slices=(z, z + 1))
    while a.ndim > 2:
        a = a[..., 0] if a.shape[-1] == 1 else a[0]
    return a


def pow2_below(step):
    """The largest power of 2 not above `step` (at least 1): the step to read at for a view that needs `step` (a little sharper)."""
    return 1 << (max(1, int(step)).bit_length() - 1)


def pow2_above(step):
    """The smallest power of 2 not below `step`."""
    return 1 << (max(1, int(math.ceil(step))) - 1).bit_length()


def read_at_step(path, index, c, z, x, y, w, h, step):
    """The part (x, y, w, h) of a channel of a plane at `step`: ceil(h / step) x ceil(w / step) pixels, each standing for the step x step
    full-resolution pixels from (x + i * step, y + j * step). Read as a whole number of steps across and down (so that nothing is
    stretched to fit: see the module's note); a last, partial step at the right or bottom edge repeats the pixel before it."""
    cols, rows = math.ceil(w / step), math.ceil(h / step)
    if step == 1:
        return _fit(read(path, index, c, z, (x, y, w, h)), rows, cols)
    fw, fh = (w // step) * step, (h // step) * step
    fw, fh = fw or w, fh or h                           # (narrower than one step: that bit, as one pixel)
    a = read(path, index, c, z, (x, y, fw, fh), (max(1, fw // step), max(1, fh // step)))
    return _fit(a, rows, cols)


def overview(path, index, c, z, width, height, step):
    """Every step-th pixel of a channel of a plane, whole: (ceil(height / step), ceil(width / step)), as the viewer's overviews are."""
    return read_at_step(path, index, c, z, 0, 0, width, height, step)


def _fit(a, rows, cols):
    """`a` cut or padded (with its last row / column) to rows x cols: slideio's scaled reads can be a pixel off at an edge."""
    if a.shape[0] > rows or a.shape[1] > cols:
        a = a[:rows, :cols]
    if a.shape[0] < rows or a.shape[1] < cols:
        a = np.pad(a, ((0, rows - a.shape[0]), (0, cols - a.shape[1])), mode="edge")
    return np.ascontiguousarray(a)


class BlockReader:
    """Rows of channel c of plane z of a scene at a step, as xlef_viewer._RowReader gives them: row(r, x0, n, step) is every step-th
    pixel of row r, n of them from column x0 (r and x0 multiples of step). Read in blocks of BLOCK x BLOCK pixels at that step, each
    once, kept until close(); several threads may ask at once (a block asked for while it is being read waits for that read)."""

    def __init__(self, path, index, c, z, width, height, dtype):
        self.path, self.index, self.c, self.z, self.width, self.height = path, index, c, z, width, height
        self.out_dtype = np.dtype(dtype).newbyteorder("=")
        self.blocks, self.lock = {}, threading.Lock()

    def _block(self, step, bx, by):
        key = (step, bx, by)
        with self.lock:
            got = self.blocks.get(key)
            if got is None:
                got = self.blocks[key] = [threading.Event(), None, None]
                mine = True
            else:
                mine = False
        if mine:
            try:
                x, y = bx * BLOCK * step, by * BLOCK * step
                w, h = min(self.width, x + BLOCK * step) - x, min(self.height, y + BLOCK * step) - y
                got[1] = read_at_step(self.path, self.index, self.c, self.z, x, y, w, h, step)
            except Exception as e:
                got[2] = e
                with self.lock:
                    self.blocks.pop(key, None)          # (asked again later: read again)
            finally:
                got[0].set()
        else:
            got[0].wait()
        if got[2] is not None:
            raise got[2]
        return got[1]

    def row(self, r, x0, n, step):
        out = np.empty(n, self.out_dtype)
        R, X = r // step, x0 // step                     # (in pixels at this step)
        by, k = divmod(R, BLOCK)
        i = 0
        while i < n:
            bx, j = divmod(X + i, BLOCK)
            block = self._block(step, bx, by)
            take = min(n - i, block.shape[1] - j)
            if take <= 0:
                raise OSError("a row past the image's edge")
            out[i:i + take] = block[k, j:j + take]
            i += take
        return out

    def cost(self, n, step):
        """About the bytes read for row(r, x0, n, step): a share of the compressed tiles (JPEG 2000: about a quarter of the pixels')."""
        return max(1, n * self.out_dtype.itemsize // 4)

    def close(self):
        with self.lock:
            self.blocks.clear()
