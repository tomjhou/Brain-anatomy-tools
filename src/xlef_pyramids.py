"""Pyramids for a Leica experiment that has none, as LAS X makes them, in a NEW project beside the experiment: <name>_pyramid.xlef.

LAS X keeps, for a big image, reduced copies of it ("pyramid levels"): 8-bit, each half the size of the one before, with a small .xlif
describing each. The viewer reads a view from the smallest level that is sharp enough, which is far quicker than from the full-resolution
TIFFs (see xlef_viewer: pyramid_level). This makes the same, so that an experiment LAS X made no pyramid for reads as quickly.

The experiment's own files are left alone. Beside <name>.xlef it writes
    <name>_pyramid.xlef
    <name>_pyramid/Metadata/            a copy of every .xlef/.xlcf/.xlif of the experiment, their references to the pixel files pointed
                                        back at the originals (the pixels are not copied), each image that gets a pyramid also naming it
    <name>_pyramid/Pyramid/             the levels' TIFFs (8-bit, uncompressed, so that rows can be read straight from them)
    <name>_pyramid/Pyramid/Metadata/    and their .xlif files
The levels are encoded the way LAS X's are: level = 255 * (original / max) ** (1 / 1.8) (see xlef_viewer.PYRAMID_GAMMA)."""

import concurrent.futures
import os
import re
import threading
import uuid
from urllib.parse import quote, unquote

import numpy as np
import tifffile

GAMMA = 1.8
SMALLEST = 512          # levels are made until the longest side of the next would be under this
MIN_IMAGE = 1024        # an image whose longest side is less than this gets no pyramid: its own overview is small already
SUFFIX = "_pyramid"
LEVEL_WORKERS = 3       # planes made at once (each: a whole plane in memory, up to ~0.5 GB)


class Cancelled(Exception):
    """The user stopped the job."""


def project_files(project_path):
    """(the new .xlef, the folder beside it for everything else) for the experiment at project_path."""
    base = os.path.splitext(project_path)[0] + SUFFIX
    return base + ".xlef", base


def level_sizes(width, height):
    """(width, height) of each pyramid level: half, a quarter, ... of the image, rounded down (as LAS X's are), while the longest side
    is at least SMALLEST (the first level always)."""
    out, k = [], 1
    while (width >> k) >= 1 and (height >> k) >= 1:
        size = (width >> k, height >> k)
        if out and max(size) < SMALLEST:
            break
        out.append(size)
        k += 1
    return out


def why_not(entry, is_online_only):
    """None if `entry` (an ImageEntry) can get a pyramid, else the reason it cannot (for the list the user is shown)."""
    if entry.read_only or entry.lif_block is not None:
        return "its data is inside a .lif file"
    if entry.problem:
        return "the viewer cannot read it"
    if entry.pyramid:
        return "it has a pyramid already"
    if max(entry.width, entry.height) < MIN_IMAGE:
        return "it is small already"
    wanted = [(z, c) for z in range(entry.nz) for c in range(len(entry.channels))]
    if not entry.channels or any(k not in entry.files for k in wanted):
        return "some of its files are missing"
    if any(is_online_only(f) for f in entry.all_files()):
        return "its files are online-only (download them first)"
    return None


def estimate(entries):
    """(planes to make, extra bytes on disk, rough seconds) for these ImageEntries: each plane gives levels of 1/4, 1/16, ... of its
    pixels, one byte each. The time is a guess (measured 2026-11: a compressed plane of 231 million pixels took 2.8 s in one thread, 84 million
    pixels/s; three at a time, say 160 million)."""
    planes = sum(e.nz * len(e.channels) for e in entries)
    pixels = sum(e.width * e.height * e.nz * len(e.channels) for e in entries)
    return planes, int(pixels / 3), pixels / 160e6


def _resolve(ref, base_dir):
    rel = unquote(ref).replace("\\", os.sep).replace("/", os.sep)
    return os.path.normpath(os.path.join(base_dir, rel))


def _ref(path, base_dir):
    """A reference to `path` from a file in base_dir, as LAS X writes them ("..%5Cmetadata%5Cx%20y.xlif")."""
    return quote(os.path.relpath(path, base_dir).replace(os.sep, "\\"), safe="")


def _read_text(path):
    with open(path, "rb") as f:
        raw = f.read()
    bom = raw.startswith(b"\xef\xbb\xbf")
    return raw.decode("utf-8-sig"), bom


def _write_text(path, text, bom=True):
    with open(path, "wb") as f:
        f.write((b"\xef\xbb\xbf" if bom else b"") + text.encode("utf-8"))


def _gamma_table(top):
    """8-bit pyramid value for each value 0..top of the image."""
    return np.rint(255 * (np.arange(top + 1, dtype=np.float64) / top) ** (1 / GAMMA)).astype(np.uint8)


def first_level(plane, size, table):
    """The pyramid's first level: 2 x 2 blocks of the plane averaged, then encoded 8-bit. size: (width, height) of the level."""
    lw, lh = size
    out = np.empty((lh, lw), np.uint8)
    for r in range(0, lh, 256):                 # (in chunks of rows: a plane is hundreds of MB, its temporaries would be several times that)
        r1 = min(lh, r + 256)
        blk = plane[2 * r:2 * r1, :2 * lw]
        s = (blk[0::2, 0::2].astype(np.uint32) + blk[0::2, 1::2] + blk[1::2, 0::2] + blk[1::2, 1::2]) >> 2
        out[r:r1] = table[s]
    return out


def next_level(level, size):
    """The level after `level` (8-bit): 2 x 2 blocks averaged. size: (width, height) of the new one."""
    lw, lh = size
    out = np.empty((lh, lw), np.uint8)
    for r in range(0, lh, 512):
        r1 = min(lh, r + 512)
        blk = level[2 * r:2 * r1, :2 * lw]
        out[r:r1] = (blk[0::2, 0::2].astype(np.uint16) + blk[0::2, 1::2] + blk[1::2, 0::2] + blk[1::2, 1::2] + 2) >> 2
    return out


def _lengths(text):
    """{DimID: (Length, Unit)} of an image's .xlif, as written there."""
    out = {}
    for m in re.finditer(r"<DimensionDescription\b[^>]*>", text):
        tag = m.group(0)
        d, n, u = (re.search(rf'\b{a}="([^"]*)"', tag) for a in ("DimID", "Length", "Unit"))
        if d and n:
            out[d.group(1)] = (n.group(1), u.group(1) if u else "m")
    return out


class Job:
    """Make the pyramids for `entries` (the experiment's ImageEntries that can have one: see why_not) and the new project.
    read_plane(entry, z, c) gives one plane as a 2-D array. progress(done, total, name) is called from the worker threads after each
    plane. cancel is a threading.Event that stops the job (Cancelled); what was written by then is deleted again."""

    def __init__(self, project_path, entries, read_plane, progress=None, cancel=None):
        self.project_path = os.path.normpath(project_path)
        self.entries = list(entries)
        self.read_plane = read_plane
        self.progress = progress
        self.cancel = cancel or threading.Event()
        self.xlef, self.folder = project_files(self.project_path)
        self.meta = os.path.join(self.folder, "Metadata")
        self.pyr = os.path.join(self.folder, "Pyramid")
        self.pyr_meta = os.path.join(self.pyr, "Metadata")
        self.made = []                  # files and folders written by this job, for taking them back
        self.done = 0
        self.lock = threading.Lock()
        self.by_path = {os.path.normcase(os.path.abspath(e.xlif_path)): e for e in self.entries}
        self.copied = {}                # normcased original path -> where its copy is
        self.names = set()

    # --- the files this job writes

    def _make_dir(self, path):
        if not os.path.isdir(path):
            parent = os.path.dirname(path)
            if parent and not os.path.isdir(parent):
                self._make_dir(parent)
            os.mkdir(path)
            self.made.append(path)

    def _note(self, path):
        with self.lock:
            self.made.append(path)

    def _save(self, path, text):
        existed = os.path.exists(path)
        _write_text(path, text)
        if not existed:
            self._note(path)

    def _remove_everything_written(self):
        for path in reversed(self.made):
            try:
                os.rmdir(path) if os.path.isdir(path) else os.remove(path)
            except OSError:
                pass

    def check_cancelled(self):
        if self.cancel.is_set():
            raise Cancelled()

    # --- the pyramid levels

    def _level_name(self, entry, k):
        return f"{os.path.splitext(os.path.basename(entry.xlif_path))[0]}_pmd_{k}"

    def _tiff(self, entry, k, z, c):
        return os.path.join(self.pyr, f"{self._level_name(entry, k)}--z{z:02d}--c{c:02d}.tif")

    def _make_plane(self, entry, z, c, table, sizes):
        """All levels of plane (z, c): read it, make the levels, write each as a TIFF."""
        self.check_cancelled()
        plane = self.read_plane(entry, z, c)
        level = first_level(plane, sizes[0], table)
        del plane
        for k, size in enumerate(sizes):
            self.check_cancelled()
            if k:
                level = next_level(level, size)
            path = self._tiff(entry, k, z, c)
            existed = os.path.exists(path)
            tmp = path + ".part"
            try:
                tifffile.imwrite(tmp, level, photometric="minisblack")      # (uncompressed: rows are read straight from the file)
                os.replace(tmp, path)
            finally:
                if os.path.exists(tmp):
                    os.remove(tmp)
            if not existed:
                self._note(path)
        with self.lock:
            self.done += 1
            done = self.done
        if self.progress is not None:
            self.progress(done, self.total, entry.name)

    def _level_xlif(self, entry, k, size, text, nxt):
        """The .xlif of level k of `entry` (text: the image's own .xlif, for its physical lengths)."""
        w, h = size
        nz, nch = entry.nz, len(entry.channels)
        lengths = _lengths(text)
        (lx, ux), (ly, uy) = lengths.get("1", ("0", "m")), lengths.get("2", ("0", "m"))
        name = self._level_name(entry, k)
        chans = "".join(
            f'          <ChannelDescription DataType="0" ChannelTag="0" Resolution="8" NameOfMeasuredQuantity="" Min="0.000000e+00" '
            f'Max="2.550000e+02" Unit="" LUTName="{_attr(ch.lut)}" IsLUTInverted="0" BytesInc="{c * w * h * nz}" BitInc="0" />\n'
            for c, ch in enumerate(entry.channels))
        dims = (f'          <DimensionDescription DimID="1" NumberOfElements="{w}" Origin="0.000000e+00" Length="{lx}" Unit="{ux}" '
                f'BitInc="0" BytesInc="1" />\n'
                f'          <DimensionDescription DimID="2" NumberOfElements="{h}" Origin="0.000000e+00" Length="{ly}" Unit="{uy}" '
                f'BitInc="0" BytesInc="{w}" />\n')
        if nz > 1:
            lz, uz = lengths.get("3", ("0", "m"))
            dims += (f'          <DimensionDescription DimID="3" NumberOfElements="{nz}" Origin="0.000000e+00" Length="{lz}" Unit="{uz}" '
                     f'BitInc="0" BytesInc="{w * h}" />\n')
        frames = "".join(
            f'      <Frame File="{_ref(self._tiff(entry, k, z, c), self.pyr_meta)}" Offset="{(c * nz + z) * w * h}" Size="{w * h}" '
            f'UUID="{uuid.uuid4()}" />\n' for c in range(nch) for z in range(nz))
        children = (f'    <Children>\n      <Reference File="{_ref(os.path.join(self.pyr_meta, nxt + ".xlif"), self.pyr_meta)}" '
                    f'UUID="{uuid.uuid4()}" />\n    </Children>\n') if nxt else "    <Children />\n"
        return (
            '<?xml version="1.0" encoding="utf-8"?>\n'
            '<LMSDataContainerHeader Version="2" ApplicationDescription="LMSIOManager" ApplicationVersion="1.15.5732.0">\n'
            f'  <Element Name="{_attr(name)}" Visibility="2" CopyOption="2" UniqueID="{uuid.uuid4()}">\n'
            '    <Data>\n      <Image TextDescription="">\n        <ImageDescription>\n          <Channels>\n'
            + chans + '          </Channels>\n          <Dimensions>\n' + dims + '          </Dimensions>\n        </ImageDescription>\n'
            '        <Attachment Name="ImagePyramid" Application="LAS AF" NumberValidStacks="1" Gamma="1.0000000000" />\n'
            '      </Image>\n    </Data>\n'
            '    <Attributes>\n      <Attribute>Pyramid</Attribute>\n      <Attribute>PreventLIFSave</Attribute>\n    </Attributes>\n'
            f'    <Memory Size="{w * h * nz * nch}" MemoryBlockID="">\n' + frames + '    </Memory>\n'
            + children + '  </Element>\n</LMSDataContainerHeader>\n')

    def _write_levels(self, entry, text):
        sizes = level_sizes(entry.width, entry.height)
        for k, size in enumerate(sizes):
            nxt = self._level_name(entry, k + 1) if k + 1 < len(sizes) else None
            self._save(os.path.join(self.pyr_meta, self._level_name(entry, k) + ".xlif"), self._level_xlif(entry, k, size, text, nxt))

    # --- the copy of the experiment's metadata

    def _new_name(self, src):
        base = os.path.basename(src)
        stem, ext = os.path.splitext(base)
        name, n = base, 1
        while name.lower() in self.names:
            n += 1
            name = f"{stem}_{n}{ext}"
        self.names.add(name.lower())
        return name

    def _copy_file(self, src, dest=None):
        """A copy of the .xlef/.xlcf/.xlif at `src` in the new Metadata folder (or at `dest`), its references pointed at the right
        files; the copies of the files it refers to are made too. Returns the copy's path."""
        key = os.path.normcase(os.path.abspath(src))
        if key in self.copied:
            return self.copied[key]
        dest = dest or os.path.join(self.meta, self._new_name(src))
        self.copied[key] = dest
        text, bom = _read_text(src)
        old_dir, new_dir = os.path.dirname(os.path.abspath(src)), os.path.dirname(dest)

        def fix(m):
            tag, ref = m.group("tag"), m.group("val")
            target = _resolve(ref, old_dir)
            if tag == "Reference" and target.lower().endswith((".xlif", ".xlcf", ".xlef")) and os.path.isfile(target):
                target = self._copy_file(target)
            return m.group("head") + _ref(target, new_dir) + '"'
        text = re.sub(r'(?P<head><(?P<tag>[A-Za-z]+)\b[^>]*?\bFile=")(?P<val>[^"]*)"', fix, text)
        entry = self.by_path.get(key)
        if entry is not None:
            self._write_levels(entry, text)
            first = _ref(os.path.join(self.pyr_meta, self._level_name(entry, 0) + ".xlif"), new_dir)
            ref = f'<Reference File="{first}" UUID="{uuid.uuid4()}" />'
            if re.search(r"<Children\s*/>", text):
                text = re.sub(r"<Children\s*/>", f"<Children>\n      {ref}\n    </Children>", text, count=1)
            elif "</Children>" in text:
                text = text.replace("</Children>", f"  {ref}\n    </Children>", 1)
            else:
                text = text.replace("</Element>", f"  <Children>\n      {ref}\n    </Children>\n  </Element>", 1)
        if dest == self.xlef:
            text = re.sub(r'(<Element\b[^>]*?\bName=")[^"]*"', lambda m: m.group(1) + _attr(os.path.basename(dest)) + '"', text, count=1)
            text = re.sub(r'(<Experiment\b[^>]*?\bPath=")[^"]*"', lambda m: m.group(1) + _attr(dest.lower()) + '"', text, count=1)
        self._save(dest, text)
        return dest

    # --- the whole job

    def run(self):
        """Make everything; the new project's .xlef path. Raises Cancelled (after taking back what was written) or the error that stopped
        it (likewise)."""
        try:
            for d in (self.folder, self.meta, self.pyr, self.pyr_meta):
                self._make_dir(d)
            jobs = []
            self.total = sum(e.nz * len(e.channels) for e in self.entries)
            with concurrent.futures.ThreadPoolExecutor(max_workers=LEVEL_WORKERS) as pool:
                for e in self.entries:
                    sizes = level_sizes(e.width, e.height)
                    table = _gamma_table(max(c.max_value for c in e.channels))
                    jobs += [pool.submit(self._make_plane, e, z, c, table, sizes) for c in range(len(e.channels)) for z in range(e.nz)]
                try:
                    for j in concurrent.futures.as_completed(jobs):
                        j.result()
                except BaseException:
                    self.cancel.set()
                    for j in jobs:
                        j.cancel()
                    raise
            self.check_cancelled()
            self._copy_file(self.project_path, dest=self.xlef)
            return self.xlef
        except BaseException:
            self.cancel.set()
            self._remove_everything_written()
            raise


def _attr(value):
    return str(value).replace("&", "&amp;").replace('"', "&quot;").replace("<", "&lt;").replace(">", "&gt;")
