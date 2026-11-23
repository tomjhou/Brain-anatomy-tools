"""A .lif (LAS X's single-file experiment) saved as an .xlef experiment, laid out the way LAS X itself writes one, so that it can be
edited (the viewer only reads a .lif) and opened again in LAS X. The .lif (and its .lifext) is only read.

Beside <name>.lif it writes the folder <name>/ (as LAS X's export, compared file by file with LAS X's own .xlef of the same experiment):
    <name>/<name>.xlef                      the experiment: a <Reference> to each of its folders and images
    <name>/Metadata/<image>.xlif            each image: its <Element> from the .lif unchanged, but for <Memory> (a <Frame> per TIFF)
    <name>/<image>--Z00--C00.tif ...        one uncompressed TIFF per channel (and z layer): "<image>--C00.tif" without layers,
                                            "<image>.tif" for a single channel; the raw pixels of the .lif, byte for byte
    <name>/<folder>/Metadata/<folder>.xlcf  a folder (LAS X's "collection", e.g. a slide), its images in <name>/<folder>/ as above
    <name>/Additional Data/                 images LAS X hides (Visibility 2: the sample overview)
    <name>/Pyramid/                         LAS X's reduced copies from the .lifext (<image>_pmd_0 ...), one TIFF per channel and
    <name>/Pyramid/Metadata/                layer, and each image's histogram (.lof); an image in a folder gets the folder's name
                                            first ("Slide1_Section1_pmd_0")
An image whose pixels are not one run of rows per channel and layer (colour, tiles, time series, floating point) is written as one
.lof file instead (LAS X's raw format: the .lif's data block with a short header in front and its XML behind), as LAS X's "LOF"
export does; the viewer and LAS X read both."""

import copy
import os
import re
import threading
import uuid
import xml.etree.ElementTree as ET
from urllib.parse import quote

HEADER_ATTRS = {"Version": "2", "ApplicationDescription": "LMSIOManager", "ApplicationVersion": "1.15.5732.0"}  # (as LAS X's files)
TIFF_DESCRIPTION = "<LMSIOManager><ProductVersion>1,15,5732,0</ProductVersion><BitShift>UnProcessed</BitShift></LMSIOManager>"
CHUNK = 32 << 20               # bytes copied at a time
TIFF_LIMIT = (1 << 32) - (16 << 20)     # a plane bigger than this doesn't fit in a classic TIFF: written as a .lof instead
BAD_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
HIDDEN_FOLDER = "Additional Data"
PYRAMID_FOLDER = "Pyramid"
META = "Metadata"


class Cancelled(Exception):
    """The user stopped the job."""


def target_paths(lif_path):
    """(the new .xlef, its folder) for the .lif at lif_path: <name>/<name>.xlef beside it; if a folder <name> is there already (without
    that .xlef in it), <name> (2), <name> (3), ..."""
    lif_path = os.path.abspath(lif_path)
    here, stem = os.path.dirname(lif_path), os.path.splitext(os.path.basename(lif_path))[0]
    n = 1
    while True:
        folder = os.path.join(here, stem if n == 1 else f"{stem} ({n})")
        xlef = os.path.join(folder, stem + ".xlef")
        if os.path.isfile(xlef) or not os.path.exists(folder) or (os.path.isdir(folder) and not os.listdir(folder)):
            return xlef, folder
        n += 1


def paths_in(lif_path, parent):
    """(the new .xlef, its folder) for the .lif at lif_path, in the folder `parent` chosen for it: <parent>/<name>/<name>.xlef, as
    beside the .lif (target_paths); <name> (2), <name> (3), ... if a folder <name> with files in it is there already. If `parent` is
    itself an empty folder called <name> (made in the dialog for it), the copy goes straight into it."""
    parent = os.path.abspath(parent)
    stem = os.path.splitext(os.path.basename(lif_path))[0]
    if os.path.basename(parent).lower() == stem.lower() and os.path.isdir(parent) and not os.listdir(parent):
        return os.path.join(parent, stem + ".xlef"), parent
    n = 1
    while True:
        folder = os.path.join(parent, stem if n == 1 else f"{stem} ({n})")
        if not os.path.exists(folder) or (os.path.isdir(folder) and not os.listdir(folder)):
            return os.path.join(folder, stem + ".xlef"), folder
        n += 1


def existing_copy(lif_path):
    """The .xlef made from this .lif earlier (target_paths), or None."""
    xlef, _ = target_paths(lif_path)
    return xlef if os.path.isfile(xlef) else None


def needed_bytes(lif_path):
    """Roughly what the copy takes on disk: the .lif and its .lifext (their pixels are copied; the metadata is small)."""
    total = os.path.getsize(lif_path)
    ext = os.path.splitext(lif_path)[0] + ".lifext"
    return total + (os.path.getsize(ext) if os.path.isfile(ext) else 0)


# --- reading the .lif / .lifext

def _u32(b):
    return int.from_bytes(b[:4], "little")


def read_header(path):
    """A .lif's (or .lifext's) XML, parsed (its namespace prefixes kept for writing it again), and where its first data block's header
    is. Reads only the start of the file."""
    with open(path, "rb") as fh:
        head = fh.read(13)
        if len(head) < 13 or _u32(head) != 0x70 or head[8:9] != b"*":
            raise ValueError(f"{os.path.basename(path)} is not a LAS X .lif file")
        n = _u32(head[9:13])
        text = fh.read(n * 2).decode("utf-16-le")
    parser = ET.XMLPullParser(events=("start-ns",))
    parser.feed(text)
    for _, (prefix, uri) in parser.read_events():
        if prefix:
            ET.register_namespace(prefix, uri)
    return ET.fromstring(text), 13 + n * 2


def scan_blocks(path, first):
    """{block id: (where its data starts, its size)} of a .lif / .lifext, from the header in front of each block (0x70, the header's
    length, '*', the size, '*', the id's length and the id), one after another."""
    blocks, size = {}, os.path.getsize(path)
    with open(path, "rb", buffering=0) as fh:
        pos = first
        while pos + 22 <= size:
            fh.seek(pos)
            head = fh.read(22)
            if len(head) < 22 or _u32(head) != 0x70:
                break
            n = _u32(head[18:22])
            block = fh.read(2 * n).decode("utf-16-le", "replace")
            start, length = pos + 22 + 2 * n, int.from_bytes(head[9:17], "little")
            blocks[block] = (start, length)
            pos = start + length
    return blocks


def _memory(el):
    """(block id, size) of an element's data ("" and 0 if it has none)."""
    m = el.find("Memory")
    if m is None:
        return "", 0
    try:
        return m.get("MemoryBlockID") or "", int(m.get("Size") or 0)
    except ValueError:
        return m.get("MemoryBlockID") or "", 0


def _is_item(el):
    """An image (or other data: an .xlif of its own), not a folder (.xlcf)."""
    return el.find("Data/Image") is not None or _memory(el)[1] > 0


def _kids(el):
    k = el.find("Children")
    return [c for c in (k if k is not None else []) if c.tag == "Element"]


def plane_layout(el):
    """How an image's data block splits into TIFFs: {"width", "height", "bytes" (per pixel), "nz", "nch", "plane" (bytes),
    "planes": [(c, z, offset in the block)], "pixel_m": (x, y) pixel size in metres or None}, channel after channel and layer after
    layer as in LAS X's .xlif; None if it doesn't (colour, tiles, time, floating point...: written as a .lof)."""
    image = el.find("Data/Image")
    desc = image.find("ImageDescription") if image is not None else None
    if desc is None:
        return None
    try:
        dims = {d.get("DimID"): d for d in desc.iter("DimensionDescription")}
        chans = list(desc.iter("ChannelDescription"))
        w, h = int(dims["1"].get("NumberOfElements")), int(dims["2"].get("NumberOfElements"))
        x_inc, y_inc = int(dims["1"].get("BytesInc")), int(dims["2"].get("BytesInc"))
        if x_inc not in (1, 2) or y_inc != x_inc * w or not chans or w <= 0 or h <= 0:
            return None
        if any(k not in ("1", "2", "3") and int(d.get("NumberOfElements") or 1) > 1 for k, d in dims.items()):
            return None
        if any((c.get("DataType") or "0") != "0" or int(c.get("Resolution") or 0) > 8 * x_inc for c in chans):
            return None
        nz = int(dims["3"].get("NumberOfElements")) if "3" in dims else 1
        z_inc = int(dims["3"].get("BytesInc")) if "3" in dims else 0
        plane = y_inc * h
        planes = [(c, z, int(cd.get("BytesInc") or 0) + z * z_inc) for c, cd in enumerate(chans) for z in range(nz)]
    except (KeyError, TypeError, ValueError):
        return None
    _, size = _memory(el)
    spans = sorted(at for _, _, at in planes)
    if spans[0] < 0 or spans[-1] + plane > size or any(b - a < plane for a, b in zip(spans, spans[1:])) or plane > TIFF_LIMIT:
        return None

    def pixel(d):          # (as LAS X works it out for its TIFFs: the length over n - 1 pixels)
        try:
            n, length = int(d.get("NumberOfElements")), float(d.get("Length"))
            scale = {"m": 1.0, "mm": 1e-3, "um": 1e-6, "µm": 1e-6}.get((d.get("Unit") or "m").strip(), None)
            return length * scale / (n - 1) if scale and n > 1 and length > 0 else None
        except (TypeError, ValueError):
            return None
    px = pixel(dims["1"])
    return {"width": w, "height": h, "bytes": x_inc, "nz": nz, "nch": len(chans), "plane": plane, "planes": planes,
            "pixel_m": (px, px) if px else None}            # (LAS X gives both directions the x pixel size)


def tiff_name(stem, c, z, nz, nch):
    """A channel's TIFF as LAS X names them: "<stem>--Z03--C01.tif", "<stem>--C01.tif" (no layers), "<stem>.tif" (one channel)."""
    parts = [stem] + ([f"Z{z:02d}"] if nz > 1 else []) + ([f"C{c:02d}"] if nch > 1 or nz > 1 else [])
    return "--".join(parts) + ".tif"


def _ref(target, base_dir):
    """A reference to `target` from a file in base_dir, as LAS X writes them: relative, lower case, backslashes, URL-encoded
    (".%5Cmetadata%5Csection%2017.xlif")."""
    rel = os.path.relpath(target, base_dir).replace("/", "\\")
    if not rel.startswith(".."):
        rel = ".\\" + rel
    return quote(rel.lower(), safe="")


def safe_name(name):
    """A file name made from an element's name (characters Windows doesn't allow become _)."""
    s = BAD_CHARS.sub("_", name or "").strip().rstrip(".")
    return s or "Image"


# --- writing

def tiff_header_and_ifd(width, height, bytes_per_pixel, pixel_m):
    """(the 8 bytes at the start of a TIFF, the directory written after the pixels) for one uncompressed plane whose rows start at byte 8,
    one row per strip, as LAS X writes them (with its description, and the resolution in pixels per cm)."""
    plane = width * height * bytes_per_pixel
    ifd_at = 8 + plane + (plane & 1)
    tags = []                          # (tag, type, count, value bytes)

    def short(tag, v):
        tags.append((tag, 3, 1, v.to_bytes(2, "little")))

    def long_(tag, values):
        tags.append((tag, 4, len(values), b"".join(v.to_bytes(4, "little") for v in values)))
    long_(256, [width])
    long_(257, [height])
    short(258, 8 * bytes_per_pixel)
    short(259, 1)
    short(262, 1)
    tags.append((270, 2, len(TIFF_DESCRIPTION) + 1, TIFF_DESCRIPTION.encode("ascii") + b"\0"))
    row = width * bytes_per_pixel
    long_(273, [8 + r * row for r in range(height)])
    short(277, 1)
    long_(278, [1])
    long_(279, [row] * height)
    if pixel_m:
        for tag, size in ((282, pixel_m[0]), (283, pixel_m[1])):
            per_cm = 0.01 / size
            num = 0xFFFFFFFF
            den = max(1, min(0xFFFFFFFF, round(num / per_cm)))
            tags.append((tag, 5, 1, num.to_bytes(4, "little") + den.to_bytes(4, "little")))
        short(296, 3)
    tags.sort(key=lambda t: t[0])
    entries, extra = b"", b""
    extra_at = ifd_at + 2 + 12 * len(tags) + 4
    for tag, typ, count, data in tags:
        if len(data) <= 4:
            value = data.ljust(4, b"\0")
        else:
            if (extra_at + len(extra)) & 1:
                extra += b"\0"
            value = (extra_at + len(extra)).to_bytes(4, "little")
            extra += data
        entries += tag.to_bytes(2, "little") + typ.to_bytes(2, "little") + count.to_bytes(4, "little") + value
    if ifd_at + 2 + len(entries) + 4 + len(extra) >= 1 << 32:
        raise ValueError("the plane is too big for a TIFF")
    head = b"II*\0" + ifd_at.to_bytes(4, "little")
    ifd = (b"\0" if plane & 1 else b"") + len(tags).to_bytes(2, "little") + entries + b"\0\0\0\0" + extra
    return head, ifd


def lof_header(size):
    """The block in front of a .lof's data (as LAS X writes it): 0x70, its length, then '*' "LMS_Object_File", '*' 1, '*' 0, '*' size."""
    name = "LMS_Object_File"
    body = (b"*" + len(name).to_bytes(4, "little") + name.encode("utf-16-le") + b"*" + (1).to_bytes(4, "little") + b"*"
            + (0).to_bytes(4, "little") + b"*" + size.to_bytes(8, "little"))
    return (0x70).to_bytes(4, "little") + len(body).to_bytes(4, "little") + body


def lof_trailer(element):
    """The block after a .lof's data: the element's XML (UTF-16), as LAS X writes it."""
    root = ET.Element("LMSDataContainerHeader", {"Version": "2"})
    root.append(copy.deepcopy(element))
    text = ET.tostring(root, encoding="unicode", short_empty_elements=True)
    body = b"*" + len(text).to_bytes(4, "little") + text.encode("utf-16-le")
    return (0x70).to_bytes(4, "little") + len(body).to_bytes(4, "little") + body


def xml_text(element):
    """An .xlif/.xlcf/.xlef's text: the element in LAS X's container, indented as LAS X does, with CRLF line ends."""
    root = ET.Element("LMSDataContainerHeader", dict(HEADER_ATTRS))
    root.append(element)
    ET.indent(root, "  ")
    text = '<?xml version="1.0" encoding="utf-8"?>\n' + ET.tostring(root, encoding="unicode")
    text = text.replace("\r\n", "\n").replace("\n", "\r\n")
    return text.replace("&#13;", "&#xD;").replace("&#10;", "&#xA;").replace("&#09;", "&#x9;")      # (line breaks in values, as LAS X writes them)


class _Source:
    """A .lif or .lifext being read: its file and where its blocks are."""

    def __init__(self, path):
        self.path = path
        self.root, first = read_header(path)
        self.blocks = scan_blocks(path, first)
        self.fh = None

    def start(self, block):
        if block not in self.blocks:
            raise ValueError(f"{os.path.basename(self.path)} has no data block {block}")
        return self.blocks[block][0]


class Job:
    """Save the .lif at lif_path as an .xlef (see the module's notes). progress(done, total, name) is called (from the thread that
    runs it) as the pixels are copied, in bytes. cancel is a threading.Event that stops it (Cancelled); what was written by then is
    deleted again. After run(): written[UniqueID] is the new .xlif/.xlcf of each of the .lif's elements (and the experiment's .xlef)."""

    def __init__(self, lif_path, progress=None, cancel=None, xlef_path=None):
        """xlef_path: where the new .xlef goes (its folder holds the rest, see paths_in); target_paths() if None."""
        self.lif_path = os.path.abspath(lif_path)
        self.progress = progress
        self.cancel = cancel or threading.Event()
        if xlef_path:
            self.xlef = os.path.abspath(xlef_path)
            self.folder = os.path.dirname(self.xlef)
        else:
            self.xlef, self.folder = target_paths(self.lif_path)
        self.made = []                  # files and folders written, for taking them back
        self.names = {}                 # normcase(folder) -> lower-case names used in it
        self.written = {}
        self.done = self.total = 0
        self.lof_images = []            # names of the images written as .lof (their layout doesn't split into TIFFs)

    # --- files and names

    def _make_dir(self, path):
        if not os.path.isdir(path):
            parent = os.path.dirname(path)
            if parent and not os.path.isdir(parent):
                self._make_dir(parent)
            os.mkdir(path)
            self.made.append(path)

    def _unique(self, folder, name):
        used = self.names.setdefault(os.path.normcase(folder), set())
        stem, n = name, 1
        while stem.lower() in used:
            n += 1
            stem = f"{name}_{n}"
        used.add(stem.lower())
        return stem

    def _save(self, path, text):
        self._make_dir(os.path.dirname(path))
        self.made.append(path)
        with open(path, "wb") as f:
            f.write(b"\xef\xbb\xbf" + text.encode("utf-8"))

    def check_cancelled(self):
        if self.cancel.is_set():
            raise Cancelled()

    def _copy(self, src, start, n, out, name):
        """n bytes of the source from `start` into the open file `out`."""
        if src.fh is None:
            src.fh = open(src.path, "rb", buffering=0)
        src.fh.seek(start)
        buf = memoryview(bytearray(min(CHUNK, max(n, 1))))
        left = n
        while left:
            self.check_cancelled()
            k = src.fh.readinto(buf[:min(len(buf), left)])
            if not k:
                raise OSError(f"{os.path.basename(src.path)} ends before the data of {name}")
            out.write(buf[:k])
            left -= k
            self.done += k
            if self.progress is not None:
                self.progress(self.done, self.total, name)

    def _write_tiff(self, path, src, start, layout, name):
        head, ifd = tiff_header_and_ifd(layout["width"], layout["height"], layout["bytes"], layout["pixel_m"])
        self.made.append(path)
        with open(path, "wb") as out:
            out.write(head)
            self._copy(src, start, layout["plane"], out, name)
            out.write(ifd)

    def _write_lof(self, path, src, start, size, element, name):
        self.made.append(path)
        with open(path, "wb") as out:
            out.write(lof_header(size))
            self._copy(src, start, size, out, name)
            out.write(lof_trailer(element))

    # --- the experiment's tree

    def _write_item(self, el, src, folder, own_prefix, folder_prefix, pyramid_of):
        """An image (or other data) element: its pixels (TIFFs or a .lof) in `folder`, its .xlif in folder/Metadata, its pyramid (from
        pyramid_of: {block id: [elements]} of the .lifext) in Pyramid/. Returns the .xlif's path."""
        stem = self._unique(folder, safe_name(own_prefix + (el.get("Name") or "")))
        meta = os.path.join(folder, META)
        xlif = os.path.join(meta, stem + ".xlif")
        self._make_dir(meta)
        block, size = _memory(el)
        new = copy.deepcopy(el)
        mem = new.find("Memory")
        if mem is None:
            mem = ET.SubElement(new, "Memory", {"Size": "0"})
        for child in list(mem):
            mem.remove(child)
        mem.set("MemoryBlockID", "")
        if size > 0 and block:
            start = src.start(block)
            layout = plane_layout(el)
            name = el.get("Name") or stem
            if layout is not None:
                for c, z, at in layout["planes"]:
                    tif = os.path.join(folder, tiff_name(stem, c, z, layout["nz"], layout["nch"]))
                    self._write_tiff(tif, src, start + at, layout, name)
                    ET.SubElement(mem, "Frame", {"File": _ref(tif, meta), "Offset": str(at), "Size": str(layout["plane"]),
                                                 "UUID": str(uuid.uuid4())})
            else:
                lof = os.path.join(folder, stem + ".lof")
                self._write_lof(lof, src, start, size, el, name)
                ET.SubElement(mem, "Block", {"File": _ref(lof, meta), "Offset": "0", "Size": str(size), "UUID": str(uuid.uuid4())})
                if el.find("Data/Image") is not None and pyramid_of is not None:
                    self.lof_images.append(name)
        refs = []
        for k in _kids(el):             # (a pyramid level names the next; an image in a .lif has none)
            refs.append((self._write_item(k, src, folder, own_prefix, folder_prefix, None), k))
        if pyramid_of is not None and block in pyramid_of:
            pyr = os.path.join(self.folder, PYRAMID_FOLDER)
            for k in pyramid_of[block]:
                refs.append((self._write_item(k, self.ext, pyr, folder_prefix, folder_prefix, None), k))
        self._set_children(new, [(p, k.get("UniqueID") or str(uuid.uuid4())) for p, k in refs], meta)
        self._save(xlif, xml_text(new))
        if el.get("UniqueID"):
            self.written[el.get("UniqueID")] = xlif
        return xlif

    def _write_folder(self, el, folder, own_file, prefix, pyramid_of):
        """A folder element (LAS X's collection, or the experiment): its children (images in `folder`, sub-folders in folders of their
        own), then its own file (own_file: its .xlcf, or the .xlef)."""
        refs = []
        for k in _kids(el):
            self.check_cancelled()
            if _is_item(k):
                where = os.path.join(folder, HIDDEN_FOLDER) if k.get("Visibility") == "2" else folder
                if where != folder:
                    self._unique(folder, HIDDEN_FOLDER)
                path = self._write_item(k, self.lif, where, "", prefix, pyramid_of)
            else:
                name = self._unique(folder, safe_name(k.get("Name")))
                sub = os.path.join(folder, name)
                path = self._write_folder(k, sub, os.path.join(sub, META, name + ".xlcf"), prefix + name + "_", pyramid_of)
            refs.append((path, k.get("UniqueID") or str(uuid.uuid4())))
        new = copy.deepcopy(el)
        mem = new.find("Memory")
        if mem is not None:
            for child in list(mem):
                mem.remove(child)
            mem.set("MemoryBlockID", "")
        self._set_children(new, refs, os.path.dirname(own_file))
        if own_file == self.xlef:
            self._as_experiment(new)
        self._save(own_file, xml_text(new))
        if el.get("UniqueID"):
            self.written[el.get("UniqueID")] = own_file
        return own_file

    @staticmethod
    def _set_children(new, refs, base_dir):
        """The element's <Children>: a <Reference> to each written file (after <Memory>, where LAS X has it)."""
        old = new.find("Children")
        index = list(new).index(old) if old is not None else len(new)
        if old is not None:
            new.remove(old)
        kids = ET.Element("Children")
        for path, uid in refs:
            ET.SubElement(kids, "Reference", {"File": _ref(path, base_dir), "UUID": uid})
        new.insert(index, kids)

    def _as_experiment(self, new):
        """The experiment's element as in LAS X's .xlef: named after it, its path, saved; not "being saved"."""
        new.set("Name", os.path.basename(self.xlef))
        exp = new.find("Data/Experiment")
        if exp is not None:
            exp.set("Path", self.xlef)
            exp.set("IsSavedFlag", "1")
        attrs = new.find("Attributes")
        if attrs is not None:
            for a in [a for a in attrs if (a.text or "") == "___Saving"]:
                attrs.remove(a)
            if not len(attrs):
                new.remove(attrs)

    def _experiment(self):
        """The .lif's experiment element; a .lif whose top is an image (or several tops) gets one made around them."""
        tops = self.lif.root.findall("Element")
        if len(tops) == 1 and tops[0].find("Data/Experiment") is not None:
            return tops[0]
        el = ET.Element("Element", {"Name": os.path.basename(self.xlef), "Visibility": "1", "CopyOption": "1",
                                    "UniqueID": str(uuid.uuid4())})
        ET.SubElement(ET.SubElement(el, "Data"), "Experiment", {"Path": self.xlef, "IsSavedFlag": "1"})
        ET.SubElement(el, "Memory", {"Size": "0", "MemoryBlockID": ""})
        kids = ET.SubElement(el, "Children")
        for t in tops:
            kids.append(t)
        return el

    def run(self):
        """Write everything; the new .xlef's path. Raises Cancelled (after taking back what was written) or the error that stopped it
        (likewise)."""
        self.lif = self.ext = None
        try:
            if os.path.exists(self.xlef):
                raise FileExistsError(f"{self.xlef} is there already")
            self.lif = _Source(self.lif_path)
            ext_path = os.path.splitext(self.lif_path)[0] + ".lifext"
            pyramid_of = {}
            if os.path.isfile(ext_path):
                try:
                    self.ext = _Source(ext_path)
                except (OSError, ValueError):
                    self.ext = None             # (an unreadable .lifext: no pyramid copied; the viewer or LAS X can make one again)
            if self.ext is not None:
                for g in self.ext.root.findall("ChildrenOf"):
                    els = [e for e in g if e.tag == "Element"]
                    if g.get("MemoryBlockID") and els:
                        pyramid_of[g.get("MemoryBlockID")] = els
            top = self._experiment()
            self.total = sum(_memory(e)[1] for e in top.iter("Element"))
            used = {_memory(e)[0] for e in top.iter("Element")}
            self.total += sum(_memory(e)[1] for b, els in pyramid_of.items() if b in used for g in els for e in g.iter("Element"))
            self._make_dir(self.folder)
            self._write_folder(top, self.folder, self.xlef, "", pyramid_of)
            return self.xlef
        except BaseException:
            self.cancel.set()
            self._close()
            for path in reversed(self.made):
                try:
                    os.rmdir(path) if os.path.isdir(path) else os.remove(path)
                except OSError:
                    pass
            raise
        finally:
            self._close()

    def _close(self):
        for src in (self.lif, self.ext):
            if src is not None and src.fh is not None:
                src.fh.close()
                src.fh = None
