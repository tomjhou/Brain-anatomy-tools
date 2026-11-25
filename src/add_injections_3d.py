"""Show the lab's injection sites (and optic fibers) in a 3D page written by the ABC atlas browser, as white balls.

    python add_injections_3d.py PAGE.html SURGERIES.csv [--radius 0.15] [--no-open]

PAGE.html is a 3D page from the ABC browser (its axes are CCF millimetres). SURGERIES.csv is the surgery viewer's "Export CSV" file
(JhouLab-surgery-viewer), so filter the list there first (by target, person, dates...) to choose which injections are drawn. The page is not copied: the sites
go in a small side file (injections-DATE-TIME.js in Brain-anatomy-tools/tmp/injection_sites, on this computer; the newest 30 are kept), and the page is opened
in a NEW browser tab with the side file named after a # in its address. The page reads it as it opens (see injection_loader.py), so each export gets its own tab
that keeps what it was opened with, a refresh included, and the tabs can be compared; close the ones you do not want. Two more legend entries come with it:
"injection sites", the balls, and "needle tracks", a white line from where each DV was measured from (the brain surface, or bregma's height) down to its ball,
at the injection angle. Optic fibers (see site_kinds) are drawn in blue, with a smaller ball and a thicker rod, as "optic fibers" and "fiber tracks".
Hovering any of them shows the animal, target, substance and coordinates. Opened without the # part, the page shows no sites.
A page written before it carried the loader is upgraded once, in place (it looks the same, and is stored compressed); see compressed_page.py.

From stereotaxic coordinates to the CCF (all approximate; there is no exact mapping):
  * Bregma is placed in the CCF where the ABC browser's Settings put it (viewer_settings.json;
    by default the IBL's estimate), with no scaling or tilt, as the browser does.
  * AP: when the row's "AP measured from" says the zero is some mm behind bregma (the
    "adjusted intersection point (~0.4 mm behind Paxinos bregma)"), that is subtracted.
  * ML: the side column picks the hemisphere; bilateral rows give two balls. A row with no
    side is drawn on the left (the lab's usual side) and says so when hovered.
  * DV: zeroed at bregma: measured down from bregma's height. Zeroed at the skull, dura or brain
    surface: measured down from the top of the CCF brain at the entry point (the CCF has no skull,
    so skull-zeroed sites land a little too deep). No zero given: the brain surface, the lab's
    usual, except for people who usually use bregma (BREGMA_DV_PEOPLE: Seung-Chan).
  * Angle: AP/ML are where the needle enters and DV is the distance along the tilted track,
    tilted in the coronal plane toward the midline (as in the lab's standard coordinates: VTA
    entered at ML 1.3, 4.8 mm at 10 degrees, ends near ML 0.5).
  * Rat rows, and rows without coordinates, are left out (and counted).

Needs numpy only (plus the CCF annotation volume from the ABC browser's atlas cache).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import struct
import sys
import time
import urllib.parse
from pathlib import Path

import numpy as np

import compressed_page
import slice_slider_scripts

SRC = Path(__file__).resolve().parent
SETTINGS_FILE = SRC / 'viewer_settings.json'                              # the ABC browser's (bregma position)
ATLAS_LOCATION_FILE = SRC / 'cache_local' / 'abc_atlas_cache_location.json'
DEFAULT_ATLAS_CACHE = SRC.parent.parent / 'abc_atlas_cache'
DEFAULT_BREGMA_AP_MM, DEFAULT_BREGMA_DV_MM, BREGMA_ML_MM = 5.40, 0.332, 5.739   # as abc_core (the IBL's bregma)
VOXEL_MM = 0.025
DEFAULT_RADIUS_MM = 0.15
BALL_COLOR = '#ffffff'
TRACK_WIDTH_PX = 3
FIBER_COLOR = '#3d8bff'      # optic fibers: blue, a smaller ball, a thicker rod
FIBER_RADIUS_SCALE = 0.5
FIBER_WIDTH_PX = 5
SIDE_FOLDER = SRC.parent / 'tmp' / 'injection_sites'      # (git ignores tmp/)
KEEP_SIDE_FILES = 30


# ---------- the CCF side ----------
def bregma_in_ccf(settings_file=SETTINGS_FILE):
    """(AP, DV, ML) of bregma in the CCF, mm, from the ABC browser's settings (its defaults if unreadable)."""
    ap, dv = DEFAULT_BREGMA_AP_MM, DEFAULT_BREGMA_DV_MM
    try:
        on_disk = json.loads(Path(settings_file).read_text(encoding='utf-8'))
        if 0.0 < float(on_disk.get('bregma_ap_mm')) < 13.5:
            ap = float(on_disk['bregma_ap_mm'])
        if -5.0 < float(on_disk.get('bregma_dv_mm')) < 8.0:
            dv = float(on_disk['bregma_dv_mm'])
    except Exception:
        pass
    return ap, dv, BREGMA_ML_MM


def find_annotation_volume(location_file=ATLAS_LOCATION_FILE):
    """The 25 um CCF annotation volume in the ABC browser's atlas cache, or None."""
    roots = []
    try:
        roots.append(Path(json.loads(Path(location_file).read_text(encoding='utf-8'))['path']))
    except Exception:
        pass
    roots.append(DEFAULT_ATLAS_CACHE)
    for root in roots:
        found = sorted((root / 'image_volumes' / 'Allen-CCF-2020').glob('*/annotation_25.nii.gz'))
        if found:
            return found[-1]
    return None


def read_nifti_volume(path):
    """The voxel array of a gzipped NIfTI-1 file as [AP, DV, ML] (the file's own axis order)."""
    with gzip.open(path, 'rb') as handle:
        raw = handle.read()
    dims = struct.unpack('<8h', raw[40:56])
    shape = tuple(dims[1:1 + dims[0]])
    code, = struct.unpack('<h', raw[70:72])
    dtype = {2: '<u1', 4: '<i2', 8: '<i4', 16: '<f4', 64: '<f8', 256: '<i1', 512: '<u2', 768: '<u4'}[code]
    offset = int(struct.unpack('<f', raw[108:112])[0])
    return np.frombuffer(raw, dtype=dtype, count=int(np.prod(shape)), offset=offset).reshape(shape, order='F')


class BrainTop:
    """The depth (CCF DV, mm) of the top of the brain at a CCF (AP, ML), from the annotation volume."""

    def __init__(self, volume):
        inside = volume != 0
        any_tissue = inside.any(axis=1)
        self.top = np.where(any_tissue, inside.argmax(axis=1), -1)   # [AP, ML] -> first DV voxel in the brain

    def __call__(self, ccf_ap, ccf_ml):
        i, k = int(round(ccf_ap / VOXEL_MM)), int(round(ccf_ml / VOXEL_MM))
        if not (0 <= i < self.top.shape[0] and 0 <= k < self.top.shape[1]) or self.top[i, k] < 0:
            return None
        return self.top[i, k] * VOXEL_MM


SKULL_SPAN_MM = 4.0       # dips in the brain's top narrower than about twice this are bridged when the skull's underside is estimated (see SkullSurface)
SKULL_STEP = 4            # the estimate is made on a grid of this many voxels (0.1 mm)


def _window_reduce(grid, offsets, reduce, fill, radius):
    """`reduce` (np.maximum or np.minimum) of `grid` shifted by each of `offsets` (a disk of cells): a dilation or an erosion; `fill` beyond the edges."""
    rows, cols = grid.shape
    padded = np.pad(grid, radius, mode='constant', constant_values=fill)
    out = padded[radius:radius + rows, radius:radius + cols].copy()
    for di, dj in offsets:
        reduce(out, padded[radius + di:radius + di + rows, radius + dj:radius + dj + cols], out=out)
    return out


class SkullSurface:
    """The depth (CCF DV, mm) of the underside of the skull, where the surgeon zeroed DV on the dura or brain surface, estimated from the top of the brain.
    The brain's top is not always right under the skull: behind the visual cortex, in front of the cerebellum, it dips about a millimetre while the skull does
    not. The estimate is the top of the brain with such dips bridged (a morphological closing of the height map by a flat disk of radius `span_mm`: dips
    narrower than that are filled to the height of their surroundings, the rest of the surface is left as it is). It is never lower than the brain's top.
    Call it like BrainTop. sheet() gives the surface for drawing."""

    def __init__(self, brain_top, span_mm=SKULL_SPAN_MM, step=SKULL_STEP):
        self.brain_top = brain_top
        self.step_mm = step * VOXEL_MM
        coarse = brain_top.top[::step, ::step]
        valid = coarse >= 0
        height = np.where(valid, -coarse * VOXEL_MM, -1e6)                       # (up is larger; outside the brain is far down)
        radius = max(1, int(round(span_mm / self.step_mm)))
        offsets = [(di, dj) for di in range(-radius, radius + 1) for dj in range(-radius, radius + 1) if di * di + dj * dj <= radius * radius]
        closed = _window_reduce(_window_reduce(height, offsets, np.maximum, -1e6, radius), offsets, np.minimum, 1e6, radius)
        self.depth = np.where(valid, -np.maximum(closed, height), np.nan)         # (CCF DV: smaller is higher)

    def __call__(self, ccf_ap, ccf_ml):
        raw = self.brain_top(ccf_ap, ccf_ml)
        if raw is None:
            return None
        i, k = int(round(ccf_ap / self.step_mm)), int(round(ccf_ml / self.step_mm))
        if 0 <= i < self.depth.shape[0] and 0 <= k < self.depth.shape[1] and np.isfinite(self.depth[i, k]):
            return min(raw, float(self.depth[i, k]))
        return raw

    def sheet(self, every=2):
        """The estimated surface as {'x': ML values, 'y': AP values, 'z': depth rows (None outside the brain)} for a Plotly surface, every `every` grid cells."""
        depth = self.depth[::every, ::every]
        step = self.step_mm * every
        z = [[None if not np.isfinite(v) else round(float(v), 3) for v in row] for row in depth]
        return {'x': [round(k * step, 3) for k in range(depth.shape[1])], 'y': [round(i * step, 3) for i in range(depth.shape[0])], 'z': z}


# ---------- the surgery side ----------
def ap_shift_mm(reference):
    """How far behind bregma the row's AP zero is (mm), from text like 'adjusted intersection point (~0.4 mm behind
    Paxinos bregma)'; 0 when it is bregma itself or not given."""
    m = re.search(r'([\d.]+)\s*mm\s*(behind|posterior|in front|anterior)', str(reference or ''), re.I)
    if not m:
        return 0.0
    return float(m.group(1)) * (1 if m.group(2).lower() in ('behind', 'posterior') else -1)


def dv_zero_kind(text):
    """'surface' for skull/dura/brain-surface zeroes, 'bregma' for bregma, None when not given."""
    t = str(text or '').strip().lower()
    if not t:
        return None
    if 'same' in t or re.search(r'\bbr(e?g(ma)?)?\b|\bbrg\b', t):
        return 'bregma'     # ('Dura/Bregma (same at this site)': the two heights agree)
    if re.search(r'dura|surf|skull|brain', t):
        return 'surface'
    return None


BREGMA_DV_PEOPLE = {'seung-chan'}    # people who usually zero DV at bregma: a blank zero in their rows means bregma


def default_dv_zero(person):
    """The DV zero of a row that gives none: the brain surface (the lab's usual), or bregma for BREGMA_DV_PEOPLE."""
    return 'bregma' if str(person or '').strip().lower() in BREGMA_DV_PEOPLE else 'surface'


def sides_of(side, ml_text=''):
    s = str(side or '').strip().lower()
    if s.startswith('bi') or re.search(r'±|\+/-', str(ml_text or '')):
        return ['Left', 'Right'], ''
    if s.startswith('r'):
        return ['Right'], ''
    if s.startswith('l'):
        return ['Left'], ''
    return ['Left'], 'side not given: drawn on the left'


def cell_values(text):
    """The numbers in a coordinate cell that may hold several sites ('-4.15\\n-4.7', '1.5/0.98', '4.5, 4.85',
    '(a) -3.1 (b) -4.23', '(+)1.4'), in order."""
    t = re.sub(r'\([a-z](\s*-\s*[a-z])?\)|\(\+\)', ' ', str(text or ''), flags=re.I)
    return [float(v) for v in re.findall(r'[-+]?\d*\.?\d+', t)]


def split_sites(r):
    """The row's sites: [(ap, ml, dv, angle or None, side or None, target, note)], one per line of a multi-site cell.
    A cell with one value is shared by every site (two depths on one track: one AP and ML, two DVs). An ML cell
    mirrored by sign ('1.63, -1.63') gives each site its own side (negative: left)."""
    aps, mls, dvs = cell_values(r.get('ap')), cell_values(r.get('ml')), cell_values(r.get('dv'))
    angles = cell_values(r.get('angle'))
    n = max(len(aps), len(mls), len(dvs))
    if n <= 1 or any(len(v) not in (1, n) for v in (aps, mls, dvs)):
        # one site, or counts that do not line up: the viewer's own single numbers
        return [(float(r['ap_mm']), float(r['ml_mm']), float(r['dv_mm']), None, None, r.get('target') or '',
                 '' if n <= 1 else 'several values that do not line up: drawn at the first')]
    targets = [t.strip() for t in str(r.get('target') or '').split('\n')]
    mirrored = len(mls) == n and min(mls) < 0 < max(mls)
    pick = lambda values, m: values[m] if len(values) == n else values[0]
    out = []
    for m in range(n):
        ml = pick(mls, m)
        out.append((pick(aps, m), abs(ml), -abs(pick(dvs, m)), angles[m] if len(angles) == n else None,
                    ('Left' if ml < 0 else 'Right') if mirrored else None,
                    targets[m] if len(targets) == n else (r.get('target') or ''), f"site {m + 1} of {n} in this row"))
    return out


def injection_sites(rows, bregma, brain_top):
    """[(x, y, z, hover, entry x, entry z, kind, zero)] in the page's axes (CCF ML, AP, DV in mm; see site_balls; kind
    'injection' or 'fiber', see site_kinds; zero 'surface' or 'bregma': what the DV was measured from) and a dict of what was left out."""
    sites, skipped = [], {'rat': 0, 'no coordinates': 0, 'other species': 0}
    for r in rows:
        species = (r.get('species') or '').strip().lower()
        if species == 'rat':
            skipped['rat'] += 1
            continue
        if species not in ('', 'mouse'):
            skipped['other species'] += 1
            continue
        try:
            row_sites = split_sites(r)
        except (KeyError, TypeError, ValueError):
            skipped['no coordinates'] += 1
            continue
        kinds = site_kinds(r, row_sites)
        for site, (kind, note) in zip(row_sites, kinds):
            if note:
                site = site[:6] + ('; '.join(n for n in (site[6], note) if n),)
            sites.extend(b[:6] + (kind,) + b[6:] for b in site_balls(r, *site, bregma=bregma, brain_top=brain_top))   # (the loader reads the kind at index 6)
    return sites, skipped


FIBER_SUBSTANCE = re.compile(r'optic(al)?\s*fib(er|re)|fib(er|re)\s*optic', re.I)


def site_kinds(r, row_sites):
    """('fiber' or 'injection', note) for each site of a row. A fiber: the substance says optic fiber, or the
    procedure is a fiber placement. An "Infusion + Fiber" row with one AP and ML and two depths (Nick's '-4.3/-4.1')
    is the injection and the fiber above it: the shallower site is the fiber. (When its AP or ML differ too, which
    site is which can't be told, and both are drawn as injections.)"""
    procedure = str(r.get('procedure') or '')
    if FIBER_SUBSTANCE.search(str(r.get('substance') or '')) or re.search(r'fib(er|re)\s*placement', procedure, re.I):
        return [('fiber', '')] * len(row_sites)
    if (re.search(r'\+\s*fib(er|re)', procedure, re.I) and len(row_sites) == 2
            and len(cell_values(r.get('ap'))) <= 1 and len(cell_values(r.get('ml'))) <= 1):
        fiber = min(range(2), key=lambda m: abs(row_sites[m][2]))
        return [('fiber', 'the fiber: the shallower depth of an Infusion + Fiber row') if m == fiber else
                ('injection', '') for m in range(2)]
    return [('injection', '')] * len(row_sites)


def one_line(text):
    return re.sub(r'\s*\n\s*', ' / ', str(text or '').strip())


def site_balls(r, ap, ml, dv, site_angle, site_side, target, site_note, bregma, brain_top):
    """The ball(s) of one site of a row, (x, y, z, hover, entry x, entry z, zero) each ('surface' or 'bregma': what the DV was measured from): one, or one per side when bilateral.
    x, y, z is the tip; the entry point (same y: tracks tilt only in the coronal plane) is where the DV was zeroed,
    on the brain surface or at bregma's height."""
    bregma_ap, bregma_dv, midline = bregma
    notes = [site_note] if site_note else []
    shift = ap_shift_mm(r.get('ap_reference'))
    ap_b = ap - shift
    if shift:
        notes.append(f"AP zero {shift:g} mm behind bregma: AP {ap_b:+.2f} from bregma")
    zero = dv_zero_kind(r.get('dv_zero'))
    if zero is None:
        zero = default_dv_zero(r.get('person'))
        notes.append(f"DV zero not given: taken as {'bregma' if zero == 'bregma' else 'the brain surface'}"
                     f"{' (as this person usually writes)' if zero == 'bregma' else ''}")
    angles = [site_angle] if site_angle is not None else (cell_values(r.get('angle')) or [0.0])
    if site_side:
        sides = [site_side]
    else:
        sides, side_note = sides_of(r.get('side'), r.get('ml'))
        if side_note:
            notes.append(side_note)
    ccf_ap = bregma_ap - ap_b
    depth = abs(dv)
    balls = []
    for n, side in enumerate(sides):
        angle = angles[min(n, len(angles) - 1)]
        sign = 1 if side == 'Right' else -1          # CCF ML (the page's left-right axis) grows toward the right
        theta = math.radians(abs(angle))
        entry_ml = midline + sign * abs(ml)
        tip_ml = midline + sign * (abs(ml) - depth * math.sin(theta))   # tilted toward the midline
        down = depth * math.cos(theta)
        side_notes = []
        if zero == 'surface':
            top = brain_top(ccf_ap, entry_ml) if brain_top else None
            if top is None:
                side_notes.append('no CCF brain at the entry point: DV taken from bregma')
                top = bregma_dv
        else:
            top = bregma_dv                          # DV zeroed at bregma: the track starts at bregma's height
        ccf_dv = top + down
        if angle:
            side_notes.append(f"{angle:g} deg: enters at ML {abs(ml):g}, tip at ML {abs(ml) - depth * math.sin(theta):.2f}")
        head = f"<b>{r.get('animal_id') or '?'}</b>  {r.get('person', '')}  {r.get('surgery_date', '')}"
        what = ' · '.join(one_line(v) for v in (target, r.get('substance'), r.get('procedure')) if v)
        coords = (f"AP {ap:g}  ML {abs(ml):g}  DV {dv:g}{'  angle ' + format(angle, 'g') if angle else ''}  ({side}"
                  f"{', DV from ' + one_line(r.get('dv_zero')) if str(r.get('dv_zero') or '').strip() else ''})")
        balls.append((round(tip_ml, 3), round(ccf_ap, 3), round(float(ccf_dv), 3),
                      '<br>'.join([head, what, coords] + notes + side_notes),
                      round(entry_ml, 3), round(float(top), 3), zero))
    return balls


def read_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


# ---------- the page ----------
def unit_sphere(subdivisions=2):
    """Vertices and faces of an icosphere of radius 1."""
    t = (1 + 5 ** 0.5) / 2
    verts = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
             (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6),
             (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10),
             (8, 6, 7), (9, 8, 1)]
    verts = [np.array(v, dtype=float) / np.linalg.norm(v) for v in verts]
    for _ in range(subdivisions):
        cache, new_faces = {}, []

        def middle(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                v = verts[a] + verts[b]
                verts.append(v / np.linalg.norm(v))
                cache[key] = len(verts) - 1
            return cache[key]
        for a, b, c in faces:
            ab, bc, ca = middle(a, b), middle(b, c), middle(c, a)
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = new_faces
    return np.round(np.array(verts), 4), np.array(faces)


class PageNeedsUpgrade(ValueError):
    """The page was written before 3D pages carried the loader that shows injection sites: compressed_page.upgrade_page adds it."""


def side_file_data(sites, radius, label):
    verts, faces = unit_sphere(2)
    return {'label': label, 'sites': sites, 'groups': [
                {'kind': 'injection', 'name': "injection sites", 'trackName': "needle tracks", 'color': BALL_COLOR,
                 'radius': radius, 'width': TRACK_WIDTH_PX},
                {'kind': 'fiber', 'name': "optic fibers", 'trackName': "fiber tracks", 'color': FIBER_COLOR,
                 'radius': radius * FIBER_RADIUS_SCALE, 'width': FIBER_WIDTH_PX}],
            'sphere': {'v': verts.ravel().tolist(), 'f': faces.ravel().tolist()}}


def clean_side_files(folder, keep=KEEP_SIDE_FILES):
    """Deletes all but the newest `keep` side files in `folder` (a tab open on a deleted one shows a note that it is gone)."""
    files = sorted(Path(folder).glob('injections-*.js'), key=lambda f: f.stat().st_mtime, reverse=True)
    for old in files[keep:] + sorted(Path(folder).glob('injections-open-*.html'), key=lambda f: f.stat().st_mtime, reverse=True)[keep:]:
        try:
            old.unlink()
        except OSError:
            pass


LANDMARK_COLOR = '#e01818'
LANDMARK_CUBE_MM = 0.15       # the cubes' edge: half the width of an injection ball
LAMBDA_BEHIND_BREGMA_MM = 4.2  # lambda's place: a typical adult mouse (approximate; the head is leveled so that lambda is at bregma's height)

# Draws the red cubes at bregma and lambda and, when there is one, the estimated skull sheet (the page's axes: x = ML, y = AP, z = DV). It is part of the side
# file, not of the page's loader, so that pages that already have the loader show it too. Waits for the plot the same way the loader does.
EXTRAS_JS = r"""
(function (e) {
  function draw() {
    var gd = document.querySelector('.plotly-graph-div');
    if (!window.Plotly || !gd || !gd._fullLayout || !gd._fullLayout.scene) return setTimeout(draw, 200);
    var traces = [], have = {}, firstOwn = gd._fullData.length;
    gd._fullData.forEach(function (t) { have[t.name] = true; });   // (an export from the ABC browser may already have the cubes, and its own spacer)
    var haveCubes = gd._fullData.some(function (t) { return /^(bregma|lambda)/.test(t.name || ''); });
    if (e.cubes.length && !haveCubes) {   // both cubes are ONE object (one legend entry, switched on and off together)
      var x = [], y = [], z = [], fi = [], fj = [], fk = [], text = [];
      e.cubes.forEach(function (b, n) {
        var h = b.edge / 2;
        [[-1, -1, -1], [-1, 1, -1], [1, 1, -1], [1, -1, -1], [-1, -1, 1], [-1, 1, 1], [1, 1, 1], [1, -1, 1]].forEach(function (c) {
          x.push(b.x + c[0] * h); y.push(b.y + c[1] * h); z.push(b.z + c[2] * h); text.push(b.text);
        });
        [7, 0, 0, 0, 4, 4, 2, 6, 4, 0, 3, 7].forEach(function (v, q) {
          fi.push(v + 8 * n); fj.push([3, 4, 1, 2, 5, 6, 5, 5, 0, 1, 2, 2][q] + 8 * n); fk.push([0, 7, 2, 3, 6, 7, 1, 2, 5, 5, 7, 6][q] + 8 * n);
        });
      });
      traces.push({type: 'mesh3d', x: x, y: y, z: z, i: fi, j: fj, k: fk, color: e.cubes[0].color, opacity: 1, flatshading: true,
        lighting: {ambient: 0.6, diffuse: 0.8, specular: 0.2, roughness: 0.6},
        name: e.cubesName, showlegend: true, text: text, hoverinfo: 'text', legendrank: 1});
    }
    if (e.discs && e.discs.points.length) {   // a small flat disc where each skull-referenced track starts (all in one mesh)
      var dx = [], dy = [], dz = [], di = [], dj = [], dk = [], dt = [], sides = 16, rad = e.discs.diameter / 2, half = e.discs.thickness / 2;
      function face(a, b, c) { di.push(a); dj.push(b); dk.push(c); }
      e.discs.points.forEach(function (pt) {
        var top = dx.length, bottom = top + sides + 1;   // vertices: the top's centre then its ring, the bottom's centre then its ring
        (half > 0 ? [-half, half] : [0]).forEach(function (dzv) {   // (no thickness: one flat face, drawn from both sides)
          dx.push(pt[0]); dy.push(pt[1]); dz.push(pt[2] + dzv); dt.push(e.discs.text);
          for (var s = 0; s < sides; s++) {
            var ang = 2 * Math.PI * s / sides;
            dx.push(pt[0] + rad * Math.cos(ang)); dy.push(pt[1] + rad * Math.sin(ang)); dz.push(pt[2] + dzv); dt.push(e.discs.text);
          }
        });
        for (var m = 0; m < sides; m++) {
          var next = (m + 1) % sides;
          face(top, top + 1 + m, top + 1 + next);                                  // the top
          if (half > 0) {
            face(bottom, bottom + 1 + next, bottom + 1 + m);                       // the bottom
            face(top + 1 + m, bottom + 1 + m, bottom + 1 + next);                 // the rim
            face(top + 1 + m, bottom + 1 + next, top + 1 + next);
          }
        }
      });
      traces.push({type: 'mesh3d', x: dx, y: dy, z: dz, i: di, j: dj, k: dk, color: e.discs.color, opacity: 1, flatshading: true,
        lighting: {ambient: 0.7, diffuse: 0.7, specular: 0.1, roughness: 0.8}, name: e.discs.name, showlegend: true, text: dt, hoverinfo: 'text'});
    }
    if (e.skull) {   // the sheet as a mesh (so that the slice cut can trim it like the other surfaces): two triangles per grid cell, none over gaps
      var s = e.skull, id = [], vx = [], vy = [], vz = [], fi = [], fj = [], fk = [];
      for (var r = 0; r < s.y.length; r++) {
        id.push([]);
        for (var c = 0; c < s.x.length; c++) {
          if (s.z[r][c] === null) { id[r].push(-1); continue; }
          id[r].push(vx.length); vx.push(s.x[c]); vy.push(s.y[r]); vz.push(s.z[r][c]);
        }
      }
      function tri(a, b, d) { if (a >= 0 && b >= 0 && d >= 0) { fi.push(a); fj.push(b); fk.push(d); } }
      for (var r2 = 0; r2 < s.y.length - 1; r2++) {
        for (var c2 = 0; c2 < s.x.length - 1; c2++) {
          var p = id[r2][c2], q = id[r2][c2 + 1], u = id[r2 + 1][c2], w = id[r2 + 1][c2 + 1];
          tri(p, q, w); tri(p, w, u);
        }
      }
      traces.push({type: 'mesh3d', x: vx, y: vy, z: vz, i: fi, j: fj, k: fk, color: e.skullColor, opacity: 0.3, flatshading: false,
        lighting: {ambient: 0.8, diffuse: 0.3, specular: 0, roughness: 1}, name: 'estimated skull underside', showlegend: true, hoverinfo: 'skip',
        legendrank: 3});
    }
    // bregma, lambda and the skull estimate are first in the legend (rank), then one blank line (a blank-named stand-in trace, which the page hides), as in
    // the page's own legend; the page's other entries have the default rank, 1000
    var spacer = '​​​';
    if (traces.length && !have[spacer]) {
      traces.push({type: 'scatter3d', mode: 'markers', x: [null], y: [null], z: [null], name: spacer, showlegend: true, hoverinfo: 'skip',
                   marker: {size: 1, color: 'rgba(0,0,0,0)'}, legendrank: 4});
    }
    Plotly.addTraces(gd, traces).then(function () { makeSliceable(gd, firstOwn); groupLegend(gd, 40); fixSwatches(gd); installSliders(); });
  }

  // The slice sliders' drag, and the diagnostics when they are on (Settings > Troubleshooting): the same as the page's own, for pages written before them. Each
  // part does nothing when the page already has it.
  function installSliders() {__SLIDER_SCRIPTS__
  }

  // Pages written before the page's own legend script cleared it can leave an entry without its icon: plotly re-uses a legend row for another entry when the
  // legend is rebuilt (as it is when these objects are added), and a row that was a surface's still has its swatch hidden by that script. Unhide the swatch of
  // every entry that is not a surface, edges or envelope (those are drawn by the script itself), after each redraw.
  function fixSwatches(gd) {
    function fix() {
      Array.prototype.forEach.call(gd.querySelectorAll('.legend g.traces'), function (group) {
        var data = group.__data__, trace = data && data[0] && data[0].trace;
        if (!trace || / (surface|edges|envelope)$/.test(trace.name || '')) return;
        var layers = group.querySelector('.layers');
        if (layers && layers.style.display === 'none') layers.style.removeProperty('display');
      });
    }
    var waiting = 0;   // (one check after a burst of redraws, not one for each)
    function soon() { clearTimeout(waiting); waiting = setTimeout(fix, 60); }
    gd.on('plotly_afterplot', soon);
    gd.on('plotly_restyle', soon);
    setTimeout(fix, 300);
  }

  // The injection sites, their tracks, the optic fibers and the skull discs are one legend group, with a small space above it. A legend group's title row has
  // a font of its own, so a blank title in a tiny font is the space. The injection traces are drawn
  // by the page's loader, perhaps a moment after these, so this looks again until they are there. Clicking an entry still switches only that entry.
  function groupLegend(gd, tries) {
    var groups = [['injections', /^(injection sites|needle tracks|optic fibers|fiber tracks|DV from the skull)/]];
    var found = groups.map(function (g) {
      var ids = [];
      gd._fullData.forEach(function (t, index) { if (g[1].test(t.name || '')) ids.push(index); });
      return ids;
    });
    if (!found[0].length && tries > 0) return setTimeout(function () { groupLegend(gd, tries - 1); }, 250);
    try {
      Plotly.relayout(gd, {'legend.groupclick': 'toggleitem'});
      groups.forEach(function (g, n) {
        if (found[n].length) Plotly.restyle(gd, {legendgroup: g[0], 'legendgrouptitle.text': ' ', 'legendgrouptitle.font.size': 4}, found[n]);
      });
    } catch (err) { /* the legend stays as it was */ }
  }

  // The page's Slice checkbox and sliders cut the objects that were in the page when it opened. These were added later (as were the injection sites the
  // loader drew), so they are cut here, with the same rule (a face is kept when its centre is inside the box; a point of a line is blanked outside it).
  // The page's slice script calls gd.__rebuildCells() after every slider move and leaves the box in gd.__sliceCut (absent when the checkbox is off).
  function makeSliceable(gd, firstOwn) {
    if (gd.__extrasSliceable) return;
    gd.__extrasSliceable = true;
    var mine = /^(injection sites|optic fibers|needle tracks|fiber tracks|DV from the skull|bregma|lambda|estimated skull)/;
    var known = {};
    var previous = gd.__rebuildCells;
    // The page calls this after every slider move, and does not wait for what it starts, and each update of the plot costs about as much as one of the page's
    // own (a fraction of a second with the brain's shells). So the cut here waits until the slider has paused for a moment (during a fast drag only the page's own
    // objects follow; these catch up when the thumb slows or stops), is done one pass at a time, and takes at most two updates of the plot per pass, one for all
    // the meshes and one for all the lines.
    var working = false, again = false, timer = 0;
    gd.__rebuildCells = function () { if (previous) previous(); clearTimeout(timer); timer = setTimeout(schedule, 150); };
    function schedule() {
      if (working) { again = true; return; }
      working = true;
      var jobs = [];
      try { jobs = cutExtras(); } catch (err) { /* the page's own cut still works */ }
      function done() { working = false; if (again) { again = false; schedule(); } }
      Promise.all(jobs).then(done, done);
    }
    function cutExtras() {
      var cut = gd.__sliceCut;
      var meshes = {i: [], j: [], k: []}, meshIds = [], lines = {x: [], y: [], z: [], text: []}, lineIds = [];
      gd._fullData.forEach(function (t, index) {
        if (!mine.test(t.name || '') || (t.type !== 'mesh3d' && t.type !== 'scatter3d')) return;
        if (index < firstOwn && /^(bregma|lambda|estimated skull)/.test(t.name)) return;   // (the page's own: its slice script cuts those)
        var r = known[index];
        if (!r) {
          r = known[index] = {mesh: t.type === 'mesh3d', x: Array.from(t.x), y: Array.from(t.y), z: Array.from(t.z), modified: false, last: -1,
                              text: t.text && t.type === 'scatter3d' ? Array.from(t.text) : null};
          if (r.mesh) { r.i = Int32Array.from(t.i); r.j = Int32Array.from(t.j); r.k = Int32Array.from(t.k); }
        }
        function putMesh(ni, nj, nk) { meshes.i.push(ni); meshes.j.push(nj); meshes.k.push(nk); meshIds.push(index); }
        function putLine(nx, ny, nz, nt) { lines.x.push(nx); lines.y.push(ny); lines.z.push(nz); lines.text.push(nt); lineIds.push(index); }
        if (!cut) {   // the checkbox is off: put back what was cut
          if (!r.modified) return;
          r.modified = false; r.last = -1;
          if (r.mesh) putMesh(r.i, r.j, r.k); else putLine(r.x, r.y, r.z, r.text);
          return;
        }
        if (t.visible !== true) { r.last = -1; return; }   // not drawn: brought up to date when it is switched on
        function inside(x, y, z) {
          return x >= cut.lo[0] && x <= cut.hi[0] && y >= cut.lo[1] && y <= cut.hi[1] && z >= cut.lo[2] && z <= cut.hi[2];
        }
        if (r.mesh) {
          var keep = [];
          for (var f = 0; f < r.i.length; f++) {
            var a = r.i[f], b = r.j[f], d = r.k[f];
            if (inside((r.x[a] + r.x[b] + r.x[d]) / 3, (r.y[a] + r.y[b] + r.y[d]) / 3, (r.z[a] + r.z[b] + r.z[d]) / 3)) keep.push(f);
          }
          if (keep.length === r.i.length && !r.modified) { r.last = keep.length; return; }
          if (keep.length === r.last) return;
          r.last = keep.length; r.modified = true;
          var ni = new Int32Array(keep.length), nj = new Int32Array(keep.length), nk = new Int32Array(keep.length);
          for (var m = 0; m < keep.length; m++) { ni[m] = r.i[keep[m]]; nj[m] = r.j[keep[m]]; nk[m] = r.k[keep[m]]; }
          putMesh(ni, nj, nk);
        } else {
          var nx = [], ny = [], nz = [], nt = r.text ? [] : null, out = 0;
          for (var q = 0; q < r.x.length; q++) {
            var gap = r.x[q] === null || r.y[q] === null || r.z[q] === null;
            var away = !gap && !inside(r.x[q], r.y[q], r.z[q]);
            if (away) out++;
            nx.push(gap || away ? null : r.x[q]); ny.push(gap || away ? null : r.y[q]); nz.push(gap || away ? null : r.z[q]);
            if (nt) nt.push(r.text[q]);
          }
          if (out === 0 && !r.modified) { r.last = 0; return; }
          if (out === r.last) return;
          r.last = out; r.modified = true;
          putLine(nx, ny, nz, nt);
        }
      });
      var jobs = [];
      if (meshIds.length) jobs.push(Plotly.restyle(gd, meshes, meshIds));
      if (lineIds.length) jobs.push(Plotly.restyle(gd, lines, lineIds));
      return jobs;
    }
  }
  draw();
})(__EXTRAS__);
"""
SKULL_SHEET_COLOR = '#c9b27c'
SKULL_DISC_COLOR = '#6b3f12'     # darker than the sheet, so that both show
SKULL_DISC_THICKNESS_MM = 0      # a flat disc (0), drawn from both sides; a number makes it a thin cylinder


def landmark_cubes(bregma):
    """The cubes at bregma and lambda, for EXTRAS_JS, from bregma_in_ccf's (AP, DV, ML). In the CCF, AP grows backward, so lambda is behind bregma at a larger AP;
    it is put at bregma's height (the usual leveling of the head) and on the midline."""
    ap, dv, ml = bregma
    lam_ap = ap + LAMBDA_BEHIND_BREGMA_MM
    return [{'name': 'bregma', 'rank': 1, 'x': round(ml, 3), 'y': round(ap, 3), 'z': round(dv, 3), 'edge': LANDMARK_CUBE_MM, 'color': LANDMARK_COLOR,
             'text': f"<b>bregma</b><br>CCF: AP {ap:g}, DV {dv:g}, ML {ml:g} mm"},
            {'name': 'lambda', 'rank': 2, 'x': round(ml, 3), 'y': round(lam_ap, 3), 'z': round(dv, 3), 'edge': LANDMARK_CUBE_MM, 'color': LANDMARK_COLOR,
             'text': f"<b>lambda</b> (approximate)<br>{LAMBDA_BEHIND_BREGMA_MM:g} mm behind bregma, at bregma's height<br>CCF: AP {lam_ap:g}, DV {dv:g}, ML {ml:g} mm"}]


def write_side_file(sites, radius, label, folder=None, bregma=None, skull=None, diagnostics=False):
    """Writes the sites as a side file: a script that calls window.__injectionSites({...}), which the page's loader defines (see injection_loader.py), and
    (given `bregma`, as bregma_in_ccf gives it) draws the red cubes at bregma and lambda, and (given `skull`, a SkullSurface) the estimated skull sheet.
    Written through a temporary file, so that a page never reads half of it. A new name every time. Returns its path."""
    folder = Path(folder or SIDE_FOLDER)
    folder.mkdir(parents=True, exist_ok=True)
    stem = 'injections-' + time.strftime('%Y%m%d-%H%M%S')
    path, n = folder / (stem + '.js'), 1
    while path.exists():
        n += 1
        path = folder / f"{stem}-{n}.js"
    temporary = path.with_name(path.name + '.part')
    text = 'window.__injectionSites(' + json.dumps(side_file_data(sites, radius, label)).replace('</', '<\\/') + ');\n'
    if bregma or skull:
        tops = [[s[4], s[1], s[5]] for s in sites if len(s) > 7 and s[7] == 'surface']      # where the tracks zeroed on the skull, dura or brain surface start
        extras = {'cubes': landmark_cubes(bregma) if bregma else [], 'cubesName': 'bregma and lambda', 'skull': skull.sheet() if skull else None, 'skullColor': SKULL_SHEET_COLOR,
                  'discs': {'points': tops, 'diameter': LANDMARK_CUBE_MM, 'thickness': SKULL_DISC_THICKNESS_MM, 'color': SKULL_DISC_COLOR,
                            'name': 'DV from the skull or dura',
                            'text': "DV measured from the skull, dura or brain surface here, not from bregma (the estimated skull underside is used)"}}
        text += EXTRAS_JS.replace('__EXTRAS__', json.dumps(extras).replace('</', '<\\/')).replace('__SLIDER_SCRIPTS__', slice_slider_scripts.page_js(diagnostics))
    temporary.write_text(text, encoding='utf-8')
    os.replace(temporary, path)
    clean_side_files(folder)
    return path


def view_address(page, side_file):
    """The address that opens `page` showing the sites of `side_file`: the page's file address, then  #sites=  and the side file's, percent-encoded."""
    return Path(page).resolve().as_uri() + '#sites=' + urllib.parse.quote(Path(side_file).resolve().as_uri(), safe='')


def make_view(page, surgeries, radius=DEFAULT_RADIUS_MM, upgrade=False, skull_span=SKULL_SPAN_MM):
    """Write the side file for the surgery list. Returns (the address to open the page at, the lines to report); raises ValueError when there is nothing to draw or
    the page is not a 3D page, and PageNeedsUpgrade when the page has no loader (with upgrade=True the page is upgraded first). `skull_span` (mm): see SkullSurface."""
    page = Path(page)
    if not compressed_page.looks_like_3d_page(page):
        raise ValueError(f"{page.name} does not look like a 3D page from the ABC browser (no plotly in it).")
    needs_upgrade = not compressed_page.has_loader(page)
    if needs_upgrade and not upgrade:
        raise PageNeedsUpgrade(f"{page.name} was made before 3D pages could show injection sites: it needs a small one-time change (it will look the same).")
    bregma = bregma_in_ccf()
    volume_path = find_annotation_volume()
    brain_top = BrainTop(read_nifti_volume(volume_path)) if volume_path else None
    report = []
    if brain_top is None:
        report.append("The CCF annotation volume was not found in the atlas cache: surface-zeroed DVs are taken from bregma.")
    skull = SkullSurface(brain_top, span_mm=skull_span) if brain_top is not None else None   # (where DV is zeroed on the dura or brain surface)
    sites, skipped = injection_sites(read_csv(surgeries), bregma, skull)
    if not sites:
        raise ValueError("No mouse injection sites with coordinates in that surgery list.")
    if needs_upgrade:
        compressed_page.upgrade_page(page)
        report.append(f"{page.name} was changed so that it can show injection sites (it looks the same, and is stored compressed).")
    left_out = ', '.join(f"{n} {why}" for why, n in skipped.items() if n)
    fibers = sum(1 for s in sites if s[6] == 'fiber')
    label = f"{time.strftime('%H:%M')}, {len(sites) - fibers} injections" + (f", {fibers} fibers" if fibers else "")
    side = write_side_file(sites, radius, label, bregma=bregma, skull=skull, diagnostics=slice_slider_scripts.diagnostics_enabled(SETTINGS_FILE))
    report.insert(0, f"{len(sites) - fibers} injection sites and {fibers} optic fibers written to {side.name}"
                  + (f" (left out: {left_out})" if left_out else ""))
    report.append(f"Bregma in the CCF: AP {bregma[0]:g}, DV {bregma[1]:g}, ML {bregma[2]:g} mm")
    if skull is not None:
        report.append(f"DV zeroed on the dura or brain surface is measured from the estimated underside of the skull (dips in the brain's top narrower "
                      f"than about {2 * skull_span:g} mm are bridged); it is drawn as a translucent sheet.")
    return view_address(page, side), report


# ---------- the window (run-add-injections.bat / .command) ----------
LAST_USED_FILE = SRC / 'cache_local' / 'add_injections_3d.json'


def load_last_used():
    try:
        value = json.loads(LAST_USED_FILE.read_text(encoding='utf-8'))
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def save_last_used(values):
    try:
        LAST_USED_FILE.parent.mkdir(parents=True, exist_ok=True)
        LAST_USED_FILE.write_text(json.dumps(values, indent=1), encoding='utf-8')
    except OSError:
        pass


def newest_export(folder=Path.home() / 'Downloads'):
    """The surgery viewer's newest Export CSV in `folder` (jhoulab-surgeries-*.csv), or None."""
    found = sorted(Path(folder).glob('jhoulab-surgeries-*.csv'), key=lambda p: p.stat().st_mtime)
    return found[-1] if found else None


def default_surgeries(last_used):
    """The newest export in Downloads when it is newer than the last list used (a fresh export is the usual next
    step), else the last list used."""
    last = last_used.get('surgeries')
    newest = newest_export()
    try:
        if newest and (not last or not Path(last).exists() or newest.stat().st_mtime > Path(last).stat().st_mtime):
            return str(newest)
    except OSError:
        pass
    return last or ''


def bring_forward(root):
    """Puts the window in the middle of the monitor the pointer is on and in front of the others. Started by the surgery
    viewer's server (a background program), it would otherwise open behind the browser: Windows only lets the program in
    front raise a window, but a moment of "always on top" is allowed."""
    try:
        import screen_area
        screen_area.center_on_home(root, vertical=1 / 3)
    except Exception:
        pass
    root.lift()
    try:
        root.attributes('-topmost', True)
        root.after(400, lambda: root.attributes('-topmost', False))
    except Exception:
        pass
    root.focus_force()


def open_in_browser(address):
    """Opens `address` (a file address with a # part) in the default browser, in a new tab; whether that worked.
    Windows drops the # part of a file address when it launches the browser (the page then opens with no sites), so the browser is given a tiny launcher
    page (a plain file path, nothing to drop) that sends it on to `address`."""
    import html
    import webbrowser
    try:
        SIDE_FOLDER.mkdir(parents=True, exist_ok=True)
        launcher = SIDE_FOLDER / f"injections-open-{time.strftime('%Y%m%d-%H%M%S')}-{os.getpid()}.html"
        launcher.write_text(
            '<!doctype html><meta charset="utf-8"><title>Opening the 3D page...</title>'
            f'<script>location.replace({json.dumps(address)});</script>'
            f'<body style="font:14px sans-serif">Opening the 3D page... if nothing happens, <a href="{html.escape(address)}">click here</a>.</body>',
            encoding='utf-8')
        return bool(webbrowser.open(launcher.resolve().as_uri(), new=2))
    except Exception:
        return False


def open_window(surgeries=None):
    """The window. `surgeries` (a CSV path) fills in the surgery list, as the surgery viewer's Show in 3D button does;
    otherwise the newest export in Downloads or the last list used."""
    import threading
    import tkinter as tk
    from tkinter import filedialog, ttk
    last = load_last_used()
    root = tk.Tk()
    root.title("Injection sites in a 3D page")
    frame = ttk.Frame(root, padding=14)
    frame.pack(fill='both', expand=True)
    frame.columnconfigure(1, weight=1)
    page_var = tk.StringVar(root, last.get('page', ''))
    csv_var = tk.StringVar(root, str(surgeries) if surgeries else default_surgeries(last))
    radius_var = tk.StringVar(root, f"{last.get('radius', DEFAULT_RADIUS_MM):g}")
    span_var = tk.StringVar(root, f"{last.get('skull_span', SKULL_SPAN_MM):g}")
    status_var = tk.StringVar(root, "")

    def browse(var, title, types):
        current = Path(var.get()) if var.get().strip() else None
        path = filedialog.askopenfilename(parent=root, title=title, filetypes=types,
                                          initialdir=str(current.parent) if current and current.parent.exists() else None)
        if path:
            var.set(path)

    ttk.Label(frame, text="3D source page").grid(row=0, column=0, sticky='w', pady=4)
    ttk.Entry(frame, textvariable=page_var, width=70).grid(row=0, column=1, sticky='ew', padx=6)
    ttk.Button(frame, text="Browse...", command=lambda: browse(
        page_var, "Choose a 3D page from the ABC browser", [("3D page", "*.html"), ("All files", "*.*")])).grid(row=0, column=2)

    ttk.Label(frame, text="Surgery list").grid(row=2, column=0, sticky='w', pady=(10, 4))
    ttk.Entry(frame, textvariable=csv_var, width=70).grid(row=2, column=1, sticky='ew', padx=6, pady=(10, 4))
    ttk.Button(frame, text="Browse...", command=lambda: browse(
        csv_var, "Choose the surgery viewer's Export CSV", [("CSV file", "*.csv"), ("All files", "*.*")])).grid(row=2, column=2, pady=(10, 4))
    ttk.Label(frame, foreground='#555555', text=(
        "The surgery viewer's current list (from its Show in 3D button)." if surgeries else
        "In the surgery viewer, filter the list to the injections you want, then click Show in 3D\n"
        "(or Export CSV: the newest export in Downloads is filled in for you).")).grid(row=3, column=1, sticky='w', padx=6)

    row = ttk.Frame(frame)
    row.grid(row=4, column=1, sticky='w', padx=6, pady=(12, 0))
    ttk.Label(row, text="Ball radius (mm)").pack(side='left')
    ttk.Spinbox(row, textvariable=radius_var, from_=0.02, to=1.0, increment=0.05, width=6).pack(side='left', padx=6)
    ttk.Label(row, text="Skull estimate span (mm)").pack(side='left', padx=(18, 0))
    ttk.Spinbox(row, textvariable=span_var, from_=0.5, to=8.0, increment=0.5, width=6).pack(side='left', padx=6)

    buttons = ttk.Frame(frame)
    buttons.grid(row=5, column=0, columnspan=3, sticky='e', pady=(14, 0))
    make = ttk.Button(buttons, text="Show in 3D")
    make.pack(side='left')
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side='left', padx=(8, 0))
    status = ttk.Label(frame, textvariable=status_var, wraplength=620, justify='left')
    status.grid(row=6, column=0, columnspan=3, sticky='w', pady=(10, 0))

    def run(upgrade=False):
        page, surgeries = page_var.get().strip().strip('"'), csv_var.get().strip().strip('"')
        try:
            radius = float(radius_var.get())
            if not 0 < radius <= 2:
                raise ValueError
        except ValueError:
            status_var.set("Enter the ball radius as a number of mm (e.g. 0.15).")
            return
        try:
            span = float(span_var.get())
            if not 0.2 <= span <= 10:
                raise ValueError
        except ValueError:
            status_var.set("Enter the skull estimate span as a number of mm between 0.2 and 10 (e.g. 4).")
            return
        for label, path in (("3D source page", page), ("surgery list", surgeries)):
            if not path or not Path(path).is_file():
                status_var.set(f"Choose the {label} first (no file at '{path}')." if path else f"Choose the {label} first.")
                return
        make.state(['disabled'])
        root.config(cursor='watch')
        status_var.set("Working..." + (" (changing the page takes a few seconds)" if upgrade else ""))
        save_last_used({'page': page, 'surgeries': surgeries, 'radius': radius, 'skull_span': span})

        def work():
            try:
                address, report = make_view(page, surgeries, radius=radius, upgrade=upgrade, skull_span=span)
                result = ('ok', address, report)
            except PageNeedsUpgrade as e:
                result = ('upgrade', None, [str(e)])
            except Exception as e:
                result = ('error', None, [str(e)])
            root.after(0, lambda: finished(*result))
        threading.Thread(target=work, daemon=True).start()

    def finished(kind, address, report):
        make.state(['!disabled'])
        root.config(cursor='')
        if kind == 'upgrade':
            from tkinter import messagebox
            if messagebox.askyesno("Change the page?", report[0] + "\n\nThe file is rewritten in place (through a temporary file, so a failure "
                                   "leaves it as it was). Change it now?", parent=root):
                run(upgrade=True)
            else:
                status_var.set("Not changed: the page can not show injection sites until it is.")
        elif kind == 'ok':
            status_var.set('\n'.join(report) + "\nOpening it in a new tab.")
            if sys.stdout:                # (kept in the log once the window has closed; pythonw may have no stdout)
                print('\n'.join(report))
            if open_in_browser(address):
                root.after(800, root.destroy)   # the page is open: the window isn't needed any more
            else:
                status_var.set('\n'.join(report) + f"\nCould not open a browser: open this address yourself:\n{address}")
        else:
            status_var.set("Could not show the sites: " + report[0])

    make.config(command=run)
    root.bind('<Return>', lambda _e: run())
    root.after_idle(lambda: bring_forward(root))
    root.mainloop()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        return open_window()
    if argv[0] == '--window':      # the surgery viewer: --window --surgeries LIST.csv
        w = argparse.ArgumentParser(prog='add_injections_3d.py --window')
        w.add_argument('--surgeries', help="the surgery list to fill in")
        return open_window(w.parse_args(argv[1:]).surgeries)
    p = argparse.ArgumentParser(description="Show the surgery viewer's injection sites in an ABC browser 3D page "
                                            "(with no arguments: a window to choose the files).")
    p.add_argument('page', help="a 3D page (.html) written by the ABC atlas browser")
    p.add_argument('surgeries', help="the surgery viewer's Export CSV file")
    p.add_argument('--radius', type=float, default=DEFAULT_RADIUS_MM, help=f"ball radius in mm (default {DEFAULT_RADIUS_MM})")
    p.add_argument('--skull-span', type=float, default=SKULL_SPAN_MM, help=f"dips in the brain's top narrower than twice this (mm) are bridged when the skull's underside is estimated (default {SKULL_SPAN_MM:g})")
    p.add_argument('--no-open', action='store_true', help="only write the side file and print the address")
    a = p.parse_args(argv)
    try:
        address, report = make_view(a.page, a.surgeries, a.radius, upgrade=True, skull_span=a.skull_span)
    except ValueError as e:
        sys.exit(str(e))
    print('\n'.join(report))
    print(address)
    if not a.no_open and not open_in_browser(address):
        print("Could not open a browser: open the address above yourself.")
    return address


if __name__ == '__main__':
    main()
