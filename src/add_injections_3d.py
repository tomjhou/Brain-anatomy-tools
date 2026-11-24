"""Add the lab's injection sites to a 3D page written by the ABC atlas browser, as white balls.

    python add_injections_3d.py PAGE.html SURGERIES.csv [-o OUT.html] [--radius 0.15]

PAGE.html is a 3D page from the ABC browser (its axes are CCF millimetres). SURGERIES.csv is the
surgery viewer's "Export CSV" file (JhouLab-surgery-viewer), so filter the list there first (by
target, person, dates...) to choose which injections are drawn. The result is a copy of the page
(PAGE_injections.html unless -o is given) with one more legend entry, "injection sites", which
shows or hides the balls; hovering a ball shows the animal, target, substance and coordinates.
Running it again on its own output replaces the balls rather than adding a second set.

From stereotaxic coordinates to the CCF (all approximate; there is no exact mapping):
  * Bregma is placed in the CCF where the ABC browser's Settings put it (viewer_settings.json;
    by default the IBL's estimate), with no scaling or tilt, as the browser does.
  * AP: when the row's "AP measured from" says the zero is some mm behind bregma (the
    "adjusted intersection point (~0.4 mm behind Paxinos bregma)"), that is subtracted.
  * ML: the side column picks the hemisphere; bilateral rows give two balls. A row with no
    side is drawn on the left (the lab's usual side) and says so when hovered.
  * DV: zeroed at bregma (or with no zero given): measured down from bregma's height. Zeroed at
    the skull, dura or brain surface: measured down from the top of the CCF brain at the entry
    point (the CCF has no skull, so skull-zeroed sites land a little too deep).
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
import re
import struct
import sys
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parent
SETTINGS_FILE = SRC / 'viewer_settings.json'                              # the ABC browser's (bregma position)
ATLAS_LOCATION_FILE = SRC / 'cache_local' / 'abc_atlas_cache_location.json'
DEFAULT_ATLAS_CACHE = SRC.parent.parent / 'abc_atlas_cache'
DEFAULT_BREGMA_AP_MM, DEFAULT_BREGMA_DV_MM, BREGMA_ML_MM = 5.40, 0.332, 5.739   # as abc_core (the IBL's bregma)
VOXEL_MM = 0.025
DEFAULT_RADIUS_MM = 0.15
BALL_COLOR = '#ffffff'
START, END = '<!-- injection-sites:start -->', '<!-- injection-sites:end -->'


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
    """[(x, y, z, hover)] in the page's axes (CCF ML, AP, DV in mm) and a dict of what was left out."""
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
        for site in row_sites:
            sites.extend(site_balls(r, *site, bregma=bregma, brain_top=brain_top))
    return sites, skipped


def one_line(text):
    return re.sub(r'\s*\n\s*', ' / ', str(text or '').strip())


def site_balls(r, ap, ml, dv, site_angle, site_side, target, site_note, bregma, brain_top):
    """The ball(s) of one site of a row, (x, y, z, hover) each: one, or one per side when bilateral."""
    bregma_ap, bregma_dv, midline = bregma
    notes = [site_note] if site_note else []
    shift = ap_shift_mm(r.get('ap_reference'))
    ap_b = ap - shift
    if shift:
        notes.append(f"AP zero {shift:g} mm behind bregma: AP {ap_b:+.2f} from bregma")
    zero = dv_zero_kind(r.get('dv_zero'))
    if zero is None:
        notes.append('DV zero not given: taken as bregma')
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
            ccf_dv = top + down
        else:
            ccf_dv = bregma_dv + down
        if angle:
            side_notes.append(f"{angle:g} deg: enters at ML {abs(ml):g}, tip at ML {abs(ml) - depth * math.sin(theta):.2f}")
        head = f"<b>{r.get('animal_id') or '?'}</b>  {r.get('person', '')}  {r.get('surgery_date', '')}"
        what = ' · '.join(one_line(v) for v in (target, r.get('substance'), r.get('procedure')) if v)
        coords = (f"AP {ap:g}  ML {abs(ml):g}  DV {dv:g}{'  angle ' + format(angle, 'g') if angle else ''}  ({side}"
                  f"{', DV from ' + one_line(r.get('dv_zero')) if zero else ''})")
        balls.append((round(tip_ml, 3), round(ccf_ap, 3), round(float(ccf_dv), 3),
                      '<br>'.join([head, what, coords] + notes + side_notes)))
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


PAGE_SCRIPT = """
<script>
(function () {
  var data = __DATA__;
  function add() {
    var gd = document.querySelector('.plotly-graph-div');
    if (!window.Plotly || !gd || !gd._fullLayout || !gd._fullLayout.scene) return setTimeout(add, 200);
    var uv = data.sphere.v, uf = data.sphere.f, n = uv.length / 3, r = data.radius;
    var x = [], y = [], z = [], i = [], j = [], k = [], text = [];
    data.sites.forEach(function (s, b) {
      for (var v = 0; v < n; v++) {
        x.push(s[0] + r * uv[3 * v]); y.push(s[1] + r * uv[3 * v + 1]); z.push(s[2] + r * uv[3 * v + 2]);
        text.push(s[3]);
      }
      for (var f = 0; f < uf.length; f += 3) { i.push(b * n + uf[f]); j.push(b * n + uf[f + 1]); k.push(b * n + uf[f + 2]); }
    });
    Plotly.addTraces(gd, {type: 'mesh3d', x: x, y: y, z: z, i: i, j: j, k: k, color: data.color, opacity: 1,
                          flatshading: false, lighting: {ambient: 0.55, diffuse: 0.8, specular: 0.3, roughness: 0.5},
                          name: data.name, showlegend: true, text: text, hoverinfo: 'text'});
    if (!gd.__hoverPassThrough) passThrough(gd);   // (newer pages do this themselves)
  }
  // Hover through the background, for pages written before the ABC browser did it (the same code as its
  // THREED_WEB_LEGEND_ROUTER_SCRIPT): objects whose trace skips hover are left out of plotly's pick pass, so the
  // translucent brain surface no longer hides the balls and cells inside it from the pointer.
  function passThrough(gd) {
    gd.__hoverPassThrough = function () {
      var scene = gd._fullLayout && gd._fullLayout.scene && gd._fullLayout.scene._scene;
      if (!scene || !scene.traces || !scene.glplot) return;
      var changed = false;
      Object.keys(scene.traces).forEach(function (uid) {
        var trace = scene.traces[uid];
        if (!trace || !trace.data || trace.data.hoverinfo !== 'skip') return;
        Object.keys(trace).forEach(function (key) {
          var object = trace[key];
          if (object && typeof object === 'object' && typeof object.drawPick === 'function') { object.drawPick = null; changed = true; }
        });
      });
      if (changed) scene.glplot.update();
    };
    ['plotly_afterplot', 'plotly_restyle'].forEach(function (name) {
      gd.on(name, function () { setTimeout(gd.__hoverPassThrough, 0); });
    });
    gd.__hoverPassThrough();
  }
  add();
})();
</script>
"""


def add_to_page(page_text, sites, radius, color=BALL_COLOR):
    """The page's HTML with the balls' script in it (replacing an earlier one)."""
    page_text = re.sub(re.escape(START) + r'.*?' + re.escape(END) + r'\n?', '', page_text, flags=re.S)
    verts, faces = unit_sphere(2)
    data = {'sites': sites, 'radius': radius, 'color': color, 'name': f"injection sites ({len(sites)})",
            'sphere': {'v': verts.ravel().tolist(), 'f': faces.ravel().tolist()}}
    script = START + PAGE_SCRIPT.replace('__DATA__', json.dumps(data).replace('</', '<\\/')) + END + '\n'
    at = page_text.rfind('</body>')
    if at < 0:
        return page_text + script
    return page_text[:at] + script + page_text[at:]


def output_path(page):
    """Where the page with the balls goes: PAGE_injections.html beside PAGE (PAGE itself if it already is one)."""
    page = Path(page)
    return page.with_name(re.sub(r'_injections$', '', page.stem) + '_injections.html')


def make_page(page, surgeries, out=None, radius=DEFAULT_RADIUS_MM):
    """Write the page with the balls. Returns (out path, the lines to report); raises ValueError when there is nothing
    to draw or the page is not a 3D page."""
    page = Path(page)
    out = Path(out) if out else output_path(page)
    text = page.read_text(encoding='utf-8')
    if 'plotly' not in text[:20000].lower():
        raise ValueError(f"{page.name} does not look like a 3D page from the ABC browser (no plotly in it).")
    bregma = bregma_in_ccf()
    volume_path = find_annotation_volume()
    brain_top = BrainTop(read_nifti_volume(volume_path)) if volume_path else None
    report = []
    if brain_top is None:
        report.append("The CCF annotation volume was not found in the atlas cache: surface-zeroed DVs are taken from bregma.")
    sites, skipped = injection_sites(read_csv(surgeries), bregma, brain_top)
    if not sites:
        raise ValueError("No mouse injection sites with coordinates in that surgery list.")
    out.write_text(add_to_page(text, sites, radius), encoding='utf-8')
    left_out = ', '.join(f"{n} {why}" for why, n in skipped.items() if n)
    report.insert(0, f"{len(sites)} balls written to {out}" + (f" (left out: {left_out})" if left_out else ""))
    report.append(f"Bregma in the CCF: AP {bregma[0]:g}, DV {bregma[1]:g}, ML {bregma[2]:g} mm")
    return out, report


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


def open_window():
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
    csv_var = tk.StringVar(root, default_surgeries(last))
    radius_var = tk.StringVar(root, f"{last.get('radius', DEFAULT_RADIUS_MM):g}")
    status_var = tk.StringVar(root, "")
    out_var = tk.StringVar(root, "")

    def show_out(*_args):
        out_var.set(f"Writes {output_path(page_var.get()).name} beside the page." if page_var.get().strip() else "")
    page_var.trace_add('write', show_out)
    show_out()

    def browse(var, title, types):
        current = Path(var.get()) if var.get().strip() else None
        path = filedialog.askopenfilename(parent=root, title=title, filetypes=types,
                                          initialdir=str(current.parent) if current and current.parent.exists() else None)
        if path:
            var.set(path)

    ttk.Label(frame, text="3D page").grid(row=0, column=0, sticky='w', pady=4)
    ttk.Entry(frame, textvariable=page_var, width=70).grid(row=0, column=1, sticky='ew', padx=6)
    ttk.Button(frame, text="Browse...", command=lambda: browse(
        page_var, "Choose a 3D page from the ABC browser", [("3D page", "*.html"), ("All files", "*.*")])).grid(row=0, column=2)
    ttk.Label(frame, textvariable=out_var, foreground='#555555').grid(row=1, column=1, sticky='w', padx=6)

    ttk.Label(frame, text="Surgery list").grid(row=2, column=0, sticky='w', pady=(10, 4))
    ttk.Entry(frame, textvariable=csv_var, width=70).grid(row=2, column=1, sticky='ew', padx=6, pady=(10, 4))
    ttk.Button(frame, text="Browse...", command=lambda: browse(
        csv_var, "Choose the surgery viewer's Export CSV", [("CSV file", "*.csv"), ("All files", "*.*")])).grid(row=2, column=2, pady=(10, 4))
    ttk.Label(frame, foreground='#555555', text=(
        "In the surgery viewer, filter the list to the injections you want, then click Export CSV.\n"
        "The newest export in Downloads is filled in for you.")).grid(row=3, column=1, sticky='w', padx=6)

    row = ttk.Frame(frame)
    row.grid(row=4, column=1, sticky='w', padx=6, pady=(12, 0))
    ttk.Label(row, text="Ball radius (mm)").pack(side='left')
    ttk.Spinbox(row, textvariable=radius_var, from_=0.02, to=1.0, increment=0.05, width=6).pack(side='left', padx=6)

    buttons = ttk.Frame(frame)
    buttons.grid(row=5, column=0, columnspan=3, sticky='e', pady=(14, 0))
    make = ttk.Button(buttons, text="Make page")
    make.pack(side='left')
    ttk.Button(buttons, text="Close", command=root.destroy).pack(side='left', padx=(8, 0))
    status = ttk.Label(frame, textvariable=status_var, wraplength=620, justify='left')
    status.grid(row=6, column=0, columnspan=3, sticky='w', pady=(10, 0))

    def run():
        page, surgeries = page_var.get().strip().strip('"'), csv_var.get().strip().strip('"')
        try:
            radius = float(radius_var.get())
            if not 0 < radius <= 2:
                raise ValueError
        except ValueError:
            status_var.set("Enter the ball radius as a number of mm (e.g. 0.15).")
            return
        for label, path in (("3D page", page), ("surgery list", surgeries)):
            if not path or not Path(path).is_file():
                status_var.set(f"Choose the {label} first (no file at '{path}')." if path else f"Choose the {label} first.")
                return
        make.state(['disabled'])
        root.config(cursor='watch')
        status_var.set("Working... (a large page takes a few seconds)")
        save_last_used({'page': page, 'surgeries': surgeries, 'radius': radius})

        def work():
            try:
                out, report = make_page(page, surgeries, radius=radius)
                result = ('ok', out, report)
            except Exception as e:
                result = ('error', None, [str(e)])
            root.after(0, lambda: finished(*result))
        threading.Thread(target=work, daemon=True).start()

    def finished(kind, out, report):
        make.state(['!disabled'])
        root.config(cursor='')
        if kind == 'ok':
            status_var.set('\n'.join(report) + "\nOpening it in your browser.")
            import webbrowser
            webbrowser.open(Path(out).resolve().as_uri())
        else:
            status_var.set("Could not make the page: " + report[0])

    make.config(command=run)
    root.bind('<Return>', lambda _e: run())
    root.mainloop()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        return open_window()
    p = argparse.ArgumentParser(description="Draw the surgery viewer's injection sites as white balls in an ABC browser 3D page "
                                            "(with no arguments: a window to choose the files).")
    p.add_argument('page', help="a 3D page (.html) written by the ABC atlas browser")
    p.add_argument('surgeries', help="the surgery viewer's Export CSV file")
    p.add_argument('-o', '--out', help="the page to write (default: PAGE_injections.html beside PAGE)")
    p.add_argument('--radius', type=float, default=DEFAULT_RADIUS_MM, help=f"ball radius in mm (default {DEFAULT_RADIUS_MM})")
    a = p.parse_args(argv)
    try:
        out, report = make_page(a.page, a.surgeries, a.out, a.radius)
    except ValueError as e:
        sys.exit(str(e))
    print('\n'.join(report))
    return out


if __name__ == '__main__':
    main()
