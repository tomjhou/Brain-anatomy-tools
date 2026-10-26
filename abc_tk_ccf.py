"""The Allen CCF (2020) atlas drawn over the UMAP viewer's section maps, as the outlines of its regions.

The cells carry their position in the CCF (the registration the atlas ships), so each section's own (x, y) can be fitted to CCF
coordinates: a low-order polynomial for each section, which follows the section's tilt and warp (a section spans a few hundred
microns of the AP axis). Every pixel of a section map is then looked up in the 25 um annotation volume, and the boundary of each
region is traced and kept as a polygon in the section's own mm coordinates, so it can be drawn crisply at any zoom."""
import cv2
import numpy as np
from scipy import ndimage

VOXEL_MM = 0.025              # the annotation volume: [AP, DV, ML], 25 um voxels, labels are the cells' parcellation_index
FIT_DEGREE = 2
RASTER_MM = 0.0125            # the grid the regions are traced on (half a voxel)
SMOOTH_SIGMA_PX = 1.6         # blurs the voxel steps out of a region's mask before its boundary is traced
MIN_CONTOUR_PX = 6            # a traced boundary shorter than this is dropped
SLICE_PITCH_MM = 0.01         # the resampled annotation (Allen's reconstructed space): 10 um voxels in the section's plane
SLICE_SPACING_MM = 0.2        # ... and one slice per section, this far apart
SLICE_SMOOTH_SIGMA_PX = 2.0   # (the same blur as above, in the finer pixels)


def poly_terms(x, y, center, scale):
    u, v = (x - center[0]) / scale, (y - center[1]) / scale
    return np.column_stack([u ** i * v ** j for i in range(FIT_DEGREE + 1) for j in range(FIT_DEGREE + 1 - i)])


def fit_section(x, y, ccf):
    """The fit of a section's CCF coordinates (n, 3: AP, DV, ML in mm) on its own (x, y): (center, scale, coefficients (terms, 3)), or None
    if too few cells have a CCF position."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(ccf).all(axis=1)
    if ok.sum() < 200:
        return None
    x, y, ccf = x[ok], y[ok], ccf[ok]
    center = (float(x.mean()), float(y.mean()))
    scale = float(max(np.ptp(x), np.ptp(y), 1e-6))
    coef, *_ = np.linalg.lstsq(poly_terms(x, y, center, scale), ccf, rcond=None)
    return center, scale, coef


def trace_regions(volume, fit, home):
    """The boundaries of the CCF regions across the section map whose view is `home` (xmin, xmax, ymin, ymax; mm): a list of (n, 2)
    float arrays of (x, y) in mm, one closed outline each."""
    center, scale, coef = fit
    x0, x1, y0, y1 = home
    xs = np.arange(x0, x1, RASTER_MM) + RASTER_MM / 2
    ys = np.arange(y0, y1, RASTER_MM) + RASTER_MM / 2
    gx, gy = np.meshgrid(xs, ys)
    ccf = poly_terms(gx.ravel(), gy.ravel(), center, scale) @ coef
    idx = np.floor(ccf / VOXEL_MM).astype(np.int64)
    inside = np.ones(len(idx), bool)
    for axis, n in enumerate(volume.shape):
        inside &= (idx[:, axis] >= 0) & (idx[:, axis] < n)
    labels = np.zeros(len(idx), np.int64)
    labels[inside] = volume[idx[inside, 0], idx[inside, 1], idx[inside, 2]]
    return trace_label_grid(labels.reshape(gy.shape), x0, y0, RASTER_MM, SMOOTH_SIGMA_PX)


def trace_label_grid(labels, x0, y0, pitch, sigma):
    """The boundaries of the regions in `labels` (rows: y, columns: x; 0 is outside the brain), as a list of (n, 2) float arrays of (x, y) in mm, one
    closed outline each. `x0`, `y0` is the corner of the grid and `pitch` the size of a pixel (mm); each region's mask is blurred by `sigma` pixels
    first, to take the steps out of its outline."""
    labels = np.asarray(labels).astype(np.int64)
    polygons = []
    pad = int(3 * sigma) + 2
    for region, box in enumerate(ndimage.find_objects(labels), start=1):          # (the list's first item is label 1)
        if box is None:
            continue
        ys, xs = box
        r0, r1 = max(ys.start - pad, 0), min(ys.stop + pad, labels.shape[0])
        c0, c1 = max(xs.start - pad, 0), min(xs.stop + pad, labels.shape[1])
        mask = (labels[r0:r1, c0:c1] == region).astype(np.uint8) * 255
        if sigma:
            mask = cv2.GaussianBlur(mask, (0, 0), sigma)
        mask = (mask > 127).astype(np.uint8)
        contours, _ = cv2.findContours(mask, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
        for c in contours:
            if len(c) < 3 or cv2.arcLength(c, True) < MIN_CONTOUR_PX:
                continue
            p = c[:, 0, :].astype(np.float64)
            polygons.append(np.column_stack([x0 + (c0 + p[:, 0] + 0.5) * pitch, y0 + (r0 + p[:, 1] + 0.5) * pitch]))
    return polygons


def slice_outlines(volume, section_z, homes):
    """{section: polygons} from the CCF annotation resampled into Allen's reconstructed space: `volume` is [x, y, slice] (10 um voxels, one slice per
    section, SLICE_SPACING_MM apart), the cells' reconstructed (x, y, z) read straight off it, with no fitting. `section_z`: {section: z in mm};
    `homes`: {section: (xmin, xmax, ymin, ymax) of its map, mm}. A section without a slice (one Allen did not register) gets an empty list."""
    out = {}
    for section, home in homes.items():
        k = int(round(section_z.get(section, np.nan) / SLICE_SPACING_MM)) if np.isfinite(section_z.get(section, np.nan)) else -1
        if not 0 <= k < volume.shape[2]:
            out[section] = []
            continue
        x0, x1, y0, y1 = home
        i0, i1 = max(int(np.floor(x0 / SLICE_PITCH_MM)), 0), min(int(np.ceil(x1 / SLICE_PITCH_MM)), volume.shape[0])
        j0, j1 = max(int(np.floor(y0 / SLICE_PITCH_MM)), 0), min(int(np.ceil(y1 / SLICE_PITCH_MM)), volume.shape[1])
        if i1 <= i0 or j1 <= j0:
            out[section] = []
            continue
        grid = np.asarray(volume[i0:i1, j0:j1, k]).T
        out[section] = trace_label_grid(grid, i0 * SLICE_PITCH_MM, j0 * SLICE_PITCH_MM, SLICE_PITCH_MM, SLICE_SMOOTH_SIGMA_PX)
    return out


def section_outlines(volume, cells, homes):
    """{section: polygons} for every section in `homes` ({section: its map's view}); `cells`: {section: (x, y, ccf (n, 3))}. A section with
    too few registered cells gets an empty list."""
    out = {}
    for section, home in homes.items():
        data = cells.get(section)
        fit = fit_section(*data) if data is not None else None
        out[section] = trace_regions(volume, fit, home) if fit is not None else []
    return out
