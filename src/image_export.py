"""Exporting a picture of an image (the XLEF viewer's Export): PNG, JPEG, TIFF, or SVG (the picture with a vector scale bar and label,
and the annotations as vector shapes and text).

The caller gives the picture a band of rows at a time, unturned (read_band(r0, r1) -> RGB uint8 [r1 - r0, width, 3]); it is turned a
number of quarter turns clockwise and mirrored left for right as the viewer shows it, the scale bar is drawn at its bottom left (as
on screen), and it is written. A TIFF is written straight into its file (a memory-mapped TIFF: any size, a slide at full resolution
too, BigTIFF past 4 GB); the others are made in memory, so are limited to MAX_MEMORY_PIXELS.

    export(path, fmt, width, height, read_band, quarter=0, flipped=False, um_per_px=None, scale_bar=False, progress=None, cancel=None)
"""
import base64
import io
import math
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

FORMATS = {"png": ".png", "jpeg": ".jpg", "tiff": ".tif", "svg": ".svg"}
MAX_MEMORY_PIXELS = 300_000_000      # PNG, JPEG and SVG are made whole in memory: 3 bytes a pixel (900 MB), and more while encoding
JPEG_MAX_SIDE = 65500                # (the JPEG format's limit)
BAND_ROWS = 512


class Cancelled(Exception):
    pass


def fits(fmt, width, height):
    """None if a width x height picture can be exported as `fmt`; else why not, in a few words."""
    if fmt == "tiff":
        return None
    if fmt == "jpeg" and max(width, height) > JPEG_MAX_SIDE:
        return f"JPEG allows at most {JPEG_MAX_SIDE:,} pixels a side"
    if width * height > MAX_MEMORY_PIXELS:
        return f"too large for {fmt.upper()} (at most {MAX_MEMORY_PIXELS // 1_000_000} million pixels): use TIFF"
    return None


def nice_length(um):
    """The largest 1/2/5 x 10^n not above um."""
    if um <= 0:
        return 0
    e = 10 ** math.floor(math.log10(um))
    return next((m * e for m in (5, 2, 1) if m * e <= um), e)


def scale_bar_geometry(width, height, um_per_px):
    """The scale bar of a width x height picture: a dict of x, y (its bar's top left), length, thickness, the label, its font size,
    and where the label's bottom middle is. About a fifth of the width, at the bottom left, as the viewer draws it."""
    um = nice_length(width / 5 * um_per_px)
    length = um / um_per_px
    unit = max(1.0, min(width, height) / 100)       # (sizes in proportion to the picture: legible small and large)
    thickness = max(3, round(unit * 0.8))
    font = max(10, round(unit * 2.6))
    x = round(unit * 2.5)
    y = height - round(unit * 2.5) - thickness
    label = f"{um / 1000:g} mm" if um >= 1000 else f"{um:g} µm"
    return {"x": x, "y": y, "length": length, "thickness": thickness, "label": label, "font": font,
            "label_x": x + length / 2, "label_y": y - max(2, round(unit * 0.6))}


def _font(size):
    for name in ("arial.ttf", "Arial.ttf", "segoeui.ttf", "Helvetica.ttc", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def draw_scale_bar(im, um_per_px, size=None, top=0):
    """The scale bar drawn on the PIL image `im` (RGB), in place: a white bar with a black edge, its length above it in white with a
    black outline (legible on any background). size, top: `im` is the part of a picture of that size (w, h) from row `top` down
    (the bar is placed and sized for the whole picture)."""
    w, h = size or im.size
    g = scale_bar_geometry(w, h, um_per_px)
    d = ImageDraw.Draw(im)
    d.rectangle([g["x"], g["y"] - top, g["x"] + g["length"], g["y"] + g["thickness"] - top], fill="white", outline="black",
                width=max(1, g["thickness"] // 4))
    d.text((g["label_x"], g["label_y"] - top), g["label"], fill="white", font=_font(g["font"]), anchor="ms",
           stroke_width=max(1, g["font"] // 12), stroke_fill="black")


def scale_bar_svg(width, height, um_per_px):
    """The scale bar as SVG elements (see scale_bar_geometry), to lay over the picture."""
    g = scale_bar_geometry(width, height, um_per_px)
    edge = max(1, g["thickness"] // 4)
    label = g["label"].replace("&", "&amp;").replace("<", "&lt;")
    return (f'  <g id="scale-bar">\n'
            f'    <rect x="{g["x"]}" y="{g["y"]}" width="{g["length"]:.2f}" height="{g["thickness"]}" fill="#ffffff" stroke="#000000" '
            f'stroke-width="{edge}"/>\n'
            f'    <text x="{g["label_x"]:.2f}" y="{g["label_y"]}" font-family="Arial, Helvetica, sans-serif" font-size="{g["font"]}" '
            f'text-anchor="middle" fill="#ffffff" stroke="#000000" stroke-width="{max(1, g["font"] / 12):.1f}" '
            f'paint-order="stroke">{label}</text>\n'
            f'  </g>\n')


# --- annotations: drawn into the picture (PNG, JPEG, TIFF) or as vector elements over it (SVG)
#
# The caller gives them as plain shapes, already in the picture's own pixels (after it is turned and mirrored: see orient_point), so
# that one description serves every format:
#     {"kind": "polygon", "points": [(x, y), ...], "color": "#rrggbb", "width": w}       (closed outline, not filled)
#     {"kind": "line", "points": [(x, y), ...], "color": ..., "width": w}                 (open)
#     {"kind": "head", "points": [tip, left, right], "color": ...}                        (a filled arrowhead)
#     {"kind": "text", "x": x, "y": y, "text": "...", "color": ..., "size": px, "bold": False, "anchor": "nw" | "w" | ...,
#      "border": False}      (anchor: which point of the text's box is at (x, y), as Tk's: n/s/e/w, combinations, or "center")

def orient_point(x, y, width, height, quarter, flipped):
    """A point of the unturned picture (width x height, x right, y down, in pixels: 0..width) -> where it is once the picture is
    turned `quarter` quarter turns clockwise and then mirrored left for right (as _oriented_band does)."""
    q = quarter % 4
    if q == 1:
        x, y = height - y, x
    elif q == 2:
        x, y = width - x, height - y
    elif q == 3:
        x, y = y, width - x
    if flipped:
        x = (height if q % 2 else width) - x
    return x, y


def _text_box(shape, font):
    """(left, top, right, bottom) of a text shape's lines, placed by its anchor, and [(line, x, y-top)] to draw them (PIL, "la")."""
    lines = shape["text"].split("\n") or [""]
    step = round(shape["size"] * 1.2)
    widths = [font.getlength(t) for t in lines]
    w, h = max(widths, default=0), step * len(lines)
    anchor = shape.get("anchor", "nw")
    x, y = shape["x"], shape["y"]
    left = x - (w if "e" in anchor else 0 if "w" in anchor else w / 2)
    top = y - (h if "s" in anchor else 0 if "n" in anchor else h / 2)
    return (left, top, left + w, top + h), [(t, left, top + i * step) for i, t in enumerate(lines)]


def _shape_bounds(shape, font=None):
    """A shape's box in the picture (left, top, right, bottom), its line width around it."""
    if shape["kind"] == "text":
        box = _text_box(shape, font or _font(max(1, round(shape["size"]))))[0]
        pad = 2 + shape["size"] / 8
        return box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    xs, ys = [p[0] for p in shape["points"]], [p[1] for p in shape["points"]]
    pad = shape.get("width", 1) + 2
    return min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad


def draw_shapes(im, shapes, left=0, top=0):
    """The shapes drawn on the PIL image `im` (RGB), in place; `im` is the part of the picture whose top-left corner is (left, top)."""
    d = ImageDraw.Draw(im)
    for s in shapes:
        col = s["color"]
        if s["kind"] == "text":
            font = _font(max(1, round(s["size"])))
            box, lines = _text_box(s, font)
            for t, x, y in lines:
                d.text((x - left, y - top), t, fill=col, font=font, anchor="la")
            if s.get("border"):
                d.rectangle([box[0] - 2 - left, box[1] - 2 - top, box[2] + 2 - left, box[3] + 2 - top], outline=col,
                            width=max(1, round(s["size"] / 12)))
            continue
        pts = [(x - left, y - top) for x, y in s["points"]]
        if s["kind"] == "head":
            d.polygon(pts, fill=col)
        else:
            w = max(1, round(s.get("width", 2)))
            d.line(pts + (pts[:1] if s["kind"] == "polygon" else []), fill=col, width=w, joint="curve")


def _draw_shapes_on(out, shapes, final_w, final_h):
    """The shapes drawn into the array `out` (a TIFF's memory map: final_h x final_w x 3), each on the part of it under the shape only,
    so that a picture of any size needn't be in memory whole."""
    for s in shapes:
        x0, y0, x1, y1 = _shape_bounds(s)
        x0, y0 = max(0, int(math.floor(x0))), max(0, int(math.floor(y0)))
        x1, y1 = min(final_w, int(math.ceil(x1)) + 1), min(final_h, int(math.ceil(y1)) + 1)
        if x1 <= x0 or y1 <= y0:
            continue
        part = Image.fromarray(np.array(out[y0:y1, x0:x1]))
        draw_shapes(part, [s], x0, y0)
        out[y0:y1, x0:x1] = np.asarray(part)


def _esc(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def shapes_svg(shapes):
    """The shapes as SVG elements (strokes, arrowheads and text that stay editable), to lay over the picture."""
    if not shapes:
        return ""
    out = ['  <g id="annotations" fill="none" stroke-linecap="round" stroke-linejoin="round">\n']
    for s in shapes:
        col = s["color"]
        if s["kind"] == "text":
            lines = s["text"].split("\n") or [""]
            size = s["size"]
            step = size * 1.2
            anchor = s.get("anchor", "nw")
            ta = "end" if "e" in anchor else "start" if "w" in anchor else "middle"
            h = step * len(lines)
            top = s["y"] - (h if "s" in anchor else 0 if "n" in anchor else h / 2)
            weight = ' font-weight="bold"' if s.get("bold") else ""
            spans = "".join(f'<tspan x="{s["x"]:.2f}" y="{top + size + i * step:.2f}">{_esc(t)}</tspan>' for i, t in enumerate(lines))
            out.append(f'    <text font-family="Arial, Helvetica, sans-serif" font-size="{size:.2f}"{weight} text-anchor="{ta}" '
                       f'fill="{col}" stroke="none">{spans}</text>\n')
            if s.get("border"):
                box = _text_box(s, _font(max(1, round(size))))[0]
                out.append(f'    <rect x="{box[0] - 2:.2f}" y="{box[1] - 2:.2f}" width="{box[2] - box[0] + 4:.2f}" '
                           f'height="{box[3] - box[1] + 4:.2f}" stroke="{col}" stroke-width="{max(1, size / 12):.2f}"/>\n')
            continue
        pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in s["points"])
        if s["kind"] == "head":
            out.append(f'    <polygon points="{pts}" fill="{col}" stroke="none"/>\n')
        else:
            tag = "polygon" if s["kind"] == "polygon" else "polyline"
            out.append(f'    <{tag} points="{pts}" stroke="{col}" stroke-width="{s.get("width", 2):.2f}"/>\n')
    out.append('  </g>\n')
    return "".join(out)


def _placement(r0, r1, width, height, quarter, flipped):
    """Where band rows r0..r1 of the unturned picture (width x height) go in the turned and mirrored one: ("rows" | "cols", a, b)."""
    q = quarter % 4
    if q == 0:
        where = ("rows", r0, r1)
    elif q == 1:
        where = ("cols", height - r1, height - r0)
    elif q == 2:
        where = ("rows", height - r1, height - r0)
    else:
        where = ("cols", r0, r1)
    if flipped and where[0] == "cols":
        final_w = height if q % 2 else width
        where = ("cols", final_w - where[2], final_w - where[1])
    return where


def _oriented_band(band, quarter, flipped):
    """A band turned `quarter` quarter turns clockwise, then mirrored left for right (see _placement for where it goes)."""
    b = np.rot90(band, -(quarter % 4))
    return b[:, ::-1] if flipped else b


def export(path, fmt, width, height, read_band, quarter=0, flipped=False, um_per_px=None, scale_bar=False, progress=None,
           cancel=None, quality=95, shapes=()):
    """Write the picture (width x height, unturned, given by read_band) to `path` as `fmt` (see FORMATS). progress(done, total) as
    bands come in; cancel() true stops it (Cancelled, the file removed). The written picture's size is (height, width) after an odd
    number of quarter turns. shapes: annotations, in the written picture's pixels (see orient_point and the shapes above draw_shapes),
    drawn into it, or for an SVG laid over it as vector elements."""
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r}")
    why = fits(fmt, width, height)
    if why:
        raise ValueError(why)
    shapes = list(shapes or ())
    if fmt == "svg" and not (scale_bar and um_per_px) and not shapes:
        raise ValueError("SVG is for a picture with a scale bar or annotations")
    q = quarter % 4
    final_w, final_h = (height, width) if q % 2 else (width, height)
    dpi = 25400 / um_per_px if um_per_px else None
    out = None
    try:
        if fmt == "tiff":
            import tifffile
            kwargs = {"photometric": "rgb", "bigtiff": final_w * final_h * 3 > 3_900_000_000}
            if um_per_px:
                kwargs.update(resolution=(1e4 / um_per_px, 1e4 / um_per_px), resolutionunit="CENTIMETER")
            out = tifffile.memmap(path, shape=(final_h, final_w, 3), dtype=np.uint8, **kwargs)
        else:
            out = np.zeros((final_h, final_w, 3), np.uint8)
        for r0 in range(0, height, BAND_ROWS):
            if cancel is not None and cancel():
                raise Cancelled()
            r1 = min(height, r0 + BAND_ROWS)
            band = read_band(r0, r1)
            kind, a, b = _placement(r0, r1, width, height, q, flipped)
            block = _oriented_band(band, q, flipped)
            if kind == "rows":
                out[a:b] = block
            else:
                out[:, a:b] = block
            if progress is not None:
                progress(r1, height)
        if cancel is not None and cancel():
            raise Cancelled()
        if fmt == "tiff":
            _draw_shapes_on(out, shapes, final_w, final_h)      # (first: the scale bar goes on top of them)
            if scale_bar and um_per_px:                 # (drawn on the bottom-left corner only: the rest stays in the file)
                g = scale_bar_geometry(final_w, final_h, um_per_px)
                top = max(0, g["label_y"] - g["font"] * 2)
                right = min(final_w, int(g["x"] + max(g["length"], g["font"] * len(g["label"])) + g["font"] * 2))
                corner = Image.fromarray(np.array(out[top:, :right]))
                draw_scale_bar(corner, um_per_px, (final_w, final_h), top)
                out[top:, :right] = np.asarray(corner)
            out.flush()
            del out
            out = None
            return path
        im = Image.fromarray(out)
        out = None
        if fmt == "svg":
            buf = io.BytesIO()
            im.save(buf, "PNG", **({"dpi": (dpi, dpi)} if dpi else {}))
            data = base64.b64encode(buf.getvalue()).decode("ascii")
            with open(path, "w", encoding="utf-8") as f:
                f.write(f'<?xml version="1.0" encoding="UTF-8"?>\n'
                        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
                        f'width="{final_w}" height="{final_h}" viewBox="0 0 {final_w} {final_h}">\n'
                        f'  <image id="image" x="0" y="0" width="{final_w}" height="{final_h}" preserveAspectRatio="none" '
                        f'xlink:href="data:image/png;base64,{data}"/>\n'
                        + shapes_svg(shapes)
                        + (scale_bar_svg(final_w, final_h, um_per_px) if scale_bar and um_per_px else "") + '</svg>\n')
            return path
        draw_shapes(im, shapes)
        if scale_bar and um_per_px:
            draw_scale_bar(im, um_per_px)
        extra = {"dpi": (dpi, dpi)} if dpi else {}
        if fmt == "png":
            im.save(path, "PNG", **extra)
        else:
            im.save(path, "JPEG", quality=quality, subsampling=0, **extra)
        return path
    except BaseException:
        out = None
        try:
            if os.path.exists(path):
                os.remove(path)
        except OSError:
            pass
        raise
