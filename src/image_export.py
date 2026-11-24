"""Exporting a picture of an image (the XLEF viewer's Export): PNG, JPEG, TIFF, or SVG (the picture with a vector scale bar and label).

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
           cancel=None, quality=95):
    """Write the picture (width x height, unturned, given by read_band) to `path` as `fmt` (see FORMATS). progress(done, total) as
    bands come in; cancel() true stops it (Cancelled, the file removed). The written picture's size is (height, width) after an odd
    number of quarter turns."""
    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r}")
    why = fits(fmt, width, height)
    if why:
        raise ValueError(why)
    if fmt == "svg" and not (scale_bar and um_per_px):
        raise ValueError("SVG is for a picture with a scale bar")
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
                        + scale_bar_svg(final_w, final_h, um_per_px) + '</svg>\n')
            return path
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
