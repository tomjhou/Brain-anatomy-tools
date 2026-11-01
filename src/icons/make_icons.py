"""Draws the launcher icons (xlef_viewer.ico, ccf_browser.ico, abc_atlas_browser.ico) that INSTALL_APPS.bat gives the shortcuts.
The .ico files are kept in the repository, so this only needs running again to change a design:  python icons/make_icons.py"""
import os
import random
from PIL import Image, ImageDraw, ImageFilter

HERE = os.path.dirname(os.path.abspath(__file__))
BIG = 1024                                   # drawn at this size, then scaled down for each size in the .ico
SIZES = [(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]


def tile(fill):
    """A rounded square filling the canvas (transparent corners)."""
    im = Image.new("RGBA", (BIG, BIG), (0, 0, 0, 0))
    ImageDraw.Draw(im).rounded_rectangle((16, 16, BIG - 16, BIG - 16), radius=180, fill=fill)
    return im


def clip(layer, mask):
    """Keeps the layer only where the mask is set."""
    out = Image.new("RGBA", layer.size, (0, 0, 0, 0))
    out.paste(layer, (0, 0), Image.composite(layer.getchannel("A"), Image.new("L", layer.size, 0), mask))
    return out


def xlef_viewer():
    """A coronal section in two flat reds, traced from a real fluorescence image (xlef_section.webp): dim tissue, the two bright striata,
    and the ventricles in black. Needs numpy and scipy (both in the .venv once the ABC browser has run)."""
    import numpy as np
    from PIL import ImageEnhance
    from scipy import ndimage

    photo = Image.open(os.path.join(HERE, "xlef_section.webp")).convert("RGB")
    # The section alone: the largest bright area, without thin things (a stray fibre, the scale bar's text) or holes (the ventricles).
    tissue = ndimage.binary_opening(np.asarray(photo).max(axis=2) > 22, iterations=4)
    labels, _ = ndimage.label(tissue)
    section = ndimage.binary_fill_holes(ndimage.binary_closing(labels == np.bincount(labels.ravel())[1:].argmax() + 1, iterations=6))
    ys, xs = np.nonzero(section)
    box = (xs.min(), ys.min(), xs.max() + 1, ys.max() + 1)
    sec = photo.crop(box)
    sec = ImageEnhance.Color(ImageEnhance.Contrast(ImageEnhance.Brightness(sec).enhance(1.9)).enhance(1.3)).enhance(1.3)
    w = int(BIG * 0.92)                                        # filling most of the tile, a little taller than the photo
    h = int(w * sec.height / sec.width * 1.12)
    if h > BIG * 0.92:
        w, h = int(w * BIG * 0.92 / h), int(BIG * 0.92)
    x, y = (BIG - w) // 2, (BIG - h) // 2
    rgb = Image.new("RGB", (BIG, BIG)); rgb.paste(sec.resize((w, h), Image.LANCZOS), (x, y))
    mask = Image.new("L", (BIG, BIG)); mask.paste(Image.fromarray((section * 255).astype(np.uint8)).crop(box).resize((w, h)), (x, y))
    mask = mask.filter(ImageFilter.GaussianBlur(10)).point(lambda v: 255 if v > 128 else 0)   # a smooth outline
    inside = np.asarray(mask) > 0

    # Bright: the brightest third of the (smoothed) tissue, keeping only the two largest areas, the striata (no bright patches of cortex).
    lum = np.asarray(rgb.filter(ImageFilter.MedianFilter(31)).filter(ImageFilter.GaussianBlur(8)).convert("L")).astype(float)
    labels, _ = ndimage.label((lum > np.percentile(lum[inside], 67)) & inside)
    sizes = np.bincount(labels.ravel()); sizes[0] = 0
    bright = ndimage.binary_fill_holes(np.isin(labels, np.argsort(sizes)[-2:]))
    bright = np.asarray(Image.fromarray((bright * 255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(6))) > 128
    # Ventricles: large dark areas well inside the section.
    raw = np.asarray(rgb.convert("L").filter(ImageFilter.GaussianBlur(3))).astype(float)
    holes = ndimage.binary_opening((raw < 25) & ndimage.binary_erosion(inside, iterations=25), iterations=3)
    labels, _ = ndimage.label(holes)
    holes = (np.bincount(labels.ravel()) > 600)[labels] & (labels > 0)

    out = np.zeros((BIG, BIG, 3), np.uint8)
    out[inside] = (95, 18, 10)
    out[bright & inside] = (240, 45, 25)
    out[holes] = (0, 0, 0)
    layer = Image.fromarray(out).filter(ImageFilter.ModeFilter(9)).convert("RGBA")
    layer.putalpha(mask)
    im = tile((0, 0, 0, 255))
    im.alpha_composite(clip(layer, im.getchannel("A")))
    return im


def mouse_brain_mask():
    """A sagittal mouse brain outline (nose to the left): olfactory bulb, cortex dome, cerebellum, brainstem."""
    m = Image.new("L", (BIG, BIG), 0)
    d = ImageDraw.Draw(m)
    d.ellipse((230, 250, 820, 700), fill=255)                  # forebrain
    d.ellipse((90, 430, 330, 640), fill=255)                   # olfactory bulb
    d.ellipse((660, 420, 930, 730), fill=255)                  # cerebellum
    d.polygon([(560, 600), (820, 640), (900, 830), (700, 820)], fill=255)   # brainstem
    return m.filter(ImageFilter.GaussianBlur(6)).point(lambda v: 255 if v > 128 else 0)


def ccf_browser():
    """The Allen CCF: a sagittal mouse brain coloured by region (in the atlas's own colours) on blue."""
    im = tile((30, 60, 120, 255))
    mask = mouse_brain_mask()
    regions = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(regions)
    d.rectangle((0, 0, BIG, BIG), fill=(255, 144, 255, 255))   # midbrain / hindbrain (pink)
    d.ellipse((380, 420, 640, 640), fill=(255, 112, 128, 255))  # thalamus
    d.ellipse((230, 470, 470, 700), fill=(152, 214, 249, 255))  # striatum
    d.ellipse((470, 330, 760, 520), fill=(126, 208, 75, 255))   # hippocampus
    d.ellipse((150, 120, 860, 470), fill=(112, 255, 113, 255))  # isocortex (a dome over the top)
    d.ellipse((60, 400, 330, 660), fill=(154, 210, 189, 255))   # olfactory bulb
    d.ellipse((660, 420, 960, 740), fill=(255, 252, 145, 255))  # cerebellum
    im.alpha_composite(clip(regions, mask))
    lines = Image.new("RGBA", im.size, (0, 0, 0, 0))           # the outline, drawn as the mask's edge
    edge = mask.filter(ImageFilter.FIND_EDGES).filter(ImageFilter.MaxFilter(13))
    lines.paste((20, 30, 60, 255), (0, 0), edge)
    im.alpha_composite(lines)
    return im


def abc_atlas_browser():
    """A UMAP: clusters of coloured cell dots on white."""
    im = tile((250, 250, 250, 255))
    dots = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(dots)
    rnd = random.Random(7)
    clusters = [((300, 330), (230, 60, 60)), ((640, 280), (40, 120, 230)), ((720, 620), (40, 170, 80)),
                ((330, 700), (240, 160, 20)), ((510, 500), (150, 70, 200))]
    for (cx, cy), colour in clusters:
        for _ in range(55):
            x, y, r = cx + rnd.gauss(0, 70), cy + rnd.gauss(0, 60), rnd.uniform(16, 30)
            d.ellipse((x - r, y - r, x + r, y + r), fill=colour + (255,))
    im.alpha_composite(clip(dots, im.getchannel("A")))
    ImageDraw.Draw(im).rounded_rectangle((16, 16, BIG - 16, BIG - 16), radius=180, outline=(150, 150, 150, 255), width=14)
    return im


if __name__ == "__main__":
    for name, draw in (("xlef_viewer", xlef_viewer), ("ccf_browser", ccf_browser), ("abc_atlas_browser", abc_atlas_browser)):
        big = draw()
        big.resize((256, 256), Image.LANCZOS).save(os.path.join(HERE, name + ".ico"), sizes=SIZES)
        print("wrote", name + ".ico")
