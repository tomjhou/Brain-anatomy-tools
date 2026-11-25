"""The ABC browser's 3D pages are stored compressed: the page is mostly numbers (Plotly itself and the surfaces' triangle lists shrink to a third), so the
file is a few lines of script around the page, deflated and written as base64. The script unpacks it when the page opens (the browser's own
DecompressionStream) and writes it out as the page. The picture is the same, the file about half the size. (A browser from before 2020 to early 2023 says so
instead of showing the page.)

This module packs and unpacks such pages, so that the ABC browser (writing them) and add_injections_3d.py (adding the injection-sites loader to pages written
before it was included: upgrade_page) work on the same format. Nothing here needs more than the standard library."""
from __future__ import annotations

import base64
import os
import re
import zlib
from pathlib import Path

import injection_loader

MARKER = 'packed-3d-page'

WRAPPER = """<!doctype html>
<html>
<head>
    <meta charset="utf-8" />
    <!-- __MARKER__ __FLAGS__ -->
    <style>html, body {height: 100%;} body {background: #111; color: #aab; font: 14px sans-serif;}</style>
</head>
<body>
<div id="unpacking" style="padding:12px">Opening the 3D page...</div>
<script>
(async function () {
  var raw = atob('__PAGE__'), bytes = new Uint8Array(raw.length);
  for (var q = 0; q < raw.length; q++) bytes[q] = raw.charCodeAt(q);
  var html = await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'))).text();
  document.open(); document.write(html); document.close();
})().catch(function (error) {
  document.getElementById('unpacking').textContent = 'This page needs a newer browser (it could not unpack itself: ' + error + '). Please update the browser.';
});
</script>
</body>
</html>
"""

OLD_INJECTION_BLOCK = re.compile(re.escape('<!-- injection-sites:start -->') + r'.*?' + re.escape('<!-- injection-sites:end -->') + r'\n?', re.S)


def pack_text(html):
    """The page `html` (text) as the compressed page (text). The first lines say whether it carries the injection-sites loader (see has_loader)."""
    flags = 'loader=1' if injection_loader.LOADER_MARK in html else 'loader=0'
    packed = base64.b64encode(zlib.compress(html.encode('utf-8'), 9)).decode('ascii')
    return WRAPPER.replace('__MARKER__', MARKER).replace('__FLAGS__', flags).replace('__PAGE__', packed)


def is_packed(head):
    """Whether `head` (the start of a file, text or bytes) is the start of a compressed page: the current kind (marked MARKER) or the first
    test version of it (no marker, but the same unpacking script: it has the 'unpacking' note)."""
    if isinstance(head, bytes):
        head = head.decode('utf-8', 'ignore')
    return MARKER in head[:2000] or 'id="unpacking"' in head[:2000]


def unpack_text(text):
    """The page inside a compressed page (text), or None if `text` is not one."""
    if not is_packed(text):
        return None
    found = re.search(r"atob\('([A-Za-z0-9+/=]*)'\)", text)
    if not found:
        return None
    return zlib.decompress(base64.b64decode(found.group(1))).decode('utf-8')


def write_atomic(path, text):
    """Writes `text` to `path` through a temporary file beside it, so that a reader never sees half a page (and a failed write leaves the old one)."""
    path = Path(path)
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(text, encoding='utf-8')
    os.replace(temporary, path)


def pack_file(path):
    """Rewrites the plain page at `path` as the compressed page; returns (size before, size after) in bytes."""
    path = Path(path)
    html = path.read_text(encoding='utf-8')
    before = path.stat().st_size
    write_atomic(path, pack_text(html))
    return before, path.stat().st_size


def read_page(path):
    """(the page as plain HTML text, whether the file was a compressed page)."""
    text = Path(path).read_text(encoding='utf-8')
    inner = unpack_text(text)
    return (inner, True) if inner is not None else (text, False)


def has_loader(path):
    """Whether the 3D page at `path` carries the injection-sites loader, found without unpacking it (the compressed page says so in its first lines; a plain
    page has the loader near its end)."""
    with open(path, 'rb') as f:
        head = f.read(4000)
        if is_packed(head):
            return b'loader=1' in head       # (the first test version has no flag: it has no loader)
        f.seek(0, os.SEEK_END)
        size = f.tell()
        f.seek(max(0, size - 400000))
        return injection_loader.LOADER_MARK.encode() in f.read()


def looks_like_3d_page(path):
    """Whether `path` is a 3D page from the ABC browser (compressed or plain)."""
    with open(path, 'rb') as f:
        head = f.read(20000)
    return is_packed(head) or b'plotly' in head.lower()


def add_loader(html):
    """The page `html` (text) with the injection-sites loader in it, in place of the older kind of injection script some pages carry."""
    html = OLD_INJECTION_BLOCK.sub('', html)
    if injection_loader.LOADER_MARK in html:
        return html
    script = '<script>' + injection_loader.LOADER_JS + '</script>\n'
    at = html.rfind('</body>')
    return html + script if at < 0 else html[:at] + script + html[at:]


def upgrade_page(path):
    """Gives the page at `path` the injection-sites loader (and stores it compressed). The page looks the same; only the file changes, through a temporary
    file, so a failure leaves the original."""
    html, _was_packed = read_page(path)
    write_atomic(path, pack_text(add_loader(html)))
