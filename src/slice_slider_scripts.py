"""The JavaScript of the 3D pages' slice sliders that is not part of the cutting itself, shared by the ABC browser (which writes it into every 3D page, in its Slice
script) and add_injections_3d.py (whose side file adds it to pages written before it existed). Standard library only.

DRAG_JS is always included: a press on a slider thumb is taken by the page, which captures the pointer and sets the slider from the pointer's position at every move.
The browser's own drag of a thumb was found to be dropped part-way on some computers (the thumb stopped following while the pointer went on).

DIAGNOSTICS_JS is included only when the Settings window's Troubleshooting tab has "Show diagnostics in 3D pages" ticked (stored as 'diagnostics_3d' in
viewer_settings.json; see diagnostics_enabled): a line under the sliders naming the graphics that draw the 3D scene, and an amber note when a thumb falls well behind
the pointer or a press on a slider is taken by something else, with the details in window.__sliderLog. With the setting off the pages contain none of it."""
import json

DRAG_JS = r'''
  // The slice sliders' drag. The browser's own drag of a slider's thumb can be dropped part-way (the thumb stops following while the pointer goes on), so this takes
  // the press itself: the pointer is captured by the slider and its value is set from the pointer's position at every move, for as long as the button is down.
  // The page's own 'input' handlers run as before.
  (function () {
    if (window.__sliderDrag) return;
    var panel = document.querySelector('.slice-box');
    if (!panel) return;
    window.__sliderDrag = true;
    var drag = null, THUMB = 7;   // (half the width of a thumb)
    Array.prototype.forEach.call(panel.querySelectorAll('input.slice-dual'), function (input) { input.style.touchAction = 'none'; });
    function valueAt(input, x, offset) {   // the slider's value for a pointer at x (the pointer keeps its place on the thumb)
      var r = input.getBoundingClientRect(), min = +input.min, span = +input.max - min;
      var frac = Math.min(1, Math.max(0, (x - offset - r.x - THUMB) / ((r.width - 2 * THUMB) || 1)));
      return min + frac * span;
    }
    document.addEventListener('pointerdown', function (e) {
      var t = e.target;
      if (!t || !t.classList || !t.classList.contains('slice-dual') || e.button !== 0) return;
      var r = t.getBoundingClientRect(), min = +t.min, span = (+t.max - min) || 1;
      var thumbX = r.x + THUMB + ((+t.value - min) / span) * (r.width - 2 * THUMB);
      drag = {input: t, offset: e.clientX - thumbX, id: e.pointerId};
      e.preventDefault();                       // (no drag by the browser: it is done here, below)
      try { t.setPointerCapture(e.pointerId); } catch (err) { /* the document-level handler below still follows the pointer */ }
    }, true);
    document.addEventListener('pointermove', function (e) {
      if (!drag || !(e.buttons & 1)) return;   // (a move with the button no longer down is ignored: the release will end the drag)
      var want = valueAt(drag.input, e.clientX, drag.offset);
      if (+drag.input.value !== want) {
        drag.input.value = want;
        drag.input.dispatchEvent(new Event('input', {bubbles: true}));
      }
    }, true);
    function finish() {
      if (!drag) return;
      try { drag.input.releasePointerCapture(drag.id); } catch (err) { /* it was not held */ }
      drag.input.dispatchEvent(new Event('change', {bubbles: true}));
      drag = null;
    }
    document.addEventListener('pointerup', finish, true);
    document.addEventListener('pointercancel', finish, true);
  })();
'''

DIAGNOSTICS_JS = r'''
  // Diagnostics (Settings > Troubleshooting). Which renderer draws the 3D scene: the sliders are only as quick as it is, and a browser can fall back to software
  // drawing (about ten times slower) after its graphics process has had trouble. Shown under the sliders.
  (function () {
    var box = document.querySelector('.slice-box');
    if (!box || document.querySelector('.renderer-note')) return;
    var name = '';
    try {
      var canvas = document.createElement('canvas'), gl = canvas.getContext('webgl'), info = gl && gl.getExtension('WEBGL_debug_renderer_info');
      name = gl ? String(gl.getParameter(info ? info.UNMASKED_RENDERER_WEBGL : gl.RENDERER)) : 'WebGL unavailable';
      var lose = gl && gl.getExtension('WEBGL_lose_context');
      if (lose) lose.loseContext();
    } catch (err) { name = 'unknown'; }
    var software = /swiftshader|llvmpipe|software|basic render/i.test(name);
    var note = document.createElement('div');
    note.className = 'renderer-note';
    note.textContent = (software ? 'Software 3D, slow: ' : '3D drawn by: ') + name;
    note.style.cssText = 'margin-top:6px;font:10px sans-serif;max-width:170px;word-wrap:break-word;color:' + (software ? '#e0a24f' : '#8b93a6') + ';';
    box.appendChild(note);
  })();
  // A drag that goes wrong: if a slider's thumb stays well behind the pointer (more than 8% of the track) for over 0.6 s while the button is held down, say so under
  // the sliders and keep the details in window.__sliderLog (F12, Console): whether mouse moves were still reaching the page and how long the page had been blocked,
  // what was under the pointer, and how the press ended. A press on a slider that something else takes is reported too.
  (function () {
    if (window.__sliderDiag) return;
    var panel = document.querySelector('.slice-box');
    if (!panel) return;
    window.__sliderDiag = true;
    window.__sliderLog = [];
    var watch = null, line = null, THUMB = 7;
    function describe(el) { return el ? (el.tagName + (el.className && el.className.baseVal === undefined && el.className ? '.' + String(el.className).split(' ')[0] : '')) : 'nothing'; }
    function report(why, e) {
      var found = {why: why, afterMs: Math.round(performance.now() - watch.t0), moves: watch.moves, inputEvents: watch.inputEvents, longestGapMs: Math.round(watch.maxGap),
                   msSinceLastMove: Math.round(performance.now() - watch.last), buttons: e ? e.buttons : null, pointerType: watch.type,
                   under: describe(document.elementFromPoint(watch.x, watch.y)), windowFocused: document.hasFocus(), offBy: Math.round(100 * watch.offBy) + '%'};
      window.__sliderLog.push(found);
      console.warn('slice slider lost the pointer', found);
      if (!line) {
        line = document.createElement('div');
        line.style.cssText = 'margin-top:6px;font:10px sans-serif;max-width:170px;word-wrap:break-word;color:#e0a24f;';
        panel.appendChild(line);
      }
      line.textContent = 'Slider fell behind the pointer: ' + found.why + '; page blocked at most ' + found.longestGapMs + ' ms; under the pointer: ' + found.under +
                         '; ' + found.moves + ' moves, ' + found.inputEvents + ' slider events (' + found.pointerType + '; details: __sliderLog in the console)';
    }
    function newWatch(input, e, offset) {
      return {t0: performance.now(), last: performance.now(), maxGap: 0, moves: 0, inputEvents: 0, input: input, type: e.pointerType, x: e.clientX, y: e.clientY,
              offSince: 0, offBy: 0, reported: false, offset: offset || 0};
    }
    document.addEventListener('pointerdown', function (e) {
      var t = e.target;
      if (!t || !t.classList || !t.classList.contains('slice-dual')) {
        // a press that lands on a slider's track but is taken by something else (a cover over the slider, or the plot itself) would not move the thumb at all
        var tracks = panel.querySelectorAll('input.slice-dual');
        for (var q = 0; q < tracks.length; q++) {
          var box = tracks[q].getBoundingClientRect();
          if (box.width > 0 && e.clientX >= box.x && e.clientX <= box.right && e.clientY >= box.y - 2 && e.clientY <= box.bottom + 2 && !(t && panel.contains(t) && t !== panel)) {
            watch = newWatch(tracks[q], e, 0);
            report('a press on a slider was taken by ' + describe(t), e);
            watch = null;
            break;
          }
        }
        return;
      }
      var r = t.getBoundingClientRect(), min = +t.min, span = (+t.max - min) || 1;
      var thumbX = r.x + THUMB + ((+t.value - min) / span) * (r.width - 2 * THUMB);
      watch = newWatch(t, e, e.clientX - thumbX);
    }, true);
    document.addEventListener('input', function (e) { if (watch && e.target === watch.input) watch.inputEvents++; }, true);
    document.addEventListener('pointermove', function (e) {
      if (!watch) return;
      var now = performance.now();
      watch.maxGap = Math.max(watch.maxGap, now - watch.last); watch.last = now; watch.moves++; watch.x = e.clientX; watch.y = e.clientY;
      var r = watch.input.getBoundingClientRect(), min = +watch.input.min, span = (+watch.input.max - min) || 1;
      var want = min + Math.min(1, Math.max(0, (e.clientX - watch.offset - r.x - THUMB) / ((r.width - 2 * THUMB) || 1))) * span;
      watch.offBy = Math.abs(+watch.input.value - want) / span;
      if (watch.offBy > 0.08) { if (!watch.offSince) watch.offSince = now; } else watch.offSince = 0;
      if (watch.offSince && now - watch.offSince > 600 && !watch.reported) { watch.reported = true; report('thumb behind the pointer', e); }
    }, true);
    window.setInterval(function () {   // (the pointer may stop moving, with the thumb still behind it)
      if (watch && watch.offSince && !watch.reported && performance.now() - watch.offSince > 600) { watch.reported = true; report('thumb behind a pointer that stopped', null); }
    }, 200);
    document.addEventListener('pointerup', function () { watch = null; }, true);
    document.addEventListener('pointercancel', function (e) { if (watch) { watch.reported = true; report('pointercancel (the browser ended the press)', e); } watch = null; }, true);
  })();
'''


def page_js(diagnostics):
    """The JavaScript to put in a page's Slice script: the drag, and the diagnostics when `diagnostics` is true."""
    return DRAG_JS + (DIAGNOSTICS_JS if diagnostics else '')


def diagnostics_enabled(settings_file):
    """Whether the Troubleshooting setting 'diagnostics_3d' is on in `settings_file` (viewer_settings.json). Off when it is missing or unreadable."""
    try:
        with open(settings_file, 'r', encoding='utf-8') as f:
            on_disk = json.load(f)
        return isinstance(on_disk, dict) and on_disk.get('diagnostics_3d') is True
    except Exception:
        return False
