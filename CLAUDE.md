# Notes for working on this repository

## Dialogs open on the monitor that triggered them

Every dialog, popup or new window must open on the same monitor as the window (or widget) that opened it, never jump to the main
monitor. Users have monitors to the left of the main one, where screen x is negative (e.g. -2560 to 0).

- Place a dialog relative to its parent window or the widget clicked (`parent.winfo_rootx() + dx`), then keep it on that monitor
  with the existing helpers: `XlefViewer._keep_on_screen(dlg)`, `tooltip.place_on_monitor` and `tooltip.monitor_area(widget, x, y)`
  (in `src/FileConverters/tooltip.py`), `screen_area.center_over(dialog, owner)` and `screen_area.center_on_home(window)`.
- Never clamp absolute screen coordinates with `max(0, x)` or `max(0, y)`, and never take `winfo_screenwidth()` /
  `winfo_screenheight()` as the screen's bounds: they are only the main monitor's. (`max(0, ...)` on an offset inside the parent
  window is fine.)
- Always pass `parent=` to `messagebox.*`, `filedialog.*`, `colorchooser.askcolor` and `simpledialog.*`; without it Windows may
  put them on the main monitor. To find calls without it, parse the code (`ast`): a line search misses calls that span lines.
- Only a dialog with no window behind it at all (an error at startup, before any window) may fall back to the main monitor.
- To test a placement, fake a widget at negative x (override its `winfo_rootx`) and check that `monitor_area` of where the dialog
  went is that widget's monitor.
