"""Where on the screens a small window (a tip) can go: the monitor a point is on, and a place beside the pointer kept on it.

Tk only knows the main monitor's size (winfo_screenwidth/height): a window on a second monitor is outside it (at negative x when the
monitor is left of the main one), so a tip kept "on the screen" by those numbers ran off the second monitor's edge or jumped to the main
one. On Windows, Windows is asked which monitor the point is on; elsewhere, the main monitor's size is the best Tk can say.

    left, top, right, bottom = work_area_at(widget, x, y)
    place_beside(tip, x, y, widget, dx=12, dy=18)

A window that is not placed over another (a progress box, a figure window, a dialog with no owner) goes on the monitor the user is working on, not the
main one: track(window) on the launcher and on the project browser remembers where they are (the last place any of them was), home_point() says it
(else the monitor the pointer is on), center_on_home(window) puts a window in the middle of that monitor, and center_over(dialog, owner) puts a dialog
over its owner, kept on the owner's monitor (a position clamped to x >= 0 threw it onto the main monitor from one to its left).
"""
import sys


def work_area_at(widget, x, y):
    """(left, top, right, bottom) of the usable part (less the taskbar) of the monitor with the point x, y (screen pixels) on it, or the
    nearest monitor's."""
    if sys.platform == "win32":
        try:
            import ctypes
            from ctypes import wintypes

            class MONITORINFO(ctypes.Structure):
                _fields_ = [("cbSize", wintypes.DWORD), ("rcMonitor", wintypes.RECT), ("rcWork", wintypes.RECT), ("dwFlags", wintypes.DWORD)]
            user32 = ctypes.windll.user32
            user32.MonitorFromPoint.restype = ctypes.c_void_p
            user32.MonitorFromPoint.argtypes = [wintypes.POINT, wintypes.DWORD]
            user32.GetMonitorInfoW.restype = wintypes.BOOL
            user32.GetMonitorInfoW.argtypes = [ctypes.c_void_p, ctypes.POINTER(MONITORINFO)]
            monitor = user32.MonitorFromPoint(wintypes.POINT(int(x), int(y)), 2)         # (MONITOR_DEFAULTTONEAREST)
            info = MONITORINFO()
            info.cbSize = ctypes.sizeof(MONITORINFO)
            if monitor and user32.GetMonitorInfoW(monitor, ctypes.byref(info)):
                r = info.rcWork
                return r.left, r.top, r.right, r.bottom
        except Exception:
            pass
    return 0, 0, widget.winfo_screenwidth(), widget.winfo_screenheight()


def place_beside(tip, x, y, widget, dx=12, dy=18, gap_up=8):
    """Puts the window `tip` below and to the right of the point x, y (the pointer), kept on its monitor: to the point's left where it
    would run off the right edge, above it where it would run off the bottom (gap_up above), and against the edge where even that doesn't
    fit."""
    tip.update_idletasks()
    w, h = tip.winfo_reqwidth(), tip.winfo_reqheight()
    left, top, right, bottom = work_area_at(widget, x, y)
    px, py = x + dx, y + dy
    if px + w > right:
        px = x - w - dx
    if py + h > bottom:
        py = y - h - gap_up
    px, py = max(left, min(px, right - w)), max(top, min(py, bottom - h))
    tip.wm_geometry(f"+{px}+{py}")


_home = {"point": None}
_pointer = {"dx": 0, "dy": 0, "watching": None}


def _note_motion(event):
    """A mouse event: the pointer's place by the event (in the coordinates the windows have) against what Windows says (GetCursorPos): the difference
    is kept, see pointer_xy."""
    try:
        px, py = event.widget.winfo_pointerxy()
        _pointer["dx"], _pointer["dy"] = event.x_root - px, event.y_root - py
    except Exception:
        pass


def watch_pointer(window):
    """Starts keeping the difference pointer_xy corrects by (once for the program): from every mouse movement over any of its windows."""
    if _pointer["watching"] is None:
        try:
            window.bind_all("<Motion>", _note_motion, add="+")
            _pointer["watching"] = window
        except Exception:
            pass


def tip_window(master):
    """A new top-level window for a tip, hidden until show_tip_window.

    A tip has no title bar or border. Tk makes such a window with wm_overrideredirect, which on Windows is a pop-up window, and with monitors of
    different scalings Windows puts a pop-up in the wrong place and at the wrong size (bigger, displaced) now and then, depending on which window has
    the focus. An ordinary window is placed and scaled right every time, so on Windows the tip is one whose frame is taken off (show_tip_window)."""
    import tkinter as tk
    window = tk.Toplevel(master)
    window.withdraw()
    if sys.platform != "win32":
        window.wm_overrideredirect(True)
    return window


def show_tip_window(window):
    """Shows a tip made by tip_window, once its contents are in (and, usually, its place set): on Windows, with its title bar and border taken off
    and not taking the focus."""
    if sys.platform == "win32":
        try:
            import ctypes
            from ctypes import wintypes
            window.update_idletasks()
            user32 = ctypes.windll.user32
            user32.GetWindowLongPtrW.restype = ctypes.c_ssize_t
            user32.GetWindowLongPtrW.argtypes = [wintypes.HWND, ctypes.c_int]
            user32.SetWindowLongPtrW.restype = ctypes.c_ssize_t
            user32.SetWindowLongPtrW.argtypes = [wintypes.HWND, ctypes.c_int, ctypes.c_ssize_t]
            hwnd = wintypes.HWND(user32.GetParent(window.winfo_id()) or window.winfo_id())
            style = user32.GetWindowLongPtrW(hwnd, -16)
            style &= ~(0xC00000 | 0x40000 | 0x80000 | 0x20000 | 0x10000 | 0x800000)      # caption, thick frame, system menu, minimize, maximize, border
            user32.SetWindowLongPtrW(hwnd, -16, style)
            exstyle = user32.GetWindowLongPtrW(hwnd, -20)
            exstyle = (exstyle | 0x80 | 0x08000000) & ~0x40000                           # tool window (no taskbar button), no activation; not an app window
            user32.SetWindowLongPtrW(hwnd, -20, exstyle)
            user32.SetWindowPos(hwnd, None, 0, 0, 0, 0, 0x0001 | 0x0002 | 0x0004 | 0x0010 | 0x0020)      # NOSIZE | NOMOVE | NOZORDER | NOACTIVATE | FRAMECHANGED
            window._frame_removed = hwnd
            window.bind("<Configure>", lambda e, w=window: fit_tip(w) if e.widget is w else None, add="+")      # (Tk sizes the window again when its contents change)
        except Exception:
            window.wm_overrideredirect(True)
    fit_tip(window)
    window.deiconify()


def fit_tip(window):
    """Makes a tip whose frame was taken off (show_tip_window) as large as its contents ask: Tk still counts the frame it had, and the window came out
    larger than its contents, with a gray margin. Call it again when the contents change size."""
    hwnd = getattr(window, "_frame_removed", None)
    if hwnd is None:
        return
    try:
        import ctypes
        from ctypes import wintypes
        want = (window.winfo_reqwidth(), window.winfo_reqheight())
        rect = wintypes.RECT()
        ctypes.windll.user32.GetWindowRect(hwnd, ctypes.byref(rect))
        if (rect.right - rect.left, rect.bottom - rect.top) == want:
            return                                                   # (as large as it should be: nothing to do, and no loop of Configure events)
        ctypes.windll.user32.SetWindowPos(hwnd, None, 0, 0, want[0], want[1], 0x0002 | 0x0004 | 0x0010)      # NOMOVE | NOZORDER | NOACTIVATE
    except Exception:
        pass


def debug_tip(tip, widget, wanted):
    """With the environment variable ABC_TIP_DEBUG set (to anything), prints what the pointer, the monitors and the tip's window say, once a tip is
    placed: for finding out where a tip that comes up in the wrong place goes wrong. `wanted`: where the tip was asked to go."""
    import os
    if not os.environ.get("ABC_TIP_DEBUG"):
        return
    try:
        tip.update_idletasks()
        lines = [f"[tip] winfo_pointerxy={widget.winfo_pointerxy()} correction=({_pointer['dx']}, {_pointer['dy']}) asked to go to {wanted}; "
                 f"the tip is at ({tip.winfo_rootx()}, {tip.winfo_rooty()}) size {tip.winfo_width()}x{tip.winfo_height()}; "
                 f"Tk screen {widget.winfo_screenwidth()}x{widget.winfo_screenheight()}, Tk dpi {widget.winfo_fpixels('1i'):.0f}"]
        if sys.platform == "win32":
            import ctypes
            from ctypes import wintypes
            user32 = ctypes.windll.user32
            pt = wintypes.POINT()
            user32.GetCursorPos(ctypes.byref(pt))
            phys = wintypes.POINT()
            got = user32.GetPhysicalCursorPos(ctypes.byref(phys))
            info = f"GetCursorPos=({pt.x}, {pt.y}) GetPhysicalCursorPos=({phys.x}, {phys.y})" if got else f"GetCursorPos=({pt.x}, {pt.y})"
            try:
                info += f" window dpi: project={user32.GetDpiForWindow(widget.winfo_toplevel().winfo_id())} tip={user32.GetDpiForWindow(tip.winfo_id())}"
            except Exception:
                pass
            lines.append("[tip] " + info + f"; work area at the pointer {work_area_at(widget, *pointer_xy(widget))}")
        print("\n".join(lines), flush=True)
    except Exception as e:
        print(f"[tip] could not report: {e}", flush=True)


def pointer_xy(widget):
    """The pointer's place on the screen, in the coordinates the windows are placed in (those of winfo_rootx, a geometry, an event's x_root).

    winfo_pointerxy asks Windows for the cursor, and with monitors of different scalings Windows can answer in other coordinates than it gives the
    windows: with the focus in a window on the monitor of another scaling, a tip put by it came up displaced, and in the other monitor's scaling
    (bigger). An event's x_root is in the windows' own coordinates, so the difference between the two, taken at every mouse event, is applied to
    what winfo_pointerxy says between events. Where they agree (one monitor, or same scalings) the difference is nil and nothing changes."""
    x, y = widget.winfo_pointerxy()
    return x + _pointer["dx"], y + _pointer["dy"]


def track(window):
    """Remembers where `window` is (the point in its middle) whenever it is moved or resized: the monitor the user is working on, for the windows that
    are not placed over another."""
    def remember(event=None):
        try:
            if event is not None and event.widget is not window:
                return
            _home["point"] = (window.winfo_rootx() + window.winfo_width() // 2, window.winfo_rooty() + window.winfo_height() // 2)
        except Exception:
            pass
    try:
        window.bind("<Configure>", remember, add="+")
        window.after_idle(remember)
    except Exception:
        pass
    watch_pointer(window)


def home_point(widget=None):
    """A point on the monitor the user is working on: where a tracked window last was, else the pointer, else None (the main monitor)."""
    if _home["point"] is not None:
        return _home["point"]
    try:
        return pointer_xy(widget) if widget is not None else None
    except Exception:
        return None


def home_area(widget):
    """work_area_at for the monitor the user is working on (see home_point)."""
    point = home_point(widget)
    return work_area_at(widget, *point) if point else (0, 0, widget.winfo_screenwidth(), widget.winfo_screenheight())


def center_on_home(window, vertical=0.5):
    """Puts `window` in the middle of the monitor the user is working on (vertical: 0.5 in the middle, 1/3 a little above it, as dialogs are),
    kept inside it."""
    window.update_idletasks()
    w, h = window.winfo_reqwidth(), window.winfo_reqheight()
    left, top, right, bottom = home_area(window)
    x = max(left, min(left + (right - left - w) // 2, right - w))
    y = max(top, min(top + int((bottom - top - h) * vertical), bottom - h))
    window.geometry(f"+{x}+{y}")


def center_over(dialog, owner, vertical=0.5):
    """Puts the dialog in the middle of `owner` (a window), kept on the owner's monitor."""
    dialog.update_idletasks()
    w, h = dialog.winfo_reqwidth(), dialog.winfo_reqheight()
    ox, oy, ow, oh = owner.winfo_rootx(), owner.winfo_rooty(), owner.winfo_width(), owner.winfo_height()
    left, top, right, bottom = work_area_at(owner, ox + ow // 2, oy + oh // 2)
    x = max(left, min(ox + (ow - w) // 2, right - w))
    y = max(top, min(oy + int((oh - h) * vertical), bottom - h))
    dialog.geometry(f"+{x}+{y}")
