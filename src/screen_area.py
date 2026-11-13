"""Where on the screens a small window (a tip) can go: the monitor a point is on, and a place beside the pointer kept on it.

Tk only knows the main monitor's size (winfo_screenwidth/height): a window on a second monitor is outside it (at negative x when the
monitor is left of the main one), so a tip kept "on the screen" by those numbers ran off the second monitor's edge or jumped to the main
one. On Windows, Windows is asked which monitor the point is on; elsewhere, the main monitor's size is the best Tk can say.

    left, top, right, bottom = work_area_at(widget, x, y)
    place_beside(tip, x, y, widget, dx=12, dy=18)
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
