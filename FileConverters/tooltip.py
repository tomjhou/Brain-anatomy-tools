"""A small note that appears when the pointer rests on a widget: shared by the alignment tool and the section loader."""
import sys
import tkinter as tk


def monitor_area(widget, x, y):
    """(left, top, right, bottom) of the work area (without the taskbar) of the monitor showing screen point (x, y), or the nearest
    one. With several monitors, winfo_screenwidth() is only the main monitor's width, and a monitor to its left has negative x."""
    if sys.platform == "win32":
        try:
            import ctypes
            from ctypes import wintypes

            class MonitorInfo(ctypes.Structure):
                _fields_ = [("cbSize", wintypes.DWORD), ("rcMonitor", wintypes.RECT), ("rcWork", wintypes.RECT),
                            ("dwFlags", wintypes.DWORD)]
            user32 = ctypes.windll.user32
            user32.MonitorFromPoint.restype = wintypes.HANDLE
            user32.MonitorFromPoint.argtypes = [wintypes.POINT, wintypes.DWORD]
            user32.GetMonitorInfoW.argtypes = [wintypes.HANDLE, ctypes.POINTER(MonitorInfo)]
            monitor = user32.MonitorFromPoint(wintypes.POINT(int(x), int(y)), 2)        # (MONITOR_DEFAULTTONEAREST)
            info = MonitorInfo()
            info.cbSize = ctypes.sizeof(info)
            if monitor and user32.GetMonitorInfoW(monitor, ctypes.byref(info)):
                r = info.rcWork
                return r.left, r.top, r.right, r.bottom
        except (OSError, AttributeError, ValueError):
            pass
    vx, vy = widget.winfo_vrootx(), widget.winfo_vrooty()
    return vx, vy, vx + widget.winfo_vrootwidth(), vy + widget.winfo_vrootheight()


def place_on_monitor(win, x, y, avoid_y=None):
    """Put the (overrideredirect) window `win` with its top left at screen point (x, y), kept on that point's monitor: moved left if it
    would stick out on the right; if it would stick out at the bottom, put above avoid_y instead (the pointer, or the widget's top)."""
    win.update_idletasks()
    w, h = win.winfo_reqwidth(), win.winfo_reqheight()
    left, top, right, bottom = monitor_area(win, x, y)
    x = max(left, min(x, right - w - 4))
    if y + h > bottom:
        y = (avoid_y if avoid_y is not None else bottom) - h - 4
    win.wm_geometry(f"+{int(x)}+{int(max(top, y))}")


class Tooltip:
    """A small note that appears when the pointer rests on a widget. `text` may be a function returning the text, called each time the
    note is shown (for a note that depends on the state of the program)."""

    def __init__(self, widget, text, delay=500, wraplength=260):
        self.widget, self.text, self.delay, self.tip, self.job = widget, text, delay, None, None
        self.wraplength = wraplength
        widget.bind("<Enter>", self.schedule, add="+")
        widget.bind("<Leave>", self.hide, add="+")
        widget.bind("<ButtonPress>", self.hide, add="+")

    def schedule(self, event=None):
        self.hide()
        self.job = self.widget.after(self.delay, self.show)

    def show(self):
        self.job = None
        text = self.text() if callable(self.text) else self.text
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        tk.Label(self.tip, text=text, justify=tk.LEFT, background="#ffffe0", foreground="black", relief=tk.SOLID, borderwidth=1,
                 wraplength=self.wraplength, padx=4, pady=2).pack()
        place_on_monitor(self.tip, self.widget.winfo_rootx(), self.widget.winfo_rooty() + self.widget.winfo_height() + 4,
                         avoid_y=self.widget.winfo_rooty())         # (below the widget, on its monitor)

    def hide(self, event=None):
        if self.job:
            self.widget.after_cancel(self.job)
            self.job = None
        if self.tip:
            self.tip.destroy()
            self.tip = None
