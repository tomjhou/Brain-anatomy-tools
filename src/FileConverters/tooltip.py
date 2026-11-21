"""A small note that appears when the pointer rests on a widget: shared by the alignment tool and the section loader."""
import os
import sys
import tkinter as tk

try:
    import screen_area           # (src: tips that come up in the right place and size with monitors of different scalings)
except ImportError:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))     # (after this folder: its modules come first)
    try:
        import screen_area
    except ImportError:
        screen_area = None


def tip_window(master):
    """A new, hidden window for a note, shown by show_tip_window once its contents are in and its place is set. On Windows it is an
    ordinary window with its frame taken off, not Tk's frameless pop-up: with monitors of different scalings, Windows put a pop-up in the
    wrong place and at the wrong size now and then, when the focus was in a window on the other monitor (see screen_area.tip_window)."""
    if screen_area is not None:
        return screen_area.tip_window(master)
    win = tk.Toplevel(master)
    win.withdraw()
    win.wm_overrideredirect(True)
    return win


def show_tip_window(win):
    """Shows a note made by tip_window (if it is hidden); once shown, call it again after its contents change, to fit it to them."""
    if str(win.state()) == "withdrawn":
        if screen_area is not None:
            screen_area.show_tip_window(win)
        else:
            win.deiconify()
    elif screen_area is not None:
        screen_area.fit_tip(win)


def watch_pointer(window):
    """Starts correcting pointer_xy (once per program: see screen_area.watch_pointer)."""
    if screen_area is not None:
        screen_area.watch_pointer(window)


def pointer_xy(widget):
    """The pointer's place on the screen in the windows' own coordinates. winfo_pointerxy can be in another monitor's scaling, while the
    focus is in a window on the monitor of another scaling (see screen_area.pointer_xy)."""
    if screen_area is not None:
        return screen_area.pointer_xy(widget)
    return widget.winfo_pointerxy()


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

    CURSOR_HEIGHT = 28       # (screen pixels below the pointer that the arrow cursor may cover: the note goes below that)

    def __init__(self, widget, text, delay=500, wraplength=260, at_pointer=False):
        self.widget, self.text, self.delay, self.tip, self.job = widget, text, delay, None, None
        self.wraplength = wraplength
        self.at_pointer = at_pointer     # (the note just below the pointer, not below the widget: for a tall widget, e.g. a slider)
        self.cursor = None       # (the widget's own cursor, while the note shows: the plain arrow instead, smaller than a hand)
        widget.bind("<Enter>", self.schedule, add="+")
        widget.bind("<Leave>", self.hide, add="+")
        widget.bind("<ButtonPress>", self.hide, add="+")

    def schedule(self, event=None):
        self.hide()
        self.job = self.widget.after(self.delay, self.show)

    def show(self):
        self.job = None
        text = self.text() if callable(self.text) else self.text
        watch_pointer(self.widget)
        self.tip = tip_window(self.widget)
        tk.Label(self.tip, text=text, justify=tk.LEFT, background="#ffffe0", foreground="black", relief=tk.SOLID, borderwidth=1,
                 wraplength=self.wraplength, padx=4, pady=2).pack()
        try:
            self.cursor = str(self.widget.cget("cursor"))
            self.widget.configure(cursor="arrow")
        except tk.TclError:
            self.cursor = None
        self.place()
        show_tip_window(self.tip)

    def place(self):
        """The note below the widget and the pointer (at_pointer: below the pointer only), on its monitor."""
        if self.at_pointer:
            x, y = pointer_xy(self.widget)
            place_on_monitor(self.tip, x, y + self.CURSOR_HEIGHT, avoid_y=y)
            return
        below = max(self.widget.winfo_rooty() + self.widget.winfo_height() + 4, pointer_xy(self.widget)[1] + self.CURSOR_HEIGHT)
        place_on_monitor(self.tip, self.widget.winfo_rootx(), below, avoid_y=self.widget.winfo_rooty())

    def hide(self, event=None):
        if self.job:
            self.widget.after_cancel(self.job)
            self.job = None
        if self.tip:
            self.tip.destroy()
            self.tip = None
        if self.cursor is not None:
            try:
                self.widget.configure(cursor=self.cursor)
            except tk.TclError:
                pass
            self.cursor = None


class HeadingTips(Tooltip):
    """Notes for some of a ttk.Treeview's column headings: tips {column key: text}, each shown when the pointer rests on that heading
    (a Treeview's headings are not widgets of their own, so a Tooltip can't be put on one)."""

    def __init__(self, tree, tips, delay=500, wraplength=320):
        self.widget, self.tips, self.delay, self.tip, self.job = tree, dict(tips), delay, None, None
        self.wraplength, self.cursor, self.text, self.column = wraplength, None, "", None
        tree.bind("<Motion>", self.moved, add="+")
        tree.bind("<Leave>", self.left, add="+")
        tree.bind("<ButtonPress>", self.hide, add="+")

    def column_at(self, x, y):
        """The key of the column whose heading is at (x, y) in the tree, or None (not on a heading)."""
        if self.widget.identify_region(x, y) != "heading":
            return None
        col = self.widget.identify_column(x)            # ("#2": the second column shown)
        try:
            return self.widget.column(col, "id")
        except tk.TclError:
            return None

    def moved(self, event):
        col = self.column_at(event.x, event.y)
        if col == self.column:
            return
        self.column = col
        self.hide()
        if col in self.tips:
            self.text = self.tips[col]
            self.job = self.widget.after(self.delay, self.show)

    def left(self, event=None):
        self.column = None
        self.hide()

    def place(self):
        """The note below the heading row, at the pointer, on its monitor."""
        x, y = pointer_xy(self.widget)
        place_on_monitor(self.tip, x, y + self.CURSOR_HEIGHT, avoid_y=y)
