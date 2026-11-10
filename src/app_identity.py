"""Gives a program run from the source its own icon and taskbar button on Windows (and its own Dock icon on macOS), so that it shows as
itself rather than as Python.

Windows groups taskbar buttons by an "app ID". Every Python program shares python.exe's, so a running viewer would sit under a Python
button beside, not on, its pinned shortcut. Each program declares its own ID here (for the process, and on each of its windows), and
install_shortcuts.ps1 stamps the same ID on the shortcuts it makes (keep APPS below and $apps there in step). Each Tk window gets the
program's icon from the icons folder.

On macOS there is no such ID: the Dock shows a running program as the app it belongs to, here Python. Each program only sets the Dock
icon to its own (from the .png in the icons folder), so the running program looks right, though still named "Python".

For a program started with pythonw.exe (no console window), or on macOS from its app made by INSTALL_APPS_MAC.command (which sets
JHOULAB_WINDOWLESS=1), output and errors would vanish: use_windowless_log sends them to a log file and shows errors in a message box."""
import os
import sys
import traceback

ICON_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "icons")
APPS = {   # program: (app ID, icon files' name: .ico for Windows, .png for macOS)
    "xlef_viewer": ("JhouLab.BrainAnatomyTools.XlefViewer", "xlef_viewer"),
    "ccf_browser": ("JhouLab.BrainAnatomyTools.CCFBrowser", "ccf_browser"),
    "abc_atlas_browser": ("JhouLab.BrainAnatomyTools.ABCAtlasBrowser", "abc_atlas_browser"),
    "section_assembler": ("JhouLab.BrainAnatomyTools.SectionAssembler", "section_assembler"),
}


def use_app_identity(program):
    """Call before the program's first window. Not for the packaged .exe: its pinned shortcut has no app ID, so declaring one would split
    the running window from the pinned button, and the .exe has its own icon."""
    if getattr(sys, "frozen", False):
        return
    app_id, icon = APPS[program]
    if sys.platform == "darwin":
        _use_dock_icon(os.path.join(ICON_DIR, icon + ".png"))
        return
    if sys.platform != "win32":
        return
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
    except (AttributeError, OSError):
        pass
    icon_path = os.path.join(ICON_DIR, icon + ".ico")
    if not os.path.exists(icon_path):
        icon_path = None
    import tkinter

    def stamp(window):
        # The process-wide ID above isn't enough on its own: a window started from the pinned shortcut still got a Python button of its
        # own until the window itself carried the ID. So each window gets it when it is shown (Tk makes the outer window only then).
        try:
            _set_window_app_id(int(window.wm_frame(), 16), app_id)
        except (tkinter.TclError, ValueError, OSError):
            pass

    def patch(cls, set_icon):
        init = cls.__init__

        def init_with_identity(self, *args, **kwargs):
            init(self, *args, **kwargs)
            if set_icon and icon_path:   # (on a Tk root, "default" covers every window in it)
                try:
                    self.iconbitmap(default=icon_path)
                except tkinter.TclError:
                    pass
            self.bind("<Map>", lambda event: stamp(self) if event.widget is self else None, add="+")
        cls.__init__ = init_with_identity
    patch(tkinter.Tk, True)          # every Tk root (some programs make several)
    patch(tkinter.Toplevel, False)   # and every other window, e.g. dialogs


def _use_dock_icon(png):
    """macOS: each Tk root sets the Dock icon (and its windows' icons) to the png."""
    if not os.path.exists(png):
        return
    import tkinter
    init = tkinter.Tk.__init__

    def init_with_icon(self, *args, **kwargs):
        init(self, *args, **kwargs)
        try:
            self._dock_icon = tkinter.PhotoImage(master=self, file=png)   # (kept: Tk drops an image nothing refers to)
            self.iconphoto(True, self._dock_icon)
        except tkinter.TclError:
            pass
    tkinter.Tk.__init__ = init_with_icon


def _set_window_app_id(hwnd, app_id):
    """Puts the app ID on one window (its System.AppUserModel.ID property, through the shell's property store for the window)."""
    import ctypes

    class GUID(ctypes.Structure):
        _fields_ = [("data1", ctypes.c_uint32), ("data2", ctypes.c_uint16), ("data3", ctypes.c_uint16), ("data4", ctypes.c_ubyte * 8)]

    class PROPERTYKEY(ctypes.Structure):
        _fields_ = [("fmtid", GUID), ("pid", ctypes.c_uint32)]

    class PROPVARIANT(ctypes.Structure):
        _fields_ = [("vt", ctypes.c_ushort), ("reserved", ctypes.c_ushort * 3), ("value", ctypes.c_void_p), ("padding", ctypes.c_void_p)]

    def guid(text):
        g = GUID()
        ctypes.windll.ole32.IIDFromString(ctypes.c_wchar_p(text), ctypes.byref(g))
        return g
    iid = guid("{886D8EEB-8CF2-4446-8D02-CDBA1DBDCF99}")   # IPropertyStore
    key = PROPERTYKEY(guid("{9F4C2855-9F79-4B39-A8D0-E1D42DE1D5F3}"), 5)   # System.AppUserModel.ID
    store = ctypes.c_void_p()
    if ctypes.windll.shell32.SHGetPropertyStoreForWindow(ctypes.c_void_p(hwnd), ctypes.byref(iid), ctypes.byref(store)) != 0:
        return
    # its methods: QueryInterface, AddRef, Release (2), GetCount, GetAt, GetValue, SetValue (6), Commit
    methods = ctypes.cast(store, ctypes.POINTER(ctypes.POINTER(ctypes.c_void_p))).contents
    set_value = ctypes.WINFUNCTYPE(ctypes.c_long, ctypes.c_void_p, ctypes.POINTER(PROPERTYKEY), ctypes.POINTER(PROPVARIANT))(methods[6])
    release = ctypes.WINFUNCTYPE(ctypes.c_ulong, ctypes.c_void_p)(methods[2])
    text = ctypes.create_unicode_buffer(app_id)
    value = PROPVARIANT(31, (0, 0, 0), ctypes.addressof(text))   # VT_LPWSTR; SetValue copies it
    try:
        set_value(store, ctypes.byref(key), ctypes.byref(value))
    finally:
        release(store)


LOG_EDITOR = "TextEdit" if sys.platform == "darwin" else "Notepad" if sys.platform == "win32" else "a text editor"


def open_log_folder(path):
    """The log's folder in Windows Explorer (or Finder), the log selected."""
    import subprocess
    path = os.path.normpath(path)
    if sys.platform == "win32":
        subprocess.Popen(f'explorer /select,"{path}"')
    elif sys.platform == "darwin":
        subprocess.Popen(["open", "-R", path])
    else:
        subprocess.Popen(["xdg-open", os.path.dirname(path)])


def open_log_in_editor(path):
    """The log in Notepad (TextEdit on macOS, the default text editor elsewhere)."""
    import subprocess
    if sys.platform == "win32":
        subprocess.Popen(["notepad.exe", path])
    elif sys.platform == "darwin":
        subprocess.Popen(["open", "-e", path])
    else:
        subprocess.Popen(["xdg-open", path])


def show_error(parent, title, message, path):
    """An error message box with buttons that open the log `path`: its folder (the log selected), or the log itself in a text editor."""
    import tkinter as tk
    from tkinter import ttk
    dlg = tk.Toplevel(parent)
    dlg.title(title)
    dlg.resizable(False, False)
    if parent is not None and parent.winfo_viewable():
        dlg.transient(parent)
    body = ttk.Frame(dlg, padding=14)
    body.pack(fill=tk.BOTH, expand=True)
    top = ttk.Frame(body)
    top.pack(fill=tk.X)
    try:
        ttk.Label(top, image="::tk::icons::error").pack(side=tk.LEFT, anchor=tk.N, padx=(0, 12))
    except tk.TclError:
        pass
    ttk.Label(top, text=message, wraplength=520, justify=tk.LEFT).pack(side=tk.LEFT, anchor=tk.W)
    bar = ttk.Frame(body)
    bar.pack(fill=tk.X, pady=(14, 0))
    ok = ttk.Button(bar, text="OK", command=dlg.destroy)
    ok.pack(side=tk.RIGHT, padx=(6, 0))
    ttk.Button(bar, text=f"Open log in {LOG_EDITOR}", command=lambda: open_log_in_editor(path)).pack(side=tk.RIGHT, padx=(6, 0))
    ttk.Button(bar, text="Open log folder", command=lambda: open_log_folder(path)).pack(side=tk.RIGHT, padx=(6, 0))
    ok.focus_set()
    dlg.bind("<Return>", lambda e: dlg.destroy())
    dlg.bind("<Escape>", lambda e: dlg.destroy())
    try:            # (centred over the program's window, on its monitor: left to itself, Windows may put it on the main monitor)
        if parent is not None and parent.winfo_viewable():
            dlg.update_idletasks()
            w, h = dlg.winfo_reqwidth(), dlg.winfo_reqheight()
            x = parent.winfo_rootx() + (parent.winfo_width() - w) // 2
            y = parent.winfo_rooty() + max(0, (parent.winfo_height() - h) // 3)
            dlg.geometry(f"+{x}+{y}")
    except tk.TclError:
        pass
    dlg.lift()
    dlg.attributes("-topmost", True)
    dlg.after(300, lambda: dlg.winfo_exists() and dlg.attributes("-topmost", False))
    dlg.grab_set()
    dlg.wait_window()
    return dlg


def use_windowless_log(program, title):
    """Started with pythonw.exe there is no console (nor from a macOS app, which says so with JHOULAB_WINDOWLESS=1): send what would have
    been printed to a log file, and show errors (uncaught ones, and the first in a Tk callback) in a message box that names the log.
    Returns the log's path, or None when there is a console."""
    if sys.stdout is not None and sys.stderr is not None and os.environ.get("JHOULAB_WINDOWLESS") != "1":
        return None
    if sys.platform == "darwin":
        folder = os.path.join(os.path.expanduser("~"), "Library", "Logs", "JhouLab")
    else:
        folder = os.path.join(os.environ.get("LOCALAPPDATA") or os.path.expanduser("~"), "JhouLab", "logs")
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, program + ".log")
    sys.stdout = sys.stderr = open(path, "w", encoding="utf-8", buffering=1)
    shown = []

    def tell(kind, value, tb):
        traceback.print_exception(kind, value, tb)
        if shown:
            return
        shown.append(True)
        try:
            import tkinter
            own = None if tkinter._default_root is not None else tkinter.Tk()
            if own is not None:
                own.withdraw()
            show_error(tkinter._default_root, title, f"{title} hit an error:\n\n{kind.__name__}: {value}\n\nThe details are in\n{path}",
                       path)
            if own is not None:
                own.destroy()
        except Exception:
            pass

    sys.excepthook = tell
    import tkinter
    tkinter.Tk.report_callback_exception = lambda self, kind, value, tb: tell(kind, value, tb)
    return path
