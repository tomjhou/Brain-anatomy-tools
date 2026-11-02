"""ColorButton: a button that shows the colors it is given (bg, fg, activebackground, ...), on macOS too.

On Windows it is tk.Button itself. On macOS, tk.Button is drawn by macOS (a light, rounded button) whatever background it is given, but
it does take the foreground: the dark theme's light text, and the white of a toggle that is on, came out on white. There a ColorButton is
a tk.Label that acts as a button: it takes command, highlights while the pointer is over it (state "active", so activebackground),
sinks while pressed, and runs its command when the button is let go over it."""
import sys
import tkinter as tk


class _LabelButton(tk.Label):
    def __init__(self, master=None, command=None, **kw):
        kw.setdefault("relief", tk.RAISED)
        kw.setdefault("bd", 2)
        self._command = command
        super().__init__(master, **kw)
        self._relief = self.cget("relief")      # (the relief it is given: shown again after a press)
        self._pressed = False
        self.bind("<Enter>", lambda e: self._hover(True), add="+")
        self.bind("<Leave>", lambda e: self._hover(False), add="+")
        self.bind("<ButtonPress-1>", self._press, add="+")
        self.bind("<ButtonRelease-1>", self._release, add="+")

    def configure(self, cnf=None, **kw):
        if isinstance(cnf, dict):
            kw = {**cnf, **kw}
        elif cnf is not None:
            return self.cget(cnf) if cnf == "command" else super().configure(cnf)
        if "command" in kw:
            self._command = kw.pop("command")
        if "relief" in kw:
            self._relief = kw["relief"]
            if self._pressed:
                kw["relief"] = tk.SUNKEN
        return super().configure(**kw) if kw else None

    config = configure

    def cget(self, key):
        return self._command if key == "command" else super().cget(key)

    __getitem__ = cget

    def __setitem__(self, key, value):
        self.configure(**{key: value})

    def invoke(self):
        if str(super().cget("state")) != tk.DISABLED and self._command is not None:
            return self._command()

    def _hover(self, inside):
        if str(super().cget("state")) != tk.DISABLED:
            super().configure(state=tk.ACTIVE if inside else tk.NORMAL)

    def _press(self, event):
        if str(super().cget("state")) != tk.DISABLED:
            self._pressed = True
            super().configure(relief=tk.SUNKEN)

    def _release(self, event):
        if not self._pressed:
            return
        self._pressed = False
        super().configure(relief=self._relief)
        if 0 <= event.x < self.winfo_width() and 0 <= event.y < self.winfo_height():
            self.invoke()


ColorButton = _LabelButton if sys.platform == "darwin" else tk.Button
