"""A message box of the app's own, in the middle of the screen, for when there is no window of the app to put it over.

Tk's native message box is placed over its parent window. With no window open yet (at start-up, say) the parent is a hidden one, and on
macOS that sits at the screen's top left corner, so the box came up partly off the left of the screen. This one measures itself and is
placed in the middle of the screen, whatever the platform. Only tkinter is needed, so it works even when other libraries are missing.

    choice = show("Title", "Message", buttons=("Yes", "No"))       # the button clicked, or None if the box was closed
"""
import tkinter as tk
from tkinter import ttk

WRAP_PX = 520


def show(title, message, buttons=("OK",)):
    """Shows `message` with the given buttons (the first is the default: Return clicks it; Escape or closing the box gives None) and waits.
    Returns the text of the button clicked, or None."""
    root = tk.Tk()
    root.withdraw()                                     # (until it is in place: no flash of it in the corner)
    root.title(title)
    root.resizable(False, False)
    try:
        root.attributes("-topmost", True)
    except tk.TclError:
        pass
    choice = tk.StringVar(root, "")
    closed = {"none": False}

    def pick(text):
        choice.set(text)

    def close():
        closed["none"] = True
        choice.set("\0")

    frame = ttk.Frame(root, padding=(20, 16, 20, 14))
    frame.pack(fill=tk.BOTH, expand=True)
    ttk.Label(frame, text=message, wraplength=WRAP_PX, justify=tk.LEFT).pack(anchor="w")
    row = ttk.Frame(frame)
    row.pack(anchor="e", pady=(16, 0))
    for k, text in enumerate(buttons):
        b = ttk.Button(row, text=text, command=lambda t=text: pick(t), default="active" if k == 0 else "normal")
        b.pack(side=tk.LEFT, padx=(0 if k == 0 else 8, 0))
        if k == 0:
            first = b
    root.bind("<Return>", lambda e: pick(buttons[0]))
    root.bind("<Escape>", lambda e: close())
    root.protocol("WM_DELETE_WINDOW", close)
    # in the middle of the screen (a little above it, as dialogs are), measured first
    root.update_idletasks()
    w, h = root.winfo_reqwidth(), root.winfo_reqheight()
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry(f"+{max((sw - w) // 2, 0)}+{max((sh - h) // 3, 0)}")
    root.deiconify()
    root.lift()
    try:
        root.focus_force()
        first.focus_set()
    except tk.TclError:
        pass
    root.wait_variable(choice)
    result = None if closed["none"] else choice.get()
    root.destroy()
    return result


def info(title, message):
    show(title, message)


def ask_yes_no(title, message):
    """True for Yes; False for No, or if the box was closed."""
    return show(title, message, buttons=("Yes", "No")) == "Yes"
