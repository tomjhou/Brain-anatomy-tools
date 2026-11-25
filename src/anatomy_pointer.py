"""Tells other programs on this computer where Brain-anatomy-tools is.

The JhouLab surgery viewer's "Show in 3D" button starts src/add_injections_3d.py from here, but the two repositories can be cloned anywhere, so the viewer
cannot guess the folder. Each Brain-anatomy-tools program therefore writes this folder's path to a small text file in a fixed place when it starts
(app_identity.use_app_identity calls register), and the viewer reads it (python/anatomy_tools.py there, which has to build the same file name). The
file is rewritten only when it is missing or names another folder, so a moved or re-cloned folder corrects it the next time any program is opened here.
Nothing has to be installed or re-run after an update.

JHOULAB_ANATOMY_POINTER names another file (tests). Standard library only; it never raises, since a program must start whether or not this worked."""
import os
import sys
from pathlib import Path

FILE_NAME = 'brain-anatomy-tools.txt'


def pointer_file():
    """Where the pointer file lives: <per-user local data>/JhouLab/brain-anatomy-tools.txt."""
    override = os.environ.get('JHOULAB_ANATOMY_POINTER')
    if override:
        return Path(override)
    if sys.platform == 'win32':
        base = Path(os.environ.get('LOCALAPPDATA') or Path.home())
    elif sys.platform == 'darwin':
        base = Path.home() / 'Library' / 'Application Support'
    else:
        base = Path(os.environ.get('XDG_CONFIG_HOME') or Path.home() / '.config')
    return base / 'JhouLab' / FILE_NAME


def tools_folder():
    """The Brain-anatomy-tools folder this copy of the code is in (the folder above src)."""
    return Path(__file__).resolve().parent.parent


def register():
    """Writes this folder's path to the pointer file, unless it is there already. Returns True when the file now names this folder."""
    try:
        if getattr(sys, 'frozen', False):    # (the packaged programs are not a clone with the scripts)
            return False
        folder = str(tools_folder())
        target = pointer_file()
        try:
            if target.read_text(encoding='utf-8').strip() == folder:
                return True
        except OSError:
            pass
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(folder + '\n', encoding='utf-8')
        return True
    except Exception:
        return False
