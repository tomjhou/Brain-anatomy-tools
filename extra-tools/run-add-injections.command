#!/bin/bash
# Opens the window that adds the surgery viewer's injection sites to an ABC browser 3D page, on macOS: the counterpart of
# run-add-injections.bat. Double-click it in Finder (the first time, right-click it and choose Open: macOS asks before running a
# script it hasn't seen). Finds Python and the private .venv-macos: see src/ensure_python.sh.
cd "$(dirname "$0")/../src" || exit 1
source "./ensure_python.sh"
if [ -z "$PYTHON_CMD" ]; then press_enter_to_close 1; fi

if ! "$PYTHON_CMD" add_injections_3d.py "$@"; then
    echo; echo "It stopped with an error (see above)."
    press_enter_to_close 1
fi
close_window_and_exit 0
