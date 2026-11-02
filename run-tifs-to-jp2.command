#!/bin/bash
# Starts the TIFF to JPEG2000 converter on macOS: the counterpart of run-tifs-to-jp2.bat. Double-click it in Finder (the first time, right-click it and
# choose Open: macOS asks before running a script it hasn't seen). Finds Python and the private .venv-macos: see src/ensure_python.sh.
cd "$(dirname "$0")/src/FileConverters" || exit 1
source "../ensure_python.sh"
if [ -z "$PYTHON_CMD" ]; then press_enter_to_close 1; fi

if ! "$PYTHON_CMD" tifs_to_jp2.py "$@"; then
    echo; echo "The TIFF to JPEG2000 converter stopped with an error (see above)."
    press_enter_to_close 1
fi
close_window_and_exit 0
