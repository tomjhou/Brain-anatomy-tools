#!/bin/bash
# Starts the Leica .xlef viewer on macOS: the counterpart of run-xlef-viewer.bat. Double-click it in Finder (the first time, right-click it and
# choose Open: macOS asks before running a script it hasn't seen). Finds Python and the private .venv-macos: see src/ensure_python.sh.
cd "$(dirname "$0")/src" || exit 1
source "./ensure_python.sh"
if [ -z "$PYTHON_CMD" ]; then read -r -p "Press Enter to close this window. " _; exit 1; fi

if ! "$PYTHON_CMD" -c "import numpy, tifffile, PIL, imagecodecs" >/dev/null 2>&1; then
    echo "Installing the libraries it needs (numpy, tifffile, pillow, imagecodecs)..."
    if ! "$PYTHON_CMD" -m pip install numpy tifffile pillow imagecodecs; then
        echo; echo "The libraries could not be installed (see above)."
        read -r -p "Press Enter to close this window. " _; exit 1
    fi
fi
# tkinterdnd2 (dragging .xlef files onto the window) is optional: the viewer works without it.
"$PYTHON_CMD" -c "import tkinterdnd2" >/dev/null 2>&1 || "$PYTHON_CMD" -m pip install tkinterdnd2 >/dev/null 2>&1

if ! "$PYTHON_CMD" xlef_viewer.py "$@"; then
    echo; echo "The Leica .xlef viewer stopped with an error (see above)."
    read -r -p "Press Enter to close this window. " _
fi
