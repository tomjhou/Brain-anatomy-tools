#!/bin/bash
# Starts the LIF converter on macOS: the counterpart of run-lifconverter.bat. Double-click it in Finder (the first time, right-click it and
# choose Open: macOS asks before running a script it hasn't seen). Finds Python and the private .venv-macos: see src/ensure_python.sh.
cd "$(dirname "$0")/src/FileConverters" || exit 1
source "../ensure_python.sh"
if [ -z "$PYTHON_CMD" ]; then read -r -p "Press Enter to close this window. " _; exit 1; fi

if ! "$PYTHON_CMD" -c "import numpy, cv2, PIL" >/dev/null 2>&1; then
    echo "Installing the libraries it needs (numpy, opencv-python, pillow)..."
    if ! "$PYTHON_CMD" -m pip install numpy opencv-python pillow; then
        echo; echo "The libraries could not be installed (see above)."
        read -r -p "Press Enter to close this window. " _; exit 1
    fi
fi

if ! "$PYTHON_CMD" LifConverterGUI.py "$@"; then
    echo; echo "The LIF converter stopped with an error (see above)."
    read -r -p "Press Enter to close this window. " _
fi
