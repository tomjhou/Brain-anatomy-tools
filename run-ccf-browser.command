#!/bin/bash
# Starts the CCF browser on macOS: the counterpart of run-ccf-browser.bat. Double-click it in Finder (the first time, right-click it and
# choose Open: macOS asks before running a script it hasn't seen). Finds Python and the private .venv-macos: see src/ensure_python.sh.
cd "$(dirname "$0")/src" || exit 1
source "./ensure_python.sh"
if [ -z "$PYTHON_CMD" ]; then press_enter_to_close 1; fi

# (CCF_browser.py installs these itself too; here so that its app can take over once they are in: see hand_off_to_app)
if ! "$PYTHON_CMD" -c "import numpy, pandas, PIL" >/dev/null 2>&1; then
    echo "Installing the libraries it needs (numpy, pandas, pillow)..."
    "$PYTHON_CMD" -m pip install --only-binary :all: numpy pandas pillow
fi
hand_off_to_app "numpy, pandas, PIL"     # (started by its app to set up: the app runs it from here: see terminal_window.sh)

if ! "$PYTHON_CMD" CCF_browser.py "$@"; then
    echo; echo "The CCF browser stopped with an error (see above)."
    press_enter_to_close 1
fi
close_window_and_exit 0
