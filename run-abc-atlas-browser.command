#!/bin/bash
# Starts the ABC atlas browser on macOS: the counterpart of run-abc-atlas-browser.bat. Double-click it in Finder (the first time, right-click it and
# choose Open: macOS asks before running a script it hasn't seen). Finds Python and the private .venv-macos: see src/ensure_python.sh.
cd "$(dirname "$0")/src" || exit 1
source "./ensure_python.sh"
if [ -z "$PYTHON_CMD" ]; then read -r -p "Press Enter to close this window. " _; exit 1; fi

if ! "$PYTHON_CMD" ABC_atlas_browser.py "$@"; then
    echo; echo "The ABC atlas browser stopped with an error (see above)."
    read -r -p "Press Enter to close this window. " _
fi
