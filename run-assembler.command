#!/bin/bash
# Starts the section assembler on macOS: the counterpart of run-assembler.bat. Double-click it in Finder (the first time, right-click it and
# choose Open: macOS asks before running a script it hasn't seen). Finds Python and the private .venv-macos: see src/ensure_python.sh.
cd "$(dirname "$0")/src/FileConverters" || exit 1
source "../ensure_python.sh"
if [ -z "$PYTHON_CMD" ]; then read -r -p "Press Enter to close this window. " _; exit 1; fi

if ! "$PYTHON_CMD" SerialAssembler.py "$@"; then
    echo; echo "The section assembler stopped with an error (see above)."
    read -r -p "Press Enter to close this window. " _
fi
