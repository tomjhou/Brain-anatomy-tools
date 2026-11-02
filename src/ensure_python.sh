# The macOS counterpart of ensure_python.bat, for the run-*.command files:  source "path/to/src/ensure_python.sh"
# Sets PYTHON_CMD to the python of a private virtual environment, the ".venv-macos" folder in the repository's top folder (not ".venv",
# so that a folder shared with Windows keeps both), built from a Python 3.10 to 3.13 found on this Mac. If there is none, it says where
# to get one and leaves PYTHON_CMD empty: unlike on Windows, installing Python on a Mac needs an administrator's password.

PYTHON_CMD=""
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_venv="$(dirname "$_here")/.venv-macos"

_acceptable() {   # a Python 3.10 to 3.13, 64-bit, with tkinter (Homebrew's Python has it only with the python-tk package)
    "$1" -c "import sys, tkinter; sys.exit(0 if (3,10) <= sys.version_info[:2] <= (3,13) and sys.maxsize > 2**32 else 1)" >/dev/null 2>&1
}

_found=""
for _v in 3.13 3.12 3.11 3.10; do      # python.org's installers, then Homebrew's, then whatever is on PATH
    for _p in "/Library/Frameworks/Python.framework/Versions/$_v/bin/python$_v" "/opt/homebrew/bin/python$_v" "/usr/local/bin/python$_v" \
              "$(command -v "python$_v" 2>/dev/null)"; do
        if [ -n "$_p" ] && [ -x "$_p" ] && _acceptable "$_p"; then _found="$_p"; break 2; fi
    done
done
if [ -z "$_found" ] && command -v python3 >/dev/null 2>&1 && _acceptable python3; then _found="$(command -v python3)"; fi

if [ -z "$_found" ]; then
    echo "No suitable Python (3.10 to 3.13, with tkinter) was found on this Mac."
    echo "Install Python 3.13 from https://www.python.org/downloads/macos/ (the macOS 64-bit universal2 installer), then run this again."
    echo "(Homebrew's Python works too, with its Tk package:  brew install python@3.13 python-tk@3.13)"
else
    # The virtual environment: rebuilt if it no longer starts (its Python was removed or upgraded) or is the wrong version.
    if [ -x "$_venv/bin/python" ] && ! _acceptable "$_venv/bin/python"; then
        echo "The virtual environment in $_venv no longer works with a suitable Python, so it will be rebuilt."
        rm -rf "$_venv"
    fi
    if [ ! -x "$_venv/bin/python" ]; then
        echo "Setting up a private Python environment in $_venv (this happens once and takes a minute)..."
        if "$_found" -m venv "$_venv"; then
            "$_venv/bin/python" -m pip install --quiet --upgrade pip >/dev/null 2>&1
        else
            echo "The virtual environment could not be created, so $_found will be used instead."
            rm -rf "$_venv"
        fi
    fi
    if [ -x "$_venv/bin/python" ]; then PYTHON_CMD="$_venv/bin/python"; else PYTHON_CMD="$_found"; fi
fi
unset _here _venv _found _v _p
