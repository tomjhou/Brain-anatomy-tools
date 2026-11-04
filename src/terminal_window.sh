# For the .command files that Terminal runs (source it; ensure_python.sh does).
#
# Terminal keeps the window of a finished script open ("[Process completed]") unless its profile says to close it, so pressing Enter at
# "Press Enter to close this window" only printed "logout ... [Process completed]". close_window_and_exit closes the window once the
# shell Terminal started for the script has ended too; a window the script was typed into still has its shell, so it stays open.
#
# hand_off_to_app: an app made by INSTALL_APPS_MAC.command runs its run-*.command in Terminal the first time, to set up Python and the
# libraries. Had the program then run from Terminal too, macOS would have counted its files as Terminal's: on a network folder it could
# list them (the file dialog is macOS's own) but not read them, with no question asked, until the app itself asked the next time. So
# once the setup is done, the run-*.command starts the app again, and the program runs as the app from the very first time.

close_window_and_exit() {   # [exit code]
    local tty
    tty="$(tty 2>/dev/null)"
    if [ "$TERM_PROGRAM" = Apple_Terminal ] && [ -z "$APPS_TEST" ] && [ "${tty#/dev/}" != "$tty" ]; then
        # (in the background, ignoring the hangup that ends the window's shell: it closes the window after this script and that shell end)
        ( trap '' HUP; sleep 1
          osascript -e "tell application \"Terminal\"
              repeat with w in windows
                  if tty of selected tab of w is \"$tty\" and (count of processes of selected tab of w) is 0 then close w
              end repeat
          end tell" ) </dev/null >/dev/null 2>&1 &
    fi
    exit "${1:-0}"
}

press_enter_to_close() {   # [exit code]
    read -r -p "Press Enter to close this window. " _
    close_window_and_exit "${1:-0}"
}

hand_off_to_app() {   # the libraries the app checks for, as an import list (as in INSTALL_APPS_MAC.command)
    # The app leaves its path in this file just before it opens the run-*.command in Terminal (an hour or two is ample for the setup).
    local note app fresh
    note="$HOME/Library/Caches/JhouLab/$(basename "$0" .command).handoff"
    [ -f "$note" ] || return 0
    app="$(cat "$note")"
    fresh="$(find "$note" -mmin -120 2>/dev/null)"
    rm -f "$note"
    [ -n "$fresh" ] && [ -d "$app" ] || return 0
    # Only when the app will find what it checks for: else it would open this window again, and again.
    case "$PYTHON_CMD" in */.venv-macos/bin/python) ;; *) return 0 ;; esac
    "$PYTHON_CMD" -c "import $1" >/dev/null 2>&1 || return 0
    echo
    echo "All set up. Starting $(basename "$app" .app) (macOS may ask it for access to your folders: click Allow)..."
    open "$app" && close_window_and_exit 0
}
