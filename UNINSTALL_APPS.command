#!/bin/bash
# macOS counterpart of UNINSTALL_APPS.bat: removes the apps INSTALL_APPS.command made, and the programs' logs and thumbnails; then asks
# (Enter = no) about the private Python environment (.venv-macos), the ABC atlas download, the browsers' caches and settings, and the
# xlef viewer's settings. Never touched: projects, images, the ABC browser's output folders, this folder's own files, Python itself.
# Apps kept in the Dock stay there as question marks: drag them out of the Dock.

REPO="$(cd "$(dirname "$0")" && pwd)"
APPS_DIR="${APPS_DIR:-$HOME/Applications}"      # (tests set these, and ANSWERS: the answers to the questions, e.g. "ynyn")
LIBRARY="${LIBRARY:-$HOME/Library}"
XLEF_SETTINGS="${XLEF_SETTINGS:-$HOME/.xlef_viewer_settings.json}"
asked=0

ask() {   # question -> 0 (yes) or 1 (no)
    local answer
    if [ -n "$ANSWERS" ]; then
        answer="${ANSWERS:$asked:1}"; asked=$((asked + 1))
        echo "$1 [y/N] $answer"
    else
        read -r -p "$1 [y/N] " answer
    fi
    [[ "$answer" == [yY]* ]]
}

size_of() { du -sh "$1" 2>/dev/null | cut -f1; }

remove() {   # path, what
    if [ -e "$1" ]; then rm -rf "$1" && echo "  Removed $2"; fi
}

echo "==== Uninstalling the Brain anatomy tools ===="
echo
# Programs running from the environment keep files in it (and in the caches) open: they must be closed first.
if pgrep -f "$REPO/.venv-macos/" >/dev/null 2>&1; then
    echo "Some of the tools are still running. Close them (and any Terminal window started by a run-*.command file), then run this again."
    [ -z "$ANSWERS" ] && read -r -p "Press Enter to close this window. " _
    exit 1
fi

echo "Removing the apps, logs and thumbnails:"
for name in "Leica xlef viewer" "CCF browser" "ABC atlas browser"; do      # (only ours: made by INSTALL_APPS.command)
    app="$APPS_DIR/$name.app"
    if grep -q "org.jhoulab.brainanatomytools" "$app/Contents/Info.plist" 2>/dev/null; then remove "$app" "$app"; fi
done
remove "$LIBRARY/Logs/JhouLab" "the logs"
remove "$LIBRARY/Caches/JhouLab/xlef_viewer_thumbs" "the xlef viewer's thumbnails"
rmdir "$LIBRARY/Caches/JhouLab" 2>/dev/null

echo
if [ -d "$REPO/.venv-macos" ]; then
    if ask "Remove the private Python environment ($REPO/.venv-macos, $(size_of "$REPO/.venv-macos"))? A run-*.command file rebuilds it (needs the internet)"; then
        remove "$REPO/.venv-macos" "the .venv-macos"
    fi
fi
# The ABC atlas download: where the browsers were told it is (or the default, beside this folder). Asked before the caches, which hold its location.
atlas="$(dirname "$REPO")/abc_atlas_cache"
where="$REPO/src/cache_local/abc_atlas_cache_location.json"
if [ -f "$where" ]; then
    saved="$(sed -n 's/.*"path": *"\(.*\)".*/\1/p' "$where" | head -1)"
    [ -n "$saved" ] && atlas="$saved"
fi
if [ -d "$atlas" ]; then
    if ask "Remove the ABC atlas download ($atlas, $(size_of "$atlas"))? Getting it back means downloading it again"; then
        remove "$atlas" "the ABC atlas download"
    fi
fi
if [ -e "$REPO/src/cache_local" ] || [ -e "$REPO/src/viewer_settings.json" ]; then
    if ask "Remove the ABC and CCF browsers' caches and settings (src/cache_local and src/viewer_settings.json, $(size_of "$REPO/src/cache_local"))? Rebuilt as needed, slowly"; then
        remove "$REPO/src/cache_local" "cache_local"
        remove "$REPO/src/viewer_settings.json" "viewer_settings.json"
    fi
fi
if [ -f "$XLEF_SETTINGS" ]; then
    if ask "Remove the xlef viewer's settings ($XLEF_SETTINGS: its recent files and preferences)?"; then
        remove "$XLEF_SETTINGS" "the xlef viewer's settings"
    fi
fi

echo
echo "Done. Still to do by hand, if you want them gone:"
echo "  - apps kept in the Dock: drag each out of the Dock"
echo "  - the program files themselves: delete the folder $REPO"
[ -z "$ANSWERS" ] && read -r -p "Press Enter to close this window. " _
exit 0
