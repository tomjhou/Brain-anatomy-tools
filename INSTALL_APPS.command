#!/bin/bash
# macOS: makes an app, with its own icon, for the xlef viewer, CCF browser and ABC atlas browser, in your Applications folder
# (~/Applications), to start them from Launchpad, Spotlight or the Dock. Double-click it in Finder (the first time, right-click it and
# choose Open). Run it again after moving this folder, or when an icon has changed. UNINSTALL_APPS.command removes them.
#
# Each app runs this folder's run-*.command file in Terminal the first time (to set up Python and the libraries, which shows its
# progress there); after that it starts the program directly, with no Terminal window, its output in ~/Library/Logs/JhouLab (the ABC
# atlas browser always opens in Terminal: it prints its progress there and asks questions there).

REPO="$(cd "$(dirname "$0")" && pwd)"
APPS_DIR="${APPS_DIR:-$HOME/Applications}"      # (tests set APPS_DIR, and APPS_TEST: see the end)
echo "==== Making apps for the Brain anatomy tools ===="
echo

make_app() {   # name, run-*.command, icon (src/icons/<icon>.icns), bundle id, program, libraries that must import first ("" = always in Terminal)
    local name="$1" command="$2" icon="$3" id="$4" program="$5" libs="$6"
    local app="$APPS_DIR/$name.app"
    rm -rf "$app"
    mkdir -p "$app/Contents/MacOS" "$app/Contents/Resources"
    cp "$REPO/src/icons/$icon.icns" "$app/Contents/Resources/$icon.icns"
    cat > "$app/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key><string>$name</string>
    <key>CFBundleDisplayName</key><string>$name</string>
    <key>CFBundleIdentifier</key><string>$id</string>
    <key>CFBundleExecutable</key><string>launcher</string>
    <key>CFBundleIconFile</key><string>$icon</string>
    <key>CFBundlePackageType</key><string>APPL</string>
    <key>CFBundleShortVersionString</key><string>1.0</string>
    <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
PLIST
    {
        echo '#!/bin/bash'
        echo "# Made by INSTALL_APPS.command: starts $name from $REPO (run INSTALL_APPS.command again if that folder moves)."
        printf 'REPO=%q\n' "$REPO"
        printf 'COMMAND="$REPO/%s"\n' "$command"
        echo 'if [ ! -f "$COMMAND" ]; then'
        echo "    osascript -e 'display alert \"$name\" message \"Its folder has moved or is gone. Run INSTALL_APPS.command in the folder where it is now.\"'"
        echo '    exit 1'
        echo 'fi'
        if [ -n "$libs" ]; then
            echo '# Set up already (the environment, with the libraries the program needs)? Then start it directly, without a Terminal window'
            echo '# (its output goes to ~/Library/Logs/JhouLab: see app_identity.py); else in Terminal, where the setup shows its progress.'
            echo 'PY="$REPO/.venv-macos/bin/python"'
            echo "if [ -x \"\$PY\" ] && \"\$PY\" -c \"import $libs\" >/dev/null 2>&1; then"
            echo "    cd \"\$REPO/src\" && JHOULAB_WINDOWLESS=1 exec \"\$PY\" $program >/dev/null 2>&1"
            echo 'fi'
        fi
        echo 'exec open -a Terminal "$COMMAND"'
    } > "$app/Contents/MacOS/launcher"
    chmod +x "$app/Contents/MacOS/launcher" "$REPO/$command"
    touch "$app"                                  # (so that Finder and the Dock read the new icon)
    echo "  Made $app"
}

mkdir -p "$APPS_DIR" || exit 1
make_app "Leica xlef viewer" run-xlef-viewer.command xlef_viewer org.jhoulab.brainanatomytools.xlefviewer xlef_viewer.py \
    "numpy, tifffile, PIL, imagecodecs"
make_app "CCF browser" run-ccf-browser.command ccf_browser org.jhoulab.brainanatomytools.ccfbrowser CCF_browser.py "numpy, pandas, PIL"
make_app "ABC atlas browser" run-abc-atlas-browser.command abc_atlas_browser org.jhoulab.brainanatomytools.abcatlasbrowser \
    ABC_atlas_browser.py ""
chmod +x "$REPO"/run-*.command "$REPO/UNINSTALL_APPS.command" 2>/dev/null

# Let Launchpad and Spotlight know about them (and their icons) now rather than later
LSREGISTER=/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister
if [ -x "$LSREGISTER" ]; then
    for app in "$APPS_DIR"/{"Leica xlef viewer","CCF browser","ABC atlas browser"}.app; do "$LSREGISTER" -f "$app" >/dev/null 2>&1; done
fi

echo
echo "Done. The apps are in $APPS_DIR (opening it now): drag them to the Dock to keep them there."
echo "They are also in Launchpad and Spotlight."
if [ -z "$APPS_TEST" ]; then                      # (tests set it: no Finder window, no waiting)
    open "$APPS_DIR"
    read -r -p "Press Enter to close this window. " _
fi
