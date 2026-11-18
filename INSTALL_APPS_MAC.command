#!/bin/bash
# macOS: makes an app, with its own icon, for the xlef viewer, CCF browser, ABC atlas browser and section assembler, in your Applications folder
# (~/Applications), to start them from Launchpad, Spotlight or the Dock. Double-click it in Finder (the first time, right-click it and
# choose Open; on macOS 15 and later that no longer works, and a downloaded folder needs its quarantine flag cleared instead: see
# README_MAC_INSTALL.md). Run it again after moving this folder, or when an icon has changed. UNINSTALL_APPS_MAC.command removes them.
#
# Each app runs this folder's run-*.command file in Terminal the first time (to set up Python and the libraries, which shows its
# progress there); after that it starts the program directly, with no Terminal window, its output in ~/Library/Logs/JhouLab (the ABC
# atlas browser always opens in Terminal: it prints its progress there and asks questions there).

REPO="${INSTALL_REPO:-$(cd "$(dirname "$0")" && pwd)}"   # (tests set INSTALL_REPO, running a copy of this from outside the folder)
source "$REPO/src/terminal_window.sh"     # (press_enter_to_close)
# Not from the Downloads folder (it gets cleared out, and the ABC atlas browser keeps its atlas download, up to about 60 GB, and its results
# next to this folder by default): there, offer to move this folder to a permanent place (a BrainTools folder in the home folder, or one
# chosen) and run this again from there. The atlas download and output folder beside it go along when they are there, on the same disk
# (an instant rename), and the records that name them are put right. (macOS lets a folder be moved while a script runs from it: bash goes
# on reading this file as it is.) Downloads is compared without case, as macOS's disks normally are, and by the real path too.
DOWNLOADS_DIR="${DOWNLOADS_DIR:-$HOME/Downloads}"     # (tests set it, and MOVE_BASE and MOVE_ANSWER: move, cancel or a folder)
lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }
is_in_downloads() {
    local here
    for here in "$1" "$(cd "$1" 2>/dev/null && pwd -P)"; do
        case "$(lower "$here")/" in "$(lower "${DOWNLOADS_DIR%/}")"/*) return 0 ;; esac
    done
    return 1
}
disk_of() { if [ "$(uname)" = Darwin ]; then stat -f '%d' "$1" 2>/dev/null; else stat -c '%d' "$1" 2>/dev/null; fi; }
tell_cancelled() {
    echo "This folder is in your Downloads folder:"
    echo "    $REPO"
    echo
    echo "Files in Downloads are often cleared out, which would remove the tools (and the ABC atlas download and your results, which are"
    echo "kept beside this folder). Please move the whole $(basename "$REPO") folder somewhere permanent, for example a BrainTools folder in"
    echo "your home folder (not an iCloud Drive folder), then run INSTALL_APPS_MAC.command again from there."
    if [ -d "$(dirname "$REPO")/abc_atlas_cache" ]; then
        echo
        echo "Move the abc_atlas_cache folder next to it as well, into the same new folder, so the ABC atlas browser finds its download there."
    fi
    echo
    echo "No apps were made."
    if [ -z "$APPS_TEST" ]; then press_enter_to_close 1; fi
    exit 1
}
replace_in_file() {   # file, old text, new text (the records are JSON: a path in them is written as it is, on a Mac)
    [ -f "$1" ] || return 0
    local content; content="$(cat "$1")"
    printf '%s\n' "${content//"$2"/"$3"}" > "$1"
}
if is_in_downloads "$REPO"; then
    name="$(basename "$REPO")"
    parent="$(dirname "$REPO")"
    base="${MOVE_BASE:-$HOME/BrainTools}"
    while :; do
        dest="$base/$name"; k=2
        while [ -e "$dest" ]; do dest="$base/$name-$k"; k=$((k + 1)); done   # (never into a folder that is there already)
        probe="$base"; while [ ! -d "$probe" ]; do probe="$(dirname "$probe")"; done
        moving=(); text_moving=""; text_staying=""
        for item in abc_atlas_cache:"the ABC atlas download" ABC_atlas_browser_output:"your ABC browser results"; do
            n="${item%%:*}"; what="${item#*:}"
            [ -d "$parent/$n" ] || continue
            if [ "$(disk_of "$parent")" = "$(disk_of "$probe")" ] && [ ! -e "$base/$n" ]; then
                moving+=("$parent/$n"); text_moving="$text_moving
    $parent/$n  ($what)"
            else
                text_staying="$text_staying
    $parent/$n  ($what)"
            fi
        done
        text="This folder is in your Downloads folder:
    $REPO

Files in Downloads are often cleared out, which would remove the tools, and the ABC atlas download (up to about 60 GB) and your results, which are kept beside them.

Move the folder here?
    $dest

The installation then goes on by itself from there, in this same Terminal window: there is nothing to run by hand."
        [ -n "$text_moving" ] && text="$text

Moved along with it (instant: nothing is copied):$text_moving"
        [ -n "$text_staying" ] && text="$text

Left where it is (offered again when you open the ABC atlas browser):$text_staying"
        if [ -n "$MOVE_ANSWER" ]; then
            case "$MOVE_ANSWER" in move|cancel) choice="$MOVE_ANSWER" ;; *) choice=choose ;; esac
        else
            choice="$(osascript - "$text" 2>/dev/null <<'APPLESCRIPT'
on run argv
    set r to display dialog (item 1 of argv) with title "The tools are in your Downloads folder" buttons {"Cancel", "Choose another place…", "Move and install"} default button "Move and install" cancel button "Cancel" with icon 2
    if button returned of r is "Move and install" then return "move"
    return "choose"
end run
APPLESCRIPT
)"
            [ -n "$choice" ] || choice=cancel                              # (Cancel, or no dialog could be shown)
        fi
        [ "$choice" = cancel ] && tell_cancelled
        [ "$choice" = move ] && break
        if [ -n "$MOVE_ANSWER" ]; then picked="$MOVE_ANSWER"; MOVE_ANSWER=move; else
            picked="$(osascript -e 'POSIX path of (choose folder with prompt "Choose where to put the tools folder (not Downloads or iCloud Drive)")' 2>/dev/null)"
        fi
        picked="${picked%/}"
        [ -n "$picked" ] || continue                                       # (the folder picker cancelled: ask again)
        if is_in_downloads "$picked" || case "$(lower "$picked")/" in "$(lower "$REPO")"/*) true ;; *) false ;; esac; then
            [ -n "$MOVE_ANSWER" ] && tell_cancelled
            osascript -e 'display alert "Choose another place" message "That is in Downloads, or in this folder itself."' >/dev/null 2>&1
            continue
        fi
        base="$picked"
    done
    if pgrep -f "$REPO/.venv-macos/bin/" >/dev/null 2>&1; then
        echo "A program from this folder is open, so it can't be moved. Close the Brain anatomy tools, then run INSTALL_APPS_MAC.command again."
        echo "Nothing was moved."
        if [ -z "$APPS_TEST" ]; then press_enter_to_close 1; fi
        exit 1
    fi
    echo "Moving the tools out of Downloads, to"
    echo "    $dest"
    mkdir -p "$base" && mv "$REPO" "$dest" || { echo "The folder could not be moved. Nothing was changed."; [ -z "$APPS_TEST" ] && press_enter_to_close 1; exit 1; }
    echo "  Moved the tools to $dest"
    for old in "${moving[@]}"; do
        if mv "$old" "$base/$(basename "$old")"; then
            echo "  Moved $old to $base/$(basename "$old")"
            # the records that name it: the copy's own, and the atlas location kept per user
            for record in "$dest/src/cache_local/abc_atlas_cache_location.json" "$dest/src/cache_local/recent_output_folders.json" \
                          "$HOME/.brain_anatomy_tools_settings.json"; do
                replace_in_file "$record" "$old" "$base/$(basename "$old")"
            done
        else
            echo "  Could not move $old; it is still there (the ABC atlas browser offers it again)."
        fi
    done
    if [ "$(lower "${parent%/}")" != "$(lower "${DOWNLOADS_DIR%/}")" ]; then    # (the folder it was unzipped into, if now empty)
        rmdir "$parent" 2>/dev/null
    fi
    echo
    echo "The folder is moved. The installation goes on from the new place, here (nothing to run by hand):"
    echo
    INSTALL_REPO="" exec bash "$dest/INSTALL_APPS_MAC.command"     # (from its new place, which it then finds for itself)
fi

APPS_DIR="${APPS_DIR:-$HOME/Applications}"      # (tests set APPS_DIR, and APPS_TEST: see the end)

# Installed already? (As check_previous_install.ps1 on Windows.) From this same folder (run again, by mistake or as a "repair"): say what
# that does and doesn't fix, and ask whether to go on. From another folder that is still there: ask whether this one is to be the
# installed copy; if so, leave a note for the programs, the first of which opened from here offers to bring over the old copy's
# settings, atlas download and results (see src/previous_install.py). An old copy in Downloads is warned about, strongly. Where the tools
# are installed from: the folder the existing apps start them from. The same folder is judged by the folder on disk (device and inode),
# not by the spelling of its path.
NOTE="${NOTE:-$HOME/.brain_anatomy_tools_switch.txt}"   # (tests set it, and INSTALL_ANSWER: yes or no instead of asking)
installed_folder() {
    local f line
    for f in "$APPS_DIR"/{"ABC atlas browser","Leica xlef viewer","CCF browser","Section assembler"}.app/Contents/MacOS/launcher; do
        [ -f "$f" ] || continue
        line="$(grep -m1 '^REPO=' "$f")" || continue
        ( eval "$line"; printf '%s' "$REPO" )              # (the line is printf %q's own quoting, written by make_app below)
        return
    done
}
folder_id() { if [ "$(uname)" = Darwin ]; then stat -f '%d:%i' "$1" 2>/dev/null; else stat -c '%d:%i' "$1" 2>/dev/null; fi; }
ask() {   # title, text, default button (Yes or No), icon (1 note, 2 caution): 0 for Yes
    if [ -n "$INSTALL_ANSWER" ]; then [ "$INSTALL_ANSWER" = yes ]; return; fi
    local answer
    answer="$(osascript - "$1" "$2" "$3" "$4" 2>/dev/null <<'APPLESCRIPT'
on run argv
    set r to display dialog (item 2 of argv) with title (item 1 of argv) buttons {"No", "Yes"} default button (item 3 of argv) with icon ((item 4 of argv) as integer)
    return button returned of r
end run
APPLESCRIPT
)"
    if [ -z "$answer" ]; then                              # (no dialog could be shown: ask here in Terminal)
        echo "$2"; echo
        read -r -p "Yes or no? " answer
    fi
    case "$(lower "$answer")" in y|yes) return 0 ;; *) return 1 ;; esac
}
stop_here() { echo "$1"; if [ -z "$APPS_TEST" ]; then press_enter_to_close 0; fi; exit 0; }
old="$(installed_folder)"
if [ -n "$old" ] && [ -d "$old" ]; then
    if [ "$(folder_id "$old")" = "$(folder_id "$REPO")" ]; then
        ask "Already installed from this folder" "The Brain anatomy tools are already installed from this folder:
    $REPO

Running the installer again only remakes the apps in $APPS_DIR. It does not repair the programs themselves. For those:
  - missing or damaged program files: the Update button in one of the programs (or download the tools again);
  - missing Python libraries, or a broken Python environment: delete the .venv-macos folder here; it is rebuilt, with the libraries, when a program next starts.

Remake the apps anyway?" No 1 || stop_here "Already installed from this folder: nothing was changed."
    else
        text="The Brain anatomy tools are already installed from another folder:
    $old

This new folder is:
    $REPO

Make this new folder the installed copy? Its apps replace the old ones, and the first time you open one of the programs from here, it offers to bring over the old copy's settings, its ABC atlas download and your results, so that nothing is downloaded again. The old copy itself is left as it is."
        icon=1
        case "$(lower "$old")/" in "$(lower "${DOWNLOADS_DIR%/}")"/*)
            text="$text

The old copy is in your Downloads folder. Files there are often cleared out, and the tools, the ABC atlas download (up to about 60 GB) and your results kept beside them would go too. Switching to this folder is strongly recommended, and so is moving the atlas download and your results out of Downloads (offered when you open a program from here)."
            icon=2 ;;
        esac
        ask "Already installed from another folder" "$text" Yes "$icon" || stop_here "Nothing was changed: the tools are still installed from $old"
        printf 'previous=%s\nnew=%s\n' "$old" "$REPO" > "$NOTE"
        echo "This folder will be the installed copy (the old one, $old, is left as it is)."
        echo
    fi
fi

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
        echo "# Made by INSTALL_APPS_MAC.command: starts $name from $REPO (run INSTALL_APPS_MAC.command again if that folder moves)."
        printf 'REPO=%q\n' "$REPO"
        printf 'COMMAND="$REPO/%s"\n' "$command"
        echo 'if [ ! -f "$COMMAND" ]; then'
        echo "    osascript -e 'display alert \"$name\" message \"Its folder has moved or is gone. Run INSTALL_APPS_MAC.command in the folder where it is now.\"'"
        echo '    exit 1'
        echo 'fi'
        if [ -n "$libs" ]; then
            echo '# Set up already (the environment, with the libraries the program needs)? Then start it directly, without a Terminal window'
            echo '# (its output goes to ~/Library/Logs/JhouLab: see app_identity.py); else in Terminal, where the setup shows its progress.'
            echo 'PY="$REPO/.venv-macos/bin/python"'
            echo "if [ -x \"\$PY\" ] && \"\$PY\" -c \"import $libs\" >/dev/null 2>&1; then"
            echo "    cd \"\$REPO/src\" && JHOULAB_WINDOWLESS=1 exec \"\$PY\" $program >/dev/null 2>&1"
            echo 'fi'
            echo '# Not yet. Being set up already (the app clicked again while Terminal was busy with it)? Then just show that window: a second'
            echo '# run would start the program from Terminal after the first had started the app, giving two windows.'
            echo 'if pgrep -f "$COMMAND" >/dev/null 2>&1; then exec open -a Terminal; fi'
            echo '# Once the setup is done, the run-*.command starts this app again (it finds this note), so that the program runs as'
            echo '# this app, not as Terminal, from the first time (see hand_off_to_app in src/terminal_window.sh).'
            echo 'mkdir -p "$HOME/Library/Caches/JhouLab"'
            echo "cd \"\$(dirname \"\$0\")/../..\" && pwd > \"\$HOME/Library/Caches/JhouLab/${command%.command}.handoff\""
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
make_app "Section assembler" run-assembler.command section_assembler org.jhoulab.brainanatomytools.sectionassembler \
    FileConverters/SerialAssembler.py ""
chmod +x "$REPO"/run-*.command "$REPO"/extra-tools/run-*.command "$REPO/UNINSTALL_APPS_MAC.command" 2>/dev/null

# A folder from a downloaded zip carries macOS's quarantine flag, and the apps start the run-*.command files in it, each of which would show "could not
# verify" until allowed one by one. Running this installer was the user's decision to trust the folder, so clear the flag on it and on the apps.
if command -v xattr >/dev/null 2>&1; then
    xattr -dr com.apple.quarantine "$REPO" 2>/dev/null
    for app in "$APPS_DIR"/{"Leica xlef viewer","CCF browser","ABC atlas browser","Section assembler"}.app; do
        [ -d "$app" ] && xattr -dr com.apple.quarantine "$app" 2>/dev/null
    done
fi

# Let Launchpad and Spotlight know about them (and their icons) now rather than later
LSREGISTER=/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister
if [ -x "$LSREGISTER" ]; then
    for app in "$APPS_DIR"/{"Leica xlef viewer","CCF browser","ABC atlas browser","Section assembler"}.app; do "$LSREGISTER" -f "$app" >/dev/null 2>&1; done
fi

echo
echo "Done. The apps are in $APPS_DIR (opening it now): drag them to the Dock to keep them there."
echo "They are also in Launchpad and Spotlight."
if [ -z "$APPS_TEST" ]; then                      # (tests set it: no Finder window, no waiting)
    open "$APPS_DIR"
    press_enter_to_close 0
fi
