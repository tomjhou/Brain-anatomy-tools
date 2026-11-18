# Installing on a Mac

Short steps for installing the Brain anatomy tools on macOS, and what to do when macOS won't open the installer. (The main
[README](README.md) describes the tools themselves.)

## 1. Python

Install **Python 3.13** from [python.org](https://www.python.org/downloads/macos/) (the macOS 64-bit universal2 installer) before the
first run. Python 3.10 to 3.13 with tkinter also works; Homebrew's Python needs its Tk package too:
`brew install python@3.13 python-tk@3.13`.

## 2. Download and unzip

On GitHub, click the green **Code** button, then **Download ZIP**, and unzip it (double-click the zip in Finder).

## 3. If macOS says it "can't verify" the file

Files from a downloaded zip carry a *quarantine* flag, and macOS won't run a script that isn't signed by Apple until you allow it. The
message says it can't verify that the file is free of malware, and the file won't open. Two ways past it:

**Option A (recommended): clear the flag on the whole folder, once.** The apps later start several `run-*.command` files from this folder, and
each would show the same warning, so clear them all at once. Open Terminal, type `xattr -dr com.apple.quarantine ` (with a space after it),
drag the unzipped folder from Finder onto the Terminal window to fill in its path, and press Return:

```
xattr -dr com.apple.quarantine /path/to/the/unzipped/folder
```

Then double-click `INSTALL_APPS_MAC.command`. Nothing is printed when it works.

**Option B: allow the one file in System Settings.**
1. Double-click `INSTALL_APPS_MAC.command`; dismiss the warning.
2. Open **System Settings > Privacy & Security**, scroll down to the message about `INSTALL_APPS_MAC.command`, and click **Open Anyway**.
3. Confirm with your password, then run the installer again and click Open.

On older versions of macOS, right-clicking the file and choosing **Open** does the same. From macOS 15 (Sequoia) on, that no longer
bypasses the check: use Option B, or Option A.

Downloading with `git clone` instead of the zip doesn't set the flag.

## 4. If it says "permission denied"

The zip may have lost the files' executable setting. Give it back to the installer and the launchers:

```
chmod +x /path/to/the/unzipped/folder/*.command /path/to/the/unzipped/folder/extra-tools/*.command
```

## 5. What the installer does

Double-click `INSTALL_APPS_MAC.command`. If the folder is in Downloads (which gets cleared out), it asks whether to move it to a permanent
place first (a `BrainTools` folder in your home folder, or one you choose), together with an ABC atlas download and results folder beside it,
and then installs from there. It makes four apps in `~/Applications`, to start from Launchpad, Spotlight or the Dock. Run it again after
moving the folder.

Once it is running, the installer also clears the quarantine flag on the whole folder and on the apps, so the apps themselves should not
show the "can't verify" warning. (If one still does, use Option A above on the folder.)

The first time each app starts, it opens a Terminal window to set up its own Python environment (`.venv-macos` in this folder) and download
the libraries it needs: a few minutes, with an internet connection. macOS may ask whether the app can access your folders (for example
Documents, or a network drive): click Allow.

To remove the apps, double-click `UNINSTALL_APPS_MAC.command`; it asks before removing anything large and never touches your images,
projects or results.

The less often used launchers (the LIF converter and the TIFF to JPEG2000 converter) are in the `extra-tools` folder.
