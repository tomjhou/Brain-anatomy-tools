"""Updating the apps from GitHub: the folder they run from is a git clone, so an update is a fetch and a fast-forward pull.

check() and pull() run git (slowly: the network), so call them from a thread; relaunch() starts the app again."""
import os
import subprocess
import sys

TIMEOUT = 120
_FLAGS = getattr(subprocess, "CREATE_NO_WINDOW", 0)     # (no console window flashing up under pythonw)


class UpdateError(Exception):
    """Something to tell the user, in a sentence."""


def _git(root, *args):
    try:
        done = subprocess.run(["git", "-C", root, *args], capture_output=True, text=True, timeout=TIMEOUT, creationflags=_FLAGS,
                              env={**os.environ, "GIT_TERMINAL_PROMPT": "0"})
    except FileNotFoundError:
        raise UpdateError("Git is not installed (or not on the PATH), so the apps cannot update themselves.\n"
                          "Install it from https://git-scm.com/downloads, or download the newest version from GitHub by hand.")
    except subprocess.TimeoutExpired:
        raise UpdateError("Git took too long to answer (is the network down?).")
    if done.returncode != 0:
        raise UpdateError((done.stderr or done.stdout).strip() or f"git {args[0]} failed.")
    return done.stdout.strip()


def repo_root(start=None):
    """The top folder of the clone this file is in; UpdateError if it is not a clone."""
    folder = os.path.dirname(os.path.abspath(start or __file__))
    try:
        return _git(folder, "rev-parse", "--show-toplevel")
    except UpdateError as e:
        if "not a git repository" in str(e).lower():
            raise UpdateError("These apps were not installed from GitHub (this folder is not a git clone), so they cannot update "
                              "themselves.\nDownload the newest version from https://github.com/tomjhou/Brain-anatomy-tools.")
        raise


def check(root):
    """Fetches, and returns the one-line summaries of the commits waiting on GitHub (newest first; empty: up to date)."""
    if not _has_upstream(root):
        raise UpdateError("The current branch does not follow a branch on GitHub, so there is nothing to update from.")
    _git(root, "fetch", "--quiet")
    log = _git(root, "log", "--format=%s", "HEAD..@{u}")
    return log.splitlines() if log else []


def _has_upstream(root):
    try:
        _git(root, "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}")
        return True
    except UpdateError:
        return False


def pull(root):
    """Brings the folder up to date (no merge commits: if the folder has changes of its own that clash, nothing changes and
    UpdateError says so). Returns git's own report."""
    try:
        return _git(root, "pull", "--ff-only", "--quiet") or "Updated."
    except UpdateError as e:
        raise UpdateError("The update could not be applied:\n" + str(e) + "\n\n(Changes made in the apps' folder itself, or commits "
                          "that are not on GitHub, can stop it.)")


def relaunch():
    """Starts this app again with the arguments it was started with (call once its window has closed)."""
    subprocess.Popen([sys.executable, os.path.abspath(sys.argv[0])] + sys.argv[1:],close_fds=True, creationflags=_FLAGS)
