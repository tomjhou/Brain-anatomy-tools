"""Updating the apps from GitHub, in either of the two ways the folder can have come from:
- a git clone: a fetch and a fast-forward pull (needs git);
- a ZIP download (GitHub's Code > Download ZIP): no git needed. GitHub writes the commit into src/version.txt when it makes the ZIP
  (export-subst in .gitattributes), so the folder knows its version; GitHub's API says what has changed since, and the update
  downloads the newest ZIP and copies its files over the folder's, removing the files the new version no longer has. Nothing else in
  the folder (the .venv, caches, settings, anything added) is touched.

repo_root() says which folder (and the way), check() and pull() use the network (slowly), so call them from a thread; relaunch()
starts the app again."""
import json
import os
import queue
import shutil
import ssl
import subprocess
import sys
import threading
import urllib.error
import urllib.request
import zipfile

TIMEOUT = 120
_FLAGS = getattr(subprocess, "CREATE_NO_WINDOW", 0)     # (no console window flashing up under pythonw)
REPO, BRANCH = "tomjhou/Brain-anatomy-tools", "main"
API = f"https://api.github.com/repos/{REPO}"
HOME = f"https://github.com/{REPO}"
VERSION_FILE = os.path.join("src", "version.txt")       # (in the top folder)
STAGING = ".update-download"                            # (in the top folder while an update is unpacked; gone afterwards)
_pending = {}                                           # top folder -> (version installed or None, version to install), from check()


class UpdateError(Exception):
    """Something to tell the user, in a sentence."""


def _git(root, *args):
    try:
        done = subprocess.run(["git", "-C", root, *args], capture_output=True, text=True, timeout=TIMEOUT, creationflags=_FLAGS,
                              env={**os.environ, "GIT_TERMINAL_PROMPT": "0"})
    except FileNotFoundError:
        raise UpdateError("Git is not installed (or not on the PATH), so this folder (a git clone) cannot update itself.\n"
                          "Install it from https://git-scm.com/downloads, or download the newest version from GitHub by hand.")
    except subprocess.TimeoutExpired:
        raise UpdateError("Git took too long to answer (is the network down?).")
    if done.returncode != 0:
        raise UpdateError((done.stderr or done.stdout).strip() or f"git {args[0]} failed.")
    return done.stdout.strip()


def _is_clone(root):
    return os.path.exists(os.path.join(root, ".git"))


def repo_root(start=None):
    """The apps' top folder (above src): a git clone, or a folder unzipped from GitHub."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(start or __file__)))
    if _is_clone(root):
        return _git(root, "rev-parse", "--show-toplevel")
    if not os.path.isdir(os.path.join(root, "src")):
        raise UpdateError(f"The apps' folder was not found (looked in {root}).")
    return root


def check(root):
    """Returns the one-line summaries of the changes waiting on GitHub (newest first; empty: up to date)."""
    return _git_check(root) if _is_clone(root) else _zip_check(root)


def pull(root):
    """Brings the folder up to date. Returns a report."""
    return _git_pull(root) if _is_clone(root) else _zip_pull(root)


def relaunch():
    """Starts this app again with the arguments it was started with (call once its window has closed)."""
    subprocess.Popen([sys.executable, os.path.abspath(sys.argv[0])] + sys.argv[1:],close_fds=True, creationflags=_FLAGS)


# --- a git clone

def _git_check(root):
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


def _git_pull(root):
    """No merge commits: if the folder has changes of its own that clash, nothing changes and UpdateError says so."""
    try:
        return _git(root, "pull", "--ff-only", "--quiet") or "Updated."
    except UpdateError as e:
        raise UpdateError("The update could not be applied:\n" + str(e) + "\n\n(Changes made in the apps' folder itself, or commits "
                          "that are not on GitHub, can stop it.)")


# --- a ZIP download

def installed_version(root):
    """The commit the folder was unzipped from (40 hex digits), or None: a ZIP from before version.txt, or one made by hand."""
    try:
        with open(os.path.join(root, VERSION_FILE), encoding="utf-8") as f:
            sha = f.readline().strip()
    except OSError:
        return None
    return sha if len(sha) == 40 and all(c in "0123456789abcdef" for c in sha) else None


def _context():
    """Certificates to check GitHub's with: certifi's when it is there (Python from python.org on a Mac has none of its own until its
    "Install Certificates" is run), else the system's."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def _get(url, what):
    request = urllib.request.Request(url, headers={"User-Agent": "Brain-anatomy-tools-updater", "Accept": "application/vnd.github+json"})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT, context=_context()) as r:
            return r.read()
    except urllib.error.HTTPError as e:
        if e.code == 403 and e.headers.get("X-RateLimit-Remaining") == "0":
            raise UpdateError("GitHub is limiting how often it can be asked (about 60 times an hour from one network). Try again later.")
        if e.code == 404:
            raise
        raise UpdateError(f"GitHub could not give {what} (error {e.code}).")
    except (urllib.error.URLError, OSError) as e:
        reason = getattr(e, "reason", e)
        raise UpdateError(f"GitHub could not be reached for {what} ({reason}). Is the computer online?")


def _api(path, what):
    return json.loads(_get(f"{API}/{path}", what))


def _zip_check(root):
    old = installed_version(root)
    if old is not None:
        try:
            data = _api(f"compare/{old}...{BRANCH}", "the list of changes")
        except urllib.error.HTTPError:          # (404: a version GitHub doesn't know, e.g. from a fork: update to the newest anyway)
            old = None
    if old is None:
        newest = _api(f"commits/{BRANCH}", "the newest version")
        _pending[root] = (None, newest["sha"])
        date = newest["commit"]["committer"]["date"][:10]
        return [f"(The version in this folder is not known, so all of it is replaced by the newest, of {date}.) "
                + newest["commit"]["message"].splitlines()[0]]
    commits = data.get("commits") or []
    if data.get("status") == "identical" or not commits:
        return []
    _pending[root] = (old, commits[-1]["sha"])
    lines = [c["commit"]["message"].splitlines()[0] for c in reversed(commits)]     # (newest first, as git log)
    more = data.get("ahead_by", len(commits)) - len(commits)                       # (the API lists 250 at most)
    return lines + [f"(and {more} older ones)"] * (more > 0)


def _zip_pull(root):
    if root not in _pending:
        _zip_check(root)
    old, new = _pending.get(root, (None, None))
    if new is None:
        return "Already up to date."
    staging = os.path.join(root, STAGING)
    shutil.rmtree(staging, ignore_errors=True)
    os.makedirs(staging)
    try:
        archive = os.path.join(staging, "update.zip")
        with open(archive, "wb") as f:
            f.write(_get(f"https://codeload.github.com/{REPO}/zip/{new}", "the new version"))
        files = _unpack(archive, os.path.join(staging, "files"))
        gone = _removed_files(old, files) if old else []
        if installed_version(os.path.join(staging, "files")) != new:
            raise UpdateError("The download from GitHub was not the version asked for, so nothing was changed.")
        for path, mode in files.items():         # (all unpacked and checked first: copied only now, so a failed download changes nothing)
            dest = os.path.join(root, *path.split("/"))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            temp = dest + ".update-new"
            shutil.copyfile(os.path.join(staging, "files", *path.split("/")), temp)
            if mode and os.name != "nt":
                os.chmod(temp, mode)
            os.replace(temp, dest)
        for path in gone:
            dest = os.path.join(root, *path.split("/"))
            if os.path.isfile(dest):
                os.remove(dest)
                _remove_empty_folders(root, os.path.dirname(dest))
    except zipfile.BadZipFile:
        raise UpdateError("The download from GitHub was damaged, so nothing was changed. Try again.")
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    _pending.pop(root, None)
    return f"Updated: {len(files)} files" + (f", {len(gone)} removed" if gone else "") + "."


def _unpack(archive, into):
    """Unpacks GitHub's ZIP (whose files are all in one folder, named after the repository and commit) without that folder.
    Returns {path (with /): Unix permissions, or 0}."""
    files = {}
    with zipfile.ZipFile(archive) as z:
        names = [n for n in z.namelist() if not n.endswith("/")]
        top = names[0].split("/")[0] + "/" if names else ""
        for info in z.infolist():
            if info.is_dir() or not info.filename.startswith(top):
                continue
            path = info.filename[len(top):]
            if not path or path.startswith("/") or ".." in path.split("/"):
                continue
            dest = os.path.join(into, *path.split("/"))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with z.open(info) as src, open(dest, "wb") as out:
                shutil.copyfileobj(src, out)
            files[path] = (info.external_attr >> 16) & 0o777
    if VERSION_FILE.replace(os.sep, "/") not in files:
        raise UpdateError("The download from GitHub did not look like these apps, so nothing was changed.")
    return files


def _removed_files(old, new_files):
    """The files the installed version had that the new one doesn't (renamed or deleted on GitHub): from the old version's list of
    files. Files the folder has that neither version has (the .venv, caches, anything added) are not in it."""
    try:
        tree = _api(f"git/trees/{old}?recursive=1", "the list of the installed version's files")
    except urllib.error.HTTPError:
        return []
    if tree.get("truncated"):                   # (a list too long for GitHub to give whole: leave old files rather than guess)
        return []
    return [t["path"] for t in tree.get("tree", []) if t.get("type") == "blob" and t["path"] not in new_files]


def _remove_empty_folders(root, folder):
    root = os.path.abspath(root)
    while os.path.abspath(folder).startswith(root + os.sep) and os.path.isdir(folder) and not os.listdir(folder):
        os.rmdir(folder)
        folder = os.path.dirname(folder)


def watch_for_update(root, show):
    """At launch: looks on GitHub for newer commits in a background thread, and calls show(summaries) (from the Tk loop) when there are
    some. Nothing waits for it: a slow or unreachable network only keeps that thread waiting (it is a daemon: it never holds the program
    open), and any failure (offline, no git, rate limit) is silent: the Check for updates button (or menu) reports those when it is used."""
    found = queue.Queue()

    def work():
        try:
            found.put(check(repo_root()))
        except Exception:
            found.put(None)

    threading.Thread(target=work, daemon=True).start()
    UpdateWatch(root, show, found).poll_later()


class UpdateWatch:
    """Waits for watch_for_update's answer from the Tk loop. An object with a method, not a function that schedules itself: such a function
    refers to itself, a reference cycle that keeps the window it was given alive after it is closed, until the cyclic collector frees it, on
    whichever thread happens to run it (the analysis's threads do, e.g. pynndescent's gc.collect()), and Tcl ends the whole program then:
    "Tcl_AsyncDelete: async handler deleted by the wrong thread". Here the window's own after-commands are the only reference, and closing
    the window drops them."""

    def __init__(self, root, show, found):
        self.root, self.show, self.found = root, show, found

    def poll_later(self):
        self.root.after(1000, self.poll)

    def poll(self):
        try:
            value = self.found.get_nowait()
        except queue.Empty:
            self.poll_later()
            return
        if value:
            self.show(value)
