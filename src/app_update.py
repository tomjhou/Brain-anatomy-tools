"""Updating the apps from GitHub, in either of the two ways the folder can have come from:
- a git clone: a fetch and a fast-forward pull (needs git: looked for on the PATH, then where its installers and GitHub Desktop put
  it). Without git, a clone can still be checked (its branch and commit read from .git itself, and GitHub's API asked what is new) but
  not updated: copying files over it would leave git seeing them all as changed, so no_pull_reason() says to install git instead;
- a ZIP download (GitHub's Code > Download ZIP): no git needed. GitHub writes the commit into src/version.txt when it makes the ZIP
  (export-subst in .gitattributes), so the folder knows its version; GitHub's API says what has changed since, and the update
  downloads the newest ZIP and copies its files over the folder's, removing the files the new version no longer has. Nothing else in
  the folder (the .venv, caches, settings, anything added) is touched.

repo_root() says which folder (and the way), check() and pull() use the network (slowly), so call them from a thread; relaunch()
starts the app again."""
import glob
import json
import os
import queue
import re
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


GIT_NEEDED = ("This folder is a git clone, and git was not found on this computer, so the update cannot be installed here.\n\n"
              "Install git from https://git-scm.com/downloads, then check for updates again.")
_found_git = []                                         # ([path or None], once looked for)


def git_exe():
    """The git program: the one on the PATH, else one where Git for Windows, GitHub Desktop or Homebrew install it (apps started from a
    shortcut or the Finder often have a PATH without it); None: not found."""
    if not _found_git:
        _found_git.append(_find_git())
    return _found_git[0]


def _find_git():
    found = shutil.which("git")
    if found:
        return found
    places = []
    if sys.platform == "win32":
        for var in ("ProgramFiles", "ProgramW6432", "ProgramFiles(x86)"):
            if os.environ.get(var):
                places.append(os.path.join(os.environ[var], "Git", "cmd", "git.exe"))
        local = os.environ.get("LOCALAPPDATA")
        if local:
            places.append(os.path.join(local, "Programs", "Git", "cmd", "git.exe"))      # (Git for Windows, installed for one user)
            desktop = glob.glob(os.path.join(local, "GitHubDesktop", "app-*", "resources", "app", "git", "cmd", "git.exe"))
            places += sorted(desktop, key=lambda p: [int(n) for n in re.findall(r"\d+", os.path.basename(p.split("resources")[0].rstrip("\\/")))],
                             reverse=True)                                              # (the newest GitHub Desktop's first)
    else:
        places += ["/opt/homebrew/bin/git", "/usr/local/bin/git",
                   "/Applications/GitHub Desktop.app/Contents/Resources/app/git/bin/git"]
    return next((p for p in places if os.path.isfile(p)), None)


class NoGit(UpdateError):
    """A git clone, and no git to update it with."""


def _git(root, *args):
    exe = git_exe()
    if exe is None:
        raise NoGit(GIT_NEEDED)
    try:
        done = subprocess.run([exe, "-C", root, *args], capture_output=True, text=True, timeout=TIMEOUT, creationflags=_FLAGS,
                              env={**os.environ, "GIT_TERMINAL_PROMPT": "0"})
    except FileNotFoundError:
        raise NoGit(GIT_NEEDED)
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
        return _git(root, "rev-parse", "--show-toplevel") if git_exe() else root
    if not os.path.isdir(os.path.join(root, "src")):
        raise UpdateError(f"The apps' folder was not found (looked in {root}).")
    return root


def check(root):
    """Returns the one-line summaries of the changes waiting on GitHub (newest first; empty: up to date)."""
    if _is_clone(root):
        return _git_check(root) if git_exe() else _clone_check_without_git(root)
    return _zip_check(root)


def no_pull_reason(root):
    """Why pull() cannot update this folder (a git clone, and no git), to say with the changes check() found instead of offering to
    download them; None: it can."""
    return GIT_NEEDED if _is_clone(root) and git_exe() is None else None


def pull(root, inside_start=False):
    """Brings the folder up to date. Returns a report. Only one copy of the apps at a time does it (instance_lock), and not while this copy is
    downloading the atlas: replacing the program files under a running download (or a second copy writing the same files) is what went wrong
    when several updates met. `inside_start`: asked by the start-up's own update step, which already holds the lock."""
    reason = no_pull_reason(root)
    if reason:
        raise NoGit(reason)
    import instance_lock
    lock = instance_lock.maintenance
    if lock.held and not inside_start:
        raise UpdateError("This copy of the app is still downloading the atlas. The update waits until that has finished (it is applied the next "
                          "time you check, or when the app starts).")
    if not lock.acquire("updating the apps"):
        holder = lock.who_holds()
        raise UpdateError("Another copy of the app is updating or downloading right now" + (f" ({holder})" if holder else "") +
                          ". Try again when that has finished.")
    try:
        return _git_pull(root) if _is_clone(root) else _zip_pull(root)
    finally:
        lock.release()


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


def _git_dirs(root):
    """(the clone's own git folder, the folder with its refs and config): the same but for a worktree (whose .git is a file naming its
    folder)."""
    own = os.path.join(root, ".git")
    if os.path.isfile(own):
        with open(own, encoding="utf-8") as f:
            line = f.read().strip()
        if not line.startswith("gitdir:"):
            raise UpdateError("This folder's .git could not be read.")
        own = os.path.normpath(os.path.join(root, line[len("gitdir:"):].strip()))
    common = own
    if os.path.isfile(os.path.join(own, "commondir")):
        with open(os.path.join(own, "commondir"), encoding="utf-8") as f:
            common = os.path.normpath(os.path.join(own, f.read().strip()))
    return own, common


def _read_ref(common, name):
    """The commit a ref (refs/heads/...) points at, from its file or packed-refs; None: neither has it."""
    try:
        with open(os.path.join(common, *name.split("/")), encoding="utf-8") as f:
            return f.read().strip()
    except OSError:
        pass
    try:
        with open(os.path.join(common, "packed-refs"), encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2 and parts[1] == name:
                    return parts[0]
    except OSError:
        pass
    return None


def _read_config(common):
    """{(section, name, key): value} from the clone's config, e.g. ("branch", "main", "remote")."""
    values, section = {}, None
    try:
        with open(os.path.join(common, "config"), encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                m = re.match(r'\[\s*([\w.-]+)(?:\s+"(.*)")?\s*\]$', line)
                if m:
                    section = (m.group(1).lower(), m.group(2))
                elif section and "=" in line and not line.startswith(("#", ";")):
                    key, value = line.split("=", 1)
                    values[section + (key.strip().lower(),)] = value.strip().strip('"')
    except OSError:
        pass
    return values


def _clone_check_without_git(root):
    """check() for a clone where git isn't found: the branch and its commit read from .git, and GitHub's API asked what its branch
    on GitHub has that this commit hasn't (git's HEAD..@{u}, as far as GitHub knows)."""
    own, common = _git_dirs(root)
    try:
        with open(os.path.join(own, "HEAD"), encoding="utf-8") as f:
            head = f.read().strip()
    except OSError:
        raise UpdateError("This folder's .git could not be read.")
    if not head.startswith("ref: refs/heads/"):
        raise UpdateError("This folder (a git clone) is not on a branch, so there is nothing to update from.")
    branch = head[len("ref: refs/heads/"):]
    sha = _read_ref(common, "refs/heads/" + branch)
    config = _read_config(common)
    remote, merge = config.get(("branch", branch, "remote")), config.get(("branch", branch, "merge"), "")
    url = config.get(("remote", remote, "url"), "")
    m = re.search(r"github\.com[:/]+([^/]+/[^/]+?)(?:\.git)?/?$", url)
    if not sha or not merge.startswith("refs/heads/") or not m:
        raise UpdateError("The current branch does not follow a branch on GitHub, so there is nothing to update from.")
    try:
        data = json.loads(_get(f"https://api.github.com/repos/{m.group(1)}/compare/{sha}...{merge[len('refs/heads/'):]}",
                               "the list of changes"))
    except urllib.error.HTTPError:
        raise UpdateError("GitHub does not know this folder's newest commit (made here, and not on GitHub?), so what is new there "
                          "cannot be worked out without git.\n\nInstall git from https://git-scm.com/downloads, then check again.")
    return _summaries(data)


def _summaries(data):
    """The one-line summaries (newest first, as git log) of the commits in GitHub's comparison."""
    commits = data.get("commits") or []
    if data.get("status") == "identical" or not commits:
        return []
    lines = [c["commit"]["message"].splitlines()[0] for c in reversed(commits)]
    more = data.get("ahead_by", len(commits)) - len(commits)                       # (the API lists 250 at most)
    return lines + [f"(and {more} older ones)"] * (more > 0)


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
    return _summaries(data)


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
    open), and any failure (offline, rate limit) is silent: the Check for updates button (or menu) reports those when it is used."""
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


# --- the setting: update automatically (one for all the apps on this computer)

SETTINGS_PATH = os.path.join(os.path.expanduser("~"), ".brain_anatomy_tools_settings.json")


def _read_settings():
    try:
        with open(SETTINGS_PATH, encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def auto_update_enabled():
    """The "Download updates automatically at start-up, then restart" setting, shared by the apps (checked unless it was cleared)."""
    return bool(_read_settings().get("auto_update", True))


def set_auto_update(on):
    """Saves the setting (best effort: a folder that cannot be written to leaves it as it was)."""
    data = _read_settings()
    data["auto_update"] = bool(on)
    try:
        with open(SETTINGS_PATH, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    except OSError:
        pass


class AutoUpdate:
    """With the setting on, at launch, once watch_for_update has found changes: downloads them without a dialog, then restarts the app. Every
    app gives it its own hooks, so that the same logic serves them all:
      busy()        -> True while the user has started an update of their own (nothing is done then);
      blocked()     -> True while a dialog has the input (the update waits: it would end the dialog's window);
      status(text)  -> shows `text` where the app says an update is available (a button, a menu), or restores the usual text for None;
      keep_open()   -> None, or why the app must not be restarted now (the atlas is downloading): it stays as it is, the new version runs at the next start;
      restart()     -> ends this run and starts the app again.
    A class with methods that are scheduled again, not a function that schedules itself: such a function is a reference cycle that keeps the window
    alive after it is closed (see UpdateWatch)."""

    def __init__(self, root, hooks):
        self.root, self.hooks = root, hooks
        self.thread = None
        self.result = queue.Queue()
        self.pulled = False

    def start(self):
        self.root.after(1000, self.step)

    def step(self):
        try:
            if not self.root.winfo_exists():
                return                                          # (the window was closed meanwhile: nothing to do)
            hooks = self.hooks
            if self.thread is None:
                if hooks["busy"]():
                    return                                      # (the user's own check or update is under way)
                if hooks["blocked"]():
                    self.root.after(1000, self.step)
                    return
                reason = no_pull_reason(repo_root())
                if reason:                                      # (a clone with no git: the menu or button says what is new, and why it isn't installed)
                    print("Updates are available, but not installed automatically: " + reason.replace("\n\n", " "))
                    return
                hooks["status"]("Updating...")
                self.thread = threading.Thread(target=self.work, daemon=True)
                self.thread.start()
            elif not self.pulled:
                try:
                    ok, value = self.result.get_nowait()
                except queue.Empty:
                    self.root.after(300, self.step)
                    return
                if not ok:
                    print(f"The automatic update could not be applied: {value}")
                    hooks["status"](None)
                    return
                self.pulled = True
            if self.pulled:
                if hooks["blocked"]():
                    self.root.after(1000, self.step)            # (a dialog came up while it was downloading: restarting waits for it)
                    return
                reason = hooks["keep_open"]()
                if reason:
                    hooks["status"]("Updated: restart to use it")
                    print(f"Updated. Not restarting now: {reason}. The new version runs the next time the app starts.")
                    return
                hooks["restart"]()
                return
            self.root.after(300, self.step)
        except Exception as e:
            print(f"The automatic update stopped: {type(e).__name__}: {e}")

    def work(self):
        try:
            self.result.put((True, pull(repo_root())))
        except UpdateError as e:
            self.result.put((False, str(e)))
        except Exception as e:
            self.result.put((False, f"{type(e).__name__}: {e}"))
