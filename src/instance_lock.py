"""A lock between copies of the apps (separate programs on one computer), for the jobs that must not run twice at once: updating the apps' own
files, and downloading the atlas. Held with the operating system's file locking, so it is let go of the moment the program ends, whatever way
it ends (a crash, the window closed, a restart): there is never a stale lock to clear.

    lock = instance_lock.Lock("maintenance")
    if lock.acquire():              # (False at once if another copy holds it; wait=seconds waits a while: a restart waits for the copy it replaces)
        try: ...
        finally: lock.release()

Within one program the lock can be taken more than once (each acquire() needs its release()): it is only other copies that are kept out.
who_holds() says what the copy holding it is doing, for the note a copy that did not get it prints."""
import os
import threading
import time

try:
    import msvcrt
except ImportError:                                             # (macOS, Linux)
    msvcrt = None
    import fcntl

FOLDER = os.path.expanduser("~")
LOCK_AT = 1 << 20


class Lock:
    def __init__(self, name, folder=None, path=None):
        self.path = path or os.path.join(folder or FOLDER, f".brain_anatomy_tools_{name}.lock")
        self._file = None
        self._count = 0
        self._guard = threading.Lock()

    def _try_lock(self, f):
        try:
            if msvcrt:
                f.seek(LOCK_AT)                                 # (a byte far past the text: Windows refuses to read a locked byte, and who_holds() reads the text)
                msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            return True
        except OSError:
            return False

    def acquire(self, doing="", wait=0.0):
        """True if this program holds the lock now (it already did, or it was free, possibly after waiting up to `wait` seconds); False if another
        copy holds it. `doing`: what this copy is going to do, kept for who_holds()."""
        deadline = time.monotonic() + wait
        while True:
            with self._guard:
                if self._count:
                    self._count += 1
                    return True
                try:
                    f = open(self.path, "a+")                  # (a+: created if need be, never emptied by a copy that fails to get the lock)
                except OSError:
                    return True                                # (a folder that cannot be written to: no lock, so no protection, but nothing stopped either)
                if self._try_lock(f):
                    self._file, self._count = f, 1
                    try:
                        f.seek(0)
                        f.truncate()
                        f.write(f"{os.getpid()} {doing}\n")
                        f.flush()
                    except OSError:
                        pass
                    return True
                f.close()
            if time.monotonic() >= deadline:
                return False
            time.sleep(0.5)

    def release(self):
        with self._guard:
            if not self._count:
                return
            self._count -= 1
            if self._count == 0 and self._file is not None:
                f, self._file = self._file, None
                try:
                    if msvcrt:
                        f.seek(LOCK_AT)
                        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
                    else:
                        fcntl.flock(f.fileno(), fcntl.LOCK_UN)
                except OSError:
                    pass
                f.close()

    @property
    def held(self):
        return self._count > 0

    def who_holds(self):
        """What the copy holding the lock wrote: 'process number what it is doing', or '' if not known."""
        try:
            with open(self.path, encoding="utf-8", errors="replace") as f:
                return f.read().strip()
        except OSError:
            return ""


maintenance = Lock("maintenance")        # updating the apps' files; downloading the atlas
