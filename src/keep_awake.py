"""Long downloads and installs that carry on when the laptop's lid is closed, and pip runs that do not hang.

    with keep_awake("downloading the atlas"):
        ...                                         # the computer is asked not to sleep, and this program not to be paused, meanwhile

On Windows a power request is made (PowerCreateRequest / PowerSetRequest) for the time inside the block: "system required" (no sleep while
it lasts) and "execution required" (this program keeps running rather than being paused: on a laptop with Modern Standby, closing the lid
with the lid action "Do nothing" still lets Windows go into standby after half a minute or so, and that pauses desktop programs; a
download then sat dead for good). It can be seen in an administrator's  powercfg /requests  while it lasts. If that is not available
(Windows before 7), SetThreadExecutionState is used instead. On macOS  caffeinate -ims -w <this program>  runs for the time inside the
block (no idle sleep, no sleep on mains power; it stops by itself if this program ends). Elsewhere nothing is done. It never raises: being
kept awake is a help, never a reason to fail. Blocks can be nested and used from several threads at once; each holds its own request.

What it cannot do: a lid set to "Sleep" (Windows), or a Mac's lid closed on battery with no external display, sleeps anyway; that is the
computer's own setting, which a program cannot override.

pip: pip_install_command(specs) gives the command line, with network timeouts and retries, and run_pip(...) runs it, kept awake, and starts it
again if it has printed nothing for a long time (a dead connection after the computer slept: pip can wait on it for ever).

Only the standard library is used, so it works before any other library is installed."""
import contextlib
import os
import subprocess
import sys
import threading
import time

# ---- staying awake ---------------------------------------------------------------------------------------------------------------------

_POWER_REQUEST_SYSTEM_REQUIRED = 1
_POWER_REQUEST_EXECUTION_REQUIRED = 3
_POWER_REQUEST_CONTEXT_VERSION = 0
_POWER_REQUEST_CONTEXT_SIMPLE_STRING = 0x1
_ES_CONTINUOUS = 0x80000000
_ES_SYSTEM_REQUIRED = 0x00000001
_thread_state = threading.local()          # (the SetThreadExecutionState fallback is per thread: how many blocks of this thread hold it)


def _windows_power_request(reason):
    """A Windows power request for `reason`, set: (kernel32, handle, [request types set]), or None if it could not be made."""
    import ctypes
    from ctypes import wintypes

    class _Detailed(ctypes.Structure):
        _fields_ = [("LocalizedReasonModule", wintypes.HMODULE), ("LocalizedReasonId", wintypes.ULONG),
                    ("ReasonStringCount", wintypes.ULONG), ("ReasonStrings", ctypes.POINTER(wintypes.LPWSTR))]

    class _Reason(ctypes.Union):
        _fields_ = [("Detailed", _Detailed), ("SimpleReasonString", wintypes.LPWSTR)]

    class _ReasonContext(ctypes.Structure):
        _fields_ = [("Version", wintypes.ULONG), ("Flags", wintypes.DWORD), ("Reason", _Reason)]

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    if not hasattr(kernel32, "PowerCreateRequest"):
        return None
    kernel32.PowerCreateRequest.restype = wintypes.HANDLE
    kernel32.PowerCreateRequest.argtypes = [ctypes.POINTER(_ReasonContext)]
    kernel32.PowerSetRequest.restype = wintypes.BOOL
    kernel32.PowerSetRequest.argtypes = [wintypes.HANDLE, ctypes.c_int]
    kernel32.PowerClearRequest.restype = wintypes.BOOL
    kernel32.PowerClearRequest.argtypes = [wintypes.HANDLE, ctypes.c_int]
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    text = ctypes.create_unicode_buffer(str(reason)[:200] or "working")
    context = _ReasonContext()
    context.Version = _POWER_REQUEST_CONTEXT_VERSION
    context.Flags = _POWER_REQUEST_CONTEXT_SIMPLE_STRING
    context.Reason.SimpleReasonString = ctypes.cast(text, wintypes.LPWSTR)
    handle = kernel32.PowerCreateRequest(ctypes.byref(context))
    if not handle or handle == wintypes.HANDLE(-1).value:
        return None
    done = [kind for kind in (_POWER_REQUEST_SYSTEM_REQUIRED, _POWER_REQUEST_EXECUTION_REQUIRED) if kernel32.PowerSetRequest(handle, kind)]
    if not done:
        kernel32.CloseHandle(handle)
        return None
    return kernel32, handle, done


def _windows_begin(reason):
    try:
        request = _windows_power_request(reason)
        if request is not None:
            return ("request",) + request
    except Exception:
        pass
    try:                                                        # (no power requests here: the older call, for this thread)
        import ctypes
        count = getattr(_thread_state, "count", 0)
        if count or ctypes.windll.kernel32.SetThreadExecutionState(_ES_CONTINUOUS | _ES_SYSTEM_REQUIRED):
            _thread_state.count = count + 1
            return ("thread",)
    except Exception:
        pass
    return None


def _windows_end(token):
    try:
        if token[0] == "request":
            _kind, kernel32, handle, done = token
            for kind in done:
                kernel32.PowerClearRequest(handle, kind)
            kernel32.CloseHandle(handle)
        elif token[0] == "thread":
            _thread_state.count = max(getattr(_thread_state, "count", 1) - 1, 0)
            if _thread_state.count == 0:
                import ctypes
                ctypes.windll.kernel32.SetThreadExecutionState(_ES_CONTINUOUS)
    except Exception:
        pass


def _mac_begin():
    try:
        # -i: no idle sleep; -m: no disk idle sleep; -s: no sleep on mains power. -w: ends by itself when this program does (a crash too).
        process = subprocess.Popen(["/usr/bin/caffeinate", "-ims", "-w", str(os.getpid())],
                                   stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return ("caffeinate", process)
    except Exception:
        return None


def _mac_end(token):
    try:
        process = token[1]
        process.terminate()
        try:
            process.wait(timeout=2)
        except Exception:
            process.kill()
    except Exception:
        pass


def begin(reason="working"):
    """Asks the computer to stay awake (see the module's text) until end(token) is called with what this returns. Never raises."""
    try:
        if sys.platform == "win32":
            return _windows_begin(reason)
        if sys.platform == "darwin":
            return _mac_begin()
    except Exception:
        pass
    return None


def end(token):
    """Lets go of what begin() asked for (nothing for None). Never raises."""
    if token is None:
        return
    try:
        if sys.platform == "win32":
            _windows_end(token)
        elif sys.platform == "darwin":
            _mac_end(token)
    except Exception:
        pass


@contextlib.contextmanager
def keep_awake(reason="working"):
    """The computer stays awake, and this program running, for the time inside the block (see the module's text). Never raises itself."""
    token = begin(reason)
    try:
        yield
    finally:
        end(token)


# ---- pip that does not hang ------------------------------------------------------------------------------------------------------------

PIP_TIMEOUT_S = 30                 # pip's own: a connection that sends nothing for this long is given up on (and tried again)
PIP_RETRIES = 10                   # ...this many times (pip waits a little longer between each)
PIP_QUIET_LIMIT_S = 10 * 60        # pip printed nothing for this long (counting only time the computer was awake): stopped and started again
PIP_ATTEMPTS = 3                   # times pip is started in all, when it has to be stopped for that


def _pip_version():
    try:
        from importlib.metadata import version
        return tuple(int(p) for p in version("pip").split(".")[:2])
    except Exception:
        return (0, 0)


def pip_install_command(specs, extra=()):
    """The command line for installing `specs` with this program's Python: with network timeouts and retries, and (pip 24.1 and later) progress
    as plain lines ("Progress <bytes> of <total>"), so that a long download is not mistaken for a stall."""
    command = [sys.executable, "-m", "pip", "install", "--timeout", str(PIP_TIMEOUT_S), "--retries", str(PIP_RETRIES)]
    if _pip_version() >= (24, 1):
        command += ["--progress-bar", "raw"]
    return command + list(extra) + list(specs)


def _size_text(n):
    return f"{n / 1e6:.1f} MB" if n >= 1e6 else f"{n / 1e3:.0f} kB"


def _stop_tree(process):
    """Stops pip and whatever it started (a build of a package from source runs in its own process)."""
    try:
        if sys.platform == "win32":
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(process.pid)], capture_output=True, timeout=30,
                           creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        else:
            import signal
            os.killpg(process.pid, signal.SIGKILL)
    except Exception:
        pass
    try:
        process.kill()
    except Exception:
        pass
    try:
        process.wait(timeout=10)
    except Exception:
        pass


def run_pip(command, on_line=print, on_idle=None, quiet_limit=PIP_QUIET_LIMIT_S, attempts=PIP_ATTEMPTS, reason="installing Python libraries"):
    """Runs a pip command (see pip_install_command) with the computer kept awake. Each line pip prints is given to on_line (pip's progress lines
    are summed up instead: one line every few seconds); on_idle(), if given, is called about ten times a second meanwhile (a window repainting
    itself). If pip prints nothing for quiet_limit seconds of the computer being awake (the time it was asleep is not counted), it is stopped and
    started again, up to `attempts` times in all. Returns (pip's exit code, the last lines it printed); exit code None if it could not be started."""
    import queue
    tail = []
    creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if sys.platform == "win32" else 0
    with keep_awake(reason):
        for attempt in range(1, max(attempts, 1) + 1):
            try:
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, text=True,
                                           bufsize=1, errors="replace", creationflags=creationflags, start_new_session=sys.platform != "win32")
            except OSError as e:
                return None, [str(e)]
            lines = queue.Queue()

            def drain(stream=process.stdout, out=lines):
                try:
                    for line in stream:
                        out.put(line)
                except Exception:
                    pass
                out.put(None)                                       # (the end of pip's output)

            threading.Thread(target=drain, daemon=True).start()
            quiet, last_tick, last_progress_note = 0.0, time.monotonic(), 0.0
            stalled = False
            while True:
                try:
                    line = lines.get(timeout=0.1)
                except queue.Empty:
                    line = ""
                now = time.monotonic()
                step, last_tick = now - last_tick, now
                if line is None:
                    break
                if line:
                    quiet = 0.0
                    text = line.rstrip()
                    if text.startswith("Progress ") and " of " in text:
                        try:
                            got, total = (int(x) for x in text[len("Progress "):].split(" of "))
                        except ValueError:
                            got, total = None, None
                        if got is not None and (now - last_progress_note > 5 or (total and got >= total)):
                            last_progress_note = now
                            on_line(f"  downloaded {_size_text(got)}" + (f" of {_size_text(total)}" if total else ""))
                    elif text:
                        tail.append(text)
                        del tail[:-40]
                        on_line(text)
                elif step < 5:                                      # (a longer step: the computer was asleep, which is not pip's silence)
                    quiet += step
                if on_idle is not None:
                    try:
                        on_idle()
                    except Exception:
                        pass
                if quiet > quiet_limit:
                    stalled = True
                    break
            if not stalled:
                code = process.wait()
                process.stdout.close()
                return code, tail
            _stop_tree(process)               # (its output pipe is left to the reading thread: closing it under that thread could block)
            message = (f"pip has printed nothing for {max(int(quiet_limit // 60), 1)} minutes (a connection that died while the computer slept?): "
                       + ("starting it again." if attempt < attempts else "giving up."))
            tail.append(message)
            on_line(message)
    return 1, tail


def _main(argv):
    """python keep_awake.py install <pip install arguments...>: pip install, as run_pip does it (for the run-*.bat / run-*.command files)."""
    if len(argv) < 2 or argv[0] != "install":
        print("Usage: python keep_awake.py install <what pip install is given>")
        return 2
    code, _tail = run_pip(pip_install_command(argv[1:]), on_line=lambda text: print(text, flush=True))
    return 1 if code is None else code


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
