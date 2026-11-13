"""Downloads of the atlas files that continue where they stopped.

abc_atlas_access writes a download straight to the file's final name, in one request, with no way to continue: a download cut short (the
program closed or crashed, the network went away) leaves a file that is there but not whole, and the next time it is asked for the library
deletes it and starts from the first byte. For the 50 GB imputed dataset that can be an hour or more lost.

install() replaces the library's S3 download method (for every cache, whoever makes it) with download_file() here:
- a big file (100 MB or more) is fetched in pieces of 32 MB, sixteen at a time, each over its own connection: one connection to the server
  carries only a few MB/s, however fast the link, and sixteen carry several times that. The finished pieces are listed in "<name>.part.json", so
  a later ask does only the pieces that were not finished (and a .part file from the one-stream version is continued: its whole pieces count);
- the file is written to "<name>.part" beside its final name, and only moved to that name, and recorded as downloaded, once it is whole and
  its hash is right: nothing else (the library, the program's checks) ever sees a half-written file under the real name;
- if it stops, the next ask continues from the end of the .part file (an HTTP range request) instead of from the start;
- a network error does not end it: it waits (5 s, then longer, up to a minute) and continues, for up to two hours without a single byte
  arriving; errors that retrying cannot mend (access denied, no such file, no room on the disk) end it with the reason;
- a file cut short by the library's own downloader (under its real name, not recorded as finished) is taken over as the beginning of the
  .part file, since the library writes in order: what was already downloaded is not lost.

Progress for a window: downloaded_bytes(path) (the .part file's size while it is downloading) and status_note(path) (a line to show when it
is waiting for the network, or checking the file, else None). The console shows tqdm's bar as the library did."""
import errno
import hashlib
import json
import os
import queue
import threading
import time
from pathlib import Path

CHUNK = 1 << 20
RETRY_FIRST_S = 5                # the wait after the first failure, doubled each time up to RETRY_MAX_S
RETRY_MAX_S = 60
GIVE_UP_S = 2 * 3600             # no byte for this long: stop (the next ask of the file continues from the same place)
WORKERS = 16                     # connections at once for a big file
SEGMENT = 32 * 1000 * 1000       # bytes in a piece: a piece is what is done again if the program stops in the middle of it
PARALLEL_MIN = 100 * 1000 * 1000  # files smaller than this are fetched in one stream
HASH_ATTEMPTS = 2                # whole downloads whose hash was wrong, before giving up
_locks = {}
_notes = {}
_active = {}                       # path -> {'bytes': ...} for a download of pieces under way in this process
_guard = threading.Lock()


def _size_text(n):
    return f"{n / 1e9:.2f} GB" if n >= 1e9 else f"{n / 1e6:.0f} MB" if n >= 1e6 else f"{n / 1e3:.0f} KB"


def part_path(path):
    path = Path(path)
    return path.with_name(path.name + ".part")


def _state_path(part):
    return Path(str(part) + ".json")


def _read_state(part):
    """The list of finished pieces kept beside a .part file by a download of pieces: {'total', 'segment', 'done': [...]}, or None."""
    try:
        with open(_state_path(part), "r", encoding="utf-8") as f:
            state = json.load(f)
        return {"total": int(state["total"]), "segment": int(state["segment"]), "done": [int(i) for i in state["done"]]}
    except Exception:
        return None


def _write_state(part, total, done):
    target = _state_path(part)
    temporary = target.with_name(target.name + ".tmp")
    try:
        with open(temporary, "w", encoding="utf-8") as f:
            json.dump({"total": total, "segment": SEGMENT, "done": done}, f)
        os.replace(temporary, target)
    except OSError:
        pass                                                    # (only a progress note: the pieces are redone if it is missing)


def downloaded_bytes(path):
    """How much of the file is on disk: while it is being downloaded, the bytes of the finished pieces and those arriving (parallel download), or the
    .part file's size (one stream); else the file's own (0 if neither is there)."""
    active = _active.get(str(path))
    if active is not None:
        return active["bytes"]
    sizes = [0]
    part = part_path(path)
    state = _read_state(part)
    if state is not None:
        seg, total = state["segment"], state["total"]
        sizes.append(sum(min((i + 1) * seg, total) - i * seg for i in state["done"]))
    else:
        try:
            sizes.append(part.stat().st_size)
        except OSError:
            pass
    try:
        sizes.append(Path(path).stat().st_size)
    except OSError:
        pass
    return max(sizes)


def status_note(path):
    """A line saying what the download of `path` is doing when it is not simply downloading ("waiting for the network...", "checking the file..."), else None."""
    return _notes.get(str(path))


def _lock_for(path):
    with _guard:
        return _locks.setdefault(str(path), threading.Lock())


def install():
    """Makes every S3 cache of abc_atlas_access download with download_file() below. Returns True if it did (False: already, or the library is not there)."""
    try:
        from abc_atlas_access.abc_atlas_cache.cloud_cache import S3CloudCache
    except Exception:
        return False
    if getattr(S3CloudCache, "_resumable_installed", False):
        return False
    S3CloudCache._download_file = download_file
    S3CloudCache._resumable_installed = True
    return True


def _fatal(error):
    """Whether retrying cannot help: the server says no (a 4xx other than a timeout or a rate limit), or the disk has no room or no permission."""
    try:
        from botocore.exceptions import ClientError
        if isinstance(error, ClientError):
            status = error.response.get("ResponseMetadata", {}).get("HTTPStatusCode", 0)
            return 400 <= status < 500 and status not in (408, 429)
    except ImportError:
        pass
    if isinstance(error, OSError) and getattr(error, "errno", None) in (errno.ENOSPC, errno.EACCES, errno.EROFS, getattr(errno, "EDQUOT", -1)):
        return True
    return isinstance(error, (MemoryError, KeyboardInterrupt))


def _pause(key, seconds, why, say):
    """Waits `seconds`, one second at a time, with the reason in status_note."""
    end = time.monotonic() + seconds
    while True:
        left = end - time.monotonic()
        if left <= 0:
            break
        _notes[key] = f"waiting for the network ({why}); trying again in {int(left) + 1} s"
        time.sleep(min(1.0, left))
    _notes.pop(key, None)


def _object_size(cache, key, say, gave_up_at, key_text):
    """The size of the object on the server, asking again while the network is down."""
    delay = RETRY_FIRST_S
    while True:
        try:
            response = cache.s3_client.list_objects_v2(Bucket=cache._bucket_name, Prefix=key)
            for item in response.get("Contents", []):
                if item["Key"] == key:
                    return int(item["Size"])
            raise FileNotFoundError(f"{key} is not on the server")
        except Exception as e:
            if _fatal(e) or isinstance(e, FileNotFoundError) or time.monotonic() > gave_up_at:
                raise
            say(f"Could not ask the server about {key.split('/')[-1]} ({type(e).__name__}); trying again in {delay} s.")
            _pause(key_text, delay, type(e).__name__, say)
            delay = min(delay * 2, RETRY_MAX_S)


def _md5(path, key_text, total, name=None):
    """The MD5 of the file. With a `name` and a big file, the console says it is checking and shows a bar (reading 50 GB back takes minutes)."""
    import tqdm
    hasher = hashlib.md5()
    done, last = 0, 0.0
    bar = None
    if name and total >= PARALLEL_MIN:
        tqdm.tqdm.write(f"Checking {name} (it is read back once to make sure it arrived whole)...")
        bar = tqdm.tqdm(desc="checking " + name, total=total, unit_scale=True, unit_divisor=1000., unit="B")
    try:
        with open(path, "rb") as f:
            while True:
                chunk = f.read(8 * CHUNK)
                if not chunk:
                    break
                hasher.update(chunk)
                done += len(chunk)
                if bar is not None:
                    bar.update(len(chunk))
                if time.monotonic() - last > 1.0:
                    last = time.monotonic()
                    _notes[key_text] = f"checking the download ({100 * done / max(total, 1):.0f}%)"
    finally:
        if bar is not None:
            bar.close()
    _notes.pop(key_text, None)
    return hasher.hexdigest()


class _NoRange(Exception):
    """The server ignored a Range header (it sent the whole object): parallel pieces cannot be used."""


def _sequential(self, bar, part, key, total, note_key, say, last_progress_box):
    """Downloads `part` from its end to `total` bytes in one stream, continuing after errors (the way for a small file, or when the server does not do ranges)."""
    have = part.stat().st_size if part.exists() else 0
    if have > total:
        part.unlink()
        have = 0
        bar.n = 0
        bar.refresh()
    delay = RETRY_FIRST_S
    last_progress_box[0] = time.monotonic()
    while have < total:
        try:
            request = {"Bucket": self._bucket_name, "Key": key}
            if have:
                request["Range"] = f"bytes={have}-"
            response = self.s3_client.get_object(**request)
            if have and "ContentRange" not in response:         # (a server that ignores the range sends it all: from the beginning)
                have = 0
                bar.n = 0
                bar.refresh()
            with open(part, "ab" if have else "wb") as out:
                for chunk in response["Body"].iter_chunks(CHUNK):
                    out.write(chunk)
                    have += len(chunk)
                    bar.update(len(chunk))
                    last_progress_box[0] = time.monotonic()
                    delay = RETRY_FIRST_S
                    _notes.pop(note_key, None)
        except Exception as e:
            if _fatal(e):
                raise
            have = part.stat().st_size if part.exists() else 0         # (what is on the disk is where it continues from)
            bar.n = have
            bar.refresh()
            if time.monotonic() - last_progress_box[0] > GIVE_UP_S:
                raise RuntimeError(f"Could not download {key}: no data for {GIVE_UP_S // 3600} hours ({type(e).__name__}: {e}). "
                                   "What was downloaded is kept; asking for the file again continues from it.")
            say(f"The download of {key.split('/')[-1]} was interrupted ({type(e).__name__}); continuing from {_size_text(have)} in {delay} s.")
            _pause(note_key, delay, type(e).__name__, say)
            delay = min(delay * 2, RETRY_MAX_S)
            continue
        have = part.stat().st_size if part.exists() else 0
        if have < total:                                                # (the stream ended early without an error: ask for the rest)
            if time.monotonic() - last_progress_box[0] > GIVE_UP_S:
                raise RuntimeError(f"Could not download {key}: the server sent no more data for {GIVE_UP_S // 3600} hours.")
            time.sleep(1)


def _client_for(cache):
    """The S3 client the pieces use: the library's own allows ten connections at a time (boto's default), so for more a sibling client with room for WORKERS
    (same bucket access, so the same anonymous or signed requests) is made. Anything else (no boto client, a failure) uses the cache's."""
    client = cache.s3_client
    try:
        if type(client).__module__.startswith("botocore"):
            import boto3
            from botocore import UNSIGNED
            from botocore.client import Config
            kwargs = {"max_pool_connections": WORKERS + 2}
            if not getattr(cache, "_auth_required", False):
                kwargs["signature_version"] = UNSIGNED
            return boto3.client("s3", config=Config(**kwargs))
    except Exception:
        pass
    return client


def _parallel(self, bar, part, key, total, note_key, say, last_progress_box, want_hash=False):
    """Downloads `part` in pieces of SEGMENT bytes, WORKERS at a time, each its own connection (one connection to the server runs at a fraction of
    what the link can carry). The pieces that are finished are listed in the .part.json file beside it, so a later ask does only the rest. Raises
    _NoRange if the server does not do ranges. With want_hash, a thread hashes the pieces in order as they land (reading them back from the disk
    cache), so that little is left to read when the last one arrives; returns the file's MD5 then (None if it could not be worked out), else None."""
    count = (total + SEGMENT - 1) // SEGMENT
    client = _client_for(self)

    def span(i):
        return i * SEGMENT, min((i + 1) * SEGMENT, total)

    state = _read_state(part)
    have = part.stat().st_size if part.exists() else 0
    done = set()
    if state is not None and state["total"] == total and state["segment"] == SEGMENT and have >= max((span(i)[1] for i in state["done"]), default=0):
        done = set(state["done"])
    elif state is None and have > 0:                           # (a file from the one-stream downloader: its whole pieces are the beginning)
        done = {i for i in range(count) if span(i)[1] <= have}
    else:
        part.unlink(missing_ok=True)                              # (a leftover that does not fit this file: start over)
        _state_path(part).unlink(missing_ok=True)
    part.touch()
    prog = {"bytes": sum(span(i)[1] - span(i)[0] for i in done)}
    bar.n = prog["bytes"]
    bar.refresh()
    lock = threading.Lock()
    cond = threading.Condition(lock)
    abort = threading.Event()
    errors = []
    pending = queue.Queue()
    for i in range(count):
        if i not in done:
            pending.put(i)
    _active[note_key] = prog
    last_progress_box[0] = time.monotonic()
    hashed = {"bytes": 0, "next": 0, "digest": None}

    def hash_in_order():
        """Hashes piece 0, 1, 2, ... as each is finished (a piece is finished when it is in `done`)."""
        try:
            md5 = hashlib.md5()
            while hashed["next"] < count:
                i = hashed["next"]
                with cond:
                    while i not in done:
                        if abort.is_set():
                            return
                        cond.wait(0.5)
                start, stop = span(i)
                with open(part, "rb") as f:
                    f.seek(start)
                    left = stop - start
                    while left > 0:
                        chunk = f.read(min(left, 8 * CHUNK))
                        if not chunk:
                            return                                  # (the file is shorter than it should be: the caller hashes it itself)
                        md5.update(chunk)
                        left -= len(chunk)
                        hashed["bytes"] += len(chunk)
                hashed["next"] = i + 1
            hashed["digest"] = md5.hexdigest()
        except Exception:
            hashed["digest"] = None

    hash_thread = threading.Thread(target=hash_in_order, daemon=True) if want_hash else None

    def save_state(force=False):
        with lock:
            _write_state(part, total, sorted(done))

    def run_segment(i):
        start, stop = span(i)
        pos, delay = start, RETRY_FIRST_S
        while pos < stop:
            if abort.is_set():
                return False
            try:
                response = client.get_object(Bucket=self._bucket_name, Key=key, Range=f"bytes={pos}-{stop - 1}")
                if "ContentRange" not in response:
                    raise _NoRange()
                with open(part, "r+b") as out:
                    out.seek(pos)
                    for chunk in response["Body"].iter_chunks(CHUNK):
                        if abort.is_set():
                            return False
                        chunk = chunk[:stop - pos]
                        out.write(chunk)
                        pos += len(chunk)
                        with lock:
                            prog["bytes"] += len(chunk)
                            bar.update(len(chunk))
                            last_progress_box[0] = time.monotonic()
                        delay = RETRY_FIRST_S
                        _notes.pop(note_key, None)
                        if pos >= stop:
                            break
            except _NoRange:
                raise
            except Exception as e:
                if _fatal(e):
                    raise
                if time.monotonic() - last_progress_box[0] > GIVE_UP_S:
                    raise RuntimeError(f"Could not download {key}: no data for {GIVE_UP_S // 3600} hours ({type(e).__name__}: {e}). "
                                       "What was downloaded is kept; asking for the file again continues from it.")
                say(f"A part of {key.split('/')[-1]} was interrupted ({type(e).__name__}); continuing in {delay} s.")
                _pause(note_key, delay, type(e).__name__, say)
                delay = min(delay * 2, RETRY_MAX_S)
                continue
            if pos < stop:                                              # (the stream ended early without an error: ask for the rest)
                time.sleep(1)
        return True

    def worker():
        while not abort.is_set():
            try:
                i = pending.get_nowait()
            except queue.Empty:
                return
            try:
                if run_segment(i):
                    with cond:
                        done.add(i)
                        cond.notify_all()
                    save_state()
            except BaseException as e:
                errors.append(e)
                abort.set()
                return

    if done:
        save_state()
        say(f"Continuing {key.split('/')[-1]} from {_size_text(prog['bytes'])} of {_size_text(total)}.")
    threads = [threading.Thread(target=worker, daemon=True) for _ in range(min(WORKERS, max(1, count - len(done))))]
    try:
        if hash_thread is not None:
            hash_thread.start()
        for t in threads:
            t.start()
        for t in threads:
            while t.is_alive():
                t.join(0.5)
    finally:
        with cond:
            abort.set()
            cond.notify_all()
        for t in threads:
            t.join()
        _active.pop(note_key, None)
    if errors:
        raise errors[0]
    if part.stat().st_size != total:
        abort.set()
        raise RuntimeError(f"Could not download {key}: the file is {part.stat().st_size} bytes, not {total}.")
    if hash_thread is None:
        abort.set()
        return None
    if hash_thread.is_alive():                                            # (the last pieces are still being hashed: show that)
        import tqdm
        name = key.split("/")[-1]
        say(f"Checking {name} (it is checked as it arrives; {_size_text(total - hashed['bytes'])} is left to check)...")
        check_bar = tqdm.tqdm(desc="checking " + name, total=total, initial=hashed["bytes"], unit_scale=True, unit_divisor=1000., unit="B")
        try:
            while hash_thread.is_alive():
                check_bar.n = hashed["bytes"]
                check_bar.refresh()
                _notes[note_key] = f"checking the download ({100 * hashed['bytes'] / max(total, 1):.0f}%)"
                hash_thread.join(0.5)
        finally:
            check_bar.close()
            _notes.pop(note_key, None)
    hash_thread.join()
    abort.set()
    return hashed["digest"]


def download_file(self, file_attributes, force_download=False, skip_hash_check=False):
    """The replacement for S3CloudCache._download_file (same arguments; True if the file was downloaded, False if it was there). See the module's text."""
    import tqdm
    local = Path(file_attributes.local_path)
    part = part_path(local)
    note_key = str(local)
    with _lock_for(local):                                      # (two threads asking for the same file: the second finds it done)
        recorded = bool(self._check_successful_download(file_attributes))
        if local.exists() and not force_download:
            return False
        if local.exists():
            if recorded:                                        # (a fresh copy of a good file was asked for)
                self._remove_download_from_list(file_attributes)
                local.unlink()
                part.unlink(missing_ok=True)
                _state_path(part).unlink(missing_ok=True)
            else:                                               # (cut short by the library's own downloader: a good beginning, continue from it)
                part.unlink(missing_ok=True)
                _state_path(part).unlink(missing_ok=True)
                os.replace(local, part)
        local.parent.mkdir(parents=True, exist_ok=True)
        key = str(file_attributes.relative_path).replace("\\", "/")
        say = lambda text: tqdm.tqdm.write(text)
        started = time.monotonic()
        total = _object_size(self, key, say, started + GIVE_UP_S, note_key)
        have = downloaded_bytes(local)
        bar = tqdm.tqdm(desc=key.split("/")[-1], total=total, initial=min(have, total), unit_scale=True, unit_divisor=1000., unit="B")
        if have and total < PARALLEL_MIN:
            say(f"Continuing {key.split('/')[-1]} from {_size_text(have)} of {_size_text(total)}.")
        last_progress = [time.monotonic()]
        wrong_hash = 0
        parallel = total >= PARALLEL_MIN
        try:
            while True:
                if parallel:
                    try:
                        digest = _parallel(self, bar, part, key, total, note_key, say, last_progress, want_hash=not skip_hash_check)
                    except _NoRange:
                        say("The server does not support downloading in parts; downloading in one stream.")
                        parallel = False
                        part.unlink(missing_ok=True)
                        _state_path(part).unlink(missing_ok=True)
                        bar.n = 0
                        bar.refresh()
                        continue
                else:
                    digest = None
                    _sequential(self, bar, part, key, total, note_key, say, last_progress)
                if not skip_hash_check:
                    if (digest or _md5(part, note_key, total, key.split("/")[-1])) != file_attributes.file_hash:
                        part.unlink(missing_ok=True)
                        _state_path(part).unlink(missing_ok=True)
                        wrong_hash += 1
                        if wrong_hash >= HASH_ATTEMPTS:
                            raise RuntimeError(f"Could not download {key}: the file's check value was wrong {wrong_hash} times.")
                        say(f"The downloaded {key.split('/')[-1]} did not pass its check; downloading it again.")
                        bar.n = 0
                        bar.refresh()
                        continue
                break
        finally:
            bar.close()
            _notes.pop(note_key, None)
        os.replace(part, local)
        _state_path(part).unlink(missing_ok=True)
        self._update_list_of_downloads(file_attributes=file_attributes)
        return True
