"""Downloads of the atlas files that continue where they stopped.

abc_atlas_access writes a download straight to the file's final name, in one request, with no way to continue: a download cut short (the
program closed or crashed, the network went away) leaves a file that is there but not whole, and the next time it is asked for the library
deletes it and starts from the first byte. For the 50 GB imputed dataset that can be an hour or more lost.

install() replaces the library's S3 download method (for every cache, whoever makes it) with download_file() here:
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
import os
import threading
import time
from pathlib import Path

CHUNK = 1 << 20
RETRY_FIRST_S = 5                # the wait after the first failure, doubled each time up to RETRY_MAX_S
RETRY_MAX_S = 60
GIVE_UP_S = 2 * 3600             # no byte for this long: stop (the next ask of the file continues from the same place)
HASH_ATTEMPTS = 2                # whole downloads whose hash was wrong, before giving up
_locks = {}
_notes = {}
_guard = threading.Lock()


def _size_text(n):
    return f"{n / 1e9:.2f} GB" if n >= 1e9 else f"{n / 1e6:.0f} MB" if n >= 1e6 else f"{n / 1e3:.0f} KB"


def part_path(path):
    path = Path(path)
    return path.with_name(path.name + ".part")


def downloaded_bytes(path):
    """How much of the file is on disk: the .part file's size while it is being downloaded, else the file's own (0 if neither is there)."""
    sizes = [0]
    for candidate in (part_path(path), Path(path)):
        try:
            sizes.append(candidate.stat().st_size)
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


def _md5(path, key_text, total):
    hasher = hashlib.md5()
    done, last = 0, 0.0
    with open(path, "rb") as f:
        while True:
            chunk = f.read(8 * CHUNK)
            if not chunk:
                break
            hasher.update(chunk)
            done += len(chunk)
            if time.monotonic() - last > 1.0:
                last = time.monotonic()
                _notes[key_text] = f"checking the download ({100 * done / max(total, 1):.0f}%)"
    _notes.pop(key_text, None)
    return hasher.hexdigest()


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
            else:                                               # (cut short by the library's own downloader: a good beginning, continue from it)
                part.unlink(missing_ok=True)
                os.replace(local, part)
        local.parent.mkdir(parents=True, exist_ok=True)
        key = str(file_attributes.relative_path).replace("\\", "/")
        say = lambda text: tqdm.tqdm.write(text)
        started = time.monotonic()
        total = _object_size(self, key, say, started + GIVE_UP_S, note_key)
        have = part.stat().st_size if part.exists() else 0
        bar = tqdm.tqdm(desc=key.split("/")[-1], total=total, initial=min(have, total), unit_scale=True, unit_divisor=1000., unit="MB")
        if have:
            say(f"Continuing {key.split('/')[-1]} from {_size_text(have)} of {_size_text(total)}.")
        wrong_hash = 0
        try:
            while True:
                have = part.stat().st_size if part.exists() else 0
                if have > total:
                    part.unlink()
                    have = 0
                    bar.n = 0
                    bar.refresh()
                delay, last_progress = RETRY_FIRST_S, time.monotonic()
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
                                last_progress = time.monotonic()
                                delay = RETRY_FIRST_S
                                _notes.pop(note_key, None)
                    except Exception as e:
                        if _fatal(e):
                            raise
                        have = part.stat().st_size if part.exists() else 0         # (what is on the disk is where it continues from)
                        bar.n = have
                        bar.refresh()
                        if time.monotonic() - last_progress > GIVE_UP_S:
                            raise RuntimeError(f"Could not download {key}: no data for {GIVE_UP_S // 3600} hours ({type(e).__name__}: {e}). "
                                               "What was downloaded is kept; asking for the file again continues from it.")
                        say(f"The download of {key.split('/')[-1]} was interrupted ({type(e).__name__}); continuing from {_size_text(have)} in {delay} s.")
                        _pause(note_key, delay, type(e).__name__, say)
                        delay = min(delay * 2, RETRY_MAX_S)
                        continue
                    have = part.stat().st_size if part.exists() else 0
                    if have < total:                                                # (the stream ended early without an error: ask for the rest)
                        if time.monotonic() - last_progress > GIVE_UP_S:
                            raise RuntimeError(f"Could not download {key}: the server sent no more data for {GIVE_UP_S // 3600} hours.")
                        time.sleep(1)
                if not skip_hash_check:
                    if _md5(part, note_key, total) != file_attributes.file_hash:
                        part.unlink(missing_ok=True)
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
        self._update_list_of_downloads(file_attributes=file_attributes)
        return True
