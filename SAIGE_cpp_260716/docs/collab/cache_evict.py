#!/usr/bin/env python3
"""cache_evict.py MODE OUT_JSON PATH [PATH ...]  -- start a timed run with a cold page cache.

PATH: files the run reads (genotype files, null-model directories -- read recursively,
symlinks followed --, variance-ratio files, sparse GRM files).
MODE:
  auto         drop_caches if passwordless sudo works, else vmtouch if installed, else fadvise
  drop_caches  sync; echo 3 | sudo -n tee /proc/sys/vm/drop_caches   (whole machine)
  vmtouch      vmtouch -e PATH...                                     (only these files)
  fadvise      per file: fsync, then posix_fadvise(POSIX_FADV_DONTNEED) (only these files, no root)
  none         nothing (the run is then NOT cold; recorded as such)
Verification: page-cache residency of every file before and after, measured with mincore(2)
on a private mapping (Python ctypes; the same numbers fincore(1) / vmtouch(8) report).
OUT_JSON records the method actually used, the bytes still resident after eviction, and the
filesystem type of each file (on Lustre / GPFS / NFS, client-side eviction may not reach the
server's cache: report the filesystem type with the results).
"""
import ctypes, ctypes.util, json, mmap, os, shutil, subprocess, sys

PAGE = os.sysconf("SC_PAGE_SIZE")
libc = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)
libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_ubyte)]


def resident_bytes(path):
    size = os.path.getsize(path)
    if size == 0:
        return 0
    with open(path, "rb") as f:
        mm = mmap.mmap(f.fileno(), size, access=mmap.ACCESS_COPY)
        try:
            buf = (ctypes.c_char * size).from_buffer(mm)
            npg = (size + PAGE - 1) // PAGE
            vec = (ctypes.c_ubyte * npg)()
            if libc.mincore(ctypes.addressof(buf), size, vec) != 0:
                return -1
            n = sum(v & 1 for v in vec)
            del buf
        finally:
            mm.close()
    return min(n * PAGE, size)


def files_of(paths):
    out = []
    for p in paths:
        if os.path.isdir(p):
            for root, _, fs in os.walk(p, followlinks=True):
                out += [os.path.join(root, f) for f in sorted(fs)]
        elif os.path.exists(p):
            out.append(p)
    return [os.path.realpath(f) for f in out if os.path.isfile(f)]


def fstype(path):
    try:
        return subprocess.run(["stat", "-f", "-c", "%T", path], capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


def sudo_ok():
    return subprocess.run(["sudo", "-n", "true"], capture_output=True).returncode == 0


def main():
    mode, outj, paths = sys.argv[1], sys.argv[2], sys.argv[3:]
    files = files_of(paths)
    before = {f: resident_bytes(f) for f in files}
    if mode == "auto":
        mode = "drop_caches" if sudo_ok() else ("vmtouch" if shutil.which("vmtouch") else "fadvise")
    err = ""
    if mode == "drop_caches":
        r = subprocess.run("sync; echo 3 | sudo -n tee /proc/sys/vm/drop_caches > /dev/null", shell=True)
        if r.returncode != 0:
            err = "drop_caches failed (no passwordless sudo?); fell back to fadvise"
            mode = "fadvise"
    if mode == "vmtouch":
        subprocess.run(["vmtouch", "-q", "-e"] + files, check=False)
    if mode == "fadvise":
        for f in files:
            fd = os.open(f, os.O_RDONLY)
            try:
                try:
                    os.fsync(fd)
                except OSError:
                    pass
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            finally:
                os.close(fd)
    after = {f: resident_bytes(f) for f in files}
    tot = sum(os.path.getsize(f) for f in files)
    rec = {
        "method": mode, "note": err, "n_files": len(files), "bytes_total": tot,
        "resident_before": sum(before.values()), "resident_after": sum(after.values()),
        "cold": mode != "none" and sum(after.values()) <= 0.01 * max(tot, 1),
        "verify": "mincore",
        "fstypes": sorted({fstype(f) for f in files}),
        "largest_files": [{"bytes": os.path.getsize(f), "resident_before": before[f], "resident_after": after[f],
                           "fstype": fstype(f), "ext": os.path.splitext(f)[1]}
                          for f in sorted(files, key=os.path.getsize, reverse=True)[:3]],
    }
    json.dump(rec, open(outj, "w"), indent=1)
    print("cache: %s, %d files, %.1f MB, resident before %.1f MB, after %.1f MB%s" % (
        mode, len(files), tot / 1e6, rec["resident_before"] / 1e6, rec["resident_after"] / 1e6,
        "" if rec["cold"] else "  (NOT cold)"))


if __name__ == "__main__":
    main()
