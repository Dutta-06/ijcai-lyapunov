#!/usr/bin/env python3
"""Data access for the repeat-imaging test set of M2 (SDSS r band, other runs of Stripe 82, against the HST HUGS catalogue): fetch, verify, record checksums.

Companion of `data/astro/fetch_second.py` (whose resumable download it uses) and `data/astro/fetch.py`. It only moves and checks data; the geometry, the selection of the frames, the
patches and the validation are in `experiments/astro_repeat.py`, the findings in `results/astro/repeat_spec.md`. Standard library only; downloads untrusted third-party files into their
own new directory each and never executes them.

Sub-commands
  list                                  every file, its group, the size recorded in checksums_repeat.json and whether it is on disk
  fetch  [--group G ...] [--only KEY ...] [--runs R ...] [--frames RUN:CAMCOL:FIELD ...]
                                        download (resumable) into $ASTRO_DATA_3/<dir>/<key>/<filename>, verify sha256 against checksums_repeat.json (for a key that is not there yet: against
                                        checksums_second.json or the M2 hashes of checksums.json when the same file was fetched before) and record a new key on its first download
                                        (size, sha256, Last-Modified of the host, UTC date)
  verify [--group G ...] [--only KEY ...] ...
                                        re-hash the files on disk and compare with checksums_repeat.json

Groups
  survey      the SDSS table of imaging runs (photoRunAll-dr17.fits, 187 kB)
  m2          the inputs of the benchmark frame: r frame, psField and photoField of run 2583 camcol 2 field 136, HUGS NGC 7089 (method 1; 107 MB)
  photofield  the photoField files (one per run and camcol) of the candidate runs, RUNS below (the output of `astro_repeat.py runs`); --runs adds runs
  frames      the r frame and the psField of every candidate frame, FRAMES below (the output of `astro_repeat.py geometry`); --frames adds frames
Default group: survey.  `--group all` is every group.

Environment
  ASTRO_DATA_3   download root, default /home/user/astro_data_3 (never inside the repository)
Network: urllib through whatever HTTPS proxy / CA bundle the environment already provides; TLS verification is never disabled; a refusal (HTTP 401/403/404) is reported and not retried
or worked around. Hosts: data.sdss.org (SDSS science archive server, DR17, rerun 301), mast.stsci.edu (Download endpoint, mast:HLSP/hugs/ngc7089/...).
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fetch_second as F2      # noqa: E402  (resumable download, sha256, the URL builders; standard library only)

CHECKSUMS = os.path.join(HERE, "checksums_repeat.json")
DEFAULT_ROOT = "/home/user/astro_data_3"

# runs of photoRunAll-dr17 whose scan circle passes within 1.5 deg of the M2 windows (the output of `astro_repeat.py runs`); photoField files for all six camcols
RUNS = [1009, 1033, 1040, 1755, 2506, 2659, 2662, 2728, 3256, 3360, 3384, 3388, 3427, 3434, 3437, 4184, 4187, 4192, 4198, 4203, 4207, 4797, 4822, 4849, 4858, 4874, 4927, 4930, 4933,
        7778, 7812]
# frames (run, camcol, field) whose photoField box, enlarged by 200 px, holds a window centre (the output of `astro_repeat.py geometry`)
FRAMES = [(1040, 2, 155), (1755, 2, 131), (1755, 2, 132), (3360, 2, 116), (3388, 2, 75), (3427, 2, 102), (3427, 2, 103), (3434, 2, 108), (4187, 2, 105), (4192, 2, 111), (4203, 2, 160),
          (4822, 5, 197), (4822, 5, 198), (4930, 2, 147), (7778, 6, 72), (7812, 6, 101)]

GROUPS = ["survey", "m2", "photofield", "frames"]


def build_items(extra_runs=(), extra_frames=()):
    """key -> (group, directory under $ASTRO_DATA_3, url, filename)."""
    it = {}
    it["photorunall"] = ("survey", "survey", F2.SAS + "/photoObj/photoRunAll-dr17.fits", "photoRunAll-dr17.fits")
    for k, u, f in (F2._frame(2583, 2, 136), F2._psfield(2583, 2, 136), F2._photofield(2583, 2), F2._hugs(7089)):
        it[k] = ("m2", "ngc7089", u, f)
    for run in sorted(set(RUNS) | set(extra_runs)):
        for cam in range(1, 7):
            k, u, f = F2._photofield(run, cam)
            it[k] = ("photofield", "survey", u, f)
    for (run, cam, fld) in sorted(set(FRAMES) | set(extra_frames)):
        for k, u, f in (F2._frame(run, cam, fld), F2._psfield(run, cam, fld)):
            it[k] = ("frames", "frames", u, f)
    return it


# ----------------------------------------------------------------------------- helpers
def root():
    return os.environ.get("ASTRO_DATA_3", DEFAULT_ROOT)


def load_sums(path=None):
    path = path or CHECKSUMS
    if not os.path.exists(path):
        return {"schema": 1, "note": ("sha256 of every file fetched by data/astro/fetch_repeat.py, measured from the hosts named in results/astro/repeat_spec.md. fetch_repeat.py verifies "
                                      "against this file and records a new entry on the first download of an unlisted key; a key that data/astro/checksums_second.json also records "
                                      "must have the same sha256 and size."), "files": {}}
    with open(path) as f:
        return json.load(f)


def save_sums(d, path=None):
    path = path or CHECKSUMS
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(d, f, indent=1, sort_keys=False)
        f.write("\n")
    os.replace(tmp, path)


def earlier_record(key, second_path=F2.CHECKSUMS):
    """(sha256, bytes) recorded for the same key by the second-cluster work (checksums_second.json) or, for the M2 files, by fetch.py (F2.M2_EXPECTED); (None, None) if none."""
    if os.path.exists(second_path):
        with open(second_path) as f:
            rec = json.load(f).get("files", {}).get(key)
        if rec and rec.get("sha256"):
            return rec["sha256"], rec.get("bytes")
    return F2.M2_EXPECTED.get(key, (None, None))


def dest_path(items, key):
    return os.path.join(root(), items[key][1], key, items[key][3])


def select(args, items):
    if args.only:
        bad = [k for k in args.only if k not in items]
        if bad:
            sys.exit("unknown key(s): %s (see `list`)" % ", ".join(bad))
        return list(args.only)
    groups = args.group or ["survey"]
    if "all" in groups:
        groups = GROUPS
    return [k for k in items if items[k][0] in groups]


def check_or_record(items, key, path, sums, lastmod=None, second_path=F2.CHECKSUMS):
    """Verify the file against checksums_repeat.json, or (first download) against the earlier record of the same key, then record it. A mismatch moves the file to *.bad and stops."""
    digest, size = F2.sha256_file(path)
    rec = sums["files"].get(key)
    earlier_digest, earlier_bytes = earlier_record(key, second_path)
    exp_digest, exp_bytes = (rec["sha256"], rec.get("bytes")) if rec and rec.get("sha256") else (earlier_digest, earlier_bytes)
    if exp_digest and (digest != exp_digest or (exp_bytes and size != exp_bytes)):
        bad = path + ".bad"
        os.replace(path, bad)
        raise SystemExit("CHECKSUM MISMATCH for %s: got %s (%d bytes), expected %s (%s bytes); moved to %s" % (key, digest, size, exp_digest, exp_bytes, bad))
    if rec and rec.get("sha256"):
        print("  ok   %s  sha256 %s...  %d bytes" % (key, digest[:16], size))
    else:
        sums["files"][key] = {"filename": items[key][3], "bytes": size, "sha256": digest, "source_last_modified": lastmod or F2.head_lastmod(items[key][2]), "first_recorded_utc": F2.utcnow()}
        if earlier_digest:
            sums["files"][key]["matches_earlier_record"] = True
        save_sums(sums)
        print("  recorded (first download) %s  sha256 %s  %d bytes%s" % (key, digest, size, "  (equal to the earlier record)" if earlier_digest else ""))
    return digest


# ----------------------------------------------------------------------------- commands
def cmd_list(args, items):
    sums = load_sums()
    print("ASTRO_DATA_3 = %s" % root())
    for k in (select(args, items) if (args.only or args.group) else list(items)):
        g, d, url, fn = items[k]
        rec = sums["files"].get(k, {})
        here = "present" if os.path.exists(dest_path(items, k)) else "absent "
        print("%-34s %-10s %s %-14s %s" % (k, g, here, ("%s B" % rec.get("bytes", "?")), url.split("?")[0] if "Download" not in url else "MAST Download API (%s)" % fn))


def cmd_fetch(args, items):
    sums = load_sums()
    for key in select(args, items):
        d = dest_path(items, key)
        print("[%s]" % key)
        lastmod = None
        if os.path.exists(d):
            print("  already on disk; verifying")
        else:
            lastmod = F2.download(items[key][2], d)
        check_or_record(items, key, d, sums, lastmod)


def cmd_verify(args, items):
    sums = load_sums()
    rc = 0
    for key in select(args, items):
        d = dest_path(items, key)
        if not os.path.exists(d):
            print("[%s] missing: %s" % (key, d))
            rc = 1
            continue
        digest, size = F2.sha256_file(d)
        rec = sums["files"].get(key, {})
        ok = digest == rec.get("sha256") and size == rec.get("bytes")
        print("[%s] %s  %s  %d bytes" % (key, "OK " if ok else "BAD", digest[:16], size))
        rc |= (not ok)
    sys.exit(rc)


def parse_frames(specs):
    out = []
    for s in specs or []:
        run, cam, fld = (int(v) for v in s.split(":"))
        out.append((run, cam, fld))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("list", "fetch", "verify"):
        p = sub.add_parser(name)
        p.add_argument("--group", action="append", choices=GROUPS + ["all"])
        p.add_argument("--only", nargs="+")
        p.add_argument("--runs", nargs="+", type=int, default=[], help="extra runs for the photofield group")
        p.add_argument("--frames", nargs="+", default=[], help="extra frames RUN:CAMCOL:FIELD for the frames group")
    args = ap.parse_args()
    items = build_items(args.runs, parse_frames(args.frames))
    {"list": cmd_list, "fetch": cmd_fetch, "verify": cmd_verify}[args.cmd](args, items)


if __name__ == "__main__":
    main()
