#!/usr/bin/env python3
"""Data access for the second-cluster test set (SDSS r band against the HST HUGS catalogue): fetch, verify, record checksums.

Companion of `data/astro/fetch.py` (the M2 / NGC 7089 patches). It only moves and checks data; the alignment, the units and the
patches are in `experiments/astro_second.py`, the findings in `results/astro/second_cluster_spec.md`.
Standard library only; downloads untrusted third-party files into their own new directory each and never executes them.

Sub-commands
  list                                  every file, its group, the size recorded in checksums_second.json and whether it is on disk
  fetch  [--group G ...] [--only KEY ...] [--photofield-runs R ...]
                                        download (resumable) into $ASTRO_DATA_2/<dir>/<key>/<filename>, verify sha256 against
                                        checksums_second.json (or against the M2 hashes of checksums.json for the M2 files), record
                                        a new key on its first download (size, sha256, Last-Modified of the host, UTC date)
  verify [--group G ...] [--only KEY ...]
                                        re-hash the files on disk and compare with checksums_second.json
  resolve                               name-resolve the cluster centres with the MAST name lookup and record them in checksums_second.json
  probe                                 ask the MAST Download endpoint (1-byte range request) which HUGS method-1 catalogues exist for the clusters of PROBED and
                                        PROBE_ABSENT, and record status and size under 'hugs_probe' in checksums_second.json

Groups
  survey      the SDSS table of imaging runs (photoRunAll-dr17.fits, 187 kB) used to find the runs whose scan strip crosses a cluster
  m2          the inputs of the validation on M2: r frame, psField and photoField of run 2583 camcol 2 field 136, HUGS NGC 7089 (method 1)
  hugs        the HUGS catalogues (method 1) of the northern clusters of the coverage survey (2.7 to 106 MB each)
  photofield  photoField files (one per run and camcol; 0.2 to 1.3 MB) of the candidate runs of the survey (--photofield-runs adds runs)
  target      r frame and psField of the frames examined as hosts of the test patches: M53 (4 patches), M92 (2), and the two that gave fewer than two (M3, M5)
Default group: target.  `--group all` is every group.

Environment
  ASTRO_DATA_2   download root, default /home/user/astro_data_2 (never inside the repository)
Network: urllib through whatever HTTPS proxy / CA bundle the environment already provides; TLS verification is never disabled; a refusal
(HTTP 401/403/404) is reported and not retried or worked around. Hosts: data.sdss.org (SDSS science archive server, DR17, rerun 301),
mast.stsci.edu (Download endpoint, mast:HLSP/hugs/<cluster>/...).
"""
import argparse
import datetime
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
CHECKSUMS = os.path.join(HERE, "checksums_second.json")
DEFAULT_ROOT = "/home/user/astro_data_2"

SAS = "https://data.sdss.org/sas/dr17/eboss"
MAST_DL = "https://mast.stsci.edu/api/v0.1/Download/file"
MAST_INVOKE = "https://mast.stsci.edu/api/v0/invoke"
HUGS = "mast:HLSP/hugs/%s/hlsp_hugs_hst_wfc3-uvis-acs-wfc_%s_multi_v1_catalog-meth1.txt"

# the sha256 of the M2 files as recorded by data/astro/fetch.py (checksums.json): the files fetched here for the validation must be the same bytes
M2_EXPECTED = {
    "sdss_frame_r_002583_2_0136": ("9129b6c6779ad8155741998e333e0f69395226aa9fc425dc9933f9f6f53ec04c", 3524824),
    "sdss_psfield_002583_2_0136": ("3d24d270c12dcd8429357a81433d0151d3a6e8739401b18273d4bfc206bcd7bd", 414720),
    "sdss_photofield_002583_2": ("052132367a8faaee18160b71c3a0beef423259df7b9218c313d564f1d5f6e69e", 650880),
    "hugs_ngc7089": ("4a96ae437122468785200570dd9121f4b677ec6ed0200e65a687babcf5b58f49", 106805615),
}

# the northern (and equatorial) clusters of the HUGS survey that have candidate SDSS runs in the survey (HUGS catalogues present on MAST: all 52 tried)
CLUSTERS = [5024, 5053, 5272, 5466, 5904, 6205, 6254, 6341, 6535, 6838, 6934, 6981, 7078, 7089]
# all clusters probed on MAST for the coverage table (a HUGS method-1 catalogue answers HTTP 206 to a range request for these)
PROBED = [1261, 1851, 2298, 2808, 3201, 4590, 4833, 5024, 5053, 5272, 5286, 5466, 5904, 5927, 5986, 6093, 6101, 6121, 6144, 6171, 6205, 6218, 6254, 6304,
          6341, 6352, 6362, 6366, 6388, 6397, 6441, 6496, 6535, 6541, 6584, 6624, 6637, 6652, 6656, 6681, 6715, 6717, 6723, 6752, 6779, 6809, 6838, 6934,
          6981, 7078, 7089, 7099]
# clusters tried with the same file name pattern whose catalogue the endpoint did not find (HTTP 404; 47 Tuc, omega Cen and a few others are filed under other names or are not in HUGS)
PROBE_ABSENT = [104, 288, 362, 1904, 2419, 4147, 5139, 5634, 6229, 6402, 6426, 6712, 7492]
# runs of photoRunAll-dr17 whose scan strip crosses one of the clusters (the output of `astro_second.py survey`); photoField files for all six camcols
SURVEY_RUNS = [1458, 1463, 1478, 1738, 1739, 1749, 2190, 2305, 2327, 2328, 2507, 2566, 3170, 3206, 3225, 3226, 4588, 4623, 4646, 4649, 4682, 4822,
               5301, 5308, 5314, 5327, 5360, 5375, 5376, 5382, 5390, 5391, 5392, 5398, 5403, 5414, 5416, 6282, 6895]
# frames of the test patches: cluster -> (run, camcol, field)
TARGETS = {"ngc6341": (4682, 6, 208), "ngc5024": (5390, 6, 162), "ngc5904": (1458, 4, 700), "ngc5272": (4646, 3, 81)}


def _frame(run, camcol, field):
    return ("sdss_frame_r_%06d_%d_%04d" % (run, camcol, field), SAS + "/photoObj/frames/301/%d/%d/frame-r-%06d-%d-%04d.fits.bz2" % (run, camcol, run, camcol, field),
            "frame-r-%06d-%d-%04d.fits.bz2" % (run, camcol, field))


def _psfield(run, camcol, field):
    return ("sdss_psfield_%06d_%d_%04d" % (run, camcol, field), SAS + "/photo/redux/301/%d/objcs/%d/psField-%06d-%d-%04d.fit" % (run, camcol, run, camcol, field),
            "psField-%06d-%d-%04d.fit" % (run, camcol, field))


def _photofield(run, camcol):
    return ("sdss_photofield_%06d_%d" % (run, camcol), SAS + "/photoObj/301/%d/photoField-%06d-%d.fits" % (run, run, camcol), "photoField-%06d-%d.fits" % (run, camcol))


def _hugs(n):
    c = "ngc%d" % n
    return ("hugs_" + c, MAST_DL + "?" + urllib.parse.urlencode({"uri": HUGS % (c, c)}), "hlsp_hugs_hst_wfc3-uvis-acs-wfc_%s_multi_v1_catalog-meth1.txt" % c)


def build_items(extra_runs=()):
    """key -> (group, directory under $ASTRO_DATA_2, url, filename)."""
    it = {}
    it["photorunall"] = ("survey", "survey", SAS + "/photoObj/photoRunAll-dr17.fits", "photoRunAll-dr17.fits")
    for k, u, f in (_frame(2583, 2, 136), _psfield(2583, 2, 136), _photofield(2583, 2), _hugs(7089)):
        it[k] = ("m2", "ngc7089", u, f)
    for n in CLUSTERS:
        if n != 7089:
            k, u, f = _hugs(n)
            it[k] = ("hugs", "ngc%d" % n, u, f)
    for run in sorted(set(SURVEY_RUNS) | set(extra_runs)):
        for cam in range(1, 7):
            k, u, f = _photofield(run, cam)
            it[k] = ("photofield", "survey", u, f)
    for c, (run, cam, fld) in TARGETS.items():
        for k, u, f in (_frame(run, cam, fld), _psfield(run, cam, fld)):
            it[k] = ("target", c, u, f)
    return it


GROUPS = ["survey", "m2", "hugs", "photofield", "target"]


# ----------------------------------------------------------------------------- helpers
def root():
    return os.environ.get("ASTRO_DATA_2", DEFAULT_ROOT)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_sums():
    if not os.path.exists(CHECKSUMS):
        return {"schema": 1, "note": ("sha256 of every file fetched by data/astro/fetch_second.py, measured from the hosts named in results/astro/second_cluster_spec.md. "
                                      "fetch_second.py verifies against this file and records a new entry on the first download of an unlisted key."), "files": {}}
    with open(CHECKSUMS) as f:
        return json.load(f)


def save_sums(d):
    tmp = CHECKSUMS + ".tmp"
    with open(tmp, "w") as f:
        json.dump(d, f, indent=1, sort_keys=False)
        f.write("\n")
    os.replace(tmp, CHECKSUMS)


def sha256_file(path, opener=open):
    h = hashlib.sha256()
    n = 0
    with opener(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
            n += len(blk)
    return h.hexdigest(), n


def select(args, items):
    if args.only:
        bad = [k for k in args.only if k not in items]
        if bad:
            sys.exit("unknown key(s): %s (see `list`)" % ", ".join(bad))
        return list(args.only)
    groups = args.group or ["target"]
    if "all" in groups:
        groups = GROUPS
    return [k for k in items if items[k][0] in groups]


def dest_path(items, key):
    return os.path.join(root(), items[key][1], key, items[key][3])


# ----------------------------------------------------------------------------- resumable download
def download(url, dest, attempts=5, timeout=60):
    """Download url to dest, resuming dest+'.part' when the server honours Range. A stream that ends before Content-Length is kept as a .part file and
    resumed on the next attempt. Returns the Last-Modified header of the response (or None)."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + ".part"
    for attempt in range(1, attempts + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        req = urllib.request.Request(url, headers={"User-Agent": "astro-fetch/1.0"})
        if have:
            req.add_header("Range", "bytes=%d-" % have)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                status = resp.status
                lastmod = resp.headers.get("Last-Modified")
                if have and status == 200:          # server ignored Range: start over
                    have = 0
                mode = "ab" if have and status == 206 else "wb"
                clen = resp.headers.get("Content-Length")
                want = have + int(clen) if clen is not None else None
                print("  GET %s (%s, resume at %d bytes, content-length %s)" % (url.split("?")[0], status, have if mode == "ab" else 0, clen))
                with open(part, mode) as out:
                    for blk in iter(lambda: resp.read(1 << 20), b""):
                        out.write(blk)
            if want is not None and os.path.getsize(part) < want:
                print("  attempt %d/%d: stream ended after %d of %d bytes; will resume" % (attempt, attempts, os.path.getsize(part), want))
                time.sleep(min(2 ** attempt, 20))
                continue
            os.replace(part, dest)
            return lastmod
        except urllib.error.HTTPError as e:
            if e.code == 416 and have:              # range not satisfiable: file may already be complete
                os.replace(part, dest)
                return None
            print("  attempt %d/%d: HTTP %s %s" % (attempt, attempts, e.code, e.reason))
            if e.code in (401, 403, 404):           # policy / missing: do not hammer, do not work around
                raise SystemExit("giving up on %s (HTTP %s); not retrying a refusal" % (url.split("?")[0], e.code))
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            print("  attempt %d/%d: %s" % (attempt, attempts, e))
        time.sleep(min(2 ** attempt, 20))
    raise SystemExit("download failed after %d attempts: %s" % (attempts, url.split("?")[0]))


def head_lastmod(url):
    """Last-Modified of a host file through a 1-byte range request (for files already on disk that are recorded for the first time)."""
    req = urllib.request.Request(url, headers={"User-Agent": "astro-fetch/1.0", "Range": "bytes=0-0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.headers.get("Last-Modified")
    except (urllib.error.URLError, OSError, TimeoutError):
        return None


def check_or_record(items, key, path, sums, lastmod=None):
    digest, size = sha256_file(path)
    rec = sums["files"].get(key)
    exp_digest, exp_bytes = (rec["sha256"], rec.get("bytes")) if rec and rec.get("sha256") else M2_EXPECTED.get(key, (None, None))
    if exp_digest:
        if digest != exp_digest or (exp_bytes and size != exp_bytes):
            bad = path + ".bad"
            os.replace(path, bad)
            raise SystemExit("CHECKSUM MISMATCH for %s: got %s (%d bytes), expected %s (%s bytes); moved to %s" % (key, digest, size, exp_digest, exp_bytes, bad))
    if rec and rec.get("sha256"):
        print("  ok   %s  sha256 %s...  %d bytes" % (key, digest[:16], size))
    else:
        sums["files"][key] = {"filename": items[key][3], "bytes": size, "sha256": digest, "source_last_modified": lastmod or head_lastmod(items[key][2]),
                              "first_recorded_utc": utcnow()}
        if exp_digest:
            sums["files"][key]["matches_data_astro_checksums_json"] = True
        save_sums(sums)
        print("  recorded (first download) %s  sha256 %s  %d bytes%s" % (key, digest, size, "  (equal to the hash recorded by fetch.py)" if exp_digest else ""))
    return digest


# ----------------------------------------------------------------------------- commands
def cmd_list(args, items):
    sums = load_sums()
    print("ASTRO_DATA_2 = %s" % root())
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
            lastmod = download(items[key][2], d)
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
        digest, size = sha256_file(d)
        rec = sums["files"].get(key, {})
        ok = digest == rec.get("sha256") and size == rec.get("bytes")
        print("[%s] %s  %s  %d bytes" % (key, "OK " if ok else "BAD", digest[:16], size))
        rc |= (not ok)
    sys.exit(rc)


def cmd_resolve(args, items):
    """MAST name lookup (SIMBAD) of the cluster centres; recorded under 'cluster_centres' with the date."""
    sums = load_sums()
    cen = sums.setdefault("cluster_centres", {})
    for n in PROBED:
        c = "ngc%d" % n
        req = json.dumps({"service": "Mast.Name.Lookup", "params": {"input": "ngc %d" % n, "format": "json"}, "format": "json"})
        url = MAST_INVOKE + "?" + urllib.parse.urlencode({"request": req})
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "astro-fetch/1.0"}), timeout=60) as resp:
            r = json.load(resp)["resolvedCoordinate"][0]
        cen[c] = {"name": " ".join(r["canonicalName"].split()), "ra": r["ra"], "dec": r["decl"], "resolver": r["resolver"], "recorded_utc": utcnow()}
        print("%-8s %-10s RA %.5f Dec %+.5f (%s)" % (c, cen[c]["name"], r["ra"], r["decl"], r["resolver"]))
        time.sleep(0.2)
    save_sums(sums)


def cmd_probe(args, items):
    """One 1-byte range request per cluster to the Download endpoint: 206 and the total size (Content-Range) if the catalogue exists, 404 if it does not."""
    sums = load_sums()
    rec = sums.setdefault("hugs_probe", {})
    for n in sorted(PROBED + PROBE_ABSENT):
        k, url, fn = _hugs(n)
        req = urllib.request.Request(url, headers={"User-Agent": "astro-fetch/1.0", "Range": "bytes=0-0"})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                total = (resp.headers.get("Content-Range") or "").split("/")[-1]
                rec["ngc%d" % n] = {"status": resp.status, "bytes": int(total) if total.isdigit() else None, "last_modified": resp.headers.get("Last-Modified")}
        except urllib.error.HTTPError as e:
            rec["ngc%d" % n] = {"status": e.code, "bytes": None, "last_modified": None}
        print("ngc%-5d HTTP %s  %s bytes" % (n, rec["ngc%d" % n]["status"], rec["ngc%d" % n]["bytes"]))
        time.sleep(0.2)
    sums["hugs_probe_recorded_utc"] = utcnow()
    save_sums(sums)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("list", "fetch", "verify"):
        p = sub.add_parser(name)
        p.add_argument("--group", action="append", choices=GROUPS + ["all"])
        p.add_argument("--only", nargs="+")
        p.add_argument("--photofield-runs", nargs="+", type=int, default=[], help="extra runs for the photofield group")
    sub.add_parser("resolve")
    sub.add_parser("probe")
    args = ap.parse_args()
    items = build_items(getattr(args, "photofield_runs", []))
    {"list": cmd_list, "fetch": cmd_fetch, "verify": cmd_verify, "resolve": cmd_resolve, "probe": cmd_probe}[args.cmd](args, items)


if __name__ == "__main__":
    main()
