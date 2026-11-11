#!/usr/bin/env python3
"""Data access for the M2 (NGC 7089) crowded-field patch benchmark: fetch, verify, extract.

This script only moves and checks data. It runs no cataloguing method.

Sub-commands
  list                      show the catalogue of files, their groups and the recorded checksums
  fetch  [--group G ...] [--only KEY ...]
                            download (resumable) into $ASTRO_DATA/<key>/<filename>, verify sha256
                            against checksums.json (or record it there on first download)
  verify [--group G ...] [--only KEY ...]
                            re-hash files already on disk and compare with checksums.json
  extract-patch [--x0 310 --y0 630 --size 100]
                            rebuild the DN / electron image of the 100x100 r-band patch and the
                            SDSS PSF stamp from the downloaded frame + psField (needs numpy,
                            scipy, astropy). Run it as `python -I data/astro/fetch.py extract-patch`
                            because it reads downloaded (untrusted) files.

Groups: core (r-band frame + psField), geometry (photoField), hst (HST catalogue, ~107 MB each),
        multiband (g and i frames), all (all of these). The separate `shipped` sub-command fetches the four small
        input files shipped in a public git repository (needs ASTRO_LION_GIT_URL).

Environment
  ASTRO_DATA            download root, default /home/user/astro_data (never inside the repository)
  ASTRO_LION_GIT_URL    https git URL of the public repository named `lion` (see
                        results/astro/data_access.md); only for the `shipped` sub-command.
Network: urllib through whatever HTTPS proxy / CA bundle the environment already provides.
TLS verification is never disabled by this script.
"""
import argparse
import bz2
import datetime
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
CHECKSUMS = os.path.join(HERE, "checksums.json")
DEFAULT_ROOT = "/home/user/astro_data"

SAS = "https://data.sdss.org/sas/dr17/eboss"
MAST_DL = "https://mast.stsci.edu/api/v0.1/Download/file"
HUGS = "mast:HLSP/hugs/ngc7089/hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc7089_multi_v1_catalog-meth%d.txt"

# key -> (group, url, filename). Sizes and sha256 live in checksums.json (single source of truth).
ITEMS = {
    "sdss_frame_r": ("core", SAS + "/photoObj/frames/301/2583/2/frame-r-002583-2-0136.fits.bz2", "frame-r-002583-2-0136.fits.bz2"),
    "sdss_psfield": ("core", SAS + "/photo/redux/301/2583/objcs/2/psField-002583-2-0136.fit", "psField-002583-2-0136.fit"),
    "sdss_photofield": ("geometry", SAS + "/photoObj/301/2583/photoField-002583-2.fits", "photoField-002583-2.fits"),
    "sdss_frame_g": ("multiband", SAS + "/photoObj/frames/301/2583/2/frame-g-002583-2-0136.fits.bz2", "frame-g-002583-2-0136.fits.bz2"),
    "sdss_frame_i": ("multiband", SAS + "/photoObj/frames/301/2583/2/frame-i-002583-2-0136.fits.bz2", "frame-i-002583-2-0136.fits.bz2"),
    "hugs_ngc7089_meth1": ("hst", MAST_DL + "?" + urllib.parse.urlencode({"uri": HUGS % 1}), "hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc7089_multi_v1_catalog-meth1.txt"),
    "hugs_ngc7089_meth2": ("hst", MAST_DL + "?" + urllib.parse.urlencode({"uri": HUGS % 2}), "hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc7089_multi_v1_catalog-meth2.txt"),
}
GROUPS = ["core", "geometry", "hst", "multiband"]


# ----------------------------------------------------------------------------- helpers
def root():
    return os.environ.get("ASTRO_DATA", DEFAULT_ROOT)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_sums():
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


def select(args):
    keys = list(ITEMS)
    if args.only:
        bad = [k for k in args.only if k not in ITEMS]
        if bad:
            sys.exit("unknown key(s): %s (see `list`)" % ", ".join(bad))
        return args.only
    groups = args.group or ["core"]
    if "all" in groups:
        groups = GROUPS
    return [k for k in keys if ITEMS[k][0] in groups]


def dest_path(key):
    return os.path.join(root(), key, ITEMS[key][2])


# ----------------------------------------------------------------------------- resumable download
def download(url, dest, attempts=5, timeout=60):
    """Download url to dest, resuming dest+'.part' when the server honours Range.
    A stream that ends before Content-Length is kept as a .part file and resumed on the next attempt."""
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
            return
        except urllib.error.HTTPError as e:
            if e.code == 416 and have:              # range not satisfiable: file may already be complete
                os.replace(part, dest)
                return
            print("  attempt %d/%d: HTTP %s %s" % (attempt, attempts, e.code, e.reason))
            if e.code in (401, 403, 404):           # policy / missing: do not hammer
                raise SystemExit("giving up on %s (HTTP %s); not retrying a refusal" % (url.split("?")[0], e.code))
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            print("  attempt %d/%d: %s" % (attempt, attempts, e))
        time.sleep(min(2 ** attempt, 20))
    raise SystemExit("download failed after %d attempts: %s" % (attempts, url.split("?")[0]))


def check_or_record(key, path, sums):
    digest, size = sha256_file(path)
    rec = sums["files"].get(key)
    if rec and rec.get("sha256"):
        if digest != rec["sha256"] or (rec.get("bytes") and size != rec["bytes"]):
            bad = path + ".bad"
            os.replace(path, bad)
            raise SystemExit("CHECKSUM MISMATCH for %s: got %s (%d bytes), expected %s (%s bytes); moved to %s"
                             % (key, digest, size, rec["sha256"], rec.get("bytes"), bad))
        print("  ok   %s  sha256 %s...  %d bytes" % (key, digest[:16], size))
    else:
        sums["files"][key] = {"filename": ITEMS[key][2], "bytes": size, "sha256": digest, "first_recorded_utc": utcnow()}
        save_sums(sums)
        print("  recorded (first download) %s  sha256 %s  %d bytes" % (key, digest, size))
    return digest


# ----------------------------------------------------------------------------- commands
def cmd_list(args):
    sums = load_sums()
    print("ASTRO_DATA = %s" % root())
    for k, (g, url, fn) in ITEMS.items():
        rec = sums["files"].get(k, {})
        here = "present" if os.path.exists(dest_path(k)) else "absent "
        print("%-20s %-9s %s %-14s %s" % (k, g, here, ("%s B" % rec.get("bytes", "?")), url.split("?")[0] if "Download" not in url else "MAST Download API"))
    print("shipped sub-command: git, pinned commit %s (needs ASTRO_LION_GIT_URL)" % sums["shipped_inputs"]["commit"])


def cmd_fetch(args):
    sums = load_sums()
    for key in select(args):
        d = dest_path(key)
        print("[%s]" % key)
        if os.path.exists(d):
            print("  already on disk; verifying")
        else:
            download(ITEMS[key][1], d)
        check_or_record(key, d, sums)
        if key == "sdss_frame_r":                   # also pin the decompressed frame
            dig, n = sha256_file(d, opener=bz2.open)
            rec = sums["derived"].get("frame_r_decompressed")
            if rec and rec["sha256"] != dig:
                raise SystemExit("decompressed frame checksum mismatch: %s vs %s" % (dig, rec["sha256"]))
            print("  decompressed frame sha256 %s...  %d bytes %s" % (dig[:16], n, "(matches record)" if rec else ""))


def cmd_verify(args):
    sums = load_sums()
    rc = 0
    for key in select(args):
        d = dest_path(key)
        if not os.path.exists(d):
            print("[%s] missing: %s" % (key, d)); rc = 1; continue
        digest, size = sha256_file(d)
        rec = sums["files"].get(key, {})
        ok = digest == rec.get("sha256") and size == rec.get("bytes")
        print("[%s] %s  %s  %d bytes" % (key, "OK " if ok else "BAD", digest[:16], size)); rc |= (not ok)
    sys.exit(rc)


def cmd_shipped(args):
    """Shallow, pinned git fetch of the four shipped input files (counts, pixel parameters, PSF, truth)."""
    sums = load_sums()["shipped_inputs"]
    url = os.environ.get("ASTRO_LION_GIT_URL")
    if not url:
        sys.exit("set ASTRO_LION_GIT_URL to the https git URL of the public repository named `lion` (see results/astro/data_access.md)")
    d = os.path.join(root(), "shipped_inputs")
    if os.path.exists(d) and os.listdir(d):
        sys.exit("%s exists and is not empty; each download must go to a new empty directory" % d)
    os.makedirs(d, exist_ok=True)
    env = dict(os.environ, GIT_LFS_SKIP_SMUDGE="1")
    for cmd in (["git", "init", "-q", d], ["git", "-C", d, "remote", "add", "origin", url],
                ["git", "-C", d, "fetch", "-q", "--depth", "1", "origin", sums["commit"]],
                ["git", "-C", d, "checkout", "-q", "FETCH_HEAD"]):
        subprocess.run(cmd, check=True, env=env)
    for rel, want in sums["files"].items():
        dig, n = sha256_file(os.path.join(d, rel))
        print("%-26s %s  %s" % (rel, "OK " if dig == want["sha256"] else "BAD", dig[:16]))
        if dig != want["sha256"]:
            sys.exit("checksum mismatch for %s" % rel)


# ----------------------------------------------------------------------------- patch extraction
def sdss_psf_recon(kl, xpos, ypos):
    """KL-basis PSF image (51x51) at 0-based (column xpos, row ypos); centre pixel = star centre.
    Same polynomial-in-position recipe as the SDSS psField data model."""
    import numpy as np
    rc = 0.001
    rnrow, rncol = int(kl["RNROW"][0]), int(kl["RNCOL"][0])
    nr, nc = int(kl["nrow_b"].max()), int(kl["ncol_b"].max())
    coeffs = np.outer(((ypos + 0.5) * rc) ** np.arange(nr), ((xpos + 0.5) * rc) ** np.arange(nc))
    img = np.zeros(len(kl["RROWS"][0]), dtype=np.float64)
    for i in range(len(kl)):
        c = np.asarray(kl[i]["c"])[: int(kl[0]["nrow_b"]), : int(kl[0]["ncol_b"])]
        img += np.asarray(kl[i]["RROWS"], dtype=np.float64) * (c * coeffs[: c.shape[0], : c.shape[1]].T).sum()
    return img[: rnrow * rncol].reshape(rncol, rnrow)


def cmd_extract(args):
    try:
        import numpy as np
        from astropy.io import fits
        from astropy.wcs import WCS
        from scipy.ndimage import map_coordinates
    except ImportError as e:
        sys.exit("extract-patch needs numpy, scipy and astropy (%s)" % e)
    frame, psfp = dest_path("sdss_frame_r"), dest_path("sdss_psfield")
    for p in (frame, psfp):
        if not os.path.exists(p):
            sys.exit("missing %s; run `fetch --group core` first" % p)
    x0, y0, n = args.x0, args.y0, args.size
    with fits.open(frame, memmap=False) as h:        # astropy reads .bz2 transparently
        hdr = h[0].header.copy()
        nmgy = h[0].data.astype(np.float64)           # (rows, cols), nanomaggies, sky-subtracted
        cal = h[1].data.astype(np.float64)            # nanomaggies per DN, one per column
        t = h[2].data
        allsky, xi, yi = t["ALLSKY"][0].astype(np.float64), t["XINTERP"][0].astype(np.float64), t["YINTERP"][0].astype(np.float64)
    with fits.open(psfp, memmap=False) as h:
        kl = h[3].data                                # r band KL eigen-images
        row6 = h[6].data
        gain, dark_var = float(row6["gain"][0][2]), float(row6["dark_variance"][0][2])
        sky_ps, skysig_ps = float(row6["sky"][0][2]), float(row6["skysig"][0][2])
    yy, xx = np.meshgrid(yi, xi, indexing="ij")
    sky = map_coordinates(allsky, [yy, xx], order=1, mode="nearest")        # DN, bilinear from ALLSKY
    dn_full = nmgy / cal[None, :] + sky                                      # bias-subtracted, flat-fielded DN
    sl = (slice(y0, y0 + n), slice(x0, x0 + n))
    dn, p_nmgy, p_sky = dn_full[sl], nmgy[sl], sky[sl]
    psf = sdss_psf_recon(kl, x0 + (n - 1) / 2.0, y0 + (n - 1) / 2.0)
    psf = np.clip(psf, 0, None)
    psf /= psf.sum()                                                          # 51x51, unit sum
    c = psf.shape[0] // 2
    psf25 = psf[c - 12: c + 13, c - 12: c + 13]
    psf25 = psf25 / psf25.sum()
    out = os.path.join(root(), "patch_r_002583_2_0136_x%d_y%d_%dx%d" % (x0, y0, n, n))
    os.makedirs(out, exist_ok=True)                    # derived products are overwritten
    w = WCS(hdr)
    w.wcs.crpix = w.wcs.crpix - np.array([x0, y0])
    ph = fits.PrimaryHDU(dn.astype(np.float32), header=w.to_header())
    ph.header["BUNIT"] = ("DN", "bias-subtracted flat-fielded counts (frame/cal + sky)")
    for k in ("RUN", "RERUN", "CAMCOL", "FRAME", "FILTER"):
        ph.header[k] = hdr[k]
    ph.header["PX0"], ph.header["PY0"] = x0, y0
    ph.header["GAINR"] = (gain, "e-/DN, psField HDU6")
    fits.HDUList([ph]).writeto(os.path.join(out, "patch_dn.fits"), overwrite=True)
    np.savez(os.path.join(out, "patch_arrays.npz"), dn=dn, electrons=dn * gain, nmgy=p_nmgy, sky_dn=p_sky,
             psf51=psf, psf25=psf25, gain=gain, dark_variance=dark_var, x0=x0, y0=y0)
    json.dump(dict(x0=x0, y0=y0, size=n, gain=gain, dark_variance=dark_var, psfield_sky=sky_ps, psfield_skysig=skysig_ps,
                   note="origin is the lower-left pixel, 0-based (column x0, row y0); array index [row, column]"),
              open(os.path.join(out, "meta.json"), "w"), indent=1)
    exp = load_sums()["derived"].get("patch_x310_y630_stats")
    st = dict(min=float(dn.min()), median=float(np.median(dn)), mean=float(dn.mean()), max=float(dn.max()))
    print("wrote", out)
    print("patch DN stats:", {k: round(v, 3) for k, v in st.items()}, "gain %.3f e-/DN, dark variance %.3f" % (gain, dark_var))
    if exp and (x0, y0, n) == (310, 630, 100):
        ok = all(abs(st[k] - exp[k]) < 1e-3 * max(1.0, abs(exp[k])) for k in st)
        print("regression check against recorded statistics:", "OK" if ok else "MISMATCH %s" % exp)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    for name in ("fetch", "verify"):
        p = sub.add_parser(name)
        p.add_argument("--group", action="append", choices=GROUPS + ["all"])
        p.add_argument("--only", nargs="+")
    sub.add_parser("shipped")
    p = sub.add_parser("extract-patch")
    p.add_argument("--x0", type=int, default=310)
    p.add_argument("--y0", type=int, default=630)
    p.add_argument("--size", type=int, default=100)
    args = ap.parse_args()
    {"list": cmd_list, "fetch": cmd_fetch, "verify": cmd_verify, "shipped": cmd_shipped, "extract-patch": cmd_extract}[args.cmd](args)


if __name__ == "__main__":
    main()
