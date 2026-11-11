# M2 patch benchmark: data access, tools and code locations

Measured on 2026-11-07 from the sandbox, through the session's HTTPS proxy (no proxy or trust setting was changed, TLS
verification was never disabled). Raw data live outside the repository under `$ASTRO_DATA` (default `/home/user/astro_data`),
one new directory per download; nothing in this repository is raw data. Scope: access and setup only, no method is tuned or
claimed. Field facts and paper statements are in `m2_field_spec.md`. All download steps are scripted in `data/astro/fetch.py`
with sha256 values in `data/astro/checksums.json`.

## 1. What is reachable and what is not

| host | result | evidence |
|---|---|---|
| data.sdss.org (SAS dr17, data-model pages) | reachable | HTTP 200, byte ranges supported; files below downloaded and hashed |
| mast.stsci.edu | reachable, with one exception | name resolver, CAOM search, product listing and `api/v0.1/Download/file` all answer. The pseudo-URIs of the form `mast:HLSP/url/cgi-bin/getdata.cgi?...` return the 20-byte body "Unable to read file" |
| arxiv.org | reachable | abstract page, PDF and source tarball of both papers |
| pypi.org | reachable | `pip download pydl` and the install of its wheel worked (crowdsource was installed from the local clone, not from PyPI) |
| github.com, git clone / fetch by commit | reachable for public repositories | shallow clones and `git fetch --depth 1 origin <sha>` worked |
| api.github.com | search works, per-repository calls refused | HTTP 403 "GitHub access to this repository is not enabled for this session" |
| archive.stsci.edu | **fails TLS** | `curl: (60) SSL certificate problem: unable to get local issuer certificate` |
| skyserver.sdss.org | **fails TLS** | `curl: (60) SSL certificate problem: self-signed certificate in certificate chain` |
| hla.stsci.edu, catalogs.mast.stsci.edu, www.astro.ufl.edu | **refused by the egress policy** | CONNECT answered 403 (proxy status page lists each host as `connect_rejected`). Plain-http to hla.stsci.edu does not resolve because only HTTPS is proxied |

I did not look for ways round any of these.

## 2. SDSS data (DR17 science archive server, rerun 301)

Base: `https://data.sdss.org/sas/dr17/eboss/`.

| key in fetch.py | path under the base | bytes | Last-Modified | sha256 |
|---|---|---|---|---|
| `sdss_frame_r` | `photoObj/frames/301/2583/2/frame-r-002583-2-0136.fits.bz2` | 3,524,824 | 2012-06-13 | `9129b6c6779ad8155741998e333e0f69395226aa9fc425dc9933f9f6f53ec04c` |
| `sdss_psfield` | `photo/redux/301/2583/objcs/2/psField-002583-2-0136.fit` | 414,720 | 2009-05-13 | `3d24d270c12dcd8429357a81433d0151d3a6e8739401b18273d4bfc206bcd7bd` |
| `sdss_photofield` | `photoObj/301/2583/photoField-002583-2.fits` | 650,880 | 2013-11-01 | `052132367a8faaee18160b71c3a0beef423259df7b9218c313d564f1d5f6e69e` |
| `sdss_frame_g` | `photoObj/frames/301/2583/2/frame-g-002583-2-0136.fits.bz2` | 3,197,457 | 2012-06-13 | `e543a3ff3d7214f503e3d2ad43f5c4dd80b97a8ff0746fb63b7d203093dcad2a` |
| `sdss_frame_i` | same directory, `frame-i-...` | 3,660,448 | 2012-06-13 | `f10681cc441903b124e7597d84f0bcbc4cf60332e9fa6cfe57def36e79569273` |

HEAD sizes of the other bands in the same directory: `frame-u` 3,156,348 bytes, `frame-z` 3,556,747 bytes.
Decompressed r frame: 12,447,360 bytes, sha256 `09ca40d8f5e7254a0c3c152427b63bd97169a48edc882f4aa18f11412e8ace9f`, byte-identical to the
r (and g, i) frame files shipped in the public code repository (section 5).

How the field was located with SDSS files only: the brief's cluster coordinates and the paper's run/camcol/field were checked against
`photoField-002583-2.fits` (one row per field; columns `RA, DEC, RAMIN, RAMAX, DECMIN, DECMAX`). Row for field 136: centre
RA 323.3429022, Dec -0.7339748; box RA 323.2609545 to 323.4247425, Dec -0.8466525 to -0.6214204. The cluster
(323.36258, -0.82325) and the patch centre (323.3357, -0.8071) are inside. Frame header: `STRIPE = 82` (a stripe with many repeat runs; I did not count
them, so RA/Dec alone does not pick run 2583 and the paper's run was taken as given). The frame WCS (TAN, 0.396 arcsec/pixel) puts the
cluster centre at (col 212.5, row 923.9) and a smoothed-image peak at (211, 925), see the figure note below.

Not in the SAS tree for this run (checked): a calibrated `fpC` image (`.../corr/2/fpC-002583-r2-0136.fit.gz` returns 404; the run directory
lists only `astrom, dr8_astrom, dr8_logs, logs, nfcalib, objcs, photo, photobot, psFangs, resolve, ssc, window`) and the raw `idR` image.
The shipped counts array (section 5) was cut from a raw image; it is the only copy I have of those exact counts.
Photo products for field 136 are empty: `photoObj-...-0136.fits` 8,640 bytes, `fpAtlas` 8,640, `psBB` (r) 8,640, `fpObjc` 37,440, `fpM` (r) 152,738.
The psField for field 136 is intact (the paper's pipeline PSF comes from it).

Frame to counts (SDSS data model for `frame` files, applied by `fetch.py extract-patch`): the image is sky-subtracted and in
nanomaggies; HDU1 is the calibration vector (nanomaggies per DN per column), HDU2 holds a 256 x 192 sky image with `XINTERP`
and `YINTERP`; DN = image / calibration + bilinear sky; electrons = DN x gain; DN error = sqrt(DN / gain + dark variance).
psField HDU3 (r) holds four 51 x 51 KL eigen-images; the PSF at (x, y) is their sum weighted by polynomials in (0.001 (y+0.5),
0.001 (x+0.5)); HDU6 carries per-filter `gain`, `dark_variance`, `sky`, `skysig`.
Documented r-band camcol-2 values: gain 4.6 e-/DN, dark variance 1.00 (data-model table); psField agrees (4.6, 1.0).

Verification figure (not committed; made by a scratch script outside the repository): frame overview with the cluster centre and the patch box,
the patch with the HST sources overlaid, and a sigma = 15 px smoothed core test (peak 0.7 arcsec from the predicted centre). Scratch copy at
`/home/user/astro_data/derived_patch_r_002583_2_0136/m2_patch_check_v2.png`.

## 3. HST catalogues

### 3.1 What MAST holds for NGC 7089

* CAOM search by position (0.1 deg around M2, collection HLSP) returns 397 records: 387 TESS light curves, 9 SKYSURF drizzled-image
  records, and one `hlsp_acsggct` record (ACS/WFC F606W, programme 10775, exposure 340 s, data rights PUBLIC, footprint RA 323.3204 to 323.4037,
  Dec -0.8649 to -0.7816). The whole `hlsp_acsggct` project is 71 F606W image records; **no photometric catalogue record**.
* The product list of that record is one 144,002,880-byte image FITS and one preview. Both URIs are the unreadable `getdata.cgi` pseudo-URIs and the
  `dataURL` points to hla.stsci.edu (refused). Guessed paths of the form `mast:HLSP/acsggct/...` return 404.
  **The ACS Survey of Galactic Globular Clusters photometric catalogue that the papers used cannot be downloaded from the reachable hosts.**
* HUGS (the later HST UV legacy survey of Galactic globular clusters) is not in CAOM near M2, but its files are served by the Download endpoint under
  `mast:HLSP/hugs/ngc7089/`. Both photometry methods were downloaded (4 to 5 s each).

### 3.2 HUGS NGC 7089 catalogues

| key | file | bytes | Last-Modified | sha256 |
|---|---|---|---|---|
| `hugs_ngc7089_meth1` | `hlsp_hugs_hst_wfc3-uvis-acs-wfc_ngc7089_multi_v1_catalog-meth1.txt` | 106,805,615 | 2018-06-29 08:54:25 | `4a96ae437122468785200570dd9121f4b677ec6ed0200e65a687babcf5b58f49` |
| `hugs_ngc7089_meth2` | same name with `meth2` | 106,805,615 | 2018-06-29 08:54:01 | `19b602780a25fb03d9eb5167aa75e02386d977bb175dd8ed7db3b80918e84dc2` |

Format: ASCII, 54 comment lines then 294,222 rows of 37 whitespace-separated columns (all rows have 37 fields). Header text: the
catalogue is built from programmes GO-10775 (the ACS survey) and GO-13297 (WFC3/UVIS), and asks users to cite two HUGS papers
(a 2015 overview and the 2018 catalogue paper). No licence text. Columns (1-based): 1-2 X, Y in a reference frame whose cluster
centre is at (5000, 5000) (pixel scale measured here: 0.0395 arcsec/pixel, X toward west, affine fit rms 0.001 arcsec; the header
does not state the scale); then for each of F275W, F336W, F438W, F606W, F814W six columns (calibrated magnitude, rms, quality-of-fit,
sharpness, number of exposures the source is found in, number well measured; `99` = saturated, magnitude `-99.9999` = not measured);
33 membership probability (`-1.0` = not available, 130,890 rows have a value); 34-35 RA, Dec in degrees, J2000, epoch 2015.0,
astrometrised with Gaia DR1; 36 identifier (string like `R0000001`); 37 detection iteration (1-5 found in F606W and F814W images,
6 F438W, 7 F336W, 8 F275W; only 1 to 6 occur in this file). The header does not state the magnitude system.

Content: method 1 has 287,653 rows with a valid F606W magnitude (method 2: 291,917); 95,672 (95,051) brighter than F606W = 22.
RA 323.327309 to 323.397115, Dec -0.859411 to -0.787322. Patch counts (box [309.5, 409.5) x [629.5, 729.5) after mapping RA/Dec with the
frame WCS): 5,414 rows, 5,394 with F606W, 1,070 brighter than 22 (method 2: 5,413, 5,407, 1,071). Read it with pandas using
`sep=r'\s+', comment='#', header=None, dtype={35: str}` (plain `numpy.loadtxt` fails on the string identifier).

### 3.3 Shipped truth file (from a public git repository, section 5)

`Data/sdss0921_true.txt`: 6,048 rows `x y flux` in patch pixel coordinates (0-based, pixel centres at integers, rows span -0.498 to 99.448 in x) and SDSS counts;
magnitude = 22.5 - 2.5 log10(flux x 0.00546689). Counts: 591 brighter than 21, 1,069 brighter than 22, 1,611 brighter than 23, 6,044
brighter than 30. It is the deepest HST-derived list I can obtain and is the table the shipped metric code evaluates against. Its origin
table (the raw HST catalogue in SDSS frame pixels that the shipped code reads as `hubble_catalog_2583-2-0136_astrans.txt`) is not in the repositories.

## 4. Environment

Pinned repository venv: untouched (numpy / scipy / matplotlib / pillow only); `fetch.py` runs in it using the standard library.

Domain venv (`scratchpad/venv_domain`, Python 3.13.18):

* The PyPI package named `crowdsource` (an unrelated web application) was **uninstalled** (`pip uninstall -y crowdsource`, version 0.1.0 removed, no leftover directory).
* The astronomy package from the read-only clone (version 0.5.9, package name on PyPI `crowdsourcephoto`, import name `crowdsource`, MIT licence,
  clone commit `e648ee9b7a5d80701788efdbe88a58de8af72910` dated 2026-10-16, working tree clean before and after) was installed with
  `pip install --no-deps` from a copy of the tree made without `.git`, so the clone stayed untouched. `--no-deps` is deliberate: `setup.py`
  demands `tensorflow >= 2`, `keras`, `scikit-image` and `guppy3`, but only `nebulosity_mask.py` (keras), `deconv_psf.py` (scikit-image) and
  `decam_proc.py` (guppy3) import them; `crowdsource_base.py` and `psf.py` need only numpy, scipy and astropy. `import crowdsource` and
  `import crowdsource.crowdsource_base` succeed and resolve to the new package.
* `pydl 1.0.1` (BSD-3, pure Python) was installed with `--no-deps` for its psField reconstruction function; `fetch.py` carries its own
  15-line implementation of the same recipe (maximum difference to pydl 1.1e-8 on the unit-sum PSF).
* photutils 3.1.1, astropy 8.0.2 were already present.

Inputs for the smoke tests (both use `fetch.py extract-patch` data): the 100 x 100 patch in DN, the psField r-band PSF at the patch centre
(51 x 51 stamp, unit sum, FWHM 2.22 pixel), per-pixel error `sqrt(DN / 4.6 + 1.0)`. Single process, `OMP_NUM_THREADS=1`, `nice`.

**photutils (DAOStarFinder + IterativePSFPhotometry, runs):** background = 1st percentile of the patch (184.98 DN), `DAOStarFinder(threshold = 5 x psField skysig = 32.11 DN,
fwhm = 2.22)`, `IterativePSFPhotometry(ImagePSF(psf, origin=(25, 25)), fit_shape=(7, 7), finder, grouper=None, maxiters=3, mode='new', aperture_radius=3)`.
DAOStarFinder alone: 98 detections, 0.10 s. IterativePSFPhotometry: 206 sources in the result table (98, 63, 45 per iteration), 176 with `flags == 0`,
0.68 s (first run 1.27 s). API as of photutils 3.1.1: `FittableImageModel` is not in the `photutils.psf` namespace any more, `ImagePSF` was used; the analytic models whose names end in `PSF` are documented as unsuitable for fitting (use the `PRF` variants).

**crowdsource (`fit_im`, runs):** inputs as in its `simple_proc.py` (`weights = 1/sigma`, zero flag image, `maxstars = 10000000`, `fewstars = 60`, `refit_psf = False`), with `ntilex = ntiley = 1` because the patch is tiny.
With the SDSS PSF (`SimplePSF`, 51 x 51 stamp): 563 catalogue rows in 3.16 s. With the Gaussian default of `simple_proc.py` (FWHM 2): 854 rows in 5.29 s.
Fields returned: `x, y, flux, passno, dx, dy, dflux, qf, rchi2, fracflux, fluxlbs, dfluxlbs, fwhm, spread_model, dspread_model, fluxiso, xiso, yiso, flags, sky`.
Problems found, none needing a trust or system change:
1. The flag image must be 64-bit. With `int32` the call stops at `crowdsource_base.py` line 242 with `OverflowError: Python integer 2147483648 out of bounds for int32` under NumPy 2.5.4.
2. `crowdsource/simple_proc.py` cannot be run as a script: it imports a module `mosaic` that is neither shipped nor declared, and it expects the older tuple return value, whereas `fit_im` now returns a dict (`stars`, `model`, `sky`, `psfs`, `iter_history`). I called the same `fit_im` directly.
3. The README documents installation and algorithm notes only; no end-to-end usage example.
4. `fit_im` prints per-iteration centroid lines even with `verbose=False`.
Larger regions (same inputs, SDSS PSF, single core): with the `ntilex = ntiley = 4` that `simple_proc.py` hard-codes, none finished: the whole frame (2048 x 1489) was stopped by me after about 13 minutes
(12 min 54 s of CPU, no iteration completed, peak RSS 2.8 GB); a 600 x 600 cutout (columns 60-659, rows 380-979, includes the core) hit a 14 minute limit; a 300 x 300 cutout (columns 210-509, rows 530-829) hit a 7 minute limit.
The same 300 x 300 cutout with `ntilex = ntiley = 1` finished in 66.9 s (10 iterations, 4,794 rows, 558 of them inside the 100 x 100 patch box, peak RSS 740 MB). With these inputs the 4 x 4 tiling is therefore the practical blocker for regions larger than the patch (cause not investigated; the Gaussian-PSF variant was not tried on large regions).
Counts above are for the run only; no statement about their quality is made.

## 5. PCAT code: where it is and how it is run

PCAT code is public on github.com, but the owner handles are the original authors' personal accounts and are left out of this file
(no-author-names rule); they are given in the hand-off report, and `fetch.py shipped` reads the URL from `ASTRO_LION_GIT_URL`.

| repository name | what it is | last commit | size | licence file |
|---|---|---|---|---|
| `lion` | the fast single-band implementation named in the follow-up ("Lion"), with the M2 inputs in `Data/`, repository history from 2017-04-27 | 2026-09-08 (HEAD `5f0f3257490a351dc9ce1745aac17b16383ff0f3`); last commit before that burst 2021-01-13 (`07be01f39dfa6995eba51bf62e6dc63550192645`) | 38 MB shallow | none at root |
| `multiband_pcat` (the repository named in the follow-up; the plain owner name and a `-zz` form of the owner name both resolve to this HEAD, checked with `git ls-remote`) | multiband extension, the M2 benchmark data directory (frames, counts, PSFs, truth, DAOPHOT photometry), metric scripts | 2022-08-15 (`fa32cda6acb29173afaa2c3eaa3f68f38304eaae`) | 174 MB shallow | none at root |
| `pcat` (general framework, MIT) | modern general metamodel sampler; its "M2" example is a "scaled synthetic analog" (8 sources, 24 x 24 pixels), **not** the benchmark | 2026-10-27 | 280 MB shallow | MIT |
| `PCAT` (2017 notebook) | proof of concept with NUTS, unrelated to M2 | 2017-04-06 | 13 MB | none |

Three more public repositories turned up in the search and were not inspected: `pcat-de` (MIT), a 2016 pipeline for probabilistic catalogues, and a C++ code for the earlier probabilistic crowded-field approach (GPL-3.0).
The 2017 paper itself says the code was not yet released ("we are currently working on computational efficiency so that the sophisticated user might be presented with code");
the follow-up says "the code used in this work, Lion 0.1, is publicly available on GitHub".

Shipped inputs (identical bytes in `lion/Data/sdss0921_*` and `multiband_pcat/Data/sdss.0921/*`; the four hashes are in `checksums.json`):
`cntp` 100 x 100 raw counts (bias included, rows = frame rows, origin (310, 630)); `pixl` (`100 100 1`, then `1044. 4.62`); `psfn` (`25 5` then a 125 x 125 array: the 25 x 25 psField PSF
at the patch centre upsampled 5 times); `true` (6,048 truth rows). The multiband repository also holds the g/i/z frames, per-band counts, pixel and PSF files,
asTrans grids for the patch (`...-0310-0630`), a 356-row DAOPHOT catalogue in patch pixels (x, y, then fluxes in counts; the third column equals the r flux from the table's r magnitude, 355 of 356 rows match the photometry table within 1 pixel), and `m2_2583.phot` (7 MB table of the SDSS-run photometry of the M2 field; the metric script reads RA, Dec, pixel x, y, r magnitude and error from columns 5, 6, 21, 22, 23, 24 as the DAOPHOT catalogue).

How they are run (read from the code, **not run**):

* `multiband_pcat/pcat.py`: Python 2 source (`print` statements, `xrange`, `raw_input`, `time.clock`; `ast.parse` under Python 3.13 fails at line 19); numpy, scipy (`KDTree`), matplotlib (TkAgg backend), astropy FITS;
  a C shared library `pcat-lion.so` built from `pcat-lion.c` (Makefile: `icc -mkl -shared -static-intel -liomp5 -fPIC -O2`; README alternative: OpenBLAS with `gcc -shared -fPIC`);
  command line `python pcat.py <dataname> <verbtype> <testpsfn> <datatype> <multiband>` with interactive prompts for the band names and a comment;
  hard-coded base paths of the authors' cluster; random seed 20170501; `nsamp` 300 or 1500 saved samples, `nloop = 1000` steps between samples (about 1.5e6 steps), margin 10 pixels, minimum flux 236 ADU (22.22 mag).
  The metric scripts in `sdss/` (`m2plots.py`, `helpers.py`) are Python 2 too.
* `lion/main.py`: parses as Python 3; imports numpy, scipy, h5py, matplotlib, astropy and `tdpy` (a third-party utility package); `blas.so` is built with `gcc -mavx -shared -fPIC -O2 blas.c -o blas.so`;
  environment variables `LION_PATH` (and `LION_DATA_PATH`); started as `python main.py <function> [args]`; `main()` defaults to `strgdata='sdss0921'`, `numbswep=100`, `numbloop=1000`; the supplied `cnfg_defa`
  reads a mock HDF5 file that is not in the repository. Its README says it is "research-oriented rather than a turnkey black-box package".
* Runtime claims (not tested here): 2017 implementation, a day of wall clock on 12 cores for the patch (about 1e9 model evaluations, 12 CPU-days); `Lion 0.1`, half an hour on one CPU for 1e6 steps.

## 6. Licences and terms stated by the hosts

| source | what is stated |
|---|---|
| SDSS SAS pages and data-model pages | no licence, policy or acknowledgement text found on the three pages I read (SAS dr17 root, `frame` and `psField` data models; zero matches); the SDSS citation and acknowledgement policy is on a different site that I did not fetch |
| MAST | the landing page is a script shell without terms text; the CAOM record of the ACS survey image lists `dataRights: PUBLIC`; the HUGS header asks for citation of two papers and states no licence |
| arXiv | both abstract pages link the non-exclusive distribution licence |
| public code repositories | `lion` and `multiband_pcat` have no licence file at the repository root, so no licence is granted for their code or data: use them for reading and local checking, and do not redistribute their files from this repository; the general `pcat` repository has an MIT licence file |
| software | crowdsource MIT (licence file in the clone); photutils and astropy BSD-3-Clause; pydl BSD 3-Clause (package metadata) |

## 7. `data/astro/fetch.py`

```
export ASTRO_DATA=/home/user/astro_data                 # default; always outside the repository
python data/astro/fetch.py list
python data/astro/fetch.py fetch --group core            # r frame + psField (3.9 MB)
python data/astro/fetch.py fetch --group core --group geometry --group multiband
python data/astro/fetch.py fetch --group hst             # two 107 MB HUGS catalogues
python data/astro/fetch.py verify --group all
python -I data/astro/fetch.py extract-patch              # needs numpy, scipy, astropy (domain venv)
ASTRO_LION_GIT_URL=<https url of the repository named lion> python data/astro/fetch.py shipped
```

Behaviour: standard library only; downloads to `$ASTRO_DATA/<key>/<filename>.part` with HTTP Range resume (a stream that ends before Content-Length is resumed, not accepted) and renames it when complete; refuses to continue
after a 401/403/404; then verifies sha256 and size against `checksums.json` and moves a bad file to `*.bad`; records a new key on first download; `shipped` uses a shallow fetch pinned to the 2021 commit and refuses a non-empty target.
Tested in this session: fresh downloads of all five SDSS files and one HUGS catalogue matched the recorded hashes; an interrupted `.part` resumed with HTTP 206 at 100,000 bytes (and a truncated stream from a local test server resumed at 40,000 of 100,000 bytes); a corrupted copy was
reported BAD (exit 1) and moved aside; `shipped` fetched and verified the four files; `extract-patch` reproduced the patch bit-for-bit (statistics regression check printed OK).

## 8. Blockers and doubts for reproducing the benchmark

1. The ACS-survey photometric catalogue used by the papers is not reachable (MAST holds only pointer records that lead to hla.stsci.edu, archive.stsci.edu fails TLS). Stand-ins that are reachable: HUGS via MAST and the shipped truth file via git. They agree with each other on positions to 0.05 pixel and on the count brighter than 21 mag (591) but not in the faint tail, and neither reproduces the paper's 1,049 exactly (1,051 to 1,070 depending on the box).
2. The patch is fixed by the follow-up's origin (310, 630) and confirmed by the shipped counts (correlation 0.99993); the paper's quoted patch-centre RA/Dec is 1.42 arcsec (3.6 pixel) away from the frame WCS and is unexplained.
3. Counts: the shipped array comes from a raw image that is not on the SAS; calibrated-frame counts differ at the percent level (mean 10.5 DN on a median of 310 DN). Gain 4.6 (psField, data model) versus 4.62 (papers' code) is a 0.4% difference.
4. The DAOPHOT comparison catalogue exists only as shipped files (I did not look for its original table), and its pixel positions sit about 1.06 pixel from the shipped HST truth in y (only 25% of bright DAOPHOT sources are within 0.75 pixel of an HST source as shipped, 94% after removing the shift); see `m2_field_spec.md` item 10.
5. PCAT: the public code that matches the benchmark is Python 2 with a compiled library and hard-coded paths (`multiband_pcat`), or a refactored Python 3 version without a working real-data entry point (`lion`); neither was run, and I did not find the 2017 implementation described in the first paper: the public code is the later Lion line.
6. Neither the bin edges of the 2017 figures nor the 1,049 selection are in the text; the shipped metric script (17 bins of 0.5 mag from 15.0 to 23.5, 0.5 pixel and 0.5 mag match, 2.5-pixel margin) is the only definition found.
7. The recent history of `lion` ends in a burst of four generic-message commits dated 2026-09-08; the data files are byte-identical before and after, which is why the 2021 commit is pinned.
