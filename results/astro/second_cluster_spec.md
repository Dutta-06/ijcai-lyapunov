# A second globular cluster as a test set: SDSS r band against HST (HUGS), M53 and M92, validated on M2

Scope: coverage survey, astrometry and photometry validated on the M2 benchmark, and six 100 x 100 pixel test patches of two other clusters in the on-disk format of the M2 patches. No cataloguing method (the
sliding inverse scale space, photutils, crowdsource) was run on any new patch; the only fits to patch pixels are the rigid-shift scans of section 4.4, which fit a PSF-convolved copy of the truth list and
detect nothing. Everything below was measured on 2026-11-07 through the session's HTTPS proxy (data.sdss.org, mast.stsci.edu; no other host was contacted, no refusal was met, TLS verification was never
disabled). Numbers are printed by `experiments/astro_second.py` (`survey`, `validate-m2`, `build`), `data/astro/fetch_second.py` (`probe`, `resolve`, `fetch`, `verify`) or are fields of
`data/astro/checksums_second.json` and of the per-cluster JSON that `build` writes; commands are in section 9. Raw data live outside the repository under `$ASTRO_DATA_2` (default `/home/user/astro_data_2`).
Numbers marked "scratch" come from exploratory scripts of this session that are not delivered. Field facts of M2 and the earlier access notes are in `m2_field_spec.md` and `data_access.md`.

## 1. Result and the caveats that matter

Test set (all patches: SDSS rerun 301, r band, camcol 6; directory names as in M2, `patch_r_<run>_<camcol>_<field>_x<x0>_y<y0>_100x100`; truth `truth_hugs_<cluster>_meth1.npz`, frame pixels):

| patch | cluster | frame (run, camcol, field) | origin (x0, y0) | F606W < 22 in the box | per pixel |
|---|---|---|---|---|---|
| `patch_r_005390_6_0162_x760_y590_100x100` | M53 (NGC 5024) | 5390, 6, 162 | (760, 590) | 402 | 0.0402 |
| `patch_r_005390_6_0162_x640_y680_100x100` | M53 | 5390, 6, 162 | (640, 680) | 671 | 0.0671 |
| `patch_r_005390_6_0162_x400_y660_100x100` | M53 | 5390, 6, 162 | (400, 660) | 950 | 0.0950 |
| `patch_r_005390_6_0162_x660_y420_100x100` | M53 | 5390, 6, 162 | (660, 420) | 1302 | 0.1302 |
| `patch_r_004682_6_0208_x1100_y220_100x100` | M92 (NGC 6341) | 4682, 6, 208 | (1100, 220) | 718 | 0.0718 |
| `patch_r_004682_6_0208_x1430_y550_100x100` | M92 | 4682, 6, 208 | (1430, 550) | 1184 | 0.1184 |

* **Primary set: M53, four patches** at 0.040 to 0.130 per pixel. **Secondary set: M92, only two patches**: the rules of section 2.4 leave two non-overlapping places in its HUGS footprint, none below 0.07 per pixel, so it does not meet the
  "at least three patches" criterion by itself. M3 (one patch) and M5 (none) were examined and rejected by the same rules (section 2.3).
* **The seeing is much worse than M2's.** The FWHM of the psField PSF at the patch centre is 3.27 to 3.40 px (M53) and 2.94 to 2.99 px (M92) against 2.12 to 2.16 px for M2 (same estimator, section 5.3), and the
  peak signal to noise of an isolated F606W = 22 source is 4.2 to 4.4 (M53), 5.2 to 5.3 (M92) and 6.3 to 6.5 (M2). F606W = 22 on M2 corresponds to F606W = 21.5 (M53) and 21.8 (M92) at equal peak signal to noise.
  The frames are a test of transfer across seeing and depth as well as across fields; a comparison of completeness at 22 mag between the sets is not like for like.
* **The mapping of the new truth is the full asTrans solution, not the header TAN of M2's `make_truth.py`**: on the M53 patches the header mapping is off by up to 0.14 px and the image prefers asTrans (section 4.3).
  Both truth lists are written (`..._meth1.npz`: asTrans, `..._meth1_header.npz`: header).
* The truth lists contain a handful of sources too few per patch (F606W not measured: 2 to 3 sources with F814W < 22 per patch, 0 or 1 with F814W < 20) and the F606W of HST-saturated stars (flag 99) is not reliable at the
  bright end; the selection avoids the worst cases (section 2.4) but does not remove the effect (section 5.4).
* The truth is HUGS epoch 2015.0, the SDSS frames are 2004 and 2005: a net offset of 0.05 to 0.19 px between truth and image is measured and not corrected (section 4.4).

## 2. Coverage survey

### 2.1 Which HUGS catalogues exist and are downloadable

`fetch_second.py probe` asks the MAST Download endpoint for `mast:HLSP/hugs/<cluster>/hlsp_hugs_hst_wfc3-uvis-acs-wfc_<cluster>_multi_v1_catalog-meth1.txt` with a 1-byte range request for 65 NGC numbers: 52 answer
HTTP 206 and give the size, 13 answer 404. The 404 are NGC 104, 288, 362, 1904, 2419, 4147, 5139, 5634, 6229, 6402, 6426, 6712, 7492 (another file name, or not a HUGS cluster; I did not investigate, so this is a possible gap
of the survey). Answering 206: NGC 1261, 1851, 2298, 2808, 3201, 4590, 4833, 5024, 5053, 5272, 5286, 5466, 5904, 5927, 5986,
6093, 6101, 6121, 6144, 6171, 6205, 6218, 6254, 6304, 6341, 6352, 6362, 6366, 6388, 6397, 6441, 6496, 6535, 6541, 6584, 6624, 6637, 6652, 6656, 6681, 6715, 6717, 6723, 6752, 6779, 6809, 6838, 6934, 6981, 7078, 7089, 7099.
The 14 catalogues I downloaded (the clusters with a candidate SDSS run, below, and M2) all parse with the column indices of the M2 file (seven of them were also checked to have 54 comment lines and 37 columns). Their header names GO-10775 (PI Sarajedini) and
GO-13297 (PI Piotto) for M53, M92 and M2, "Photometric Method: 1", positions "J2000, epoch 2015.0, astrometrised using Gaia DR1".

### 2.2 Which have SDSS r-band imaging over their footprint (SDSS tables only, SkyServer not used)

Tables used, all on the SAS (`https://data.sdss.org/sas/dr17/eboss/`): `photoObj/photoRunAll-dr17.fits` (765 runs, all rerun 301, flavour science: node, inclination, `MU_REF`, `FIELD_REF`, first and last field of each run) and
`photoObj/301/<run>/photoField-<run6>-<camcol>.fits` (one row per field: great-circle box `MU_START`..`NU_END`, r-band `PSF_WIDTH`, `PSF_NSTAR`, `PSP_STATUS`, `SATURATION_LEVEL`). The 352 MB `resolve/2013-07-29/window_flist.fits` exists
but was not needed. Method (`experiments/astro_second.py survey`): the cluster centre (MAST name lookup, SIMBAD, recorded by `fetch_second.py resolve`) is converted to the great-circle coordinates of every run;
a run is a candidate when the cluster lies within 1.5 deg of its scan circle and inside its field range; for every candidate run the six camcol photoField tables give the field boxes, the fraction of the F606W < 22 HUGS
sources inside the field box, the r-band seeing and the number of 100 x 100 windows (10 px grid, 0.03 to 0.15 sources brighter than 22 per pixel, 30 px from the box edge, inside the HUGS footprint) in it. The best field
of a cluster is the one with the largest footprint fraction, then the smallest seeing. The frame WCS then confirms the geometry: the cluster centres fall at (x, y) = (556.9, 547.7) in the M53 frame and (1330.3, 454.8) in the M92 frame.

`seeing` = r-band `PSF_WIDTH` of the field (arcsec; M2's field 136: 1.048). `PSP_STATUS` is "the maximum value of status over all 5 filters" (the data-model page says no more); `PSF_NSTAR(r)` is the number of PSF stars of the field.

| NGC | name | candidate runs | best frame (run, camcol, field) | seeing | footprint fraction | windows | density range of windows | PSP_STATUS / PSF_NSTAR(r) |
|---|---|---|---|---|---|---|---|---|
| 5024 | M53 | 7 | 5390, 6, 162 | 1.317 | 1.00 | 929 | 0.03 to 0.15 | 0 / 14 |
| 5053 | | 4 | 5390, 5, 168 | 1.204 | 0.95 | 0 | none in the range | 0 / 29 |
| 5272 | M3 | 2 | 4646, 3, 81 | 1.322 | 0.94 | 199 | 0.07 to 0.15 | 2 / 5 |
| 5466 | | 4 | 4623, 1, 324 | 1.646 | 1.00 | 71 | 0.03 only | 0 / 28 |
| 5904 | M5 | 4 | 1458, 4, 700 | 1.135 | 1.00 | 581 | 0.06 to 0.15 | 3 / 0 |
| 6205 | M13 | 5 | 3226, 5, 127 | 1.575 | 0.69 | 34 | 0.10 to 0.15 | 2 / 4 |
| 6254 | M10 | 3 | 5391, 5, 185 | 1.268 | 0.50 | 286 | 0.06 to 0.15 | 2 / 5 |
| 6341 | M92 | 6 | 4682, 6, 208 | 1.345 | 1.00 | 716 | 0.06 to 0.15 | 0 / 16 |
| 6535 | | 2 | 1463, 3, 61 | 1.558 | 1.00 | 87 | 0.03 to 0.04 | 3 / 0 |
| 6838 | M71 | 2 | 6895, 3, 52 | 1.355 | 0.52 | 82 | 0.03 to 0.04 | 3 / 0 |
| 6934 | | 2 | 1738, 1, 29 | 1.027 | 0.79 | 452 | 0.03 to 0.15 | 2 / 1 |
| 6981 | M72 | 1 | none: no field box holds 2 % of the footprint (the cluster falls in a camcol gap) | | | | | |
| 7078 | M15 | 4 | 1749, 6, 141 | 1.744 | 1.00 | 338 | 0.07 to 0.15 | 2 / 3 |
| 7089 | M2 | 32 | the reference (benchmark frame 2583, 2, 136); not part of the selection | 1.048 | | | | 2 / 12 |

The other 38 HUGS clusters answering 206 have no candidate run (the scan strips of rerun 301 do not cross their position): NGC 1261, 1851, 2298, 2808, 3201, 4590, 4833, 5286, 5927, 5986, 6093, 6101, 6121, 6144, 6171, 6218,
6304, 6352, 6362, 6366, 6388, 6397, 6441, 6496, 6541, 6584, 6624, 6637, 6652, 6656, 6681, 6715, 6717, 6723, 6752, 6779, 6809, 7099. Northern clusters with SDSS imaging over (nearly) the whole HUGS footprint in one
frame: M53, M92, M5, M15 and NGC 5466, NGC 6535 (fraction 1.00); NGC 5053 (0.95) and M3 (0.94, the cluster centre near the camcol edge); M13 only 0.69 in its best frame.

### 2.3 What the placement rules leave (`build`, rules of section 2.4; origins on the 10 px grid, after each rule in turn)

| cluster, frame | density range | frame margin | HUGS footprint | saturation | bright star | lost star | patches |
|---|---|---|---|---|---|---|---|
| M53, 5390-6-162 | 1737 | 1737 | 934 | 786 | 246 | 206 | **4** |
| M92, 4682-6-208 | 1857 | 1857 | 711 | 417 | 78 | 37 | **2** |
| M3, 4646-3-81 (dry run) | 1311 | 1219 | 209 | 162 | 10 | 5 | 1 (0.128 per pixel) |
| M5, 1458-4-700 (dry run) | 1873 | 1873 | 570 | 227 | 0 | 0 | 0 |

Candidates per density class [0.03, 0.05), [0.05, 0.08), [0.08, 0.11), [0.11, 0.15]: M53 41, 99, 42, 24; M92 0, 5, 30, 2 (the M92 candidates form two groups, near (1100 to 1180, 220) and (1430 to 1500, 490 to 550),
each smaller than one patch, so at most two non-overlapping patches). M5's closer cluster has a bright star (F606W < 14.5) within 10 px of the box of every origin that survives the saturation rule.

### 2.4 Choice of cluster and of patches

M53 first because it is the only frame that leaves four patches in four density classes, with the whole footprint inside one frame (rows 220 to 882, columns 224 to 907 of the 2048 x 1489 frame), a measured PSF
(`PSP_STATUS` 0, 14 stars) and `QUALITY` 1; M92 second for the same reasons (footprint columns 1050 to 1614, rows 151 to 745; 16 PSF stars). M5 has the best seeing of the frames that cover their whole footprint (1.135) but its PSF is extrapolated (status 3, no star in
the field) and no patch survives. The rules (module docstring of `astro_second.py`; origins on a 10 px grid of integer pixels, box [x0 - 0.5, x0 + 99.5)^2, density = truth sources with F606W < 22 in the box / 10^4):
(1) 40 px from the frame edges; (2) inside the HUGS footprint: every 20 x 20 px cell of the box grown by 20 px holds a truth source; (3) no pixel above 0.8 x the photoField r `SATURATION_LEVEL` (52,900 DN for both
frames) within 10 px of the box; (4) no truth source brighter than F606W = 14.5 within 10 px; (5) no HUGS source with F814W < 17 whose F606W is not measured within 10 px; (6) density 0.03 to 0.15 per pixel; one patch per
class [0.03, 0.05), [0.05, 0.08), [0.08, 0.11), [0.11, 0.15], the candidate closest to the class target (0.04, 0.065, 0.095, 0.13), no overlap with a patch already chosen, the scarcest class first, ties by (y0, x0).
The current M2 test patches span 0.103 to 0.231 per pixel (section 5.2); the range 0.03 to 0.15 is the brief's.

History, for honesty: the first selection used rules (1), (2), (3), (6). Before any method was run, inspection of the PSF-convolved truth model against the data (the machinery of section 4.4) showed two defects of the truth in
chosen patches, and rules (4) and (5) were added and the selection rerun from scratch; those two rules were not adjusted afterwards. The truth mapping was also changed from the header TAN to asTrans after the rigid-shift scan
(section 4.3); the patches were selected again with the asTrans truth. The defects (exploratory scratch scripts, not delivered): a M92 patch contained a
HUGS source F606W = 12.18 flagged saturated (`nf606` = 99) at which the SDSS image peaks at only 6,237 DN, so its magnitude is not a measurement; a M53 patch contained a bright star (F814W = 14.39, F606W not measured, flag 99) that
is in the image and not in the truth list. M13, M5, M3 were not tried with different thresholds. No detector output was used at any point.

## 3. Inputs

All files are fetched by `fetch_second.py fetch`, sha256 and size recorded in `data/astro/checksums_second.json` on the first download and verified afterwards (`verify --group all`: 260 files OK). The four M2 files were fetched
again and are byte-identical (same sha256 and size) to those recorded by `fetch.py` in `checksums.json` (flag `matches_data_astro_checksums_json`). Each download has its own directory `$ASTRO_DATA_2/<cluster>/<key>/`.
Resume was tested against the real host (a 100,000-byte `.part` of `photoRunAll-dr17.fits` resumed with HTTP 206 at byte 100,000 and verified) and a corrupted copy was reported as a checksum mismatch and moved to `.bad`; the same
two behaviours are in `tests/test_astro_second.py` against a local server.

| key | what | bytes | Last-Modified of the host | sha256 |
|---|---|---|---|---|
| `photorunall` | survey table | 187,200 | 06 Apr 2016 04:39:39 | `ff2bb31f3f5bbfda7b2494f11c07f633053da7cc9c540c454635b81c99b830f0` |
| `sdss_frame_r_005390_6_0162` | M53 r frame | 3,273,246 | 24 Jun 2012 08:36:41 | `7e0ef80e37bb3d99664a2a30d9124d961b2c08ad80a8f13cdf508b87387e959e` |
| `sdss_psfield_005390_6_0162` | M53 psField | 391,680 | 05 Jun 2009 21:32:17 | `d0abe86685835ae19508fdbb09c97f1fc8f97e2579981932f0d858ee7caf5a0b` |
| `sdss_photofield_005390_6` | M53 photoField (camcol 6) | 538,560 | 31 Oct 2013 13:30:50 | `6c230d7f3d7c6ce7537693a1370a873e82130ed73200856ec8a25ddf87716015` |
| `hugs_ngc5024` | HUGS NGC 5024 (M53) method 1 | 94,152,161 | 29 Jun 2018 08:12:18 | `a00427f80528109795f94f59ca94aaa603147ef1c5493befa3655ca40779a257` |
| `sdss_frame_r_004682_6_0208` | M92 r frame | 3,614,530 | 21 Jun 2012 06:31:50 | `b293fc2998b370da7c52c64ee6c154f2aac3690f75d5d36d4d9f5f2ce9fcb2e5` |
| `sdss_psfield_004682_6_0208` | M92 psField | 455,040 | 04 Jun 2009 23:42:56 | `ab488c25fb55e954e2326bd61f1b2c636dc970e62d65167940dff3f41f67a06c` |
| `sdss_photofield_004682_6` | M92 photoField (camcol 6) | 792,000 | 31 Oct 2013 23:27:13 | `58e9d0253fb7b5e8527660156dd9e6e0f20c2a757b8a208efba93b9504b79572` |
| `hugs_ngc6341` | HUGS NGC 6341 (M92) method 1 | 57,436,526 | 29 Jun 2018 08:27:52 | `85e300a56de95e6b6deeb091cc70661ec4b16d99c9da8bc8874208a478e17a9c` |
| `sdss_frame_r_002583_2_0136` | M2 r frame (validation) | 3,524,824 | 13 Jun 2012 12:12:49 | `9129b6c6779ad8155741998e333e0f69395226aa9fc425dc9933f9f6f53ec04c` |
| `sdss_psfield_002583_2_0136` | M2 psField | 414,720 | 13 May 2009 05:43:44 | `3d24d270c12dcd8429357a81433d0151d3a6e8739401b18273d4bfc206bcd7bd` |
| `sdss_photofield_002583_2` | M2 photoField (camcol 2) | 650,880 | 01 Nov 2013 02:48:21 | `052132367a8faaee18160b71c3a0beef423259df7b9218c313d564f1d5f6e69e` |
| `hugs_ngc7089` | HUGS NGC 7089 (M2) method 1 | 106,805,615 | 29 Jun 2018 08:54:25 | `4a96ae437122468785200570dd9121f4b677ec6ed0200e65a687babcf5b58f49` |
| `sdss_frame_r_004646_3_0081` | M3 r frame (examined) | 3,133,162 | 20 Jun 2012 19:10:01 | `ff13054dcd440c1ca1d8f34bfeb6357f9197b20290ac3269a5b174a813c387ed` |
| `sdss_psfield_004646_3_0081` | M3 psField | 504,000 | 31 May 2009 21:34:34 | `4d96dcffeebee8b18ad078ce5b911370da29f5b1bc74bd15e1510b761653d840` |
| `sdss_frame_r_001458_4_0700` | M5 r frame (examined) | 3,675,609 | 11 Jun 2012 04:56:23 | `43138e7c06071489de77b5db66ed9517cd40f2b057162abec31e2dcec966f5ca` |
| `sdss_psfield_001458_4_0700` | M5 psField | 786,240 | 22 May 2009 09:01:36 | `7953cbfb97e510df4e34213b33cf7a7015f3844303fcf584005013814173bb4c` |

The other 234 photoField files (39 candidate runs x 6 camcols), the other HUGS catalogues of the survey (13 files, 2.7 to 105.7 MB) and the 52 name-lookup centres are recorded in the same JSON (`files`, `cluster_centres`, `hugs_probe`).

Frames (r band, rerun 301 of the DR17 SAS). Header and photoField values:

| | M2 (reference) | M53 | M92 |
|---|---|---|---|
| run, camcol, field; stripe, strip | 2583, 2, 136; 82 S | 5390, 6, 162; 17 N | 4682, 6, 208; 39 S |
| `DATE-OBS` (HUGS: epoch 2015.0) | 2001-09-20 (13.28 yr before) | 2005-06-06 (9.57 yr) | 2004-06-15 (10.55 yr) |
| `NMGY` (nanomaggies per count) | 0.00546689 | 0.00534669 | 0.00545525 |
| gain (e-/DN), dark variance | 4.600, 1.0000 | 4.895, 0.9025 | 4.895, 0.9025 |
| `SATURATION_LEVEL` r (DN) | 50,631 | 52,900 | 52,900 |
| photoField r: `PSF_WIDTH`, `SKY`, `AIRMASS`, `QUALITY` | 1.048, 5.218, 1.335, 1 | 1.317, 3.735, 1.110, 1 | 1.345, 3.025, 1.087, 1 |
| `PSP_STATUS`, `PSF_NSTAR`(r) | 2, 12 | 0, 14 | 0, 16 |

## 4. Astrometry and photometry, validated on M2

### 4.0 What "the shipped M2 truth" is

The brief calls `truth_hugs_meth1.npz` a shipped file. It is not: `data/astro/make_truth.py` builds it from the HUGS NGC 7089 catalogue with `astropy.wcs.WCS(frame header).all_world2pix(ra, dec, 0)`. A comparison with it checks my code path
and conventions, not the astrometry. The independent shipped M2 truth is `Data/sdss0921_true.txt` of the public repository (6,048 rows, ACS survey catalogue, patch pixels of the patch (310, 630), SDSS counts), which I use in
test [3]. The shipped files are read from `$ASTRO_SHIPPED` (here `/home/user/astro_data/fetchtest_4/shipped_inputs`; the repository URL is withheld from the repository by the no-author-names rule of `data_access.md`).

### 4.1 The two mappings

* Header: the frame's TAN solution (`CRVAL`, `CRPIX`, `CD`); pixels 0-based with the pixel centre at the integer (x = u + CRPIX1 - 1). Pure numpy (`tan_world2pix`), equal to astropy's.
* asTrans: HDU 3 of the frame, the SDSS data model (`PHOTO_REDUX/RERUN/RUN/astrom/asTrans`): rowm = row + dRow(col), colm = col + dCol(col) (cubics in the column), mu = a + b rowm + c colm, nu = d + e rowm + f colm, then
  (mu, nu) to (RA, Dec) for the great circle (node, inclination). The SDSS convention puts pixel centres at half-integers: row = y + 0.5, col = x + 0.5. This convention was found by checking the asTrans position of the header's
  reference pixel against `CRVAL`: 0.002 and 0.003 arcsec apart with the half-pixel, 0.20 arcsec without (M2 frame, scratch check; the data model says the header WCS is corrected to the asTrans solution). The colour term (`csRow`, `csCol`: r - i)
  is off in the truth: it moves positions by 0.0075 px (r - i = 0.3) and 0.025 px (r - i = 1.0) in x (printed in [2]).

### 4.2 Validation on M2 (`astro_second.py validate-m2`, M2 frame `2583-2-136`, 287,653 HUGS rows with a valid F606W)

| test | result |
|---|---|
| [1] header mapping against `truth_hugs_meth1.npz` (all rows; F606W < 20: 27,653 rows) | median 1.4e-10 px, 90th percentile 2.5e-10 px, maximum 3.8e-10 px; `m606`, `e606`, `nf606` identical (maximum magnitude difference 0): code path and conventions reproduced |
| [2] asTrans against header, benchmark patch box, F606W < 20 (233 rows) | median 0.0097 px, 90th percentile 0.0199 px, maximum 0.0224 px; signed median (asTrans - header) (+0.0064, -0.0075) px |
| [2] asTrans against header, whole HUGS footprint inside the frame (272,963 rows) | median 0.0555 px, 90th percentile 0.1404 px, maximum 0.1905 px (largest at small x) |
| [3] independent shipped truth, sources brighter than 20 (234): matched to a HUGS source within 1.5 px and 0.5 mag | 233 (99.6 %) |
| [3] |HUGS - shipped|, header mapping | **median 0.0519 px, 90th percentile 0.0597 px**, max 0.136 px; median (dx, dy) = (-0.0464, +0.0236) px, scatter 0.0040, 0.0054 px |
| [3] |HUGS - shipped|, asTrans mapping | **median 0.0435 px, 90th percentile 0.0570 px**, max 0.135 px; median (dx, dy) = (-0.0403, +0.0162) px, scatter 0.0136, 0.0046 px |
| [3] magnitude, HUGS F606W - shipped (flux at zero point 22.5, 0.00546689), F606W < 20 | median +0.0047 mag, scatter 0.0050 mag, 90th percentile of the absolute difference 0.0163 mag |
| [3] the same for all 1,069 shipped sources brighter than 22 | 1,065 matched (99.6 %); median 0.0517 px (header), 0.0419 px (asTrans); dm median +0.0030, scatter 0.0062 |
| [4] units: `extract_patch` on the M2 frame and psField against the patch `x310_y630` of `fetch.py` | `dn`, `electrons`, `nmgy`, `sky_dn`, `psf51`, `psf25`, `gain`, `dark_variance`, `x0`, `y0` identical bit for bit, `meta.json` identical; `write_patch` output: `patch_dn.fits` data identical, 36 header cards equal, one extra card (`NMGY`), npz arrays and dtypes identical |
| [4] PSF against `pydl.photoop.image.sdss_psf_recon` at the patch centre | maximum difference 1.06e-8 (unit-sum PSF) |
| [4] shipped raw counts (bias 1044 removed) against `dn` | correlation 0.99993, median difference -7.87 DN (the numbers of `m2_field_spec.md`: the shipped array is built from raw counts) |

The acceptance criterion of the brief (the new truth must reproduce the shipped M2 truth to better than 0.15 px in the median) is met by a factor 3 (0.052 px, header; 0.044 px, asTrans) against the independent list, and exactly
against the derived file. The position difference to the independent list is a small offset (dx = -0.04 to -0.05 px, dy = +0.02 px) with 0.004 to 0.014 px scatter. The scatter in dx is smaller for the header mapping (0.0040) than for asTrans
(0.0136): the shipped list varies with position like the header's linear solution, not like asTrans (the two mappings differ by +-0.02 px across the patch), so test [3] cannot tell the two mappings apart; it
shows that both reproduce an independent HST-to-SDSS mapping to 0.05 px.

### 4.3 Which mapping the new truth uses

The two mappings differ by at most 0.03 px in the M2 benchmark patch and in the M92 patches but by up to 0.14 px in the M53 patches (smaller x). The rigid-shift scan of section 4.4 tells which one the image prefers. Shift to
add to the truth positions (px) and chi^2 at zero shift, per patch (10,000 pixels):

| patch | asTrans - header, median |d| (signed median dx, dy) | image shift, asTrans mapping | chi^2 at zero shift | image shift, header mapping | chi^2 at zero shift |
|---|---|---|---|---|---|
| M2 (310, 630) | 0.0099 (+0.006, -0.008) | (+0.092, +0.032) | 73,419 | (+0.099, +0.025) | 76,603 |
| M53 (760, 590) | 0.0460 (-0.045, +0.008) | (-0.043, +0.015) | 43,640 | (-0.088, +0.023) | 50,527 |
| M53 (640, 680) | 0.0747 (-0.074, +0.010) | (-0.048, +0.043) | 110,556 | (-0.123, +0.052) | 138,957 |
| M53 (400, 660) | 0.1259 (-0.125, +0.017) | (-0.037, +0.062) | 205,087 | (-0.162, +0.080) | 307,402 |
| M53 (660, 420) | 0.0723 (-0.071, +0.012) | (-0.044, +0.030) | 95,012 | (-0.112, +0.041) | 123,630 |
| M92 (1100, 220) | 0.0133 (+0.013, +0.000) | (+0.134, -0.079) | 272,337 | (+0.147, -0.077) | 282,100 |
| M92 (1430, 550) | 0.0110 (-0.008, +0.008) | (+0.182, -0.063) | 354,271 | (+0.175, -0.049) | 345,995 |

On M53 the header mapping leaves a position-dependent offset (dx from -0.088 to -0.162 px) which the asTrans mapping removes (dx from -0.037 to -0.048 px, the same in the four patches) and the chi^2 at zero shift is lower in all
four patches; on M2 and M92 the mappings are indistinguishable (chi^2 differences of a few percent, one in each direction). So the truth lists use asTrans; the header-mapped lists are written next to them. This is a
departure from `make_truth.py`; it moves the M2-style truth by less than 0.02 px in M2's own patch and by 0.05 to 0.13 px (median) in the M53 patches.

### 4.4 The image-based check (`rigid_shift`), and what it does not do

The truth sources (F606W < 25, relative flux 10^(-0.4 m), positions of the 124 x 124 px area around the patch) are rendered with the psField PSF of the patch centre (bicubic interpolation of the 51 x 51 stamp), the
model is shifted rigidly on a grid of 0.05 px over +-0.4 px, scale and a sky plane are fitted by weighted least squares (variance DN / gain + dark variance), and the shift at the chi^2 minimum (parabola through the grid) is
reported. The reduced chi^2 at the best shift is 4 to 26 (the models are imperfect: F606W is not r, saturated stars have unreliable magnitudes); the shift comes from the shape of chi^2 against the shift. This is a measurement of the
net offset between the truth and the image (astrometry plus the proper motion between the epochs: one mas/yr over the epoch difference of the table in section 3 is 0.034 px for M2, 0.027 px for M92, 0.024 px for M53, arithmetic only,
no cluster proper motion was looked up); nothing is applied. On the benchmark patch it reproduces the image-based offset of `m2_field_spec.md` (+0.10, +0.02 px).
Checks that bright truth stars sit on bright pixels: in every new patch all scored F606W < 18 truth stars (6 to 29 per patch) have a maximum DN in the 3 x 3 pixels around them of at least twice the patch median (median ratio 6.7 to 19.8).

## 5. The patches

### 5.1 Counts (`build`; the truth box is [x0 - 0.5, x0 + 99.5)^2, `scored` = at least 2.5 px from the edges, the rule of `astro_common`)

| patch | F606W < 22 | per pixel | < 21 | < 23 | all valid F606W | scored, < 22 | brightest truth F606W | with the header mapping, < 22 |
|---|---|---|---|---|---|---|---|---|
| M53 (760, 590) | 402 | 0.0402 | 169 | 785 | 4,527 | 365 | 15.12 | 403 |
| M53 (640, 680) | 671 | 0.0671 | 288 | 1,180 | 6,427 | 594 | 15.39 | 671 |
| M53 (400, 660) | 950 | 0.0950 | 412 | 1,703 | 8,024 | 858 | 14.76 | 950 |
| M53 (660, 420) | 1,302 | 0.1302 | 574 | 2,317 | 10,004 | 1,163 | 15.25 | 1,303 |
| M92 (1100, 220) | 718 | 0.0718 | 423 | 1,115 | 3,141 | 638 | 14.75 | 718 |
| M92 (1430, 550) | 1,184 | 0.1184 | 741 | 1,767 | 4,460 | 1,052 | 14.51 | 1,184 |

Truth lists (arrays `px`, `py`, `m606`, `e606`, `nf606`, float64, frame pixels, every HUGS row with a valid F606W; same format as `truth_hugs_meth1.npz`; read them with `astro_common.truth_for_patch(truth, x0, y0)`):

| file | rows | bytes | sha256 |
|---|---|---|---|
| `truth_hugs_ngc5024_meth1.npz` (asTrans) | 256,593 of 259,364 | 10,264,936 | `7a31394d4f7523703f7085cb1315761b941e4f38fd1e77805411713ce374881e` |
| `truth_hugs_ngc5024_meth1_header.npz` (header TAN) | 256,593 | 10,264,936 | `61ad0b7742ed51bfbe04abfc89157c51611f50999bfd951f8443cb8cb71dce86` |
| `truth_hugs_ngc6341_meth1.npz` (asTrans) | 155,859 of 158,219 | 6,235,576 | `6ae0284af6b9b6c58af8fc422739621a3ff26ff2e2342b4365263e8156375f11` |
| `truth_hugs_ngc6341_meth1_header.npz` (header TAN) | 155,859 | 6,235,576 | `068cbdc2490a63b2b86314684c76fb9cff6435d7c597ef1f79d14e8ccdf731f8` |

### 5.2 Crowding, signal and quality, with the M2 patches by the same functions

PSF area = pi (FWHM / 2)^2; "per PSF area" = density x PSF area; blend = fraction of the F606W < 22 truth sources with another one closer than 2.22 px (the M2 FWHM of `m2_field_spec.md`) and than the patch's own FWHM;
peak S/N = peak of an isolated F606W = 22 source above the sky over sqrt((sky + peak) / gain + dark); `F606W(M2 S/N)` = the magnitude with the peak S/N of the M2 benchmark patch at 22.0 (6.33).

| patch | density | FWHM (px) | per PSF area | blend (2.22 px / own FWHM) | peak S/N at 22 | F606W (M2 S/N) | sky (DN) | DN max | DN median | 20 x 20 blocks, F606W < 24 (min / median / max) | F606W missing, F814W < 20 / < 22 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| M2 (110, 630), development | 0.1130 | 2.15 | 0.41 | 0.86 | 6.35 | | 149.2 | 31,125 | 339 | 15 / 101 / 160 | |
| M2 (210, 630), development | 0.1388 | 2.16 | 0.51 | 0.87 | 6.34 | | 149.2 | 46,860 | 380 | 57 / 106 / 194 | |
| M2 (310, 630), test, benchmark | 0.1070 | 2.16 | 0.39 | 0.81 | 6.33 | 22.00 | 149.2 | 26,861 | 310 | 61 / 90 / 161 | |
| M2 (310, 730), test | 0.2309 | 2.15 | 0.84 | 0.96 | 6.37 | | 149.1 | 35,390 | 566 | 109 / 177 / 356 | |
| M2 (360, 1040), test | 0.1027 | 2.13 | 0.37 | 0.79 | 6.47 | | 149.0 | 47,035 | 310 | 57 / 86 / 145 | |
| M2 (100, 1100), test | 0.1623 | 2.12 | 0.57 | 0.91 | 6.52 | | 149.0 | 38,789 | 446 | 79 / 137 / 208 | |
| M2 (410, 630), replaced | 0.0171 | 2.16 | 0.06 | 0.68 | 6.32 | | 149.2 | 29,724 | 234 | 0 / 0 / 92 | |
| M2 (410, 730), replaced | 0.0606 | 2.15 | 0.22 | 0.86 | 6.35 | | 149.1 | 24,895 | 299 | 0 / 28 / 148 | |
| M53 (760, 590) | 0.0402 | 3.36 | 0.36 | 0.53 / 0.79 | 4.22 | 21.51 | 109.5 | 14,890 | 156 | 28 / 47 / 82 | 0 / 2 |
| M53 (640, 680) | 0.0671 | 3.40 | 0.61 | 0.66 / 0.89 | 4.16 | 21.50 | 109.5 | 13,084 | 193 | 32 / 68 / 155 | 1 / 2 |
| M53 (400, 660) | 0.0950 | 3.38 | 0.85 | 0.75 / 0.94 | 4.18 | 21.50 | 109.5 | 21,934 | 235 | 50 / 104 / 183 | 1 / 2 |
| M53 (660, 420) | 0.1302 | 3.27 | 1.09 | 0.86 / 0.98 | 4.42 | 21.57 | 109.6 | 12,899 | 290 | 74 / 136 / 262 | 0 / 2 |
| M92 (1100, 220) | 0.0718 | 2.94 | 0.49 | 0.69 / 0.85 | 5.31 | 21.78 | 83.8 | 17,171 | 222 | 47 / 67 / 115 | 1 / 3 |
| M92 (1430, 550) | 0.1184 | 2.99 | 0.83 | 0.83 / 0.95 | 5.17 | 21.75 | 84.0 | 21,650 | 323 | 70 / 110 / 165 | 0 / 2 |

The M2 rows are the current protocol of `exp_astro_m2.py` (read on 2026-11-07 at 22:04): two development patches, four test patches (310, 630), (310, 730), (360, 1040), (100, 1100); (410, 630) and (410, 730), the first choice of test
patches, were replaced after a referee round found them mostly outside the HUGS footprint, which the block counts reproduce (22 and 13 of the 25 blocks hold fewer than 50 sources of F606W < 24, and (410, 630) has blocks with none). In
PSF units the M53 patches span 0.36 to 1.09 per PSF area (M2 test patches: 0.37 to 0.84; the densest M53 patch is above M2's densest), at a lower signal to noise and a broader PSF. "F606W missing" = sources with a measured F814W but
no F606W, absent from the truth.

Footprint check by blocks (the M2 protocol requires every 20 x 20 block of a patch to hold at least 50 sources of F606W < 24): every block of the six new patches holds sources, the smallest block counts are 28 to 74 (M53) and 47 and 70
(M92) against medians of 47 to 136 and 67 and 110, so no block is a hole (ratio smallest to median 0.47 to 0.70). The absolute threshold of 50 is calibrated on M2's density and does not apply to M53, whose F606W < 24 counts are about
half of M2's in its sparsest patches (13 of the 25 blocks of the patch (760, 590) and 5 of (640, 680) are below 50). An exploratory run (scratch script) over all 206 + 37 candidates of rule (2) found a ratio of the smallest block to the median
block of at least 0.34 for every one, and requiring 0.3 or 0.4 would not have changed the selection.

### 5.3 PSF, units, noise

* PSF: psField (r) reconstruction at the patch centre (x0 + 49.5, y0 + 49.5), 51 x 51, negative pixels clipped, unit sum (the recipe of `fetch.py`; against `pydl` at the same positions the maximum difference is 5.6e-9 to 6.7e-9 for
  the four M53 patches and 3.6e-9, 7.5e-9 for the two M92 patches). FWHM from the azimuthal mean about the centre (32 directions, cubic interpolation; a Gaussian of sigma 1.5 px gives 3.535 against 3.532): 2.12 to 2.16 px for the
  M2 patches (`m2_field_spec.md` quotes 2.22 px with another estimator), 3.27 to 3.40 px for M53, 2.94 and 2.99 px for M92; peak of the unit-sum PSF 0.1422 (M2), 0.074 to 0.079 (M53), 0.086 and 0.088 (M92).
* Units are those of M2: `dn` = calibrated frame / calibration vector + bilinearly interpolated sky (bias-subtracted, flat-fielded counts), `electrons` = `dn` x gain, noise variance `dn / gain + dark_variance`. Gain and dark variance
  are the psField values (4.895 and 0.9025 for camcol 6). The sky noise per pixel is sqrt(sky / gain + dark) = 5.8 DN (M2), 4.8 DN (M53), 4.3 DN (M92) (arithmetic from the table values).
* The flux to magnitude conversion of `astro_common` (`counts_to_mag`: 22.5 - 2.5 log10(flux x 0.00546689)) uses M2's `NMGY`; each frame has its own (section 3). With it, a source's magnitude from DN counts differs from the
  frame's own calibration by -0.0241 mag (M53) and -0.0023 mag (M92); `NMGY` is stored in the FITS header of each new patch (the only addition to the M2 format; the npz keys and `meta.json` are exactly the M2 ones).

### 5.4 Truth quality, for the transfer

* F606W is "used as r" (the M2 convention, no colour term). In exploratory fits (a scratch script, not delivered) the best-fit flux scale of the truth model differed from 1 by up to about 10 % in individual patches and the
  residual maps show star-by-star residuals at the bright end; I did not quantify the colour term.
* The HUGS F606W of HST-saturated stars (`nf606` = 99: 309 of the 155,859 M92 rows, 493 of the 256,593 M53 rows, 988 of the 287,653 M2 rows, mostly brighter than F606W = 17; scratch count) can be wrong by several magnitudes (section 2.4).
* Sources whose F606W is not measured but whose F814W is are not in the truth list (section 5.2: 2 to 3 per patch brighter than F814W = 22).
* No fpM mask, cosmic-ray or satellite-trail check was made on the new patches (as for M2).

### 5.5 Files (nothing here is in the repository; `$ASTRO_DATA_2`)

Each patch directory holds exactly the three files of the M2 patches. `patch_arrays.npz` is `np.savez` of `dn, electrons, nmgy, sky_dn` (float64, 100 x 100), `psf51`, `psf25`, `gain`, `dark_variance`, `x0`, `y0`; the npz bytes are
reproducible here (two builds gave identical files); the `content` hash is the sha256 of names, dtypes, shapes and bytes of the arrays, independent of the zip container.

| patch | `patch_arrays.npz` bytes, sha256 | content hash | `patch_dn.fits` bytes, sha256 | `meta.json` bytes, sha256 |
|---|---|---|---|---|
| `005390_6_0162_x760_y590` | 348,286 `335d00455ca3a201d8b4268d93736a50daa9e5390a521504ee48a9f836641cad` | `b9f800c4f5ab544604b9203af5b427488993d811c5b71fd26ec6b056e28f5eb6` | 46,080 `676d0b9a8dab703f1d9250a62168faed543febe316ee321a924c75827776e9df` | 270 `a9332b5e125c712a1f437b7a83fe482f5cff91318b3f513132898b30bdac6ac3` |
| `005390_6_0162_x640_y680` | 348,286 `f8cf1405e1b2c5d59aa67fe3a864649e7b2614a7a189ed6414508b6b41166382` | `ddd6c13287629173ae813e86cd21b300943864cd185409916a92b684ef3595de` | 46,080 `cb83530f02d5c7ee5a3f5c56f55783999e3574f75f4b21ac64ff032df689375d` | 270 `8f47a6a7d7705f0749dc492bacd9d0d22f40be7581018681cbbd3e42d0096a24` |
| `005390_6_0162_x400_y660` | 348,286 `5af5d8e03a077e30a5f3c2e88c4c259160a6008620885775dfb604aa8ad211ab` | `49c93591e279fef05db400247840405ca28365cdc7609d4a686469a2fbcd4d1f` | 46,080 `58b6de947413bdaeada8bf3aca7a616464694c4e31e35d7edefc6063b9f29e73` | 270 `80a4bc73790c9b5f6cd5c3e13fcd009c927a632e7ba4cc944e53bda752769aad` |
| `005390_6_0162_x660_y420` | 348,286 `3593e8e350d328f887014d8d4bb0dc360f471bcdadf40db0ce787495b782706a` | `543d89de69ef034e1ccd9c17355a6b66a249f01159bc0efc9c3ed411f8f98d20` | 46,080 `c933daeb08a5b2b11573428060e20e5aea7209474a77e0b4fb2ce87f40371185` | 270 `c6fe9fa2046b3b1ebb61c9ffe4ba52ff82bd79a5e32655d9b11530ef42008365` |
| `004682_6_0208_x1100_y220` | 348,286 `6357f69a72670311e93d69dc4a31360f2bba5ddaf2a045ffd43d9995d5a466be` | `ef25d334619a162f8753fa572e2cf3a2ed1b90e4e198e4439f6ac942336d40de` | 46,080 `9d576e3605a3f0e8fd944bd4fb7ffe1bf56066db1236c561729f95916a6075a2` | 271 `57d563c92a5a5208da448adca4decdb87048fb948281c6ca88bf4f97ba329030` |
| `004682_6_0208_x1430_y550` | 348,286 `c505cde1f2240685397bdcc563d1a0377fee29de29b95d74ca4de35bfadc93d2` | `3f1cbb4684b8b14f75a024ab1fdd2868573fc0501aefe031b7c6ada6b83ed836` | 46,080 `f620504fa9f8ff498cc11424e03f05aba99962259cd606ab5f5889086a04f6cc` | 271 `53f578d2ba6796148f29faf339d0a027ebb39328c5d6a72fce3c7feeceab1bcc` |

(Directory names are `patch_r_` + the first column + `_100x100`. The machine-readable table, with the per-patch numbers of sections 4.3, 5.1 and 5.2, is `$ASTRO_DATA_2/second_cluster_<cluster>_patches.json`.)

## 6. What differs from M2, for a transfer test

| | M2 | M53 | M92 |
|---|---|---|---|
| patches | 2 development, 4 test | 4 | 2 |
| density range, per pixel (F606W < 22) | 0.103 to 0.231 (test), 0.113, 0.139 (development) | 0.040 to 0.130 | 0.072 to 0.118 |
| PSF FWHM (px) | 2.12 to 2.16 | 3.27 to 3.40 | 2.94 to 2.99 |
| per PSF area | 0.37 to 0.84 (test) | 0.36 to 1.09 | 0.49 to 0.83 |
| sky noise per pixel (DN), gain | 5.8, 4.60 | 4.8, 4.895 | 4.3, 4.895 |
| peak S/N of a F606W = 22 source | 6.3 to 6.5 | 4.2 to 4.4 | 5.2 to 5.3 |
| epoch offset truth - image (yr) | 13.28 | 9.57 | 10.55 |
| measured net offset truth - image (px, asTrans mapping) | (+0.09, +0.03) | (-0.04 to -0.05, +0.02 to +0.06) | (+0.13 to +0.18, -0.06 to -0.08) |
| truth mapping | header TAN (`make_truth.py`) | asTrans | asTrans |

The configuration of the M2 method and the tuning of the baselines were not touched. Any parameter given in pixels (the exclusion radius, the high-pass widths) or tied to the noise (the 5-sigma rules) refers to the M2 PSF and
noise and will act differently on frames with a 36 to 57 % wider PSF; that is the transfer being tested, and the PSF is given to the method as in M2 (the psField reconstruction).

## 7. Licences and terms stated by the hosts (what I read; as in `data_access.md`)

| source | what is stated |
|---|---|
| SDSS SAS directory listings (dr17/eboss, resolve, photoObj, photo/redux, sweeps) and the data-model pages of `frame`, `asTrans`, `psField`, `photoField` | no licence, copyright, acknowledgement or citation text (zero matches in the eleven pages read: seven directory listings and four data-model pages); the SDSS citation and acknowledgement policy is on another site, not fetched |
| MAST Download endpoint, HUGS method-1 catalogues | no licence text; the file header asks users to cite Nardiello et al. (2018) and Piotto et al. 2015, AJ, 149, 91 ("When you use this catalog in any published work, please cite these papers"); the data are HLSP products of GO-10775 and GO-13297 |
| `fetch_second.py`, `astro_second.py`, the tests, this file | new files of this repository; third-party data are never stored in the repository |

## 8. Blocked, uncertain, not done

* Blocked: nothing. Only data.sdss.org and mast.stsci.edu were used (the name lookup is `https://mast.stsci.edu/api/v0/invoke`, service `Mast.Name.Lookup`; the `api/v0.1/resolve` and `.../Tools/name/lookup` paths answer 404). SkyServer, hla.stsci.edu,
  archive.stsci.edu were not contacted.
* Uncertain: (a) the meaning of `PSP_STATUS` values 2 and 3 (not defined on the data-model page I read; I use `PSF_NSTAR(r)` instead); (b) the magnitudes of HST-saturated stars; (c) the colour term between F606W and r;
  (d) the proper-motion share of the measured offsets of section 4.4 (no cluster proper motion was fetched); (e) whether M13 (0.69 of its footprint in the best frame) could host patches with relaxed rules: not tried; (f) the HUGS 404 names (section 2.1).
* Not done: a patch set of M5 or M3 (the rules leave none, one); a seeing-matched set (the frames with seeing near M2's are M5 and NGC 6934; the first fails the rules, the second covers 0.79 of its footprint with `PSF_NSTAR` 1);
  proper-motion or colour corrections; photometry in other bands; a check of the psField PSF against stars of the new frames.

## 9. Commands and tests

```
export ASTRO_DATA_2=/home/user/astro_data_2 ASTRO_DATA=/home/user/astro_data ASTRO_SHIPPED=/home/user/astro_data/fetchtest_4/shipped_inputs
DOMAIN=/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/venv_domain/bin/python        # numpy, scipy, astropy (pydl optional)
python3 -I data/astro/fetch_second.py fetch --group survey --group m2 --group hugs --group photofield --group target     # 260 files, sha256 verified or recorded
python3 -I data/astro/fetch_second.py resolve ; python3 -I data/astro/fetch_second.py probe                              # cluster centres, HUGS availability
python3 -I data/astro/fetch_second.py verify --group all                                                                  # 260 files OK
OMP_NUM_THREADS=1 nice -n 10 $DOMAIN -I experiments/astro_second.py survey                      # section 2.2 (35 s)
OMP_NUM_THREADS=1 nice -n 10 $DOMAIN -I experiments/astro_second.py validate-m2                 # section 4.2, 4.3 (M2 rows), 5.2 (M2 rows) (10 s)
OMP_NUM_THREADS=1 nice -n 10 $DOMAIN -I experiments/astro_second.py build --cluster ngc5024     # M53: truth lists, 4 patches, sections 2.3, 4.3, 5.1, 5.2
OMP_NUM_THREADS=1 nice -n 10 $DOMAIN -I experiments/astro_second.py build --cluster ngc6341     # M92: 2 patches
OMP_NUM_THREADS=1 nice -n 10 $DOMAIN -I experiments/astro_second.py build --cluster ngc5904 --dry-run     # M5 (0 patches); ngc5272: M3 (1)
PYTHONPATH=<domain site-packages> venv/bin/python -m pytest tests/test_astro_second.py -q       # astropy for the data tests
```

Everything is deterministic: a second `build` gave identical patch files, JSON and truth files. Tests (`tests/test_astro_second.py`, 24 tests: 23 functions, one parametrised over the two clusters): made-up TAN header
against astropy, asTrans round trip, convention and colour terms, agreement of asTrans with TAN for an equatorial scan, sky interpolation and DN conversion, the psField recipe, patch keys and slices, patch files, HUGS reader,
truth edge rules against `astro_common`, density windows, blend fraction, footprint, block counts, candidate rules and selection, run geometry, S/N, FWHM, rigid shift, `fetch_second` URLs and hashes, resumable download and checksum refusal against
a local server, and data tests that are skipped when the data are absent (the M2 validation numbers of section 4.2; the reported counts of the six patches of section 5.1 and their rebuild from the downloaded frames).
With the data and astropy: 24 passed (24 s). In the repository venv, without astropy: 19 passed, 5 skipped (the astropy comparison, the patch-file writer and the three data tests). A mutation check (12 deliberate bugs in a copy of the
module: CRPIX offset, closed truth box, integer pixel centres in asTrans, calibration multiplied, overlap rule removed, bright-star rule removed, untransposed PSF polynomial, swapped patch axes, wrong shift sign, wrong magnitude
limit, no footprint collar, swapped RA/Dec columns) was caught by a test in all 12 cases (scratch script, not delivered).
