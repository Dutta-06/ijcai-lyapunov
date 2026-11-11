# M2 (NGC 7089) crowded-field patch: field specification

Scope: data access and setup only. Nothing here tunes, runs or claims anything about a cataloguing method.
Facts come from (a) the text of two arXiv papers, quoted, and (b) commands run on 2026-11-07 against the hosts listed in
`data_access.md`. "Not confirmed" marks what I could not establish. The papers are cited by arXiv number only:

| tag | paper | version read | sha256 of the PDF read |
|---|---|---|---|
| P17 | arXiv:1703.01303, probabilistic cataloguing of an SDSS r-band image of M2 | v2 | `0ad4554a2c9f27f190d069fb3e01228195493b99b3501d4ecc5507e84f796f97` |
| MB | arXiv:1907.04929, multiband follow-up (same cutout, plus i and g) | v2 | `b4cb52b9760f1e73cd08de4592fc626b68b97f87c90db3f6172b12fcca039541` |

I also read the code and input files that the follow-up points to (public git repositories, see `data_access.md`
section 5). They contain the exact cutout origin, the metric code and the input arrays, and they settle several points
the papers leave open. Where I rely on them I say "shipped code" or "shipped file".

## 1. The target object

* M2 is "a rich and compact globular cluster containing about 150,000 stars located 11.5 kpc away", core radius 0.34 arcmin
  (1.1 pc), half-light radius 1.08 arcmin (3.6 pc) (P17, quoting the 2010 edition of the standard cluster catalogue).
* Centre used here: the MAST name resolver (SIMBAD) returns RA 323.36258 deg, Dec -0.82325 deg
  (21h33m27.02s, -00d49m23.7s, J2000). The brief's value (about 21h33m27s, -00d49m24s) is confirmed.
* The 2017 paper's goal: a catalogue ensemble of a crowded image ("~0.1 sources per pixel brighter than 22nd magnitude
  in F606W"), judged against a deeper HST catalogue of the same region.

## 2. SDSS field and patch

| item | what the papers say | what I verified |
|---|---|---|
| run / camcol / field / band | P17: "SDSS run 2583, field 136, camcol 2, in the r band". MB: "run 2583, camcol 2, field 136" | frame header: `RUN=2583, CAMCOL=2, FRAME=136, FILTER=r`, `STRIPE=82` |
| rerun | not stated in either paper | DR17 serves this frame under rerun 301. The frame files shipped with the public code are byte-identical to the rerun-301 files (cmp, r, g and i). Other reruns not checked |
| no SDSS catalogue | P17: "There is no SDSS catalog for this field, because ... Photo ... timed out" | `photoObj-002583-2-0136.fits` is 8,640 bytes (empty); `photoField-002583-2.fits` row for field 136: `NOBJECTS=0, NSTARS=0`. Neighbour 135 is 7,490,880 bytes; 137 is also 8,640 bytes. The r-band frame (3,524,824 bytes) and the psField exist anyway |
| patch size | "100 x 100 pixel cutout"; MB: "roughly 0.43 arcmin^2" | (100 x 0.396 arcsec)^2 = 0.436 arcmin^2 |
| patch origin | MB: "pixel coordinate (x0,y0) = (310, 630) in the frame file defining the lower left corner" | x0 = column, y0 = row, 0-based: patch = `frame[630:730, 310:410]` (frame is 2048 columns x 1489 rows, 0.396 arcsec/pixel). See "orientation test" below |
| centre of patch | MB: "R.A. = 21h 33m 20.5s, decl. = -0d 48' 26.5\"" | frame-header WCS gives 21h33m20.57s, -00d48m25.49s for the centre pixel (col 359.5, row 679.5). Differs from the quoted value by 1.42 arcsec = 3.6 pixel. **Not resolved**, see section 8 |
| distance from cluster centre | P17: "approximately 2 arcminutes from the center" | 1.88 arcmin (patch centre to name-resolver centre, header WCS) |
| cluster core in the frame | not stated | name-resolver centre maps to (col 212.5, row 923.9). A sigma = 15 px Gaussian-smoothed copy of the frame peaks at (col 211, row 925): 0.7 arcsec from the prediction. The patch sits at pixel offset (+147, -244) from the core |
| sky orientation | not stated | header CD matrix: Dec increases with column, RA increases with row (north = +column, east = +row) |
| exposure | P17: "the 34.9 s SDSS exposure" | frame header `EXPTIME = 53.907456` s, `DATE-OBS = 2001-09-20`. The paper's 34.9 s is probably a typo; not resolved |
| saturation | P17: "away from the heavily saturated core" | patch DN max 26,861 (frame max 51,135, 323 frame pixels above 3e4, none of the patch pixels) |

Orientation test (settles the origin convention): the 100 x 100 counts array shipped with the public code, after
subtracting its bias (1044), correlates with `frame[630:730, 310:410]` at 0.99993 with no rotation and no flip. The other
seven orientations give correlations between -0.03 and +0.10. Median pixel difference (shipped minus mine) is -7.9 DN,
mean -10.5 DN, rms 12.5 DN: the shipped array is built from raw counts, mine from the calibrated frame, so a percent-level
flat-field difference is expected.

Patch statistics from the calibrated frame (DN = nanomaggy image / calibration vector + interpolated sky; run
`python -I data/astro/fetch.py extract-patch`): min 168.996, median 309.978, mean 442.763, max 26,861.107 DN.
Shipped counts minus 1044: min 166, median 302, mean 432.22, max 27,006.

## 3. HST catalogue

What the papers say:

* P17: "the ACS Globular Cluster survey cataloged M2 using imaging from the Advanced Camera for Surveys on HST ... five
  340 s HST exposures vs the 34.9 s SDSS exposure". The survey is the ACS Survey of Galactic Globular Clusters (HST programme
  GO-10775); both papers cite a 2007 survey paper for "the HST catalog".
* P17: "The HST catalog identifies 1,049 sources brighter than 22nd magnitude in the F606W band in this patch."
  F606W "is much broader than SDSS r, but centered at roughly the same wavelength, and F606W magnitudes for main sequence stars
  agree with SDSS r magnitudes to a few tenths of a magnitude."
* MB: "HST detects 6051 sources down to F606W ~ 30, 1428 of which have F606W < 23"; truth treated as reliable "down to
  F606W ~ 23"; "The region we use has been observed by HST in the ACS F606W and F814W bands. We compare with F606W."
* P17 and MB give no catalogue file name, no magnitude system (Vega or AB) and no description of how HST positions were mapped
  to SDSS pixels. The follow-up only says an SDSS astrometric transformation file (asTrans) was used for band-to-band mapping.
  The shipped metric code reads a file named `hubble_catalog_2583-2-0136_astrans.txt` (x, y in frame pixels, F606W, F814W), which suggests
  an asTrans mapping for the HST positions too; that file is not in the public repositories.

What exists and what I can reach (details in `data_access.md`):

| catalogue | status | F606W<22 in the patch | F606W<23 | all valid F606W |
|---|---|---|---|---|
| ACS survey photometric catalogue itself | **not reachable**: MAST only holds a pointer record that leads to hla.stsci.edu (refused by the egress policy); archive.stsci.edu fails TLS | n/a | n/a | n/a |
| HUGS catalogue, method 1 (MAST, same GO-10775 images plus GO-13297) | downloaded, sha256 pinned | 1,070 | 1,625 | 5,394 |
| HUGS catalogue, method 2 | downloaded, sha256 pinned | 1,071 | 1,623 | 5,407 |
| truth file shipped with the public code (6,048 rows, x y flux in SDSS counts) | downloaded by git, sha256 pinned | 1,069 | 1,611 | 6,044 (mag<30) |
| paper | | 1,049 | 1,428 (MB) | 6,051 (MB) |

Counting boxes: HUGS rows use the continuous pixel box [309.5, 409.5) x [629.5, 729.5) after mapping RA/Dec with the
frame-header WCS. Shipped-truth variants (F606W<22): all rows 1,069; 0 <= x,y < 100: 1,059; 0.5 <= x,y < 99.5: 1,051;
1 < x,y < 99 (the selection in the shipped helper `get_hubble`): 1,032; 2.5 < x,y < 96.5 (the margin in the shipped metric
script): 949. For F606W<23 the 2.5-pixel-margin selection gives 1,429 (follow-up: 1,428). With that metric-script selection (F606W and F814W both valid, more than 2.5 pixels inside the edges) HUGS method 1 gives 951 brighter
than 22 and 1,445 brighter than 23 (method 2: 952 and 1,444), the shipped truth 949 and 1,429, the follow-up 1,428 brighter than 23.
So 1,049 is reproduced to within 0.2% by the half-pixel box (1,051) but **I could not find the selection that gives exactly 1,049**.

Magnitude conversion (shipped code, not in the papers): `mag = 22.5 - 2.5 log10(flux_counts x nmgy_per_count)` with
`nmgy_per_count = 0.00546689` (the frame header `NMGY`). The F606W magnitude is used directly as an r magnitude, no colour term.
With this conversion the shipped truth spans 14.999 to 30.502 mag. The shipped file lacks F814W.

Agreement between the two independent HST-derived lists (HUGS method 1 mapped with the frame-header WCS, versus the shipped
truth): for the 1,069 shipped sources brighter than 22, 99.6% have a HUGS neighbour within 1.0 pixel and 99.3% within 0.3
pixel; median offset (HUGS minus shipped) dx = -0.046, dy = +0.023 pixel, rms 0.024 / 0.036 pixel; median magnitude difference
0.003 (scatter 0.20). Counts brighter than 21 are identical (591 and 591). They differ at the faint end (HUGS 5,394 sources
in the box versus 6,048 shipped; in the magnitude bins [26,27) and [27,28) HUGS has 353 and 21, the shipped file 852 and 358),
so the shipped truth is deeper than HUGS.
The number of F606W exposures per patch source in HUGS is 1 to 5 (99 flags saturated stars), consistent with five exposures;
the 340 s per exposure matches the exposure time on the MAST record of the F606W HLSP image.

Image-based check of the HST-to-SDSS alignment (a measurement, not tuning): a PSF-convolved model of the HUGS F606W sources
(relative flux 10^(-0.4 m), free scale and constant sky) fitted to the patch by scanning a rigid shift gives the best offset
dx = +0.10, dy = +0.02 pixel (0.04 arcsec); chi^2 falls from 97,335 (zero shift) to 59,282 (reduced chi^2 5.9, imperfect model).
So HST positions mapped with the frame-header WCS already agree with the image to about a tenth of a pixel, which is well
inside the 0.75 / 0.5 pixel match radii.

## 4. How completeness and false discovery rate were computed

From P17 (section 3.3-3.4):

* Match: a catalogue source matches an HST source if it is "within 0.75 pixels of the HST position (cp. the full width at half
  maximum of the SDSS point spread function, 2.2 pixels) and its r magnitude is within 0.5 mag of the 606W magnitude". The
  criterion is called loose because F606W is wider than r and because of faint blended HST neighbours. MB tightens the
  position criterion to 0.5 pixel (same 0.5 mag) and estimates spurious associations at about 3% at r = 22 and below 1% for r < 20.
* Completeness: fraction of HST sources matched, plotted against HST F606W magnitude. For an ensemble of M samples a match in
  one sample counts 1/M (average over samples). For the condensed catalogue a source with prevalence p counts p.
* FDR = FP / (FP + TP): a true positive is a catalogue source with an HST match, a false positive one without; plotted against
  SDSS r magnitude; per sample 1/M weights as above.
* No bin edges in the text; the figures run from 15 to 24 mag. **Shipped metric script**: 17 bins of 0.5 mag between 15.0 and 23.5
  (`minr, maxr = 15.0, 23.5`, `nbins = 17`), `dr = 0.5` pixel, `dmag = 0.5` mag (the 0.75 pixel value is the labelling radius for
  the condensed catalogue), a kd-tree ball query followed by the magnitude test; HST sources kept only if both F606W and F814W
  are positive and they lie more than 2.5 pixels inside the patch edges (`hwhm = 2.5`). FDR is 1 - (matched / total) per catalogue
  magnitude bin; completeness is matched / total per HST magnitude bin.
* Reported reference results (claims of the papers, not checked by me): the DAOPHOT catalogue (356 sources in the region per MB)
  starts losing completeness at 18th mag and is incomplete by 22nd; its FDR peaks at 47% at 18th mag and stays above 10%
  fainter; the ensemble is "over 1 magnitude deeper at every level of completeness" and has a lower FDR than DAOPHOT brighter
  than 20th mag; MB reports FDR between 20% and 30% for the 2017 ensemble between mag 18 and 20.

## 5. SDSS PSF, sky and gain inputs

| input | paper statement | shipped files / DR17 files (verified) |
|---|---|---|
| PSF | P17: "the pixel-convolved point spread function extracted by the standard SDSS pipeline for the center of our image"; assumed constant over the patch and "taken as known". MB: PSF from psField files, upsampled with a Lanczos kernel | Shipped PSF (25 x 25 stamp, 5x oversampled) equals the psField KL reconstruction at 0-based (col 359.5, row 679.5), unit sum over 25 x 25: correlation 1.000000 at the 25 x 25 sampling points, max residual 1.5e-6 after normalisation |
| PSF width | P17: FWHM 2.2 pixels | reconstructed 51 x 51 PSF: FWHM 2.22 px (radial-profile half maximum), peak 0.1422 of unit flux; psField r: `psf_sigma1 = 0.8802, psf_sigma2 = 2.3084 px, psf_b = 0.0588`, `psf_width = 1.0476` (arcsec per the data model) |
| sky | P17: uniform sky `I_sky` is a free parameter, log-uniform prior. MB: background fixed in the main runs; mock values `(b_r, b_i, b_g) = (179, 310, 115)` ADU "based off empirical background estimates from SDSS run 2583, camcol 4 [sic], field 136" | shipped analysis code lists 179 / 310 / 95 for r / i / g (g differs from the paper's 115); psField r sky 150.75 DN, sky sigma 6.42 DN; interpolated ALLSKY mean over the patch 149.18 DN; lowest patch pixel 169.0 DN |
| gain | P17: not stated (data are "photoelectron counts", Gaussian noise with variance equal to the expected count). MB: "bias and gain ... obtained from the opECalib ... file", mock "gain set to 4.62 across all bands" | shipped r pixel file: bias 1044, gain 4.62 (g: 1113, 4.2; i: 1177, 4.39). psField r gain 4.6, dark variance 1.0; the SDSS data model table for camcol 2 r also gives 4.6 and 1.00. 4.62 vs 4.6 is a 0.4% difference |
| counts | "photoelectron counts k_lm" | shipped code reads the raw counts file, subtracts the bias (1044) and keeps ADU units, with variance `= (counts - bias) / gain` and inverse-variance weights from the observed image, i.e. the Poisson variance in DN is DN / gain (the papers' "photoelectrons" differs by the factor gain). The calibrated frame has no bias; `fetch.py extract-patch` rebuilds DN as `frame / cal + sky` (SDSS data model: electrons = DN x gain, error in DN = sqrt(DN / gain + dark variance)) |
| noise model | log-likelihood `-(k - lambda)^2 / (2 lambda)` with the `-0.5 log(lambda)` term dropped | same in MB |

## 6. Model and run settings stated in P17 (context only, not run)

Flux power-law index 2 between `f_min` and a free `f_max`; `f_min` set at the SDSS 95% completeness limit of a sparse field,
"approximately 22 mag in r band" (shipped code: 236 ADU, which is 22.22 mag); exponential prior `pi(N) ~ exp(-3N/2)`; at most 3,000
sources; four chains for the convergence statistic; labelling radius 0.75 pixel; flux degradation-factor cut of 8; about 1,100
sources inferred; sampler of that paper: diffusive nested sampling; cost: "a day of wall clock time" on 12 cores, 1 core-ms per sample
(MB restates it as about 1e9 model evaluations over 12 CPU-days, versus half an hour on one CPU for the later implementation).

## 7. Files needed to rebuild the benchmark from reachable hosts

| need | file | where |
|---|---|---|
| r-band image | `frame-r-002583-2-0136.fits.bz2` (rerun 301) | data.sdss.org |
| PSF, gain, sky | `psField-002583-2-0136.fit` | data.sdss.org |
| HST reference | HUGS catalogue (method 1 or 2), or the shipped truth file | mast.stsci.edu / public git |
| exact original inputs | four small shipped files (counts, pixel parameters, PSF, truth) | public git, `fetch.py shipped` |

## 8. Doubts and unconfirmed items

1. Patch centre: the quoted centre (21h33m20.5s, -00d48m26.5s) is 1.42 arcsec (3.6 pixel) from the header-WCS centre of the
   box (310, 630). The box itself is confirmed by the pixel comparison above, so the quoted centre probably came from a different
   astrometric solution; the HST positions mapped with the frame-header WCS agree with the image to 0.1 pixel.
2. The number 1,049 is not reproduced exactly (1,051 to 1,070 depending on the box); the selection is unstated.
3. SDSS exposure 34.9 s (paper) versus 53.907456 s (header).
4. Camcol "4" in the MB mock-background sentence, and the g-band background 115 (paper) versus 95 (shipped code).
5. Gain is not given in P17; 4.62 comes from MB and the shipped files, the psField and data model say 4.6.
6. Magnitude system of F606W (Vega or AB) is not stated in either paper or in the HUGS file header; the shipped code simply
   maps F606W magnitudes to SDSS counts with the r nanomaggy zero point.
7. The original ACS-survey photometric catalogue that the papers used cannot be downloaded from the reachable hosts; HUGS and the
   shipped truth file are stand-ins that agree on positions to 0.05 pixel and on the count brighter than 21 mag, but are not the same table.
8. Magnitude-bin edges are not in P17; the 15.0-23.5 / 0.5 mag binning comes from shipped code written for MB.
9. The DAOPHOT comparison catalogue exists only as shipped files. The 356-row `truth/...-tru.txt` file in the multiband repository
   is it in patch pixels (x, y, then fluxes in counts): its third column equals the r-band flux computed from the r magnitude of the
   shipped photometry table (ratio 1.000), and 355 of its 356 rows have a table row within 1 pixel when the table's RA/Dec are mapped
   with the frame WCS (median distance 0.87 pixel, so the table's sky positions sit about 0.9 pixel from the frame-WCS positions).
   The meaning of columns 4 and 5 (other-band fluxes) was not verified.
10. DAOPHOT baseline positions carry an offset. Relative to the shipped HST truth, the 356-row file is offset by (HST minus DAOPHOT)
    dx = -0.04, dy = -1.06 pixel (peak of the pair-offset histogram, 165 of the 241 DAOPHOT sources brighter than 20 mag sit in the peak).
    As shipped, only 11.6% (0.5 pixel) and 25.3% (0.75 pixel) of those 241 sources have an HST source brighter than 22 within the radius;
    after removing the shift, 87.1% and 94.2%. Using the table's RA/Dec through the frame WCS instead gives a different offset (-0.63, -0.41).
    The table's own y column is about 1361.7 above the frame row (the shipped metric script subtracts 1992 = 630 + 1362; the shipped preparation script uses 630 on the same column). How the papers handled
    this is not established: the DAOPHOT evaluation lines of the shipped metric script are commented out. Any DAOPHOT baseline rebuilt from
    these files needs an explicit, documented position convention.
