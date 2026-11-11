# Analysis plan for the held-out repeat-frame test (declared before the sliding inverse scale space is run on any repeat patch)

Status at the time of writing (2026-11-08). The six repeat patches of M2 exist and are validated (`repeat_spec.md`; frames `m2r1`, `m2r2`, `m2r3`: other SDSS runs of the same four test windows,
psField FWHM 2.28, 3.05 and 3.52 px against 2.12 to 2.16 px for the M2 frame of the frozen test). The catalogues of photutils and crowdsource on them were written once at the settings of
`freeze.json` (`baseline_catalogs_repeat.npz`) and were never scored against the truth. **The sliding ISS has not been run on any of the six patches, and no catalogue of any method has been scored
against their truth.** This file, the script `experiments/exp_astro_repeat.py` and the freeze record `freeze_repeat.json` (hashes of this file, of the code and of the inputs) are committed
before the run; the run refuses to start unless the record matches, and it never overwrites a result.

## 1. What the test is for

Two claims of the project rest on the astronomy results and were made after the data were seen:

* (T) the configuration chosen on two development patches of M2 (sigma_b 3 px, r_ex 1.2 px, kappa 2; the a priori default 5, 1.2, 5; the BIC stop) beats tuned crowdsource and tuned photutils on four
  test windows of the same frame (frozen test, `m2.md`); on M53 and M92, imaged with worse seeing, it did not (small gain on M53, loss on M92);
* (K) the failure is explained by the kernel: the isotropic three-Gaussian mixture fitted to the psField stamp (`imaging.fit_psf_parts`) misfits the elliptical, lumpy stamps of M53 and M92, and a
  kernel of separable Gaussians with own widths and offsets chosen from the stamp alone by `psf_general.fit_adequate` (rms goal 2e-3 of the peak; all six M2 stamps keep the isotropic fit, bit for bit)
  repairs a large part of it (`psf.md`). K is POST HOC: the family, the rule and its goal were fixed after M53 and M92 had been scored, on the same patches.

The repeat patches are images, PSFs and noise realisations that neither the configuration, nor the baselines' tuning, nor the kernel rule has seen. They are NOT independent sky: the same four
windows, the same HUGS stars as the M2 test patches (the truth is shared, with its defects: HST-saturated stars, stars without an F606W). Six patches, three frames: the unit of replication for
seeing is the frame.

## 2. Methods, all frozen

* The configuration, the a priori default and the tuned baselines are read from `freeze.json` (sigma_b 3, r_ex 1.2, kappa 2; default 5, 1.2, 5; photutils `photutils_ext_all_local_i6`, crowdsource
  `crowdsource|snr5|ff0|qf0|spinf`) and checked against it; nothing is re-tuned on any repeat patch. The configuration is in pixels, as it was on M53 and M92.
* ISS with the isotropic kernel: `exp_astro_m2.run_path` unchanged. ISS with the adequate kernel: `exp_astro_psf.run_path_psf` with `psf_general.fit_adequate` of the patch's stamp. Where the rule keeps
  the isotropic fit the two are the same program on the same input; the run is done once and used for both rows (the table says so).
* Both kernels are run at the tuned configuration (kappa 2), at the a priori default (kappa 5) and at the a priori default with the BIC stop (key 0), exactly the three readouts of the frozen test.
* Baselines: default and tuned, from `baseline_catalogs_repeat.npz`.
* Scoring: `astro_common.score` (match at 0.5 px and 0.5 mag, one to one; truth within 2.5 px of the edge not scored), completeness, purity and F at mag < 22; 95 % percentile intervals of a block
  bootstrap (20 x 20 px blocks, 2000 resamples, seed 0, stratified within each patch, the same resamples for every row of a table, so differences are paired). The offset between truth and image
  (up to 0.26 px in y for m2r2, `repeat_spec.md` 6.2) is NOT corrected for any method.
* Further baselines may be added to the tables after the main run, under three conditions: they are chosen from a declared family on the two M2 development patches only, their catalogues are written to a file of their own before they are scored, and the report labels them as added later; they cannot change P1 to P3, which concern the ISS (a family of baselines given the ISS's high-pass front end is planned: referee round R3, finding F7).
* Exploratory, labelled as such: the configuration with its lengths scaled by the seeing, sigma_b and r_ex multiplied by s = FWHM(stamp of the patch) / FWHM_ref, FWHM_ref = the mean FWHM of the stamps of the
  two M2 development patches (`astro_second.psf_fwhm`), kappa unchanged, adequate kernel, tuned configuration and a priori default. It is a question ("do the pixel-unit constants need to follow the seeing?"),
  not a result to be reported as a method: the row is in the supplement, whatever it shows, with this sentence.

## 3. Outcomes declared before the run

Primary (the wording of the paper follows the rule, whatever the result):

* **P1 (the kernel rule on held-out frames).** On the patches where the rule picks a general kernel, the pooled paired difference dF = F(adequate) - F(isotropic) at the tuned configuration (kappa 2).
  Positive 95 % interval: "confirmed on held-out frames"; interval containing 0: "not confirmed"; negative: "contradicted". The same for the a priori default and the BIC stop (secondary). On patches where the
  rule keeps the isotropic kernel dF is 0 by construction and these patches do not enter P1.
* **P2 (competitiveness at the frozen configuration).** dF of the adequate-kernel ISS (kappa 2) against tuned crowdsource and against tuned photutils, per patch (win / tie / loss by the paired interval) and
  pooled, for all six patches and for each frame. The same for the a priori default and the BIC stop against the baselines' defaults. The record is stated as it is: if the ISS loses on a frame, the paper says so.
* **P3 (the diagnosis, label-free).** The robust standard deviation of the standardised residual of the catalogue (1.4826 MAD), isotropic against adequate kernel, per patch. If the isotropic kernel's value is
  not above the adequate kernel's on the patches where the rule picks a general kernel, the diagnosis (K) is not supported on these frames.

Reported beside F in every table (referee round IS1b-A: F at mag < 22 does not see faint false positives): the purity over all magnitudes (matched detections over all detections) and the number of unmatched detections.

Descriptive (no claim rests on them): F of each method against the seeing of the frame, with the M2 results of the same windows (frozen, `m2_raw.npz`) beside the repeat results; the entrants, the BIC block,
the run time; the stars brighter than 19.5 mag matched; the stamp fits of the repeat patches (rms of the isotropic fit, kind chosen, ellipticity).

## 4. What would change the project's statements

* If P1 is not confirmed, the statement in the paper that the kernel explains the transfer failure is reduced to "post hoc on M53 and M92, not confirmed on repeat frames of M2".
* If P2 shows losses to a baseline on a frame, the paper's astronomy claim is restricted to the seeing range in which the ISS wins (M2's), and the repeat frames are reported as a limitation.
* Nothing in the method, the kernel rule or the configuration is changed after the run. A bug found afterwards is corrected only with the first output kept (`repeat_raw.npz` is write-once; a rerun
  needs the first file removed in a commit that says why), and the report states both.

## 5. Known limits (declared now)

* The windows are the M2 test windows: shared truth and shared sky with the frozen M2 test; two of the four windows per frame, not the same two; the selection rule found no window passing all five rules
  (a bright HST-saturated star lies in or within 10 px of every window) and was refined twice before building (`repeat_spec.md` 1, 2, 2b).
* Seeing 2.28 px (m2r1) is a same-seeing control, not a step; the ladder has two broader frames (3.05, 3.52 px).
* Epochs: the truth is HUGS 2015.0, the frames 6 to 13 years earlier; no proper motion or colour term. Magnitudes from counts use M2's fixed zero point (0.03 to 0.08 mag offset between frames).
* m2r2 is rotated by 1.7 degrees against M2 and has the largest astrometric offset (+0.24 to +0.26 px in y).
