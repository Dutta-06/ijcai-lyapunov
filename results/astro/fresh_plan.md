# Analysis plan for the fresh-cluster test (declared before any catalogue of the baselines was built, before the ISS was run, before anything was scored on a fresh patch)

**Declaration time: 2026-11-08T22:08:48Z** (`date -u`, taken immediately before this file was written). Author: the agent that carries out the coordinator's brief for a fresh cluster-level test. This file, the rule file `fresh_rule.md`, the
scripts `experiments/astro_fresh.py`, `experiments/exp_astro_fresh_baselines.py`, `experiments/exp_astro_fresh.py` and the freeze record `freeze_fresh.json` (hashes of this file, of the code and of the inputs) are meant to be committed
before the ISS run; the run refuses to start unless the record matches, and it never overwrites a result. The run script's functions are written AFTER this file and tested only on planted worlds (never on a fresh patch);
the freeze record is taken after the baseline catalogues exist and before the run, so a bug found in the script before the freeze is corrected before it, and the report says so if it happens.

## 0. State at the time of writing (what exists, what does not)

* **Exists:** the two clusters, their frames and the three patches (`fresh_rule.md` with its amendments A1 to A3, `fresh_survey.txt`, `fresh_fields.json`, `fresh_selection.txt`, `fresh_selection_A1.txt`, `fresh_patches.json`),
  the truth lists (`truth_hugs_<cluster>_meth1.npz`, HUGS method 1, the full asTrans mapping), the label-free stamp fits (below) and the descriptive rigid offset of the truth mapping against each patch (`fresh_patches.json`, `shifts`).
* **Does not exist:** any catalogue of any method on a fresh patch (the baselines' catalogues are built by `exp_astro_fresh_baselines.py` after this file), any ISS run, any score, any comparison of a detection with a truth source of
  a fresh cluster. The builder of the patches read the frame's pixel maximum (saturation rule) and fitted a PSF-convolved copy of the truth list to the pixels (rigid offset, descriptive); no detector output enters either.
* **What I know about the clusters from the earlier work of the project** (stated, because it can bias): nothing about M3 or M15 beyond the survey tables and the placement dry runs of `second_cluster_spec.md`. The frozen
  results of the other clusters are known to me (`baselines_hp2.md`: the ISS's lead over the front-end crowdsource is +0.118 on M2, -0.058 on M53 and -0.215 on M92 for the frozen isotropic kernel, +0.049 and -0.060 with the stamp kernel,
  +0.054 on the repeat patches).
* **Label-free facts known now:** the stamp of M3's patch is fitted to rms 1.02e-3 of its peak by the isotropic mixture, so the kernel rule keeps the isotropic kernel there (the isotropic and the rule's rows are the same program on the same input);
  the stamps of M15's patches have isotropic rms 2.83e-3 (above the goal 2e-3), so the rule picks the general family (3 parts, rms 0.96e-3, ellipticity 0.055). Stamp FWHM: 2.81 px (M3), 3.18 px (M15).

## 1. What the test is for

The astronomy evidence of the paper is three clusters, one SDSS frame each (M2: the configuration is tuned on two of its patches and wins on its four test patches and, with the repeat frames, at other seeing; M53 and M92: transfer, where the
frozen isotropic kernel loses to crowdsource, and the stamp-chosen kernel was found after the fact). Every choice of the project (sigma_b 3 px, r_ex 1.2 px, kappa 2; the a priori default; the kernel rule; the front-end baselines
selected on the M2 development patches) was made with those three clusters on the table. **This test asks one question: on clusters that no choice has seen, does the frozen pipeline keep its lead over the strongest baselines?**

The **unit of replication is the cluster**. There are **n = 2 clusters** (M3 with ONE patch, M15 with TWO patches of one frame). Patches of a cluster are windows of one frame: not independent, and the intervals below are conditional on them.
With n = 2 a result can contradict a general statement ("the ISS beats crowdsource on any cluster") by one loss; it cannot establish one by two wins. The report says so beside every pooled number.

## 2. Data (all fixed)

| cluster | frame (SDSS rerun 301, r) | stamp FWHM (frame centre / patch) | patches (x0, y0): F606W < 22 per pixel | HUGS | rigid offset of the truth (descriptive) |
|---|---|---|---|---|---|
| M3 (NGC 5272) | run 4646, camcol 3, field 81 (2004-05-22) | 2.81 / 2.81 px | (1870, 820): 0.128 | method 1, 72 MB | (+0.14, -0.03) px |
| M15 (NGC 7078) | run 2507, camcol 6, field 39 (2001-08-23) | 3.18 / 3.18 px | (640, 1180): 0.078; (920, 1180): 0.134 | method 1, 106 MB | (+0.36, -0.06) and (+0.38, -0.12) px |

How the clusters were chosen: `fresh_rule.md` (rule written before it was executed, executed at 22:00 UTC, three amendments A1 to A3 declared as each failure appeared, all before any catalogue existed; each amendment is stated there with what I knew).
In short: of the 11 HUGS clusters left in the SDSS coverage survey after M2, M53 and M92, three have a frame with a measured PSF (PSP_STATUS <= 2, at least one PSF star), a stamp FWHM in 2.3 to 3.6 px and at least 100 survey
windows of moderate density: M10, M3 and M15. M10 leaves no patch in any of its two frames (stars brighter than F606W = 14.5 within 10 px of every admissible origin), M3 one, M15 two (in its third passing frame; the first two leave none).
**Consequence: two clusters, three patches, one cluster with a single patch; the nearby dense clusters that the patch rules remove, and the far sparse ones that the density rule removes, are NOT in this test (a limit, section 6).**
Patch rules: those of M53 and M92 unchanged (`astro_second.py`; 100 x 100 px; 40 px from the frame edge; inside the HUGS footprint; no pixel above 0.8 x the photoField saturation within 10 px; no truth source brighter than F606W = 14.5 and
no HUGS source with F814W < 17 and no F606W within 10 px; density 0.03 to 0.15 of F606W < 22 per pixel; one patch per density class, no overlap) plus the patch-centre stamp FWHM in [2.3, 3.6] px. Truth: HUGS F606W used as r, zero point 22.5
and the M2 frame's 0.00546689 nanomaggies per count for every frame (the frames' own values differ by 0.03 mag: M3 0.00530650, M15 0.00526951), positions by the full asTrans solution, no proper motion, no colour term.

## 3. Methods, all frozen

* **Configuration.** From `freeze.json` (the M2 development patches): sigma_b 3 px, r_ex 1.2 px, kappa 2 ("tuned"); the a priori default (5, 1.2, 5); the BIC readout (key 0 of a path: the block of the smallest BIC along it; reported for the a priori default path, "default with the BIC stop", and stored for the tuned path). Nothing is re-tuned; the lengths are in pixels (not scaled with the seeing).
* **ISS rows.** (i) *isotropic kernel*: `exp_astro_m2.run_path` unchanged (the frozen M2-test program); (ii) *rule kernel*: `exp_astro_psf.run_path_psf` with `psf_general.fit_adequate` of the patch's stamp (rms goal 2e-3 of the peak, at most 10 parts); where the rule
  keeps the isotropic fit (M3) the two are one program on one input, run once and used for both rows (the table says so). Each at the tuned configuration (kappa 2), the a priori default (kappa 5) and the default with the BIC stop.
* **Guarded path.** EVERY ISS run calls `imaging_guard.sliding_path_guarded` (the two-strike guard against the enter/exit limit cycle of `imaging.sliding_path`, `cyclefix.md`) through a thin wrapper that replaces `imaging.sliding_path`
  in the worker process, as `exp_astro_cyclefix.py` does; the unguarded path is never used. For each run the raw stores the cap flag (entry cap 2500), the entry count, the guard activations (barred lattice points with their entry time)
  and the exit/re-entry counts; a run that reaches the cap is flagged "capped" in every table and is never pooled silently (it is pooled with its flag, stated in the table caption).
* **Baselines (no re-selection).** From the catalogues written by `exp_astro_fresh_baselines.py` with the functions of the earlier scripts, imported not copied: photutils default (`photutils`) and tuned (`photutils_ext_all_local_i6`);
  crowdsource default (unfiltered) and tuned (post-filter `crowdsource|snr5|ff0|qf0|spinf`) - the two tunings frozen in `freeze.json`; **and the front-end baselines exactly as selected by round 2 of the baselines given the ISS's front end
  (`baselines_hp2.md`, selected on the two M2 development patches): photutils `phu2_hp3_all_i6_b1`, crowdsource `cs2_hp2_t3_pm`** (read from the verified selection `exp_astro_hp2.verify_selection2`). These are the baselines
  "crowdsource tuned" and "crowdsource + high-pass" of the primary outcomes; photutils likewise.
* **Scoring.** `astro_common.score` (match at 0.5 px and 0.5 mag, one to one; truth within 2.5 px of the edge not scored), completeness, purity and F at mag < 22; 95 % percentile intervals of a block bootstrap (20 x 20 px blocks, 2000 resamples,
  seed 0, stratified within each patch, the same resamples for every row of a table, so differences are paired). The offset between truth and image is NOT corrected in the primary scoring (as for M53, M92 and the repeat patches).
  **Secondary scoring ("registration-corrected truth"):** the truth positions of a patch shifted by the rigid offset in `fresh_patches.json` (fitted to the image and the truth before any catalogue existed), the same catalogues, the same
  bootstrap; it never decides a verdict.

## 4. Outcomes declared before the run

**Primary** (per cluster; a verdict is *win* if the 95 % paired interval of dF is above 0, *loss* if below, *tie* if it contains 0; dF = F of the first minus F of the second, patches of the cluster pooled with the stratified bootstrap):

* **O1** (the claim as written in the paper): dF(ISS tuned (kappa 2), rule kernel; crowdsource tuned).
* **O2** (the strict test): dF(ISS tuned (kappa 2), rule kernel; crowdsource with the ISS's front end, `cs2_hp2_t3_pm`).

The ISS row of O1 and O2 is the **rule-kernel row** (identical to the isotropic row on M3). The same two differences for the **isotropic-kernel row** (the program frozen on M2; O1i, O2i) are reported beside them; if a verdict differs between the
two rows the report says so and the paper's statement is the rule-kernel row's with the isotropic one named.

What each result means (the wording of the paper follows this, whatever the result):

| O2 over the 2 clusters | reading |
|---|---|
| win in both | the lead over the strongest baseline we have holds at two clusters never seen by any choice, at stamp FWHM 2.8 to 3.2 px; stated with n = 2 and with M3 a single window |
| win in one, tie or loss in the other | the lead is not general: it is cluster-dependent; the paper gives the per-cluster verdicts next to the earlier ones (wins at M2 and its repeat frames, loss at M92) and no general sentence |
| tie or loss in both | the paper's claim is restricted to the M2 field and seeing; against baselines that have the ISS's front end the ISS has no demonstrated advantage outside it |
| (added to any row) a loss at a cluster | that cluster is reported as a loss in the statement of limits, whatever the other cluster shows |

O1 is read the same way against the old tuned crowdsource; a win at O1 together with a tie/loss at O2 means "the lead over the baseline as originally tuned survives, the lead over the baseline given the front end does not".

**Secondary** (reported, no claim rests on them): O1/O2 against photutils (tuned, and the round-2 front-end photutils); the a priori default and the BIC readouts against the baselines' defaults and the front-end baselines; **S1** the pooled
paired dF = F(rule kernel) - F(isotropic) at kappa 2, a priori default and BIC, on the patches where the rule picks a general kernel (M15's two) with the wording of the repeat plan (positive interval "confirmed on a held-out cluster", containing 0
"not confirmed", negative "contradicted"; M3's patch enters with dF = 0 by construction and is left out); **S2** the registration-corrected truth versions of O1, O2; **S3** both clusters pooled (a descriptive row: M3 enters with one patch, M15 with two, so the pooled rows
are weighted by patches; no verdict rests on it).

**Descriptive:** the stamp fits (isotropic rms, kernel chosen, parts, ellipticity), the residual scale of the catalogues (1.4826 MAD of the standardised residual, isotropic and rule kernel), completeness by magnitude bin of the main rows, purity over all
magnitudes and the unmatched detections beside every F, entrants, BIC block, cap flags, guard activations, run times.

## 5. What would change the project's statements

* O2 in the reading of the table above; nothing in the method, the configuration, the kernel rule or the baselines changes after the run. A bug found afterwards is corrected only with the first output kept (`fresh_raw.npz` is write-once; a rerun needs the
  first file removed in a commit that says why), and the report states both.
* If a path hits the entry cap or the guard acts, the report says so beside the affected rows; no row is dropped for it.
* If the verdicts of the primary scoring and the registration-corrected scoring differ, both are stated; the primary decides.

## 6. Known limits (declared now)

* **n = 2 clusters, 3 patches; M3 has ONE patch** (25 blocks of one 100 x 100 px window). The cluster-level intervals are conditional on the windows. This is a smaller test than the three-cluster record it adds to.
* **Selection of the clusters.** Only M3 and M15 satisfied the written rule at all; the rule needed three amendments (the single-patch cluster; frames other than the best-footprint frame; the next passing frame when one leaves no patch), each declared before
  the next step and with what was known (`fresh_rule.md`). The patch rules (bright stars within 10 px, HUGS footprint, density) remove the nearby dense clusters (M10, M5, M13 and most of M15) and the density rule the far sparse ones (NGC 5053, 5466):
  **the test is of the clusters that these rules let through, not of an unbiased sample of globular clusters**. The universe is the 11 clusters of the coverage survey (the 38 HUGS clusters with no SDSS run and the 13 whose file name the
  endpoint did not find are outside it; the latter are southern or not in HUGS, `second_cluster_spec.md` 2.1).
* **Frame choice by windows.** M15's frame is the third passing frame of the cluster (the first two leave no patch). Its footprint fraction is 0.42 (the patches lie in the part of the HUGS field it covers).
* **Offsets.** The net offset of the truth against the image is 0.14 px (M3) and 0.36 to 0.38 px (M15) in x; the match tolerance is 0.5 px, so the primary scoring penalises every method for it, the ISS and the baselines alike; S2 removes it. The epochs differ by
  10.6 (M3) and 13.4 (M15) years; the proper motion of the cluster is part of the offset, the scatter of the stars' own proper motions is not corrected.
* **Truth defects of HUGS** (stars saturated in HST, stars without an F606W magnitude, blends below the HST resolution) are shared by every method; magnitudes from counts use M2's zero point for every frame (0.03 mag).
* **The ISS's constants are not scaled with the seeing** (the exploratory scaled run of the repeat test is not repeated); the front-end baselines are the round-2 winners of two M2 development patches (a family of 8 + 6 variants, not exhaustive;
  the options not in the family are listed in `baselines_hp2_plan.md` section 7).
* **Guard limitation:** the guard removes repeated exits at one entry time; a cycle that advances tau is not touched (`imaging_guard.py`). The astronomy paths use `order='dual'`.
* The windows of M15 are two windows of one frame: they share the PSF model, the sky and the calibration. Time since the truth: see above.

## 7. Run protocol

    python -I experiments/exp_astro_fresh_baselines.py            # the baselines' catalogues, once (domain environment; no truth is read)
    python experiments/exp_astro_fresh.py --freeze                # freeze_fresh.json: hashes of this file, the code, the registry, the patches, the truth, the catalogues
    nice -n 10 python experiments/exp_astro_fresh.py --run        # the ISS on the three patches, once each, guarded path; refuses unless the freeze record matches; write-once raw results/astro/fresh_raw.npz
    python experiments/exp_astro_fresh.py --from-raw              # regenerates results/astro/fresh.md byte for byte from the raw

One heavy process at a time, `nice -n 10`, `OMP_NUM_THREADS=1`. The raw is `results/astro/fresh_raw.npz` (write-once); the report is regenerated from the raw alone (the raw stores the SHA-256 of this plan and of the freeze record, the cluster titles, the
offsets and every score; `--from-raw` reads nothing else).

## 8. History of this file

* Declared at the time on the first line. One edit before the freeze record was written and before anything was run (the sentence closing section 7, which had named other inputs of the report than the raw); nothing else changed since the declaration.
