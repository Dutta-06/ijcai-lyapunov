# Referee report, round BHP: the astronomy baselines given the ISS's front end and a wider pipeline search

Scope: `results/astro/baselines_hp_plan.md`, `experiments/exp_astro_baselines_hp.py`, `experiments/exp_astro_hp.py`, `results/astro/baseline_catalogs_hp.npz`, `baselines_hp_raw.npz`, `baselines_hp.md`, the two test files, the two mutation specifications (all committed at HEAD b3df74e). Nothing in the repository was modified (`git diff HEAD` is empty for these files and for every hash-frozen file). My work is in `/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_BHP_work/` (`myscore.py` is my scorer; `probe_cs.json`, `probe_phu.json` hold the development probes). Cost: about 35 CPU minutes at nice 10, one busy process at a time, except that two numpy scripts of under 10 s each ran while a single photutils probe was running in the background. No network.

## Verdict

No blocker. The chain of custody is consistent, the numbers reproduce (my own scorer matches every stored total and block of 112 catalogue scores; the report regenerates byte for byte; two catalogues re-run bit for bit), and the front end is consistent with the data operator. What does not hold as the coordinator intends to write it:

* "Photutils is not helped" is unsupported as a statement about photutils. The family omitted `xy_bounds`. On the development patches the selected variant with `xy_bounds=1` scores F 0.263 against 0.192 (+0.071), above the old tuned photutils (0.236). (Finding 1.)
* The crowdsource "wider pipeline search" was screened by catalogue size. On the development patches the number-of-passes options add +0.017, and +0.030 with a 2 px high-pass and threshold 3. The 44 % is a floor, not the share of a converged baseline. (Findings 2, 3.)
* "9 wins, 3 ties, 4 losses of 16" is not what the report says (8, 4, 4). (Finding 4.)
* The M53 and M92 sentences hold only as "ISS as frozen, pooled, one frame each". (Finding 5.)

## Findings

### 1. [major] The photutils family is handicapped by omitted options; "photutils is not helped" is not established

Location: plan lines 70-72 ("no local background, no xy bounds"), `experiments/exp_astro_baselines_hp.py:153` (`IterativePSFPhotometry(...)` without `xy_bounds`), `PHOT_HP` lines 55-63; plan section 9 (lines 146-157).

Evidence, development patches only (110,630) and (210,630), pooled F at m<22 with my scorer, same rule as the selection:

| variant | F | note |
|---|---|---|
| selected `phu_hp3_all_i6` | 0.1915 | c 0.119, p 0.482, 421 unmatched |
| same, `xy_bounds=1` | **0.2628** | **+0.071**; c 0.170, p 0.581, 302 unmatched, 150 CPU s for both patches (stored variant: 178 s) |
| `new` i30, relaxed cuts, `xy_bounds=1` | 0.2479 | +0.056 |
| `new` i30, relaxed cuts, `xy_bounds=2` | 0.2271 | +0.036 |
| `new` i30, default cuts, `xy_bounds=1` | 0.2113 | +0.020 |
| old tuned photutils (raw data, local background) | 0.236 | report line 30 |

My wrapper reproduces the stored catalogues exactly without the option (`all` i6 on (210,630): 432 sources, equal; `new` i10 on (110,630): equal), so the difference is the option. The stored catalogues show why: unbounded fits run away. `phu_hp3_all_i6` has 16-19 of about 400 sources per development patch outside the patch box (9 and 2 of them farther than 10 px), x from -5013.9 to 2233.5, flux up to 4.27e6 DN (the brightest star of the patch is about 3e5); `phu_hp5_all_i6` reaches 1.8e14 DN. The old tuned catalogues have the same defect (1.4e12 DN). Few such sources fall in the scored area (4 over the 16 test patches), so they do not move F directly; the harm is drift onto neighbours (completeness and purity).

Effect on the claim: the report's "photutils HP minus old tuned" (-0.045 [-0.057, -0.033] on M2) compares the old family with the best of a family that cannot bound its fits. Other probes did not rescue photutils: fit shape 5, 9, 11 (`new` i30): 0.1727, 0.1691, 0.1651 against 0.1729; 60 iterations: 0.1728; high-pass 2 or 4 px: 0.1646, 0.1724; relaxed cuts without bounds: 0.1977. A second pass on the residual of the stored `all` catalogue was attempted but is invalid for (110,630): the 4.27e6 DN source at (31.0, 1.4) corrupts the residual (std 8187 DN against a frame std of 767), so I discard those numbers.

Cheapest fix: write "photutils, in the declared configuration (unbounded fits, `all` mode limited to 6 iterations), is not helped (0.160 against 0.205 on M2)", and state that `xy_bounds=1` raises the development F of the selected variant from 0.19 to 0.26. A clean answer needs a new declared round with `xy_bounds` in the family (150 CPU s for the two development patches), selected on the development patches, then the 16 test patches. I did not score anything new on a test patch.

### 2. [major] Crowdsource knobs were screened by catalogue size; the selected variant is below its development plateau

Location: plan section 3 check 8 (lines 61-66, "Only the detection threshold changes the catalogue by more than a few per cent; it is the knob searched"), section 5 (lines 110-112), `CS_HP` (`exp_astro_baselines_hp.py:67-75`), `cs_fit` (line 177: `miniter`, `maxiter`, `blendthreshu`, `titer_thresh` left at the package defaults).

Evidence, development patches only (selected `cs_hp3_t5` = 0.4176, reproduced by my wrapper; sb = high-pass width in px, thr = threshold):

| option | F | change |
|---|---|---|
| miniter 8, maxiter 20 (fewstars 20, same value) | 0.4345 | +0.017 (1589 sources against 1564) |
| miniter 12, maxiter 30, fewstars 10 | 0.4341 | +0.017 |
| fewstars 20 alone | 0.4265 | +0.009 |
| sb 2 thr 3 | 0.4342 | +0.017 |
| **sb 2 thr 3 + miniter 8, maxiter 20, fewstars 20** | **0.4478** | **+0.030** |
| sb 2 thr 5 + same passes | 0.4435 | +0.026 |
| thr 4 / 6 / 7 at sb 3 | 0.4119 / 0.3908 / 0.3767 | -0.006 / -0.027 / -0.041 |
| sb 2.5 / 3.5 / 4 at thr 5 | 0.3979 / 0.3972 / 0.3987 | about -0.02 |
| blendthreshu 0.5 / 1 / 4 | 0.4123 / 0.4086 / 0.4030 | lower |
| titer_thresh 1 / 3 | 0.4002 / 0.4049 | lower |
| fewstars 200 | 0.3543 | -0.063 |
| psfderiv False | 0.2044 | -0.21 |
| high-passed stamp cut to 25 / 31 / 41 px | 0.4132 / 0.4218 / 0.4185 | flat |
| flux convention total flux instead of the central 19 x 19 | 0.4132 | -0.004 |
| old raw-count weights instead of the ISS weights | 0.3989 | -0.019 |
| significance cut flux/sigma >= 2 ... 6 on the catalogue | 0.4176 ... 0.4172 | 0.000 |

So the plan's inference from catalogue size to F fails: +1.6 % more sources moved F by +0.017. The development F of the declared 16 variants is also a best-of-many number: neighbouring settings fluctuate by 0.01-0.02 (sb 2.5 is below both sb 2 and sb 3), and the selected `cs_hp3_t5` beats `cs_hp3_t3` by 0.002 (0.4176 against 0.4156), a coin flip. The +0.030 is the maximum of about 45 probes and partly selection optimism; the +0.017 from more passes is consistent across settings.

Consequence: the "search" gain on M2 (+0.018 [0.006, 0.030], selected minus reference) is exactly the width change sigma_b 5 to 3, because the declared threshold search returned the package default 5. The ISS's own tuned width entered the baseline family. The family is a front-end and threshold family, not "a wider pipeline search", and the remaining 0.146 is an upper bound on what a tuned baseline leaves.

Cheapest fix: do not call it a pipeline search; say that the options not searched (passes, high-pass below 3 px) gain +0.02 to +0.03 on the development patches. Or add miniter, maxiter, fewstars and sb 2 to a declared round.

### 3. [major, wording] "The front end accounts for about 44 % of the gap": the number is right, the attribution is not

Evidence (my bootstrap, M2 test patches): gap ISS minus old tuned crowdsource 0.260; closed by the selected crowdsource 0.115, share 44.1 % [38.7, 49.7]. By the front end alone with the default pipeline (`cs_hp5_t5` minus the default) +0.097, share 37.2 % [32.3, 42.1]. On the six repeat patches the shares are 53.9 % [47.5, 59.9] and 41.4 %. On M53 the ISS-minus-old gap is only +0.031, so the share exceeds 100 %.

Three problems with "accounts for":

* It is one counterfactual (give the high-pass to crowdsource). The reverse counterfactual gives a different number: on the benchmark patch (310,630) the ISS falls from F 0.653 to 0.260 without the high-pass (`V|x310_y630|nohp`, below old crowdsource 0.377), a loss of 0.393, whereas crowdsource gains 0.121 (0.377 to 0.498). The ISS minus crowdsource sign flips with the high-pass (-0.117 without, +0.155 with). The ISS's machinery depends on the front end far more than crowdsource does, so no ordering-free "share" exists.
* "Front end" in the 44 % bundles the high-pass of image and stamp, the ISS noise weights (worth about 0.02 on the development patches: finding 2) and the choice of the ISS's own tuned width 3 px.
* At the same 5 px front end the ISS (a priori default, kappa 5) leads crowdsource (`cs_hp5_t5`) by +0.121 [0.106, 0.135] on M2 (0.468 against 0.347) and +0.062 [0.050, 0.074] on the repeat patches. The ISS-specific ingredients are worth about as much as at 3 px.

Cheapest fix, suggested sentence: "Giving crowdsource the ISS's front end (high-pass of image and stamp, the ISS weights) with a threshold and width search lifts its F on the four M2 test patches from 0.251 to 0.365, closing 0.114 (44 %, 39-50 %) of the ISS's 0.260 lead; at least that much, since options that were not searched add 0.02-0.03 on the development patches."

Where the remaining 0.146 sits (M2 test patches pooled, stored detections): completeness by truth magnitude, ISS / crowdsource-with-front-end / old crowdsource: 15-19 mag 0.90 / 0.86 / 0.82; 19-20 0.80 / 0.58 / 0.38; 20-21 0.52 / 0.29 / 0.13; 21-21.5 0.30 / 0.10 / 0.03; 21.5-22 0.16 / 0.04 / 0.01. It is depth in the 19.5-22 mag range (completeness +0.183, purity only +0.041), and it is a position-level effect, not a magnitude-matching artefact (position-only completeness at 20-21 mag: 0.70 against 0.43). It is not explained by what the plan lists as possibly unfair to the baselines: the ISS significance clean-up (adding the analogous flux/sigma cut to crowdsource changes dev F by at most 0.001, because crowdsource's own 5 sigma detection already enforces it); the threshold (threshold 3 gives -0.002); the border (finding in Verified). The ISS's PSF being a three-Gaussian fit while the baselines get the stamp itself favours the baselines. What remains is the joint sparse path, its time stopping rule and exclusion zone, none of which is ablated against a baseline here.

### 4. [minor] The per-patch tally in the commit message and the brief (9 wins, 3 ties, 4 losses) is not in the report

Location: commit message of b3df74e; the 16 rows "tuned ISS minus crowdsource selected" of `results/astro/baselines_hp.md`.

Parsing the 16 per-patch verdicts of the report gives **8 wins, 4 ties, 4 losses** (M53: loss, tie [-0.043 [-0.080, +0.001]], tie, loss; M92: two losses; M2 four wins; repeat four wins, two ties). My own bootstrap (different seed) gives 8 / 3 / 5 because the M53 patch (640,680) has an upper end of -0.001. Photutils: 14 / 1 / 1 in both. The M53 verdicts sit on the margin of the interval. The 9 and 3 are a hand count and wrong by one.

Cheapest fix: "8 of 16 patches (all four M2, four of six repeat), level on 4, behind on 4 (two M53, two M92)", or add a script-generated count to the report.

### 5. [minor] Scope of the M53, M92 and repeat statements

* M53: pooled ISS minus crowdsource-with-front-end -0.027 [-0.047, -0.008] (I reproduce -0.027 [-0.048, -0.008]). Per patch -0.063, -0.043, +0.021, -0.044; with the patch as unit the mean is -0.032, t-interval [-0.090, +0.026], which contains 0. The "ISS as frozen" is kappa 2 from M2, a poor operating point there (purity 0.27). The stored a priori default ISS (sigma_b 5, kappa 5) is level with crowdsource: 0.310 against 0.291, difference +0.019 [-0.003, +0.039]. The post hoc stamp-kernel ISS is ahead by +0.079 [0.059, 0.099].
* M92: two windows of one frame, -0.247 and -0.173; the default ISS also loses (-0.173 [-0.201, -0.147]); the post hoc kernel row loses by 0.046 [0.018, 0.077]. A two-window cluster has no patch-level interval at all.
* Repeat: +0.075 [0.064, 0.087] pools six patches that are four sky windows (windows (310,630) and (360,1040) appear twice) sharing truth with the M2 windows. Per patch the differences are +0.133, +0.123, -0.021, +0.013, +0.085, +0.084. By sky window (mean of the M2 patch and its repeat patches) the mean is +0.103, t-interval [+0.049, +0.157]; repeat patches only, by window, +0.078 [+0.024, +0.132]. The sign survives; the width of the report's interval understates the uncertainty by a factor of about two.
* "The ISS still leads on M2": supported (patch-level t-interval [+0.094, +0.187]; interior >10 px +0.133 [0.114, 0.153]).

Cheapest fix: wording: "pooled over four windows of one frame", "the ISS at its M2-frozen operating point", "intervals are conditional on the patches".

### 6. [minor] The chronology rests on file times and logs, with no external witness

* The plan was never committed or pushed before the scoring; the repeat experiment's `freeze_repeat.json` had that. The original 10:02 text is not preserved (the amendment overwrote the file, mtime 10:35:43); the surviving `hp_work/run_dev.sh` (10:07:44), which still lists `phu_hp3_all_i10 phu_hp5_all_i10`, corroborates amendment 1.
* The selection has no timestamp of its own (`save_part` rewrites the whole raw on every part; zip member times are all 1980-01-01). The stored plan hash is self-attested. The time of the freeze is only bracketed: last development catalogue 10:41:37; the two test-stage logs print the selected tags in their first line and end at 10:45:16 (2.4 CPU min) and 10:53:32 (9.7 CPU min), so the selection existed by about 10:43.
* `atime` is useless as evidence on reading `m2.md`, `repeat_raw.npz`, `repeat.md` (relatime: only the first read after a write is recorded; they show 09:36:12, 10:17:02, 10:17:12, the coordinator's repeat run).

I found no sign of access to ISS results before the selection (see Verified). Cheapest fix: for the next declared round, commit the plan and push before the first scoring and record the commit hash with the selection; for this one, say in the report that the witness is file times, logs and the stored hash.

### 7. [minor] Housekeeping

* `experiments/exp_astro_hp.py:401`: `'... 95 %% intervals:\n'` is not a format string, so "95 %%" prints literally in the 16 per-patch headings of `baselines_hp.md` (lines 69, 95, 121, ...). Use a single %.
* `tools/mutants/wp_astrohp_own.py` skips the id F11 (F10 then F12, lines 77-79). The author's result files (`scratchpad/mut_hp_rep/res_*.jsonl`) show F11 ("frozen_rows: the post hoc row for a patch the psf raw lacks") surviving in an earlier run ("25 passed"); it was dropped without a note, contrary to HANDOFF.md ("the docstring of each states ... what survives and why"). The survivor looks equivalent (the final filter in `frozen_rows` drops absent rows), but it should be recorded.
* `exp_astro_hp.py` was edited at 11:55:55, after the results (10:53:53). The report generations kept in `hp_work/` differ from the final only in the pipe-to-slash change of old tag names (50 lines), and the stored scores reproduce under my scorer, so I count it presentation-only.

## Verified

1. **Frozen files untouched.** `git diff HEAD` is empty for the 11 BHP files, `freeze.json`, `freeze_repeat.json`, `m2_raw.npz`, `psf_raw.npz`, `repeat_raw.npz`, `imaging.py`, `astro_common.py`. HEAD adds only new files. `exp_astro_repeat.check_freeze()` passes; `exp_astro_m2.check_freeze(raw)` returns the record. `git status` lists only other studies' untracked files and a modified `HANDOFF.md`.
2. **Hashes.** SHA-256 of the plan equals `SEL|plan_sha256` (944ad6f8...); `dev_hash` recomputed equals the stored one; `verify_selection` passes.
3. **Catalogues scored are the ones stored.** `baseline_catalogs_hp.npz` (768 arrays) equals the union of the author's six partial files, array for array; its development part equals `dev_merged_backup.npz`. My scorer (`myscore.py`, written from the `astro_common` docstrings, not importing `score`) applied to the stored catalogues gives totals and 20 x 20 blocks identical to the raw for all 80 test (5 rows x 16 patches) and 32 development scores.
4. **Selection.** Only `DEV|` keys enter `choose_hp`/`dev_F`; `part_test` and `tags_of_test_stage` call `verify_selection` first; ISS keys are read only in `load_frozen`/`frozen_rows`/`patch_rows` (called by `report`) and `part_old` reads names from `freeze.json`. My scorer on the dev catalogues selects the same variants with the same rule: photutils `phu_hp3_all_i6` (0.1915; next 0.1729), crowdsource `cs_hp3_t5` (0.4176; next 0.4156). The tie rule is identical to `exp_astro_m2.choose_baselines`. The author's exploration scripts (09:46-10:00) contain no truth path. First report (10:54:29) differs from the final only by the pipe escape.
5. **Numbers.** With my own stratified block bootstrap (4000 resamples, seed 4711): M2 ISS 0.511 [0.488, 0.535], selected crowdsource 0.365, difference +0.146 [0.127, 0.164]; photutils HP 0.160, difference +0.351; front end alone +0.097 [0.081, 0.113]; search +0.018 [0.006, 0.031]. M53 -0.027, post hoc +0.079; M92 -0.201, post hoc -0.046; repeat +0.075 [0.063, 0.088] and +0.085 [0.071, 0.098]; m2r2 -0.005 [-0.023, +0.014] tie. All agree with the report to the third decimal; the stored ISS detections re-scored by me give the report's pooled F for the ISS (0.511, 0.263, 0.122, 0.324); all 955 frozen block arrays sum to their stored totals.
6. **Report regenerates byte for byte**: `HP_OUT=<my dir> exp_astro_hp.py --from-raw`, `cmp` identical.
7. **Re-runs bit for bit** (domain interpreter): `cs_hp3_t5` on both development patches (752 and 812 sources), `cs_default` on (110,630), `phu_hp3_all_i6` on test patch (310,630) (405 sources; 142.9 s CPU against 143.6 s stored).
8. **Front end.** Planted stars with the real psField stamp of (110,630), rendered by a Fourier shift (independent of the fitters' splines), noise-free, 6 sub-pixel draws per distance, the author's runners: interior (12 and 25 px) photutils flux +0.3 to +0.6 % (max 0.9 %), crowdsource +0.1 to +0.3 % (after the central-19 flux convention), positions within 0.005 px, no ghost. At 3 px from the border photutils -0.8 % (3 px high-pass) / -2.0 % (5 px) and +0.033 / +0.020 px along the border axis, crowdsource +0.4 % / -0.9 % and +0.037 / +0.025 px; at 6 px at most 1.3 %, 0.012 px. This agrees with the plan's -1.1 % / -1.7 % at 3 px. On real catalogues, F with truth and detections restricted to more than 10 px from the border (M2 / M53 / M92 / repeat): ISS 0.507 / 0.258 / 0.116 / 0.323 (all 0.511 / 0.263 / 0.122 / 0.324); crowdsource `cs_hp3_t5` 0.373 / 0.304 / 0.328 / 0.255 (all 0.365 / 0.291 / 0.323 / 0.248); `cs_default` 0.257 / 0.239 / 0.214 / 0.165. The ISS lead over `cs_hp3_t5` becomes +0.133 [0.114, 0.153] / -0.046 [-0.069, -0.022] / -0.212 / +0.067 [0.054, 0.081]. Border-band to interior F ratio on M2: ISS 1.03, `cs_hp3_t5` 0.93, `cs_default` 0.93. The model's interior stamp operator therefore costs the baselines at most 0.01-0.02 near the edge and favours no conclusion: no verdict changes, and the ISS's exact border operator is worth about 0.013 of the 0.146 on M2.
9. **Plan claims.** DAOStarFinder threshold 3 and 5 give bit-identical catalogues (new, 10 it, 3 px, (110,630)); crowdsource `nskyx=nskyy=3` raises `ValueError` as stated; every `psf51` sums to 1.0000 and the central-19 convention shifts crowdsource fluxes by 1-3 % (0.01-0.03 mag), with no effect on F.
10. **Tests.** 54 passed, 9 skipped in the repository venv (26 s); the 9 skipped (domain) tests pass with the domain interpreter (8 s): 63 of 63. The tests are strong on the selection rule (tie after rounding, pooling, hash refusals) and on the front end.
11. **Mutation specifications.** The author's result files show 108 of 108 (`exp_astro_hp.py`) and 110 of 110 (`exp_astro_baselines_hp.py`) mutants killed by failed assertions, none by import error, control passing. Coverage of the decisive lines is good: tie and rounding (S01-S03), eligibility (S05), freeze (V01-V10), test routing (T01-T11), weights (E01-E06), high-pass and stamp (E07-E18, P01-P04, C14-C19), threshold (C04, C12, C13), flux convention (C16, C17), field swap (P26, C01). Not covered: options that are absent from the code (finding 1 and 2), and the bootstrap itself (hash-frozen `exp_astro_m2.py`, covered by earlier specifications). I did not re-run the mutants.

## Could not check

* Whether the author opened `m2.md`, `repeat_raw.npz`, `repeat.md` or the ISS keys before the selection: only the code paths, the logs and the exploration scripts can be read, and they show no access; a negative cannot be proved (finding 6).
* The exact moment of the selection (bracket 10:41:37 to about 10:43).
* Whether the development probes of findings 1 and 2 carry over to the test patches (not allowed to score them); the +0.030 and +0.071 are optimistic by selection, the first more than the second.
* Photutils `all` mode beyond 6 iterations: the bounded `all` i10 run used 7.4 CPU minutes on the first development patch and was stopped by process id. The grouper radius and the fit shape in `all` mode were not probed either.
* An explicit second crowdsource pass on the residual image (the passes options of finding 2 stand in for it), and DAOStarFinder options other than threshold and cuts.
* The mutation suites were not re-run; the HUGS truth, the patch extraction and the frozen ISS code are outside this round.
* `exp_astro_hp.py` as it was at the time of scoring (edited later); only its outputs, which reproduce, could be compared.
