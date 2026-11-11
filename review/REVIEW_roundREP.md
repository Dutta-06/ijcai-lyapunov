# Referee report, round REP: the held-out repeat-frame test of the sliding inverse scale space on SDSS images of M2

Reviewer: senior, adversarial, fair. Reviewed HEAD 851e5a7 (branch claude/lyapunov-context-setup-yg51ew). Nothing in the repository was modified: bytecode writing and the pytest cache were off, and `git status` shows only other agents' untracked files. My scripts and logs are in `/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_REP_work/`. The main ones are `myscore.py` (own scorer), `rescore_all.py`, `pooled.py`, `p1_fragility.py`, `shift_boot.py`, `rerun_iss.py`, `cap_check2.py`, `data_check.py`, `star_check.py`, `baseline_rerun.py`.

## Summary

No blocker. The chain of custody holds, the numbers reproduce, and the ISS reproduces bit for bit.

- **Chronology.** Plan, script and freeze were committed and pushed before the raw result existed. The freeze record matches the code and inputs, both at HEAD and at f05ca8d.
- **Numbers.** My own scorer reproduces 60 stored rows exactly, including the 25-block arrays. A re-run of the frozen ISS on m2r2 (133, 281) reproduces the stored catalogue bit for bit. The baselines re-run bit for bit on all six patches. The report regenerates byte for byte.
- **P2 holds.** The ISS beats the tuned baselines on F at m < 22 on all six patches. This survives alignment of the truth, different match tolerances and depth cuts down to m < 21, with one counter-example at the bright end (finding 4).
- **P1 "confirmed on held-out frames" does not stand as worded.** It is one frame and two patches. Its sign depends on the 0.25 px truth offset that the plan left uncorrected, and one of its secondary verdicts is a Monte Carlo accident (finding 1).
- **The "entry cap" runs are a defect, not benign censoring.** They are an enter/exit limit cycle of a single lattice point in `imaging.sliding_path`. It invalidates 3 of 12 a priori default runs here. The same cap flag also appears on 2 of 2 earlier M92 default runs and on 6 of 12 development-grid paths that `choose` ranked as if valid (finding 2).
- **The inference to K is weaker than the plan assumes** (finding 3).

## Findings

**1. major. P1 "confirmed on held-out frames" is one frame and two patches, and its sign depends on the uncorrected truth offset.**
Location: `results/astro/repeat.md:49-54` and `:39-41`; `experiments/exp_astro_repeat.py:284-297`; `results/astro/repeat_plan.md:33,43-45,58`.
Evidence (my scorer, which matches the report's counts exactly):
- (a) **One frame, one patch carries the gain.** The rule picks a general kernel on m2r2 only. Its isotropic stamp rms is 4.4e-3 and 5.6e-3. All other held-out stamps are 1.35e-3 to 1.87e-3, and the six M2 stamps are 1.0e-3. Per patch, dF(tuned) is -0.004 [-0.044, +0.029] on (133, 281) and +0.076 [+0.032, +0.118] on (195, 689).
  - On the sigma_b = 3 path, (133, 281) is negative at every readout (kappa 2 to 8 and BIC: -0.004 to -0.027). (195, 689) is positive at every one (+0.036 to +0.076).
  - The pooled [0.006, 0.062] resamples blocks inside two fixed patches. With two patches of one frame there is no between-patch or between-frame variance (patch effects: mean 0.036, sd 0.057).
  - Within that conditional frame it is stable: block size 10/20/25/50 px gives lower ends +0.012/+0.005/+0.007/+0.005, and 200 seeds give +0.003 to +0.008.
- (b) **The frame where the rule acts has the largest truth offset in the project.** The m2r2 offset is (-0.09 to -0.01, +0.24 to +0.26) px (`repeat_spec.md:355-378`). The others are at most 0.09 (M2), 0.14 (m2r1), 0.04 (m2r3), 0.18 (M53) and 0.06 (M92).
  - I re-measured it with template fits of 6 bright isolated stars using the patches' own PSF stamps: median (-0.02, +0.23) px.
  - The two kernels' unit positions differ systematically on m2r2. The median adq minus iso over flux-consistent unit pairs is (+0.06, -0.09) and (+0.06, -0.11) px. Against the unshifted truth the isotropic run sits at median y offset +0.21/+0.27 px and the general-kernel run at +0.12/+0.15 px, so the 0.5 px match radius costs the isotropic run more.
  - Re-scoring with the truth moved by alpha times the rigid-scan offsets of the spec gives pooled dF = +0.035 (alpha 0), +0.012 (0.5), +0.001 (0.75) and -0.015 [-0.043, +0.012] (1.0; per patch -0.042 and +0.013).
  - With offsets from my own star fits (-0.018, +0.227) the result is -0.007 [-0.036, +0.021]. With only +0.15 px in y it is +0.005 [-0.023, +0.033]. With the truth moved and per-method median magnitude offsets removed it is -0.002.
  - The cause of the 0.1 px kernel difference is not identified. A noise-free star built from the stamp is recovered by both kernels to within 0.013 px, so it is not a kernel-origin convention.
- (c) **The a priori default readout is not confirmed.** The report prints "0.028 [0.000, 0.056], confirmed". The seed-0 lower end is +0.00001. Over 200 other seeds (B = 2000) it averages -0.0014 and is positive in only 7 %. B = 40000 gives [-0.0013, +0.0573]. By the plan's own rule this is "not confirmed". Both secondary readouts also compare a complete general-kernel path on (133, 281) (1545 entries) with an isotropic path that hit the cap (finding 2).
- (d) **The all-six row is diluted.** The "All 6 patches" row (`repeat.md:39`, win 0.009 [0.001, 0.016]) pools four patches whose dF is 0 by construction. It narrows the interval and carries the same offset sensitivity.

Cheapest fix:
- Report P1 as "positive under the declared scoring on one held-out frame (patches -0.004 and +0.076); not robust to the 0.25 px truth offset (-0.015 [-0.043, +0.012] after alignment)".
- Add the alpha table to `repeat.md` as a labelled sensitivity analysis beside P1.
- State the two secondary readouts as not confirmed, or print B = 40000 intervals.
- Drop or label the all-six dF row. Say "one held-out frame", not "held-out frames".

**2. major. The "entry cap" runs are an enter/exit limit cycle of one lattice point in `sliding_path`, not censoring.**
Location: `src/core/imaging.py:316, 336-340` (entry), `:391-395` (exit resets `U` to `rho` and `taken` to False), `:401` (`n_ent` counts re-entries); `experiments/exp_astro_repeat.py:76, 358-368`; `results/astro/repeat.md:170-192`; `experiments/exp_astro_m2.py:221-231` (`choose`).
Evidence:
- **The re-run.** I re-ran m2r3 (110, 501) at the a priori default. It gave the same 2500 entries (160 s on the loaded machine). There were 155 distinct entries before tau = 2.1897. Then one lattice point, (110, 54), entered and exited 2345 times at that same tau until the cap. The run logged 4850 events: 2500 enters and 2350 exits.
- **What the readouts are.** The stop times for kappa 8 to 2 are 12.1 to 48.2, so the time rule never acts. Every kappa readout is the final block. The BIC block is 202 of 2501 with 150 units, against 1100 to 1900 blocks for uncapped default paths. Completeness on the two m2r3 default rows is 0.133 and 0.046.
- **Probable mechanism.** Exit is decided on the correlation at the slid centre (`g <= vtol*vscale`). Re-entry is decided on the correlation at the lattice start (`v > vtol*vscale`, `dt = (rho-U)/v = 0`). These can disagree forever.
- **Scope of the verification.** I re-ran only (110, 501). The other two capped runs show the same flag and were not re-run: (133, 281) isotropic default (BIC 508/2501) and (370, 441) (465/2501).
- **The same cap flag elsewhere (not re-run).**
  - In `m2_raw.npz`, 6 of 12 development-grid paths carry it: hp 8 at both r_ex, and hp 5 at r_ex 0.8, on both development patches. `choose` ranks them as valid rows (F 0.36 to 0.52 against 0.562 for the frozen configuration). The uncapped hp 5, r_ex 1.2 gives 0.555, so whether (3, 1.2, 2) is the grid optimum is not established.
  - Both M92 a priori default runs carry it (BIC 630/2501 and 293/2501).

Consequences:
- The default and BIC rows for both m2r3 patches, and the isotropic arm on (133, 281), are not the stated procedure. They sit in the pooled default tables (`repeat.md:18-31`), in P1's secondary lines, in P2's default pairs (for example the (110, 501) "tie +0.008") and in "default beats crowdsource default on 5 of 6".
- The report's "(a capped run is censored)" (`exp_astro_repeat.py:358`) understates this. The runs stop after about 150 distinct entries.
- Plan line 60 prescribes the handling: keep the first output, fix, rerun, report both.

Cheapest fix:
- Count distinct lattice points toward `max_entries`, and abort and flag when the same id re-enters at the same tau more than N times (or keep an exited point blocked until the residual changes). Add a regression test, since no existing test exercises this.
- Mark every capped row invalid in all tables now, and exclude them from pooled default rows and P1/P2 secondaries.
- Rerun the three runs, the M92 defaults and the six development-grid paths as labelled additions.

**3. major (wording). What the held-out frames say about K and about the M53/M92 failure.**
Location: `results/astro/repeat_plan.md:15-17, 56-58`; `repeat.md:49-54`.
Evidence:
- **Dose-response.** Kernel gain against stamp misfit:

  | patches | isotropic stamp rms | ellipticity | tuned gain |
  |---|---|---|---|
  | M2 | 1.0e-3 | 0.003-0.005 | 0 by construction |
  | m2r2 (133, 281) | 4.4e-3 | 0.04 | -0.004 |
  | m2r2 (195, 689) | 5.5e-3 | 0.09 | +0.076 |
  | M53 (four patches, post hoc) | 6.5e-3 to 8.2e-3 | 0.09 to 0.12 | +0.07 to +0.13 |
  | M92 (1100, 220), (1430, 550) (post hoc) | 7.3e-3, 8.4e-3 | 0.10, 0.13 | +0.07, +0.21 |

  The held-out gain is two to six times smaller than on the post hoc patches and not robust (finding 1). It adds two points to a post hoc dose-response and is no independent confirmation.
- **Confound.** Every frame the rule treats as non-isotropic is a camcol-6 frame: M53 5390-6-162, M92 4682-6-208 and m2r2 7778-6-72. Every round-stamp frame is camcol 2: 2583-2-136, 4187-2-105 and 3388-2-75. Stamp shape cannot be separated from camcol 6 (gain 4.895, dark 0.9025, saturation 52900, stripe 83) with these data.
- **What the result says.** Seeing of 3.0 to 3.5 px does not by itself reproduce the M92 loss (-0.104 against tuned crowdsource) or the M53 near tie (+0.031). The isotropic-kernel ISS beats tuned crowdsource on both m2r2 patches (+0.073, +0.102) with elliptical stamps and on both m2r3 patches (+0.181, +0.117).
- **What it does not say.** It is silent on what else differs in M53 and M92 (field, crowding, truth depth, background, offsets). It does not exercise the rule between stamp rms 1.9e-3 and 4.4e-3. The m2r1 stamps (1.83e-3 and 1.87e-3) sit 7 to 9 % under the 2e-3 goal.
- The plan's reduction (line 58) should be the default wording.

Cheapest fix: state K as "post hoc on M53/M92; one held-out elliptical-stamp frame gives a gain that does not survive alignment of the truth; rule action is confounded with camcol 6; no held-out frame reproduces a loss to a baseline".

**4. minor. P2 is true as declared, but its scope is narrower than the sentence.**
Location: `repeat.md:18-31`; `experiments/exp_astro_repeat.py:267-271`; `repeat_plan.md:34`.
Evidence:
- **Robust.** The ISS beats tuned crowdsource and tuned photutils on every patch, by paired interval.

  | level | vs crowdsource (tuned) | vs photutils (tuned) |
  |---|---|---|
  | pooled | +0.173 | +0.214 |
  | per patch | +0.069 to +0.241 | +0.134 to +0.277 |
  | per frame | +0.235, +0.122, +0.143 | +0.265, +0.200, +0.162 |

  - The F ratio ISS/crowdsource is 2.2, 1.7 and 2.8 on m2r1, m2r2 and m2r3 (2.0 on M2), so there is no erosion with seeing.
  - After aligning the truth, pooled margins are +0.170 and +0.207 and the minimum patch margin is +0.083.
  - Match tolerance 0.35 px/0.5 mag, 0.5/0.3 and 0.75/0.75 all keep the ordering on every patch.
  - In every 1-mag bin from 15 to 22, ISS completeness is higher (0.77/0.57/0.31/0.11 against crowdsource 0.66/0.21/0.06/0.01), with purity comparable or higher.
  - Removing per-method magnitude offsets does not change the ordering.
- **Bright-end exception.** At a depth cut of m < 20, crowdsource's F exceeds the ISS's on (133, 281): 0.458 against 0.345 (isotropic 0.431). At m < 21 and fainter the ISS wins on every patch.
- **Purity.** The ISS is less pure than photutils at m < 22 (d purity -0.064 [-0.098, -0.028]). Over all magnitudes it has 3923 unmatched detections against 1165 (crowdsource) and 735 (photutils). The win is completeness.
- **"Tuned" is thin.**
  - "Tuned crowdsource" equals the default minus five detections: F 0.160 in both rows.
  - The tuned family is post-filters of the output, not the detector's controls.
  - The baselines get no high-pass front end. The plan lists that family as planned (R3 F7), but the report header does not say so.
- **Not independent trials.** The six patches are four sky windows, and the HUGS truth is shared with M2 and across frames. Windows (310, 630) and (360, 1040) each appear in three frames. The replication unit is the exposure, three of them, two at broad seeing. This supports "holds on three further exposures of the same four windows". It does not support other fields.
- **Plan item missing from the tables.** The plan asks for the default and BIC readouts against "the baselines' defaults". `PAIRS_MAIN` has crowdsource-default pairs only, with no photutils-default pairs.

Cheapest fix: add one scope sentence (the baselines' tuning family, no high-pass front end, shared windows and truth). Add the photutils-default pairs. Report the m < 20 counter-example.

**5. minor. P3 is a declared primary outcome with no verdict, and by the plan's literal rule it is mixed.**
Location: `repeat.md:157-168`; `experiments/exp_astro_repeat.py:349-357`; `repeat_plan.md:48-50`.
Evidence:
- At the tuned configuration, isotropic against general is 1.20 against 1.05 on (133, 281) and 1.22 against 1.17 on (195, 689), so the isotropic value is above on both.
- At the default, (195, 689) is 1.50 against 1.55, so the isotropic value is not above.
- P3 does not track P1 patch by patch. (133, 281) has the largest scale improvement (default 2.32 to 1.43) and no F gain. (195, 689) has the F gain (+0.076) and almost no scale change.

Cheapest fix: print the P3 verdict per readout and state the mismatch with P1.

**6. minor. "Scaling the lengths by the seeing does not help (0.328 against 0.333)" holds for the tuned configuration only.**
Location: `repeat.md:204-216`; `experiments/exp_astro_repeat.py:110-116, 387-392`.
Evidence:
- **Tuned configuration.** Frame-level F is 0.433 vs 0.434 (m2r1), 0.307 vs 0.303 (m2r2) and 0.224 vs 0.203 (m2r3). The pooled table hides that m2r3 is worse. Kappa was not re-tuned for the scaled lengths.
- **A priori default.** Scaling raises pooled F from 0.297 to 0.305. On m2r3 it raises F from 0.134 to 0.184, and on (110, 501) from 0.084 to 0.174. It removes all three cap hits.

Cheapest fix: restrict the sentence to the tuned configuration and give frame-level rows. Cross-reference finding 2.

**7. minor. The freeze does not cover the decision code, and files were edited while the run was going.**
Location: `experiments/exp_astro_repeat.py:51-52` (FILES), `:139-145` (`frozen_code_hash`).
Evidence:
- FILES omits `exp_astro_repeat.py`, `exp_astro_m2.py` and `exp_astro_psf.py`. The code hash covers named functions only. `p1_text`, `method_table`, `boot_*`, `verdict` and the seed are unhashed.
- No harm found. `exp_astro_repeat.py` is unchanged since f05ca8d. `exp_astro_m2.py` (09:32) and `exp_astro_psf.py` (09:33) were edited during the run, but only in docstrings and `report()`. Their `boot_*`, `verdict` and `cpf` are untouched. The run's modules had been imported at 09:18.

Cheapest fix: put git blob ids of every imported repository file in the freeze record.

**8. minor. Part of the chronology evidence lies outside git, and the run log has no timestamps.**
Location: `results/astro/repeat_spec.md:40, 445-449, 405-420`; `repeat_plan.md:5`.
Evidence:
- **Pre-git order.** The declared-rule-first order is evidenced only by coordinator files outside git: `repeat_spec_v1.time` = 07:58:49Z, `rule_section_v1.sha256` = 484c03db, and survey logs at 08:09/08:14/08:16. I checked the hash against the scratch rule text (match). The committed section 1 differs from it by one blank line. The committed spec contains neither the hash nor the time. The first git commit containing the rule is 4b51b24 (09:14:16).
- **Run start.** `repeat_run.log` has no timestamps. The start is bounded only by the freeze file mtime 09:17:57, commit 09:18:04 and push 09:18:05 (origin reflog), with the raw file at 10:16:43.
- **"Never scored."** The plan says the baselines "were never scored against the truth". Spec §7 printed position-only counts against the truth. This was disclosed in the spec and is weaker than F.
- **Cannot be verified from git.** The R3'/R5' refinements being made "after survey tables, before any pixel result" is supported by the scratch logs and the explore scripts e1-e9 (no pixel statistics beyond the saturation maximum). It cannot be verified from git.

Cheapest fix: add the rule's sha256 and time to the spec, timestamp the start and end lines of the run log, and say "position-only counts" in the plan.

**9. minor. The report omits some items the plan lists.**
- The "stars brighter than 19.5 mag matched" descriptive item is absent.
- There is no P3 verdict (finding 5).
- Photutils-default pairs are missing (finding 4).

**10. minor. Data audit: nothing found that favours the ISS in the comparison with the baselines; four conventions matter.**
- **PSF stamp against the frame's stars** (6 to 14 isolated stars per frame, stamp stretch fit).

  | frame | width stretch | flux ratio |
  |---|---|---|
  | m2r1 | 0.983 | 0.94 |
  | m2r2 | 0.991 | 0.97 |
  | M2 | 0.999 | 1.01 |
  | m2r3 | 0.926 (mad 0.03, n = 8): stamp about 7 % wider than the stars | 1.00 |

  The reported seeing labels are stamp FWHMs. The adequacy rule checks the isotropic fit only against the stamp, so it cannot see a stamp-versus-star mismatch.
- **Zero point.** M2's NMGY is used for every frame, so detected magnitudes are 0.08, 0.03 and 0.07 mag brighter than the frames' own calibration (m2r1, m2r2, m2r3). Method-specific biases of the median (m_det - m_true), true mag 17 to 20, are larger than that:

  | patch | ISS adq | ISS iso | crowdsource | photutils |
  |---|---|---|---|---|
  | m2r2 (133, 281) | +0.23 | +0.14 | +0.12 | +0.07 |
  | m2r3 (370, 441) | -0.19 (adq = iso) | -0.19 | -0.13 | -0.15 |

  Removing them leaves P2 unchanged and moves P1 (finding 1).
- **Same-seeing control.** m2r1 (+6 % FWHM) loses 0.08 and 0.06 in ISS F against M2 on the same windows (0.653 to 0.572, 0.417 to 0.355). That is larger than the P1 effect. Cross-frame F is not controlled. The m2r1 truth offset is +0.12 to +0.14 px in y, against 0.00 to 0.05 px in y for M2. The M2 truth is header-mapped and the repeat truths are asTrans-mapped (differences at most 0.09 px on the M2 windows).
- **Proper motion.** None is applied (disclosed). By arithmetic, 4 mas/yr over 6 to 13 yr is 0.06 to 0.13 px, coherent. It does not explain m2r2's 0.25 px.
- **Saturated stars.** HST-saturated stars (nf606 = 99) are 0.5 to 1.9 % of truth below 22, and unmatched detections are not concentrated near them (for example 42 % within 10 px against a 41 % area share).

## Verified

**Chronology and integrity.**
- Commit order by `git log` and reflog:
  - 4b51b24 and b1a665a (data, plan) at 09:14:16, pushed 09:14:19 (origin reflog).
  - f05ca8d (script and freeze) at 09:18:04, pushed 09:18:05.
  - 851e5a7 (raw and report) at 10:17:56.
  - 851e5a7 touches only `repeat.md` and `repeat_raw.npz`.
  - The plan and `freeze_repeat.json` are unchanged since commit.
- The files in the freeze record (code files, plan, M2 freeze, checksums, baseline catalogue) are byte-identical between f05ca8d and HEAD.
- `check_freeze()` passes at HEAD and on a `git archive` snapshot of f05ca8d. `frozen_code` fbd56f73... is identical. The nine data-file hashes and `fwhm_ref` match.
- The baseline catalogue file (mtime 08:24:26) predates the script (09:15), the freeze (09:17:57) and the run.
- The baselines re-run bit for bit: all 14 arrays on all six patches, with photutils and crowdsource in the domain venv, which closes the spec §9 "not checked".
- File mtimes and `psf_raw.npz` give no sign of the kernel rule being evaluated on a repeat stamp before the freeze. `psf_general.py` was last modified 00:39 on Nov 8, before any repeat frame was fetched (07:59).
- The log shows a single "6 patches to run, 0 exist", so there was one invocation.
- The sha256 of all frames, psFields, the M2 frame and the HUGS catalogue match `checksums_repeat.json`.

**Numbers.**
- My scorer (no import of `astro_common`; max-weight matching per connected component, cross-checked by Hopcroft-Karp cardinality) reproduces 60 rows exactly, totals and 25-block arrays. These are ISS kappa 2 and 5 for iso, adq and scaled, plus four baseline rows, over six patches, with baselines read through the cs_dets filter.
- All 348 stored `blk22`/`tot` pairs are internally consistent. The truth side is identical across methods on each patch.
- Pooled F matches the report: 0.333, 0.324, 0.160 and 0.119. My bootstrap (B = 5000, stratified by patch) gives adq minus crowdsource +0.173 [0.158, 0.188], adq minus photutils +0.214 [0.198, 0.229] and P1 [+0.0061, +0.0615]. The repo procedure at seed 0 gives [0.00604, 0.06164].
- `repeat.md` regenerates byte for byte from `repeat_raw.npz` (sha256 7b9406...).
- ISS re-run: m2r2 (133, 281), tuned configuration with the rule's general kernel.
  - The kappa-2 catalogue is bitwise equal (1209 units).
  - Entries 1253, BIC block 1013 of 1254 and rsd 1.0547 are identical.
  - All nine readouts (kappa 2 to 8 and BIC) agree with the stored totals under my scorer.
- Aliasing is legitimate: `X.run_path` and `PS.run_path_psf` with an isotropic fit are bitwise identical on a 50 x 50 crop, for all kappa and BIC.
- In `exp_astro_repeat.py` I found no bug in the scaled-lengths run, the P1 selection, the pooling or the ladder's window matching. The ladder values match my recount.

**Data (domain venv, astropy).**
- Patches rebuilt from the frames with my own bilinear sky interpolation: dn, sky_dn and nmgy are bit-exact for all six patches. Gain and dark match the psField. PSF stamps agree with pydl to 3e-9 to 1.8e-8.
- HUGS truth through the frame header WCS against the stored lists: m2r1 within 0.03 px, m2r3 0.006 and 0.08 px, and m2r2 0.15 to 0.17 px, as spec 6.2 states. The M2 truth file equals the header mapping exactly. The truth rows are the same HUGS rows in all frames.
- Orientation and windows:
  - m2r2 is rotated 1.68 deg and the others at most 0.02 deg.
  - The repeat patch centres are 0.32 to 0.55 px from the M2 window centres (header WCS).
  - There is no residual rotation or scale trend across the patches (slopes at most 0.16 px per 100 px; a missed 1.7 deg would be 1.7).
- Astrometric offsets confirmed by star fits (truth to image):
  - m2r1 (+0.115, +0.133)
  - m2r2 (-0.018, +0.227)
  - m2r3 (+0.019, +0.026)
  - M2 (+0.097, +0.019)

**Tests.** `tests/test_exp_astro_repeat.py`: 18 passed. `tests/test_astro_repeat.py`: 20 passed with the domain site-packages, and 19 passed with 1 skipped in the repository venv alone.

## Could not check

- The exact start time of the run (the log has no timestamps; only bounds, finding 8). I could not confirm the push reached the remote; the reflog is local.
- The scratch-directory evidence for the pre-git order (07:58 to 08:24) is coordinator-written. I could read it but not authenticate it.
- I re-ran one ISS patch (all kappa) and one capped path. I did not re-run the other five patches. The BIC-readout detections are not stored, so BIC rows beyond (133, 281) adq tuned can be checked only for internal consistency.
- The limit-cycle mechanism is verified on (110, 501) only. The other two capped runs, the six development-grid paths and the two M92 default runs carry the same cap flag and were not re-run.
- Whether the post hoc M53/M92 kernel gains have the same truth-offset sensitivity as P1 was not tested.
- I did not re-implement asTrans. I used the header WCS and image-based star fits. The star fits agree with the spec that asTrans is closer than the header mapping for m2r2 (residual 0.23 px against about 0.33 px).
- Proper motion (the cluster value was not looked up) and `fpM`, cosmic-ray and satellite-trail masks were not checked.
- Bitwise reproducibility under other BLAS thread counts was not tested. My re-runs used one thread.
