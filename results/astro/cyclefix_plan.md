# Plan for the guarded re-runs of the astronomy paths that the enter/exit limit cycle of `imaging.sliding_path` made invalid (declared before any of them was re-run)

Declared: **2026-11-08T11:28:54Z** (`date -u`, taken immediately before this file was written). Author: the agent that carries out the coordinator's brief `agent_cyclefix_brief.md`. Nothing in this file was
decided after a guarded result of any published path was seen; the state at the time of writing is in section 0. The experiment script `experiments/exp_astro_cyclefix.py` records the sha256 of this file in
its raw (`meta|plan_sha256`) when it runs, and the report prints it.

## 0. State at the time of writing

* **No published astronomy path has been run with the guard to its end.** Guarded or diagnostic runs that exist: the coordinator's first check of `src/core/imaging_guard.py` (m2r3 (110, 501), a priori default, `max_entries`
  400) and this agent's reproduction of the stall without the guard (same patch and configuration, `max_entries` 200: the lattice point (110, 54) enters and exits at tau = 2.1896895598963453 until the cap, 250 events, 8 s).
  No result of this plan exists.
* The files that the freeze records `results/astro/freeze.json` and `freeze_repeat.json` hash, and the raw results, are as follows (sha256 at declaration; they are read-only for this work and are checked again at the
  end; `exp_astro_m2.check_freeze(m2_raw)` and `exp_astro_repeat.check_freeze()` both pass now):

  | file | sha256 |
  |---|---|
  | `src/core/imaging.py` | 23ec17ea1e85f43161b4d6262d5555c37d9ddccfd6c3e13d880ba3f51f5e5be6 |
  | `src/core/psf_general.py` | 5e253b9303271f4350d83a8c2dced6dd329567ddbcef3bc3b0ca31db6443dab9 |
  | `src/core/imaging_guard.py` (new, the coordinator's guarded copy) | 5c9cf8b81781d766981552fc1fcc05c533b9c1477377e92da643357e75789ed5 |
  | `experiments/astro_common.py` | a25d698f21adfc0a15647477715bb64238c981900470310eabc8dbbd86e5d1da |
  | `experiments/astro_second.py` | 3688e6b5f424c1de650606fe82d64daa67f1b55e3d4c5f7378c53bbfcb6e953f |
  | `experiments/astro_repeat.py` | fc1186670adcc824a76c7bc69edcbc54da8809bf6636b78deba37267d8148e81 |
  | `experiments/exp_astro_baselines.py` | 4c5e7451247d73857b77b7d521c83ae682f87fc58ea9c0ec98d0177c04061b32 |
  | `experiments/exp_astro_m2.py` | 8acf3c302e48f7ef3e1c93c98af3c120f954e72a7fac5a0ece788a34a3b7de74 |
  | `experiments/exp_astro_psf.py` | 51e3fb651c528518536a7bc4d71e698a9c8af0653e1fd63dce46c97e6bda7718 |
  | `experiments/exp_astro_repeat.py` | f22922d054f839112a1d67885475f3910e0bd77b808405d1d4c8dca3228a882f |
  | `results/astro/m2_raw.npz` | b77dedee026a4a3b45816a51b1f2807921a849ab72a3072fdff20551fa94e90b |
  | `results/astro/psf_raw.npz` | 089c4ca289cd347c44a70792451bc1c7f86f6eed5def1ab864064b53da9de144 |
  | `results/astro/repeat_raw.npz` | 0f54188187d766c94f7c1611d4e38879abcf6528cc8a02d78b10aaeaf8ea6071 |
  | `results/astro/baseline_catalogs.npz` | 071f5b8f873c84eb35ce02c56f427ddb98447e7d8b164be33b212506d4649b00 |
  | `results/astro/baseline_catalogs_second.npz` | 2eed2f5fcc6d9e3550a8ffdef8854e90897993e68e925dfb4bcf395a32f90a01 |
  | `results/astro/baseline_catalogs_repeat.npz` | 5d52dd166233bc1c742d02e3b84a8470eb4e0ed71d129cc67aebf22a45abf092 |
  | `results/astro/repeat_plan.md` | b409c4ec06cbae513c891df1c04ecd66ebdbdf1e7ace934f8db1b28a51fefbdb |
  | `results/astro/freeze.json` | 1fc355bb0b28b44a7d03f8f303a466406e7451d8a14435c4015d360f0b2bd6f7 |
  | `results/astro/freeze_repeat.json` | 2e64698a48382c8b0939ae2beb267ba12c2f669494a961cece14b4850ef9be1b |
  | `results/astro/m2.md`, `psf.md`, `repeat.md` | 34ecf2e2228e64a0395f384272f939c1a5c9b689cdf7f04e020cef24cb2cb684, 4f89296229f93178abfa63aa8efce0b906c1c270515410f5a6eeaf2999e5488b, 7b9406439cc26fa5d1c70e62a6e2e454279eb2d2a5b2ea4780c53e0e4a9bb944 |

* Environment: Python 3.13.18, numpy 2.5.4, scipy 1.18.1 (the repository venv and the domain venv are the same versions); `OMP_NUM_THREADS=1`, `nice -n 10`, `PYTHONDONTWRITEBYTECODE=1`, at most two heavy processes of this
  work at a time on the shared 4-core machine.

## 1. The defect and the guard (semantics fixed here)

`imaging.sliding_path` lets a lattice kernel enter when its dual reaches the wall; the entry test uses the correlation at the kernel's *lattice start* (`v`). A unit leaves (`gone`) when its coefficient is 0 and the
correlation at its *slid centre* (`g`) is not positive; its dual then goes back to the wall (`U = rho`), so a kernel with `v > 0` re-enters at `dt = 0`, at the same entry time. If the block fit of the crowded
cluster returns its coefficient 0 each time, the path loops at one tau until `max_entries` (2500); the time rule (`tau_max`) never acts because tau does not advance. Such a run stops after about 150 distinct
entries (repeat patch m2r3 (110, 501), a priori default: lattice point (110, 54), 2345 enter/exit turns at tau = 2.18969) whatever its stopping rule says. The cap flag of the raw files (`bic = [block, hist length, 1]`,
hist length 2501) is the signature of the cycle, not a harmless censoring.

The guard, `imaging_guard.sliding_path_guarded` (default `cycle_guard=True`): a lattice point that has *exited twice at the same entry time* (within 1e-12 relative) is barred from entering for the rest of the path; the
result has the extra key `barred = [(tau, (b, a))]`. A first exit followed by a successful re-entry is untouched, so a path without a double exit at one entry time is the unguarded path bit for bit;
`cycle_guard=False` is `imaging.sliding_path` bit for bit. In `order='greedy'` tau is the step count and advances at every entry, so no two exits share an entry time and the guard never acts (a limitation, stated in the report;
the astronomy paths and the SMLM `fit_frame` use `order='dual'`; the one greedy astronomy path is the design variant `V|x310_y630|greedy_snr`).

The re-runs call `exp_astro_m2.run_path`, `exp_astro_psf.run_path_psf` and the repeat code of the frozen modules unchanged; the only difference is that `imaging.sliding_path` is replaced, in the worker process, by a
thin wrapper that calls `sliding_path_guarded` with its defaults and keeps `barred` and the event counts in a side channel (the frozen functions drop the path's result). Every call site used looks the function up at call time
(`IM.sliding_path` in `exp_astro_m2.run_path` and `exp_astro_psf.run_path_psf`; `sliding_path` in `imaging.fit_frame`).

## 2. Which paths are affected (scan, then list)

The list is derived from the raw files by scanning every `bic` key whose third element (cap flag) is 1. Expected, and declared here: 11 distinct path runs.

* `m2_raw.npz`, development grid: `D|x110_y630|hp5|rex0.8`, `D|x110_y630|hp8|rex0.8`, `D|x110_y630|hp8|rex1.2`, `D|x210_y630|hp5|rex0.8`, `D|x210_y630|hp8|rex0.8`, `D|x210_y630|hp8|rex1.2` (6);
* `m2_raw.npz`, M92 a priori default: `T2|m92|x1100_y220|default`, `T2|m92|x1430_y550|default` (2);
* `repeat_raw.npz`, a priori default, isotropic kernel: `R|m2r2|x133_y281|iso|default`, `R|m2r3|x370_y441|iso|default`, `R|m2r3|x110_y501|iso|default` (3; the two m2r3 runs are also stored as `adq`, aliased: the rule keeps the
  isotropic kernel there, `kernel|same` = 1);
* `psf_raw.npz`: none.

**Added by the scan beyond the brief's criterion (found while writing this plan, from the stored entry counts):** the design variant `V|x310_y630|norex` in `m2_raw.npz` (no exclusion zone, M2 test patch (310, 630), tuned
configuration) stored 2500 entries (`V|x310_y630|norex|secs = [1606 s, 2500]`), i.e. it reached the cap, but the variants part stores no `bic` key, so the `bic` scan cannot see it. Whether that is a cycle or a long path is not
known; it is re-run with the guard in stage 1 as a twelfth path (labelled as an addition), because with no exclusion zone, double exits are plausible. The scan code also lists every stored path whose entry count is the cap, whatever
key stores it.

## 3. Stage 1: the affected paths are re-run with the guard (item 1(a))

Each of the 11 (+ the variant) is re-run exactly as the original run: same patch, same configuration (sigma_b, r_ex), same kernel (the isotropic mixture of `imaging.fit_psf_parts` through `exp_astro_m2.run_path`; the
variant through `run_path` with the overrides of `exp_astro_m2.VARIANTS`), same weights (Poisson-plus-dark from the smoothed data), same lattice (h 0.5, margin 2), same fit radius 4, `tol` 1e-6, one pass,
stopping time 1/(2 x noise scale) (kappa 2, the smallest of the grid), `max_entries` 2500, same readouts (kappa 2 to 8 and key 0, the BIC block) and same scoring (`astro_common.score` against the patch's truth).
The original rows are not touched; the new rows are stored under a namespace of their own (`G|` + the original key) in the new write-once raw `results/astro/cyclefix_raw.npz`, and the report shows both side by side:
entries, hist length (BIC block / blocks), the cap flag, the run time, the guard activations (entry time and lattice point of each barred point), the numbers of exits and of distinct lattice points that entered, the
largest number of entries of one lattice point (a value of 5 or more is flagged as possible residual cycling that the two-strike guard does not remove), and F, completeness, purity at the reported readouts.

Reading rules for stage 1 (declared):

* A guarded re-run that again reaches 2500 entries **without** any guard activation is a legitimate censoring of a long path; it is reported as such with its entry times (the stored `taus`). If it reaches 2500 entries **with**
  activations, the report says that the guard acted and that the cap was reached all the same, and shows the entry times; its rows are flagged as capped in every table, as the as-run table flags them.
* A re-run whose as-run cap was a legitimate censoring (no activation expected) must reproduce the as-run row bit for bit; if it does not, the difference is reported.
* The corrected rows are labelled **"guarded re-run (post hoc correction of a defect found after the runs)"**; the as-run rows are labelled "as run (entry cap hit: an enter/exit cycle, invalid)". Nothing is deleted or overwritten.

## 4. Stage 2: the corrected development table and the choice (item 1(c))

The 12 development paths (two development patches x sigma_b 3, 5, 8 x r_ex 0.8, 1.2) are taken with the 6 capped paths replaced by their guarded re-runs, and the 8 kappa readouts of each are scored with
`exp_astro_m2.choose` unchanged (largest pooled F at mag < 22 over (sigma_b, r_ex, kappa) of the two development patches' counts added; ties on `round(F, 6)` go to the larger kappa; capped rows are
treated exactly as `choose` treats them, that is, as valid rows). The report prints the as-run table and the corrected table side by side and the argmax of each. Declared outcomes:

* **(A)** the argmax is the frozen configuration (sigma_b 3, r_ex 1.2, kappa 2): nothing else changes; the report says so (and states the corrected pooled F of the frozen point and of the best other point).
* **(B)** another configuration wins: the report says so with its pooled F. Only then, as a *labelled post hoc addition* (no claim that it is held out), that configuration is run once with the guard on **all 16 test patches** (the brief says
  14; the test patches of M2 (4), M53 (4), M92 (2) and the six repeat patches are 16) and scored beside the as-run results with the paired differences to the unchanged baselines; the run uses the same
  functions and the kernel as the as-run runs (isotropic on M2, M53, M92; both kernels on the repeat patches where the rule picks a general kernel).
* Descriptive robustness line, not the outcome rule: the argmax over the rows of paths that did not reach the cap (the guarded re-runs that are capped again excluded), and the corrected F of the frozen point.

## 5. Stage 3: the corrected rows that the paper uses (item 1(d))

Computed by the script from the raw alone, with the frozen scoring, matching and bootstrap functions (`astro_common.score`, `exp_astro_m2.boot_*`, `boot_cpf`, `boot_diff`, `verdict`, `fmt_iv`, and the repeat report's helpers
`pooled_blocks`, `method_table`, `diff_row`, `p1_text`'s method) and the same block bootstrap (20 x 20 px blocks, 2000 resamples, seed 0, stratified within each patch for pooled rows):

1. the M92 a priori default and BIC rows, pooled and per patch, with the paired differences to tuned/default crowdsource and photutils as the M2/M53/M92 report prints them (the report section of `exp_astro_m2.report` for
   M92 is regenerated from the raw with the guarded default rows in place of the capped ones and printed in full beside the as-run section);
2. the repeat frames' default and BIC rows, the stamp kernel and the isotropic kernel, pooled over the six patches and by frame, with the paired differences to the baselines' defaults as `repeat.md` prints them
   (`PAIRS_MAIN`); the photutils-default pairs, which `repeat.md` does not print, are added and labelled as added after the plan;
3. the secondary lines of P1 (a priori default and BIC stop) of `repeat.md` recomputed with the corrected isotropic arm on m2r2 (133, 281): the pooled dF of the chosen against the isotropic kernel over the two patches
   where the rule picks a general kernel, 95 % intervals of the same block bootstrap with 2000 resamples (seed 0) and, beside them, with 40000 resamples (seed 0), with the wording of the plan of the repeat test
   (a positive interval: "confirmed on held-out frames", an interval containing 0: "not confirmed", a negative one: "contradicted"). The as-run lines are reproduced by the same function as a check against `repeat.md`.

Extras computed from the same re-runs and labelled as added after the plan: the residual scale (P3) of the corrected isotropic default catalogues of the repeat patches, and the residual scale of the guarded M92 default catalogues
(`psf.md` uses the capped catalogues in its "frozen" column).

All of them carry the label "guarded re-run (post hoc correction of a defect found after the runs)".

## 6. Stage 4: check mode (item 1(b))

Every other published astronomy path is re-run with the guard, with its original configuration, and compared with its stored result. There are 66 such paths:

| family | paths | count |
|---|---|---|
| development grid, not capped | `D|x110_y630|hp3|rex0.8`, `hp3|rex1.2`, `hp5|rex1.2`, `D|x210_y630|hp3|rex0.8`, `hp3|rex1.2`, `hp5|rex1.2` | 6 |
| M2 test patches | the 4 patches x tuned, a priori default | 8 |
| M53 test patches | the 4 patches x tuned, a priori default | 8 |
| M92, tuned | the 2 patches | 2 |
| design variants of the path on the M2 patch (310, 630) | `greedy_snr`, `nohp`, `snr`, `nocleanup` (`norex` is in stage 1) | 4 |
| kernel study (`psf_raw`, adequate kernel) | the 6 patches of M53 and M92 x tuned, a priori default | 12 |
| repeat patches | isotropic, tuned (6); isotropic, default (3 uncapped: (316, 124), (316, 224), (195, 689)); chosen kernel, tuned and default on the two patches where the rule picks a general kernel (4); lengths scaled by the seeing, tuned and default (12) | 25 |

(the aliased `adq` runs on the patches where the rule keeps the isotropic kernel are the same runs as the `iso` ones and are checked once). Order of priority, if the CPU time does not allow all of them: the tuned and default runs
on the 16 test patches and the 6 uncapped development paths first, then the kernel study, then the repeat patches' tuned and scaled runs, then the variants; what was not checked is listed in the report.

**Identity rule (declared).** A path is *identical* when its stored entry count (`secs[1]`), its `bic` triple, the scores of every one of its readouts (every array `exp_astro_m2._pack` stores, for kappa 2 to 8 and key 0),
the catalogue at the reported readout (`det`, where stored) and the stored residual scales (`rsd`) are equal to the stored ones **exactly** (`numpy.array_equal`, NaN equal to NaN). The run time is not compared. The number of
guard activations is recorded for each path. A path with at least one activation had a double exit at one entry time: its unguarded path took a further entry at that time, so it is *expected to differ* (its entry
count differs) even if its detections do not; a path without activation must be identical, because the guard cannot act. Any path that differs is reported with F (and completeness, purity) before and after at its
reported readouts, the entry counts, and the activations.

## 7. Item 5: the SMLM spot check

`experiments/exp_cycle_scan.py` (its own raw and report, `results/smlm/cycle_scan_raw.npz`, `cycle_scan.md`): the first 40 test frames of each of the eight sets (the four 2016 sets, the four 2013 sets) at the tuned
(E-optimal) operating point and at the a priori default, fitted as the experiments fit them (the same functions and the same calibration), run with the guard: the number of frames with a double exit at one entry time (read
from `barred`) per set and operating point, and, for every frame, the detections of the guarded and the unguarded run compared (the unguarded run is made for every frame that has an activation and, as a control of the
equivalence claim, for a sample of frames without one). Reported: activations, and any change of detections (frames, emitters). The E/J tables are not recomputed.

## 8. What is reported, and what is not claimed

The report `results/astro/cyclefix.md` (regenerated byte for byte from the raw by `--from-raw`) states: the defect and the guard; the as-run and guarded rows side by side; the corrected development table with the
outcome (A) or (B); the corrected rows of section 5; the check-mode table (identical or not, per path, with the activations); every limitation. Not claimed: that any frozen row, report or freeze record is replaced (the
coordinator decides which of `m2.md`, `psf.md`, `repeat.md` and the paper's numbers are regenerated, and how the first outputs are kept); that the guard is the only possible repair (it removes the repeated entry at one
entry time, not its cause, the disagreement of the entry test at the lattice start with the exit test at the slid centre); that the guarded rows of a patch are held out (the defect was found by the referee after the
results had been seen, so every corrected row is a post hoc correction and says so); any change of the SMLM tables. Anything the report adds that is not in this plan is labelled "added after the plan".

## 9. Amendments (each before the first guarded result of any published path; nothing above was changed)

* **Amendment 1, 2026-11-08T11:37:13Z.** The check-mode table of section 6 lists 6 + 8 + 8 + 2 + 4 + 12 + 25 = **65** paths, not 66 as the sentence above the table says (an arithmetic slip; the table is right). The script's inventory
  of the stored path runs gives 77 in all: 12 in stage 1 (the 11 flagged by `bic` and `V|x310_y630|norex`) and 65 in check mode.
* **Amendment 2, 2026-11-08T11:37:13Z.** Documentation-only edit of `src/core/imaging_guard.py`, no change of code or semantics: the module docstring now lists the one line that is not a guard line but must change with them (the line that
  adds the exclusion zone to `blocked`), and states the limit of the rule (it keys on the entry time, so in `order='greedy'`, where tau advances by 1 at every entry, it never acts; only repeated entries at one entry time are removed); the
  function docstring documents the result key `barred`. sha256 of the file after the edit: `8f67aff1c8b9e7ff2d5102ebb132956fd9f4323f42b2a206b89b464e98242ae0` (the table in section 0 has the hash of the file before).
