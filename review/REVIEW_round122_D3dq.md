# Referee report D3dq, round 122: Proposition D‴.61 as revised after round 121

## Scope

**What I reviewed.** Commit 04d9ead on claude/vigilant-clarke-hf939b, exported with `git archive 04d9ead | tar -x` into W/repo. W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref122, and W/repo/venv is a symlink to the repository's venv (python-flint 0.9.0).

**Hashes of the files reviewed** [own]:
- experiments/gh_lastann.py: 474ba8c9…
- results/chaos/gh_peano.md: e7ac9779…7b00
- results/chaos/gh_peano_raw.npz: 184cb64c…
- Round 121's patched file (ref121/fix/lastann.py, copied to W/fix/lastann_fix.py): 66422a18…, the hash round 121 reports.

**The repository moved during the round.** Other sessions advanced HEAD to f346010. Nothing of theirs is reviewed here.

**How I ran things.**
- I used read-only git only and wrote nothing to the repository or to ref121.
- Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, with GP_OUT set to W/out, pytest's `--basetemp` inside W, and at most two heavy processes of mine at a time.
- I killed no process; all of mine exited on their own.

**Read in full.**
- review/drafts/D61.md.
- gh_lastann.py, line by line.
- In exp_gh_peano.py: `part_VA`, `va_lines`, `vc_lines` and what they call.
- results/chaos/gh_peano.md, sections VC and VA.
- tests/test_gh_peano.py: the VA test, the series test and the VC label test.
- review/REVIEW_round120_D3do.md and review/REVIEW_round121_D3dp.md.
- The Round 120 and Round 121 entries of review/REVIEW_LOG.md.
- tools/mutants.py and tools/mutants/round121_D61_own.py.
- `git diff fbf237e 04d9ead` for these files, plus HANDOFF, PAPERS and drafts/README.
- I also looked at D58.md where (iv)'s upgrade touches it.

**What I reused from round 121.** ref121/att/ind.py, copied unchanged to W/scripts/ind.py, and the patched file, read-only.

## Summary

**Design.** This commit adopts round 121's NaN-safe `handoff_range` patch verbatim; only the comment changed (`diff` against ref121/fix/lastann.py) [own]. It adds:
- `outer_box`: the connector is always certified, and the certified (lo, hi) is returned;
- `handoff_check`, a containment check at 15 η per accepted η-box, plus an area-sum tiling check;
- `soundness`, which samples the forms below their box bounds;
- `bb`'s keep list;
- the asserted B0 check in `validate`;
- `arb(hi) − arb(lo)` in `mvbound`;
- a |σ| = 0 guard in `seg_bound`;
- a rewired `main` and `part_VA`.

On that basis it upgrades three claims:
- D61 (iii)(b) to the whole open disc |η| < 2, uniformly on compact subsets ("conjecture P proved, computer-assisted");
- (iv) to every frame with r < 1;
- D‴.58's sharp constant to every σ*/σ.

**Verdicts.**
1. **Reproducibility is exact.**
   - `--from-raw` regenerates the md byte for byte.
   - 42/42 tests pass (22.0 s).
   - `part_VA()` recomputes equal to the raw's 'VA' (2.9 s).
2. **The adopted certificate is round 121's patch plus changes that only record, check or tighten.**
   - It accepts exactly the patched file's boxes and segments.
   - The only difference is one form-B box whose float bound differs in the last digit, from the exact width.
   - It is sound.
   - The new checks (`handoff_check`, `soundness`) print true statements. But they do not test the whole joint, and they print the same statements for unsound variants (D3dq-1, D3dq-2).
3. **Independently, it is a proof on 1.995 ≤ |η| < 2** (computer-assisted; it assumes Arb/python-flint is correct and the hand-derived identities of rounds 120–121).
   - Round 121's ind.py on the adopted file reproduces round 121's numbers exactly.
   - My adversarial check at 1200 bits (1493 η, 1.62 M path points, every η-box's corners and edges, θ = 0, θ = π/2, ε down to 10⁻⁵⁰) finds 0 defects.
   - My structural and rigorous joint checks pass.
4. **D61's upgraded claims follow.**
   - (iii)(b) holds on |η| < 2, uniformly on compact subsets.
   - (iv) holds for every r < 1: c_need² = 4r ln(1/r)/(1 − r) < 4.
   - The analytic steps still rest on earlier rounds: (a); the Λ⁺ contour bound and the band near |σ| = 1; Duhamel/Laplace; route (2)'s 0-side and chord.
   - No overreach in substance. There are wording issues, and one cross-document inconsistency with D58.md (D3dq-3).
5. **Mutants.**
   - **Results.** Of 45 new mutants, 21 are killed and 24 survive; none of the survivors changes the raw or the md. The shipped L01–L08 are killed and L09 survives, as logged.
   - **What does the killing.** Every unsound mutant that the new checks miss (H10, O02, O05, T01, V02) is killed only by the test's exact raw equality.

**Findings.**
1. **D3dq-1 (minor).** `handoff_check` checks the connector that `outer_box` reports, not the segments it certified. It does not check y_c ∈ [lo, hi] or the polyline's start. Three unsound joints pass it with "0 outside":
   - O02 and O05 certify only part of the connector;
   - H10 starts the polyline 10⁻³ above the connector.

   `soundness` misses V02 (the mean value width dropped). In the adopted code the joint is right by construction (my J1, J2).
2. **D3dq-2 (minor).** The tiling check is an area sum. It is redundant with `failed == 0`, and T01 fools it: duplicating each lower θ-half in place of the upper keeps the area, and "tiled" stays True.
3. **D3dq-3 (minor, text).** REVIEW_LOG's D3dp-6 row lists "D‴.58's tiers (every rule with a frame, the sharp constant at every σ*/σ)" as applied. D58.md at 04d9ead is unchanged and still says t < 0.9900 (σ*/σ < 199), "conjectural beyond" (lines 14, 25, 171–192, 287, 310). Meanwhile HANDOFF l. 34 and PAPERS already state the upgrade, and HANDOFF l. 184 schedules D58 for round 123.
4. **D3dq-4 (cosmetic, D61/md wording).**
   - (iii)'s header still reads "proved on most of the disc".
   - "The true value below the bound" is the form's value fA/fB, not F_η's.
   - "1319 outside the annulus" includes 77 boxes skipped for X > δ₀.
   - "The joint checked at every η-box" means 15 sampled η per box.
   - "Arb's correctness assumed" should also name python-flint and CPython's float/nextafter.
   - The forms are "exact" only up to two one-sided simplifications.
5. **D3dq-5 (cosmetic, code and tools).**
   - `va_lines`'s conclusion ignores `tiled`.
   - "128 bits" is a literal.
   - The series test uses `overlaps`, not `contains`.
   - The shipped harness's md column is vacuous for `--exp gh_lastann`: the control L00 is reported "md: true".
   - HANDOFF l. 14 says "63 + 39".
   - The log says "9 mutants and a control" twice and calls L09 "equivalent" where it is only "equal on the data".
6. **D3dq-6 (minor).** Of this round's new code, 24 survivors (list below):
   - two unsound in principle but equal on the data: H07 (the NaN branch's lo = f1) and H08 (y_c's clamp removed);
   - the rest weaken checks: the top edge at the float π/2, tolerances 10⁻³ and 0.5, route 'd' unchecked, and others.

   The soundness plant (bounds lowered by 1.0) is too coarse to catch S03.

## What I checked and how

**(1) Reproducibility** [own]
- **The md.** `GP_OUT=W/out exp_gh_peano.py --from-raw`, then `cmp` against the committed md: identical (sha256 e7ac9779…).
- **The tests.** `pytest -q -p no:cacheprovider --basetemp W/pt/base1 tests/test_gh_peano.py`: 42 passed in 22.00 s (W/logs/pytest_gh.log).
- **Part VA.** W/va_recompute.py: `GP.part_VA() == list(raw['VA'])` is True, in 2.93 s.
- **Numbers in Numerics (VA) and the md.** Each one is in the raw [raw]:
  - inner boxes: 40, 2011, 320 (634 skipped), 769 (1319 skipped);
  - O: 76 η-boxes and 5193 segments, worst −1.4993e-7; end −1.3314;
  - validate: 1479 points, A_vs_def 9.3e-25, B_vs_centred 7.6e-23, A0 1.5e-59, B0 2.1e-59;
  - handoff: 1140 points, 0 failed, area_rel 0;
  - soundness: 80, 5358, 640 and 1538 points, worsts −7.8e-3, −2.9e-6, −5.9e-2, −2.4e-2.

**(2) The adopted file against round 121's patch** [own]
- **The patch.** `diff ref121/fix/lastann.py gh_lastann.py`: `handoff_range` is identical except for its comment. Every other hunk is one of the changes listed in the header.
- **Same acceptance** (W/scripts/cmp_fix.py: both files run, every accepted box's bounds and every `seg_bound` call recorded):
  - A-J, A-FC and B-J are identical box for box and bound for bound;
  - B-F has 769 = 769 boxes, one bound differing by 1 ulp (−0.049904430874173504 against −0.04990443087417351), the exact arb width against the float one;
  - O makes 6091 identical `seg_bound` calls, and its record differs only by the new `R1` field.
  - So on the data nothing is accepted that the patched file rejects.
- **Each change on its own** [own, reading]:
  - `arb(hi) − arb(lo)` encloses the exact width. The float version could round either way, so this is rigorous where that one is not; strictly, it may accept where float rounding up made the patched file reject, which is sound.
  - The `seg_bound` guard returns ∞ where the patched file set ∇C = 0 (non-rigorous). It is never reached.
  - The always-certified connector adds a degenerate segment when lo = hi, which can only reject more. No accepted box has lo = hi, so it is a no-op on the data.
  - `bb`'s keep and `recs` only record.
  - B0's assertion raises if fB0 < fB(·, 0) at a sample.
  - `main` orchestrates the same pieces. No test runs it; `all --eps0 0.11` still crashes in `validate` (ValueError in `Lam`, |v| > 1/2: D3dp-5's residue). The pieces alone at ε₀ = 0.11 pass: 54/2533/413/1189 boxes, O 220/16083 with worst −4.23e-9, handoff 3300 points with 0 outside and tiled, soundness 0 failed. This matches round 121 (W/cli).
- **The special cases Y = 0** [own, reading]:
  - On 0 ≤ θ ≤ π/2, Im A = ω²sin 2θ ≥ 0 and Re A < 0 for η ≠ 2, so κ ≥ 0 and 0 ≤ Y = min(κδ₀, θ) ≤ θ.
  - At θ = 0 this gives Y = 0. At the true π/2, sin π = 0 gives κ = 0 and Y = 0. Both are exact.
  - The NaN branch's [0, θ_hi] is therefore a valid enclosure.
  - `exact0` fires only at arb(f1) == 0 and at the top box's true π/2. Mutant L06 shows it is needed: without it the straddling ball fails 3 points.
- **Trying to fool the checks** (mutants, see (5); each of the following leaves the check at "0 failed"):
  - O02 certifies the connector only from y_c up;
  - O05 certifies it only up to y_c;
  - H10 starts the polyline at hi + 10⁻³;
  - T01 duplicates each split's lower θ-half, keeping the area and `tiled` True;
  - V02 sets the mean value width to 0, an unsound bound, and `soundness` still reports 0 above it.
  - All five are killed only by `va == list(raw['VA'])`.
  - `soundness` evaluates fA/fB themselves. These equal E/X² on the path, up to the one-sided (Re σ − 1)₊ in A and the dropped −P₊² in B, so it tests the box machinery (closed form, mean value step, dual enclosure), not the identities, which `validate` covers.
  - The tiling check sums float areas with tolerance 10⁻¹². It detects a missing box (the test's plant), but not overlap plus a gap of equal area.

**(3) Independent proof check** [own]
- **Round 121's ind.py, unchanged, on the adopted file** (W/logs/ind.log), with output identical to round 121's ind_fix.log:
  - 444 η, 437 720 points at 400 bits, 0 defects;
  - max E = −6.25e-47; inner true − certified ≤ −3.05e-16; 0 misses;
  - outer max −3.98e-5; 5369 outer segment calls × 10 points with true − bound ≤ −2.19e-5.
- **My variant ind_corners.py, at 1200 bits** (W/logs/ind_corners.log). Every η is checked against every accepted η-box that contains it, not just the first. The η are:
  - for each of the 76 accepted η-boxes, R ∈ {the lower edge, 2 − 10⁻¹⁵, 2 − 10⁻⁹, the middle} × θ ∈ {f1, f2 (the true π/2 at the top), the middle, f1 + 10⁻¹², f2 − 10⁻¹²};
  - ε ∈ {10⁻⁵⁰, 10⁻³⁰, …, 0.0049} × θ ∈ {0, 10⁻³⁰⁰, 10⁻¹⁰⁰, …, 0.0122718, 0.0122719, 0.02} (the NaN box and its split at π/256);
  - θ = π/2 − {10⁻³⁰⁰, …, 10⁻³} and the true π/2.

  Result: 1493 η and 1 619 161 points, 0 defects, max E = −6.25e-152, inner true − certified ≤ −3.05e-51, outer max −3.84e-5.
- **My joint check joint122.py** (each η-box's certified segments recorded by wrapping `seg_bound`):
  - J1: in every accepted η-box the certified segments tile the connector [δ₀ + i lo, δ₀ + i hi] and every polyline segment exactly; y_c ∈ [lo, hi]; the polyline starts at δ₀ + iy_c; its end lies within END_R in balls. 0 failures in 76 boxes.
  - J2: Y(η) ⊂ [lo, hi] is proved in Arb on 19 456 sub-boxes. I used my own real formula κ = ω²sin 2θ/(1 − ω²cos 2θ) and the monotonicity of min, with [0, θ_hi] where the denominator's ball meets 0, and the true π/2 at the top. 0 failures.
  - J3: the accepted η-boxes have pairwise disjoint interiors, lie in the sector and fill it exactly (area_rel 0).
- **Joints and edges, by reading.**
  - b → J (λ = 0 at z = t) → F (X = t, u = 4 inside A-FC's u ∈ [0, 4]) → C (`branchesA` certifies both branches when undecided; F holds iff α > b̂) → δ₀ + iY(η) → the connector → the polyline → the end piece. The end bound holds for every R ≤ 2.
  - Every cap point has X ≥ θ tan θ ≥ θ², so q ≤ 1 and form A covers the cap.
  - Form B's skipped boxes are 1242 + 634 with ε > ε₀ and 77 with X > δ₀, the latter covered by O (W, a classification run).
  - The edges 1.995f/ε₀ and the true π/2 are as round 121 found.
  - Route (2)'s hypotheses: t = ε/4 < ε/(2 − ε); the chord lies in the disc; Re(φ″u²) = 8Re A/|η|² < 0; on the real segment E ≤ log(1 − x²) + x².
  - The reflections leave E and the z-path unchanged.

**(4) The changed claims, against the md, the raw and the mathematics**
- **(iii)(b) uniformly on compact K ⊂ {|η| < 2}.** On K, ε ≥ ε_K > 0, so X ≥ t ≥ ε_K/4 along the inner path. The per-box bounds then give E ≤ X²·bd ('v' route) or E ≤ X³·sup f′ ('d' route) uniformly; the outer bound is ≤ −1.5e-7 and the end piece −1.33. The chord is uniform by Laplace with |φ″| bounded below on K, and finitely many boxes are involved. So the margin is uniform, and the band argument (h(r₀) below the margin, round 108) applies. Correct [own].
- **(iv) for every r < 1.** c_need² = 4 ln u/(u − 1) = 4r ln(1/r)/(1 − r) < 4 ⇔ ln x < x − 1 at x = 1/r > 1. The mass needs |η|² ≈ 4r < 4. The tail bound uses |ε̃_K|²e^{−Ax²−By²} ≤ 4e^{−(u−1)|ξ|²/2}, which I re-derived (B − 1 = γy − ½). Every frame has r < 1: γx > ½, and γy > ½ ⇔ ρ < 1. At the half-squeeze rule γx = γy = 1/(2t), so r = t < 1 for every σ*/σ. Correct [own].
- **"Conjecture P proved, computer-assisted".** Justified, given D (|η| ≤ 1.8927), V (to 1.995) and VA, all in ball arithmetic, plus the earlier rounds' analytic steps:
  - (a);
  - the Λ⁺ contour bound with ρ₀ → 1 and the band r < |σ| < 1;
  - Duhamel and Laplace on the chord;
  - route (2)'s 0-side and chord;
  - (v)'s frames.

  These are stated in D61 and were not re-refereed here.
- **Wording issues:** D3dq-4.
- **D‴.58.** D61 (iv) and the Status claim D‴.58 (i)'s sharp constant "for every σ*/σ", but D58.md at this commit still states the 0.9900 tier (D3dq-3).

**(5) Tests and mutants** [own]
- **Harness.** W/scripts/mut122.py runs each mutant in a slot copy of W/repo, then:
  - (a) `pytest -x` on test_gh_peano.py;
  - (b) `part_VA()` against the committed raw;
  - (c) the md rendered with the recomputed VA, against the committed md.
- **Spec.** W/scripts/spec122.py: the control Q00 and 45 mutants of the new or changed code. Results are in W/mut/res_spec122_{a,b}.jsonl.
- **Killed: 21.**
  - Some fire the checks themselves: H04 (handoff 216 failures), C06 (1) and P01 (soundness of B with fA, 566+ failures).
  - The rest fall to the raw pin: H01, H02, H03, H05, H10, O02, O04, O05, C07, C08, K01, K03, S04, S07, V02, V03, P02.
  - W04 (the md literal) falls to the md-from-raw test.
- **Survivors: 24, none changing the raw or the md.**
  - Equivalent on the data: H06 (assert), H09, O01, O03 (the old `hi > lo` guard), G01 (the guard is never reached, and NaN fails anyway), V01 (float width), B01, W02 and W05 (redundant: `ok` already covers the joint's and soundness' `failed`), W03.
  - Unsound in principle, equal on the data: H07, H08.
  - Weakened checks: C01 (top at the float π/2), C02, C03 (tiling tolerance 10⁻³), C04 (Re A unchecked), C05, S01 (midpoint), S02 and S06/K02 (route 'd' only against v), S03 (tolerance 0.5), S05, W01 ("NOT exactly" never printed).
- **Shipped spec round121_D61_own.py, through my harness** (W/mut/res_spec122b_c.jsonl):
  - L01–L08 are killed and L09 survives, as REVIEW_LOG says;
  - T01 (the tiling fooled) is killed only by the pin; T02 is equal on the data, since no R-split occurs.
- **The shipped tools/mutants.py with `--exp gh_lastann`** reports `"md": true` for the control L00 and for L09 alike. `gh_lastann.py --from-raw` writes no md, so that column is vacuous (W/mut/shipped/res_a.jsonl).
- **Which killed mutants change a published result** (i.e. the md):
  - md: H04, H10, O02, O05, C06, C08, K01, K03, S04, V02, V03, P01, P02, W04, L02, L06, L07, L08, T01;
  - the raw only: H01, H02, H03, H05, O04, C07, S07.

**(6) Round 121's and round 120's issues.** See the section below.

**(7) Literature.** Web search; arxiv.org was not fetchable from this session, so I relied on search snippets. Details are under Sources.

## Verdicts per part

| part | verdict |
|---|---|
| (i), (ii), (iii)(a), (iii)(c), (v) | Unchanged; not re-refereed. |
| (iii)(b) on \|η\| < 2, uniformly on compact subsets | **Proved, computer-assisted:** D + V + VA in ball arithmetic, plus earlier rounds' analytic steps. VA as adopted is a proof on 1.995 ≤ \|η\| < 2 [own: ind, ind_corners, J1–J3, cmp_fix]. Only (iii)'s header is stale (D3dq-4). |
| (iv) for every r < 1 | **Correct.** The c² argument checks; every frame has r < 1; the half-squeeze constant holds at every σ*/σ. |
| Status, summary | Correct in substance. Wording issues (D3dq-4); D58.md is behind (D3dq-3). |
| Numerics (VA), md section VA | **Every number matches the raw.** The "checks" overstate what `handoff_check`, `soundness` and the tiling check exercise (D3dq-1, D3dq-2, D3dq-4). |
| Tests | 42 pass. Strong pins (raw equality) and good plants (NaN connector, short connector, missing box, a lowered bound). No plant for an under-certified joint or a fake tiling; the series test uses `overlaps`. |

## Status of the round-121 issues

**D3dp-1 (the NaN handoff): resolved.**
- The patch is verbatim, and the connector is always certified.
- O gives 76 η-boxes and 5193 segments. The NaN box is split at π/256 with connector [0, 0.0122718463…].
- ind.py shows 0 defects, and J1 and J2 pass [own].

**D3dp-2 (the joint unchecked): partly resolved.**
- `handoff_check` exists and is printed, and plants exist (L02, L03, L05–L07 killed).
- It does not test the certified segments, y_c or the polyline's start, and its tiling is an area sum (D3dq-1, D3dq-2).

**D3dp-3 (the vacuous controls): resolved** for the repository's two controls, which now run beside their positive twins. Round 121's m ∈ {0, 0.5, 0.9, 1.1, 2} caps are not ported (minor).

**D3dp-4 (under-asserting checks): partly resolved.**
- Done: B0 asserted; true ≤ bound asserted at samples; series against Arb.
- Deferred, as logged: difference quotients against the dual enclosure, `I.__mul__`, `D.pos`, and `mvbound` on synthetic f.
- Round 121's survivors C04, C07, G01, G03, G04 and G08 still rest on reading.

**D3dp-5 (floats): mostly resolved.**
- Done: the arb width, the ∞ guard, and θ_hi from the true π/2.
- Not done: `validate` at ε₀ = 0.11 still crashes in the CLI [own]; the 1.8900000000000001 edge is unstated. Both are logged as not run.

**D3dp-6 (text): resolved for D61, the md, HANDOFF and PAPERS; not for D58.md,** although the log lists it (D3dq-3).

**Round 120's residue.**
- D3do-2, D3do-4 and D3do-5 are fixed.
- D3do-3: VA at ε₀ = 0.11 is still not run in the repository. I ran the pieces, and they pass.
- D3do-6 is fixed except:
  - D61's Numerics (VC) still lacks rounds 115 and 117's numbers;
  - VC_MID still has no point on the cap's far piece: (3, 5, 3) has P < 0.

  Both are still deferred.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dq-1 | minor | **The new checks do not check the joint they report.** `handoff_check` trusts `outer_box`'s returned (lo, hi). It never checks that the certified segments cover [lo, hi], that y_c ∈ [lo, hi], or that the polyline starts at δ₀ + iy_c. Unsound O02/O05 (connector partly certified) and H10 (polyline 10⁻³ above the connector) pass with 0 outside. `soundness` passes V02 (mean value width 0). All four are killed only by `va == raw`. The adopted code is right by construction: J1 shows 0 failures in 76 boxes, and J2 proves containment on 19 456 sub-boxes [own: mut/res_spec122_*.jsonl, scripts/joint122.py]. | gh_lastann.py: have `outer_box` return its certified segments (or their union per original segment). In `handoff_check`, assert that they tile [lo, hi] and each polyline segment, that y_c ∈ [lo, hi] and that pts[0] = δ₀ + iy_c. Add J2's rigorous sub-box containment. test_gh_peano.py: plants O02, H10 and V02 that must make the checks fail, not just the pin. |
| D3dq-2 | minor | **The tiling check is an area sum.** It is redundant with O's `failed == 0`, and T01 (each lower θ-half duplicated in place of the upper) leaves `tiled` True with the area exact [own: res_spec122b_c.jsonl]. `va_lines`'s conclusion also ignores `tiled` (W01 survives). | `handoff_check`: also check pairwise disjoint interiors and containment in the sector (as J3), which with the area makes a real tiling check. `va_lines`: require `tiled` for the conclusion. Plant T01. |
| D3dq-3 | minor | **D‴.58 is behind D61.** D61 (iv) and the Status, HANDOFF l. 34 and PAPERS claim D‴.58 (i)'s constant for every σ*/σ and t(v) for every rule with a frame. D58.md at 04d9ead still says t < 0.9900 (σ*/σ < 199), "conjectural beyond" (ll. 14, 25, 171–192, 287, 310). REVIEW_LOG's D3dp-6 row lists "D‴.58's tiers" as applied; HANDOFF l. 184 defers it to round 123 [md]. | Either update D58.md's tiers, or mark them deferred in REVIEW_LOG (Round 121, D3dp-6) and keep HANDOFF and PAPERS' D‴.58 sentences conditional until round 123. |
| D3dq-4 | cosmetic | **Wording in D61 and the md:** | D61 (iii) header, Status, Numerics (VA); `va_lines`. |
| | | • (iii)'s header still says "proved on most of the disc". | "proved on the open disc, computer-assisted beyond \|η\| < 1 and the axes". |
| | | • "The true value below the certified bound" is fA/fB at the point (forms, one-sided where they simplify), not F_η's value. | "the forms' value". |
| | | • "(1319 outside the annulus)": 77 of them are skipped for X > δ₀. | "outside the inner pieces' domain". |
| | | • "The joint checked at every η-box" means 15 sampled η per box. | Say so, and note that the enclosure itself is by construction. |
| | | • "Arb's correctness assumed" omits python-flint and CPython's float/`math.nextafter`. | Name them as assumed too. |
| | | • "Two exact scaled forms" are exact up to (Re σ − 1)₊ and the dropped −P₊² (upper bounds, exact on the path). | Say so. |
| D3dq-5 | cosmetic | **Code and tools:** | gh_lastann.py / exp_gh_peano.py / tests / tools / HANDOFF / REVIEW_LOG as listed. |
| | | • `va_lines` prints "at {128} bits", a literal. | Use LA.PREC, stored in the raw. |
| | | • The series test asserts `overlaps`. | Assert `contains`. |
| | | • tools/mutants.py with `--exp gh_lastann` reports `md: true` for the control L00 (W/mut/shipped/res_a.jsonl). | Run part VA instead, as W/scripts/va_md.py does. |
| | | • HANDOFF l. 14 says "63 + 39 pass". | 42. |
| | | • The log repeats "9 mutants and a control" and calls L09 "equivalent". | Delete the repeat; say "equal on the data". |
| | | • The CLI `all --eps0 0.11` crashes in `validate`. | Clip validate's X to δ₀. |
| D3dq-6 | minor | **Survivors in the new code: 24 of 45**, none changing the raw or the md. | |
| | | • Unsound in principle, equal on the data: H07 (NaN branch lo = f1), H08 (y_c's clamp removed). These rest on reading, and are correct as written. | J1's y_c assertion kills H08; a NaN-box plant with f1 > 0 kills H07. |
| | | • Weakened checks: C01, C02, C03, C04, C05, S01, S02, S03, S05, S06, K02. | A soundness plant lowering the bounds by 10⁻⁶ (S03); a route-'d' plant (S02, S06, K02); a check at the top box's true π/2 that fails at the float one (C01). |
| | | • Still deferred from round 121: the dual-enclosure, `I`, `D.pos` and `mvbound` unit tests. | As round 121 proposed. |

## Sources

**Repository** [md][raw]: D61.md, gh_lastann.py, exp_gh_peano.py, gh_peano.md and its raw, test_gh_peano.py, REVIEW_round120_D3do.md, REVIEW_round121_D3dp.md, REVIEW_LOG.md (Rounds 120–121), D58.md, HANDOFF.md, PAPERS.md, tools/mutants*, all at 04d9ead. Also ref121/fix/lastann.py and ref121/att/ind.py.

**My scripts and logs** [own], all under W:
- scripts: va_recompute.py, scripts/ind.py (round 121's, unchanged), scripts/ind_corners.py, scripts/joint122.py, scripts/cmp_fix.py, scripts/mut122.py, scripts/va_md.py, scripts/spec122.py, scripts/spec122b.py;
- logs: logs/{pytest_gh, va_recompute, ind, ind_corners, mut_a, mut_b, mut_c}.log;
- mutant results: mut/res_spec122_a.jsonl, mut/res_spec122_b.jsonl, mut/res_spec122b_c.jsonl, mut/shipped/res_a.jsonl;
- CLI runs: cli/eps011*.log.

**Literature, web search** [snippet]:
- Johansson, "Arb: efficient arbitrary-precision midpoint-radius interval arithmetic", IEEE Trans. Comput. 66(8), 1281–1292 (2017), DOI 10.1109/TC.2017.2690633, arXiv:1611.02831. This confirms D61's citation; ARITH 24 was the talk.
- Rump, "Verification methods: rigorous results using floating-point arithmetic", Acta Numerica 19 (2010), 287–449.
- Gómez-Serrano, "Computer-assisted proofs in PDE: a survey", arXiv:1810.00745 (2018); interval arithmetic as the standard tool.
- IEEE Std 1788-2015 (drafts; GNU Octave interval package documentation): NaI ("Not an Interval") is distinct from the empty set. This is the standard's form of round 121's lesson that a NaN is not an empty interval.
- de Boor, "Divided differences", arXiv:math/0502036: B-splines as the Peano kernels of divided differences; Hermite–Genocchi.
- arXiv:2312.07940: root-exponential convergence of Gauss–Hermite for analytic integrands; it cites Donaldson–Elliott's contour representation of the remainder.
- Saad (1992), SIAM J. Numer. Anal.: an a posteriori Krylov error formula (Thm 5.1), seen through secondary sources (Jia–Lv, arXiv:1307.7219, and a PMC article), not the original.
- Bui-Thanh–Villatoro–Krishnanunni, arXiv:2602.18634 (2026): error kernels of two-point Hermite quadrature are Legendre polynomials. This is not Gauss–Hermite and does not overlap D61 (i).
- Pei Liu, arXiv:2509.23532: asymptotic errors of Gauss–Legendre.
- A Numer. Math. paper on the asymptotics of Peano kernels for Filippi, Pólya and Clenshaw–Curtis rules (DOI 10.1007/BF01400357) [snippet: record only].
- I found no source stating the Gauss–Hermite Peano kernel as a law, or the limit of its Laplace transform. D61's "new, as far as these sources show" survives my search, which is not a proof of novelty.

**From memory** [memory]: Tucker, *Validated Numerics* (2011); Moore, *Interval Analysis* (1966).

## Lesson

1. **A check must look at what was certified, not at what was reported.** `handoff_check` verifies the interval that `outer_box` returns, so an `outer_box` that certifies less than it returns passes. The joint check should rebuild the cover from the certified pieces, and an area sum is not a tiling without disjointness.
2. **When a check is introduced, plant a mutant it must catch, not only one the pins catch.** Here every unsound joint was stopped by the exact-raw pin while the new checks printed "0 outside". That is good regression protection, but it does not show the checks are what the md says.
3. **Upgrade every document a tier lives in, in the same commit or with an explicit deferral.** D61, HANDOFF and PAPERS moved to "every σ*/σ" while D58.md, the proposition the constant belongs to, still says σ*/σ < 199.
