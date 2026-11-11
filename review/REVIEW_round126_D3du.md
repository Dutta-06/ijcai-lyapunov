# Referee report D3du, round 126: Proposition D‴.61 as revised after round 122

## Scope

**What I reviewed.** Commit 5f56b63 on claude/vigilant-clarke-hf939b. I exported it with `git archive 5f56b63 | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref126. W/repo/venv is a symlink to the repository's venv (Python 3.13.18, python-flint 0.9.0).

**Hashes (sha256 prefixes) [own]:**
- gh_lastann.py 359f668beb469ae6
- exp_gh_peano.py e649ed303404d3c0
- test_gh_peano.py a5c4c0ef62ba40a5
- D61.md 004b37b6778cd90e
- gh_peano.md 8ed172468bfa5c8f
- gh_peano_raw.npz 4c4c4a4d314e8fbf
- round122_D61_own.py 0c9feb674a24f3ee

**How I ran things.**
- Read-only git only. The repository moved during the round (HEAD is now 18d5b41, with other sessions' edits to D58, width_mismatch and noise_resonance). I reviewed none of that.
- Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, with GP_OUT inside W and at most two heavy processes.
- I killed one process: PID 22348, my own `topedge.py y45` run, after checking /proc/22348/cmdline. Every other process of mine exited on its own.

**Read in full:**
- D61.md;
- gh_lastann.py;
- `part_VA` and `va_lines`;
- gh_peano.md section VA;
- the VA test and the series test;
- round122_D61_own.py;
- REVIEW_round122_D3dq.md and the Round 122 log entry;
- `git diff 04d9ead 5f56b63` for these files.

I also looked at the parts of HANDOFF and D58.md that bear on round 122's issues.

## Summary

**Design.** The revision does the following:
- `outer_box` returns the pieces it certified, in order, together with lo, hi, y_c and the polyline.
- `handoff_check` splits the pieces into a connector prefix and a polyline.
  - `_chain` must take the connector from δ₀+i·lo to δ₀+i·hi, and the polyline from its start to its end through every vertex.
  - y_c must lie in [lo, hi], and the polyline must start at δ₀+i·y_c.
- `_y_enclosed` proves Y(η) = min(κδ₀, arg η) ∈ [lo, hi] on adaptive sub-boxes in Arb.
- The tiling is now containment plus pairwise disjoint interiors plus the area to 10⁻¹².
- `soundness` adds four corners at the expansion variable's upper end.
- `va_lines`' conclusion requires `tiled` and `joints`.
- A 'meta' record stores the precision.
- The series test now uses `contains`.
- `validate` skips (and counts) points outside the forms' domain.
- D61 and the md are reworded.

**Verdicts.**
1. **Reproducibility is exact** [own].
   - `--from-raw` regenerates the md byte for byte (sha 8ed17246…).
   - 42/42 tests pass (40.5 s).
   - `part_VA()` recomputes equal to the raw's VA (6.9 s).
   - Every VA number in D61's Numerics and in the md is in the raw: 484 sub-boxes, 1140 η, 240/16074/1920/4614 points, worst −2.4253e-13, O 76/5193 with worst −1.4993e-7, end −1.3314, 7.6e-23 at 1479 points, 2.1e-59 [raw].
2. **`_y_enclosed` is a rigorous proof of what it claims**, given that the box lies in the sector (R ≤ 2, 0 ≤ θ ≤ π/2), which `inside` checks.
   - Every clamp and exact endpoint is justified (details in "What I checked").
   - Independently, at 200 bits on a 41×41 grid per box, the true Y range lies in [lo, hi] for all 76 boxes. The lower margin is exactly 0 on bottom edges where lo = f1 = Y (tight, and handled by θ's exact ends); the upper margin is at least 5.3e-12.
   - 484 sub-boxes, depth ≤ 6, one fallback sub-box (the corner η = 2) [own: joint126].
3. **`_chain` and the split cannot be fooled into accepting an unsound joint.**
   - My plants of reordered pieces, a connector after the polyline, a connector overshooting hi, lo = hi = y_c, the connector's last piece missing, the whole connector dropped (round 120's defect) and a dropped middle half-piece are all rejected.
   - A detour off the polyline through two certified pieces is accepted. That is sound: a chain of certified segments is a path.
   - The joint's last link, the polyline's end inside |1 − z| < END_R where the end piece takes over, is **not** rebuilt (D3du-2).
4. **The tiling check does not imply a tiling with float boxes.**
   - A 10⁻¹⁴ gap passes with `tiled` True.
   - The tiling does hold: the exact rational area equals the sector's (D3du-1).
5. **Soundness corners.** They sit at the upper end of the expansion variable as intended.
   - In the junction pieces u is degenerate, so the 4 corners collapse to 2 (A-J 240 → 160 distinct, B-J 1920 → 1280), and λ's upper end is never taken.
   - The lower end, where the excess is closer to 0 (−1.5e-13), is never sampled (D3du-5).
6. **Wording** is accurate on "proved" against "sampled". The one overreach: D61's Status says "the η-boxes tiling the sector", which the check does not establish (D3du-1). HANDOFF §5 is stale (D3du-6).
7. **Mutants.**
   - Of 58 new single-point mutants, 24 are killed and 34 survive; no survivor changes the raw or the md.
   - 20 survivors are unsound in principle: they leave the shipped data unchanged but would let a wrong joint, enclosure or tiling through.
   - With the raw pin removed, the unsound `_y_enclosed` mutants Y09, Y12, Y14 and Y20 pass all tests with "0 failures", so only the pinned count 484 stops them.
   - The shipped P01–P09 are all killed. Round 121's C01 still survives (D3du-3, D3du-4).

**No major issue. The certificate as shipped is sound.**

**Findings.**
1. **D3du-1 (minor).** The tiling check (area to 10⁻¹²) passes a 10⁻¹⁴ θ-gap. The tiling holds exactly, but by construction and by my rational check, not by the shipped check.
2. **D3du-2 (minor).** `handoff_check` does not rebuild the last link. A polyline cut halfway (0.53 from z = 1) passes with 0 failures, and O01 (outer_box's END_R test relaxed to < 1) survives everything.
3. **D3du-3 (minor).** The new checks' failure paths are under-planted. Q02, Q03, Q04, Q05, Y06, Y08, Y19, H03, H07, T01–T05, T07, T10 and W01 survive. Y09, Y12, Y14 and Y20 are stopped only by the pinned sub-box count.
4. **D3du-4 (minor).** The true π/2 at the top edge is load-bearing but untested. Y04 makes `_y_enclosed` "prove" lo = 10⁻²⁰ on the top box, although Y = 0 at the true π/2. Y04, Y05 and C01 survive; the log says C01 is killed.
5. **D3du-5 (cosmetic).** The soundness corners are duplicated in J pieces and miss λ's upper end and the expansion variable's lower end. The printed counts include duplicates.
6. **D3du-6 (cosmetic).** Three small points:
   - `_y_enclosed` does not assert its preconditions; on a box to |η| = 2.001 it returns a "proof".
   - The docstring says "covering the polyline" where the check is a chain through its vertices.
   - The adaptive bisection has no budget.
   - Also: HANDOFF §5 still says "two exact scaled forms", "joint checked at 1140 η" and "7616 sampled points".

## What I checked and how

**(1) Reproducibility** [own]
- **The md.** `exp_gh_peano.py --from-raw` with GP_OUT=W/out, then `cmp` against the committed md: identical.
- **The tests.** `pytest -q -p no:cacheprovider --basetemp W/pt/base1 tests/test_gh_peano.py`: 42 passed in 40.48 s (W/logs/pytest_gh.log).
- **Part VA.** W/scripts/va_recompute.py: equal, in 6.94 s (W/logs/va_recompute.log).
- **The CLI at ε₀ = 0.11.** `gh_lastann.py validate --eps0 0.11` now runs (2 points skipped). Its B_vs_centred is 1.1e-12, against 7.6e-23 at ε₀ = 0.005. My vB.py shows every A and B deviation lies within the balls' radii (ratio ≤ 0.973) with no inconsistency: it is series truncation carried in the radius, not a defect.

**(2) `_y_enclosed` and the joint** [own, reading]
- **The formula.** κ = Im A/|Re A| with A = ω²e^{2iθ} − 1. Then |Re A| = 1 − ω²cos 2θ = (1 − ω²) + 2ω²sin²θ, which is correct and free of cancellation.
- **The clamps.**
  - On the true box, R ≤ 2 gives 1 − ω² ≥ 0, and 0 ≤ θ ≤ π/2 gives sin 2θ ≥ 0 and κ ≥ 0.
  - The balls do overshoot the box: `hull(0, x)` has a lower end of about −2⁻³⁰·x (Arb's 30-bit radius) [own].
  - Applying `pos` to an enclosure of a quantity whose true values are ≥ 0 still encloses those values. So the clamps on 1 − ω², sin 2θ and k_lo are valid.
  - They are valid only under R2 ≤ 2 and 0 ≤ f1 ≤ f2 ≤ π/2. `_y_enclosed` does not assert this; `inside` does (T3 plant: a box to 2.001 gives `_y_enclosed` → 4 and `inside` False, so `tiled` is False and there is no conclusion).
- **The bounds on Y.**
  - Y(η) = min(κδ₀, θ) ≥ min(inf κδ₀, θ₁) ≥ min(k_lo, θ₁), and Y ≤ min(k_hi, θ₂). Both hold by the monotonicity of min.
  - `lower()`/`upper()` return exact points (radius 0, checked), so the comparisons are exact.
- **The fallback.** Where `den > 0` fails, [0, θ₂] is valid because 0 ≤ Y ≤ θ. It forces lo = 0, which the corner box has.
  - A NaN K also falls through to [0, t_hi], which is also valid.
  - At η = 2 itself κ is undefined; that point is outside the open disc.
- **The top edge.** It is taken at the true π/2 in both Tb and t_hi, and upper sub-boxes inherit t2 = math.pi/2 exactly.
- **Measured** (W/scripts/joint126.py):
  - 484 sub-boxes; depth histogram {0: 19, 1: 171, 2: 214, 3: 51, 4: 18, 5: 7, 6: 4}; 1 fallback sub-box.
  - The true Y range lies within [lo, hi] on all 76 boxes. The top box has lo = 0 and Y ∈ [2.8e-63, 0.00123].
- **The certified pieces against the polyline (J1).**
  - 133 connector pieces and 5060 polyline pieces; every polyline segment is reached; the end lies at most 0.00995 from 1.
  - 403 pieces are off the exact polyline by float midpoint rounding (cross product ≤ 5.8e-19). This is harmless, since the proof needs a chain of certified segments, not the polyline itself.
- **Plants against `handoff_check`** (joint126, plants2):

  | plant | result |
  |---|---|
  | detour off the polyline (both pieces certified, −5.7e-3 and −6.6e-3) | accepted (sound) |
  | swapped pieces | 'polyline cover' |
  | connector after the polyline | rejected |
  | lo = hi = y_c | 'Y not enclosed' plus samples |
  | lo → y_c → hi with a shared vertical piece | rejected (conservatively) |
  | connector overshooting hi | rejected |
  | connector's last piece missing | 'connector cover' |
  | whole connector dropped | 'connector cover' |
  | a middle half-piece dropped | 'polyline cover' |
  | **P7** (one vertex trimmed, end 0.0209 > END_R from 1) | **0 failures** |
  | **P8** (cut halfway, 0.53 from 1) | **0 failures** |

**(3) The tiling**
- On the data, `inside` and `disjoint` are exact float comparisons. Plant T1 (one box's f1 raised by 10⁻¹⁴) gives inside True, disjoint True and area_rel 6.4e-15, so `tiled` is True with an uncovered sliver.
- The exact check closes this: in rationals, the sum of (b − a)(f2 − f1) equals (2 − 1.995f)·(π/2)f exactly (J3). Containment plus disjoint interiors plus equal measure, with a closed union, then gives union = rectangle.
- The top boxes extend to the true π/2 in `seg_bound` and `_y_enclosed`.
- V's band ends at the float 1.995 and VA's R1 is the same float, so there is no radial gap.

**(4) Soundness corners** (joint126, bj_corners.py; mutants S01 and S02)
- **Distinct points.** A-J 240 → 160 distinct and B-J 1920 → 1280, because u = (UJ, UJ) makes the u-bit idle. A-FC and B-F are all distinct. No J corner has λ at its upper end.
- **λ's upper end (S01).** It changes only B-J's worst, from −0.0107 to −0.0592, so it is not the tight side.
- **The lower end (S02).** Corners at the lower end give A-FC's largest excess as −1.5e-13 instead of −2.4e-13.
- **Dense grid.** My dense grid (3×3 × 5 points including both ends) on the shipped A-FC gives −6.06e-14.
- **V02 (the mean value width dropped).** It is still not flagged by `soundness`. On V02's own accepted boxes, a dense check finds no violation either: A-FC 108 720 points, worst −6.1e-14; B-F 31 239; B-J 8334. So V02 is unsound in principle but true at every sampled point, and no sampling check can catch it; the pins (box and point counts) do.

**(5) Wording** [md]: see D3du-1 and D3du-6. The rest of what D61 and the md say matches the code:
- "proved inside it in balls on 484 sub-boxes (and sampled at 1140 η)";
- "chaining the connector … the polyline through every vertex";
- "the forms' value", "outside its domain";
- the assumed correctness of Arb, python-flint and CPython's float arithmetic and `math.nextafter`.

**(6) Mutants** [own]
- **Harness.** W/scripts/mut126.py. Each mutant runs in a slot copy and is judged three ways:
  - (a) `pytest -x` on `-k "part_va or md_is_the_report or report_and_main"`;
  - (b) `part_VA()` against the committed raw;
  - (c) the md rendered with the recomputed VA against the committed md.
- **Other runs.** mut126b.py records the failing test line. mut126c.py reruns with the raw pin `assert va == list(raw['VA'])` removed.
- **Spec.** W/scripts/spec126.py: the control Q00 and 58 mutants. W/scripts/spec126b.py: C01 and the shipped P00–P09.
- **Killed, 24.**
  - By plants (raw unchanged): H05, H06, H08, H09, Q01, T06, Y07, Y11.
  - With the raw changed: O02, S01, S02, T09, V02a, V02b, W03, Y02, Y09, Y10, Y12–Y16, Y20.
  - **With the pin removed**, Y09, Y12, Y14, Y20, S01 and S02 survive. Y13 and Y15 fall to the short-connector plant (l. 1292). Y02, Y10, Y16 and O02 fall to `failed == 0`, T09 to `tiled`, V02a/b to the soundness point counts, and W03 to `'joints' not in O`.
- **Survivors, 34; none changes the raw or the md.**
  - **Unsound in principle, equal on the data (20):**
    - Q02 (the connector chain's end unchecked: the O05 shape passes);
    - Q03 and Q04 (links unchecked);
    - Q05 (an empty chain passes: round 120's dropped connector);
    - H03 (the polyline chain starts from its own first piece);
    - H07 (the joints flag constant);
    - O01 (END_R relaxed);
    - T01–T04 (containment weakened), T05 (area to 1e-3), T07 (`tiled` without `inside`), T10 (disjointness only between neighbours in the list);
    - W01 (the conclusion without joints);
    - Y04 and Y05 (the top edge at the float π/2);
    - Y06 (the fallback's lower end θ₁);
    - Y08 (no min with θ₁);
    - Y19 (`_y_enclosed`'s None not reported).
  - **Sound weakenings or equivalent (14):** H01, H02, H04, O03 (the y_c check now stands behind it), S03, S04, S05 (tolerance 1e-6, below the 1e-5 plant), T08 (no R-split on the data), W02, W04 (ε₀ = 0.11 only), Y01, Y03, Y17 (looser, still valid enclosures), Y18 (equivalent).
- **Shipped P01–P09:** all killed. P01–P03 and P07–P09 by plants; P04, P05 and P06 by the raw pin first. **C01 survives.**
- **Which kills change published results.**
  - The md: O02, S02, T09, V02a, V02b, Y02, Y09, Y10, Y12–Y16, Y20, P04, P05, P06. For the `_y_enclosed` mutants this is only through the printed count "484 sub-boxes".
  - The raw only: S01, W03.
- **The top-edge plant** (W/scripts/topedge.py). On the top box with lo ∈ {1e-30, 1e-20, 1e-18}:
  - the shipped code returns None, which is correct since Y = 0 at the true π/2;
  - **Y04 returns 12.7 M, 12.8 M and 25 M sub-boxes, i.e. a false "proof"**, after minutes of bisection;
  - Y05 returns None.
  - The combined Y04 + Y05 run was stopped by me; Y04 alone suffices.

## Verdicts per part

| part | verdict |
|---|---|
| (i), (ii), (iii)(a), (iii)(c), (iv), (v) | Unchanged this round; not re-refereed. Round 122's checks of (iv) for every r < 1 stand. |
| (iii)(b) on \|η\| < 2 (VA) | **Sound as shipped.** `_y_enclosed` is rigorous; the joint's pieces chain; the tiling holds exactly (rational area). The checks are weaker than stated in two places (D3du-1, D3du-2), and their failure paths are under-planted (D3du-3, D3du-4). |
| Numerics (VA), md section VA | Every number matches the raw. The "seeded points and corners" counts include duplicates (D3du-5). |
| Status | Correct, except that "the η-boxes tiling the sector" is not what the check proves (D3du-1). |
| Tests | 42 pass; strong pins; good plants for O02, H10, y_c, the shortcut, the short connector, T01, S03 and route 'd'. Missing plants: the connector's end, links, the empty chain, the polyline's end, a 1e-14 gap, the true π/2, `_y_enclosed`'s None and its lo side. |

## Status of the round-122 issues

**D3dq-1 (the joint from what was reported): resolved in substance.**
- `outer_box` returns its certified pieces; `handoff_check` rebuilds the connector and polyline chains; Y(η) is proved inside the connector in balls.
- Plants exist for O02, H10, y_c and the shortcut; the shipped P01, P02, P04, P08 and P09 are killed.
- Residue:
  - the end link is not rebuilt (D3du-2);
  - the O05 shape, links and the dropped connector are caught by the code but unplanted (Q02, Q03, Q05 survive);
  - V02 is still not flagged by `soundness`, but it is true at 100k+ dense points, so it is caught by the counts, as is inherent.

**D3dq-2 (the tiling): mostly resolved.**
- Containment and disjointness are added, and `va_lines` requires `tiled` and `joints`.
- Residue: the area is checked only to 10⁻¹² (D3du-1); T01–T05, T07, T10 and W01 survive.

**D3dq-3 (D58 behind): resolved** at 783d55b. D58.md now gives the earlier tiers as "before VA" history [md].

**D3dq-4 (wording): resolved in D61 and the md.** HANDOFF §5 is still stale (D3du-6).

**D3dq-5 (code and tools): resolved.**
- The 'meta' precision record, the series test with `contains` against a 400-bit reference, `validate` running at ε₀ = 0.11 [own], and HANDOFF "63 + 42" are all in.
- The md column is documented as vacuous for `--exp gh_lastann`.

**D3dq-6 (survivors): partly resolved.**
- H07: the equivalence argument holds (f1 = 0 is the only case on this grid), and `_y_enclosed` would now catch a wrong lo anyway.
- H08: caught (P09).
- The S03 tolerance plant at 1e-5 is in; S05 at 1e-6 survives.
- The route-'d' plant is in.
- **C01 is listed as killed but survives**, and its rigorous counterparts Y04/Y05 survive too (D3du-4).
- The unit tests for G01, G03, G04 and G08 are still deferred, as logged.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3du-1 | minor | **The tiling check does not imply a tiling.** It checks containment, disjoint interiors and \|area/sector − 1\| < 10⁻¹². T1 (a 10⁻¹⁴ θ-gap) gives `tiled` True with area_rel 6.4e-15, and T05 (tolerance 1e-3) survives. The tiling does hold: the exact rational area equals (2 − 1.995f)(π/2)f (J3) [own: joint126]. D61's Status says "the η-boxes tiling the sector" [md]. | gh_lastann.py `handoff_check`: compare `sum(Fraction(b)-Fraction(a))*(Fraction(f2)-Fraction(f1))` with `(Fraction(2.0)-Fraction(R1))*Fraction(math.pi/2)` exactly. test_gh_peano.py: plant T1. D61 Status: "inside the sector, disjoint, of exactly its area (rationals)". |
| D3du-2 | minor | **The polyline → end piece link is not rebuilt.** `handoff_check` never checks \|1 − pts[−1]\| < END_R. P7 (end 0.0209 from 1) and P8 (cut halfway, 0.53 from 1) pass with 0 failures, and O01 (outer_box's END_R test relaxed to < 1) survives tests, raw and md [own: joint126, mut]. | `handoff_check`: assert `abs(1 - acb(pts[-1])) < END_R` in balls. Plant P8. |
| D3du-3 | minor | **The new checks' failure paths are under-planted.** Unsound in principle and surviving: Q02, Q03/Q04, Q05 (an empty chain, i.e. round 120's dropped connector), H03, H07, T01–T04, T07, T10, W01, Y06, Y08, Y19. With the raw pin removed, Y09, Y12, Y14 and Y20 also pass, with `handoff_check` reporting 0 failures. The shipped code catches the O05 shape, the dropped connector and a dropped middle piece (plants2), so these are test gaps [own: res_spec126_*, res_nopin_*]. | test_gh_peano.py, on the in-memory Oj, plant: the connector's last piece dropped; the whole connector dropped; a middle polyline half-piece dropped. For `_y_enclosed`: lo raised above inf Y (lo side); a connector short only in the upper-right quarter; `_y_enclosed` monkeypatched to return None, expecting `failed > 0`. Also a box outside the sector (`inside`) and `joints=False` in `va_lines`. |
| D3du-4 | minor | **The true π/2 is load-bearing but untested.** With Y04 (Tb at the float π/2), `_y_enclosed` "proves" lo = 10⁻²⁰ on the top box (12.8 M sub-boxes), although Y = 0 at the true π/2; the shipped code returns None. Y04, Y05 and round 121's C01 survive, contrary to the round-122 log's "C01 … killed" [own: topedge, spec126b]. | test: `assert LA._y_enclosed(*top_box[:4], 1e-20, hi, DELTA0) is None`, with a sub-box budget (D3du-6), so that a regression fails fast. REVIEW_LOG: correct the C01 line. |
| D3du-5 | cosmetic | **Soundness corners.** In A-J and B-J, u is degenerate, so the 4 corners are 2 distinct points (240 → 160, 1920 → 1280) and λ's upper end is never sampled. No piece samples the expansion variable's lower end, where A-FC's excess is closer to 0 (S02: −1.5e-13; my dense grid: −6.1e-14 against the printed −2.4e-13). The printed counts include duplicates [own]. | `soundness`: corners over the non-degenerate coordinates (u or λ, and q/p) at both ends of the expansion variable. Print the distinct counts in D61 and the md. |
| D3du-6 | cosmetic | **Small code and text points:** | |
| | | • `_y_enclosed`'s clamps need R2 ≤ 2 and 0 ≤ f1 ≤ f2 ≤ π/2, which only `inside` checks (T3: it returns 4 on a box to 2.001). | Assert the preconditions in `_y_enclosed`. |
| | | • The bisection has no budget (25 M sub-boxes in a Y04 variant). | Add a sub-box budget that returns None. |
| | | • The `handoff_check` docstring says "covering … the polyline"; the check is a chain through its vertices (detours accepted, which is sound). | Say "a chain of certified pieces through its vertices". |
| | | • HANDOFF §5 still says "two exact scaled forms", "the joint checked at 1140 η" and "the true value below the bound at 7616 sampled points" [md]. | Align it with D61's round-122 wording and counts. |

## Sources

**Repository** [md][raw] at 5f56b63: D61.md, gh_lastann.py, exp_gh_peano.py (`part_VA`, `va_lines`), gh_peano.md and gh_peano_raw.npz, test_gh_peano.py, tools/mutants/round122_D61_own.py, REVIEW_round122_D3dq.md, REVIEW_LOG.md (Round 122), HANDOFF.md §5 and D58.md (for D3dq-3 and D3dq-4).

**My scripts and logs** [own], all under W:
- scripts/va_recompute.py, scripts/joint126.py (J1–J3, Y's true range, the plants, the corners), scripts/plants2.py, scripts/vB.py, scripts/dense_v02.py, scripts/topedge.py, scripts/bj_corners.py;
- the harness: scripts/mut126.py, scripts/mut126b.py, scripts/mut126c.py (raw pin removed), scripts/evalva.py;
- the specs: scripts/spec126.py, scripts/spec126b.py;
- logs: logs/{pytest_gh.log, va_recompute.log, joint126.log, joint126.json, cli_validate011.log, dense_v02a.log, dense_v02b.log, mut_a–f.log};
- mutant results: mut/res_spec126_{a,b,c}.jsonl, mut/res_spec126b_d.jsonl, mut/res_nopin_spec126_{e,f}.jsonl.

**Literature** [memory]:
- Moore, *Interval Analysis* (1966);
- Tucker, *Validated Numerics* (2011);
- Rump, "Verification methods", Acta Numerica 19 (2010): enclosures are used only through monotone operations on valid bounds, and a clamp to a known sign is valid on an enclosure of a quantity of that sign;
- IEEE 1788-2015: NaI is distinct from the empty interval.

No new web search this round.

## Lesson

1. **When the claim is exact, the check must be exact.** A tiling claimed of float boxes is decidable exactly, in rationals. A 10⁻¹² tolerance turns "tiled" into "nearly tiled", and a 10⁻¹⁴ gap passes.
2. **Rebuild every link of a joint, not just the new ones.** The connector and polyline are now rebuilt from what was certified, but the hand-off to the end piece still rests on `outer_box`'s word.
3. **Plant the failure path of a rigorous sub-proof as well as its success.** Six unsound variants of `_y_enclosed` printed "0 failures". Four of them were stopped only because the md happens to print the sub-box count, and the variant that drops the true π/2 "proved" a false enclosure on 12.8 M sub-boxes.
