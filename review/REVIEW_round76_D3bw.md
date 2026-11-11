# Referee report D3bw, round 76: Proposition D‴.60 as revised after round 69

## Scope

- **Commit:** eeef69c, exported with `git archive eeef69c | tar -x` to W/repo. W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref76/D3bw_work (repo, out, scripts, mut, logs, build).
  - HEAD has since moved to 6764e41. `git diff eeef69c HEAD` is empty over every file reviewed.
- **Files read:**
  - review/drafts/D60.md (222 lines);
  - experiments/exp_noise_resonance.py (1723 lines), and what it uses from exp_closed_transit.py, exp_dead_zone.py (`fitted_state`, `noise_of`), exp_cohesion_leak.py (`load`/`save`) and src/core/closed_loop.py/.c;
  - results/chaos/noise_resonance.md (494 lines) and noise_resonance_raw.npz;
  - tests/test_noise_resonance.py (47 tests);
  - review/REVIEW_round69_D3bp.md, and the "Round 69" entry of review/REVIEW_LOG.md (lines 2881–2921).
- **I modified nothing in the repository** and ran no git command that writes. The C extension was built into W/build (`LCK_BUILD`), `NR_OUT` pointed into W/out, and the mutants ran in private trees under W/mut.
- **Run settings:** `nice -n 10`, OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1, pytest `-p no:cacheprovider`. At most two heavy processes of mine ran at a time, and none is left running.
- **Marks:** every number below is marked (md), (raw) or (mine). Web search confirmed two sources (Anderson 1985; Kokotović–Riedle–Praly 1985); the rest are [memory].

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (md5 61f003dab65d545068d56c5df53c3bc7).
  - The 47 tests pass in 95 s.
  - Twelve cells rerun with the script's own job functions equal the raw in every field:
    - `c7_job` ×4: s = 0; 21.3 rad/s at 0.02; draw 7's window from period 125 at 0.02; at k_e = 20, the window from period 50 at 0.08;
    - `edge_job` ×3: kernels 1 and 3; WFR;
    - `c3_job` ×2: C6's draw 3 raised by 0.03; C4's −n twin of draw 7;
    - `c_job` ×2 and `w2_job` ×1.
  - The raw's C records predate the fields `c_end`, `w_left` and `pre`. Every field they do hold is identical.
- **Every number in D60.md is printed by the md.** Of 508 distinct numbers, 22 are not printed verbatim. They are bracket midpoints (0.0744, the C7 single-tone edges), roundings (1150, 1540, 210, 310, −2.6·10⁻⁴, 1.085), 146.5 = 29.3/0.2, and years.
- **(i) is right.**
  - The sensitivity |S|, |S|² − 1 = K(2ω² − K)/|K − ω² + ik_eω|² and |S(√K)| all re-derive. So do |S| = 1.414, 1.385, 1.240, 1.110 at the bound (mine).
  - The plant's slope does cancel in the ratio.
  - Bode's integral applies:
    - L = K/(s(s + k_e)) has relative degree 2;
    - it has no pole in the open right half-plane (the integrator sits on the axis);
    - the closed loop s² + k_es + K is stable.
  - Freezing K is the quasi-static approximation, and A measures its error: 0.999–1.082 at ω ≥ 10 (md).
- **(ii) is right.**
  - The Floquet algebra re-derives; one Re is missing (D3bw-10).
  - B5 at s = 0.00125 shows the averaged loop is the constant part's response at O(s²): 0.9979–1.0054 (md).
- **C4 and W2 are like for like.** The rise's excess over s²⟨w₂⟩ is steady over the ten periods before the leave and falls further back (mine). So the excess is a growing higher order, not a first-passage selection at the leave.

  | periods before the leave (s = 0.02) | 0 | 5 | 10 | 80 |
  |---|---|---|---|---|
  | median of rise / s²⟨w₂⟩ | 1.25 | 1.27 | 1.26 | 1.08 |

- **The noise-free edge is real.**
  - On a fine grid it is monotone: the kernel returns from raises ≤ 0.0730 and is ejected from 0.0735 (mine).
  - Its trends with m, k_e and the metric explain the first jumps.
  - s* is calibrated: scaled by 0.9 or 1.1, the consistent count drops from 33 to 28 or 31 (mine).
- **"Single tones never jump" holds.**
  - It also holds in ten further runs at 7, 11.5, 18.9, 25 and 29.3 rad/s, s = 0.08 and 0.12. The largest distance is ≤ 0.42 there, although the first coefficient reaches 2.3 (mine).
  - At 10 rad/s and s = 0.04 the edge is 0.70–0.82 under four seeds, against a drift of +0.47 … +0.68 (mine; raw).

**What does not hold (majors first)**

1. **C7's "the least window edge decides the jump" is not a test (D3bw-1).**
   - It compares a second-order proxy at period 300 with the least of eight windows taken from periods 0–185.
   - Using each run's own coefficient it scores 11 of 14, not 13 of 14.
   - With 30 windows instead of eight, the verdict at m = 4, s = 0.08 flips.
   - Out of sample it scores 8 of 12, against 6 of 12 for a constant guess.
   - Done time-resolved and like for like, the same windows do much better on my runs.
2. **Under noise the "edge" is not a threshold (D3bw-2).** Ejection is not monotone in the coefficient: there are ejecting bands with gaps between them. This, not a threshold, explains C6's draw 3.
3. **The (k_e/v₁)^{1/2} law fails outside its fitted window and leaves out the orbit's speed (D3bw-3).**
4. **Minors:**
   - E·s² is clean for 13.1 rad/s only (D3bw-4);
   - the anatomy's quasi-static error and its "downward leg" are not borne out (D3bw-5);
   - s*'s framing, and its draw-specific part carries no information (D3bw-6);
   - 27 behaviour-changing mutants survive (D3bw-7);
   - md residue (D3bw-8); sources and novelty (D3bw-9); cosmetics (D3bw-10).

## What I checked and how

- **Reproduction:**
  - `repro.py` reruns the cells with the job functions and compares them with the raw field by field.
  - `num_check.py` looks up every D60 number in the md.
  - By hand: the algebra of (i) and (ii), and the sign argument in (iii).
- **C7 like for like** (`drift_level.py`, `c7_analysis.py`, `windows*.py`):
  - C's loop under +n and −n for C7's band cases and for draws 1–3, giving the period-mean coefficients;
  - draw 7's windows tiled every 10 periods over periods 0–290, at k_e = 20 (s = 0.04, 0.06) and at m = 4 (0.06, 0.08);
  - draws 1–3 with the draft's eight windows, at the same configurations and levels;
  - draws 1–3 at m = 4 tiled;
  - draw 7's windows around its leaves.
- **Thresholds** (`scans*.py`): `edge_run` on fine grids of the raise Δ, for:
  - the noise-free loop, 21.3 rad/s and 13.1 rad/s;
  - the band's windows, and draw 3 near its leave;
  - k_e = 2.5 and 5;
  - the k_e = 20, s = 0.08 windows at Δ = 0.
- **Single tones** (`tones.py`):
  - `c7_job` under −n at 10 and 13.1 rad/s;
  - 10 rad/s at s = 0.04 (four seeds) and at 0.06; 13.1 rad/s at 0.04; 21.3 rad/s at 0.025;
  - `c_job` at 7–29.3 rad/s.
- **Scaling** (`edge_ext.py`, `edge_fit.py`):
  - `edge_job` at m = 1.5 and 24; at k_e = 2.5 and 80; under WFR at k_e = 20 and 40; kernel 3 at three more settings; and with CT.OM = 0.1 and 0.4;
  - the power law refitted with standard errors.
- **Anatomy** (`anatomy*.py`): dense traces of c(t) and ẽ(t), 4000 points a period, on both sides of the edge, with the raise applied at x_m = 0 and at x_m = −1.
- **s\*, C4 and C6** (`star_analysis.py`, `c4_selection.py`, `c6_check.py`):
  - s* scaled by λ;
  - the draws permuted within each configuration;
  - rise/s²⟨w₂⟩ in the periods before the leave;
  - raised against un-raised coefficients.
- **Mutants** (`mut/driver.py`):
  - 50 single-point mutants of the new code, plus 16 replicas of round 69's survivors;
  - each applied in a private tree, run against the targeted tests, then `--from-raw` and `cmp`;
  - targeted survivors that are not equivalent were rerun against all 47 tests. The rest call functions that only their targeted test reaches (checked with grep), so the targeted run already covers them.

## Verdicts per part

| claim (D60.md line) | verdict | evidence |
|---|---|---|
| (i) sensitivity, the slope cancelling, \|S\|² − 1, \|S(√K)\|, Bode (18–25, 88–92) | verified | algebra; relative degree 2, integrator on the axis, stable loop; quadrature 3·10⁻¹² (md) |
| (i) bound on m; A2's factors (23–25) | verified | 1.495/1.248/1.370 → 1.225/1.095/1.152 (md) |
| (ii) ⟨z⟩_j = Σ_v(Pv)a_j + Re(e^{iΩjT}H) and its bound (34–40) | verified, Re missing | D3bw-10 |
| (ii) averaged loop = constant part's response at O(s²); B5 (43–49, 166–171) | verified | 0.9979–1.0054 at s = 0.00125 (md) |
| (iii) noise-free edge 0.0731, 0.0744 (51–53, 203–206) | verified | bit-identical; monotone on a fine grid (mine) |
| (iii) repulsion, "ẽ of the sign of −w_i" (53–54, 105) | at the onset only | D3bw-5 |
| (iii) carried out "on its downward leg" (54–57, 206–208) | only at hi, from the draft's start | D3bw-5 |
| (iii) edge ∝ (k_e/v₁)^{1/2} "across m, k_e and the metric" (57–60, 134–135) | interpolation only | D3bw-3 |
| (iii) C4: the rise is the even part at the level, 1.02–1.59 × s²⟨w₂⟩ (61–64, 200–203) | verified | twin bit-identical; the ratio is steady before the leave (mine) |
| (iii) s*: 33 of 52 (64–67, 210–212) | numbers verified, overstated | D3bw-6 |
| (iii) single tones raise the edge by E·s² (67–72, 216–218) | for 13.1 rad/s only | D3bw-4 |
| (iii) single tones never jump (72, 186–187) | verified | ten more runs (mine) |
| (iii) the least window edge decides the jump: 3 of 3, 10 of 11 (72–78, 218–222) | not established | D3bw-1, D3bw-2 |
| (iii) C6: "the drift sets when … the fluctuation when" (78–80) | contradicted by draw 3 | D3bw-2 |
| Status (108–119) | mostly accurate | omits the non-monotone edge and the law's range |
| Classical / New (121–138) | partly | D3bw-9 |

## Status of the round-69 issues

| id | status | what remains |
|---|---|---|
| D3bp-1 | fixed | C4 and W2 are like for like and reproduce; the excess is a higher order, not selection (mine) |
| D3bp-2 | partly | The noise-free edge, its anatomy at hi and its trends are in and right. The noise criterion (C7), the scaling law and E·s² are not established (D3bw-1 to D3bw-4) |
| D3bp-3 | fixed | only a missing Re (D3bw-10) |
| D3bp-4 | fixed | — |
| D3bp-5 | fixed | — |
| D3bp-6 | fixed | — |
| D3bp-7 | fixed | — |
| D3bp-8 | fixed | — |
| D3bp-9 | fixed | — |
| D3bp-10 | partly | 14 of 16 replicas of round 69's survivors are now killed. Two survive: C3c (`c3_job` stores w as C) and AR1 (`a2_rows`' K = 0 reference without f̂′). The new code adds 27 survivors (D3bw-7) |
| D3bp-11 | fixed | the bursting analogy is overstated (D3bw-9) |
| D3bp-12 | fixed | — |

## Issues

| id | sev. | issue | fix |
|---|---|---|---|
| D3bw-1 | major | **C7's band verdict — "the least window edge decides the jump", 3 of 3 and 10 of 11 (D60:72–78, 218–222; `c7_lines` 1557–1581) — is not a test.**<br>• **Not like for like.** The "drift" is the second-order proxy s²⟨w₂⟩ at period 300, although C4 shows that the even part at the level exceeds it by 2–59% (md). On the same C runs, each run's own period-mean coefficient at period 300 is 1.25× the proxy at m = 4, s = 0.06; 1.38× at 0.08; 1.41× at k_e = 20, s = 0.04; and 0.58× under WFR at 0.24 (mine).<br>• **With the run's own coefficient the score is 11 of 14.** Three runs exceed their least window and never jump: m = 4 at 0.06 (0.180 against 0.179), m = 4 at 0.08 (0.354 against 0.293), and k_e = 20 at 0.04 (0.187 against 0.168). The WFR "exception" disappears (0.487 against 0.655) (mine).<br>• **Wrong time.** The eight windows cover periods 0–185, 27% of the horizon. Draw 7's jumps at period 147 (s = 0.03) and at period 41 (k_e = 20, s = 0.08) fall in no window (md; raw).<br>• **The statistic depends on the window count.** With 30 windows the least edge falls to 0.237 at m = 4, s = 0.08, below even the proxy 0.256, so the verdict flips. At k_e = 20 it falls to 0.158 (s = 0.04) and 0.309 (s = 0.06) (mine).<br>• **Out of sample** (draws 1–3; m = 4 at s = 0.06 and 0.08; k_e = 20 at 0.04 and 0.06; the draft's rule) it scores 8 of 12 and predicts 3 of the 6 jumps. The noise-free edge, or either constant guess, scores 6 of 12 (mine).<br>• **At k_e = 20, s = 0.08 no drift is needed.** In the window from period 50 the kernel moves only 0.33 away at Δ = 0, a leave, not a jump. The window from period 35 ejects it from the exact fit, 1.11 away within 2 periods (mine), so the draw's jump at period 41 is pure noise.<br>• **Done time-resolved and like for like, the windows predict.** I compared each run's own coefficient at a window's start with that window's edge (mine):<br>&nbsp;&nbsp;– it flags no window in the 4 tiled runs that did not jump;<br>&nbsp;&nbsp;– it flags the window holding the leave in 5 of 5 jumps: windows 123–128 for draw 7's leave at 133, and 137–143 for its leave at 147 (near-leave windows only); for draws 3, 2 and 1 at m = 4, windows 150, 250 and 260 against leaves at 157, 256 and 263;<br>&nbsp;&nbsp;– it raises one false alarm (draw 1 at period 220). | • Replace the verdict by the time-resolved comparison with each run's own coefficient (or the ±n half-sum at the level), using windows that tile the horizon and an edge defined by a jump.<br>• Test it on draws it was not built on, against a null rule.<br>• Until then, drop "the edge decides the jump" and "its least window deciding the jump" from (iii), the Status and New. |
| D3bw-2 | major | **Under noise the edge is not a threshold (D60:61, 67–80, 136–137). `bracket_edge` assumes it is monotone.**<br>• **Fine scans of the raise find ejecting bands with gaps (mine):**<br>&nbsp;&nbsp;– 21.3 rad/s, s = 0.02: Δ = 0.20 ejects, 0.21–0.26 do not, from 0.27 on all do. The md's edge, 0.198, is the isolated band, and its own recheck prints "0.25: no";<br>&nbsp;&nbsp;– draw 7, window from period 125 at 0.02: 0.0675 ejects, 0.07–0.0725 do not;<br>&nbsp;&nbsp;– draw 3 at 0.02, windows from periods 147 and 150: 0.085–0.115 eject, 0.125–0.15 do not, 0.16 and above do;<br>&nbsp;&nbsp;– without noise the edge is monotone.<br>• **These bands, not "the drift sets when … the fluctuation when" (D60:80), explain C6's draw 3** (raw):<br>&nbsp;&nbsp;– at period 157 the un-raised coefficient, +0.112, lies in the ejecting band, and that run leaves there;<br>&nbsp;&nbsp;– the raised runs' coefficients, +0.132 … +0.153, lie in the gap, and those runs stay until 166–192;<br>&nbsp;&nbsp;– yet the raised coefficients were higher at every earlier period, by at least 0.019.<br>&nbsp;&nbsp;A time-resolved rule keyed to the lowest band would also misjudge these runs.<br>• **Leave versus jump.** The edge is a leave threshold (0.3 away), while the verdict counts jumps (0.6 away). At k_e = 20, s = 0.06, draw 7 leaves at period 41 (0.314) and is scored "did not jump" (raw). | • Report each window's ejecting set in Δ from fine scans, not a bisected edge.<br>• Define the edge by a jump.<br>• State that under noise ejection is not monotone in the coefficient.<br>• Rewrite C6's conclusion. |
| D3bw-3 | major | **The scaling ≈ (k_e/v₁)^{1/2} "across m, k_e and the metric" (D60:57–60, 134–135; `edge_lines` 1469–1486) is an interpolation inside its fitted window.**<br>• **The fit** (the raw's 15 rows; mine): α = 0.537 ± 0.025 and β = 0.555 ± 0.028. The residuals are curved: +0.14 and +0.11 at the corners, −0.02 … −0.05 in the centre.<br>• **Against the law 0.0809·(k_e/v₁)^{1/2} outside the window (mine):**<br>&nbsp;&nbsp;– k_e = 2.5 at m = 12: 0.0594, the same as at k_e = 5, i.e. 1.53× the law;<br>&nbsp;&nbsp;– WFR at k_e = 20: 0.581, 1.61× the law;<br>&nbsp;&nbsp;– WFR at k_e = 40: no ejection up to Δ = 1, where the law gives 0.51;<br>&nbsp;&nbsp;– m = 24: 0.90×; k_e = 80: 1.13×; m = 1.5: 1.04×.<br>• **Kernel 3:** Δ_e(v₃/k_e)^{1/2} is 0.0706 (raw) and 0.0715–0.0777 (mine), below the quoted range 0.0756–0.1017, which is built from kernel 1 only.<br>• **The orbit's speed:** at m = 12, k_e = 10, Δ_e = 0.0431, 0.0731 and 0.1253 at ω_ref = 0.1, 0.2 and 0.4, i.e. ∝ ω_ref^0.76–0.78; also 0.68–0.74 at k_e = 20 and 0.72–0.73 at m = 3 (mine). The law has no ω_ref, σ or w in it. | • State the exponents with their errors and the window they come from: m = 3–12, k_e = 5–40, ω_ref = 0.2, one target.<br>• Print the deviations outside that window.<br>• Drop "across the metric" and the law from New, or derive a scaling that includes ω_ref and test it. |
| D3bw-4 | minor | **Single tones: "raises it by E·s²" (D60:67–72, 136, 216–218).**<br>• **Only 13.1 rad/s is clean.** The +n and −n edges agree, 0.0670/0.0664 at s = 0.0025 up to 0.4070/0.4070 at 0.02, and E = 797 and 834 (mine).<br>• **10 rad/s (50ω_ref, commensurate):** the edges under +n and −n differ (md; mine). So the shift has an odd part of about ±3s, and the even part's E is 692, 734 and 500 at s = 0.005, 0.01 and 0.02, not 100–419.<br><br>&nbsp;&nbsp;\| s \| 0.0025 \| 0.005 \| 0.01 \| 0.02 \|<br>&nbsp;&nbsp;\|---\|---\|---\|---\|---\|<br>&nbsp;&nbsp;\| +n (md) \| 0.0795 \| 0.0758 \| 0.1087 \| 0.2410 \|<br>&nbsp;&nbsp;\| −n (mine) \| 0.0648 \| 0.1054 \| 0.1847 \| 0.3055 \|<br><br>• **15 rad/s:** E falls from 1537 to 1149 (md).<br>• **21.3 rad/s:** E = 209 and 311 at s = 0.01 and 0.02 (md), then 0.773 at 0.025, i.e. E = 1119 (mine), and > 1 at 0.03 (md).<br>• **"Never jump" holds, but by a modest margin.** At 10 rad/s, s = 0.04, the edge exceeds the drift by 18–32% (0.702–0.824 against +0.474 … +0.678), not "beyond its drift's reach". At 13.1 rad/s, s = 0.02, the edge equals the second-order drift (0.407 against 0.404, md). | • Run ±n and report the even part.<br>• Give E per tone and level rather than as a law.<br>• Measure the edge at C2's levels, 0.04–0.08.<br>• Soften New. |
| D3bw-5 | minor | **The anatomy (D60:53–57, 105–106, 206–208).**<br>• **"ẽ ≈ −Δφ_i/k_e near the feature."** The measured ratio ẽ/(−Δ·sgn(w)·φ/k_e) has median 0.27–0.30 and takes both signs during the passage (mine). The centre's fast response cancels most of the excess.<br>• **"On its downward leg."** Below the edge the excursions run ahead of the orbit point in both directions: kernel 1 moves −0.068 on the way down and +0.041 on the way up; kernel 3 moves −0.087 and +0.049 (mine). Started at x_m = −1, kernel 1 raised by 0.1 is carried upward on its first passage, to +1.36 (mine). Under noise, kernel 1 ends at +1.39 … +1.41 in 3 of 11 C3 draws at s = 0.02 and 8 of 11 at 0.03 (raw). Why the downward excursions are larger is not explained.<br>• **"0.05 off" as the moment of detachment.** Kernel 3 at Δ = 0.075 reaches +0.051 at t/T = 0.113 and comes back (mine). The md's 0.406 for kernel 3 relies on 64 samples a period missing that crossing. | • Present the sign argument as describing the onset only.<br>• Say that the ejection's direction depends on the passage, and that under noise it goes both ways.<br>• Mark detachment by a non-returning criterion on a fine time grid. |
| D3bw-6 | minor | **s\* (D60:64–67, 210–212).**<br>• "33 of 52" leaves the 17 runs that never jumped out of the denominator; it is 33 of 69 testable cases.<br>• The draw-specific part carries no information: with the draws' ⟨w₂⟩ permuted within each configuration, the consistent count averages 32.95 over 200 permutations (mine).<br>• The "consistent" windows are 0.20–0.50 of s* wide.<br>• "Early … where the edge is small" includes k_e = 40 (0.83×).<br>• The list of the 17 omits the 25–75 band at m = 12 (md). | • Report 33 of 69, the permutation null and the window widths.<br>• Say that the prediction is per configuration, not per draw. |
| D3bw-7 | minor | **27 behaviour-changing mutants of the new code survive the tests; 15 are killed and 8 are equivalent.** Mutants of the following survive:<br>• `c7_lines`, all 8: the drift taken at period 1; jumped = left; any seed; every configuration's base set to C7_CFG's edge; the second-least window; the median verdict replaced by the least; the WFR W2 lookup broken; the rise counted only from s > 0.0025;<br>• `edge_lines`, 5: kernel 3's s* computed from kernel 1's ⟨w₂⟩, or from its jumps; the midpoint replaced by hi; WFR's v₁; kernel 3 added to the fit;<br>• `edge_run`, 3: the wrong kernel raised; the 0.05 marker; the centre "after";<br>• `edge_job`: the anatomy taken at lo;<br>• `c7_job`: no rechecks; every window at period 0;<br>• `bracket_edge`: an extra recheck;<br>• `b5_lines`: the summary over the last mark only;<br>• `c3_lines`, 2; `c4_lines`, 2 (a twin that leaves at p; r² for r);<br>• `w2_cases` with 10 draws;<br>• `bode_integral` with log₁₀ (trivial).<br>Nineteen of them change the md on the raw. The C7 test's `': yes' in md` (tests:828) passes for any md. | Add tests:<br>• `c7_lines` on synthetic windows with a time-varying W2, two configurations with different noise-free edges, a run that leaves without jumping, and other seeds;<br>• `edge_lines` with kernel-3 and WFR rows and wide brackets;<br>• `edge_run` and `edge_job` on a fake loop (which kernel is raised; the anatomy at hi);<br>• the noise and the rechecks `c7_job` passes on;<br>• edge cases of `c3_lines` and `c4_lines`; the length of `w2_cases`; `c3_job`'s C; `a2_rows`' jit0. |
| D3bw-8 | minor | **md residue and presentation.**<br>• `c3_lines` (1148–1161) still prints "The first kernel's drift against (ii)'s even part at its saturation" (md:352). That is the comparison D3bp-1 showed is not like for like, and it has no caveat.<br>• `c7_lines` drops the '>1' windows from the median: at k_e = 20, s = 0.06 the median is 0.614, and 0.658 with those windows counted as +∞.<br>• The single-tone "rise per s²" list drops the '>1' rows (15 and 21.3 rad/s at s = 0.03), which contradict E·s².<br>• C7 prints brackets, while D60 quotes their midpoints. | • Drop the C3 paragraph or caveat it.<br>• Count '>1' windows as above every finite edge.<br>• Keep the '>1' rows in the list.<br>• Print the midpoints. |
| D3bw-9 | minor | **Sources and novelty (D60:121–138).**<br>• New lists the criterion and the law of D3bw-1 to D3bw-3, neither of which is established.<br>• "The phenomenology … is bursting's": Anderson's bursting is intermittent (a drift, a burst, then a return; checked). Here the kernel leaves for good and ends beyond ±1.<br>• Accurate [memory]: Bode with an integrator on the axis; KRP 1985 (checked); Egardt; Ioannou–Kokotović. Hoppensteadt and Khalil §10.1 are rightly confined to the formed block, because 36 multipliers are neutral. | • Restrict New to what is derived or measured.<br>• Say "a drift, then a jump, as in bursting but without return". |
| D3bw-10 | cosmetic | • The mode sum in ⟨z⟩_j (D60:38) is complex and needs Re.<br>• "One jump at 0.32" under WFR at k_e = 10 is kernel 3's jump (D60:193).<br>• `raised`' test compares magnitudes only, so a flip of kernel 2's sign passes (tests:601). | wording; a sign check in the test |

## Sources

- H. W. Bode, *Network Analysis and Feedback Amplifier Design*, 1945; J. S. Freudenberg, D. P. Looze, IEEE TAC 30(6):555–565, 1985 — [memory]. Integrators in Bode's integral: Seron, Braslavsky, Goodwin, *Fundamental Limitations in Filtering and Control*, 1997 — [memory]. The search gave background only: [arXiv 1902.11302](https://ar5iv.labs.arxiv.org/html/1902.11302).
- B. D. O. Anderson, Automatica 21(3):247–258, 1985 — checked against the [ANU repository](https://openresearch-repository.anu.edu.au/entities/publication/b9a832b6-dae1-41f8-bf77-66c2fe5de1df). Bursting as intermittent: [arXiv 2108.11336](https://arxiv.org/pdf/2108.11336), [Espana–Praly 1993](https://cas.minesparis.psl.eu/~praly/Telechargement/Journaux/1993-SIAM-Espana-Praly.pdf).
- P. V. Kokotović, B. D. Riedle, L. Praly, Systems & Control Letters 6(1):7–14, 1985 — checked: [PDF](https://cas.minesparis.psl.eu/~praly/Telechargement/Journaux/1985-SCL-Kokotovic-Riedle-Praly.pdf).
- B. Egardt, *Stability of Adaptive Controllers*, LNCIS 20, 1979; P. A. Ioannou, P. V. Kokotović, Automatica 20(5):583–594, 1984; F. C. Hoppensteadt, Trans. AMS 123:521–535, 1966; H. K. Khalil, *Nonlinear Systems*, 3rd ed., §10.1 — [memory].

## Lesson

A threshold measured under noise has to be shown to be a threshold before it is used as one. Scanning the raise finely shows ejecting bands with gaps between them, so a bisected "edge" is just the lower boundary of whichever band the grid happens to hit. A criterion that is meant to decide events also has to be tested at the events' time and level, and on draws it was not built on. The least of eight windows, set against a second-order proxy at period 300, agrees 13 of 14 times on its own runs. With each run's own coefficient it agrees 11 of 14 times, and on new draws 8 of 12. Yet the same windows, compared with each run's own coefficient at each window's start, flag the windows that contain the leaves.
