# Referee report D3bj, round 63: Proposition D‴.60 as revised after round 56

## Scope

- **Commit reviewed:** ada0e1c ("D‴.60 round 56: ejection as a first-passage event, the period means as Floquet sums, the revised proposition").
  - Other sessions moved HEAD (now f3061bd) and left uncommitted changes to unrelated files.
  - `git diff ada0e1c HEAD` is empty over every file below, and the md and raw md5s are unchanged.
- **Files read:**
  - review/drafts/D60.md (110 lines);
  - experiments/exp_noise_resonance.py and its imports: exp_closed_transit.py as CT, exp_dead_zone.py as DZ, exp_cohesion_leak.py (load/save), src/core/closed_loop.py, closed_loop.c, src/core/mirror_free.py;
  - results/chaos/noise_resonance.md and noise_resonance_raw.npz;
  - tests/test_noise_resonance.py;
  - review/REVIEW_round56_D3bc.md and the "Round 56" entry of review/REVIEW_LOG.md (lines 2286–2327);
  - PARADIGM.md D‴.31 (line 3808) and D‴.33 (line 4164). D‴.49 is not in PARADIGM.md at this commit, so I read review/drafts/D49.md.
- closed_loop.c is unchanged since round 56 (17d6217), so round 56's probe scripts still apply. I re-ran the one I rely on.
- **I modified nothing in the repository.** All work was in a `git archive ada0e1c` copy under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref63/D3bj_work:
  - W/repo: the copy; W/out: the raw copy and outputs; W/scripts: my scripts; W/mut: the mutants.
- **Run settings:** every run used `nice -n 10`, OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1, LCK_BUILD=W/build and pytest `-p no:cacheprovider`.
- **Load:** at most two heavy processes of mine at a time (one Pool(2), or two mutant drivers). Waits ran only in background shells, and no process of mine is left.

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (md5 94b5f10966f067843a0aee13ab7d7e1b).
  - The 24 tests pass in 28.0 s.
  - Four cells rerun with the experiment's own job functions are bit-identical to raw:
    - B4 at 21.4 and with the 16 tones, s = 0.00125, every field;
    - C at m = 12, k_e = 10, s = 0.015, seed 2 (left 256, −1, −1);
    - C2 #45 (pair 18.9/22.4, seed 7).
- **(i) re-derives.** The oscillator law, K(t), √K_i = (Gm)^{1/2}|w_i|/(σe^{1/2}) = 15.758, 10.505, 13.132 (md 15.76, 10.51, 13.13; at m = 3: 7.88, 5.25, 6.57) and the bound 100/20.693 = 4.833 all check.
- **(ii)'s expansion and parity are exact.** For single tones, the B3 raw shows a cubic odd remainder (×7.5–8.2 per doubling of s) and an even-part error ∝ s².
- **The first-order "rotating increments" are real, not an integrator artefact.**
  - scipy's DOP853 (`second_order`'s first order) reproduces the C loop's B4 at 21.4: ⟨w₁⟩ = −0.00112, −0.00896, −0.0355 at periods 1, 3, 10, against −0.00112, −0.0089, −0.0353. It also matches at 21.3.
  - The pattern holds at six more tones I ran:
    - a constant per-period centre mean at 10.0, 15.0, 25.0 and 29.4 (50, 75, 125, 147·ω_ref);
    - a period-4 cycle at 15.05 (+0.248, +0.431, −0.248, −0.430);
    - alternation at 15.1 (±0.286).
- **Single-tone ⟨θ₂⟩ relaxes with the loop's slow Floquet modes.**
  - The monodromy multipliers are 0.99400, 0.99245, 0.99237 (τ = 166, 132, 131 periods).
  - B4's even-part increments shrink ×0.685 (kernels 1, 3) and ×0.740 (kernel 2) per 50 periods, which is exactly ρ⁵⁰.
- **C's cells are robust to seeds.** Seven more seeds (4–6, 8–11) reproduce them:
  - m = 12, k_e = 10: kernel 1 leaves at s = 0.015 in 1 of 7 seeds, and at 0.02 and 0.03 in 7 of 7; kernel 2 never leaves up to 0.03.
  - m = 3, k_e = 5: kernel 1 leaves at 0.06 in 7 of 7, and nothing leaves at 0.04.
- **WFR at k_e = 10 never ejects up to s = 0.32** (largest distance 0.19).
- **Every number in the draft is printed by the md.** A few ranges are pooled by hand from printed values.

**What is wrong (majors first)**

1. **Tone pairs eject (D3bj-1).** "Single tones and tone pairs at the same levels do not eject" is false for pairs.
   - C2 ran pairs only at s = 0.02, and only with the two difference classes round 56 had found inert.
   - Round 56's ejecting pair still ejects at s = 0.02 (period 218, distance 0.90).
   - C2's own pairs eject at s = 0.08 in 2 of 8 runs.
2. **"Where they add up secularly" is false (D3bj-2).** The linearised loop is exponentially stable, so commensurate period means saturate.
   - At 21.4 = 107ω_ref the first kernel's odd part per s is −0.467 at period 300 and −0.5215 at period 1200.
   - The statement's bound and the proof's bound differ by a factor 2.
3. **13.1's averaging failure is commensurability, not the resonance crossing (D3bj-3).**
   - At 13.05 and 13.15, where K₁(t) crosses ω² just as at 13.1, the averaged loop matches the exact second order to 0.986 and 1.004. At 13.1 the ratio is 0.666.
   - A third of 13.1's "resonant" drift is the 2ω part adding coherently, because 2·13.1 = 131ω_ref.
4. **The 16 tones' "saturation" is a 300-period snapshot (D3bj-4).** Over 1200 periods the first kernel's ⟨w₂⟩ moves from the quoted 225–234 to 254–285, and the centres' even part never settles.
5. **Resonance is not isolated, and the bound is not a threshold (D3bj-5).**
   - Moving the band to 25–75 is mostly roll-off: ×2.2 of the ×2.9.
   - Both the resonance factor and the ejection rate rise smoothly with m through 4.83.
   - C2's single-seed m = 4.8 row is not representative.
6. **Smaller issues:**
   - "New" still claims a classical condition, and one citation is imprecise (D3bj-6).
   - C2 is single-seed, though the Status says "four noise seeds per cell" (D3bj-7).
   - Ejection timing (D3bj-8) and the WFR margin at k_e = 5 (D3bj-9).
   - The 0.3 cut counts drifts as leaves (D3bj-10).
   - B3's odd-part sentence is misleading (D3bj-11).
   - 24 of 51 mutants survive, 21 of them non-equivalent, all in untested new code (D3bj-12).
   - Two cosmetic items (D3bj-13, -14).

## What I checked and how

All scripts are in W/scripts; logs and outputs are in W and W/out.

**Reproducibility**
- I copied the raw to W/out, ran `NR_OUT=W/out venv/bin/python experiments/exp_noise_resonance.py --from-raw`, and compared with `cmp`: identical.
- `pytest tests/test_noise_resonance.py`: 24 passed.
- `rerun_cells.py` ran `b4_job` (21.4 and the 16 tones), `c_job` (C, and C2_CASES[45]) and compared with raw: maximum relative difference 0.0.

**Derivations**
- I re-derived (i) and the proof of (ii) by hand.
- `floquet.py`: the monodromy of the 7-dimensional loop linearised at the fit (e, w₁₋₃, c₁₋₃; DOP853, rtol 10⁻¹¹), for the paradigm's metric at m = 12 and 3 and for WFR.
- The 18 unformed kernels are frozen: μ = 8.3·10⁻¹⁰ and v = 0 at α = 0.
- `first_order_check.py`: `second_order`'s first order against B4's odd part.
- `where_dc.py`: where in the period the first order's net change of w₁ accrues at 21.4. It is −0.0039 per period per unit s, about 60% on kernel 1's downward pass and 40% on its upward pass. This is the non-adiabatic response to the fast K(t) sweep.

**Long horizon and other tones** (the experiment's `sim_means`, ±s·n, s = 0.00125)
- `long_b4.py` and `long_analyse.py`: 1200 periods at 21.4 and with the 16 tones.
- `odd_tones.py`: 100 periods at 10.0, 15.0, 15.05, 15.1, 25.0, 29.4.
- `avg_commens.py`: `second_order` and `envelope_means` (s = 0.005) at 13.05 and 13.15 over 10 periods, against raw's 13.1.
- `even_err.py`: the loop's even part at 13.05 and 13.15.

**Ejection controls** (`controls.py`: `NR.simulate` exactly as in `c_job`, 300 periods, cut 0.3, m = 12, k_e = 10 unless stated; 149 runs)
- Round 56's pair (seed 7's 22.5019/18.9015 with the draw's phases) at s = 0.02, 0.04, 0.08.
- C2's two pairs under four seeds at 0.04 and 0.08.
- One tone at 13.1, 15.7, 21.3 at 0.06 and 0.08.
- Integer-multiple tones 10.0, 15.0, 16.0, and the draw's 16.0033 with its phase, at 0.02, 0.04, 0.08.
- WFR at k_e = 5 and 10, s = 0.12–0.32.
- Seeds 4–11 at (m = 12, k_e = 10, s = 0.0075–0.03) and at (m = 3, k_e = 5, s = 0.04–0.08).
- m = 3.5, 4.0, 4.5, 4.83, 5.5 × s = 0.06, 0.08 × seeds 7, 1, 2, 3.
- `inspect_C.py`: per-seed tables and leave-time statistics from raw.

**Quasi-static analyses** (A's `transfer`)
- `resonance_vs_m.py`, `band_split.py`, `peaks_qs.py`.

**Mutants**
- `mutants.py` defines 51 single-point mutants:
  - `b4_job` 8; `c_noise`/`c_job` 8; `part_C`/`part_C2`/`part_B4`/C2_CASES 4; `c_cells` 5;
  - report code (B4, C, C2, B3 centres, B column, A per metric) 16;
  - round 56's targets 10 (`envelope_rhs` ē terms V2–V6, `b3_job` T1–T4, `b2_jobs` T5).
- `mutdriver.py` applies each mutant in a private copy, runs the test file (`-x`) and `--from-raw` (byte comparison). `mut_summary.py` tallies the results.

**Literature**
- Web search only, so sources are marked [snippet] or [memory]. www2.eecs.berkeley.edu could not be fetched from this session.

## Verdicts per part

| claim (D60.md line) | verdict | evidence |
|---|---|---|
| (i) first order = the linearisation along the orbit (8) | verified | `variational` is the exact linearisation (unformed kernels frozen); A: loop over linearisation 0.9899–1.0072 (md:30) |
| (i) oscillator law, δw_i = (μ_iφ_i/K)y, δc_i = (v_i∂_cφ_i/K)y (12) | verified (quasi-static derivation) | ẏ = Kd₁ + orbit terms, ḋ₁ = ṅ − f′n − y − k_ed₁; A: 0.999–1.082 at ω ≥ 10 |
| (i) √K_i = (Gm)^{1/2}\|w_i\|/(σe^{1/2}) (15–16) | verified | K_i max 248.35, 110.39, 172.47 → 15.76, 10.51, 13.13 |
| (i) bound m ≲ ω_lo²/(G max_tΣw_i²(∂_cφ_i)²) (17) | algebra verified (4.833) | its meaning: D3bj-5 |
| (i) "a tone near the top of the sweep is amplified" (14); A "peaks near each kernel's √K_i" (79) | roughly | quasi-static peaks 14.16, 9.16, 11.64 = 0.90, 0.87, 0.89·√K_i (D3bj-13) |
| (ii) even part s²⟨θ₂⟩ + O(s⁴), odd part s⟨θ₁⟩ + O(s³) on finite horizons (19–22) | proved | analyticity; the start does not depend on n; B3 raw scaling as above |
| (ii) rotating increments, alternation at half-integer tones (22–25) | verified at 9 tones | scipy = C loop; 15.05 period-4 cycle, 15.1 alternation |
| (ii) "bounded, by about 1/\|sin(ω_kT/2)\|" (24); proof "1/(2\|sin(ω_kT/2)\|)" (50) | inconsistent; both ignore the damping | D3bj-2 |
| (ii) "unless ω_k is a multiple of ω_ref, where they add up secularly" (24–25) | refuted | multipliers ≤ 0.99400; 21.4 saturates at −0.5215/s (D3bj-2) |
| (ii) the averaged loop = the constant part's response (25–26) | verified (structure) | `envelope_rhs`; round 56's split |
| (ii) averaging domain (26–27); proof (53–54); Status "domain derived" (60) | partly refuted, not derived | 13.05/13.15 against 13.1; 29.3 (D3bj-3) |
| (ii) ⟨θ₂⟩ a bounded offset reached with the slow mode's τ (28–29); "saturation numerical (300 periods)" (60) | verified for single tones; for the 16 tones bounded, but the quoted values are wrong | Floquet τ = B4's ratios; 1200-period runs (D3bj-4) |
| (iii) ejection is a jump (31) | verified | C leaves: largest distance median 0.93, 90% above 0.76 |
| (iii) "first-passage …, after hundreds of periods" (31) | partly | leaves from period 7; near threshold none before period 114 (D3bj-8) |
| (iii) "single tones and tone pairs at the same levels do not eject" (32) | refuted for pairs; marginal for single tones | D3bj-1 |
| (iii) level falls with Gmw (32) | verified | C (m = 3 against 12), m-scan, extra seeds |
| (iii) band covering √K_max lowers the level, band above raises it (32–33) | not isolated from roll-off | D3bj-5 |
| (iii) level rises with k_e (33–34) | verified for the levels | C, m = 12: 0.0075, 0.015, 0.06, 0.08 (1/4) at k_e = 5, 10, 20, 40; the rectification attribution has no D‴.60 data (D3bj-14) |
| (iii) WFR: no ejection at the levels tried (34–35) | verified, small margin at k_e = 5 | k_e = 5 ejects from s = 0.12 (D3bj-9) |
| (iii) the bound or WFR "removes the resonance", not the ejection (35–36) | second half verified; first half only nominal | D3bj-5 |
| Numerics (73–109) | all numbers match the md | ranges at lines 85, 88, 96–98 pooled from printed values |
| Status (56–62) | mostly accurate | "domain derived" (D3bj-3); "four noise seeds per cell" (D3bj-7) |
| New (68–71) | partly overstated | D3bj-6 |

## Status of the round-56 issues

| id | status | what remains |
|---|---|---|
| D3bc-1 (major) | partly fixed | (iii) is now numerical and C2 adds controls. The pair control ran only the inert difference classes, at one level, and pairs do eject (D3bj-1). The band control is mostly roll-off and the bound is not a threshold (D3bj-5) |
| D3bc-2 (major) | fixed for C | 4 seeds, fractions, '*', medians and the 300-period horizon are stated; 7 more seeds confirm the cells. C2's 41 non-pair runs are single-seed (D3bj-7) |
| D3bc-3 (major) | partly fixed | Krylov–Bogoliubov and "slow beats" are gone, and the constant-part / rotating-increment split is right. The K-crossing condition is refuted, 13.1's failure is commensurability, and "derived" is unsupported (D3bj-3) |
| D3bc-4 (major) | partly fixed | Rotating increments, alternation and rotation are shown and verified. "Secular" is false, and statement and proof differ by a factor 2 (D3bj-2) |
| D3bc-5 (major) | partly fixed | The offset-with-τ picture holds for single tones. τ is not stated, and the 16 tones' figures are a 300-period snapshot (D3bj-4) |
| D3bc-6 (major) | partly fixed | The paper's own and the classical sources are now cited. "The commensurability condition" is still claimed new, and "Riedle–Kokotović" is imprecise (D3bj-6) |
| D3bc-7 (minor) | fixed in wording | No D‴.60 data on the sign reversal (D3bj-14) |
| D3bc-8 (minor) | fixed | C2's WFR rows; the margin is 1.5× at k_e = 5 (D3bj-9) |
| D3bc-9 (minor) | fixed | Status, D60.md:58–59 |
| D3bc-10 (minor) | fixed | The peaks sit about 10% below √K_i (D3bj-13) |
| D3bc-11 (minor) | partly fixed | 9 of my 10 replicas of round 56's targets are killed (V2–V6, T1, T2, T4, T5). T3, an offset in `b3_job`'s centres' even part, survives. 21 non-equivalent mutants of the new code survive (D3bj-12) |
| D3bc-12 (minor) | fixed | Every draft number is printed by the md |
| D3bc-13 (minor) | fixed | The md:36 column is against the even part |
| D3bc-14 (cosmetic) | fixed, one residue | "Three sources" and the bound 4.83 are printed. The odd-part sentence's contrast is misleading (D3bj-11) |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bj-1 | major | **(iii)'s "single tones and tone pairs at the same levels do not eject" and "multi-tone" (D60.md:31–32; C2, D60.md:105–106) are not supported.**<br>• C2 ran pairs only at s = 0.02, and only with the differences round 56 found inert (exactly 18ω = 3.6 and 17.5ω = 3.5) and fresh phases (exp_noise_resonance.py:62).<br>• At ada0e1c, round 56's ejecting pair still ejects kernel 1 at s = 0.02 (period 218, distance 0.90) and at 0.08 (period 273). That pair is seed 7's 22.5019/18.9015 with the draw's phases: Δω = 18ω_ref + 4·10⁻⁴, a 500-period beat against the orbit's 18th harmonic.<br>• C2's own pairs at s = 0.08: 18.9/22.5 with seed 7 ejects (period 124, distance 0.95); 18.9/22.4 with seed 3 crosses (period 69, distance 0.358); four more reach 0.22–0.297. At 0.04 none leaves (≤ 0.062).<br>• The 16 tones eject kernel 1 at 0.02 and 0.08 under 4 of 4 seeds.<br>• Single tones never jump. But at s = 0.08 three of seven cross the cut: 13.1 (period 150, distance 0.312), 15.0 (period 132, 0.356) and the draw's 16.0033 (period 57, 0.365). The others reach 0.23–0.27 (10.0, 15.7, 16.0, 21.3). | Rerun C2's pairs at 0.02–0.08 with near-commensurate differences and the draw's phases. Restate (iii): pairs whose difference is near a harmonic of ω_ref eject at the 16 tones' levels; drop "multi-tone". Report jumps apart from the cut (D3bj-10) |
| D3bj-2 | major | **(ii)'s "bounded, by about 1/\|sin(ω_kT/2)\|, unless ω_k is a multiple of ω_ref, where they add up secularly" (D60.md:23–25) is false.** The proof's partial-sum bound (D60.md:48–51) treats the propagated past as the identity.<br>• The linearised loop is exponentially stable: multipliers 0.99400, 0.99245, 0.99237 (τ = 166, 132, 131 periods); the rest ≤ 0.202.<br>• So commensurate increments saturate at \|δ\|/(1 − ρ). At 21.4 = 107ω_ref the first kernel's odd part per s is −0.467, −0.510, −0.519, −0.5215 at periods 300, 500, 700, 1200 (my run). Increments shrink ×0.47 per 100 periods, which is ρ¹⁰⁰.<br>• The md's own B4 shows the slowdown already: increments −0.078, −0.053, −0.036, −0.025 per 50 periods from period 100 to 300.<br>• The statement's 1/\|sin(ω_kT/2)\| and the proof's \|1 − e^{iω_kT}\|⁻¹ = 1/(2\|sin(ω_kT/2)\|) differ by a factor 2. The first is the maximum of the partial sums; the second is their eventual amplitude. | In (ii) and its proof, on each slow mode m_{j+1} = ρm_j + e^{iω_kjT}δ, so \|m_j\| ≤ \|δ\|(1 + ρ^j)/\|1 − ρe^{iω_kT}\|. At commensurability this is an offset ≈ \|δ\|τ/T, reached with τ, like ⟨θ₂⟩. Print the multipliers. Note that the increment δ falls steeply with ω: the per-period centre mean at commensurate tones is 1.31, 0.504, 0.033, 0.0082, 0.0009 per s at 10.0, 15.0, 21.4, 25.0, 29.4 |
| D3bj-3 | major | **The averaging domain (D60.md:26–27), the proof's sentence (D60.md:53–54), "domain derived" (D60.md:60) and "the largest at 13.1 …, near its resonance" (D60.md:84) misattribute the failure.**<br>• At 13.05 and 13.15, K₁(t) crosses ω² as at 13.1 (K₁ max 248 against ω² = 170–173), but 2ω = 130.5ω_ref and 131.5ω_ref. There the averaged loop over the exact second order is, for kernel 1, 0.955 / 0.986 at periods 3 / 10 (13.05) and 1.001 / 1.004 (13.15), against 0.663 / 0.666 at 13.1.<br>• The exact ⟨w₂⟩ at period 10 is 51.3, 76.0, 50.4 for 13.05, 13.1, 13.15. So 13.1's averaging failure, and a third of its "resonant" drift, is the 2ω part adding coherently because 2·13.1 = 131ω_ref. The K-crossing condition is not needed, and the sweep sentence's conclusion fails where its premise holds.<br>• Commensurability of 2ω_k is not sufficient either: 29.3 (2ω = 293ω_ref) averages to 1.0003–1.0019 (md:80), and so does WFR at 13.1, because the increments there are small.<br>• "After the orbit's modulation" is undefined.<br>• (The O(s⁴) part is resonance-related: the loop's even part over s²⟨w₂⟩ at s = 0.005 is 0.889 at 13.05 against 0.792 at 13.1.) | State the averaged loop's error as Σ over the dropped parts of \|δ_Ω\|/\|1 − ρe^{iΩT}\|, with the increments numerical. Drop the K-crossing condition. Mark the domain "numerical". Run every commensurate test tone with its non-commensurate neighbours, and disclose that 13.1, 15.7, 21.3 and 29.3 are half-integer multiples |
| D3bj-4 | major | **⟨θ₂⟩'s saturation (D60.md:28–29, 60, 90–93) holds for single tones, not as quoted for the 16 tones.**<br>• Single tones: verified; at 21.4 over 1200 periods ⟨w₂⟩ goes 82.3 → 93.2, −47.7 → −57.1, 60.0 → 66.7.<br>• 16 tones: "reaches +225 at period 200 and stays between +225 and +234" (D60.md:92–93) is a 300-period snapshot. The first kernel's ⟨w₂⟩ is 258, 257, 281, 259, 254, 264, 254, 285, 269 at periods 400–1200 (range 236–290 over periods 601–1200), and the second's goes −77 → −91.<br>• The centres' even part never settles: first kernel −171 to +12, third +2 to +132 over periods 601–1200. B3 quotes the period-30 snapshot −45, +9.9, +22 (D60.md:86).<br>• No time constant is stated. | Quote τ = 131–166 periods. Extend B4 to ≥ 5τ. For the 16 tones give bands over a long horizon, not period-300 values |
| D3bj-5 | major | **(iii)'s resonance attribution (D60.md:32–36) and C2's bound rows (D60.md:106–108) are not isolated.**<br>(a) The band control: moving the 16 tones from 10–30 to 25–75 cuts the quasi-static centre jitter per unit rms (m = 12) by 2.92, 2.61, 2.79. Of this, 2.24, 2.21, 2.22 is the K-independent roll-off and only 1.30, 1.18, 1.26 is resonance. So most of the ×3 ejection level (seed 7: 0.02 → 0.06) is roll-off.<br>(b) The bound is no threshold. The resonance factor (response over its K = 0 value; kernel 1; the 16 tones) is 1.09, 1.14, 1.19, 1.23, 1.28, 1.36, 1.50 at m = 2, 3, 4, 4.83, 6, 8, 12. With k_e = 10 the half-power width is k_e, and the band's lowest tone still gains ×1.42 at the bound's K peak.<br>(c) Ejection is smooth in m: at s = 0.06 kernel 1 leaves in 0, 0, 1, 2, 3, 4 of 4 seeds at m = 3, 3.5, 4.0, 4.5, 4.83, 5.5. C2's single-seed "m = 4.8 ejects only at s = 0.08" flips at m = 4.83 (seed 7 at 0.06, period 260). | Say the bound "reduces" the resonance, with the factor. Split the band effect into roll-off and resonance. Report the m-dependence over seeds |
| D3bj-6 | minor | **Novelty and sources (D60.md:63–71).**<br>• "With the commensurability condition of its period means" is classical:<br>  – a T-periodic linear system's forced response near a Floquet multiplier scales as 1/\|e^{iωT} − ρ\| [memory];<br>  – averaging with periodic signals [snippet];<br>  – the multi-frequency resonance relation k·ω = 0 [snippet].<br>• "Riedle–Kokotović" for drift under a disturbance (D60.md:65) points to the disturbance-free slow-adaptation papers. Disturbance-driven drift by averaging is Riedle–Cyr–Kokotović 1984 and Fu–Sastry 1987 [snippet].<br>• Related: LMS "non-Wiener" mean-weight shifts under narrowband interference [snippet]. | Limit "New" to the instance (Gmw as the stiffness, the exact second order with centres, the ejection data). Cite the sources |
| D3bj-7 | minor | **Status "four noise seeds per cell" (D60.md:61).** C2's 41 non-pair runs use seed 7 only (exp_noise_resonance.py:58–61), and one of them flips under a 0.6% change of m (D3bj-5c). | Run C2 over SEEDS_C, or mark the single-seed rows |
| D3bj-8 | minor | **"After hundreds of periods" and "first-passage" (D60.md:31).**<br>• C's leaves start at period 7, and 45 of 180 come before period 100; C2's m = 8, s = 0.08 leave is at period 41.<br>• Near threshold the times are not memoryless: at m = 12, k_e = 10, s = 0.02–0.03, all 22 kernel-1 ejections (11 seeds) fall in periods 114–274, after about one slow-mode τ. This suggests a build-up on the slow modes, which "first-passage" leaves unexplained. | Report the leave-time distribution per level against τ |
| D3bj-9 | minor | **WFR margin (D60.md:34–35, 108–109).** No ejection to 0.08 is confirmed, and at k_e = 10 none to 0.32. But at k_e = 5 WFR ejects from s = 0.12 (kernel 1, period 151), and all three kernels leave at 0.16 (periods 75–77). | State the margin |
| D3bj-10 | minor | **The 0.3 cut (D60.md:100) counts slow excursions as leaves.**<br>• 7 of C's 180 leaves have a largest distance of 0.31–0.49, e.g. m = 12, k_e = 20, s = 0.06, seed 7: 0.314 at period 41.<br>• Non-leaves reach 0.296.<br>• Single tones reach 0.23–0.37 at s = 0.08 (D3bj-1). | Report jumps (largest distance > 0.6) separately |
| D3bj-11 | minor | **B3's asymptotics (D60.md:86–87; md:80).** The sentence says the odd part follows s⟨w₁⟩ + O(s³) for single tones but departs for the 16 tones (2.78·10⁻⁴ against 4.80·10⁻⁴).<br>• Single tones depart more at s = 0.005: 21.37's first kernel shows 1.93·10⁻⁵ against 8.06·10⁻⁶ (×2.39), and 13.1 ×1.31 (md:71–72). Their remainder is cubic.<br>• The 16 tones differ by a non-cubic remainder (×14.5, ×20.6 per doubling), which the md does not print. Their first kernel's even error is non-monotone in s (0.30%, 1.19%, 0.19%), against the caption's "falls as s²". | Print the remainders' scaling and reword |
| D3bj-12 | minor | **Tests and mutants.** Of 51 mutants, the tests kill 16, `--from-raw` kills 11 more (all in report code), and 24 survive both.<br>• 3 survivors are (near-)equivalent: Cc5, R5, R10.<br>• The other 21 change results. They are:<br>  – all 8 in `b4_job`: the −n phase moved by π/2, /s instead of /s², no ½, marks off by one, the wrong kernel column, an offset centre, another seed's phase, the WFR metric;<br>  – both in `c_noise`;<br>  – all 6 in `c_job`, including the metric ignored, so C2's WFR rows would run the paradigm's metric, and k_e ignored;<br>  – all 4 drivers (`part_C` with one seed, `part_C2` without options, `part_B4` with one s, C2's WFR case under the paradigm's metric);<br>  – T3, an offset in `b3_job`'s centres' even part.<br>• Nothing tests `b4_job`, `c_noise` or `c_job`'s options, and `--from-raw` cannot see them. The short-horizon test does not check even_c (tests:394–408). | Add small tests: `b4_job` on 2 periods (as for `b3_job`); `c_noise`/`c_job` options reaching `simulate` (seed, band, freqs, metric, k_e); C2's labels; an even_c check |
| D3bj-13 | cosmetic | **"Peaks near each kernel's √K_i" (D60.md:14, 79).** The quasi-static peaks sit at 0.87–0.90·√K_i. Kernel 3's grid maximum is at 10.7 (1.040), not at 13.9 (0.976), which is nearest √K₃ = 13.13. | Say "about 10% below √K_i" |
| D3bj-14 | cosmetic | **Two wording points.**<br>• "D‴.33(iv)'s rectification, whose sign k_e reverses" (D60.md:34) has no D‴.60 data: C reports levels, and the sign reversal is round 56's run.<br>• "The slow coefficient mode" (D60.md:28): there are three, with τ = 131, 132, 166 periods. | Add a k_e row to B3, or cite the referee run; use the plural |

## Sources

- Widrow et al., "Adaptive noise cancelling: principles and applications", Proc. IEEE 63(12):1692–1716, 1975 — [snippet] https://ieeexplore.ieee.org/document/1451965
- J. R. Glover, IEEE Trans. ASSP 25(6):484–491, 1977 — [snippet] for the citation; the notch width set by the gain is [memory].
- T. E. Gibson, A. M. Annaswamy, E. Lavretsky, IEEE Access 2013, DOI 10.1109/ACCESS.2013.2284005 — [snippet] https://doaj.org/article/95064d112d7e497c9c7a41f5477e341b, https://ar5iv.labs.arxiv.org/html/1304.7278
- B. D. Riedle, B. Cyr, P. V. Kokotović, "Disturbance instabilities in an adaptive system", IEEE TAC 29:822–824, 1984 — [snippet]
- P. V. Kokotović, B. D. Riedle, L. Praly, Systems & Control Letters 1985 — [snippet] https://cas.minesparis.psl.eu/~praly/Telechargement/Journaux/1985-SCL-Kokotovic-Riedle-Praly.pdf
- L.-C. Fu, S. Sastry, "Slow-drift instability in model reference adaptive systems — an averaging analysis", Int. J. Control 45(2):503–527, 1987 — [snippet] https://scholars.lib.ntu.edu.tw/handle/123456789/332237 (the ERL M86/25 version was not fetchable here).
- Anderson et al., *Stability of Adaptive Systems*, MIT Press 1986 — [snippet] https://find.library.upatras.gr/Record/83452; content [memory].
- Kosut–Anderson–Mareels, IEEE TAC 1987 — [snippet] https://researchportalplus.anu.edu.au/en/publications/stability-theory-for-adaptive-systems-method-of-averaging-and-per/
- Bodson et al., Systems & Control Letters 1986 — [snippet] https://openresearch-repository.anu.edu.au/entities/publication/5bd59d65-27db-4a7d-8895-10b8359b9016/full
- A. I. Neishtadt, Russian Math. Surveys 2014 (with Arnold 1965 and Neishtadt 1976; resonance relation k·ω = 0) — [snippet] https://www.mathnet.ru/eng/rm9603
- "Non-Wiener effects in narrowband interference mitigation using adaptive transversal equalizers", Virginia Tech (Butterweck expansion) — [snippet] https://vtechworks.lib.vt.edu/handle/10919/26772
- From memory only: Yakubovich–Starzhinskii 1975; Wereley 1991 (harmonic transfer functions); Sanders–Verhulst–Murdock 2007; Bogoliubov–Mitropolsky; Rohrs et al. 1985; Sastry–Bodson 1989; Kapitza; Blekhman; Meerkov 1980.

## Lesson

The revision answered round 56 with two new dichotomies, and one cheap check the draft did not make settles each:
- "Bounded unless commensurate, where secular": the monodromy's spectrum (0.992–0.994 per period) makes every coherent sum saturate with the same 130–170-period constant as ⟨θ₂⟩.
- "Not single tones or pairs": rerunning the pair round 56 had already seen eject shows that it still ejects.

A control is informative where the previous round found the effect, not where it found none. A test tone at a special frequency (13.1 = 65.5ω_ref, so 2ω = 131ω_ref) needs its ordinary neighbours. Here, 0.05 rad/s away, they move a third of the "resonant" drift and the whole averaging failure from the resonance crossing to commensurability.
