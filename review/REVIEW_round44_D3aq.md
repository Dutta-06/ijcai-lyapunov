# Referee report D3aq (round 44): Proposition D‴.49, revised after round 38 (D3ak)

The revision computes the right averaged object and sweeps the knobs, but four sentences fall to one cheap check each: the low-gain
mechanism (the state never leaves the range; the kernels are ejected by the fast active block within a period), the tie between the
leak's hold and its bias (the repository's own desired regressor breaks it), the near-range thresholds (a doubled horizon, a flipped
sign), and "a smaller mass forms fewer kernels" (false at low gain, where shrinking ε makes it worse).

## What I checked and how

**Setup**
- I worked only in a scratch copy: `SCRATCH/repo`, made by tar of the working tree at HEAD 3c821aa, venv excluded.
- I ran only read-only git commands, and nothing in the repository was modified.
- Every process used the original venv's Python with OMP_NUM_THREADS=1 and `nice -n 15`, and I ran at most one heavy job at a time. None
  of my processes remain.
- Scripts are in `SCRATCH/work` and logs in `SCRATCH/logs` (`noise_runs.jsonl`, `avg_runs.jsonl`, `mut44.txt`, `q1–q8.log`), where SCRATCH
  = the session's scratchpad `review44/`.

**Objects reviewed**
- These are byte-identical to e7ed133: `closed_loop.c`, `closed_loop.py`, `exp_closed_transit.py`, `test_closed_transit_loop.py`,
  `closed_transit.md` and `closed_transit_raw.npz`.
- `FrozenOrbit`, `frozen_orbit_rhs` and `OrbitBumps` are unchanged; the uncommitted diff to `rich_transit.py` touches only `sliding_iss`.
- `D49.md` differs from HEAD only by one deleted "NUMERICS" line.

**Tests**
- `test_closed_transit_loop.py`: 4/4 pass (229 s).
- `test_rich_transit.py`: 13/13 pass.

**Numbers**
- I checked every number in (ii)–(iv) and in the numerics paragraph against the md.
- Where the md lacks a number, I checked it against the raw npz, which I dumped in full: entries, centres, coefficients and displacements
  of every A–E2 row.

**Independent numerics.** I used the repository's C integrator and `FrozenOrbit`, plus my variant `cl_var.c`. The variant is
`closed_loop.c` with three additions: a continuous dead zone, a 1st-order low-pass of the measured state and reference used by the
adaptation, and the true error. It reproduces the md's s = 0.01 row to 4 digits at rtol 10⁻⁷, and my k_e = 2 trace reproduces the B row's
six entries exactly.
- **(a) Low-gain traces.** Real loop at k_e = 2 for ε = 10⁻⁴, 10⁻⁵, 10⁻⁶, 10⁻⁷, with dense output (state range per period, path of each
  expelled kernel).
- **(b) Threshold.** f×2 at k_e = 4, 6, 8; f×½ at k_e = 2, 3 (ε = 10⁻⁵).
- **(c) Averaging check.** The real loop's one-period averages ⟨eφ_i(x)⟩ and ⟨e∂_cφ_i⟩ (2048 samples) against `FrozenOrbit.drive`: at τ = 140
  and 200 (k_e = 10), and at τ = 20, 50, 66 (k_e = 2).
- **(d) Active block.** Jacobian of the averaged active block at the fit (k_e = 2…10, amplitudes ×½, 1, 2).
- **(e) Low-gain frozen orbit.** `FrozenOrbit` at k_e = 2, 3, 4, which is below (i)'s hypothesis (Lip f = 4.56).
- **(f) Leading order and asymmetry.** Leading-order error for k_e = 5…160, and the antisymmetric part of the exact drive's Jacobian.
- **(g) Count at small masses.** k_e = 20 and ω = 1 at ε = 10⁻⁷ and 10⁻⁸.
- **(h) Part C variants.** All |w| listed; the in-range bump's sign flipped; 1.6σ and 2.9σ run to τ = 2·10⁴.
- **(i) Part D.** Fixed centres run to τ = 3·10⁴, compared with the grid's least-squares fit; WFR at 10⁻⁸.
- **(j) Shape.** A feature of width 0.2.
- **(k) Noise** (Part E setting: ε = 10⁻⁴, τ = 600): the desired regressor at s = 0, 0.01, 0.02; the mirror leak at λΛ = 5·10⁻³ and 10⁻²
  (s = 0.02), and at 3·10⁻³ with seed 11 (s = 0.01); the filtered adaptation (5 rad/s) at s = 0.01 and 0.02; a dead zone of 0.03 at
  s = 0.01; s = 0.005 with seed 11.

**Mutation testing**
- 11 of the 15 planned mutants of the current code were run (`mutate44.py`), each against the tests that could kill it.
- The other 4 (C3–C6), plus W1, survive by construction; the reason is given under Tests.

**Literature**
- Web search confirmed: Sadegh–Horowitz, IJRR 9(3) 1990 (desired-trajectory adaptation, which removes noise correlation); Peterson–
  Narendra, IEEE TAC 27 (1982), 1161–1168 (dead zone); Riedle–Kokotović, "Integral manifolds of slow adaptation" (Riedle's 1986 thesis).
- From memory: Anderson et al. 1986 (nonlinear averaging along the frozen system).

## Verdicts

### (i) The averaged loop
- **Contraction.** Correct as a sufficient condition. At k_e = 2, 3, 4 (below Lip f) the Picard iteration still converges (137, 52, 41
  iterations), and the orbit is unique and attracting: five starts converge, and the Floquet exponent is −3.39 at k_e = 2 although
  sup r′ = 3.96 > k_e. Those B rows therefore lie outside (i)'s hypothesis, which the md never says.
- **Leading-order error.** Confirmed O(k_e⁻²): k_e²·max|g_lead − g| = 0.106–0.119 for k_e = 5…160 (relative error 16% at k_e = 10).
- **"Not a gradient."** Confirmed: the antisymmetric part of ∂g/∂w is 12.6% (w = 0) and 4.9% (mid-transit) of the symmetric part.
- **Averaging for the inactive kernels while the active block is fast: not justified.** The duals integrate eφ_j along the actual orbit,
  which is the frozen orbit only if the active block is frozen over the period. The proof sketch's "small for the inactive kernels (their
  masses)" addresses the wrong variables.
- **Measured at k_e = 10.** Active centres swing 0.041–0.056 per period. The five largest inactive drives match `FrozenOrbit` to ≤ 4·10⁻³
  relative; kernels next to the active block deviate by up to 4% (τ = 140) and 13% (τ = 200) of the largest inactive drive.
- **Measured at k_e = 2.** Swings are 0.20, then 0.37 and 0.30. The real error differs from the frozen orbit's by 0.059–0.11, i.e. 17–40%
  of max|e|; the competitors' drives are off by 0.9–6.5%, the size of the near-tie margins.
- **Rates of the averaged active block at the fit.** k_e = 10: −2.70, −2.08, −0.98, i.e. 4.9–13.5ω (the draft's "5–13" is reproduced);
  k_e = 2: −8.9, −7.3, −2.8, faster than the error filter itself.

### (ii) The real loop
- **Base-point numbers.** Right: entries within 0.363 at ε = 10⁻⁴…10⁻⁶.
- **"Same coefficients" at 10⁻⁴: false.** At the same τ = 600 (Part E's noise-free row against A's averaged row) the pair splits are
  −0.5351/−0.0649 vs −0.5176/−0.0824, and 0.7199/0.0306 vs 0.7035/0.0464; kernel 16 sits at 0.6314 in the real loop vs 0.6101 averaged:
  the averaged pair has merged, the real one has not. Part A also compares the real loop at τ ≈ 337 with the averaged loop at τ = 600.
- **Over which knobs (raw npz; the B table omits the averaged entries).** Entries agree within 0.13–0.39 at k_e = 4–10 and a = 0.6, 0.8;
  they differ by 0.54 at ω = 0.05, 1.6 at ω = 1, and 2.5 at k_e = 20; kernel sets differ at k_e = 2 and 3; k_e = 40 compares horizons of
  600 vs 720 (both mid-transit); the quasi-periodic reference has no averaged loop. So the Status line "agreement over gains … reference
  shapes" overclaims.
- **"One kernel per feature at small masses."** Supported where averaging holds, but only by my runs: ε = 10⁻⁷ and 10⁻⁸ give three
  kernels at k_e = 20 (kernels 5, 11, 16) and at ω = 1 (5, 10, 17). There is no such run in the md, the margins are never computed, and
  D‴.51 (an open-loop result) is not shown to transfer.
- **No mass restores it at low gain.** At k_e = 2 the real loop forms 5, 6, 6, 6 kernels at ε = 10⁻⁴, 10⁻⁵, 10⁻⁶, 10⁻⁷, with 1, 3, 3, 4
  expelled. At 10⁻⁷ no kernel is left at the feature 0.07 by τ = 82.
- **Unstated shape hypothesis.** "At the feature's centre and amplitude, then tracks exactly" assumes f is a sum of kernel-width
  Gaussians. A middle feature of width 0.2 is formed by two kernels (at −0.054 and 0.185), with orbit-error rms 1.5·10⁻³ (ε = 10⁻⁶,
  τ = 600). With a 1.6σ edge feature, the in-range feature at 0.07 ends as three coincident kernels with w = −1.62, −0.20, +1.22.

### (iii)(a) Low gain
- **The mechanism is false.** At k_e = 2 (every ε) the state stays in [−1.027, 1.0025]; the per-period maximum is ≤ 1.0025. The frozen
  orbit even stays in [−0.60, 1.00]. The stranded kernels were formed inside the range, at the features: kernel 17 at 0.61 (τ = 41.1), 10
  at 0.07, 16 at 0.61. They were ejected while the state was inside: kernel 17, at 0.63 with w = 0.81, went 0.578 → 1.345 within the
  single period at τ = 63.6, while the entrant (kernel 10) swung by 0.40–0.53 per period. The exact averaged loop, which contains the
  orbit's error exactly, expels nothing. The cause is the active block's within-period dynamics.
- **The threshold belongs to this instance.** f×2 expels at k_e = 4 (kernel 17 to 1.50, w = 1.56) but not at 6 or 8; f×½ forms exactly
  three kernels at k_e = 2 and 3. The ratio λ_c/k_e is not a sharp criterion either: 2.4 expels (×1, k_e = 3), 2.9 does not (×2, k_e = 6).
- **Quoted numbers.** |x| ≈ 1.47–1.51 and |w| = 0.51–0.82 are right, but they appear only in the raw npz.

### (iii)(b) Noise and leak
**What does not hold up:**
- "Keeps the kernels in the range up to 1% only at λΛ = 3·10⁻³" is contradicted by the md's own E2 row: λΛ = 10⁻³ at s = 0.01 pushes out
  0 kernels.
- "Tied" holds only within the one-parameter mirror-leak family, where it is trivially true.
- The repository's own desired regressor (D‴.31(iv), D‴.33; already a flag in `closed_loop.c`, unused) breaks it. At s = 0.01 it forms
  exactly three kernels at (−0.527, 0.058, 0.648) with w = 0.93, −0.61, 0.765.

**Remedy comparison at s = 0.01:**

| remedy | kernels formed / expelled | in-range fit rms | true error |
|---|---|---|---|
| desired regressor | 3 / 0 | 0.059 | 0.031 |
| mirror leak, λΛ = 3·10⁻³ (seed 11) | 3 / 0 | 0.26 | 0.058 |
| filtered adaptation, 5 rad/s | 9 / 5 | 0.011 | 0.023 |
| dead zone 0.03 | 3 / 0, never forms the bump at 0.07 | 0.20 | 0.063 |

**What holds up:**
- **"Every leak fails at 2%"** rested on a single leak strength. It survives my runs: λΛ = 5·10⁻³ leaves one kernel, which travelled from
  −0.6 to 0.63 (w = 0.16); λΛ = 10⁻² forms nothing; the desired regressor (12 kernels, 11 expelled) and the filtered adaptation (14
  kernels, all expelled) also fail.
- **Seed robustness:** seed 11 reproduces seed 7 at s = 0.005 and at s = 0.01 with λΛ = 3·10⁻³.

**Misstated numbers:**
- "0.67 and 0.8 of the largest" are max|w|; the fractions are 0.74 and 0.89.
- "Halved": the amplitudes are 0.48, 0.27 and 0.44 of the true ones, and the centres move by 0.055 and 0.097 (0.8σ).
- "Five at the wrong places": the final centres are 0.61 (three kernels), 0.155 and 0.07, so the −0.52 feature is the one missing.
- The reported errors are measured errors (x̃ − x_m). For example, 0.101 measured is 0.0965 true.

### (iii)(c) Features near the range
- **Horizon artefact.** The 1.6σ bump is formed at τ = 5344, just past Part C's horizon. By τ = 2·10⁴ the fit is exact (rms 7.6·10⁻⁷), so
  the "compensating pair" was a transient.
- **Instance artefact.** With the bump at 0.07 flipped to +0.6: the 0.92σ bump is the first entrant (τ = 49.7, against 1236.7); the 1.6σ
  bump forms at τ = 59.1 with a clean fit.
- **2.9σ.** Not formed by τ = 2·10⁴ (with the flip, not by 5000), so this part holds.
- **0.92σ.** The fit hides a fourth kernel with w = 0.0085, which is 0.94 of the formation threshold.
- **Averaged loop only.** Part C no longer has a real-loop run.

### (iv) Transit time and price
- **τ₁/k_e numbers.** Right.
- **"A smaller mass forms fewer kernels."** Fails at k_e = 2 (5, 6, 6, 6 kernels).
- **Fixed centres.** "Leave the grid's residual" is false: the six-kernel state is the least-squares fit on those six grid kernels (orbit
  error 1.228·10⁻³, reproduced); the full grid's residual is 8.2·10⁻⁵ (orbit error 7.4·10⁻⁶); kernels keep entering at τ = 9713, 20156,
  22795, 22813 and 27765: 11 kernels by τ = 3·10⁴, orbit error 9.2·10⁻⁴.
- **WFR at 10⁻⁸.** A sub-threshold kernel with w = −0.0058 sits at 0.07.

### Tests
**Killed:** O5 (OrbitBumps mean dropped — a round-38 survivor, now killed); F1 (FrozenOrbit mean dropped); F2 (drive's regressor evaluated
at x_m); F6 (mobility with twice the mass); C8 (leak sign flipped) — killed only by a 900 s hang, not by an assertion.

**Survive:** F4: `frozen_orbit_rhs`'s WFR branch (Part D); F5: `frozen_orbit_rhs`'s fixed branch (Part D); C1: C's WFR (unused by the
experiment); C2: f̂ evaluated at the true state (Part E); C7: entry times not interpolated (every real-loop entry).

**Survive by construction** (no test has noise, a₂ ≠ 0, regx, reads e_last, or calls tones): C3: the error without the noise; C4: the
desired-regressor flag (regx) ignored; C5: the quasi-periodic reference's second tone dropped from ẋ_m; C6: e_last taken over three
periods instead of one; W1: `tones` rms off by √2; every function of `exp_closed_transit.py` — no test imports it.

### Round-38 issues

| issue | status | note |
|---|---|---|
| D3ak-1 | fixed | except at low gain |
| D3ak-2 | partly | |
| D3ak-3 | partly | mechanism wrong; threshold instance-specific; (iv) unqualified |
| D3ak-4 | partly | |
| D3ak-5 | partly | |
| D3ak-6 | partly | |
| D3ak-7 | partly | |
| D3ak-8 | fixed | |
| D3ak-9 | partly | ε swept, horizon not |
| D3ak-10 | partly | Part C e_last and the error history still missing |
| D3ak-11 | partly | "the others left in place" in New; the real displacement 0.205 at k_e = 3 is not in the md |
| D3ak-12 | fixed | |
| D3ak-13 | fixed | but the numerics run outside the hypothesis |
| D3ak-14 | partly | |
| D3ak-15 | not | deferred to D‴.52 |
| D3ak-16 | fixed | τ₁ notation remains |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3aq-1 | major | (iii)(a)'s mechanism is false. The state stays in [−1.027, 1.0025] while kernels formed at the features are ejected within one period (17: 0.578 → 1.345). The exact averaged loop expels nothing. | (iii)(a): restate as an averaging failure of the fast active block. Part B: report the transit-time swing and the state's range. |
| D3aq-2 | major | The averaging's small parameter is not ε. Expulsion persists and worsens as ε shrinks (1, 3, 3, 4 kernels expelled at 10⁻⁴…10⁻⁷). The threshold moves with amplitude (×2 expels at k_e = 4; ×½ is clean at 2). The active rates exceed k_e there, and the proof sketch averages the wrong variables. | Proof sketch: average with the fast subsystem (e plus the active block), state the rest condition. (ii)/(iii)(a): give the condition, or say "this instance, k_e ≥ 4". |
| D3aq-3 | major | "Tied" is refuted by the desired regressor at 1%: fit rms 0.059 vs 0.26, true error 0.031 vs 0.058. "Only at 3·10⁻³" is contradicted by the md's own E2 row. "Every leak fails at 2%" rested on one strength. | (iii)(b), Status: restrict "tied" to the mirror leak. Part E3: add the DR, filter, dead-zone and leak 5·10⁻³ and 10⁻² rows. Cite D‴.31 and D‴.33. |
| D3aq-4 | major | (iii)(c)'s thresholds are horizon- and sign-specific. The 1.6σ bump forms at τ = 5344 and the fit becomes exact; with the neighbour's sign flipped, 0.92σ enters first (49.7) and 1.6σ at 59.1. | (iii)(c): state it in terms of the drive along the orbit with the horizon given; Part C: longer horizon and a sign sweep. |
| D3aq-5 (D3ak-7 cont.) | major | Surviving mutants: F4, F5, C1, C2, C7, plus by construction C3–C6, W1 and all experiment functions. | Tests: C loop with noise, WFR and a₂ against Python; e_last; entry interpolation; frozen_orbit_rhs WFR and fixed branches; the experiment's `averaged`, `crossings` and `report`. |
| D3aq-6 | minor | Coefficients at 10⁻⁴ differ at equal τ (0.016–0.018; kernel 16 at 0.631 vs 0.610). Part A compares τ = 337 with τ = 600. | (ii); Part A at equal τ; claim agreement only for sums over coincident kernels. |
| D3aq-7 | minor | The B table lacks averaged entries, coefficients and displacements, and compares unequal horizons (k_e = 40: 600 vs 720). Entry agreement over the knobs is up to 2.5 in τ, not 0.4. The draft quotes numbers that exist only in the raw npz. The Status overclaims. | B table; Status. |
| D3aq-8 | minor | (ii) needs f to be a sum of kernel-width Gaussians: the width-0.2 feature takes two kernels and leaves a residual of 1.5·10⁻³. | Setting / (ii): state the hypothesis. |
| D3aq-9 | minor | "A smaller mass forms fewer kernels" fails at k_e = 2. | (iv) |
| D3aq-10 | minor | "The grid's residual" is a τ = 3000 plateau of the six-kernel fit; 11 kernels by τ = 3·10⁴. | (iv), numerics D. |
| D3aq-11 | minor | The averaged loop is used at k_e = 2–4, below Lip f, and the md caption treats convergence of the iteration as validity. | md caption, (i). |
| D3aq-12 | minor | Noise numbers misstated: "0.67 and 0.8 of the largest", "halved", "wrong places", and measured vs true error. max\|e\| hides partial fits. | Numerics E: add the true error and the in-range fit rms. |
| D3aq-13 | minor | Neither the swing nor the active rates are in the md; the "swing" column is measured after the fit (5·10⁻⁴). | Part B. |
| D3aq-14 | minor | Counts depend on the 1% threshold (sub-threshold kernels at 0.0085 and −0.0058). | Report the sub-threshold mass. |
| D3aq-15 | minor | "New: the exact averaged loop" is classical frozen-system averaging; Sadegh–Horowitz and Peterson–Narendra are uncited. | Status. |
| D3aq-16 | minor | Part C is averaged-loop only; the e_last and error history D3ak-10 asked for are still missing. | Part C, (iv). |
| D3aq-17 | cosmetic | τ₁ means two things; "within 0.003" is 0.0032; the quasi-periodic averaged loop is not "undefined"; "the two loops alike" is applied to a case with no averaged loop. | (ii), (iv). |

## Lesson

The revision computed the right averaged object and swept the knobs. But it read each sweep through a mechanism it had not traced, and
turned the edge of each sweep — one leak family, one horizon, one sign pattern of f, one amplitude — into a general statement. Four cheap
checks each overturn a sentence: one traced event (kernel 17 leaves while the state stays inside [−1, 1]); one remedy from another family
(the repository's own desired regressor); one doubled horizon (the 1.6σ feature forms at τ = 5344); one flipped sign (the 0.92σ feature
becomes the first entrant). For an averaging claim, find the small parameter and test it where it fails. Here it is the active block's
speed against k_e and ω, not the mass, and shrinking ε makes low gain worse, not better.
