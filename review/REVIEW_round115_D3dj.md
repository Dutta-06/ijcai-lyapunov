# Referee report D3dj, round 115: Proposition D‴.61 as revised after round 110, with part VC

## Scope

**Commit and setup**
- I reviewed 5be79c1 on claude/vigilant-clarke-hf939b.
- W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref115.
- I exported the commit with `git archive 5be79c1 | tar -x -C W/repo`. The sha256 of exp_gh_peano.py, test_gh_peano.py, D61.md, gh_peano.md and gh_peano_raw.npz all match `git show 5be79c1:…`.
- W/repo/venv is a symlink to the repository's venv (python-flint 0.9.0, numpy 2.5.4, scipy 1.18.1).
- GP_OUT and WM_OUT pointed at W/out, which holds copies of gh_peano_raw.npz and width_mismatch_raw.npz.
- Every job ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`. Pytest always ran with `--basetemp` inside W.
- My mutant runner, W/tools_ref/mutants_bt.py, is tools/mutants.py with REPO fixed to W/repo, a per-slot basetemp and WM_OUT added.
- At most two heavy processes of mine ran at a time.
- I used only read-only git commands, and `git status` is clean. I touched no process of anyone else. The only thing I stopped was my own monitor task, with TaskStop.

**Read in full**
- review/drafts/D61.md.
- experiments/exp_gh_peano.py:
  - every hunk since 88ef79f: `c_need_ball`, `d_lines`, `v_lines`, part VC (`corner_psi`, `corner_floor`, `corner_phi_dom`, `corner_sector`, `part_VC`, `vc_lines`), `report`'s DB merge, KEYS, `main`;
  - the functions they call: `cert_E`, `gap_closed`, `bound_sup`, `box_bound_ball`, `dom_bb`, `r_star`, `ratio_floor`.
- results/chaos/gh_peano.md.
- tests/test_gh_peano.py, the round-110 and part-VC tests line by line.
- tools/mutants/round110_D61_own.py and partVC_D61_own.py.
- review/REVIEW_round110_D3de.md and the Round 110 entry of review/REVIEW_LOG.md.
- PAPERS.md Paper 18 and HANDOFF.md §§1, 5, 7, 8, 10.
- `git diff 88ef79f 5be79c1` of these files.

## Summary

**Design.** Round 110's fixes are in:
- `c_need_ball`: the radius each printed threshold r* needs, rechecked in Arb at r's 4-decimal value;
- `d_lines` and `v_lines`: no fallback to the angle grid's radius, a guard when no branch and bound is verified, and "both … in ball arithmetic" only when the ball bound is negative;
- `report` keeps D's own bb when it is in balls;
- a dense-supremum test of `box_bound_ball` on 61 boxes;
- 13 new mutants.

Part VC is new. With ε = 2 − |η|, θ = arg η = εφ and z = ε(x + iy):
- `cert_E` = ε³Ψ_φ + o(ε³), where Ψ_φ = −(x² − y²) − 4φxy + (x − ½)₊[(y + φ)² − (x − ½)²/3];
- the valley floor is y = min(2φx, φ);
- φ_dom = 1.4397;
- part VC prints the ratio of the true bound to ε³Ψ_φ along the floor, and the straight segment's sector over ε.

D61's Status calls VC "a candidate route, its remainder not yet bounded uniformly". The md says it becomes "a candidate for P on 1.995 < |η| < 2, once the expansion's remainder is bounded uniformly in (ε, φ) and the floor joined to a path to η".

**Verdicts**
- **Reproducibility is exact.**
  - The md regenerates byte for byte (sha256 6394a131…71fec6ae6b).
  - 38 tests pass in 17.00 s.
  - `part_VC()` reruns in 0.016 s, identical to the raw.
  - My sample of D, all of DB and V bands 0–1 with the ball recheck are identical to the raw.
- **Part VC's mathematics is correct as a pointwise limit.**
  - I re-derived Ψ_φ.
  - The floor, its values, φ_dom as the root of 8φ³ − 12φ² + 1 = 0 (equivalently (φ − 1)²(2φ + 1) = ¾) and the −¼ at φ = 0 all check out.
  - The straight segment of VC is exactly D's domination gap, and the printed ratios and sectors match the raw.
- **"Not a proof" is accurate.** What VC says remains is not. Two things are wrong with it:
  - The remainder is not uniform in φ. The relative correction near the saddle is −ε(2φ² + O(1)).
  - On the true bound, VC's floor turns positive at θ ≈ √(ε/2), where a third limit problem (θ ~ √ε) governs. So no bound "uniform in (ε, φ)" can exist, and the floor is not the path there.
- **Round 110's fixes work.** All the shipped mutants they target are killed on assertions, not crashes. Part VC's own computation is never recomputed in a test, though: 7 non-equivalent mutants of `part_VC` and `corner_sector` survive.
- **No major issue.** No proved claim changes.

**Findings**
1. **D3dj-1 (minor; the main finding).** VC's account of what remains understates it.
   - Measured coefficients: (E/ε³ − Ψ_φ)/(εΨ_φ) → −0.249, −1.002, −2.652, −8.709, −50.63, −200.26 at φ = 0, 0.5, 1, 2, 5, 10 (x = 0.05).
   - At ε = 10⁻³, the true bound on VC's floor is positive from φ ≈ 25 (θ = 0.025: E = +1.54·10⁻⁹, while min_y E = −2.78·10⁻⁹).
   - At θ = √εψ: E/ε² → Ψ̃_ψ, a different limit (formula below), which I verified to 6 digits.
   - The route therefore needs at least four parameter charts, not one uniform remainder bound.
2. **D3dj-2 (minor). Tests.**
   - No test recomputes `part_VC` or `corner_sector`. The survivors C01, C02, C03, C05, C06, C07 and C09 each change the raw and the md's VC line on a rerun.
   - C14 survives the 10⁻³ tolerance even though it moves φ_dom to 1.4390.
   - The planted value 2.5 in the round-110 guard test cannot tell the radius apart from 2 (R04–R06 survive).
3. **D3dj-3 (cosmetic).** "On the straight segment Ψ's maximum is −0.2500" is the maximum beyond the disc (x ≥ ½). On the whole segment it is 0, at the saddle. `corner_phi_dom`'s docstring has the same slip.
4. **D3dj-4 (cosmetic).** Floating point cannot resolve the ε³ margins below ε ≈ 10⁻⁴.
   - `corner_sector` at |η| = 1.99999 gives lo/ε = 1.2485.
   - `cert_E`/ε³ is off by up to 1.2·10⁻⁴ (relative) at ε = 10⁻⁴.
   - This affects no printed number, but it bounds any extension of VC.
5. **D3dj-5 (cosmetic).** HANDOFF still has stale lines:
   - §7: "not yet in any md or draft", "Put … as a part (VC)", "A probe (not in any md yet)", "Re-reviews due: … D‴.61 (after round 90)";
   - §10: "test_gh_peano.py 35" (38 now), "Last commits", "Referee ids";
   - §5: "cover the overlap with the straight segment (φ < φ_dom)", which is not needed.
6. **D3dj-6 (cosmetic). Odds and ends.**
   - `vc_lines` prints VC_X and groups the rows by VC_EPS from the module's constants, not from the raw.
   - A zero Ψ at a floor point would divide by zero (my C04 crashes).
   - The dense-supremum test's tolerance is 10⁻¹², against measured float noise of 1.2·10⁻¹⁶.
   - Round 108's acceptance threshold (`u < 0` in `dom_bb`) is not tested by a planted positive ball bound (my R10 survives).

## What I checked and how

**(1) Reproducibility** [own][md][raw]
- **The md.** `--from-raw` into W/out matches committed gh_peano.md byte for byte (cmp; both sha256 6394a13157033431702d9b9a1d6b4227f4b38fcaced9898a80488d71fec6ae6b).
- **Tests.** `pytest -q -p no:cacheprovider --basetemp=W/pt tests/test_gh_peano.py`: 38 passed in 17.00 s (W/logs/pytest.txt).
- **Part VC** (W/own/rerun_vc.py). `GP.part_VC()` gives the same 12 rows × 7 points as the raw (x, y, E, ψ all identical), φ_dom 1.4396926208478287, gap0 −0.25, and sector lo 0.014842092994976668, 0.007307117029307617, 0.0014439590610928607.
- **The printed ratios are verbatim.** For example, at ε = 0.01 and φ = 5 the range is 0.4958–1.1036, printed as "0.496–1.104". D61's "(at ε = 0.01 0.913–1.053 for φ ≤ 2)" is the union over φ ≤ 2.
- **Sample of D, DB and V** (W/own/sample_dbv.py, W/logs/sample_dbv.log).
  - D: the axes' gaps [−2.5021917871015376·10⁻¹⁰, −5.002502190265545·10⁻⁷] and `dom_sector()` (0.08605763030080119, 0.6175580488484276, Re 1.0051079…–1.0999465…) are identical to the raw. `dom_radius` at 13 of the 361 angles is identical too.
  - DB: `part_DB()` is identical at both radii, (4935, −3.281257330289833·10⁻⁷, 0, −3.2812573349693787·10⁻⁷) and (402664, −1.9432233600014115·10⁻¹¹, 0, −1.9432017877494493·10⁻¹¹), in 26.5 s.
  - V: bands 0 and 1 with ball are dict-identical (9.3 s and 46.4 s).

**(2) Part VC's mathematics** [own]
- **The centred form.** `cert_E` equals F_η(σ) − F_η(η/2), built from the definitions (log|η − σ| + Re(η − σ)²/2 + log|σ| + Re σ²/2 + Λ⁺ outside the unit disc). On 3000 random (η, z) with 1 < |η| < 2, the largest relative difference is 6.9·10⁻¹⁶, computed in Arb (W/own/blowup.py).
- **Derivation of Ψ_φ** (by hand).
  - log|1 − z²| + Re(cz²) = Re((c − 1)z²) + Re(log(1 − z²) + z²).
  - c − 1 = ε(−1 + 2iφ) + ε²(¼ − 2φ² − 2iφ) + O(ε³), which gives ε³[−(x² − y²) − 4φxy] + O(ε⁴); the log term is −ε⁴Re w⁴/2 + ….
  - Write δ = |σ| − 1 = ε(x − ½) + ε²(y² − x)/2 + O(ε³) and ρ′ = Re σ − 1 = ε(x − ½) + ε²(−φ²/2 − φy − x/2) + O(ε³). Then δ − ρ′ = ε²(y + φ)²/2 + O(ε³).
  - C = h(1 + δ)1_{δ>0} − ρ′₊², with h(1 + δ) = δ² − δ³/3 + δ⁴/4 − …. Since t ↦ t₊² is C¹ with a Lipschitz derivative, this gives C = ε³(x − ½)₊[(y + φ)² − (x − ½)₊²/3] + O(ε⁴).
  - The O(ε⁴) also covers the O(ε)-wide layer around x = ½, where δ > 0 ≥ ρ′; there C = O(ε⁴). So the limit is right pointwise, and uniform on compact sets of (w, φ).
- **The floor.**
  - Ψ_φ is convex in y, with y²-coefficient 1 inside the disc and 1 + (x − ½) outside.
  - The minimiser is y = 2φx inside. Outside, ∂_yΨ = 0 gives y(1 + u) = φ(1 + u), so y = φ.
  - The values are −(1 + 4φ²)x² and −x² − φ² − (x − ½)³/3. The md is correct.
  - From b = (x_b, 0), x_b < ½, the vertical segment up to the floor is negative by convexity (≤ −x_b²). So "joined to b" is trivial in the limit.
- **φ_dom.**
  - The maximum of −(u + ½)² + φ²u − u³/3 over u ≥ 0 is at u = φ − 1 for φ > 1. Its value is (2/3)φ³ − φ² + 1/12, which vanishes when 8φ³ − 12φ² + 1 = 0, i.e. (φ − 1)²(2φ + 1) = ¾.
  - Exact bisection gives 1.4396926207859084. The script's value is 6.2·10⁻¹¹ above it, and it prints 1.4397 correctly. At φ ≤ 1 the maximum is −¼ at u = 0.
- **VC's straight segment is D's domination gap** (W/own/vc_more.py). The maximum of `cert_E` over the real z-segment beyond the disc equals `gap_closed(R, θ)`. On 32 (R, θ) pairs with a 10⁴-point grid the difference is at most 1.4·10⁻⁹; at |η| = 1.999, θ = 0 the two values are −2.502189·10⁻¹⁰ and −2.502192·10⁻¹⁰.
- **The sector over ε** (`corner_sector`):

  | \|η\| | 1.95 | 1.99 | 1.995 | 1.999 | 1.9995 | 1.9999 | 1.99999 |
  |---|---|---|---|---|---|---|---|
  | lo/ε | 1.72115 | 1.48421 | 1.46142 | 1.44396 | 1.44182 | 1.44023 | 1.24847 |
  | (lo/ε − φ_dom)/ε | 5.63 | 4.45 | 4.35 | 4.27 | 4.26 | 5.40 | −19122 |

  So θ_dom = φ_dom·ε + ≈4.3ε² down to ε = 5·10⁻⁴. Below that, float noise takes over (D3dj-4). The printed 1.4842, 1.4614 and 1.4440 are right.

**(3) Is the remainder uniform, and what a proof on 1.995 < |η| < 2 needs** [own]
- **The remainder is O(ε) at fixed (w, φ)** (Arb at 200 bits, ε = 10⁻⁴, 10⁻⁶, 10⁻⁸; W/logs/blowup.log). (E/ε³ − Ψ_φ)/ε converges, but its size relative to Ψ_φ grows like φ²:

  | φ | 0 | 0.5 | 1 | 2 | 5 | 10 |
  |---|---|---|---|---|---|---|
  | x = 0.05 | −0.249 | −1.002 | −2.652 | −8.709 | −50.63 | −200.26 |
  | x = 0.25 | −0.219 | −1.062 | −2.694 | −8.425 | −47.84 | −188.47 |

  - Near the saddle E/(ε³Ψ_φ) = 1 − ε(2φ² + O(1)) + ….
  - Far out the relative correction is O(ε|w|), i.e. O(|z|): +4.97 at x = 8, φ = 0.
  - The expansion is therefore valid only for θ ≪ √ε and |z| ≪ 1. The 0.496 that the md prints at ε = 0.01, φ = 5 is exactly 1 − 2εφ².
- **The floor fails on the true bound beyond θ ≈ √(ε/2)** (W/logs/floor_large_phi.log; min_y over y ∈ [0, 3·floor]):

  | ε | φ (θ) | ratio at x = 0.05 | true E on the floor | min_y E |
  |---|---|---|---|---|
  | 10⁻³ | 20 (0.020) | +0.202 | −8.07·10⁻¹⁰ | −2.22·10⁻⁹ |
  | 10⁻³ | 25 (0.025) | −0.247 | +1.54·10⁻⁹ | −2.78·10⁻⁹ |
  | 10⁻³ | 30 (0.030) | −0.795 | +7.16·10⁻⁹ | −3.21·10⁻⁹ |
  | 10⁻⁴ | 50 (0.005) | +0.501 | −1.25·10⁻¹¹ | −1.67·10⁻¹¹ |
  | 10⁻⁴ | 100 (0.010) | −0.995 | +9.95·10⁻¹¹ | −3.33·10⁻¹¹ |

  - At ε = 10⁻³, φ = 25 the floor is positive at x = 0.25 and 0.45 too (+2.69·10⁻⁸, +3.99·10⁻⁸), while min_y E stays negative.
  - These θ lie inside the domination's violated sector at |η| = 1.999 (at least [0.004, 0.764]; D).
  - The true valley does exist there: the exact quadratic's floor, Im z = (Im A/|Re A|)Re z with A = η²/4 − 1, gives −|A|²/|Re A|·(Re z)². VC's floor, which takes A ≈ ε(−1 + 2iφ), is the wrong path there.
- **A distinguished regime at θ ~ √ε.** Put θ = √εψ and z = εx + i√εy.
  - Then Re A ≈ −ε(1 + 2ψ²), Im A ≈ 2√εψ, |σ| − 1 ≈ εD and Re σ − 1 ≈ εP, with D = x − ½ + y²/2 and P = x − ½ − ψ²/2 − ψy. This gives
    E/ε² → Ψ̃_ψ(x, y) = (1 + 2ψ²)y² − 4ψxy − y⁴/2 + D₊² − P₊².
  - Checked in Arb (W/logs/mid_regime.log). At ψ = 1, (x, y) = (1, 0.5): E/ε² = −0.890684, −0.890626, −0.890625, −0.890625 at ε = 10⁻⁴ … 10⁻¹⁰, against Ψ̃ = −0.890625. At ψ = 3, (5, 3): 31.513697 → 31.500000, against 31.5.
  - Matching between the regimes:
    - Ψ_φ(x, φỹ)/φ² (φ → ∞) and Ψ̃_ψ(x, ψỹ)/ψ² (ψ → 0) share the limit Ψ_∞(x, ỹ) = ỹ² − 4xỹ + (x − ½)₊(ỹ + 1)². Its floor is ỹ = min(2x, 1), with values −4x² and −1; this governs ε ≪ θ ≪ √ε.
    - Ψ̃_ψ(x, ŷ/ψ) (ψ → ∞) tends to 2ŷ² − 4xŷ + (x − ½)₊². Its floor is ŷ = x, with value −2x² + (x − ½)₊². That is HANDOFF §7's fixed-θ limit as θ → 0, valid for θ ≫ √ε.
  - So the parameter corner (ε, θ) = (0, 0) needs a weighted blow-up with weights (2, 1) for (ε, θ), plus a second blow-up at ψ = 0 (the φ chart). VC covers one chart of four.
  - In the Ψ_∞ and Ψ̃ charts the real segment through b sits on the zero level at leading order (Ψ_∞(x, 0) = Ψ̃_ψ(x, 0) = 0 inside the disc). There b's margin is the lower-order −εX², so the path must leave y = 0 at once: Ψ̃ < 0 for 0 < y < 4ψx/(1 + 2ψ²).
- **The rest of a proof along route (2).**
  - What is already proved for every η: 0 → a → b on the real z-segment inside the disc, bound ≤ log(1 − x²) + x²Re c ≤ −(1 − Re c)x² − x⁴/2 < 0.
  - So the saddle's degeneration φ″ = 2 − 8/η² → 0 is harmless on this route. It enters only condition (1)'s Taylor ratio q, and Laplace's main term is needed only on compact subsets of |η| < 2.
  - No overlap with the straight segment is needed: at φ = 0 the floor is the straight segment, and φ_dom only explains why D fails.
  - In z there are three pieces:
    - an inner chart, |w| ≤ W;
    - a matching chart, εW ≤ |z| ≤ δ₀, where E ≈ −X³/3. Certify E/X³ with s = ε/X ∈ [0, 1/W] as an interval;
    - an outer chart, |z| ≥ δ₀, by perturbing η = 2e^{iθ} (at θ = 0, log(1 − X) + X + X²/2 ≤ −X³/3), plus V's end piece.
- **Ball arithmetic with ε ∈ [0, ε₀] as an interval works for the φ chart.** My prototype is W/own/chart1_ball.py.
  - It uses this cancellation-free form of E/ε³:
    - Re(αw²) with α = (1 − ε/2)²·2iφ·E₁ − 1 + ε/4, where E₁ = (e^{iτ} − 1)/(iτ) = sinc τ + i(τ/2)sinc²(τ/2) and τ = 2εφ;
    - the log term −εRe(w⁴G(ε²w²)), with G(t) = Σ_{k≥2}t^{k−2}/k plus a tail bound;
    - D = (2x + ε|w|²)/(|1 + εw| + 1) − ½|1 + εw|;
    - P = (1 − ε/2)[−εφ²·½sinc²(εφ/2) + x cos εφ − y sin εφ] − ½;
    - I = (1 − ε/2)[φ sinc(εφ)(1 + εx) + y cos εφ] and J = I²/(|σ| + Re σ);
    - C/ε³ = D₊³κ(εD) + μ(D₊ + P₊), with κ(t) = Σ(−1)^{k+1}t^k/(k+3) and μ = (D₊ − P₊)/ε = min(D₊/ε, J) ∈ [0, J]; μ = J when P > 0.
  - It agrees with the centred form in Arb to 1.0·10⁻¹⁵ (relative). It equals `corner_psi` at ε = 0 to 3.6·10⁻¹⁵.
  - Ball enclosures with ε ∈ [0, 0.005], on a w-box of half-width 0.01 around the floor:

    | φ | x | enclosure, ε ∈ [0, 0.005] | Ψ at the centre | enclosure, ε ∈ [0, 0.001] |
    |---|---|---|---|---|
    | [1, 1.01] | [0.99, 1.01] | [−2.3103, −1.7877] | −2.0517 | |
    | [1, 1.01] | [0.49, 0.51] | [−1.3701, −1.1079] | | |
    | [1, 1.01] | [0.44, 0.46] | [−1.0968, −0.9301] | | |
    | [5, 5.01] | [0.99, 1.01] | [−40.55, −13.00] | | |
    | [10, 10.01] | [0.44, 0.46] | [−117.63, +35.41] | | [−88.62, −60.78] |
    | [10, 10.01] | [0.99, 1.01] | [−412.40, +64.95] | | [−145.99, −61.14] |

  - So the φ chart closes in balls for φ ≲ c·ε₀^{−1/2}, which is about 7 at ε₀ = 0.005. Each of the other charts needs its own cancellation-free form with its small parameters as intervals that include 0: (ε/θ, θ/√ε), (√ε, ψ) and (ε/θ², θ). One path per parameter box suffices, since the valley condition is pointwise in η.
  - Precedent: a thesis on rigorous numerics for a singular-perturbation problem covers ε ∈ (0, ε₀] by factoring ε out before the interval check [snippet]. Blow-up charts for degenerate fast-slow points [snippet].
- **Float limit.** The float `cert_E`/ε³ differs from Arb by up to 1.19·10⁻⁴ (relative) at ε = 10⁻⁴ (W/logs/chart1_ball_check.log). Any extension of VC to smaller ε must be done in balls.

**(4) Every new or changed claim against the md and raw** [md][raw]
- **Verbatim and earned:**
  - D61's (VC) numerics: 1.000–1.005 … 0.949–1.010; 0.913–1.053; −0.2500; φ_dom = 1.4397; 1.4842, 1.4614, 1.4440;
  - the Status's Ψ_φ and floor;
  - the c_need floors 1.89269820, 1.99497695, 0.99982566. In Arb: c_need(0.8052) = 1.892698192521…, c_need(0.99) = 1.994976941718…, c_need(0.0966) = 0.999825651279…, so the 8-decimal ceilings are right;
  - the new sources sentence.
- **Wording issues:**
  - "uniformly in (ε, φ)" (md), "its remainder not yet bounded uniformly" (Status), "to be certified uniformly in ε" (PAPERS Missing) and "the fixed-θ regime, φ → ∞, has its own scale ε²" (HANDOFF §5) all miss the θ ~ √ε regime (D3dj-1);
  - "Ψ's maximum on the straight segment" (D3dj-3).
- **Round 110's PAPERS and D61 edits are in:** the Claim is qualified (r < 0.99), Missing says |ρ| ≥ 1 (14 of 36), and Johansson, Tucker and Hales are cited.

**(5) Round 110's code** [own]
- **`c_need_ball`** forms the decimal exactly as arb(int)/10⁴ (a ball) and takes `.upper()` and a ceiling. It is sound.
- **The guards in `d_lines` and `v_lines`** and **`report`'s merge** read correctly.
- **The dense-supremum test.** On its 61 boxes, 51 have the maximiser on the quadratic piece, and min(u − dense) = −1.2·10⁻¹⁶ (W/own/dense_test.py).

**(6) Tests and mutants** (W/mut/: log_ship.txt, log_own.txt, ship_*/res_a.jsonl, own/res_b.jsonl; W/logs/semantic.log, vc_recompute.log)
- **Shipped mutants.**
  - partVC: V01–V05 are all killed.
  - round110: Z01, Z02 and Z04–Z13 are killed. Z03 survives, recorded as sound.
  - I re-ran 15 of the kills with the targeted test. Each fails on an assertion, not a crash; for example Z01's is "0.3490333 <= 0.3321243 + 1e-12".
- **My 27 mutants** (W/mut/spec115.py):

  | class | mutants | result |
  |---|---|---|
  | `part_VC` and `corner_sector`, non-equivalent | C01 θ sign, C02 z conjugated, C03 floor → y = 0, C05 \|η\| = 2 − ε/2, C06 top 0.01, C07 bisection reversed, C09 VC_R | survive; each changes the raw and the md's VC line on a rerun (vc_recompute.log) |
  | `part_VC`, non-equivalent | C04 ψ at φ = 0 | survives; a rerun would crash (ZeroDivisionError) |
  | equivalent | C08, C15 | survive |
  | `corner_psi`, `corner_floor` | C10–C13 | killed |
  | `corner_phi_dom` | C14 (cubic over 3.1, φ_dom = 1.438968) | survives; the test allows 10⁻³ |
  | `c_need_ball` | R01 binary float, R02 mid | survive; equivalent at 8 decimals |
  | `c_need_ball` | R03 floor | killed |
  | guards | R04–R06 ("inside" / "below 1" compared with 2) | survive; the planted 2.5 cannot tell |
  | `report` | R07 (any for all) | survives; equivalent on the data |
  | `d_lines` | R08 (r0 not floored) | survives; 0.09665 prints 0.0966 either way |
  | `box_bound_ball` | R09 (no rounding up) | killed by the 1.78 run's exact pin |
  | `dom_bb` | R10 (ball accepted below 10⁻⁹; round 108's code) | survives; no planted positive ball bound |
  | `d_lines` | R11 (no min with the grid's radius) | survives; mathematically equivalent |

**(7) Literature** [snippet] unless marked. arxiv.org could not be fetched in this session, so arXiv items are from search snippets only.
- Coalescence of a saddle with a singularity:
  - Chester–Friedman–Ursell;
  - Bleistein 1966 (DLMF §2.3; arXiv:2204.09305);
  - uniform asymptotics as a stationary point approaches an endpoint (arXiv:1707.07927);
  - two saddles coalescing with two branch points (arXiv:1610.04472);
  - a numerical method for oscillatory integrals with coalescing saddle points (arXiv:1806.06549, title only).
- These are uniform expansions. None of them is a validated valley inequality, which is all the claim needs.
- Validated numerics with a small parameter, and blow-up: rigorous numerics treating ε ∈ (0, ε₀] by factoring ε out (arXiv:1909.06207); the blow-up survey for fast-slow systems (arXiv:1901.01402); quasi-homogeneous weights (arXiv:1710.01629); rigorous slow-manifold enclosures by interval analysis (Guckenheimer–Johnson–Meerkamp); saddle-type blow-up with computer-assisted proofs (arXiv:2103.12390). Krupa–Szmolyan 2001 [memory].
- I found no validated certificate of a degenerate saddle–branch-point valley.

## Verdicts per part

- **(i), (ii), (iii)(a), (iii)(c), (v):** unchanged; proved as before (P_j's polynomial form is still a sketch).
- **(iii)(b)**
  - **|η| ≤ 1.995:** proved, computer-assisted, both steps in ball arithmetic. DB and V bands 0–1 recompute identically; c_need's floors are rigorous.
  - **1.995 < |η| < 2:** conjecture P.
- **Part VC:** correct as a pointwise limit (Ψ_φ, the floor, φ_dom, the ratios and the sector are verified). "Not a proof" is accurate. Its description of what a proof needs is incomplete (D3dj-1).
- **(iv):** unchanged: r < 0.8052 and r < 0.9900 (σ*/σ < 9.26, 199.00), with the thresholds now rigorous in Arb.

## Status of the round-110 issues

| id | status |
|---|---|
| D3de-1 (tests) | **Resolved.** Z01, Z02 and Z04–Z13 are killed on assertions (checked); part_DB stores the rigorous bound. Residue: the planted guard value is not discriminating (R04–R06; D3dj-2). |
| D3de-2 (HANDOFF) | **Largely resolved.** The status row, §8 and §10's state paragraph are fixed. Still stale: §10 "test_gh_peano.py 35" (38), "Last commits", "Referee ids"; §7's VC lines (D3dj-5). |
| D3de-3 (§7's θ_dom) | **Resolved.** 0.08606, 0.01484, 0.00731, 0.00144 and φ_dom; my values agree. |
| D3de-4 (floors, fallback) | **Resolved.** `c_need_ball` is printed with its comparison; there is no grid fallback; r0 is floored. |
| D3de-5 (PAPERS) | **Resolved.** The Claim is qualified; Missing says \|ρ\| ≥ 1 (14 of 36). |
| D3de-6 (sources, log, report) | **Resolved.** Johansson, Tucker and Hales are cited; the log carries the `part_D` note; `report` keeps a ball-checked D bb (planted 4936 test). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dj-1 | minor | **VC's remaining gap is misdescribed.** The blow-up's remainder is relative −ε(2φ² + O(1)) near the saddle (coefficients −0.249 … −200.26 at φ = 0 … 10) and O(ε\|w\|) far out. It cannot be bounded uniformly in φ. On the true bound VC's floor is positive from θ ≈ √(ε/2): ε = 10⁻³, φ = 25 gives E = +1.54·10⁻⁹ against min_y −2.78·10⁻⁹; ε = 10⁻⁴, φ = 100 gives +9.95·10⁻¹¹. A distinct limit governs θ = √εψ: E/ε² → Ψ̃_ψ = (1 + 2ψ²)y² − 4ψxy − y⁴/2 + D₊² − P₊², verified to 6 digits. The chain Ψ_φ → Ψ_∞ → Ψ̃_ψ → fixed-θ needs four parameter charts [own] (W/logs/blowup.log, floor_large_phi.log, mid_regime.log). | `vc_lines`' last sentence, D61 Status and "What (b) still lacks", PAPERS Missing, HANDOFF §§5, 7: "Ψ_φ governs θ ≪ √ε only (relative remainder ≈ −2εφ²); the route needs the regimes ε ≪ θ ≪ √ε (Ψ_∞), θ ~ √ε (Ψ̃_ψ) and θ ≫ √ε (fixed θ), with their matchings, and inner/matching/outer charts in z". Part VC: print the ratio at φ = 10, 20, 30 for ε = 10⁻³ and the coefficient (E/ε³ − Ψ)/(εΨ). Also state that φ″ → 0 is harmless on route (2) (the chord is the real segment) and that no overlap with the straight segment is needed. A certificate: one ball computation per chart, its small parameters as intervals including 0, with cancellation-free forms (my prototype, W/own/chart1_ball.py). |
| D3dj-2 | minor | **Part VC's computation is untested.** The test reads raw['VC'] but never calls `part_VC()` or `corner_sector()`. C01, C02, C03, C05, C06, C07 and C09 survive, and each changes the raw and the md's VC line on a rerun; C04 would crash a rerun. C14 moves φ_dom to 1.438968 and passes the 10⁻³ tolerance. R04–R06 survive because the planted `c_need_ball` = 2.5 exceeds every candidate threshold [own]. | tests/test_gh_peano.py: `assert GP.part_VC() == raw['VC']` (0.02 s); `abs(pd − root of 8φ³ − 12φ² + 1) < 1e-9`; plant 1.9999 for `v_lines`, 1.8928 for `d_lines` and 1.5 for the "below 1" check. tools/mutants/partVC_D61_own.py: add C01–C09 and C14. |
| D3dj-3 | cosmetic | "On the straight segment (y = 0) Ψ's maximum is −0.2500 at φ = 0" (md; D61 Numerics "Ψ's maximum −0.2500"): this is the maximum over x ≥ ½. On the whole segment the maximum is 0, at the saddle. `corner_phi_dom`'s docstring says "(and −x² below ½), negative iff φ < φ_dom" [md]. | "beyond the disc (x ≥ ½)" in `vc_lines`, D61 and the docstring. |
| D3dj-4 | cosmetic | Floating point cannot resolve the ε³ margins below ε ≈ 10⁻⁴. `corner_sector`: lo/ε = 1.44023 at 1.9999 (correction 5.40 against 4.26) and 1.24847 at 1.99999. `cert_E`/ε³ is off by up to 1.19·10⁻⁴ (relative) at ε = 10⁻⁴ [own]. No printed number is affected (VC_R ≤ 1.999). | `corner_sector`'s docstring and VC: valid for ε ≥ 10⁻³ in floating point; any extension in Arb (the cancellation-free form). |
| D3dj-5 | cosmetic | HANDOFF is stale. §7: "(2026-11-01, not yet in any md or draft)", "Put the blow-up … as a part (VC)", "A probe (not in any md yet)", "Re-reviews due: … D‴.61 (after round 90)". §10: "test_gh_peano.py 35" (38), "Last commits: b98e652 …", "Referee ids: … next D3cz = 105". §5: "cover the overlap with the straight segment (φ < φ_dom)" is not needed [md]. | HANDOFF §§5, 7, 10. |
| D3dj-6 | cosmetic | `vc_lines` prints VC_X and groups the rows by VC_EPS from the module's constants, not the raw. A zero Ψ at a floor point divides by zero. The dense-supremum test's tolerance is 10⁻¹², against measured noise of 1.2·10⁻¹⁶. Round 108's `u < 0` acceptance is untested (R10 survives) [own]. | `vc_lines`: take x and ε from the raw's rows and skip Ψ = 0. Test: tolerance 10⁻¹⁴, and a planted box with a ball bound in (0, 10⁻⁹) rejected by `dom_bb`. |

## Sources

- [Method of Chester–Friedman–Ursell](https://en.wikipedia.org/wiki/Method_of_Chester%E2%80%93Friedman%E2%80%93Ursell) [snippet]
- [DLMF §2.3 (Bleistein 1966; Bleistein–Handelsman)](https://dlmf.nist.gov/2.3) [snippet]
- [Recent Advances in Asymptotic Analysis (arXiv:2204.09305)](https://arxiv.org/pdf/2204.09305) [snippet]
- [Uniform asymptotics as a stationary point approaches an endpoint (arXiv:1707.07927)](https://arxiv.org/pdf/1707.07927) [snippet]
- [Uniform asymptotic expansions: two saddles coalescing with two branch points (arXiv:1610.04472)](https://arxiv.org/abs/1610.04472v2) [snippet]
- [A numerical method for oscillatory integrals with coalescing saddle points (arXiv:1806.06549)](https://www.arxiv.org/pdf/1806.06549) [snippet: title only]
- [Rigorous numerics for a singular perturbation problem (arXiv:1909.06207)](https://arxiv.org/pdf/1909.06207) [snippet]
- [A survey on the blow-up method for fast-slow systems (arXiv:1901.01402)](https://arxiv.org/pdf/1901.01402) [snippet]
- [Stabilization of slow-fast control systems: the non-hyperbolic case (arXiv:1710.01629)](https://arxiv.org/pdf/1710.01629) [snippet]
- [Guckenheimer–Johnson–Meerkamp, rigorous enclosures of slow manifolds](https://pi.math.cornell.edu/~pmeerkamp/slowmfdenclosures.pdf) [snippet]
- [Saddle-type blow-up solutions with computer-assisted proofs (arXiv:2103.12390)](https://arxiv.org/pdf/2103.12390) [snippet]
- Krupa–Szmolyan, SIAM J. Math. Anal. 2001 [memory]
- Johansson, Arb, IEEE Trans. Comput. 66 (2017), arXiv:1611.02831 [memory; as cited in D61]

## Lesson

1. **Compute the next-order coefficient as a function of the blow-up parameters before calling a remainder "to be bounded uniformly".** Here (E/ε³ − Ψ_φ)/(εΨ_φ) ≈ −2φ² gives away a second distinguished limit at θ ~ √ε. A sample stopping at φ = 5 cannot show it.
2. **Evaluate the proposed path on the true function past the sampled range.** φ = 25 at ε = 10⁻³ turns VC's floor positive, while the true valley stays negative.
3. **A test that reads a part's raw without recomputing the part pins only the report.** When the part is cheap (part VC takes 0.02 s), recompute it in the test.
4. **A planted value in a guard test must lie between the true threshold and the plausible wrong ones.** Otherwise a wrong comparison survives.

My files are all under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref115:
- own/: rerun_vc.py, blowup.py, floor_large_phi.py, vc_more.py, mid_regime.py, chart1_ball.py, chart1_dbg.py, vc_recompute.py, dense_test.py, sample_dbv.py, semantic.py;
- logs/: pytest.txt, rerun_vc.log, blowup.log, floor_large_phi.log, vc_more.log, mid_regime.log, chart1_ball.log, chart1_ball_check.log, vc_recompute.log, sample_dbv.log, semantic.log;
- mut/: spec115.py, log_ship.txt, log_own.txt, ship_partVC_D61_own/res_a.jsonl, ship_round110_D61_own/res_a.jsonl, own/res_b.jsonl;
- tools_ref/mutants_bt.py;
- out/: the regenerated md and the raw copies.
