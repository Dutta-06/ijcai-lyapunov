# Referee report D3x (round 25): Proposition D‴.35, the split of a class as an S_n-symmetric transcritical bifurcation

**Scope.** PARADIGM §5d: Proposition D‴.35 (statement (i)–(v), proof sketch, Status, *Numerics* (A)–(E)). Also the §0 and §7 rows "| D‴.35 |", the
§8 bullet "Done (Proposition D‴.35)", `results/RESULTS_PARADIGM.md` (Proposition D‴.35), LITERATURE_TRANSPORT [S63] and the audit bullet
"Proposition D‴.35", PAPERS.md (paper 13) and README.md where they mention D‴.35, `experiments/exp_split_bifurcation.py` →
`results/chaos/split_bifurcation.md` (+ raw), and `tests/test_split_bifurcation.py`. Context read: Proposition D‴.32 (iii)–(iv) and Numerics (G),
Proposition D‴.34, `src/core/cohesion_leak.py`, `src/core/mass_transfer.py`, and REVIEW_LOG rounds 23–24.

Branch `claude/festive-bardeen-eopej2` at 44fe0f9 (`git log --oneline -5`: 44fe0f9, dc2787b, 9d122dd, 364fab8, 1ddff8c). The tree was
extracted with `git archive 44fe0f9` into `ref_D3x/tree`, and everything below ran there. Uncommitted edits by the concurrent session in the
shared tree were ignored. No repository file was modified, and no git command touching the index or the tree was run. Python ran with
`PYTHONDONTWRITEBYTECODE=1`, `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`, with at most two processes at a time. pytest ran with
`-p no:cacheprovider`.

**Code used.** My instance is `scripts/inst.py`, the independent implementation written by the round-23 referee (D3u) and imported by no
repository module. It has the saturated network, analytic ∇Φ̃ and ∇²Φ̃, the averaged learning drift, both kernels, analytic Jacobians and
single linkage. My own additions are in `scripts/common.py`:
- ϑ^b(q) by Newton continuation in ε, which matches the repository's BFGS approximant, as saved by the round-23 referee, to 5.7·10⁻¹¹;
- the gauge transform;
- the mass-free force;
- c₃, both from the analytic Hessian differenced and from a 6-point stencil of ⟨rΦ̃⟩;
- a midpoint integrator.

Every confirmed number below comes from these scripts, not from the code under review. Only `tests_check.py` imports the repository, to run
its tests and in-memory mutants of its code.

Line numbers are those at 44fe0f9:
- PARADIGM: the statement l. 4380–4387, (i) l. 4389–4395, (ii) l. 4397–4401, (iii) l. 4403–4408, (iv) l. 4410–4412, (v) l. 4414–4418, the
  proof sketch l. 4420–4425, the Status l. 4427–4436, the *Numerics* l. 4438–4460 ((A) l. 4439–4442, (B) l. 4442–4446, (C) l. 4446–4452,
  (D) l. 4452–4458, (E) l. 4458–4460).
- The §0 row l. 145, the §7 row l. 6901, the §8 bullet l. 7165–7172.
- RESULTS_PARADIGM l. 917–930. LITERATURE_TRANSPORT [S63] l. 83 and the audit bullet l. 529–537. PAPERS.md l. 226. README l. 75–77.
- `exp_split_bifurcation.py`: `coefficients` l. 63–72, `b_job` l. 97–133, `reduced_crossing` l. 136–150, `c_job` l. 153–185 (the 0.1
  criterion l. 169–174), `MassFree` l. 188–194, `d_job` l. 197–219.
- `tests/test_split_bifurcation.py` l. 29–84 (test 3 l. 64–84, its σ at l. 68).

**Scripts.** Scratch scripts and their `.out` files are in
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3x/scripts/`:

| script | what it does |
|---|---|
| `inst.py`, `common.py` | the independent instance (see above), the approximant, the gauge transform, the mass-free force, the midpoint rule |
| `coef.py` → `coef.out` | (A): ϑ^b(q), λ⁺, v and c₃ for each class. (i): the even part of each drift term along random zero-sum ξ (learning; compact and Gaussian cohesion; mass-free leak) against (c₃/2)(ξ² − S₂/n). The mean term. The full Jacobian's spectrum at the threshold, at μ = 5·10⁻⁵, at σ₀/16, and for the mass-free leak at σ₀/32 |
| `branches.py` → `branches.out` | (ii), (B): the split states q = 1–4 by Newton on my full drift with my analytic Jacobian, at μ = +5·10⁻⁷ and −5·10⁻⁷: positions, eigenvalues on the separations along v, unstable count |
| `q5.py`, `q5b.py`, `cubic.py` → `.out` | the reduced equation on Fix(S_q × S_{n−q}) with every other coordinate slaved: the effective cubic coefficient, split into the cohesion's direct term, the learning's direct term and the slaving; the S₅ × S₅ states by Newton and their spectra |
| `partC.py` → `C_<μ>_<seed>.out` | (C): all 12 runs reproduced from their seeds with my drift (midpoint rule, dt = 0.25); the departures at 0.1 from the median; exits from the kernel's range; the reduced equation (DOP853), the linear equation alone and the flipped-c₃ equation; the end state's margin decomposed |
| `preds12.py` → `preds12.out` | the predictions of the three reduced equations for all 12 seeds; the first particle on the side sign(c₃) under each |
| `endgeom.py` → `endgeom.out` | the end states of the 12 runs: the departed particles' whitened distances to the main piece and to each other |
| `massfree_static.py` → `.out` | (v): the mass-free leak's linearisation at the lift; the sum of the forces; the asymmetry of G⁻¹J, against the cohesion force's |
| `partD.py` → `D_mf_0_{a,b,c}.out` | (D): the mass-free leak, seed 0, Part G's schedule reproduced from `default_rng(0)`, then extended to σ₀/128 and σ₀/256. Each stage prints the classes, their margins, the fit error and the largest drift |
| `tests_check.py` → `tests_check.out` | the tests' regime; test 3's jitter under the three reduced equations; seven mutants of the code under review patched in memory |

Also in `ref_D3x/`: `tests.out` (3 passed, 37 s), `from_raw.out`, `md5_before.out` and `md5_after.out`. `--from-raw` rewrites
`split_bifurcation.md` byte-identically: md5 `370eaf78c899d90c6aa5c6c57e072698`, the same as `git show 44fe0f9:results/chaos/split_bifurcation.md`.

## Verdicts

| object | verdict |
|---|---|
| (i) the reduced equation | **Holds to second order.** The even part of the full learning drift along v is (c₃/2)(ξ_i² − S₂/n) to 3.2·10⁻⁸ against the term's scale 5.2·10⁻⁴; the error falls 4.0× when ε is halved. The mean term is common (the constant −2.13·10⁻⁴), so it drops with the projection. The compact cohesion, the Gaussian cohesion and the mass-free leak have no even part (≤ 1.5·10⁻¹⁴). λ⁺ is simple, and the other blocks have nonzero rates at the threshold. Caveat (ID-7): the "near μ = 0" regime does not cover Part C at σ₀/16, where three stable modes are slower than μ. At third order the class mean's slaving is the largest term. |
| (ii) the branches | **Holds for 1 ≤ q < n/2, for both signs of μ.** Every one is a saddle, and the eigenvalues converge as stated; all of (B)'s numbers reproduce. The list of equilibria near 0 is incomplete for even n (ID-6). The S₅ × S₅ states exist at μ < 0, at ξ = ±1.21·10⁻² (μ = −5·10⁻⁷); they are saddles too, so "no small split is stable" survives. |
| (iii) order, time, "alone" | **Order and time hold; "alone" fails in 3 of 6 runs at σ₀/16 (ID-3).** The order holds on the centre manifold, and it follows from S_n-equivariance alone. The reduced quadratic equation gives the first departure's particle in 12 of 12 runs; the linear equation alone gets 6 of 12, the flipped c₃ 2 of 12. Its time is within 0.01%–1.15% of the measured one, against 0.5%–22.4% for the linear equation. In every run where a piece left the main piece's range, the top-ranked particle was in it. At σ₀/16, three of six runs shed a **bound pair**, the two particles 0.01–0.27 apart and both 1.06–1.19 from the main piece. The 12-of-12 "first on the side sign(c₃)v" statistic is a tautology (ID-5). The "+ side" explanation is wrong (ID-4). |
| (iv) the cascade | **The margin's accounting is wrong (ID-2).** The remaining particles' amplitude rises by n/(n − 1) (the D‴.28 gauge), so the residual curvature rises by λ⁺/(n − 1). That is 2.2× the stiffness drop σ/m at σ₀/16 and 1.35× at μ = 5·10⁻⁵. Stopping the cascade needs a drop of more than μ + σ₀/m, not σ₀/m. In all 12 runs both terms are swamped by the residual's change (+1.8·10⁻⁴ to +3.4·10⁻⁴), so the margin turns positive within the stage. |
| (v) the mass-free leak | **The blocks hold. "Sheds once per crossing" and "the mass factor drives the cascade" are refuted (ID-1).** The blocks reproduce to the digit: separations at −σ, class means 0. The forces do not sum to zero, and G⁻¹J is asymmetric (0.22 against 3·10⁻⁹ for the cohesion). Extended to σ₀/128 and σ₀/256, seed 0 sheds several particles per stage. It ends at 9 classes with the fit error 0.0028, Part G's end. The 4 classes at σ₀/64 reflect where the schedule stops. |
| Numerics (A) | All confirmed (my code): λ⁺, v, c₃ to the printed digits, and the thresholds 0.0424, 0.1213, 0.0460. |
| Numerics (B) | All confirmed: positions 0.31%, 0.44%, 0.83%, 3.0% from the prediction; eigenvalues +1.25108μ (×8), −1.02305μ, and so on; 8, 7, 6, 5 unstable directions. |
| Numerics (C) | All 12 runs reproduced to the md's digits: departures, times, clusters, main pieces and margins. The readings overstate (ID-3, ID-4, ID-5). |
| Numerics (D) | The mass-free run of seed 0 reproduced at every stage (4 classes, centres −1.50, −0.44, +0.07, +1.35, fit 0.0149). The reading is refuted (ID-1). The replicator runs were not reproduced. |
| Numerics (E) | The reading is not supported against γ = 0.1 (ID-9). |
| tests | 3 pass. Test 3 runs at a **negative** cohesion strength (ID-8). Caught: the flipped-c₃, doubled-c₃ and t*-without-the-2 mutants, and a quadratic cohesion term at ≥ 8% of \|c₃\|/2. Not caught: a quadratic cohesion term at 2%. `reduced_crossing` is untested. |
| §0, §7, §8, RESULTS, PAPERS, README | They carry ID-1 to ID-3 and ID-5 ("alone", "the mass-dependent stiffness makes the cascade", "the target's count", "12 of 12"). |
| novelty, [S63] | Classical as the audit says. Part of the "new" list is generic or elementary, and precedents are missing (ID-10). |

## Confirmed issues

### ID-1 (major): the mass-free leak does not shed once per crossing; its 4 classes reflect where the schedule stops

**Where.** (v) l. 4414–4418; *Numerics* (D) l. 4455–4458; the Status l. 4430–4431 and 4435–4436; the §0 row l. 145 ("a mass-free leak … sheds
once per crossing — the target's count, not its fit", "the cascade through the mass factor"); the §7 row l. 6901; the §8 bullet l. 7169–7171;
RESULTS l. 926–927; PAPERS l. 226 and README l. 76–77 ("the class's mass-dependent stiffness makes the cascade"); the audit l. 536.

**Claim.** "(v) *(the mass factor drives the cascade)* … a departure leaves the class's stiffness at σ …: a class sheds once per crossing of its
threshold". (D): "at 4 classes in 3 of 3 runs, the target's count … one particle leaving per split".

**Evidence** (`partD.py`, seed 0; the first seven stages reproduce the md exactly).

1. **The schedule shifts by the factor w_k, and that confounds the comparison.** The mass-free leak's separation stiffness is σ, the cohesion
   leak's σw_k = σ/3. Its thresholds are therefore λ⁺_k/σ₀ = 0.0141, 0.0404, 0.0153 against 0.0424, 0.1213, 0.0460 (`coef.out`). The halving
   schedule takes the cohesion leak three stages past the middle class's threshold (σ₀/16, /32, /64) and past the outer classes' at σ₀/32 and
   /64. It takes the mass-free leak only two stages past the middle class's (σ₀/32, /64) and none past the outer ones.
2. **At σ₀/64 the outer class is just crossing.** The end-of-stage margins at σ₀/64 are +1.38·10⁻⁵ (the middle class's 9-particle piece) and
   −5.9·10⁻⁶ (class 3); the largest drift is still 1.4·10⁻⁵.
3. **Extending the schedule.** At σ₀/128 the mass-free leak goes to **9 classes** (fit 0.0081). The middle piece sheds three particles in one
   stage (7 + 1 + 1 + 1), and class 3 splits into 5 + 3 + 1 + 1. At σ₀/256 it stays at **9 classes** with the fit error **0.00278**: singletons
   at −0.28 to −0.23 and +1.68 to +1.77, the middle piece at +0.21, class 3's piece at +1.14. This is Part G's end under the cohesion leak (8–10
   classes, 0.0028–0.0031, singletons at −0.29 to −0.23 and +1.69 to +1.76).

**Conclusion.** Without the mass factor the cascade has the same extent, one to two stages later, as the factor-3 shift predicts. "Sheds once
per crossing", "one particle leaving per split" and "the target's count" describe the point where the schedule was stopped, not the law.

**Fix.**
- Restate (v): the mass-free leak has the blocks −(σI − R^G_k) on the separations, so its thresholds are σ = λ⁺_k, lower than the cohesion
  leak's by the factor w_k. A departure leaves the stiffness at σ, but the residual curvature still changes (ID-2), so nothing limits the
  number of departures per crossing. Remove "(the mass factor drives the cascade)" and "sheds once per crossing".
- Part D: extend the mass-free schedule to σ₀/256, or run it at the matched stiffness σ = w_kσ_coh, and report the counts there.
- Replace "4 in 3 of 3 … the target's count" in the §0 row, §8, RESULTS, PAPERS and README with "4 classes at σ₀/64 under the shifted
  thresholds, 9 at σ₀/256 (seed 0)".

### ID-2 (major): (iv)'s margin after a departure omits the gauge rise of the residual curvature, and its condition is off by μ

**Where.** (iv) l. 4410–4412; the proof sketch l. 4424 ("Proposition D‴.32(iv) at the new configuration"); the §0 row ("after a departure the
class's stiffness σ₀w_k has dropped by σ₀/m, which drives the cascade"); §8 l. 7169–7170; PAPERS l. 226; README l. 76–77.

**Claim.** "the leak's stiffness falls with the class's particle fraction, so the class sheds again unless the departed particle lowers the
residual curvature by more than σ₀/m".

**Evidence.**

1. **The arithmetic.** Before the departure the margin is σw − λ⁺ = −μ < 0. After it, the margin is σ(w − 1/m) − λ⁺′, which is positive only if
   λ⁺ − λ⁺′ > μ + σ/m. A drop of more than σ/m leaves the margin negative whenever μ > 0.
2. **The gauge.** The class's weighted amplitude w_ka_k is carried by n − 1 particles, so each particle's amplitude, and with it R^G_k, scales by
   n/(n − 1): this is Proposition D‴.28's q_k/w_k law. The curvature therefore **rises** by λ⁺/(n − 1) ≈ σw/(n − 1) ≥ σ/m, as large as the
   stiffness drop or larger.
3. **Measured on all 12 Part C runs** (`partC.py`; λ⁺ at the main piece's mean with the end residual). The margin change splits into three
   terms:

   | run set | stiffness drop | gauge rise (computed: residual unchanged, amplitude ×10/n_main) | residual and atom |
   |---|---|---|---|
   | σ₀/16, one piece of 1 | −1.47·10⁻⁵ | −3.18·10⁻⁵ | +1.93·10⁻⁴ |
   | σ₀/16, main piece of 8 | −2.95·10⁻⁵ | −7.15·10⁻⁵ | +2.9·10⁻⁴ to +3.4·10⁻⁴ |
   | μ = 5·10⁻⁵ | −2.36·10⁻⁵ | −3.18·10⁻⁵ | +1.8·10⁻⁴ to +2.7·10⁻⁴ |

   The gauge rise is 1.35× (μ = 5·10⁻⁵) to 2.4× (σ₀/16) the stiffness drop. Both are an order of magnitude below the residual's change, which
   is why the margin is positive at the end of every run (7.5·10⁻⁶ to 1.6·10⁻⁴, matching the md).
4. **What actually drives the cascade.** Within a stage the class does not cascade. What restarts it in Part G is the halving of σ, together
   with the residual's evolution as the classes move toward the fourth bump, which ID-1 shows proceeds without the mass factor.

**Fix.**
- Restate (iv): the margin after k departures is σ(w_k − k/m) − λ⁺(r′, ϑ′), with λ⁺ evaluated with the remaining particles' raised amplitude
  (the gauge factor n/(n − k) when the residual is unchanged). The class sheds again unless λ⁺ falls by more than μ + kσ/m.
- Report the three-term decomposition (C), and state that on the instance the residual's change dominates and stops the cascade within a stage.
- Remove "which drives the cascade" and "the mass-dependent stiffness makes the cascade" from §0, §8, PAPERS and README, or replace it with
  "the stiffness drop and the gauge rise of the curvature both favour further departures; the residual's change opposes them".

### ID-3 (major): "alone" and "sheds single particles" fail in half of the runs at Part G's first split; "leave" means 0.1 from the median, well inside the range

**Where.** (iii) l. 4406–4407 ("alone unless the next particle's ζ₀ is within a factor e^{μΔs} of its own", "A split sheds single particles");
(C) l. 4450–4451 ("One to four particles leave within 10⁵ slow units"); §8 l. 7165 and 7168 ("why a split sheds single particles", "the extreme
particle of the jitter leaves first, alone"); PAPERS l. 226 and README l. 76 ("leaves first and alone"); the code's criterion at
`exp_split_bifurcation.py` l. 169–174.

**Evidence** (`partC.py` and `endgeom.py`; all 12 runs reproduce the md).

1. **The criterion.** A particle counts as "leaving" when it is 0.1 from the median of ξ. That is far inside the kernel's range R = 1, and
   inside the distance R/√7 = 0.38 where the compact pair force peaks.
2. **Who left the main piece's range by s = 10⁵**, with the departed particles' whitened distances to the main piece:
   - σ₀/16:

     | seed | what left the range | distances |
     |---|---|---|
     | 0 | particle 3 alone | 1.94 |
     | 1 | the **pair** 8 + 2 (0.27 apart) | 1.06, 1.12 |
     | 2 | the **pair** 3 + 5 (0.01 apart) | 1.19, 1.19 |
     | 3 | the **pair** 0 + 4 (0.06 apart) | 1.10, 1.13 |
     | 4 | particle 1 alone | 1.52 |
     | 5 | particle 8 alone | 1.55 |

   - μ = 5·10⁻⁵: one particle alone in seeds 0, 1 and 5 (1.23–1.24 from the main piece). In seeds 2–4 nothing has left the range yet: the
     departed particles are 0.28–0.76 from the main piece.
3. **The md mixes criteria.** For σ₀/16 seed 2 it lists four "departures" (2, 5, 3, 8). Particles 2 and 8 never left the range; 5 and 3 left
   as a doublet. The top-ranked particle is in every departing piece.
4. **When pairs form.** The pairs form where the top two ζ₀ are within 1.2%–11% of each other (seed 2: 7.82 against 7.73·10⁻⁴; seed 3: 8.21
   against 7.85; seed 1: 10.8 against 9.72). The single departures had ratios of 1.29–4.7. "Within a factor e^{μΔs}" leaves Δs undefined, so
   the stated exception has no content.

**Fix.**
- (iii): "the top-ranked particle leaves first along sign(c₃)v; when the next particle's ζ₀ is within a few per cent of it (on the instance
  ≤ 11%), the two leave together as a bound pair". Define the window by the time the first one takes to cross the range.
- (C): count departures from the kernel's range and report the pairs (3 of 6 at σ₀/16).
- §8, PAPERS, README: replace "alone" and "sheds single particles" by "sheds the extreme particle, alone or with a near-tied neighbour".
- Part G's "single particles" may include such pairs split by the 0.25 linkage (σ₀/16 seed 1's pair is 0.27 apart and counts as two
  singletons). Report the pieces at the range.

### ID-4 (minor): the "+ side" first departures are explained wrongly

**Where.** (C) l. 4449–4450: "At σ₀/16 two runs flag first a particle on the other [side] — the remaining particles' common drift −(c₃/2n)Σξ²
measured against the median, which the reduced equation reproduces." The word "side" is missing in the text.

**Evidence.**
1. A drift common to all particles shifts the median by the same amount, so it cancels exactly in ξ_i − median and cannot flag anything.
2. The **linear** reduced equation (c₃ = 0, no common term at all) predicts the same + side particles first: σ₀/16 seed 2, particle 2 at
   2.85·10⁴; seed 4, particle 8 at 2.99·10⁴ (`preds12.out`). The flipped-c₃ equation does too.
3. The cause is the departure level. At σ₀/16 the level 0.1 lies below 2μ/|c₃| = 0.221, the excursion at which the quadratic term equals the
   linear one. A particle with the largest |ξ_i(0) − median| on the other side therefore reaches 0.1 first by linear growth.
4. Those particles never leave the class: at s = 10⁵, particles 8 and 0 are still in seed 4's main piece (single linkage at 0.25), and particle 2 is still in seed 2's.

**Fix.** "At σ₀/16 the level 0.1 lies below 2μ/|c₃| = 0.22, so the first excursion to 0.1 is decided by linear growth: in two runs it is a
particle on the other side, which the linear equation predicts as well, and which does not leave the class; the particle that leaves the range
is the top-ranked one on the side sign(c₃)v." Report the departure at the kernel's range, or at a level above 2μ/|c₃|.

### ID-5 (minor): "first on the side sign(c₃)v is the top-ranked one in 12 of 12" is a tautology; §8's "predicted in 12 of 12" is false as worded

**Where.** (C) l. 4447; the §0 row ("the first particle to leave on the side sign(c₃)v the predicted one in 12 of 12 runs"); RESULTS l. 922;
§8 l. 7168–7169 ("the extreme particle of the jitter leaves first, alone, on the side the quadratic term picks — predicted in 12 of 12 runs").

**Evidence.**
1. The particle with the largest sign(c₃)ξ_i is by definition the farthest from the median on that side (12 of 12 in `preds12.out`). Any order-
   preserving law makes it the first on that side: the linear equation alone gives the same particle in 12 of 12, so the statistic says nothing
   about c₃.
2. The discriminating facts are the first departure overall and its time. The quadratic equation gets the particle in 12 of 12, the linear one
   in 6 of 12, the flipped c₃ in 2 of 12. The quadratic times are within 0.01%–1.15%; the linear ones are off by −10.5% to +22.4%.
3. The first departure overall is on the side sign(c₃)v in 10 of 12 runs, not 12.

**Fix.**
- Report the discriminating statistics: the particle in 12/12 against 6/12 (linear) and 2/12 (flipped), and the times within 1.2% against up
  to 22%.
- Drop the tautological count, or label it "order preservation".
- §8: "on the side the quadratic term picks in 10 of 12 runs (the other two: linear growth to the 0.1 level, ID-4)".

### ID-6 (minor): (ii)'s list of the equilibria near 0 is incomplete for even n; "leaves every neighbourhood" is too strong

**Where.** (ii) l. 4397–4401; the proof sketch l. 4421–4423.

**Evidence.**
1. On Fix(S_{n/2} × S_{n/2}) the quadratic term vanishes identically, so the q = n/2 states are cubic (pitchfork) states at |ξ| ~ √|μ|. The
   instance has n = 10.
2. With everything else slaved (`q5b.py`, `cubic.py`), the reduced equation is dξ/ds = μξ + κξ³ with κ = +3.43·10⁻³; the value is stable to
   10⁻⁴ over t = 2.5·10⁻³ to 2·10⁻². So the states exist for μ < 0.
3. Newton on the full drift finds them at ξ = ±1.2067·10⁻² (μ = −5·10⁻⁷) and ±2.409·10⁻² (μ = −2·10⁻⁶), with the drift 5.5·10⁻¹⁷ and
   3.9·10⁻¹⁷. Each is a saddle with 5 unstable directions: 4 inside a group (+31.5|μ|) and the radial one (+1.98|μ|).
4. At μ > 0 there is no such state.
5. (ii)'s other claims hold for both signs (`branches.out`). At μ = −5·10⁻⁷ the q = 1–4 states are saddles with 1, 2, 3, 4 unstable directions,
   for example q = 1: −1.249|μ| (×8) and +0.977|μ|.
6. "The class leaves every neighbourhood of the replicated state" holds generically, outside the saddles' stable manifolds. For a fixed
   neighbourhood it holds as μ → 0⁺, since the O(μ) states lie in every neighbourhood.

**Fix.**
- (ii): "for 1 ≤ q < n/2 …; for even n the S_{n/2} × S_{n/2} states are cubic, at |ξ| ~ √(|μ|/|κ|), and on the instance (κ > 0, μ < 0) they are
  saddles too".
- "leaves every fixed neighbourhood for almost every initial condition".

### ID-7 (minor): the centre-manifold regime does not cover Part C at σ₀/16; the third order is mostly the class mean's slaving, not the kernel's concavity

**Where.** (i) l. 4389–4395; (iii) l. 4407–4408 ("the kernel's concavity … acts at third order and in the same sense"); (C) l. 4446–4447
(σ₀/16 presented as a check of (iii)); the Status l. 4427–4429.

**Evidence.**
1. **No spectral gap at σ₀/16** (`coef.out`, full Jacobian at the lift). The centre eigenvalue is +1.387·10⁻⁴ (×9). The next ones are
   −3.89·10⁻⁵ (class 3's separations), −4.74·10⁻⁵ (class 1's) and −1.474·10⁻⁴ (class 2's separations along its amplitude direction): three
   stable modes slower than or as slow as μ. At μ = 5·10⁻⁵ the gap is 1.28·10⁻⁴ = 2.6μ.
2. **Why Part C still agrees.** The stable modes couple to ξ only at third order or through common terms, and the jitter is small. That, not
   the centre manifold, explains the agreement at σ₀/16.
3. **The third order** (`cubic.py`, q = 1 at the threshold, in the group parameter t). The effective cubic coefficient is +0.155, made of:
   - the cohesion's concavity, +0.086;
   - the residual's direct cubic, −0.588;
   - the slaving, +0.657, mainly the class mean, which moves by 24t² ≈ 0.30ξ_A² along v.

   "In the same sense" holds for the net cubic on this instance. But the kernel is one of three cubic terms, and the two larger ones nearly
   cancel.

**Fix.**
- (i)/(Status): "valid for |μ| small against the other modes' rates (on the instance |μ| ≪ 1.8·10⁻⁴); at σ₀/16 three stable modes are slower
  than μ, and the reduced equation's agreement there reflects the weak coupling".
- (iii): "at third order the class mean's slaving and the residual's own cubic term (opposite, nearly cancelling) join the kernel's concavity".

### ID-8 (minor): test 3 runs at a negative cohesion strength; `reduced_crossing` is untested

**Where.** `tests/test_split_bifurcation.py` l. 68 (`sigma = (lam - 4e-4) / w[1]`); RESULTS l. 929–930; the §7 row.

**Evidence** (`tests_check.out`).

1. **The test's regime.** At the counts (4, 4, 4), λ⁺ = 2.8607·10⁻⁴ < 4·10⁻⁴, so σ = −3.418·10⁻⁴ = −0.048σ₀: a repulsion. Every class then has
   two unstable separation directions. For class 2, σw − eig R^G = (+1.87, −1.14, −4.0)·10⁻⁴; for class 1, (+0.13, −1.14, −2.14)·10⁻⁴. The
   test therefore checks a law outside the proposition.
2. **What it discriminates.** Its jitter does discriminate the quadratic term: the linear equation predicts particle 0 on the + side, and the
   quadratic one particle 2 on the − side, which is what happens. But the two times differ by only 0.7% (1.345 against 1.354·10⁴).
3. **Mutants** (patched in memory). Caught:
   - c₃'s sign flipped (all three tests fail);
   - c₃ doubled, a missing ½ (tests 1 and 2 fail; test 3 passes);
   - t* without the factor 2 (test 2 fails);
   - an even, quadratic cohesion term at 8% or 40% of |c₃|/2 (test 1 fails).

   Not caught: the same term at 2%. That is acceptable.
4. **`reduced_crossing`.** It produces all of Part C's predictions and is covered by no test.

**Fix.**
- Test 3 at μ = 1·10⁻⁴ (σ > 0), with an assertion that σ > 0 and that the other blocks are stable.
- A test of `reduced_crossing` against the closed-form Riccati time in the one-particle limit.
- A test of `MassFree`'s blocks (−σ on the separations, 0 on the class means).

### ID-9 (minor): Part E's reading is not supported by its three seeds

**Where.** (E) l. 4458–4460 ("a faster mass transfer shortens the cascade"); RESULTS l. 927–928; §8 "the replicator's rate against the
cascade's".

**Evidence** (from the md).
1. γ = 1 ends at 5, 6 and 8 classes, γ = 0.1 at 6, 6 and 7: the same mean, 6.33, with overlapping ranges.
2. Only against the cohesion leak alone (8, 9, 10) are both shorter.
3. Part G's runs draw the jitter from `default_rng(700 + seed)` (`exp_cohesion_leak.g_job`), Parts D and E from `default_rng(seed)`, so the
   "seeds 0–2" comparison is not paired.
4. γ = 1, seed 0 goes from 6 classes at σ₀/32 to 5 at σ₀/64: a class merges or is absorbed, which the reading does not mention.

**Fix.** "With the replicator (γ = 0.1 or 1) the cascade ends at 5–8 classes against 8–10 without; the rate γ makes no difference in three
seeds". Pair the seeds, or say they are not paired.

### ID-10 (minor): novelty: part of the "new" list is generic or elementary, and precedents are missing

**Where.** The Status l. 4431–4436; the §0 row; LITERATURE_TRANSPORT [S63] l. 83 and the audit l. 529–537.

**Evidence.**
1. **Order preservation is generic.** It holds for every S_n-equivariant flow on R^n: Fix(transposition) is invariant, so no two coordinates
   cross. It is not specific to the cohesion leak.
2. **The escape time is elementary.** It is the Riccati solution of dζ/ds = μζ + aζ², checked here: s = (1/μ)ln[ζ₁(μ + aζ₀)/(ζ₀(μ + aζ₁))].
3. **Missing from [S63] and the audit:**
   - Stewart, Elmhirst and Cohen, *Symmetry-breaking as an origin of species*, in *Bifurcation, Symmetry and Patterns* (Birkhäuser, 2003),
     pp. 3–54: speciation as an all-to-all coupled S_N-symmetric jump bifurcation. Details checked by a web search; the [S63] entry cites only
     the Golubitsky–Stewart book.
   - Cohen and Stewart's earlier polymorphism paper, from memory.
   - The Dias–Stewart cubic truncation, whose branch stabilities Elmhirst computes; seen in a citing abstract.
   - Brea et al. 2019 (arXiv:1907.02911, permutation saddles of duplicated units) and Şimşek et al. 2021 (ICML, symmetry-induced critical points
     of overparameterised networks): the replicated class is such a point. Both from memory.
   - Rose's phase transitions in deterministic annealing ([S46], [S60]): a cluster's split at a critical parameter set by its covariance.
     Already in the bibliography, but not in D‴.35's audit.
4. **What was not seen.** Ihrig–Golubitsky and Chossat–Lauterbach–Melbourne were not seen, as the [S63] entry says. arxiv.org was unreachable
   from this machine, so I could not check Elmhirst's and Arjevani–Field's stability statements. That at fixed nonzero quadratic coefficient no
   small branch is stable agrees with this instance: every state found is a saddle.

**Fix.**
- Restrict "new" to the identification of the quadratic coefficient as the residual's third derivative (the cohesion and the mean term
  entering only at third order and in the common direction), and to the instance's quantification. Call order preservation a consequence of
  equivariance and the escape time the Riccati solution.
- Add the precedents above to [S63] and the audit, with what was checked.

### ID-11 (cosmetic)

- (v) l. 4414: "σ/ψ̄_i times the cohesion force" should be "σ/(σ₀ψ̄_i) times". The cohesion force is σ₀ψ̄_i(m_i − ϑ_i).
- (iii) l. 4405–4406: the closed form is for ζ "up to the common term", but the departure is measured against the median. For the gap D between
  the escapee and the median, the quadratic coefficient is (|c₃|/2)(ζ + ζ_med)/D ≈ (|c₃|/2)(n − 2)/n when the others are near the median: 20%
  smaller for n = 10. Say that the formula is the leading behaviour, not the measured time.
- (iii) l. 4406: Δs is undefined (see ID-3), and the sentence ends without a period before "A split sheds".
- (C) l. 4449: "on the other — the remaining" is missing "side".
- (D) l. 4456 says "the split delayed to σ₀/32 (its threshold λ⁺ = 0.040σ₀ …)". It should say that the outer classes' thresholds (0.014σ₀,
  0.015σ₀) are not reached by σ₀/64 (ID-1).

## Independent confirmations (what was recomputed and agrees)

**Coefficients (A)** (`coef.out`):
- ϑ^b(q) matches the repository's to 5.7·10⁻¹¹; the drift at the lift is 2.0·10⁻¹⁷ and the fit error 0.016953.
- λ⁺ = 9.9999·10⁻⁵, 2.8607·10⁻⁴, 1.0845·10⁻⁴; the other eigenvalues (−1.27, −3.01, −1.74)·10⁻⁴ and 0 (the amplitude direction).
- v = (0, 0.7484, −0.6633), (0, 0.7159, −0.6982), (0, 0.7849, −0.6196).
- c₃ = 1.41876·10⁻⁵, −1.257488·10⁻³, 1.196565·10⁻⁴; the analytic Hessian differenced at 10⁻³ and 10⁻⁴ agrees to 1.3·10⁻⁵ relative, and the
  6-point stencil gives the same magnitude.
- The thresholds are 0.0424, 0.1213 and 0.0460σ₀.

**(i)'s quadratic coefficient** (`coef.out`):
- The learning drift's even part along v matches (c₃/2)(ξ_i² − S₂/n) to 3.2·10⁻⁸ and 7.9·10⁻⁹ at ε = 10⁻² and 5·10⁻³ (scale 5.2·10⁻⁴; the
  remainder quartic).
- Before the projection the even part differs from (c₃/2)ξ_i² by a constant common to the particles, −2.13·10⁻⁴: the mean term.
- The compact and Gaussian cohesion and the mass-free leak have no even part (≤ 1.5·10⁻¹⁴).
- At the threshold the full Jacobian has a 9-fold zero eigenvalue (≤ 5·10⁻¹⁴); the next one is −1.776·10⁻⁴.

**(ii)/(B)** (`branches.out`, my Newton and analytic Jacobian):
- At μ = 5·10⁻⁷, q = 1–4: positions +0.31%, +0.44%, +0.83%, +3.02% from t*.
- The eigenvalues on the separations along v are (−1.02305; +1.25108 ×8), (−1.67266, −1.06834; +1.67077 ×7), (−2.51785 ×2, −1.19908;
  +2.51652 ×6) and (−5.14261 ×3, −1.94179; +5.14443 ×5) times μ; the full Jacobian has 8, 7, 6, 5 unstable directions. All as in the md.
- At μ = −5·10⁻⁷ every state is a saddle, with q unstable directions.

**(iii) order and time** (`preds12.out`, `C_*.out`): all 12 runs reproduce the md:
- the departure orders and sides, and the first departure's time (2.686·10⁴, 3.093·10⁴, 3.18·10⁴, 3.349·10⁴, 3.30·10⁴, 2.68·10⁴; 6.61·10⁴,
  7.67·10⁴, 8.61·10⁴, 8.414·10⁴, 8.837·10⁴, 6.623·10⁴);
- the clusters at the end, the fit errors and the main pieces' sizes and margins;
- the reduced equation's particle in 12 of 12, its time within 1.15%.

The closed-form Riccati time was re-derived and agrees.

**(v)** (`massfree_static.out`):
- At the lift the mass-free force's Jacobian is exactly −σ on every separation, residual ≤ 1.1·10⁻¹⁹, and maps the class means to ≤ 1.9·10⁻¹⁵.
- The forces sum to zero at the lift, not at spread clouds: (−2.6, 1.3, −0.4)·10⁻⁵ against forces of 6·10⁻⁵, and −6.9·10⁻⁶ for three
  particles on a line.
- ‖G⁻¹J − (G⁻¹J)ᵀ‖/‖G⁻¹J‖ = 0.22–0.23 (the cohesion force's: 3–4·10⁻⁹), so the leak is not a G-gradient field.
- The thresholds are 0.0141, 0.0404, 0.0153σ₀.

**(D) mass-free, seed 0** (`D_mf_0_*.out`):
- Three classes through σ₀/16; 4 classes at σ₀/32 (fit 0.0156) and σ₀/64 (fit 0.0149, centres −1.50, −0.44, +0.07, +1.35, sizes 10, 1, 9, 10),
  as in the md.
- The escapee sits at −0.44 to −0.49, on the side sign(c₃)v (centre lower, width larger: s = 0.32–0.38 against the class's 0.37). It does not
  move toward the bump at +0.8. At σ₀/128 and σ₀/256 the bump is covered by the middle class moving to +0.15–+0.21 and class 3's piece to +1.14.
  The escapee would not reach the bump given more time: the bump is fitted from the shoulders.

**The tests and the md:** 3 passed in 37 s, and the md regenerates byte-identically from the raw file.

## Not checked

- Parts D and E with the replicator (γ = 0.1 and 1) were not rerun; their numbers are taken from the md. The mass-free leak's seeds 1 and 2
  were not rerun.
- Part G of D‴.32 (the cohesion leak's own schedule) was not rerun here; it was reproduced in round 23.
- Elmhirst's and Arjevani–Field's theorems, Ihrig–Golubitsky and Chossat–Lauterbach–Melbourne were not read: arxiv.org was unreachable.
  Brea et al., Şimşek et al. and Cohen–Stewart are cited from memory.
- The Gaussian kernel's split, with the inter-class bias making the lift non-stationary, was not checked beyond the vanishing of its even part.
- The real loop (CL) with the split was not checked.
