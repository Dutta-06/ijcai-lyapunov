# Referee D3ax, round 51: Proposition D‴.49 (the rich transit in the closed loop, as revised after rounds 38 and 44)

**Summary.** The numbers are sound.
- The md regenerates from the raw file with 0 differing lines.
- Every number in the statement and in the Numerics matches the md and the raw, except four values rounded twice and one misquote.
- An independent integration reproduces the traced ejection to every printed digit.
- k_e = 4 stays clean at ε = 10⁻⁶ and 10⁻⁷.
- The desired regressor's result at 1% noise survives a second noise seed.
- At ε = 10⁻⁶ and 10⁻⁷ the exact averaged loop forms exactly one kernel per feature at k_e = 20 and at ω = 1.
- All 20 mutants of the simulation and averaging code are killed (one only by a hang).

What fails is the explanatory layer added in round 44:
- **The "rest condition"** (the active block's rates below k_e) is not the condition the averaging argument needs, and it does not separate failures from successes in the experiment's own raw file. Neither does the within-period swing.
- **Part C's "drive at its nearest kernel"** reports the in-range neighbour's tail, which has the opposite sign to the outside feature.
- **The experiment's own code is still untested.** All 5 mutants there survive, and each one changes the md.

Two of round 44's noise misstatements are still in the text, and (ii) and (iv) each contain a wrong sentence.

## 1. What I checked and how

**Objects**
- review/drafts/D49.md at 4bded0d, unchanged at HEAD 95c3afe.
- exp_closed_transit.py, rich_transit.py, test_closed_transit_loop.py, closed_transit.md and closed_transit_raw.npz are also unchanged since 4bded0d.
- closed_loop.c and closed_loop.py were changed after 4bded0d by 753eb06 (D‴.57 round 50). I read the diff: it only adds options `etrue`, `cbox` and `al0`, whose defaults reproduce the 4bded0d equations. The md and raw predate that change.

**Setup**
- I worked only in a scratch copy at `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref51`: `repo/` plus the mutation trees `mt/` and `mt2/`. Scripts are in `work/` and logs in `logs/`.
- Settings: LCK_BUILD=ref51/build, OMP_NUM_THREADS=1, `nice -n 15`, at most two of my processes at a time, `-B`/PYTHONDONTWRITEBYTECODE, `-p no:cacheprovider`.
- I ran only read-only git commands in the repository, and its working tree stayed clean.
- results/chaos/dead_zone.md is currently empty in the working tree: another session is rerunning exp_dead_zone.py, and the file is marked assume-unchanged. I read D‴.57's Part C from HEAD with `git show`.

**md against raw**
- I ran `exp_closed_transit.py --from-raw` with CT_OUT pointing at a copy: 0 differing lines.
- I dumped every row of the raw: A–E3, B2's trace and rates, C's checkpoints, and D's drive zeros.
- I checked about 150 numbers in the statement, proof sketch and Numerics against both.

**Tests**
- `pytest -q -p no:cacheprovider tests/test_closed_transit_loop.py tests/test_closed_loop.py`: 9 passed in 987 s wall (434 s CPU).

**Mutants**
- 25 single replacements:
  - 11 in closed_loop.c;
  - 2 in closed_loop.py;
  - 7 in rich_transit.py (FrozenOrbit, OrbitBumps, frozen_orbit_rhs, curve_mobility, transit_scale);
  - 5 in exp_closed_transit.py (crossings, report, lip_f, active_rates, c_job).
- Each mutant ran against the tests that import the mutated module, kill-likely tests first, with `-x`.
- For the experiment's mutants I measured the effect on the md: `--from-raw` for the report path, or a recomputation of the affected quantity.

**Targeted numerics** (each under 10 CPU-minutes)
- (a) An independent integration of the traced ejection: scipy LSODA on (x, α, log β, c), with no C code, at k_e = 2, ε = 10⁻⁴, τ ≤ 70 (606 s).
- (b) The repository's `real` at k_e = 4 for ε = 10⁻⁶ and 10⁻⁷, to B's τ = 139.25 (144 s and 484 s).
- (c) The desired regressor at s = 0.01 with noise seed 11 (297 s).
- (d) The exact averaged loop at ε = 10⁻⁶ and 10⁻⁷ for k_e = 20 and ω = 1 (6–8 s each).
- (e) The formed centres' swing in every period of τ ∈ [84, 100] at a = 0.8 (41 s).
- (f) Part C's drives:
  - signed, and split into the outside bump's own part and the rest (exact and leading order);
  - the same drives at the stored τ = 10³ state.
- (g) The symmetry of the leading-order and exact drives' Jacobians, and the eigenvectors of the active block's linearisation.
- (h) Lip f, sup f′, and the unadapted frozen orbit.
- (i) B's k_e = 2 row with a 10× looser error control.

**Literature**
- Web searches (arxiv.org and ams.org are not reachable) for:
  - Sanner–Slotine, Gorinevsky, Wang–Hill, Kingravi et al., Fabri–Kadirkamanathan, Boffi–Slotine;
  - Rohrs et al., Riedle–Kokotović;
  - convergent dynamics (Pavlov et al.);
  - saddle-to-saddle dynamics;
  - small-initialisation or growing-network adaptive control.

## 2. Verdicts

### (i) The exact averaged loop

**Contraction and periodic orbit: correct.**
- ‖L‖_∞ = 1/k_e, the factor is Lip(r)/k_e, and convergence is at rate k_e − Lip(r).
- This is Demidovich's convergence for a scalar system. The one-sided condition k_e > sup r′ already gives a unique, exponentially stable periodic orbit (quasi-periodic for the quasi-periodic reference).
- Here sup f′ = 4.549 and max|f′| = 4.558; the md prints 4.56 and the experiment's docstring says 4.55.
- The unadapted frozen orbit at k_e = 10 has max|e_p| = 0.0894 and rms 0.0370.

**Leading order and its additive error O(|r|²/(k_e²σ)): correct.**

**"Not a gradient": true for the exact drive only.**
- For x_m = sin ωt, and for any reference symmetric under t ↦ t₀ − t, the odd part of the filter L drops out of ⟨φ_i, L r⟩.
- The leading-order drive is therefore ⟨φ_i, L_s r⟩ with L_s = Re H = k_e/(k_e² + ω²), which is the gradient of ½⟨r, L_s r⟩.
- Measured:
  - the antisymmetric part of the leading-order Jacobian is 2.3·10⁻¹⁶ of its symmetric part;
  - g_lead equals the drive computed with Re H to 3.7·10⁻¹⁶;
  - the exact drive's antisymmetric part is 0.126 of the symmetric part at w = 0, and 0.088 with kernel 5 formed.
- So the leading-order loop is D‴.48's transit for a weighted fit, not merely "formally". The non-gradient part is exactly the O(1/k_e²) correction that decides the near ties.

**The active block: the stated condition is not the one the argument needs, and it is not a criterion in the data.** See D3ax-1.
- The formed amplitudes are slow. The active block's slowest modes at the fit are −7.1·10⁻⁴, −6.4·10⁻⁴ and −4.0·10⁻⁴ at k_e = 10, and −3.0·10⁻³ to −1.6·10⁻³ at k_e = 2. Their eigenvectors lie at least 99.4% in α. These are the md's "3 null — the gauge of coinciding kernels".
- Only the centres are fast: rates 0.98–2.70, i.e. 4.9–13.5ω at k_e = 10.
- Classical averaging needs the centres to be slow against the period. The proof's step "converge to its periodic state within a period, i.e. the active block's rates below the error loop's" (D49.md:60–61) does not follow.
- The real-minus-averaged gap in entry times does not shrink with ε for the later entries (k_e = 10):
  - first entry: −0.008, −0.0015, +0.0013 at ε = 10⁻⁴, 10⁻⁵, 10⁻⁶;
  - second entry: −0.036, −0.005, −0.001;
  - third and later entries: +0.18 to +0.29, then +0.32, then +0.36.

### (ii) The real loop

**The numbers hold.**
- The entrants are the same in both loops for every row with k_e ≥ 4.
- Entry times agree within 0.363 at ε = 10⁻⁵ and 10⁻⁶ (at most 0.14%).
- They agree within 0.54 over k_e = 4–10, ω = 0.05 and a = 0.6, 0.8, and within 2.48 at k_e = 20.
- The coefficients agree within 9.4·10⁻⁴; the sums are −0.6000 and 0.7505.

**Two sentences are wrong.**
- **"Five kernels at k_e = 20 and ω = 1"** (D49.md:26): ω = 1 forms four, both in the md and in the Numerics' own (B).
- **"The same coefficients summed over coinciding kernels … a gauge"** (D49.md:23, 86–87): at 0.61 the real loop's two kernels do not coincide. Only the 0.07 pairs are a gauge (Δc = 1.3·10⁻⁵). The real loop lags the averaged loop's merge:

  | case | real loop: centres at 0.61 | averaged loop: centres at 0.61 | real in-range fit rms | averaged orbit-error rms |
  |---|---|---|---|---|
  | ε = 10⁻⁴ | 0.6314, 0.6091 (0.19σ apart) | 0.6101, 0.6100 | 1.4·10⁻⁴ | 1.6·10⁻⁷ |
  | k_e = 20 | 0.6092, 0.6208 | 0.6100, 0.6101 | — | — |

**"One kernel per feature at small masses" is supported, but only outside the md.**
- The exact averaged loop forms three kernels at k_e = 20 (kernels 5, 16, 11) and at ω = 1 (kernels 5, 17, 10), at both ε = 10⁻⁶ and 10⁻⁷.
- At k_e = 20 the kernel identity changes with the mass: kernels 17 and 10, with their partners, at 10⁻⁵; kernels 16 and 11 below that.
- Gaps the draft leaves open:
  - the real loop is not run at these masses;
  - the closed loop's margins are not computed;
  - no separation hypothesis is stated (D3ak-15 is still open, and D‴.52 is not cited).

### (iii)(a) Low gain

**The counterexample holds.**
- My independent integration reproduces B2(b) to every printed digit:
  - kernel 17 forms at τ = 41.118 at 0.696;
  - it is out at τ = 65.839, having moved from 0.601 to 1.238 in the period before;
  - its |w| is 0.816;
  - the state stays in [−0.9945, 0.9950];
  - the other formed kernels swing by 0.415.
- The other entries also match: 9.918, 63.284, 68.609 and 69.337.
- B's k_e = 2 row is unchanged at the md's precision under a 10× looser error control.

**"Clean for k_e ≥ 4" survives a smaller mass.**
- Three kernels form at the features at ε = 10⁻⁶ (entries 50.71, 66.62, 107.20) and at 10⁻⁷ (50.74, 66.61, 107.19).
- None is pushed out. The transit swing is 0.111 and 0.125; the final error is 9.8·10⁻⁹ and 1.8·10⁻⁹.

**Neither the rates nor the swing separate ejecting runs from clean ones.** The table below covers every relevant row in the raw file; rates/k_e is the fastest active rate over k_e.

| run | kernels ejected | max swing | median swing | fastest rate / k_e |
|---|---|---|---|---|
| k_e = 2 | 3 | 0.419 | 0.179 | 4.47 |
| k_e = 3 | 3 | 0.176 | 0.098 | 2.41 |
| k_e = 4 | 0 | 0.126 | 0.105 | 1.49 |
| a = 0.6 | 0 | 0.181 | 0.056 | 0.44 |
| a = 0.8 | 0 | 0.330 | 0.047 | 0.43 |
| f × 2, k_e = 4 | 1 | 0.153 | 0.123 | — |
| f × 2, k_e = 6 | 0 | 0.144 | 0.094 | 2.71 |
| f × 2, k_e = 8 | 0 | 0.097 | 0.067 | 1.45 |
| f × ½, k_e = 2 | 0 | 0.217 | 0.107 | 1.12 |

- My per-period run at a = 0.8 shows the entrant (kernel 5, |w| ≈ 0.08) swinging by 0.357 within one period near τ = 88.1. The 90th percentile over periods is 0.30, and nothing is ejected.
- "Kernel 17, formed at 0.70 with |w| = 0.82" (D49.md:33) merges two moments: the formation at τ = 41.1, and the coefficient one period before the ejection at τ ≈ 65.8.

### (iii)(b) Measurement noise

**The numbers hold, and a second seed confirms the desired regressor at 1%.** With seed 11:
- kernels [4, 11, 16], none out, at centres (−0.544, 0.067, 0.598);
- true-error rms 0.0077 (seed 7: 0.0071);
- in-range fit rms 0.050 (seed 7: 0.058).

The free loop at s = 0.01 ends at 0.0968 (max) and 0.0372 (rms), at or above the unadapted level of 0.0894 and 0.0370.

**Two misstatements flagged by D3aq-12 remain.**
- **"Halves the amplitudes"** (D49.md:40, 114): the mirror leak's coefficients are 0.433, −0.160 and 0.331, i.e. 0.48, 0.27 and 0.44 of the true amplitudes. Two centres moved by 0.054 and 0.097 (0.8σ).
- **"Five at the wrong places"** (D49.md:116–117), for λΛ = 10⁻³ at s = 0.01:
  - three kernels sit at 0.611–0.621, i.e. at the 0.61 feature;
  - one sits at 0.070 with w = −0.029;
  - one sits at 0.155;
  - the −0.52 feature is lost, and the three kernels at 0.61 started at −0.7, −0.6 and −0.3.

**Not reported: the "kept" kernels oscillate within each period.** The fit columns are therefore snapshots.
- Mirror leak at λΛ = 3·10⁻³: end-of-run swing 0.244, 0.242 and 0.247 at s = 0, 0.005 and 0.01.
- Desired regressor at 1%: 0.098 (seed 7) and 0.121 (seed 11).

**The dead-zone citation and "every remedy fails at 2%" need qualifying.**
- "The dead zones of D‴.57 (seven kernels and twelve)" are its desired-regressor variants. With the plain regressor x̃, both dead zones end with 16 kernels, all pushed out.
- "At 2% every remedy tried fails": D‴.57's row combining the desired regressor with the mirror leak keeps three kernels at 2% with none out (fit 0.259, true-error rms 0.020). That is the outcome the text counts as "keeps three kernels" for the leak at 1%.
- D‴.57's Part C is being regenerated with round-50 code, so these citations should be re-read after that rerun.

### (iii)(c) Features near the range

**The formation facts hold.**
- Formed at 0.92σ (τ = 1237) and at 1.6σ (τ = 5344).
- With the sign flip, formed early.
- Not formed at 2.9σ or 3.9σ within τ = 2·10⁴.

**The explanatory column measures the wrong thing.** See D3ax-3.
- "The drive at its nearest kernel before any entry" (exp_closed_transit.py:320) is |g| of the net drive. With the given signs it is negative in all four rows: it is the tail of the in-range bump (−0.6 at 0.07), not the outside feature's pull.

  | distance | exact net drive | leading order: outside bump | leading order: the rest |
  |---|---|---|---|
  | 1.6σ | −0.00203 ("0.188 of the largest") | +0.00129 | −0.00313 |
  | 0.92σ | −3.3·10⁻⁴ | +0.00328 | −0.00313 |

- At 0.92σ the leading-order net drive is +1.4·10⁻⁴, so the two parts nearly cancel and the sign flips between the leading order and the exact drive.
- As a result the column is non-monotone in distance: 0.0305, 0.188, 0.0087, 0.0006.
- The bump is formed by kernel 9, not kernel 8, after the in-range transit. At the stored τ = 10³ state, the drives at kernels 8 and 9, as fractions of the initial largest drive, are positive and ordered like the formation times:
  - 0.92σ: +0.053 and +0.109;
  - 1.6σ: +0.020 and +0.034.

**"Beyond ≈ 3σ … the in-range fit is clean" fails with the sign flip at 2.9σ.**
- Kernel 12 enters at τ = 1.51·10⁴ at 0.120 (w = 0.058), and kernel 11 sits at 0.063 (w = 0.548).
- The orbit's error at τ = 2·10⁴ is 1.30·10⁻⁴, so the Numerics' "2.0·10⁻⁴ … at both signs" holds only to τ = 10⁴.

**The "cluster that compensates on the way" is the end state.** At τ = 2·10⁴ there are still five kernels, three within 0.0055 of 0.07 with coefficients −1.62, −0.20 and +1.22.

**The md contradicts itself on the flipped rows.**
- Its "entry near it" column says "none" while "fitted" says "yes".
- The bump was formed by kernel 11, which started at 0.1, entered at τ = 49.7 or 59.1, and slid 0.26–0.34 out of the range.
- `first` (exp_closed_transit.py:318) only looks at kernels started within 0.15 of the bump.

### (iv) The price

- The scalings hold: τ_{e,1}/k_e = 4.93, 11.86, and 12.65–13.86. Part D holds.
- **"t ≈ k_eΛ(ε)τ_{e,1}" (D49.md:49) is off by a factor k_e.** With τ_{e,1} the first entry's τ, the first entry is at t = Λτ_{e,1} ≈ 14k_eΛ, and the transit ends about 1.9 times later.

### Numerics against the md

All transcribed numbers match, except these:

| draft says | should be | raw value |
|---|---|---|
| 1.1·10⁻⁷ (D49.md:84) | 1.0·10⁻⁷ | 1.049·10⁻⁷ |
| 3.4·10⁻³ (D49.md:104) | 3.3·10⁻³ | 3.345·10⁻³ |
| 0.059 (D49.md:39, 113) | 0.058 | 0.05850 |
| 0.41 (D49.md:115) | 0.40 | 0.40496 |
| 2.0·10⁻⁴ at both signs (D49.md:106) | 1.30·10⁻⁴ for the flipped row at τ = 2·10⁴ | — |

The first four come from rounding the md's three-digit values a second time. The draft also hand-computes several numbers that the md does not print: τ_{e,1}/k_e, the entry-time differences 0.4, 0.6 and 2.5, the "7%", and the pair sums. They all check.

### Round-38 and round-44 fixes

**Present and correct:**
- Round 38: D3ak-1, 3, 5, 6, 8, 9, 12, 13, 16.
- Round 44: D3aq-1, 7, 8, 9, 10, 11, 13, 15, 16.

**Present but incomplete or wrong:**
- D3ak-2 and D3aq-2: the rest condition (D3ax-1).
- D3ak-4: (ii) gives the wrong count at ω = 1.
- D3ak-7 and D3aq-5: the experiment's code is untested (D3ax-2).
- D3ak-10: the "unadapted level" the text invokes is not printed in the md; e_first is in the raw only.
- D3ak-11: "the others left in place" remains in the Status' "New" (D49.md:73). The md shows:
  - inactive displacements of 0.119 (Part A) and 0.205 (real loop, k_e = 3);
  - in the averaged runs at k_e = 3 and 4, inactive kernel 6 sliding 0.12 onto the −0.52 feature with w = 1.4·10⁻³ and 4.2·10⁻⁴.
- D3ak-14 and D3aq-15: some classical sources are still missing (D3ax-10).
- D3ak-15: close features are not treated; still open.
- D3aq-3: "every remedy fails at 2%", and the D‴.57 rows cited are the desired-regressor variants.
- D3aq-4: the drive variable is mismeasured (D3ax-3).
- D3aq-6: the pair at 0.61 does not coincide in the real loop.
- D3aq-12: "halved" and "wrong places" remain.
- D3aq-14: "below 0.009 in every row" is true by definition. The 0.92σ row hides a fourth kernel at 0.61 with w = 0.0085, which is 0.94 of the threshold.
- D3aq-17: τ_{e,1} was introduced, but (iv)'s formula using it is wrong.

### Tests and mutants

**Simulation and averaging code: 20/20 killed (19 by an assertion).**
- R1–R7 (rich_transit), P1–P2 (closed_loop.py) and C1–C5, C7–C10, C12 (closed_loop.c) were killed by an assertion.
- C6 (the leak's sign flipped) was killed only by a hang. I stopped it after 778 s; the test it hangs in takes 160 s on baseline. Round 44 saw the same.
- Two kills are narrow:
  - C10 (10× looser error control) fails only the centres' 10⁻⁶ bound, at 1.09·10⁻⁶. It changes no md number beyond the fourth decimal of one state range.
  - R6 (leading-order filter 10% off) fails only the truncation-gap bound, 0.0099 against 0.0082.

**exp_closed_transit.py: 0/5 killed.** No test imports it, and test_dead_zone uses only its constants, `law` and `th_abs`. Each surviving mutant changes the md:

| mutant | change | effect on the md |
|---|---|---|
| X1 | entries not interpolated | A's "132.9, 157.7, 253.8" becomes "133, 157.8, 254" |
| X2 | the rate cut in the report | 12 lines of B(b) change; the "null" counts become 0 (4 at k_e = 20) |
| X3 | lip_f made one-sided | 4.56 becomes 4.55, twice |
| X4 | the factor 2 dropped in active_rates | every quoted rate doubles (k_e = 10: −5.39 instead of −2.70) |
| X5 | c_job reads the drive one kernel inward | C's column at 0.92σ reads 0.175 instead of 0.0305 |

**Test docstring.** test_frozen_orbit_is_the_periodic_solution announces a "halving r halves e_p" check that its body does not perform (tests/test_closed_transit_loop.py:78).

### Literature

**Classical, confirmed by search unless noted:**
- Frozen-system averaging: Anderson et al.; Bodson–Sastry–Anderson et al.; Sastry–Bodson. Cited.
- Desired compensation: Sadegh–Horowitz. Cited.
- Persistence of excitation along periodic orbits and deterministic learning: Wang–Hill (including ICCA 2005); Kurdila et al. Cited.
- Growing RBF networks:
  - Platt, MRAN and GAP-RBF, as cited by the draft;
  - Fabri–Kadirkamanathan, IEEE TNN 7(5):1151–1167, 1996;
  - Kingravi–Chowdhary–Vela–Johnson, IEEE TNNLS 23(7):1130–1141, 2012.
  - All cited.
- Mirror-descent adaptation laws whose geometry can favour sparsity: Boffi–Slotine, Neural Computation 33(3), 2021. Cited.

**Classical but not cited in D49:**
- Sanner–Slotine, "Gaussian networks for direct adaptive control", IEEE TNN 3:837–863, 1992. This is Gaussian kernels on a lattice in a Lyapunov adaptive controller, i.e. this controller without the free law and the moving centres.
- Demidovich's convergent dynamics (Pavlov–Pogromsky–van de Wouw–Nijmeijer, Systems & Control Letters 52, 2004), for (i).
- The slow-adaptation conditions and their failure, for (iii)(a):
  - Kokotović–Riedle–Praly 1985, and Riedle's 1986 thesis "Integral manifolds of slow adaptation";
  - Rohrs–Valavani–Athans–Stein, IEEE TAC 30(9), 1985.
- Saddle-to-saddle and incremental learning at vanishing initialisation (Jacot et al. 2021; Pesme–Flammarion 2023; Berthier 2023). "Forms its kernels one at a time from tiny masses" is that phenomenon; the repository cites these for D‴.48, but D49's "New" does not.
- Gorinevsky's 1995 persistence-of-excitation result for RBF identification could not be confirmed by search.

**New, as far as the searches show:**
- the closed-loop realisation of the small-mass transit as an implicit growth rule;
- the frozen orbit's O(1/k_e²) non-gradient correction deciding the near ties;
- the low-gain ejection of formed kernels by the centres' motion within a period.

All three are numerical observations on one instance.

## 3. Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ax-1 | major | The "rest condition" (the active block's rates below k_e) is not the averaging condition, and it is not a criterion in the data. The formed amplitudes are slow, but the centres are 4.9–13.5ω. The real-vs-averaged gap for later entries does not shrink with ε. In the raw, f × 2 at k_e = 6 is clean at 2.71k_e while k_e = 3 ejects at 2.41k_e. The swing fails as a criterion too (0.33 and 0.181 clean, 0.176 ejects). | (i), proof, Status: average the inactive duals and the formed amplitudes classically; state the centres' within-period motion as an unbounded error of the averaged loop, measured by the swing; drop "rates below k_e" as a condition, or call it a heuristic and list its counterexamples. (iii)(a): keep "clean for k_e ≥ 4 in this instance" without the ratio. md: print B2's rates and the median swing beside each outcome. (ii): state the ε-independent gap. |
| D3ax-2 | major | No test touches exp_closed_transit.py (D3ak-7 and D3aq-5 only partly fixed). All 5 mutants there survive and change the md, including the rates and the Part C drive the draft uses as mechanisms. | Tests: crossings against event-located entries on a short run; active_rates against a linear right-hand side with known eigenvalues; c_job's drive (signed) against FrozenOrbit; lip_f against A·e^{−1/2}/σ for one bump; report's null count against a synthetic raw; a short b_job end to end. |
| D3ax-3 | major | The quantity behind (iii)(c) and D3aq-4's fix is mismeasured: the drive at the nearest kernel is the neighbour's tail with the opposite sign (1.6σ: −0.00203 = +0.00129 − 0.00313). "Beyond ≈ 3σ the in-range fit is clean" fails with the flip at 2.9σ. The cluster persists. "Entry near it: none" contradicts "fitted: yes". | md C: the signed drive and its split; the drive at the forming kernels after the in-range transit (τ = 10³); the forming kernel wherever it started. (iii)(c) and Numerics (C): restate. |
| D3ax-4 | minor | (ii) says "five kernels at k_e = 20 and ω = 1" (four at ω = 1). It calls the 0.61 pair "coinciding … a gauge" when the real pair is 0.19σ apart, with a 1.4·10⁻⁴ fit error against 1.6·10⁻⁷. | (ii), Numerics (A): "five at k_e = 20, four at ω = 1"; "the real loop's pair at 0.61 has not merged by τ = 600". |
| D3ax-5 | minor | "t ≈ k_eΛ(ε)τ_{e,1}" is off by a factor k_e. | (iv): t_{e,1} = Λτ_{e,1} with τ_{e,1} ≈ 14k_e; the transit ends at about 1.9 t_{e,1}. |
| D3ax-6 | minor | Noise: "halves" and "five at the wrong places" are still in the text; the in-period swings of the "kept" kernels (0.24 with the leak, 0.10–0.12 with the desired regressor) are unreported; the D‴.57 rows cited are its desired-regressor variants; "every remedy fails at 2%" ignores the desired regressor + leak row. | (iii)(b) and Numerics (E): exact fractions; "the −0.52 feature lost, three kernels at 0.61"; an end-swing column in E; name D‴.57's variants; qualify "fails". |
| D3ax-7 | minor | "N null — the gauge of coinciding kernels" mislabels the formed amplitudes' slow modes. | exp_closed_transit.py:447–458: report the amplitude modes and the gauge modes separately. |
| D3ax-8 | minor | (i) understates the leading order, which is an exact gradient for reversible references; the hypothesis could be one-sided. | (i), proof: add both points. |
| D3ax-9 | minor | "One kernel per feature at small masses" has no support in the md where 10⁻⁵ gives 5 or 4 kernels; it lacks the margins and a separation hypothesis. | md B: averaged rows at 10⁻⁶ and 10⁻⁷ for k_e = 20 and ω = 1 (seconds each); (ii): margins beyond the active block's drive error, and well-separated features (cite D‴.52). |
| D3ax-10 | minor | Novelty: "the others left in place" is contradicted by the md; classical sources are missing. | Status: remove "left in place"; cite Sanner–Slotine, Demidovich/Pavlov et al., Kokotović–Riedle–Praly/Riedle, Rohrs et al., Jacot et al./Pesme–Flammarion/Berthier. |
| D3ax-11 | cosmetic | Twice-rounded numbers (1.1·10⁻⁷, 3.4·10⁻³, 0.059, 0.41); "2.0·10⁻⁴ at both signs"; hand-computed derived numbers; the "formed at 0.70 with \|w\| = 0.82" conflation; Status "k_e ≥ 4" against (ii)'s k_e > 4.56; the docstring's 4.55; the test docstring's halving; the tautological sub-threshold sentence. | Numerics, Status, docstrings. |

**Evidence**
- **D3ax-1:**
  - Where it is stated: D49.md:17–18, 30–31, 37, 59–61, 65.
  - Eigenvector α-shares of the slow modes ≥ 0.994 (check2).
  - Entry gaps by ε as in §2 (i), from the raw.
  - The rate and swing table in §2 (iii)(a) comes from the raw's B and B2 rows (their `rates` field).
  - a = 0.8 per period: 0.357 within one period, median 0.056, 90th percentile 0.30.
- **D3ax-2:**
  - 25 mutants; C/P/R 20/20 killed (C6 only by a hang at 778 s).
  - X1–X5 survive; their md diffs were produced with `--from-raw` or by recomputation (`logs/mutants_summary.txt`).
- **D3ax-3:**
  - exp_closed_transit.py:318 and :320 (`first`, `drive_near = |g[j]|`); D49.md:44–47, 103–107.
  - Signed split from check3/check1: 0.92σ exact −3.3·10⁻⁴, leading order +0.00328 − 0.00313; 3.9σ: +2.6·10⁻⁹ against −3.6·10⁻⁶.
  - At τ = 10³ the drives at kernels 8 and 9 are +0.053/+0.109 (0.92σ) and +0.020/+0.034 (1.6σ).
  - Flipped 2.9σ: entries 11 at 64.8, 17 at 92.5, 12 at 15089.6; errors 1.95·10⁻⁴ → 1.30·10⁻⁴.
- **D3ax-4:**
  - md B: ω = 1 kernels [5, 10, 11, 17].
  - Raw A at ε = 10⁻⁴: real centres 0.6314/0.6091, fit_rms 1.39·10⁻⁴; averaged 0.6101/0.6100, rms 1.6·10⁻⁷. B at k_e = 20: 0.6092/0.6208.
- **D3ax-5:** D49.md:49. At k_e = 10, τ_{e,1} = 132.8 while the last entry is at 254.1.
- **D3ax-6:**
  - Raw E: (0.433, −0.160, 0.331) at (−0.472, 0.170, 0.603); E2 λ = 10⁻³, s = 0.01: centres 0.611, 0.621, 0.611, 0.155, 0.070.
  - Raw `swing` 0.244/0.242/0.247; E3 regressor 0.098; seed 11: 0.121.
  - D‴.57 md (HEAD) Part C: x̃ dead zones 16/16; x_m + leak at s = 0.02: 3 kernels, 0 out, fit 0.259, rms 0.0203.
- **D3ax-7:** X2 regenerates the md with "0 null" everywhere except k_e = 20 (4) — that is, the gauge alone.
- **D3ax-8:** check2: 2.3·10⁻¹⁶ and 3.7·10⁻¹⁶; exact 0.126 and 0.088; sup f′ = 4.549 against max|f′| = 4.558.
- **D3ax-9:** averaged loop at 10⁻⁶ and 10⁻⁷: k_e = 20 → [5, 11, 16]; ω = 1 → [5, 10, 17]; fits exact.
- **D3ax-10:** D49.md:72–74; raw B k_e = 3 and 4, kernel 6 at −0.5216 and −0.5216 with w = 1.4·10⁻³ and 4.2·10⁻⁴.
- **D3ax-11:** the table in §2 "Numerics against the md"; D49.md:33, 66, 117; exp_closed_transit.py:8; tests/test_closed_transit_loop.py:78.

## 4. Lesson

The round-44 revision answered every issue with a number, and the numbers are right. The md regenerates exactly, and the one traced event reproduces in an independent integration. What failed is the layer that explains them. Each indicator was read off the rows that suggested it and never checked against every row the experiment already held:
- the active rates against k_e;
- the within-period swing;
- the drive at the nearest kernel.

The f × 2 rates (2.71k_e, clean), the amplitude rows (swing 0.33, clean) and the signs of Part C's drives were all in the raw file. When a proposition names a quantity that separates failure from success, it should tabulate that quantity for every run beside the outcome, keep its sign, and split it into its sources. An averaging claim should also be tested for the dependence of its error on the small parameter. Here the first entries converge as ε shrinks, but the later ones keep a gap of about 0.3 at every ε. That is the signature of an error carried by the active centres, not by the mass, and the condition that makes it harmless is still unknown.
