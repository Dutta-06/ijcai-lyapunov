# Referee report D3bd, round 57: Proposition D‴.59 (the transit with moving widths)

**Scope.** I reviewed commit d0ad068 from a scratch copy made with `git archive`:
- `review/drafts/D59.md`;
- `experiments/exp_moving_width.py`;
- `src/core/moving_width.py`;
- `src/core/rich_transit_nd.py`: `sliding_iss_nd`, `critical_rhs_nd`, `critical_jac_nd`, `loop_nd`, `polish_block_nd`;
- `src/core/rich_transit.py`: `critical_w`, `_alpha_of`, `relax_chunk`;
- `results/chaos/moving_width.md` and `moving_width_raw.npz`;
- `tests/test_moving_width.py`.

Other sessions moved HEAD to cfdf29c during the review. Among these files, the only change is round 54's change to the 1-D `polish_block`/`sliding_iss`, which the moving-width path never calls. So everything below also holds at HEAD.

I modified nothing in the repository; the working tree is clean. All runs were in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref57/d3bd_work/`, with scripts in `scripts/` and logs in `out/`. I used `nice -n 10`, `OMP_NUM_THREADS=1`, and never more than two processes at a time.

## Summary

**What holds.**
- **Reproducibility is exact.**
  - `--from-raw` into `MW_OUT` is byte-identical to the committed md.
  - Fresh reruns are bit-identical to the raw: A at σ*/σ₀ = √2 and 4 (both paths) and 0.25 (moving), D at d = 0 and 0.6, and G at Δ ≥ 0.4, q = 9 (the flow and both loops).
  - 6 tests pass in 11 s.
- **(i) holds in the continuum.** J, its gradient, the unique critical point and the limits re-derive. The convergence argument is right up to minor gaps.
- **The dipole asymptotics are stronger than the draft claims.** On all six G targets:
  - the moment p is constant and equals the optimal dipole fit's p to 6 digits;
  - |Δc| and |Δσ| each fall as s^−1/3, so the exponent is not an artefact of the distance metric;
  - rms − rms_inf falls as s^−0.6667;
  - μ = Gβa′w′ = 0.1|w| exactly, already at |w| = 0.57;
  - the prefactor closes: v³ = 3(4Gc/A)K·s, with K = −⟨r_dip, ∇³φ:p³⟩/24, holds to 4–5 digits.
- **The six C budget stops are the dipole.** I verified this by instrumenting C's own path. E's single unstopped target is the same dipole.
- **(iii) holds per draw.** On 18 of 20 E targets, the stopped fit is the least-squares fit started at the truth; the errors are equal to 6 digits.

**What does not hold.**
- **(ii)(a)'s moment bookkeeping misses the mechanism the data use.** The pair's width difference absorbs the third moment: σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), matched to 4–5 digits on all three separated pairs. Without it, the draft's own unequal-share example (Δ ≥ 0.4, q = 1, residual 2.5·10⁻¹³) would sit at 7·10⁻¹¹.
- **The pairs are not D‴.54 near ties.** They come from sequential entries.
- **(ii)(b)'s prefactor is off by a factor 2,** and κ comes from the third-order Taylor term only.
- **"Lands … an exact fit … a frozen cancelling pair" is not shown.** At τ = 400, only 1 of the 6 G loops has landed. The other five are dipoles still in transit, and they are still moving at τ = 1600–3200. Which feature a dipole lands on depends on ε.
- **"A target of one sign cannot make them" is false.** Two close positive features give negative kernels in 10 of 10 variants. They also produce non-dipole budget stalls, multi-kernel cancellations with |w| ≈ 3.2, and sub-sample widths. (ii)'s dichotomy holds only on C's target distribution.
- **(iii) omits both of its exceptions.** F's target 2 loop is not a gauge copy of the path.
- **C reports the state before the stalled block.**
- **Tests:** 25 of 55 single-replacement mutants survive both the tests and `--from-raw`.
- **Misquotes:** two numbers are misquoted.

## What I checked and how

**1. Reproducibility.**
- `--from-raw` compared with `cmp`.
- `scripts/rerun_AD.py`: `a_job(√2)`, `a_job(4)`, the moving path at 0.25, `d_job(0.6)`, `d_job(0)`.
- `scripts/rerun_G.py`: `g_job` at Δ ≥ 0.4, q = 9.
- All max |diff| = 0.
- pytest with `-p no:cacheprovider`: 6 passed.

**2. The raw data** (`scripts/rawlook.py`, `scripts/rawlook2.py`).
- Pairs: their positions, widths and coefficients.
- Stopped paths: their stored state.
- E, F and G in full.

**3. (i).**
- Re-derived J, ∂_cJ, ∂_σJ, the critical point and the limits by hand.
- Computed J_N on the actual 321 samples, over c ∈ [−3, 3] × σ ∈ [10⁻³, 10²] at ratios 0.25, 1.031 and 4, together with its σ → 0 and σ → ∞ limits (`scripts/jn_samples.py`).

**4. (ii)(a)** (`scripts/pairs3.py`, `scripts/c_pairs_hist.py`).
- For each stored pair, compared the second and third moments against the common width and the width split.
- Recomputed the residual with the widths equalised, and re-fitted with the pair's widths tied.
- Traced the path histories of C's three pair targets.

**5. (ii)(b).**
- `scripts/dipole_flow.py`: G's flow with the full state at 17 log-spaced s; the moment; the |Δc| and |Δσ| slopes; the optimal dipole fit by least squares (the infimum rms_inf).
- `scripts/dipole_law.py`: K = v²(E − E_inf) at s = 10⁹ and 10¹¹; v³/s; K from the 2nd- and 3rd-order Taylor terms at the optimal dipole; μ(w)/|w| at |w| = 0.57–4.
- I tried to continue past s = 10¹¹. LSODA failed on restart ("Unexpected istate"), and Radau failed in both absolute and log time ("Required step size is less than spacing between numbers"). So 10¹¹ is the reachable horizon.

**6. C's stalls.**
- `scripts/c_stall.py` instruments C's own path: every chunk, every polish, and the state when the budget runs out.
  - Δ ≥ 0.4, q = 9 with the full 10⁶ budget (362 s).
  - The other five targets with a 4·10⁴ budget.
- `scripts/polish_status.py`: the polish's `least_squares` status.

**7. The real loop.**
- `scripts/g_loop.py`: τ up to 1600 at ε = 10⁻⁸ and up to 3200 at ε = 10⁻⁶. At each τ I recorded the formed kernels, the closest pair's distance, its coefficients and centres, and the residual. I also took τ = 400 for the other five targets.
- `scripts/f_check.py`: F's target 2 at 1.3× and 6× the last entry.

**8. (iii).**
- `scripts/e_check.py` reproduces E's targets and noise draws from the seeds. Per draw, it fits the 3-kernel moving-width model by least squares starting from the truth, then compares errors.
- `scripts/e12.py` traces the unstopped target.

**9. One-sign targets.**
- `scripts/onesign.py`: 24 random targets, C's generator with |A|, seed 591.
- `scripts/onesign_adv.py`, `scripts/onesign_detail.py`, `scripts/onesign_asym.py`: close positive pairs at ±0.08 … ±0.2, symmetric and perturbed.
- `scripts/c_stall_custom.py`: the instrumented one-sign stall.

**10. Mutants** (`scripts/mutants.py`, `scripts/run_mutants.py`, log `out/mutants.log`).
- 55 single replacements in `pairs_of`, `PAIR_TOL`, `n_formed`, `outcome`, `dipole_exponents`, `g_job`, the report's C, D and G sections, and part B's check.
- Each one run against the tests and against `--from-raw`.

**11. Numbers and literature.**
- Every number in the numerics paragraph was compared with the md and the raw.
- Literature from web-search snippets; items from memory are marked [memory]. ftp1.cs.wisc.edu, which hosts de Boor's γ-polynomials paper, is not reachable from this session, so I relied on its snippet.

## Verdicts per part

**(i): true in the continuum for one feature; the proof is essentially complete, but the claim's scope is overstated.**
- The formula, the gradient (∂_cJ = (2d/v)E₁, ∂_σJ = −(1/σ − 2σ/v + 2σd²/v²)E₁), the unique critical point and the limits all check.
- The flow is ẇ = −μ(w)∂_wE, θ̇ = −Gm∂_θE, with μ > 0 for w ≠ 0. Łojasiewicz applies once the sublevel set is shown to be compact in (w, c, σ) with |w| bounded away from 0. That follows from E ≥ J, as Lemma D‴.53a(ii) shows for fixed widths; the proof skips this step.
- E has no critical points at w = 0, since ∂_wE = −⟨f*, φ⟩ ≠ 0. The "w = 0 critical points" in the proof are equilibria of the degenerate metric.
- **On the samples the argument does not transfer.** J_N's σ → ∞ limit is 0.47–0.97 of ½⟨f²⟩ and its σ → 0 limit at a sample is 0.81–0.99, so the sublevel sets just below ½‖f*‖² are not compact. J_N also has sub-sample flat valleys at σ ≈ 1–1.7·10⁻³, at 0.64–0.98 of ½⟨f²⟩. A's nine exact fits are therefore the evidence on the samples.
- **B adds nothing about the samples.** It is the continuum closed form, and it is scale-invariant: on B's grids, J/σ* is the same landscape at all nine ratios.
- **"One kernel per isolated feature at any width" is proved only for a single feature.** With several features, 9 of 80 C targets do not get one kernel per feature.

**(ii)(a): the pairs are real, but the bookkeeping is incomplete and the mechanism is misnamed.**
- In each pair the net coefficient equals the feature's amplitude (q = 1: 0.778261 vs 0.778260; q = 7: −0.429651 vs −0.429652; D 0.6: 0.800001).
- The second moment checks: for q = 1, the common width's shift σ* − σ̄ = 4.103·10⁻⁷ equals ½(w₁w₂/W²)Δc²/σ.
- The third moment is absorbed by the width split, which the draft omits:

| pair | Δc/σ | Δσ observed | Δσ predicted | residual stored | residual, widths equalised |
|---|---|---|---|---|---|
| C q = 1 | 3.5·10⁻³ | +7.9453·10⁻⁸ | +7.9449·10⁻⁸ | 2.5·10⁻¹³ | 7.05·10⁻¹¹ |
| C q = 7 (cancelling) | — | −1.3404·10⁻⁷ | −1.3405·10⁻⁷ | 1.5·10⁻¹³ | 4.34·10⁻¹⁰ |
| D d = 0.6 | — | −3.0498·10⁻⁷ | −3.0516·10⁻⁷ | 1.9·10⁻¹³ | 1.96·10⁻¹⁰ |

- For the cancelling pair, equal shares are impossible, so the width split is the only way to fit. A common-width re-fit stays at 1.0·10⁻¹¹ (q = 7) and 1.3·10⁻¹⁰ (D 0.6).
- The code's own comment at `exp_moving_width.py:53-54` ("O(Δc⁴)") is right; the draft's proof is not.
- **The C pairs are not near ties.**
  - q = 2: entries at τ = 14.141 and 18.491, Δτ/τ = 0.31. After block 2, kernel 10 sits at (0.077, 0.166) and kernel 11 at (0.066, 0.515). They coalesce only at block 4.
  - q = 1: Δτ/τ = 0.41. Kernel 13 sits between two features at (0.212, 0.095) until the last feature is fitted.
  - q = 7: a positive entrant at τ = 105.4 lands on the −0.43 feature.
- "The split is a gauge" holds only at coincidence; q = 2 (Δc = 2.8·10⁻⁸) is the only one.

**(ii)(b): the dipole is confirmed and the exponent is right; the prefactor is off by 2.**

| target | v³/s (measured) | 3μ₀K | draft's law (6μ₀K) | \|p_σ\|/\|p\| |
|---|---|---|---|---|
| Δ ≥ 0.4, q = 9 | 7.6794·10⁻⁴ | 7.6793·10⁻⁴ | 1.536·10⁻³ | 0.96 |
| Δ ≥ 0.4, q = 18 | 5.3390·10⁻⁴ | 5.3389·10⁻⁴ | 1.068·10⁻³ | 0.99 |
| Δ ≥ 0.6, q = 12 | 6.3211·10⁻⁵ | 6.3168·10⁻⁵ | 1.263·10⁻⁴ | 0.96 |
| Δ ≥ 0.6, q = 20 | 3.5087·10⁻⁵ | 3.5087·10⁻⁵ | 7.018·10⁻⁵ | 0.73 |
| Δ ≥ 0.6, q = 27 | 1.8328·10⁻⁴ | 1.8327·10⁻⁴ | 3.666·10⁻⁴ | 0.99 |
| Δ ≥ 0.6, q = 29 | 4.0340·10⁻⁴ | 4.0340·10⁻⁴ | 8.068·10⁻⁴ | 0.99 |

- **The moment is constant.** The log-slope of |p| is below 10⁻⁵ by s = 10⁹. For q = 9, p_flow = (−0.100307, 0.351904) equals p_opt.
- **Each component scales on its own.** Over 3.2·10¹⁰ → 10¹¹ the slopes are |Δc| −0.3333 and |Δσ| −0.3333, with |w| at +0.3330.
- **The residual approaches its infimum as s^−2/3.** The slope of rms − rms_inf is −0.6667.
- **K is constant.** K = v²(E − E_inf) agrees between s = 10⁹ and 10¹¹ to 4 digits.
- **K is entirely third order.** K₃ = K to 4–5 digits, while the second-order term is 10⁻⁹–10⁻¹². It vanishes at the optimal dipole by stationarity in θ̄, so "κ … from the second- and third-order terms" is wrong in part.
- **The factor 2.** The correct law is dv/ds = μK/v³, which gives v³ = 3μ₀Ks and matches the table to 4–5 digits. The draft's "2μκp²/w³" applies an individual coefficient's mobility to the half-difference's gradient, which gives 6μ₀K.
- **The mobility.** μ/(μ₀|w|) = 1.0000 at |w| = 0.57, 0.6, 0.9 and 4.0, as the draft says.
- **The dipoles are width dipoles.** Five of the six are 96–99% in σ: a difference of Gaussians of coalescing widths whose limit is the σ-derivative, a Mexican hat. All five are on alternating-sign triples (+−+ or −+−).

**The C stalls are dipoles, but the path stalls at its polished state, not on G's flow.**
- G's flow reaches s = 10¹¹ in 8,115–12,363 RHS evaluations.
- C's path behaves differently. Its polish after the stalled block's first 10⁵ chunk (2,900–4,400 trf evaluations, stopped on xtol/ftol with status 2–3) is accepted because the block has no seed-level kernels.
- The polished states are far beyond the flow's reach. |w| is 3910, 3964, 1821, 1431, 3002 and 3956, with separations 1.9–6.1·10⁻⁴ and rms equal to the dipole infimum to about 10⁻⁹. By the s^1/3 law, that corresponds to s ≈ 10¹⁴.
- For q = 9, the next chunk alone spent 994,226 evaluations without visibly moving the state.
- The polished |w| is numerically arbitrary. A replication whose start differs at about 10⁻¹⁶ (α rebuilt by `_alpha_of`, as `g_job` does) lands at |w| = 3486 instead of 3910.
- C's table and its "Not exact" list report the state before the block, with the entrant at w = 1.0·10⁻⁷: "2 kernels, residual 2.7e-01". The C table's max residual (2.7·10⁻¹, 2.2·10⁻¹) and its formed-kernel minimum of 1 are these pre-block values.
- **Other stalls.** C has none: all six budget stops are dipoles and nothing is "not fitted". E's q = 12 is the same dipole: its second block polishes to w = ±1821 and the budget runs out at block rms 0.216 against 1.1s = 0.0011. A one-sign close pair, however, gives a **non-dipole** budget stall: its 4th block's polish is rejected, and LSODA spends 139,516 evaluations without moving a seed-level entrant (w ≈ −2·10⁻⁷, max|w| 1.02). So `outcome`'s label "stopped by the budget" (glossed in its docstring as "a block that never settles: G") is not the dipole in general.

**The real loop's landing: not shown; "exact" and "frozen" are contradicted.**
- At τ = 400 and ε = 10⁻⁸, the "coinciding pairs" have distances d = 4.6·10⁻⁴ (q = 9), 2.1·10⁻² (q = 12), 2.1·10⁻² (q = 18), 3.0·10⁻² (q = 20), 7.5·10⁻³ (q = 27) and 1.9·10⁻² (q = 29).
- Their midpoints sit 0.008σ* (q = 9) and 0.09–0.21σ* (the other five) from the feature they "share", with residuals 2.2·10⁻⁷–1.4·10⁻⁶, 22–140 times the draft's exact threshold.
- PAIR_TOL = 0.05 is ten times the largest C/D pair distance (0.0049). It counts dipoles still in transit as coinciding pairs.
- q = 18 does not freeze. From τ = 400 to 800 to 1600 at ε = 10⁻⁸, d goes 0.0213 → 0.0184 → 0.0161 and the residual 1.38·10⁻⁶ → 7.6·10⁻⁷ → 4.5·10⁻⁷.
- At ε = 10⁻⁶, q = 18's dipole lands on a *different* feature: c = 0.541 rather than −0.017. At τ = 3200 it has d = 0.018 and residual 8.8·10⁻⁷.
- At ε = 10⁻⁶, q = 20 forms no dipole at all. Instead the −0.507 feature splits into a same-sign pair (c = 0.627/0.690 at τ = 3200), with residual 1.2·10⁻⁵ at τ = 400 and 1.7·10⁻⁶ at τ = 3200.
- So τ = 400 is not long enough, and the configuration being approached is reached only algebraically, if at all.

**(iii): true per draw on 18 of 20 targets; both exceptions are omitted.**
- E, per draw: on 18 targets the path's error equals the oracle fit's to 6 digits.
- Target 2 (the coinciding pair) stops with 4 kernels, at 1.049× the oracle fit at s = 10⁻³ and 1.267× at s = 10⁻².
- Target 12 never stops: it is the dipole.
- The draft compares the median with s(3K/N)^1/2, which is (E err²)^1/2, so the median is the wrong statistic. The rms of the 19 errors is 1.648·10⁻⁴ and 1.750·10⁻³, against 1.674·10⁻⁴ and 1.674·10⁻³. The s = 10⁻³ median of 1.42·10⁻⁴ is χ²₉ skew, not a gain on the oracle.
- F's target 2 is not a gauge of the path's pair.
  - At 1.3× the last entry: the loop's pair has Δc = 1.5·10⁻³, Δσ = 5.5·10⁻⁴ and residual 5.2·10⁻⁷, against Δc = 2.8·10⁻⁸ and residual 1.0·10⁻¹⁵ on the path.
  - At 6×: Δc = 2.2·10⁻³, Δσ = 3.9·10⁻⁴, the split frozen at 0.329/0.489, residual still 5.0·10⁻⁷, while the inactive kernels' coefficients are 5.4·10⁻⁶.
- On the dipole targets the path never ends, so "the loop forms the path's kernels" cannot apply there.

**One-sign targets: the claim is false and (ii)'s dichotomy is unqualified.**
- 24 of 24 random one-sign targets at Δ ≥ 0.4 end exact with no negative entrant.
- Two positive features at ±0.08 tell a different story (Δ = 0.16 = 1.33σ₀ = 2.7–4σ*).
  - σ* = 0.04, equal amplitudes: six symmetric tie pairs of negative kernels enter (down to −0.43). The path ends with a null triple at (0, 0.399): +0.109, −0.055, −0.055.
  - σ* = 0.06, equal amplitudes: 17 entries. The end has w = +3.24 against six negative kernels (−0.06 … −0.61) and five kernels at σ = 1–9·10⁻⁴, residual 4.3·10⁻¹⁰. `outcome` calls this "not fitted".
  - Eight perturbations (amplitudes 0.99/0.95/0.9, a shifted centre, a changed width): all eight give negative kernels (−0.33 … −0.64). Six are stopped by a 1.5·10⁵ budget and two are "not fitted" (|w| 3.1–3.3, σ_min 5·10⁻⁵ and 2.4·10⁻⁴).
  - With amplitudes (1, 0.6), every variant is exact.

**Numerics paragraph.** Every number matches the md or raw except:
- line 76: "1.8–5.9" should be 1.8–6.0 (the raw is 5.9512);
- line 75: "0.035–0.073" should be 0.035–0.072 (the raw is 0.072467, double-rounded through the md's 7.25e-02).

Two phrases are not measurements:
- "fixed widths: 4 kernels" (line 72) is the MAX_ENT_FIXED cap.
- "forms … with a cancelling coinciding pair" (line 77) depends on PAIR_TOL; see above.

**Tests and mutants.**
- 16 of 55 mutants are killed by the tests and 25 by `--from-raw`.
- 25 survive both. Two of those (P3, RD1) are effectively equivalent.
- Survivors in full:
  - all five `g_job` mutants (no test runs it);
  - N1, N2, N4 in `n_formed` (it has no unit test);
  - in `outcome`: O2 (EXACT_TOL → 10⁻³), O3 (the pair's rms condition dropped), O4 (pairs among all kernels), O5 (formed by the max amplitude), O6 (exact rms 10⁻⁸ → 10⁻⁶);
  - P6 (/S0 in place of /sm) and P7 (PAIR_TOL 0.05 → 0.5);
  - RD2 and RD4 in the D report;
  - in part B: B2 (2 points), B3 (absolute error), B4 (last point only), B5 (points near the centre), B7 (non-strict minima), B9 (c-grid ±σ*).
- `dipole_exponents` (X1–X5) and most of G's report are killed.

**Literature and novelty.**
- **Present and correct:** Golub–Pereyra (SIAM J. Numer. Anal. 1973; Inverse Problems 2003 review); Braess's book (1986).
- **Missing: the non-closedness and runaway precedents.**
  - de Boor (1969) on γ-polynomials: the family is not closed for n > 1 and is closed by adjoining limits of divided differences. That covers the dipole *and* higher confluent terms.
  - Kammler (thesis 1971; 1976 papers) on existence of best approximations.
  - Petersen–Raslan–Voigtlaender (FoCM 2021): fixed-size network classes are not closed in L^p, and approximating a closure point forces the parameters to explode. This is the modern statement of the dipole runaway.
- **Missing: the σ-dipole.** Marr–Hildreth: a difference of Gaussians approximates the Laplacian of Gaussian.
- **Missing: free-width fitting and mixture degeneracies.**
  - Moving-width kernel learning: Poggio–Girosi HyperBF (1990), with moving centres and adjustable widths by gradient descent.
  - CL-OMPR (Keriven–Bourrier–Gribonval–Pérez, 2016–18): greedy atom selection plus joint continuous optimisation, a sliding-Frank–Wolfe variant [the optimisation of the variances is from memory]. It is the closest algorithmic analogue of the moving-width transit.
  - ADCG (Boyd–Schiebinger–Recht).
  - Chizat's conic particle gradient descent (Math. Program. 2022), where θ = (c, σ) is just a particle position.
  - Mixture-likelihood degeneracies: Kiefer–Wolfowitz 1956; Day 1969; Biernacki–Chrétien 2003, where EM is attracted to degenerate solutions.
- **Mis-cited: "Duval–Peyré's two-spike limits".** Their thin-grid result is that each spike is recovered as a *pair of neighbouring spikes* (Duval–Peyré 2017). That is the analogue of (ii)(a)'s split, not of the dipole. [memory: in the BLASSO the TV norm excludes |w| → ∞, so there the dipole cannot be a minimiser; here it arises because each block relaxes to the unpenalised fit.]
- **"New" is too broad.** The one-kernel landscape is elementary, and adaptive widths are standard. What is new is the free-law runaway law v³ = 3(4Gc/A)K·s, the pair's width split, and the ε-dependent landing.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bd-1 | major | "Dipoles need opposite-sign kernels: a target of one sign cannot make them" (D59.md:30) is false. Two positive narrow features at ±0.08 (σ* = 0.04, 0.06; Δ = 1.33σ₀) give negative kernels in 10 of 10 variants (down to −0.64): a null triple; w = +3.24 cancelled by six negatives with five kernels at σ = 1–9·10⁻⁴ ("not fitted"); and, under perturbation, 6 budget stops that are not dipoles. (ii)'s "On separated features the path ends exact or in one of two configurations" (line 15) never defines "separated" and is only observed on C's distribution (71 exact, 3 pairs, 6 dipoles of 80). | (ii): define "separated" (C's Δ ≥ 0.4, widths ≤ 0.3) and state the dichotomy as observed on that distribution; replace line 30 with the data (24/24 random one-sign targets clean; close one-sign pairs are not). Add a part H with the ±0.08 pair and its eight perturbations. `outcome`: a class for multi-kernel cancellations and sub-sample widths. |
| D3bd-2 | major | "Lands on a feature as a cancelling coinciding pair — an exact fit" (lines 28–30), "a frozen cancelling pair" (line 61) and "forms four kernels on all six with a cancelling coinciding pair" (line 77) are not supported. At τ = 400, ε = 10⁻⁸, only q = 9 has landed (d = 4.6·10⁻⁴, rms 1.1·10⁻¹¹). The other five are pairs in transit: d = 0.0075–0.030, midpoints 0.09–0.21σ* off the feature, rms 2.2·10⁻⁷–1.4·10⁻⁶. q = 18 is still moving at τ = 1600 (d 0.021 → 0.016, rms 1.4·10⁻⁶ → 4.5·10⁻⁷). At ε = 10⁻⁶ the q = 18 dipole lands on another feature and q = 20 forms no dipole. PAIR_TOL = 0.05 (`exp_moving_width.py:53`) is ten times the largest C/D pair distance (0.0049). | (ii)(b) and Status: "transient; it migrates toward a feature, the residual decaying algebraically; which feature depends on ε"; drop "exact" and "frozen". `g_job`: record d, the midpoint's offset in σ* units and the residual at τ = 400, 1600, 6400; tighten PAIR_TOL to the C/D scale (≤ 0.01) or report d. |
| D3bd-3 | major | (ii)(a)'s bookkeeping (lines 16–18, 42–44) says only the second moment is absorbed and that the third vanishes for equal shares, implying O(Δc³) errors for unequal shares. The data use a width split, σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), matched to 4–5 digits on q = 1, q = 7 and D 0.6. Equalising the widths raises the residuals from 1.5–2.5·10⁻¹³ to 7·10⁻¹¹–4·10⁻¹⁰. So the draft's own quoted case (q = 1, unequal shares, 2.5·10⁻¹³) contradicts its proof. | Proof (a): the common width absorbs ½(w₁w₂/W²)Δc²/σ and the width split absorbs the third moment, so any split (and a cancelling pair, where equal shares are impossible) fits to O(Δc⁴). C/D report: print Δc, Δσ and the prediction per pair. Add a test on a synthetic pair. Status: "checked on the stored pairs". |
| D3bd-4 | major | (iii) (lines 32–34, 80–83) omits its exceptions. E's target 2 stops with 4 kernels at 1.05× and 1.27× the oracle; target 12 never stops (dipole: polish to ±1821, out of budget at block rms 0.216 ≫ 1.1s). "Reaches the oracle" compares medians with (E err²)^1/2; the rms of the errors is 1.648·10⁻⁴ and 1.750·10⁻³. F's target 2 is not a gauge of the path: the loop's pair has Δc 1.5 → 2.2·10⁻³ and Δσ 3.9–5.5·10⁻⁴, residual stuck at 5·10⁻⁷ up to 6× the last entry, against 1·10⁻¹⁵ on the path. | (iii): "on 18 of 20 targets the stopped fit is the least-squares fit started at the truth (equal to 6 digits); a coinciding pair adds a noise kernel; a dipole never reaches the discrepancy". `e_job`: the per-draw oracle fit and `oob`; the md reports rms as well as median. F: restrict to paths that end, report the loop's residual and the pair's Δθ, extend τ. |
| D3bd-5 | minor | The dipole law's prefactor is off by 2, and κ's origin is misstated (lines 26–28, 46–50). Measured: v³/s = 3(4Gc/A)K with K = v²(E − E_inf) = −⟨r_dip, ∇³φ:p³⟩/24, to 4–5 digits on all six (e.g. 7.6794·10⁻⁴ vs 7.6793·10⁻⁴). The draft's 2μκp²/w³ gives 6μ₀K. The second-order term vanishes at the optimal dipole (K₂ ≈ 10⁻⁹). | Proof (b): v = (w₁ − w₂)/2, v̇ = −(μ/2)∂_vE, so dv/ds = μK/v³, v³ = 3μ₀Ks and E − E_inf = K/v² ∝ s^−2/3. G: print K, v³/s, 3μ₀K, the component slopes and the residual's slope (−0.6667). Note the dipoles are 96–99% in σ (a difference of Gaussians → the σ-derivative). |
| D3bd-6 | minor | "Same-sign (a split: D‴.54's near tie)" (line 16) is not what happened. The C pairs come from sequential entries (Δτ/τ = 0.31, 0.41): a compromise kernel collapses onto an already fitted feature at the last block. q = 7's cancelling pair is a late opposite-sign entrant. "The split is a gauge" (line 18) is exact only at Δθ = 0, which holds only for Δ ≥ 0.6, q = 2. | (ii)(a): describe the mechanism as observed; "a gauge at coincidence; a little apart, the width split ties the split to Δc". Keep "near tie" only for exactly symmetric targets (D d = 0). |
| D3bd-7 | minor | Scope of (i): the title, line 12 "one kernel per feature" and line 59 "hence one kernel per isolated feature at any width" are proved only for one feature on ℝ. "Isolated" is undefined, and "hence" is a non-sequitur for several features (9/80 C targets). On the samples, J_N's σ → ∞ limit is 0.47–0.97 and its σ → 0 limit 0.81–0.99 of ½⟨f²⟩, so the sublevel sets are not compact; there are also sub-sample valleys at σ ≈ 10⁻³. B (lines 68–69, 52–53) checks the continuum closed form, a single landscape up to scale at all nine ratios, not the samples. | Title and (i): "one Gaussian feature, one kernel, at any width (the continuum)". On the samples add D‴.53a(ii)'s sublevel condition including the σ-boundaries. B: J_N on the samples (as in `scripts/jn_samples.py`) or drop "B" from the evidence. |
| D3bd-8 | minor | Proof gaps and wording. (i): compactness in (w, c, σ) and |w| bounded away from 0 (E ≥ J) are skipped; "critical points with w = 0" are the degenerate metric's equilibria, not critical points of E. (ii)(b), line 47: "the mobility Gmw_i against the gradient w_i⟨r∇φ_i⟩" is wrong (it is Gm against ∂_θE = −w_i⟨r∇_θφ_i⟩). Lines 18–21: no minimiser does not imply the flow runs off. "Its closure holds the dipoles": it also holds higher confluent terms (k coalescing kernels give (k−1)-th derivatives). | Proof (i): add the two bounds as in Lemma D‴.53a(ii). Fix line 47. Lines 18–21: "where the flow approaches that infimum it runs off". Mention the higher-order limits (de Boor). |
| D3bd-9 | minor | C's stalls are reported at the wrong state and the classification is circular. `match` returns the state before the stalled block (`rich_transit_nd.py:304-307` breaks before updating), so the md shows "2 kernels, residual 2.7e-01" and a table maximum of 2.7·10⁻¹. The real stall is the polish's run-off to |w| = 1.4–4.0·10³, which is accepted (no seed-level kernels, `rich_transit_nd.py:284-291`) and numerically arbitrary (3910 vs 3486 from a 10⁻¹⁶ change); the budget then burns without progress (994,226 evaluations in one chunk). `outcome` equates "out of budget" with the dipole (`exp_moving_width.py:119-127`), which fails on the one-sign stall. Line 75's "flowed from their entry" is not the path's trajectory. | Record the stalled block's state (z when the budget ran out) and report its |w|, separation and rms against rms_inf. Port round 54's polish status/budget to `polish_block_nd`, and record or reject a polish that multiplies max|w| by more than ~10. `outcome`: classify the dipole by an opposite-sign pair with |w| ≫ max|A| and small d. G's text: "the path's own block was polished to |w| ≈ 10³, beyond the flow's reach at s = 10¹¹ (LSODA and Radau fail just beyond it)". |
| D3bd-10 | minor | 25 of 55 mutants survive both the tests and `--from-raw` (`out/mutants.log`): all of `g_job`; `n_formed` N1, N2, N4; `outcome` O2–O6; `pairs_of` P6 and PAIR_TOL = 0.5 (P7); D report RD2, RD4; part B B2, B3, B4, B5, B7, B9. | Tests: `n_formed` (min amplitude, abs, threshold); `outcome` at its thresholds (err 2·10⁻⁶, rms 2·10⁻⁸, a pair among unformed kernels); `pairs_of` at PAIR_TOL ± ε and with the width term; `g_job` on a cheap synthetic +−+ dipole with small S_G; part B with one bad gradient point among 20; D's formed count with an idle entrant. |
| D3bd-11 | minor | Literature and novelty (lines 55–61): see the verdict above for the missing precedents, the mis-cited Duval–Peyré result and the over-broad "New". | Status: add the missing references; re-cite Duval–Peyré for the split; narrow "New" to the runaway law, the width split and the ε-dependent landing. |
| D3bd-12 | cosmetic | Misquotes and docstrings. Line 76 "1.8–5.9" → 1.8–6.0; line 75 "0.035–0.073" → 0.035–0.072; line 72 "fixed widths: 4 kernels" is the entry cap. `moving_width.py:2` and `tests/test_moving_width.py:1` "(planned)"; `moving_width.py:12-13` "converges … from every start" (false for a wrong-sign start: α cannot cross 0); `exp_moving_width.py:122` "a cancelling dipole" → "a cancelling pair"; the md header "as many kernels as features" vs `outcome`'s "formed kernels"; the comment at `exp_moving_width.py:45-46` "~2·10³" (actual 1.5–8.8·10³ in A, median 8.3·10³ in C). | Edit as stated. |

## Lesson

A classification of failure modes drawn from one target distribution describes that distribution, not the dynamics. Two positive features closer than the starting width broke "a target of one sign cannot make them" in every variant tried.

A formal asymptotic law should be closed quantitatively before it is called derived. The constant K = lim w²(E − E_inf) and the optimal dipole fit predict the prefactor; checking it exposed a factor 2 and showed κ is purely third-order.

A moments argument should be checked against the stored states. The raw contained the width split that absorbs the third moment to five digits, and it is the only thing that explains the draft's own quoted residual.

Finally, report the state where a budget actually ran out. Do not compare transients against a tolerance set ten times above the observed cases. And do not call a configuration "frozen" without a second, later time: here the polish had run the block to |w| ≈ 10³ before the budget was spent, and the "landed" pairs were still moving at four times the reported horizon.
