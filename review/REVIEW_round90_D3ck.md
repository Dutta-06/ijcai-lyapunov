# Referee report D3ck, round 90: Proposition D‴.61 as drafted

## Scope

I reviewed commit da37d04 on branch claude/vigilant-clarke-hf939b, exported with `git archive` into W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref90. Everything ran there with GP_OUT and WM_OUT pointing into W, under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 15`, one heavy process at a time.

- **Read in full:** review/drafts/D61.md, experiments/exp_gh_peano.py, results/chaos/gh_peano.md, tests/test_gh_peano.py.
- **Read in part:** D58's (i), (ii), the Status's A5/A6 paragraph and its Fock sketch with t(v); and `gh_nodes`, `he_all`, `frame_series`, `gen_ratio`, `a6_error`, `a6_scalings`.
- **Environments:** repository code ran in the repository's venv. My independent code ran in a private venv with mpmath (W/venv2), because the repository's venv has no mpmath.
- **Repository untouched:** `git status` is clean. During the review the branch moved to 6f030c6 (commits 1bafe68 and 6f030c6, by someone else). That change adds `test_the_md_is_the_report_of_the_raw` to tests/test_gh_peano.py. This report concerns da37d04.

## Summary

**Design.** The draft has five parts:
- (i) the Gauss–Hermite (GH) Peano kernel of order 2K is K!/(2K)! times the law g_K of a Dirichlet(1, …, 1) average S of (X, z₁, z₁, …, z_K, z_K), with X ~ He_K²φ/K!;
- (ii) the defects D(m) as weighted Dyck paths reaching height K, P_1 = K² + K + 1 and the exact variance;
- (iii) conjecture P: ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2};
- (iv) a Fock-space lower bound, and the Laplace-method consequences of (iii) for D‴.58's frame series;
- (v) each scaled rule's own Mehler frame and its constant.

**Verdicts.**
- Reproducibility is exact.
- Every number in the draft matches the md and the raw verbatim.
- The mathematics of (i), (ii), (v) and the Laplace closed forms is right; I re-derived it.
- Three things are wrong in what is claimed:
  - the reach of (iv)/(v): "every scaled rule";
  - the value of the one unconditional result of (iv);
  - the logical status of (iv), which the draft undersells.
- The tests let 12 result-changing mutants through.

**Findings**
- **D3ck-1 (major).** "Every scaled rule falls at t(v)" and "the sharp asymptotics of every scaled rule" overreach.
  - The frames exist only for v ∈ (τ²/σ*², τ²), i.e. ρ ∈ (−1, 1). As stated (0 < ρ < 1) that is 12 of A6's 36 rules; with ρ < 0 and the ρ = 0 limit, 22.
  - The other 14 have no frame, among them moment matching (ρ = 1, D‴.58's rate λ).
- **D3ck-2 (major).** The "unconditional lower bound" lim inf S_K^{1/(2K)} ≥ 1/(2(γx + ½)) equals ρv′/(2(1 + ρ)).
  - That is strictly weaker, for every 0 < ρ < 1, than the trivial bound ρv′/2 from the series' own first term (D_{2K}(v′) = v′^K K!).
  - Example: at the half-squeeze rule it gives t/(1 + t) where D‴.58 proved t.
- **D3ck-3 (major).** (iv) needs no extra uniformity, but it does need (iii) for complex η.
  - The tails |ξ| ≥ c√K are negligible unconditionally, since |ε̃_K(ξ)| ≤ 1 + e^{(Im ξ)²/2}.
  - Conjecture (iii), read for complex η and locally uniformly on |η| < 2, then implies (iv)'s closed forms.
  - The draft states η ambiguously, gives evidence on the two axes only, and calls (iv) "a sketch with a uniformity it does not state".
  - The proof also never shows that r = 1/(2 min γ) equals t(v).
- **D3ck-4 (minor).** (ii)'s "polynomial of degree 2j, leading coefficient 1/j!" is provable in a few lines: P_j(K) = E_X h_{2j}(X, z₁, z₁, …, z_K, z_K), with h the complete homogeneous symmetric polynomial. I verified the identity exactly. M's "exact" cumulants rest on the interpolated polynomials; my direct path counts confirm them exactly up to K = 5120.
- **D3ck-5 (minor).** (v)'s proof does not justify integrating Mehler's series term by term against N(0, τ²). Crude domination fails for ρ ≥ σ²/τ², e.g. three of the τ/σ = 2 frames. A Parseval-frame argument fixes it for all ρ ∈ (−1, 1) and the ρ = 0 limit.
- **D3ck-6 (minor).** The Stieltjes identity holds only for w ∉ ℝ. Its Laplace representation converges only for Re w = 0, so the saddle point "as w reaches 2√K" is formal. "Within reach" is speculative.
- **D3ck-7 (minor).** "Within 0.01 at 11 of 12 rules at K = 160" is against the per-angle numerical prediction. Against (iv)'s closed forms, which the draft calls "the sharp asymptotics", it is 10 of 12, and it concerns e_K², not e_K. In S, "the closed forms over the per-angle prediction … up to 1.10577" is mislabelled.
- **D3ck-8 (minor).** "Pointwise convergence" on the real axis is not shown at K ≤ 160. My run to K = 1280 supports it.
- **D3ck-9 (minor).** Of 67 single-point mutants, 22 survive the tests and 12 of those would change the published md on rerun. They include the cusp's ₂F₁ factor and a `--from-raw` that writes nothing.
- **D3ck-10 (minor).** S's "precision check" column is vacuous (two double-precision logs compared). R's "to 7.1·10⁻³⁹" is rounding at 80 digits.
- **D3ck-11 (minor).** Novelty: (i) and (ii) are packagings of classical facts (Gauss error formula, Curry–Schoenberg/Hermite–Genocchi, B-splines as Dirichlet averages, Flajolet–Viennot paths), and Flajolet–Viennot is not cited. Conjecture P and the exact variance were not found in the literature I could reach.
- **D3ck-12 (cosmetic).** "0 < ε̃_K ≤ 1 on ℝ" is false at ξ = 0, and the upper bound is strict.
- **D3ck-13 (cosmetic).** "In cumulants" is an implication, not an equivalent, and c_1 is proved by (ii).
- **D3ck-14 (cosmetic).** P_j means different things in D‴.58 and D‴.61. Frames S's v′ = 1.6 and 2.5 at ρ = t are synthetic, not any rule's own frame. (v) does not say which v gives ρ < 0.

## What I checked and how

**Reproducibility** [own]
- `--from-raw` regenerates gh_peano.md byte for byte (sha256 eab5800d…).
- tests/test_gh_peano.py: 10 passed in 0.97 s.
- The following reruns in a fresh GP_OUT, with one worker, are bit-identical to the raw (`repr` equality):
  - P and M fully, in under 1 s;
  - L jobs (160, real), 32 s, and (20, imag);
  - S frames 0 (τ/σ = ½ cusp), 6 (1, v′ = 2.5) and 10 (2, v′ = 2.5);
  - all 12 R jobs, 0.1 s.

**Claims against the md and the raw** [md][raw]
- Every number in (P), (M), (L), (S) and (R) appears verbatim. Examples: 0.4805194805, −0.249878, the 2s(2K) − s(K) values, 2.9756/2.1569, 0.99786–0.99996, 0.99814–1.05221, v′ = 6.910477, ρ = 0.102458, 7.1·10⁻³⁹.

**Independent code** [own] (in W/own; none of it uses the repository's functions)
- **Path counts (cumulants).** A band-pruned weighted-path DP gives D(2K + 2j) directly for j ≤ 6 at every K of M_K (1.3 s at K = 5120). Cumulants come from a power-series logarithm.
  - var, scaled and nxt are bit-identical to the raw at all 10 K.
  - The fitted polynomials P_j equal the direct counts at every K up to 5120.
- **ĝ_K without nodes.** ĝ_K(ξ) = ((2K)!/K!)Σ_j D(2K + 2j)ξ^{2j}/(2K + 2j)!, from exact capped-DP integers (mpmath at up to 2.5× the cancellation in digits).
  - All 8 raw L jobs are reproduced bit for bit.
  - Extended:
    - K = 320: real axis at s = 1.9 gives 1.88943; deviation on s ≤ 1.5 is 2.95·10⁻²; imaginary axis 4.66·10⁻⁴.
    - K = 640, 1280, real axis: s = 1.9 gives 1.65243, 1.45668; s = 1.5 gives 1.01578, 1.00822 (deviation ratios 0.557, 0.535, 0.520, tending to 1/K).
    - Off the axes, at η = s·e^{iθ} with θ = π/8, π/4, 3π/8 and s ≤ 1.8, |ratio − 1| falls about as 1/K. At θ = π/8, s = 1.8: 0.261, 0.168, 0.087, 0.030, 0.0122 for K = 20 … 320.
- **S_K exactly.** D_n(v′) = Σ_m c_{n,m}v′^{m/2}D(m) in integer arithmetic, without nodes. log S agrees with the raw at all 12 frames for K ≤ 40, to the 15 digits printed.
  - My quadrature in the angle (scipy `quad`) reproduces the raw `pred` and `simple` to ≤ 5.7·10⁻¹⁴ in the log.
  - The full 2-D conjecture-P integral, without the radial Laplace step, over `pred` is 1.00004–1.00873 at K = 10 and 1.00000–1.00062 at K = 80. So S's deviations are conjecture P's own finite-K error.
- **Fock integral.** A 2-D grid integral of the true |ε̃_2|²e^{−Ax²−By²} equals e_K²/C at ρ = ¼, v′ = 3 and at ρ = −½, v′ = −¾ (ratio 1.000000).
- **Frame constant (v).** I wrote my own closed-form Gaussian L² error and my own series (120 digits, K = 2, 3, 5).
  - It holds at ρ = ¼ (deviation 2–7·10⁻⁴⁷).
  - It holds at ρ = −0.46 … −0.94 (deviation ≤ 7·10⁻⁴⁵, the series tolerance) for: τ/σ = 2 with v = 1 and with 0.5v_hs; τ/σ = 1 with Mercer-zero and with 0.8v_hs; Mercer-zero at τ/σ = ½ and ¼.
  - It holds at the Fock scaling (ρ = 0) via the limit S_K = Σ_n(v/2)ⁿD(n)²/n! with C = σ*v/τ², to 10⁻⁶¹–10⁻⁶⁹.
- **(i)'s law by Monte Carlo.** Sample X by rejection and λ ~ Dirichlet (4·10⁶ samples; 1.6·10⁷ for a second seed at K = 5).
  - E S², E S⁴, E S⁶ against the path counts, and E e^{ξS} at ξ = 1, 2, 1 + i, 3i against the value from the nodes: K = 2, 3 within 2.3σ; K = 5 within 1.7σ on rerun (the first seed's 3.5σ was one shared fluctuation, mean S 2.6σ off 0).
  - The Stieltjes identity at w = 1 + 2i, 0.5 + 3i, −2 + 1.5i holds within 1.6σ, with h_K by `quad`. The classical (E − Q_K)[1/(w − z)] = h_K/He_K holds to 10⁻¹⁶.
- **Algebra re-derived by hand:**
  - weight factorisation e^{−Ax²−By²} = |e^{−(A−B)ξ²/4}|²e^{−(A+B)|ξ|²/2}, α = 1/(ρv′), kernel constant (α/π) = (πρv′)⁻¹;
  - x² = 2K/(A + 1) and the Stirling limit 1/(2(A + 1));
  - the off-cusp closed form r^{2K}(1 ∓ r)/(ρv′(2γ(γ′ − γ))^{1/2}) and the cusp form (πK)^{1/2}r^{2K}₂F₁(−½, −½; 1; r²)·2r/(ρv′);
  - the branch assignment: γx < γy ⇔ v < v_hs ⇔ rate attained at τ²/v, mass on the real axis, factor 1 − r;
  - 1/(2γx) = h(τ²/v) and 1/(2γy) = h(vσ*²/τ²), so r = t(v);
  - A = v/(τ² − v) and B = vσ*²/(vσ*² − τ²), smooth through ρ = 0;
  - the constant of (v); the P_1 path count; the saddle τ² − 2ω²τ + 4ω² = 0, coalescing at ω = 2, η = 2.
- **Further identities, checked numerically:**
  - P_j(K) = E_X h_{2j}(X, doubled nodes): exact for j ≤ 6 at K ∈ {1, 2, 3, 5, 8, 13}.
  - Duhamel's (E − Q_K)[e^{ξz}] = (1/(K−1)!)∫₀^ξ(ξ − s)^K e^{(ξ−s)²/2}Q_K[He_{K−1}(z)e^{sz}]ds: 10⁻¹²–10⁻¹⁶ for K = 1, 3, 5, 8.
- **Mutants.** 67 single-point mutants in a private copy (W/mut). Each is run against pytest -x and against `--from-raw` (md byte comparison). Each also has the cheap parts recomputed (P, M, L at (20, ·) and (40, imag), S frames 0/6/9, all of R) and compared with the raw. Results are in W/mut/res.jsonl.
- **Literature.** Web searches; arxiv.org, ams.org and ar5iv cannot be fetched from here, and neither could viennot.org.

## Verdicts per part

- **(i): proved and correct** [own derivation, Monte Carlo].
  - The conditions suffice: E|f^{(2K)}(S)| < ∞ forces E|f(Z)| < ∞ through f = Hf + He_K²f[·].
  - Complex f is handled by real and imaginary parts.
  - ĝ_K is entire and ĝ_K ≥ 1 on ℝ.
  - Wrong only at ξ = 0 (D3ck-12).
  - The ingredients are classical; the novelty is the packaging (D3ck-11).
- **(ii): partly proved.**
  - Proved: the path sums, integer-valuedness (each path reaching K carries K!) and P_1 = K² + K + 1 (the 2K − 1 paths, weight K²K!).
  - Checked only: polynomial form, degree 2j and leading coefficient 1/j!. Provable (D3ck-4).
  - Var S = (K² + K + 1)/((K + 1)(2K + 1)) is proved.
- **(iii): a conjecture, well supported numerically**, on both axes, off the axes, and on the real axis to K = 1280 [own].
  - The Stieltjes identity is proved for w ∉ ℝ [own MC].
  - The saddle-point account is formal (D3ck-6).
  - "Locally uniformly" needs to say complex η (D3ck-3).
- **(iv): the lower bound is correct but dominated** (D3ck-2).
  - The closed forms are correct as consequences of (iii) for complex η.
  - That implication is provable along the lines in D3ck-3; the draft calls it a sketch.
  - The scope is overstated (D3ck-1).
- **(v): correct**, at ρ > 0, ρ < 0 and the ρ = 0 limit [own]. The proof's interchange step needs one line (D3ck-5).

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ck-1 | major | (iv) "so every scaled rule falls at … t(v)", and New's "the sharp asymptotics of every scaled rule". A frame needs ρ ∈ (−1, 1), i.e. v ∈ (τ²/σ*², τ²). Of A6's 36 rules: 12 have 0 < ρ < 1 (R's set), 5 have −1 < ρ < 0, 5 have ρ = 0, 9 have ρ ≥ 1 (including moment matching), 5 have ρ ≤ −1 [own, raw A6]. The Fock integral holds verbatim for ρ < 0, with ρv′ = v(1 − ρ²)/2 > 0, A, B > 0 [own grid check], and at ρ = 0 as a limit [own, 10⁻⁶¹]. On the 10 rules with ρ ≤ 0, conjecture P's prediction gives e_K²/(C·pred) at K = 160 of: 1.00022, 1.00051 (τ/σ = ¼); 1.00108, 1.00225 (½); 1.00637, 1.01167, 1.01261 (1); 1.06510 twice, 1.19578 (2). All fall with K, all on the real branch [own]. "(R) D‴.58's 12 rules with a frame" contradicts (v)'s own range of ρ ≠ 0 in (−1, 1), which gives 17. | State (iv)–(v) for v ∈ (τ²/σ*², τ²) with ρ ∈ (−1, 1) and the ρ = 0 limit; drop "every". Extend `rule_frame`, `r_rules` and `frame_series` (√v′ < 0) to ρ ≤ 0 and add those rules to R. Note that moment matching and the rules with \|ρ\| ≥ 1 are outside reach. D61 (iv), (v), New; exp_gh_peano R. |
| D3ck-2 | major | The Fock bound is 1/(2(A + 1)) = ρv′/(2(1 + ρ)). The series' first term gives S_K ≥ ρ^{2K}v′^{2K}K!²/(2K)!, i.e. rate ≥ ρv′/2, larger for every 0 < ρ < 1. Half-squeeze, τ/σ = ¼, ½, 1, 2: 0.0149, 0.0528, 0.1464, 0.2764 against t = 0.0152, 0.0557, 0.1716, 0.3820, which D‴.58 proved [own]. Only for ρ < 0 does it beat the first term (τ/σ = 2, 0.5v_hs: 0.388 against 0.102, rate 0.635). | Replace it by the first-term bound, or state the Fock bound only for ρ < 0 where it is informative. Remove "an unconditional lower bound" from Status and New unless it improves on the trivial bound. D61 (iv), Status. |
| D3ck-3 | major | (iv) proof: "the tails beyond \|η\| ≥ 2 … what (iii) must be strengthened to". By (i), 0 ≤ ε̃_K(x) on ℝ and the weights are positive, so \|ε̃_K(ξ)\| ≤ 1 + e^{y²/2}. Hence \|ε̃_K\|²e^{−Ax²−By²} ≤ 4e^{−(u−1)\|ξ\|²/2}, with u = 2 min γ = 1/r > 1. The region \|ξ\| ≥ c√K is o(r^{2K}) once c² > 4 ln u/(u − 1), which is < 4. The Laplace mass sits at \|ξ\|² ≤ 4Kr < c²K, since ln u/(u − 1) ≥ 1/u. So (iii) for complex η, locally uniformly on \|η\| < 2, gives (iv)'s closed forms rigorously; no further uniformity is needed [own derivation]. The Fock integral needs every angle, but (iii) does not say η is complex, and L checks two axes only (my off-axis run supports it). The proof never shows r = t(v): 1/(2γx) = h(τ²/v), 1/(2γy) = h(vσ*²/τ²) [own]. | Restate (iii) for complex η on compact subsets of \|η\| < 2. Restate (iv) as proved conditional on (iii), with the tail bound and the two h-identities in the proof. Add an off-axis L-row, e.g. θ = π/8, π/4, 3π/8. D61 (iii), (iv), Status; exp_gh_peano L. |
| D3ck-4 | minor | (ii) states "polynomial of degree 2j … leading coefficient 1/j!" while the proof says "checked …, not proved". A proof follows from f[x₀, …, x_{2K}] = h_{2j}(x₀, …) for f = z^{2K+2j}: P_j(K) = E_X h_{2j}(X, z₁, z₁, …). Here h of the doubled nodes is a polynomial in the power sums. These are polynomials in K by Newton's identities, since e_{2i} = (−1)^i(2i)!/(i!2^i)·C(K, 2i). E_X X^a = E[He_K²Z^a]/K! is a polynomial by the three-term recurrence. The degree is ≤ 2j because odd power sums vanish. The top term is p₂^j/(2^j j!) with p₂ = 2K(K − 1), giving 1/j!. The same identity gives P_1 = (2K + 1) + K(K − 1). The identity is exact for j ≤ 6 at K ∈ {1, 2, 3, 5, 8, 13} [own]. M evaluates the interpolated polynomials at up to K = 5120, 250 times the check range; direct path counts agree exactly at all 10 K [own]. | Prove polynomiality via h_{2j}, or mark it checked in (ii)'s statement. Compute M from direct counts (my band DP takes 1.3 s at K = 5120), or say that M uses the fitted polynomials. D61 (ii), Status; exp_gh_peano part_M. |
| D3ck-5 | minor | (v) proof: the Mehler sum is integrated against N(0, τ²) term by term without justification. The absolute-series bound is e^{−(1−2ρ)(y²+y′²)/(4σ²)+\|yy′\|/(2σ²)}, integrable only for ρ < σ²/τ². It fails, for instance, at τ/σ = 2 for ρ = 0.38, 0.5, 0.94 [own]. D‴.58's own proof of (i) addressed this. | Mehler's formula is k = Σe_n⊗e_n pointwise, so {e_n} is a Parseval frame of the RKHS (Paulsen–Raghupathi, Thm 2.10 [memory]). Then ‖∫k dμ‖² = Σ\|∫e_n dμ\|² for every finite signed μ, for ρ ∈ (−1, 1) (real e_n for ρ < 0) and for the ρ = 0 Taylor frame. D61 (v) proof. |
| D3ck-6 | minor | (iii) and Status. (a) S has full support, so E(1 − S/w)^{−(2K+1)} and h_K exist only for w ∉ ℝ; the cut of h_K is all of ℝ. "Plancherel–Rotach … of h_K outside [−2√K, 2√K]" must mean off the axis. (b) On ℝ, ĝ_K(x) ~ ((2K)!/K!)x^{−2K}e^{x²/2}, because ε̃_K → 1. So E(1 − S/w)^{−(2K+1)} = ∫t^{2K}e^{−t}ĝ_K(t/w)dt/(2K)! converges only for Re w = 0, where it sees ĝ_K on the imaginary axis only. The saddle point "as w reaches 2√K" (real w) is formal [own]. (c) "Within reach" is therefore speculative; the hard part is the real axis near η = 2. | Say w ∉ ℝ, mark the saddle-point account formal, and replace "looks within reach" by "a possible route". Optionally add the exact Duhamel representation above (checked to 10⁻¹²–10⁻¹⁶ [own]) as a second route. D61 (iii), Status. |
| D3ck-7 | minor | (R) and New: "within 0.01 … at 11 of 12" uses `log_s_pred`, the per-angle quadrature, not (iv)'s closed forms. With the closed forms at K = 160: 0.99899–1.06064, 10 of 12; outside are τ/σ = ¼ Mercer frame (1.02997) and 2, 0.8v_hs (1.06064) [own]. The ratios are of e_K²; New says "of the exact errors". S: "the closed forms over the per-angle prediction near the cusp up to 1.10577". The 1.10577 is S over the closed form [md column]; pred/closed is 1.1038, and the column's maximum is 1.17622 (τ/σ = 2, v′ = 2.5). | Name the prediction used, add an e_K²/(C·closed) column, say e_K², and fix the S sentence. D61 Numerics (S), (R), New; exp_gh_peano r_lines. |
| D3ck-8 | minor | L: "pointwise convergence, slower toward s = 2". At s = 1.9 the successive deviation ratios are 0.869, 0.838, 0.804 [md], not yet a rate. My extension continues to fall: 1.889, 1.652, 1.457 at K = 320, 640, 1280. On s ≤ 1.5 the ratios are 0.557, 0.535, 0.520, tending to 1/K [own]. S's "(as 1/K)" for the imaginary branch fails at τ/σ = ½, v′ = 1.6 (−0.046, −0.070) [md]. | "Decreasing at every s, consistent with pointwise convergence". Optionally add K = 320–1280 real-axis rows (the node-free series is cheap). D61 (L), (S). |
| D3ck-9 | minor | Mutants: 67, killed by the tests 45, survivors 22.<br>• No change to raw or md: M6 (algebraically identical), S20, S22, S24, R12, O3.<br>• Raw changes below the printed digits: S15 (4·10⁻¹⁰), S18 and S19 (stored `prec` only), S21 (S_THETA = 64: 1.1·10⁻⁸ in log pred).<br>• Would change the published md on rerun (12): P7 (P_KCHK = 14), M7 (M_K end 4096), S12 (cusp ₂F₁ → 1: S/closed 1.00153, 1.00124, 0.99987 become 1.00231, 1.00862, 1.03668), S16 and S17 (frame values), S25 (frame filter), R6 (ρ < 0 admitted; my check stopped at the changed rule count, not run to completion), R7 (r_job at 50 digits), R13 (`dv` min), R14 (R_K = (2,)), O5 (S summary order), O6 (`--from-raw` writes nothing: test_report_and_main passes because the md from `--only P,M` is still there) [own].<br>S12 survives because the cusp test compares at t = 0.2 with tolerance 0.02, and ₂F₁(0.04) − 1 ≈ 0.01. The later test in 6f030c6 (md == report(raw)) catches report-only changes, but not the compute survivors or O6. | Add tests: the cusp branch of `log_s_simple` against hyp_half exactly; `s_frames` values; S_THETA, P_KCHK, M_K, R_K and r_job's precision pinned; `r_lines` with two rules and two K; s_lines with two axes; `--from-raw` on a deleted md. tests/test_gh_peano.py. |
| D3ck-10 | minor | S's "precision check" compares float(ln S) at two precisions and is 0.0e+00 everywhere [md]. S_CHECK = 0 and a check at the same precision are both equivalent on the data (S20, S24). R's frame-constant deviations, up to 7.1·10⁻³⁹ at 80 digits, are the rounding of a cancelling sum, not a 39-digit identity. | Report the relative difference of the two Decimals; describe the constant's check as exact up to the 80-digit rounding. exp_gh_peano s_job, s_lines; D61 (R). |
| D3ck-11 | minor | Novelty and literature. (i) follows from the Gauss error R[f] = E[ω²f[x, x₁, x₁, …]], plus divided differences as B-spline integrals [snippet: de Boor], plus B-splines as Dirichlet averages [snippet: Karlin–Micchelli–Rinott 1986 record, Massopust survey]. Not found as a stated result in my searches; an explicit packaging. (ii)'s path sums and their truncation are Flajolet 1980 / Viennot: convergents count paths of bounded height [memory; snippets confirm the framework only], uncited. The exact variance was not found [searches] and is elementary. Conjecture P was not found. Nearest: Q_K[e^{ξz}] = e₁ᵀe^{ξJ_K}e₁ is K-step Lanczos for e₁ᵀe^{ξJ}e₁, whose superlinear onset at m ~ (‖A‖τ)^{1/2} is Hochbruck–Lubich's [snippet]; Derevianko–Plonka bound GH errors for exponential-type integrands, not sharply [snippet]. Sharp constants: Karvonen–Oates–Girolami's are non-matching [snippet]; D‴.61's are conditional on (iii). | Cite Flajolet/Viennot and the B-spline–Dirichlet literature; downgrade (i)'s "new" to "explicit form not found". Keep conjecture P and the exact variance as "not found". Mark the sharp asymptotics conditional on (iii) and limited per D3ck-1. D61 Status/New. |
| D3ck-12 | cosmetic | "On ℝ, 0 < ε̃_K ≤ 1": ε̃_K(0) = 0, and ε̃_K < 1 because Q_K[e^{xz}] > 0. | "0 ≤ ε̃_K < 1, zero only at 0". D61 (i). |
| D3ck-13 | cosmetic | "In cumulants, …" reads as equivalent; it is implied by (iii) but does not imply it. K(κ₂ − ½) → −¼ is proved by (ii); the exact value is (5K + 1)K/(4(2K² + 3K + 1)) − ¼ … → +5/8 for K times the difference [own]. | "Implies, in cumulants"; mark c₁ as proved. D61 (iii). |
| D3ck-14 | cosmetic | D‴.58's P_j (v′ = 2 defects, degree j) and D‴.61's P_j (z^m defects, degree 2j) collide. S's v′ = 1.6 and 2.5 frames at ρ = t are not any rule's own frame. (v) does not say that ρ < 0 means v ∈ (τ²/σ*², Fock). | Rename one P_j, call those S frames synthetic, add the ρ < 0 range. D61, gh_peano.md S header. |

## Sources

Web searches only; arxiv.org, ams.org and ar5iv cannot be fetched from here, so everything below is [snippet].
- [de Boor, Divided differences (arXiv:math/0502036)](https://arxiv.org/pdf/math/0502036): the normalized B-spline is the Peano kernel of the divided difference.
- [On the moments and distribution of discrete Choquet integrals (arXiv:0811.0468)](https://arxiv.org/pdf/0811.0468): Hermite–Genocchi with the B-spline Peano form.
- [Karlin–Micchelli–Rinott, Multivariate splines: a probabilistic perspective (JMVA 20, 1986), record](https://ideas.repec.org/a/eee/jmvana/v20y1986i1p69-90.html) and [HUJI record](https://cris.huji.ac.il/en/publications/multivariate-splines-a-probabilistic-perspective/): existence only, content unverified.
- [Massopust, Fractional operators, Dirichlet averages, and splines (arXiv:1309.0240)](https://arxiv.org/pdf/1309.0240)
- [Convergence analysis of Hermite approximations for analytic functions (arXiv:2312.07940)](https://arxiv.org/pdf/2312.07940)
- [Derevianko–Plonka, Differential approximation of the Gaussian by short cosine sums (arXiv:2307.13587)](https://arxiv.org/pdf/2307.13587)
- [Hermite quadrature error via the Peano kernel (arXiv:2602.18634)](https://arxiv.org/pdf/2602.18634)
- [Karvonen–Särkkä, Gaussian kernel quadrature at scaled Gauss–Hermite nodes (arXiv:1803.09532)](https://arxiv.org/pdf/1803.09532)
- [Karvonen–Oates–Girolami, Integration in RKHS of Gaussian kernels (Math. Comp. 2021)](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/)
- [Simoncini, slides citing Hochbruck–Lubich's Krylov bounds for exp](https://www.dm.unibo.it/~simoncin/berlin06.pdf)
- [Multiple orthogonal polynomials … branched continued fractions (arXiv:2204.11528)](https://arxiv.org/pdf/2204.11528) and [Gessel-style path decompositions of J-fraction convergents (Waterloo)](https://uwaterloo.ca/math/sites/default/files/uploads/documents/gjjcta1986-1.pdf): Flajolet's framework.
- [MathWorld, Hermite–Gauss quadrature (error term)](https://mathworld.wolfram.com/Hermite-GaussQuadrature.html)

From memory [memory], not seen here: Flajolet, Discrete Math. 32 (1980); Viennot (1983); Carlson, B-splines, hypergeometric functions and Dirichlet averages (1991); Paulsen–Raghupathi, An Introduction to the Theory of RKHS (2016), Thm 2.10; Deift–Kriecherbauer–McLaughlin–Venakides–Zhou (1999); Brass–Petras (2011).

## Lesson

1. A proved bound is not automatically a result. Before an unconditional bound is listed, compare it with the trivial bound from the very series it bounds. Here the first term beats the Fock bound in the whole stated range.
2. A claim made "for every rule" should be swept against the whole rule set the claim names, not the subset the code admits. Here that is A6's 36 rules against R's 12. The two sweeps found both the overreach (14 rules without a frame) and an extension the draft missed (10 more rules, ρ ≤ 0, where the prediction works).
3. A conditional result should be stated with exactly the hypothesis it uses. Here, once η is read as complex, the "missing uniformity" is just what (iii) already says.

Files: my independent checks are in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref90/own; the mutant spec and results are in …/ref90/mut (spec.py, res.jsonl, run.log).
