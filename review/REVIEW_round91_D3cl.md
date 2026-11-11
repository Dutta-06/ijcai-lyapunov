# Referee report D3cl, round 91: Proposition D‴.61 as revised after round 90

## Scope

I reviewed commit a02bbba on branch claude/vigilant-clarke-hf939b.

- **Setup.** I exported the commit with `git archive` into W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref91/repo. GP_OUT = W/gp_out and WM_OUT = W/wm_out each hold a copy of their raw. Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 15`, one heavy process of mine at a time. Once, two of mine overlapped; I paused the second with SIGSTOP on its verified PID until the first had finished.
- **Read in full:** review/drafts/D61.md, experiments/exp_gh_peano.py, results/chaos/gh_peano.md, tests/test_gh_peano.py, review/REVIEW_round90_D3ck.md, the Round 90 entry of review/REVIEW_LOG.md, and `git diff da37d04 a02bbba` for these files. I also read tools/mutants.py and the two round-90 mutant specs.
- **Environments.** Repository code ran in the repository's venv. My own code needs mpmath, which that venv lacks, so it ran in the round-90 referee's venv (ref90/venv2). I reused only that interpreter; I wrote all my code new, in W/own.
- **Repository untouched.** My only git commands were archive, diff, log, status and rev-parse. During the review the branch moved to 78a3a04 (D‴.58's round-88 report, not by me; no D‴.61 file changed), and the working tree gained uncommitted D‴.58 changes, also not mine.
- **Web.** arxiv.org, ams.org and ar5iv cannot be fetched from here, and neither can dlmf.nist.gov (it is outside the session's allowed domains).

## Summary

**Design.** The revision keeps (i)–(v) and adds the following:
- frames for every ρ ∈ (−1, 1), with the ρ = 0 limit (`frame_series_any`), and a census of 22 of A6's 36 rules;
- the first-term and Fock lower bounds;
- (iv) proved conditionally on (iii) for complex η;
- a Parseval-frame proof of (v), and P_j via h_{2j};
- part LC (the transform off the axes);
- a reduction of conjecture P to a conjecture Q on the corner ratio c_K(u) = Q_K[e^{uz}He_{K−1}]/(u^{K−1}e^{u²/2}). The route is: Duhamel's formula, Laplace's method on the segment u = ξs, then steepest descent for Q. Part Q supplies the numerics.

**Verdicts.**
- Reproducibility is exact: the md byte for byte, 16 of 16 tests, and every job of the raw recomputed bit-identically.
- Every number in the draft is printed by the md, verbatim.
- All 14 round-90 issues are addressed. D3ck-9 is only partly closed, because part Q is new and its tests leave result-changing mutants alive.
- In the reduction, these are all correct; I re-derived each and checked it numerically:
  - the Duhamel identity and the contour formula;
  - g′(s) = (1 − 2s)(1/(s(1 − s)) − Re η²);
  - the Laplace asymptotics (E − Q_K)[e^{ξz}] ~ (K!/(2K)!)ξ^{2K}e^{ξ²/4}(1 − ξ²/(4K))^{−1/2}c_K(ξ/2), with every constant;
  - the saddle at the Joukowski point, the exterior prefactor (1 − μ²)^{−1/2}, and the Gaussian factor.
- The reduction's one stated proviso is false as worded (D3cl-1).
- Apart from that, the draft undersells what the reduction yields (D3cl-2):
  - Conjecture Q follows in a few lines from an exact identity on a circle plus the classical uniform exterior asymptotics.
  - The corrected Laplace step then gives P for |η| < 1.89 at every angle. That covers everything (iv) uses at every frame of S and R; the worst frame needs |η| ≤ 1.78.
  - Only in a sector near |η| = 2 does the straight segment fail. There c_K's actual growth reaches the bound.

**Findings**
- **D3cl-1 (major).** The reduction's proviso, "provided c_K is not exponentially large away from s = ½", is false.
  - The segment u = ξs leaves Q's disc at s = 1/|η|. Beyond that point c_K is exponentially large off the imaginary axis: (1/K)log|c_K(ξs)| reaches +0.4865 at η = 1.95e^{0.8i}, s = 0.95, K = 160.
  - The condition actually needed is domination of the integrand by its value at s = ½. It holds for |η| < 1.89, and on both axes up to 2, but not everywhere.
- **D3cl-2 (major).** Q is provable, and the reduction then proves P on |η| < 1.89.
  - The key is an identity: Re f_μ(μ) − Re f_μ(μe^{iψ}) = (1 − cos ψ)(1 − |μ|²cos(2 arg μ + ψ)).
  - With it, (iv)'s sharp asymptotics at all 22 rules and all 12 S frames are proved, modulo the classical exterior asymptotics and one scalar inequality that I checked on a grid.
  - P itself stays open only in a sector of 1.89 < |η| < 2.
- **D3cl-3 (minor).** No script reproduces the reduction: neither the Duhamel identity nor the Laplace step. My checks: the Laplace-step ratio tends to 1 as 0.10/K (η = 1) and 0.074/K (η = 1.5i); near |η| = 1.9 off-axis it tends to 1 slowly; on the real axis at η = 1.9 it drifts 0.946 → 1.097 for K ≤ 160.
- **D3cl-4 (minor).** The Duhamel identity is a theorem (variation of constants), not a referee's check. It is also the classical Krylov error representation, whose a posteriori estimate is exactly the corner entry (Saad 1992). The draft cites neither.
- **D3cl-5 (minor).** Mutants. Of my 35, the tests kill 24. Four survivors change published results, all in part Q: QJ4, QJ5, QP1 and QM1. The round-90 referee's spec has a third equivalent survivor, O3, that the log does not account for. Several spec anchors are stale.
- **D3cl-6 (minor).** c_K has a Dirichlet form: c_K(u) = e^{−u²/2}E e^{uT}, with T the uniform-simplex average of the K nodes and Var T = (K − 1)/(K + 1) exactly. This is (i)'s mechanism, and the draft does not state it.
- **D3cl-7 (cosmetic).** Bookkeeping in the Laplace step: s^{−1}, 1/(K − 1)!, du = ξ ds, and the complex exponent's own stationary point G′(½) = 0.
- **D3cl-8 (cosmetic).** The circle identity explains why Q (and so P) converges slowest on the real axis. On the imaginary axis Q's limit continues beyond the disc, but not far.

## What I checked and how

**Reproducibility** [own]
- `--from-raw` regenerates gh_peano.md byte for byte (sha256 4c524b49…).
- tests/test_gh_peano.py: 16 passed in 0.83 s.
- With the script's own job functions and one worker, every job of the raw is bit-identical to the raw (`repr` equality, `secs` excluded):
  - P and M;
  - all 22 R jobs, among them the 5 with ρ < 0 and the 5 with ρ = 0;
  - all 12 S frames, 12 Q jobs, 12 LC jobs and 8 L jobs.
- Script: W/scripts/rerun.py; logs in W/logs/rerun*.log.

**Claims against the md and the raw** [md][raw]
- Every number in Numerics (P), (M), (L), (Q), (LC), (S) and (R) appears verbatim in the md. Examples: Q's 5.03·10⁻³ … 6.03·10⁻⁴ (0.488, 0.494, 0.497) and 8.81·10⁻¹ … 3.39·10⁻¹ (0.747, 0.728, 0.708); LC's 3.39·10⁻¹ … 3.73·10⁻² (0.703, 0.582, 0.269); S's 0.99786–0.99996 and 1.14137 … 1.09843; R's 0.99814–1.19578, 16 of 22, the six outliers and 1.83295.
- The census is right: 5 rules with ρ < 0, 5 with ρ = 0 (τ/σ = 2's Fock and Mercer-zero rules are the same rule, as the draft says), 12 with ρ > 0, and 14 with no frame.
- The six rules outside 0.01 all have their mass on the real axis, and all fall with K.

**Independent mathematics and code** [own, W/own]
- **Duhamel identity and contour formula, re-derived.**
  - Split J = (J_K ⊕ J′) + √K(e_{K−1}e_Kᵀ + e_Ke_{K−1}ᵀ). Use [e^{tJ}]_{0,K} = t^K e^{t²/2}/√(K!) and [e^{uJ_K}]_{K−1,0} = Q_K[e^{uz}He_{K−1}]/√((K − 1)!); since √K/√(K!(K − 1)!) = 1/(K − 1)!, this gives exactly the draft's formula.
  - From w_kHe_{K−1}(z_k) = (K − 1)!/He_K′(z_k): Q_K[e^{uz}He_{K−1}] = (K − 1)!(2πi)⁻¹∮e^{uz}/He_K. By Hermite–Genocchi on the K nodes this equals u^{K−1}E e^{uT}, which gives D3cl-6.
  - Numerically, by quadrature at K = 12, η = 1.5e^{iπ/5}, in 80 digits: 5.5·10⁻⁷⁶.
- **c_K without nodes.** I computed c_K from the Taylor series of e^{uJ_K}e₀ in mpmath, at 3K + 80 digits.
  - All 60 of Q's deviations are reproduced to ≤ 1.2·10⁻¹² relative; the residue is the repository's float conversions.
  - A second method, the node sum Σ e^{uz_k}/He_K′(z_k) with mpmath eigenvalues, agrees to 10⁻¹²²–10⁻²⁰³ off the real axis at K = 20, 40.
- **Laplace step.**
  - Derivation: (E − Q_K)[e^{ξz}] = (ξ^{2K}/(K − 1)!)∫₀¹s⁻¹e^{KG(s)}c_K(ξs)ds, with G(s) = log(s(1 − s)) + (η²/2)((1 − s)² + s²). Then G′(½) = 0 for every η, G(½) = −log 4 + η²/4 and G″(½) = −2(4 − η²). Since 4^{−K}√(πK)C(2K, K) → 1, this gives the draft's constant exactly.
  - The draft's formula against exact (E − Q_K)[e^{ξz}], node-free, at K = 20, 40, 80, 160, as |exact/prediction − 1|:

    | η | values at K = 20, 40, 80, 160 | K·deviation |
    |---|---|---|
    | 1 | 4.97·10⁻³ … 6.46·10⁻⁴ | 0.099–0.106 |
    | 1.5i | 3.57·10⁻³ … 4.60·10⁻⁴ | 0.071–0.074 |
    | 1.9e^{iπ/8} | 0.110, 0.081, 0.058, 0.035 | |
    | 1.9e^{0.45i} | 0.092, 0.067, 0.047, 0.021 | |
    | 1.95e^{0.5i} | 0.099, 0.081, 0.062, 0.017 | |
    | 1.9 (real), ratio itself | 0.946, 1.018, 1.069, 1.097 | |

  - The same run reproduces L and LC entries without nodes: on the real axis at s = 1.9, 1.976, 1.716, 1.438, 1.157 (the md's 2.9756 … 2.1569, minus 1); and LC's θ = π/8, s = 1.9 entries 3.39·10⁻¹ … 3.73·10⁻².
- **Q's steepest descent.**
  - Take f_μ(ζ) = μ(ζ + 1/ζ) + log ζ − ζ²/2 and F(μ) = μ²/2 + 1 + log μ. The saddle is ζ = μ, i.e. z = √K(μ + 1/μ), with f″ = (1 − μ²)/μ².
  - With Stirling, (K − 1)!K^{1−K}e^K ~ √(2πK), and the prefactor, c_K → 1 − μ²: the draft's account is right.
  - The identity Re f_μ(μ) − Re f_μ(μe^{iψ}) = (1 − cos ψ)(1 − |μ|²cos(2 arg μ + ψ)) holds to 2.7·10⁻¹⁵ at 10⁵ random points. It makes ζ = μ the strict maximum on the circle |ζ| = |μ|.
  - Exterior Plancherel–Rotach (PR) asymptotics, from He_n(z) = E(z + iZ)^n by steepest descent, give He_K(√Kx) ~ K^{K/2}ζ^{−K}e^{Kζ²/2}(1 − ζ²)^{−1/2}. Relative error at K = 20, 80, 320, 1280:

    | ζ | relative error | K·error at 1280 |
    |---|---|---|
    | 0.5 | 5.1·10⁻³ … 8.9·10⁻⁵ | 0.11 |
    | 0.5i | 1.2·10⁻³ … 1.8·10⁻⁵ | 0.023 |
    | 0.7e^{iπ/4} | 7.0·10⁻³ … 1.0·10⁻⁴ | 0.13 |
    | 0.9 | 0.30, 0.20, 0.11, 0.043 | 55 (pre-asymptotic) |

    With the prefactor (1 − ζ²)^{+1/2} instead, the ratio at ζ = ½ is 1.333 = 1/(1 − ζ²).
- **c_K outside the disc**, node-free, as (1/K)log|c_K(√Kμ)| at K = 40, 80, 160:
  - On the imaginary axis it is O(1)/K: c_K(i√Kt)/(1 + t²) − 1 = 4.2·10⁻³, 2.1·10⁻³, 1.0·10⁻³, 5.1·10⁻⁴ at t = 1.5 (K = 20 … 160). But at t = 2.5, c_K ≈ 10^{12}, 10^{25}, 10^{49}, 10^{100}.
  - On the real axis, at μ = 1.5: −0.146, −0.103, −0.076 (decaying).
  - Off the axes, the slopes of log|c_K| match the bound Λ⁺ defined in D3cl-1:

    | μ | slope, K 80→160 | slope, K 160→320 | Λ⁺ |
    |---|---|---|---|
    | 1.17e^{0.5i} | 0.0143 | 0.0267 | 0.0267 |
    | 1.3e^{0.8i} | 0.0843 | 0.0806 | 0.0826 |

- **The integrand along the segment**, with exact c_K, at K = 80, 160:
  - At η = 1.95e^{0.8i}, c_K's rate rises to +0.47 (K = 80) and +0.49 (K = 160) at s = 0.95, yet the total stays ≤ −0.0103 for s ≥ 0.55.
  - At η = 1.95e^{0.5i}, s = 0.6, the c_K-free part is −0.0214 (−0.0203 as K → ∞), while c_K's rate there is 0.0127 (K = 160) and 0.0197 (K = 320), heading to 0.0267. So the modulus at s = 0.6 will overtake that at s = ½.
  - On both axes at η = 1.9 the maximum stays at s = ½.
- **Grid scan of the domination inequality** sup_{s ≥ 1/|η|}[g(s) + Λ⁺(ηs)] < g(½) (W/own/laplace_gap2.py):
  - It holds at every angle for |η| < 1.8927.
  - It is violated at |η| = 1.95 for arg η ∈ (0.087, 0.611), and at |η| = 1.999 for (0.009, 0.759).
  - On both axes it holds up to 2: on the real axis by T′ = (1 − aw)²/w ≥ 0, on the imaginary axis by concavity — both by my hand derivation.
- **Radius (iv) needs.** (iv) needs P on |η| ≤ c with c² > 4 ln u/(u − 1), u = 1/r.
  - The worst frame is τ/σ = 2, 0.5·v_hs or 2·v_hs (r = 0.6345), needing c > 1.777. S's worst frame needs 1.690. The Laplace mass sits at |η| ≤ 1.593.
  - Hence every frame with r < 0.80 is covered.
- **(ii).** P_j(K) = E_X h_{2j}(X, z₁, z₁, …): exact to 7·10⁻⁵⁶ against the md's P_1–P_3 at K = 2, 3, 5, 8. I also re-derived the degree and top-term argument (p₂ = X² + 2K(K − 1)).
- **(iv) and (v), re-derived:**
  - the tail bound 4e^{−(u−1)|ξ|²/2} and the window 4/u ≤ 4 ln u/(u − 1) < c² < 4;
  - the radial integral ½Γ(2K + 1)γ(θ)^{−2K−1}, the off-cusp and cusp closed forms, and (2π)⁻¹∫|1 − re^{iθ}| = ₂F₁(−½, −½; 1; r²);
  - γx < γy ⇔ v < v_hs, and r = t(v);
  - the first-term bound a/2 and the Fock bound a/(2(1 + ρ));
  - the Parseval-frame step (H̃_n real for ρ < 0, the Taylor frame at ρ = 0) and the constant.

**Mutants** [own]
- **Method.** I ran the repository's tools/mutants.py from the export (its venv symlinked in), with pytest -x plus `--from-raw` and a byte comparison of the md.
- **Round-90 own spec.** The tests kill 28. Three survive (FA1, FA4, RR1); recomputing R under each leaves the raw identical, so they are equivalent, as the log says. Three anchors are gone: SL1 was rewritten, and RP1 and KY1 broke when part Q changed `report` and `KEYS` (my QR1 and the KEYS pin cover them).
- **Round-90 referee spec.** The tests kill 50, 14 anchors are gone (code rewritten), and 3 survive. M6 and S22 are equivalent, as logged. O3 (p_lines' zero filter; no coefficient of P_1…P_6 is 0) is equivalent too, but the log does not account for it.
- **My spec**, W/mut/spec91.py: 35 single-point mutants of part Q and of round-90 code the specs leave out. The tests kill 24.
- **Survivors of my spec, recomputed** with W/mut/recompute.py (patched private copy, the parts rerun, compared with the raw and the md):
  - Change the published md on rerun:
    - QJ4: q_job uses He_K instead of He_{K−1}.
    - QJ5: q_job takes the real part of c/(1 − μ²) − 1. At K = 20, π/4 the row becomes 1.61·10⁻⁵, 1.28·10⁻³, 8.96·10⁻³, …
    - QP1: part_Q runs LC's angles.
    - QM1: `main` runs part_LC for Q. By inspection, `--only Q` would save LC's output under Q, and q_lines would then crash on `mu`.
  - Change the raw only, not the md:
    - SJ1: the stored `check`.
    - LJ2: the stored LC `prec` at π/8 and 3π/8.
    - LJ3: the sign of the stored `arg`; the md prints |arg|.
  - Equivalent on the data: QP2 (LC_K = Q_K), FA8, FA9, RR2.
  - Cause: q_job is tested only with corner_ratio monkeypatched to a fake whose deviation is real; part_Q and main's Q wiring are never called.

**Literature** (web searches). I found no source giving sharp asymptotics of the corner entry of e^{uJ_K}, or of the Gauss–Hermite error of e^{ξz}, for ξ ~ √K. Details are in D3cl-4 and Sources.

## Verdicts per part

- **(i): proved.** Unchanged since round 90; the corrected 0 ≤ ε̃_K < 1 is right.
- **(ii): proved,** except that P_j's polynomial form is a sketch. The sketch is sound, and I checked its identity exactly [own].
- **(iii) P: numerically well supported.**
  - The reduction is correct in its exact parts.
  - Its proviso is wrong (D3cl-1).
  - Corrected and combined with D3cl-2, it proves P uniformly on |η| ≤ c for every c < 1.89, modulo the scalar inequality, which I checked on a grid only.
  - P stays open in a sector near |η| = 2.
  - The Stieltjes account is correctly marked formal.
- **Q:** the draft calls it a conjecture. By D3cl-2 it is a theorem modulo the classical uniform exterior asymptotics, with relative error O(1/K), uniformly on compact subsets of the disc (the maximum principle covers μ → 0). That rate matches part Q's 1/K at arg μ = π/2 [md].
- **(iv): correct as stated.**
  - The lower bounds are proved.
  - The asymptotics follow from (iii) on |η| ≤ c(r); I re-derived the tail bound and the choice of c, Laplace's method, r = t(v) and the branch.
  - With D3cl-2 they hold, modulo, at all 22 rules and all 12 S frames.
- **(v): proved,** for ρ ∈ (−1, 1) and at the ρ = 0 limit; the Parseval argument is right.

## Status of the round-90 issues

| id | status |
|---|---|
| D3ck-1 | Resolved: frames on ρ ∈ (−1, 1) and the ρ = 0 limit; 22 rules in R, 14 named without a frame; "every" now qualified [md] |
| D3ck-2 | Resolved: first-term bound a/2; Fock bound a/(2(1 + ρ)) better only for ρ < 0 |
| D3ck-3 | Resolved: (iii) for complex η; (iv) proved conditionally, with the tail bound and r = t(v); LC added |
| D3ck-4 | Resolved: the h_{2j} sketch (identity exact [own]); M's use of the fitted polynomials stated |
| D3ck-5 | Resolved: Parseval frame |
| D3ck-6 | Resolved: w ∉ ℝ, the saddle marked formal, "a second route" |
| D3ck-7 | Resolved: both predictions printed (16 and 12 of 22), e_K², the S sentence fixed |
| D3ck-8 | Resolved |
| D3ck-9 | Partly resolved: all of round 90's still-applicable result-changing survivors are now killed, but part Q's new code adds 4 (D3cl-5) |
| D3ck-10 | Resolved: digits agreed 261–428; S20 and S24 killed |
| D3ck-11 | Resolved, but the reduction's novelty now needs the Krylov credit (D3cl-4) |
| D3ck-12, D3ck-13, D3ck-14 | Resolved |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cl-1 | major | (iii): "Q gives P … provided c_K is not exponentially large away from s = ½". The segment μ = ηs leaves Q's disc at s = 1/\|η\|, and beyond it c_K is exponentially large off the imaginary axis. Pushing the contour to the unit circle gives the upper rate Λ⁺(μ) = M(Re μ) − ½ − log\|μ\| − Re(μ²)/2, where M(a) = a² for \|a\| ≤ 1 and 2\|a\| − 1 otherwise; the rate is attained (slopes 0.0267 against 0.0267 at 1.17e^{0.5i}, 0.0806 against 0.0826 at 1.3e^{0.8i}, to K = 320). Along η = 1.95e^{0.8i}, (1/K)log\|c_K\| reaches +0.49, yet the integrand peaks at s = ½. What is needed is domination, sup_{s ≥ 1/\|η\|}[g(s) + Λ⁺(ηs)] < g(½). It holds for \|η\| < 1.8927 and on both axes up to 2, and fails in a sector (arg η ∈ (0.087, 0.611) at \|η\| = 1.95). There the actual integrand overtakes s = ½ as K grows (η = 1.95e^{0.5i}, s = 0.6: −0.0203 + 0.0267 > 0) [own]. | Replace the proviso by the domination condition. Say that Q covers only s < 1/\|η\|, state the region where the condition is checked, and mark P in the sector as open (a deformed path, or another argument, would be needed). D61 (iii), Status. |
| D3cl-2 | major | The draft undersells Q, P and (iv). (a) Q: on the ellipse x = ζ + 1/ζ with \|ζ\| = \|μ\|, which encloses the nodes for large K, the exterior asymptotics hold uniformly, and Re f_μ(μ) − Re f_μ(μe^{iψ}) = (1 − cos ψ)(1 − \|μ\|²cos(2 arg μ + ψ)) ≥ (1 − cos ψ)(1 − \|μ\|²). Laplace's method on the circle then gives c_K(√Kμ) = (1 − μ²)(1 + O(1/K)), uniformly on compact subsets of the disc; the maximum principle covers small \|μ\|. (b) With D3cl-1's condition: Q on \|μ\| ≤ r₀, the r₀-circle bound (→ Λ⁺) beyond, and Laplace's method with complex phase (G′(½) = 0) together give P uniformly on \|η\| ≤ c for c < 1.89. (c) (iv) needs only c > (4 ln u/(u − 1))^{1/2}, at most 1.777 over R and S (r < 0.80 suffices). So (iv)'s sharp asymptotics at all 22 rules are proved modulo the classical exterior asymptotics (Plancherel–Rotach [memory]; complex-plane forms and error bounds [snippets]) and the scalar inequality, which I checked on a grid only [own]. | Recast Q as a proposition with this proof. State P as proved on \|η\| < 1.89, subject to an analytic proof of the inequality or with "checked on a grid", and conjectural only in the sector. Restate (iv) as unconditional for frames with r < 0.80, modulo the cited asymptotics. D61 (iii), (iv), Status, New. |
| D3cl-3 | minor | The repository computes nothing of the reduction; the Status points to "checked by referee round 90 … not reproduced". My node-free checks are listed under "What I checked". The Laplace-step ratio tends to 1 as 0.10/K (η = 1) and 0.074/K (1.5i), and slowly at 1.9e^{iπ/8}, 1.9e^{0.45i} and 1.95e^{0.5i} (0.035, 0.021, 0.017 at K = 160). It drifts to 1.097 at η = 1.9 real: the Laplace window (K(4 − η²))^{−1/2} ≈ 0.09 still contains the edge crossing s = 0.526 [own]. | Add a part D: the Duhamel identity by quadrature, and the Laplace step against exact (E − Q_K)[e^{ξz}], with node-free c_K (cheap) and tests. exp_gh_peano, D61 Numerics, Status. |
| D3cl-4 | minor | The Duhamel identity is variation of constants (a theorem), and it is the Krylov error representation of the K-step Lanczos approximation of e^{ξJ}e₀. Saad's 1992 a posteriori estimate h_{m+1,m}\|e_mᵀe^{H_m}e₁\| is exactly the corner entry [snippet]. The Krylov bounds (Hochbruck–Lubich, Druskin–Knizhnerman, Faber-type) are bounds, not sharp asymptotics [snippet]. | Call it proved, and cite Saad 1992 and the Krylov literature. In New, claim the Laplace analysis and Q, not the identity. D61 (iii), Status, New. |
| D3cl-5 | minor | Of my 35 mutants of part Q and of uncovered round-90 code, 4 survivors change published results: QJ4, QJ5, QP1 (md changes on rerun, verified) and QM1 (raw corrupted, by inspection). The cause: q_job is tested only with corner_ratio monkeypatched to a fake with a real deviation; part_Q and main's Q branch are untested. SJ1, LJ2 and LJ3 change the raw only; QP2, FA8, FA9 and RR2 are equivalent. Referee survivor O3 is equivalent but missing from the log; anchors SL1, RP1, KY1 and 14 of the referee spec are stale [own]. | Test q_job unpatched at K = 2 (closed form sinh(u)/(ue^{u²/2})) with a complex c; part_Q's job list; and main `--only Q` with a patched part_Q. Add a round-91 spec file and refresh the stale anchors. tests/test_gh_peano.py, tools/mutants/. |
| D3cl-6 | minor | c_K(u) = e^{−u²/2}E e^{uT}, where T = Σλ_kz_k, λ is uniform on the simplex, and Var T = (K − 1)/(K + 1) exactly. So Q reads E e^{√KμT}e^{−Kμ²/2} → 1 − μ², parallel to P for S, and c_K > 0 on ℝ. It is consistent at second order: K(Var T − 1)/2 → −1, the μ² coefficient of log(1 − μ²) [own]. | State it beside the contour formula. D61 (iii). |
| D3cl-7 | cosmetic | "the integrand's modulus without c_K is e^{Kg(s)} up to \|ξ\|^{2K−1}". The s-integrand is (ξ^{2K}/(K − 1)!)s⁻¹e^{KG(s)}c_K(ξs). Laplace's method with complex phase needs G′(½) = 0, which holds, not only that Re G is maximal. "e^{ξJ} − e^{ξJ_K}" mixes the two matrix sizes. | One line each. D61 (iii). |
| D3cl-8 | cosmetic | The circle identity's curvature is ψ²(1 − r²)/2 on the real axis against ψ²(1 + r²)/2 on the imaginary axis, and the PR error is large near the edge (K·error ≈ 55 at ζ = 0.9). That is why Q (and P) converge slowest on the real axis. On the imaginary axis c_K → 1 + t² beyond the disc (t = 1.5: 5.1·10⁻⁴ at K = 160), but at t = 2.5, c_K ≈ 10^{100} at K = 160 [own]. | Optional remark. D61 Status. |

## Sources

Web searches only; everything below is [snippet] unless marked otherwise.
- [Error bounds for the Krylov subspace methods for computations of matrix exponentials (arXiv:1603.07358)](https://arxiv.org/pdf/1603.07358): bounds the corner-type quantity e_kᵀe^{itT_k}e₁ with Faber polynomials; bounds only.
- [KIOPS (arXiv:1804.05126)](https://arxiv.org/pdf/1804.05126): attributes the h_{m+1,m}·e_mᵀφ(τH_m)e₁ error estimate to Saad 1992, Thm 5.1.
- [Simoncini, slides (Berlin 2006)](https://www.dm.unibo.it/~simoncin/berlin06.pdf): Hochbruck–Lubich's 10e^{−m²/(5ρ)} bound and the corner estimate.
- [Computable upper error bounds for Krylov approximations to matrix exponentials (PMC7039864)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7039864/)
- [Druskin–Knizhnerman 1991, MathNet record](https://www.mathnet.ru/eng/zvmmf3050): record only.
- [Low-rank updates of matrix functions (arXiv:1707.03045)](https://arxiv.org/pdf/1707.03045)
- [Wang–Wong, Asymptotics of orthogonal polynomials via recurrence relations (arXiv:1101.4371)](https://www.arxiv-vanity.com/papers/1101.4371/): monic Hermite asymptotics in the complex plane away from the oscillatory interval.
- [Shi–Nemes–Wang–Wong, Error bounds for the asymptotic expansions of the Hermite polynomials (arXiv:2109.01054)](https://real.mtak.hu/191450/1/2109.01054v2.pdf): outer-interval error bounds.
- [Dominici, Asymptotic analysis of the Hermite polynomials (arXiv:math/0601078)](https://arxiv.org/pdf/math.CA/0601078)
- [Gil–Segura–Temme, nodes and weights of Gauss–Hermite (arXiv:1709.09656)](https://arxiv.org/pdf/1709.09656)
- [Convergence analysis of Hermite approximations for analytic functions (arXiv:2312.07940)](https://arxiv.org/pdf/2312.07940)
- [On asymptotic properties of certain B-splines … theta-like functions (arXiv:2601.06643)](https://arxiv.org/html/2601.06643v1): Mellin transforms of B-splines with Hermite-zero knots; tangential.

From memory [memory], not seen here: Saad, SIAM J. Numer. Anal. 29 (1992); Hochbruck–Lubich (1997); Plancherel–Rotach (1929); Szegő, Orthogonal Polynomials, Thm 8.22.9; Deift–Kriecherbauer–McLaughlin–Venakides–Zhou (1999); Paulsen–Raghupathi, Thm 2.10.

## Lesson

1. **Test a proviso along the whole path it constrains.** The segment leaves the disc where Q lives. There "not exponentially large" is false, while the weaker domination it stands for holds for |η| < 1.89 and fails in a sector. That difference decides whether the reduction proves (iv) at every rule or proves nothing.
2. **Look for the global estimate before calling a steepest descent a conjecture.** One exact identity on the contour turned "the saddle's rigorous treatment" into a short proof.
3. **A test that monkeypatches the core function leaves the wrapper's arithmetic untested.** Three of the four result-changing survivors sit in q_job and part_Q, around a patched corner_ratio.
4. **For my own checks: compare raw parts element by element.** An ndarray against a list first gave me a spurious "raw changed".

Files: my code is in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref91/own; the reruns are in …/ref91/scripts/rerun.py; the mutant spec, recompute script and results (spec91.py, recompute.py, res_*.jsonl) are in …/ref91/mut; the logs are in …/ref91/logs.
