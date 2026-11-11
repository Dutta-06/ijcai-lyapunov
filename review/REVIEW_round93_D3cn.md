# Referee report D3cn, round 93: Proposition D‴.61 as revised after round 91

## Scope

I reviewed commit 363855d on branch claude/vigilant-clarke-hf939b.

- **Setup.** The commit was exported with `git archive` into W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref93/repo.
  - GP_OUT = W/gp_out and WM_OUT = W/wm_out each hold a copy of their raw.
  - Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, with at most two heavy processes of mine at a time. Pytest ran with `--basetemp` inside W.
  - Repository code ran in the repository's venv. My own code (W/own, all new) ran in …/ref90/venv2 for mpmath.
  - For tools/mutants.py I symlinked the venv into the export.
- **Read in full:** review/drafts/D61.md, experiments/exp_gh_peano.py, results/chaos/gh_peano.md, tests/test_gh_peano.py, review/REVIEW_round91_D3cl.md, the Round 91 entry of review/REVIEW_LOG.md, tools/mutants/round91_D61_own.py, round91_D61_referee.py and tools/mutants.py. I also read `git diff a02bbba 363855d` for these files. Within that range, commit 363855d touches only D‴.58's files.
- **Repository untouched.** My only git commands were archive, diff, log, show, status and rev-parse.
  - During the review HEAD moved to 89e9e28 (D‴.60, not mine; no D‴.61 file changed).
  - The working tree has uncommitted moving_width changes, not mine.
  - Another session's mutants.py run on exp_noise_resonance was visible; I did not touch it.
- **Web.** arxiv.org does not resolve from this session (ENOTFOUND). I used search snippets only.

## Summary

**Design.** The revision applies round 91. Three conjectures become theorems:
- **(iii)(a): the corner ratio's limit, now proved.** c_K(√Kμ) = (1 − μ²)(1 + O(1/K)). The proof combines an exact representation He_K(√K(ζ + 1/ζ)) = K^{K/2}ζ^{−K}e^{Kζ²/2}I_K(ζ) with a Gaussian bound, the circle identity, Laplace's method and the maximum principle.
- **(iii)(b): the transform's limit on the open set Ω,** defined by a domination condition with the bound Λ⁺. It is proved analytically on both axes and grid-checked on |η| < 1.8927. A sector near |η| = 2 stays conjectural.
- **(iii)(c): the cumulants.**
- **(iv): the sharp asymptotics for every frame with r < r* = 0.8053.** That covers all 22 of A6's rules with a frame and the half-squeeze rule for σ*/σ < 9.27.
- **New material:** the Dirichlet form of c_K, Duhamel's identity credited to Krylov/Saad, and a part D. Part D covers quadrature of the identity, the Laplace step, the domination grid, the radius each frame needs, and c_K beyond the disc on the imaginary axis with a heuristic "alias rate".

**Verdicts.**
- **Reproducibility is exact.** The md regenerates byte for byte. All 22 tests pass. Every job of the raw is recomputed bit-identically, including all of part D.
- **Every number in the draft is in the md, verbatim.**
- **The mathematics of (a), (b) on the axes, (c) and (iv) is correct.** I re-derived every step. No major issue.
- **The domination condition holds on the whole disc |η| ≤ 1.8927.** It reduces to a closed form whose supremum over s is explicit, and a floating-point branch and bound confirms it. The grid check can therefore become a computer-assisted proof.
- **The new heuristic alias rate is wrong**, and so is its threshold 1.8637.
- **"Λ⁺ is attained" holds only where |Re σ| > 1.**
- **All 8 round-91 issues are resolved.**
- **Mutants:** of my 32 mutants of round 91's new code, 27 are killed. No survivor changes the md.

**Findings**
- **D3cn-1 (minor).** The alias rate (11/18)y² − log y − 3/2 evaluates the residue sum's modulus at a real stationary-phase point. But the modulus e^{−Kx²/4} varies on the same exponential scale, so the true saddle is complex.
  - The first alias is the second saddle t = i(1/ζ − ζ) of I_K's integrand. Its exponent is f₁(ζ) = σ(ζ + 1/ζ) + 3 log ζ − ζ² + ζ^{−2}/2, with f₁′ ∝ (ζ² − 1)(2ζ² − σζ − 1). The saddles (iy ± √(8 − y²))/4 lie on |ζ| = 1/√2.
  - The resulting rate is λ₁(y) = (5/8)y² − log y − ½ − (3/2)log 2, which vanishes at y = 1.8589, not 1.8637.
  - Exact c_K at 13 consecutive K, fitted by Prony's method, give λ₁ to 5 digits and the oscillation frequency to 5 digits. The two-point slopes the draft cites oscillate and cannot discriminate [own].
- **D3cn-2 (minor).** "Its rate is Λ⁺, attained there (referee round 91)" holds only where |Re σ| > 1, i.e. at the edge saddles ζ = ±1.
  - For |Re σ| < 1, Λ⁺ is the circle's maximum, not a saddle, and the growth is the first alias's. Examples: 0.08064 against Λ⁺ = 0.08264 at 1.3e^{0.8i}; 0.80361 against 0.80685 at 2e^{1.2i}; 0.267 against 0.807 at 2i.
  - In the violated sector at |η| = 1.95 the maximiser has Re(ηs*) ∈ [1.009, 1.099], so the failure there is genuine, as claimed [own].
- **D3cn-3 (minor).** The domination condition has a closed form, and the grid can be replaced by a computer-assisted proof.
  - The closed form: gap = log w + (p/2)(w² − ½) + M(R(1 − w)cos φ) − ½ − log R + log 4, with w = 1 − s ∈ (0, 1 − 1/R], η = Re^{iφ} and p = Re η². Its supremum over w is explicit: endpoints, the kink, and the roots of two quadratics.
  - A floating-point branch and bound confirms the condition on 1 < |η| ≤ 1.8927 (816,881 boxes) and on |η| ≤ 1.78 (7,631 boxes).
  - The quoted margins −9.14·10⁻³, −3.29·10⁻³, −9.71·10⁻⁴ are attained on the real axis. There the gap is proved, and for R ≤ 1.85 it decreases with arg η. So the margins say nothing about off-axis robustness [own].
- **D3cn-4 (minor).** Claims overreach in places and undersell in others.
  - "New" says (iv) "proves both of D‴.58's conjectures where r < 0.8053" without the grid proviso. D‴.58 (i)'s statement does the same; its Status has the proviso.
  - Conversely, the Status undersells: (b) on |η| < 1 and (c) need no grid at all.
- **D3cn-5 (cosmetic).** The proof of (b) leaves the band r < |ηs| < 1 (s just below 1/|η|) to a continuity argument it does not state. It also omits why Ω is open with a uniform gap on compacts. One sentence each closes both.
- **D3cn-6 (cosmetic).** Several rounded values are slightly optimistic.
  - On a fine grid with the exact supremum, the least radius is 1.8927327 at arg η = 0.3432; the grid's value 1.8927345 is 1.8·10⁻⁶ high.
  - Hence r* = 0.805257 and (1 + r*)/(1 − r*) = 9.2699, so "r* = 0.8053" and "σ*/σ < 9.27" round up.
  - The violated sectors are [0.0864, 0.6173] at |η| = 1.95 and [0.0016, 0.7662] at 1.999, slightly wider than the grid's [0.087, 0.615] and [0.004, 0.764].
- **D3cn-7 (cosmetic).** Small statements:
  - The ellipse encloses the nodes for every K, not only large K, and I_K is zero-free on 0 < |ζ| < 1.
  - "(so c_K > 0 on ℝ)" follows from the Dirichlet form, not from Var T.
  - (iv)/(v) write v ∈ (τ²/σ*², τ²) (σ = 1) beside formulas in σ.
  - The md and the script still call Q "conjecture Q".
  - "(1/K)log|c_K| from K = 80 to 160" is a difference quotient of a sign-changing sequence.
- **D3cn-8 (cosmetic).** Three of my mutants change the raw only (the md is identical on recompute):
  - E07 (401 points in s);
  - E08 (qx at 3K + 80 digits: only the stored `prec` changes, the values are identical);
  - E10 (40 radii in `dom_radius`).

## What I checked and how

**Reproducibility** [own]
- `--from-raw` regenerates gh_peano.md byte for byte (sha256 60606f73…).
- tests/test_gh_peano.py: 22 passed in 1.58 s.
- With the script's own job functions (W/scripts/rerun.py; logs in W/logs/rerun_*.log), every job of the raw is bit-identical (`repr` equality, `secs` excluded):
  - P and M;
  - L (8 jobs), S (12), R (22, and the 14 without a frame), LC (12), Q (12);
  - part D: all 20 `d_job`s, the 5 `duhamel_job`s, the 4 `qx_job`s and `part_dom`.

**Claims against the md and the raw** [md][raw][own]
- W/own/numcheck.py extracts 180 numeric tokens from the draft. All are in the md verbatim, except arXiv numbers and "Theorem 2.10".
- I spot-checked the semantics of these:
  - 1.8927 at 0.345;
  - 1.7773 / 0.6345, 1.6904 / 0.5279, 0.8053, 9.27;
  - the step's K·deviation 0.103 … 0.044 and ratios 0.488 … 0.518;
  - 8.811 and 15.540 (pre-asymptotic);
  - the quadrature 3.1 … 8.7·10⁻¹⁵;
  - beyond the disc +2.7, +5.1, +6.0·10⁻⁴;
  - the slopes +0.051 … +1.480.

**(iii)(a), re-derived** [own]
- **He_K(x) = E(x + iZ)^K.** Gaussian integration by parts gives the recurrence.
- **The contour shift.** Substitute t = s + iζ (a horizontal line, on which the entire integrand decays uniformly). Then ζ + 1/ζ + it = ζ^{−1}(1 + iζs), which gives exactly K^{K/2}ζ^{−K}e^{Kζ²/2}I_K(ζ).
  - Checked to 10⁻³⁹ at 40 digits, at K = 8 and 40 and three values of ζ.
- **The Gaussian bound.** It is log(1 + x) ≤ x with x = −2Im(ζ)t + |ζ|²t², and |e^{−iKζt}| = e^{K Im(ζ)t}. Together these bound the integrand by e^{−K(1 − |ζ|²)t²/2}.
  - Checked at 2·10⁵ random points: max(lhs − rhs) = −4.1·10⁻¹².
- **The limit of I_K.** h″(0) = ζ² − 1, the odd terms vanish, so I_K → (1 − ζ²)^{−1/2} with O(1/K) uniformly on |ζ| ≤ r. At ζ = 0.5 + 0.2i the deviation is 1.41·10⁻² at K = 8 and 2.99·10⁻³ at K = 40.
- **The ellipse and I_K's zeros.** ‖J_K‖ ≤ 2√(K − 1) < 2√K < √K(ρ + 1/ρ), so the ellipse encloses the nodes for every K. Since ζ + 1/ζ is real with |ζ| < 1 only for real ζ, I_K has no zeros on 0 < |ζ| < 1.
- **The contour formula.** dz = √K(1 − ζ^{−2})dζ, and z runs counterclockwise when ζ runs clockwise. The normalisation (K − 1)!K^{(1−K)/2}(√Kμ)^{1−K}e^{−Kμ²/2} is right.
- **The stationary points.** f_μ′ ∝ (ζ² − 1)(μ − ζ), so ζ = μ is the only stationary point inside the disc. f_μ(μ) and f_μ″(μ) are as stated.
- **The circle identity.** It holds to 2.2·10⁻¹⁵ at 10⁵ random points.
- **Laplace's method on the circle.** In ψ the exponent is −K(1 − μ²)ψ²/2 with Re(1 − μ²) ≥ 1 − |μ|² > 0. The principal roots cancel against I_K's. The result is (1 − μ²)·K!K^{−K}e^K(2πK)^{−1/2} → 1 − μ².
- **The maximum principle** applies because c_K is entire (it has the Dirichlet form).

**(iii)(b), re-derived** [own]
- **The phase.** G′ = (1 − 2s)(1/(s(1 − s)) − η²); G′(½) = 0, G(½) = −log 4 + η²/4, G″(½) = −2(4 − η²) with Re(−G″) > 0, and the principal root.
- **The constant.** 2(4 − η²)^{−1/2}(1 − η²/4) = (1 − η²/4)^{1/2}.
- **The bound beyond the disc.**
  - The contour on |ζ| = ρ₀ is exact for any σ. For fixed ρ₀, 1/I_K is bounded on that circle. The prefactors give −Re f_σ(σ).
  - max_θ(2a cos θ − ½cos 2θ) = M(a) + ½, and the convergence as ρ₀ → 1 is uniform on compacts. For real σ the bound decreases in ρ₀, so ρ₀ → 1 is optimal there.
  - Inside the disc Λ⁺(σ) = |σ|²/2 − ½ − log|σ| ≥ 0, which tends to 0 at |σ| = 1. This is what the band argument of D3cn-5 needs.
- **Real axis.** The total is log w + w²/2 + 2(η − w) − 3/2 − 2 log η, with derivative (w − 1)²/w. I verified the gap, its vanishing to second order at η = 2, and the third derivative, positive on (1, 2] (positive up to η = 2^{1/3}/(2^{1/3} − 1) ≈ 4.85).
- **Imaginary axis.** The total is log w − w²/2 − ½ − 2 log y, with derivative 1/w − w > 0. I verified the gap, its vanishing to first order, and the second derivative, negative on (1, 2].

**The domination condition** [own]
- **Code.** W/own/dom_exact.py, dom_profile.py, dom_where.py, dom_bb.py.
- **The closed form of D3cn-3.**
  - It agrees with `dom_gap` on the same s-grid to 3.1·10⁻¹⁵.
  - The exact supremum exceeds the 4001-point grid supremum by at most 5.3·10⁻⁹.
  - On 4001 angles the least radius is 1.8927327 at 0.3432; refined, 1.892732715 at 0.34318.
  - The violated sets are contiguous.
- **Margins on |η| ≤ 1.7, 1.7773, 1.78, 1.85:** −9.137·10⁻³, −3.426·10⁻³, −3.294·10⁻³, −9.708·10⁻⁴, all at φ = 0. For R ≤ 1.85 the gap decreases with arg η. At R = 1.8927 a bump near φ ≈ 0.31–0.37 reaches −4.9·10⁻⁵.
- **Robustness.** The grid's least radius changes by only 4·10⁻⁶ with 401 points in s instead of 4001 (mutant E07).
- **Branch and bound.** On each box [R₁, R₂] × [φ₁, φ₂] I bound the gap by log w + (p_lo/2)(w² − ½) + M(R₂(1 − w)cos φ₁) − ½ − log R₁ + log 4 on w ≤ 1 − 1/R₂, which is valid because w² < ½ and M is increasing. Its supremum over w is exact. With a slack of 10⁻¹² (floating point, not interval arithmetic), no box stays unresolved up to |η| = 1.7773, 1.78, 1.85, 1.89 and 1.8927.
- **Where the violation sits.** At |η| = 1.95, φ = 0.087 … 0.611, the maximiser has Re(ηs*) = 1.023–1.099 (the edge regime). At 1.999, φ = 0.698, it has 0.961.

**(iv), re-derived** [own]
- The tail bound: |ε̃_K| ≤ 1 + e^{y²/2}, min(A, B − 1) = (u − 1)/2, and c² > 4 ln u/(u − 1) > 4r (because 1 − 1/u < ln u).
- At the cusp the angular spread sits at |η| = 2√r < c, and off the cusp the angle concentrates. The model's tail beyond c²K is exponentially small relative to its peak. Nothing needs the limit outside |η| ≤ c.
- Both closed forms (off the cusp and at it) and the half-squeeze constant (1 − t²)^{1/4}(πK)^{1/4}t^K ₂F₁^{1/2} re-derived.
- r* = 0.805260 from the grid radius (0.805257 refined), giving σ*/σ < 9.27009 (9.26993 refined). c_need(0.6345) = 1.77733 and c_need(0.5279) = 1.69039; the worst rule's Laplace mass sits at |η| = 1.593.

**The Dirichlet form.** Var T = Σz²/(K(K + 1)) from the Dirichlet covariances, with Σz² = K(K − 1). The second-order term K(Var T − 1)/2 → −1 is consistent with log(1 − μ²) [own].

**Beyond the disc** [own]
- **Code.** W/own/alias_exact.py, alias_prony.py, corner_complex.py.
- **Exact c_K from integers.** c_K(u) = (K − 1)!e^{−u²/2}Σ_n u^{n−K+1}N_n/n!, where N_n = (A^n e₀)_{K−1} and A is multiplication by z in the monic basis modulo He_K. A has integer entries, so the N_n are exact. This reproduces the raw's qx values; for example log|c_160(2.5i√160)| = 233.548, against the raw's 2.68·10¹⁰¹.
- **Prony fits** over 13 consecutive K at K ≈ 160, 320, 640:

  | y | Prony rate | λ₁ | heuristic | frequency (= Im(f₁(ζ_s) − f_σ(σ)) mod 2π) |
  |---|---|---|---|---|
  | 1.9 | 0.07470, 0.07468, 0.07468 | 0.07468 | 0.06426 | 2.13177 |
  | 2 | 0.26716, 0.26714, 0.26713 | 0.26713 | 0.25130 | 2.28540 |
  | 2.2 | 0.69684, 0.69683, 0.69682 | 0.69682 | 0.66932 | 2.56931 |
  | 2.5 | 1.45026, 1.45024, 1.45024 | 1.45024 | 1.40315 | 2.92170 |

- **Plain slopes** over the same windows vary 0.037–0.158 at y = 1.9 and 0.632–0.878 at y = 2.2.
- **Off the axes**, the slopes for K = 160 → 1280 against the candidate rates:

  | σ | slopes | candidate rates |
  |---|---|---|
  | 1.3e^{0.8i} | 0.08065, 0.08061, 0.08064 | first alias 0.08065; Λ⁺ 0.08264; edge 0.07375 |
  | 2e^{1.2i} | 0.80361 (three times) | first alias 0.80361; Λ⁺ 0.80685 |
  | 1.17e^{0.5i} | 0.0267 | Λ⁺ = edge = 0.02673 |

**Mutants** [own]
- **My spec**, W/mut/spec93.py: 32 single-point mutants of round 91's part D that its specs leave out. The tests kill 27.
- **Survivors, recomputed** with W/mut/recompute93.py:
  - E07, E08 and E10 change the raw only; the md is identical.
  - E17 (r_star's upper end 0.9) and E27 (|Re μ| without the absolute value; the angles lie in [0, π/2]) are equivalent on the data.
- **Why the report-level mutants die.** `test_the_md_is_the_report_of_the_raw` kills every report-level mutant that changes the md.
- **Earlier specs, rerun:**

  | spec | killed | survivors | stale anchors |
  |---|---|---|---|
  | round91_D61_own | 29 of 31 | D17, D31 (equivalent, as logged) | none |
  | round91_D61_referee | 31 of 35 | QP2, FA8, FA9, RR2 (equivalent, as logged); QJ4, QJ5, QP1, QM1, SJ1, LJ2, LJ3 now killed, as the log says | none |
  | round90 own | 28 | FA1, FA4, RR1 | 3 |
  | round90 referee | 50 | M6, S22, O3 | 14 |

  The round-90 results match the log.

**Literature** (web search; arxiv.org not fetchable)
- Ye bounds the corner entry e_kᵀe^{itT_k}e₁ a priori, with Faber polynomials [snippet].
- Beckermann–Kuijlaars give superlinear Krylov convergence by constrained equilibrium problems (nth-root rates), and Beckermann–Güttel treat rational Arnoldi [snippets].
- arXiv:2312.07940 gives root-exponential Gauss–Hermite bounds by contour integrals [snippet].
- de Boor: B-splines as Peano kernels of divided differences, and Hermite–Genocchi [snippet].
- I found no sharp asymptotics of the corner entry of e^{uJ_K}, of the Gauss–Hermite error of e^{ξz} at ξ ~ √K, or of Q, P or the domination argument. The Krylov literature gives bounds, not constants. Saad 1992's a posteriori estimate is the corner entry [memory].

## Verdicts per part

- **(i): proved** (unchanged).
- **(ii): proved,** except that P_j's polynomial form remains a sound sketch (unchanged).
- **(iii)(a): proved.** Every step is right: the representation, the bound, the ellipse, the contour formula and its orientation, the circle identity, Laplace's method with Re(1 − μ²) > 0, and the maximum principle. Two remarks only: the ellipse works for every K, and I_K is zero-free (D3cn-7).
- **(iii)(b): proved** in three pieces:
  - on |η| < 1, unconditionally;
  - on (−2, 2) and i(−2, 2) (both axis proofs redone and correct);
  - on |η| ≤ 1.8927, modulo the domination inequality, which my branch and bound confirms on the whole disc (D3cn-3).

  The band below s = 1/|η| and Ω's openness need one sentence each (D3cn-5). "Λ⁺ attained" is overstated (D3cn-2). The sector stays open; its failure at |η| = 1.95 is genuine for the straight segment.
- **(iii)(c): proved.** It needs only a disc |η| < 1.
- **The Dirichlet form and Var T: proved.**
- **Heuristics.**
  - The alias rate and the threshold 1.8637 are wrong (D3cn-1).
  - The Stieltjes identity is proved, and its saddle-point use is correctly marked formal.
- **(iv): correct.**
  - Nothing beyond |η| ≤ c is used: the tails go by the crude bound, and the cusp's spread sits at |η| = 2√r.
  - r* and 9.27 round slightly up (D3cn-6).
  - At all 22 rules it holds modulo the grid, or by the branch and bound. For |η| ≤ 1.78 the binding constraint is the real axis, which is proved.
- **(v): proved** (unchanged).

## Status of the round-91 issues

| id | status |
|---|---|
| D3cl-1 | Resolved. The domination condition is stated, the sector is open, and the failure is genuine at \|η\| = 1.95. Residue: Λ⁺ is "attained" only for \|Re σ\| > 1 (D3cn-2). |
| D3cl-2 | Resolved. Q is proved, P holds on Ω, and (iv) holds for r < 0.8053. Residue: New omits the grid proviso (D3cn-4). |
| D3cl-3 | Resolved. Part D exists and is bit-identical on rerun. |
| D3cl-4 | Resolved. Saad 1992 [memory] and the Krylov literature are cited; New claims the analysis, not the identity. |
| D3cl-5 | Resolved. All result-changing survivors are killed, and the spec is kept and rerun as logged. |
| D3cl-6 | Resolved, with a test. |
| D3cl-7 | Resolved. |
| D3cl-8 | Resolved in the Status. The continuation beyond the disc that it prompted carries a wrong rate (D3cn-1). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cn-1 | minor | The alias rate (11/18)y² − log y − 3/2 uses a real stationary-phase point under an exponentially varying modulus. The first alias's exponent is f₁(ζ) = σ(ζ + 1/ζ) + 3 log ζ − ζ² + ζ^{−2}/2, with f₁′ ∝ (ζ² − 1)(2ζ² − σζ − 1). Its saddles lie on \|ζ\| = 1/√2, giving λ₁ = (5/8)y² − log y − ½ − (3/2)log 2 with zero 1.8589. Prony fits on 13 consecutive K at K ≈ 160–640 give λ₁ to 5 digits (0.07468, 0.26713, 0.69682, 1.45024 at y = 1.9, 2, 2.2, 2.5) against the heuristic's 0.064, 0.251, 0.669, 1.403, and match the oscillation frequency. Two-point slopes oscillate (0.037–0.158 over windows at y = 1.9), so they are no evidence [own]. | Replace the rate and the threshold 1.8637 → 1.8589. Add a consecutive-K fit to part D (the integer series is cheap: K = 160 takes under 1 s). Withdraw "matched by the slopes". D61 (iii), Numerics (D), Status, New; `alias_rate` and its test pin. |
| D3cn-2 | minor | "Its rate is Λ⁺, attained there" holds only where \|Re σ\| > 1 (the edge saddles ζ = ±1; 1.17e^{0.5i}: 0.0267 = Λ⁺). For \|Re σ\| < 1 the growth is the first alias's, below Λ⁺: 0.08064 against 0.08264 at 1.3e^{0.8i}; 0.80361 against 0.80685 at 2e^{1.2i}; 0.267 against 0.807 at 2i. Also "exponentially large off the imaginary axis" is imprecise: c_K decays near the real axis. In the violated sector at 1.95 the maximiser has Re(ηs*) ∈ [1.009, 1.099], so the failure there is genuine [own]. | "Λ⁺ is a bound, attained where \|Re σ\| > 1." The sector's failure is genuine where the maximiser has \|Re(ηs*)\| > 1. D61 (iii)(b) text, Status. |
| D3cn-3 | minor | The condition is gap = log w + (p/2)(w² − ½) + M(R(1 − w)cos φ) − ½ − log R + log 4 on w ∈ (0, 1 − 1/R], with its supremum over w in closed form. A branch and bound (slack 10⁻¹²) confirms it on 1 < \|η\| ≤ 1.8927 (816,881 boxes) and on ≤ 1.78 (7,631 boxes). The quoted margins are the real axis's, at φ = 0, where the gap is proved; for R ≤ 1.85 the gap decreases with arg η [own]. | State the closed form. Verify it with interval arithmetic (mpmath.iv or arb) in a `dom_verify`, and replace "checked on a grid" by a computer-assisted proof, at least for (iv)'s needed radius 1.78. Say where the margins sit. D61 (iii)(b), Status; exp_gh_peano part D. |
| D3cn-4 | minor | New: "it proves both of D‴.58's conjectures where r < 0.8053", without the grid proviso; D‴.58 (i)'s statement likewise (its Status has the proviso). Conversely, the Status leaves out that (b) on \|η\| < 1 and (c) are unconditional [md]. | Add "modulo the grid-checked condition", or drop it after D3cn-3. Mark (b) on \|η\| < 1 and (c) as unconditional. D61 New, Status; D58 (i). |
| D3cn-5 | cosmetic | The proof splits into "\|ηs\| ≤ r by (a); beyond, the ρ₀-contour". The band r < \|ηs\| < 1 lies outside the stated sup (s ≥ 1/\|η\|). Ω's openness and the uniform gap on compacts are also unstated [own]. | One sentence: continuity of g + Λ⁺ at s = 1/\|η\|, with Λ⁺ = \|σ\|²/2 − ½ − log\|σ\| → 0 inside the disc; and the sup is continuous in η and tends to −∞ as \|η\| ↓ 1. D61 (iii)(b) proof. |
| D3cn-6 | cosmetic | The exact-sup least radius is 1.8927327 at 0.3432 (the grid's value is 1.8·10⁻⁶ high). That gives r* = 0.805257 and σ*/σ < 9.2699, so 0.8053 and 9.27 round up. The violated sectors are [0.0864, 0.6173] and [0.0016, 0.7662], against the grid's [0.087, 0.615] and [0.004, 0.764] [own]. | Write r < 0.805 and σ*/σ < 9.26, or the refined values, and call the sectors "at least". D61 (iii), (iv), Status, New. |
| D3cn-7 | cosmetic | The ellipse encloses the nodes for every K (‖J_K‖ ≤ 2√(K − 1)). I_K is zero-free on 0 < \|ζ\| < 1, which the division by I_K needs. "(so c_K > 0 on ℝ)" follows from the Dirichlet form. The interval (τ²/σ*², τ²) is written in σ = 1 units beside σ-formulas. The md and the script still call Q "conjecture Q". "(1/K)log\|c_K\| from K = 80 to 160" is a difference quotient [md][own]. | One phrase each. D61 (iii), (iv), (v); exp_gh_peano docstrings and headings. |
| D3cn-8 | cosmetic | Raw-only survivors: E07 (401 points in s), E08 (qx at 3K + 80: only the stored `prec` changes, so 4K + 80 is not needed), E10 (40 radii). E17 and E27 are equivalent on the data [own]. | Pin qx's `prec` and part_dom's resolution, as `test_round91_stored_fields` does for LC. tests/test_gh_peano.py. |

## Sources

Web searches only. Everything below is [snippet] unless marked otherwise.
- [Ye, Error bounds for the Krylov subspace methods for computations of matrix exponentials (arXiv:1603.07358)](https://arxiv.org/pdf/1603.07358): a priori bounds on the corner entry e_kᵀe^{itT_k}e₁ via Faber polynomials.
- [Ye, report M106373 (UKY)](https://www.ms.uky.edu/~qye/reports/M106373.pdf): Lanczos convergence and the decay of the exponential of a Hessenberg matrix.
- [Beckermann–Kuijlaars, ETNA 14 (2002) 1–19](https://etna.math.kent.edu/vol.14.2002/pp1-19.dir/abstr1-19.pdf): superlinear convergence via constrained equilibrium problems.
- [Beckermann–Güttel, rational Arnoldi (Numer. Math.)](https://pro.univ-lille.fr/fileadmin/user_upload/pages_pros/bernhard_beckermann/abs/revised_Beckermann_Guettel_Numer_Math.pdf)
- [Convergence analysis of Hermite approximations for analytic functions (arXiv:2312.07940)](https://arxiv.org/pdf/2312.07940): contour integrals; root-exponential Gauss–Hermite convergence.
- [de Boor, Divided differences (arXiv:math/0502036)](https://arxiv.org/pdf/math/0502036): B-splines as Peano kernels, and Hermite–Genocchi.
- [Gil–Segura–Temme, nodes and weights of Gauss–Hermite (arXiv:1709.09656)](https://arxiv.org/pdf/1709.09656)

From memory [memory], not seen here: Saad, SIAM J. Numer. Anal. 29 (1992); Plancherel–Rotach (1929); Olver, *Asymptotics and Special Functions* (Laplace's method with complex phase).

## Lesson

1. **A stationary phase is not a saddle when the modulus varies on the same exponential scale.** The alias rate came from the modulus at a real stationary-phase point. The true saddle sat off the circle, at |ζ| = 1/√2.
2. **Measure growth rates by fits over consecutive K, not two-point slopes.** The sequence changed sign from one K to the next. Prony fits over 13 consecutive K gave the rate and the frequency to 5 digits, where the slopes scattered by ±0.1.
3. **Reduce an inequality to closed form before gridding it.** The supremum over s was explicit. A two-parameter branch and bound then confirmed the domination condition on the whole disc |η| ≤ 1.8927. It also showed that the quoted margins were the real axis's, which is the proved part.

Files are in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref93:
- own/: dom_exact.py, dom_bb.py, dom_profile.py, dom_where.py, alias_exact.py, alias_prony.py, corner_complex.py, numcheck.py;
- scripts/rerun.py;
- mut/: spec93.py, recompute93.py, and the result files under w93, w91o, w91r, w90o, w90r;
- logs/.
