# Referee report D3ea, round 132: Proposition D‴.62 (the rules without a Mehler frame), first review

## Scope

I reviewed commit 8ed0f31 (branch claude/vigilant-clarke-hf939b), exported with `git archive` into W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref132, with the repository's venv symlinked as W/repo/venv. I wrote nothing into the repository. The working tree's uncommitted changes, made by other sessions, were never read.

I read in full: review/drafts/D62.md, review/RESEARCH_note_frameless.md, experiments/exp_frameless.py, results/chaos/frameless.md, tests/test_frameless.py and tools/mutants/d62_initial_own.py. From D61.md I read (i)–(v), its proof and its Status, in particular the axis gaps and the openness and uniformity argument. From D58.md I read (i), (ii), (v), its Status and A6. I also read PAPERS.md §18, the HANDOFF entries and the round-129 report's c₂ passages.

The proof was re-derived step by step by hand. The numerics are my own code: W/own/own_err.py, check_sample.py, check_repr.py, new_rules.py, fit_new.py, check_F1.py, rich_stored.py and kw_threshold.py, with logs alongside in W/own. I ran 42 mutants of my own plus a control (W/own/mutants_ref132.py, results in W/mut/res_a.jsonl). I ran one heavy process at a time and killed only one PID of my own, after checking its cmdline.

## Summary

**Design.** I checked:
- reproducibility: the md from the raw, the 23 tests, and a sample of S;
- every step of (a)–(e) by hand: the change of variables, the Fubini cases, D‴.61 (iii)(b) on the axes (both gap functions recomputed), Step 2's exponents, Lemmas A, T and B, the order of choosing c′ and δ, the region off the box, moment matching, ρ = −1, the formal s₁, and the identity in (e);
- independent numerics, with my own exact-error code: nodes are certified roots of He_K from flint's `fmpz_poly.complex_roots`, and the error is assembled in the z-variables as E⊗E, E⊗Q_K and Q_K⊗Q_K of F. This shares no code with `err2` or `gh_nodes`.
- my own trapezoidal check of the representation, including frame rules and ρ = −1;
- exact errors at six rules the repository does not contain, a numerical test of assumption (F), and a literature search.

**Verdicts.**
- (a) correct.
- (b) correct and complete: I found no gap.
- (c) the two sector asymptotics and their positivity are correct; the transient is honestly labelled formal.
- (d) the formal derivation is algebraically right, and s₁ is confirmed numerically to 10⁻⁸–10⁻³⁰ at new rules.
- (e) correct. The claim of "no computer-assisted step" is correct.

There are no major findings.

**Findings.**
- **D3ea-1 (minor, literature/novelty).** Kuo–Woźniakowski 2012 (as quoted in snippets) claims a convergence threshold 2αγ² = 1 for the Gauss–Hermite rule, with exponential divergence beyond it. Under the natural mapping, D‴.62 (b) at moment matching proves convergence on both sides of that threshold. D‴.62's own cheap bound has exactly that threshold. D62.md does not engage with this.
- **D3ea-2 (minor).** Five semantic mutant survivors would change published R or T numbers on a rerun. A sixth (E30) survives only for the md-byte test's reasons and changes nothing printed.
- **D3ea-3 (minor).** The case ρ = −1 is proved but nowhere checked numerically. My own checks of it pass.
- **D3ea-4 (minor).** Two "What it changes" proposals are slightly over-stated: the half-squeeze gain formula's caveat, and the "degenerate cases" wording for |ρ| > 1.
- **D3ea-5 (minor).** The sentence saying the fitted s₂ are "what (d)'s (1 − t)^{−3} predicts" is not supported.
- **D3ea-6 (cosmetic).** "Round 129's formula (D‴.58's Status)" points to text that is not in D58.md at this commit.
- **D3ea-7 (cosmetic).** The ratio in (e)'s proof is inverted (harmless, since it equals 1).
- **D3ea-8 (cosmetic).** Two wording and caveat points in the Numerics paragraph.

## What I checked and how

**(1) Reproducibility.**
- `exp_frameless.py --from-raw` with FL_OUT in W produced frameless.md byte for byte (sha256 43a5c492…162ec, identical to the committed file) [own].
- 23 of 23 tests pass, in 11.5 s [own].
- S sample: 17 (rule, K) points, including K = 400 at four rules and K = 360 at two, all 14 rules covered. My certified enclosures (relative radius ≤ 10⁻¹⁰¹) agree with the stored 60-digit strings within 3.4·10⁻⁶⁰, i.e. the stored rounding [own]. This also confirms, at those points, that S's radius omits the error of the nodes harmlessly (D3ea-8).
- A first pass with too small a margin had my own radius above 10⁻¹²; I reran it with a margin of 300 + K bits.

**(2) The proof, by hand.**

*(a)*
- e_K² = σ*(μ − Q)^{⊗2}[κ]: I checked μ⊗μ[κ] = 1/σ* and the cross term's √2(1 + σ*²)^{−1/2} factor.
- The change of variables: g = (v/τ²)^{1/2}e^{cz²}, so c(z² + z′²) − v(z − z′)²/4 = c(z + z′)²/2 + c_b(z − z′)²/2. The terms are finite because c and c_b are below ½.
- The Gaussian representation holds for either sign of the variance.
- Fubini–Tonelli: c_b < c, so the three cases are exhaustive. At a frame the a-expectation e^{c(z+z′)²/2} is integrable under each of the four product measures.
- The frame's constants match D‴.61: A = 1/(2c) − 1 = v/(τ² − v), B = 1 + 1/(2|c_b|) = vσ*²/(vσ*² − τ²), and (2π)^{−1}(c|c_b|)^{−1/2} = (1 − ρ²)^{1/2}/(πa), using c·c_b = −(1 − ρ²)v²/16.
- (i²)^{2K} = 1. The argument holds for any finite atomic Q_K.

*(b), Step 1.*
- Stirling: (K!/(2K)!)K^K ~ 2^{−1/2}(e/4)^K, and (K!/(2K)!)²(4K − 1)!! = C(4K, 2K)/(4^K C(2K, K)) = 2^{−1/2}(1 + 1/(16K) + …).
- D‴.61 (iii)(b) on the axes, recomputed:
  - Real axis: on s ≥ 1/η, g_η + Λ⁺ has s-derivative −η(w − 1)²/w ≤ 0, so its supremum is at s = 1/η. The gap there vanishes at η = 2 together with its first two derivatives, and its third derivative 2/(η − 1)³ − 4/η³ is positive on (1, 2], so the gap is negative on (1, 2).
  - Imaginary axis: Λ⁺(iys) has the same analytic form along the whole segment. The derivative is y(w − 1/w) < 0 for w ≤ y − 1 < 1. The gap vanishes with its first derivative at y = 2 and is concave on (1, 2].
  - The gap tends to −∞ as |η| ↓ 1, so the margin is uniform on compact subsets of each interval.
  - The inputs are by hand: (iii)(a) (Plancherel–Rotach type with proof, the circle identity, the maximum principle), Duhamel's identity (round 127, rechecked in round 131), and the contour bound on |ζ| = ρ₀ < 1. That contour representation is valid for every σ, since the integrand is entire in u.
  - None of D, V, VA or the c_need constants enters. "No computer-assisted step" is correct.
- Relative form: on the real axis the limit (1 − y²/4)^{1/2} is at least (1 − c′²/4)^{1/2} > 0, and on the imaginary axis (1 + y²/4)^{1/2} ≥ 1, so the additive o(1) becomes a relative (1 + o(1)). The factor ξ^{2K} is exact, so the vanishing near y = 0 is harmless; ĝ_K(0) = 1.
- Step 1 (iii): |φ_K(iw)| ≤ 2 needs only positive weights; φ_K ≥ 0 and ε̃_K < 1 come from D‴.61 (i).

*Step 2.*
- The Jacobian K/2 gives the density (K/4π)(|c||c_b|)^{−1/2}e^{−KQ}. In p and p′ the exponent Φ is as stated, with ½ + 1/(2|c|) = 1/(2|s_a|).
- ρ < −1: Q − (y₁² + y₂²)/2 has the c̃ = c/(1 − 2c) form, with c̃/(1 + c̃) = s_a, and the tail rate is 0 after the same subtraction.

*Lemma A.* Φ decreases in u on w ≥ u, its maximum is at w = 4|s_b|, and the value there is 2 log|s_b|. Correct.

*Lemma T.*
- α > |β|, and the sector is {y₁y₂ ≥ 0}.
- R1: both cases of the sign of β give Q ≥ 2/k_a on that sector (2/k_b on the other), with equality only at the corners, and log(1 + 1/k) < 1/k.
- R2: the minimiser |β/α||y₂| < 2, so each half-line's minimum of Q is at its end. The end keeps its sector, and L(±2) = 0. The ridge points are interior with |y_i| = 2|s|^{1/2} < 2, so the closed edge is strictly below 2 log|s|. Upper semicontinuity holds, with −∞ at y₂ = 0.
- Breaking attempts: the margin tends to 0 as t(v) → 1 (2.5·10⁻⁶ at t = 0.998), but strictness for each rule is all the proof needs; no uniformity in v is claimed. The cusp is a frame rule, so Lemma A's |s_a| ≠ |s_b| always holds without a frame. At ϖ → 0 nothing degenerates for fixed ϖ > 0, and ρ = ±1 is treated separately in one dimension.
- I found no parameters at which the tail supremum reaches 2 log t.

*Step 3.*
- The box is |A| + |B| < c′√K. M_K ≥ 0, so the uniform factor passes through the expectation.
- The tilt is exact: N = ((1 − c)(1 − c_b))^{−1/2} and σ² = |s|.
- Lemma B: the size-bias is exact (Ỹ² ~ χ²_{4K+1}). (1 − ϖX²/Ỹ²)^{2K} lies in [0, 1] and tends to e^{−ϖX²/2}; m_K → 1 − s_d; the indicators tend to 1 because 2σ_d < c′. The limit is bounded convergence in probability, and E e^{−ϖX²/2} = (1 + ϖ)^{−1/2}.
- The region off the box, including where one coordinate is in the box and the other out, is covered by L̃. L̃ ≤ L⁺ − L(c′) holds in all three ranges of |y|, and the bound e^{KΦ̃} ≤ e^{(K−1)S}e^{Φ̃} is valid.
- Suprema of the coercive upper-semicontinuous Φ⁺ over the decreasing closed sets converge to the supremum over their intersection.
- Order of choices: δ comes from Lemma T, then c′ ∈ (2t^{1/2}, 2) is taken close to 2. This is consistent.
- Collecting the constant reproduces C(v), with 1 − c = (τ² + v)/(2τ²) and 1 − c_b = (τ² + vσ*²)/(2τ²).
- Moment matching: erfc(y) ≤ e^{−y²} gives P(|B| ≥ c′√K) ≤ e^{−c′²K/τ²}. The condition c′² > 2τ² ln(1/λ) = 4 ln(1 + u)/u lies in (4/(1 + u), 4).
- ρ = −1: c = τ²/(2σ*²), and −2/c̃ < 2 log s_a.
- Nothing is circular. The Laplace mass sits on the compact set |y| ≤ c′ where the o(1) is uniform, and the box boundary is fixed in the scaled variables.

*(c).* The sector forms hold verbatim. C_o > 0, since 1 − s_o > 0 on both branches. I re-derived the expansion to first order in q^K: (S(K − 1)/S(K))^{1/2} = 1 + O(K⁻²), so B₁/B₀ = s₁′ − s₁, and B₀ = 0.23832 by hand [own].

*(d).*
- I recomputed every term: Eζ = 1/(4K), Eζ² = 1/(2K), ω = 1/(1 + ϖ), E W² = 3ω². The ± cross terms cancel at first order and leave f₀″s²ϖW/K at second order. The prose "the cross term's square is s²ϖWỸ²" is true for σ_dσ_oXỸ without the 2.
- The ±3ϖ²ω²/16 terms cancel, s²f₀″/2 + s²f₀′² = 0, sf₀′/2 cancels 2F₁'s last term, and the ϖ terms reduce to −(ϖω/8)(1 + s)²/(1 − s)². The displayed s₁ follows.
- At ϖ = 0 the one-dimensional Gamma computation gives the same s₁. The rationals −13/125, −11/64, −187/2000 and 8/4913 check, and c₂ = s₁/2 follows.

*(e).* C₆₁/C = 2((τ² − v)(vσ*² − τ²))^{1/2}/(vτ²(1 − ρ²)^{1/2}) = 1 by the identity. The draft prints the reciprocal (D3ea-7). The duality check u ↦ −u works.

**(3) Numerics, independent of the repository's code.**
- *Representation (a)*, by my own trapezoid in the Gaussians' own coordinates (float64), against my own exact e_K². Relative differences are 4·10⁻¹⁶ to 7.8·10⁻¹⁴ at:
  - moment matching at τ/σ = 1 (K = 6) and 2 (K = 10);
  - ρ > 1 at (¼, 2v_hs) and (1, 2v_hs);
  - ρ < −1 at (½, 0.5v_hs), (1, 0.5v_hs) and the new (0.7, 0.3v_hs);
  - ρ = −1 at τ = 1, v = ½ (2.0·10⁻¹⁴);
  - four frame rules, (1, 0.8v_hs), (2, 1.25v_hs), (3, 3v_hs) and (1, Fock), which confirms D‴.61 (v)'s identity [own].
- *C(v) at large K*, own exact errors at K = 30 … 300 for six rules not in the repository, fitted by least squares with the limit free:

  | rule | ρ | t(v) | \|L − 1\| | s₁ against formal (relative) |
  |---|---|---|---|---|
  | τ = 1, v = ½ | −1 | 1/3 | 1.4·10⁻¹⁰ | 13/4 within 1.4·10⁻⁸ |
  | moment matching, τ = 3 | 1 | 9/11 | 10⁻²⁶ | −0.048875 within 5·10⁻²³ |
  | (3, v = 10.8) | 1.037 | 11/13 | 10⁻²⁶ | within 6·10⁻²³ |
  | (1, 4v_hs) | 2.29 | 0.700 | 1.2·10⁻³³ | within 1.7·10⁻³⁰; g₀ = C_o/C within 1.2·10⁻⁷ (q = 0.466) |
  | (0.7, 0.5v_hs) | −4.88 | 0.419 | 1.7·10⁻⁹ | within 1.1·10⁻⁷ |
  | frame (0.7, 1.15v_hs) | 0.749 | 0.168 | 1.8·10⁻²⁷ | round 129's s₁ within 2.8·10⁻²⁴ |

  [own]
- *Stored S data*, my own exact interpolation in 1/K (not least squares): s₁ is reproduced to 10⁻²⁷–10⁻³⁹ at the nine resolvable fast-branch rules, and to 2·10⁻¹¹ … 9·10⁻⁹ on the slow branch, better than the md's 3.5·10⁻⁷ at (1, 0.5v_hs) [own].
- *Assumption (F)*: K[log ĝ_K(√Kη) − Kη²/4 − ½ log(1 − η²/4)] at K = 60 … 320, extrapolated:
  - imaginary axis: equals F₁(η²/4) to 8 digits at y = 1, 1.5 and 1.9;
  - real axis: 0.837960 against 0.837963 at η = 1;
  - pre-asymptotic at η = 1.5 and 1.8 (10.9 against 11.0, 143 against 217), the known slow convergence near η = 2. Nothing contradicts (F) [own].
- *Frame rules not covered.* At a frame, a is real and b imaginary. No real measure on ℝ² can have the moment generating function e^{q} with q indefinite, and no measure on (iℝ)² can either (|characteristic function| ≤ 1). So the Status claim, restricted as it is to Gaussian measures, is right, and C(v) with ϖ < 0 does match at frames, through D‴.61 (iv).

**(5) Mutants.** 42 single-point mutants beyond d62_initial_own.py, plus the control.
- 33 are killed, 26 of them also changing the md.
- Every md-changing mutant is killed by the md-equals-report test, so the question is the computing code.
- Nine survive. Five of them would change published numbers on a rerun:

  | mutant | change | published numbers affected |
  |---|---|---|
  | E20 | r_2d's b perturbed by a relative 2⁻¹⁰⁰ | R's two-dimensional "at most 1.0·10⁻⁴⁶" becomes about 10⁻³⁰ |
  | E21 | dM from .mid() instead of .upper() | probably none (1-digit print) |
  | E32 | t_grid's L⁺ raised by 10⁻³ on the imaginary branch | T's grid margins for the five ρ > 1 rows, and "grid's maximum never above 2 log t(v)" |
  | E33 | sweep upper branch starts at 1.01σ* | sweep minima |
  | E40 | sweep lower end 0.03 | sweep minima, incl. "the smallest 2.5·10⁻⁶ … v/v_hs = 0.02" |
  | E36 | tail_b mask y₁y₂ < 0 | none (unprinted field) |

- Equivalent on the data: E30 (fit midpoints not taken exact), E35 (R1 < R2 in sector b) and E38 (the middle node not zeroed: Newton converges to 0 anyway).
- Cause: T's grid is recomputed only at T_RULES[:2], both ρ < −1; the sweep's (sa, sb) and grid are never recomputed from (τ, w); R's two-dimensional code is checked only live at K = 3 against a 10⁻¹² tolerance.
- One test-structure point: md-changing mutants die by byte comparison, not by semantic tests.

**(6) Literature** (web search only; arxiv.org, ar5iv and springer.com are blocked in this session):
- Kuo–Woźniakowski 2012 (BIT 52, 425–436), by search summaries of its abstract [snippet]:
  - an explicit formula for the Gauss–Hermite rule's error;
  - exponential convergence for 2αγ² < 1, none at 2αγ² = 1, exponential divergence beyond;
  - α is the Gaussian weight's variance, the kernel is exp(−γ²(x − y)²).
- Karvonen–Oates–Girolami (arXiv:2004.12654, Math. Comp. 2021) [snippets]: the rule converges exponentially if ℓ > α, and "the matching lower bound claimed in Kuo and Woźniakowski (2012) is erroneous as pointed out by Kuo et al. (2017)". A search summary of Karvonen–Särkkä (arXiv:1803.09532) restates the condition as "converges if and only if ε² < 1/2" [snippet]; a second, targeted search did not confirm the "if and only if" wording.
- Sriperumbudur et al. 2010 (JMLR 11): a COLT 2008 version defines the metric through the kernel's Fourier spectrum and characteristic functions [snippet]. At moment matching, (a) is exactly that identity (B = τW, W ~ N(0, ½)). For general v it adds a second Gaussian that absorbs the reweighting e^{cz²}.
- Dao–De Sa–Ré (arXiv:1709.02605) is deterministic quadrature in the frequency domain for kernel features, the opposite direction, as the draft says [snippet].
- I found no exact rate or sharp constant for any Gauss–Hermite rule's worst-case error in a Gaussian-kernel space, and no uniform asymptotics of Gauss–Hermite errors on e^{iωx} [snippets].
- Mapping: in D‴.58's variables the kernel is e^{−(y − y′)²/(4σ²)} and the weight N(0, τ²), so 2αγ² = τ²/(2σ²). D‴.62's moment-matching rate is λ = 2αγ²/(1 + 2αγ²) < 1 for every length-scale. At the threshold (τ² = 2, λ = ½) and beyond it (τ² = 8, λ = 0.8, 2αγ² = 4), my own exact errors to K = 300 give e_K²/(Cλ^{2K}) − 1 = −4.7·10⁻⁴ and −1.8·10⁻⁴, with K(·) tending to the formal s₁ (−61/432 and −0.054098) [own].
- D‴.62's own cheap bound e_K² ≲ 2^{−1/2}σ*(τ²/2)^{2K} has rate exactly 2αγ², which plausibly explains the published threshold.

## Verdicts per part

- **(a)** Correct, elementary, holds for every rule. The frame identity is confirmed [own].
- **(b)** Correct; the proof is complete, and given D‴.61 (i) and (iii)(b) on the axes it is free of computer assistance, as claimed. On the slow branch, (d) makes the o(1) visibly slow: s₁ ≈ 2t(1 + t)/(1 − t)³, so K ≫ s₁ is needed. That is not a defect.
- **(c)** The sector asymptotics, the positivity for large K (for every K when ρ < −1) and B₀ are correct. The transient is formal and labelled so.
- **(d)** Formal and labelled so. The algebra is verified, and the numbers are confirmed at 6 new rules and at the 14 stored ones. (F) is supported numerically on the imaginary axis to 8 digits.
- **(e)** Correct, apart from the cosmetic D3ea-7.
- **Status.** Accurate. The phrase "nothing new for the frame rules" is correctly restricted to Gaussian representing measures.
- **What it changes.** Correct apart from D3ea-4. The 36 = 14 + 18 + 4 count checks, and the Paper 18 quotations match PAPERS.md §18.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ea-1 | minor | **The KW threshold is not engaged.** Kuo–Woźniakowski 2012's abstract (snippets) claims convergence only for 2αγ² < 1 and divergence beyond; Karvonen–Särkkä reportedly restate it as "iff ε² < 1/2". With α = τ², γ² = 1/(4σ²), D‴.62 (b) proves the classical rule (moment matching) converges for every 2αγ², at λ = 2αγ²/(1 + 2αγ²). Own exact errors at 2αγ² = 1 and 4 confirm it [own]. (a)'s cheap bound reproduces the threshold. The divergence half is impossible anyway for positive weights summing to 1 (D‴.58's √2), and KOG report KSW 2017 found KW's lower bound erroneous [snippet]. D62's Literature, and Paper 18's proposed Novelty "moment-matching constant", rest on not having read KW's explicit error formula. | D62 Literature and "What it changes" (Paper 18 Novelty): add the mapping and the statement, conditional on checking KW 2012, KSW 2017 and KS 2019 in full: "D‴.62 (b) at ρ = 1 shows the classical rule converges at every length-scale; the threshold is the rate of the bound from \|ĝ_K(iB)\| ≤ 1." Check whether KW's explicit formula already yields λ. |
| D3ea-2 | minor | **Semantic survivors.** E20, E32, E33 and E40 would change published R and T numbers on a rerun without failing a test; E21 most likely changes nothing printed. T's grid is recomputed only for two ρ < −1 rules; the sweep's (sa, sb) and grids are never recomputed from (τ, w); R's two-dimensional values are only checked live at K = 3 against 10⁻¹². Equivalent survivors: E30, E35, E38; harmless: E36. | tests/test_frameless.py: recompute `t_rule` for one ρ > 1 rule; recompute `t_sv` from (τ, w) at a sample of sweep rules and recompute one sweep grid; run `r_2d` at small K and M against a pinned stored value to 10⁻⁴⁰. Record E30, E35, E38 as equivalent in a round-132 spec. |
| D3ea-3 | minor | **ρ = −1 is proved but untested.** No R, S or A6 rule has ρ = −1, and the 'rho=-1' paths (`params`, `C_other`) never run. Own checks: the representation at K = 6 to 2·10⁻¹⁴; C = 0.3849…, \|L − 1\| = 1.4·10⁻¹⁰, s₁ = 13/4 within 1.4·10⁻⁸ at K ≤ 300 [own]. | exp_frameless.py: add (τ/σ = 1, v = τ²/σ*²) to R (one-dimensional, rigorous) and to S; one test pinning C and s₁ there. |
| D3ea-4 | minor | **Two "What it changes" items over-state.** (i) The half-squeeze gain formula for D‴.58 (ii) is a quotient: its numerator is D‴.62 (by hand), but its denominator is D‴.61 (iv)'s cusp constant, computer-assisted for t ≥ 0.0966. The proposal omits this. (ii) D‴.61 (iv)'s new last sentence calls the \|ρ\| > 1 rules the frame picture's "degenerate cases". No frame exists there; they are the other sign case of (a), not a degeneration (only ρ → ±1 is a limit). | D62 "What it changes": add "(the half-squeeze factor by D‴.61 (iv), computer-assisted for t ≥ 0.0966)"; reword to "ρ = ±1 are the frame picture's limits; for \|ρ\| > 1 the representation (a) replaces it, on one axis's square". |
| D3ea-5 | minor | **The s₂ growth claim is unsupported.** "the fitted s₂ = −30.7, −54.2, −263 … what (d)'s (1 − t)^{−3} predicts": (d) derives s₁ only. The ratio −263/−30.7 = 8.6 does not follow (1 − t)^{−3} (ratio 1.95), and the rules also differ in ϖ. | D62 Numerics (F): "large higher coefficients on the slow branch, as the growth of s₁ suggests", or drop. |
| D3ea-6 | cosmetic | **A citation points to absent text.** "round 129's formula (D‴.58's Status)" in (d), and "D‴.58's Status" in the note: D58.md at 8ed0f31 contains no c₂ formula (round 129 is still being applied, HANDOFF). | Cite review/REVIEW_round129_D3dx.md (its c₂ formula, lines 69 and 178) until D58 is updated. |
| D3ea-7 | cosmetic | **The ratio in (e) is inverted.** "its ratio to C(v) is vτ²(1 − ρ²)^{1/2}/(2(…)^{1/2})" is C/C₆₁; C₆₁/C is the reciprocal. Both equal 1. | D62 (e) proof: "C(v)'s ratio to it". |
| D3ea-8 | cosmetic | **Two wording and caveat points.** (i) Numerics (F) says "with q^K terms where the second ridge is resolvable"; the code adds them where q^{K_min} > 10⁻²⁸. (ii) S's "largest relative arb radius" excludes the node error (midpoints taken exact); my certified recomputation shows it is harmless. Also, Step 1 (iii)'s \|φ_K(iw)\| ≤ 2 needs only positive weights, not D‴.61 (i). | D62 Numerics and Step 1 (iii): adjust the wording; add "(the nodes' midpoints taken exact; confirmed with certified nodes, round 132)". |

## Sources

- [Kuo–Woźniakowski, Gauss–Hermite quadratures for functions from Hilbert spaces with Gaussian reproducing kernels, BIT 52 (2012)](https://link.springer.com/article/10.1007/s10543-011-0358-9) and [its ACM DL record](https://dl.acm.org/doi/10.1007/s10543-011-0358-9) [snippet: abstract via search summary]
- [Karvonen–Oates–Girolami, arXiv:2004.12654](https://arxiv.org/pdf/2004.12654), [its HTML version](https://arxiv.org/html/2004.12654v2), [Math. Comp. 2021 page](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/) [snippets]
- [Kuo–Sloan–Woźniakowski 2017, preprint](https://web.maths.unsw.edu.au/~fkuo/pubs/preprint/Kuo_Sloan_Wozniakowski_2017.pdf) [snippet: record only]
- [Karvonen–Särkkä, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532) and [BIT 2019](https://link.springer.com/article/10.1007/s10543-019-00758-3) [snippets]
- [Sriperumbudur et al., JMLR 11 (2010)](https://www.jmlr.org/beta/papers/v11/sriperumbudur10a.html), [arXiv:0907.5309](https://arxiv.org/pdf/0907.5309) [snippets]
- [Dao–De Sa–Ré, arXiv:1709.02605](https://arxiv.org/pdf/1709.02605), [NeurIPS 2017](https://proceedings.neurips.cc/paper/2017/hash/62f91ce9b820a491ee78c108636db089-Abstract.html) [snippets]
- [Gnewuch–Ritter–Rüßmann, arXiv:2412.05368](https://arxiv.org/pdf/2412.05368) [snippet; nothing on exact Gauss–Hermite rates seen]
- Searches on uniform large-n asymptotics of Gauss–Hermite errors for oscillatory integrands found none (e.g. [arXiv:1806.06549](https://arxiv.org/pdf/1806.06549), [arXiv:1212.1293](https://arxiv.org/pdf/1212.1293)) [snippets]
- Repository [md, raw]: frameless.md, frameless_raw.npz, D61.md, D58.md, PAPERS.md §18, HANDOFF.md, REVIEW_round129_D3dx.md

## Lesson

The new theorem survived a step-by-step attack: its correctness rests on the one-variable axis gaps and on bookkeeping I could redo by hand.

The two real gaps are elsewhere:
- **The literature.** The theorem quietly overturns a published convergence threshold for the classical Gauss–Hermite rule and does not say so.
- **The tests.** They pin the md byte for byte but recompute the T and R computations too narrowly to catch semantic changes in them.

Branches proved but never exercised (ρ = −1) need one numerical point each.
