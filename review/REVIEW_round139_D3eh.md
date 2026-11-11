# Referee report D3eh, round 139: the next-order term — (F) to every order and c₂ = s₁/2 (research note), first review

## Scope

I reviewed `review/RESEARCH_note_c2.md` at commit 9398242, read-only. The export is in W/repo, with W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref139`.

I read in full:
- the note;
- `review/drafts/D61.md`: (i), (ii), (iii) and the proofs of (i)–(iii);
- `review/drafts/D62.md`;
- `manuscript/paper18/paper18.tex` §§3–4: Lemma 4.1, Lemmas 4.2–4.9, Definition 4.5, Theorems 1–4 and Corollary 3;
- `review/REVIEW_round129_D3dx.md` §5.

I also read the note's scripts in `review/research_c2/` (lib.py, laplace.py, check_exact.py, check_ridge.py and their .out files). I did not run them.

All my own code is in `W/own/`:
- `alg.py`: exact-rational algebra;
- `exact.py` and `run_cg.py`: c_K and ĝ_K without nodes;
- `nu0.py`: c_K near ν = 0;
- `eK.py` and `run_e.py`: e_K without nodes and without A6;
- `k4.py`: exact cumulants from my own path DP.

Each run used `OMP_NUM_THREADS=1 nice -n 15`, one heavy process at a time.

## Summary

**Verdict: the proof stands.** I tried to break Lemma L, Proposition A (i)–(iii), Proposition B, Corollary B′ and the correction to D‴.62 (d), and found no mathematical error.
- I re-derived every hand computation independently: the Wick formulas (W1) and (W2), b₁, the six-term γ₁ table, G₁ = F₁, the Hessian, the third and fourth derivatives and the log-amplitude at the ridge, T1–T5 on both axes, the Stirling constants, the prefactor C(v), and c₂ and c₃ from log S.
- My own exact numerics use methods that share nothing with the note's node-based code. They confirm γ₁, F₁, F₂ (as identified), s₁, the predicted s₂, c₂ = s₁/2 and c₃ to 10⁻¹⁵ … 10⁻²³ on the fast branches. On the slow branch they agree to within its known pre-asymptotic accuracy.

**Recommendation: accept, with minor revisions, before D‴.62 (d), D‴.58 or Paper 18 cite it.**

Findings:
- **D3eh-1 (minor).** Several cross-references to Paper 18 are off by one.
  - "Lemma 4.5 (i)/(ii)" should be Lemma 4.6 (the path lemma). Definition 4.5 is the valley condition; the lemma and definition counters are shared.
  - "Theorem 3" should be Theorem 2 (the transform's limit). Corollary 3 shares the theorem counter.
  - These appear in Proposition A (iii) and in §4 ("Lemma 4.5 (ii) and Theorem 3 gain the expansion").
- **D3eh-2 (minor).** Lemma L is classical and is presented without literature.
  - Laplace's method to every order with a uniform remainder, and coefficients by Wick/Feynman-diagram bookkeeping, are standard. (W2) has the structure of the classical 1/n correction term.
  - Proposition A (i) is the exterior Plancherel–Rotach expansion. Its coefficients and error bounds have been published (Shi–Nemes–Wang–Wong).
  - Lemma L should be stated as "standard, written out with the parameter explicit", with citations.
- **D3eh-3 (minor).** The correction to D‴.62 (d) is right, but should be phrased precisely.
  - (F) with remainder O(K⁻²), plus a 2-D Laplace remainder bound, gives only S(K) = 1 + s₁/K + O(K⁻²). That yields "no K⁻¹ term" but not c₂.
  - The minimal sufficient hypothesis is (F) with an F₂ term and remainder O(K⁻³), which is n = 3 of Proposition A.
  - D‴.62 (d) was labelled formal, and its last line already assumed the s₂ term. So "the statement's hypothesis is too weak" is more accurate than "wrong".
- **D3eh-4 (cosmetic).** Lemma L's statement should say three things:
  - for d ≥ 2, h_p and a_p may be complex as long as H_p is real positive definite;
  - L_m uses only Taylor degrees ≤ 2m+2 of h and ≤ 2m of a, so it does not depend on n;
  - (W2) assumes diagonal G.
- **D3eh-5 (cosmetic).** The step at ν = 0 in Proposition A (ii) is correct.
  - It should say that equality on 0 < |ν| ≤ r plus continuity at ν = 0 is enough.
  - It should also say that K₀ depends on r, through I_K being zero-free on |ζ| ≤ r.
- **D3eh-6 (minor).** "Now explained" for D‴.62's fitted slow-branch s₂ overstates the precision.
  - The prediction −263.0676 differs from D‴.62's fit −263.05374 by 1.4·10⁻², and from my extraction (K ≤ 360) −263.018 by 5·10⁻².
  - That is consistent, given the very large higher coefficients there, but agreement should be quoted as about 10⁻⁴ relative on the slow branch.
- **D3eh-7 (cosmetic).** Reproducibility details.
  - The Files section names `W/NOTE.md` and `W/scripts/`; the repository copy is `review/research_c2/`.
  - `check_ridge.py` hard-codes REPO = `/home/user/…` and imports `experiments/exp_frameless`.

## What I checked and how

**1. Lemma L** [own, by hand]
- **The bookkeeping.** A monomial K^k x^ν integrates to (2π/K)^{d/2} det^{−1/2} K^{k−|ν|/2} E[X^ν]. Since |ν| ≥ 3k, the K^{−m} term has k ≤ 2m and |ν| = 2(k+m). The highest single-factor degree needed is 2m+2 (for k = 1), so the Taylor truncations 2n+2 (of h) and 2n (of a) are more than enough.
- **The window K^{−1/3}.** There |KT| ≤ CK|x|³ ≤ C, so e^{KT} is bounded. Each remainder term integrates to O(K^{−d/2−n−1/2}):
  - |x|^{2n+1} gives K^{−n−1/2−d/2};
  - K|x|^{2n+3} gives K^{−n−1/2−d/2};
  - K^{2n+1}|x|^{6n+3} gives K^{−n−1/2−d/2}.
- **The tails.** Outside the window the integrand is O(e^{−δK^{1/3}}). Extending the Gaussian integrals to ℝ^d needs Re H ≥ 2δ (or H ≥ 2δI), which is assumed. Odd monomials vanish on the symmetric window.
- **(W1).** E X² = 1/(2α), E X⁴ = 3/(4α²), E X⁶ = 15/(8α³) give (W1) exactly.
- **(W2).** I recounted the Wick pairings: 3 pairings for the a₁·h₃ and h₄ terms; 9 + 6 of the 15 six-point pairings for h₃². This gives the coefficients ½, ⅛, ⅛ and 1/12. Correct.
- **Hypotheses where the lemma is applied:**
  - (i): α = (1−ζ²)/2, and Re ≤ −(1−r²)x²/2 holds on all of ℝ. I checked the bound ½log|1+iζx|² ≤ −Im ζ·x + |ζ|²x²/2.
  - (ii): Re φ_ν ≤ −(2/π²)(1−r²)θ² on [−π, π], and Re H = 1 − Re ν² ≥ 4δ.
  - (iii), the chord: Re χ_η ≤ −(1 − Re y)z².
  - Proposition B: a real Hessian diag((1+ϖ)/(ϖb²), 2/b²) at (0, ±2b), with no cross term (Ψ_pq is odd in p).
  - The Cauchy-estimate derivative bounds hold in every case, since the data are analytic on fixed neighbourhoods.
- The hypotheses are met in every application.

**2. I_K and c_K** [own, by hand; exact rationals]
- **h₃ and h₄.** I re-expanded log(1+iζx), which gives h₃ = −iζ³/3 and h₄ = −ζ⁴/4, and hence b₁.
- **The θ-integral.** I re-derived it from Paper 18's (4.2) with ζ = νe^{iθ}, the circle clockwise.
  - The prefactor is (K−1)!K^{1−K}e^K·ν, and the amplitude is (e^{−iθ} − ν²e^{iθ})/(2π·I_K).
  - φ_ν's Taylor coefficients −(1−u)/2, i(1+3u)/6 and (1−7u)/24 check out.
  - At ν = 0 the integral is K!K^{−K}·[w^{K−1}]e^{Kw} = 1 = c_K(0).
- **The amplitude derivatives.** By my own differentiation:
  - A₁/A₀ = −i(1+2u)q;
  - A₂/A₀ = −½ + uq − uq² − u²q²/2, which equals the note's form.
- **The six-term table.** Each entry and the sum u + u² are exact at 40 random rationals (`alg.py`).
- **Evenness.** c_K is even (T is symmetric), so γ_m is analytic in u = ν² on |u| < 1.
- **The maximum-principle step.** It is indeed unnecessary.
  - Both sides are analytic, since I_K(ζ)ζ^{−K}·e^{Kζ²/2}… is built from the polynomial He_K(√K(ζ+1/ζ))ζ^K.
  - They agree on 0 < |ν| ≤ r, so continuity at ν = 0 closes the gap.

**3. Passage to (F)** [own, by hand]
- **The chord data.** In Paper 18's (4.4)/(4.3) frame:
  - phase −(1−y)z² + h with h = log(1−z²) + z², so h₃ = 0 and h₄ = −½;
  - amplitude a = 1/(1+z) − y(1+z), so a₀ = 1−y, a₁ = −(1+y), a₂ = 1;
  - Stirling for K·C(2K,K)4^{−K}: (K/π)^{1/2}(1 − 1/(8K) + …).
- **G₁.** These give G₁ = γ₁ + 1/(8(1−y)²) − 1/8, which is round 129's F₁ (exact at 40 rationals).
- **Uniformity.** It is local in η, by Paper 18's scaled paths C_η = (η/η₀)C₀ with the sub-chord in a fixed compact subset of the open unit disc, and then by compactness.
- **Off the sub-chord.** Lemma 4.4, and Lemma 4.3 near ν = 0, make that part exponentially small. Leading-order bounds suffice there, so the "unchanged" claim is right.
- **The branch.** ĝ_K is real and positive on both axes for K ≥ K₀; elsewhere the principal log of the bracket, which is near 1, is a legitimate definition.

**4. Proposition B** [own, by hand; exact rationals]
- **The change of variables.** I redid it: (K/2π)(|c||c_b|)^{−1/2}, ½ + 1/(2|c|) = 1/(2|s_a|), and Ψ_max = 2log(4t) − 2.
- **The prefactor.** K^{2K}(K!/(2K)!)² = ½(e/4)^{2K}exp(1/(12K) + O(K⁻³)). With the two symmetric points this gives ½·2·t^{2K}. The leading constant is exactly C(v):
  - det H = 2(1+ϖ)/(ϖb⁴);
  - ϖb⁴/(|c||c_b|) = 1/((1−c)(1−c_b)).
- **The derivatives at (0, 2b).**
  - h_qqq = h_ppq = 1/b³, and h_pppp = h_ppqq = h_qqqq = −3/(2b⁴), from 2log(p'−p) + 2log(p'+p).
  - ℓ_q = b/(1+b²) and ℓ_pp = ℓ_qq = (1−b²)/(2(1+b²)²).
  - On the real axis: ℓ_q = −b/(1−b²) and ℓ″ = −(1+b²)/(2(1−b²)²).
- **The sum.** T1 + T2 + (T3+T4+T5) + 2F₁(s) + 1/12 = D‴.62's s₁, on both axes, at 60 random (s, ϖ). By hand, the k-terms collect to −(k/8)(1+s)²/(1−s)².
- **ρ = ±1.** The 1-D phase 4log p' − p'²/(2λ) has the same h_qqq and the same amplitude, so s₁ is the k = 0 value. The moment-matching s₁ values −13/125, −11/64 and −187/2000 reproduce.
- **Steps 1–2.** They are sound to every order:
  - off the box: D‴.62's Lemma T margin;
  - inside the box away from U_±: an upper semicontinuous Ψ on a compact closure, with unique maximisers (Lemma 4.7).
  - The second ridge lies in the Step 2 region, so J_o = O(q^K) is beyond every power of 1/K.
- **Corollary B′.** I checked it by exact series in x = 1/K: the ratio series is 0·x + (s₁/2)x² + ((s₁ + 2s₂ − s₁²)/2)x³.
  - For c₂ with an O(K⁻³) remainder, S needs n = 3; for c₃, n = 4.

**5. The correction to D‴.62 (d)** [own]
- With S = 1 + s₁/K + ε(K) and ε = O(K⁻²), the term ε(K−1) − ε(K) is O(K⁻²) in general, so c₂ is undetermined. The correction is valid (D3eh-3).

**6. Independent numerics** [own, python-flint balls]
- **Methods.** None of them uses nodes or the A6 closed form.
  - **Golub–Welsch functional calculus.** Q_K[f p_j] = f(J_K)_{j0}, with e^{xJ_K}e₀ by a Taylor series in balls:
    - c_K(u)·u^{K−1}e^{u²/2} = √(K−1)!·(e^{uJ_K})_{K−1,0};
    - φ_K(ξ) = e^{ξ²/2} − (e^{ξJ_K})₀₀.
  - **A new series for e_K.** e_K² = σ*(v/τ²)·Σ_{k even}(v/2)^k/k!·d_k², with d_k = (E−Q_K)[z^k e^{αz²}] and α = c − v/4. It is absolutely convergent because c < ½. I validated it against the stored exact e_K: all 30 compared digits agree at K = 10, 20, 40 for (½, MM), (½, 2v_hs), (1, 0.5v_hs) and (1, ρ = −1).
- **c_K/(1−ν²) at K = 40 … 640, Neville in 1/K.** Agreement with γ₁ (in brackets, my γ₂ against the note's engine):

  | ν | γ₁ agreement |
  |---|---|
  | 0.5i | 1.1·10⁻¹⁵ |
  | 0.9i | 6.3·10⁻¹⁸ |
  | −0.4+0.55i | 3.8·10⁻⁹ (γ₂ = −0.6660432+0.0648864i, engine −0.6660431+0.0648863i) |
  | 0.5 | 1.8·10⁻⁷ |
  | 0.05, 0.05i, 0.03+0.04i (near the ν = 0 claim) | ≤ 1.0·10⁻¹⁵ |

- **ĝ_K at K = 40 … 640.** Extrapolated F₁ and F₂ against the closed forms:

  | η | F₁ | F₂ |
  |---|---|---|
  | i | 1.8·10⁻¹⁸ | 1.7·10⁻¹⁷ |
  | 1.5i | 8.7·10⁻¹⁹ | 7.3·10⁻¹⁸ |
  | 1.9i | 5.4·10⁻¹⁹ | 1.4·10⁻¹⁸ |
  | 0.6+0.5i | 2.1·10⁻¹³ | 4.3·10⁻¹² |
  | 0.3+1.2i | 2.1·10⁻¹⁶ | 2.3·10⁻¹⁵ |
  | 1 (real axis, slow) | 9.7·10⁻⁹ | 6.1·10⁻⁷ |

- **Exact cumulants (`k4.py`).** My own monic path DP confirms P₁ and Paper 18's P₂ at K = 3, 5, 8. The 1/K² coefficients of κ₂K and κ₄K² give the y and y² coefficients of F₂ as −13/8 and −37/2, matching the closed form. The y² coefficient of F₁ is 35/8.
- **e_K at consecutive K, 9 points from 120 to 360.**

  | rule | s₁ extrapolated − closed | s₂ against the note's prediction | c₂ − s₁/2 | c₃ − (s₁+2s₂−s₁²)/2 |
  |---|---|---|---|---|
  | (½, MM) | 1.1·10⁻²⁰ | 0.0839040000 exact to print | 5.0·10⁻²¹ | 3.4·10⁻¹⁹ |
  | (½, 2v_hs), ρ>1, ϖ = 0.74 (tests the k-terms) | 6.2·10⁻²² | 2.4·10⁻¹¹ | 3.1·10⁻²¹ | 4·10⁻²¹ |
  | (¼, 1.25v_hs), ρ>1 | 5.9·10⁻²² | 1.2·10⁻¹¹ | 3.7·10⁻²³ | 6·10⁻²¹ |
  | (1, ρ = −1), slow | 2.2·10⁻⁷ | 2.3·10⁻⁵ against −507/16 | 7.9·10⁻⁷ | — |
  | (1, 0.5v_hs), slow | 2.1·10⁻⁴ | 4.9·10⁻² against −263.0676 | 6.1·10⁻⁴ | — |

  The two slow-branch rows are limited by their large higher coefficients, as the note says.

**7. Literature** [snippet; arXiv abstract pages read]
- I read the arXiv abstracts (via curl) of Kirwin (2010, arXiv:0810.1700), Nemes (2013, arXiv:1207.5222), Shi–Nemes–Wang–Wong (arXiv:2109.01054), Ogden (arXiv:1808.06341) and Tang–Reid (arXiv:2107.10885).
- WebFetch could not reach arxiv.org or dlmf.nist.gov, so DLMF 18.15 and 12.10 are [memory].
- No source I saw gives γ₁, F₁, F₂ or s₁ for this problem.

## Verdicts per part

| part | verdict |
|---|---|
| Lemma L (all orders, uniform, d ≥ 1, the Wick formulas) | Correct. Classical; cite it (D3eh-2, D3eh-4) |
| A (i) I_K, b₁ | Correct; matches the Plancherel–Rotach exterior expansion |
| A (ii) c_K, γ₁ and the six-term table, uniform on \|ν\| ≤ r, ν = 0 regular, no maximum principle | Correct (D3eh-5 cosmetic) |
| A (iii) (F) to every order on Ω′ ⊃ {\|η\|<1} ∪ both axes, F₁ = round 129's | Correct (D3eh-1 references) |
| A (iv) F₂ | Correctly labelled numerical. Confirmed numerically at 6 points and its first two Taylor coefficients exactly |
| B: off-box tails, the second ridge, 2-D Laplace, T1–T5 sum to s₁, both axes, ρ = ±1 | Correct |
| B′: c₂ = s₁/2, c₃ | Correct; n = 3 (resp. 4) suffices |
| Correction to D‴.62 (d) | Correct; rephrase (D3eh-3) |
| Remark F (frames) | A sketch, correctly labelled; the formal p → ip continuation is consistent. Not reviewed as a proof |
| Numerics claims | Reproduced independently; the slow-branch wording should be softened (D3eh-6) |

## Issues

| id | sev. | issue and evidence | proposed fix |
|---|---|---|---|
| D3eh-1 | minor | Paper 18 cross-references are off by one. The compiled PDF numbers Definition 4.5 (valley), Lemma 4.6 (path), Theorem 2 (limit) and Corollary 3 [read] | Write Lemma 4.6 (i)/(ii) and Theorem 2 in A (iii) and §4 |
| D3eh-2 | minor | Lemma L and A (i) are standard, but no literature is cited [snippet] | Cite Olver ch. 3, Wong ch. IX, Bleistein–Handelsman ch. 8, Kirwin 2010, Nemes 2013, Shun–McCullagh 1995, and Shi–Nemes–Wang–Wong for A (i) |
| D3eh-3 | minor | The D‴.62 (d) correction reads as "wrong", while (d) was labelled formal and implicitly assumed the s₂ term | State: "(F) with O(K⁻²) gives only 'no K⁻¹ term'; c₂ needs (F) with F₂ and O(K⁻³), i.e. Proposition A at n = 3" |
| D3eh-4 | cosmetic | Lemma L does not say h and a may be complex for d ≥ 2, that L_m does not depend on n, or that (W2) assumes diagonal G | Add one sentence |
| D3eh-5 | cosmetic | At ν = 0: identity theorem or continuity, and K₀(r) not stated | Add one clause |
| D3eh-6 | minor | Slow-branch s₂ "now explained": −263.0676 against the fit −263.05374 (1.4·10⁻²) and my −263.018 [own, repo] | Quote agreement of about 10⁻⁴ relative on the slow branch |
| D3eh-7 | cosmetic | The Files section's W/ paths differ from the repository's, and REPO is hard-coded in check_ridge.py | Point to `review/research_c2/` and parametrise REPO |

## Sources

- [read] D61.md, D62.md, paper18.tex and paper18.pdf (numbering via pdftotext), REVIEW_round129_D3dx.md §5, the note and its scripts.
- [repo] `results/chaos/frameless_raw.npz`, used only to validate my e_K implementation.
- [snippet, arXiv abstract pages read]:
  - [Kirwin, Higher Asymptotics of Laplace's Approximation, Asympt. Anal. 70 (2010), arXiv:0810.1700](https://arxiv.org/abs/0810.1700);
  - [Nemes, An explicit formula for the coefficients in Laplace's method, Constr. Approx. 38 (2013), arXiv:1207.5222](https://arxiv.org/abs/1207.5222);
  - [Shi–Nemes–Wang–Wong, Error bounds for the asymptotic expansions of the Hermite polynomials, Proc. RSE A, arXiv:2109.01054](https://arxiv.org/abs/2109.01054): Plancherel–Rotach expansions with recursive coefficients;
  - [Ogden, arXiv:1808.06341](https://arxiv.org/abs/1808.06341);
  - [Tang–Reid, arXiv:2107.10885](https://arxiv.org/abs/2107.10885).
- [snippet] Shun–McCullagh (JRSS B 1995), known only through the secondary papers above; [Daniels 1954 (Ann. Math. Statist.)](https://www.projecteuclid.org/euclid.aoms/1177728652) and uniform saddlepoint work (Adv. Appl. Prob. 1988), for Edgeworth/saddlepoint context.
- [memory] Olver, Asymptotics and Special Functions, ch. 3–4; Wong, Asymptotic Approximations of Integrals, ch. IX; Bleistein–Handelsman, ch. 8; DLMF §18.15 and §12.10. WebFetch could not reach arxiv.org or dlmf.nist.gov.

## Lesson

A proof that redoes classical Laplace bookkeeping is best attacked at its hand algebra and its hypotheses-at-application, not at the lemma. Everything here held up, including the one non-classical step: removing the maximum principle via the θ-integral at ν = 0.

The strongest independent test was a node-free exact e_K: a positive series of squared defects d_k = (E−Q_K)[z^k e^{αz²}]. It reproduced the predicted s₂ and c₃ to 10⁻¹¹ … 10⁻²¹. That test exercises F₂ and the order-2 ridge Laplace at once.

The remaining weaknesses are presentation:
- cross-references off by one;
- missing classical citations;
- the exact hypothesis that c₂ needs: (F) through F₂ with remainder O(K⁻³).
