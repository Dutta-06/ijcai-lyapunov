# Referee report D3ed, round 135: the Paper 18 manuscript, first review

Manuscript: `manuscript/paper18/paper18.tex` (610 lines), commit 835d9a0, branch claude/vigilant-clarke-hf939b. I exported it read-only to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref135/repo`. My own scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref135/own` (exact.py, t1.py–t4.py).

## Summary of the submission

The paper studies scaled Gauss–Hermite rules Q_{K,v}. Each is the K-point Gauss–Hermite rule of N(0,v), applied to f·dμ/dN(0,v), with μ = N(0,τ²). The error is the relative worst-case error e_K(v) in the RKHS of the Gaussian kernel e^{-(y-y')²/(4σ²)}. The headline (Theorem A) has four parts:
- **Rate.** e_K(v)^{1/K} → t(v) = tanh((r+|ln(v/v_hs)|)/2), with r = ln(σ*/σ).
- **Constant off the cusp.** For v ≠ v_hs, e_K ~ C(v)^{1/2} t(v)^K, with C(v) explicit.
- **The cusp.** At v_hs, e_K ~ (1−t²)^{1/4}(πK)^{1/4} t^K ₂F₁(−½,−½;1;t²)^{1/2}.
- **Classical rule.** At moment matching, e_K ~ 2^{−1/4}(1+λ)^{3/4} λ^K for every length-scale. Per the paper, this corrects the divergence threshold reported for Kuo–Woźniakowski 2012.

The proof has five ingredients:
1. Theorem 1: the Peano kernel of order 2K, written as K!/(2K)! times the law of a Dirichlet average S. Equivalently S = sT+(1−s)Z with s ~ Beta(K,K+1).
2. A Duhamel identity for the Jacobi matrix (Lemma 4.1).
3. A corner-ratio limit (Theorem 2) and a crude bound beyond the disc (Lemma 4.2).
4. A "valley path" Laplace lemma (Lemma 4.4), giving γ̂_K(√K η)e^{−Kη²/4} → (1−η²/4)^{1/2} on |η|<2 (Theorem 3). This is by hand on |η|<1 and on the two axes, and by ball-arithmetic certificates D and VA elsewhere.
5. A Gaussian representation e_K² = σ*(v/τ²) E[E_K(a+b)E_K(a−b)] (Proposition 2.2). Laplace's method in the plane then gives the frame rules (Corollary 4), and a one-axis ridge analysis gives the rules without a frame (Theorem 5).

The rest of the paper is two appendices (A: the certificates; B: heuristics), a numerical section and several formal remarks.

## Overall assessment (recommendation: major revision)

The mathematics is, as far as I could check, correct. I re-derived every displayed identity in the main-line proofs: Theorems 1, 2, 3 (by-hand part), 5, Corollary 4, Proposition 2.2 and Lemmas 2.1, 4.1, 4.2, 4.4. I found no error and no circularity. The constants of Corollary 4 (both off the cusp and at the cusp) and of Theorem 5 at moment matching match exact errors that I computed independently. Certificates D and VA (ε₀ = 0.11) re-run at this commit and reproduce every count in Propositions A.3 and A.9.

The results are interesting and, as far as I can tell, new in their sharp form. The main candidates are the exact rate t(v) for all v, the cusp K^{1/4}, the sharp constants, and convergence of the classical rule for every 2αγ².

I still cannot recommend acceptance in the present form, for three reasons:
- **(i) The literature comparison is unfinished.** The authors say so themselves (lines 85, 87, 212–213). The closest competitor (Karvonen–Oates–Girolami 2021) proves two-sided exponential bounds for scaled Gauss–Hermite rules, with N^{±1/4} prefactors. The novelty claim must be checked against that paper's explicit constants and scaling.
- **(ii) Some analytic steps are asserted rather than proved.** These are the branch/orientation in Lemma 4.4, the uniformity of the Laplace steps in their parameters, and the main-term tail outside the disc in Corollary 4. Each is routine but should be written out.
- **(iii) The presentation is far too compressed for SINUM or Math. Comp.** The notation is heavily overloaded, several proofs are single run-on paragraphs with embedded pseudo-lemmas, and internal-process material (CHECK comments, the "Verification" acknowledgement, a citation to a repository md table) is still in the text.

The computer-assisted part is unusually well documented for a first submission. But it is not checkable from the text alone without the program, and propose and certify run in one program.

## Major comments

**M1. Novelty against Karvonen–Oates–Girolami 2021 (KOG21) and Kuo–Woźniakowski 2012 (KW12) is not established (lines 83–88, 210–214).**
- What I could see of KOG21 (abstract, via search snippets) says: two classes of rules; no constraint on the length-scale; for the scaled Gauss–Hermite class, lower and upper bounds of the forms exp(−c₁N^{1/d})N^{1/(4d)} and exp(−c₂N^{1/d})N^{−1/(4d)}, with c₁ > c₂; and the matching lower bound claimed in KW12 is erroneous, as Kuo et al. 2017 pointed out.
- So exponential convergence of a scaled Gauss–Hermite rule at every length-scale, with a two-sided bound, is known. The N^{1/4} is reminiscent of the cusp's K^{1/4}.
- The paper must identify KOG21's scaling as one of the v's (KOG21 reports observed rates r ≈ 0.528 and 0.203, per a snippet). It must compare c₁, c₂ with 2 ln(1/t(v)), and state exactly what is new: presumably the exact rate for all v, the optimal (half-squeeze) scaling, the sharp constants and the cusp law.
- The statement "converges for every 2αγ²" for the classical rule must be compared with what KSW17 and KOG21 already say ("converges exponentially if ℓ > α").
- Remark 3.7 ("conditional", "plausibly explains the published threshold") is not acceptable as written. KW12, KSW17, KS19 and KOG21 must be read in full, and the remark rewritten as a definite comparison.

**M2. Lemma 4.4(ii): branch, orientation and uniformity (lines 279–298).**
- *The branch.* The definition allows any chord direction u with Re(φ''u²) < 0, and "the same direction class" (line 293) is undefined. The sign of the Gaussian integral, and hence the principal root, is fixed only by the claim that "the orientation … is that of the straight path". For the paths actually used it is rigorous, because both the straight path and VA's real segment have u = η/|η|, which gives exactly η/(8−2η²)^{1/2} with the principal root since Re(8−2η²) > 0.
- *Fix for the branch.* Either restrict Definition 4.3 to u = η/|η|, or add the argument: the limit is analytic on Ω′, it equals ±(1−η²/4)^{1/2} locally, and Ω′ is connected and contains (0,1).
- *Uniformity.* Laplace's method is applied with a phase φ_η and a domain that depend on η, and with an amplitude carrying a uniform (1+O(1/K)) factor. Olver's theorem is for a fixed phase. Please write the three-line uniform estimate along the chord (a quadratic upper bound for Re φ_η), as is done implicitly in Theorem 2.
- *References.* BHNOD18 (arbitrary-order saddles) is not needed. The Olver theorem number is flagged unverified (line 299).

**M3. Corollary 4 (lines 317–330): the disc-to-plane step is asserted.**
- I checked the proof in full. The tails |E_K| ≤ 2e^{x²/2}, the threshold c⋆² > 4 ln ϰ/(ϰ−1) ∈ (4t, 4), the radial Γ-integral, the angular Laplace step with width O(K^{-1/2}), the cusp average (2π)^{-1}∫|1−te^{2iθ}| = ₂F₁(−½,−½;1;t²), and the identification with C(v) via 2γ(γ′−γ) = (t−t′)/(2t²t′) and (1+ϖ)^{-1/2} = t^{1/2}/(t−t′)^{1/2} are all correct.
- *The gap.* The main term is approximated on |ξ| ≤ c⋆√K only, but then integrated over the whole plane in polar coordinates. Its part beyond the disc must be shown negligible. It is: a Γ(2K+1) tail beyond c⋆²K/(2t) > 2K(1+δ), times the polynomial factor |1−ξ²/4K|. The text says only "with its mass at …".
- *Line 380.* The claim that C(v) "agrees with the constant of the proof of Corollary 4 to 5.9·10⁻²⁹⁸" refers, in frameless.md, to |C/C₆₁ − 1|, the constant of D61 (iv). D61 (iv) is the same planar Laplace formula, S_K ~ r^{2K}(1−r)/(a(2γ_x(γ_y−γ_x))^{1/2}). So the manuscript's route is not really different from the drafts' route. Since the identity is algebraic and is already written out on line 330, replace the numerical claim with the algebra.
- *Non-uniformity.* The off-cusp statement is pointwise in v and is not uniform as v → v_hs (C(v) → ∞ there). Say so in Theorem A and in Corollary 4.

**M4. Theorem 5 (lines 332–353): correct, but buried.**
- I verified every step:
  - Stirling;
  - the rates L(y) and Φ;
  - Lemma A (the maximiser at w = 4|s_b|);
  - Lemma T: R1 via Q ≥ 2/k on each sector and −2/k < 2 log(k/(1+k)); R2 via convexity of Q(·,y₂), with minimiser of modulus < 2;
  - the tilt, with N = ((1−c)(1−c_b))^{-1/2};
  - Lemma B: size-biasing Y^{4K}, (1−ϖX²/Ỹ²)^{2K} → e^{−ϖX²/2}, m_K → 1−s_d;
  - the upper-semicontinuity limit c′ ↑ 2;
  - the ρ < −1 case via c̃ = c/(1−2c), with c̃/(1+c̃) = s_a;
  - the ρ = ±1 tails.
- Presentation fixes:
  - "Lemma A/T/B" should be stated as numbered lemmas outside the proof. The names clash with Theorem A and Appendix A.
  - For ρ < −1, "verbatim" hides that Lemma A's maximiser becomes (±2|s_a|^{1/2}, 0).
  - Q is used for the quadratic form while Q_K is the rule.

**M5. The computer-assisted part (Appendix A): adequate trusted base, weak independent checkability.**
- What is good:
  - The trusted base (Arb via python-flint 0.9.0, exact double-to-ball conversion, a few binary64 steps with explicit outward rounding, Fraction for the tiling) is reasonable.
  - The stated inequalities (Lemmas A.1–A.8, A.10) are correct as far as I checked by hand: A.1, A.2, A.4, A.5(ii)–(iii), A.6, A.7 and the end piece.
  - The "propose only / certify" list is clear.
- What is missing:
  - The certify step re-runs inside the same program that proposes. A referee cannot check a certificate without re-running the proposals.
  - For Math. Comp., the authors should store the accepted boxes, polylines and connectors, with a hash, and ship a short standalone verifier that reads only the stored data and applies Lemmas A.2, A.4, A.6–A.8 and A.10.
- Margins the reader should see:
  - D's largest rigorous bound is −1.94·10⁻¹¹.
  - In VA, the inner piece A-FC's worst value is −8.8·10⁻¹⁴ (in the raw output, not in the text).
  - The handoff's min_margin is 0.0, presumably the θ = 0 edge where lo = Y = 0. Please explain this.
- The threshold 0.0966 is an artefact of using the unit disc. State this honestly.
- Remark 3.2's "D alone covers t(v) < 0.8052" holds by a slack of only 7·10⁻⁶ in c² (3.582306 vs 1.8927² = 3.582313). That is fine, but it is worth stating that this comparison is exact.

**M6. Structure and length.**
- The core (Sections 2–4) is roughly 15 dense pages. The formal Remarks 3.5–3.6 (the next order, two ridges), Appendix B and much of Section 5 (10–37-digit fit residuals) are not needed for the theorems. Move them to supplementary material, or cut them.
- Move Theorem 2 and Lemma 4.2 into the proof section as lemmas.
- Remove the "[proved]" / "[computer-assisted]" tags from headings and replace them with a sentence in the introduction.
- Remark 3.1 ("partly a sketch") should either be proved or become an open problem.

**M7. Notation is overloaded across the paper.** The worst cases are:
- μ: the measure, and the corner variable of Theorem 2;
- Q: the rules, the quadratic form of Theorem 5, and the doubles P, Q of Lemma A.6;
- λ: the classical rate, the Dirichlet weights, λ ∈ [0,1] in Lemma A.3, and λ₁(y) in Appendix B;
- t: tanh(r/2), t(v), the integration variable of I_K, t′ and t_b;
- H: the Hermite interpolant, H(ϱ), and the half-length in A.6;
- κ / ϰ: the kernel, the cumulants, VA's slope, and 1/t(v);
- s / S: the Beta variable, s_a and s_b, s₁, s = u+2q, the random variable S, the series S(w), and S_K;
- c: the variance parameter, c_K, c_η, c⋆, c′, and a disc radius;
- also M, N, A, B, C, F, G, L, p, q, u, w, ℓ.

A notation table and systematic renaming are needed.

## Minor comments

1. L36–39, L429–433: placeholders. L85, 87, 213, 299, 539 and most bibitems carry `% CHECK` / `% NOT IN MD` comments. Remove them after checking.
2. L47: the abstract is long and formula-heavy. Shorten it and drop ₂F₁ and γ̂_K from it.
3. L68: the "kernel network" motivation is unsupported. Cite it or drop it.
4. L77–79: Theorem A should say "for fixed v". The constant is not uniform near v_hs.
5. L71 and L151: λ means both tanh r and the Dirichlet weights λ_i. Rename the weights.
6. L134–138: define C_o next to C, and justify the name ϖ.
7. L141: a "centred Gaussian of variance c" with c < 0 is a formal device. Call it that.
8. L160: rename μ (it clashes with μ = N(0,τ²)).
9. L164 and L307: j = 1 also holds, since K(Var S − ½) → −¼. Either include it or say why it is excluded. Line 307 should name the disc radius explicitly.
10. L186: the "---;" punctuation; "the cross-check V" is used before it is defined.
11. L189–191 and L357: "Mercer-zero" is given only by a formula, and "the Mercer frame's rule" is never defined.
12. L211: "any rule with positive weights summing to 1" refers to Q_K, not to Q_{K,v}. Reword.
13. L220: say explicitly that Q_K f = Q_K Hf because Hf interpolates f at the nodes.
14. L245: this long inline Tonelli formula will overflow. Make it a display.
15. L253–263, L310–315, L317–330, L332–353: these proofs sit outside proof environments, with a bare `\qed`. Use `\begin{proof}[Proof of …]`.
16. L255–257: the integration variable t clashes with t. The "Laplace … uniformly" step should cite a uniform version or give the bound.
17. L293: "direction class" is undefined (see M2).
18. L326: add the Γ-tail bound (see M3).
19. L330: t and t′ are reused for t(v) and |s_o|.
20. L357: `\cite[A6]{repo}` cites a table in an md file. Put the rule list in the paper.
21. L380: replace the 1000-bit identity claim with the algebra (see M3).
22. L430: delete the "Verification … no correctness issue remained open" paragraph. It is not appropriate in a journal.
23. L466 and L517: report the tightest margins (see M5).
24. L539–540: the runtimes come from D61, not from an md. I measured 15 s (DB) and 7 s (VA011), consistent with the stated 18 s and 7 s.
25. Bibliography:
    - KOG21: Math. Comp. 90 (2021) 2209–2233 (pages missing in the bibitem).
    - KS19: pages missing.
    - BHNOD18: journal version missing.
    - KW12 (BIT 52 (2012) 425–436) and KSW17 (Math. Comp. 86 (2017) 829–853): the details match what I saw.
26. L250: good that the infinite-matrix Duhamel form is "not claimed as new". Do the same for Theorem 1's law. It is an immediate combination of the classical Gauss error formula E[He_K² f[Z, doubled nodes]] with Hermite–Genocchi; say so explicitly.

## Verification I did

**[own] Exact errors.**
- *Method.* Independent arb code (`own/exact.py`): He_K roots from flint, refined by Newton at 600–900 bits; Christoffel weights; e_K² from closed-form Gaussian inner products.
- *Frame rules (Corollary 4 off the cusp).* R = e_K²/(C(v)t(v)^{2K}).

  | τ | v | C(v) | K = 20 | K = 40 | K = 80 | md value |
  |---|---|---|---|---|---|---|
  | 1 | 1.25 v_hs | 1.12130609853 | K(R−1) = −0.16042 | −0.16069 | −0.16081 | Table 2: 160(R−1) = −0.160858, s₁ = −0.16091 |
  | 1 | 0.8 v_hs | 0.507445213475 | K(R−1) = 1.718 | 1.899 | 2.022 | — |
  | 2 | 0.5 v_hs | 0.168086336244 | R−1 = 0.820 | 0.531 | 0.329 | Table 2: 0.19480 at K = 160 (slow real branch) |

  In all three, R → 1 at rate 1/K.
- *The cusp.* Ratio e_K² / [(1−t²)^{1/2}(πK)^{1/2}t^{2K}₂F₁(−½,−½;1;t²)] minus 1:

  | τ | K = 20 | K = 40 | K = 80 |
  |---|---|---|---|
  | 1 | 5.07·10⁻³ | 2.50·10⁻³ | 1.24·10⁻³ |
  | 2 | 3.7·10⁻⁴ | −9.7·10⁻⁵ | −1.3·10⁻⁴ |

  At τ = 1 this is 1/K convergence, with K(ratio − 1) ≈ 0.10.
- *Cusp, against the Mehler series.* e_K² over the first Mehler term (1−t²)^{1/2}t^{2K}4^K(K!)²/(2K)! is 1.00518 at K = 10 and 1.00617 at K = 20, tending to ₂F₁ = 1.00738. This is consistent with Remark 3.4's positivity.
- *Moment matching, τ = 1 (Theorem 5).* R−1 = −0.0010737678 and 160(R−1) = −0.17180285 at K = 160, against Table 2's −0.0010738 and −0.171803, and s₁ = −11/64.

**[own] Closed forms.**
- Var S at K = 5120 is 0.4999511957.
- Remark 3.2's root at t = 0.0966 is 0.99982565, and 4t ln(1/t)/(1−t) increases on (0,1).
- The D threshold: 3.582306 < 1.8927². The V threshold: 3.97993 < 1.995². σ*/σ = 9.26 gives t = 0.80507; σ*/σ = 199 gives t = 0.99.
- Axis gaps: the real-axis gap G(η) is < 0 on (1,2), with its maximum −2.5·10⁻¹³ near η = 2. The imaginary-axis gap is < 0, with second derivative ≤ −1 on (1,2].
- s₁ at moment matching is −13/125, −11/64, −187/2000, −61/432 and −631/11664; at ρ = −1 it is 13/4.
- The Fock, Mercer-zero and moment-matching rate formulas, and λ/t = 1 + sech r, hold at τ = ¼, ½, 1, 2.
- Every C(v) and ρ in Tables 2–3 matches, as do C² = 27/16, 729/250 and 4/27.
- The end-piece bound is −1.33142.
- The cheap bound gives 1.2247 and 10^{361.5626} = 3.653·10³⁶¹.
- The alias zero is at 1.85891, with alias values 0.07468, 0.26713, 0.69682, 1.45024.

**[own] By hand.** Every algebraic step listed in M2–M4, plus:
- Lemma 2.1 (i)–(iii), including c·c_b = −(1−ρ²)v²/16;
- the Dyck-path count P₁ = K²+K+1, and P₂(1) = 15, P₂(2) = 52;
- the Duhamel telescoping, its Tonelli bound and the s-form exponent;
- the He_K contour shift and the |1+iζt| bound;
- the circle identity (1−cos ψ)(1−|μ|²cos(2 arg μ + ψ));
- the Lemma 4.2 limit M(a) + ½;
- Lemmas A.1, A.3 (closed forms) and A.10 (|∇Λ|² identity) of Appendix A.

**[own] Certificates re-run at this commit.**
- `--only DB`: 402664 boxes, 0 open, ball worst −1.9432·10⁻¹¹.
- `--only VA011`: 54 / 2533 / 413 (549 skipped) / 1189 (1174 skipped) boxes; outer piece 220 η-boxes and 16083 segments, worst −4.228·10⁻⁹; end piece −1.33142; handoff on 1306 sub-boxes, tiled.
- All match Propositions A.3 and A.9.

**[md]** About 110 manuscript decimals, grepped in gh_peano.md, frameless.md and width_mismatch.md: all found. I checked the context of the key ones: 0.0966, 0.8052, 0.6345, 6.1·10⁻⁵, 3.2·10⁻⁷, 5.9·10⁻²⁹⁸ (which is |C/C₆₁ − 1|, see M3), and the certificate counts.

**[snippet]** Literature, as listed below.

**[memory]** The Christoffel-number formula, Mehler's formula, and the Beta/Dirichlet moments, each re-derived.

## Literature

Only search snippets were available. arxiv.org, ams.org, link.springer.com and tskarvone.github.io did not resolve in this session; they would need to be allowed in the network settings to read the full texts.

- **KOG21**, Math. Comp. 90 (2021) 2209–2233, arXiv 2004.12654. Abstract (snippets): scaled Gauss–Hermite rules give lower and upper bounds exp(−c₁N^{1/d})N^{1/(4d)} and exp(−c₂N^{1/d})N^{−1/(4d)}, c₁ > c₂, with no length-scale restriction. A nested construction gives exp(−c₃N^{1/(2d)}). KW12's matching lower bound is described as erroneous, per KSW17. This is the main overlap (see M1).
- **KS19** (BIT 59, arXiv 1803.09532). Scaled Gauss–Hermite nodes with Mercer truncation. The snippet confirms "the Gauss–Hermite rule converges in this space if and only if ε² < 1/2". KOG21 notes that KS19 could not prove convergence of its rules.
- **KW12**, BIT 52 (2012) 425–436, DOI 10.1007/s10543-011-0358-9. I could not see the abstract.
- **KSW17**, Math. Comp. 86 (2017) 829–853 (analytic functions, Gaussian kernels, nonincreasing shape parameters).
- **Possibly missing:**
  - Al-Jarrah, Publ. Inst. Math. (Beograd) 33(47) (1983) 17–22: a Gauss–Hermite error estimate for entire functions under a growth condition.
  - arXiv 2312.07940 (2023): root-exponential convergence of Gauss–Hermite and Hermite approximation for analytic functions.
  - Recent Gaussian-kernel integration work: arXiv 2412.05368 (multi- and infinite-variate integration and L²-approximation with Gaussian kernels), and Karvonen's arXiv 2209.12473 (almost matching approximation bounds).
  - These should at least be checked.
- **Chawla–Jain 1968** (Math. Comp. 22, no. 101) is correctly cited as the saddle-point precedent.
- **Peano kernels.** Nonnegativity of the Gauss rule's Peano kernel of order 2K, the Hermite–Genocchi formula, and B-splines as Dirichlet averages (de Boor; Karlin–Micchelli–Rinott; Carlson) are classical. My searches found no prior explicit "law of S" form, but Theorem 1 is an elementary combination of known facts and should be presented as such.
- **Novelty, so far as I can judge.** The exact rate t(v) for all v, the half-squeeze optimum, the cusp K^{1/4} law with ₂F₁, the sharp constants, and convergence of the classical rule for every 2αγ² appear new. This is subject to M1.

## Issues

| id | sev. | issue | proposed fix and where |
|---|---|---|---|
| I1 | major | Literature comparison unfinished; KOG21 already has two-sided exponential bounds for scaled Gauss–Hermite rules; KW12 not read | Read KW12, KSW17, KS19, KOG21; map KOG21's scaling to v and c₁, c₂ to 2 ln(1/t(v)); rewrite §1.3 and Remark 3.7 (L83–88, L210–214) |
| I2 | major | Branch/orientation and "direction class" in Lemma 4.4 not proved; Laplace step with a parameter-dependent phase cited for a fixed phase | Restrict to u = η/|η| or add the connectivity argument; write the uniform chord estimate (L279–298) |
| I3 | major | Corollary 4: main term beyond the disc asserted; the 5.9·10⁻²⁹⁸ check misattributed (it is D61's C₆₁) and redundant | Add the Γ-tail bound (L326); replace the numerical claim with the algebra (L330, L380) |
| I4 | major | Certificates not independently checkable: propose and certify in one program; tight margins (−1.9·10⁻¹¹, A-FC −8.8·10⁻¹⁴, handoff margin 0) not reported | Store certificate data with a hash and ship a standalone verifier; report margins (App. A, L448, L466, L517) |
| I5 | major | Over-compressed presentation; pseudo-lemmas inside proofs; formal material mixed with theorems | Restructure §3–4; move Remarks 3.5–3.6, App. B and most of §5 to a supplement; state Lemmas A/T/B as numbered lemmas (L332–353) |
| I6 | major | Notation overload (μ, Q, λ, t, H, κ, s/S, c, …) | Rename systematically; add a notation table |
| I7 | minor | Theorem A / Corollary 4 do not say the constant is pointwise in v (it blows up at the cusp) | Add "for fixed v" and a sentence on non-uniformity (L77–79, L167–169) |
| I8 | minor | ρ < −1 "verbatim" hides the change of maximiser | One sentence (L349) |
| I9 | minor | Remark 3.1 is "partly a sketch" | Prove it or move it to open problems (L182) |
| I10 | minor | Internal artefacts: CHECK comments, placeholders, "Verification" acknowledgement, `\cite[A6]{repo}` | Remove or replace (L36–39, L85, 87, 213, 299, 357, 430, 539; bibitems) |
| I11 | minor | Undefined terms: "Mercer frame's rule", "direction class", Mercer-zero scaling given only by a formula | Define them (L189–191, L293, L357) |
| I12 | minor | LaTeX: long inline formula (L245), bare `\qed` outside proof environments (L263, 315, 330, 353), "---;" (L186), likely overfull tables and lines (§5, App. A) | Display formulas, use proof environments, fix punctuation, compile and fix overfull boxes |
| I13 | minor | Missing bib details (KOG21 pages 2209–2233, KS19 pages, BHNOD18 journal); unverified Olver theorem number | Complete and verify the bibliography |
| I14 | minor | j = 1 cumulant omitted; cumulant disc radius unspecified | Include K(Var S − ½) → −¼; name the radius (L164, L307) |
| I15 | minor | Theorem 1's novelty not delimited | Say it combines the classical Gauss error formula with Hermite–Genocchi (L88, L153) |
| I16 | minor | Runtimes not from an md | Keep them, but cite their source; my measurements are 15 s (DB) and 7 s (VA011) |

## Lesson

A proof can be correct and still not be a paper. Here every main-line identity checks out, by hand and against exact errors, including the planar-Laplace constant and the cusp. Yet the paper's standing depends on three things the drafts never finished: reading the closest prior work (KOG21's two-sided bounds), writing out the "routine" analytic steps (branch, uniformity, tails), and handing the reader a certificate they can check without re-running the search. In-text numerical attributions should name the check that was actually run. Here the 5.9·10⁻²⁹⁸ was |C/C₆₁ − 1|, not the manuscript's own constant, and an algebraic identity should be shown by algebra.
