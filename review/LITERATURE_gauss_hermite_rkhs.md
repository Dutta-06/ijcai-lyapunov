# Literature notes: Gauss–Hermite rules in Gaussian-kernel spaces (for Paper 18)

These are reading notes for the comparison in `manuscript/paper18/paper18.tex` §1.3 and Remark 3.4 (`rem:kw`), and for issue I1 of `review/REVIEW_round135_D3ed.md`. They were written on 2026-11-02, after round 135 (not a referee round).

The arXiv texts were downloaded into the session scratch directory and read in full from their LaTeX source (KOG21, KS19) or from `pdftotext` output (the others). Nothing from them was executed. Section, equation, figure and page numbers below refer to the arXiv versions named.

Throughout, the notation is the manuscript's:
- kernel κ(y, y′) = e^{−(y−y′)²/(4σ²)}, measure μ = N(0, τ²), σ*² = σ² + τ², r = ln(σ*/σ);
- the scaled rule Q_{K,v} is the K-point Gauss–Hermite rule of N(0, v) applied to f·dμ/dN(0, v);
- e_K(v) is the relative worst-case error; the absolute one is wce = (σ/σ*)^{1/2} e_K;
- λ = tanh r = τ²/(2σ² + τ²), and t(v) = tanh((r + |ln(v/v_hs)|)/2) with v_hs = τ²σ/σ*;
- ρ(v) = 1 + 2σ²/τ² − 2σ²/v is the Mehler-frame parameter.

## Status of each source

| Source | Read? | Where from |
|---|---|---|
| Karvonen–Oates–Girolami, *Integration in reproducing kernel Hilbert spaces of Gaussian kernels*, Math. Comp. 90 (2021), no. 331, 2209–2233; arXiv:2004.12654 **v2** (31 Mar 2021, 24 pp.) | **in full** | arXiv source tarball (`article.tex`) and PDF |
| Karvonen–Särkkä, *Gaussian kernel quadrature at scaled Gauss–Hermite nodes*, BIT 59 (2019), no. 4, 877–902; arXiv:1803.09532 **v4** (2 May 2019, 25 pp.) | **in full** | arXiv source tarball and PDF |
| Kuo–Woźniakowski, *Gauss–Hermite quadratures for functions from Hilbert spaces with Gaussian reproducing kernels*, BIT 52 (2012), no. 2, 425–436, DOI 10.1007/s10543-011-0358-9 | **not read** | not on arXiv; link.springer.com, UNSW pages (web.maths.unsw.edu.au) not reachable |
| Kuo–Sloan–Woźniakowski, *Multivariate integration for analytic functions with Gaussian kernels*, Math. Comp. 86 (2017), no. 304, 829–853 | **not read** | not on arXiv; ams.org, UNSW not reachable |
| Gnewuch–Ritter–Rüßmann, *Multi- and infinite-variate integration and L²-approximation on Hilbert spaces with Gaussian kernels*, arXiv:2412.05368 **v2** (5 Dec 2025) | in full | arXiv PDF |
| Wang–Zhang, *Convergence analysis of Hermite approximations for analytic functions*, arXiv:2312.07940 **v3** (11 Jun 2025) | in full (skimmed for overlap) | arXiv PDF |
| Karvonen–Suzuki, *Approximation in Hilbert spaces of the Gaussian and related analytic kernels*, arXiv:2209.12473 **v4** (2 Oct 2025) | in full (skimmed for overlap) | arXiv PDF |
| Al-Jarrah, Publ. Inst. Math. (Beograd) 33(47) (1983) 17–22 | **not read** | elib.mi.sanu.ac.rs not reachable |

The arXiv search (by title and by the authors Woźniakowski and F. Y. Kuo) finds neither KW12 nor KSW17. Everything said below about these two papers comes from what KOG21, KS19 and GRR24 report.

Bibliographic details were cross-checked against the reference lists of KOG21, GRR24, Karvonen–Suzuki and Wang–Zhang. They confirm:
- KW12: BIT 52(2), 425–436;
- KSW17: Math. Comp. 86(304), 829–853;
- KS19: BIT 59(4), 877–902;
- KOG21: 2209–2233 (Karvonen–Suzuki; GRR24 prints 2209–2223, a typo);
- Irrgeher–Kritzer–Leobacher–Pillichshammer, J. Complexity 31(3) (2015), 380–404;
- Fasshauer–McCourt, SIAM J. Sci. Comput. 34(2) (2012), A737–A762;
- Barrett, Comput. J. 3(4) (1961), 272–277;
- Paulsen–Raghupathi (2016), Cambridge Stud. Adv. Math. 152.

## 1. Karvonen–Oates–Girolami 2021 (KOG21)

### Setting

- §1, eq. (1.1)–(1.2):
  - the weight is N(0, α²), **α a standard deviation**;
  - the kernel is K_ℓ(x, y) = exp(−(x−y)²/(2ℓ²));
  - the worst-case error is absolute.
- In our notation: α = τ and ℓ² = 2σ².
- Hence α²/ℓ² = τ²/(2σ²) (= 2αγ² of KW12) and α²/(α²+ℓ²) = λ.
- The norm of the integration representer is (1 + 2α²/ℓ²)^{−1/4} = (σ/σ*)^{1/2} (Lemma 2.8). This is our wce(0).

### The rule (§2.2, p. 5, eqs. (2.3)–(2.5))

The rule is Q^GH_{α,ℓ,n}(f) = (β/α)·Q^GH_{β,n}(f·e^{x²/(2ℓ²)}), with β² = α²ℓ²/(α²+ℓ²). Here Q^GH_{β,n} is the n-point Gauss–Hermite rule of N(0, β²).

Since 1/β² − 1/α² = 1/ℓ², we have (β/α)e^{x²/(2ℓ²)} = dN(0,α²)/dN(0,β²). So **KOG21's rule is exactly Q_{K,v} with v = β² = 2σ²τ²/(2σ²+τ²)**. In the manuscript's units (σ = 1) this is v = 2τ²/(1+σ*²), the **Fock scaling** of Remark 3.2:
- ρ(v) = 0, the Taylor frame;
- u_v = ln(v/v_hs) = −ln cosh r;
- t(v) = λ/(2−λ).

KOG21 build it as the generalised Gaussian rule (they cite Barrow 1978) for the orthonormal basis φ_{ℓ,m} = x^m e^{−x²/(2ℓ²)}/(ℓ^m √m!) of H (eq. (1.7)), m < 2n. This basis is the ρ = 0 Mehler frame.

### The bounds (Theorem 2.5, eq. (2.9), p. 8)

C_n (ℓ/√(α²+ℓ²)) (α²/(2(α²+ℓ²)))^n n^{1/4} ≤ e(Q^GH_{α,ℓ,n}) < π^{−1/4} (ℓ/√(α²+ℓ²)) (α²/(α²+ℓ²))^n n^{−1/4},

where:
- C_n = 2^n n!/√((2n)!) · n^{−1/4}, strictly decreasing to π^{1/4} (Lemma 2.2, eq. (2.8), p. 6);
- ℓ/√(α²+ℓ²) = β/α.

In our notation: lower bound ∝ (λ/2)^K K^{1/4}, upper bound ∝ λ^K K^{−1/4}.

Where each factor comes from:
- *Lower bound.* It is the error of the single basis function φ_{ℓ,2n} (Lemma 2.3, p. 7), from the Gauss–Hermite remainder in mean-value form, eq. (2.2): f^{(2n)}(ξ)α^{2n}n!/(2n)!.
- *Upper bound.* It sums the termwise bound |I(φ_{ℓ,2q}) − Q(φ_{ℓ,2q})| ≤ C_q^{−1}(β/α)(β²/ℓ²)^q q^{−1/4} (Lemma 2.4, p. 7), which uses positivity of the weights only.
- *The n^{±1/4}.* Both are the Stirling factor 2^n n!/√((2n)!) ~ (πn)^{1/4}.

### Their constants c₁, c₂ (abstract; Remark 2.12, p. 13)

For d = 1, c₁ = ln(2/λ) and c₂ = ln(1/λ); the d-variate tensor product has N^{1/d} in place of N.

Our Theorem A gives e_K(v) = C(v)^{1/2}(λ/(2−λ))^K(1+o(1)) for this v (a frame rule off the cusp). So:
- the exact constant is ln((2−λ)/λ), strictly between c₂ and c₁;
- the "2 ln(1/t(v))" of the question is the e_K² exponent (KOG21's c's refer to wce, not wce²);
- no power of K occurs: the true K^{1/4} appears only at v_hs.

Numerical comparison:
- Exact errors are from `experiments/exp_frameless.py`'s `err2`, imported read-only. The script is `scratchpad/lit/scripts/kog_cmp.py`.
- "lo/e" and "hi/e" are KOG21's bounds (in relative form, divided by wce(0)) over the exact e_K.

| (α, ℓ) of KOG21 | τ/σ | λ | t(v) = λ/(2−λ) | C(v) | c₁ = ln(2/λ) | c₂ = ln(1/λ) | ln(1/t) | 2 ln(1/t) |
|---|---|---|---|---|---|---|---|---|
| (1.5, 1), Fig. 1 left | 2.12132 | 9/13 = 0.692308 | 9/17 = 0.529412 | 0.356892 | 1.060872 | 0.367725 | 0.635989 | 1.271978 |
| (1, 1.4), Fig. 1 right | 1.01015 | 25/74 = 0.337838 | 25/123 = 0.203252 | 1.000723 | 1.778336 | 1.085189 | 1.593309 | 3.186617 |
| (1, √2) | 1 | 1/3 | 0.2 | 1.011929 | 1.791759 | 1.098612 | 1.609438 | 3.218876 |
| (2, √2) | 2 | 2/3 | 0.5 | 0.395285 | 1.098612 | 0.405465 | 0.693147 | 1.386294 |

| case | K | e_K | e_K/(C^{1/2}t^K) | lo/e | hi/e | −ln(e_K)/K |
|---|---|---|---|---|---|---|
| (1.5, 1) | 10 | 1.311346e-03 | 1.269150 | 3.8e-02 | 6.9e+00 | 0.66367 |
| (1.5, 1) | 20 | 2.118901e-06 | 1.185688 | 6.9e-04 | 9.1e+01 | 0.65323 |
| (1.5, 1) | 40 | 5.992884e-12 | 1.121048 | 1.8e-07 | 1.7e+04 | 0.64601 |
| (1, 1.4) | 10 | 1.288485e-07 | 1.070467 | 3.4e-01 | 6.2e+01 | 1.58646 |
| (1, 1.4) | 20 | 1.503247e-14 | 1.037948 | 6.5e-02 | 8.6e+03 | 1.59143 |
| (1, 1.4) | 40 | 2.139448e-28 | 1.020349 | 1.9e-03 | 1.9e+08 | 1.59280 |
| τ = 1 | 40 | 1.128257e-28 | 1.020078 | 2.2e-03 | 2.1e+08 | 1.60879 |
| τ = 2 | 40 | 6.300306e-13 | 1.101809 | 3.8e-07 | 3.7e+04 | 0.70233 |

- The ratio e_K/(C^{1/2}t^K) tends to 1 at rate 1/K, as Corollary 3 states.
- C(v) at τ = 1 and τ = 2 equals frameless.md's Fock row: 1.01192885125 and 0.395284707521.
- KOG21's Fig. 1 caption (p. 9) reports the empirical rates r ≈ 0.528 (left) and r ≈ 0.203 (right). These are t(v) = 9/17 = 0.5294 and 25/123 = 0.2033, within the accuracy of a finite-n fit; our −ln(e_K)/K at K = 40 corresponds to 0.524 and 0.2033.

### Exact rate, optimal scaling, constant

KOG21 have none of these:
- Remark 2.6 (p. 8) attributes the exponential gap to the rough Lemma 2.4 and expects improvement.
- §4 (Conclusions, p. 21), first bullet, lists this as open.
- The scaling is chosen by the exactness requirement, not optimised; for τ > 0 it is never v_hs.

### The classical rule in KOG21

- *§1, p. 2.* KW12 and KSW17 "consider the n-point Gauss–Hermite rule … and prove that e(Q^GH_{α,n}) ≤ b_n(α²/ℓ²)^n with (b_n) a decreasing sequence converging to 2^{−1/4}. That is, the Gauss–Hermite rule converges with an exponential rate if ℓ > α."
- *Footnote 1, p. 2.* "Note that the matching lower bound claimed in Kuo and Woźniakowski (2012) is erroneous as pointed out by Kuo et al. (2017, p. 830)."
- *Continuing, p. 2.* "Their potential non-convergence when ℓ ≤ α is perhaps not surprising because these rules are not adapted to the RKHS". So the case ℓ ≤ α, i.e. 2αγ² ≥ 1, is left open.
- *§2.3, p. 9.* KW12 and KSW17 "analyse integration in H(K_ℓ) under the constraint that ℓ > α".
- *KSW17's lower bound for the minimal error.* KOG21 note that KSW17 Theorem 4.1 does not need ℓ > α. It is super-exponentially small (eq. (2.12)–(2.13)), and gives Theorem 2.7.
- *Fig. 1 caption, p. 9.* "Note that the earlier results for the standard Gauss–Hermite rule in [KW12] and [KSW17] do not apply in the left panel because α > ℓ … Interestingly, the upper bound in (2.9) appears to be almost exact for the *standard* Gauss–Hermite rule." The rate of (2.9) is λ. So **KOG21 observed numerically that the classical rule decays at about λ^n**, including at (α, ℓ) = (1.5, 1), where 2αγ² = 9/4 lies beyond the threshold. They did not prove it.

Our exact errors of the classical rule (`scratchpad/lit/scripts/mm_cmp.py`):
- e_K/(2^{−1/4}(1+λ)^{3/4}λ^K) at K = 40 is 0.997853 (τ/σ = 1), 0.998926 (τ/σ = 1.5√2, KOG21's left panel) and 0.998232 (τ/σ = √2, 2αγ² = 1).

## 2. Karvonen–Särkkä 2019 (KS19)

### Parametrisation (§2.1, eqs. (2.1)–(2.5), p. 4)

- μ = N(0, 1) and the kernel is exp(−(x−y)²/(2ℓ²)), with ε = 1/(√2 ℓ), so the kernel is e^{−ε²(x−y)²}.
- They use the Mercer eigenpairs of Fasshauer–McCourt 2012 in L²(μ_α), with μ_α of variance 1/(2α²):
  - β = (1 + (2ε/α)²)^{1/4} and δ² = α²(β²−1)/2 (eq. (2.3));
  - eigenvalues (2.4);
  - eigenfunctions φ_n^α(x) = √(β/n!) e^{−δ²x²} H_n(√2αβx) (eq. (2.5)).
- After the dilation y ↦ y/τ: ε² = τ²/(4σ²) = γ²α of KW12. Hence ε² < 1/2 ⟺ 2αγ² < 1 ⟺ τ² < 2σ².

### The rule (Theorem 2.2, eq. (2.10)–(2.11), p. 8)

- The nodes are x̃_n = x_n^GH/(√2αβ); with α = 1/√2, x̃_n = x_n^GH/β, β⁴ = 1 + 8ε².
- These are the zeros of φ_N^α, i.e. the nodes of our **Mercer-zero scaling** v = τ²(1+2τ²/σ²)^{−1/2}: with τ = σ = 1, β⁴ = 3 = β_M⁴.
- The weights w̃ (eq. (2.11)) are those of the unique rule exact for the **first N eigenfunctions** only. They contain a polynomial correction Σ_m (2α²β²/(1+2δ²) − 1)^m H_{2m}(x^GH)/(2^m m!).
- So **KS19's rule is not a Q_{K,v}**: same nodes as the Mercer-zero Q_{K,v}, different weights.

### Convergence (§3)

- Theorem 3.1 (p. 12) is a generic bound e(Q_N) ≤ (1 + C₁W_N)C₂η^{M_N}, with η = √λ_M e^{1/β²} (λ_M the Mercer eigenvalue ratio) and W_N ≥ Σ|w_n|.
- Theorem 3.2 (p. 14) gives e = O(η^N) for the KS19 rule *under the assumption* sup_N Σ|w̃_n| < ∞, observed numerically (§5.3).
- The rates are conservative (§5.4, p. 22): empirical c ≈ 0.21 and 0.98 at ℓ = 0.2 and 1, against the theoretical 0.00033 and 0.054.
  - Side observation, not used: the Mercer-zero Q_{K,v} rate ln(1/q), q = (β_M²−1)/(β_M²+1), is 0.1997 and 0.9624 at these ℓ (τ/σ = √2/ℓ). This is close to KS19's empirical rates for their rule.
- **End of §3, p. 14:** the generalised Gaussian rule Q*_N with positive weights, exact for φ_0, …, φ_{2N−1} (Barrow 1978), "would establish an exponential convergence result with the 'double rate'". It is not analysed.
  - This rule is **our Mercer frame's rule**. Q_{K,v} is exact for e^{−δ²x²}·(polynomials of degree < 2K) iff 1/v = 1/τ² + 2δ², i.e. v = 2τ²/(β_M²+1), and ρ(v) = q.

### The "iff ε² < 1/2" sentence (§3, p. 11)

"When it comes to the RKHS of the Gaussian kernel … Kuo and Woźniakowski [19] have analysed convergence of the Gauss–Hermite quadrature rule. Unfortunately, it turns out that the Gauss–Hermite rule converges in this space if and only if ε² < 1/2."

- This is a statement about the **classical** rule (v = τ²), based on KW12: ε² < 1/2 ⟺ 2αγ² < 1 ⟺ τ < √2σ.
- It is **false** by our Theorem 4(b): the classical rule converges at every τ/σ, at rate λ. Table 4 covers τ/σ = √2 (ε² = 1/2) and 2√2.
- It rests on the lower bound of KW12 that KOG21 (footnote 1) and KSW17 (p. 830) call erroneous.
- KS19 observe nothing numerically against it. Their Fig. 5.4 plots the classical rule only to N = 20 or 60, at ℓ = 0.2–1.2; where ℓ < 1, λ = 1/(1+ℓ²) is close to 1.

## 3. Kuo–Woźniakowski 2012 and Kuo–Sloan–Woźniakowski 2017 (not read)

These papers are not reachable from this session, so no explicit error formula of KW12 could be evaluated against our exact e_K. What is known, all secondhand:

- *Upper bound (KOG21 §1).* For the classical rule, e ≤ b_n(α²/ℓ²)^n = b_n(2αγ²)^n, with b_n decreasing to 2^{−1/4}. The rate is 2αγ² = τ²/(2σ²), so it converges iff 2αγ² < 1.
- *KW12's claimed lower bound (KOG21 footnote 1).* It "matches" the upper bound and is erroneous per KSW17 p. 830. KS19's "iff ε² < 1/2" relies on it.
- *KSW17's exponential upper bounds (GRR24 p. 3).* They hold "under the additional assumption max σ_j < 1/2" (GRR's shape parameter, kernel e^{−σ²(x−y)²}, μ₀ = N(0, 1)). KSW17 Theorem 1.1 is quoted under related conditions (GRR24 §4.3, p. 18).

Comparison with the manuscript:
- **The manuscript's elementary bound** (Remark 3.3 since round 148, 3.4 before; from |ĝ_K(iB)| ≤ 1) is, for the absolute wce, wce(Q_{K,τ²}) ≤ b′_K(2αγ²)^K with b′_K = (K!/(2K)!)((4K−1)!!)^{1/2}.
  - b′_1 = √3/2, and (b′_{K+1}/b′_K)² = (16K²+16K+3)/(16K²+16K+4) < 1.
  - So b′_K decreases to 2^{−1/4}; numerically b′_160 = 0.841061 against 2^{−1/4} = 0.840896.
  - This has exactly the rate and limit constant of the KW12/KSW17 bound as KOG21 report it. Whether b_n = b′_n could not be checked.
- **Against a matching lower bound.** Since λ = x/(1+x) < x for x = 2αγ², Theorem 4(b) contradicts *any* lower bound c_n x^n with c_n^{1/n} → 1, for every α, γ, not only for 2αγ² ≥ 1.
- **Does KW12's explicit formula already give λ or the constant?** Unknown (not read). It is kept as a `% CHECK` in the manuscript and an open item in D62.

## 4. Gnewuch–Ritter–Rüßmann 2024 (GRR24, arXiv:2412.05368 v2)

### Spaces

- *Gaussian space.* Kernel ℓ_σ(x, y) = exp(−σ²(x−y)²) (§2.1), μ₀ = N(0, 1). In ours, σ_GRR² = τ²/(4σ²).
- *Hermite space.* Kernel k_β = Σ β^ν h_ν h_ν, with h_ν = He_ν/√ν! (§2.2).

### Transference (Theorem 4.3, p. 14; Remark 4.4, p. 16; Corollary 4.5)

- With 1 − β = 1/(1 + 2σ_GRR²), i.e. β = 2σ_GRR²/(1+2σ_GRR²) = λ, there is an explicit bijection A ↦ B of n-point rules with e(A, ℓ_σ) = e₀(ℓ_σ)·e(B, k_β).
- Worked out (σ = 1, general τ), the Gauss–Hermite rule B of N(0, 1) on H(k_β) corresponds to the scaled rule with:
  - nodes z_k/e, where e = ((1+4σ_GRR²)/(1+2σ_GRR²))^{1/2};
  - weights (τ/c)w_k exp(σ_GRR² z_k²/c²).
  - That is Q_{K,v} with v/τ² = (σ²+σ*²)/(2σ*²), ρ(v) = λ and u_v = +ln cosh r, the mirror image of the Fock scaling. Its t(v) is again λ/(2−λ).
- Equivalently, in the manuscript's Mehler remark at ρ = λ:
  - ℓ² = v, so v′ = 1;
  - σ*(1−ρ²)^{1/2}v/τ² = 1;
  - hence e_K(v)² = Σ_{n≥2K} λ^n Q_K[He_n]²/n!.
- Checked numerically (`scratchpad/lit/scripts/grr_cmp.py`) at K = 4, 8, 12 and τ/σ = 1, 2: the Hermite-space series equals `err2` to 13 digits (e.g. 1.061960669145e-11 at τ = 1, K = 8).

### Consequence

- Theorem A gives the **exact rate β/(2−β) and constant C(v)^{1/2} for the classical Gauss–Hermite rule in the Hermite space H(k_β)**.
- GRR24 §4.4 (p. 25) quote the bound e(B_n, k_β) ≤ Cβ^n from Irrgeher et al. 2015, Proposition 1 (not read here).
- GRR24 Theorem 4.6 (p. 17) gives upper bounds via KOG21 and the minimal-error lower bound C₁(σ)(σ²/(1+2σ²))^{2n}(n+1)^{−2}, i.e. rate (λ/2)² per point. This is compatible with ours ((λ/2)² < t_hs).

There is no overlap with the manuscript's other claims.

## 5. The other works the referee lists

- **Wang–Zhang, arXiv:2312.07940.** Root-exponential convergence of Hermite projection and interpolation, and of Gauss–Hermite quadrature (Theorem 5.1, Remark 5.2, p. 22), for functions analytic in a strip S_ρ with conditions at infinity. This is a different function class (the Gaussian RKHS consists of restrictions of entire functions) and is not the worst-case setting. No overlap.
- **Karvonen–Suzuki, arXiv:2209.12473.** Approximation (not integration against a Gaussian on ℝ) on [−1, 1] for weighted power-series and stationary kernels; for the Gaussian kernel, minimal errors ≍ (ε/2)^n(n!)^{−1/2} up to sub-exponential factors. No overlap.
- **Al-Jarrah 1983.** Not reachable; not cited (the `% CHECK` remains).

## 6. Novelty claims of the manuscript, checked against what was read

| Claim | Status |
|---|---|
| Exponential convergence of some scaled GH rule at every length-scale, two-sided exponential bounds | **known** (KOG21, Thm 2.5), for the Fock scaling |
| Exact rate t(v) for every v | **new**: KOG21 only bound the Fock rule (rates λ/2 and λ); GRR24 + IKLP15 give the upper rate λ for the ρ = λ rule; KS19 analyse no Q_{K,v} |
| Half-squeeze optimum v_hs and the cusp | **new**: no paper optimises v |
| Cusp law (1−t²)^{1/4}(πK)^{1/4}t^K ₂F₁^{1/2} | **new**; KOG21's n^{±1/4} share the Stirling origin, 2^n n!/√((2n)!) ~ (πn)^{1/4}, but are artefacts there (the Fock rule has no power of K) |
| Sharp constants C(v)^{1/2} | **new** |
| Classical rule converges at every length-scale | **new as a proof**: KOG21 leave it open, KS19 assert the opposite; KW12/KSW17 unread (CHECK) |
| Classical rate λ | **observed numerically** in KOG21 Fig. 1 (caption); proved here, with the constant 2^{−1/4}(1+λ)^{3/4} — the manuscript now says so |
| Theorem 1 (Peano kernel as a law) | not in these papers; KOG21 use only the mean-value remainder (2.2) |
| Corner ratio, transform limit | not in these papers |
| Exact rate β/(2−β) of GH in the Hermite space H(k_β) | **new consequence** (via GRR24's transference, equivalently the Mehler remark at ρ = λ), added to §1.3 |

The manuscript's novelty claim that needed correcting was (iv). It presented the classical rule's rate λ as new outright, but the rate had been observed numerically in KOG21 Fig. 1. §1.3 now says this. The other claims stand.

The named rules also gain their sources:
- the Fock scaling is KOG21's rule;
- the Mercer-zero nodes are KS19's nodes;
- the Mercer frame's rule is the generalised Gaussian rule named at the end of KS19 §3.
