<!-- Research note (2026-11-02): assumption (F) on the axes to every order and c₂ = s₁/2 as a theorem for the rules without a frame.
Written by a research subagent in the session's scratchpad (research_c2/); its scripts are copied to review/research_c2/. Refereed in
round 139 (review/REVIEW_round139_D3eh.md, D3eh): accepted with minor revisions, applied (see the Review history at the end). -->

# Assumption (F) proved — to every order — and the next term c₂ = s₁/2 as a theorem

Research note, 2026-11-02. Written in the working directory `W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/research_c2`;
the repository copy of its scripts and outputs is `review/research_c2/` (Files, at the end). The repository was read, not modified. Notation as in D‴.61, D‴.62 and Paper 18: y = η²/4, u = ν², s = s_d (signed), ϖ = s_o/s_d,
k = ϖ/(1 + ϖ), t = t(v) = |s_d|.

**Main result.** (F) holds, and with every order:

  log ĝ_K(√Kη) − Kη²/4 = ½ log(1 − y) + F₁(y)/K + F₂(y)/K² + … + F_{n−1}(y)/K^{n−1} + O(K^{−n})

for each n, uniformly on compact subsets of an open set containing the disc |η| < 1 and both axes' intervals (−2, 2) and i(−2, 2). Every
step there is by hand. The same holds on compact subsets of the whole disc |η| < 2, by hand since round 148 (Ω′ is the whole disc:
review/RESEARCH_note_annulus.md, Theorem A; until then with the certificates D and VA). F₁ is round 129's
F₁, now derived by hand. F₂(y) = −y(13 + 70y + 69y² + 24y³)/(8(1 − y)⁶) is identified numerically, not derived by hand.

Consequently, for every scaled rule without a frame (ρ ∉ (−1, 1), moment matching and ρ = −1 included), e_K² = C(v)t(v)^{2K}(1 + s₁/K +
s₂/K² + … + O(K^{−n})) with D‴.62 (d)'s s₁, and

  **e_{K−1}/e_K = t(v)^{−1}(1 + s₁/(2K²) + c₃/K³ + O(K⁻⁴)),  c₂ = s₁/2,  c₃ = (s₁ + 2s₂ − s₁²)/2,**

with no K⁻¹ term. This is a theorem, by hand, with no computer-assisted step.

One correction to D‴.62 (d), the Paper 18 supplement S1 and D‴.58, which state c₂ = s₁/2 "under (F)": the hypothesis is too weak. (F)
as stated, with remainder O(K⁻²), together with a 2-D Laplace remainder bound, gives only e_K² = Ct^{2K}(1 + s₁/K + O(K⁻²)): with S(K) = 1
+ s₁/K + ε(K), ε = O(K⁻²), the term ε(K − 1) − ε(K) of the ratio is O(K⁻²) in general. That gives "no K⁻¹ term" in the ratio, not c₂.
The minimal sufficient hypothesis is (F) with an F₂ term and remainder O(K⁻³), i.e. Proposition A below at n = 3; then the K⁻² term of
e_K² exists with an O(K⁻³) remainder. (D‴.62 (d) was labelled formal, and its last line already assumed that s₂ term.)

---

## 1. Statements and proofs

Status labels: **[proved]** means by hand, complete. **[routine]** means proved, with a standard step indicated rather than written out
to the last ε. **[numerical]** means numerical evidence only. **[sketch]** means the outline of a proof.

### Lemma L (Laplace's method to every order, uniformly in a parameter) [proved; standard, written out with the parameter explicit]

Lemma L is classical: Laplace's method to every order with a uniform remainder, the coefficients by Wick's (Feynman-diagram)
bookkeeping. See Olver, *Asymptotics and Special Functions* (Academic Press, 1974), ch. 3 [memory]; Wong, *Asymptotic Approximations of
Integrals* (Academic Press, 1989), ch. IX [memory]; Bleistein–Handelsman, *Asymptotic Expansions of Integrals* (1975), ch. 8 [memory];
Kirwin, *Higher asymptotics of Laplace's approximation*, Asymptot. Anal. 70 (2010), no. 3–4, 231–248, arXiv:0810.1700; Nemes, *An
explicit formula for the coefficients in Laplace's method*, Constr. Approx. 38 (2013), no. 3, 471–487, arXiv:1207.5222. (W2) has the
structure of the classical 1/n correction term of Laplace approximations in several variables (Shun–McCullagh, J. R. Stat. Soc. B, 1995;
CHECK: volume and pages, known only through secondary sources). It is written out here only to make the uniformity in the parameter p
explicit.

Let d ≥ 1, n ≥ 1, x₀ > 0, δ > 0, M < ∞ and B = {x ∈ ℝ^d : |x| ≤ x₀}. For p in an index set 𝒫 let:

- H_p be a symmetric d × d matrix: complex if d = 1, with Re H_p ≥ 2δ, and real positive definite, H_p ≥ 2δI, if d ≥ 2; and |H_p| ≤ M.
- h_p ∈ C^{2n+3}(B) and a_p ∈ C^{2n+1}(B), with h_p(0) = 0, ∇h_p(0) = 0, ∇²h_p(0) = 0, and all derivatives up to these orders bounded by M.
- Re(−½xᵀH_px + h_p(x)) ≤ −δ|x|² on B.
- b_{p,K} measurable with |b_{p,K}| ≤ ε_K.

h_p and a_p may be complex-valued, also for d ≥ 2, as long as H_p is real positive definite there.

Then, uniformly in p,

  ∫_B e^{K(−½xᵀH_px + h_p(x))} a_p(x)(1 + b_{p,K}(x)) dx = (2π/K)^{d/2} det(H_p)^{−1/2} [ Σ_{m<n} L_m(p) K^{−m} + O(K^{−n} + ε_K) ],

with the principal root for d = 1. The coefficients are given by Wick's rule, with X a centred Gaussian of covariance G = H_p^{−1}:

  L_m = Σ_{k ≤ 2m} E[ a_p^{(≤2n)}(X)·(h_p^{(3…2n+2)}(X))^k / k! ]_{homogeneous degree 2(k+m)},

where a^{(≤2n)} and h^{(3…2n+2)} are Taylor polynomials. So L₀ = a(0). Only the Taylor degrees ≤ 2m + 2 of h and ≤ 2m of a enter L_m
(k = 1 needs h's degree 2m + 2), so L_m does not depend on n. For d = 1 with phase −αx² + h (H = 2α) and Taylor coefficients
a_j and h_j,

  L₁ = a₂/(2α) + 3a₁h₃/(4α²) + a₀(3h₄/(4α²) + 15h₃²/(16α³)).            (W1)

For diagonal G (as in every application below) and log-amplitude ℓ = log a, with h_{ijk…} the partial derivatives at 0,

  L₁/a₀ = ½Σ_i G_ii(ℓ_ii + ℓ_i²) + ½Σ_i ℓ_iG_ii Σ_k h_{ikk}G_kk + ⅛Σ_{i,k} h_{iikk}G_iiG_kk
          + ⅛Σ_k G_kk(Σ_i h_{iik}G_ii)² + (1/12)Σ_{ijk} h_{ijk}² G_iiG_jjG_kk.        (W2)

*Proof.* This is Paper 18's Lemma 4.1 carried to every order, with the same three regions.

1. The b-part is at most Mε_K∫e^{−Kδ|x|²} = O(ε_K K^{−d/2}).
2. Put x_K = K^{−1/3}. On x_K ≤ |x| ≤ x₀ the integrand is O(e^{−δK^{1/3}}).
3. On |x| ≤ x_K write h = T + r, with T the Taylor polynomial of degrees 3 … 2n+2 and |r| ≤ C|x|^{2n+3}, and a = A + r_a, with A of
   degree ≤ 2n and |r_a| ≤ C|x|^{2n+1}. Then |KT| ≤ CK|x|³ ≤ C, so:
   - e^{Kh} = e^{KT}(1 + O(K|x|^{2n+3})), and
   - e^{KT} = Σ_{k≤2n} (KT)^k/k! + O(K^{2n+1}|x|^{3(2n+1)}).

   Expanding the product with A gives monomials K^k x^ν with |ν| ≥ 3k. Integrated against e^{−KxᵀHx/2} over ℝ^d, each one contributes
   (2π/K)^{d/2} det^{−1/2} · K^{k−|ν|/2} · E[X^ν]. The tails beyond x_K are O(e^{−δK^{1/3}}). For d = 1 with complex α the Gaussian
   moments follow by analytic continuation in α, as in Lemma 4.1.

   The monomials with k − |ν|/2 = −m < −n + 1, and every remainder, are O(K^{−n}) after normalisation. The remainders are bounded by
   C∫e^{−Kδ|x|²}(|x|^{2n+1} + K|x|^{2n+3} + K^{2n+1}|x|^{6n+3}) dx = O(K^{−d/2−n−1/2}). Finally |(2π/K)^{d/2} det^{−1/2}| ≥ cK^{−d/2}.

The constants depend on n, d, x₀, δ and M only. ∎ (W1) and (W2) are the m = 1 terms written out; (W2) assumes G diagonal. The engine
`review/research_c2/laplace.py` implements this formula for every m.

All data below are analytic and bounded on a fixed complex neighbourhood of B, uniformly in the parameter. Cauchy's estimates then give
the derivative bounds. The K-dependent amplitudes are of the form a⁽⁰⁾ + a⁽¹⁾/K + … + a⁽ⁿ⁻¹⁾/K^{n−1} times (1 + O(K^{−n})), with analytic
a⁽ʲ⁾. Apply the lemma to each a⁽ʲ⁾ and put the remainder into b.

### Proposition A (the transform to every order)

**(i) I_K [proved].** He_K(√K(ζ + 1/ζ)) = K^{K/2}ζ^{−K}e^{Kζ²/2}I_K(ζ) (Paper 18, Lemma 4.3). For n ≥ 1 and r < 1,
uniformly on |ζ| ≤ r,

  (1 − ζ²)^{1/2} I_K(ζ) = 1 + Σ_{m=1}^{n−1} b_m(ζ)K^{−m} + O(K^{−n}),  b₁(ζ) = −3ζ⁴/(4(1 − ζ²)²) − 5ζ⁶/(6(1 − ζ²)³),

with b_m analytic on |ζ| < 1.

*Proof.* As in Lemma 4.3, with Lemma L instead of Lemma 4.1:
- α = (1 − ζ²)/2 and h(x) = log(1 + iζx) − iζx − ζ²x²/2, analytic and bounded for |x| ≤ 1 and |ζ| ≤ r.
- Re(−αx² + h) ≤ −(1 − r²)x²/2, and the region |x| ≥ ½ is O(e^{−K(1−r²)/8}).
- h₃ = −iζ³/3, h₄ = −ζ⁴/4 and a ≡ 1, so (W1) gives b₁.

The b_m are polynomials in the Taylor coefficients of h, which are polynomials in ζ, times powers of 1/α. ∎

(i) is the exterior Plancherel–Rotach expansion of He_K, in this parametrisation. Its coefficients and error bounds are published: Shi,
Nemes, Wang and Wong, *Error bounds for the asymptotic expansions of the Hermite polynomials*, Proc. Roy. Soc. Edinburgh Sect. A,
arXiv:2109.01054 (CHECK: volume, year and pages; the arXiv page gives no journal reference). The proof is given here for the uniformity on
|ζ| ≤ r in the form used below.

**(ii) The corner ratio [proved].** For n ≥ 1 and r < 1, uniformly on |ν| ≤ r (ν = 0 included),

  c_K(√Kν) = (1 − ν²)(1 + Σ_{m=1}^{n−1} γ_m(ν²)K^{−m} + O(K^{−n})),  **γ₁(u) = u(1 + u)/(1 − u)³**,

with γ_m analytic on the unit disc.

*Proof.* In Paper 18's (4.2), substitute ζ = νe^{iθ} and use the clockwise orientation. The prefactor ν^{1−K}e^{−Kν²/2}e^{Kf_ν(ν)} = νe^K
cancels the 1/ν of the amplitude, so

  c_K(√Kν) = (K!K^{−K}e^K) · (1/2π) ∫_{−π}^{π} e^{Kφ_ν(θ)} (e^{−iθ} − ν²e^{iθ}) I_K(νe^{iθ})^{−1} dθ,
  φ_ν(θ) = ν²(e^{iθ} − 1 − (e^{2iθ} − 1)/2) + e^{−iθ} − 1 + iθ = −(1 − ν²)θ²/2 + i(1 + 3ν²)θ³/6 + (1 − 7ν²)θ⁴/24 + O(θ⁵).

This holds for 0 < |ν| ≤ r and K ≥ K₀, since by (i) I_K has no zero on |ζ| ≤ r; so K₀ depends on r, through I_K being zero-free on |ζ| ≤
r. Both sides are continuous at ν = 0, where the right side equals 1 = c_K(0) (coefficient extraction). Equality on 0 < |ν| ≤ r plus
continuity at ν = 0 is enough: the formula holds on the whole disc |ν| ≤ r, and **the maximum-principle step of Lemma 4.3 is not
needed.**

Lemma L applies with x₀ = π:
- Re φ_ν = −(1 − cos θ)(1 − |ν|²cos(2 arg ν + θ)) ≤ −(2/π²)(1 − r²)θ², from Lemma 4.3.
- The amplitude is A(θ) = (e^{−iθ} − ν²e^{iθ})(1 − ν²e^{2iθ})^{1/2}, times 1 + Σ b̃_m(νe^{iθ})K^{−m} + O(K^{−n}), where b̃ comes from
  inverting (i). It is analytic on a strip |Im θ| < w with re^{w} < 1.
- Stirling: K!K^{−K}e^K = (2πK)^{1/2}(1 + 1/(12K) + 1/(288K²) + …).

For γ₁, use (W1) with α = (1 − u)/2, h₃ = i(1 + 3u)/6, h₄ = (1 − 7u)/24 and A₀ = (1 − u)^{3/2}. Write q = 1/(1 − u). Then:
- A₁/A₀ = −i(1 + 2u)q.
- A₂/A₀ = −½ + uq + u²q²/2 − u(1 + u)q².
- The K⁻¹ part of the amplitude at θ = 0 is −b₁(ν).

The six contributions, each times (1 − u)³, are:

| contribution | value |
|---|---|
| A₂/(2αA₀) | −½ + u − 2u² |
| 3A₁h₃/(4α²A₀) | (1 + 5u + 6u²)/2 |
| 3h₄/(4α²) | (1 − 8u + 7u²)/8 |
| 15h₃²/(16α³) | −5(1 + 6u + 9u²)/24 |
| −b₁ | 3u²(1 − u)/4 + 5u³/6 |
| Stirling 1/12 | (1 − u)³/12 |

They sum to u + u². ∎

**(iii) The transform [proved].** Let Ω′ be the valley-condition set of Paper 18, Lemma 4.6 (ii) (the valley condition is its Definition
4.5). By hand it contains an open set Ω_h ⊃ {|η| < 1} ∪ (−2, 2) ∪ i(−2, 2): the domination condition holds there and is open (Lemma 4.6
(i), and the axis gaps in the proof of Theorem 2).
Ω′ is the whole disc |η| < 2, by hand since referee round 148 (Paper 18's Proposition 4.8, the research note
review/RESEARCH_note_annulus.md's Theorem A; until then with certificates D and VA). The proof below uses Ω′ only.
For every compact 𝒦 ⊂ Ω′ and n ≥ 1, uniformly on 𝒦,

  ĝ_K(√Kη)e^{−Kη²/4} = (1 − y)^{1/2}(1 + Σ_{m=1}^{n−1} G_m(y)K^{−m} + O(K^{−n})),  y = η²/4,

with G_m analytic on |y| < 1. Taking logarithms (principal branch of the bracket; ĝ_K > 0 on both axes for K ≥ K₀),

  log ĝ_K(√Kη) − Kη²/4 = ½ log(1 − y) + Σ F_m(y)K^{−m} + O(K^{−n}),  **F₁ = G₁ = round 129's F₁**, F₂ = G₂ − G₁²/2.

*Proof.* Repeat the proof of Lemma 4.6 (ii) verbatim, with Lemma L in place of Lemma 4.1 on the sub-chord.
- The phase is χ_η(z) = −α_ηz² + h(z), with α_η = 1 − y and h = log(1 − z²) + z² independent of η.
- The amplitude is a_η(z)(1 + Σ γ_m(y(1 + z)²)K^{−m} + O(K^{−n})), with a_η(z) = (1 − y(1 + z)²)/(1 + z). This uses (ii) at ν = (η/2)(1 +
  z), which lies in a fixed compact subset of the unit disc.
- Off the sub-chord the integral is exponentially smaller, unchanged from the proof of Lemma 4.6 (ii).
- Stirling: K·C(2K, K)·4^{−K} = (K/π)^{1/2}(1 − 1/(8K) + 1/(128K²) + …).

G₁ comes from (W1) with a₀ = 1 − y, a₁ = −(1 + y), a₂ = 1, h₃ = 0 and h₄ = −½:

  G₁ = 1/(2(1 − y)²) − 3/(8(1 − y)²) + γ₁(y) − 1/8 = γ₁(y) + 1/(8(1 − y)²) − 1/8.

This is round 129's F₁, because y/(8(1 − y)²) + y/(8(1 − y)) = (2y − y²)/(8(1 − y)²) = 1/(8(1 − y)²) − 1/8. The G_m are functions of y
alone, since c_K is even and only ν² = y(1 + z)² enters. ∎

So (F) holds on Ω_h — in particular on both axes' intervals — with all orders and by hand.

*[After round 154 (review/REVIEW_round154_D3ew.md, D3ew-1).]* The proof of Lemma 4.6 (ii) gives local uniformity only near η₀ ≠ 0:
(4.4) has η^{−2K}, and the scaled path (η/η₀)𝒞₀ divides by η₀; at η = 0 the valley condition is vacuous (Paper 18's Definition 4.5
now takes 0 < |η| < 2, so Ω′ is the punctured disc). Uniformity on compact subsets of the whole disc, η = 0 included, follows by the
maximum principle: K^n[ĝ_K(√Kη)e^{−Kη²/4}(1 − y)^{−1/2} − 1 − Σ_{m<n} G_m(y)K^{−m}] is analytic on |η| < 2 and bounded on each circle
|η| = r < 2 uniformly in K (there the expansion holds uniformly), hence on |η| ≤ r. Nothing else changes; no number changes. Paper 18's
supplement §S1 now states Propositions A and B and Corollary B′ with their proofs (Lemma L, (i)–(iii), the γ₁ table, G₁, the Wick terms
T1–T5 and s₁), re-checked independently in floating point at random complex points in that revision.

**(iv) F₂ [numerical].** F₂(y) = −y(13 + 70y + 69y² + 24y³)/(8(1 − y)⁶). The evidence:
- The engine's F₂ (the exact Laplace coefficient, evaluated in 300-bit balls) times (1 − y)⁶ at 14 points y = −j/16, interpolated to
  degree 13: coefficients 0, −13/8, −35/4, −69/8, −3, then 0 to 10⁻⁶⁹, and these rationals reproduce the engine to 10⁻⁷⁶.
- The closed form against the engine at 6 more points, complex points included: within 10⁻⁸⁶.
- The closed form against the exact cumulants' 1/K² coefficients e_j for j ≤ 10: equal in rationals (§3).

A degree bound for F₂(1 − y)⁶, which would turn this into a proof, was not written out. F₂'s closed form is not needed for c₂: c₂ needs
only that the F₂ term exists with an O(K⁻³) remainder, which (iii) at n = 3 proves; the closed form enters only the values of s₂ and c₃.

### Proposition B (every rule without a frame: the full expansion and s₁) [proved]

Let v be fixed with ρ(v) ∉ (−1, 1). For every n ≥ 1,

  e_K² = C(v)t(v)^{2K}(1 + Σ_{m=1}^{n−1} s_mK^{−m} + O(K^{−n})),
  s₁ = 1/16 + 2s(1 + s)/(1 − s)³ + s/(4(1 − s)²) − ϖ(1 + s)²/(8(1 + ϖ)(1 − s)²),  s = s_d,

which is D‴.62 (d)'s s₁, and also round 129's (at ϱ = −s, √q = −ϖ).

*Proof.* Take ρ > 1 first; the integrand is in D‴.62 (a) and Paper 18 Theorem 4 (a). Put p = A/√K and p′ = B/√K. With
H_K(y) := ĝ_K(i√Ky)e^{Ky²/4}(1 + y²/4)^{−1/2} = exp(Σ F_m(−y²/4)K^{−m} + O(K^{−n})) on |y| ≤ c′ < 2 (Proposition A on the imaginary axis),

  e_K² = σ*(v/τ²)(K!/(2K)!)²K^{2K}(K/2π)(|c||c_b|)^{−1/2} ∫∫ e^{KΨ} m H_K(p + p′)H_K(p − p′) dp dp′,
  Ψ = 2 log(p′² − p²) − p²/(2|s_a|) − p′²/(2|s_b|),  m = ((1 + (p + p′)²/4)(1 + (p − p′)²/4))^{1/2}.

1. **Off the box** |p| + |p′| < c′ (c′ close to 2): exponentially smaller than t^{2K}. This is D‴.62 Step 3 / Paper 18 Lemma 4.10 (4.8 until round 148), unchanged.
2. **Inside the box, away from neighbourhoods U_± of the two dominant points (0, ±2b), b² = |s_d|:** Ψ ≤ Ψ_max − ε. This follows from
   D‴.62 Lemma A (the unique maximisers), upper semicontinuity, and Ψ = −∞ on the diagonals. With |mH_KH_K| ≤ C on the box and ∫e^{Ψ} <
   ∞, this part is O(e^{(K−1)(Ψ_max−ε)}): exponentially smaller. It contains the second ridge, so J_o is beyond every power of 1/K.
3. **On U_±:** Lemma L with d = 2.
   - The Hessian is diag((1 + ϖ)/(ϖb²), 2/b²), so G = (kb², b²/2), k = ϖ/(1 + ϖ).
   - The amplitude is m·exp((F₁(−y₁²/4) + F₁(−y₂²/4))/K + …), with H_K's remainder in b.
   - Ψ_max = 2 log(4b²) − 2, and K^{2K}(K!/(2K)!)² = ½(e/4)^{2K} exp(1/(12K) + O(K⁻³)), so the prefactor is ½t^{2K}(1 + 1/(12K) + …).
     The K-powers cancel: (K/2π)·(2π/K) = 1.
   - The leading constant is C(v), as in D‴.62 (b).

For s₁, take derivatives of Ψ at (0, 2b), in p and q = p′ − 2b:
- h_{qqq} = h_{ppq} = 1/b³.
- h_{pppp} = h_{ppqq} = h_{qqqq} = −3/(2b⁴).
- The odd-in-p derivatives vanish.

The log-amplitude: ℓ_q = b/(1 + b²), ℓ_pp = ℓ_qq = (1 − b²)/(2(1 + b²)²), and ℓ_p = ℓ_pq = 0. The terms of (W2), in s = −b², are:

| term of (W2) | value |
|---|---|
| T1 = ½ΣG(ℓ″ + ℓ′²) | −sk(1 + s)/(4(1 − s)²) − s(1 + s)/(8(1 − s)²) + s²/(4(1 − s)²) |
| T2 | −s(k + ½)/(4(1 − s)) |
| T3 + T4 + T5 | −k/8 − 1/48 |

Add 2F₁(s), from the amplitude's K⁻¹ part, and 1/12 from Stirling:
- The k-terms collect to −(k/8)(1 + s)²/(1 − s)².
- The rest give 1/16 − s/(4(1 − s)) + 2F₁(s) = 1/16 + 2s(1 + s)/(1 − s)³ + s/(4(1 − s)²).

**ρ < −1** (the real axis): Proposition A on (−2, 2), and Paper 18's c̃ absorbs e^{+K(y₁²+y₂²)/4}. Then Ψ = 2 log(p² − p′²) − p²/(2s_a) −
p′²/(2s_b), maximised at (±2b, 0), b² = s_a, with m = ((1 − y₁²/4)(1 − y₂²/4))^{1/2} and F₁(+y²/4). Compared with ρ > 1:
- G and h are the same in b.
- ℓ_q = −b/(1 − b²) and ℓ″ = −(1 + b²)/(2(1 − b²)²).
- So every term of (W2) is the same function of s = +b²: the same s₁.

**ρ = ±1:** one-dimensional, in the dominant variable only. The terms are those above with k = 0, by (W1). ∎

The engine re-derives s₁ in exact rationals at 24 random (b, ϖ), on both axes and including ϖ = 0: 24/24 equal to the closed form. ✓

**Corollary B′ (the ratio per node) [proved].** Write S(K) = e_K²/(Ct^{2K}) = 1 + s₁/K + s₂/K² + s₃/K³ + O(K⁻⁴) and log S = Σλ_mK^{−m}.
Then t·e_{K−1}/e_K = exp(½[log S(K − 1) − log S(K)]) = 1 + s₁/(2K²) + (s₁ + 2s₂ − s₁²)/(2K³) + O(K⁻⁴). This has an asymptotic expansion to
every order in 1/K, with no K⁻¹ term. Only n = 3 of Proposition A is needed for c₂ with an O(K⁻³) remainder.

**Remark F (frame rules off the cusp) [sketch].** The same proof applies to D‴.61 (iv)'s Fock integral E[φ_K(ξ)φ_K(ξ̄)], ξ = √K(p + ip′):
Ψ = 2 log(p² + p′²) − p²/(2s_a) − p′²/(2|s_b|), with m = |1 − η²/4| and the F₁ term 2Re F₁(η²/4).
- The mass sits at one axis point, off the cusp. Proposition A is needed only in a complex neighbourhood of it, which is by hand.
- Exponential smallness elsewhere needs D‴.61 (iii)(b) on |η| ≤ c, so the result needs D‴.61 (iv)'s disc: by hand since round 148 (until
  then computer-assisted for r ≥ 0.0966).
- The Wick terms are those of (W2) with the non-dominant direction rotated by i. This sends h_{ii…} → ±h, G_ii → −G_ii in pairs, so the
  same s₁(s, ϖ) results, now with ϖ < 0. The engine confirms this at 8 frame rules (§3), e.g. s₁ = −0.088, c₂ = −11/250 at (2, Mercer
  frame), round 129's value.

This would make D‴.58's c₂ statement a theorem off the cusp for the frame rules too, at D‴.61 (iv)'s tier. The cusp v_hs needs a separate,
degenerate, angular expansion and is not covered.

---

## 2. What remains open, and why

1. **No explicit constants.** The O(K^{−n}) constants depend on Ω's margin, Lemma T's gap δ and the box choice c′. These are proved
   positive but not quantified. So there is still no K-explicit bound (Paper 18 open problem (4); (5) until round 148): c₂ is an asymptotic statement.
2. **The transient B₀q^K** (D‴.62 (c)) stays formal. The second ridge is O(|s_o|^{2K}), beyond every power of 1/K. Proposition B's series
   does not determine it, and B₀q^K "as a theorem" needs J_d to relative accuracy o(q^K) — a resurgence-type statement about a Gevrey-1
   series (round 129 (f)). Numerically it is visible at the outlier (τ/σ = ¼, 2v_hs): there the extrapolated c₂ agrees only to 7.5·10⁻⁹,
   against 10⁻²⁰ elsewhere on the imaginary branch.
3. **Uniformity in v.**
   - Not near ρ = ±1: the p-curvature (1 + ϖ)/(ϖb²) → ∞ as ϖ → 0, so the 2-D expansion is not uniform onto the 1-D one.
   - Not near the cusp v_hs, where ϖ → −1.
   - Not as t(v) → 1.
4. **F₂ in closed form** is identified, not proved: an exact-rational (Q(i)(y)) run of the engine together with a degree bound would
   settle it. The general F_m (denominators (1 − y)^{3m}?) is open. Neither is needed for c₂.
5. **Frames (Remark F)** are written as a sketch. Writing them out is routine; they inherit D‴.61 (iv)'s disc, by hand since round 148 (until then
   computer-assisted for r ≥ 0.0966).
6. **P_j's polynomiality** (Paper 18 open problem (1); (2) until round 148) is still unproved. It is not used here; it is checked again for j ≤ 10
   (`check_cumulants.py`). Route (2) of the brief, cumulants for all j, is therefore not needed: the Laplace route (1) gives the
   expansion directly and uniformly, beyond the cumulant series' disc |η| < 2.

---

## 3. Numerical checks (scripts and outputs in `review/research_c2/`)

Every run used `OMP_NUM_THREADS=1 nice -n 15 venv/bin/python <script>`, one at a time; `check_exact.py` took 34 s, the others under 1 s.
The exact values are python-flint balls, with the nodes' midpoints taken exact (as in exp_frameless.py). The extrapolations are Neville in
1/K and are not rigorous.

| script | what | result |
|---|---|---|
| `check_algebra.py` A | F₁ = γ₁ + 1/(8(1 − y)²) − 1/8 | exact at 50 random rationals |
| `check_algebra.py` B | Wick engine (independent of the hand Taylor coefficients) vs b₁, γ₁, F₁ at complex points | ≤ 2.4·10⁻⁸⁹ |
| `check_algebra.py` C | ridge s₁ by the engine in exact rationals vs D‴.62 (d) | 24/24 exact; −13/125 at MM τ = ½ |
| `check_exact.py` I | exact I_K, K = 100 … 3200, vs b₁ (and the engine's b₂) | b₁ to 7·10⁻²⁰ (ζ = 0.5i), 2·10⁻¹⁷ (0.3 + 0.4i), 2·10⁻¹³ (0.5); b₂ to 2·10⁻²⁰ … 3·10⁻¹² |
| `check_exact.py` C | exact c_K(√Kν), K = 40 … 960, vs γ₁, γ₂ | γ₁ to 7·10⁻¹⁶ (ν = 0.5i), 1·10⁻¹⁶ (0.95i), 7·10⁻¹¹ (−0.4 + 0.55i), 7·10⁻⁸ (0.5); γ₂ to 2·10⁻¹⁵ … 7·10⁻⁶ |
| `check_exact.py` G | exact ĝ_K(√Kη), K = 40 … 960: extrapolated D_K vs F₁, and K(D_K − F₁) vs the engine's F₂ | imaginary axis y = 0.5, 1, 1.5, 1.9: F₁ to ≤ 6.5·10⁻¹⁸, **F₂ to ≤ 2·10⁻¹⁷**; complex 0.6 + 0.5i, 0.3 + 1.2i: F₁ to 2·10⁻¹³, 3·10⁻¹⁶, F₂ to 3·10⁻¹², 3·10⁻¹⁵; real η = 1: F₁ to 5·10⁻⁹, F₂ to 3·10⁻⁷ |
| `fit_F2.py` | F₂'s closed form | coefficients rational to 10⁻⁷⁶; 6 more points to 10⁻⁸⁶ |
| `check_cumulants.py` | exact cumulants (integer path sums → P_j, j ≤ 10, polynomial and checked at 5 extra K) | c_j, d_j (F₁) and e_j (F₂) all equal in rationals, j = 1 … 10 |
| `check_ridge.py` | stored exact e_K (frameless_raw.npz S, 17 rules) | below |
| `check_s2.py` | s₂ predicted (engine, order K⁻², with F₂) vs s₂ extracted from the exact errors | below |
| `check_frame.py` | Remark F's s₁ at 8 frame rules vs s₁(s, ϖ < 0) | ≤ 8.6·10⁻⁸⁸ |

One caveat on the rows above:
- On the real axis near the edge (η = 1.5; ν = 0.8; ζ = 0.8) the sequences are still pre-asymptotic at these K, with large higher
  coefficients (e.g. γ₂(0.64) = −4778). They are consistent (F₂(9/16) predicted −786.89, extrapolated −785.16) but not sharp. This is the
  known slow approach near the branch point.

**Theorem B against the exact errors** (Neville through 9 points K_max/2 … K_max):

| rules | c₂ = s₁/2 | s₂ predicted vs extracted |
|---|---|---|
| imaginary branch (4 MM, 2 MM beyond Kuo–Woźniakowski, 4 with ρ > 1 off the outlier) | ≤ 4.4·10⁻²⁰ | ≤ 1.2·10⁻¹⁹ |
| outlier (¼, 2v_hs), q^{K/2} = 9·10⁻¹⁵ | 7.5·10⁻⁹ | 8.6·10⁻⁸ |
| slow branch ρ ≤ −1, 6 rules | 1.7·10⁻¹² … 2.1·10⁻⁴ | 9.8·10⁻¹² … 1.7·10⁻² |

The predicted s₂ include, e.g., 0.0839040000 (MM, τ = ½), −30.7153839 (¼, 0.5v_hs) and −263.0676 (1, 0.5v_hs). The last two are consistent
with D‴.62's fitted −30.715381 and −263.05374, but on the slow branch only to about 10⁻⁴ relative: −263.0676 differs from the fit by
1.4·10⁻², and from referee round 139's own extraction (K ≤ 360), −263.018, by 5·10⁻². That is the precision the very large higher
coefficients there allow, not a sharp test.

The prediction of s₂ is the sharpest test of the all-orders statement. It uses F₂, the 2-D Laplace coefficient of order 2 and Stirling, and
none of the data.

---

## 4. How it would enter

**D‴.62 (d).**
- Replace "Assume (F) … Then, formally" by a reference to Proposition A, which gives (F) to every order on an open set containing the axes'
  intervals, by hand.
- State Proposition B and Corollary B′: e_K² = Ct^{2K}(1 + s₁/K + s₂/K² + …) for every rule without a frame, and e_{K−1}/e_K =
  t^{−1}(1 + s₁/(2K²) + c₃/K³ + …).
- Status (d): "proved (by hand), given D‴.61 (i), Duhamel's identity and (iii)(a)'s representations"; the remainder bound is Lemma L plus
  Steps 1–2 of the proof.
- The size-biasing computation of (d) can stay as a second evaluation of the same coefficient. The rigorous route is Laplace in (p, p′),
  with the Wick terms T1–T5 above.
- Add the correction: (d)'s hypothesis is too weak — (F) with remainder O(K⁻²) gives only "no K⁻¹ term"; c₂ needs (F) with an F₂
  term and remainder O(K⁻³), i.e. Proposition A at n = 3.
- The "What stays open" list loses "c₂ … rigorously, which needs (F)". It keeps the transient Bq^K and K-explicit bounds.
- Numerics (G) gains the s₂ prediction and F₂.

**D‴.58 (ii), c₂.**
- "formal at next order — D‴.61 and D‴.62 prove only the leading limit; the next-order uniformity and the Laplace remainder are not proved"
  becomes: proved for the 14 rules without a frame, by hand (D‴.62 (d) / Proposition B).
- For the frame rules off the cusp it becomes proved at D‴.61 (iv)'s tier (Remark F, once written out); the cusp stays formal.
- The slow branch's large s₁ (and s₂ = −30.7 … −263) becomes a theorem about the coefficients, which explains "K^{−1.52 …}" rigorously.

**Paper 18.**
- Lemma 4.1 → Lemma L, all orders.
- Lemma 4.3 (corner ratio) gains "with a full expansion, γ₁ = u(1 + u)/(1 − u)³", and drops the maximum-principle step (the θ-integral
  is regular at ν = 0).
- Lemma 4.6 (ii) and Theorem 2 gain the expansion, G₁ = F₁ — by hand on the axes and |η| < 1, computer-assisted on the disc (by hand on the whole disc
  since round 148).
- Theorem 4 gains (e): the full expansion and c₂ = ς₁/2.
- Supplement S1 moves into the paper as a theorem, with the correction above, and keeps the transient (S2) formal.
- Open problem (4) shrinks to the transient and the cusp/uniformity. (5) is unchanged. (Numbered (3) and (4) since round 148.)
- The numerical section can quote the s₂ prediction (agreement ≤ 10⁻¹⁹ on the imaginary branch).
- Table 1's ς₁ entry: "formal" → "Theorem 4 (e)".

---

## Files

The repository copy (written as W/NOTE.md, W/scripts/ and W/*.out in the working directory):
- `review/RESEARCH_note_c2.md` (this note)
- `review/research_c2/lib.py`: nodes, closed forms, Neville
- `review/research_c2/laplace.py`: Wick engine
- `review/research_c2/check_algebra.py`, `check_exact.py`, `check_ridge.py`, `check_s2.py`, `check_frame.py`, `check_cumulants.py`,
  `fit_F2.py`
- `review/research_c2/check_*.out`, `fit_F2.out`: the outputs quoted above

`check_ridge.py` (and `check_s2.py`, `check_frame.py` through it) reads `results/chaos/frameless_raw.npz` and imports
`experiments/exp_frameless.py` from the repository root, taken from the environment variable `LCK_REPO` if set and otherwise from the
script's own location (two directories up).

---

## Review history

- Round 139 (review/REVIEW_round139_D3eh.md, D3eh, 2026-11-02): **accepted with minor revisions**; no mathematical error found, the
  hand algebra and the numerics reproduced independently. Changed: Paper 18 cross-references corrected to Lemma 4.6 (i)/(ii), Definition
  4.5 and Theorem 2 (D3eh-1); Lemma L labelled standard with literature, and Proposition A (i) cited as the exterior Plancherel–Rotach
  expansion (D3eh-2); the correction to D‴.62 (d) rephrased — its hypothesis is too weak, (F) with O(K⁻²) gives only "no K⁻¹ term", c₂
  needs Proposition A at n = 3 (D3eh-3); Lemma L's statement says that h, a may be complex for d ≥ 2, that L_m does not depend on n, and
  that (W2) assumes diagonal G (D3eh-4); the ν = 0 step and K₀'s dependence on r stated (D3eh-5); the slow-branch s₂ agreement quoted as
  about 10⁻⁴ relative instead of "now explained" (D3eh-6); the Files section points to `review/research_c2/` and `check_ridge.py`'s
  REPO is parametrised (D3eh-7). No result, proof or number changed; Lemma L's label now says it is standard.
- Round 148 (review/REVIEW_round148_D3eq.md, D3eq, 2026-11-05; a review of the annulus note review/RESEARCH_note_annulus.md, not of
  this one): Ω′ is the whole open disc by hand. Propagated here: Proposition A's statement and (iii) — the expansion to every order on
  compact subsets of the whole disc is by hand, as (iii)'s proof uses Ω′ only —, Remark F's dependence on D‴.61 (iv) (by hand; Remark F
  stays a sketch, not refereed as a proof), §2 item 5, §5's Paper 18 line; Paper 18's renumbering (Lemma 4.8 → 4.10, open problems
  (2), (5) → (1), (4)). Nothing else changed; no number changed.
- Round 154 (review/REVIEW_round154_D3ew.md, D3ew, 2026-11-05; a review of Paper 18, not of this note): the gap at η = 0 in the
  uniformity of Proposition A (iii) (D3ew-1), closed by the maximum principle (a dated note after (iii)'s proof). Paper 18's supplement
  §S1 now contains Propositions A and B and Corollary B′ with proofs (D3ew-5), so the paper no longer cites this note. No number changed.
