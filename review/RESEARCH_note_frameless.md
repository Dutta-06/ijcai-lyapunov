<!-- Research note (2026-11-02), not refereed: the scaled Gauss–Hermite rules without a Mehler frame. Written by a research subagent in
the session's scratchpad (research62/, with scripts/ and out/ referenced below as W/…); kept here as provenance. It is to be superseded
by review/drafts/D62.md (Proposition D‴.62), whose numbers come from its md; a referee round on D‴.62 is due before the claim is
used elsewhere. -->

# research62: the scaled Gauss–Hermite rules without a Mehler frame

**Result.** D‴.58's conjectured rate t(v) = tanh((ln(σ*/σ) + |ln(v/v_hs)|)/2) holds for every scaled rule without a frame:
moment matching (ρ = 1), ρ > 1 and ρ ≤ −1. Each such rule also has a sharp constant. The proof uses D‴.61 (i) and D‴.61 (iii)(b)
only on the two axes' intervals (−2, 2) and i(−2, 2). D‴.61 proves (iii)(b) there by hand, so this step has **no
computer-assisted part**.

The mechanism is a Bochner (Gaussian) representation of e_K². Without a frame, the representing Gaussian lives on **one axis** of
the ξ-plane: the imaginary axis for v ≥ τ², the real axis for v ≤ τ²σ²/σ*². With a frame it lives off the axes, in D‴.61's Fock
plane.

- **Moment matching.** e_K² is a one-dimensional integral of the square of the Peano law's characteristic function. The Laplace
  mass sits at η = 2iλ^{1/2}, inside i(−2, 2). The result is ‖f* − f_K‖/‖f*‖ ~ 2^{−1/4}(1 + λ)^{3/4}λ^K, with λ = tanh r.
- **ρ ∉ [−1, 1].** e_K² is a two-dimensional Gaussian integral of a product φ_K(α)φ_K(β), with α and β on one axis. It has two
  separate ridges, with rates t(v) and the other branch |h|. Their sector integrals are both positive and add: this proves round
  113's two-axis picture. The outlier's transient constant comes out as B₀ = 0.2383158459, against round 129's fitted 0.238316.
- **The next term c₂ (formal).** It follows from the same representation by a short Laplace computation. That computation
  reproduces round 129's formula, with its signed √q equal to −s_o/s_d.

All of this is checked against exact errors: the referee's 95-digit scans, my own scans at consecutive K up to 200–400 for the
other 11 frameless rules (python-flint arb), and A6's stored 26-digit errors.

Labels used below: **[proved]**, **[formal]** (an asymptotic expansion or identification not proved), **[numerical]**.
Sources: **[snippet]** means web-search snippet only; **[memory]** means recalled, not seen this session.

---

## 1. Proposition draft — D‴.61 (vi), *the rules without a frame: a Bochner representation on one axis*

**Setting.** As in D‴.61 (v) and D‴.58, with σ = 1:

- kernels of width σ, a feature of width σ*, τ² = σ*² − 1, r = ln σ*, λ = tanh r = τ²/(1 + σ*²);
- the Gauss–Hermite rule of N(0, v) for the reweighted law N(0, τ²), with nodes y_k = √v z_k;
- e_K = ‖f* − f_K‖/‖f*‖, and h(y) = (y − 1)/(y + 1);
- φ_K(ξ) = (E − Q_K)[e^{ξz}] = e^{ξ²/2}ε̃_K(ξ) = (K!/(2K)!)ξ^{2K}ĝ_K(ξ), where ĝ_K(ξ) = E e^{ξS} is the transform of D‴.61 (i)'s
  Peano law g_K.

**The parameters.** Put

- c = (1 − v/τ²)/2 and c_b = c − v/2 = (τ² − vσ*²)/(2τ²);
- s_a = c/(1 − c) = h(τ²/v) and s_b = c_b/(1 − c_b) = −h(vσ*²/τ²), both signed.

D‴.61 (v)'s ratio satisfies 1 − ρ = 2(τ² − v)/(vτ²) and 1 + ρ = 2(vσ*² − τ²)/(vτ²). So **|ρ| < 1 ⟺ c·c_b < 0 ⟺ s_a·s_b < 0**,
and:

- the frame rules are exactly those with s_a s_b < 0;
- moment matching is c = s_a = 0 (ρ = 1);
- ρ = −1 is c_b = s_b = 0;
- ρ > 1 means s_a, s_b < 0; ρ < −1 means s_a, s_b > 0.

**The rate and the two ridges.** D‴.58's rate is t(v) = max(|s_a|, |s_b|). Off the cusp v_hs, |s_a| ≠ |s_b|. Let s_d be the one
of larger modulus (so t(v) = |s_d|), s_o the other, and ϖ = s_o/s_d ∈ (−1, 1). Then ϖ < 0 at a frame and 0 ≤ ϖ < 1 without one.
Without a frame, t(v) = h(vσ*²/τ²) if ρ ≥ 1 (the fast branch) and h(τ²/v) if ρ ≤ −1 (the slow branch).

**(a) *(the Bochner representation; every v > 0)*.** Let a and b be independent centred Gaussians with variances c and c_b. A
negative variance c means a = i|c|^{1/2}X with X ~ N(0, 1), and likewise for b. Then

  e_K² = σ*(v/τ²)·E[φ_K(a + b)φ_K(a − b)].

- **Frame (c·c_b < 0).** α = a + b and β = a − b are complex conjugates (a real, b imaginary), and this is D‴.61's Fock integral
  E|φ_K(ξ)|², from its preamble and (v); the constant agrees by (e).
- **No frame.** α and β lie on one axis.
  - **v ≥ τ² (ρ ≥ 1).** With A ~ N(0, |c|) and B ~ N(0, |c_b|) independent,
    e_K² = σ*(v/τ²)(K!/(2K)!)²·E[(B² − A²)^{2K}ĝ_K(i(A + B))ĝ_K(i(A − B))].
    - At moment matching A ≡ 0 and **e_K² = σ*(K!/(2K)!)²·E[B^{4K}ĝ_K(iB)²]** with B ~ N(0, τ²/2). Here ĝ_K(iB) = E cos(BS) is
      the Peano law's characteristic function: real, and at most 1 in modulus.
  - **v ≤ τ²/σ*² (ρ ≤ −1).** With a ~ N(0, c) and b ~ N(0, c_b) real,
    e_K² = σ*(v/τ²)(K!/(2K)!)²·E[(a² − b²)^{2K}ĝ_K(a + b)ĝ_K(a − b)]. The integrand is ≥ 0.

**(b) *(the rate and the sharp constant)*.** For every rule without a frame (ρ ∉ (−1, 1)):

  e_K² = C(v)·t(v)^{2K}(1 + o(1)),  C(v) = √2·σ*·v·(1 − s_d) / [(τ² + v)(τ² + vσ*²)(1 + ϖ)]^{1/2}.

Hence:

- e_K ~ C(v)^{1/2}t(v)^K, lim e_K^{1/K} = t(v), and e_{K−1}/e_K → 1/t(v). This is D‴.58's conjecture: t(v) = h(vσ*²/τ²) for
  ρ ≥ 1 and h(τ²/v) for ρ ≤ −1, with its sharp constant.
- **Moment matching** (s_d = −λ, ϖ = 0): C = (1 + λ)^{3/2}/√2 = 2σ*³/(1 + σ*²)^{3/2}, so
  **‖f* − f_K‖/‖f*‖ ~ 2^{−1/4}(1 + λ)^{3/4}λ^K** with no power of K.
- **The ingredients.** The proof uses D‴.61 (i)'s bounds 0 ≤ φ_K(x) ≤ e^{x²/2} and |φ_K(iw)| ≤ 2. It uses D‴.61 (iii)(b) only on
  the intervals (−2, 2) and i(−2, 2) of the two axes, where D‴.61 proves it without a computer. The Laplace mass sits at
  |η| = 2t(v)^{1/2} < 2.

**(c) *(two ridges whose contributions add)*.** Let ρ ∉ [−1, 1]. In the real coordinates of (a), split the expectation by the two
sectors, |b| > |a| and |a| > |b|. This gives e_K² = J_d(K) + J_o(K), where J_d is the integral over the dominant variable's
sector, with:

- J_d ~ C(v)t(v)^{2K}, and J_o ~ C_o(v)|s_o|^{2K};
- C_o(v) = √2σ*v(1 − s_o)/[(τ² + v)(τ² + vσ*²)(1 + 1/ϖ)]^{1/2} > 0, that is C(v) with d and o exchanged;
- |s_o| = min(|h(τ²/v)|, |h(vσ*²/τ²)|), the formula's other branch.

Both contributions are positive and their sum is e_K². This is round 113's "two arguments' contributions to e_K² add", now a
theorem. (For ρ ≤ −1 the integrand itself is nonnegative.)

*Formally* — granting J_d to relative accuracy o(q^K) — the ratio is e_{K−1}/e_K = R_K^{(d)} + B q^K(1 + B₁/K + …), with

- q = ϖ² = (h_other/t(v))²;
- B = (C_o/C)(1 − q)/(2t(v)q);
- B₁/B₀ = s₁′ − s₁, where s₁′ is (d)'s s₁ with d and o exchanged.

At the outlier (τ/σ = ¼, v = 2v_hs) this gives B = 0.2383158459 and B₁/B₀ = −0.0050412.

**(d) *(the next term; formal)*.** Assume D‴.61 (iii)(b) has a uniform next-order term on the axes' intervals,
log ĝ_K(√Kη) − Kη²/4 = ½log(1 − η²/4) + F₁(η²/4)/K + O(K⁻²), with round 129's

  F₁(y) = y(1 + y)/(1 − y)³ + y/(8(1 − y)²) + y/(8(1 − y)).

Then each ridge's integral is C t^{2K}(1 + s₁/K + O(K⁻²)), with s = s_d (signed) and

  s₁ = 1/16 + 2s(1 + s)/(1 − s)³ + s/(4(1 − s)²) − ϖ(1 + s)²/(8(1 + ϖ)(1 − s)²).

So e_{K−1}/e_K = t(v)^{−1}(1 + c₂/K² + O(K⁻³)), with c₂ = s₁/2 and no K⁻¹ term.

- This is round 129's formula, with ϱ = −s and its signed √q = −ϖ. In the frame region ϖ < 0, and the same algebra applies there.
- At moment matching ϖ = 0 and s₁ = 1/16 − λ/(4(1 + λ)²) − 2λ(1 − λ)/(1 + λ)³. That is −13/125, −11/64 and −187/2000 at τ/σ =
  ½, 1, 2, so c₂ = −13/250, −11/128 and −187/4000.
- On the real branch s = t(v) → 1 gives s₁ ∝ (1 − t)^{−3} → ∞; on the imaginary branch s = −t stays bounded. This is why the slow
  branch converges slowly.

**(e) *(one constant for all rules but the cusp)*.** For ϖ < 0, (b)'s C(v) is D‴.61 (iv)'s constant at every frame off the cusp.
The identity is algebraic: 4(τ² − v)(vσ*² − τ²) = τ⁴(1 − ρ²)v². So e_K ~ C(v)^{1/2}t(v)^K is the sharp asymptotics of every
scaled rule except v_hs.

**Continuity with the frames.** As ρ ↑ 1 (v ↑ τ²), D‴.61's weight e^{−Ax²}, with A = v/(τ² − v) → ∞, concentrates its Fock
integral on the imaginary axis, and the limit is (a)'s one-dimensional integral. C(v) is continuous at ϖ = 0, and likewise at
ρ ↓ −1, where B → ∞ and the integral concentrates on the real axis. So moment matching is the boundary case of the frame formula,
as suggested.

The rate's symmetry v ↦ v_hs²/v is the exact duality (s_a, s_b) ↦ (−s_b, −s_a), which exchanges the real and imaginary axes. It
maps the ρ > 1 rules onto the ρ < −1 ones with the same t and q, though not the same C.

### Proof

**(a)**

- **The error.** By D‴.61 (v)'s proof, or directly from the referee's closed form (arberr.py, a6_error),
  e_K² = σ*(μ − Q)^{⊗2}[κ]. Here κ(y, y′) = e^{−(y − y′)²/4}, μ = N(0, τ²), Q[f] = Σ_k w_k r(y_k)f(y_k), and
  r = dN(0, τ²)/dN(0, v). Expanding gives 1, the cross term and the Gram term of a6_error, using μ⊗μ[κ] = 1/σ* and
  ∫κ(·, y_k)dμ = √2(1 + σ*²)^{−1/2}e^{−y_k²/(2(1 + σ*²))}.
- **The change of variables.** Put y = √v z. Then μ[f] = E[g(Z)f(√vZ)] and Q[f] = Q_K[g·f(√v·)], with
  g(z) = r(√vz) = (v/τ²)^{1/2}e^{cz²}. Hence (μ − Q)^{⊗2}[κ] = (E − Q_K)^{⊗2}[F], with
  F(z, z′) = (v/τ²)exp(c(z + z′)²/2 + c_b(z − z′)²/2), because c(z² + z′²) = c(z + z′)²/2 + c(z − z′)²/2.
  Each of the four terms of (E − Q_K)^{⊗2} is finite, since c < ½ and c_b < ½.
- **The Gaussian representation.** For real p, E e^{pa} = e^{cp²/2} in both cases (c ≥ 0, or a = i|c|^{1/2}X), and likewise for b.
  So F(z, z′) = (v/τ²)E e^{(a + b)z + (a − b)z′}.
- **Fubini.** Fubini–Tonelli applies to each of the four terms:
  - if c, c_b ≤ 0, the integrand has modulus 1;
  - if c, c_b ≥ 0, it is positive with expectation F/(v/τ²), which is E⊗E-integrable;
  - at a frame its modulus is e^{a(z + z′)}, with a real.
  Finally (E − Q_K)^{⊗2}[e^{αz + βz′}] = φ_K(α)φ_K(β), and (iw)^{2K}(iw′)^{2K} = (ww′)^{2K} gives the displayed forms. ∎

**(b)** I write the proof for ρ > 1, then say what changes for moment matching, for ρ < −1 and for ρ = −1.

*Step 1 (the inputs).*

- (i) **Stirling.** (K!/(2K)!)K^K = 2^{−1/2}(e/4)^K(1 + O(1/K)), and (K!/(2K)!)²(4K − 1)!! → 2^{−1/2}. In fact it is
  2^{−1/2}(1 + 1/(16K) + …).
- (ii) **D‴.61 (iii)(b) on the axes.** For each c′ < 2, uniformly on |y| ≤ c′:
  - (−1)^Kφ_K(i√Ky) = 2^{−1/2}e^{KL(y)}(1 + y²/4)^{1/2}(1 + o(1));
  - φ_K(√Ky) = 2^{−1/2}e^{K(L(y) + y²/2)}(1 − y²/4)^{1/2}(1 + o(1)).

  Here L(y) = 1 + log(y²/4) − y²/4. It is ≤ 0, increasing on (0, 2], zero only at |y| = 2, and decreasing beyond.
- (iii) **D‴.61 (i)'s bounds.** |φ_K(iw)| ≤ e^{−w²/2} + Q_K[1] ≤ 2, and 0 ≤ φ_K(x) ≤ e^{x²/2}, since 0 ≤ ε̃_K < 1 on ℝ.

*Step 2 (the exponent).* Write A + B = √K y₁ and A − B = √K y₂. The law of (A + B, A − B) has density ∝ K e^{−KQ(y)}, with
Q(y) = (y₁ + y₂)²/(8|c|) + (y₁ − y₂)²/(8|c_b|). On the box |y_i| ≤ c′, the integrand's exponential rate per K is
Φ(y) = L(y₁) + L(y₂) − Q(y). In the tilted coordinates a = (y₁ + y₂)/2 and b = (y₁ − y₂)/2,

  Φ = 2 + 2log(|a² − b²|/4) − a²/(2|s_a|) − b²/(2|s_b|),

since 1/(2|c|) + ½ = 1/(2|s_a|).

- **For ρ < −1**, the factor e^{K(y₁² + y₂²)/2} of Step 1 (ii) is absorbed. This replaces Q by
  Q̃ = Q − (y₁² + y₂²)/2 = (y₁ + y₂)²/(8c̃) + (y₁ − y₂)²/(8c̃_b), with c̃ = c/(1 − 2c) > 0 and c̃_b = c_b/(1 − 2c_b) > 0. It gives
  the same Φ with s_a = c̃/(1 + c̃) and s_b = c̃_b/(1 + c̃_b). In the tail, the bound φ_K ≤ e^{x²/2} becomes rate 0 after the same
  subtraction. So everything below is literally the ρ > 1 argument with (|c|, |c_b|) replaced by (c̃, c̃_b).

*Lemma A (the maximum).* Φ ≤ 2 log t(v), with equality only at the dominant ridge's points (a, b) = (0, ±2|s_b|^{1/2}) or
(±2|s_a|^{1/2}, 0). On the closed sector {|a| ≥ |b|}, Φ ≤ 2 log|s_a| with equality only at (±2|s_a|^{1/2}, 0); on {|b| ≥ |a|}
likewise with s_b. These points have |y₁| = |y₂| = 2|s|^{1/2} < 2.

*Proof.* With u = a² and w = b², on w ≥ u the function Φ decreases in u. So Φ ≤ 2 + 2log(w/4) − w/(2|s_b|) ≤ 2 log|s_b|, with
equality iff u = 0 and w = 4|s_b|. The other sector is symmetric. ∎

*Lemma T (the tails).* Let L⁺(y) = L(min(|y|, 2)), that is L on |y| ≤ 2 and 0 beyond, and Φ⁺ = L⁺(y₁) + L⁺(y₂) − Q. Then
sup{Φ⁺ : max|y_i| ≥ 2} < 2 log t(v), and the same holds on each closed sector with 2 log|s| in place of 2 log t(v).

*Proof.* Write Q = α(y₁² + y₂²) + 2βy₁y₂, with α = (1/k_a + 1/k_b)/8 > |β| and β = (1/k_a − 1/k_b)/8, where (k_a, k_b) =
(|c|, |c_b|). The sector {|a| ≥ |b|} is {y₁y₂ ≥ 0}.

- (R1) **Both |y_i| ≥ 2.** Here Φ⁺ = −Q. On the sector y₁y₂ ≥ 0, Q ≥ 8(α + β) = 2/k_a: use 2βy₁y₂ ≥ 8β if β ≥ 0, and
  2βy₁y₂ ≥ β(y₁² + y₂²) if β < 0. Likewise Q ≥ 2/k_b on y₁y₂ ≤ 0. And −2/k < 2 log(k/(1 + k)), because log(1 + 1/k) < 1/k.
- (R2) **|y₁| ≥ 2 > |y₂|** (R3 is symmetric). Here Φ⁺ = L(y₂) − Q(y₁, y₂). Q(·, y₂) is convex with minimiser −(β/α)y₂, of
  modulus below 2. So over each half-line {y₁ ≥ 2} or {y₁ ≤ −2} it is least at the endpoint ±2. Since L(±2) = 0, this gives
  Φ⁺ ≤ max_{ε=±1}Φ(2ε, y₂), restricted to the ε of the sector when one is imposed. That is a value of Φ on the boundary of the box
  [−2, 2]², inside the same closed sector. There Φ < 2 log t, or < 2 log|s|, by Lemma A, because the maximisers are interior.
  The supremum over the closed boundary segments is attained, as Φ is upper semicontinuous with Φ = −∞ at y₂ = 0. So it is
  strictly smaller. ∎

*Step 3 (assembly).* Fix c′ ∈ (2t(v)^{1/2}, 2), and let the box be max(|A + B|, |A − B|) = |A| + |B| < c′√K.

- **The box.** By Step 1, φ_K(α)φ_K(β) = M_K·(1 + o(1)) uniformly on the box, where
  M_K = (K!/(2K)!)²(B² − A²)^{2K}e^{−(A² + B²)/2}m_K ≥ 0 and m_K = (1 + (A + B)²/4K)^{1/2}(1 + (A − B)²/4K)^{1/2}.
- **The tilt.** E[e^{−(A² + B²)/2}G] = N·E[G(σ_AX, σ_BY)], with X, Y iid N(0, 1), σ² = |s|, and N = ((1 − c)(1 − c_b))^{−1/2}.
- **Lemma B (the ridge integral).** Size-bias Y by Y^{4K}: E[Y^{4K}G(Y)] = (4K − 1)!!·E G(Ỹ), with Ỹ² ~ χ²_{4K+1} and
  Ỹ²/4K → 1 in probability. On the sector {|b| > |a|} = {Y² > pX²}, p = σ_A²/σ_B² = ϖ:
  - the factor (1 − pX²/Ỹ²)^{2K} lies in [0, 1] and tends to e^{−pX²/2} in probability;
  - m_K lies in [1, 1 + c′²/4] on the box and tends to 1 + σ_B² = 1 − s_d;
  - the box's indicator tends to 1, since 2σ_B < c′.

  Dominated convergence then gives E[M_K; box, |B| > |A|] = (K!/(2K)!)²(4K − 1)!!·N|s_d|^{2K}(1 − s_d)(1 + ϖ)^{−1/2}(1 + o(1)).
  On the other sector it gives the same with o and 1/ϖ: this is (c)'s C_o. The tilt and the size-biasing are exact; only the
  limit uses dominated convergence.
- **The tails.** Off the box, Step 1 (ii)–(iii) give |φ_K(α)φ_K(β)| ≤ C·e^{K(L̃(y₁) + L̃(y₂))}, with L̃ = L on |y| < c′ and 0
  beyond. So E[|φφ|; off the box] ≤ C′K·e^{(K−1)S}·∫e^{Φ̃}. Here Φ̃ = L̃ + L̃ − Q is coercive with ∫e^{Φ̃} < ∞, and
  S = sup_{max|y_i| ≥ c′}Φ̃.
  - Since L is increasing on (0, 2], L̃ ≤ L⁺ − L(c′). So S ≤ sup_{max|y_i| ≥ c′}Φ⁺ − 2L(c′).
  - As c′ ↑ 2, −L(c′) → 0. The suprema over the decreasing closed sets tend to the supremum over {max|y_i| ≥ 2}, because Φ⁺ is
    upper semicontinuous and coercive. By Lemma T that supremum is below 2 log t.
  - So for c′ near 2 the part off the box is O(e^{K(2log t − δ)}) = o(t^{2K}).
- **Collecting.** e_K² = σ*(v/τ²)N·2^{−1/2}t^{2K}(1 − s_d)(1 + ϖ)^{−1/2}(1 + o(1)), and σ*(v/τ²)N/√2 =
  √2σ*v/((τ² + v)(τ² + vσ*²))^{1/2}, since 1 − c = (τ² + v)/(2τ²) and 1 − c_b = (τ² + vσ*²)/(2τ²). This is (b).

**The other cases.**

- **Moment matching.** The integral is one-dimensional. The tail bound is −c′²/(2|c_b|) → −2/|c_b| < 2 log λ = 2 log σ_B² (that is
  Lemma T's R1), and the Gamma integral is exact: ∫W^{4K}e^{−γW²}dW = Γ(2K + ½)γ^{−2K−½}, with γ = 1/(2λ). The tail condition
  in closed form is c′² > 4 ln(1 + u)/u with u = 2/τ². It is always below 4, and above the mass point's 4λ = 4/(1 + u).
- **ρ = −1.** Likewise one-dimensional, on the real axis.
- **ρ < −1.** As said in Step 2. In addition, φ_K ≥ 0 makes the integrand nonnegative, so no o(1) can change a sign. ∎

**(c)** This is (b)'s argument run on each sector, with Lemma A, T and B in their sector forms. ∎ The ratio's transient is formal.
Writing e_K² = C t^{2K}(1 + s₁/K + …) + C_o|s_o|^{2K}(1 + s₁′/K + …) and expanding gives B and B₁/B₀ = s₁′ − s₁. A proof of the
transient would need J_d to relative accuracy o(q^K), beyond any power-series statement.

**(d)** [formal]

- **Set-up.** On the dominant ridge write Y² = 4K(1 + ζ). Then Eζ = 1/(4K) + O(K⁻²) and Eζ² = 1/(2K) + O(K⁻²). Take
  W = X² ~ χ²₁ tilted by e^{−ϖW/2}, with E W = ω = 1/(1 + ϖ) and E W² = 3ω². Write α², β² = s(Y² + ϖX²) ± 2ab, with
  (ab)² = s²ϖX²Y².
- **The integrand.** It is s^{2K}Y^{4K}(1 − ϖX²/Y²)^{2K}·exp(f(α²/4K) + f(β²/4K)), where f(u) = ½log(1 − u) + F₁(u)/K.
- **Expansion.** Expand to O(1/K), writing α²/4K = s + δ_±, δ_± = sζ + (sϖX² ± 2ab)/(4K):
  - the ±ab terms cancel at first order;
  - their squares give ½f″(s)·(2ab/4K)²·2 = f″(s)s²ϖX²(Y²/4K)/K ≈ f″(s)s²ϖW/K;
  - log(1 − ϖX²/Y²)^{2K} = −(ϖW/2)(1 − ζ + ζ²) − ϖ²W²/(16K) + …;
  - the ζ-terms are averaged to second order.

  With f′ = −1/(2(1 − s)) and f″ = −1/(2(1 − s)²), and adding Stirling's 1/16, this gives
  s₁ = 1/16 − ϖω/8 + sf′/2 + s²f″/2 + s²f′² + ϖωsf′ + ϖωs²f″ + 2F₁(s).
- **Simplification.** Here s²f″/2 + s²f′² = 0, and the expression simplifies to the displayed formula.
- **Checks.** At ϖ = 0, s = −λ, the one-dimensional Gamma computation gives the same. Substituting s = −ϱ and ϖ = −√q gives round
  129's s₁ term by term.
- **Why it extends.** Every term is a polynomial in Gaussian moments times rational functions of (s, ϖ), so the formal expansion is
  one analytic expression across ϖ = 0. This is why round 129's frame-derived c₂ also fits the frameless rules. ∎

**Status.**

- (a) proved (elementary).
- (b) proved, given D‴.61 (i) and (iii)(b) on the axes' intervals, both of which D‴.61 proves by hand. No computer-assisted step
  is involved: none of D‴.61's branch and bound, V or VA is needed, because the integrand never leaves the axes.
- (c) the two sector asymptotics proved; the transient in the ratio formal.
- (d) formal. It needs (iii)(b) to next order uniformly on the axes' intervals, and round 129's F₁. F₁ comes from the cumulants'
  1/K coefficients d_j = −c_j·j(8j² + j + 1)/4, exact for j ≤ 22 (round 129).
- (e) proved (algebra).

Numerically everything is checked to the resolution of the exact errors (§3).

---

## 2. What remains open, and why

1. **c₂ and every higher order, rigorously.** Two things are needed.
   - A uniform next-order statement for ĝ_K on the axes' intervals, assumption (F) in (d). D‴.61 (iii) gets its limit by two
     Laplace steps:
     - (a) the corner ratio on the circle |ζ| = |μ|, an analytic phase with a unique maximum;
     - (b) Duhamel's s-integral at s = ½, where the parts beyond the unit disc are exponentially small.

     Carrying both to O(K⁻²) looks routine, but it is not done here. Then F₁ must be identified: either by proving round 129's
     cumulant formula d_j for all j, or by computing the next order of the two Laplace steps directly.
   - The next order of the ridge integrals, which is standard once (F) holds. The tails are exponentially small.
2. **The transient Bq^K in the ratio, as a theorem.** (c) proves that the two sector integrals have their own asymptotics and add.
   Seeing the subdominant one in e_{K−1}/e_K requires J_d to relative accuracy o(q^K). That is a resurgence-type statement about
   the 1/K series. Round 129 (f) found that series Gevrey-1, with its Borel singularity at the other axis's exponent.
3. **Nothing new for the frame rules.** Their representing Gaussian is genuinely off the axes. No Gaussian measure on (ℝ or iℝ)²
   has the required covariance, because the quadratic form's eigenvalues c and c_b have opposite signs. So D‴.61 (iv)'s use of
   the full disc, computer-assisted beyond r = 0.0966, is not avoidable by this route.
4. **Non-asymptotic bounds.** The representation gives cheap unconditional bounds: for example e_K² ≤ 4σ*v/τ² for ρ ≥ 1, since
   |φ_K(iw)| ≤ 2, and at moment matching e_K² ≤ σ*(K!/(2K)!)²(4K − 1)!!(τ²/2)^{2K} ~ 2^{−1/2}σ*(τ²/2)^{2K} from |ĝ_K(iB)| ≤ 1. These are not sharp: the
   rate τ²/2 exceeds λ. A K-explicit sharp bound would need a K-explicit version of (iii)(b) on the axes.
5. **Not needed, and not addressed.** ĝ_K beyond |η| = 2 on the imaginary axis (D‴.61's alias heuristic λ₁). The Laplace mass of
   every frameless rule sits at |η| = 2t(v)^{1/2} < 2, and outside the box the trivial bound |φ_K(iw)| ≤ 2 suffices (Lemma T).
   So the alias rate plays no role, and an upper bound beyond 2i was not needed.

---

## 3. Numerical checks

Scripts are in W/scripts/ and outputs in W/out/, with W = this note's directory. Everything ran under `OMP_NUM_THREADS=1 nice -n 15`
with the repository's venv (python-flint 0.9; the venv has no mpmath), one heavy process at a time.

**(R) The representation itself** (`check_repr.py` → `out/check_repr.txt`).

| case | rules | method | agreement with `err2`'s closed form |
|---|---|---|---|
| moment matching, 1-D | τ/σ = 1, K = 6; τ/σ = 2, K = 10; τ/σ = ¼, K = 5 | rigorous `acb.integral` | relative ≤ 2.4·10⁻⁸⁹ |
| ρ > 1 and ρ < −1, 2-D | six rules, K = 6–8 | 120² and 160² Gauss–Hermite product rules | relative ≤ 10⁻⁴⁶ |

In the 2-D rows the change from M = 120 to 160 is ≤ 7·10⁻³⁴. The 2-D quadrature is not rigorous; it is a check of the normalisation
and of the Fubini step.

**(C) The constant C(v), (b)** (`check_constants.py`, `fit_all.py` → `out/check_constants.txt`, `out/fit_all.txt`). Fit
e_K²/(C t^{2K}) to L + Σ s_jK^{−j} (+ q^K(g₀ + …) where a second ridge exists), on consecutive K, with the limit L free. The fits
are linear least squares in arb, with the data's midpoints taken as exact.

- **Sources.** The referee's 95-digit scans are scan_1_mm, scan_2_mm and scan_q_2vhs; my own scans, with the same formula and arb
  radii ≤ 10⁻⁷¹ relative, are out/scan_τ_rule.txt.
- **Frame-rule control.** At the five scanned frame rules, C(v) equals D‴.61 (iv)'s constant identically: difference 0 at 2400 bits.

| τ/σ | rule | ρ | kind | t(v) | C(v) | K | \|L − 1\| (3 limit-free fits) | formal s₁ | fitted s₁ (L = 1) | rel. diff |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | moment matching | 1.0000 | MM | 0.03030303030 | 0.7394902765824 | 2…200 (own) | 2.8e-32 … 7.9e-24 | 0.00162833299410 | 0.00162833299410 | 2.3e-27 |
| 0.5 | moment matching | 1.0000 | MM | 0.1111111111 | 0.8281733249999 | 2…300 (own) | 6.7e-33 … 2.2e-24 | -0.104000000000 | -0.104000000000 | 9.1e-30 |
| 1 | moment matching | 1.0000 | MM | 0.3333333333 | 1.088662107904 | 2…400 (referee) | 1.8e-35 … 3.6e-26 | -0.171875000000 | -0.171875000000 | 2.1e-32 |
| 2 | moment matching | 1.0000 | MM | 0.6666666667 | 1.521451548625 | 2…400 (referee) | 1.3e-37 … 1.9e-27 | -0.0935000000000 | -0.0935000000000 | 3.2e-34 |
| 0.25 | 2·v_hs | 16.508 | rho>1 | 0.3467367306 | 0.9156974364610 | 2…400 (referee) | 1.3e-25 … 7.2e-21 | -0.184873969084 | -0.184873969084 | 8.5e-19 |
| 0.25 | 1.25·v_hs | 6.6121 | rho>1 | 0.1260538456 | 0.6663247689700 | 2…240 (own) | 4.2e-30 … 6.2e-24 | -0.149237084205 | -0.149237084205 | 2.0e-23 |
| 0.5 | 1.25·v_hs | 1.8446 | rho>1 | 0.1658124885 | 0.7911996228148 | 2…300 (own) | 3.7e-34 … 1.7e-26 | -0.158691259843 | -0.158691259843 | 4.1e-31 |
| 0.5 | 2·v_hs | 4.5279 | rho>1 | 0.3819660113 | 0.9862320002593 | 2…360 (own) | 1.9e-35 … 3.2e-26 | -0.177022048219 | -0.177022048219 | 9.9e-33 |
| 1 | 2·v_hs | 1.5858 | rho>1 | 0.4775922501 | 1.179105345634 | 2…400 (own) | 1.6e-37 … 4.2e-27 | -0.150996760566 | -0.150996760566 | 2.4e-34 |
| 0.25 | 0.5·v_hs | -32.970 | rho<-1 | 0.3467367306 | 0.2220892500962 | 2…300 (own) | 3.0e-12 … 9.4e-11 | 3.36076057505 | 3.36076057412 | 2.7e-10 |
| 0.25 | 0.8·v_hs | -8.2311 | rho<-1 | 0.1260538456 | 0.4137151854575 | 2…240 (own) | 5.9e-21 … 4.4e-16 | 0.439276729183 | 0.439276729183 | 2.4e-18 |
| 0.5 | 0.5·v_hs | -8.8885 | rho<-1 | 0.3819660113 | 0.2205281794165 | 2…360 (own) | 9.6e-12 … 8.5e-10 | 4.51872070453 | 4.51872070058 | 8.7e-10 |
| 0.5 | 0.8·v_hs | -2.1803 | rho<-1 | 0.1658124885 | 0.4529090919388 | 2…300 (own) | 8.3e-19 … 3.2e-15 | 0.726673566907 | 0.726673566907 | 2.6e-16 |
| 1 | 0.5·v_hs | -2.6569 | rho<-1 | 0.4775922501 | 0.2084383464078 | 2…400 (own) | 6.4e-09 … 7.6e-08 | 10.1351971970 | 10.1351936620 | 3.5e-07 |

- **Columns.** "K" is the range of consecutive K, with its source. "\|L − 1\|" spans three limit-free fits: [K_max/2, K_max] with 10
  powers, [K_max/2 − 40, K_max] with 14, and [max(60, K_max/4), K_max] with 18; the q^K-terms are included where resolvable. The
  smallest of the three is the agreement of the limit with the proved C(v), the largest is the most conservative. "fitted s₁" is
  from the third fit with L = 1 fixed.
- **Moment matching and ρ > 1.** e_K²/(C(v)t^{2K}) → 1 to 21–37 digits. The formal s₁ is matched to 19–34 digits.
- **ρ < −1, the slow branch.** Limits within 10⁻¹⁵–10⁻²¹ at the two rules with small t(v), and 10⁻¹⁰–10⁻¹² at the two with
  t(v) ≈ 0.35–0.38. At τ/σ = 1, 0.5·v_hs (t(v) = 0.478, s₁ = 10.1) they are within 6·10⁻⁹–8·10⁻⁸. There s₁ is matched to
  3.5·10⁻⁷, and to 10⁻¹⁰–10⁻¹⁸ at the others. The slow branch's large higher coefficients — fitted s₂ = −30.7, −54.2 and −263 at the three rules with t(v) = 0.35–0.48,
  against −0.94 and −2.10 at the two with small t(v) — are what (d)'s (1 − t)^{−3} predicts.
- **Precision of my scans.** Their largest arb relative radius is 5.6·10⁻⁷², at τ/σ = 1, 0.5·v_hs, K = 400; elsewhere ≤ 3·10⁻⁷⁹.
- **Cross-check against A6.** My scans reproduce A6's stored 26-digit strings, for example τ/σ = 1, 0.5·v_hs at K = 19, 20.

**(S) The formal s₁, (d).** The last three columns of the table above give it: the L = 1 fit's s₁ against
s₁ = 1/16 + 2s(1 + s)/(1 − s)³ + s/(4(1 − s)²) − ϖ(1 + s)²/(8(1 + ϖ)(1 − s)²). My formula equals round 129's with ϱ = −s and
√q = −ϖ to 3·10⁻⁵⁸ at 20 random (s, ϖ), and equals the unsimplified expansion of (d)'s proof as well.

At moment matching s₁ is exactly rational: −11/64 and −187/2000 to 32–34 digits, and −13/125 from A6's raw. The next coefficients
fitted are s₂ = 93/8192 (= 0.0113525390625) at τ/σ = 1 and −0.002692875 at τ/σ = 2.

**(B) The second ridge, (c)** (`fit_transient.py`, `fit_all.py`).

| τ/σ | rule | kind | q | C_o/C_d (model) | B₀ (model) | s₁′ (formal) | fits: \|g₀/(C_o/C_d) − 1\|, g₁/g₀ |
|---|---|---|---|---|---|---|---|
| 0.25 | 2·v_hs | rho>1 | 0.85062954 (q^400 = 7.9e-29) | 0.941149256337 | 0.23831584594 | -0.1899151843 | [100,400] n=18,m=3: 2.1e-10, -0.18991525; [50,400] n=22,m=4: 8.5e-13, -0.18991518 |
| 0.25 | 0.5·v_hs | rho<-1 | 0.85062954 (q^300 = 8.4e-22) | 0.999971084858 | 0.25321058631 | 2.672647476 | [75,300] n=18,m=3: 2.5e-04, 2.6098237; [37,300] n=22,m=4: 2.3e-04, 2.6305507 |
| 0.25 | 0.8·v_hs | rho<-1 | 0.58142937 (q^240 = 3.0e-57) | 0.903132182413 | 2.5789149274 | 0.2729516532 | [60,240] n=18,m=3: 6.9e-05, 0.28199941; [30,240] n=22,m=4: 3.1e-07, 0.27297881 |
| 0.25 | 1.25·v_hs | rho>1 | 0.58142937 (q^240 = 3.0e-57) | 0.850006759918 | 2.4272140493 | -0.1376659156 | [60,240] n=18,m=3: 1.8e-09, -0.13766623; [30,240] n=22,m=4: 1.4e-10, -0.13766594 |
| 0.5 | 0.5·v_hs | rho<-1 | 0.54839384 (q^360 = 1.2e-94) | 0.998538787999 | 1.0764094747 | 1.937925580 | [90,360] n=18,m=3: 2.0e+13, -230.79466; [45,360] n=22,m=4: 2.1e+07, -144.82733 |
| 0.5 | 0.8·v_hs | rho<-1 | 0.11295718 (q^300 = 7.5e-285) | 0.656238713893 | 15.539797559 | 0.1009336240 | [31,300] n=18,m=3: 1.7e+19, -90.633324; [20,300] n=22,m=4: 5.0e+11, -69.453346 |
| 0.5 | 1.25·v_hs | rho>1 | 0.11295718 (q^300 = 7.5e-285) | 0.524990971115 | 12.431838047 | -0.1142879188 | [31,300] n=18,m=3: 4.7e+05, -90.294855; [20,300] n=22,m=4: 6.3e-04, -0.15557695 |
| 0.5 | 2·v_hs | rho>1 | 0.54839384 (q^360 = 1.2e-94) | 0.798831030399 | 0.86112757974 | -0.1950736387 | [90,360] n=18,m=3: 9.9e-10, -0.19507396; [45,360] n=22,m=4: 2.6e-10, -0.19507369 |
| 1 | 0.5·v_hs | rho<-1 | 0.12905734 (q^400 = 0.0e+00) | 0.950474300152 | 6.7152155473 | 0.6481812161 | [33,400] n=18,m=3: 1.7e+28, -96.306817; [20,400] n=22,m=4: 2.3e+20, -67.065722 |
| 1 | 2·v_hs | rho>1 | 0.12905734 (q^400 = 0.0e+00) | 0.475237150076 | 3.3576077736 | -0.1915080866 | [33,400] n=18,m=3: 1.6e+04, -96.038712; [20,400] n=22,m=4: 7.1e-07, -0.19155479 |

The fit is e_K²/(C t^{2K}) − 1 = Σ_{j≤n}s_jK^{−j} + q^K Σ_{j≤m}g_jK^{−j}. Its g₀ and g₁/g₀ are compared with the model's C_o/C_d
and the formal s₁′.

- **Resolved.** The second ridge's leading coefficient C_o/C_d — proved, (c) — and its formal 1/K coefficient s₁′ are matched at
  five rules:

  | rule | C_o/C_d agreement | s₁′ agreement |
  |---|---|---|
  | the outlier (¼, 2·v_hs) | 8.5·10⁻¹³ | 2·10⁻⁸ |
  | (¼, 1.25·v_hs) | 1.4·10⁻¹⁰ | 2·10⁻⁷ |
  | (½, 2·v_hs) | 2.6·10⁻¹⁰ | 3·10⁻⁷ |
  | (1, 2·v_hs) | 7·10⁻⁷ | 2·10⁻⁴ |
  | (¼, 0.8·v_hs), slow branch | 3·10⁻⁷ | 10⁻⁴ |

- **Roughly resolved.** Within 2–6·10⁻⁴ at (¼, 0.5·v_hs) and (½, 1.25·v_hs).
- **Not resolvable** at three slow-branch rules: the fitted g₀ is off by orders of magnitude. That is a limit of the data and the
  fit, not a failure of the model.
  - (½, 0.5·v_hs): q = 0.55, but the power series' coefficients are large (s₂ = −54).
  - (½, 0.8·v_hs): q = 0.11, so q^K falls below the series' truncation error by K ≈ 20.
  - (1, 0.5·v_hs): both, with s₂ = −263 and q = 0.13.
- **The outlier's B₀.** The model's 0.2383158459 is round 129's fitted 0.238316–0.238317 (§6 there). Its B₁/B₀ = s₁′ − s₁ =
  −0.00504 is round 129's measured −0.005 … −0.006.

**(A6) A6's stored 26-digit errors** (`check_raw.py` → `out/check_raw.txt`; 12 points per rule, K − 1 and K at K = 20 … 160).
The raw was read from the working tree, read-only. Another session had uncommitted changes to width_mismatch_raw.npz when I read
it. The A6 records I used agree with my independent scans where they overlap, for example τ/σ = 1, 0.5·v_hs at K = 19, 20 to all
26 digits.

- **The 14 frameless rules, limit free, 7 powers.**
  - Moment matching: |L − 1| ≤ 1.9·10⁻¹⁴, with s₁ equal to the formal value to 8 digits at all four.
  - ρ > 1: |L − 1| from 4.7·10⁻¹⁵ to 9·10⁻¹⁰ at four rules; 2.2·10⁻⁶ at the outlier, whose transient is not modelled in this fit.
  - ρ < −1, the slow branch: |L − 1| from 9.6·10⁻¹⁰ to 3.3·10⁻⁴. These are pre-asymptotic at K ≤ 160, s₁ being as large as 10.1, and (C)'s
    consecutive-K scans resolve them.
- **The 22 frame rules.** C(v) with ϖ < 0 is identically D‴.61 (iv)'s constant (an algebraic identity; equal to 0 at 2400 bits at
  the five scanned frame rules), so A6's frame-rule comparison is D‴.61's (R). The cusp (v = v_hs) is excluded, because the formula
  is singular there (1 + ϖ → 0).

**(T) The tail inequality, Lemma T** (`check_tail.py` → `out/check_tail.txt`; floating point, a sanity check of a proved lemma).

- At A6's 10 two-dimensional frameless rules, sup over {max|y_i| ≥ 2} of Φ⁺ lies 0.48–5.37 below 2 log t(v).
- The global maximum is at the interior ridge.
- On a sweep τ/σ ∈ [0.05, 20], v/v_hs ∈ [0.02, 50] (frameless), the margin stays positive. It is smallest, 1.4·10⁻⁵, as t(v) → 1.

**The model's numbers for A6's 14 frameless rules** (`table.py` → `out/table.txt`). R₁₆₀ is A6's e₁₆₀²/(C t^{320}); 160(R₁₆₀ − 1)
is to be compared with s₁.

| τ/σ | rule | ρ | kind | t(v) | C(v) | s₁ | c₂ | \|s_o\| | q | C_o/C_d | B₀ | B₁/B₀ | A6: R₁₆₀ − 1 | 160(R₁₆₀ − 1) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | moment matching | 1.00000 | MM | 0.03030303030 | 0.739490276582 | 0.001628332994 | 0.0008141664970 | — | — | — | — | — | 1.2370e-5 | 0.00197926 |
| 0.25 | 0.5*v_hs | -32.9697 | rho<-1 | 0.3467367306 | 0.222089250096 | 3.360760575 | 1.680380288 | 0.31979383 | 0.85062954 | 0.9999710849 | 0.2532105863 | -0.688113 | 0.019954 | 3.19265 |
| 0.25 | 0.8*v_hs | -8.23106 | rho<-1 | 0.1260538456 | 0.413715185457 | 0.4392767292 | 0.2196383646 | 0.096117968 | 0.58142937 | 0.9031321824 | 2.578914927 | -0.166325 | 0.0027098 | 0.433562 |
| 0.25 | 1.25*v_hs | 6.61212 | rho>1 | 0.1260538456 | 0.666324768970 | -0.1492370842 | -0.07461854210 | 0.096117968 | 0.58142937 | 0.8500067599 | 2.427214049 | 0.0115712 | -0.00092955 | -0.148727 |
| 0.25 | 2*v_hs | 16.5076 | rho>1 | 0.3467367306 | 0.915697436461 | -0.1848739691 | -0.09243698454 | 0.31979383 | 0.85062954 | 0.9411492563 | 0.2383158459 | -0.00504122 | -0.0011532 | -0.184510 |
| 0.5 | moment matching | 1.00000 | MM | 0.1111111111 | 0.828173325000 | -0.1040000000 | -0.05200000000 | — | — | — | — | — | -0.00064672 | -0.103475 |
| 0.5 | 0.5*v_hs | -8.88854 | rho<-1 | 0.3819660113 | 0.220528179417 | 4.518720705 | 2.259360352 | 0.28285965 | 0.54839384 | 0.9985387880 | 1.076409475 | -2.58080 | 0.026457 | 4.23312 |
| 0.5 | 0.8*v_hs | -2.18034 | rho<-1 | 0.1658124885 | 0.452909091939 | 0.7266735669 | 0.3633367835 | 0.055728090 | 0.11295718 | 0.6562387139 | 15.53979756 | -0.625740 | 0.0044626 | 0.714012 |
| 0.5 | 1.25*v_hs | 1.84458 | rho>1 | 0.1658124885 | 0.791199622815 | -0.1586912598 | -0.07934562992 | 0.055728090 | 0.11295718 | 0.5249909711 | 12.43183805 | 0.0444033 | -0.00098903 | -0.158246 |
| 0.5 | 2*v_hs | 4.52786 | rho>1 | 0.3819660113 | 0.986232000259 | -0.1770220482 | -0.08851102411 | 0.28285965 | 0.54839384 | 0.7988310304 | 0.8611275797 | -0.0180516 | -0.0011045 | -0.176716 |
| 1 | moment matching | 1.00000 | MM | 0.3333333333 | 1.08866210790 | -0.1718750000 | -0.08593750000 | — | — | — | — | — | -0.0010738 | -0.171803 |
| 1 | 0.5*v_hs | -2.65685 | rho<-1 | 0.4775922501 | 0.208438346408 | 10.13519720 | 5.067598599 | 0.17157288 | 0.12905734 | 0.9504743002 | 6.715215547 | -9.48702 | 0.055943 | 8.95083 |
| 1 | 2*v_hs | 1.58579 | rho>1 | 0.4775922501 | 1.17910534563 | -0.1509967606 | -0.07549838028 | 0.17157288 | 0.12905734 | 0.4752371501 | 3.357607774 | -0.0405113 | -0.00094277 | -0.150843 |
| 2 | moment matching | 1.00000 | MM | 0.6666666667 | 1.52145154863 | -0.09350000000 | -0.04675000000 | — | — | — | — | — | -0.00058447 | -0.0935156 |

---

## 4. How it would enter D‴.58, D‴.61 and Paper 18

- **D‴.61.** Add a part (vi) as drafted in §1, or make it a short Proposition D‴.62 resting on D‴.61 (i) and (iii)(b).
  - (iv)'s last sentence, "Moment matching (ρ = 1) and the rules with |ρ| ≥ 1 are outside this frame picture", becomes: they are
    the frame picture's degenerate cases. The representing Gaussian, off the axes for a frame, collapses onto one axis (ρ = ±1) or
    becomes a real Gaussian on one axis's square (|ρ| > 1).
  - C(v) is one formula for all of them (e).
  - Status: (vi) is proved without any computer-assisted step.
- **D‴.58 (ii).**
  - "conjectural for the others, moment matching (ρ = 1) and |ρ| > 1" becomes proved, by D‴.61 (vi), for every v.
  - Moment matching's λ, so far "numerical (A3)", becomes proved, with ‖f* − f_K‖/‖f*‖ ~ 2^{−1/4}(1 + λ)^{3/4}λ^K.
  - The half-squeeze rule's gain over moment matching, "(λ/t)^K up to a power of K", becomes exactly
    (λ/t)^K(1 + λ)^{3/4}2^{−1/4}/[(1 − t²)^{1/4}(πK)^{1/4}₂F₁(−½, −½; 1; t²)^{1/2}](1 + o(1)).
- **D‴.58's Status and A6.**
  - "at moment matching and the 10 rules with |ρ| > 1 the curve rests on A6's numbers alone" is withdrawn. All 36 of A6's rules now
    have proved limits, the 14 frameless ones without computer assistance.
  - The fits' brackets and biases on the slow branch (K^{−1.52…−1.98}) are calibrated by (d)'s s₁. The outlier's "geometric
    transient of the formula's other branch" becomes (c): two sector integrals that add, with B₀ = 0.2383158459 and
    B₁/B₀ = s₁′ − s₁ predicted and matched.
- **Paper 18.**
  - Claim: "every scaled Gauss–Hermite rule … falls at t(v) with sharp constants" — no longer "conjectural for the rules without a
    frame". The constant is C(v)^{1/2} off the cusp.
  - "Missing": delete "the rules with |ρ| ≥ 1 … in the frame picture".
  - "Novelty": add the axis representation, the moment-matching constant 2^{−1/4}(1 + λ)^{3/4}, and the two-ridge decomposition.
  - The rate t(v) for every v is then computer-assisted only for frame rules with r ≥ 0.0966.
  - New "Missing" items: the next order c₂ rigorously, which needs (F); and the transient as a theorem.
  - The existing "Missing" item "the real branch's slower convergence near the branch point η = 2" gets a formal explanation:
    s₁ ∝ 2t(1 + t)/(1 − t)³ on the real branch, against a bounded s₁ on the imaginary one. For example s₁ = 43.98 at the frame
    rule τ/σ = 2, 0.5·v_hs, fitted 43.86 to K = 400.
- **Literature.**
  - The MMD of a shift-invariant kernel as a weighted L² distance of characteristic functions is classical
    (Sriperumbudur–Gretton–Fukumizu–Schölkopf–Lanckriet 2010) [memory]. The Bochner step of (a) is that identity, specialised.
  - Bochner-based quadrature for kernel features (Dao–De Sa–Ré, arXiv:1709.02605) [snippet] is the opposite direction.
  - The search found no exact rate or sharp constant for the Gauss–Hermite rule's worst-case error in a Gaussian-kernel space
    [snippets: Karvonen–Oates–Girolami arXiv:2004.12654 and its Newcastle eprint record — two-sided bounds e^{−c₁N}N^{1/4} ≲ · ≲
    e^{−c₂N}N^{−1/4}; Karvonen–Särkkä arXiv:1803.09532].
