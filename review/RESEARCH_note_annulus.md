<!-- Research note (2026-11-05): the valley condition of D‴.61 (iii)(b) / Paper 18 Definition 4.5 on the last annulus
1.995 ≤ |η| < 2 — and on the whole annulus 1 < |η| < 2 — proved by hand, along a new path. Written by a research subagent; its
scripts and outputs are in review/research_annulus/. Refereed in round 148 (review/REVIEW_round148_D3eq.md, D3eq): the proof stands,
minor and editorial findings applied (see the Review history at the end). -->

# The valley condition on the last annulus, by hand — and on the whole annulus 1 < |η| < 2

Research note, 2026-11-05. The scripts and their outputs are in `review/research_annulus/` (Files, at the end). When the note
was written the repository was read, not modified (D61.md, PARADIGM.md, the manuscripts and the experiments untouched); the
changes it suggests were made after round 148 (section 5). Notation as in D‴.61 (iii)(b) and Paper 18 §4: η = Re^{iθ},
ε = 2 − R, ω = R/2, F₀ = F_η(η/2), M(a) = a² (|a| ≤ 1), 2|a| − 1 (|a| > 1), H(r) = (r² − 1)/2 − log r, g_η(s) = log(s(1 − s)) +
(Re η²/2)((1 − s)² + s²).

**Main result.** For every η with 1 < |η| < 2 the valley condition holds along the three-segment path

  Π_η:  Ch = [0, η/|η|]  (the chord, to the unit circle),
        V  = [η/|η|, ν₁]  (vertical: Re ν fixed),
        H  = [ν₁, η]      (horizontal: Im ν = Im η),      ν₁ = Re(η)/|η| + i Im η,

and the maximum of F_η − F₀ on V ∪ H is attained at η/|η| and equals

  m(η) = g_η(1/R) − g_η(½) = log(1 − ε²/R²) + (ε²/4)cos 2θ ≤ −ε²(4 − R²cos 2θ)/(4R²) < 0.

The proof is by hand, with no computation: on the unit circle and beyond it F_η = log|η − ν| + ½Re(η − ν)² + M(Re ν) − ½ (that is
the definition of Λ⁺), and this function is strictly monotone along V and along H. In particular the last annulus 1.995 ≤ |η| < 2
(part VA's range) is covered by hand, and so is everything part D (the branch and bound to 1.8927) and part V certify: with the
straight path on |η| ≤ 1, **Ω′ is the whole open disc |η| < 2 by hand, and Paper 18's Theorem 2 (D‴.61 (iii)(b)) needs no
computer-assisted step.** The route asked for — a remainder bound uniform on VA's (u, q) square — was not needed and was not
carried out (section 3).

Status labels: **[proved]** by hand, complete; **[numerical]** numerical evidence only.

---

## 0. What had to be shown

*The condition* (Paper 18, Definition 4.5; D‴.61 (iii)(b), condition (2)). For |η| < 2, F_η(ν) = log|η − ν| + ½Re(η − ν)² +
log|ν| + ½Re ν² + Λ(ν), with Λ = 0 on |ν| ≤ 1 and Λ = Λ⁺ = M(Re ν) − ½ − log|ν| − ½Re ν² beyond (Paper 18, Lemma 4.4: the bound of
(1/K)log|c_K(√Kν)|, "a bound beyond the unit disc"). η satisfies the valley condition if some piecewise C¹ path from 0 to η contains the chord {(η/2)(1 + z) :
|z| ≤ x₀} (x₀ ∈ (0, ½], the chord in the open unit disc) and F_η(ν) < F₀ on it except at ν = η/2. In D61's saddle coordinates σ =
(η/2)(1 + z) the bound is E(z; η) = F_η(σ) − F₀ = log|1 − z²| + Re(η²z²/4) + C(σ), C = 0 for |σ| ≤ 1 and H(|σ|) − (|Re σ| − 1)₊²
beyond; D61's form (2): it suffices that a path from b = (1 + t)η/2, some 0 < t < 2/|η| − 1, to η keeps E < 0, the chord
[(1 − t)η/2, (1 + t)η/2] and 0 → a lying on the straight segment. By the symmetries F_{η̄}(ν̄) = F_{−η}(−ν) = F_η(ν) the sector
0 ≤ θ ≤ π/2 suffices. Part VA certifies this on 2 − ε₀ ≤ |η| < 2 (ε₀ = 0.005: the last annulus; ε₀ = 0.11: VA011).

*VA's path and its square* (D61 A.3, VC). VA uses part VC's family: the real segment z ∈ [−1, t], t = ε/4; J: X = t, 0 ≤ Y ≤
min(κt, θ); F: Y = κX (κ = Im A/|Re A|, A = η²/4 − 1, the quadratic floor); C: Y = θ, up to X = δ₀ = 0.05; then, per η-box, a
connector, a float descent polyline to |1 − z| < 0.02 and V's end piece. On J, F, C, in the weight-0 coordinates (u, q, ρ) =
(ε/X, θ²/X, X) (and (u, p, ϑ) = (ε/X, X/θ², θ²) for q large), E/X² = fA(u, q, ρ) exactly (Lemma VA1), and fA(u, q, 0) is the
weight-4 limit in closed form (Lemma VA3). In (u, q) the four blow-up limits are the faces of one square: the rays q/u = ψ² give
Ψ̃_ψ (θ = √εψ), the face q = 0 gives Ψ_φ (θ = εφ) and Ψ_∞, the corner u = q = 0 the matching region (E ≈ −εθ² − θ⁴/2 − X³/3), and
q → ∞ (p → 0) the fixed-θ limit; each limit has a proved negative path (Lemma VC). The limit vanishes identically on the face q = 0
and at the junction's λ = 0, where VA's mean value form needs sup ∂_ρfA < 0. A hand proof along VA's path therefore needed (a) a
two-term expansion of fA with an explicit remainder, uniform on the closed square including the faces where the first term
vanishes, and (b) a hand treatment of the outer piece X ≥ δ₀, which VA certifies with 5193 segments over 76 η-boxes (VA011: 16083
over 220), the largest bound only −1.5·10⁻⁷. That is what the handoff asked for. Neither is needed below.

## 1. The proof

**Lemma 1 (F on and beyond the unit circle) [proved].** For |ν| ≥ 1,

  F_η(ν) = G_η(ν) := ℓ(η − ν) + M(Re ν) − ½,   ℓ(τ) = log|τ| + ½Re τ².

*Proof.* For |ν| > 1 substitute Λ⁺: log|ν| + ½Re ν² + Λ⁺(ν) = M(Re ν) − ½. On |ν| = 1, Λ = 0 and Λ⁺(ν) = (Re ν)² − ½ − ½(2(Re ν)² −
1) = 0, so the same identity holds. ∎ (Inside the disc F_η = G_η − H(|ν|) ≤ G_η; this is not used. D61 writes the same function as
log|1 − z²| + Re(cz²) + C; the identity C = Λ⁺ = M(a) − ½ − log|σ| − ½Re σ² is D61 A's own.) So beyond the circle F_η depends on ν
only through τ = η − ν and Re ν, and M is C¹, even and nondecreasing in |a|, with 0 ≤ M′(a) ≤ 2 for a ≥ 0.

**Lemma 2 (vertical lines) [proved].** Fix p ∈ ℝ, write ν = p + ir and τ = η − ν = a + ib (a = Re η − p fixed, b = Im η − r). Then
G_η(ν) = ½log(a² + b²) + ½(a² − b²) + M(p) − ½ and

  ∂G_η/∂b = b(1 − a² − b²)/(a² + b²).

So on a vertical segment contained in {|η − ν| < 1}, traversed towards Im ν = Im η (|b| decreasing), G_η strictly decreases.
*Proof.* Differentiate; for 0 < a² + b² < 1 the derivative has the sign of b, and |b| decreasing keeps a² + b² < 1. ∎

**Lemma 3 (the horizontal line through η) [proved].** On ν = p + i Im η with p ≥ 0 and a = Re η − p ∈ (0, 1),

  dG_η/dp = −(a + 1/a) + M′(p) = −(a − 1)²/a − (2 − M′(p)) < 0     (M′(p) = 2 min(p, 1)).

*Proof.* τ = a is real, so G_η = log a + ½a² + M(p) − ½; differentiate; (a − 1)²/a > 0 for 0 < a < 1, and M′(p) ≤ 2. ∎

**Theorem A (the valley condition on 1 < |η| < 2) [proved].** Let 1 < R = |η| < 2. Along Π_η, F_η(ν) < F₀ for every ν ≠ η/2; the
maximum over V ∪ H is m(η) < 0 above, attained only at η/|η|; Π_η contains the chord {(η/2)(1 + z) : |z| ≤ x₀} with x₀ = min(½,
(2/R − 1)/2), which lies in the open unit disc; and Π_η satisfies D61's form (2) for every 0 < t < 2/R − 1. Hence every such η
satisfies the valley condition.

*Proof.* By the symmetries (which map Π_η to Π_η̄ and Π_{−η}: the formula for ν₁ is equivariant) let 0 ≤ θ ≤ π/2.

(Ch) For ν = sη, 0 < s ≤ 1/R, |ν| ≤ 1, so Λ = 0 and F_η(sη) = g_η(s), with g_η′(s) = (1 − 2s)(1/(s(1 − s)) − Re η²) and 1/(s(1 −
s)) ≥ 4 > R² ≥ Re η²: g_η strictly increases on (0, ½] and strictly decreases on [½, 1); F_η → −∞ at s → 0. So F_η < F₀ on Ch
except at s = ½, and at s = 1/R, with d = 1/R − ½ = ε/(2R), 4s(1 − s) = 1 − 4d² and (1 − s)² + s² − ½ = 2d²,
F_η(η/|η|) − F₀ = log(1 − ε²/R²) + Re(η²)d² = m(η); log(1 − x) < −x gives m(η) < −ε²/R² + (ε²/4)cos 2θ = −ε²(4 − R²cos 2θ)/(4R²) < 0.

(V) ν = cos θ + ir, sin θ ≤ r ≤ R sin θ. Here |ν|² = cos²θ + r² ≥ 1, so F_η = G_η (Lemma 1), and η − ν = a + ib with a = (R −
1)cos θ, b = R sin θ − r ∈ [0, (R − 1)sin θ], so a² + b² ≤ (R − 1)² < 1. By Lemma 2, F_η(ν) ≤ F_η(η/|η|) = F₀ + m(η), strictly for
r > sin θ. (θ = 0: V is the point 1.)

(H) ν = p + iR sin θ, cos θ ≤ p ≤ R cos θ. Here |ν|² ≥ cos²θ + sin²θ = 1, so F_η = G_η; a = R cos θ − p ∈ [0, (R − 1)cos θ] ⊂ [0,
1) and p ≥ 0. By Lemma 3, F_η(ν) < F_η(ν₁) ≤ F₀ + m(η) for cos θ < p < R cos θ, and F_η(η) = −∞. (θ = π/2: H is the point η.)

The chord {|z| ≤ x₀} is {sη : |s − ½| ≤ x₀/2} ⊂ Ch, and R(1 + x₀)/2 < 1 since x₀ < 2/R − 1. For D61's form (2), Π_η from b = (1 +
t)η/2 runs along Ch (F_η = g_η decreasing there), then V and H. ∎

*Explicit margins on Π_η.* Paper 18's Lemma 4.6 (ii) needs, for η near η₀, a uniform gap off the sub-chord {|z| < x₁} along the
path it uses, which is the scaled path (η/η₀)Π_{η₀}, not Π_η. That gap needs nothing beyond the lemma as written: its own
compactness argument applies unchanged (F_η is jointly continuous in (η, ν) off ν ∈ {0, η}, and Λ is continuous across the circle,
Λ⁺ = H(1) = 0 there). The paragraph below is a bound on Π_η itself, an aside not needed by the lemma: by the monotonicity of g_η on
Ch and the bound m(η) on V ∪ H, F_η − F₀ ≤ max(log(1 − x₁²) + x₁²Re(η²)/4, m(η)) on Π_η off the sub-chord — continuous in η and
negative on 1 < |η| < 2, where m is defined (Π_η meets the circle at η/|η| only). (Referee round 148, D3eq-2.)

**Corollary B [proved, given Paper 18's Lemmas 4.4–4.6 as they stand].** Ω′ = {|η| < 2}: for |η| ≤ 1 the straight path (the
domination condition is void), for 1 < |η| < 2 Theorem A. So ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2} uniformly on compact subsets of
the open disc |η| < 2 by hand — Paper 18's Theorem 2 and D‴.61 (iii)(b) — and everything downstream of Ω′ (D‴.61 (iv) for every
r < 1, Paper 18's corollary for frame rules (cor:frame) for every frame, the half-squeeze rule at every σ*/σ) loses its computer-assisted step. The
certificates D, V and VA become independent cross-checks. Nothing here touches the rest of Theorem 2's proof (Lemmas 4.3–4.6),
which was not re-verified.

## 2. Why the corner analysis is not needed

VA's path leaves the real axis at b, a distance ε/4 from the saddle, and crosses the circle |σ| = 1 at z ~ ε, where C (≈ δ₊² −
(Re σ − 1)₊²) is of order ε² and cancels to order ε³; that cancellation, at four different scalings of θ against ε, is the
corner's four regimes. Π_η instead stays on the straight chord — where F_η is exactly g_η and needs no expansion — up to the
circle itself, z = 2/R − 1 ≈ ε/2, and beyond the circle follows two lines on which F_η is monotone. There is no saddle-point
structure left beyond the circle once ν moves only in Im ν (at fixed Re ν, M is constant and ℓ is increasing in |Im τ| on |τ| < 1)
and then only in Re ν (along τ real, ℓ′ = a + 1/a ≥ 2 ≥ M′). The cost is only that the margin m(η) is of order ε²sin²θ + ε³ (it
must vanish as |η| → 2: the saddle reaches the circle).

*Consistency with D61's blow-up limits [numerical; the closed forms by hand].* In VC's charts Π_η becomes:

- Ψ_φ chart (θ = εφ, z = ε(x + iy)): y = 0 (x ≤ ½), then x = ½ (0 ≤ y ≤ φ), then y = φ (x ≥ ½). Along it Ψ_φ = −x², then −¼ −
  y(2φ − y), then −x² − φ² − (x − ½)³/3 (on y = φ the terms −4φ²x and (x − ½)·4φ² cancel to −2φ²): negative for every φ, whereas
  the straight segment holds only for φ < φ_dom = 1.4397.
- Ψ̃_ψ chart (θ = √εψ, z = εx + i√εy): y = 0 (x ≤ ½), then the segment x = ½ + ψy (0 ≤ y ≤ ψ), then y = ψ (x ≥ ½ + ψ²). On the
  segment P = −ψ²/2 < 0 and Ψ̃_ψ = −ψ²[t(2 − t) + ψ²t²(1 − t/2)²], t = y/ψ; on y = ψ it is Lemma VC's cap, from −ψ² − ψ⁴/4 to −ψ² −
  ψ⁴/2. (On y = 0 the limit is 0 for x ≤ ½: that face is degenerate, but there F_η = g_η exactly.)

`check_algebra.py` compares E/ε³ and E/ε² along Π_η at ε = 10⁻⁴⁰ with these closed forms (agreement to 10⁻³³). In z, V has slope
cot θ (direction ie^{−iθ}), steeper than VA's floor κ = ω²sin 2θ/((1 − ω²) + 2ω²sin²θ) < cot θ; the two agree in the fixed-θ limit,
where VA's floor is vertical in σ and its cap Im z = θ is, to first order, the line Im σ = Im η. So Π_η is VA's route in the
fixed-θ regime, made exact, and differs near the corner by leaving the chord at the circle instead of at b. The (u, q) square does
not arise.

## 3. What remains

- **For the valley condition: nothing.** The last annulus 1.995 ≤ |η| < 2, and the whole annulus 1 < |η| < 2, are proved by hand
  (Theorem A).
- **Not done, and no longer needed:** the hand proof along VA's own path (a remainder bound for fA, fB uniform on the (u, q)
  square, and the outer piece X ≥ δ₀). VA stays a valid computer-assisted proof of the same statement along a different path.
- **To be checked by a referee** (all four checked in round 148, D3eq: confirmed). (1) That Paper 18's F_η beyond the circle is Lemma 1's G_η — it is the definition of Λ⁺ (Lemma
  4.4), and D61's C is the same function (checked symbolically and in balls, section 4). (2) Lemmas 2 and 3 (one derivative each)
  and the geometry of V and H (outside the closed disc, inside |η − ν| < 1). (3) That Lemma 4.6 (ii) asks nothing of the path
  beyond Definition 4.5 — in particular Lemma 4.4's bound holds on every compact subset of ℂ, so a path leaving the disc far from
  η/2 is allowed (it is: the integrand of the path integral in Lemma 4.6's proof, (η − ν)^Kν^{K−1}e^{K((η−ν)²+ν²)/2}c_K(√Kν),
  is entire, and D61's form (2) is phrased for an arbitrary path from b to η).
  (4) Whether any statement downstream used VA's specific path rather than Ω′ (I found none: D61 (iii)(b), (iv) and Paper 18's
  Theorem 2 and its frame corollary use only Ω′).
- The proof is short enough to be suspicious of. I looked for the catch and did not find one: the identity F = G is exact beyond
  the circle, the monotonicities are exact, and the numerical checks below find no positive value on 282 064 path points, including
  |η| = 2 − 10⁻¹², while the same evaluator finds the straight path positive in the sector where the domination condition fails.

## 4. Numerical checks (`review/research_annulus/`)

None of these is part of the proof.

**`check_cvh.py`** (python-flint 0.9.0, Arb at 320 bits; output `check_cvh.out`, 15 s):

1. Lemma 1: |F_η − G_η| ≤ 3.2·10⁻⁹⁴ at 2000 random points beyond the circle (F_η from Paper 18's definition with Λ), and
   |G_η − F_η − H(|ν|)| ≤ 3.6·10⁻⁹⁵ inside.
2. E = F_η − F₀ < 0, as a rigorous ball inequality, at all 282 064 sample points of Π_η (chord points clustered geometrically at
   s = ½, V and H points clustered at their starts, the last H point at relative distance 10⁻³⁰ from η): |η| = 2 − ε for ε = 0.005, 0.002, 10⁻³, 10⁻⁴,
   10⁻⁶, 10⁻⁹, 10⁻¹², and |η| = 1.0001, 1.01, 1.1, 1.3, 1.5, 1.7, 1.8, 1.8927, 1.9, 1.95, 1.98, 1.99; 33–35 angles in [0, π/2]
   (including ε/10, ε, √ε, π/2 − 10⁻⁹) and 7 in the other quadrants (the same formula for ν₁). E computed twice — from Paper 18's
   F_η and from D61's centred form log|1 − z²| + Re(cz²) + C — and, on every fourth angle, a third time with the repository's
   `gh_lastann.E_centred_safe` (read-only import); all agree. 0 failures. The largest values: on V ∪ H exactly m(η) (e.g.
   −3.139·10⁻⁸ at ε = 0.005, θ = 0; −2.501·10⁻³⁷ at ε = 10⁻¹²), max(E − m) on V ∪ H ≤ 2·10⁻⁹² (rounding), as Theorem A says.
3. Control: on the straight path [0, η] the same evaluator finds E > 0 (rigorously, in balls) at |η| = 1.95, arg η = 0.3 and
   0.097; |η| = 1.98, arg 0.2; |η| = 1.995, arg 0.1 — the sector where the domination condition fails (D61: at |η| = 1.95 at least
   0.087 … 0.615).
4. A computer proof for 197 values of η that does not use Lemmas 2 and 3: V and H covered by sub-intervals, on each the smaller of
   the interval enclosure of G_η and its centred form G_η(mid) + rad·|G_η′(interval)| below F₀ (200 bits), at ε = 0.005, 0.002,
   10⁻³, 10⁻⁴, 10⁻⁶, 10⁻⁹ and 32–33 angles each including θ = 0 and π/2: none failed (up to 114 530 pieces for H at ε = 10⁻⁹, θ = 0).
5. Lemmas 2 and 3's derivative formulas against central differences of F_η at 600 random points: agreement to 4·10⁻⁵⁰, signs as
   claimed.

**`check_algebra.py`** (sympy 1.13.3, mpmath 1.3.0; output `check_algebra.out`): g_η(1/R) − g_η(½) = m(η) symbolically; the
derivatives of Lemmas 2 and 3 symbolically; the blow-up limits along Π_η of section 2 at ε = 10⁻⁴⁰ (φ = 0.1 … 25, ψ = 0.02 … 20).

## 5. How it would enter (suggestions only; nothing was edited)

(Incomplete, D3eq-4: the full list of downstream statements is the "Consequences" table of review/REVIEW_round148_D3eq.md, applied
after round 148.)


- **Paper 18.** §4, proof of Theorem 2: the "By hand" paragraph covers the whole disc with Theorem A (three lemmas, half a page);
  the "Computer-assisted" paragraph and Appendix A (D, VA, V) become an optional cross-check or go; the abstract's "with a
  computer-assisted proof in ball arithmetic elsewhere", the introduction's "One ingredient is computer-assisted …" and the
  keywords change accordingly; the frame corollary's "computer-assisted for t ≥ 0.0966" becomes by hand.
- **D‴.61.** (iii)(b): Ω′ is the open disc by Theorem A; D, V, VA, VC become history and cross-checks; Status and Appendix A
  likewise. Lemma VC and the four regimes remain correct statements about VA's path.
- **HANDOFF.md** "Open mathematics": the hand proof of the last annulus is done, and of 1 < |η| < 2; §8 item 1's "computer-assisted"
  becomes by hand.
- **RESEARCH_note_c2.md**, Proposition A (iii): "With certificates D and VA, Ω′ is the whole disc" becomes by hand (that proof uses
  Ω′ only; its owner should confirm).

## Files

- `review/RESEARCH_note_annulus.md` — this note.
- `review/research_annulus/check_cvh.py`, `check_cvh.out` — section 4, items 1–5. Run from the repository:
  `OMP_NUM_THREADS=1 venv/bin/python review/research_annulus/check_cvh.py`.
- `review/research_annulus/check_algebra.py`, `check_algebra.out` — section 4, last paragraph. Needs sympy and mpmath, which the
  venv lacks: `PYTHONPATH=<dir with sympy 1.13.3 and mpmath 1.3.0> venv/bin/python review/research_annulus/check_algebra.py`.

## Review history

- **Round 148** (2026-11-05; review/REVIEW_round148_D3eq.md, D3eq): the proof stands — no catch found; Theorem A, Lemmas 1–3 and the
  geometry of Π_η re-derived by hand and checked independently (mpmath at 60 digits, 506 + 444 values of η, F_η from Paper 18's
  definition), the note's scripts rerun with identical outputs. Changed here: the margins paragraph says that it bounds Π_η itself
  and that Lemma 4.6 (ii)'s compactness argument covers the scaled path the lemma uses, and "negative" is restricted to
  1 < |η| < 2 (D3eq-2); section 5 points to the report's full list (D3eq-4); section 3's four points marked checked; the scratch path
  in the header replaced by `review/research_annulus/`. No statement, proof step or number changed. D3eq-1 and -3 (Paper 18's
  sentence after Definition 4.5; the axes as θ = 0, π/2 of Theorem A, where m(η) equals Paper 18's two axis gaps) and D3eq-5 (the
  scripts' printed "computer-assisted") were applied downstream: Paper 18, D‴.61, D‴.58, D‴.62, the c₂ note, PAPERS.md, HANDOFF.md,
  `exp_gh_peano.py`, `exp_width_mismatch.py` and the certificates' verifier.
