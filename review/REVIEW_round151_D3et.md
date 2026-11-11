# Referee report D3et, round 151: the Paper 15 manuscript, third review

## Scope
This is the third review of `manuscript/paper15/paper15.tex` and `supplement.tex` at commit c5de42f, with the committed PDFs, `README.md` ("Revision after round 145"), the Round 145 entry of `REVIEW_LOG.md` and `review/REVIEW_round145_D3en.md`. The main object is the new **Lemma 2.6** (l. 188–190), proved in Supplement S6 (S l. 263–284), and its use in Theorem 3.4(b).

- **Workspace.** The commit was exported to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref151/repo`, with the venv symlinked. The repository was not modified.
- **Line numbers.** They refer to `paper15.tex`; "S" marks `supplement.tex`.
- **Evidence tags.** [read] = in the files; [md] = in the result reports; [own] = my own computation; [paper] = an arXiv text I downloaded and read at the cited passage; [memory] = background knowledge I did not check.
- **Out of scope.** Bibliographic placeholders and `% CHECK` items.

## Lemma 2.6 (verdict, step by step)

**Verdict: correct.** The proof is complete for what the lemma states, and the hypotheses are verified correctly for both flows. Every step reproduces numerically [own].

What remains is presentation, plus one place where a downstream proof uses more regularity than the lemma gives (D3et-1).

**(i) The linearisation at the limit.**
- "Semisimple 0" is read from "ker DG(p) = T_pΓ′ and all *other* eigenvalues, counted with algebraic multiplicity, in Re < 0". That reading is right, but the statement should say it outright (D3et-3).
- P_p is the Riesz projection of the 0-group. It is continuous on 𝒦, because at every p ∈ 𝒦 exactly d eigenvalues are 0, so the 0-group cannot split.
- The bounds |e^{sA_p}| ≤ C and |R_p(s)| ≤ Ce^{−νs} hold uniformly on the compact 𝒦. This needs ν below the normal gap; "we decrease ν where needed" covers it, and also covers Jordan blocks among the nonzero eigenvalues.

**(ii) The first variation.**
- I re-derived the variation-of-constants formula.
- Gronwall gives a bound on Y that is uniform in x.
- The identity Y − Y_∞ = R_p(t) + ∫₀ᵗ R_p(t−u)ΔY du − P_p∫_t^∞ ΔY du is exact. It gives O((1+t)e^{−νt}).
- The bound |Δ| ≤ C|x(t) − p| needs only DG to be Lipschitz near 𝒦 and bounded elsewhere: trivial.

**(iii) The second variation.**
- The equation Ẇ = DG W + D²G[Y, Y] is right.
- The bound |b − b_∞| ≤ C(1+t)e^{−νt} uses only that D²G is Lipschitz.
- The key claim holds. For a C² surface σ in Γ′, differentiating DG(σ)∂₁σ ≡ 0 in s₂ gives D²G(p)[a, a′] = −DG(p)∂₁∂₂σ ∈ range DG(p). This needs only G ∈ C² and Γ′ ∈ C².
- Since Y_∞ maps into T_pΓ′, P_p b_∞ = 0, and the constant forcing lies in the decaying subspace.
- The splitting into ∫e^{(t−u)A_p}[ΔW + b − b_∞] + ∫₀ᵗ R_p(s)b_∞ ds is exact. The limit W_∞ is as stated.
- The three error terms are O(t²e^{−νt}), O((1+t)e^{−νt}) and O(e^{−νt}), so the rate O((1+t)²e^{−νt}) is right. On the paper's own augmented flow the t² factor actually appears (Verification).

**(iv) The C² conclusion.**
- Uniform convergence of Dφ_t and D²φ_t on 𝒩 gives Φ ∈ C², by the standard theorem. All constants are uniform because they depend only on 𝒦, C₀, ν and bounds on the compact set containing the orbits.
- The statement "φ_t are C³" is true but not needed.
- The proof uses only that D²G is locally Lipschitz (G ∈ C^{2,1}). That is the hypothesis SSP attribute to [Kat91, Cor. 5.1] [paper, SSP p. 13: "L three times differentiable with locally Lipschitz third derivative … Φ ∈ C²(U) [Kat91, Corollary 5.1]", via "[Fal83, Theorem 5.1]"].

**Hypotheses, (a) the flow of F.**
- F ∈ C³ (f̂ ∈ C⁴, M smooth) and Γ ∈ C³.
- A is diagonalisable (similar to −M^{1/2}HM^{1/2}), so 0 is semisimple.
- The Polyak–Łojasiewicz-type bounds hold: c d² ≤ J ≤ Cd² and c d ≤ |∇J|_M ≤ Cd, from Taylor's formula at the nearest point of Γ with H positive definite on the normal space. They give J̇ ≤ −κJ, a trajectory length O(d(z₀)) and the rate e^{−κt/2}.
- The limits lie in a compact 𝒦 ⊂ Γ₀ for 𝒩 thin enough. All correct.

**Hypotheses, (b) the augmented flow.**
- The field (F + ⟨ζU⟩, (Λη_x)_x) is C³, with U = mM∇f̂ ∈ C³.
- The linearisation at (p, 0) is block upper triangular, with diagonal blocks A and I⊗Λ. Its kernel is T_pΓ × {0}. The algebraic multiplicity of 0 equals A's, because Λ is invertible, so 0 is semisimple.
- J̇ = −|∇J|²_M + ∇J·⟨ζU⟩ ≤ −½|∇J|²_M + C|η|² by Young's inequality in the M/M⁻¹ pairing. Hence J̇ ≤ −(κ/2)J + C|η₀|²e^{−2ν_Λt}, and then J ≤ C(J(0) + |η₀|²)e^{−2νt} with 2ν < min(κ/2, 2ν_Λ). Then |ż| ≤ C(d(z₀) + |η₀|)e^{−νt}.
- The bootstrap ("as in (a)") that keeps the orbit where the Polyak–Łojasiewicz bounds hold is routine. All correct.

**Comparison with the sources.**
- **LWA [paper].** Lemma B.1 is one line: "Applying Theorem 5.1 of Falconer (1983) with f(·) = φ(·,1) suffices". The surrounding hypotheses are a gradient flow, L ∈ C³, Γ ∈ C², rank ∇²L maximal (Assumptions 3.1–3.2).
- **SSP [paper].** They attribute C² to [Kat91, Cor. 5.1] for L ∈ C³ with locally Lipschitz third derivative, and the attracting neighbourhood to [Kat91, Prop. 3.5].
- **Parsons–Rogers [paper].** §II (procedure step 1: λ₁ = … = λ_m = 0, Re λ_{m+1..n} < 0) and App. B ("f twice continuously differentiable … globally attractive manifold of equilibria … satisfies the conditions of Katzenberger") support Remark 2.3's description.
- **What the time-one map adds.** Falconer applied to φ(·, 1) needs nothing gradient-specific. The paper's direct proof is the continuous-time version, and it needs only the same C^{2,1} data.
- **Classical background [memory].** The result is the smoothness of the asymptotic phase (the stable foliation) of a normally attracting manifold of equilibria (Fenichel 1974/1977; Aulbach, LNM 1058, 1984). A sentence saying the result is classical would help (D3et-5).
- **Not read.** Katzenberger 1991 and Falconer 1983 (projecteuclid.org connection refused by the proxy). A web search for "Kat91 Cor. 5.1" returned nothing relevant.

**Use in Theorem 3.4(b).** Lemma 2.6 delivers exactly the C² that Katzenberger's drift ½∂²Φ̃[Σ̃] needs. Equation (2.1) and Lemma 3.2 are obtained by differentiating identities with a vanishing factor at points of Γ. They hold with Φ only C², by a one-sided limit argument (or by differentiating Φ∘φ_t = Φ twice and then in t at t = 0).

One gap remains: the expansion at l. 288 claims a remainder O(|η₀|³) (D3et-1).

## Status of the round-145 issues

| D3en id | resolved? | evidence |
|---|---|---|
| 1 (SSP paragraph) | Yes | l. 101 checked against arXiv:2404.12293v1 [paper]: Thm 3.15 (for −∇L), Thm 4.1 (½Δ_ηL̂(w,0), time scale ασ²), Thm 4.7, Cor 5.10, Cor 5.11 ("effectively the same result as [LWA21, Cor. 5.2]"), §5.1.3 (L(w+η), i.i.d.), §6.2 (correlation across the coordinates of one draw), Rem 4.9 (Itô). All accurate. One normalisation nit (D3et-6). |
| 2 (Prop 5.7 hypothesis) | Yes | l. 478, 492: every N columns independent; the dimension formula follows; (a) re-derived. |
| 3 (τ_c = 100; (c) asymptotic) | Yes | l. 518. Every K2 number matches `colour_selection.md` l. 257–260 [md]. 9.8 / 13.9 / 15.5 are differences of the md medians +5.18, +15, +28.9, +44.4. Remainder in closed form (l. 448, 486). |
| 4 (O(h) at fixed intensity) | Yes | l. 377–381: (1−e^{−x})/(1+e^{−x}) = tanh(x/2), checked. |
| 5 (Thm 3.4(b) hypotheses) | Yes, with D3et-1 to D3et-3 | Lemma 2.6 and S6 are correct; Table 1 and Remark 2.3 updated. |
| 6 (credit for no-diffusion) | Yes | l. 245. |
| 7 (temporal-spectrum overclaims) | Yes | l. 62, 103, 272; abstract. |
| 8 (abstract: N-atom fits) | Yes | "some fits … stay local minima". |
| 9 (Lemma 3.3 on the spectrum) | Yes | l. 258–263. |
| 10 (isotropic noise) | Yes | l. 103 matches LWA Cor 5.1 [paper]: −½∂Φ∇ln\|∇²L\|₊. |
| 11 (Table 1 rows) | Yes | l. 80–81. |
| 12 (supplement status table) | Yes | S l. 43–58 (Table S1). |
| 13 (abstract length) | Yes | 246 words by my count (formulas as one word). |
| 14 (Heun at hλ ≈ 1) | Yes | l. 402; S l. 158. |
| 15 (notation overloads) | Yes for the listed symbols | No leftover Π, 𝖻, 𝖼 or ν_j/ν_k. Lemma 2.6 introduces new clashes (D3et-4). |
| 16 (Remark S1.4) | Yes | S l. 103. |
| 17 (SSP's non-degenerate case) | Yes | l. 101. |

## New findings

1. **D3et-1. The O(|η₀|³) remainder in Theorem 3.4(b)'s proof is not supplied by Lemma 2.6** (l. 286–288, p. 12).
   - ψ = ∂Φ·U is only C¹ when Φ is C². Then ψ(z(t)) − ∇ψ(z₀)·(z(t) − z₀) is o(|η₀|), not O(|η₀|²).
   - So Φ̃(z₀, η₀) − z₀ = ∫⟨ζ∇ψ(z₀)·δz⟩ + o(|η₀|²).
   - An O(|η₀|³) remainder needs Φ ∈ C^{2,1}. With F only C³ (f̂ ∈ C⁴), the lemma's method gives at most C² for Φ: a third variation would need D³G Lipschitz.
   - The conclusion is unaffected. Φ̃ ∈ C² makes its second-order Taylor polynomial unique, so o(|η₀|²) suffices both for "no linear term" and for identifying ½∂²_ηηΦ̃[η₀, η₀].
   - Fix: write o(|η₀|²) (and "+ o(|η₀|²)" in l. 286's sentence), or assume f̂ ∈ C⁵ and say why.
2. **D3et-2. Remark 2.3, l. 147: "Under this spectral condition the limit map is C² (Lemma 2.6)" omits a hypothesis.**
   - Lemma 2.6 also assumes a neighbourhood that is attracted at a uniform exponential rate, with orbits in a compact set.
   - The pointwise spectral condition implies that only through normal-hyperbolicity theory, which the paper neither proves nor cites. The paper verifies it only for its two flows, through J (S6).
3. **D3et-3. Lemma 2.6's statement (l. 189).**
   - Say explicitly "0 is a semisimple eigenvalue of DG(p) with eigenspace T_pΓ′", since the proof's first line relies on it.
   - The last sentence refers to equation (3.1), which carries the noise and is defined later. Say "the noiseless augmented field of (3.1)".
   - Optionally, note that G ∈ C² with locally Lipschitz D²G suffices. That is all the proof uses, and it matches SSP's description of [Kat91, Cor. 5.1].
4. **D3et-4. Notation clashes introduced by §2.4 and S6.**
   - G (the vector field) against G(a), G̃ in Table 2 and Theorem 4.2.
   - x (a point) and x(t) against the node index x, used in the same proof of Theorem 3.4(b) (η_x, ζ_x).
   - The lemma's ν and the section's 𝒦 are not in Table 2.
   - Fix: rename the field (for example V) and the point (for example y), or list them in Table 2.
5. **D3et-5. Attribution sentence (l. 192).**
   - "This is what Falconer's theorem … gives" is asserted about a text the authors say they have not read.
   - Reword to what was read: LWA obtain C² for gradient flows from [Fal83, Thm 5.1] applied to φ(·, 1); SSP attribute it to [Kat91, Cor. 5.1] under L ∈ C³ with locally Lipschitz third derivative.
   - Optionally, add that the statement is the classical smoothness of the asymptotic phase of a normally attracting manifold of equilibria (Fenichel; Aulbach) [memory; to be checked before citing].
6. **D3et-6. SSP normalisation (l. 101).**
   - "The deterministic flow of Δ_wL/(2N)" is in SSP's normalisation, L = N⁻¹Σ(f − y)² (their (86)) [paper].
   - The sentence sits in a paper where L = (2N)⁻¹Σr² (l. 123). In that notation it reads Δ_wL/N.
   - Add "(in their normalisation)". This is cosmetic: same flow up to a constant.
7. **Nothing else found wrong** in the fresh read of:
   - Theorem 3.4's statement, including the covariance and normalisation of (b) and the cross-check of Σ_η;
   - Table 1's statuses;
   - the SSP paragraph;
   - Prop. 5.7 with its corrected hypothesis (re-derived);
   - Prop. 5.4's exact remainder;
   - §6.3's numbers.

## Verification I did

- **Script.** `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref151/scripts/lemma26_check.py`, with second-order forward-mode jets in `jet.py` giving exact DG and D²G.
  - Integration: DOP853, rtol 10⁻¹², atol 10⁻¹⁴.
  - D²φ_t comes from the joint first and second variational equations. D²Φ also comes from central finite differences of the integrated limit map.
  - Output: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref151/lemma26_out.txt`.
- **Results [own]:**

| example | normal spectrum | max\|D²φ_T − FD\| (ε = 10⁻³) | \|P_pD²G[a,a′]\|, tangent / random | \|W(t) − W_T\|·e^{μt}/(1+t)² |
|---|---|---|---|---|
| E1: non-gradient field on R³, curve of equilibria {x₂ = ½sin x₁, x₃ = 0.3x₁²}, two starts | −1.06 ± 2.13i; −1.22 ± 2.02i | 2.6·10⁻⁸, 3.6·10⁻⁸ (\|D²Φ\| ≈ 0.45, 0.53) | 3.7·10⁻¹⁶, 2.4·10⁻¹⁶ / 7.5·10⁻², 7.1·10⁻² | bounded, 0.01–0.04, for t = 4–28 (then the 10⁻¹³ floor) |
| E2: Jordan normal block, exact Φ = (x₁, ½sin x₁, 0.3x₁²) | −1, −1 (defective) | 1.7·10⁻⁸; against the exact D²Φ: 1.8·10⁻¹⁴ | 0 / 0 (range DG holds all values of G here) | 0.03–0.11 (t·e^{−t} behaviour) |
| E3: the paper's noiseless augmented flow (3.1), on D3en's model f̂₁ = z₁ + 0.3z₂z₃ + 0.2z₁², f̂₂ = z₂ + 0.4 sin z₁ z₃ + 0.1z₃², ϖ = ½, OU τ_c = 2 (R⁵, Γ×{0} a curve) | −0.430, −0.5 (×2), −0.679 | 4.0·10⁻⁸ (\|D²Φ̃\| ≈ 0.33) | 8.1·10⁻¹⁷ / 8.3·10⁻² | ≈ 0.008–0.011 for t = 12–48: the (1+t)² factor of S6(iii) is attained |

  - In all cases |Y_T − P_pY_T| ≤ 10⁻¹⁴: DΦ maps into T_p.
  - Continuity (C²): along a segment of 11 points in E1, D²Φ changes by 0.026–0.029 per step of 0.1, smoothly, with no jumps.
- **Rate.** D²φ_t converges at e^{−μt} times a polynomial of degree at most 2, μ being the normal gap. The S6 rate is confirmed and essentially sharp. The "decrease ν" step costs nothing in practice.
- **Literature [paper].** The arXiv PDFs were downloaded to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref151/dl/` (own directories, pdftotext only, nothing executed):
  - SSP 2404.12293v1: §3.3 (Def 3.7, the Kat91/Fal83 attributions), Thm 3.15, Thm 4.1, Thm 4.7, Rem 4.9, §5.1.3, Cor 5.10, Cor 5.11, §6.2;
  - LWA 2110.06914: Assumptions 3.1, 3.2, Lemma B.1, Thm B.7, Cor 5.1;
  - Parsons–Rogers 1510.07031v2: §II, App. B.
- **Build [own].** In a scratch copy, three passes each: 0 errors, 0 undefined references, 0 overfull boxes; 28 + 13 pp. The text of both PDFs is identical to the committed ones (pdftotext diff). The labels resolve as the supplement hard-codes them (Lemma 2.6, (3.1), Remark 2.4, Prop 5.5, Remark 4.5, Assumption 2.5).
- **Number check [own].** `numcheck15.py`: 407 tokens; 3 not in the md (5.4, 15.5, 16.1), each derived visibly in §6.3. I checked the §6.3 K2 figures against the md table by eye.
- **Not checked.** Katzenberger 1991 and Falconer 1983 (not reachable). The journal version of SSP.

## Overall assessment
**Minor revision (close to accept).**

- Lemma 2.6 is correct.
- Its proof is a clean continuous-time version of the Falconer and Katzenberger regularity results.
- Its hypotheses are verified correctly for both flows, and it discharges the last "assumed" in Theorem 3.4(b).
- All 17 round-145 items are resolved, and the SSP paragraph is now accurate against the full text.
- What remains is local:
  - one remainder order in Theorem 3.4(b)'s proof that over-claims the regularity delivered (harmless to the result);
  - two hypotheses to make explicit (in Remark 2.3 and in the lemma's statement);
  - notation clashes from the new section;
  - an attribution wording;
  - a normalisation nit.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3et-1 | minor | Thm 3.4(b) proof: "+O(\|η₀\|³)" needs ψ ∈ C^{1,1} (Φ ∈ C^{2,1}). Lemma 2.6 gives Φ ∈ C² only (F ∈ C³), so the remainder is o(\|η₀\|²). That suffices, since Φ̃ ∈ C² has a unique second-order Taylor polynomial. | l. 286–288: write o(\|η₀\|²), adding one line on why it suffices; or assume f̂ ∈ C⁵ and extend the lemma. |
| D3et-2 | minor | Remark 2.3 says the spectral condition alone gives C² "by Lemma 2.6". The lemma also needs uniform exponential attraction of a neighbourhood, with orbits in a compact set. | l. 147: "under this spectral condition and exponential attraction of a neighbourhood (verified for our two flows, Supplement S6)". |
| D3et-3 | minor | The lemma's statement leaves the semisimplicity of 0 implicit (the proof uses it). It cites (3.1), which includes the noise, by forward reference. It assumes C³ where C^{2,1} is used. | l. 189: "0 is a semisimple eigenvalue with eigenspace T_pΓ′"; "the noiseless field of (3.1)"; optionally "G ∈ C² with D²G locally Lipschitz"; S l. 281: drop "φ_t are C³". |
| D3et-4 | minor | Notation: G (field) vs G(a), G̃ (Table 2); x (point) vs node x (same proof as η_x); ν and 𝒦 missing from Table 2. | §2.4, S6: rename (e.g. V, y); Table 2. |
| D3et-5 | minor (wording) | "This is what Falconer's theorem … gives" is said of an unread text. The result is classical (normal hyperbolicity) [memory]. | l. 192: attribute via LWA Lemma B.1 and SSP p. 13 (Kat91 Cor 5.1); optionally cite Fenichel/Aulbach after checking them. |
| D3et-6 | trivial | "Δ_wL/(2N)" uses SSP's L = N⁻¹Σ(f−y)² (their (86)); the paper's L = (2N)⁻¹Σr² (l. 123). | l. 101: "(in their normalisation)". |

## Sources
- [read] `manuscript/paper15/paper15.tex`, `supplement.tex`, `README.md` (Revision after round 145; Mathematical checks done for round 145); `review/REVIEW_round145_D3en.md`; `review/REVIEW_LOG.md` (Round 145); `numcheck15.py`.
- [md] `results/chaos/colour_selection.md`, Part K2 (l. 245–260).
- [own] Scripts `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref151/scripts/lemma26_check.py` and `jet.py`. Output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref151/lemma26_out.txt`. Build in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref151/build/`.
- [paper] Read at the cited passages:
  - [Shalova–Schlichting–Peletier, arXiv:2404.12293v1](https://arxiv.org/abs/2404.12293)
  - [Li–Wang–Arora, arXiv:2110.06914](https://arxiv.org/abs/2110.06914)
  - [Parsons–Rogers, arXiv:1510.07031v2](https://arxiv.org/abs/1510.07031)
- Not reachable: Katzenberger 1991 and Falconer 1983 (projecteuclid.org connection refused). A web search for Kat91 Cor. 5.1 found nothing relevant.
- [memory] Fenichel (1974, 1977) and Aulbach (LNM 1058, 1984): smoothness of the asymptotic phase near normally attracting manifolds of equilibria. Unchecked.

## Lesson
When a regularity lemma is added to discharge an assumption, check every downstream Taylor remainder against the regularity it actually delivers. A C² limit map makes ∂Φ·U only C¹, so a remainder written as O(|η|³) silently needs one more (Lipschitz) derivative, although o(|η|²) is all the argument uses.

A sentence that summarises an unread source ("this is what Falconer's theorem gives") should be phrased as what the read sources say about it.
