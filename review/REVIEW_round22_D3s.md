# Referee report D3s (round 22) — PARADIGM §5d: Lemma D‴.26a (statement (i)–(v), proof, Status, *Numerics*), the places it upgrades — Corollary D‴.26(c) and its Status, Proposition D‴.30(ii) and its Status — the §0 rows "| D‴.26 |", "| D‴.26a |", "| D‴.30 |", the §7 rows "| D‴.26 |", "| D‴.26a |", "| D‴.30 |", the §8 bullets "Done (Corollary D‴.26)", "Done (Lemma D‴.26a)", "Done (Proposition D‴.30)", `results/RESULTS_PARADIGM.md` (Lemma D‴.26a), LITERATURE_TRANSPORT [S59] and the audit bullet, PAPERS.md and README.md where they cite the lemma, `src/core/adaptive_floquet.py` (`CollapsedVariational.chang_grid`, `chang`, `slow_generator`), `src/core/measurement_noise.py` (`noise_decoupling`, `decoupled_covariance`), `experiments/exp_adaptive_floquet.py` Part D → `results/chaos/adaptive_floquet.md` (+ raw), `tests/test_adaptive_floquet.py` (the three lemma tests)

Branch `claude/festive-bardeen-eopej2` at 3e8f8e5, HEAD when the review started (`git log --oneline -5`: 3e8f8e5, 9204ea7, ac7795c,
637a051, 91658b9).

During the review the concurrent session committed a6ed653 (Proposition D‴.32, part 1) and e325236 (manuscript part 4d:
`manuscript/sections/thmDppp.tex`, `refs.tex`). It also left uncommitted edits to `src/core/adaptive_floquet.py` and
`src/core/measurement_noise.py` (57 insertions, 17 deletions). None of this is reviewed here. Everything below is about 3e8f8e5's tree,
which was extracted with `git archive 3e8f8e5` into `ref_D3s/tree` before those changes and is unaffected by them.

Every script imports nothing from the repository except where marked "cross-check" (the repository's object compared against mine).

Line numbers are 3e8f8e5's:
- PARADIGM: the lemma l. 3182–3226 (hypotheses l. 3183, (★) l. 3186, the identification l. 3188–3191, ε₀ l. 3191–3192, (i) l. 3194–3198,
  (ii) l. 3200–3202, (iii) l. 3204–3208, (iv) l. 3210–3214, (v) l. 3216–3226); the proof l. 3228–3267; the Status l. 3269–3281; the
  *Numerics* l. 3283–3295. Corollary D‴.26: (c) l. 3083–3097, "Hence" l. 3099–3103, the proof l. 3104–3139, the Status l. 3140–3157.
  Proposition D‴.30: (ii) l. 3606–3617 (the floor l. 3610), the proof of (ii) l. 3636–3647 (the citation l. 3642), the Status l. 3657–3672.
  §0 rows l. 130, 131, 135; §7 rows l. 6286, 6287, 6291; §8 bullets l. 6535–6539, 6540–6545, 6564–6570.
- RESULTS l. 642–653. LITERATURE_TRANSPORT: [S59] l. 79, the audit bullet l. 409–416. PAPERS.md l. 105, 114, 180–181, 191–192. README l. 57.
- Code: `adaptive_floquet.py` l. 24–26 (docstring), 254–301; `measurement_noise.py` l. 482–525; `exp_adaptive_floquet.py` l. 204–248 (Part D),
  348–373 (its report), 394–401 (`--add-D`); `tests/test_adaptive_floquet.py` l. 154–205.

No repository file was modified, and no git command touching the index or the tree was run. Python ran with `PYTHONDONTWRITEBYTECODE=1`
and one thread (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`); pytest with `-p no:cacheprovider` on the scratch tree. `--from-raw` on the
scratch tree rewrites `adaptive_floquet.md` byte-identically (md5 `f8356ea1cfe33fa36e48840fd5e217d2`, equal to
`git show 3e8f8e5:results/chaos/adaptive_floquet.md`).

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3s/scripts/`:

| script | what it does |
|---|---|
| `inst.py` (module) | independent instance: the exact representation at w = (6, 5, 4)/15, analytic Ĵ(x) and its x-derivatives to order 4 (Hermite polynomials), K(t) = Ĵ(x_m(t)) and its t-derivatives to order 4 (Faà di Bruno), the plant's slope a = f′(x_m); (★)'s K and Λ at η = 0 (K = Ĵ, Λ = σ₀I), finite η and ∞ |
| `check_inst.py` | the instance against finite differences; h_j against the md on both references; Ĵ against the repository's `whitened_jacobian` (cross-check) |
| `lin.py` (module) | Gauss–Legendre propagation (4 stages, order 8) of linear periodic systems; the right slow subspace (P = EY⁻¹/ε, forward at +ε, backward at −ε) and the left slow rows (backward), S from them; the periodic covariance by the exponential midpoint rule with Van Loan's block exponential per step and a Kronecker solve of the discrete Lyapunov equation |
| `check_partD.py` | Part D's deterministic columns at both η, k_e = 12.5–400: the full monodromy, the slow subspace and slow monodromy at ±ε, the parity per mode and in the md's normalisation, the first-order deviation, the mismatch against its first-order coefficient, |P − P₂|/ε³ (also with P₂'s Λ term flipped or its cubic term dropped) and |P − P₍₃₎|/ε⁴, the Magnus parameter (sup and integral forms), Liouville's fast sum |
| `check_reversal.py` | the parity's two halves as identities (P[K](−t; −ε) = P[K̃](t; ε), Φ_K̃(ε) = Φ_K(−ε)⁻¹, the full multipliers); the repository's `chang` against the propagated P (cross-check) |
| `check_parity_own.py` | the parity on references and gains of my own: n = 2 with a random trigonometric K; Λ symmetric; K = [J, 0] at finite η ≠ σ (σ = 0.5, η = 2); complex pairs (an indefinite D); a control without the symmetry; a mismatch; the truncations with an O(1) Λ |
| `check_magnus.py` | Ā₁ of (iii) (Magnus formula, own quadrature) and c₁ against the exact slow generator at ±ε, k_e = 100–800 |
| `check_noise.py` (+ `check_noise_hi.py`, `check_noise_conv.py`) | (v) at η = σ₀, ∞, 0: S from the left slow subspace and its equation, b against b₀, the periodic covariance at t = 0 (Richardson over N_c, 2N_c), its transform to (v, Y′), the md's columns, the fast jitter, the period-averaged floor against D‴.30(ii)'s formula and its range over the period, the decomposition of Cov(Y′) − Σ̄; k_e to 1600; the resolution study |
| `check_noise2.py` | η = 0 under the lemma's literal identification; the pointwise prediction of the floor over the period; b's first-order term; the cross term of Cov(Y) − Cov(Y′) |
| `check_tests.py` | the three tests' asserted quantities; 12 mutants (of the tests' formulas, or of the code patched in memory) |
| `check_misc.py` | the D-symmetry with a mismatch in the error channel; P₃ against K⃛; the mismatch coefficient per mode |
| `check_jitter_next.py` | the next term of Cov(Y) − Cov(Y′) − (ε/2)KᵀK at t = 0, against −(ε²/2)(d(KᵀK)/dt + 2aKᵀK), k_e = 100–400 |

Also in `ref_D3s/`: `tests.out` (3 passed, 18.83 s), `from_raw.out`, `md5_before.out`, `md5_after.out`. Part D was not rerun with the
committed code (`--add-D`); it was recomputed by independent code and every printed number reproduces.

## Verdicts

| object | verdict |
|---|---|
| **(i) the slow subspace** | **holds.** Re-derived: 𝒢 (∫₀^∞e^{−σ}g(t − εσ)dσ, causal for ε > 0, anti-causal for ε < 0, ‖𝒢‖ ≤ 1, 𝒢g − g = −ε𝒢ġ); the self-map ‖𝒢F(P) + K‖ ≤ \|ε\|‖K̇‖ + ε²β ≤ r = 2\|ε\|(‖K̇‖ + β) (it needs only \|ε\| ≤ 2; r ≤ 1 gives ‖P‖ ≤ ‖K‖ + 1, which is what β uses); the Lipschitz constant ε²(\|Λ\| + 2‖K‖(‖K‖ + 1)); the recursion (P₂ = −(K̈ + KKᵀK + KΛ), the Λ term's sign right; P₃ = K⃛ + 2K̇KᵀK + KK̇ᵀK + 2KKᵀK̇ + 2K̇Λ); ρ = O(ε⁴) needs K⁗ (C⁴ is the right hypothesis); the block-triangular form including the εṖY term. Uniqueness: the Lipschitz estimate needs only ‖P‖, ‖P′‖ ≤ ‖K‖ + 1, so P is unique even in ‖P + K‖ ≤ 1, not globally (the quadratic equation has other periodic solutions of size ε⁻²); "unique with ‖P + K‖ = O(ε)" is right, and the parity's "by uniqueness" needs only the ball (P̃ lies in K̃'s ball of the same radius). Numerically: the repository's P equals the propagated one to ≤ 8.2·10⁻¹⁵, |P − P₍₃₎|/ε⁴ is bounded. The md's and the test's truncation checks cannot see P₂'s cubic and Λ terms (D3s-6). |
| **(ii) the slow multipliers** | **holds.** The logarithmic norm of −ε⁻¹(I + ε²PKᵀ) is ≤ −ε⁻¹ + ε‖P‖‖K‖, so the fast block's monodromy is ≤ e^{−T₀/(2ε)} once ε²‖P‖‖K‖ ≤ ½; the slow multipliers are ≥ e^{−εT₀‖A‖}; the q slow ones are the Y block's. Liouville: the k_e⁻² term is zero (⟨tr K̇Kᵀ⟩ = 0), and the k_e⁻³ coefficient is C₃ = ⟨\|K̇\|²⟩ − ⟨‖KKᵀ‖²_F⟩ − ⟨tr KΛKᵀ⟩ = 10.742, measured 10.734 at k_e = 400. The fast multipliers (e^{−k_eT₀}) are below double precision, so their bound is only analytic. |
| **(iii) the slow generator** | **holds.** The truncation, Grönwall (O(ε⁵) for the monodromies), the analyticity of log Φ₍₃₎(ε)/(εT₀) and the division by ε are right. Ā₁ with the commutator coefficient 1/(2T₀) is confirmed against the exact generator (the symmetric difference converges to it like ε²). Two imprecise sentences (D3s-9): the Magnus series "of εA collected by powers of ε" (A is not analytic in ε; the proof uses εA₍₃₎), and "conjugate … up to O(ε⁴)" (D‴.26's averaging is carried to Ā₂). Moan–Niesen: the statement's εT₀‖A‖_∞ < π (sup of the 2-norm) is sufficient for the theorem's ∫₀^{T₀}‖εA‖₂ < π. The code's "Magnus parameter" is the sup form over the 512 samples: 0.480 at k_e = 12.5 (the integral form 0.219). |
| **(iv) the parity** | **holds as an identity, and more.** Re-proved: under Mᵀ = DMD⁻¹ the system with M_ε(−t) has monodromy D⁻¹Ψ(T₀)ᵀD, so the same multipliers; D = diag(−ε²I, I, (σ₀/η)I) is right for (★) at finite η (diag(−ε²I, I) at η = ∞, and for η = 0 without ζ); P[K](−t; −ε) is K̃'s P at ε; Φ_K̃(ε) = Φ_K(−ε)⁻¹ exactly as matrices, not only up to conjugacy; the principal logarithm of the inverse is minus the logarithm. The identity gives the whole complex spectrum of Ā(ε) even, not only the simple real rates (D3s-7). Confirmed on the instance (per mode ≤ 3.4·10⁻¹⁰) and on six cases of my own, including complex pairs; a control without the symmetry breaks it at O(ε). The real error channel is εẇ = −(1 − εa)w − KY with a = ∂_x(f − f̂), the lemma's sign. The mismatch breaks the parity at first order with the predicted coefficient 2εℓ_j⟨aKᵀK⟩r_j. "μ_j¹ = μ_j³ = 0" follows, and is consistent with c_j (measured c₁/λ₁ = −108.16 against the formula's value from the exact Ā₁, Ā₂, and the md's −108.2). "Which holds when a enters only the noise" is misleading (D3s-7). ε₀ also depends on the gap (D3s-4). |
| **(v) the slow noise** | **holds asymptotically.** Re-derived: S's equation and its anti-causality (independently, S from the left slow subspace satisfies it to ≤ 10⁻¹⁴); b = −(K̇ + aK)ᵀ + O(ε), with the next term −ε(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ; the moments' orders and the averaging argument for Z = Σ′ − Σ̄; Cov(Y) = Σ̄ + (εN₀/2)KᵀK + O(ε²N₀). The noise model is D‴.30's: the input −(k_e + a)ν into the error, Ĵᵀν into the law, ζ noise-free in both. The floor is the period average of tr(H Cov y), not of tr(Ĵ(t)ᵀĴ(t)Cov y(t)), whose fast part would be ⟨\|Ĵ\|⁴⟩. Four gaps. The floor holds only as a period average (D3s-2). η = 0 is outside the identification, where Ā₀ is singular (D3s-3). The O(ε²N₀) remainder exceeds the jitter below k_e ≈ 50 (D3s-2). "Exactly … −ε(K̇ + aK)ᵀν" is overstated (D3s-5). |
| Proof | **correct.** Imprecisions: D3s-7, D3s-9. |
| *Status* | **accurate on what is proved**, with D3s-3 and D3s-4. "The slow noise identified exactly as the high-passed −ε(K̇ + aK)ᵀν" is overstated (D3s-5). The novelty audit omits the stochastic and time-varying Chang precedents (D3s-10). |
| *Numerics* | **every number reproduces** with independent code: to all printed digits for the truncation, the first-order deviation, the mismatch, the Magnus parameter and the whole noise table. "Halving per doubling" holds (2.08, 2.06, 2.02, 2.00, 2.00), and so does "its constant set by the regressor's third derivative" (max\|P₃\| = 3341, max\|K⃛\| = 3346). Two readings do not hold. "Dominated by the weakest mode's rate correction, O(ε²) on the instance" is wrong (D3s-1): the slow noise's correction is 2.5 times the rate correction and of the opposite sign, and beyond k_e = 400 the deviation is O(ε). The normalisations are relative to the largest rate (D3s-8). |
| Corollary D‴.26(c) and its Status | **the upgrade holds** for simple eigenvalues under ∂_xf̂ = ∂_xf. Not stated where "proved" is claimed: K ∈ C⁴ (x_m ∈ C⁴), H_L(w) ≻ 0 for the stability "Hence", and a threshold that depends on the spectral gap (D3s-4). The local exponential stability of the real loop then follows from the multipliers (Floquet–Lyapunov; the field is C¹ near the lift because the cells are constant there). |
| Proposition D‴.30(ii) and its Status | **proved for the linearised loop under white noise**, for the period-averaged floor at η ∈ (0, ∞], with the relative remainder O(1/k_e). η = 0 needs K = Ĵ, Λ = σ₀I (D3s-3). The pointwise "+O(k_e⁻²)" fails at first order through the jitter's oscillation (D3s-2). |
| §0 rows, §7 rows, §8 bullets, RESULTS, PAPERS, README | **numbers hold**; wording D3s-5, D3s-8, D3s-11. |
| LITERATURE [S59], audit bullet | **bibliographic details right** (Moan–Niesen FoCM 8(3):291–301, 2008, Theorem 3, the constant π sharp; Chang's DOI 10.1137/0503050 exists). The audit omits precedents (D3s-10). |
| `adaptive_floquet.py` (`chang_grid`, `chang`, `slow_generator`) | **holds.** `chang` equals the propagated slow subspace to ≤ 8.2·10⁻¹⁵ for both signs of ε at both η. The slow block's rates are the full monodromy's. The module docstring's wording: D3s-11. |
| `measurement_noise.py` (`noise_decoupling`, `decoupled_covariance`) | **holds.** b, Cov(v), E[vY′ᵀ], Cov(Y′) and the floor at t = 0 are reproduced to the printed digits by independent code. `decoupled_covariance` does not apply at η = 0: the ζ block has a zero rate (D3s-3). |
| `exp_adaptive_floquet.py` Part D, md, raw | **holds**: the regeneration is byte-identical. job_D's "even" and "first" are normalised by the largest rate (D3s-8), and the floor is evaluated only at t = 0, where the jitter is 0.07% of it (D3s-2). |
| `tests/test_adaptive_floquet.py` (the three lemma tests) | **pass (3/3, 18.83 s); partly vacuous** (D3s-6). Two P₂ mutants pass, and so do Σ̄ without the slope and the slope dropped everywhere. The parity test is strong: it catches a generator even to 8·10⁻⁹. |

## Confirmed issues

**D3s-1 (minor; *Numerics* l. 3293–3294 "‖Cov(Y′) − Σ̄‖/‖Σ̄‖ from 0.36 to 0.0011 (dominated by the weakest mode's rate correction, O(ε²) on
the instance)"; RESULTS l. 650; §0 row D‴.26a l. 131 "the noise's moments at their stated orders").** Both halves of the reading are wrong
(`check_noise.py`, `check_noise_hi.py`, `check_noise_conv.py`; η = σ₀).
- *The mechanism.* To k_e = 400 the largest deviation sits on the weakest mode: the entry (5, 5), also Σ̄'s largest. Its cause is the
  slow noise's attenuation, ⟨bbᵀ⟩ < Q, not the rate correction. Σ̄ was recomputed with ⟨bbᵀ⟩ in place of Q̄ ("noise only") and with the
  exact slow generator log Φ/(εT₀) in place of Ā₀ ("rate only"). Below, the weakest mode's variance (y block, H's eigenbasis) is relative
  to Σ̄:

  | k_e | exact Cov(Y′) | noise only | rate only |
  |---|---|---|---|
  | 12.5 | −34.9% | −56.4% | +47.2% |
  | 25 | −17.7% | −28.6% | +15.1% |
  | 50 | −6.11% | −10.0% | +4.29% |
  | 100 | −1.70% | −2.79% | +1.12% |
  | 200 | −0.436% | −0.719% | +0.284% |
  | 400 | −0.110% | −0.181% | +0.0713% |

  The noise's correction is 2.5 times the rate's and has the opposite sign. The exact variance lies below Σ̄; the rate correction
  alone would raise it. In max-norm the noise part is 0.562 … 1.80·10⁻³ and the rate part 0.448 … 6.90·10⁻⁴. The mechanism is the
  filter's attenuation of the regressor's harmonics, b = −(K̇ + aK)ᵀ − ε(K̈ + …)ᵀ + … (D3s-5). It acts on the noise as on the rate:
  about −290ε² against +114ε² on this mode (−c₁/λ₁ = 108), from the k_e = 400 values.
- *The order.* Beyond k_e = 400 the deviation is O(ε). ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ = 1.109·10⁻³, 3.290·10⁻⁴, 1.652·10⁻⁴ at k_e = 400, 800,
  1600; after 3.89 and 3.97 the successive ratios drop to 3.37 and 1.99. The maximum moves off the weakest mode, to the entry (0, 1).
  - *The (0, 1) entry is O(ε) throughout:* +1.283·10⁻³, +6.526·10⁻⁴, +3.289·10⁻⁴, +1.650·10⁻⁴ at k_e = 200, 400, 800, 1600
    (relative to max|Σ̄|).
  - *It is resolved.* Richardson over (N_c, 2N_c) and over (2N_c, 4N_c) gives 3.28960·10⁻⁴ and 3.28931·10⁻⁴ at k_e = 800, and
    1.65185·10⁻⁴ and 1.65044·10⁻⁴ at 1600 (`check_noise_conv.py`). The weakest mode's own entry is not resolved beyond k_e = 400,
    so it enters only as a bound, below the (0, 1) entry.
  - *The averaged corrections* (noise and rate together) leave a remainder of 9.67·10⁻⁴, 4.90·10⁻⁴, 2.46·10⁻⁴ (k_e = 100–400),
    which halves per doubling.
  - This is the O(ε) relative term that the lemma's Cov(Y′) = Σ̄ + O(ε²N₀) allows. The weakest mode's large O(ε²) terms hide it to
    k_e = 400. That is the lesson of D3q-1: an order read off a limited range.

*Fix:* "(dominated to k_e = 400 by the weakest mode's slow-noise attenuation ⟨bbᵀ⟩ < Q — 2.5 times its rate correction and of the
opposite sign — and O(ε) beyond, 3.3·10⁻⁴ and 1.7·10⁻⁴ at k_e = 800, 1600, as the lemma's bound allows)". Add k_e = 800, 1600 to Part D,
or print the decomposition.

**D3s-2 (minor; (v) l. 3224–3226 "hence Cov(Y)(t) = Σ̄ + (εN₀/2)K(t)ᵀK(t) + O(ε²N₀) — the slow part plus the fast jitter — and
Proposition D‴.30(ii)'s floor … holds"; D‴.30(ii) l. 3610 "E‖f̂ − f‖²_ρ = (N₀/2k_e)[…] + O(k_e⁻²)"; the Status's "Checked" l. 3281;
*Numerics* l. 3294–3295; the md Part D's floor columns).** The fast part of the floor is a period average, and the jitter is not tested.
- *The floor is a period average.* The proof (l. 3266) proves the period average of tr(H Cov y). The statement (v) and D‴.30(ii) state
  E‖f̂ − f‖²_ρ without it. At time t the fast part is (N₀/2k_e)tr(HĴ(t)ᵀĴ(t)): an O(k_e⁻¹) oscillation, of the same order as the floor.
  - *Exact pointwise floor* (`check_noise.py`, `check_noise_hi.py`; η = σ₀): tr(H Cov y(t)) over the period lies in [0.9937, 1.0098] of
    D‴.30(ii)'s value at k_e = 400, [0.9945, 1.0099] at 800 and [0.9948, 1.0099] at 1600. The width does not shrink.
  - *Period average:* 0.99952, 0.99990, 0.99999 at the same k_e.
  - *Predicted pointwise floor:* tr(H(Σ̄ + (ε/2)Ĵ(t)ᵀĴ(t))) ranges over [0.99481, 1.00987] of the formula (`check_noise2.py`). The
    jitter's share is 0.52% on average, 1.51% at its maximum and 0.07% at t = 0.
  - Pointwise, the floor's error is O(1/k_e) (±1% of the floor here), not O(k_e⁻²). The averaged deficit is 0.77% at k_e = 100 and
    0.048% at 400.
- *The jitter is not tested.* Part D compares the floor at t = 0 only, where the jitter is 0.07% of it ("fast share 0.001"), so that
  column tests the slow part. Compared directly, Cov(Y) − Cov(Y′) differs from (ε/2)K(0)ᵀK(0) by 4.60, 2.89, 1.07, 0.315, 0.162, 0.081
  (relative, k_e = 12.5–400), and by 0.041 and 0.020 at 800 and 1600. Below k_e ≈ 50 the O(ε²N₀) remainder is larger than the "fast
  jitter" it accompanies.
  - *The remainder's leading part* is −(ε²N₀/2)(d(KᵀK)/dt + 2aKᵀK). The cross term −ε(SC + CᵀSᵀ) is ≈ −ε²N₀(d(KᵀK)/dt + 2aKᵀK)
    (3.556·10⁻⁵ against 3.565·10⁻⁵ at k_e = 400, `check_noise2.py`). The jitter's own correction, ε²SVSᵀ − (εN₀/2)KᵀK, cancels half
    of it.
  - *Check.* In max-norm the remainder is 1.779·10⁻⁵ against 1.783·10⁻⁵ at k_e = 400. The entrywise difference falls like ε³:
    2.66·10⁻⁴, 3.45·10⁻⁵, 4.35·10⁻⁶ at k_e = 100, 200, 400 (`check_jitter_next.py`).
  - *Its size.* Relative to the jitter it is of order ε|d(KᵀK)/dt|/|KᵀK|: it carries the regressor's frequency.
- *D‴.30's own functional is never computed in Part D.* The period-averaged floor is 0.7587, 0.9044, 0.9710, 0.99228, 0.99804, 0.99952
  of the formula at k_e = 12.5–400 (η = σ₀).

*Fix:* in (v) write "the period-averaged floor", and in D‴.30(ii) "E‖f̂ − f‖²_ρ averaged over a period (at time t the fast part is
(N₀/2k_e)tr(HĴ(t)ᵀĴ(t)))". Add the next term of Cov(Y) − Cov(Y′) − (εN₀/2)KᵀK, −(ε²N₀/2)(d(KᵀK)/dt + 2aKᵀK). In Part D print the period
average, and Cov(Y) − Cov(Y′) against the jitter.

**D3s-3 (minor; the identification l. 3188–3191; D‴.30's Status l. 3658 "(ii) proved for the linearised loop under white noise
(Lemma D‴.26a(v) …"; the §0 row D‴.30 l. 135; PAPERS l. 181; `measurement_noise.decoupled_covariance`).** The identification lists
K = [Ĵ, 0], Λ = [[σ₀I, −σ₀I], [−ηI, ηI]] and the case η = ∞. It omits η = 0, the fixed prior, whose filter factor h_j/(h_j + σ₀) is part of
the floor D‴.30(ii) now calls proved.
- *The listed Λ at η = 0.* The ζ block has a zero rate: Ā₀ has 9 eigenvalues with zero real part (`check_noise2.py` (a)). (v)'s
  "If Ā₀ is Hurwitz" fails, and Σ̄ does not exist.
- *The right identification.* With K = Ĵ, Y = y, Λ = σ₀I the lemma applies: Ā₀ = −(H + σ₀I), margin 7.5·10⁻³. The period-averaged floor
  is 0.90327, 0.99219, 0.99952 of the η = 0 formula at k_e = 25, 100, 400 (`check_noise.py zero`).
- *The code.* `decoupled_covariance` has the same gap: through `CollapsedVariational` it builds a ζ block at every finite η, including 0.

*Fix:* add "(η = 0, the fixed prior: K = Ĵ, Y = y, Λ = σ₀I)" to the identification.

**D3s-4 (minor; l. 3191–3192 "There is ε₀ > 0, depending on the C⁴ norms of K, on |Λ|, ‖a‖_{C¹} and T₀"; Corollary D‴.26(c)
l. 3083–3103 and its Status l. 3140).** The stated dependence is right for (i)–(iii) and incomplete for (iv) and (v).
- *(iv).* The continued eigenvalue stays simple only while the perturbation is small against the gap of −λ_j in the spectrum of Ā₀,
  and c_j contains Σ_{i≠j}(ℓ_jĀ₁r_i)(ℓ_iĀ₁r_j)/(λ_i − λ_j). So ε₀ and the O(ε⁴) constant depend on the gap and on the eigenvector's
  condition, which no C⁴ bound on K controls: two eigenvalues can merge at fixed norms. On the instance (η = σ₀) the smallest gaps are
  1.42·10⁻⁴ and 1.82·10⁻⁴, against ε·max|Ā₁| = 1.37·10⁻²/k_e: the two are equal at k_e ≈ 100 (`check_misc.py`, `check_magnus.py`).
- *(v).* Σ̄ scales like the inverse Hurwitz margin, and so does the constant of the relative remainder. Here λ₁ = 2.15·10⁻⁴ against
  λ_max = 0.194, and the weakest mode's variance is 35% off at k_e = 12.5.
- *Corollary D‴.26(c).* Its "Hence … exponentially stable" uses λ_1^− > 0, i.e. H_L(w) ≻ 0 (Ā₀ Hurwitz). The lemma needs K ∈ C⁴, i.e.
  x_m ∈ C⁴. Neither is stated where (c) is called proved.

*Fix:* "ε₀ and the constants depending on …, and, for (iv), on the distance of −λ_j from the rest of the spectrum of Ā₀ and the condition of
its eigenvector, for (v) on Ā₀'s stability margin". Add "x_m ∈ C⁴ and H_L(w) ≻ 0" to D‴.26(c).

**D3s-5 (minor; the Status l. 3277 "the slow noise identified exactly as the high-passed −ε(K̇ + aK)ᵀν"; (v) l. 3221 "driven by the
high-passed noise alone"; the §0 row l. 131 "Chang's second transformation decouples the slow coordinates from the noise's fast part
exactly, their noise the high-passed −ε(K̇ + aK)ᵀν"; the §8 bullet l. 6542–6543; LITERATURE l. 414–415; D‴.30's Status l. 3658–3659 "their
noise the high-passed one").** The exact slow noise is εb(t)ν, with b = ε⁻¹[Kᵀ − ε⁻¹(1 + εa)S − εSPKᵀ]. −(K̇ + aK)ᵀ is only its leading term.
- *Size.* On the instance |b − b₀|/|b₀| = 62.7%, 40.3%, 22.5%, 11.7%, 5.90%, 2.96% at k_e = 12.5–400 (the md: 63%–3.0%).
- *Next term.* b₁ = −(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ, derived here. |b − b₀ − εb₁|/ε² = 2630, 3272, 3337 at k_e = 25, 100, 400, so it is bounded
  (`check_noise2.py`). The K̈ in it attenuates ⟨bbᵀ⟩, which is D3s-1's dominant correction.
- *What is decoupled.* The drift is decoupled, not the noise: Y′ and v are driven by the same W, and E[vY′ᵀ] = −N₀bᵀ ≠ 0.

*Fix:* "the slow noise identified exactly as εb(t)ν, b = −(K̇ + aK)ᵀ − ε(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ + O(ε²) (the high-passed noise at leading
order)", and "decouples the slow coordinates' drift from the fast ones".

**D3s-6 (minor; `tests/test_adaptive_floquet.py` l. 154–205; the Status's "Checked" l. 3278–3281; the md's truncation column).** The tests'
asserted quantities (`check_tests.py`, at (2, 1, 1)/4):

| test | quantity | measured | tolerance |
|---|---|---|---|
| chang | invariance residual (k_e = 25, 50) | 2.2·10⁻¹⁵, 1.2·10⁻¹⁵ | < 10⁻¹² |
| chang | slow block vs full | 1.1·10⁻¹¹, 6.8·10⁻¹¹ | < 10⁻⁸ |
| chang | \|P − P₂\| ratio k_e = 25 → 50 | 6.819 | (6, 10) |
| parity | \|r(ε) − r(−ε)\|/max r (η = σ₀, ∞) | 7.2·10⁻¹⁴, 4.4·10⁻¹⁴ | < 10⁻¹⁰ |
| parity | first-order deviation | 3.42·10⁻², 3.53·10⁻² | > 10⁻² |
| parity | mismatch | 1.33·10⁻², 1.38·10⁻² | > 10⁻² |
| noise | coupling (k_e = 50, 100) | 2.1·10⁻¹⁷, 6.0·10⁻¹⁸ | < 10⁻¹⁴ |
| noise | \|b − b₀\|/\|b₀\|: ratio, value at k_e = 100 | 1.922, 0.117 | (1.7, 2.3), < 0.15 |
| noise | Cov(v)2ε³ − 1; E[vY′ᵀ] + bᵀ; Cov(Y′) − Σ̄ | 2.07·10⁻³; 1.66·10⁻²; 1.72·10⁻² | 10⁻²; 3·10⁻²; 3·10⁻² |

Mutants that pass (patched in memory):
- P₂ with the Λ term's sign flipped: the ratio stays 6.819.
- P₂ without KKᵀK: 6.822.
- Σ̄ built from Q without the slope, ⟨K̇ᵀK̇⟩: Cov(Y′) 2.20·10⁻², under 3·10⁻².
- The slope dropped everywhere (FPRIME = 0): every assertion passes. The noise test compares the code's quantities with each other,
  never b₀ with an independent a.

The md's |P − P₂|/ε³ column is equally blind. Flipping the Λ term changes it by (6–7)·10⁻⁵ (relative), dropping KKᵀK by 0.1–0.2%
(`check_partD.py`): the O(ε³) term (K⃛, ~3.3·10³) dominates both at every k_e of the table. P₃ and the recursion beyond k = 2 are checked
nowhere.

Mutants that are caught:
- P₂ without K̈ (ratio 3.56).
- chang's F with +Λ (slow vs full 2.1·10⁻⁵).
- The causal filter at negative ε (parity 4.9·10⁻⁴).
- The second-order averaged generator in place of the slow block (parity 8.0·10⁻⁹, caught by a factor 80).
- S by the causal filter (coupling 1.2·10⁻²).
- b without (1 + εa) (ratio 1.51).
- The noise input −k_eν in the exact covariance (Cov(Y′) 3.63·10⁻², a margin of 1.2).

*Fix:*
- Assert that |P − P₍₃₎|/ε⁴ is bounded.
- Add a case with an O(1) Λ. In my case B the stated P₂ stays at 33–53, while the flipped sign grows 155 → 922 as ε halves
  (`check_parity_own.py`).
- Test b against b₀ + εb₁ at O(ε²).
- Test Q (and so Σ̄) against an independent quadrature with an analytic Ĵ̇ and a.
- Replace the 3·10⁻² tolerances by the O(ε) scaling between two k_e.

**D3s-7 (cosmetic; (iv) l. 3212–3214; the proof l. 3248–3257; Corollary D‴.26(c) l. 3089–3092).**
- *"Which holds when a enters only the noise"* (l. 3213). The D-symmetry persists with a mismatch in the error channel: with the plant's
  slope as the mismatch, max|Mᵀ − DMD⁻¹|/|M| = 1.6·10⁻¹⁶ (`check_misc.py` (a)). What a breaks is the invariance of P's equation under
  (t, ε) ↦ (−t, −ε), through the εaP term (l. 3256); Corollary D‴.26's proof says as much ("it persists with a ≢ 0").
  *Fix:* "(which persists with a mismatch in the error channel; the parity needs a ≡ 0 there, below)".
- *Scope.* The proof gives more than (iv) states: eig Ā(−ε) = eig Ā(ε) as multisets of complex numbers.
  - Φ_K̃(ε) = Φ_K(−ε)⁻¹ holds exactly (5.6·10⁻¹⁶), and K̃'s slow multipliers are K's (7.3·10⁻¹⁵, `check_reversal.py`).
  - The complex eigenvalues of log Φ/(εT₀) at ±ε agree to ≤ 5.2·10⁻¹³ on the instance. They agree to 1.8·10⁻¹⁴ in a case with complex
    pairs and an indefinite D (`check_parity_own.py` case D).
  - Because (iv) states only the simple real rates, it does not cover the Corollary's "at a degenerate one only the rates are even"
    (l. 3090) as written.
  *Fix:* state "the spectrum of Ā(ε) is even".
- *A shorter route.* The full system at −ε is the time reversal of K̃'s at +ε, so with the D-symmetry its slow multipliers are the
  inverses of K's at +ε. This needs no P: 6.4·10⁻¹⁵ (`check_reversal.py`).

**D3s-8 (cosmetic; *Numerics* l. 3287–3289; the §0 row l. 131 "(their first-order deviation up to 9.5%)"; the §7 row l. 6287; RESULTS
l. 647–648; job_D l. 233–235; the Status l. 3278).**
- *The normalisation.* "The rates at ±ε agree to ≤ 1.6·10⁻¹² (relative)" and "their deviation from the first-order prediction is
  9.5·10⁻²–1.4·10⁻⁴" are max|Δ|/max|ν|, relative to the largest rate.
- *Per mode.* The ±ε agreement is ≤ 3.4·10⁻¹⁰. The slowest mode's first-order deviation is 29.3% (η = σ₀) and 29.6% (η = ∞) at
  k_e = 12.5, and 12.3% at 25. "(Up to 9.5%)" understates it by a factor of three; the parity's contrast is 10⁹.
- *η = ∞.* Its values (9.8·10⁻²–1.5·10⁻⁴; with the mismatch 1.1·10⁻²–3.25·10⁻⁴) are not quoted.
- *The Status.* "(on the asymmetric reference at η = σ₀, k_e = 12.5–400 …)" merges Part D (η = σ₀ and ∞, w = (6, 5, 4)/15) with the tests
  ((2, 1, 1)/4, k_e = 25–100).

*Fix:* quote per-mode values for both η, and describe Part D and the tests separately.

**D3s-9 (cosmetic; (iii) l. 3205–3208; the proof l. 3241–3245).**
- *"The higher coefficients those of the Magnus series of εA collected by powers of ε."* A depends on ε through P, which has only an
  asymptotic expansion. Collecting by powers is legitimate for εA₍₃₎, a polynomial in ε, which is what the proof uses.
- *"Conjugate to Ā(ε) up to O(ε⁴)."* Corollary D‴.26 carries the near-identity averaging only to Ā₂ (the slaving to L₃), so its
  generator is conjugate up to O(ε³).
- *"Its Taylor coefficients are obtained by … or the near-identity averaging."* The averaging gives a conjugate generator (by
  I + εU(0) + …), whose coefficients differ from the Taylor coefficients by commutators with Ā₀.
- *On the instance.* The Magnus Ā₁ (max entry 1.37·10⁻²) is not the md's averaging Ā₁ (5.6·10⁻³); both have ℓ_jĀ₁r_j = 0. The Magnus Ā₁
  is the commutator term alone: ⟨KᵀK̇⟩ = 2.5·10⁻¹⁶, an exact differential in one dimension, even on the asymmetric reference.

*Fix:* write "of εA₍₃₎", "up to O(ε³) with D‴.26's Ā₂", and "conjugate generators, with the same eigenvalues".

**D3s-10 (minor; the Status l. 3273–3277; the §0 row's last column l. 131; the audit bullet l. 409–416; [S59] l. 79).** The audit omits
the closest precedents.
- *The noise part (v).*
  - Chang's transformation of linear singularly perturbed systems with white-noise inputs, with the slow subsystem's transformed noise
    and the expansion of its covariance: Kokotović–Khalil–O'Reilly ([S9]), Ch. 4 "Stochastic linear filtering and control" and App. A
    (singularly perturbed systems driven by white noise).
  - A. H. Haddad and P. V. Kokotović, "Stochastic control of linear singularly perturbed systems", IEEE TAC 22(5):815–821, 1977.
  - Stochastic slow manifolds and their covariance: N. Berglund and B. Gentz, J. Differential Equations 191(1):1–54, 2003, and
    *Noise-Induced Phenomena in Slow-Fast Dynamical Systems*, Springer, 2006.
- *The time-varying slaving (i).*
  - KKO Ch. 5 "Linear time-varying systems".
  - X. Yang and J. J. Zhu, "Chang transformation for decoupling of singularly-perturbed linear slowly time-varying systems", 51st IEEE
    CDC, 2012: the differential Riccati and Sylvester equations of the transformation are solved by a contraction mapping, the
    structure of (i) and (v). The same authors' "A generalization of Chang transformation for linear time-varying systems" is at the
    49th CDC, 2010.
  - The periodic case is a specialisation: a bounded solution of a periodic equation that is unique in a ball is periodic.
- *The parity.* The identity Φ_K(−ε) = Φ_K̃(ε)⁻¹ is the adjoint relation Φ*_h = Φ_{−h}⁻¹ of symmetric methods ([S56]). "New: the parity
  as an identity at negative ε" is that argument applied to the slow block.

The verdict "classical tools; a rigorous foundation" is fair. What remains specific is the application: the exact slaving of this loop,
the cancellation Kᵀ − (1 + εa)(Kᵀ + εK̇ᵀ) = −ε(K̇ + aK)ᵀ that makes the slow noise high-passed, and through them the proofs of D‴.26(c) and
D‴.30(ii). The bibliographic details above were checked by a web search. The content is from memory and the abstracts; Haddad–Kokotović
and Yang–Zhu come only from reference lists and an abstract snippet. *Fix:* cite them in [S59] and in the audit, and restrict "new"
accordingly.

**D3s-11 (cosmetic; wording).**
- *README l. 57.* "Perfect tracking with an adaptive prior is exponentially stable in the real loop" should say locally, for each m:
  the basin is not proved, as the §8 bullet says. "An expansion in 1/k_e with only even powers" holds in slow time; the real-time
  exponents (1/k_e)[λ_j + c_j/k_e² + …] have odd powers.
- *`adaptive_floquet.py` l. 25–26.* "Whose eigenvalues are exactly the slow Floquet exponents (in slow time)": the eigenvalues' real
  parts are minus the rates, as `slow_generator`'s own docstring says.
- *Corollary D‴.26(c) l. 3090–3092.* It still explains the parity by the count "every term of Ā_k has the parity (−1)^k". After the
  lemma the explanation is the identity Ā[K̃](ε) = Ā[K](−ε).

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.** `--from-raw` rewrites the md byte-identically (md5 `f8356ea1…`). The three lemma tests pass in 18.83 s
   (8.85, 4.98 and 4.61 s).
2. **The instance** (`check_inst.py`).
   - f̂ − f = 0.0 on [−4, 4].
   - Ĵ against central differences of f̂ in the unconstrained parameters: 2.3·10⁻¹⁰ (relative). Its t-derivatives to order 4 against
     differences: ≤ 2.0·10⁻⁶.
   - h_j = 4.3720·10⁻⁴ … 0.18694 on the asymmetric reference and 3.0798·10⁻⁴ … 0.18336 on the symmetric one, the md's values.
   - Ĵ equals the repository's to 2.2·10⁻¹⁶ (cross-check).
3. **(i)** (`check_partD.py`, `check_reversal.py`, `check_parity_own.py`, `check_misc.py`).
   - The propagated slow subspace (Gauss–Legendre, order 8, 4096 steps) gives the repository's Fourier-contraction P to ≤ 8.2·10⁻¹⁵
     (relative), at ε = ±1/12.5, ±1/50, ±1/400 and both η.
   - |P − P₂|/ε³ on the 512-grid is 1849.7804, 2621.0622, 3074.0166, 3268.2155, 3305.7286, 3332.2482 (η = σ₀) and 1849.8753 … 3332.4406
     (η = ∞): the md's values to all printed digits.
   - |P − P₍₃₎|/ε⁴ = 32126, 45347, 53614, 56934, 57960, 58240 (k_e = 12.5–400) is bounded, so P₃ and the recursion at k = 3 hold.
   - max|P₃| = 3341.2 and max|K⃛| = 3346.3 (max|P₃ − K⃛| = 22.7): the md's "its constant set by the regressor's third derivative"
     holds.
   - On my cases (n = 2, q = 3; Λ symmetric; K = [J, 0] at σ = 0.5, η = 2), |P − P₍₂₎|/ε³ → 30.0, 32.8, 10.1 and |P − P₍₃₎|/ε⁴ → 128,
     164, 31, both bounded. With P₂'s Λ term flipped the ratio grows like 1/ε: 155, 252, 471, 922 and 24, 36, 64, 123 as ε halves. The
     sign −KΛ is right.
4. **(ii)** (`check_partD.py`).
   - The slow block's rates along my subspace equal my full monodromy's to ≤ 5.5·10⁻¹⁰.
   - Liouville: (fast sum − (k_e − tr H/k_e))k_e³ = 7.01, 9.24, 10.28, 10.62, 10.71, 10.73 → C₃ = 10.742 at η = σ₀ (10.746 at η = ∞).
   - The fast sum from the monodromy equals n/ε + ε⟨tr PKᵀ⟩ from P to 10⁻⁹.
5. **(iii)** (`check_magnus.py`, `check_partD.py`).
   - Ā₁ by the Magnus formula (own quadrature) against (Ā(ε) − Ā(−ε))/(2ε): 6.97·10⁻², 1.77·10⁻², 4.45·10⁻³, 1.12·10⁻³ (relative,
     k_e = 100–800), falling by 4 per doubling (the ε²Ā₃ term). ℓ_jĀ₁r_j ≤ 1.4·10⁻¹⁶.
   - The Magnus parameter (sup form) is 0.4800, 0.2470, 0.1246, 0.0624, 0.0312, 0.0156 (η = σ₀) and 0.4764 … 0.0155 (η = ∞), the md's
     values. The integral form ∫‖εA‖₂ is 0.2187 at k_e = 12.5.
6. **(iv)** (`check_partD.py`, `check_reversal.py`, `check_parity_own.py`, `check_magnus.py`).
   - *The identities.* P[K](−t; −ε) = P[K̃](t; ε) to ≤ 6.5·10⁻¹⁵. Φ_K̃(ε) = Φ_K(−ε)⁻¹ to ≤ 5.6·10⁻¹⁶. The slow multipliers of K and K̃
     at +ε agree to ≤ 7.3·10⁻¹⁵. The full system at −ε has the inverses of K's slow multipliers at +ε, to ≤ 6.4·10⁻¹⁵.
   - *The instance.* max|μ(ε) − μ(−ε)|/max μ ≤ 2.7·10⁻¹², and ≤ 3.4·10⁻¹⁰ per mode. The first-order deviation in the md's
     normalisation is 9.448·10⁻² … 1.423·10⁻⁴ (σ₀) and 9.838·10⁻² … 1.481·10⁻⁴ (∞): the md's values to all printed digits.
   - *My own cases* (absolute differences): 2.5·10⁻¹⁴ (n = 2), 3.1·10⁻¹⁴ (Λ symmetric), 2.5·10⁻¹⁴ (finite η ≠ σ), 1.8·10⁻¹⁴ (complex
     pairs, rates and complex eigenvalues). Their first-order deviations at ε = 0.1 are 0.23, 0.27, 0.017, 0.013. The control without
     the symmetry gives 3.6·10⁻³, 1.8·10⁻³, 9.2·10⁻⁴ (O(ε)), so the test can fail.
   - *The mismatch.* The md's column reproduces to all printed digits: 1.077·10⁻² … 3.121·10⁻⁴, and 1.129·10⁻² … 3.248·10⁻⁴ at η = ∞.
     The ratios per doubling are 2.08, 2.06, 2.02, 2.00, 2.00. Per mode, (μ(ε) − μ(−ε))/(2ε) → ℓ_j⟨aKᵀK⟩r_j: at k_e = 400,
     +1.207·10⁻² against +1.206·10⁻², −1.212·10⁻² against −1.213·10⁻², +7.956·10⁻⁴ against +7.960·10⁻⁴. My case F gives 0.18,
     0.070, 0.033.
   - *c₁.* At k_e = 800, (μ₁ − λ₁)/ε² = −2.32786·10⁻², against −2.32785·10⁻² from the perturbation formula with the exact Ā₁ and Ā₂.
     c₁/λ₁ = −108.16 (σ₀) and −109.89 (∞); the md has −108.2. The residual falls by 3.85, then 4.0: there is no ε³ term.
7. **(v)** (`check_noise.py`, `check_noise_hi.py`, `check_noise_conv.py`, `check_noise2.py`).
   - *S.* S from the left slow subspace satisfies S − εṠ = εKᵀ − ε²(SPKᵀ + AS) to ≤ 8.4·10⁻¹⁴ (relative, spectral derivative).
     |S − εKᵀ − ε²K̇ᵀ|/ε³ → 166, bounded.
   - *The md's noise table*, to all printed digits:
     - Cov(v)2ε³ − 1 = −9.6004·10⁻² … −1.8060·10⁻⁴.
     - ‖E[vY′ᵀ] + bᵀ‖/‖b‖ = 0.15426 … 4.8248·10⁻³.
     - ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ = 0.35808 … 1.1093·10⁻³.
     - The floor at t = 0: 3.23656·10⁻¹ against 4.28924·10⁻¹ (0.7546) … 1.33969·10⁻² against 1.34039·10⁻² (0.9995).
     - |b − b₀|/|b₀| = 0.6271 … 0.02957 (the md's 0.6271 … 0.0295; my maximum is over 4096 points, not 512).
   - *Resolution.* The Richardson step change is ≤ 1.6·10⁻⁵ (k_e = 12.5) … 2.6·10⁻⁹ (1600); the propagated covariance is periodic to
     ≤ 5.4·10⁻¹⁴. The resolution study at k_e = 200–1600 (N_c, 2N_c, 4N_c) resolves the max-norm deviation to ≤ 0.1%.
   - *The period-averaged floor of D‴.30(ii):*
     - η = σ₀, k_e = 12.5–1600: 0.75870, 0.90440, 0.97096, 0.99228, 0.99804, 0.99952, 0.99990, 0.99999;
     - η = ∞, k_e = 12.5–400: 0.76642, 0.90727, 0.97182, 0.99251, 0.99810, 0.99954;
     - η = 0 with K = Ĵ, Λ = σ₀I, k_e = 25, 100, 400: 0.90327, 0.99219, 0.99952.
   - *b's first-order term.* b₁ = −(K̈ + KKᵀK + KΛᵀ + aK̇)ᵀ: |b − b₀ − εb₁|/ε² = 2630, 3272, 3337, bounded.
   - *The cross term.* At t = 0, Cov(Y) − Cov(Y′)'s cross term is 3.556·10⁻⁵, against 3.565·10⁻⁵ from −ε²(d(KᵀK)/dt + 2aKᵀK) at
     k_e = 400.
   - *The net remainder.* Cov(Y) − Cov(Y′) − (ε/2)KᵀK is 1.779·10⁻⁵ against −(ε²/2)(d(KᵀK)/dt + 2aKᵀK)'s 1.783·10⁻⁵, with an
     entrywise difference of O(ε³) (`check_jitter_next.py`).
8. **The tests** (`check_tests.py`): the asserted quantities and the mutants of D3s-6.
9. **Bibliography** (web search).
   - Moan–Niesen: FoCM 8(3):291–301, 2008, DOI 10.1007/s10208-007-9010-0. Their Theorem 3 gives convergence for ∫₀ᵗ‖A‖₂ < π, the
     constant sharp.
   - Chang: SIAM J. Math. Anal. 3(3), 1972, DOI 10.1137/0503050; the record exists, the pages were not seen.
   - KKO's chapters 4–5 and Appendix A: titles from library records.
   - Haddad–Kokotović, IEEE TAC 22(5):815–821, 1977: from reference lists.
   - Yang–Zhu: CDC 2010 and 2012.
   - Berglund–Gentz: JDE 191(1):1–54, 2003, DOI 10.1016/S0022-0396(03)00020-2, and the 2006 book.

## Not checked

- **The constants.** ε₀ and the remainders were not computed; the lemma does not claim them.
- **Part D with the committed code.** `--add-D` was not rerun. Part D was recomputed by independent code, every printed number
  reproduces, and the md was regenerated from the raw record.
- **The other tests.** The four Corollary D‴.26 tests of `tests/test_adaptive_floquet.py` were not run (`-k` selected the three lemma
  tests).
- **The weakest mode at high k_e.** Its covariance entry is not resolved by my solver at k_e = 800 and 1600. Only its bound (below the
  resolved (0, 1) entry) is used in D3s-1.
- **The floor at high k_e.** At k_e = 800 and 1600 it is quoted to five decimals only. Its first-order coefficient on this instance
  was not resolved: the deficit falls faster than ε² between k_e = 400 and 1600. That suggests a first-order term of the opposite
  sign, as in D3q-1, but it is not established.
- **The noise part at n > 1.** The noise part was checked only on the instance (n = 1). The parity and the truncations were also
  checked at n = 2 (own case A); the noise was not.
- **The §8 "Open: the evenness of the noise floor without the slope".** Not examined.
- **The fast multipliers.** Their bound e^{−T₀/(2ε)} lies below double precision, so it was checked analytically only.
- **The real nonlinear loop.** The lifts of Corollary D‴.26 and the real loop of D‴.30 are outside the lemma and were not simulated.
- **The precedents of D3s-10.** Their content is from memory and abstracts; no full text was read (KKO Ch. 4–5, Haddad–Kokotović 1977,
  Yang–Zhu 2012).
- **The concurrent session.** Its commits a6ed653 and e325236, and its uncommitted edits to the two modules, are outside the material.
  If those edits change `chang`, `slow_generator`, `noise_decoupling` or `decoupled_covariance`, this report does not cover them.
