# Referee report D3ee, round 136: the Paper 18 manuscript, second review

Manuscript: `manuscript/paper18/paper18.tex` (985 lines) and `paper18.pdf` (35 pp.), commit f4825c8, branch claude/vigilant-clarke-hf939b. I exported it read-only to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref136/repo`.
- My own scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref136/own`: `exact.py`, `t1.py`, `t2.py`, `forms.py`, `plant.py`.
- The literature sources I read are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref136/lit`.

## Summary

The paper determines the exact asymptotics of the relative worst-case error e_K(v) of scaled Gauss–Hermite rules Q_{K,v} in the Gaussian RKHS. The main result is Theorem A:
- the rate is t(v) = tanh((r+|ln(v/v_hs)|)/2);
- off the cusp, the constant is C(v)^{1/2};
- at the cusp v_hs the law is (1−t²)^{1/4}(πK)^{1/4}t^K ₂F₁^{1/2};
- the classical rule satisfies e_K ~ 2^{−1/4}(1+λ)^{3/4}λ^K at every length-scale.

The revision makes six main changes:
1. **A uniform Laplace lemma (Lemma 4.1).** It is used for I_K, for the circle in the corner ratio, for the chord and for the angular integral.
2. **A rewritten path lemma (Lemma 4.6).** The chord is fixed along η/|η|, and the phase splits exactly as ψ_η = −α_ηz² + (log(1−z²)+z²).
3. **Corollary 3 with a Gamma-tail bound.** The constant is identified by algebra.
4. **Theorem 4's internal steps as Lemmas 4.7–4.9.**
5. **A rewritten §1.3, after reading KOG21, KS19 and GRR24.** It identifies KOG21's rule as the Fock scaling, states that its constants c₁ = ln(2/λ) and c₂ = ln(1/λ) bracket the exact ln((2−λ)/λ), says that KS19's "only if" is false, and adds a consequence through the GRR24 transference (the exact rate of Gauss–Hermite in the Hermite space of base β).
6. **Stored certificates and a standalone verifier** (`tools/verify_certificates.py`).

## Overall assessment (recommendation: minor revision)

The first review's major points are resolved in substance.
- **Mathematics.** I re-checked every new or rewritten argument step by step and found no error. The one exception is a sign typo in the proof of Lemma 4.3, and the proof is correct with the sign fixed.
- **Literature.** I read KOG21, KS19 and GRR24 myself. I confirmed the mapping of KOG21's rule to the Fock-scaled Q_{K,v} numerically, to rounding of the inputs, and KOG21's empirical rates 0.528 and 0.203 are t(v) = 9/17 and 25/123. KOG21's constants are as stated. KS19's statement is indeed contradicted by Theorem 4(b), which I confirmed with independent exact errors at 2αγ² = 1 and 4. The GRR24 consequence holds: GRR24 relate β and σ for integration through (g20), not (g2), which gives β = λ, and their Remark 4.4 maps Q_K exactly onto Q_{K,v} at ρ(v) = λ.
- **Verifier.** It runs (exit 0, 14 s), rejects a corruption I planted myself, and checks what it claims.

What remains is not mathematical:
- **(a) Three works are still cited secondhand.** KW12 and KSW17 are cited as KOG21 reports them. IKLP15 is marked "not read", although it is on arXiv (1403.5102; I read its §6.1).
- **(b) The appendix describes the wrong checker.** Appendix A and the code-availability statement describe the proposing program's certify step, not the verifier, which is now the natural proof of record. Several details differ: the tails in Arb rather than in doubles × (1+10⁻⁷), and a direct cover of the parameter domain in rationals rather than Lemma A.8 (v)–(vii). The verifier also labels lemmas with the old internal numbering ("Theorem 3", "VA1–VA9").
- **(c) Editorial residue.** 18 `% CHECK` and 6 `% NOT IN MD` comments remain, along with the placeholders. Some notation still clashes (β, ψ, ℓ, h, M, δ). The paper is long (35 pp.).

None of these affects correctness. A journal would accept after a minor revision that settles (a)–(c).

## Status of the first review's issues

| id | round-135 issue | status | how / what remains |
|---|---|---|---|
| I1 | literature unfinished; KOG21 overlap | **resolved in substance** | §1.3 and Remark 3.4 rewritten from KOG21 (v2), KS19 (v4) and GRR24 (v2), all read in full. I confirmed each mapping (see Verification). The novelty claim was correctly narrowed: the rate λ was observed numerically in KOG21's Fig. 1. **Remains:** KW12 and KSW17 are cited only as reported by KOG21 and KS19, IKLP15 is marked unread although it is on arXiv, and Al-Jarrah was not reached. |
| I2 | path lemma: branch, direction class, uniformity | **resolved** | Definition 4.5 fixes the chord along η/|η|. On the chord the integrand is exactly (η/2)^{2K}e^{Kη²/4}e^{Kψ_η}a_η dz with a_η(0) = α_η, so only the principal root of α_η enters (through Lemma 4.1). Openness under γ_η = (η/η₀)γ₀ is correct. |
| I3 | Corollary 3: tail beyond the disc; 5.9·10⁻²⁹⁸ | **resolved** | Γ(2K+1, X_K) and Γ(2K+2, X_K) with Chernoff are correct. The constant is now shown by algebra, which I re-did, and the numerical claim is gone. |
| I4 | certificates not independently checkable; margins | **resolved, with documentation gaps** | Stored data with SHA-256, a standalone verifier, and the margins reported. Gaps: Appendix A and the code statement describe the proposer, not the verifier; the lemma labels do not match; A.4 does not list the verifier command (M2). |
| I5 | structure | **largely resolved** | Proof environments are used, Lemmas 4.7–4.9 are stated, and the supplementary Appendix B is split off. The paper is still 35 pp.: move Appendix B and V to a separate supplement (m14). |
| I6 | notation overload | **largely resolved** | A notation table and systematic renaming are in. Remaining clashes are listed in m9. |
| I7 | non-uniformity in v | **resolved** | L83, L89, L241, L260. |
| I8 | ρ < −1 "verbatim" | **resolved** | L618 names the maximiser (±2s_a^{1/2}, 0), and Lemma 4.9 covers both signs. |
| I9 | Remark 3.1 "partly a sketch" | **resolved** | Now Open problem (2), with the proved facts stated after Theorem 1. |
| I10 | internal artefacts | **partly resolved** | The Verification paragraph and `\cite[A6]{repo}` are removed. Remaining: 18 `% CHECK` and 6 `% NOT IN MD` comments, and placeholders at L42–45, L711, L714 and L971. |
| I11 | undefined terms | **resolved** | Mercer-zero and the Mercer frame's rule are defined (Remark 3.2), and I checked both against KS19. "Direction class" is gone. |
| I12 | LaTeX | **resolved** | My compile: 35 pp., 0 overfull boxes, and only the hyperref bookmark warnings. |
| I13 | bibliography | **partly resolved** | KOG21, KS19 and GRR24 are complete. Remaining: Olver's chapter and ten "from memory" items (m12), and IKLP15 should get its arXiv number. |
| I14 | j = 1 cumulant; disc radius | **resolved** | Theorem 2 includes j = 1. The proof uses the circle \|η\| = ½, and K(Var S − ½) = −K(K−1)/(2(K+1)(2K+1)) is checked. |
| I15 | Theorem 1's novelty | **resolved** | §1.3, L128. |
| I16 | runtimes | **resolved** | Sourced, and my 15 s and 7 s are quoted. The verifier's 14 s should be added (m7). |

## Major comments

None that block acceptance. The following two points should be settled in the revision.

**M1. Literature that is still secondhand (L113, L300, L947, L961–963).**
- **Mapping.** KOG21's footnote 1 and §1 report KW12/KSW17's bound as b_n(α²/ℓ²)^n, with convergence only for ℓ > α. That is strong evidence that neither paper contains λ or convergence beyond the threshold. Still, claim (iv) and Remark 3.4 cite KW12 for its claimed lower bound, and KSW17 for page 830 and the form of b_K, without the authors having seen either. Before acceptance, obtain both (BIT 52 (2012) 425–436; Math. Comp. 86 (2017) 829–853) through a library and replace the CHECK comments with direct statements.
- **IKLP15 is on arXiv (1403.5102).** Its Proposition 1 (§6.1, p. 15) gives e²(A_n, K_{a,b,ω}) ≤ ω^{a(2n)^b}√(8π)/(1−ω²). With b = 1 and ω^a = β this is the rate β, via Cramér's bound. This confirms that claim (v), the exact rate β/(2−β), is not in IKLP15. Cite it directly and remove the CHECK at L947.
- **Minimal errors.** Optionally, mention in "Other work" (L125) KOG21's Theorem 2.7, the minimal-error bounds derived from KSW17's Theorem 4.1. GRR24 (§4.2) state that their Theorem 4.6 improves on it.

**M2. Make the standalone verifier the certificate of record, and describe it accurately (L714, L720–849, L867–884).**

The verifier is sound as far as I read and tested it (details under Verification). But the text does not yet match it:
- **(i) The code statement overstates.** It says the verifier "re-applies Lemmas A.2, A.4–A.11 and A.12". For Lemma A.8 (domains and cover) this is not so. The verifier does not use the root grids, clips or float-midpoint bisection of A.8 (v)–(vii). It proves directly, in rationals, that the stored boxes cover the closed domain of each piece. For A-FC that domain is {0 ≤ u ≤ 4, 0 ≤ q ≤ q_A, 0 ≤ ϱ ≤ δ₀, uϱ ≤ ε₀}, and similarly for A-J, B-J and B-F. This is a simpler and stronger statement. Restate Lemma A.8 as (i)–(iii), the analytic domain facts, plus "the stored boxes cover these closed domains", and drop (iv)–(vii) or mark them as propose-only.
- **(ii) Other passages describe the proposer.** Lemma A.5's last sentence (tails "in doubles, multiplied by 1+10⁻⁷") and trusted-base item (3) describe the proposer. The verifier computes the tails in Arb. If the verifier is the proof of record, update both, or say explicitly which program each paragraph describes.
- **(iii) The verifier transcribes the forms.** Its f_A, f_B, f_A0 and f_B0 are transcriptions of the proposer's expressions, with the same variable names and the same ordering. "Shares no code" is true at the level of imports, but the algebra is shared. That is acceptable, because Lemmas A.4 and A.6 carry the algebra and the verifier's own docstring says it does not check them. Say so in the text: "it re-implements the stated bounds; the algebra of Lemmas A.4 and A.6 is proved in the text". Also add the forms-against-definition check (from the test suite) to the verifier's output as an informational line. I ran this check independently: the forms equal E/X², with E from (A.1), at 4000 random path points.
- **(iv) The labels do not match.** The verifier's output, docstring and `results/certificates/README.md` speak of "Theorem 3", "Lemmas VA1–VA9" and "D61". The paper has Theorem 2 and Lemmas A.2, A.4–A.12. Add a mapping, or print the paper's numbers.
- **(v) A.4 does not mention the verifier.** "One command per certificate" (A.4) lists only the proposers. Put `venv/bin/python -I tools/verify_certificates.py --dir results/certificates` (14 s) first, as the check a reader should run.
- **(vi) Stale build notes.** `manuscript/paper18/README.md` still says the verifier "has not been written yet".

## Minor comments

1. **L381 (p. 14), proof of Lemma 4.3: sign typo.** It should read h(x) = log(1+iζx) − iζx **−** ζ²x²/2. With "+" the identity −αx² + h = log(1+iζx) − iζx − x²/2 fails, and h″(0) = 2ζ² ≠ 0, which violates Lemma 4.1's hypothesis.
2. **L516 (p. 18): Lemma 4.1's hypothesis is not met literally.** Corollary 3 applies Lemma 4.1 "with the O(K^{−1/2}) above", but the lemma assumes |b_{p,K}| ≤ M/K. Either state Lemma 4.1 with |b_{p,K}| ≤ ε_K and the error O(1/K + ε_K), as the proof gives verbatim, or add one sentence. Also say that the second peak at θ = π is handled by γ(θ+π) = γ(θ), and that the rest of the circle is exponentially smaller.
3. **L631 (p. 22): the deviation is relative.** "Within 6.1·10⁻⁵ of 1/t(v)" is a relative deviation, as the md's "limit over 1/t(v), minus 1" shows. At τ/σ = 2, v = 0.5v_hs, the absolute gap in Table 2 is 1.576110 − 1.576014 = 9.6·10⁻⁵. Write "relative deviation". The same applies to 3.2·10⁻⁷.
4. **L65 vs L414: ψ means two things.** ψ_v = dμ/dN(0,v) and ψ_η is the chord phase. Rename one of them.
5. **L102, L115, L132, L269: β means four things.** The Beta variable, KOG21's β (a standard deviation), KS19's β with β_M, and GRR24's base parameter β. Use β only for one of them, and for instance b_K or β_G for the others.
6. **L280 vs L96 vs L330: ℓ means three things.** Mehler's ℓ² = 2ρ/(1−ρ²), KOG21's length-scale, and the Dyck-path bump level.
7. **L867–884: add the verifier command and runtime** (see M2(v)).
8. **L125: make the minimal-error rate explicit.** Write GRR24's bound as C(λ/2)^{2K}(K+1)^{−2}; "rate (λ/2)² per point" is easy to misread.
9. **Further clashes.** M is the constant of Lemma 4.1, the function M(a) of Lemma 4.4, and the bound in Lemma 4.6. δ is the δ of Lemma 4.1, δ = R_⋆²/(4t) − 1 in Corollary 3, δ = |ν|−1 in A.4, and δ₀. h is h_p in Lemma 4.1, h(z) in Lemma 4.6, h_{2j}, and the half-length in A.9. γ is the path, γ_x and γ_y, and KW12's γ. Either localise these explicitly or rename.
10. **L188 / Proposition 2.2:** add "proved in §4.6".
11. **L290: wording.** "Divergence is excluded … by the elementary bound e_K² ≤ 4σ_*v/τ²" says only that the errors are bounded. Say "the errors stay bounded"; convergence is Theorem 4(b).
12. **Bibliography.** Resolve the remaining CHECKs: Olver, *Asymptotics and Special Functions*, Ch. 3 (Laplace's method), which is correct as far as I recall. From memory, the following details look right: Johansson, IEEE TC 66(8) (2017) 1281–1292; Jawecki–Auzinger–Koch, BIT 60 (2020) 157–197; and Carlson, de Boor, Flajolet, Gautschi–Varga, KMR86, PR29 and SGFSL10 as given. They still need checking against the sources.
13. **L292: tone.** "This disproves the statement of [KS19, §3]" is correct. It would be fairer to add that the statement is a one-line remark in KS19 that rests on KW12's lower bound, which KSW17 had already called erroneous.
14. **Length.** Move Appendix B, and V's cross-check paragraph (L853–865), to a separate supplement file, as the README already allows. Shorten §5 to Tables 2 and 4 and one sentence on the transform.
15. **A figure would help.** For instance t(v) against ln(v/v_hs), with the cusp, and e_K(v)^{1/K} at K = 40. There are no figures now.
16. **Tensor products.** A referee will ask about d > 1. One sentence would do: KOG21's tensor-product argument (their Theorem 2.10) combined with Theorem A gives the upper rate t(v)^{N^{1/d}} for the d-fold product of Q_{K,v}. Open problem (6) could then be narrowed to sharp constants.
17. **L42–45, L711, L714, L971: placeholders**, and the 24 internal comments, before submission.

## Verification I did

**[own] Mathematics, by hand, step by step.**
- *Lemma 4.1.* Taylor with |Kh_p| ≤ 2M on |x| ≤ K^{−1/3}; the odd part vanishes; the error terms are O(K^{−3/2}); and the principal root by analytic continuation. All correct.
- *Lemma 4.3.* He_K(z) = E(z+iZ)^K and the shift x ↦ x + iζ. The bound ½log|1+iζx|² ≤ −Im ζ·x + |ζ|²x²/2. This is where the h-sign typo is (m1). On the circle: f_ν(ν) = 1+ν²/2+log ν, f_ν″(ν) = (1−ν²)/ν², the identity (1−cos θ)(1−|ν|²cos(2arg ν+θ)), the clockwise orientation, and the assembly with Stirling to (1−ν²)(1+O(1/K)).
- *Lemma 4.4.* max_θ(2a cos θ − ½cos 2θ) = M(a)+½, and Λ⁺ = H(|ν|) on |ν| ≤ 1.
- *Lemma 4.6.* The chord identity (4.4); ψ_η = −α_ηz² + h with h independent of η; a_η(0) = α_η; the constant (2K)!/(K!(K−1)!)4^{−K}(π/K)^{1/2} → 1; and the openness argument.
- *Theorem 2 (by hand).* On the real axis Γ_ℝ and its derivatives vanish at 2 and Γ_ℝ‴ > 0. On the imaginary axis the gap is concave, with second derivative ≤ −1. The cumulant coefficients are −(2j−1)!/4^j.
- *Corollary 3.*
  - the tails, and 4t < f(t) < 4;
  - Chernoff, Γ(m,X) ≤ Γ(m)(X/m)^m e^{m−X}, with X_K = 2K(1+δ);
  - the radial Lipschitz step and the angular Laplace step;
  - 2γ_x(γ_y−γ_x) = (t−t_o)/(2t²t_o), (c_a|c_b|)^{−1/2} = (tt_o)^{−1/2}[(1−c_a)(1−c_b)]^{−1/2}, and the result C(v);
  - the cusp: (c_a|c_b|)^{−1/2}t_hs = (1−t_hs²)^{1/2}, and the ₂F₁ by Parseval.
- *Lemmas 4.7–4.9 and Theorem 4.*
  - Φ in y_± and its maximum 2log|s_b|;
  - R1 via Θ ≥ 2/k on each sector; R2 via convexity of Θ;
  - size-biasing to χ²_{4K+1}; the tilt N_t;
  - (K!/(2K)!)²(4K−1)!! → 2^{−1/2};
  - the ρ < −1 change of variables c̃ = c/(1−2c);
  - the tails at ρ = ±1, with 4λ < R_b² < 4;
  - (b′_{K+1}/b′_K)² = (16K²+16K+3)/(16K²+16K+4).

**[own] Exact errors** (`own/exact.py`: independent Arb code, numpy start values, Newton, Christoffel weights, closed-form Gaussian kernel sums).
- *Classical rule at 2αγ² = 1 and 4.* C(v) = 1.299038105677 and 1.707629936491. At K = 100, y_K − 1 = −1.41245·10⁻³ and −5.40949·10⁻⁴, matching Table 4. K(y−1) is −0.141245 and −0.054095, against ς₁ = −0.1412037 and −0.0540981.
- *KOG21's rule, computed in KOG21's own units, against our Fock-scaled e_K.* At (α̃,ℓ) = (3/2,1) and (1,7/5), K = 20–61, the ratio −1 is ≤ 2·10⁻¹⁴, which is the rounding of the double-precision inputs. The ratios e₆₁/e₆₀ are 0.528891 and 0.203209 (KOG21 report 0.528 and 0.203), with t(v) = 9/17 and 25/123.
- *Hermite space.* Σw_iw_j k_β(z_i,z_j) − 1, with Mehler's closed form for k_β, against e_K(v) at ρ(v) = β. At β = 1/3 and 0.6, K = 10 and 30, the ratio −1 is ≤ 10⁻¹⁴, and y_K → 1 at the rate β/(2−β).

**[own] Remark 3.1, in Arb.**
- f(0.0966)^{1/2} = 0.999825651, and f(t) = 1 at t = 0.0966496.
- f(0.8052) = 3.58230645 < 1.8927² = 3.58231329.
- f(0.99) = 3.979933 < 1.995² = 3.980025.
- The stored doubles: fl(1.8927) > 1.8927 and fl(1.89) < 1.89, so the seam holds for the real radii.

**[own] Verifier.**
- `(cd results/certificates && sha256sum -c SHA256SUMS)`: all four files OK.
- `venv/bin/python -I tools/verify_certificates.py --dir results/certificates` exits 0 in 14 s:
  - D: 402664 leaves, worst −1.9432·10⁻¹¹, tiled;
  - inner pieces: 54 / 2533 / 413 / 1189 boxes; A-FC worst −8.781·10⁻¹⁴; covers proved;
  - outer piece: 16083 segments, worst −4.228·10⁻⁹;
  - 1306 handoff sub-boxes; 220 η-boxes tiled; end piece −1.3314.
- Code read in full for:
  - `d_bound`, against Lemma A.2: the corners, the stable root formula, the undecided sign of q_a, M taken at the upper end;
  - `tiles`: exact dyadic areas and a sweep for overlaps;
  - `fA`, `fB`, `fA0`, `fB0` and `mean_value`, against Lemmas A.4, A.6 and A.7;
  - `cover`: a rational decomposition of the domain;
  - `E_mid` and `seg_bound`, against Lemma A.9, where I re-derived ∂_s, ∂_R, ∂_θ and ∇Λ;
  - `y_enclosed`, `chain` and `verify_meta`.
- The verifier imports nothing from `experiments/`. Its form code is a transcription of the proposer's (M2(iii)).
- *Forms against the definition* (`own/forms.py`): f_A on F, C and J, and f_B on F with the dropped term added back, against (A.1) at 4000 random path points: equal to ~10⁻¹².
- *Tests:* `pytest -k "round135 and (planted or lemmas)"`: 2 passed (27 s).
- *My own planted corruption* (`own/plant.py`): one B-F box shrunk by 10⁻³ of its p-width, then re-hashed. The verifier rejects it with "[FAIL] VA5 B-F … uncovered piece", NOT VERIFIED.

**[own] Build.** I compiled with pdflatex: 35 pp., 0 overfull boxes. I rendered pp. 1, 7, 23 and 27; the typesetting, Table 1 and Tables 2–3 are fine. Lemma numbering in the PDF: A.1–A.14, with Propositions A.3 and A.13.

**[md]** I grepped 58 manuscript decimals (Tables 2–4, §5's transform, Remark 3.1, the counts and margins of D, VA and V, Appendix B) in gh_peano.md, frameless.md and width_mismatch.md: all found. The context of 6.1·10⁻⁵ shows it is a relative deviation (m3).

**[read]**
- **KOG21**, arXiv:2004.12654 v2, LaTeX source plus Fig. 1. I checked: §1, p. 2 and footnote 1; §2.2's rule with β² = α²ℓ²/(α²+ℓ²); Lemmas 2.2–2.4; Theorem 2.5 and eq. (2.9); Remark 2.6; Remark 2.12; §4; and Fig. 1 with (α,ℓ) = (1.5,1) and (1,1.4) and r ≈ 0.528 and 0.203.
- **KS19**, arXiv:1803.09532 v4, source plus PDF. I checked: eq. (2.5); Theorem 2.2 and eq. (2.10), x̃ = x^GH/β with β⁴ = 1+8ε² at α = 1/√2; §3, p. 11's "if and only if ε² < 1/2"; Theorems 3.1–3.2; and the end of §3 on the generalised Gaussian rule with 2N eigenfunctions. That rule is Q_{K,v} at v = 2τ²/(β_M²+1), ρ(v) = q, as the manuscript says.
- **GRR24**, arXiv:2412.05368 v2 (the latest version), source. I read §2.5 (the isometry (g2), Theorem 2.? from GHHR21) and §4 in full: Lemmas 4.1–4.2, Theorem 4.3 with (g20) 1−β = 1/(1+2σ²), Remark 4.4, Corollary 4.5, Theorem 4.6, and §4.4's IKLP15 Prop. 1 Cβ^n.
  - With σ_G² = τ²/4, (g20) gives β = λ.
  - Remark 4.4's map sends Q_K to nodes z_k/e with e² = (1+4σ²)/(1+2σ²), i.e. v = τ²(2+τ²)/(2σ_*²), and to weights w_k·dμ/dN(0,v). That is exactly the manuscript's Q_{K,v} at ρ(v) = λ.
  - The prefactor (1+4σ²)^{−1/4} is the initial error.
- **IKLP15**, arXiv:1403.5102: §6.1 only (Proposition 1).

**[snippet]** An arXiv site search for later work on scaled or classical Gauss–Hermite in Gaussian or Hermite spaces found nothing beyond KOG21 and KS19 (export.arxiv.org was blocked; the arxiv.org search worked).

**[memory]** The bibliographic details in m12, and Olver Ch. 3.

## Issues

| id | sev. | issue | proposed fix |
|---|---|---|---|
| J1 | minor (should fix) | KW12 and KSW17 are still cited secondhand for specific content (p. 830, b_K, the claimed lower bound) | Obtain both through a library and replace the CHECKs at L113 and L300 with direct statements |
| J2 | minor | IKLP15 is marked "not read" but is on arXiv (1403.5102). Its Prop. 1 is rate β (e² ≤ ω^{a(2n)^b}√(8π)/(1−ω²)), which confirms claim (v) | Cite it directly with the arXiv number and drop the CHECK at L947; optionally mention KOG21's Theorem 2.7 at L125 |
| J3 | minor (should fix) | The code statement says the verifier re-applies A.8, but it proves the cover directly; A.5 and the trusted base describe the proposer's double-precision tails | Make the verifier the proof of record and restate A.8 as "domain facts (i)–(iii) + stored boxes cover the closed domains"; update A.5 and trusted-base item (3) (L714, L730, L783, L798–807) |
| J4 | minor | The verifier's form code is a transcription of the proposer's; the algebra is shared and not checked by the verifier | Say so in Appendix A; print the forms-against-(A.1) sample check as an informational line |
| J5 | minor | Label mismatch: the verifier and its README use "Theorem 3", "VA1–VA9" and "D61"; the paper uses Theorem 2 and A.2, A.4–A.12; the manuscript README says the verifier is "not written yet" | Align the labels or add a mapping table; fix the READMEs |
| J6 | minor | A.4 lists only the proposers' commands | Add the verifier command (14 s) first |
| J7 | typo | Proof of Lemma 4.3 (L381): h(x) = … **+**ζ²x²/2 should be **−** | Fix the sign |
| J8 | minor | Lemma 4.1 is applied with an O(K^{−1/2}) perturbation in Corollary 3 (L516); the hypothesis is \|b\| ≤ M/K | State Lemma 4.1 with \|b\| ≤ ε_K and the error O(1/K+ε_K); add the θ = π peak |
| J9 | minor | "Within 6.1·10⁻⁵ of 1/t(v)" (and 3.2·10⁻⁷) are relative deviations (L631) | Write "relative" |
| J10 | minor | Notation clashes: ψ, β, ℓ, h, M, δ, γ (L65/414, 102/132/115/269, 96/280/330, …) | Rename, or localise explicitly |
| J11 | editorial | 18 `% CHECK` and 6 `% NOT IN MD` comments, and placeholders (L42–45, L711, L714, L971) | Resolve and remove; archive the code with a DOI |
| J12 | editorial | 35 pp.; no figure | Split off Appendix B and V as a supplement, trim §5, add a t(v) figure; add one sentence on d > 1 via KOG21's Theorem 2.10 |

## Lesson

A standalone verifier is only as independent as its formulas. Re-typing the proposer's expressions removes the shared-code risk but not the shared-algebra risk. The algebra has to be carried by the paper, the forms checked against the definition, and the document has to describe the checker that is actually of record, not the one it replaced. A second lesson concerns "not available": IKLP15 was marked unreadable, yet it has been on arXiv since 2014. A search by title or author before giving up would have settled claim (v) in the first pass.
