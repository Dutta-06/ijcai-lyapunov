# Referee report D3eq, round 148: the valley condition on the annulus by hand (research note)

## Scope
I reviewed review/RESEARCH_note_annulus.md and review/research_annulus/ (check_cvh.py, check_algebra.py and their outputs) at commit 29eb4bc. I worked on an export at /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref148/repo, with the venv symlinked. The repository was not modified.

Sources read:
- Paper 18 §4: Lemmas 4.3 and 4.4, Definition 4.5, Lemma 4.6 (i)/(ii), Theorem 2 and its proof, Remark rem:tiers, the frame corollary, Open problems (1), Appendix A headings.
- D61 (iii)(b), including forms (1)/(2), F_η, C = Λ⁺, M and H; its (iv); Appendix A.
- The V/VC/VA history in D61's text and the REVIEW_LOG.
- The downstream files: D58, D62, RESEARCH_note_c2, HANDOFF and PAPERS.

## Summary
**Verdict: the proof stands. I found no catch.** Theorem A, with Lemmas 1–3 and the geometry of Π_η, is correct and complete. Together with Paper 18's Lemmas 4.3–4.6 as they stand, it gives Ω′ = {|η| < 2} by hand. So Paper 18's Theorem 2 and D‴.61 (iii)(b) no longer need certificates D, V or VA.

The whole argument rests on one fact: beyond the closed unit disc Paper 18's bound is exactly F_η = ℓ(η − ν) + M(Re ν) − ½. This bound
- does not change along a vertical line (M is fixed), where ℓ is monotone in |Im τ| while |τ| < 1, and
- falls along the horizontal line through η, because ℓ′ = a + 1/a > 2 ≥ M′.

Earlier rounds worked in the saddle coordinate z and never wrote the bound in this form.

Findings, all minor or editorial:
- **D3eq-1 (minor):** Paper 18's sentence after Definition 4.5 ("all the paths used have this form … with (|η|/2)(1 + x_b) < 1") does not fit Π_η. Its first piece ends on the circle, where (|η|/2)(1 + x_b) = 1. This is harmless, because Definition 4.5 asks only that the sub-chord lie in the open disc, but the sentence must be rewritten.
- **D3eq-2 (minor/editorial):** The note's "Explicit margins" paragraph gives a gap on Π_η itself. Lemma 4.6 (ii) as written uses the scaled path (η/η₀)Π_{η₀}, not Π_η. The lemma's own compactness argument still covers it, so nothing is missing. The paragraph should say which path it bounds. "Negative on |η| < 2" should read 1 < |η| < 2, because m is undefined for |η| ≤ 1.
- **D3eq-3 (editorial, an improvement):** At θ = 0 and θ = π/2, Π_η is the straight path. There m(η) equals Paper 18's two axis gaps exactly. On the real axis m = log 4 + log(R − 1) − 2 log R + R²/4 − R + 1 = Γ_ℝ. On the imaginary axis the R²/4 − R + 1 becomes −R²/4 + R − 1. Theorem A's bound log(1 − x) < −x therefore replaces the third-derivative and concavity arguments in Theorem 2's "By hand" paragraph. The axes become a special case and no longer need separate treatment.
- **D3eq-4 (minor):** The consequences list in the note's §5 is incomplete. The full list is in "Consequences" below.
- **D3eq-5 (informational):** Generated outputs and code strings that say "computer-assisted" will be stale until they are regenerated: results/chaos/gh_peano.md, results/chaos/width_mismatch.md, and the printed lines in experiments/exp_gh_peano.py and exp_width_mismatch.py.

## What I checked and how
1. **By hand.** I re-derived every step: Lemma 1 (including the M branch change at |Re ν| = 1), Lemmas 2 and 3, the geometry of V and H, m(η) and its sign, the symmetries, and each requirement of Definition 4.5, Lemma 4.6 (ii) and D61 form (2). I also re-read Lemma 4.4's proof to confirm its bound holds on every compact set, including far outside the disc. Details are under "Verdicts per step".

2. **Independent numerics (my own code, mpmath at 60 digits, F_η taken straight from Paper 18's definition with Λ).** Scripts are in /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref148/mine/.
   - **path_check.py:** 506 values of η with ε = 10⁻¹², 10⁻⁹, 10⁻⁶, 10⁻⁴, 10⁻³, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.99 and 0.999 (so R from 1.001 to 2 − 10⁻¹²). The 32 angles per radius include 0, ε/10, ε, √ε/10, √ε, 3√ε, 23 interior angles, π/2 − 10⁻⁹, π/2 − ε and π/2.
     - Each segment was sampled at 300 uniform points plus geometric clustering at both ends.
     - Result: 0 failures off the sub-chord. On V ∪ H the maximum equals m(η) to relative 10⁻⁴⁰ in every case, attained at the start of V.
     - Example values of m: −2.5·10⁻³⁷ at ε = 10⁻¹², θ = 0; −5.0·10⁻²⁵ at ε = 10⁻¹², θ = π/2; −3.139·10⁻⁸ at ε = 0.005, θ = 0. These agree with the note.
   - **control.py, other quadrants:** 444 values of η over all four quadrants, with no symmetry used and the path built from ν₁ = Re η/|η| + i Im η. All points have |ν| ≥ 1, and the maximum is at most m with 0 failures.
   - **control.py, straight-path control:** on [0, η] beyond the circle the maximum of E is
     - positive: +4.3·10⁻³ (1.95, θ = 0.3), +4.1·10⁻⁵ (1.95, 0.097), +2.5·10⁻³ (1.98, 0.2), +4.8·10⁻⁴ (1.995, 0.1), +4.3·10⁻⁴ (1.9, 0.3);
     - negative: −3.6·10⁻⁵ (1.95, 0.05), −5.6·10⁻⁴ (1.95, 0.7), −1.1·10⁻⁴ (1.8927, 0.3).
     
     This matches D61's failing sector (0.087 … 0.615 at |η| = 1.95) and D's radius 1.8927.
   - **finiteK.py (a sanity check of the input Lemma 4.4 on the new region; not needed for the proof):** I computed c_K(√Kν) exactly from Gauss–Hermite nodes refined by Newton's method at 450 digits, for K = 40, 80 and 120, at points of V and H for η = 1.95e^{0.3i}, 1.9e^{1.3i}, 1.98e^{0.05i} and 1.5e^{0.8i}.
     - Inside the disc the values match 1 − ν².
     - Everywhere (1/K)log|c_K| stays below Λ(ν). For example at |ν| = 1.894 it is 0.566 against Λ = 0.655, and at |ν| = 1.491 it is −0.080 against Λ = −0.026 (M's linear branch).
     - The finite-K integrand rate relative to η/2 tends to the asymptotic F − F₀ with O(1/K)-type drift.
     
     Output: /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref148/mine/finiteK.out.

3. **The note's scripts.**
   - check_cvh.py: run with OMP_NUM_THREADS=1 and nice, 15.5 s. Its output is byte-identical to the committed check_cvh.out (all 5 checks passed; 197 values of η covered rigorously).
   - check_algebra.py: run with sympy 1.13.3 and mpmath 1.3.0 installed into W/pylib, not into the venv. Its output is identical to the committed one: m simplifies to 0 difference, the derivatives check, and the blow-up charts agree to ≤ 2·10⁻³³.

## Verdicts per step
- **Lemma 1: correct.** For |ν| > 1, log|ν| + ½Re ν² + Λ⁺ = M(Re ν) − ½ by definition. On |ν| = 1, Λ⁺ = M(a) − a² = 0 because |a| ≤ 1, so Λ is continuous and the identity extends to the circle. The formula holds for Re ν > 1 as well: M's branch enters only through M(Re ν), and M is C¹ with |M′| ≤ 2. D61's C = H(|σ|) − (|Re σ| − 1)₊² is the same function, which I checked by hand. Inside the disc, F = G − H(|ν|) also checks.
- **Lemma 2: correct.** With a, b and the derivative as defined in the note, ∂G/∂b = b(1 − a² − b²)/(a² + b²). On V, a = (R − 1)cos θ and b ∈ [0, (R − 1)sin θ], so a² + b² ≤ (R − 1)² < 1. G strictly decreases towards Im ν = Im η, and at θ = π/2 it tends to −∞ at η.
- **Lemma 3: correct.** dG/dp = −(a + 1/a) + M′(p) = −(a − 1)²/a − (2 − M′(p)) < 0 for 0 < a < 1. This holds for every p, since |M′| ≤ 2 everywhere; at p = 1, M′ = 2 from both sides.
- **Geometry: correct.** On V, |ν|² = cos²θ + r² > 1 for r > sin θ. On H, |ν|² ≥ cos²θ + R²sin²θ > 1 for θ > 0. At θ = 0, H is [1, R] on the real axis. Π_η meets the circle only at η/|η|. The symmetries F_{η̄}(ν̄) = F_{−η}(−ν) = F_η(ν) hold (M is even; |ν| and Re ν² are invariant), ν₁ is equivariant, and the direct check over all four quadrants confirms this.
- **m(η) and its sign: correct.** At s = 1/R, 4s(1 − s) = 4(R − 1)/R² = 1 − ε²/R², and the quadratic term is (ε²/4)cos 2θ. Since 0 < ε²/R² < 1 for R > 1, log(1 − x) < −x gives m < −ε²(4 − R²cos 2θ)/(4R²) < 0, including at cos 2θ = ±1. The ends of the range behave as they should: m → −∞ as R ↓ 1 and m → 0⁻ as R ↑ 2. Theorem 2 needs only compact subsets of the open disc.
- **Definition 4.5: met.** The path is piecewise C¹ from 0 to η. The chord with x₀ = min(½, (2/R − 1)/2) lies in the open disc. F < F₀ everywhere except at η/2, including F = −∞ at 0 and at η. The condition is about values along the path only. Laplace's method uses only the chord's line, and Re(φ″u²) = (2R²cos 2θ − 8)/R² < 0. No transversal or neighbourhood condition enters. The integrand is entire, so the path integral is well defined.
- **D61 form (2): met.** From b = (1 + t)η/2 with any 0 < t < 2/R − 1, the path runs along Ch, where g_η decreases, then along V and H.
- **Lemma 4.6 (ii): applies.** Lemma 4.4 is stated and proved for every compact subset of ℂ: the contour bound is uniform on compact subsets of ℂ∖{0}, and Lemma 4.3 handles a neighbourhood of 0. The scaled paths (η/η₀)Π_{η₀} lie in a fixed compact set, F_η is jointly continuous off {0, η}, and Λ is continuous across the circle. Nothing in the lemma depends on VA's or V's specific path.
- **Section 2 (blow-up charts):** consistent according to check_algebra. It is not part of the proof and I did not re-derive it independently.

## Consequences for the downstream statements
Everything below rests on Ω′ only. None of it used VA's specific path.

| Where | Current text | Change to |
|---|---|---|
| Paper 18, keywords (l.51) | include "computer-assisted proof, ball arithmetic" | drop or keep as optional |
| Paper 18, abstract (l.56) | "with a computer-assisted proof in ball arithmetic elsewhere" | by hand on the whole disc |
| Paper 18, introduction (l.97) | "One ingredient is computer-assisted …" | delete or rewrite |
| Paper 18, after Theorem 2 (l.247) | by hand on the unit disc and the axes, certificates on the rest | by hand on the whole disc |
| Paper 18, frame corollary (l.254) and its proof (l.557) | "R⋆ < 1 if t(v) < 0.0966" | no longer needed |
| Paper 18, Remark rem:tiers (l.271–273) | the whole remark | obsolete |
| Paper 18, Remark named rules (l.283) | "computer-assisted for t_hs ≥ 0.0966" | by hand |
| Paper 18, proof of Theorem 2 (l.475–491) | "By hand" paragraph, axes argued separately; "Computer-assisted" paragraph | "By hand" uses Theorem A, axes as θ = 0, π/2 (D3eq-3); "Computer-assisted" goes |
| Paper 18, sentence after Definition 4.5 (l.440) | see D3eq-1 | include Π_η |
| Paper 18, Open problem (1) | analytic proof of the last annulus | solved |
| Paper 18, notation table (l.222), Appendix A, app:cmd, Code and data availability | the verifier as "certificate of record" | optional cross-check or removed |
| Paper 18, supplement S3 and S5, manuscript/paper18/README.md | — | supplement S3 and S5 become informational; README updated to match |
| D61 | (iii) heading and lines ≈48–58; (iv) lines ≈86–99 (the r < 0.9900 / r < 1 tiers); Status ≈211–248 and 303–316; Appendix A preamble | Ω′ = disc by hand; D, V, VA and VC become history and cross-checks. Lemma VC and the four regimes stay true statements about VA's path. |
| D58 | lines 13, 26, 29, 174, 190, 201, 257, 273, 340, 365, 369 ("computer-assisted for t ≥ 0.0966" and similar) | by hand |
| D62 | lines 162–177 (including "D‴.61 (iv)'s use of the whole disc, computer-assisted beyond r = 0.0966, is not avoided by this route"), 224, 242, 250 | by hand |
| RESEARCH_note_c2 | Proposition A (iii) (l.16, l.166: "With certificates D and VA …"), Remark F (l.259, 284), l.360 | by hand; the proof repeats Lemma 4.6 (ii) verbatim and uses Ω′ only |
| HANDOFF.md | lines 212, 223, 276, 299, 314, 324, 350 | by hand |
| PAPERS.md | lines 425, 432, 443, 472 (the "Missing" item), 476, 479 | by hand |
| tests/test_gh_peano.py, tools/verify_certificates.py, results/certificates/README.md | describe the certificates | can stay as cross-checks; their wording should not claim the theorem depends on them |

## Issues
| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3eq-1 | minor | Paper 18 l.440: "All the paths used have this form … whose first piece is the real segment z ∈ [−1, x_b] with (|η|/2)(1 + x_b) < 1". Π_η's first piece ends at the circle, where (|η|/2)(1 + x_b) = 1. Definition 4.5 needs only the sub-chord in the open disc, so this is harmless but the sentence is wrong for Π_η. | Rewrite: "the straight path for \|η\| ≤ 1 and Π_η for 1 < \|η\| < 2, whose first piece is the chord to η/\|η\|". Paper 18 §4.4 and D61 (iii)(b). |
| D3eq-2 | minor | The note's "Explicit margins" paragraph bounds the gap on Π_η, but Lemma 4.6 (ii) uses (η/η₀)C₀. "Negative on \|η\| < 2" overstates: m is defined for 1 < \|η\| < 2 only. | State that the lemma's compactness argument applies unchanged to the scaled path, or restate the lemma for an η-dependent family with a continuous gap; restrict the range. In the note, and in Paper 18 if the paragraph is carried over. |
| D3eq-3 | editorial | m(η) at θ = 0 and θ = π/2 equals Paper 18's real-axis gap Γ_ℝ and its imaginary-axis gap exactly (checked by hand), and Π_η is the straight path there. | Replace the two axis paragraphs in Theorem 2's proof with Theorem A. Paper 18 §4.4 and D61 (iii)(b) proof. |
| D3eq-4 | minor | The note's §5 omits D58, D62, PAPERS.md, Remark rem:tiers, Remark named rules, the corollary's proof (l.557), Open problem (1), the notation table and the supplement. | Use the list in "Consequences" above. |
| D3eq-5 | info | Generated md outputs and printed strings in exp_gh_peano.py and exp_width_mismatch.py say "computer-assisted". | Regenerate after the text changes. |

## Sources
- Note and scripts: review/RESEARCH_note_annulus.md and review/research_annulus/check_cvh.py, check_cvh.out, check_algebra.py, check_algebra.out.
- Paper 18: manuscript/paper18/paper18.tex, l.229–275 and l.344–495; supplement.tex (S3, S5).
- review/drafts/D61.md: (iii) l.23–80, proof l.155–200, Appendix A l.402ff.
- Other files: review/REVIEW_LOG.md, review/RESEARCH_note_c2.md (l.10–26, 160–185, 255–262), review/drafts/D58.md, review/drafts/D62.md, HANDOFF.md, PAPERS.md.
- My scripts: /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref148/mine/path_check.py, control.py, finiteK.py and finiteK.out.
- My reruns of the note's scripts: /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref148/check_cvh.mine.out and check_algebra.mine.out.

## Lesson
About thirty rounds of certificates (D, V, VC, VA) were built in the saddle coordinate z, where Λ⁺ looks like a delicate correction to an ε³ cancellation. Written in the variables τ = η − ν and Re ν, the bound beyond the circle separates into ℓ(τ) + M(Re ν) − ½. A vertical-then-horizontal path makes it monotone, and the hard corner at η = 2 disappears. Before investing in a computer-assisted proof of an inequality along a chosen path, write the bound in the coordinates where it separates, and ask whether a different path avoids the hard region altogether.
