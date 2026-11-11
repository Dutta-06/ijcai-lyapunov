# Referee report D3ew, round 154: the Paper 18 manuscript after the by-hand annulus

## Scope
I reviewed `manuscript/paper18/paper18.tex` (28 pp.), `supplement.tex` (12 pp.), both committed PDFs and `manuscript/paper18/README.md` at commit b83e915. I worked on an export in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref154/repo` with the venv symlinked; the repository was not modified.

Background read:
- `review/RESEARCH_note_annulus.md` and its review, round 148 (D3eq);
- `review/RESEARCH_note_c2.md` and its review, round 139 (D3eh);
- the round-135 and round-136 reports, and their issue tables I1–I16 and J1–J12;
- `git diff 638fb1b b83e915` on `manuscript/paper18/`, i.e. everything changed since round 136.

Venue assumed: SINUM, Constr. Approx. or FoCM. Placeholders and `% CHECK` items are out of scope.

## The new proof (verdict per step)
All of these I checked by hand, and the claims about F_η also numerically (see Verification).

- **Lemma 4.7 (i), F_η = G_η on |ν| ≥ 1: correct.** For |ν| > 1 it follows from the definition of Λ⁺. On |ν| = 1, Λ⁺ = M(Re ν) − (Re ν)² = 0, so Λ is continuous across the circle.
- **Lemma 4.7 (ii): correct.** ∂G/∂b = b/(a²+b²) − b = b(1−a²−b²)/(a²+b²), so G strictly decreases towards Im ν = Im η inside {0 < |η−ν| < 1}.
- **Lemma 4.7 (iii): correct.** On b = 0, ∂G/∂Re ν = −(a+1/a) + M′(Re ν) = −(1−a)²/a − (2−M′) < 0. The hypothesis Re ν ≥ 0 is not needed, since |M′| ≤ 2 everywhere; it is harmless.
- **Proposition 4.8, first segment: correct.** On ν = xη with x ≤ 1/R, Λ = 0. g_η is unimodal at ½, because g′ = (1−2x)(1/(x(1−x)) − Re η²) and Re η² < 4.
  - With d = 1/R − ½ one gets 4x(1−x) = 1 − 4d², and the quadratic term is d²R²cos 2θ.
  - So m(η) = log(1 − (2−R)²/R²) + ((2−R)²/4)cos 2θ, as printed.
  - The bound via log(1−y) < −y and its sign (R²cos 2θ < 4) are right.
- **Vertical segment: correct.** a = (R−1)cos θ and b ∈ [0, (R−1)sin θ], so a² + b² ≤ (R−1)² < 1 and |ν| ≥ 1. Lemma 4.7 (ii) applies, and the endpoints θ = 0 and θ = π/2 are handled.
- **Horizontal segment: correct.** |ν|² ≥ cos²θ + R²sin²θ ≥ 1, b = 0 and a ∈ [0, (R−1)cos θ]; Lemma 4.7 (iii) applies. The maximum over the vertical and horizontal segments is attained only at η/R.
- **Chord: correct.** x₀ = min(½, (2/R−1)/2) < 2/R − 1, so the chord lies on the first segment and in the open unit disc.
- **Quadrant reduction: correct.** F_{η̄}(ν̄) = F_{−η}(−ν) = F_η(ν) because M is even, and ν₁ = Re η/R + i Im η is equivariant.
- **Definition 4.5 and the sentence after it: correct.** D3eq-1 is fixed: "[0, η] for |η| ≤ 1, Π_η for 1 < |η| < 2, first piece [0, η/|η|]".
- **Lemma 4.6 (i): correct.** Pointwise strict and sup-strict are equivalent here, because g → −∞ as x → 1 and the value at x = 1/R is negative.
- **Lemma 4.6 (ii), including its use of the scaled path (η/η₀)C₀: correct for η₀ ≠ 0.**
  - (η, ν₀) ↦ F_η((η/η₀)ν₀) − F_η(η/2) is jointly upper semicontinuous into [−∞, ∞). On the compact set C₀ minus the open sub-chord it is negative, which gives ε₀ uniformly.
  - On the sub-chord the integrand is exactly (η/2)^{2K}e^{Kη²/4}e^{Kχ_η}a_η dz, with a_η(0) = α_η. Lemma 4.1's hypotheses hold uniformly, and the constant (2K)!/(K!(K−1)!)·4^{−K}(π/K)^{1/2} → 1 is right.
  - The scaled path is not Π_η, and Proposition 4.8's monotonicity is neither needed nor claimed for it; the compactness argument suffices (D3eq-2 is consistent).
  - **Gap at η₀ = 0: see D3ew-1.**
- **Proof of Theorem 2.**
  - "The disc": correct, up to D3ew-1.
  - "The axes" as θ = 0 and π/2: Π_η is the straight path and m(η) is the domination gap; the formula log 4 + log(R−1) − 2 log R ± (R²/4 − R + 1) and both bounds are correct.
  - "Cumulants": unchanged; it needs uniformity on |η| ≤ ½, including η = 0 (D3ew-1).
- **Corollary 3 with R⋆ < 2 only: correct.** 4t ≤ f(t) = 4t ln(1/t)/(1−t) < 4 gives Δ⋆ > 0, and the Γ-tail Chernoff bound is unchanged. The 0.0966 clause is gone with nothing depending on it.
- **Theorem 4 (d): correct.** Only the axes are used, and they are the cases θ = 0 and π/2.

## Consistency of the "by hand" upgrade
| Place | Status |
|---|---|
| Keywords and MSC | "computer-assisted proof, ball arithmetic" and 65G30 removed ✓ |
| Abstract | by hand on the whole disc ✓ (wording: D3ew-2) |
| Introduction | ✓ (history wording: D3ew-3) |
| Statement of Theorem 2 and the sentence after it | ✓ |
| Corollary 3 and its proof | 0.0966 gone ✓ |
| Remark on named rules | ✓ |
| Old Remark rem:tiers | removed ✓ |
| Open problems | old (1) removed and renumbered ✓ |
| Notation table | ✓ |
| Code statement | ✓ |
| Supplement intro, Appendix A preamble, S3, S5 | all say "cross-check, which no proof uses" ✓ |

Nothing in the paper still depends on D, V or VA. In particular, nothing uses the domination condition off the axes:
- it is defined (Definition 4.5);
- it is void on |η| ≤ 1;
- it is used on the axes only, where Proposition 4.8 gives it by hand with gap m(η).

Corollary 3 uses Ω′ = disc via Π_η, not via D. The c₂ note's Proposition A (iii) uses Ω′ only. So D is not needed anywhere, and nothing needs to be stated about it.

Supplement S3 correctly records that the domination condition fails on a sector at |η| = 1.95, which is why the deformed path is needed.

Remaining wording issues: D3ew-2 (the "monotone exponential rate") and D3ew-3 (revision-history text in the PDFs, a stale "(by hand)" at l.631, the code statement's dangling "and V on …").

## c₂
Every c₂ = ς₁/2 statement agrees with RESEARCH_note_c2.md and its round-139 review:
- the notation table's ς₁ row: proved for |ρ| ≥ 1, formal at a frame;
- §5 and Table 3: "proved for these rules without a frame";
- open problem (3): needs (F) with an F₂ term and an O(K⁻³) remainder; F₂'s closed form is numerical only;
- supplement §S1 and §S4: the D3eh-3 correction, that (F) with O(K⁻²) gives only "no K⁻¹ term", is stated.

I recomputed ς₁ and C² by hand:
- −61/432, −631/11664 and 13/4;
- 27/16, 729/250 and 4/27;
- −13/125, −11/64 and −187/2000;
- ς₁ = −0.18487 for τ = ¼, v = 2v_hs (Table S1).

All match. The c₂ note's cross-references to the paper (Lemmas 4.1, 4.3, 4.6, 4.10, Theorem 4, open problems (1), (3), (4)) match the compiled numbering.

The only issue is presentational (D3ew-5): "proved" rests on a markdown file in `review/` of the code release.

## Status of the round-135/136 issues
- **Still resolved:** I2, I3, I7, I8, I9, I11, I14 (up to D3ew-1's η = 0 point), I15, J2, J7, J8, J9, J12.
- **Still resolved, now about a cross-check:** I4, J3, J4, J5, J6. I re-ran the verifier: VERIFIED in 13 s.
- **I5 / J12:** further improved; 28 pp.
- **I12:** my build has 0 overfull boxes and no undefined references.
- **I6 / J10:** the restructuring introduced new clashes (D3ew-4).
- **I10 / J11:** placeholders remain (out of scope), but there are new revision-history artefacts in the PDFs (D3ew-3).
- **I1 / J1, I13:** `% CHECK` items, out of scope.

## New findings
1. **D3ew-1 (minor; l.445, l.452, l.511, l.519; pp. 16, 18–19).** Lemma 4.6 (ii) proves local uniformity only near η₀ ≠ 0: (4.4) has η^{−2K}, and the scaled path (η/η₀)C₀ divides by η₀. At η = 0 the valley condition is vacuous: the chord is {0}. Yet Theorem 2 claims uniformity on compact subsets of the whole disc, and the Cumulants step needs the closed disc |η| ≤ ½. The fix is one sentence: ĝ_K(√Kη)e^{−Kη²/4} − (1 − η²/4)^{1/2} is analytic on |η| < 2, so uniform convergence on a circle |η| = r gives it on |η| ≤ r (maximum principle, as already used in Lemma 4.3). Alternatively, define Ω′ inside 0 < |η| < 2. The same remark applies to the c₂ note's Proposition A (iii). This gap predates the restructuring and was not raised in rounds 135–148.
2. **D3ew-2 (minor; abstract l.57, p.1; l.248, p.9; notation table l.220).**
   - The text says "a three-segment path on which the exponential rate of the integrand is monotone". On the first segment F_η rises to the saddle and then falls, so it is not monotone.
   - F_η is also not the rate but an upper bound for it beyond the unit circle. Λ⁺ strictly exceeds the true rate there: on the imaginary axis Λ⁺ = H(|ν|) > 0, while c_K ≈ 1 + y² (supplement S3).
   - Suggested wording: "…along which an upper bound for the integrand's exponential rate decreases away from the saddle". The table entry should read "the rate bound".
3. **D3ew-3 (minor/editorial; supplement p.1 l.50 and Appendix A l.129; paper l.98, l.248, l.631, l.746).** Revision history and stale qualifiers are visible in the PDFs:
   - the supplement: "until referee round 148" and "This appendix was the paper's Appendix A until referee round 148";
   - the paper's introduction: "(certificates …, from an earlier version of the proof)";
   - the code statement: "An earlier version of the proof of Theorem 2 covered …", ending in a dangling "and V on 1.8927 ≤ |η| ≤ 1.995", which reads as if the verifier checks V; it does not;
   - Theorem 4's proof: "Theorem 2 on the axes (by hand)", now redundant.
4. **D3ew-4 (editorial; l.221, l.502, l.537, l.560).** The new notation clashes with existing symbols:
   - Π_η (the path) against Π (the constant in Corollary 3's proof);
   - m(η) (the margin) against m(ξ) (the integrand) and m_K^± (Lemma 4.11). The table admits these clashes instead of removing them;
   - r = ln(σ*/σ) is reused as the vertical coordinate in Proposition 4.8's proof (l.502);
   - a, b in Lemma 4.7 (also the Gaussians of Proposition 2.2), and p on the horizontal segment;
   - the stale parenthetical "(…precisions, not paths)" at l.560: no γ-paths remain.
5. **D3ew-5 (minor; Table 1's ς₁ row, §5 before Table 3, open problem (3) l.738, supplement §S1 and §S4).**
   - "Proved [c2note]" rests on `review/RESEARCH_note_c2.md`, an unpublished markdown file in the code release. A journal will not accept this as the source of a proof.
   - Either put Propositions A and B and Corollary B′ in the supplement as a proved section, or write "established in a companion note (not part of this paper)".
   - Open problem (3) also carries a proof sketch, which belongs in the supplement.
6. **D3ew-6 (cosmetic; README.md l.18, l.24, l.44).**
   - "Build of round 136" heads a table that the next line says is the build of round 148.
   - "§S5 … Lemma A.14" should be Lemma S5.1.
   - "Lemmas A.1–A.12, Propositions A.3 and A.13": A.3 is a proposition.
   - The supplement's subsection numbers A.1–A.4 coincide with Lemmas A.1–A.4, so "the commands (A.4)" is ambiguous.
7. **D3ew-7 (editorial, optional; l.513).** The "axes" paragraph of Theorem 2's proof is no longer part of the proof; the disc paragraph covers the axes. Say that it records the gaps used for Theorem 4 (d), or shorten it.

## Verification I did
- **Build.** pdflatex ×2–3 on both files: 28 and 12 pp., 0 overfull boxes, no undefined references. Text identical to the committed PDFs (pdftotext diff empty). I checked all cross-references against the `.aux` files, including the supplement's references to the paper's numbers (Theorems 1, 2, 4 (c), Lemmas 4.3, 4.6, Definition 4.5, Proposition 4.8, Section 5) and Appendix A's numbering.
- **F_η along Π_η** (my own code, Arb at 300 bits, F_η taken from the paper's definitions with Λ via M), in `…/ref154/mine/path_check.py` and `path_check.out`:
  - 224 values of η: R from 1 + 10⁻¹² to 2 − 10⁻¹², 14 angles including 0, 10⁻¹², π/2 − 10⁻¹² and π/2; about 480 points per segment, clustered at the ends;
  - 0 failures. The first segment is unimodal, and F < F₀ off the sub-chord;
  - on the vertical and horizontal segments F strictly decreases, with maximum equal to m(η) at η/R to 10⁻⁷⁹ (e.g. m = −2.5·10⁻³⁷ at R = 2 − 10⁻¹², θ = 0; −5.0·10⁻²⁵ at θ = π/2).
- **All four quadrants without symmetry, and the axis formulas** (`quad_control.py`, `quad_control.out`):
  - 65 cases, 0 failures, and m < bound < 0 in each;
  - the axis formulas and both axis bounds agree with (4.5);
  - straight-path control: the maximum is +4.3·10⁻³ at (1.95, 0.3) and +4.1·10⁻⁵ at (1.95, 0.097), but negative at (1.8927, 0.3), (1.95, 0.05) and (1.95, 0.7). This matches S3's failing sector and D's radius.
- **End to end, where the new proof matters** (`finiteK_thm2.py` and `_big.py`; exact Gauss–Hermite rules in Arb, up to 9000 bits): |ĝ_K(√Kη)e^{−Kη²/4}/(1−η²/4)^{1/2} − 1| at K = 20 … 640:
  - 1.5e^{0.5i}: 0.081 → 0.0022;
  - 1.9e^{0.3i}: 0.47 → 0.014 (a point where the straight path fails);
  - 1.95e^{0.3i}: 0.53 → 0.020;
  - 1.95e^{1.2i}: 0.0070 → 0.00022;
  - 1.98e^{0.2i}: 0.89 → 0.22, slow near the real axis as S4 says.
- **Verifier.** `tools/verify_certificates.py` gives VERIFIED in 13 s, exit 0, all bounds as stated. Its docstring and `results/certificates/README.md` call it a cross-check.
- **Number check.** There is no LaTeX number-check script (`tools/numcheck.py` works on markdown drafts), so I wrote `…/ref154/mine/texnum.py`.
  - Of 76 distinct decimal tokens in the paper and 114 in the supplement, all are verbatim in `gh_peano.md`, `frameless.md` or `width_mismatch.md`. The exceptions are the marked NOT-IN-MD tokens (0.528, 0.203, 8.8·10⁻¹⁴, 7.3·10⁻¹³, 5.8·10⁻¹³), arXiv ids, cross-reference numbers and a column width.
  - I also sampled integer counts (402664, 160111, 16083, 1306, …) and the exact fractions in Table 3; all are present.
  - The mds' generated text is updated to "by hand since round 148".

## Overall assessment
**Minor revision.** The new proof of Theorem 2 is correct and complete apart from the one-line gap at η = 0 (D3ew-1), and it is a real simplification. The paper is now fully by hand with no residual dependence on D, V or VA, and its c₂ statements are accurate. What remains is wording (D3ew-2), history artefacts (D3ew-3), notation (D3ew-4), the status of the companion note (D3ew-5), and the README.

## Issues
| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ew-1 | minor | Lemma 4.6 (ii) gives local uniformity only near η₀ ≠ 0: (4.4) has η^{−2K}, and (η/η₀)C₀ divides by η₀. At η = 0 the valley condition is vacuous. Theorem 2 and its Cumulants step need uniformity on \|η\| ≤ ½, including 0. Not raised before. | One sentence in Theorem 2's proof (l.511): maximum principle for the analytic difference on \|η\| ≤ r; or define Ω′ in 0 < \|η\| < 2. Same in the c₂ note, Proposition A (iii). |
| D3ew-2 | minor | "a three-segment path on which the exponential rate of the integrand is monotone" (l.57, l.248): F rises to the saddle on the first segment, and F_η is an upper bound for the rate (Λ⁺ > true rate off the disc, S3). Table: "F_η … the rate" (l.220). | "…along which an upper bound for the integrand's exponential rate decreases away from the saddle"; in the table, "the rate bound". |
| D3ew-3 | minor/editorial | Revision history in the PDFs: supplement l.50 and l.129 ("until referee round 148"); paper l.98 ("from an earlier version of the proof"); l.746 ("An earlier version…", dangling "and V on …"); stale "(by hand)" at l.631. | Rephrase as "independent computer-assisted checks of the valley condition along other paths (supplement, Appendix A; V in §S5, checked by its proposing program)"; drop "(by hand)" at l.631. |
| D3ew-4 | editorial | New notation clashes: Π_η / Π (l.537); m(η) / m(ξ) / m_K^±; r = ln σ* reused at l.502; a, b and p in Lemma 4.7 and Proposition 4.8; stale "(…precisions, not paths)" at l.560. | Rename Π → N_Π, m(ξ) → q(ξ), r → s at l.502; drop the parenthetical; update the table (l.221). |
| D3ew-5 | minor | "Proved [c2note]" (Table 1, §5, open problem (3), S1, S4) cites a markdown file in `review/` of the code release; open problem (3) carries a proof sketch. | Put the note's Propositions A and B and Corollary B′ into the supplement as a proved section, or phrase as "established in a companion note, not part of this paper"; shorten open problem (3). |
| D3ew-6 | cosmetic | README: "Build of round 136" header over round-148 numbers; "Lemma A.14" should be S5.1; "Lemmas A.1–A.12" includes Proposition A.3; supplement subsections A.1–A.4 coincide with Lemmas A.1–A.4 ("the commands (A.4)"). | Fix the README lines 18, 24, 44; write "Appendix A.4 (Runtimes and commands)" or un-number those subsections. |
| D3ew-7 | editorial (optional) | The "axes" paragraph (l.513) is no longer part of Theorem 2's proof. | Label it as the record of the axis gaps for Theorem 4 (d), or shorten it. |

## Sources
- **Paper and supplement:** `manuscript/paper18/paper18.tex` (l.57, 98, 220–221, 248, 425–520, 535–600, 631, 738, 746), `supplement.tex` (l.50–129), `README.md`, and the committed PDFs.
- **Review files:** `review/RESEARCH_note_annulus.md`, `review/REVIEW_round148_D3eq.md`, `review/RESEARCH_note_c2.md` (l.1–80, 160–185, 360–404), `review/REVIEW_round139_D3eh.md`, `review/REVIEW_round135_D3ed.md`, `review/REVIEW_round136_D3ee.md`.
- **Results and tools:** `results/chaos/gh_peano.md`, `frameless.md`, `width_mismatch.md`; `tools/verify_certificates.py`, `results/certificates/README.md`, `tools/numcheck.py`.
- **Diff:** `git diff 638fb1b b83e915 -- manuscript/paper18/`.
- **My files**, in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref154/`:
  - `mine/path_check.py` and `.out`;
  - `mine/quad_control.py` and `.out`;
  - `mine/finiteK_thm2.py` and `.out`, `mine/finiteK_thm2_big.py` and `.out`;
  - `mine/texnum.py`;
  - `mine/verify.out`;
  - `build/` (my compile, `.aux` and `.log` files).

## Lesson
When a hard region of a proof is replaced by a clean argument, the old scaffolding leaves two kinds of residue:
- **Text residue.** The new proof's own wording ("monotone", "the rate") is written by people who know which parts are bounds. And "by hand", "earlier version" and "round 148" are history that readers of the new paper cannot place.
- **Degenerate points.** The old proof was never stress-tested at them, because the effort went into the hard region. Here that is η = 0, where the scaling argument silently divides by η₀.

After a big simplification, re-read the statement's full range, including its trivial points, as carefully as the new lemma.
