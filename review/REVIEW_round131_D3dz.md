# Referee report D3dz, round 131: Proposition D‴.61 as revised after round 127

## Scope

**What I reviewed.** Commit a15cc7e on claude/vigilant-clarke-hf939b.
- I exported it with `git archive a15cc7e | tar -x` into W/repo, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref131. W/repo/venv is a symlink to the repository's venv (Python 3.13.18, python-flint 0.9.0; mpmath is not installed, so I used Arb).
- `git diff --stat eb6e7d3 a15cc7e` is empty for every D‴.61 file, so a15cc7e carries exactly round 127's application (eb6e7d3).
- The repository's HEAD moved to aa96ad3 during the round because of other sessions. I reviewed only a15cc7e.
- One correction to the brief: the (iii) header "proved on the open disc — computer-assisted beyond |η| < 1 and the axes" was already in 31487a6; eb6e7d3 did not change it [own: `git show`].

**Hashes (sha256 prefixes) [own]:**
- gh_peano.md 6df1ab16d26158f1
- gh_peano_raw.npz 643069e12657d319
- exp_gh_peano.py 4667a29ba0512540
- gh_lastann.py f81fcf3eb545eefa
- test_gh_peano.py bfe45fd30b5f02a6
- D61.md f0f9f25e5365f515
- round127_D61_own.py 90440efa6d3b50fe

**How I ran things.**
- Read-only git only.
- Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with GP_OUT = W/out.
- Never more than one heavy process at a time. I killed no process; every process of mine exited on its own.

**Read in full:**
- D61.md;
- REVIEW_round127_D3dv.md and the Round 127 entry of REVIEW_LOG.md;
- `git diff 31487a6 eb6e7d3` for D61.md, exp_gh_peano.py, gh_lastann.py, test_gh_peano.py, gh_peano.md and the new round127_D61_own.py;
- in exp_gh_peano.py: `va_seams`, `va_lines`, `report`, `v_lines`, `part_VA`, `part_DB`, `cert_box`, `cert_band`, `cert_seg_parts`, `cert_E`, `cert_end_ball`, `dom_sector`;
- in gh_lastann.py: `I`, `D`, `pos`, `mvbound`, `piece_outer`, `piece_end`, `handoff_range`, `handoff_check` (head), `soundness`, `validate` (head), `main`;
- tests: `test_part_va_the_last_annulus_certificate`, `test_part_va_series_against_arb`, `test_round127_duhamel_identity_moment_by_moment`;
- D58.md where it uses D‴.61, and the HANDOFF lines on D‴.61.

## Summary

**Design.** Round 127 made these changes:
- a proof of Duhamel's identity for the infinite Hermite Jacobi matrix, by telescoping on banded matrices, placed before "(a) He_K(x) = E(x + iZ)^K";
- `va_seams`: the whole-disc sentence is printed only if D is resolved in balls, V is contiguous and passes in balls, and VA's R1 equals V's radius;
- seam plants and threshold plants (P7, the 2⁻¹⁴ shift, the clamps on every side, the budget);
- distinct versus point–branch counts in `soundness`;
- citations (O'Sullivan 2019; Bennett–Howls–Nemes–Olde Daalhuis);
- a moment-by-moment test of the identity.

**Verdicts.**
1. **Reproducibility is exact** [own].
   - `--from-raw` regenerates the md byte for byte.
   - 43/43 tests pass (30.3 s).
   - The shipped mutants reproduce: S01 and S03–S06 are killed, and S02 survives (equivalent, as logged).
2. **Duhamel's proof paragraph is correct for the infinite matrix.** I re-derived every step.
   - Arb, 256 bits, rigorous integration: it holds at 9 points (K = 1…40, both axes, the sector, |η| up to 1.99), relative difference ≤ 2.1·10⁻⁵³.
   - Its s-form holds at the same points.
   - It holds exactly in integers for K ≤ 7, m ≤ 3K + 5 (126 identities, 0 mismatches).
   - The s-form gives (ii)'s variance exactly.
   - Only tidy-ups remain (D3dz-7). The citation of the finite-matrix precedent should change (D3dz-6).
3. **Nothing else in (iii)–(v) is load-bearing and unproved,** apart from the three certificates, which are labelled computer-assisted.
   - Several classical inputs are cited by name or [memory] only (D3dz-6).
   - The labels are accurate, with two wording exceptions (D3dz-8).
4. **`va_seams` is sound on the shipped data and stops the four planted seams.**
   - Of my 25 plants on the raw, it stops every one aimed at the seams it checks. That includes a 1-ulp gap between V's bands.
   - It still passes:
     - a V band without its end piece in balls (default `ball_end = −1.0`);
     - a VA raw missing its inner pieces, `validate`, soundness or meta records;
     - a mismatch of meta.ε₀ against O.R1.
   - Five single-point mutants of it survive and are unsound in principle (D3dz-5). None changes a published result.
5. **Mutants.**
   - Mine (Y01–Y18, full test file): 7 killed, 11 survive. Of the survivors, 5 are unsound in principle, 1 (Y06) only weakens the "in ball arithmetic" claim, and 5 are equivalent.
   - Round 127's gh_lastann mutants rerun on the current code: its unsound survivors Z02–Z06 are now killed, and Z07 too. Z08 (sound), Z09 and Z13 (equivalent) and Z10 (the known tolerance survivor S05) survive.
6. **The deferred residue does not affect correctness.**
   - My own unit checks of `I`, `D.pos`, difference quotients and `mvbound` find 0 failures.
   - **VA at ε₀ = 0.11 passes in 6.5 s** [own], so D and VA alone cover |η| < 2 and V is redundant (D3dz-3).
   - The standalone CLI `gh_lastann.py all` crashes at the outer piece (D3dz-4).
7. **Journal readiness.** The mathematics is ready. The paper is not: two publication-major items (D3dz-1, D3dz-2).

**No correctness issue was found.**

**Findings:**
1. **D3dz-1 (major for publication; not a correctness defect).** The inequalities the computer verifies appear only in code docstrings, not in the text: D's closed-form supremum and box bound, V's segment bound, VA's forms A and B, their closed forms f₀ and the mean-value form, the handoff and the tiling.
2. **D3dz-2 (major for publication).** Structure and length.
   - The statement alone is 2434 words, with 13 referee-round mentions.
   - Theorem, heuristics, superseded numerical routes and verification history are interleaved.
   - Two analytic routes are used where one suffices.
3. **D3dz-3 (minor).** V is redundant: VA at ε₀ = 0.11 (R1 = 1.89 < 1.8927) passes every piece [own]. This has been deferred since round 121.
4. **D3dz-4 (minor).** The documented re-run path `gh_lastann.py all` crashes: `piece_outer` logs `json.dumps` of complex joints (l. 832). Separately, at ε₀ = 0.11 `validate`'s B agreement is 1.1·10⁻¹², so the tests' 10⁻²⁰ bar is specific to ε₀ = 0.005.
5. **D3dz-5 (minor).** Residual gaps in `va_seams` and `va_lines`, and unsound-in-principle survivors Y01, Y05, Y07, Y08, Y09.
6. **D3dz-6 (minor).** Load-bearing classical inputs are cited without a precise reference.
   - Golub–Welsch has no reference at all (D3dv-5(f) residue).
   - The complex-phase Laplace method is cited without a theorem number.
   - Papadakis (Paulsen–Raghupathi Thm 2.10) is [memory].
   - Mehler, the Fock reproducing kernel and the Gauss weight formula are uncited.
   - The finite-matrix Duhamel precedent is cited as Saad 1992 [memory], but the identity proper is the defect-based integral representation (Jawecki–Auzinger–Koch) [snippet].
7. **D3dz-7 (cosmetic).** Tidy-ups in the Duhamel paragraph, and an undersold s-form.
8. **D3dz-8 (cosmetic).** Wording:
   - the md's "the valley condition holds … on the whole open disc";
   - "(D, rounded outward)" for float bisections;
   - the md's VA header lacks round 127;
   - D58's narrower assumption sentence;
   - `va_seams` crashes on the raw's numpy array.

## What I checked and how

**(1) Reproducibility** [own]
- `GP_OUT=W/out exp_gh_peano.py --from-raw`, then `cmp` against results/chaos/gh_peano.md: identical.
- `pytest -q -p no:cacheprovider tests/test_gh_peano.py`: 43 passed in 30.32 s (W/logs/pytest_gh.log).
- Shipped mutants (`tools/mutants.py --spec tools/mutants/round127_D61_own.py -k "part_va or md_is_the_report"`): S00 is clean; S01, S03, S04, S05 and S06 are killed; S02 survives (W/logs/mut_shipped.log).

**(2) Duhamel's identity, independently** [own: re-derivation; W/scripts/duhamel_arb.py, W/logs/duhamel_arb.log]

Step by step:
- **The telescoping.** Banded infinite matrices form an associative algebra, since every entry of a product is a finite sum. So J^m − J₀^m = Σ_{j<m}J^{m−1−j}VJ₀^j holds entrywise, with no operator theory needed.
- **The entries.** (J^m)_{00} = ⟨p₀, z^m p₀⟩ = E Z^m. J₀ is block diagonal, so (J₀^m)_{00} = (J_K^m)_{00} = Q_K z^m (Golub–Welsch).
- **The coupling.** J₀^je₀ is supported on levels 0…K−1, so V acting on it gives √K·x_{K−1}e_K. Hence (E − Q_K)z^m = √KΣ_j(J^{m−1−j})_{0K}(J_K^j)_{K−1,0}.
- **The two corners.**
  - (J^n)_{0K} = E[Z^n p_K] = E[Z^nHe_K]/√K!.
  - (J_K^j)_{K−1,0} = Σ_k w_k p_{K−1}(z_k)z_k^j (from the eigenvectors (p_i(z_k)√w_k)_i of J_K) = Q_K[z^jHe_{K−1}]/√(K−1)!.
  - √K/(√K!·√(K−1)!) = 1/(K−1)!. Correct.
- **The passage from J^m to e^{ξJ}.**
  - The proof never forms the operator e^{ξJ}; it sums scalar series.
  - Σ_m ξ^m E Z^m/m! = E e^{ξZ} by dominated convergence (E e^{|ξZ|} < ∞).
  - The double series Σ_{n,j}|a_n||b_j|·|ξ|^{n+j+1}/(n+j+1)! equals ∫_0^{|ξ|}[Σ_n|a_n|(|ξ|−t)^n/n!][Σ_j|b_j|t^j/j!]dt ≤ |ξ|·E[e^{|ξ||Z|}|He_K|]·Q_K[e^{|ξ||z|}|He_{K−1}|] < ∞. So Fubini–Tonelli allows the interchange along u = tξ/|ξ|. The text names only the first factor; the second is a finite sum.
  - If the operator reading ("the Lanczos approximation of e^{ξJ}e₁") is wanted: ‖J^me₀‖² = E Z^{2m} = (2m − 1)!! ≤ 2^m m!. So Σ|ξ|^m‖J^me₀‖/m! < ∞, and e₀ is an entire vector, with e^{ξJ}e₀ = (ξ^ne^{ξ²/2}/√n!)_n ∈ ℓ².
- **The generating functions.** E[e^{wZ}He_K] = w^Ke^{w²/2}, from e^{wz − w²/2} = Σ He_n w^n/n! and orthogonality. The Beta integral is ∫_0^ξ(ξ − u)^n u^j du = n!j!ξ^{n+j+1}/(n+j+1)!.
- **c_K.** c_K(u) = Q_K[e^{uz}He_{K−1}]/(u^{K−1}e^{u²/2}). Its Dirichlet form u^{K−1}E e^{uT} follows from w_kHe_{K−1}(z_k) = (K−1)!/He_K′(z_k) and Hermite–Genocchi.
- **The change of variables.** u = ξs gives (ξ²((1−s)² + s²))/2 − ξ²/4 = ξ²(1 − 2s)²/4 and (2K)!/(K!(K−1)!) = K·C(2K, K). Correct.
- **The s-form.** S =_d sT + (1 − s)Z with s ~ Beta(K, K + 1), since 1/B(K, K + 1) = K·C(2K, K), and both sides have entire MGFs. Correct.
- **New consequence.** Var S = E s²·Var T + E(1 − s)² = (K² + K + 1)/((K + 1)(2K + 1)) exactly, which is (ii)'s variance in one line. Checked exactly for K ≤ 12 (W/logs/sform_var.log).

Numerics in Arb at 256 bits:
- **Method.** Nodes from `fmpz_poly(He_K).complex_roots()`. LHS = e^{ξ²/2} − Σw_ke^{ξz_k}. RHS by `acb.integral` (rigorous) along [0, ξ]. The s-form by `acb.integral` on [0, 1].

| K | ξ | rel. diff. LHS vs RHS (upper bound) |
|---|---|---|
| 1 | 0.7 + 0.3i | 6.7·10⁻⁷⁶ |
| 3 | 1.5 − 2i | 2.4·10⁻⁷⁵ |
| 5 | 2 + i | 3.8·10⁻⁷⁴ |
| 8 | 3.1 + 0.4i | 9.2·10⁻⁷³ |
| 12 | √12·1.9e^{0.3i} | 3.4·10⁻⁷⁰ |
| 20 | √20·1.95e^{0.35i} (the sector) | 5.8·10⁻⁶⁶ |
| 20 | i√20·1.8 | 2.2·10⁻⁷³ |
| 30 | √30·1.99e^{0.2i} | 2.8·10⁻⁵⁸ |
| 40 | √40·1.5 | 2.1·10⁻⁵³ |

- Every pair of balls overlaps.
- The s-form against the exact ĝ_Ke^{−ξ²/4}: ≤ 2.1·10⁻⁵³.
- c_K by Q_K against the Dirichlet divided difference: ≤ 4.8·10⁻⁶³.
- The exact integer check of the telescoped moment identity in the monic basis: 0 mismatches for K ≤ 7, m ≤ 3K + 5.

**(3) What in (iii)–(v) is proved, cited, computer-assisted or heuristic** [md, own reading]

| step | status | label in D61 |
|---|---|---|
| c_K's Dirichlet form; Var T = (K − 1)/(K + 1) | proved; uses the Gauss–Hermite weight formula (uncited) and Hermite–Genocchi (de Boor [snippet]) | accurate |
| Q_K f = e₀ᵀf(J_K)e₀ | cited by name only ("Golub–Welsch"; no reference) | — |
| Duhamel's identity and its s-form | proved in the text (round 127) | accurate |
| (a): He_K = E(x + iZ)^K, the shifted representation, the Gaussian bound, I_K ≠ 0 on the ellipse, the circle identity, the maximum principle | proved | accurate |
| (a): the Laplace / saddle-point method with complex phase, parameter-dependent amplitude and uniform O(1/K) | cited (O'Sullivan 2019; Bennett et al. [snippets]), no theorem number | "proved" — acceptable given a precise theorem |
| (b) Ω: the contour bound → Λ⁺, the band near \|σ\| = 1, openness, both axes | proved | accurate |
| (b) Ω on 1 < \|η\| ≤ 1.8927 | computer-assisted (D, in balls); the closed form and the box-bound lemma are stated only in code | labelled |
| (b) Ω′: Cauchy deformation, the chord by Laplace, the 0-side | proved | accurate |
| (b) Ω′: path b → η on 1.8927–1.995 and 1.995–2 | computer-assisted (V, VA); V's segment-bound inequality and VA's forms are stated only in code | labelled |
| uniformity on compact sets | proved (finitely many boxes, margins) | accurate |
| the violated sector at 1.95 and the maximiser's Re(ηs) | numerical (float bisection of the closed form), not used | "(D, rounded outward)" suggests rigour (D3dz-8) |
| Λ⁺ attained; the imaginary-axis alias rate; the η = 2 Stieltjes saddle | heuristic | labelled ("conjecturally", "heuristic", "formal") |
| VC's four regimes | heuristic support for VA's design; Lemma VC is proved but not in the chain | labelled |
| (c) | proved (Cauchy's estimates) | accurate |
| (iv): lower bounds | proved; the Fock reproducing kernel is classical and uncited | accurate |
| (iv): asymptotics | proved given (iii)(b) on \|η\| ≤ c | accurate |
| (iv): c_need thresholds | computer-assisted in balls (trivial) | labelled |
| (iv): t(v) ≤ 0.6345 over A6's rules | numerical, far from r* | fine |
| (v): Mehler's formula | classical, uncited | — |
| (v): Parseval frame | cited, Paulsen–Raghupathi Thm 2.10 [memory]; secondary sources call it Papadakis' theorem [snippet] | — |
| (v): the rest | proved | accurate |

**(4) `va_seams` and the conclusion** [own: W/scripts/seam_plants.py, W/logs/seam_plants.log]

Behaves correctly:
- **Stopped (no whole-disc sentence), 8 plants:**
  - V without `ball_worst`;
  - V's `ball_end` NaN;
  - V's middle band starting 1 ulp above 1.98;
  - D's 1.8927 run open;
  - D's `ball_worst` > 0;
  - D in floats only (DB not substituted);
  - O failed = 1;
  - no V and no D.
- **Sound variants still printed, correctly, 2 plants:** an overlap 1 ulp past 1.995, and an extra passing band 1.5–1.995.

Gaps — **still printed, 15 plants:**
- **Realistic:**
  - V bands without `ball_end` (the default −1.0, while `ball_worst` defaults to 0.0 and `v_lines` uses 0);
  - VA without A-J, or without all four inner pieces;
  - VA without `validate`, without the soundness records, or without meta;
  - meta.ε₀ = 0.004 with O.R1 = 1.995;
  - `validate` disagreement 1.0.
  - These are guarded only by `test_part_va_…`'s pins, not where the conclusion is printed (D3dv-2's own lesson).
- **Records the code cannot produce:**
  - zero-box bands or pieces;
  - handoff sub_boxes = 0;
  - worst > 0 with failed = 0;
  - an end bound > 0 with ok = True.

Wording:
- "So the valley condition holds … on the whole open disc |η| < 2" is true only through Ω ⊆ Ω′: domination is the valley condition along the straight path, whose chord lies along η/|η|, and Re(φ″η²)/|η|² = (2Re η² − 8)/|η|² < 0. D61 asserts "Ω′ ⊇ Ω" but never writes this one-line argument, and D's branch and bound certifies domination, not the valley condition.

Robustness:
- `va_seams(raw['VA'], …)` raises ValueError (`VA or []` on a numpy array). It is reached only through `va_lines`' `list(VA)`, so it is harmless.

**(5) Mutants** [own: W/scripts/spec131_gp.py, W/mut/res_gp.jsonl; full test file]

| id | mutant | result |
|---|---|---|
| Y01 | D's `ball_worst` default −1.0 | **survives, unsound in principle** (D in floats accepted) |
| Y02 | D's `open` default 0 | survives, equivalent (key always present) |
| Y03 | R_bb = min | killed |
| Y04 | V's float `end` unchecked | survives, equivalent (`ball_end` checked) |
| Y05 | V's `ball_worst` unchecked | **survives, unsound in principle** ("both in ball arithmetic" false) |
| Y06 | V's `ball_end` unchecked | survives; the float end stays checked, but the ball claim is weakened (the shipped default is already this lax) |
| Y07 | V's contiguity dropped | **survives, unsound in principle** (a middle gap) |
| Y08 | contiguity with a 10⁻⁶ tolerance | **survives, unsound in principle** |
| Y09 | the seam R_v == R1 within 10⁻⁶ | **survives, unsound in principle** (the test's plant is 10⁻³ away) |
| Y10 | R_v ≥ R1 | survives, sound (overlap) |
| Y11 | R_v ≤ R1 | killed |
| Y12 | R_v = b['R2'] | survives, equivalent (sorted, contiguous) |
| Y13 | R1 from meta.ε₀ | killed (by the 1.996 plant's expected wording) |
| Y14 | `if seams` → `if True` | killed |
| Y15 | the printed V radius taken from VA's R1 | survives, equivalent (printed only when equal) |
| Y16 | no D passed | killed |
| Y17 | the end piece's `ok` dropped (round 127's G01) | killed |
| Y18 | evaluations printed as distinct | killed |

- No survivor changes the md or a published number.
- Round 127's gh_lastann mutants rerun (W/mut/res_la127.jsonl): Z01–Z07, Z11, Z12, Z14 and Z15 are killed.
  - Z08 (default budget 50) survives; sound, since the budget only fails boxes, and the plant passes budget = 4 explicitly.
  - Z09 (area ≥) survives; equivalent given inside and disjoint.
  - Z10 (soundness tolerance) survives; round 122's S05.
  - Z13 (the end link in floats) survives; equivalent up to rounding.

**(6) The deferred residue** [own: W/scripts/machinery.py, W/logs/machinery.log; W/scripts/va011.py, W/logs/va011b.log]
- **`I.__mul__`, `sq` and `**3`:** checked on 4000 sign-mixed float intervals against exact rational products; 0 failures.
- **`D.pos`:** checked inside (x − c)₊² and (x − c)₊³ on straddling balls; 20 000 value and derivative enclosures; 0 failures.
- **Difference quotients:** of fA and fB on 2070 points of the pieces' own accepted boxes (A-FC, B-F, A-J, B-J), each quotient inside the dual enclosure; 0 failures.
- **`mvbound`:** on 11 synthetic f with known f′ and maximum, including tight boundaries (bound −0.04 against a true maximum of −0.04, and 0.02 against 0.02), route 'd', and a NaN derivative. Sound in every case, and as expected wherever the case was decidable.
- **Verdict on G01, G03, G04, G08:** current correctness holds; the unit tests matter only as regression guards.
- **VA at ε₀ = 0.11 (R1 = 1.89).**
  - The CLI `gh_lastann.py all --eps0 0.11` crashes at `outer`: TypeError, a complex number is not JSON serializable (W/logs/va011.err). My driver (W/scripts/va011.py) calls the pieces the way `part_VA` does.
  - validate: A 1.4·10⁻²⁴, B 1.1·10⁻¹², A0/B0 ≤ 2.1·10⁻⁵⁹ (2 cases skipped).
  - Inner pieces: A-J 54, A-FC 2533, B-J 413 (549 skipped), B-F 1189 (1174 skipped); 0 failed.
  - O: 220 η-boxes and 16 083 segments, 0 failed, the largest bound −4.2·10⁻⁹. The end piece −1.3314.
  - handoff: 0 failures on 1306 sub-boxes; area exact, inside, disjoint, tiled.
  - Soundness: 324, 22 258, 3982 and 11 634 points, 0 above the bound.
  - 6.5 s in all.
  - With D in balls to 1.8927, this covers 1 < |η| < 2 without V.
- **D61's Numerics (VC) and a VC_MID point on the cap's far piece:** VC is not in the proof's chain (VA certifies directly), so these do not affect correctness.

**(7)–(8)** See Journal readiness and Sources. arxiv.org and ncbi.nlm.nih.gov could not be fetched from this session, so for those I saw search snippets only.

## Verdicts per part

| part | verdict |
|---|---|
| (i), (ii), (iii)(c), (v) | Unchanged; checked only where (iii)–(iv) use them. The s-form reproves (ii)'s variance [own]. |
| (iii) Duhamel's identity | **Proved as written** for the infinite matrix. Cosmetic tidy-ups (D3dz-7); the precedent's citation should change (D3dz-6). |
| (iii)(a) | Proved, modulo a precisely cited complex Laplace theorem. |
| (iii)(b) Ω (D) | Sound; computer-assisted in balls. Its box-bound lemma is not in the text (D3dz-1). |
| (iii)(b) Ω′ (V, VA) | Sound. V is now redundant (D3dz-3). |
| `va_seams` / `va_lines` | Sound on the data, and the four planted seams are stopped. Gaps and unsound-in-principle survivors remain (D3dz-5). |
| (iv) | Correct given (iii)(b); labelled. |
| md, D61 wording | Every number I checked is in the raw. Wording points in D3dz-8. |
| Tests | 43 pass. Survivors are listed above. |

## Status of the round-127 issues

- **D3dv-1 (Duhamel unproved): resolved.** The proof is in the text, the Status is updated, and a test checks 60-digit moments for K ≤ 6. Residue: the citation (D3dz-6) and tidy-ups (D3dz-7). Round 127's suggested test of the Beta-mixture moments for K ≤ 8 was replaced by an equivalent moment test.
- **D3dv-2 (seams unchecked): resolved for the four planted seams and the end piece.** G01 (my Y17) is killed. Residue: D3dz-5.
- **D3dv-3 (plants far from thresholds): resolved.** Z02–Z07 are now killed; Z08 is sound.
- **D3dv-4 (distinct counts): resolved.** 18 326 distinct points and 24 770 evaluations; HANDOFF has 29 056 and 35 500.
- **D3dv-5 (wording): resolved for (a)–(e).** (f) is half done: the Laplace method is cited, Golub–Welsch is still not.
- **D3dv-6 (HANDOFF): resolved.** The parts list includes VC and VA (l. 62), the assumed arithmetic is stated (l. 325), and the counts are corrected (l. 206).

## Journal readiness

As a referee for SIAM J. Numer. Anal. or Math. Comp., I would require the following before acceptance.

1. **Split the Proposition.**
   - **Theorem 1:** the Peano kernel as the law of S, with the s-form S = sT + (1 − s)Z, and the moments. The s-form is undersold: it is the cleanest structural fact here and gives Var S in one line.
   - **Theorem 2:** the corner ratio (a).
   - **Theorem 3:** ĝ_K(√Kη)e^{−Kη²/4} → (1 − η²/4)^{1/2} uniformly on compact subsets of |η| < 2, computer-assisted beyond |η| < 1 and the axes.
   - **Corollary 4:** (iv)–(v), and D‴.58's rates and sharp constants.
   - **Size.** The statements are 2434 words now; aim for under 400. Remove every "referee round" from the mathematics (42 occurrences, 13 of them in the statement) into a verification-history note.
2. **One route, two certificates.**
   - Domination is the valley condition on the straight path. State a single "path lemma": a path 0 → η through η/2, with a chord inside the disc along a descent direction and F_η < F_η(η/2) off η/2, gives the limit. Then D's branch and bound is the straight-path verification, and the separate Laplace argument for Ω shrinks to a remark.
   - Delete sufficient condition (1), the numerically checked grid of 57 points, which is superseded.
   - With VA at ε₀ = 0.11 [own], drop V, or keep it as an independent cross-check in an appendix. The proof then needs D on 1 < |η| ≤ 1.8927 and VA on 1.89 ≤ |η| < 2.
3. **The computer-assisted appendix.** For each certificate:
   - the exact inequality verified per box: D's closed-form supremum over s and the box bound's monotonicities; VA's forms A and B as identities, f₀ in closed form, the mean-value form and its domains; the connector and handoff; the rational tiling;
   - proofs that the implemented bounds are upper bounds, including the series tails;
   - which steps only propose (floats choosing boxes and paths) and which certify (balls, rationals);
   - the trusted base: Arb's ball operations and elementary functions, python-flint 0.9.0's exact float→arb conversion, IEEE-754 doubles with `math.nextafter` where floats feed exact inputs, `fractions`; and that the programs implement the stated inequalities, which the appendix makes checkable;
   - runtimes: DB about 15 s, VA about 5–7 s; V 6 + 29 + 267 s if kept (round 127's timings);
   - a pinned environment (tools/requirements-pinned.txt exists), an archived code and data release with a DOI, and one command per certificate. The standalone `gh_lastann.py all` must work (D3dz-4).
   - Model the presentation on Rump's Acta Numerica survey [snippet] and on Arb-based work in SIAM journals (Johansson–Mezzarobba [snippet]).
4. **Proof details to add.**
   - Ω ⊆ Ω′ in one line.
   - A precise theorem for Laplace's method with complex phase and K-dependent analytic amplitude.
   - Golub–Welsch with a reference.
   - Tonelli for Duhamel's double series.
   - Probably an O(1/K) rate on compact sets. Laplace's method with (a)'s uniform O(1/K), and exponentially small path pieces, gives it; it is undersold as o(1).
5. **Overreach and underselling.**
   - "Conjecture P is proved" is fine with "computer-assisted".
   - The md's "the valley condition holds on the whole open disc" should read "the limit holds".
   - "(D, rounded outward)" should read "numerically" for the float bisections.
   - "New, as far as these sources show" rests on snippets; say "not found in …".
6. **D‴.58.**
   - Cite Theorem 3 and Corollary 4 by number. State the result once, for every t < 1, computer-assisted, and move the t < 0.0966 / 0.8052 / 0.9900 tiers to a remark.
   - Align its assumption sentence ("Arb's correctness assumed", ll. 16 and 176) with D61's full list.
   - Keep moment matching and |ρ| ≥ 1 explicitly conjectural, as now.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dz-1 | major (publication) | **The verified inequalities are not in the text.** D's box bound, V's `cert_seg_parts` bound (\|g′(m)\| + ℓ sup\|g″\| with \|1 − z²\| ≥ \|1 − m²\| − ℓ(2\|m\| + ℓ); L_η via h′), VA's forms A/B, f₀, the mean-value form, the handoff and the tiling are stated only in docstrings [md; own reading]. | D61: an appendix with one lemma (statement and proof) per bound, pseudo-code, the trusted base and one command per certificate. |
| D3dz-2 | major (publication) | **Length and structure.** Statement 2434 words, proof 2110 words, 42 "referee round" mentions; route (1) numerical but in the statement; heuristics (alias, η = 2, VC) interleaved; Ω and Ω′ argued separately [own: word counts]. | Restructure as in Journal readiness items 1–2: Theorems 1–3, Corollary 4, Appendix A (certificates), Appendix B (heuristics and numerics). |
| D3dz-3 | minor | **V is redundant.** VA at ε₀ = 0.11: inner pieces 54/2533/413/1189, O 220/16 083 with worst −4.2·10⁻⁹, end −1.3314, handoff 0 failures on 1306 sub-boxes, tiled, soundness 0 above; 6.5 s; R1 = 1.89 < 1.8927 [own: va011b.log]. Deferred since round 121. | exp_gh_peano: a part VA011 (or VA's ε₀ = 0.11) with its seam to D (R1 ≤ R_bb). D61: D + VA as the proof, V as a cross-check. Tests: pin the counts. |
| D3dz-4 | minor | **The documented CLI crashes.** `gh_lastann.py all` fails at `outer`: `log(json.dumps(out))` on complex joints (l. 832) [own: va011.err]. Also, at ε₀ = 0.11, `validate`'s B_vs_centred is 1.1·10⁻¹², so the 10⁻²⁰ bar is ε₀-specific [own]. | gh_lastann `piece_outer`: dump without `joints` (or with `default=str`). A CLI smoke test. Document `validate`'s tolerance as a function of ε₀ (the series tails). |
| D3dz-5 | minor | **Seam-check gaps.** `ball_end` defaults to −1.0 (lax); `ok` checks only the records present (a VA raw without its inner pieces, validate, soundness or meta still prints); meta.ε₀ = 2 − O.R1 unchecked; validate not thresholded [own: seam_plants.log]. Unsound-in-principle survivors Y01, Y05, Y07, Y08, Y09 [own: res_gp.jsonl]. No published change. | `va_seams`/`va_lines`: default `ball_end` to 0.0; require the set {validate, A-J, A-FC, B-J, B-F, O, end, handoff, sound ×4, meta} with boxes > 0; require O.R1 == 2 − meta.eps0; threshold validate (< 10⁻²⁰ at ε₀ = 0.005). Tests: plants for a V band without `ball_worst`, a middle band missing, V ending at nextafter(1.995, 0), D's bb in floats, and a missing A-J. |
| D3dz-6 | minor | **Imprecise citations of load-bearing classical inputs.** Golub–Welsch named only; the complex Laplace method without a theorem number; Papadakis [memory] (secondary snippets: Paulsen–Raghupathi Thm 2.10); Mehler, the Fock reproducing kernel and the Gauss weight formula uncited. The finite-matrix Duhamel form is cited as Saad 1992 [memory]; snippets show Saad's Thm 5.1 and Prop. 2.1 but not the corner-entry statement, while Jawecki–Auzinger–Koch write the defect-based integral representation L_m(t)v = ∫E(t − s)D_m(s)v ds [snippet]. | D61 Classical: Golub–Welsch 1969; a theorem of O'Sullivan 2019 (or Olver's textbook [memory]); Paulsen–Raghupathi Thm 2.10 (Papadakis); Jawecki–Auzinger–Koch, BIT (2020), as the finite-matrix identity; Saad 1992 as the a posteriori estimate. |
| D3dz-7 | cosmetic | **Duhamel paragraph tidy-ups.** Tonelli is named for the E-factor only; e₁ in (iii) but e₀ in the proof; "V" denotes both the coupling and part V; the operator reading is unexplained (e₀ is an entire vector, ‖J^me₀‖² = (2m − 1)!!). The s-form is undersold: it gives Var S in one line [own]. | D61 (iii) proof: one sentence on the double series; rename the coupling (e.g. Δ); one sentence on e^{ξJ}e₀. Move the s-form into (i) with the one-line variance. |
| D3dz-8 | cosmetic | **Wording.** The md's "the valley condition holds … on the whole open disc" needs Ω ⊆ Ω′, which is asserted but not argued; "(D, rounded outward)" for float bisections (`dom_sector`); the md's VA header still says "(referee rounds 120–122, 126)"; D58's "Arb's correctness assumed" (ll. 16, 176); `va_seams` raises on the raw's numpy array (`VA or []`) [md; own]. | `va_lines`: "the limit (iii)(b) holds on the whole open disc", plus the one-line inclusion in D61. "Numerically" for the sector. The header to include 127. D58: D61's assumption list. `va_seams`: `list(VA) if VA is not None else []`. |

## Sources

**Repository** [md][raw] at a15cc7e:
- review/drafts/D61.md, review/drafts/D58.md;
- review/REVIEW_round127_D3dv.md; review/REVIEW_LOG.md (Round 127); review/REVIEW_round121_D3dp.md (G01–G08);
- experiments/exp_gh_peano.py, experiments/gh_lastann.py;
- results/chaos/gh_peano.md and gh_peano_raw.npz;
- tests/test_gh_peano.py;
- tools/mutants.py, tools/mutants/round127_D61_own.py;
- HANDOFF.md;
- round 127's mutant spec (/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref127/scripts/spec127_la.py, read and copied, not modified).

**My scripts and logs** [own]
- All under /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref131.
- Scripts in scripts/: duhamel_arb.py, seam_plants.py, spec131_gp.py, spec127_la_rerun.py, machinery.py, va011.py.
- Logs in logs/: pytest_gh.log, duhamel_arb.log, sform_var.log, seam_plants.log, mut_shipped.log, mut_gp.log, mut_la127.log, machinery.log, va011.log, va011.err, va011b.log.
- Mutant results in mut/: res_shipped.jsonl, res_gp.jsonl, res_la127.jsonl.

**Literature** (web search; snippets only, since arxiv.org and ncbi.nlm.nih.gov could not be fetched)
- **Krylov and Duhamel error formulas:**
  - Saad, SIAM J. Numer. Anal. 29 (1992) 209–228, https://epubs.siam.org/doi/10.1137/0729014 [snippet; secondary sources cite its Thm 5.1 and Prop. 2.1].
  - Jawecki–Auzinger–Koch, Computable upper error bounds for Krylov approximations to matrix exponentials and associated φ-functions, BIT (2020), https://link.springer.com/10.1007/s10543-019-00771-6, arXiv:1809.03369 [snippet: the defect-based integral representation; the defect a multiple of the next Krylov vector].
  - Jawecki, arXiv:2001.11922 [snippet].
  - Hochbruck–Lubich, SIAM J. Numer. Anal. 34 (1997) 1911–1925, https://epubs.siam.org/doi/10.1137/S0036142995280572 [snippet].
- **Gauss quadrature and Lanczos:**
  - Golub–Meurant, Matrices, Moments and Quadrature (Princeton), https://oeaw.ac.at/resources/Record/9781400833887 [snippet].
  - https://par.nsf.gov/servlets/purl/10299648 [snippet: Gauss-quadrature error estimates for Lanczos f(A)v].
  - https://www.imacm.uni-wuppertal.de/fileadmin/imacm/preprints/2015/imacm_15_21.pdf [snippet: the Lanczos error for Stieltjes functions].
- **Ball arithmetic and computer-assisted proofs:**
  - Rump, Verification methods, Acta Numer. 19 (2010) 287–449 [snippet].
  - Johansson–Mezzarobba, SIAM J. Sci. Comput. (2018), arXiv:1802.03948 [snippet].
  - Gómez-Serrano, SeMA J. 76 (2019) 459–484, arXiv:1810.00745 [snippet].
  - FLINT/Arb overview, https://flintlib.org/doc/overview.html [snippet].
- **Papadakis' theorem as Paulsen–Raghupathi Thm 2.10:** arXiv:2102.09585, arXiv:2310.18078 [snippets].
- **Not found:** the Gauss–Hermite Peano kernel's Beta-mixture (s-) form; Duhamel's identity for the infinite Jacobi matrix; conjecture P's limit; a comparable computer-assisted saddle-point limit in a numerical-analysis journal.

## Lesson

1. **Plant at the threshold, also for the joins.** A seam plant 10⁻³ away cannot tell exact equality from a 10⁻⁶ tolerance; the shipped check catches a 1-ulp gap, but no test makes it.
2. **Check completeness where the conclusion is printed.** "No recorded failure" is not "every required piece present and passing".
3. **Run the cheap redundant certificate.** VA at ε₀ = 0.11 takes 6.5 s and makes the most expensive certificate, V, unnecessary. A journal version should be built from the fewest certificates whose mathematics is written out in the text.
