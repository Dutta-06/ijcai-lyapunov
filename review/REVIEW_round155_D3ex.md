# Referee report D3ex, round 155: acceptance check of Papers 16 and 14

**Scope.** I checked commit b83e915, exported to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref155/repo` with the venv symlinked. The repository was not modified. Line numbers refer to the `.tex` files at that commit, and "S" marks the supplement. Bibliographic placeholders and `% CHECK` items are out of scope.

## Paper 16 (`manuscript/paper16/`, last report D3em, round 144)

### Status of the round-144 issues

| D3em id | resolved? | evidence |
|---|---|---|
| 1 (plain gradient flow) | yes | §1.1 l. 66 now says "the gains G and γ, a diagonally preconditioned gradient flow". |
| 2 (exponential loss of mass) | yes | l. 66 now gives the exponential rate below the critical curve with a correlated residual, and algebraic decay on the critical curve. This matches Thm 4.1(c). |
| 3 (core profile in §7.5) | yes | The 8 mismatches are now attributed to the self-consistent coefficient. The profile's errors (median 0.0259, maximum 204) are quoted from md E. |
| 4 ("direct check") | yes | §7.5 now says "consistent with … at this one starting mass (a consistency check …)". The β⁰ sequence is listed under README Open items. |
| 5 (PVRF) | yes | §1.3 describes θ⊙θ, the h-weighted Lasso, the penalty ∝ step × noise, and the O(√δ) fluctuation. The open problem (4) cites them. The \|φ_i\|²-weighted basis pursuit is stated as untested in §1.3, §7.4 and §8(4). |
| 6 (Kunin et al.) | yes | §1.3 credits App. F, and identifies (c) as stochastic mirror descent with the invariance of GLSS18/AH19. The novelty claim is narrowed. |
| 7 (supplement history) | yes | The history is gone and the literal `\n` is fixed. The disclaimer now says "(S1)–(S5) rest on runs reported in S1–S5". |
| 8 (Remark 4.7, ε²log(1/ε)) | yes, re-derived | See below. |
| 9 (φ̄'s norm) | yes | Prop. 6.2 uses \|φ\|_∞. |
| 10 (the O(t^{−1/2}) line) | yes, re-derived | See below. |
| 11 (README order, Fig. 3 caption) | yes | The README gives O(ε^{1/2}/(μ(1−μ))). The caption says "asymptotic slope … the data's local slope is still larger". |

**Prop. 4.9 (l. 382–397) [own].**
- Under the rank condition, Φ_T⁺Φ_T = I and Φw − y = Φ_T(w_T − w*_T) + Φ_{Tᶜ}w_{Tᶜ}.
- The off-T term is bounded by \|Φ_{Tᶜ}w_{Tᶜ}\| ≤ max\|φ_i\|·Σ_{i∉T}w_i, using w ≥ 0. Also Σ_{i∉T}w_i ≤ \|v\|\|Φw − y\|/μ, because Σ h_iw_i = vᵀ(Φw − y) and vᵀy = hᵀw* = 0.
- Finally \|Φw − y\|² = 2NJ/m ≤ 2NB(w*, w(0))/t by Thm 4.1(a), since J ≤ mB/t and w* ∈ 𝒮.
- So the stated bound and "at least as fast as t^{−1/2}" are correct.
- The counterexample still checks: h = Φᵀ(−1, 1) = (0, 0, 1.5).

**Remark 4.7 (l. 360–362) [own].**
- Free core: dw/dη = mD ≈ Gε²r²/m.
- Branches: log\|w_j\| = −L + s_j(Φᵀz)_j/ρ_j + O(1), with L = log(1/ε).
- |S| = N: z = Lξ + O(1), so the free coefficients are ε²L(Φᵀξ)_j, and the off-support coefficients are ε^{μ_j}, as in Prop. 4.6(b).
- |S| < N:
  - P_⊥Φw = 0 forces the active off-support atoms to ≍ ε²L, that is s_j(Φᵀz)_j = −ρ_jL + ρ_j log L + O(1).
  - Dividing the dual by ε²L² gives the barrier QP, which is strictly convex because the free columns span R⁸. It is strictly feasible by Gordan's alternative.
  - The relative correction is O(log L/L). All of this checks.
- My own SLSQP solution of the QP on S3's design gives (G/m)Σ\|Φᵀζ*\| = **230.2205**, with active set {7, 11, 13, 15, 19} (slacks ≈ 10⁻⁶; atoms 3 and 17 inactive). Both match md G.
- md G's deficits fit (c₁log L + c₂)/L, with c₁ ≈ 1.6 and c₂ ≈ 10, from L = 18 to L = 295. The exponents 1.958–1.962 at L = 294.7 against 2 − log L/L = 1.981 are consistent with an O(1/L) term.
- One flaw remains in the stated QP: D3ex-1.

### New findings
1. **D3ex-1 (minor), l. 361: the QP omits the free masses' weights.**
   - The free dual terms are (Gε²r_j²/2m)(Φᵀz)_j². The remark itself writes w_j ≈ Gε²r_j²η_j/m.
   - So ζ* must minimise Σ_free r_j²(Φᵀζ)_j², and the constant is (G/m)Σ r_j²\|Φᵀζ*\|_j.
   - The text has the unweighted sum. This is correct only for equal masses, which S3's run uses (`np.full(DC, eps)` in `part_G`, and `mixed_qp` is unweighted), so no number changes.
2. **D3ex-2 (cosmetic), l. 361 and §8(2): "At a vertex with fewer than N atoms the approach is still a power of ε, with an exponent we do not compute."**
   - For the mixed start with |S| < N, the new derivation now gives the approach: the support refits against outputs of size ε^{2−o(1)}.
   - md G's distances (1.96·10⁻⁶ and 3.98·10⁻¹³) follow the active atoms' ε^{1.5–1.6} at the ε resolved.
   - The sentence and Open problem (2) should be scoped to all-branch starts at degenerate vertices.

### Verification
- **[raw→md]** I copied `mirror_branch_raw.npz` into the scratchpad `mb_out/` directory and ran `MB_OUT=… exp_mirror_branch.py --from-raw`. `cmp` with `results/chaos/mirror_branch.md` is IDENTICAL.
- **[own run]** `pytest tests/test_mirror_branch.py`: 5 passed.
- **[own run]** `numcheck16.py results/chaos paper16.tex supplement.tex`: 302 tokens, 1 not found (10.4, derived in the text), as the README states.
- **[own run] Build:** pdflatex ×3 for the paper and ×2 for the supplement. Result 26 + 8 pp., 0 errors, 0 overfull boxes, no undefined references. `pdftotext` is identical to the committed PDFs. Open problems end on p. 23 and the availability statement is on p. 24, as the README says.
- **[own]** `qp.py`: the independent QP above.

**Recommendation: accept**, after the one-line correction D3ex-1 and optionally D3ex-2. Both are editorial and need no further review round.

## Paper 14 (`manuscript/paper14/`, last report D3ev, round 153, then Part M)

### Status of the round-153 issues

| D3ev id | resolved? | evidence |
|---|---|---|
| 1 (tautological χ²) | yes | §9.5 and S2.4 now call the χ² "bookkeeping". **[md J]** gives the per-state bias ratio at Δt over Δt/2 as 1.9939–2.0001. Without the Itô term the ratio is 1.001–1.008. This is a real test of the generator, because b_k uses the closed-form E_n[Υ]. |
| 2 (Euler bias over the long runs) | yes | Part K's decomposition (md K) is quoted correctly in S2.4: (1), (2), 12.3, 3.8 and 16.1 all check. Part M then supersedes K's reading. |
| 3 (Prop. 7.3) | yes, re-derived | See below. |
| 4 (wall layer) | yes | l. 498, Derived st. 7.2 and Prop. 7.4(a) are qualified. One residue remains: D3ex-3. |
| 5 (Prop. 7.4) | yes | τ_k > 0 is moved into (a). (b) is restricted to κ_F > 1 and its target is named in §S1.20. (c) and (d) are moved to §S1.21, with (d) stated for the chain. Table 1 l. 112 is reworded. |
| 6 ("analytic centre") | yes | No occurrence is left; W♯ is called the D-optimal fit or the weighted-barrier maximiser. |
| 7 (any chart) | yes, re-derived | See below. |
| 8 (abstract, title) | yes | The new title; the layer-law exponents; isotropic only; "drift only at order γ²τ". |
| 9 (coloured-noise errors) | yes | S2.4: 0.496 ± 0.123, 0.454 ± 0.082, 0.237 ± 0.124 (+0.2, +1.1, −0.4). |
| 10 (numbers) | yes | 0.0100 of 0.0220; "at most"; Part I's s.e. with b's estimation (md L). |
| 11 (notation) | yes | 𝓑 at l. 448, θ₀, d^S, Υ_+. |
| 12 (literature) | yes | Independent versus correlated shocks; VW16 and HJ95 (both `% CHECK`). |
| 13 (length) | yes | 30 + 26 pp. The body ends on p. 27. |

**Re-derivations [own].**
- **Prop. 7.3.**
  - Υ̇_k = Υ′α̇_k = (β_kϱ_k/G_α)(−G_α∂_α(J + V)) = β_k(u_k − ϖ_k).
  - Ẇ_k = γW_k(u_k − ΣW_lu_l) with W_k/β_k = 1/m gives the stated form.
  - Integration by parts without an Itô term, since W is absolutely continuous.
  - In stationarity E[Υ_+Υ̇_+] = E[Υ_kΥ̇_k] = 0. Because the amplitudes carry no noise, Γ(Υ_l, Υ_k) = 0 and E[Υ_lΥ̇_k] is antisymmetric, hence zero for a reversible law. It is O(γ²τ).
  - "Absolutely continuous (C¹ for continuous noise)" and "dΥ_k/dt integrable" are now in place.
- **Remark 7.5(ii).**
  - u_k = −(a_kψ′/β_k)∂_{α̃}J_fit.
  - With Υ̃′ = β_kϱ̃_k/G̃, ∂_{α̃}ϱ̃ = (1 + aψ″/ψ′)/β.
  - At the fit, (1/m)∂_W[log W − log ψ′(A/W)] = (1 + aψ″/ψ′)/(mW), which equals −(τ/2)∂ℓ_ψ by Remark 4.3's ℓ_ψ.
  - The tanh chart gives (1 + s²)/(βa′) = ϱ′. The physical chart gives τ_k/β_k. All correct.
- **Part M's offset.**
  - From Ẇ = (γ/m)(I − W1ᵀ)Υ̇, the quantity W − (γ/m)(I − W1ᵀ)Υ is conserved up to O(γ²). Hence W − Y = (1) + M(W − Y), with dΥ^fit/dW = −h²/(G_αW), h = a/a′.
  - Check: a_k = 3A_k = (3, 2.1, 2.4) and m = 3 give −diag M = 0.3135, 0.0560, 0.0938, as printed.
  - The O(γ²) remainder, of order (γ/m)Υ_+\|W − Y\| ≈ 0.1 × 10⁻⁴, matches the reported ≤ 1.4·10⁻⁵.
- **Arithmetic of Part M.** Richardson 2r(Δt/2) − r(Δt), the ratios, every z, the γ factors 1.44 and 1.61 against 1.43 and 1.52, the four-term split of (e), (2) as a rate (−6.04·10⁻⁵/(7.5·10⁵·10⁻⁸) = −0.00805), and Y − 1/3 − (2) all check by hand.

### New findings
3. **D3ex-3 (minor), Prop. 7.4(a) l. 520 and §S1.20 l. 155: the start condition "from every start … (in P_B, with the box)" needs the walls, outside Derived st. 7.2's scope.**
   - [own, `pb_face.py`] On the instance with levels (2, 1, 0.5), the ascent direction of 𝓕_τ points out of P_B on part of the face W₃ = W_{B,3} = 0.2358, namely W₁ ∈ [0.2947, 0.317].
   - For equal levels no face point is outward, so Prop. 6.3(d) is fine.
   - Trajectories from those starts are held by the box's term ϖ. That term is exactly what l. 498 and Derived st. 7.2 now exclude, since it is up to 5.65 times the noise drive at the wall.
   - Fix: "from every start whose path stays in the interior of P_B at a distance from the walls (e.g. a superlevel set of 𝓕_τ inside it); near a wall as in Prop. 6.3(c), (d), at the leading order there".
4. **D3ex-4 (minor), §9.5 l. 603 and S2.4 l. 296: the γ = 3·10⁻³ results carry an unquoted Richardson remainder.**
   - At γ = 3·10⁻³ the step effect is 4.3 times the drift: the Δt chain gives −0.04844 against Richardson's −0.00912. With two steps the O(Δt²) term is not controlled.
   - Part M(d)'s closed forms for the chains' laws are linear in Δt only to about 1%. The symmetric part is −0.00180 and −0.00089, against −0.00090 for a linear law.
   - So the remainder is about 0.0004 per γ²τ, about 4% of the drift. That is larger than the quoted ±0.02 and ±0.04 on the γ factors 1.44 and 1.61.
   - The second class's 1.61 ± 0.04 against 1.52 is already about 2.3σ. The research note says "somewhat above"; the paper says only "the finite-γ factor of Katzenberger's limit".
   - At γ = 10⁻² the effect is about 1% and the conclusions there stand.
   - Fix: state the remainder and the second class's excess, or add a Δt/4 chain.
5. **D3ex-5 (cosmetic): the fixed-mass evidence is described two ways.**
   - l. 515 says "at fixed masses the drift has the closed form's size".
   - l. 603 says "a shortfall as at fixed masses (0.86 ± 0.08, 0.93 ± 0.10)".
   - The fixed-mass values are z = +1.8 and −0.8 from the closed form, and 1.08 ± 0.20 at τ = 3·10⁻⁵.
   - Use one wording, for example "consistent with the fixed-mass ratios, which themselves lie within 1.8σ of the closed form".
6. **No error found** in §9.5, S2.4 or Table S5 against Part M and the research note. Every decimal printed there matches md M.

### Verification
- **[raw→md]** I copied `gauge_noise_raw.npz` into the scratchpad `gn_out/` directory and ran `GN_OUT=… exp_gauge_noise.py --from-raw`. `cmp` with `results/chaos/gauge_noise.md` is IDENTICAL.
- **[own, from raw]** `richM.py` recomputes Part M(c) from raw key `M`, giving the same values as md M(c):
  - γ = 10⁻²: Richardson −0.00633 ± 0.00065, +0.00538 ± 0.00059; ratio 0.788, 0.814; z = +2.59, −2.08; step effect −0.00730, +0.01071.
  - γ = 3·10⁻³: ratio 0.794, 0.859.
  - Linear fits over the three times give the same rates.
- **[own run]** `pytest tests/test_gauge_noise.py`: 40 passed in 83 s, including 8 Part M tests.
- **[own run]** `numcheck14.py`: 1409 tokens, 7 not found (1.54 ×2, 0.631, 2.14·10⁻⁴, 2.25, 7.32·10⁻⁴, 8.39·10⁻⁴), each derived in the text, as the README states.
- **[own run] Build:** pdflatex, paper and supplement alternately, 3 times each. Result 30 + 26 pp., 0 errors, 0 overfull boxes, no undefined references, no "??". `pdftotext` is identical to the committed PDFs.

**Recommendation: minor revision.** D3ex-3 and D3ex-4 are text changes; a Δt/4 chain is optional. Accept once they are made, with no further review round. The O(γ²τ) result is correctly derived and reproduced, and its remaining shortfall is honestly marked open.

## Issues

| id | paper | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|---|
| D3ex-1 | 16 | minor | Remark 4.7's QP drops the r_j² weights of the free dual terms (Gε²r_j²/2m)(Φᵀz)_j². It is correct only for equal masses; S3 uses r ≡ 1, so no number changes [own] | l. 361: minimise Σ_free r_j²(Φᵀζ)_j²; the constant (G/m)Σ r_j²\|Φᵀζ*\|_j; README D3em-8 derivation likewise |
| D3ex-2 | 16 | cosmetic | "exponent we do not compute" at \|S\| < N is now superseded for mixed starts (ε^{2−o(1)}, md G) | l. 361 and §8(2): scope the sentence to all-branch starts |
| D3ex-3 | 14 | minor | Prop. 7.4(a) "every start in P_B": for levels (2, 1, 0.5) the ascent points out of P_B on W₃ = 0.2358, W₁ ∈ [0.2947, 0.317] [own `pb_face.py`]. The wall dynamics there are ϖ ≠ 0, outside Derived st. 7.2 | l. 520 and §S1.20: restrict to paths at a distance from the walls, or invoke Prop. 6.3(c), (d) explicitly |
| D3ex-4 | 14 | minor | At γ = 3·10⁻³ the step effect is 4.3× the drift. The two-step Richardson remainder is about 1% of it (M(d) closed forms), about 4% of the drift, which exceeds the quoted γ-factor errors. The second class's 1.61 ± 0.04 against 1.52 is not mentioned | l. 603, S2.4 l. 296: quote the systematic and the second-class excess, or add Δt/4 |
| D3ex-5 | 14 | cosmetic | "has the closed form's size" (l. 515) against "a shortfall as at fixed masses" (l. 603) | Use one wording with the z values |

## Sources
- **[read]**
  - `paper16.tex` (§1.1, §1.3, §2, §4 in full, §7.4–7.5, §8), `supplement.tex` S3 and S6 header, README "Revision after round 144";
  - the 6dfcafc diff, `exp_mirror_branch.py` (`part_G`, `mixed_qp`, `main`);
  - `paper14.tex` (§7 in full, §9.5, Table 1, Remark 4.3, Prop. 6.3, the title), `supplement.tex` (§S1.19–S1.21, §S2.4 and Table S5, the coloured-noise paragraph);
  - README "Revision after round 153" and "The O(γ²τ) drift";
  - `RESEARCH_note_gamma2_drift.md`, the 9d08c80..b83e915 diff, `exp_gauge_noise.py` (`main`, `m_loop` and the rate summary);
  - `REVIEW_round144_D3em.md` and `REVIEW_round153_D3ev.md`.
- **[md]** `mirror_branch.md` C and G; `gauge_noise.md` J, K, L and M.
- **[raw]** `mirror_branch_raw.npz`; `gauge_noise_raw.npz`, key M.
- **[own]**
  - Scripts in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref155/own/`: `pb_face.py`, `qp.py`, `richM.py`, `peekM.py` and `tests_gn.out`.
  - The regenerated md files in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref155/mb_out/` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref155/gn_out/`.
  - The builds in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref155/build/`.
- The first version of `pb_flow.py` (a slow time-stepping loop, my own background task) was stopped with TaskStop and replaced by the face test.

## Lesson
When a leading-order statement is restricted to the interior of a constraint set, its convergence claim "from every start in the set" also needs the path to stay in that interior. A check of the flow's direction on the set's faces settles this in seconds. Likewise, a two-step Richardson extrapolation is only as good as the linearity of the step effect it removes. When that effect is several times the signal, quote the remainder, measured with a third step or bounded by a closed form, before quoting percent-level ratios.
