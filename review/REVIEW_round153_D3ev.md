# Referee report D3ev, round 153: the Paper 14 manuscript, third review

## Scope

This is the third review of `manuscript/paper14/paper14.tex` and `supplement.tex` at commit f498c54 (branch claude/vigilant-clarke-hf939b), as revised after round 150. I also read the two PDFs, `README.md` ("Revision after round 150" and "General-noise section") and `numcheck14.py`. I read these against `review/REVIEW_round150_D3es.md` and the README's account of the revision.

The main new material, refereed here for the first time:
- **Section 7**, "White noise of any covariance in the unconstrained chart" (ll. 470–531, pp. 19–21): Props. 7.1, 7.3 and 7.4, Derived st. 7.2 and Remark 7.5;
- **§9.5** (ll. 602–607, p. 27);
- **the supplement's** S1.17–S1.21 (ll. 134–163) and S2.4 (ll. 275–300, Table S5);
- **the data:** `results/chaos/gauge_noise.md` (Parts A–I) with its raw file, and the code behind it (`exp_gauge_noise.py`, `src/core/gauge_noise.py`).

I worked on an export at `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref153/repo`, with the venv symlinked. The repository was not modified. Line numbers are those of `paper14.tex` unless marked S (`supplement.tex`). Venue: SIADS or MMS.

## Section 7 (verdicts per result)

**Prop. 7.1 (the generator identity): correct.**
- **(a) Existence.** I re-derived LΛ for Λ = e^{ε|Θ|²}: Λ(−2εΘ·G∇(J+V) + ε trΣ + 2ε²ΘᵀΣΘ).
  - G∇V = σ_B(Θ − Π_BΘ) holds because G is diagonal and the box is a product, so the G⁻¹-projection is the coordinatewise clip.
  - −G∇J is bounded: the saturating maps have bounded derivatives, s ≥ s_min, ρ has compact support and the cohesion's force is bounded.
  - So LΛ ≤ −cΛ + C″ for ε < σ_B/λ_max(Σ), and existence follows by Krylov–Bogoliubov, with ∫Λ ≤ C″/c for every stationary law.
  - Uniqueness for Σ ≻ 0 follows from ellipticity and a locally Lipschitz drift. The README is right to drop the source's "hypoelliptic", which alone does not give uniqueness.
- **(b) The identity.** Υ′ = (a_max/2G_α)sinh(2α/a_max) = βϱ/G_α. The gauge identity u_k = −ϱ_k∂_{α_k}J_fit holds at every state, and LΥ/β_k = u_k − ϖ_k + τ_kϱ′_k. Υ, Υ′ and Υ″ grow like e^{2|α|/a_max}, and ϖ like |α|e^{2|α|/a_max}; all are integrable under Gaussian moments, so Itô's formula and the martingale property hold. (b) is exact.
- **(c) The box term's sign.** ∂_αV has the sign of α outside 𝓑, and ϱ_k has the sign of a_k, which is the sign of α_k. Correct.
- **(d) The algebra.** ϱ′ = (a′² − aa″)/(βa′²) with a″ = −2aa′/a_max² gives (1 + s²)/(β(1 − s²)), s = a/a_max. At the fit this is m⁻¹∂_W log(W − q/W), so −τ_kϱ′_k = −∂_{β_k}𝓕_τ. Correct.
- **Consistency with Derived st. 4.4.** With τ_k = τ, 𝓕_τ = (τ/2)ℓ, and Derived st. 7.2(b)'s flow equals (4.1): −γW(τ/(2m))(∂_Wℓ − …).

**Derived st. 7.2: the statuses are appropriate.**
- (a) "formal" is, if anything, conservative. Given the exact identity, the leading order needs only E[ϱ′_k] → ϱ′_k(fit), together with exit terms (quasi-stationary law) that are not quantified.
- The hypotheses are right: compact parts, times short against the exit, a large Péclet number, and unique ergodicity for a degenerate Σ.
- One gap in scope: with the box, "ϖ_k = 0 in the well" (l. 498) fails near the walls of P_B (D3ev-4).

**Prop. 7.3 (the freezing): correct as an identity; the statement is slightly imprecise.**
- The pathwise identity, the masses' equation (W_k/β_k = 1/m), the integration by parts without an Itô term, and D = 0 at every γ all check. For D = 0: Tᵀgξ = Σ_k t_{α,k}ξ_{α_k}/G_α, because the gauges have δc = δs = 0 and g's position block is diagonal.
- Υ_k's integrability is assumed in (b), but the proof (S l. 149) needs dΥ_k/dt = β_k(u_k − ϖ_k) integrable as well. ϖ_k grows faster than Υ_k, so the stated hypothesis is too weak.
- For non-Gaussian (jump) noise on the centres, α_k is only absolutely continuous, not C¹. The identity holds almost everywhere, so only the wording changes.
- (a)'s second term, (γ/m)∫Υ_Σ dW_k, is not bounded in time. Its stationary mean per unit time is (γ/m)²·(an antisymmetric correlation E[Υ_lΥ̇_k], l ≠ k), which is nonzero for a non-reversible positions' law. That term is the O(γ²τ) drift of Part C. "Freezes … at any γ" should say: no diffusion at any γ, and a drift O(γ²τ) instead of O(γτ) (D3ev-3).

**Prop. 7.4 (the correction at the classes' levels): correct; four precision points (D3ev-5).**
- (a) On 𝓜 the masses form gives κ_F∂_β𝓕_τ, and LaSalle gives convergence to W♯_τ. I recomputed W♯_τ for all four level vectors by root-finding on the Lagrange condition.
- (b) The levels max(τ_k, τ_B) keep the net weights κ_F·max(τ_k, τ_B) − τ_k positive only for κ_F > 1, which (b) does not state. At κ_F = 1 the weights of the classes with τ_k ≥ τ_B vanish.
- Also in (b): the target is the maximiser for those weights, which is neither W♯_τ nor W♯ in general, and it must lie in P_B. Further, (b) allows τ_k = 0, which contradicts the proposition's head hypothesis "every τ_k > 0".
- (c) The plug-in cancellation of the noise's term is exact in every stationary law, wherever the floor does not bind. Correct.
- (d) "Without bias from Δα − bΔt" holds for the Euler–Maruyama chain. For a sampled continuous diffusion an O(Δt²) cross term remains. The proof is about the chain, but the statement does not say so.

**Remark 7.5: all three formulas are correct.**
- **(i) Coloured noise.**
  - My own linear-response derivation: E[δΘ(t)η(t)ᵀ] = ∫₀^∞ e^{−GHs}e^{−s/τ_c}ds Σ/(2τ_c) = ½(I + τ_cGH)⁻¹Σ. Since E[η] = 0, E[ϱη_α] ≈ ϱ′E[δα η_α], so the level is [(I + τ_cGH)⁻¹Σ]_{αα}/(2G_α), with the white limit at τ_c → 0. The next order is O(τ²): the odd Gaussian moments vanish.
  - I built H independently: f̂'s gradient, ρ by 8192-point quadrature of 2 sin t, and the instance's maps.
  - By two routes (G^{1/2}-symmetrised eigen-decomposition, and quadrature of the time integral) I get 0.4686, 0.3667, 0.2806 at τ_c = 100, and also 0.9854/0.9806/0.9691 (τ_c = 1), 0.8779/0.8395/0.7675 (10) and 0.1070/0.0638/0.0508 (1000). These are exactly Part H's values.
  - The rates of GH run from 3.059·10⁻⁴ to 0.1952, as in the md.
- **(ii) Stratonovich noise.** The added Stratonovich drift c_α enters LΥ/β as ϱc_α/G_α. Its sign and its order τ are correct.
- **(iii) The physical chart.** ψ = a²/(2βG_a) gives E[u_k] = −τ_k/β_k plus boundary terms, so the selection descends Σ_kτ_k log W_k. This is consistent with Remark 4.3(ii).
- **A strengthening the authors can take (D3ev-7).** The identity works in any chart α̃ = ψ(a):
  - u_k = −(a_kψ′(a_k)/β_k)∂_{α̃}J_fit;
  - the generator gives E[u_k] = −(τ/β_k)E[1 + a_kψ″/ψ′] + ϖ-terms;
  - at the fit this equals −(τ/2)∂_{β_k}ℓ_ψ, exactly Remark 4.3's ℓ_ψ.
  
  So Remark 4.3's "selection in another chart" no longer needs the Laplace step.

**Consistency with Remark 4.3 and Derived st. 4.4.** Consistent throughout. Prop. 7.1(d) (the centres' and widths' noise does not enter) agrees with Remark 4.3's "the chart of the centres and widths does not enter".

## Status of the round-150 issues

| D3es id | resolved? | evidence | what remains |
|---|---|---|---|
| 1 | yes | §9.4 l. 598: (1.04 ± 0.07)ΔJ_𝒪 at ς = ½, 0.91–1.17 over ς, ς = 0.36 ± 0.25, ς = 1 at 2.57 s.e., χ² = 2.56 on 5 dof; Table 1 l. 106 | none |
| 2 | partly | abstract l. 64 now gives −0.045 ± 0.133, 0.652 ± 0.242 and "not a sharp test" | "consistent with it" refers to δ^{κ_F−1}, whose exponents are 0 and 1, while −0.093 and 0.823 are the layer law's (D3ev-8) |
| 3 | yes | "Its divergence there, not its direction, is the dithered chart's Jacobian" | none |
| 4 | yes | Prop. 6.2 l. 440 (start point); S1.15 gradient-path proof (checked: ½∫\|ẏ + ∇J\|² ≥ 2∫∇J·ẏ); second minimum 6.397τ; exit points (S l. 352) | none |
| 5 | yes | l. 587: jackknife over shared seeds, 0.697 ± 0.214, z = 3.25 | none |
| 6 | yes | −0.062 and +0.776 as a second reference; ticks in Fig. 3(a) | none |
| 7 | yes | S2.6: 0.253, 0.123, 0.203 (κ_F = 1); 0.179, 0.087, 0.143 (κ_F = 2) | none |
| 8 | yes | SSP Cor. 5.7 at l. 124; §6.2 in Open problem (1); DLLN19 Prop. 13 at ll. 596/598 | DLLN19's [H1]–[H3] (∂_nJ > 0 on ∂𝒪) not checked, as the README says |
| 9 | yes | S2.6: first exits 2135.5–3107.25; lag fits 0.931–0.946 | none |
| 10 | yes | 3.74·10⁶–1.27·10⁷; 6.91–9.66·10⁵; 0.65–1.02; γ ≲ 10⁻⁴ | none |
| 11 | yes | the abstract is 249 words by my count (math as one word) | at the limit; any addition must be offset |
| 12 | mostly | 𝓑, ς, N_e, t_{α,k} (S l. 63), numbered Remark 4.3, Fig. 3 ticks | "vanishes inside $B$" (l. 448) still uses B for the box; new clashes from Section 7 (D3ev-11) |

## New findings

### Major

None.

### Minor

1. **Part D's χ² test is tautological with respect to the identity (§9.5 l. 603; S2.4 "The identity").**
   - The code defines lpsi = u − g·nonfit + τ_k·gp (the identity's right side) and bias = (E_nΥ_{n+1} − Υ_n)/Δt − lpsi (`exp_gauge_noise.py` ll. 210–214). So E[LΥ] + E[𝓔] ≡ (Υ_end − Υ_start)/T − (the martingale's mean), whatever lpsi is.
   - Hence the reported "χ² = 24.8 for 18 values" only tests that a sum of martingale increments averages to zero. It cannot detect an error in the identity. The supplement half-says this ("by construction"); the paper presents it as the test.
   - Likewise, Part I's "the gauge identity splits each run's log-masses exactly (largest residual 3.4·10⁻¹⁴)" (S2.4) is exact because the chain's bias is defined as the remainder (l. 637).
   - The real evidence is elsewhere:
     - the per-state bias is 0.17–1.9% of τ_kE[ϱ′_k], i.e. O(Δt), but only one Δt was run;
     - the size of the box's term;
     - the smallness of Part I's remainder (≤ 4.5·10⁻⁵ in W₁).

2. **The Euler chain's O(Δt) bias is not carried into the long-run comparisons (S2.4, Table S5).**
   - In Part I the chain's bias moves W₁ by −4.4·10⁻⁵ in 2·10⁴ time units. Linearly extrapolated to Part E's second half (mean time 7.5·10⁵), that is about −1.6·10⁻³, comparable to the s.e. of the noisy variants (0.0018 in W₁).
   - The "centres, widths" row has |z| = 21.8, and the "+κ_F = 1" row 8.3, at s.e. < 10⁻⁴. S2.4 calls these "deviations of that order".
   - My decomposition of Part E's "nonlinear only" records (`freeze.py`), second-half means against 1/3:

     | component | W₁ | W₂ | W₃ |
     |---|---|---|---|
     | observed | −2.0·10⁻⁴ | +2.4·10⁻⁴ | −0.4·10⁻⁴ |
     | Prop. 7.3(a)'s bounded term | −0.4·10⁻⁴ | +0.8·10⁻⁴ | −0.4·10⁻⁴ |
     | Katzenberger's drift (Part C) × time | −0.6·10⁻⁴ | +0.5·10⁻⁴ | +0.1·10⁻⁴ |
     | Part D's centres/widths chain bias, scaled ∝ τ to 10⁻⁴ | about −1.0·10⁻⁴ | +1.8·10⁻⁴ | −0.8·10⁻⁴ |

   - The bounded term tracks the per-record fluctuations with correlation 0.997.
   - The scaled chain bias has the right sign and order to account for the rest. This is my estimate; Part D's τ = 10⁻² run is contaminated by the box (class 2 outside it 37% of the time).
   - Accumulating the chain's bias along Part E's runs, as Part I already does, or halving Δt for the two deterministic-amplitude variants would turn "deviations of that order" into an explained term.

3. **Prop. 7.3's statement (ll. 507–512).**
   - (b) needs ϖ_k (equivalently dΥ_k/dt), not only Υ_k, to be integrable.
   - Replace "C¹" by "absolutely continuous" for jump noise.
   - (a)'s second term carries the O(γ²τ) secular drift, so "freezes … at any γ" means no diffusion at any γ and no drift at order γτ.
   - The abstract's "the gauge freezes" deserves the same qualifier.

4. **"ϖ_k = 0 in the well" (l. 498) fails near the walls of P_B.** Part I's fixed-mass table gives the box's term as a multiple of the noise's drive:

   | W₁ | 0.2947 (the wall) | 0.3050 | 0.3100 | 0.3125 | 0.3150 | 0.3200 |
   |---|---|---|---|---|---|---|
   | box term / noise drive | 5.650 | 1.090 | 0.237 | 0.122 | 0.014 | 0 |

   - W_{B,1} = 0.2947 is only 0.039 from 1/3, so the layer of about 0.02 is half of that distance.
   - With the box, Derived st. 7.2 and Prop. 7.4(a) need 𝒦 compact in the interior of P_B, at a distance large against the positions' spread.
   - On the instance W♯_{τ,3} = 0.2567 is only 0.021 from W_{B,3} = 0.2358.

5. **Prop. 7.4 (ll. 517–525).**
   - (b) omits κ_F > 1.
   - (b)'s target is unnamed: it is the maximiser for the weights κ_Fmax(τ_k, τ_B) − τ_k, and it must lie in P_B.
   - (b) breaks the head hypothesis "every τ_k > 0"; move that hypothesis into (a).
   - (d) is a statement about the Euler chain; say so.
   - Table 1's "(a), (b) proved at leading order" (l. 112) should read "the flow's convergence proved; the leading-order fitness as in Derived st. 7.2(a) (formal)".

6. **"Analytic centre" is a misnomer (l. 78, l. 479; l. 128 for W♯).** log(W_k − W²_min,k/W_k) is not the log of an affine constraint. My comparison (`centre.py`):

   | levels τ_k/τ | W♯_τ (barrier maximiser) | τ-weighted analytic centre of 𝒲 |
   |---|---|---|
   | (2, 1, 0.5) | (0.4609, 0.2824, 0.2567) | (0.4643, 0.2821, 0.2536) |
   | (1, 1, 1) | (0.3765, 0.2988, 0.3247) | (0.375, 0.300, 0.325) |

   Call W♯_τ "the maximiser of the weighted barrier (the D-optimal fit for equal levels)".

7. **Remark 4.3 can be upgraded with Section 7's identity** (see Remark 7.5 above). One line in Remark 7.5(iii), and Table 1's row for Remark 4.3 then reads "the selection in another chart: as Derived st. 7.2(a)".

8. **Abstract and title (l. 64, l. 53).**
   - "Consistent with it" mixes δ^{κ_F−1} with the layer law's numbers (D3es-2 remainder).
   - "Several particles per class obey the same laws, slowed by K/m" holds for the isotropic dither only (Remark 8.3). After Section 7 it reads as covering any noise.
   - "The gauge freezes": see finding 3.
   - The title is a four-item list. "The dither's confinement" misleads: §6 shows the dithered positions have no stationary law, and the confinement belongs to the correction and to the dead zone.

9. **The coloured-noise simulation (l. 607).** "0.496, 0.454, 0.237 (linear response: 0.469, 0.367, 0.281)" carries no errors.
   - The delta-method errors (independent runs) are about 0.12, 0.08 and 0.13, so the agreement is 0.2, 1.1 and 0.4σ [own].
   - The quoted +3.2/+5.1/+4.5 test only that the coloured fitness differs from the white, not the agreement.

10. **Small number issues.**
    - "0.0220 from it, of which 0.0103" (l. 605) compares class 1's distance with class 3's shift. From the raw, the zero at γ = 10⁻² is (0.45083, 0.28211, 0.26706), so class 1's shift is 0.0100 of its 0.0220.
    - "Below 1.1·10⁻²²" is the maximum ("at most"). "Reads 0.0005 of τ_k" is also the largest of 0.0005/0.0004/0.0003.
    - Part I's regressed s.e. omits b̂'s uncertainty (the m̄²/S_xx term). That adds about 9% for the 36 seeds at 2·10⁴ and is negligible for the 192 seeds, which the paper uses.

11. **Notation.**
    - "Vanishes inside $B$" (l. 448).
    - ε in Prop. 7.1(a) clashes with the leak strength ε of Table 2.
    - c for the Stratonovich drift (Remark 7.5(ii)) clashes with the centre c.
    - The subscript in Υ_Σ clashes with the covariance Σ.

12. **Literature (§1.3 l. 132).**
    - In Fudenberg–Harris and Imhof, "each strategy's own payoff noise" holds for independent shocks; with correlated shocks the covariances enter. The sharper contrast is that here correlations cancel exactly at leading order (Prop. 7.1(d)).
    - Optional additions [memory]: Volpe–Wehr (Rep. Prog. Phys. 79 (2016) 053901) for convention- and chart-induced drifts (Remark 4.3, Remark 7.5(ii)–(iii)), and Hänggi–Jung (Adv. Chem. Phys. 89 (1995)) for the filtered response to coloured noise.
    - The citations to configurational temperature (Rugh; Jepps–Ayton–Evans), the blue-moon mean force with its Jacobian term (Carter–Ciccotti–Hynes–Kapral; den Otter), Fudenberg–Harris, Imhof, Hofbauer–Imhof and SSP are apt [memory].

13. **Length.**
    - 31 + 22 pp. The body ends on p. 28 (Section 7 on pp. 19–21, §9.5 on p. 27).
    - The general-noise material is a genuine gain: it explains why the dither selects, as a generator identity. But the paper now carries two stories.
    - Suggested moves to the supplement: Prop. 7.4(c), (d), Remark 7.5(ii) and §9.5's Part I narrative. Shorten the title.

## Verification I did

All marked [own]; scripts in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref153/own/`.

- **Regeneration.** I copied `gauge_noise_raw.npz` into `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref153/gn/`. `GN_OUT=<that dir> exp_gauge_noise.py --from-raw` ran in 0.8 s. `cmp` with `results/chaos/gauge_noise.md` is byte-identical.
- **Part I from the raw** (`partI.py`, my own z, two-sided p, Holm step-down, and least-squares regression; no repository functions used):
  - Holm-adjusted p for the 192 seeds: 0.3108, 0.2151, 0.4964 (minimum 0.2151); for the 36 seeds: 0.0363; for the 12 seeds: 0.0704.
  - Martingale-regressed W₁ for the 192 seeds: 0.3299 ± 0.0002, 0.3282 ± 0.0005, 0.3265 ± 0.0009, i.e. 0.90 ± 0.06, 0.84 ± 0.08 and 1.06 ± 0.14 of Katzenberger's displacement. The correlations are 0.978, 0.925 and 0.801, and the martingale's contribution +0.0050 ± 0.0034 at 2·10⁴. All match the md.
  - Note: the 192 seeds' regressed ratios lie 1.7σ and 2σ below 1 at 5000 and 10⁴. Correlated over time, they still hint at a slower early displacement of about 10%.
- **Coloured level** (`coloured.py`): my own H and two routes; Part H's level ratios are reproduced to four digits at τ_c = 1, 10, 100 and 1000.
- **λ_min(H)** (`lmin.py`): 0.00120049 at W = 1/3 and 0.00120001 at W♯, so the "Laplace parameter 0.083 at W♯" (l. 200) stands.
- **W♯_τ against the analytic centre** (`centre.py`): finding 6.
- **Prop. 7.3(a) on Part E's records** (`freeze.py`): finding 2.
- **Part F's zeros from the raw** (`peekF.py`): finding 10.
- **Section 7 by hand:**
  - Υ′ = βϱ/G_α and LΥ/β;
  - ϱ′ and its value at the fit;
  - the Lyapunov inequality;
  - the box's sign;
  - the freezing identity and its integration by parts;
  - D = 0;
  - the max(τ_k, τ_B) weights;
  - the quadratic-variation bias (Euler chain);
  - the coloured linear response;
  - the Stratonovich term;
  - the physical chart, and its extension to any chart ψ.
- **Code reading:** `exp_gauge_noise.py` Part D (ll. 180–228), Part I (ll. 600–640, 1020–1170), `src/core/gauge_noise.py`.
- **Tests:** `tests/test_gauge_noise.py`, 24 passed (60 s).
- **Build:** three alternating pdflatex passes in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref153/build/` give 31 + 22 pp., 0 warnings, no undefined or multiply defined references, no overfull boxes, no "??".
- **`numcheck14.py`:** 1277 tokens, 7 not found (1.54, 0.631, 2.14·10⁻⁴, 2.25, 7.32·10⁻⁴, 8.39·10⁻⁴), as the README lists.
- **§9.5 against the md:** every decimal checked against Parts A–I. All are printed there (with the rounding and component remarks of finding 10).
- **Abstract:** 249 words (`wc.py`).

## Overall assessment

**Recommendation: minor revision** (accept after minor changes), for SIADS or MMS.

Section 7 is mathematically sound and is the paper's best addition. The generator identity is exact and nicely proved. It explains the dither's selection without Laplace's method, isolates each amplitude's own level, and gives a clean freezing result and a correction at the classes' levels. Every formula I re-derived holds, including Remark 7.5(i)'s filtered level, reproduced from my own Hessian. The md regenerates byte for byte, and Part I's statistics recompute from the raw. All round-150 issues are resolved except a half-fixed abstract sentence and one notation residue.

What remains is presentation:
- the statement-level precision of Props. 7.3 and 7.4;
- the box's boundary layer, where ϖ ≠ 0 in the well;
- Part D's tautological χ² test;
- the Euler chain's bias, which is unaccounted for in the long-run comparisons and the likely source of the |z| = 21.8;
- "analytic centre";
- the abstract's last sentences and the title;
- length.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ev-1 | minor | Part D's χ² = 24.8 and Part I's "exact split (3.4·10⁻¹⁴)" hold by construction (bias := remainder; code ll. 210–214, 637); they cannot test the identity | §9.5 l. 603, S2.4: present the per-state bias (0.17–1.9% of τ_kE[ϱ′_k]) as the check, with a Δt-halving run for its O(Δt) scaling; label the χ² as a check of the averaging |
| D3ev-2 | minor | Euler bias not in Part E/F. Part I's −4.4·10⁻⁵ by 2·10⁴ extrapolates to about 1.6·10⁻³; the "nonlinear only" offset (−2.0, +2.4, −0.4)·10⁻⁴ (\|z\| 21.8) is bounded term (−0.4, +0.8, −0.4) + Katzenberger (−0.6, +0.5, +0.1) + scaled chain bias of about (−1.0, +1.8, −0.8), all ·10⁻⁴ [own] | S2.4, Table S5: accumulate the chain's bias along Part E (as Part I does) or halve Δt; replace "deviations of that order" |
| D3ev-3 | minor | Prop. 7.3: (b) needs ϖ_k / dΥ_k/dt integrable (the proof uses it); "C¹" is false for jump noise; (a)'s integral carries the O(γ²τ) drift (antisymmetric correlation (γ/m)²E[Υ_lΥ̇_k]) | ll. 507–512; abstract "the gauge freezes": "at leading order (no diffusion at any γ; drift O(γ²τ))" |
| D3ev-4 | minor | "ϖ_k = 0 in the well" (l. 498) fails within about 0.02 of P_B's walls: box term / noise drive = 5.65, 1.09, 0.24, 0.01 at W₁ = 0.2947, 0.305, 0.31, 0.315 (md I) | l. 498, Derived st. 7.2, Prop. 7.4(a): 𝒦 compact in the interior of P_B, at a distance large against the spread; quote the layer |
| D3ev-5 | minor | Prop. 7.4: κ_F > 1 missing in (b); (b)'s target unnamed and its membership of P_B not required; (b) contradicts "every τ_k > 0"; (d) holds for the Euler chain only; Table 1 "proved at leading order" overstates | ll. 517–525; Table 1 l. 112 |
| D3ev-6 | minor | "Analytic centre" for W♯_τ (and W♯): the barrier log(W − q/W) is not affine-log; the weighted analytic centre is (0.4643, 0.2821, 0.2536) against (0.4609, 0.2824, 0.2567) [own] | ll. 78, 128, 479: "weighted-barrier maximiser" |
| D3ev-7 | minor | The identity extends to any chart ψ: E[u_k] = −(τ/β_k)E[1 + aψ″/ψ′] + ϖ-terms, which equals −(τ/2)∂_{β_k}ℓ_ψ at the fit | Remark 7.5(iii) one line; Remark 4.3 and Table 1 row |
| D3ev-8 | minor | Abstract: "consistent with it" against the layer law's −0.093, 0.823 (D3es-2 remainder); "Several particles per class obey the same laws" holds for the isotropic dither only; title long, and "the dither's confinement" misleads | l. 64; l. 53 |
| D3ev-9 | minor | Coloured simulation: ratios without errors; delta-method errors about 0.12, 0.08, 0.13, so 0.2/1.1/0.4σ from linear response [own] | l. 607, S2.4 |
| D3ev-10 | cosmetic | 0.0220 (class 1) "of which 0.0103" (class 3; class 1's shift is 0.0100 from raw); "below 1.1·10⁻²²" and "0.0005" are maxima; Part I's regressed s.e. omits b̂'s uncertainty (+9% at 36 seeds) | l. 605, l. 603, l. 607; S2.4 |
| D3ev-11 | cosmetic | "vanishes inside $B$"; ε (Prop. 7.1) against the leak strength ε; c (Stratonovich drift) against the centre c; Υ_Σ against Σ | l. 448, l. 484, l. 530, l. 507; Table 2 |
| D3ev-12 | minor | §1.3: FH/Imhof's "own variance" holds for independent shocks only; state the contrast (correlations cancel here); optionally Volpe–Wehr 2016 and Hänggi–Jung 1995 | l. 132, l. 126; Remark 7.5 |
| D3ev-13 | minor | 31 + 22 pp., body to p. 28 | Prop. 7.4(c), (d), Remark 7.5(ii) and §9.5's Part I paragraph to the supplement; shorter title |

## Sources

- **[read]**
  - `paper14.tex` (all 742 lines);
  - `supplement.tex`: S1.5–S1.7, S1.15, S1.16–S1.21 and S2.4 in full, and S2.6 (exits, caps);
  - README "Number check", "Revision after round 150" and "General-noise section";
  - `review/REVIEW_round150_D3es.md`;
  - `exp_gauge_noise.py` (header, Part D ll. 180–228, Part I ll. 600–640 and 1020–1170);
  - `src/core/gauge_noise.py`;
  - `exp_cohesion_leak.py` (instance).
- **[md]/[raw]** `gauge_noise.md` A–I; `gauge_noise_raw.npz` (keys I, E_rec, F).
- **[own]** `keys.py`, `peekI.py`, `peekE.py`, `peekF.py`, `partI.py`, `freeze.py`, `coloured.py`, `lmin.py`, `centre.py`, `wc.py` and `tests.out`, in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref153/own/`; the regenerated md in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref153/gn/gauge_noise.md`; the build in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref153/build/`.
- **[memory]** Rugh 1997; Jepps–Ayton–Evans 2000; Carter–Ciccotti–Hynes–Kapral 1989; den Otter 2000; Fudenberg–Harris 1992; Imhof 2005; Hofbauer–Imhof 2009; Volpe–Wehr 2016; Hänggi–Jung 1995. No downloads this round.

## Lesson

When a check defines one term as "what is left over", the identity it then verifies is exact by construction. Its χ² or its 10⁻¹⁴ residual only tests the bookkeeping. The evidence lies in the remainder itself: its size, and how it scales with the step. Here that remainder, the Euler chain's O(Δt) bias, is small per step but accumulates over 10⁶ time units, to the size of the unexplained |z| = 21.8. Before calling an exact identity "confirmed in simulation", ask which computed quantity would have changed had the identity been false, and track the leftover term in every long-run comparison.
