# Referee report: Corollary E.0′.3 and Lemma E.0″ (PARADIGM.md §6), with the code behind them

Scope: `PARADIGM.md` lines 1622–1686 (Corollary E.0′.3) and 1688–1811 (Lemma E.0″, Remark E.0″.1, Numerics, Checked);
`src/core/excitation.py` lines 522–763 (`dither_covariance_bound`, `exact_stationary_covariance`, `exp_average`,
`ExpAverageLyapunov`, `lyapunov_exponent_rate`); the last two tests of `tests/test_excitation.py`;
`experiments/exp_explicit_rate.py::dither_table` (table 6) and `experiments/exp_pe_general.py` → `results/explicit_rate/general_pe.md`
(the regenerated table of 22:03, read after the run finished). Lemma E.0′ and the (LTV) setting of Lemma E.0 were read first.
Scratch scripts and raw outputs: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_E7/`
(`check_E0p3*.py` → `out_E0p3*.txt`; `check_E0dp_A…E.py` → `out_E0dp_*.txt`). No repository file was modified.

## VERDICTS

| item | verdict |
|---|---|
| Corollary E.0′.3 (statement + proof) | **holds**. Itô correction, block form of 𝐇, martingale step, Grönwall, the three bounds and the explicit ν̄ are all correct; the averaged covariance C̄ = kX, M̄X + XM̄ = ΣΣᵀ has the right sign and factor. One numerical law quoted in its Numerics paragraph is wrong (issue 1). |
| Lemma E.0″ (statement + proof) | **holds**. The exponential-average identities, Ḣ, the order-ε² residual −ε²ζᵀN̄₂ζ (sign checked), every per-direction bound in the time-local norm, the carried-over formulas for c₁, c₂, d_w, r_ζ (no hidden (T/4) bound), the PE equivalence (both directions and the identity) are correct. Independent certificate confirms it; no counterexample in nine regressor families. Its Status/Numerics text has several inconsistencies with its own table (issues 5–9). |
| Code | **holds with fixes** (all minor): a vacuous assertion in the E.0′.3 test, a docstring/code mismatch in `ExpAverageLyapunov.constants(explicit=True)`, cosmetic docstring typos. `exact_stationary_covariance`, `dither_covariance_bound`, `exp_average`, `ExpAverageLyapunov`, `lyapunov_exponent_rate` all agree with independent computations (numbers below). The version of `exp_pe_general.py` I read at the start had a spline-extrapolation bug (Lorenz orbit shorter than the padded grid for T = 20); it was fixed at 22:00:13, before the published table (22:03:40) was generated, and I verified the table's T = 20 row reproduces exactly with the fixed script. |

Both tests pass (`pytest -k "E0prime_3 or E0dprime"`: 2 passed in 35.7 s).

## ISSUES

### Corollary E.0′.3

1. **[minor — wrong formula in the text]** `PARADIGM.md:1679` ("The exact E‖e‖² is ‖Σ‖_F²/k² to four digits on both instances") and the
   same claim in `experiments/exp_explicit_rate.py:284` (table.md reading: "the exact E‖e‖² is of order ‖Σ‖²_F/k²"). The table's own
   numbers refute it: on the sinusoid ‖Σ‖_F² = 4·10⁻³ and E‖e‖² = 1.000·10⁻⁴, 4.000·10⁻⁵, 1.000·10⁻⁵ at k = 20, 50, 200, i.e.
   **‖Σ‖_F²/(2k)**, not ‖Σ‖_F²/k² (= 10⁻⁵ at k = 20, off by a factor k/2 = 10, 25, 100). On the ideal, ‖Σ‖_F² = 15·2·10⁻³ = 0.03 and
   E‖e‖² = 5.000·10⁻⁵, 3.750·10⁻⁵, 2.500·10⁻⁵ = 0.03/(2k) again. The law is general: e ≈ −(1/k)Ψζ on the slow manifold gives
   E‖e‖² ≈ k⁻² tr(M̄ C̄) = k⁻¹ tr(M̄X) = ‖Σ‖_F²/(2k), because tr(M̄X) = ½ tr(M̄X + XM̄) = ½‖Σ‖_F². I confirmed the ratio
   E‖e‖²/(‖Σ‖_F²/(2k)) = 1.00000 in all 14 (instance, k) cases I ran (six regressors, Σ isotropic, diagonal, rank-one, random 2×3,
   non-commuting with M̄). *Fix:* replace "‖Σ‖_F²/k²" by "‖Σ‖_F²/(2k) (= tr(M̄X)/k, the slow manifold e ≈ −εΨζ carrying the
   stationary spread)" in both places; the parenthetical explanation is right, the exponent is not.

2. **[minor — observation on the "averaging is exact to O(k⁻²)" claims]** `PARADIGM.md:1675–1677` and the test
   (`tests/test_excitation.py:318`, "within 2/k²") hold on the instances quoted (I reproduce the deviation exactly: exact/averaged − 1 =
   1.250/k² for every Σ tested on the sinusoid, 3.3·10⁻⁴ on the ideal). But the O(k⁻²) agreement is a property of ΣΣᵀ commuting with
   M̄ (more precisely of C̄ being diagonal in M̄'s eigenbasis), not of the slow regime: with Ψ = (sin t, 0.3 cos t) (κ(M̄) = 11.1) and
   Σ = [[0.05, 0.02],[0, 0.03]] the deviation is +0.236/k, +0.227/k, +0.224/k at k = 50, 200, 1000 (O(1/k)); with Ψ = (sin t, sin(t + 0.4))
   (κ = 24) and Σ = diag(0.05, 0.01) it is −0.247/k at k = 200; with the same anisotropic M̄ and Σ ∝ I it is back to 0.34/k². The
   first-order averaging correction is a commutator [⟨ΨᵀΨ̇⟩, C̄] (traceless), whose Lyapunov-equation response has non-zero trace when
   M̄ is anisotropic. The bound still dominates in every case (ratios 5.13/4.59, 4.30/3.82, 4.13/3.67 for δ = 0.2/0.1 — the trace
   price (q/λ_min)/Σ 1/λ_i ≈ 3.7 times 1/(1 − δ) and the O(1/k) factors; 2.65/2.36 on the κ = 24 regressor at k = 200). *Fix:* qualify
   "essentially exact" (line 1676) by "up to O(k⁻²) when ΣΣᵀ commutes with M̄, O(1/k) in general"; the corollary itself is unaffected.

3. **[cosmetic — proof wording]** `PARADIGM.md:1655–1656`: "the linear SDE has finite moments of every order" is not available under the
   hypothesis "every initial law with finite second moment"; the martingale property only needs E∫₀ᵗ‖𝐇(s)x_s‖² ds < ∞, which follows
   from the second moments (bounded coefficients). *Fix:* "finite second moments, which propagate under the linear SDE".

4. **[minor — vacuous test assertion]** `tests/test_excitation.py:316`: `assert b['bound_w'] >= ex['e_mean'] * 0.5` with the comment
   "E‖w‖² ≥ ½E‖e‖² − ‖e₂‖² terms" compares the bound on E‖w‖² with an unrelated quantity; in fact E‖w‖² ≪ E‖e‖² (exact E‖w‖² =
   tr([I, −L]C[I, −L]ᵀ) averaged over the period = 1.28·10⁻⁷ vs E‖e‖² = 1.0·10⁻⁴ at k = 20; 1.25·10⁻¹⁰ at k = 200), so the assertion
   passes only because the bound (5.2·10⁻³, 4.9·10⁻²) is 4·10⁴–4·10⁸ times larger than either. *Fix:* assert `bound_w ≥ max_t E‖w_t‖²`
   with the exact value from C(t) and L(t) (cheap: `quadratic_forms` gives L); it holds with a factor ≥ 4·10⁴, which is the factor
   1/β the text already calls "coarse".

### Lemma E.0″

5. **[minor — Status contradicts the table]** `PARADIGM.md:1732`: "the conditions hold from k = 20 on" is stated for both regressors.
   For the Lorenz regressor they hold only from k = 100 on (general_pe.md: rows k = 20, 50 fail (A) and/or (C) for all three windows:
   at T = 5, k = 50, δ = 0.2 only "ABD"; at T = 10, 20, k = 50 only "BCD"); the Numerics paragraph itself says "from k = 100 on"
   (line 1797). `PARADIGM.md:1734`: "within a factor 4–5 on the chaotic one" holds only for the best window T = 20 at k = 100–200
   (5.0/4.3, 4.2/3.7; 3.9/3.5 at k = 400); at T = 5 the ratio is 14–138 (c_ζ = 0.07 at k = 100 because δ₂ = 0.69) and at T = 10 it is
   7.5–11. *Fix:* "from k = 20 on (quasi-periodic) and k = 100 on (chaotic)"; "within a factor 4–5 at its best window (T = 20)".

6. **[minor — wrong range]** `PARADIGM.md:1791`: "the explicit closed-form constants satisfy the conditions from k = 50 on and give
   rates 10–25% below the measured-constant ones". From the table (and my recomputation): the deficit is 24/21% (T = 2π), 42/37% (4π),
   79/70% (8π) at k = 50 (δ = 0.2/0.1), 6/5, 10/9, 19/17% at k = 100, 1/1, 2/2, 5/4% at k = 200 and ≤ 1% at k = 400. *Fix:* "between
   1% (k = 400) and 79% (T = 8π, k = 50), the deficit decreasing like 1/k".

7. **[minor — stale numbers]** `PARADIGM.md:1796–1797`: "‖N̄₂‖ = 1.05, 0.30, 0.17" for the Lorenz windows T = 5, 10, 20; the table has
   1.054, 0.341, 0.422 (and my reconstruction of the T = 20 row gives 0.422). Note that ‖N̄₂‖ is not monotone in T here, which the
   Remark's "n₂ ∝ T" (line 1762) and the Numerics' "‖N̄₂‖ (which grow with T)" do not reflect for this regressor. *Fix:* quote 1.05, 0.34,
   0.42 and soften "grow with T" to "typically grow with T".

8. **[minor — wrong number]** `PARADIGM.md:1789`: "(T = 2π: 3.7 → 1.9 and 3.1 → 1.9 from k = 20 to 200.)" The table gives 3.69 → 2.16
   (δ = 0.2) and 3.10 → 1.92 (δ = 0.1) at k = 200 (2.13/1.89 at k = 400). *Fix:* "3.7 → 2.2 and 3.1 → 1.9".

9. **[cosmetic]** `PARADIGM.md:1776–1777`: "the exact decay rate is λ_∞/k = 0.5/k to 0.5%": at k = 20 it is 0.02485 vs 0.02500, i.e. 0.6%
   (0.5%, 0.6%, 0.4%, 0.8% at k = 50, 100, 200, 400 from the table's 5-digit values). *Fix:* "to 1%".

10. **[cosmetic — statement wording]** `PARADIGM.md:1699`: "(any; the long-time means of M̄ and F̄ minimise the norms below)". The mean
    minimises the mean-square deviation, not the sup-norms p₁ = 2T sup‖M̄ − M̂‖, p₂, π₁, π₂ (the Chebyshev centre does). *Fix:* "a natural
    choice; the means approximately minimise…" or drop the clause. The proof only uses that M̂, F̂ are constant and symmetric
    (constant for Ṗ₁ = 2(M − M̄), Ṡ = F − F̄; symmetric so that H is symmetric) — nothing more is needed, and nothing more is assumed.

11. **[minor — finite horizon of the numerics, to be stated]** The lemma's constants are suprema/infima over t ∈ [0, ∞); the code takes
    them over [0, 300] (`ExpAverageLyapunov`, `N_valid`) and the certificate is checked there. For Ψ = (sin t, cos √2 t) this is harmless:
    with closed-form exponential averages, inf λ_min(M̄_T) = 0.39068, sup λ_max = 0.61572, ρ = 4.0000, ρ_DD = 10.000, p₁ = 2.908, π₁ = 4.396
    are the same to 4–5 digits on [0, 300], [0, 10³], [0, 3·10³] and [0, 10⁴] (the orbit is dense on the torus within 300 units). For the
    Lorenz regressor it is not exact: on [0, 1500] (same orbit, T = 20) λ_T stays 0.0506 (the minimum is at t = 120), but p₁ 5.31 → 5.67,
    p₂ 19.4 → 20.7, π₁ 23.5 → 25.2, π₂ 49.4 → 51.9, ‖N̄₂‖ 0.422 → 0.481, δ₂(k = 100) 0.083 → 0.095, and the guaranteed rate at k = 100 drops
    3% (3.049·10⁻⁴ → 2.960·10⁻⁴); the certificate still holds (sup V̇/V = −9.72·10⁻⁴ ≤ −γ = −5.92·10⁻⁴). The finite-time "exact" rate
    at k = 100 is 0.00152 over [0, 200] and 0.00156 over [0, 1000]. *Fix:* say in the Numerics paragraph that the Lorenz constants are
    [0, 300] estimates (a few per cent below the longer-horizon values) and that the certificate is a statement on [0, 300].

12. **[cosmetic]** `PARADIGM.md:1803`: "The regressor traversed at standard Lorenz speed has ρ_DD of order 10⁶–10⁷" refers to the earlier
    experiment (whose table I saw before the regeneration) that used *different components*, Ψ = (tanh(x₁/8), tanh(x₂/8), (x₃ − 25)/20),
    with ρ_DD = 1.5–2.1·10⁷; the sentence suggests the same components at speed 1, which was not computed. Either rerun or say so.

### Code

13. **[minor — docstring/code mismatch]** `src/core/excitation.py:690` says `explicit=True` uses "p₂ = ‖M − M̄‖ + T‖F̄ − F̂‖" (the
    statement's bound), but line 698 computes `p2 = PM**2 + T*(sup‖F̄‖ + ‖F̂‖)`, which is a looser (still valid) bound. Also the explicit
    branch keeps the *measured* `p1 = p1_meas` (correctly: p₁ := sup‖P₁‖ has no closed form) — the docstring should say so.

14. **[cosmetic]** `src/core/excitation.py:565`: "Ẅ := Ẇ = AW + WAᵀ + Σ̃Σ̃ᵀ" — a stray Ẅ. The returned `C0` is the covariance after
    propagation to t = period (equal to C₀ to 10⁻¹⁵ relative), fine but the name suggests the start of the period.

15. **[fixed during the review — record it]** The `exp_pe_general.py` I read at the start (mtime before 22:00) integrated the Lorenz orbit to
    T_lor = (T_END + 200)/5 + 1 = 101 Lorenz time units, while `ExpAverageLyapunov` samples Ψ on [0, T_END + 12T] = [0, 540] → Lorenz time
    108 for T = 20, so the `CubicSpline` *extrapolated* over 7 Lorenz time units. With that version I get garbage at T = 20 (Λ = 2233,
    p₂ = 8.9·10¹⁴, ‖N̄₂‖ = 2.4·10¹³, every condition failing) — the exponential average pulls the e⁻¹⁰-weighted extrapolation blow-up
    into [0, 300]. The script was changed at 22:00:13 to T_lor = (T_END + 12·20 + 20)/5 + 1 = 113 and the table was regenerated at 22:03:40;
    importing the current module, the T = 20 row reproduces exactly (λ = 0.0506, Λ = 0.2424, ρ = 23.78, ρ_D = 158.70, ρ_DD = 3118.95,
    p₁ = 5.31, p₂ = 19.40, π₁ = 23.53, π₂ = 49.41, ‖N̄₂‖ = 0.422; rates 0.00030/0.00018/0.00010 at k = 100/200/400, certificate OK), and
    the same numbers come out with a pad (T_pad = 200) that stays inside the data. *Suggested hardening:* make `ExpAverageLyapunov` (or the
    regressor) raise when Ψ is evaluated outside the tabulated range (`CubicSpline(..., extrapolate=False)` would turn it into NaNs that
    `eigvalsh` reports).

No major issue was found in either result or in the code.

## INDEPENDENT CONFIRMATIONS

### Corollary E.0′.3 (`check_E0p3.py`, `check_E0p3_b.py`, `check_E0p3_c.py`)

- **Derivation re-done.** dV = [½xᵀ𝐇̇x + xᵀ𝐇Ax + ½tr(Σ̃ᵀ𝐇Σ̃)]dt + xᵀ𝐇Σ̃ dB; with 𝐇 = [[βI, −βL],[−βLᵀ, H_ζ + βLᵀL]] (from
  V = ½ζᵀH_ζζ + (β/2)‖e − Lζ‖², matching Corollary E.0′.2's proof) Σ̃ᵀ𝐇Σ̃ = Σᵀ(I + εP₁ + ε²P₂ + βLᵀL)Σ; tr(ΣᵀAΣ) ≤ ‖A‖‖Σ‖_F² gives
  ν̄ ≤ ½‖Σ‖_F²(1 + εp₁ + ε²p₂ + βL_M²); Grönwall gives E V(t) ≤ e^{−γt}EV(0) + ν̄/γ; V ≥ ½(1 − εp₁ − ε²p₂)‖ζ‖² and V ≥ (β/2)‖w‖² give
  the two moment bounds with exactly the stated constants; the slow-regime formula k‖Σ‖_F² hi(hi + βL_M²)/(2c_ζ lo λ_min) is
  reproduced by `bound_zeta` (0.21599 at k = 20, δ = 0.2 by hand). For ζ̇ = −(1/k)M̄ζ + ΣḂ the stationary covariance satisfies
  (1/k)(M̄C̄ + C̄M̄) = ΣΣᵀ, i.e. C̄ = kX with M̄X + XM̄ = ΣΣᵀ (sign and factor as stated); tr(M̄X) = ½‖Σ‖_F² always.
- **`exact_stationary_covariance` vs an independent solve_ivp (DOP853, rtol 10⁻¹⁰) implementation** (Φ, W over the period → discrete
  Lyapunov → dense covariance ODE): max relative difference of C₀ 1.3·10⁻¹² (k = 20) to 2.3·10⁻¹¹ (k = 1000); periodicity
  ‖C(T) − C₀‖/‖C₀‖ ≤ 8·10⁻¹⁵; period means and maxima identical to 7 digits. **Euler–Maruyama Monte Carlo** (2000 paths, dt = 0.5/k,
  horizon 30k/λ_min) on the sinusoid with Σ = diag(0.05, 0.01): E‖ζ‖² ratio to exact 1.0011 (k = 20), 1.0008 (k = 50); E‖e‖² 1.0034, 1.0013.
- **Itô correction checked numerically through stationarity**: mean_t E[V̇_det](t) = mean_t tr(N_sym(t)C(t)) equals −mean_t ν(t) =
  −½ mean_t tr(ΣᵀH_ζζ(t)Σ) to 7 digits in all 24 (instance, k, δ) cases (e.g. −2.085898·10⁻³ vs −2.085898·10⁻³ at k = 20); omitting the
  βLᵀL block or the factor ½ would show up at the 8% / 100% level. Pointwise E[V̇_det](t) ≤ −γ E V(t) and mean EV ≤ mean ν/γ hold in all cases.
- **Bound dominates the exact periodic covariance (mean and max over the period) with non-isotropic Σ and non-isotropic M̄** (bound/max, δ = 0.2 / 0.1):
  sinusoid, Σ = diag(0.05, 0.01): 2.63/2.47 (k = 20), 1.56/1.40 (50), 1.31/1.17 (200); Σ = [0.05; 0.01] (q′ = 1): 2.63/2.46, 1.31/1.17 (k = 20, 200);
  random 2×3 Σ: 1.56/1.40 (k = 50) — all falling toward 1/(1 − δ) as the corollary says for isotropic M̄, independently of Σ.
  Ψ = (sin t, 0.3 cos t), κ = 11.1, Σ non-commuting: 5.12/4.59 (k = 50), 4.29/3.82 (200), 4.13/3.67 (1000); same Ψ with Σ ∝ I: 2.41 (k = 200, δ = 0.2;
  the trace price (2/0.045)/(1/0.5 + 1/0.045)/(1 − δ) = 2.29). Ψ = (sin t, sin(t + 0.4)), κ = 24.3, non-diagonal M̄: 2.65/2.36 (k = 200).
  n = 2, Ψ = [[sin t, cos t],[cos 2t, ½ sin t]], Σ = [[0.05, 0.02],[0, 0.03]]: 3.10/2.80 (k = 50), 1.95/1.74 (200). `bound_V` ≥ max_t E V(t)
  (= ½tr(𝐇C)) in every case (ratios 1.15–4.8); `bound_w` ≥ exact E‖w‖² by ≥ 4·10⁴. `nu_grid` equals my sup_t ν(t) to 5 digits.
- **Exact vs averaged**: exact/averaged − 1 = 1.250/k² on the sinusoid for every Σ (isotropic, diagonal, rank-one, random); O(1/k) for
  anisotropic M̄ with non-commuting ΣΣᵀ (issue 2). E‖e‖² = ‖Σ‖_F²/(2k) to 6 digits everywhere (issue 1).

### Lemma E.0″ (`check_E0dp_A.py` … `check_E0dp_E.py`)

- **Derivation re-done.** T d/dt⟨X⟩_T = ⟨X⟩_T − X (from ⟨X⟩_T(t) = (e^{t/T}/T)∫_t^∞ e^{−r/T}X); Ṗ₁ = 2(M − M̄), Ṡ = F − F̄ for constant M̂, F̂;
  Ḣ = 2ε(M − M̄) + ε²(F − F̄ − Ṁ + Ṃ̄) exactly; order ε: −ζᵀMζ + ζᵀ(M − M̄)ζ = −ζᵀM̄ζ; order ε²: −½ζᵀFζ + ½ζᵀṀζ + ½ζᵀ(F − F̄ − Ṁ + Ṃ̄)ζ =
  ½ζᵀ(Ṃ̄ − F̄)ζ = −ζᵀN̄₂ζ with N̄₂ = ½(F̄ − Ṃ̄) — signs as stated; R_ζ unchanged (products of (I + εP₁ + ε²P₂) with (−εM + ε²ΨᵀΨ̇ − ε³M²)).
  Every per-direction bound of (c) uses only ‖·‖_{M̄(t)} at the current t plus operator norms at other times: ‖P₁(t)v‖ = 2T‖(M̄(t) − M̂)v‖ ≤
  2T(√Λ + ‖M̂‖/√λ)‖v‖_{M̄(t)}; ‖F(r)v‖ ≤ 2p₁Ψ_M²‖v‖ ≤ (2p₁Ψ_M²/√λ)‖v‖_{M̄(t)} for every r; ‖S(t)v‖ ≤ T(‖F̄(t)v‖ + ‖F̂v‖); hence π₂, p₂, n₂ as
  stated. The r_ζ, c₁, c₂, d_w formulas of Lemma E.0′ use only √ρ, √ρ_D, √ρ_DD, Ψ_M, Ψ_D, π₁, π₂ (I checked each of the 11 terms of r_ζ and
  c₂ and the three of c₁ against the terms of R_ζ, ẇ, ζᵀHΨᵀw); the only (T/4)-based quantities of E.0′ (its π₁, π₂, p₁, p₂) are the ones the
  statement redefines, so no periodicity survives. (d): ε²|ζᵀN̄₂ζ| ≤ (ε²n₂/λ)ζᵀM̄ζ = εδ₂ζᵀM̄ζ; Young step and (A), (B) as in E.0′. PE
  equivalence: both directions verified; the identity ∫_t^{t+T′}e^{−(r−t)/T}M = T[M̄(t) − e^{−T′/T}M̄(t + T′)] checked numerically to 4·10⁻¹⁴
  (Simpson vs closed form), e^{−1}α/T = 0.137 ≤ λ_T = 0.391 on the quasi-periodic regressor.
- **Closed-form exponential averages** (⟨cos ωt⟩_T = [cos ωt − ωT sin ωt]/(1 + ω²T²), etc.) for Ψ = (sin t, cos √2 t), T = 4π: agree with
  Simpson quadrature to 7·10⁻¹² and with a backward solve_ivp of Ẏ = (Y − M)/T to 1.7·10⁻⁹. **Code's `exp_average`**: M̄ error 5.8·10⁻⁶
  (O(dt²), dt = 0.005), F̄ error ≤ 2.2·10⁻³ (0.4% of ‖F̄‖; median 1.2·10⁻⁶, largest at t → 300), P₁ 1.4·10⁻⁴, P₂ 2.6·10⁻² (0.4%), N̄₂ 1.1·10⁻³
  (0.4%). Consequently the constants agree to ≤ 0.03%: λ 0.39068 (both), Λ 0.61572/0.61573, ρ 3.9999/4.0000, ρ_D 6.0885 (both), ρ_DD 10.0001
  (both), p₁ 2.9234/2.9235, p₂ 8.9238/8.9258, π₁ 4.4187/4.4188, π₂ 12.712/12.715, n₂ 0.3512 (both); γ to 5 digits at k = 20, 50, 200, δ = 0.2, 0.1
  (e.g. 0.019977 vs 0.019976). The statement's closed-form bounds are loose on this regressor: π₁ ≤ 39.9 (measured 4.42), n₂ ≤ 6.0 (0.35),
  p₂ ≤ 9.45 (8.92) — consistent with the explicit constants failing (A)/(C) at k = 20 in the table.
- **Independent pointwise certificate** (own V from the closed-form M̄ and ODE F̄, own constants, 400 random off-grid (t, e, ζ) per case with
  ‖e‖/‖ζ‖ ∈ {10⁻³, 1, 10³}): sup over all states of V̇/V at those times = −0.03394 vs −γ = −0.01998 (k = 20, δ = 0.2); −0.03400 vs −0.02332
  (20, 0.1); −0.01466 vs −0.01115 (50, 0.2); −0.01469 vs −0.01262 (50, 0.1); −0.003822 vs −0.003057 (200, 0.2); −0.003838 vs −0.003442 (200, 0.1).
  Analytic V̇ (with the lemma's Ḣ) vs a central difference of V along an RK4 step of (LTV): relative difference 2.7·10⁻⁶ (k = 20), 1.7·10⁻⁵
  (50), 2.7·10⁻⁴ (200, cancellation-limited) — a sign error in any ε² term of Ḣ would show at the 10⁻² level.
- **Code internals**: `quadratic_forms`' Hdot vs central differences of its own H on the grid: 8.4·10⁻⁶ relative (k = 50). `lyapunov_exponent_rate`
  vs solve_ivp of Φ (DOP853, 10⁻¹⁰): 0.024848 vs 0.024848 (k = 20, T = 200), 0.024924 vs 0.024924 (T = 1000), 0.002487 vs 0.002487 (k = 200);
  renormalisation bookkeeping and sign (decay rate positive) correct.
- **Attempts to break the lemma** with the code's constants and certificate, measured and explicit, grids of 400 and 4000 times, k ∈ {20, 50,
  100, 200, 400, 1000, 3000} (10⁴ for Lorenz), δ ∈ {0.2, 0.1}, ≈ 250 admissible rows: **no violation**.
  (1) barely PE Ψ = (sin t, sin t + 0.15 cos √2 t): λ_T = 0.0051 (4π); (C) fails at k ≤ 50 (δ₂ = 2.25 at k = 20), holds from k = 100; OK.
  (2) amplitude-modulated Ψ = (sin t, (1 + 0.8 sin 0.07t) cos √2 t): ‖N̄₂‖ = 15.5 (4π), 10.3 (8π), p₂ = 280/412; (D) fails at k = 20, (C) up to
  k = 100 (4π) / 50 (8π); OK from k = 200 / 100 (c_ζ = 0.089 at k = 200, δ = 0.2: margin 0.088, still OK).
  (3) near-resonant beat Ψ = (sin t, cos 1.02t): λ_T = 0.0146, p₂ = 280, π₂ = 677, ‖N̄₂‖ = 14.1 (4π; 0.052, 887, 1210, 23 at 8π); (C) fails up to
  k = 1000 (4π) / 400 (8π); OK at k = 3000 / 1000, 3000.
  (4) constant component Ψ = (0.7 + 0.5 sin t, cos √2 t): OK from k = 20. (5) bursty Ψ = (sin t, (0.2 + 0.8 sin² 0.15t) cos √2 t): OK from k = 20
  (c_ζ = 0.196 at k = 20). (6) n = 2 matrix regressor Ψ = [[sin t, cos √2 t],[cos √3 t, ½ sin t]]: OK from k = 20. (7) q = 3, Ψ = (sin t, cos √2 t,
  sin √3 t): OK from k = 20 (explicit (A) fails at k = 20, δ = 0.1). (8) Lorenz (experiment regressor) T = 5, 10: conditions from k = 100
  (as in the table), OK up to k = 10⁴ with margins (−γ)/(−sup) = c_ζ/hi from 0.07 to 0.90; T = 20: see issue 15 (reproduced exactly with the
  current script; OK; also OK on [0, 1500]).
  In every admissible row sup V̇/V = −2λ_T/k·(1 + O(ε)) and (−γ)/(−sup V̇/V) = c_ζ/(1 + εp₁ + ε²p₂) to three digits: the certificate's slack
  is exactly the lemma's Young-step factor, nothing is hidden in the constants.
- **The published table** (`general_pe.md`, 22:03): all quasi-periodic numbers I recomputed (λ_T, ρ's, p's, π's, ‖N̄₂‖, rates, sup V̇/V) agree;
  the T = 20 Lorenz row reproduces to all printed digits; the claims "1.3–2.5" (Status), λ_∞/λ_T = 1.69, 1.28, 1.13, the Lorenz exact rates and
  the T = 20 ratios 5.0/4.3, 4.2/3.7 are correct; issues 5–9 list the sentences that do not match it.

## Summary for the authors

Both results stand as stated; no step of either proof is wrong or unjustified. What needs fixing is text: one incorrect law in the
E.0′.3 numerics (E‖e‖² = ‖Σ‖_F²/(2k), not ‖Σ‖_F²/k², in PARADIGM.md:1679 and exp_explicit_rate.py:284), the Status line of E.0″ (k = 20
vs k = 100 for the Lorenz regressor; "4–5" only at T = 20), four stale or wrong numbers in E.0″'s Numerics paragraph (lines 1777, 1789,
1791, 1797), the "minimise the norms" clause, a statement that the Lorenz constants are [0, 300] estimates, a vacuous test assertion
(test_excitation.py:316) and a docstring/code mismatch (excitation.py:690/698). The O(k⁻²) agreement of the exact and averaged
stationary covariances should be qualified (it is O(1/k) when ΣΣᵀ does not commute with M̄).
