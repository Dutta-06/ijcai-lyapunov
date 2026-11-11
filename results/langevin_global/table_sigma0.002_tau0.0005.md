# Langevin adaptive control: dither and global kernel placement (scalar plant, m = 32, 3 seeds)

σ = 0.002 (leak toward ϑ⁰ = 0 in the G⁻¹ metric), k_e = 10, G = diag(4.0, 1.0, 1.0), τ₁ = 0.0005, T = 900 s. L̄ = occupation loss ½⟨(f − f̂)²⟩_ρ at T; RMS over the last 10 periods; W₂(collapsed, spread) between the final measures in (a, c, s) for the same seed.

Reference: unregularised offline optimum of L̄ over 32 kernels in the same box (L-BFGS, 5 restarts): **7.63e-06**. Guaranteed rate of Theorem H(ii) at τ₁ (Holley–Stroock lower bound on the averaged dynamics): 2τρ_τ ≥ e^(-8401) s⁻¹, i.e. vacuous at this temperature; the table does not test that rate.

| dither | start | L̄(μ_T) (mean over seeds) | tracking RMS | W₂(final collapsed, final spread) |
|---|---|---|---|---|
| none | collapsed | 1.66e-01 | 0.0576 | 2.896 |
| none | spread | 6.94e-03 | 0.0112 | 2.896 |
| constant | collapsed | 4.80e-02 | 0.0333 | 2.094 |
| constant | spread | 9.19e-03 | 0.0127 | 2.094 |
| annealed | collapsed | 5.21e-02 | 0.0338 | 2.162 |
| annealed | spread | 7.82e-03 | 0.0118 | 2.162 |

Reading (computed from the table, not asserted):
- dither = none: collapsed/spread loss ratio 23.96 — the starts have NOT equalised within the horizon.
- dither = constant: collapsed/spread loss ratio 5.23 — the starts have NOT equalised within the horizon.
- dither = annealed: collapsed/spread loss ratio 6.66 — the starts have NOT equalised within the horizon.
