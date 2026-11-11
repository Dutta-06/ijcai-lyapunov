# Langevin adaptive control: dither and global kernel placement (scalar plant, m = 32, 3 seeds)

σ = 0.02 (leak toward ϑ⁰ = 0 in the G⁻¹ metric), k_e = 10, G = diag(4.0, 1.0, 1.0), τ₁ = 0.004, T = 600 s. L̄ = occupation loss ½⟨(f − f̂)²⟩_ρ at T; RMS over the last 10 periods; W₂(collapsed, spread) between the final measures in (a, c, s) for the same seed.

Reference: unregularised offline optimum of L̄ over 32 kernels in the same box (L-BFGS, 5 restarts): **7.63e-06**. Guaranteed rate of Theorem H(ii) at τ₁ (Holley–Stroock lower bound on the averaged dynamics): 2τρ_τ ≥ e^(-1053) s⁻¹, i.e. vacuous at this temperature; the table does not test that rate.

| dither | start | L̄(μ_T) (mean over seeds) | tracking RMS | W₂(final collapsed, final spread) |
|---|---|---|---|---|
| none | collapsed | 2.96e-01 | 0.0769 | 1.691 |
| none | spread | 9.55e-02 | 0.0427 | 1.691 |
| constant | collapsed | 1.10e-01 | 0.0465 | 0.403 |
| constant | spread | 1.06e-01 | 0.0452 | 0.403 |
| annealed | collapsed | 1.06e-01 | 0.0461 | 0.310 |
| annealed | spread | 1.03e-01 | 0.0443 | 0.310 |

Reading (computed from the table, not asserted):
- dither = none: collapsed/spread loss ratio 3.10 — the starts have NOT equalised within the horizon.
- dither = constant: collapsed/spread loss ratio 1.03 — both starts reach the same loss (within 20%).
- dither = annealed: collapsed/spread loss ratio 1.03 — both starts reach the same loss (within 20%).
