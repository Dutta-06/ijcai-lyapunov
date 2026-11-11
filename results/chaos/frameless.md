# Proposition D‴.62: the scaled Gauss–Hermite rules without a Mehler frame (`exp_frameless.py`)

D‴.58's scaled rule (σ = 1): the K-point Gauss–Hermite rule of N(0, v) for the mixing law reweighted to N(0, τ²), e_K = ‖f* − f_K‖/‖f*‖; c = (1 − v/τ²)/2, c_b = c − v/2, s_a = c/(1 − c), s_b = c_b/(1 − c_b), s_d the larger in modulus, s_o the other, ϖ = s_o/s_d, t(v) = |s_d|, C(v) = √2σ*v(1 − s_d)/[(τ² + v)(τ² + vσ*²)(1 + ϖ)]^{1/2}. Labels: **[proved]** (D‴.62), **[formal]** (an expansion not proved), **[numerical]**. Sources: research note, 2026-11-02 (review/RESEARCH_note_frameless.md); the exact-error method: referee round 129's arberr.py.

## R. The representation (a) at small K

e_K² = σ*(v/τ²)E[φ_K(a + b)φ_K(a − b)] **[proved]**, checked against the exact closed form (`err2`) at the same nodes — the identity holds for every rule, so the comparison is exact up to the integration. Moment matching: σ*E[φ_K(iB)²], B ~ N(0, τ²/2), by python-flint's `acb.integral` on |W| ≤ (200τ²)^{1/2} + 5 at 320 bits, the truncation's bound 4e^{−L²/τ²} added to the radius (rigorous):

| τ/σ | K | σ*E[φ_K(iB)²] | the closed form | relative difference (upper bound) | the balls overlap | truncation bound |
|---|---|---|---|---|---|---|
| 1 | 6 | 1.99073239819049107763656818383e-6 | 1.99073239819049107763656818383e-6 | 2.4e-89 | yes | 2.9e-159 |
| 2 | 10 | 0.000453268446733988969426045540901 | 0.000453268446733988969426045540901 | 2.1e-91 | yes | 2.1e-120 |
| 0.25 | 5 | 4.83932381521276022058502924205e-16 | 4.83932381521276022058502924205e-16 | 9.4e-80 | yes | 2.2e-506 |

ρ > 1 (a = iA, b = iB: φ_K on the imaginary axis) and ρ < −1 (a, b real): the M × M Gauss–Hermite product rule of N(0, 1)² at 320 bits (not rigorous; a check of the normalisation and of the Fubini step):

| τ/σ | rule | kind | K | M = 120 | M = 160 | the closed form | relative difference (M = 160) | change M = 120 → 160 |
|---|---|---|---|---|---|---|---|---|
| 0.25 | 2·v_hs | rho>1 | 6 | 3.64133142323250094722979029342e-6 | 3.64133142323250094722979029342e-6 | 3.64133142323250094722979029342e-6 | 9.3e-69 | 2.3e-50 |
| 1 | 2·v_hs | rho>1 | 8 | 8.47993902009990696491448201162e-6 | 8.47993902009990696491448201162e-6 | 8.47993902009990696491448201162e-6 | 7.8e-47 | 5.3e-34 |
| 0.5 | 1.25·v_hs | rho>1 | 6 | 3.33387650115369551349969038648e-10 | 3.33387650115369551349969038648e-10 | 3.33387650115369551349969038648e-10 | 2.3e-70 | 2.3e-70 |
| 0.25 | 0.5·v_hs | rho<-1 | 6 | 1.19923876238414481883835396278e-6 | 1.19923876238414481883835396278e-6 | 1.19923876238414481883835396278e-6 | 7.5e-66 | 3.5e-50 |
| 1 | 0.5·v_hs | rho<-1 | 8 | 2.38458788970181488001774739917e-6 | 2.38458788970181488001774739917e-6 | 2.38458788970181488001774739917e-6 | 1.0e-46 | 7.0e-34 |
| 0.5 | 0.8·v_hs | rho<-1 | 6 | 2.12812254666580889639359625742e-10 | 2.12812254666580889639359625742e-10 | 2.12812254666580889639359625742e-10 | 2.3e-67 | 2.3e-67 |

Largest relative difference: one-dimensional 9.4e-80, two-dimensional 1.0e-46; the largest change M = 120 → 160 7.0e-34.

ρ = −1 (b ≡ 0, a ~ N(0, c) real, c = τ²/(2σ*²)): σ*(v/τ²)E[φ_K(a)²] = E[φ_K(a)²]/σ* by `acb.integral` on |x| ≤ (200/β)^{1/2} + 5, β = 1/(2c) − 1 = 1/τ², at 320 bits, the truncation's bound (2cβ)^{−1/2}e^{−βL²} (from 0 ≤ φ_K(x) ≤ e^{x²/2}, D‴.61 (i)) added to the radius (rigorous):

| τ/σ | K | E[φ_K(a)²]/σ* | the closed form | relative difference (upper bound) | the balls overlap | truncation bound |
|---|---|---|---|---|---|---|
| 1 | 6 | 9.32876217406905042144538008886e-7 | 9.32876217406905042144538008886e-7 | 7.5e-89 | yes | 1.0e-159 |

Pinned for the tests (`r_2d` with small product rules, values to 60 digits; not accuracy checks): τ/σ = 0.5, 1.25·v_hs, K = 3, M = 50: 1.57174928393024411529925578509511765330211000704578272574518e-5; τ/σ = 1, 2·v_hs, K = 4, M = 32: 0.00307800847557830583261338820479002406014244547325926808950291; τ/σ = 1, 0.5·v_hs, K = 4, M = 32: 0.00103529567185696350108105966027561166687922231074875088923787.

## S. The exact errors at consecutive K **[numerical]**

e_K for K = 2 … K_max (python-flint arb at 2K·log₂(1/t(v)) + 450 bits; referee round 129's arberr.py: the nodes Newton iterates from LAPACK's, independent of A6's decimal code), stored to 60 digits; the radius is the arithmetic's at those nodes; e_{K_max} recomputed with 128 more bits:

| τ/σ | rule | kind | t(v) | K | bits | largest relative radius of e_K | e_{K_max} with more bits: relative difference | seconds |
|---|---|---|---|---|---|---|---|---|
| 0.25 | moment matching | MM | 0.03030303030 | 2…200 | 470…2467 | 1.7e-112 | the same 60 digits | 34.1 |
| 0.5 | moment matching | MM | 0.1111111111 | 2…300 | 462…2351 | 9.3e-103 | the same 60 digits | 109.6 |
| 1 | moment matching | MM | 0.3333333333 | 2…400 | 456…1717 | 1.8e-96 | the same 60 digits | 160.3 |
| 2 | moment matching | MM | 0.6666666667 | 2…400 | 452…917 | 1.5e-101 | the same 60 digits | 71.2 |
| 0.25 | 1.25·v_hs | rho>1 | 0.1260538456 | 2…240 | 461…1884 | 6.8e-111 | the same 60 digits | 39.0 |
| 0.25 | 2·v_hs | rho>1 | 0.3467367306 | 2…400 | 456…1672 | 3.0e-105 | the same 60 digits | 154.0 |
| 0.5 | 1.25·v_hs | rho>1 | 0.1658124885 | 2…300 | 460…2005 | 7.9e-105 | the same 60 digits | 83.7 |
| 0.5 | 2·v_hs | rho>1 | 0.3819660113 | 2…360 | 455…1449 | 4.5e-108 | the same 60 digits | 94.7 |
| 1 | 2·v_hs | rho>1 | 0.4775922501 | 2…400 | 454…1302 | 4.9e-104 | the same 60 digits | 109.2 |
| 0.25 | 0.5·v_hs | rho<-1 | 0.3467367306 | 2…300 | 456…1366 | 2.1e-88 | the same 60 digits | 50.7 |
| 0.25 | 0.8·v_hs | rho<-1 | 0.1260538456 | 2…240 | 461…1884 | 2.1e-104 | the same 60 digits | 39.6 |
| 0.5 | 0.5·v_hs | rho<-1 | 0.3819660113 | 2…360 | 455…1449 | 2.9e-79 | the same 60 digits | 104.7 |
| 0.5 | 0.8·v_hs | rho<-1 | 0.1658124885 | 2…300 | 460…2005 | 8.9e-97 | the same 60 digits | 83.5 |
| 1 | 0.5·v_hs | rho<-1 | 0.4775922501 | 2…400 | 454…1302 | 5.6e-72 | the same 60 digits | 112.7 |
| 1 | v = τ²/σ*² (ρ = −1) | rho=-1 | 0.3333333333 | 2…300 | 456…1400 | 1.8e-93 | the same 60 digits | 51.0 |
| √2 | moment matching | MM | 0.5000000000 | 2…300 | 454…1050 | 2.3e-108 | the same 60 digits | 37.8 |
| 2√2 | moment matching | MM | 0.8000000000 | 2…300 | 451…643 | 2.5e-111 | the same 60 digits | 23.5 |

14 rules, 4586 values, the largest relative radius 5.6e-72; 1247 s in all.
Beyond A6's rules (F5): 3 rules, 897 values, the largest relative radius 1.8e-93; 112 s.
S against A6's stored 26-digit errors where they overlap (K − 1, K at K in A6_K): the largest |e_S/e_A6 − 1| 4.2e-26 over 168 values.

## F. The constants and the fits

### F1. The model's numbers

t(v) and C(v) **[proved]** ((b): e_K² ~ C(v)t(v)^{2K}); |s_o|, q = ϖ² and C_o/C **[proved]** ((c): J_o ~ C_o|s_o|^{2K}); s₁ and c₂ = s₁/2 **[proved]** ((d), by hand since referee round 139: review/RESEARCH_note_c2.md); s₁′ (s₁ with d and o exchanged), B₀ = (C_o/C)(1 − q)/(2tq) and B₁/B₀ = s₁′ − s₁ **[formal]** ((c), the transient):

| τ/σ | rule | ρ | kind | t(v) | C(v) | s₁ | c₂ | \|s_o\| | q | C_o/C | s₁′ | B₀ | B₁/B₀ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | moment matching | 1.00000 | MM | 0.03030303030 | 0.7394902765824 | 0.00162833299410 | 0.000814166497049 | — | — | — | — | — | — |
| 0.5 | moment matching | 1.00000 | MM | 0.1111111111 | 0.8281733249999 | -0.104000000000 | -0.0520000000000 | — | — | — | — | — | — |
| 1 | moment matching | 1.00000 | MM | 0.3333333333 | 1.088662107904 | -0.171875000000 | -0.0859375000000 | — | — | — | — | — | — |
| 2 | moment matching | 1.00000 | MM | 0.6666666667 | 1.521451548625 | -0.0935000000000 | -0.0467500000000 | — | — | — | — | — | — |
| 0.25 | 1.25·v_hs | 6.61212 | rho>1 | 0.1260538456 | 0.6663247689700 | -0.149237084205 | -0.0746185421027 | 0.09611796798 | 0.5814293663 | 0.850006759918 | -0.1376659156 | 2.427214049 | 0.011571169 |
| 0.25 | 2·v_hs | 16.5076 | rho>1 | 0.3467367306 | 0.9156974364610 | -0.184873969084 | -0.0924369845420 | 0.3197938296 | 0.8506295388 | 0.941149256337 | -0.1899151843 | 0.2383158459 | -0.0050412152 |
| 0.5 | 1.25·v_hs | 1.84458 | rho>1 | 0.1658124885 | 0.7911996228148 | -0.158691259843 | -0.0793456299216 | 0.05572809000 | 0.1129571802 | 0.524990971115 | -0.1142879188 | 12.43183805 | 0.044403341 |
| 0.5 | 2·v_hs | 4.52786 | rho>1 | 0.3819660113 | 0.9862320002593 | -0.177022048219 | -0.0885110241094 | 0.2828596527 | 0.5483938411 | 0.798831030399 | -0.1950736387 | 0.8611275797 | -0.018051590 |
| 1 | 2·v_hs | 1.58579 | rho>1 | 0.4775922501 | 1.179105345634 | -0.150996760566 | -0.0754983802829 | 0.1715728753 | 0.1290573422 | 0.475237150076 | -0.1915080866 | 3.357607774 | -0.040511326 |
| 0.25 | 0.5·v_hs | -32.9697 | rho<-1 | 0.3467367306 | 0.2220892500962 | 3.36076057505 | 1.68038028752 | 0.3197938296 | 0.8506295388 | 0.999971084858 | 2.672647476 | 0.2532105863 | -0.68811310 |
| 0.25 | 0.8·v_hs | -8.23106 | rho<-1 | 0.1260538456 | 0.4137151854575 | 0.439276729183 | 0.219638364592 | 0.09611796798 | 0.5814293663 | 0.903132182413 | 0.2729516532 | 2.578914927 | -0.16632508 |
| 0.5 | 0.5·v_hs | -8.88854 | rho<-1 | 0.3819660113 | 0.2205281794165 | 4.51872070453 | 2.25936035227 | 0.2828596527 | 0.5483938411 | 0.998538787999 | 1.937925580 | 1.076409475 | -2.5807951 |
| 0.5 | 0.8·v_hs | -2.18034 | rho<-1 | 0.1658124885 | 0.4529090919388 | 0.726673566907 | 0.363336783453 | 0.05572809000 | 0.1129571802 | 0.656238713893 | 0.1009336240 | 15.53979756 | -0.62573994 |
| 1 | 0.5·v_hs | -2.65685 | rho<-1 | 0.4775922501 | 0.2084383464078 | 10.1351971970 | 5.06759859850 | 0.1715728753 | 0.1290573422 | 0.950474300152 | 0.6481812161 | 6.715215547 | -9.4870160 |

At moment matching C = (1 + λ)^{3/2}/√2 and s₁ = 1/16 − λ/(4(1 + λ)²) − 2λ(1 − λ)/(1 + λ)³, rational with τ² (`mm_C2`, `mm_s1`): τ/σ = 0.25: C² = 19652/35937, s₁ = 8/4913, c₂ = 4/4913; τ/σ = 0.5: C² = 500/729, s₁ = -13/125, c₂ = -13/250; τ/σ = 1: C² = 32/27, s₁ = -11/64, c₂ = -11/128; τ/σ = 2: C² = 125/54, s₁ = -187/2000, c₂ = -187/4000.

### F2. The limit and s₁ from the exact errors **[numerical]**

Least squares (F_PREC = 1000 bits, the data's midpoints exact) of e_K²/(C(v)t(v)^{2K}) on L + Σ_{j ≤ n}s_jK^{−j}, plus q^K(g₀ + … + g_3K^{−3}) where a second ridge exists and q^{K_min} > 1e-28, on three windows: [K_max//2, K_max] with 10 powers, [K_max//2 − 40, K_max] with 14 powers, [max(60, K_max//4), K_max] with 18 powers. |L − 1| is the limit's distance from the proved C(v) (its smallest and largest over the windows); s₁ and s₂ are the last window's with L = 1 fixed:

| τ/σ | rule | kind | K | \|L − 1\| (three limit-free fits) | s₁ of (d) (proved) | fitted s₁ (L = 1) | relative difference | fitted s₂ |
|---|---|---|---|---|---|---|---|---|
| 0.25 | moment matching | MM | 2…200 | 2.8e-32 … 7.9e-24 | 0.00162833299410 | 0.00162833299410 | 2.3e-27 | 0.056381755 |
| 0.5 | moment matching | MM | 2…300 | 6.7e-33 … 2.2e-24 | -0.104000000000 | -0.104000000000 | 9.1e-30 | 0.083904000 |
| 1 | moment matching | MM | 2…400 | 1.8e-35 … 3.6e-26 | -0.171875000000 | -0.171875000000 | 2.1e-32 | 0.011352539 |
| 2 | moment matching | MM | 2…400 | 1.3e-37 … 1.9e-27 | -0.0935000000000 | -0.0935000000000 | 3.2e-34 | -0.0026928750 |
| 0.25 | 1.25·v_hs | rho>1 | 2…240 | 4.2e-30 … 6.2e-24 | -0.149237084205 | -0.149237084205 | 2.0e-23 | 0.081639441 |
| 0.25 | 2·v_hs | rho>1 | 2…400 | 1.3e-25 … 7.2e-21 | -0.184873969084 | -0.184873969084 | 8.5e-19 | 0.058036330 |
| 0.5 | 1.25·v_hs | rho>1 | 2…300 | 3.7e-34 … 1.7e-26 | -0.158691259843 | -0.158691259843 | 4.1e-31 | 0.071234706 |
| 0.5 | 2·v_hs | rho>1 | 2…360 | 1.9e-35 … 3.2e-26 | -0.177022048219 | -0.177022048219 | 9.9e-33 | 0.048818316 |
| 1 | 2·v_hs | rho>1 | 2…400 | 1.6e-37 … 4.2e-27 | -0.150996760566 | -0.150996760566 | 2.4e-34 | 0.024330574 |
| 0.25 | 0.5·v_hs | rho<-1 | 2…300 | 3.0e-12 … 9.4e-11 | 3.36076057505 | 3.36076057412 | 2.7e-10 | -30.715381 |
| 0.25 | 0.8·v_hs | rho<-1 | 2…240 | 5.9e-21 … 4.4e-16 | 0.439276729183 | 0.439276729183 | 2.4e-18 | -0.93601943 |
| 0.5 | 0.5·v_hs | rho<-1 | 2…360 | 9.6e-12 … 8.5e-10 | 4.51872070453 | 4.51872070058 | 8.7e-10 | -54.233955 |
| 0.5 | 0.8·v_hs | rho<-1 | 2…300 | 8.3e-19 … 3.2e-15 | 0.726673566907 | 0.726673566907 | 2.6e-16 | -2.0994000 |
| 1 | 0.5·v_hs | rho<-1 | 2…400 | 6.4e-09 … 7.6e-08 | 10.1351971970 | 10.1351936620 | 3.5e-07 | -263.05374 |

- moment matching (4 rules): the best window's |L − 1| from 1.3e-37 to 2.8e-32 (every window: at most 7.9e-24); s₁ within 3.2e-34 … 2.3e-27 of (d)'s proved value.
- ρ > 1 (5 rules): the best window's |L − 1| from 1.6e-37 to 1.3e-25 (every window: at most 7.2e-21); s₁ within 2.4e-34 … 8.5e-19 of (d)'s proved value.
- ρ < −1 (5 rules): the best window's |L − 1| from 5.9e-21 to 6.4e-09 (every window: at most 7.6e-08); s₁ within 2.4e-18 … 3.5e-07 of (d)'s proved value.

### F3. The second ridge **[numerical]**

Least squares of e_K²/(C t^{2K}) − 1 on Σ_{j ≤ n}s_jK^{−j} + q^K Σ_{j ≤ m}g_jK^{−j}: g₀ against C_o/C (**[proved]** by (c)), g₁/g₀ against the formal s₁′ (J_o's own 1/K coefficient) and g₁/g₀ − s₁ (the fit's) against s₁′ − s₁ = B₁/B₀ (the ratio's transient, (c)), on [k₀, K_max] and [max(20, k₀//2), K_max], k₀ = max(20, min(K_max//4, 30/log₁₀(1/q))). Resolved when the better window's |g₀/(C_o/C) − 1| < 1e-05, roughly below 1e-02, otherwise not resolvable; in brackets the q^K term's size at the window's start against the last fitted power's, (C_o/C)q^{K_min} vs |s_n|K_min^{−n}:

| τ/σ | rule | kind | q (q^{K_max}) | C_o/C | s₁′ | s₁′ − s₁ | fits: [K_min, K_max] n, m: \|g₀/(C_o/C) − 1\|, g₁/g₀, g₁/g₀ − s₁ | verdict |
|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.5·v_hs | rho<-1 | 0.85062954 (8.4e-22) | 0.999971084858 | 2.672647476 | -0.68811310 | [75, 300] 18, 3: 2.5e-04, 2.6098237, -0.750937 (5.4e-06 vs 6.4e-07); [37, 300] 22, 4: 2.3e-04, 2.6305507, -0.730210 (2.5e-03 vs 1.1e-03) | roughly resolved |
| 0.25 | 0.8·v_hs | rho<-1 | 0.58142937 (3.0e-57) | 0.903132182413 | 0.2729516532 | -0.16632508 | [60, 240] 18, 3: 6.9e-05, 0.28199941, -0.157277 (6.7e-15 vs 1.5e-14); [30, 240] 22, 4: 3.1e-07, 0.27297881, -0.166298 (7.8e-08 vs 3.0e-11) | resolved |
| 0.25 | 1.25·v_hs | rho>1 | 0.58142937 (3.0e-57) | 0.850006759918 | -0.1376659156 | 0.011571169 | [60, 240] 18, 3: 1.8e-09, -0.13766623, 0.0115709 (6.3e-15 vs 4.3e-20); [30, 240] 22, 4: 1.4e-10, -0.13766594, 0.0115711 (7.3e-08 vs 7.1e-15) | resolved |
| 0.25 | 2·v_hs | rho>1 | 0.85062954 (7.9e-29) | 0.941149256337 | -0.1899151843 | -0.0050412152 | [100, 400] 18, 3: 2.1e-10, -0.18991525, -0.00504128 (8.9e-08 vs 2.1e-15); [50, 400] 22, 4: 8.5e-13, -0.18991518, -0.00504122 (2.9e-04 vs 5.0e-13) | resolved |
| 0.5 | 0.5·v_hs | rho<-1 | 0.54839384 (1.2e-94) | 0.998538787999 | 1.937925580 | -2.5807951 | [90, 360] 18, 3: 2.0e+13, -230.79466, -235.313 (3.3e-24 vs 4.6e-06); [45, 360] 22, 4: 2.1e+07, -144.82733, -149.346 (1.8e-12 vs 7.4e-02) | not resolvable: at the better window's start (C_o/C)q^{45} = 1.8e-12, below the last fitted power |s_22|·45^{−22} = 7.4e-02; fitted s₂ = -54.23 |
| 0.5 | 0.8·v_hs | rho<-1 | 0.11295718 (7.5e-285) | 0.656238713893 | 0.1009336240 | -0.62573994 | [31, 300] 18, 3: 1.7e+19, -90.633324, -91.3600 (2.9e-30 vs 3.5e-08); [20, 300] 22, 4: 5.0e+11, -69.453346, -70.1800 (7.5e-20 vs 4.3e-05) | not resolvable: at the better window's start (C_o/C)q^{20} = 7.5e-20, below the last fitted power |s_22|·20^{−22} = 4.3e-05; fitted s₂ = -2.099 |
| 0.5 | 1.25·v_hs | rho>1 | 0.11295718 (7.5e-285) | 0.524990971115 | -0.1142879188 | 0.044403341 | [31, 300] 18, 3: 4.7e+05, -90.294855, -90.1362 (2.3e-30 vs 3.4e-21); [20, 300] 22, 4: 6.3e-04, -0.15557695, 0.00311431 (6.0e-20 vs 4.2e-19) | roughly resolved |
| 0.5 | 2·v_hs | rho>1 | 0.54839384 (1.2e-94) | 0.798831030399 | -0.1950736387 | -0.018051590 | [90, 360] 18, 3: 9.9e-10, -0.19507396, -0.0180519 (2.6e-24 vs 4.0e-30); [45, 360] 22, 4: 2.6e-10, -0.19507369, -0.0180516 (1.5e-12 vs 2.4e-19) | resolved |
| 1 | 0.5·v_hs | rho<-1 | 0.12905734 (2.1e-356) | 0.950474300152 | 0.6481812161 | -9.4870160 | [33, 400] 18, 3: 1.7e+28, -96.306817, -106.442 (4.3e-30 vs 2.9e+01); [20, 400] 22, 4: 2.3e+20, -67.065722, -77.2009 (1.6e-18 vs 3.5e+05) | not resolvable: at the better window's start (C_o/C)q^{20} = 1.6e-18, below the last fitted power |s_22|·20^{−22} = 3.5e+05; fitted s₂ = -263.1 |
| 1 | 2·v_hs | rho>1 | 0.12905734 (2.1e-356) | 0.475237150076 | -0.1915080866 | -0.040511326 | [33, 400] 18, 3: 1.6e+04, -96.038712, -95.8877 (2.2e-30 vs 2.4e-24); [20, 400] 22, 4: 7.1e-07, -0.19155479, -0.0405580 (7.8e-19 vs 4.8e-21) | resolved |

- resolved: 5 — (0.25, 1.25·v_hs) 1.4e-10, (0.25, 2·v_hs) 8.5e-13, (0.5, 2·v_hs) 2.6e-10, (1, 2·v_hs) 7.1e-07, (0.25, 0.8·v_hs) 3.1e-07
- roughly resolved: 2 — (0.5, 1.25·v_hs) 6.3e-04, (0.25, 0.5·v_hs) 2.3e-04
- not resolvable: 3 — (0.5, 0.5·v_hs) 2.1e+07, (0.5, 0.8·v_hs) 5.0e+11, (1, 0.5·v_hs) 2.3e+20

### F4. A6's stored errors at K = 160, all 36 rules

R = e_K²/(C(v)t(v)^{2K}) at K = 160 from A6's 26-digit e_K (this raw's A, a copy of results/chaos/width_mismatch_raw.npz's A6): R − 1 and 160(R − 1) against s₁ ((d)'s, proved for the rules without a frame; at a frame s₁ is round 129's formula, ϖ < 0, formal). At a frame C(v) is D‴.61 (iv)'s constant (e), the identity's |C/C₆₁ − 1| at 1000 bits **[proved; checked]**; the cusp (ϖ = −1) is excluded:

| τ/σ | rule | ρ | kind | t(v) | C(v) | \|C/C₆₁ − 1\| | R − 1 | 160(R − 1) | s₁ |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | moment matching | 1.00000 | MM | 0.03030303030 | 0.739490276582 | — | 1.2370e-05 | 0.00197926 | 0.0016283330 |
| 0.25 | Fock | 0 | frame | 0.01538461538 | 4.02832393418 | 2.3e-298 | 0.031732 | 5.07718 | 4.4195709 |
| 0.25 | Mercer-zero | -0.941125 | frame | 0.02943725152 | 0.686588790185 | 5.8e-298 | 0.00087664 | 0.140263 | 0.14086893 |
| 0.25 | Mercer frame | 0.0294373 | frame | 0.01537148874 | 4.27832487141 | 5.9e-298 | 0.029970 | 4.79528 | 4.0852706 |
| 0.25 | 0.5·v_hs | -32.9697 | rho<-1 | 0.3467367306 | 0.222089250096 | — | 0.019954 | 3.19265 | 3.3607606 |
| 0.25 | 0.8·v_hs | -8.23106 | rho<-1 | 0.1260538456 | 0.413715185457 | — | 0.0027098 | 0.433562 | 0.43927673 |
| 0.25 | 1·v_hs | 0.0151550 | cusp | 0.01515499506 | excluded | | | | |
| 0.25 | 1.25·v_hs | 6.61212 | rho>1 | 0.1260538456 | 0.666324768970 | — | -0.00092955 | -0.148727 | -0.14923708 |
| 0.25 | 2·v_hs | 16.5076 | rho>1 | 0.3467367306 | 0.915697436461 | — | -0.0011532 | -0.184510 | -0.18487397 |
| 0.5 | moment matching | 1.00000 | MM | 0.1111111111 | 0.828173325000 | — | -0.00064672 | -0.103475 | -0.10400000 |
| 0.5 | Fock | 0 | frame | 0.05882352941 | 2.04169842653 | 4.1e-299 | 0.010041 | 1.60652 | 1.5732422 |
| 0.5 | Mercer-zero | -0.797959 | frame | 0.1010205144 | 0.638909832024 | 5.1e-299 | 0.0026021 | 0.416343 | 0.42100605 |
| 0.5 | Mercer frame | 0.101021 | frame | 0.05825819935 | 2.54135527658 | 6.2e-299 | 0.0065383 | 1.04613 | 0.99742575 |
| 0.5 | 0.5·v_hs | -8.88854 | rho<-1 | 0.3819660113 | 0.220528179417 | — | 0.026457 | 4.23312 | 4.5187207 |
| 0.5 | 0.8·v_hs | -2.18034 | rho<-1 | 0.1658124885 | 0.452909091939 | — | 0.0044626 | 0.714012 | 0.72667357 |
| 0.5 | 1·v_hs | 0.0557281 | cusp | 0.05572809000 | excluded | | | | |
| 0.5 | 1.25·v_hs | 1.84458 | rho>1 | 0.1658124885 | 0.791199622815 | — | -0.00098903 | -0.158246 | -0.15869126 |
| 0.5 | 2·v_hs | 4.52786 | rho>1 | 0.3819660113 | 0.986232000259 | — | -0.0011045 | -0.176716 | -0.17702205 |
| 1 | moment matching | 1.00000 | MM | 0.3333333333 | 1.08866210790 | — | -0.0010738 | -0.171803 | -0.17187500 |
| 1 | Fock | 0 | frame | 0.2000000000 | 1.01192885125 | 8.0e-300 | 0.010824 | 1.73181 | 1.7812500 |
| 1 | Mercer-zero | -0.464102 | frame | 0.2679491924 | 0.535898384862 | 9.2e-300 | 0.012358 | 1.97720 | 2.0568103 |
| 1 | Mercer frame | 0.267949 | frame | 0.1883451609 | 1.99261596067 | 1.5e-299 | 0.00075083 | 0.120133 | 0.11509415 |
| 1 | 0.5·v_hs | -2.65685 | rho<-1 | 0.4775922501 | 0.208438346408 | — | 0.055943 | 8.95083 | 10.135197 |
| 1 | 0.8·v_hs | -0.535534 | frame | 0.2773958090 | 0.507445213475 | 1.2e-299 | 0.013105 | 2.09679 | 2.1851564 |
| 1 | 1·v_hs | 0.171573 | cusp | 0.1715728753 | excluded | | | | |
| 1 | 1.25·v_hs | 0.737258 | frame | 0.2773958090 | 1.12130609853 | 1.2e-299 | -0.0010054 | -0.160858 | -0.16090608 |
| 1 | 2·v_hs | 1.58579 | rho>1 | 0.4775922501 | 1.17910534563 | — | -0.00094277 | -0.150843 | -0.15099676 |
| 2 | moment matching | 1.00000 | MM | 0.6666666667 | 1.52145154863 | — | -0.00058447 | -0.0935156 | -0.093500000 |
| 2 | Fock | 0 | frame | 0.5000000000 | 0.395284707521 | 5.5e-300 | 0.071292 | 11.4067 | 13.687500 |
| 2 | Mercer-zero | 0 | frame | 0.5000000000 | 0.395284707521 | 5.5e-300 | 0.071292 | 11.4067 | 13.687500 |
| 2 | Mercer frame | 0.500000 | frame | 0.4285714286 | 2.09121444203 | 3.8e-300 | -0.00055745 | -0.0891922 | -0.088000000 |
| 2 | 0.5·v_hs | -0.736068 | frame | 0.6345120047 | 0.168086336244 | 1.5e-299 | 0.19480 | 31.1678 | 43.976004 |
| 2 | 0.8·v_hs | 0.102458 | frame | 0.4730025771 | 0.483257921375 | 9.5e-300 | 0.060640 | 9.70245 | 11.461789 |
| 2 | 1·v_hs | 0.381966 | cusp | 0.3819660113 | excluded | | | | |
| 2 | 1.25·v_hs | 0.605573 | frame | 0.4730025771 | 1.68843369197 | 1.1e-299 | -0.00078141 | -0.125025 | -0.12418625 |
| 2 | 2·v_hs | 0.940983 | frame | 0.6345120047 | 1.50340989572 | 1.5e-299 | -0.00064098 | -0.102557 | -0.10248645 |

14 rules without a frame, 18 frame rules off the cusp (the identity's largest |C/C₆₁ − 1| 5.9e-298, |r₆₁ − t(v)| 2.0e-300), 4 at the cusp (excluded). Without a frame |R − 1| from 1.2e-05 to 5.6e-02, and |160(R − 1)/s₁ − 1| from 1.7e-04 to 2.2e-01 (the O(1/K) remainder s₂/K).

### F5. Beyond A6's rules: ρ = −1, and moment matching at and beyond Kuo–Woźniakowski's threshold **[numerical]**

S_NEW's exact errors against (b): y_K = e_K²/(C(v)t(v)^{2K}), y_K − 1 at a few K, K_max(y_{K_max} − 1) against (d)'s proved s₁, and the limit-free fits of F2. Kuo–Woźniakowski's parametrisation (the kernel e^{−γ²(x − x′)²}, the weight N(0, α)) is α = τ², γ² = 1/(4σ²) for D‴.58's kernel e^{−(y − y′)²/(4σ²)} and N(0, τ²), so 2αγ² = τ²/(2σ²) and moment matching's rate is λ = 2αγ²/(1 + 2αγ²); ρ = −1 (v = τ²/σ*²) is moment matching's dual (u ↦ −u), t(v) = λ with s_d = +λ, on the real axis:

| τ/σ | rule | kind | τ²/(2σ²) = 2αγ² | t(v) | C(v) | C², s₁ (exact) | K | y_K − 1 | K_max(y − 1) | \|L − 1\| (three fits) | fitted s₁ (L = 1): relative difference |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | v = τ²/σ*² (ρ = −1) | rho=-1 | 1/2 | 0.3333333333 | 0.3849001794598 | 4/27, 13/4 | 2…300 | K = 100: 2.9901e-02; K = 200: 1.5540e-02; K = 300: 1.0507e-02 | 3.1520539 (s₁ = 3.2500000) | 8.9e-12 … 4.2e-10 | 3.24999999717: 8.7e-10 |
| √2 | moment matching | MM | 1 | 0.5000000000 | 1.299038105677 | 27/16, -61/432 | 2…300 | K = 100: -1.4125e-03; K = 200: -7.0613e-04; K = 300: -4.7073e-04 | -0.14121829 (s₁ = -0.14120370) | 2.8e-34 … 2.0e-27 | -0.141203703704: 3.2e-31 |
| 2√2 | moment matching | MM | 4 | 0.8000000000 | 1.707629936491 | 729/250, -631/11664 | 2…300 | K = 100: -5.4095e-04; K = 200: -2.7049e-04; K = 300: -1.8033e-04 | -0.054097525 (s₁ = -0.054098080) | 1.9e-35 … 1.4e-25 | -0.0540980795610: 5.8e-32 |

At moment matching the cheap bound e_K² ≤ σ*(K!/(2K)!)²(4K − 1)!!(τ²/2)^{2K} ~ 2^{−1/2}σ*(τ²/2)^{2K} (from |ĝ_K(iB)| ≤ 1; its rate τ²/2 = 2αγ² is Kuo–Woźniakowski's threshold) at K_max: τ/σ = √2 (2αγ² = 1): the bound 1.225, e²_{K_max} = 3.129e-181, λ = 0.500000; τ/σ = 2√2 (2αγ² = 4): the bound 3.653e+361, e²_{K_max} = 1.220e-58, λ = 0.800000.

## T. Lemma T's margin (floating point; a sanity check of a proved lemma) **[numerical]**

sup of Φ⁺ = L⁺(y₁) + L⁺(y₂) − Q over the tail {max|y_i| ≥ 2} against 2 log t(v): on a grid in the representation's own variables (Q with (|c|, |c_b|), the real axis unsubtracted; 1601² points on [−7, 7]²) and by the proof's reduction (R1: −2/k at the corner; R2: the largest Φ(2, y₂) on the edge, 200001 points); the sector forms against 2 log|s| of the sector's own variable:

| τ/σ | rule | kind | t(v) | 2 log t(v) | grid max − 2 log t | at \|y\| (2t^{1/2}) | margin (grid) | margin (reduction) | R1, R2 | sector margins (d, o) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.5·v_hs | rho<-1 | 0.3467367306 | -2.118379 | -3.7e-05 | 1.1813 (1.1777) | 0.8522 | 0.8450 | -3.7681, -2.9633 | 0.8450, 0.9635 |
| 0.25 | 0.8·v_hs | rho<-1 | 0.1260538456 | -4.142092 | -1.4e-05 | 0.7087 (0.7101) | 5.3662 | 5.3373 | -13.8662, -9.4794 | 5.3373, 6.4383 |
| 0.5 | 0.5·v_hs | rho<-1 | 0.3819660113 | -1.924847 | -1.4e-05 | 1.2337 (1.2361) | 0.7238 | 0.7171 | -3.2361, -2.6419 | 0.7171, 1.1652 |
| 0.5 | 0.8·v_hs | rho<-1 | 0.1658124885 | -3.593795 | -2.6e-06 | 0.8137 (0.8144) | 4.5975 | 4.5703 | -10.0618, -8.1641 | 4.5703, 10.0444 |
| 1 | 0.5·v_hs | rho<-1 | 0.4775922501 | -1.477996 | -2.4e-07 | 1.3825 (1.3822) | 0.4820 | 0.4764 | -2.1877, -1.9544 | 0.4764, 2.3152 |
| 0.25 | 1.25·v_hs | rho>1 | 0.1260538456 | -4.142092 | -1.4e-05 | 0.7088 (0.7101) | 5.3662 | 5.3373 | -13.8662, -9.4794 | 5.3373, 6.4383 |
| 0.25 | 2·v_hs | rho>1 | 0.3467367306 | -2.118379 | -3.7e-05 | 1.1813 (1.1777) | 0.8522 | 0.8450 | -3.7681, -2.9633 | 0.8450, 0.9635 |
| 0.5 | 1.25·v_hs | rho>1 | 0.1658124885 | -3.593795 | -2.6e-06 | 0.8138 (0.8144) | 4.5975 | 4.5703 | -10.0618, -8.1641 | 4.5703, 10.0444 |
| 0.5 | 2·v_hs | rho>1 | 0.3819660113 | -1.924847 | -1.4e-05 | 1.2338 (1.2361) | 0.7238 | 0.7171 | -3.2361, -2.6419 | 0.7171, 1.1652 |
| 1 | 2·v_hs | rho>1 | 0.4775922501 | -1.477996 | -2.4e-07 | 1.3825 (1.3822) | 0.4820 | 0.4764 | -2.1877, -1.9544 | 0.4764, 2.3152 |

The sweep: τ/σ in [0.05, 20] (25 geometric), v/v_hs in 8 geometric values of each frameless range ((0.02, 0.999/σ*) and (1.001σ*, 50)), 400 rules, grids of 801² points on [−9, 9]²: the smallest margin 2.5e-06 by the reduction (τ/σ = 20, v/v_hs = 0.02, t(v) = 0.998004) and 1.4e-05 on the grid; the smallest sector margin 2.5e-06; every margin positive: yes; the grid's global maximum never above 2 log t(v) (Lemma A), within 6.1e-08 of it at the closest.

## G. Assumption (F) on the axes **[numerical]**

D_K(η) = K[log ĝ_K(√Kη) − Kη²/4 − ½log(1 − η²/4)] on the real axis and K[log ĝ_K(i√Ky) + Ky²/4 − ½log(1 + y²/4)] on the imaginary (η = iy), at K = 60, 80, 120, 160, 240, 320 (50 digits stored; arb at 2.5K + 200 bits, `gh_nodes`), extrapolated to K = ∞ by the polynomial in 1/K through all the points and through the last four, against round 129's F₁(η²/4) ((d)'s assumption (F): D_K → F₁). Near η = 2 on the real axis the approach is slow (pre-asymptotic, not a contradiction):

| axis | η | η²/4 | D_K at K = 320 | extrapolated (6 points) | extrapolated (last 4) | F₁(η²/4) | \|extrapolated (6) − F₁\| |
|---|---|---|---|---|---|---|---|
| real | 1 | 0.250000 | 0.81945373 | 0.8379604885 | 0.8379418337 | 0.8379629630 | 2.5e-06 |
| real | 1.5 | 0.562500 | 9.2955010 | 10.93218668 | 10.86268767 | 11.02368805 | 9.2e-02 |
| real | 1.8 | 0.810000 | 79.177158 | 143.3139223 | 134.4628359 | 217.0859637 | 7.4e+01 |
| imag | 1 | -0.250000 | -0.14101409 | -0.1410000000 | -0.1410000001 | -0.1410000000 | 1.9e-13 |
| imag | 1.5 | -0.562500 | -0.13844517 | -0.1383120000 | -0.1383120000 | -0.1383120000 | 1.9e-14 |
| imag | 1.9 | -0.902500 | -0.10332964 | -0.1032433426 | -0.1032433426 | -0.1032433426 | 2.1e-14 |

G computed in 0.3 s.

Numbers are produced by `experiments/exp_frameless.py`; nothing is edited by hand.
