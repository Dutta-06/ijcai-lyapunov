# Referee report D3dp, round 121: independent verification of the last-annulus certificate (D‴.61)

## Scope

**What I reviewed.**
- The certificate in ref120/own (lastann.py, sha256 78de9795…, identical to round 120's file; selfcheck.py, coverage.py, vs_repo.py, README.md), copied unchanged into W/own.
- The repository at the HEAD when I started, 1ad8abb on claude/vigilant-clarke-hf939b, exported with `git archive` into W/repo. Its exp_gh_peano.py (435ce048…) and D61.md (30dc054b…) are the files round 120 reviewed.
- W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref121.

**The repository moved during the round.** Another session committed fbf237e, "part VA … adopted (pending verification)". Its experiments/gh_lastann.py is byte-identical to lastann.py apart from the header (`diff`) [own]. So every finding below applies to the adopted file as well. D61 at fbf237e (lines 176–177) keeps 1.995 ≤ |η| < 2 a conjecture "until that verification" [md].

**Read in full.**
- review/REVIEW_round120_D3do.md.
- review/drafts/D61.md: (iii)(b), route (2) and its hypotheses, the proof's Ω′ paragraph, (iv), Status, Numerics (V).
- In exp_gh_peano.py: `cert_E`, `cert_seg_parts`, `cert_seg_ball`, `cert_poly`, `cert_box*`, `cert_end*`, `cert_band`, `part_V`, `_cert_E_ball`.
- lastann.py line by line, plus selfcheck.py, coverage.py, vs_repo.py and the README.

**How I ran things.**
- Everything ran under `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with the repository venv (python-flint 0.9.0), with at most two heavy processes of mine at a time.
- I used read-only git only and wrote nothing to the repository or to ref120.
- I killed no process. I stopped four of my own background wait loops with the task tool; they had polled `pgrep -f mutants.py` and matched their own command lines. Other sessions' processes, such as a pytest run of test_moving_width, were left alone.

**My files** (all under W):
- att/: ind.py (independent path check), ident.py and ident4.py (identities), mutants.py, ctl_truecap.py and ctl2.py (negative controls), nanbox.py and nanbox2.py, hilo.py.
- fix/lastann.py: the patched certificate.
- logs/ and mut/{a,b}/results.jsonl.

## Summary

**Design.**
- **The path.** For each η = (2 − ε)e^{iθ} the path is:
  - b = (1 + ε/4)η/2;
  - J: vertically up to the floor Im z = κRe z, with κ = Im A/|Re A| and A = η²/4 − 1;
  - F and C: the floor capped at Im z = θ, out to Re z = δ₀ = 0.05;
  - O: a vertical connector at Re z = δ₀, then a descent polyline per η-box to |1 − z| < 0.02;
  - E: part V's end piece to z = 1.
- **J, F and C** are certified in two exact rescalings of E/X²:
  - A: (u, q, ρ) = (ε/X, θ²/X, X);
  - B: (u, p, ϑ) = (ε/X, X/θ², θ²).

  The box bound is the mean value form f(ρ) ≤ f(0) + ρ·sup f′. f(0) is a closed form in endpoint intervals; f′ comes from forward-mode differentiation in Arb.
- **O** uses a three-variable mean value form in (s, |η|, arg η) over each η-box.

**Verdicts.**
1. **As submitted (= fbf237e's gh_lastann.py), the certificate is not a proof on the whole annulus.**
   - In the η-box [1.995, 2] × [0, π/128], Re A's ball is [±6.2e-3]. κ is then NaN, `handoff_range` returns lo = hi = NaN, and `hi > lo` is False, so the connector is silently dropped.
   - For every η in that box (and its reflections) the certified inner path ends at δ₀ + iY(η), with Y(η) between 0 and 0.0245. The certified polyline starts at δ₀ + 0.01227i. Nothing certifies the stretch in between (D3dp-1).
   - The same happens at ε₀ = 0.11, in the box [1.993125, 2] × [0, π/512].
2. **With a 9-line patch to `handoff_range` (W/fix/lastann.py), I find the certificate a proof.** It is computer-assisted and relies on Arb being correct.
   - Patched run: O certifies 76 η-boxes with 5193 segments, worst bound −1.4993e-7, none failed. At ε₀ = 0.11: 220 η-boxes, 16083 segments, worst −4.23e-9.
   - Every identity re-derives (by hand, and numerically to 1e-84 or better with my own code).
   - Every rigour-relevant line reads correct, apart from D3dp-1 and cosmetics.
   - My independent check at 400 bits passes:
     - 444 η (144 adversarial, 300 random) with 437 720 path points: all E < 0, the largest −6.25e-47;
     - every inner point lies in an accepted box of its branch, with true − certified ≤ −3.05e-16 in E/X² units;
     - every handoff lies in its connector;
     - 5369 accepted outer segments × 10 points: true − bound ≤ −2.19e-5.
3. **Consequences.** Once D3dp-1 is fixed:
   - (iii)(b) holds on the open disc |η| < 2 (computer-assisted), with V and D's branch and bound, or with VA at ε₀ = 0.11 and D;
   - (iv) holds for every frame with r < 1;
   - D‴.58's half-squeeze sharp constant holds for every σ*/σ, and t(v) for every rule with a Mehler frame.

   |η| = 2 is not needed, and the uniformity D61 needs follows (see Consequences).
4. **The self-checks are weaker than reported.**
   - Coverage.py never tests the outer handoff (D3dp-2).
   - The negative controls reject the true path too (D3dp-3).
   - Validate's B0 check and the soundness sampling do not assert what they print (D3dp-4).

**Findings.**
1. **D3dp-1 (major, soundness).** A NaN handoff interval drops the connector of one η-box, so η with |arg η| ≤ π/128 near ±2 are uncertified. Fixed constructively.
2. **D3dp-2 (minor).** No check exercises the outer piece's handoff. "Coverage complete" (round 120) rests on inner-piece sampling only. "Floats only choose the boxes and the outer polyline" (README) is false at `hi > lo`.
3. **D3dp-3 (minor).** The negative-control harness is vacuous. Its closed form for the wrong caps drops −P₀₊², so it rejects the true cap m = 1 as well (0 boxes accepted, 20001 failures, all at the corner). With a corrected harness, m = 1 passes (2011 boxes) and m ∈ {0, 0.5, 0.9, 1.1, 2} fail genuinely.
4. **D3dp-4 (minor).** Validate's B0 comparison only records the deviation and never asserts; the soundness sampling asserts only E < 0, not true ≤ bound. An optimistic, unsound fB0 (mutant C06) passes `lastann.py all`.

   Nine mutants survive every check. One is equivalent on the domain (C02). Five could be unsound in principle (C04, C07, G01, G03, G04); G01, G03 and G04 are equivalent on the data. S05 and G12 change nothing at these arguments and are killed only by my series unit test; G08 survives every check. For these lines, soundness rests on reading the code, which I did; they are correct as written.
5. **D3dp-5 (cosmetic).**
   - `(hi − lo)` in `mvbound` is a float subtraction with no directed rounding (34 of 3185 accepted uses are inexact; effect ≤ 2.1e-16 of the margin).
   - `seg_bound` sets ∇C = 0 when |σ|'s ball contains 0 (never reached).
   - `handoff_range`'s `hi` uses the float f2 rather than the true π/2.
   - `validate` crashes under `--eps0 0.11` (it samples F points with X > δ₀), so round 120's ε₀ = 0.11 numbers come from the pieces only.
   - At ε₀ = 0.11 the outer piece's lower edge is 1.8900000000000001, not 1.89 (harmless: it overlaps 1.8927).
6. **D3dp-6 (minor, text).** Round 120's counts for O (75 η-boxes, 5126 segments) and its "coverage complete" become 76 and 5193 with the fix. D61 at fbf237e must keep "a conjecture there" until the fix is in.

## The identities

All of these I re-derived by hand. "Numerically" means my own code at 300 bits, not lastann's [own].

- **Log term.** With v = |1 − z²|² − 1 = −2Re z² + |z|⁴ and Λ(v) = (v − log(1 + v))/v²:
  - log|1 − z²| = ½log(1 + v) = v/2 − v²Λ/2;
  - adding Re z² gives log|1 − z²| + Re z² = |z|⁴/2 − v²Λ(v)/2;
  - Λ = Σ(−v)^k/(k + 2), tail ≤ r^K/((K + 2)(1 − r)), derivative tail ≤ r^{K−1}/(1 − r);
  - numerically 4.2e-84.
- **h.** h(1 + δ) = δ + δ²/2 − log(1 + δ) = δ²(1 + δH₁(δ)), with H₁ = Σ(−1)^{k+1}δ^k/(k + 3); tails analogous; numerically 1.7e-88.
- **A.** With ω = 1 − ε/2:
  - |Re A| = 1 − ω²cos 2θ = ε(1 − ε/4) + 2ω²sin²θ, which is ρα;
  - Im A = ω²sin 2θ = θ·2ω²sinc 2θ, which is θb̂;
  - S(w) = sinc√w has tail ratio r/((2K + 2)(2K + 3)), with the derivative's ratio analogous; Cv likewise;
  - numerically 4.6e-87 and 2.3e-90.
- **σ.**
  - |σ| − 1 = (ω²|1 + z|² − 1)/(ω|1 + z| + 1), with ω²|1 + z|² − 1 = −ε(1 − ε/4) + ω²X(2 + |z|²/X); numerically 9.6e-88;
  - Re σ − 1 = −ε/2 + ω[X − (1 + X)θ²Cv(θ²) − Yθ sinc θ]; numerically 2.7e-87.
- **The floor.** With Y = κX, Re(Az²) = −|Re A|(X² − Y²) − 2Im A·XY = −X²|A|²/|Re A|; numerically 2.9e-90.
- **Form A.** With Y = θm, Y²/X² = qm²/ρ and θY/X = qm:
  - Q = qm(αm − 2b̂) − αρ;
  - |z|²/X = ρ + qm², and v/X = −2(ρ − qm²) + ρ(|z|²/X)²;
  - D̂ and P̂ as stated;
  - C = X²[D̂₊²(1 + ρD̂₊H₁(ρD̂₊)) − P̂₊²]. This is exact inside the disc, where D̂₊ = P̂₊ = 0, and an upper bound outside, because (Re σ − 1)₊ ≤ (|Re σ| − 1)₊.
- **The branch values of m.** m = b̂/α on the floor, 1 on the cap, and λb̂/α on the junction (u = 4 because X = t = ε/4).
- **Form B** (X = pϑ, Y/X = m̃/θ): Q = m̃(α_Bm̃ − 2b̂) − ϑα_B, |z|²/X = p(ϑ + m̃²), and v/X = −2p(ϑ − m̃²) + X(|z|²/X)². Dropping −P₊² gives an upper bound. Floor iff b̂p < α_B; since sinc 2θ/sinc²θ = θ/tan θ ≤ 1, every q > 1 is on the floor.
- **Closed forms at ρ = 0** (ε = 0, ω = 1, S = 1, Cv = ½, α = s = u + 2q, b̂ = 2):
  - floor: −4q/s − 8q²/s⁴ + D₀₊² − P₀₊², with D₀ = 1 − u/2 + 2q/s² and D₀ − P₀ = q(2/s² + ½ + 2/s);
  - cap: q(s − 4) − q²/2 + C₀, with D₀ − P₀ = 2q; for P₀ ≥ 0 this equals −uq − q²/2 exactly;
  - junction: −(4q/s)λ(2 − λ) − 8λ⁴q²/s⁴ + D₀₊², and D₀₊ = 0 at u = 4;
  - fB0 is fA0 with q = 1/p (checked algebraically).

  Against the true E/X² at ρ = 1e-40, computed at 800 bits from F_η's definition (211 points), the deviation is 5.6e-16 (my reference is in double). lastann's fA0 bounds the limit from above at every point.
- **Convexity bound.** D₊² − P₊² ≤ 2D₊(D − P) is the tangent inequality of the convex x ↦ x₊², whose derivative is 2x₊. 10⁵ samples give a maximum of 0.
- **The ε = 0 face** reproduces round 120: −q²/2, then 7q²/4 − 3q + 1, then −(1 − 1/(2q))²; it is negative for q > 0.
- **The centred form.** F_η(σ) − F_η(η/2) = log|1 − z²| + Re(cz²) + C, with C = h − (|Re σ| − 1)₊² = Λ⁺ (derived).
- **Series enclosures.** At 23 arguments, including the extremes w = π² (for S(4w)) and |v|, |δ| = 0.4999, Arb's sin, cos and log values and derivatives lie inside the series' value and derivative enclosures.
- **Earlier comparisons, reproduced.** Validate: A 9.3e-25, B 7.6e-23 (series truncation), A0 1.5e-59, B0 2.1e-59. vs_repo: 339 points, 2.26e-23.

## The code's rigour

**Comparisons.** Every acceptance decision in the inner pieces is an Arb comparison. Float decisions are:
- box splitting and the polyline (harmless);
- `adm`'s clipping, rounded outward by a factor 1 + 1e-12;
- `seg_bound`'s returned `up(b) < 0`, a rigorous upper end;
- one place that matters: `handoff_range`/`outer_box`. Floats derived from κ's ball go through Python's min and max, which do not propagate NaN in a defined way, and then `hi > lo` decides whether the connector exists (D3dp-1).

**The defect's mechanism** [own: W/att/nanbox.py; logs/ind.log, 264 η flagged, θ ≤ 0.0224]. Of the 75 accepted η-boxes, only [1.995, 2] × [0, π/128] has a NaN handoff. Its 77 polyline segments were accepted without a connector.

**The fix** (W/fix/lastann.py, sha256 66422a18…):

```
kl, kh = float((kap*delta0).lower()), float((kap*delta0).upper()); th_hi = up(Tb)
if not (isfinite(kl) and isfinite(kh)): lo, hi = 0.0, th_hi      # 0 ≤ min(κδ₀, θ) ≤ θ_hi
else: lo = min(nextafter(kl, -inf), f1); hi = max(min(nextafter(kh, inf), th_hi), 0.0); lo = max(min(lo, hi), 0.0)
assert isfinite(lo) and isfinite(hi) and lo <= hi
```

**Endpoint versus midpoint-radius.**
- Closed forms use the `I` class. Its endpoints are always exact arbs from `.lower()`/`.upper()`, so `I.pos`'s tests are decisive. Its `mul`, `div`, `sq` and `pow` are correct enclosures.
- The forms and derivatives use balls, where only width, not sign, is at stake.
- `pos()` on balls is rigorous.

**NaN.** NaN derivatives fail (`dpos is None`), NaN v0 and fl fail, a NaN series radius fails, and an inf κ is handled correctly. The only escape is the one above.

**Mean value form.**
- With lo = 0, sup f′ is taken over the whole [0, hi] (dual on hull(0, hi)).
- With lo > 0, f(lo) is bounded by min(direct value, f0 + lo·sup f′[0, lo]), and the rest by (hi − lo)·max(sup f′[lo, hi], 0).
- The 'd' route needs v0 ≤ 0 and sup f′ < 0, and the path has ρ ≥ t > 0.
- The MVT is applied at fixed (u, q), so the curved path in (u, q, ρ) does not matter.
- fA is C¹ in ρ: x₊² and x₊³ are C¹, and `D.pos`'s derivative hull(0, x′) encloses 1{x > 0}x′.
- Removable singularities are absorbed by the entire or analytic series. The only divisions are by α (an exception, hence failure, at the corner; the corner region is cap-only, which never divides) and by ω|1 + z| + 1 ≥ 1.
- 880 difference quotients all lie inside the dual enclosures.

**Branch choice.** `branchesA` compares α and b̂ in balls over ρ ∈ [0, ρ_b]; undecided boxes certify both branches, and at θ = 0 the two branches coincide. Junction below the cap: in A, α ≥ 4 − ε > 2 ≥ b̂; in B, `Bfloor_ok` checks it per box.

**Skipped boxes** lie outside {uρ ≤ ε₀, ρ ≤ δ₀} (A) or {pϑ ≤ δ₀, upϑ ≤ ε₀, ϑ ≤ π²/4} (B), with ρ ≤ ε₀/u_lo and ϑ ≤ min(δ₀/p_lo, ε₀/(u_lo p_lo)) rounded outward.

**Outer.** `seg_bound` is correct:
- ∂_s = Re(g′H) + Re(conj G·(η/2)H);
- ∂_R = Re(z²Re^{2iθ}/2) + Re(conj G·e^{iθ}(1 + z)/2);
- ∂_θ = Re(2icz²) + Re(conj G·iσ);
- G = (r − 1/r)₊σ/r − 2sgn(a)(|a| − 1)₊ is valid on both sides of |σ| = 1, and E is C¹ away from z = ±1, which is checked;
- rT includes the true π/2;
- η-boxes and segments tile exactly (shared float midpoints).

**End piece.** It runs from |1 − z_n| < 0.02 (checked in balls) to z = 1. The bound is log 0.02 + log 2.02 + 1.0404 + h(2.02) = −1.3314, valid for |η| ≤ 2 because h increases on [1, ∞) and C ≤ h.

**Cosmetics** (D3dp-5): float `hi − lo`; `seg_bound`'s `else acb(0)`; the float f2 in `hi`.

## Coverage

**Parameters.**
- Path points with X ≤ δ₀ have u ∈ [ε/δ₀, 4], q = θ²/X, and uρ = ε ≤ ε₀.
- q ≤ 2.5 is form A; q ≥ 2 (p ≤ 0.5) is form B. Every q > 1 is on the floor, so B needs only F and J.
- The junction has q = 4θ²/ε and p = 1/q.
- The root grids tile [0, 4] × [0, 2.5] × [0, δ₀] and [0, 4] × [0, 0.5] × [0, π²/4]. The B&B leaves no gaps.

**Edges.**
- 1.995f − (2 − 0.005f) = 1.07e-16 > 0, and the outer piece's R1 equals 1.995f exactly, part V's last edge Rb.
- The true π/2 is covered (ϑ ≤ up(π²/4); `th_ball`).
- ε → 0 is covered (u = 0 in A; |η| = 2 included in O); θ → 0 is covered (q = 0).

**Sampling.** Coverage.py reproduces 28 000/28 000 inner points and 9401 points with true − bound ≤ −7.4e-6. My ind.py adds every outer handoff and polyline. Only that check found D3dp-1; after the fix it finds 0 defects.

**Continuity.** b → J → F → C reaches δ₀ + iY(η). The connector sub-segment then leads to δ₀ + iy_c, followed by the polyline and the end piece. With the fix every piece is certified on its closed η-box.

**Route (2)'s hypotheses.**
- t = ε/4 < ε/(2 − ε).
- The chord is in the open unit disc: (1 − ε/2)(1 + ε/4) < 1.
- Re(φ″u²) = 8Re A/|η|² < 0.
- On the chord and the 0-side (real z ∈ (−1, t]), C = 0 and E ≤ log(1 − x²) + x² < 0 except at 0.
- Ω′'s conditions are met as D61 states them.
- The reflections η ↦ η̄ and η ↦ −η leave E and the z-path unchanged (C depends on |σ| and |Re σ| only).

## Attacks

**Runs.**
- `lastann.py all` reproduces round 120's line exactly: A-J 40, A-FC 2011, B-J 320 (634 skipped), B-F 769 (1319 skipped), O 75/5126 with worst −1.4993e-7, end −1.3314.
- Selfcheck: 1500, 1200 and 128/400 reproduce. Coverage and vs_repo reproduce.
- At ε₀ = 0.11 (pieces only): 54/2533/413/1189 boxes and 220/16082 for O.

**Mutants (37, on the patched file).**
- Kill criteria: K1 validate; K2 a piece fails; K3 my ind.py.
- 28 killed:
  - F01–F10, B01–B03, S01–S04 and C01, C03, C05 by K1. Several of these (F02, F04, F05, F08, F09, F10, B02, B03) pass the branch and bound itself, so validation carries the weight.
  - C06, G02, G05, G06, G07, G09 and G11 by K3; G10 by K2.
- 9 survive:
  - C02 is equivalent on the domain (P₀ < 0 on A's floor).
  - C04, C07, G01, G03, G04 and G08 are non-equivalent machinery or closed-form changes.
  - S05 and G12 are near-equivalent at these arguments; my series unit test (ident.py (6)) kills both, while G08 survives (6) and (7).

**Adversarial η.** I took ε ∈ {1e-15, 1e-12, 1e-9, 1e-6, 1e-4, 1e-3, 0.004999, 2 − 1.995f}, crossed with θ ∈ {0, 1e-15, 1e-12, εφ (φ = 0.5, 1.44, 5, 20), √εψ (ψ = 0.1, 1/√2, 1, 3), the floor/cap switch at δ₀, 0.05, 0.3, √(QA·δ₀), 1, π/2 − 1e-15, π/2 (Arb's)}, plus 300 random η. The results are in Verdict 2; the closest inner margin is at ε = 1e-15.

**Negative controls, corrected** (logs/ctl2.log, att/ctl_truecap.py): the true cap passes; wrong caps fail at the corner, where f(0) ≈ q(m − 1)² + qu(m² − 2m − 1)/2 > 0 for small u.

## Consequences

With D3dp-1 fixed:
- **(iii)(b) holds on |η| < 2** (computer-assisted). The pieces are D (|η| ≤ 1.8927), V (to 1.995) and VA (1.995 ≤ |η| < 2); alternatively D and VA at ε₀ = 0.11 (patched: all pass).
- **Uniformity on compact K ⊂ {|η| < 2}.** Since ε ≥ ε_K > 0 and the inner path has X ≥ ε_K/4:
  - E ≤ X³·max(sup f′) < 0 on the 'd' route, and E ≤ X²·max(bd) < 0 on the 'v' route, over finitely many boxes;
  - the outer bound is ≤ −1.5e-7 per η-box;
  - the paths have bounded length, and the chord ends vary continuously.

  So D61's "open, margins continuous" holds with these per-box paths.
- **(iv) holds for every r < 1.** c_need(r)² = 4r ln(1/r)/(1 − r) < 4 because ln x < x − 1, and the mass bound 4r < c² holds, so a disc |η| ≤ c < 2 suffices. Hence D‴.58's half-squeeze sharp constant holds for every σ*/σ (t < 1), and t(v) for every rule with a frame.

**Still conditional or not done:**
- the correctness of Arb/python-flint;
- D61's analytic steps behind route (2) (the contour bound by Λ⁺, the band near |σ| = 1), refereed in rounds 99, 104 and 108 and not re-refereed here;
- D‴.58's rules without a frame (|ρ| ≥ 1).

|η| = 2 is not needed.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dp-1 | major | In the η-box [1.995, 2] × [0, π/128], Re A's ball [±6.2e-3] makes κ NaN, `handoff_range` returns NaN, and `hi > lo` is False, so there is no connector. η with \|arg η\| ≤ π/128 near ±2 have an uncertified gap between δ₀ + iY(η) and δ₀ + 0.0123i. The same occurs at ε₀ = 0.11 ([1.993125, 2] × [0, π/512]) [own: nanbox.py, ind.log]. | In experiments/gh_lastann.py (fbf237e) `handoff_range` (l. 737–739): the patch above (finite check, else [0, θ_hi]; θ_hi from `th_ball`; assert). In `outer_box` (l. 765): fail the box instead of `if hi > lo else []`. Rerun: O 76/5193, worst −1.4993e-7; at 0.11, 220/16083. Update the raw, md and tests. |
| D3dp-2 | minor | No check covers the outer handoff. Coverage.py samples J/F/C only. Round 120's "coverage complete" and the README's "floats only choose…" are wrong at `hi > lo` [own]. | Test per accepted η-box: lo ≤ hi finite, and Y(η) ∈ [lo, hi] at the corners, θ = 0, the true π/2 and random η (W/att/ind.py). Plant mutant G06 (no connector) in tools/mutants. |
| D3dp-3 | minor | The negative controls are vacuous: `generic_f0` drops −P₀₊², so m = 1 is rejected too (0 boxes accepted, 20001 failed at the corner) [own: ctl_truecap.py]. | Run each control together with m = 1 in the same harness, using the closed form with −P₀₊² and the convexity bound. Then m = 1 gives 2011/0; m = 0, 0.5, 0.9, 1.1, 2 fail [own: ctl2.log]. |
| D3dp-4 | minor | Self-checks under-assert. Validate's B0 never asserts, and the samplers do not assert true ≤ bound, so C06 (an unsound fB0) passes `all`. Survivors C04, C07, G01, G03, G04 and G08 rest on reading [own: mut/*/results.jsonl]. | Assert B0 like A0, and assert worst true − bound ≤ 0. Unit tests: the series against Arb's functions (kills S05, G12), difference quotients against the dual enclosure, `I.__mul__` on sign-mixed intervals, `D.pos` on straddling balls, and `mvbound` on synthetic f with known sign of f′ (aimed at G01, G03, G04, G08). |
| D3dp-5 | cosmetic | Float `hi − lo` (34/3185 inexact, ≤ 2.1e-16 of the margin) [own: hilo.py]; `seg_bound`'s `else acb(0)`; float f2 in `hi`; validate crashes at `--eps0 0.11`; R1 at 0.11 is 1.8900000000000001. | Use `arb(hi) − arb(lo)`; return inf when \|σ\|'s ball contains 0; restrict validate to X ≤ δ₀; state the edge as 2 − ε₀ (float). |
| D3dp-6 | minor | Text: O's counts (75/5126) and "coverage complete" in D61 at fbf237e, the md, HANDOFF and round 120's report predate the fix [md]. | After the fix: D61 Status "certified (part VA)", without the conditional; (iii)(b) on \|η\| < 2, (iv) for r < 1, D‴.58's constant for every σ*/σ; update the numbers. |

## Sources

- D61.md, REVIEW_round120_D3do.md, HANDOFF.md, exp_gh_peano.py at 1ad8abb, and gh_lastann.py, D61.md and test_gh_peano.py at fbf237e [md].
- ref120/own and ref120/logs [own].
- Johansson, Arb, IEEE TC 66 (2017), arXiv:1611.02831 [memory].
- Moore, Interval Analysis (1966) [memory]; Rall, Automatic Differentiation (1981) [memory]; Tucker, Validated Numerics (2011) [memory].
- Sterbenz, Floating-Point Computation (1974): exact subtraction when y/2 ≤ x ≤ 2y [memory].
- IEEE 754 NaN comparison semantics [memory].
- No web search.

## Lesson

1. **A NaN is not an empty interval.** Every float that a ball arithmetic certificate takes out of a ball and lets decide the certificate's shape must be checked finite. That covers whether a segment exists, an interval's ends and a loop's range. The round-120 lesson "never let a NaN pass as 0" applies equally to "never let a NaN pass as nothing to check".
2. **Sample the joints, not only the pieces.** Coverage sampling inside each piece cannot see a hole between pieces. Each handoff needs its own test, at the boxes' corners and degenerate edges (θ = 0, |η| = 2), where balls of A = η²/4 − 1 contain 0.
3. **Run a negative control with its positive twin in the same harness, and make every self-check assert what it reports.** A rejection that the true path would also get proves nothing. An unasserted number is not a check.

**Files:** W/fix/lastann.py (the patch); W/att/ind.py, ident.py, mutants.py, ctl2.py, nanbox.py; W/logs/ind_fix.log, ind.log, ident.log, ctl2.log, run_eps011_*.log; W/mut/a, W/mut/b.
