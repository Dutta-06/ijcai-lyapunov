# Referee report D3bc, round 56: Proposition D‴.60 (measurement noise and the formed kernels)

**Scope.** I reviewed these files at 17d6217:
- `review/drafts/D60.md`
- `experiments/exp_noise_resonance.py`, and what it imports (`exp_closed_transit`, `exp_dead_zone`, `src/core/closed_loop.{py,c}`, `src/core/mirror_free`)
- `results/chaos/noise_resonance.md` and `noise_resonance_raw.npz`
- `tests/test_noise_resonance.py`
- `PARADIGM.md`, for D‴.31, D‴.33 and D‴.49

Other sessions moved HEAD to 726ca94 during the review. `git diff --stat 17d6217 HEAD` over all of these files is empty. I modified nothing in the repository. All work ran in a `git archive` copy in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref56/D3bc/` (S below), with `nice -n 10`, `OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1` and pytest `-p no:cacheprovider`. I ran at most two processes of my own at a time, except for about a minute when a third was running; I stopped it by its PID.

## Summary

**What holds**

- **Reproducibility is exact.**
  - `--from-raw` gives the committed md byte for byte.
  - The 18 tests pass in 16.7 s.
  - Fresh reruns are bit-identical to raw: an A row (paradigm, ω = 13.9), and a full B3 cell (paradigm, 13.1: `second_order` plus the ±n runs at s = 0.00125).
- **(i) is correct, and classical.**
  - The oscillator law re-derives. `variational` is the exact linearisation (up to O(ε²) from the 18 unformed kernels).
  - √K_max goes 7.881 → 15.759 from m = 3 to 12, a ratio of 2.000.
  - μ = 0.09, 0.06, 0.075 at both m, so μ ≈ 0.1|w| and does not depend on m.
- **`second_order` is right term by term,** and the ±n parity argument is exact.
- **The even part converges to s²⟨θ₂⟩ with an error ∝ s².** Better than the draft claims: at s ≤ 0.0025 it stays accurate to 300 periods.
- **The drift's resonance is real.** At m = 12 the coefficient drift is ⟨w₂⟩ = 221 at 13.1 against 9.6 at 41.3. The m = 12 / m = 3 ratio falls from 20 (13.1) to 1.1 (41.3) and 0.9 (59.3).
- **The resonance-sweep part of (ii)'s averaging story checks out in an exact split.** The dropped 2ω part is 31% of ⟨w₂⟩ for the kernel whose K(t) crosses 13.1² transversally, and 3% for the kernel whose K(t) only touches 13.1² at its maximum.
- **The numbers are quoted correctly.** I found no misquote in the numerics paragraph.

**What does not hold**

1. **(iii)'s mechanism and design rule are not supported.**
   - A kernel leaving its feature is a sudden ejection driven by multi-tone noise, not the resonant drift that (ii) analyses.
   - Single resonant tones at twice the rms never eject.
   - Moving the 16 tones far above √K_max raises m = 12's threshold from 0.02 to 0.08. It still ejects there, while m = 3 never does.
   - So keeping √K_max below the band is not sufficient, and the centres' mobility matters on its own.
2. **Part C's thresholds are single-seed first-passage numbers.** Other seeds shift them by one grid step, and they are not monotone in the level.
3. **(ii)'s condition for the averaging's failure is mis-identified.** It is the commensurability of the tone differences with the reference's harmonics, not "slow beats". The proof's Krylov–Bogoliubov appeal does not apply.
4. **The odd part is not "a bounded offset left by the noise's start".** It alternates sign every period at half-integer tones. At a commensurate tone it grows secularly, beyond the even part.
5. **"⟨θ₂⟩ grows with time" is false beyond about 30 periods.** ⟨w₂⟩ saturates (τ ≈ 130 periods), and ⟨c₂⟩ is constant from period 3.
6. **"New" is overstated.** The four sources are in the paper's own D‴.31(iv)(c) and D‴.33.
7. **Further problems:**
   - The k_e effect is attributed only to damping, but much of it is the noise rectification's sign change (D‴.33(iv)).
   - The WFR threshold claim has no data in D‴.60.
   - Several numbers are hand-computed.
   - 19 of 42 mutants survive both the tests and `--from-raw`.

## What I checked and how

Scripts are in S.

**Reproducibility**
- `NR_OUT=S/out … --from-raw` and `cmp`: identical (md5 f1829f01…).
- pytest: 18 passed.
- `rerun_cells.py`:
  - A (paradigm, 13.9): all nine fields bit-identical.
  - B3 (paradigm, 13.1): w1, c1, w2, c2 at the marks and the s = 0.00125 rows bit-identical (maximum relative difference 0.0).

**Derivations**
- (i) and `second_order` by hand, including each sign and the regressor point.
- The curves' μ(w), v(w) and their derivatives.
- Per-kernel √K_max,i = √(Gm)|w_i|/(σ√e):
  - m = 12: 15.76, 10.51, 13.13
  - m = 3: 7.88, 5.25, 6.57
- The design rule's bound: m ≲ 100/20.70 = 4.83 (4.93 with the actual lowest tone, 10.105).

**Raw** (`inspect_b3.py`): growth across the marks 3, 10, 30, the scaling of the even and odd parts, and the cubic remainder.

**Controls with the C loop** (`controls.py`; the even and odd parts of the period means under ±n)
- **E3:** m = 3 and 12 at nine tones that are half-integer multiples of ω_ref (5.1 … 59.3), s = 0.00125.
- **KE:** k_e = 5, 20, 40 at 13.1 and 41.3.
- **E5:** a commensurate tone 21.4 = 107ω_ref against 21.3, over 60 periods.
- **E2:** 300 periods at 21.37 and with the 16 tones, at s = 0.00125 and 0.0025.
- **E4:** the 16 tones one at a time and together (both metrics), plus five tone pairs (WFR).

**Averaging error** (`split2w.py`): the exact second order split into the response to the quadratic forcing's tone average (DC) and to its 2ω part, from the complex first order. DC + 2ω equals raw w2 to every printed digit.

**Part C** (`trace_C.py`, `seeds_C.py`, `single_C.py`, `eject_probe{,2,3,4}.py`, `farband_C.py`, `wfr_C.py`):
- a dense trace of an ejection;
- seeds 1–3;
- single tones;
- tone pairs (commensurate, detuned, rephased);
- a band moved ×2.5;
- the WFR metric.

**Tests** (`mutants.py`, `mut_impact.py`): 42 single-replacement mutants, each run against the tests and `--from-raw`. For the survivors I measured the numerical effect.

**Literature:** web search; items from memory are marked [memory].

## Verdicts per part

**(i): correct, classical.**
- The law y″ + k_eẏ + Ky = K(ṅ − f̂′n) and δw_i, δc_i follow quasi-statically.
- The A comparison is fair: the same filter, and the start transient decays within about 1 time unit of 126.
- Under the paradigm's metric the centre term dominates K: μφ² ≤ 0.09 against K_max = 248.4.
- The draft never states that the resonances are per kernel: 15.76, 10.51, 13.13. That is why kernels 2 and 3 peak at 10.7, not 13.9 (A table, md:13–14).
- The "amplification" is modest. At 37.1, the paradigm/WFR jitter ratio is 12.2, 7.5, 9.6. That is the mobility ratio v = 10.8, 7.2, 9.0. Resonance adds only ×1.24–1.38 at 13.9 (14.9, 9.0, 12.4 against 10.8, 7.2, 9.0).

**(ii): the algebra holds; three claims do not.**
- **The algebra.** Every term of `second_order` (exp_noise_resonance.py:282–318) matches the expansion. That includes ½f″(e₁² − d₁²) = −½f″n(2e₁ + n), using ∂_x∂_cφ = −∂²_cφ.
- **Parity.** The start x = x_m does not depend on n, so ±n flips the odd orders exactly.
- **Horizon.** The expansion is regular on finite horizons, and in fact uniform to 300 periods at small s. E2: the even parts per s² at s = 0.0025 and 0.00125 agree within 0.31% (21.37) and 1.7% (16 tones) at every mark up to period 300.
- **The drift saturates.** ⟨w₂⟩ (21.37, kernel 1) is 18.1, 32.9, 48.1, 61.6, 70.8, 77.2, 81.5 at periods 30, 60, 100, 150, 200, 250, 300. That is geometric with ratio 0.685 per 50 periods: τ ≈ 130 periods, asymptote ≈ 91, five times the quoted period-30 value. ⟨c₂⟩ is set by period 3 (−13.8 → −15.9 at period 300). It is a bounded second-order offset, not the unbounded parameter drift of Ioannou–Sun.
- **The odd part.** For a tone at a half-integer multiple of ω_ref the centres' period mean alternates sign every period. Per unit s it is ±0.309 for kernel 1 at 13.1 (raw c1 at marks 3, 10, 30: −0.308, +0.309, +0.309), and ±0.0298 at 21.3. At ω = 21.4 = 107ω_ref, ⟨w₁⟩/s grows −0.0011, −0.035, −0.103, −0.188 at periods 1, 10, 30, 60. At s ≤ 0.006 that exceeds the even part.
- **The averaged system.**
  - "Exact along the orbit" is right for the envelopes (the test agrees to 10⁻⁷). The mean is averaged.
  - Krylov–Bogoliubov does not apply: the mean state has fast modes (k_e = 10, √K up to 15.76).
  - What is exact is this: the averaged system's second order is the response to the DC parts, and its error is the period-mean response to the dropped 2ω_k and ω_j ± ω_k parts.
- **The 16-tone failure is entirely the cross terms.** Under WFR, the averaged loop equals the sum of the 16 single-tone exact second orders to four digits: ratios 0.9998, 1.0000, 0.9997. The full second order differs from that sum by the cross terms, which account for the 0.8745.
- **Which pairs carry the cross terms.** Measured as a share of the pair's diagonal:

  | pair difference | relation to the reference's harmonics | cross term, share of diagonal |
  |---|---|---|
  | 3.6004 | 18ω_ref + 4·10⁻⁴ | 56%, −29%, −45% |
  | 2.0029 | 10ω_ref + 0.003 | 38–55% |
  | 0.057, the slowest beat | — | 3% |
  | 1.4992, 3.2981 | about (l + ½)ω_ref | ≤ 0.3% |

  So "slow beats" (D60.md:23, 44) is the wrong criterion.
- **The resonance sweep is supported.** At 13.1 the dropped 2ω part is 31.2%, −4.1%, 2.9% of ⟨w₂⟩(10). Kernel 1's K₁ crosses 13.1² on its flanks; kernel 3's K₃ only reaches 13.13² at its maximum.
- **Two of the listed terms are numerically empty.** The "mobilities' correlation" (μ′, v′) and the blur's w–c term change ⟨w₂⟩ by ≤ 5·10⁻⁴ (mutants S5, S7, S9).

**(iii): the qualitative m = 3 / m = 12 gap is robust; its cause and the design rule are not.**
- **Robust across seeds.**
  - At m = 12, k_e = 10, s = 0.02, kernel 1 leaves under all four seeds (periods 133, 274, 164, 157).
  - At m = 3, k_e = 5, s = 0.04, nothing leaves under any seed.
- **The mechanism is a sudden ejection.**
  - Seed 7, m = 12, k_e = 10, s = 0.02: |Δc| stays ≤ 0.023 for periods 1–131. It then jumps to 0.31 and 0.82 and ends at c = −1.43.
  - At k_e = 5, s = 0.015, kernels 1 and 3 both jump at period 92 to about +1.38.
- **The ejection is not the drift.**
  - Single tones at 13.1, 15.7 and 21.3, at s = 0.02 and 0.04, never eject: max |Δc| ≤ 0.069. Yet w₁ ends at 0.96–1.60.
  - The 16 tones eject at 0.02. So do 5 near-resonant tones (period 116), and the 16 tones shifted by +0.05 or +0.1, or with redrawn phases.
  - A pair (18.90, 22.50, Δω = 18ω_ref + 4·10⁻⁴) with the seed's phases ejects at s = 0.02 (period 218). It does not with other phases, nor with Δω = 3.5004, 3.6 exactly, or 3.7004.
- **The m-dependence is not only resonance.**
  - Moving the 16 tones to 25.3–69.9 rad/s (≥ 1.6√K_max at m = 12) raises m = 12's threshold from 0.02 to 0.08 (kernels 1 and 3 leave at periods 207 and 203).
  - At the same levels m = 3 never leaves (max |Δc| 0.022 against 0.14 at s = 0.08).
  - A pair at 37.8 / 41.402 ejects at m = 12 at s = 0.04, but not at 0.08. At m = 3 it never ejects up to 0.08.
  - So both resonance and mobility matter, and "keep √K_max below the band" is not sufficient. The bound m ≲ 4.8 is never stated, and no m between 3 and 12 is tested.
- **k_e.** Off resonance the drift falls with k_e through the rectification factor (k_e² − ω²)/(ω² + k_e²) (D‴.31(iv)(c), D‴.33(iv)), not through damping.
  - m = 3, ω = 41.3: ⟨w₂⟩ = 12.65, 8.61, 3.86, 0.38 at k_e = 5, 10, 20, 40.
  - m = 3, ω = 13.1: +10.9 → +0.46 → −1.57; the sign reverses.
  - Part C's raw shows the reversal too: m = 3, k_e = 40, s = 0.08 ends with w = 0.869, −0.579, 0.722.
- **WFR.** The threshold claim has no D‴.60 data. My runs support it: WFR at k_e = 5 and 10, s ≤ 0.08, max |Δc| ≤ 0.023.

**Status (D60.md:46–58).**
- Correct for (i).
- For (ii), "⟨θ₂⟩ grows with time" (line 48) is refuted and the odd part's characterisation is wrong.
- (iii) should say: one target, one law, one seed, one horizon (300 periods) and one cut (0.3).
- The second order cannot reach (iii)'s levels. Its relative error at 13.1 is 1.2%, 4.8%, 18.8% at s = 0.00125, 0.0025, 0.005 (∝ s²), so about 75% at s = 0.01.

**Numerics (D60.md:60–82).** Every quoted number matches the md. Several are hand-assembled and not printed by the script:
- 1.0%, 0.58%, 0.12% (line 46);
- 248 (line 37);
- 8% and 3% (line 47);
- the B2 ranges and the WFR ⟨w₂⟩ range;
- "the centres' even part likewise" (line 69). The md prints the centres only at s = 0.005; raw confirms ∝ s² (0.020, 0.080, 0.30 for the middle kernel at 13.1).

Two quotations mislead:
- "Peaks at ω = 13.9" holds for kernel 1 only.
- The 16-tone ⟨c₂⟩ "−45, +9.9, +22" is a snapshot. Kernel 1 reads −14, −10, −45 at periods 3, 10, 30, and between −11 and −98 up to period 300.

The odd-part sentence (lines 73–74) hides two things:
- For the 16 tones the remainder scales by 14.5 and 20.6 per doubling, not 8.
- The single tones sit at half-integer multiples of ω_ref (13.1, 29.3) or 0.03 from one (21.37), and this is not stated.

**Tests.** Of 42 mutants:
- The tests kill 17.
- `--from-raw` kills 9, all of them in the report code (L1–L7, W3, W4).
- **19 survive both:**
  - all five `b3_job` mutants (B1 the "−n" phase moved by π/2; B2 the odd part ×2; B3 the marks off by one; B4 `noise_of(0.5)`, every B3 ratio ×4; B5 an offset in the centres' even part);
  - both `b2_jobs` mutants (W1, W2);
  - `envelope_rhs`'s ē-dependent terms (V2, V3, V4, V5, V6, V9, V11; all tests sit at ē = 0 except those for 'ctl' and 'cross');
  - near-equivalent mutants (S5, S7, S9 change ⟨w₂⟩ by ≤ 5·10⁻⁴; S11 max_step = ∞ and S12 rtol = 10⁻⁶ change it by < 10⁻⁴).
- The one second-order test (21.37, 2 periods, 1% tolerance) cannot see the small terms. Nothing calls `b3_job`.

**Literature.**
- **"Classical" is right but incomplete.** It omits:
  - the paper's own **D‴.31(iv)(c)** (PARADIGM.md:3854–3866), which has exactly the four sources: g = k_eE[e], the compensation's curvature ½D²f̂[Σ], the slope–parameter term, and the correlation ∇²Φ̃E[ε_mδϑ];
  - **D‴.33(i)–(iv)** (4164–4210): the drift S_k + C_k, the threshold, and the sign change of the factor (B − k_e − 2a)/(k_e + B);
  - the adaptive noise canceller as a notch filter whose bandwidth is set by the gain (Widrow et al. 1975; Glover 1977);
  - MRAC's oscillations growing with the speed of adaptation and damped by the closed-loop reference model's gain (Gibson–Annaswamy–Lavretsky 2013), which is the oscillator of (i) with k_e as damping;
  - multi-frequency averaging and passage through resonance (Sanders–Verhulst–Murdock; Neishtadt).
- **Riedle–Kokotović** already gives frequency-dependent slow drift under persistent excitation at low gain.
- **Mean-field particle flows under noise** [memory]: mean-field Langevin (Mei–Montanari–Nguyen 2018; Chizat 2022) and WFR particle flows (Chizat's conic particle gradient descent; Rotskoff et al. 2019). There, the metric sets each particle's mobility and so its noise response.
- **What is new:**
  - the exact (non-averaged) second order along the orbit for tones, checked against the loop's even part;
  - the observation that the paradigm's metric gives the centres the gain Gmw.
- "Exact second-order drift … with its four sources" and "the averaged reduction's failure at that resonance" are not new.

Sources:
- [Kokotović–Riedle–Praly 1985 (SCL)](https://cas.minesparis.psl.eu/~praly/Telechargement/Journaux/1985-SCL-Kokotovic-Riedle-Praly.pdf)
- [UC Berkeley ERL-86-25](https://www2.eecs.berkeley.edu/Pubs/TechRpts/1986/ERL-86-25.pdf)
- [Kosut–Anderson–Mareels](https://researchportalplus.anu.edu.au/en/publications/stability-theory-for-adaptive-systems-method-of-averaging-and-per/)
- [Widrow et al. 1975](https://ieeexplore.ieee.org/document/1451965)
- [Gibson et al. 2013](https://doaj.org/article/95064d112d7e497c9c7a41f5477e341b), also [arXiv 1304.7278](https://ar5iv.labs.arxiv.org/html/1304.7278)
- [Rohrs et al. (LIDS P-1240)](https://dspace.mit.edu/bitstream/handle/1721.1/2820/P-1240-15681714.pdf)
- [Sanders–Verhulst–Murdock (table of contents)](https://psnz.umt.edu.my/seal/Record/978-0-387-48918-6/TOC)
- [Meerkov 1980 (cited in)](https://link.springer.com/10.1007/s11071-019-05361-3)
- [DCAL origin (Yao)](https://engineering.purdue.edu/~byao/Papers/JDSMC_09_DCARC.pdf)
- From memory: Anderson et al. (MIT Press 1986); Sastry–Bodson 1989; Ioannou–Sun 1996.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bc-1 | major | (iii)'s causal chain (resonance → drift → leaving) and the design rule are unsupported. Leaving is a sudden ejection: \|Δc\| ≤ 0.023 for 131 periods, then 0.82. Single resonant tones at 2× the rms never eject (max \|Δc\| ≤ 0.069, w₁ to 1.60). The 16 tones, 5 near-resonant tones, or one pair with Δω ≈ 18ω_ref (phase-dependent) do. With the band ≥ 1.6√K_max, m = 12 still ejects at s = 0.08 while m = 3 never does (max \|Δc\| 0.022 against 0.14). A far pair (37.8/41.4) ejects at m = 12 (s = 0.04), never at m = 3 (D60.md:25–29) | Restate (iii) as numerical. Ejection is a multi-tone first-passage event. Resonance lowers m = 12's threshold (0.08 → 0.02 when the band covers √K_max), and the mobility Gmw separates m = 12 from m = 3 even off resonance. The rule is necessary-ish, not sufficient. Add Part C controls: single tones, the band ×2.5, m between 3 and 12 (`c_job` options) |
| D3bc-2 | major | Thresholds are single-seed, horizon-bound, and not monotone. With seeds 2 and 3, kernel 1 leaves at 0.015 (m = 12, k_e = 10; reported 0.02). Kernel 3 leaves at 0.02 under seeds 1 and 2 (reported 0.025). m = 3, k_e = 5, s = 0.06: kernel 3 leaves only under seed 7. Kernel 2 (m = 12, k_e = 5) leaves at 0.015–0.04 and 0.08, not at 0.06. "0.08" at k_e = 20 rests on a leave at period 299 of 300 (md:86–95; `thresholds`, exp_noise_resonance.py:405) | At least 4 seeds per cell; report the fraction left and the median period; say "within 300 periods"; flag non-monotone cells |
| D3bc-3 | major | The averaging-failure condition is wrong, and the KB appeal is invalid. Under WFR the 16-tone error is all cross terms: averaged/Σ single tones = 0.9998, 1.0000, 0.9997. Cross terms are 38–56% of the pair diagonal where Δω ≈ lω_ref (3.6004, 2.0029), 3% for the slowest beat (0.057), ≤ 0.3% for half-integer differences. The mean state is not slow (k_e, √K ≈ the tones) (D60.md:22–23, 41–44) | State: the averaged second order is the response to the DC parts; its error is the period-mean response to the 2ω_k and ω_j ± ω_k parts. Large where those are near multiples of ω_ref after the orbit's modulation, or where K(t) crosses ω_k² transversally (31% against 3% at 13.1). Drop "slow beats" and KB |
| D3bc-4 | major | "Odd part … a bounded offset left by the noise's start" is false. Half-integer tones: the centres' period means alternate ±0.309·s (13.1, kernel 1) and ±0.030·s (21.3). Commensurate 21.4: ⟨w₁⟩/s = −0.0011 → −0.188 over 60 periods, above the even part at s ≤ 0.006. The 16 tones (16.0033 ≈ 80ω_ref + 0.0033): a 60-period beat (−0.004, 0.028, 0.096, −0.022, 0.062) (D60.md:20–21) | Per-period Floquet increments rotating by ω_kT; bounded by about 1/\|sin(ω_kT/2)\|, secular at commensurability. Add the condition to (ii) and the Status; disclose the tone choice |
| D3bc-5 | major | "⟨θ₂⟩ grows with time" fails beyond 30 periods. ⟨w₂⟩ saturates (τ ≈ 130 periods, asymptote ≈ 91 ≈ 5× the quoted 18.1); ⟨c₂⟩ is constant from period 3. The 16-tone ⟨c₂⟩ fluctuates between −11 and −98. "Drift" is not Ioannou–Sun's (D60.md:17–21, 48, 70–73) | "A bounded second-order offset reached with the slow coefficient mode's time constant." Quote the saturated offsets (period 300) and the 16-tone ranges |
| D3bc-6 | major | Novelty overstated. The four sources are D‴.31(iv)(c) (PARADIGM.md:3854–3866) and D‴.33(i) (4169); the k_e sign change is D‴.33(iv) (4196). The oscillator with gain-set frequency and k_e damping is classical (notch canceller; MRAC with closed-loop reference models) (D60.md:50–58) | Cite D‴.31(iv)(c), D‴.33, Widrow/Glover, Gibson–Annaswamy–Lavretsky, Riedle–Kokotović and passage through resonance. Limit "New" to the exact non-averaged second order for tones and the paradigm's mobility Gmw |
| D3bc-7 | minor | "Rises with k_e, the oscillator's damping" is a partial attribution. Off resonance (m = 3, 41.3), ⟨w₂⟩ = 12.65, 8.61, 3.86, 0.38 at k_e = 5–40 tracks the rectification factor (−0.97 … −0.03). At 13.1 the sign reverses (+10.9 → −1.57), visible in Part C's raw (m = 3, k_e = 40: w_end 0.869, −0.579, 0.722) (D60.md:28) | Name both effects and cite D‴.33(iv); "raise k_e" reverses, does not cancel, the coefficient drift |
| D3bc-8 | minor | The WFR threshold claim has no D‴.60 data; Part C is paradigm-only (my runs: WFR at k_e = 5 and 10, s ≤ 0.08, max \|Δc\| ≤ 0.023) (D60.md:26–28; `part_C`, :447) | Add WFR rows to C |
| D3bc-9 | minor | (ii) cannot reach (iii)'s levels: the relative error ∝ s² (18.8% at s = 0.005, 13.1), about 75% at 0.01 (D60.md:48–49) | Say so in the Status |
| D3bc-10 | minor | Per-kernel resonances (15.76, 10.51, 13.13) are not stated, and "peaks at ω = 13.9" holds for kernel 1 only. The metric's jitter effect is mostly mobility (×10.8, 7.2, 9.0 at 37.1), with resonance adding ×1.24–1.38 (D60.md:14–15, 65–67) | Print the per-kernel √K_max (`orbit_terms`); reword |
| D3bc-11 | minor | 19 of 42 mutants survive both checks: all of `b3_job` (B4 makes the ratios ×4), `b2_jobs`, and the averaged loop's ē terms (tests at ē = 0) (tests:91–122, 280–296) | Small tests: `b3_job` with NP_B = 2 and S_B3 = (0.002, 0.004) monkeypatched; `b2_jobs` labels; corr and blur at ē ≠ 0; the regressor point in Qw and Qc |
| D3bc-12 | minor | Hand-computed or unprinted numbers: 1.0%, 0.58%, 0.12%, 248, 8%, 3%, the B2 and WFR ranges, and "likewise" for the centres (D60.md:37, 46–47, 69–72) | Print the per-metric ranges, K_max, and the centres' ratios at every s |
| D3bc-13 | minor | The md's B column "averaged / loop" includes the loop's odd part: WFR 13.1, kernel 2 reads 1.353, but 1.001 against the even part (md:45; report :535) | Compare with the even part, or drop the column |
| D3bc-14 | cosmetic | The "mobilities' correlation" and the blur's w–c term are numerically empty (≤ 5·10⁻⁴); the design bound (≈ 4.8) is not evaluated; the odd-part remainder for the 16 tones is not cubic (×14.5, ×20.6 per doubling) (D60.md:20, 29, 73) | Say "three effective sources"; give m ≲ 4.8; reword the odd-part sentence |

## Lesson

A second order verified at small s over 30 periods cannot carry a threshold measured at 2–16 times that level over 300 periods. The cheap controls that separate mechanisms should come before the causal sentence:
- one tone against many at equal rms;
- a band moved off the resonance;
- other seeds and phases;
- a commensurate against a half-integer tone.

With the C loop these take seconds. Here they show that the ejection is a multi-tone event and that the mobility matters off resonance. In a loop driven along a periodic orbit, "fast" tones interact with the orbit's harmonics: commensurability, not slowness, decides what averaging drops, whether the odd part stays bounded, and which tone pairs push a kernel out.
