# POST HOC: the sliding methods without the early stop of the path

Made after the test result of the pre-registered run (`blobcount.md`, verdict KILL). The path of `ours` and `greedy` was cut short by a stop that was fixed in the plan and not tuned (a path stops a few entries after an entrant of matched-filter significance below 2); on crowded patches the cut path ended several units short of the true count (see `blobcount_posthoc.md`, column "MAE of the last block"). Here the same 24-configuration grids are run with that stop removed (the entry cap of 120 stays), re-tuned on the development patches with the same objective, frozen, and run once on the test patches. Everything here is POST HOC: it was not pre-registered and cannot support a claim without new held-out data.

## Development (objective: mean of the base and crowded MAE)

| method | best index | objective | MAE base | MAE crowded | bias base | bias crowded |
|---|---|---|---|---|---|---|
| ours_nc | 21 | 4.291 | 4.175 | 4.407 | +3.325 | -4.185 |
| greedy_nc | 9 | 4.062 | 4.458 | 3.667 | +3.858 | -1.815 |

## Test (the frozen configurations of the variants)

| method | stratum | count MAE | relative error | bias | F1 |
|---|---|---|---|---|---|
| ours (pre-registered) | all | 4.833 | 1.998 | +4.040 | 0.577 |
| greedy (pre-registered) | all | 4.663 | 1.882 | +3.750 | 0.572 |
| ours_nc | all | 4.763 | 1.938 | +3.597 | 0.588 |
| greedy_nc | all | 4.737 | 1.906 | +3.823 | 0.570 |
| log (best classical) | all | 2.830 | 0.733 | +2.277 | 0.637 |
| ours (pre-registered) | crowded | 7.853 | 0.349 | -6.029 | 0.650 |
| greedy (pre-registered) | crowded | 7.559 | 0.337 | -6.118 | 0.707 |
| ours_nc | crowded | 8.279 | 0.371 | -7.721 | 0.658 |
| greedy_nc | crowded | 7.559 | 0.337 | -6.118 | 0.707 |
| plm (best classical) | crowded | 2.809 | 0.139 | -0.338 | 0.825 |

Paired bootstrap over images of the difference in count MAE (variant minus other):

| comparison | crowded | all |
|---|---|---|
| ours_nc - log | +5.162 [+1.852, +8.714] | +1.933 [+1.010, +2.777] |
| ours_nc - dog | +5.309 [+1.809, +9.167] | +0.770 [-0.743, +2.117] |
| ours_nc - plm | +5.471 [+2.034, +9.342] | +0.937 [-0.443, +2.137] |
| ours_nc - ws | +5.294 [+1.833, +9.350] | +1.370 [+0.360, +2.290] |
| greedy_nc - log | +4.441 [+1.508, +7.737] | +1.907 [+1.063, +2.683] |
| greedy_nc - dog | +4.588 [+1.485, +8.215] | +0.743 [-0.723, +2.060] |
| greedy_nc - plm | +4.750 [+1.662, +8.407] | +0.910 [-0.403, +2.067] |
| greedy_nc - ws | +4.574 [+1.492, +8.431] | +1.343 [+0.393, +2.217] |
| ours_nc - greedy_nc | +0.721 [-0.265, +1.529] | +0.027 [-0.437, +0.490] |

The pre-registered rule with `ours_nc` and `greedy_nc` in place of `ours` and `greedy` (best classical: crowded plm, all log): C1 not met, C2 not met, C3 not met, rule outcome **KILL** (POST HOC, not the verdict).

Crowded MAE: ours_nc 8.279, best classical 2.809, greedy_nc 7.559; all-patch MAE: ours_nc 4.763, best classical 2.830. Share of paths whose BIC block is the last one: ours_nc 0.096, greedy_nc 0.032.

* ours_nc (grid index 21, development objective 4.291): `{"path": {"box": 1.0, "extra": 6, "fit": 3.0, "h": 1.0, "hi": 4.0, "hp": 0.0, "lo": 1.0, "margin": 1.0, "max_entries": 120, "rex": 0.8, "rtrunc": 4.0, "s0": 1.8, "snr_min": null}, "stop": {"cal": 1, "clean": 1.0, "lam": 64.0}}`
* greedy_nc (grid index 9, development objective 4.062): `{"path": {"box": 1.0, "extra": 6, "fit": 3.0, "h": 1.0, "hi": 3.0, "hp": 0.0, "lo": 0.7, "margin": 1.0, "max_entries": 120, "rex": 0.8, "rtrunc": 4.0, "s0": 1.3, "snr_min": null}, "stop": {"cal": 1, "clean": 1.0, "lam": 64.0}}`
