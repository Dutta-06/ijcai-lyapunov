# Reference: pretrained instance-segmentation models on the same test patches (run after the test result)

Pretrained models with their default parameters, nothing tuned or chosen on any patch, centroids of the instances scored exactly like the detectors of `blobcount.md`. Versions: cellpose 3.1.1.2, stardist 0.9.2. These are references, not targets: such models are reportedly trained on data that include nuclei images of this challenge (no check was made of any overlap with the test images), so their numbers are probably optimistic. Generated from `results/ml/blobcount_ref_raw.npz` and `results/ml/blobcount_raw.npz` (`--from-raw`).

## all

| method | count MAE | relative error | bias (detected - true) | detection F1 (6 px) |
|---|---|---|---|---|
| ours | 4.833 | 1.998 | +4.040 | 0.577 |
| log (best classical) | 2.830 | 0.733 | +2.277 | 0.637 |
| cellpose | 2.037 | 0.434 | -0.903 | 0.729 |
| stardist | 0.997 | 0.232 | -0.150 | 0.841 |

## crowded

| method | count MAE | relative error | bias (detected - true) | detection F1 (6 px) |
|---|---|---|---|---|
| ours | 7.853 | 0.349 | -6.029 | 0.650 |
| plm (best classical) | 2.809 | 0.139 | -0.338 | 0.825 |
| cellpose | 9.882 | 0.428 | -9.588 | 0.619 |
| stardist | 3.279 | 0.150 | -1.574 | 0.871 |

## very crowded

| method | count MAE | relative error | bias (detected - true) | detection F1 (6 px) |
|---|---|---|---|---|
| ours | 21.091 | 0.626 | -21.091 | 0.469 |
| plm (best classical) | 5.000 | 0.152 | -1.182 | 0.821 |
| cellpose | 27.091 | 0.802 | -27.091 | 0.264 |
| stardist | 7.909 | 0.237 | -4.091 | 0.842 |
