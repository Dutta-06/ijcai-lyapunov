import numpy as np
from scipy import stats
import os

print("--- Statistical Significance Test (Steady State Error) ---")

# Load steady-state data (last 40% of the trajectory, same as box plot)
def load_ss(path):
    e = np.load(path)
    return e[int(0.6 * len(e)):]

proposed = load_ss("results/proposed_error.npy")
exact_gp = load_ss("results/exact_gp_error.npy")
sparse_gp = load_ss("results/sparse_gp_fixed_error.npy")

print(f"Proposed mean:  {proposed.mean():.4f}")
print(f"Exact GP mean:  {exact_gp.mean():.4f}")
print(f"Sparse GP mean: {sparse_gp.mean():.4f}\n")

# The data is highly non-normal (bounded at 0, log-distributed peaks).
# We use the non-parametric Wilcoxon rank-sum test (Mann-Whitney U)
# to compare two independent samples, or Wilcoxon signed-rank test if paired.
# Since they are simulated over the same time steps, they are paired.
print("Using Wilcoxon signed-rank test (paired non-parametric):")

stat1, p1 = stats.wilcoxon(proposed, exact_gp)
print(f"Proposed vs Exact GP: p-value = {p1:.4e}")
if p1 < 0.05:
    print("  -> Statistically DIFFERENT (reject null)")
else:
    print("  -> Statistically INSIGNIFICANT difference (fail to reject null)")

stat2, p2 = stats.wilcoxon(proposed, sparse_gp)
print(f"Proposed vs Sparse GP: p-value = {p2:.4e}")
if p2 < 0.05:
    print("  -> Statistically DIFFERENT (reject null)")
else:
    print("  -> Statistically INSIGNIFICANT difference (fail to reject null)")
