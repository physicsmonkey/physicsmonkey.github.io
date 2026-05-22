# Calculation 1: O(N) variability of the alpha=2 weight-4 sum

This folder contains a numerical test of how much
$$
\sum_{i<j<k<l} J_{ijkl}^4
$$
can be changed by an orthogonal rotation of a random antisymmetric SYK-like
4-tensor $J$.

The script `rotate_l4.py` samples the independent tensor components as standard
Gaussians. The overall SYK normalization is not included because only ratios
under rotation are reported. It then performs a Jacobi-style coordinate search
over Givens rotations, separately trying to minimize and maximize the quartic
sum. The search preserves the $O(N)$ constraint exactly at every step.

Run from this folder with, for example,

```bash
../../../sparse-spin-syk-teleport/venv/bin/python rotate_l4.py
```

The default run writes:

- `outputs/givens_search_results.csv`: one row per disorder realization.
- `outputs/summary.json`: means and standard deviations grouped by $N$.
- `outputs/ratio_plot.png`: a compact plot of the main ratios.

Status: implemented and run once with the default parameters
`N = 8, 10, 12, 14`, 8 disorder realizations per $N$, and 8 sweeps of the
Givens search. The optimized minimum was consistently below the initial value
but remained well above the equal-component lower bound, while the optimized
maximum showed a larger upward variation.

Agent contribution: Codex (GPT-5) created this calculation folder, wrote and ran
`rotate_l4.py`, generated the output files, and edited the main entry `ngm.md`
to summarize the results on May 22, 2026.
