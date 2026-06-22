# Robust (zero-insensitive) crossover

This folder supports a follow-up to the numerical section of `../cvsd.md`. It
addresses why the `calculation-1` "exact first crossing" fails to converge for
`K >= 13` and replaces it with a crossover diagnostic that is well defined.

## Why the first crossing is not a good observable

The disconnected weight is `A1 = |B|^2` with
`B(gamma) = <0|e^{-sH}|4K>`. Two facts, both checked here:

1. **`B` has genuine zeros.** `B` is not sign definite; competing multi-step
   paths cancel and `B` vanishes at a discrete set of `gamma` (verified at
   60-digit precision). At those points `log A1 -> -inf`, so the first crossing
   of `A2 e^{-K delta}/A1` is set by the nearest zero, not by the smooth
   competition of the two weights.

2. **Double precision cannot reach the crossover region for `K >= 12`.** In the
   relevant region `|B|` falls below ~`1e-21` and the value returned by a
   double-precision eigendecomposition is deterministic round-off noise — it
   reproduces across `nmax` (same banded round-off) yet is wrong in magnitude
   and sometimes in sign. Increasing `nmax` does not help because it does not
   raise the floating-point floor. This is the real reason the
   `calculation-1` first crossings wander for `K >= 13`.

## What this calculation does

`robust_crossover.py` (NumPy/SciPy + mpmath) computes `A1` and
`A2 = <4K|e^{-2sH}|4K>` at `s=1` in **extended precision** (mpmath `eigsy`,
40 digits), so the small-`gamma` weights are the true values. It then defines
the crossover from the **smooth envelope** of `A1`: writing
`B = E(gamma) * oscillation`, the local maxima of `|B|` trace the envelope `E`
(the zeros are where the oscillation vanishes). It interpolates `log A1` through
its peaks and crosses it with the smooth `log(A2 e^{-K delta})` to get a single
`gamma_*(K)`.

Cross-checks, all reaching the same qualitative conclusion:

- `window`  -- a multiplicative running average of `A1` in `gamma`.
- `s_avg`   -- average of `A1(s), A2(s)` over `s in [0.7,1.3]` (double
  precision). The `B(gamma,s)` zeros are isolated in `s`, so the `s`-integral is
  smooth. Because `A1(s), A2(s)` grow as `s` decreases, this integral is
  weighted toward small `s` and reports a crossover at a somewhat *larger*
  `gamma` than the true `s=1` value; it is used only to confirm the trend.
- `first`   -- the original `calculation-1` double-precision first crossing,
  for contrast.

Defaults: `omega = s = delta = 1`, `K = 10..18`, `nmax = 90` for the
extended-precision scan, with cutoff-convergence spot checks at `K=14`
(`nmax = 80,120,160`) and `K=16,18` (`nmax = 90,140`). The envelope crossover
is grid-resolution limited: refining the `gamma` grid moves it by about one
percent, so the crossover is reliable to roughly two significant figures
(`gamma_* ~ 0.039`), not to the digits printed in the raw output CSVs.

## Outputs

- `outputs/crossover_summary.csv` -- `gamma_*` per `K` for all four diagnostics.
- `outputs/convergence_check.csv` -- envelope `gamma_*` vs `nmax` at `K=14`.
- `outputs/diagnostic_K15.csv`    -- the `log A1`, envelope, and
  `log(A2 e^{-K delta})` curves used in the illustration figure.
- `figures/robust_crossover_vs_K.svg` -- `gamma_*(K)` for the four diagnostics
  and the large-`K` line.
- `figures/diagnostic_K15.svg`        -- the envelope construction at `K=15`.

## Status

Self-contained; run `python robust_crossover.py` from this folder. The headline
finding is that the envelope crossover is cutoff-converged and, unlike the
`calculation-1` first crossing, **does not drift downward** with `K` — it sits
near `gamma_* ~ 0.04` over `K = 10..18`. See
`../reports/report-2026-06-22-claude-opus-4-8.md` (addendum) for discussion.
