# Diagonal-resummed perturbation theory

This folder tests an analytic improvement over the bare perturbative estimate in
`../cvsd.md`. The idea is to treat the diagonal part of the quartic interaction
exactly and then keep the leading monotone `a^4` lowering process into each
`4`-spaced Fock sector.

For the state `|4q>`, the diagonal energy is
$$
E_q=4\omega q+\frac{\gamma}{K}(96q^2+24q+3),
$$
because
$$
\langle n|(a+a^\dagger)^4|n\rangle=6n^2+6n+3.
$$
The common-factor-stripped amplitude to end in `|4(K-j)>` is approximated by
$$
d_j =
\left(\frac{\gamma}{K}\right)^j
\sqrt{\frac{(4K)!}{(4K-4j)!}}\,
[E_{K-j},\ldots,E_K]\,e^{-s x},
$$
where the bracket denotes a divided difference of the function `exp(-s x)`.
The approximation then compares
$$
A_1^{\rm diag}=|d_K|^2,\qquad
A_2^{\rm diag}=\sum_{j=0}^K |d_j|^2.
$$

The script `diagonal_resummation.py` evaluates this formula at
`omega=s=delta=1`, finds the crossing of
`A2_diag exp(-K delta) / A1_diag`, and compares it to the bare perturbative
finite-`K` prediction and the refined envelope data from `calculation-2`.

## Outputs

- `outputs/diagonal_resummed_summary.csv` -- diagonal-resummed and bare
  perturbative crossovers for `K=1..40,50,60,80`.
- `outputs/comparison_to_refined.csv` -- comparison against the `calculation-2`
  envelope data where available.
- `figures/diagonal_resummed_vs_K.svg` -- main comparison plot.
- `figures/diagonal_resummed_K15_rates.svg` -- the two diagonal-resummed weights
  versus `gamma` at `K=15`.

## Status

Self-contained. Run `python diagonal_resummation.py` from this folder or from
the repository root. This approximation captures a large correction to the bare
large-`K` result but appears to undershoot the refined envelope crossover,
suggesting that the diagonal self-energy is important but not the whole finite
`K` story.
