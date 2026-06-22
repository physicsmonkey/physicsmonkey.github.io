# Finite-K oscillator numerics

This folder supports the numerical calculation in `../cvsd.md`.

`finite_k_numerics.py` computes the finite-cutoff oscillator quantities
$$A_1=|\langle 0|e^{-sH}(a^\dagger)^{4K}|0\rangle|^2$$
and
$$A_2=\langle 0|a^{4K}e^{-2sH}(a^\dagger)^{4K}|0\rangle$$
for
$$H=\omega a^\dagger a+\frac{\gamma}{K}(a+a^\dagger)^4.$$
The calculation restricts to the even Fock sector, diagonalizes the truncated
Hamiltonian using NumPy/SciPy, and finds the first crossing of
`A2 exp(-K delta) / A1` as `gamma` is increased.

The default run uses `omega=s=delta=1`, `K=1,...,16`, and checks cutoff
convergence at `nmax=80,120,160,220`. Outputs are written to `outputs/`.
The main comparison figure is `figures/finite_k_crossings.svg`, and a cutoff
check is in `figures/cutoff_convergence.svg`. The file
`figures/k10_options_vs_gamma.svg` plots the two competing
common-factor-stripped weights for `K=10`; the underlying data are in
`outputs/k10_options.csv`.
