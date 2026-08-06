# Calculation 4: normalization-fixed same-side correlators

This calculation uses the even parent string
$$
A=\{1,2,3,4\}
$$
to test the same-replica conditional propagators without introducing an odd
background. Define
$$
X_A=\operatorname{Tr}(e^{-\beta H}\mu_A)
$$
and
$$
X^{(11)}_{A;i}(\tau,\tau')
=\operatorname{Tr}\!\left[
e^{-\beta H}T_\tau\psi_i(\tau)\psi_i(\tau')\mu_A
\right].
$$
The ED observable is the annealed, moment-normalized ratio
$$
\mathcal C^{(s)}_{11}(\tau,\tau')
=\frac{\mathbb E_J[
X^{(11)}_{A;i}(\tau,\tau')X_A]}
{\mathbb E_J[X_A^2]},
$$
averaged over $i\in A$ for the P sector and over $i\notin A$ for the AP
sector. This is the direct two-replica moment derivative, so its absolute
normalization is predicted:
$$
\mathcal C^{\mathrm P}_{11}\stackrel{?}{=}G^{\mathrm P}_{11},
\qquad
\mathcal C^{\mathrm{AP}}_{11}\stackrel{?}{=}G^{\mathrm{AP}}_{11}.
$$
No amplitude is fitted.

More explicitly, introduce a bilocal source $K_{11}$ for the probe in the
first trace and write
$$
M_A[K_{11}]=\mathbb E_J[X_A[K_{11}]X_A].
$$
Then, with the source convention used in the saddle equations,
$$
\left.\frac{\delta\log M_A}{\delta K_{11}(\tau,\tau')}\right|_{K=0}
=\frac{\mathbb E_J[X^{(11)}_{A;i}(\tau,\tau')X_A]}
{\mathbb E_J[X_A^2]}.
$$
This fixes both the denominator and the absolute coefficient; it is not a
normalization chosen after examining the curves.

The matrix discretization uses the negative one-sided contact convention,
so the diagonal prediction is $-1/2$. Off the diagonal, ordinary fermionic
time ordering fixes the antisymmetric continuation.

`compare_same_side.py` evaluates the three-operator thermal trace exactly
using parity blocks, spectral propagation, and sparse Majorana conjugation.
Delete-one-disorder-sample jackknife errors include the correlated numerator
and denominator.

The saddle uses time lattices with $L_\tau=90$ and $180$, both exactly aligned
with the nine ED times. Since the forward-difference derivative has a leading
$1/L_\tau$ error, the reported saddle surface is the parameter-free
Richardson estimate
$$
G_{\rm cont}=2G_{180}-G_{90}.
$$
This extrapolation matters mainly for the P zero-mode sector: its
off-diagonal ED residual falls from $2.37\%$ at $L_\tau=90$ and $1.18\%$ at
$L_\tau=180$ to $0.0529\%$. Repeating the estimate with $L_\tau=180,360$
gives $0.0534\%$, so the remaining difference is not a time-lattice artifact
at this resolution.

The latter check is reproduced by `check_path_extrapolation.py`; its scalar
convergence record is retained in
`outputs/path_extrapolation_180_360.json`.

## Result

For $N=20$, $\beta J=0.5$, $W=4$, 24 disorder samples, and nine time points:

| comparison | off-diagonal residual | off-diagonal cosine | contact: ED / saddle |
|---|---:|---:|---:|
| $i\in A$ vs. $G^{\rm P}_{11}$ | $5.29\times10^{-4}$ | $0.999999945$ | $-0.500000/-0.499992$ |
| $i\notin A$ vs. $G^{\rm AP}_{11}$ | $1.56\times10^{-3}$ | $0.999999905$ | $-0.500000/-0.5000001$ |

There is no fitted scale, offset, or time shift in either row. Thus the
same-side experiment fixes the normalization as well as the shape and gives
a direct empirical identification of the inside probe with the P conditional
propagator and the outside probe with the AP conditional propagator for even
$|A|$. The small P difference is statistically visible (the norm of the
difference is about four times the norm of the jackknife-error surface), so it
is a resolved finite-$N$ discrepancy rather than disorder noise. Its scaling
with $N$ has not been tested.

Run from this directory with

```bash
python3 compare_same_side.py --outdir outputs
```

Outputs:

- `outputs/summary.json`: parameters and comparison metrics, including both
  finite time lattices.
- `outputs/same_side_data.npz`: ED samples, normalized surfaces, jackknife
  errors, saddle solutions, and continuum estimates.
- `outputs/same_side_comparison.png`: side-by-side surfaces and differences.

Validation: the spectral trace evaluator was compared with a direct dense
matrix trace at $N=6$; the maximum discrepancy was
$7.2\times10^{-16}$. The reproducible check is `validate_same_side.py`, and
its scalar output is retained in `validation.json`.

Status: implemented and run by Codex (GPT-5) on July 29, 2026. The journal
entry itself was not modified at that time. Reproducible validation and the
180/360 convergence record were added by Codex (GPT 5.6 Sol) on August 3,
2026; the entry was updated with their results.
