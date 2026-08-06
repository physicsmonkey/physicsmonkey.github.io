# Calculation 3: an inside probe versus the P propagator

This calculation tests the periodic sector using
$$
A=\{1,2,3,4,5\},\qquad i=5.
$$
Since $i\in A$, this is the sector expected to correspond to the conditional
periodic propagator $G^{\mathrm P}_{12}$. At coincident time,
$$
\psi_5\mu_A=\frac{1}{\sqrt{2}}\mu_{\{1,2,3,4\}}
$$
up to the fixed string convention, so the net operator has weight four and
couples directly to the $q=4$ Hamiltonian at first order in $\beta J$.

`compare_inside_p.py` performs a new 24-sample, $N=20$ ED run and compares both
the fixed $i=5$ probe and the symmetry-improved average over all five
$i\in A$ with $G^{\mathrm P}_{12}$. It solves connected saddles at two
relative weights:

- $w=5/20=0.25$, the literal weight of the background $\mu_A$;
- $w=4/20=0.20$, a finite-$N$ sensitivity check based on the coincident net
  string weight.

Run from this directory with

```bash
python3 compare_inside_p.py --outdir outputs
```

The output contains raw and origin-normalized figures for both path weights,
all ED and saddle arrays in `outputs/inside_p_data.npz`, and scalar diagnostics
in `outputs/summary.json`.

## Initial 24-sample result

For the primary $w=5/20$ comparison:

- The fixed $i=5$ probe has raw cosine similarity $0.99999997$ and raw
  relative residual $2.61\times10^{-4}$. Its constant component is well
  reproduced, but 24 samples are insufficient to give a stable error estimate
  for its much smaller time-dependent component.
- Averaging over the five equivalent $i\in A$ probes gives origin-normalized
  residual $8.86\times10^{-5}$. More stringently, after removing the dominant
  time-independent component, the shape correlation is
  $0.999201\pm0.000203$ and the relative residual is
  $0.03997\pm0.00491$, with jackknife errors over disorder samples.
- The $w=4/20$ sensitivity run is comparably good: for the inside average its
  centered correlation is $0.99904$ and centered residual is $0.0439$.

The raw comparison is dominated by the expected periodic constant, while the
normalized figures expose the sub-percent time variation. Both components
agree, and the centered agreement is stable only after using the
symmetry-improved average over $i\in A$.

## 100-sample run

The higher-statistics run is

```bash
python3 compare_inside_p.py --samples 100 --outdir outputs_100
```

For the symmetry-improved average over all five $i\in A$:

| path weight | origin-normalized residual | centered correlation | centered residual |
| --- | ---: | ---: | ---: |
| $w=5/20$ | $0.00067\pm0.00129$ | $0.999174\pm0.000123$ | $0.04063\pm0.00303$ |
| $w=4/20$ | $0.00145\pm0.00130$ | $0.998951\pm0.000316$ | $0.04579\pm0.00689$ |

Errors are delete-one-disorder-sample jackknife estimates. The fixed $i=5$
probe is also stable at 100 samples: its centered residual is
$0.03985\pm0.00553$ for $w=5/20$.

Both path weights describe the P-sector ED surface very well. The literal
background weight $w=5/20$ is mildly favored by the central values, but its
advantage over the common net weight $w=4/20$ is not statistically
significant. Reporting both cleanly exposes the finite-$N$ ambiguity.

## Annealed moment normalization

As in calculation 2, the odd background has
$X_A=\operatorname{Tr}(e^{-\beta H}\mu_A)=0$. For an inside probe the
nonzero collision string is instead
$$
B_i=A\setminus\{i\},\qquad
\psi_i\mu_A=c_i\mu_{B_i},\qquad c_i^2=+\frac12.
$$
The no-fit comparison is
$$
\mathcal C_{\mathrm{in}}(\tau,\tau')
=\frac{\frac1W\sum_{i\in A}
\mathbb E_J[Y_{Ai}(\tau)Y_{Ai}(\tau')]}
{\frac1W\sum_{i\in A}\mathbb E_J[X_{B_i}^2]}
\stackrel{?}{=}
\frac12\frac{G^{\mathrm P}_{12}(\tau,\tau')}
{G^{\mathrm P}_{12}(0,0)}.
$$
This uses unnormalized traces inside both disorder averages, rather than the
sample-by-sample thermal normalization in $\zeta_{Ai}=Y_{Ai}/Z$.

`compare_inside_p_moment_normalized.py` replays the existing 100-sample
ensemble, verifies the collision identity with relative error
$7.8\times10^{-15}$, and uses the aligned continuum estimate
$2G_{180}-G_{90}$. The collision value is $+1/2$ on both sides.

| path weight and probe | residual away from collision | cosine similarity | difference/error norm |
| --- | ---: | ---: | ---: |
| $w=5/20$, average over $i\in A$ | $0.000825$ | $0.999999969$ | $0.631$ |
| $w=5/20$, fixed $i=5$ | $0.001886$ | $0.999999829$ | $0.475$ |
| $w=4/20$, average over $i\in A$ | $0.001622$ | $0.999999859$ | $1.240$ |
| $w=4/20$, fixed $i=5$ | $0.001088$ | $0.999999954$ | $0.274$ |

The literal-weight, symmetry-averaged comparison agrees at the
$8.3\times10^{-4}$ level without a fitted scale, and its full difference is
smaller than the jackknife-error norm. This promotes the earlier P-sector
shape agreement to an absolute moment-normalized prediction.

Run with

```bash
python3 compare_inside_p_moment_normalized.py
```

The outputs are in `outputs_100_moment_normalized/`.

Status: implemented by Codex (GPT-5) on July 29, 2026. The journal entry itself
was not modified.
