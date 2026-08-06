# Calculation 2: outside fermions versus the AP propagator

This calculation refines calculation 1 by resolving both sides of the proposed
dictionary into boundary-condition sectors. It compares
$$
C_{\mathrm{out}}(\tau,\tau')
=\frac{1}{N-W}\sum_{i\notin A}\mathbb E_J[
\zeta_{Ai}(\tau)\zeta_{Ai}(\tau')]
$$
with the conditional antiperiodic propagator $G^{\mathrm{AP}}_{12}$.

At the self-consistent connected saddle,
$$
G^{\mathrm P}
=-\left(D_{\mathrm P}-\Delta\tau^2\Sigma\right)^{-T},
\qquad
G^{\mathrm{AP}}
=-\left(D_{\mathrm{AP}}-\Delta\tau^2\Sigma\right)^{-T},
$$
and the original saddle field is recovered as
$G=wG^{\mathrm P}+(1-w)G^{\mathrm{AP}}$. The conditional propagators share
the same $\Sigma=J^2G^3$ and are not independently extremized.

`compare_outside_ap.py` reads the raw $\zeta$ samples and converged saddle from
calculation 1, reconstructs both conditional propagators, checks their weighted
recombination, and makes the sector-resolved comparison. Run it from this
directory with

```bash
python3 compare_outside_ap.py --outdir outputs
```

To recompute the conditional AP object from a connected saddle at a specified
relative weight, use for example

```bash
python3 compare_outside_ap.py --path-weight 0.20 --outdir outputs_w0p20
```

## Initial 24-sample result

For $N=20$, $A=\{1,2,3\}$, $\beta J=0.5$, 24 disorder samples, nine time
points, and an 80-slice saddle:

- A raw one-amplitude fit of $G^{\mathrm{AP}}_{12}$ to $C_{\mathrm{out}}$ has
  cosine similarity $0.99822$ and relative residual $0.0596$.
- Normalizing both objects at $\tau=\tau'=0$ gives cosine similarity $0.99822$
  and relative residual $0.0602$.
- After subtracting the two time averages, the shape correlation is $0.99938$
  and the relative residual is $0.0351$.
- The weighted P/AP reconstruction of the saved saddle has relative error
  $9.97\times10^{-7}$, consistent with the saddle iteration tolerance.

The figure is `outputs/outside_ap_comparison.png`, numerical arrays are in
`outputs/outside_ap_data.npz`, and scalar diagnostics are in
`outputs/summary.json`.

This sector-resolved comparison is substantially stronger than the original
comparison with the weighted $G_{12}$. It shows that the time-independent
mismatch comes from the periodic sector, while the AP object directly captures
the ED correlator for $i\notin A$.

## 100-sample run and common-weight check

The higher-statistics ED source data are in `ed_inputs_100/`. They were
generated with

```bash
python3 ../calculation-1/compare_ed_path_integral.py \
  --samples 100 --outdir ed_inputs_100
```

The literal-background and common-net-weight comparisons are respectively

```bash
python3 compare_outside_ap.py \
  --input ed_inputs_100/comparison_data.npz \
  --outdir outputs_100_w0p15
python3 compare_outside_ap.py \
  --input ed_inputs_100/comparison_data.npz \
  --path-weight 0.20 --outdir outputs_100_w0p20
```

The results are

| path weight | raw residual | origin-normalized residual | centered correlation | centered residual |
| --- | ---: | ---: | ---: | ---: |
| $w=3/20$ | $0.05962$ | $0.06025$ | $0.999384$ | $0.03511$ |
| $w=4/20$ | $0.06325$ | $0.06400$ | $0.999337$ | $0.03640$ |

The central values are unchanged from 24 samples. Jackknife errors on the
shape metrics are only a few parts in $10^6$, because at this temperature the
normalized ED time dependence is nearly sample-independent; disorder mainly
changes its overall amplitude. The remaining few-percent difference is
therefore a small systematic shape mismatch rather than sampling noise.

Using the common net weight $w=4/20$ remains a useful finite-$N$ control, but
the literal background weight $w=3/20$ gives a slightly better AP comparison.

## Annealed moment normalization

There is a small but essential adjustment to the literal proposal to divide
by $\mathbb E[X_A^2]$. Here $|A|=3$ is odd, so
$X_A=\operatorname{Tr}(e^{-\beta H}\mu_A)$ vanishes identically. For each
outside probe the collision instead produces the even string
$$
B_i=A\cup\{i\},\qquad
\psi_i\mu_A=c_i\mu_{B_i},\qquad c_i^2=-\frac12.
$$
The normalization-fixed observable and prediction are therefore
$$
\mathcal C_{\mathrm{out}}(\tau,\tau')
=\frac{\frac1{N-W}\sum_{i\notin A}
\mathbb E_J[Y_{Ai}(\tau)Y_{Ai}(\tau')]}
{\frac1{N-W}\sum_{i\notin A}\mathbb E_J[X_{B_i}^2]},
\qquad
\mathcal C_{\mathrm{out}}\stackrel{?}{=}
-\frac12\frac{G^{\mathrm{AP}}_{12}(\tau,\tau')}
{G^{\mathrm{AP}}_{12}(0,0)},
$$
where $Y_{Ai}=\operatorname{Tr}(e^{-\beta H}\psi_i(\tau)\mu_A)$ and
$X_{B_i}=\operatorname{Tr}(e^{-\beta H}\mu_{B_i})$ are unnormalized traces.
There is no fitted amplitude.

`compare_outside_ap_moment_normalized.py` replays the saved 100-sample
Hamiltonian ensemble to recover the partition functions and $X_{B_i}$. The
shared machinery lives in `cross_replica_moment_normalization.py` in this
folder; it is also loaded by the calculation-3 and calculation-5
moment-normalized comparisons.
At the collision it verifies the saved trace data and the exact value
$\mathcal C_{\mathrm{out}}(0,0)=-1/2$ to machine precision. The conditional
saddle is sampled on aligned $L_\tau=90,180$ grids and extrapolated as
$2G_{180}-G_{90}$.

| path weight | residual away from collision | cosine similarity | maximum difference |
| --- | ---: | ---: | ---: |
| $w=3/20$ | $0.05070$ | $0.998854$ | $0.01487$ |
| $w=4/20$ | $0.05516$ | $0.998645$ | $0.01619$ |

Thus moment normalization fixes the sign and absolute scale and modestly
improves the previous $6.02\%$ origin-normalized residual to $5.07\%$ for the
literal weight. The remaining smooth difference is far larger than the
jackknife error and is a genuine finite-$N$ shape discrepancy in this test;
identifying it specifically as a $1/N$ correction requires a controlled size
sequence.

Run with

```bash
python3 compare_outside_ap_moment_normalized.py
```

The outputs are in `outputs_100_moment_normalized/`.

Status: implemented and run by Codex (GPT-5) on July 29, 2026. The journal
entry itself was not modified.
