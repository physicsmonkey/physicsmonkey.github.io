# Note: periodic versus antiperiodic sectors

The `WeightedReplicas` saddle equation can be resolved before averaging over
fermion boundary conditions. At a fixed self-consistent self-energy,
$$
G^{\mathrm P}
=-\left(D_{\mathrm P}-\Delta\tau^2\Sigma\right)^{-T},
\qquad
G^{\mathrm{AP}}
=-\left(D_{\mathrm{AP}}-\Delta\tau^2\Sigma\right)^{-T},
$$
and
$$
G=wG^{\mathrm P}+(1-w)G^{\mathrm{AP}},
\qquad
\Sigma=J^2G^3.
$$
Thus $G^{\mathrm P}$ and $G^{\mathrm{AP}}$ are conditional propagators that
share the saddle self-energy; they are not independent saddle solutions.

For the replica-connected saddle used in calculation 1
($\beta J=0.5$, $w=3/20$, and 80 time slices), the off-diagonal blocks are

| object | time average | centered RMS |
| --- | ---: | ---: |
| $G^{\mathrm P}_{12}$ | $5.8641$ | $0.01238$ |
| $G^{\mathrm{AP}}_{12}$ | $1.68\times10^{-4}$ | $0.01413$ |

After weighting, the periodic and antiperiodic contributions to the constant
part are respectively
$$
w\,\overline{G^{\mathrm P}_{12}}=0.8796,\qquad
(1-w)\,\overline{G^{\mathrm{AP}}_{12}}=4.95\times10^{-4}.
$$
The large constant that spoiled the original comparison therefore comes
almost entirely from the periodic sector.

There is a corresponding exact-diagonalization split,
$$
C_{\mathrm{in}}(\tau,\tau')
=\frac{1}{W}\sum_{i\in A}\mathbb E_J[
\zeta_{Ai}(\tau)\zeta_{Ai}(\tau')],
$$
$$
C_{\mathrm{out}}(\tau,\tau')
=\frac{1}{N-W}\sum_{i\notin A}\mathbb E_J[
\zeta_{Ai}(\tau)\zeta_{Ai}(\tau')],
$$
with $C_{\mathrm{ED}}=wC_{\mathrm{in}}+(1-w)C_{\mathrm{out}}$. The boundary
condition dictionary is $i\in A\leftrightarrow\mathrm P$ and
$i\notin A\leftrightarrow\mathrm{AP}$.

The saved ED data show that $C_{\mathrm{out}}$ has centered RMS
$2.13\times10^{-6}$, whereas $C_{\mathrm{in}}$ has centered RMS only
$1.24\times10^{-12}$. This large hierarchy is natural at small $W$: for
$i\notin A$, $\psi_i\mu_A$ is a weight-4 string and couples directly to the
$q=4$ Hamiltonian at first order in $\beta J$; for $i\in A$ it is a weight-2
string and first appears at higher perturbative order.

Most importantly, $C_{\mathrm{out}}$ agrees directly with
$G^{\mathrm{AP}}_{12}$. A one-amplitude raw fit has cosine similarity
$0.99822$ and relative residual $0.0596$; after subtracting time averages the
shape correlation is $0.99938$ and the residual is $0.0351$. Calculation 2
records this comparison separately.

For clarity, the exactly replica-diagonal saddle discussed after the original
run has $G_{12}=0$ in both conditional sectors. The numbers above resolve the
nonzero replica-connected saddle used in calculation 1.

Status: written by Codex (GPT-5) on July 29, 2026, from the saved calculation-1
data and a direct reconstruction of the conditional propagators.
