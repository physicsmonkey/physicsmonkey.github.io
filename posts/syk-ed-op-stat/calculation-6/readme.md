# Calculation 6: exact P/AP boundary-update relation

This calculation turns the proposed Woodbury explanation of the collision
product into an exact finite-time-lattice identity.

Let

$$
A_s=D^s-\Delta\tau^2\Sigma,
\qquad G^s=-A_s^{-T},
$$

for the two-replica conditional kernels. The forward-difference derivatives
differ only at the temporal corner. If $E$ selects the first time site on
each replica and $F$ the last, then

$$
A_{\mathrm{AP}}=A_{\mathrm P}+2EF^T.
$$

Woodbury therefore gives

$$
G^{\mathrm{AP}}
=G^{\mathrm P}
+2G^{\mathrm P}F(1-2E^TG^{\mathrm P}F)^{-1}E^TG^{\mathrm P}.
$$

For the $2\times2$ replica boundary matrices

$$
C_{\mathrm P}=E^TG^{\mathrm P}E,\qquad
K_{\mathrm P}=E^TG^{\mathrm P}F,\qquad
C_{\mathrm{AP}}=E^TG^{\mathrm{AP}}E,
$$

this reduces to the exact identity

$$
C_{\mathrm{AP}}=(1-2K_{\mathrm P})^{-1}C_{\mathrm P}.
$$

Write the replica-circulant matrices as

$$
C_{\mathrm P}=\begin{pmatrix}k&-m\\m&k\end{pmatrix},\qquad
K_{\mathrm P}=\begin{pmatrix}a&-b\\b&a\end{pmatrix}.
$$

The AP off-diagonal collision is then exactly

$$
m_{\mathrm{AP}}
=\frac{(1-2a)m+2bk}{(1-2a)^2+4b^2}.
$$

In the continuum, the canonical one-sided contact and periodic boundary
condition give

$$
k\to-\frac12,\qquad a\to+\frac12,\qquad b\to m.
$$

Consequently, on the replica-connected branch ($m\ne0$),

$$
m_{\mathrm{AP}}=-\frac{1}{4m},\qquad
G^{\mathrm P}_{12}(0,0)G^{\mathrm{AP}}_{12}(0,0)=-\frac14.
$$

Thus the product is an exact continuum boundary-condition identity, not only
the leading weak-coupling zero-mode result. The zero-mode argument remains
useful because it explains why the P collision is large and the AP collision
small. The small deviations from $-1/4$ after two-grid Richardson
extrapolation are residual time-discretization errors.

`derive_boundary_woodbury.py` checks the finite-lattice matrix identity using
the saved calculation-4 and calculation-5 saddles at
$(N,\beta J,w)=(20,0.5,0.2)$ and $(24,0.5,1/3),(24,1,1/3)$. The maximum
Woodbury error is at the level of floating-point inversion error. It also
records the approach of $k,a,b$ to their continuum boundary values and the
Richardson-extrapolated collision products.

Run from this directory with

```bash
python3 derive_boundary_woodbury.py
```

The numerical record is `outputs/woodbury_check.json`.

Status: derived, implemented, and checked by Codex (GPT-5.6 Sol) on
August 3, 2026. The entry was edited with the result.
