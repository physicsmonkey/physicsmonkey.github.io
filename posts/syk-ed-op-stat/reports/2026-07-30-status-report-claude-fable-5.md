# Status report on `syk-ed-op-stat.md`

> Claude Fable 5 (Anthropic), July 30, 2026.
> Prepared at the author's request as a status assessment with emphasis on the
> meaning of the conditional P and AP propagators. Per `journal.md`, this report
> does not modify the entry or the calculations. The supporting scan in this
> folder (`scan_quarter_product.py`, `quarter_product_scan.json`) was run
> independently for this report.

## 1. Verification of quoted results

Every quantitative claim in the entry's "Testing the proposal" section was
checked against the saved outputs in calculations 1–3 and matches exactly:

- AP test (calculation 2, 100 samples, moment-normalized): residual $0.05070$,
  cosine $0.998854$ at $w=3/20$; residual $0.05516$ at $w=4/20$; ED collision
  value $-1/2$ to machine precision.
- P test (calculation 3): residual $8.25\times10^{-4}$, cosine $0.999999969$
  at $w=5/20$; $1.62\times10^{-3}$ at $w=4/20$; difference norm below the
  jackknife-error norm ($0.631$); collision $+1/2$ exact.
- The figure paths referenced by the entry exist.

Two provenance gaps:

1. The collision values $G^{\mathrm P}_{12}(0,0)=4.75227$,
   $G^{\mathrm{AP}}_{12}(0,0)=-0.052608$, and the product
   $-0.250005\simeq-1/4$ appear in the entry but are produced by no script and
   recorded in no output in any calculation folder. The scan in this folder
   independently reproduces all three (at $\beta J=0.5$, $w=0.2$:
   $4.752266$, $-0.052608$, $-0.250005$).
2. Calculation 4 — arguably the strongest evidence, a parameter-free
   absolute-normalization test of the same-side objects
   $\mathcal C^{(s)}_{11}$ against $G^{\mathrm{P/AP}}_{11}$ with residuals
   $5.3\times10^{-4}$ (P) and $1.6\times10^{-3}$ (AP), with the normalization
   *derived* from the source-derivative identity
   $\delta\log M_A/\delta K_{11}$ — is not mentioned in the entry.

## 2. Independent scan of the $-1/4$ product

`scan_quarter_product.py` solves the connected two-replica saddle over
$\beta J\in\{0.25,0.5,1,2\}$, $w\in\{0.1,0.2,0.3\}$ (same solver as
calculation 1; collision values Richardson-extrapolated from aligned
$L_\tau=90,180$ lattices). Results in `quarter_product_scan.json`:

| $\beta J$ | $w$ | $G^{\mathrm P}_{12}(0,0)$ | $G^{\mathrm{AP}}_{12}(0,0)$ | product $+1/4$ |
|---|---|---:|---:|---:|
| 0.25 | 0.1 / 0.2 / 0.3 | 11.267 / 6.696 / 4.938 | −0.0222 / −0.0373 / −0.0506 | $\approx-2\times10^{-6}$ |
| 0.50 | 0.2 / 0.3 | 4.752266 / 3.5010 | −0.052608 / −0.071410 | $-5\times10^{-6}$ |
| 1.00 | 0.3 | 2.5007 | −0.099990 | $-4.3\times10^{-5}$ |

At the remaining scan points ($w=0.1$ for $\beta J\ge0.5$; $w\le0.2$ at
$\beta J=1$; all of $\beta J=2$) the plain damped iteration collapses to the
replica-diagonal branch $G_{12}=0$, so the connected branch was not reached;
testing there requires continuation in $w$ or $\beta J$.

Conclusion: the $-1/4$ product is genuine on the reachable connected branch,
but the deviation grows steadily with coupling (about two orders of magnitude
from $\beta J=0.25$ to $1$). It is a leading-order identity, not an exact one.

## 3. Mechanism: a periodic-sector Grassmann zero mode

The P kernel $D^{\mathrm P}-\Sigma$ is singular as $\Sigma\to0$: each periodic
replica carries a fermionic zero mode $\theta_r$, and the only term pairing
$(\theta_1,\theta_2)$ is the cross-replica self-energy. Integrating the pair
out gives, at leading order in $\Sigma_{12}$,
$$
G^{\mathrm P}_{12}(0,0)\simeq\frac{1}{\iint\Sigma_{12}\,d\tau\,d\tau'},
$$
while the AP off-diagonal has no zero mode and starts at first order,
$$
G^{\mathrm{AP}}_{12}(0,0)
\simeq\big[G^{\mathrm{AP}}_{d}\star\Sigma_{12}\star G^{\mathrm{AP}}_{d}\big](0,0)
\simeq-\frac{\beta^2}{4}\,\bar\sigma_{12},
$$
using $\int_0^\beta G_d(\tau,0)\,d\tau=\beta/2$ for the nearly free diagonal
AP propagator and $\bar\sigma_{12}$ the double-time average of $\Sigma_{12}$.
Both formulas were checked against the converged saddles in the scan: the
zero-mode formula reproduces $G^{\mathrm P}_{12}(0,0)$ to 0.2–1.6% across all
points, and the first-order AP formula matches to $\sim10^{-4}$ relative
(e.g. predicted $-0.052614$ vs. actual $-0.052608$). The product cancels
$\Sigma_{12}$ exactly at this order:
$\bigl(-\tfrac{\beta^2\bar\sigma}{4}\bigr)\cdot\tfrac{1}{\beta^2\bar\sigma}
=-\tfrac14$. Consequently the entry's conjectured normalization has a concrete
identity:
$$
\Lambda=2G^{\mathrm P}_{12}(0,0)\simeq\frac{2}{\iint\Sigma_{12}},
$$
the inverse zero-mode pairing amplitude between the two twisted (periodic)
replicas — a boundary-state normalization in exactly the sense the entry
speculates.

The product is far more accurate than the individual leading-order formulas
(deviation $10^{-6}$ where the zero-mode formula errs at $10^{-2}$), i.e. the
subleading corrections cancel. A likely structural reason: $D^{\mathrm{AP}}$
and $D^{\mathrm P}$ differ by a rank-one corner update per replica, so the
Woodbury identity expresses the collision values of $G^{\mathrm{AP}}$ exactly
in terms of a $2\times2$ replica matrix of corner values of $G^{\mathrm P}$.
Writing that matrix as $\begin{pmatrix}k&m\\-m&k\end{pmatrix}$ with
$m\sim G^{\mathrm P}_{12}(0,0)$ large, the expansion gives
$G^{\mathrm{AP}}_{12}(0,0)=-c/G^{\mathrm P}_{12}(0,0)$ with $c\to1/4$ as the
equal-time contact $k\to-1/2$, and the deviation from $-1/4$ controlled by the
$O(\Sigma)$ shift of the contact — consistent with the observed growth in
$\beta J$. This appears to be the "determinant/Pfaffian/monodromy relation"
the entry asks for, and should be provable in a page.

## 4. Why inside $\leftrightarrow$ P and outside $\leftrightarrow$ AP

The dictionary follows from a one-line operator identity. For
$Y_{Ai}(\tau)=\operatorname{Tr}(e^{-\beta H}\psi_i(\tau)\mu_A)$,
$$
Y_{Ai}(\beta)=\operatorname{Tr}(\psi_i e^{-\beta H}\mu_A)
=\operatorname{Tr}(e^{-\beta H}\mu_A\psi_i),
\qquad
\mu_A\psi_i=
\begin{cases}
(-1)^{|A|}\,\psi_i\mu_A,&i\notin A,\\
(-1)^{|A|-1}\,\psi_i\mu_A,&i\in A,
\end{cases}
$$
because $\psi_i$ anticommutes with every distinct Majorana in $\mu_A$ but
commutes through its own factor. For odd $|A|$ this gives
$Y(\beta)=-Y(0)$ outside (AP) and $Y(\beta)=+Y(0)$ inside (P).

A tempting but incorrect argument runs: (i) the thermal trace gives AP
boundary conditions; (ii) transporting $\psi_i$ ($i\notin A$) around the
circle past the odd string $\mu_A$ gives an extra $(-1)$; hence P for outside
probes. The flaw is in (i): KMS antiperiodicity is not a property of the trace
by itself. A lone fermion insertion in a trace is $\beta$-periodic — cyclicity
carries no sign. The familiar AP of $\langle T\psi(\tau)\psi(0)\rangle$ is
itself the crossing sign past the *partner* insertion $\psi(0)$. In the
cross-replica object $\mathbb E_J[Y_{Ai}(\tau)Y_{Ai}(\tau')]$ the partner
lives on the other replica circle and is never crossed, so the only monodromy
source is $\mu_A$: $(-1)^{|A|}$ outside, $(-1)^{|A|-1}$ inside. Starting from
an AP baseline and then adding the $\mu_A$ crossing double-counts.

In the flavor-factorized Gaussian (fixed-$\Sigma$) description the same rule
reads: spectator flavors decouple by Wick's theorem, so reordering
$\psi_i(\tau)$ past the other-flavor content of $\mu_A$ costs a
$\tau$-independent overall sign that cannot affect the boundary condition;
only same-flavor insertions produce $\tau$-dependent signs. Each $\psi_i$
insertion at a point is a $\mathbb Z_2$ twist for flavor $i$'s kernel:
an odd number of $\psi_i$ insertions on the circle flips AP to P. Outside
probes see zero $\psi_i$ insertions in $\mu_A$ (AP); inside probes see one
(P). The count is uniform across the experiments: for the even-$|A|$
same-side objects of calculation 4 the probe pair $\psi_i(\tau)\psi_i(\tau')$
sits on one circle, contributing one crossing itself, again giving
$i\in A\to$ P, $i\notin A\to$ AP.

## 5. A caution on the entry's framing

The P-sector agreement is described as "considerably stronger"
($8\times10^{-4}$ vs. $5\times10^{-2}$). Numerically true, but unequal: the P
surface is dominated by the normalization-protected zero-mode constant, which
the moment normalization fixes by construction. After centering, the P shape
residual is $\sim\!4\%$ (calculation 3: $0.0406\pm0.0030$) — the same scale as
the AP sector's $5\%$. The honest summary is that both sectors carry a
comparable few-percent finite-$N$ shape correction; the P sector hides it
under a large protected constant.

## 6. Textual and structural items

- Line "sum over the remaining $N-W-1$ fermions": should be $N-W+1$
  (check: $N=20$, $W=4$ gives $17$ probes, as used in calculation 2).
- "the specific example of $A=\{1,2,3,5\}$": the example set up in the text
  was $A=\{1,2,3,4\}$.
- `\sum_{ i in (A \cup 5)}`: `in` should be `\in`.
- "symmetries that constraint" → "constrain"; "there are strong discrete
  since every odd weight vanishes" is missing a noun.
- Header placeholders (Confidence, Reports, Acknowledgements) are unfilled.
- Byline says "Codex with GPT 5.6"; the section footer and calculation
  readmes say "Codex (GPT-5)".

## 7. Suggested next steps

1. Add calculation 4 to the entry; its source-derivative normalization
   argument partially supplies the missing derivation of the dictionary.
2. Record provenance for the collision values and the $-1/4$ product
   (the scan in this folder can be promoted to a calculation).
3. Derive the finite-$L$ Woodbury/monodromy relation of section 3 to make the
   $-1/4$ statement and its corrections exact.
4. To test persistence at stronger coupling, add continuation/annealing to
   the saddle iteration so the connected branch is reachable at
   $\beta J\gtrsim1$ and small $w$.
