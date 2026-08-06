# Calculation 7: scan of the P/AP collision product and zero-mode formulas

This calculation scans the connected two-replica saddle over
$\beta J\in\{0.25,0.5,1,2\}$ and $w\in\{0.1,0.2,0.3\}$ and tests three
statements at each point, with collision values Richardson-extrapolated from
aligned $L_\tau=90,180$ lattices:

1. the collision product
   $G^{\mathrm P}_{12}(0,0)\,G^{\mathrm{AP}}_{12}(0,0)$ against $-1/4$;
2. the periodic zero-mode prediction
   $G^{\mathrm P}_{12}(0,0)\simeq1/\iint\Sigma_{12}\,d\tau\,d\tau'$;
3. the first-order antiperiodic prediction
   $G^{\mathrm{AP}}_{12}(0,0)\simeq
   [G^{\mathrm{AP}}_{d}\star\Sigma_{12}\star G^{\mathrm{AP}}_{d}](0,0)$,
   where $G^{\mathrm{AP}}_{d}=-(D_{\mathrm{AP}}-\Delta\tau^2\Sigma_{11})^{-T}$
   is the diagonal-only AP propagator.

`scan_quarter_product.py` reuses the calculation-1 saddle solver. Run from
this directory with

```bash
python3 scan_quarter_product.py
```

The numerical record is `outputs/quarter_product_scan.json`.

## Result

On the reachable connected branch the product deviates from $-1/4$ by
$-2\times10^{-6}$ at $\beta J=0.25$, $-5\times10^{-6}$ at $\beta J=0.5$, and
$-4.3\times10^{-5}$ at $\beta J=1$, essentially independently of $w$. The
zero-mode formula reproduces $G^{\mathrm P}_{12}(0,0)$ to $0.2$–$1.6\%$, and
the first-order AP convolution matches to $4\times10^{-6}$–$1.2\times10^{-3}$
relative (for example $-0.052614$ predicted versus $-0.052608$ at
$\beta J=0.5$, $w=0.2$). The further simplification
$-\beta^2\bar\sigma_{12}/4$ used in the entry's heuristic derivation is only
accurate at the half-percent level.

At $w=0.1$ for $\beta J\ge0.5$, at $w\le0.2$ for $\beta J=1$, and at all of
$\beta J=2$, the plain damped iteration collapses to the replica-diagonal
branch $G_{12}=0$; those points are recorded with vanishing collision values
and were not part of the tested set. Reaching them requires continuation in
$w$ or $\beta J$.

The two-grid residuals from $-1/4$ were later shown by
[calculation 6](../calculation-6/readme.md) to be pure time-discretization
error: the product is an exact continuum identity on the connected branch.

Status: implemented and run by Claude (Fable 5, Anthropic) on July 30, 2026
in support of its [status report](../reports/2026-07-30-status-report-claude-fable-5.md),
where the scan was first recorded. Promoted from `reports/` to this
calculation folder by Claude (Fable 5) on August 4, 2026, at the author's
request; the script's output path was adjusted to `outputs/` and its
calculation-1 reference was already relative. The scan itself was not rerun
during the move; the archived JSON is the July 30 record.
