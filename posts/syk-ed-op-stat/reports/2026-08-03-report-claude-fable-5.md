# Report on `syk-ed-op-stat.md` (August 3 revision)

> Reviewer: Claude (Fable 5, Anthropic)<br>
> Date: August 3, 2026<br>
> Scope: section-by-section correctness of the entry as revised on August 3,
> with emphasis on the material added after the August 3 Codex review
> (calculation 6, the cross-replica extension of calculation 5, the
> centered-shape analysis, and the layered Outlook). Per `journal.md`, the
> entry and calculations were not modified. Rerunning the deterministic
> validators (`validate_same_side.py`, `validate_implementations.py`,
> `derive_boundary_woodbury.py`, the calculation-1 self-test) regenerated
> their output files with values identical to the archived ones.<br>
> Verdict: **Sound. Every quantitative claim in the entry reproduces from the
> archived outputs; the new exact boundary-update identity is independently
> verified; and an independent from-scratch ED implementation confirms the
> same-side dictionary end-to-end at N=10.**

## Summary verdict

The entry's three-layer conclusion is stated accurately and each layer is
supported at the appropriate level of rigor:

1. The inside/P and outside/AP boundary-condition assignment follows from the
   monodromy identity, which I verified both by hand and directly on the
   operator matrices of an independent implementation.
2. The same-side $G_{11}$ identification is normalization-fixed and precise.
   Beyond confirming the archived $N=20$ and $N=24$ results, I reimplemented
   the entire ED side from scratch (own Jordan–Wigner construction, own
   Hamiltonian sampler, own trace evaluator) and reproduced the correspondence
   at $N=10$: contacts $-1/2$ to machine precision and off-diagonal residuals
   $8.2\times10^{-4}$ (P) and $8.0\times10^{-4}$ (AP) against the conditional
   saddles at $w=0.4$, with 40 disorder samples. Nothing of the agreement
   depends on the original code base.
3. The cross-replica $G_{12}$ identification is correctly framed as a
   prescription. The collision normalization's P/AP reciprocity is now an
   exact continuum identity via calculation 6, which I verified independently
   (see below); what remains underived is the physical step of comparing the
   canonically normalized kernels to the ED moments, and the entry says so.
   The AP shape discrepancy ($3.3\%$–$6.4\%$ centered, growing with
   $\beta J$) is real, coherent, and honestly reported as undiagnosed.

I found no errors of substance. The issues listed at the end are one broken
LaTeX fragment, one presentational ambiguity, and provenance housekeeping.

## Checks performed

**Archived-output verification.** Every number in the entry's tables and prose
was checked against the corresponding `summary.json`,
`shape_diagnostics.json`, `woodbury_check.json`, `validation.json`,
`path_extrapolation_180_360.json`, and `quarter_product_scan.json`. All match,
including the $N=24$ cross-replica residuals and cosines, the centered-shape
table with its jackknife errors, the denominator effective sample sizes
(7.25/7.32), the same-side error-norm ratios ("about four times" is 4.07,
"about 1.6" is 1.61), the P-sector lattice convergence sequence
($2.37\%\to1.18\%\to0.0529\%$, with the 180/360 repeat $0.0534\%$ now
archived), and the scan-based accuracy claims for the zero-mode formulas.

**Reruns.** The calculation-1 self-test, the calculation-4 dense-trace
validator (max error $7.2\times10^{-16}$), the calculation-5 validators
(compact Hamiltonian exact; streaming $3.0\times10^{-16}$; cross-replica
evaluator $9.5\times10^{-16}$), and the full calculation-6 script all rerun
successfully and reproduce their archived records.

**Independent derivation checks.**

- *Monodromy.* Verified by hand from cyclicity and the Clifford algebra, and
  numerically on independently constructed operators: $\mu_A\psi_i=-\psi_i\mu_A$
  for $i\notin A$ and $+\psi_i\mu_A$ for $i\in A$ with odd $|A|\in\{3,5\}$.
- *Collision coefficients.* $c_i^2=-1/2$ outside and $+1/2$ inside recomputed
  by hand for the phase convention in the entry.
- *Calculation 6.* The corner claim $A_{\mathrm{AP}}=A_{\mathrm P}+2EF^T$ was
  checked directly against `discrete_derivative` (exact); the Woodbury update
  for $G$, the reduction $C_{\mathrm{AP}}=(1-2K_{\mathrm P})^{-1}C_{\mathrm P}$,
  and the explicit $m_{\mathrm{AP}}$ formula were verified with random
  replica-circulant $\Sigma$ (errors at $10^{-15}$) and, for the $2\times2$
  algebra, with exact rational arithmetic. The continuum substitution
  $k\to-1/2$, $a\to+1/2$, $b\to m$ gives $m_{\mathrm{AP}}=-1/(4m)$ exactly.
  The archived Richardson-extrapolated boundary values approach these limits
  at the $10^{-5}$–$10^{-6}$ level, supporting the substitution. This
  completes, and upgrades to "derived," what my July 30 report could only
  argue structurally.
- *Zero-mode formulas.* Recomputed $1/\iint\Sigma_{12}$ directly from the
  saved calculation-4 saddle: $4.7256$ versus $G^{\mathrm P}_{12}(0,0)=4.7523$
  ($0.56\%$), consistent with the entry's quoted $0.2$–$1.6\%$ range.

**Independent end-to-end ED test.** The from-scratch $N=10$ same-side check
described above. This closes the main residual worry about any single shared
implementation: the ED convention chain (Majorana normalization, string
phase, time ordering, contact convention, moment normalization) reproduces
the saddle prediction with none of the original code on the ED side.

**Script audits.** `compare_large_cross_replica.py`,
`analyze_cross_replica_shapes.py`, `cross_replica_moment_normalization.py`,
and `derive_boundary_woodbury.py` were read in full. The collision
coefficients are computed from the actual Clifford action with a
state-independence assertion rather than hard-coded signs; the collision
replay guards against ensemble mismatch; the jackknife is the standard
delete-one ratio estimate; the centered-shape analysis correctly jackknifes
the global scalar metrics rather than pointwise surfaces. No defects found.

## Section-by-section notes

**Setup and proposal.** Unchanged in substance since the August 3 Codex
review, whose assessment I share. The transition from the sample-normalized
$\zeta$ covariance to ratios of disorder-averaged unnormalized traces remains
clearly flagged.

**Testing the proposal.** All numbers verified. The closing paragraph now
correctly distinguishes the derived reciprocity of the two kernels from the
still-conjectural use of the collision-normalized kernels as comparators.
This is the right epistemic split.

**Mechanisms.** The monodromy subsection is correct, including the
double-counting warning, and the flavor-factorized restatement is a clean way
to see that spectator flavors cannot affect a boundary condition. The
zero-mode subsection is correct with one attribution nuance (issue 2 below).
The new exact-relation passage is correct and independently verified; its
final claim that two-grid residuals from $-1/4$ are pure time-discretization
error is now justified by the finite-lattice identity holding to
$4\times10^{-14}$ at every grid.

**Same-side correlations.** Verified in full, including the explicit source
convention added in response to the earlier audit request. The $N=10$
independent test gives additional confidence that the absolute normalization
is a consequence of the algebra, not of a shared convention hiding in one
code base.

**Pushing the numerics.** All tables verified against the archived summaries,
including the literal-weight variants mentioned in prose. The centered-shape
table matches `shape_diagnostics.json` exactly, and its interpretation — the
P time-dependent signal unresolved at 16 samples, the AP discrepancy resolved
and growing with coupling — is what the data support. The stated caveats
(descriptive error norms, rank-deficient covariance, entangled changes of
$N$, $w$, and background) are appropriate.

**Outlook.** The layered summary now matches the evidence layer by layer.
The attribution of the $\beta J$ worsening to neglected $1/N$ corrections is
labeled as expectation rather than result, which resolves my earlier framing
concern and the Codex report's recommendation 2.

## On the previous reports

Both prior reports hold up. Every recommendation of the August 3 Codex
report has been implemented except the genuine fixed-$(w,\beta J)$ size
sequence (its recommendation 5), which remains the most important open
numerical task. The three provenance gaps it listed are closed:
`validate_same_side.py`/`validation.json` and the 180/360 record now exist
and rerun, and `scan_quarter_product.py` now uses a relative path. My July 30
report's Woodbury speculation is superseded by calculation 6; its suggestion
to promote the scan out of `reports/` (also Codex's) is still open — the scan
remains the one calculation-grade artifact living in the reports folder, cited
by the entry.

## Issues found (none affect conclusions)

1. **Broken LaTeX in the entry.** In the continuum-limit display of the exact
   relation, `k\to-\frac12,qquad a\to+\frac12,qquad b\to m` is missing the
   backslashes on both `\qquad`s, so "qquad" renders as literal text.
2. **Attribution of the first-order AP accuracy.** The quoted relative errors
   of $4\times10^{-6}$–$1.2\times10^{-3}$ belong to the full convolution
   $[G_d\star\Sigma_{12}\star G_d](0,0)$ computed with the interacting
   diagonal propagator (this is what the scan evaluates). The further
   simplification $-\beta^2\bar\sigma_{12}/4$ displayed in the same equation
   chain is only $\sim0.5\%$ accurate at $\beta J=0.5$ (I recomputed
   $-\iint\Sigma_{12}/4=-0.05290$ versus $-0.052608$). The text reads as if
   the quoted precision applies to the final formula; a clarifying clause
   would fix it.
3. **Rounding nit.** The dense-trace check is $7.2\times10^{-16}$
   (7.2227), quoted as $7.3\times10^{-16}$ in the entry and calculation-4
   readme.
4. **Provenance: uncommitted work.** Calculation 6, the entire cross-replica
   extension of calculation 5 (script, outputs, diagnostics, figures), the
   new validators, the August 3 Codex report, and the revised entry are all
   untracked or modified-uncommitted in git. Given how much of the entry's
   evidence now lives in these files, committing them is the single most
   useful housekeeping step.
5. **Statistical wording (minor).** "Its difference is smaller than the norm
   of the pointwise jackknife-error surface" for the P-$G_{12}$ rows is a
   descriptive norm comparison, as the entry elsewhere acknowledges; since
   the centered analysis shows the P shape is simply unresolved, the
   uncentered statement carries little inferential weight and could be
   phrased as such at first occurrence rather than only in the centered
   paragraph that follows.

## Suggested next steps

1. Commit the untracked calculation folders and outputs (issue 4).
2. Fix the `qquad` typo and the first-order-accuracy attribution (issues
   1–2).
3. The open scientific question is unchanged: a fixed-$(w,\beta J)$ size
   sequence for the cross-replica AP component — e.g. $N=16,20,24$ at
   $w=1/4$, reusing the streaming evaluator — to test whether the coherent
   AP discrepancy scales like $1/N$. My $N=10$ same-side run (residuals
   $\sim8\times10^{-4}$, comparable to $N=20$ and $24$) incidentally shows
   the same-side sector is uniformly precise down to very small $N$, which
   sharpens the puzzle of why only the cross-replica AP object deviates:
   the discrepancy is plausibly in the external-leg/normalization step
   rather than in finite-$N$ corrections to the kernels themselves. That
   observation slightly favors deriving the cross-replica comparator from a
   source derivative (as was done for the same-side case) over brute-force
   size scaling, if a derivation can be found.
4. Promote `reports/scan_quarter_product.py` and its JSON to a numbered
   calculation folder, since the entry cites its results directly.
