# Report on `syk-ed-op-stat.md`

> Reviewer: Codex (GPT-5.6 Sol)<br>
> Date: August 3, 2026<br>
> Scope: section-by-section correctness of the entry, calculations 1--5, and
> the previous Fable 5 report, following `journal.md`.<br>
> Verdict: **Substantially sound, with the cross-replica AP identification
> still quantitatively unresolved.**

## Summary verdict

The entry's central qualitative result is well supported: finite-$N$ probe
traces distinguish the conditional periodic and antiperiodic propagators of
the two-replica saddle, with inside probes mapping to P and outside probes to
AP. The boundary-condition assignment follows from an exact operator
monodromy identity, not merely from the numerics. The same-side comparisons
are particularly persuasive because their normalization follows from a
moment derivative and the ED data agree directly with the unrescaled
$G^{\mathrm P/AP}_{11}$ surfaces at the $5\times10^{-4}$--$3\times10^{-3}$
level in all six tests reported at $N=20$ and 24.

The main qualification concerns the replica-off-diagonal AP component. The
collision-normalized formula reproduces its sign-changing shape, but it is not
an equality at the tested sizes: its residual is $5.07\%$ at $N=20$,
$\beta J=0.5$, and $5.93\%$ and $11.74\%$ at $N=24$, $\beta J=0.5$ and 1.
These differences are smooth and far larger than the disorder jackknife
errors. Calling them finite-$N$ corrections is plausible, but is not yet a
diagnosis: there is no fixed-$(w,\beta J)$ size sequence, and the same-side AP
component is much more accurate at the same $N$. The precise normalized
$G^{\mathrm{AP}}_{12}$ prescription should therefore remain explicitly
conjectural.

I found no consequential sign, phase, variance, trace-contraction, or saved
metric error in the code I checked. I independently recomputed the headline
residuals and collision values from the archived `.npz` arrays. The
calculation-1 self-test passed; all 11 Python files parsed; the compact
Hamiltonian validator returned zero error; and the newer cross-replica trace
evaluator agreed with the original evaluator to $9.48\times10^{-16}$ at small
$N$. I did not rerun the expensive production disorder ensembles.

## Relation to the previous report

I agree with the July 30 Fable 5 report. Its two most important additions are
now incorporated into the entry: the monodromy derivation of inside/P versus
outside/AP, and the periodic-zero-mode explanation of the near-$-1/4$
collision product. Its request to include calculation 4 has also been
addressed. The present report mainly assesses those additions and the new
$N=24$ calculation-5 results.

## Section-by-section assessment

### Setup and proposed dictionary

The Majorana conventions are consistent. With $\gamma_i=\sqrt2\psi_i$, the
code uses

$$
\mu_A=i^{|A|(|A|-1)/2}\prod_{i\in A}\gamma_i,
$$

which is Hermitian and squares to one. For $q=4$, the Hamiltonian coefficient
in the code is $J_I/4$ multiplying $\mu_I$, exactly matching
$i^2J_I\psi_{i_1}\cdots\psi_{i_4}$, and the sampled variance is
$\mathbb E[J_I^2]=3!J^2/N^3$.

The parity discussion is correct. An odd clump $A$ is needed for a single
probe trace to be parity even, while thermal one-point functions of odd
strings vanish. The stronger $q=4$ selection rule that only weights divisible
by four have nonzero one-point functions is also correct for these
conventions, although the text could name the relevant antiunitary/Clifford
symmetry rather than leaving it as an example. Strictly, the allowed list also
contains the identity at $W=0$.

The finite-$N$ constructions around an even target string are algebraically
sound. Removing one fermion gives an odd $A_-$ and $N-W+1$ outside probes;
adding one gives an odd $A_+$ and $W+1$ inside probes. The collision algebra

$$
\psi_i\mu_A=c_i\mu_{A\triangle\{i\}},\qquad
c_i^2=-\tfrac12\ (i\notin A),\quad c_i^2=+\tfrac12\ (i\in A)
$$

is correct for odd $|A|$.

One conceptual transition should remain prominent: the initial observable
uses the sample-normalized thermal trace $\zeta=Y/Z$, whereas the successful
normalization-fixed tests use a ratio of disorder averages of unnormalized
traces. These are different ensemble observables at finite $N$. The entry now
states this clearly; it should not later summarize the result as verification
of the original $\mathbb E[\zeta\zeta]$ proposal without that qualification.

### Testing the proposal

All quoted calculation-2 and calculation-3 numbers match the saved summaries
and arrays. In particular:

- The $N=20$ outside/AP residuals are $0.050704$ at $w=3/20$ and $0.055165$
  at $w=4/20$, with exact ED collision value $-1/2$.
- The symmetry-averaged inside/P residuals are $8.249\times10^{-4}$ at
  $w=5/20$ and $1.622\times10^{-3}$ at $w=4/20$, with exact collision value
  $+1/2$.
- The quoted cosine similarities, maximum differences, and comparison with
  the jackknife-error norm also match.

The entry correctly says that dividing by $G_{12}(0,0)$ is an additional
prescription rather than a derived normalization. This is the key limitation
of the $G_{12}$ test. The excellent uncentered P residual is also dominated by
the normalization-protected periodic constant. The later warning that its
centered shape residual is about $4\%$, comparable to the AP shape residual,
is essential and prevents the earlier phrase "considerably stronger" from
being overread.

### Boundary conditions and zero-mode mechanism

The monodromy derivation is correct. Cyclicity gives

$$
Y_{Ai}(\beta)=\operatorname{Tr}(e^{-\beta H}\mu_A\psi_i),
$$

and commuting the probe through an odd $\mu_A$ yields a minus sign outside
$A$ and a plus sign inside it. The explanation of why one should not begin
with an independent "thermal AP sign" is also right: the usual fermionic KMS
sign comes from crossing the other fermionic insertion, which is on the other
replica in this observable.

The periodic-zero-mode explanation is physically and algebraically
reasonable. The saved scan reproduces the stated continuum-extrapolated
collision values, including

$$
G^{\mathrm P}_{12}(0,0)=4.752266,\qquad
G^{\mathrm{AP}}_{12}(0,0)=-0.052608,
$$

and their product $-0.250005$ at $\beta J=0.5,w=0.2$. Across the connected
points in the scan, the leading zero-mode estimate for P is accurate to
roughly $0.2$--$1.6\%$, the first-order AP estimate to
$4\times10^{-6}$--$1.2\times10^{-3}$ relative, and the product is closer to
$-1/4$ than either approximation separately. The entry appropriately calls
the product a leading-order relation and notes its growing violation with
coupling.

The Woodbury paragraph is a plausible structural explanation, not yet a
derivation. In particular, the coefficient approaching $1/4$ and the claim
that its correction is controlled by the contact shift have not been worked
out explicitly in the supporting record. Phrases such as "likely structural
reason" and "appears to be" are therefore the right confidence level. The
statement that $\Lambda$ has acquired a "concrete identity" should likewise
be read as the leading zero-mode identification
$\Lambda\simeq2/(\iint\Sigma_{12})$, not an exact identity.

### Same-side correlations

This is the strongest part of the entry. Differentiating the annealed moment
$M_A[K]$ gives a ratio with denominator $\mathbb E[X_A^2]$, so no empirical
amplitude, offset, time shift, or division by a collision value is introduced.
The $-1/2$ one-sided contact agrees with the Majorana algebra and the saddle
convention.

I independently recomputed the archived comparison metrics. At $N=20$ they
are exactly the quoted $5.287\times10^{-4}$ P and $1.556\times10^{-3}$ AP
off-diagonal residuals, with cosines $0.999999945$ and $0.999999905$. The
coarse- and fine-grid values also confirm that the P sector has a large
first-order lattice artifact removed by the aligned Richardson estimate. The
reported repeat using $L_\tau=180,360$ is not saved in an output or generated
by the current script, so that particular $0.0534\%$ number has a provenance
gap even though it is consistent with the displayed convergence pattern.

For completeness, the source-derivative equation would benefit from an
explicit one-line definition of the bilocal source term, including its sign
and any factor of $1/2$. The phrase "with the same source convention as in the
saddle equations" is enough for an internal note but makes the absolute
normalization derivation harder to audit independently. The contact identity
and numerical comparison strongly support the convention actually used.

### Pushing the numerics

The $N=24$ same-side table reproduces exactly from the saved arrays. The four
off-diagonal residuals are $1.523\times10^{-3}$, $3.372\times10^{-4}$,
$2.801\times10^{-3}$, and $5.212\times10^{-4}$ in the order printed. The
contacts, cosines, denominator effective sample sizes, and convergence flags
also match. These are useful robustness checks, but the entry correctly warns
that denominator ESS values near 7 out of 16 make them unsuitable for a
precision finite-size study.

The new $N=24$ cross-replica metrics also reproduce exactly. The trace
evaluator's independent small-$N$ validation passes at
$9.48\times10^{-16}$. The results should be interpreted asymmetrically:

- The normalized P surfaces have sub-percent full-surface residuals, but—as
  at $N=20$—a centered diagnostic is needed before concluding that their small
  time-dependent part has comparable precision.
- The AP residuals of $5.93\%$ and $11.74\%$ are decisive discrepancies from
  the proposed equality, despite good cosines. The difference/error norms of
  about 36 are descriptive rather than chi-square significances because the
  time points are highly correlated.

Likewise, saying the P cases are "within the disorder uncertainty" is a little
stronger than the norm diagnostic establishes. "Not resolved by the quoted
jackknife-error norm" would be statistically safer without a covariance-aware
test.

The comparison from $N=20,W=4$ to $N=24,W=8$ changes $N$, $w$, the neighboring
odd background, and the ensemble of collision strings together. It shows that
the AP discrepancy does not disappear under that combined change, but it is
not a finite-$N$ scaling test. A clean diagnosis needs at least two or three
sizes at fixed $w$ and $\beta J$, or an analytic correction to the external-leg
normalization.

### Outlook

The overall claim of a successful P/AP dictionary is fair if it refers to the
boundary conditions, qualitative $G_{12}$ shapes, and especially the
normalization-fixed $G_{11}$ tests. The statement that the worsening with
$\beta J$ is what one expects from neglected $1/N$ corrections is reasonable
intuition but not demonstrated by the available data. The same behavior could
also reflect the still-undemonstrated cross-replica normalization
prescription. The sentence about not expecting much more from a careful
finite-$N$ extrapolation is unclear and seems premature because no such
extrapolation has been performed.

There is also a typo: "satisyfying" should be "satisfying." Since the entry
was extended on August 3, its header date could be changed from July 30 to a
date range or supplemented by a last-updated date.

## Reproducibility and provenance

The calculation folders follow `journal.md`: each has a useful `readme.md`,
the production arrays and summaries are present, and the main entry links to
the relevant figures. I preserved the entry and calculations unchanged.

The following checks passed during this review:

- calculation-1 Clifford/trace/saddle self-test;
- syntax parsing of all 11 Python files;
- compact versus original Hamiltonian construction, maximum error zero;
- cross-replica energy-basis versus original trace evaluator, maximum error
  $9.48\times10^{-16}$;
- direct recomputation of all entry tables from the saved arrays.

Three minor provenance issues remain:

1. The claimed calculation-4 dense-trace error $1.7\times10^{-15}$ and the
   calculation-5 streaming comparison error $4.2\times10^{-16}$ are described
   in readmes, but the corresponding validation routines or outputs are not
   retained in the current files.
2. The $L_\tau=180,360$ calculation-4 check is likewise not archived.
3. `reports/scan_quarter_product.py` hard-codes the absolute path to the
   calculation-1 module, which makes the report calculation nonportable.

None of these gaps changes the numerical conclusions, but retaining the short
validators and the 360-grid scalar output would strengthen the research
record.

## Recommended revisions and next tests

1. State the main conclusion in layers: the P/AP monodromy assignment is
   exact; the raw same-side identification is numerically precise; the
   collision-normalized cross-replica identification is approximate and has
   an unresolved AP discrepancy.
2. Replace definitive references to an AP "finite-$N$ correction" with
   "finite-$N$ discrepancy" until a fixed-parameter size sequence or analytic
   correction identifies its origin.
3. Add centered-shape metrics for the $N=24$ P-$G_{12}$ tests and, if possible,
   a covariance-aware comparison rather than only the norm of pointwise
   jackknife errors.
4. Complete the finite-lattice Woodbury calculation. It is the shortest route
   to separating an exact P/AP kernel relation from the weak-coupling
   $-1/4$ approximation.
5. For a genuine finite-size test, compare several $N$ at fixed $w$,
   $\beta J$, time grid, and collision-string convention. This is the evidence
   needed to support the Outlook's attribution to $1/N$ effects.

With those qualifications, the entry is a strong and unusually transparent
piece of numerical evidence for the operator meaning of conditional
multireplica $G\Sigma$ variables.
