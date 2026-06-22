# Report on `cvsd.md` -- "A Toy Model of Wormholes"

> Reviewer: Codex (GPT-5)<br>
> Date: June 22, 2026<br>
> Scope: correctness of the entry, `calculation-1`, `calculation-2`, and the existing Opus 4.8 report, following `journal.md`.<br>
> Prior report considered: `reports/report-2026-06-22-claude-opus-4-8.md`.

## Summary

The analytic part of the entry looks correct. I rechecked the definitions of
`A_1`, `A_2`, the Cauchy-Schwarz inequality, the leading perturbative
calculation, and the large-`K` saddle formulas. I found no algebraic error in
the displayed perturbative result or in the large-`K` transition condition.

The numerical story has improved substantially relative to the original
first-crossing calculation. The refined `calculation-2` correctly identifies two
problems with the original "exact first crossing": genuine zeros of the
transition amplitude and a double-precision floor in the small-amplitude region.
The envelope diagnostic is a reasonable way to ask the smoother question the
toy model is really meant to address.

The main caution is that the refined envelope crossover is still a heuristic
diagnostic, not a derived observable, and its numerical precision is overstated
in a few places. The reported plateau near `gamma_* ~ 0.04` is plausible and
supported by the saved outputs, but the digits quoted for `gamma_*` should not
be read literally. The grid spacing, local-maximum envelope construction, and
limited cutoff checks are not consistent with 8- or 9-digit accuracy.

## Main Findings

### 1. Analytic calculation is internally consistent

The perturbative calculation in the entry checks out:

- With `H_0 = omega a^\dagger a`, there is no zero-point factor in the leading
  expressions. This matches the formulas in the entry.
- The leading contribution to
  `B(s)=<0|exp(-sH)(a^\dagger)^(4K)|0>` requires exactly `K` insertions of the
  `a^4` term, giving
  \[
  A_1(s)=\left[\frac{(4K)!}{K!}\right]^2[gI(s)]^{2K}
  \]
  at leading order.
- The per-sector expression for `A_2` is correctly described as "leading
  process in each Fock-number sector," not as the fixed-order leading expansion
  of `A_2`.
- The ratio
  \[
  \frac{A_2}{A_1}\simeq
  \sum_{\ell=0}^K \frac{1}{(4\ell)!}
  \left[\frac{K!}{(K-\ell)!}\right]^2 y^{2\ell}
  \]
  follows from the change of variables `ell = K-j`.
- The large-`K` rate function and saddle equation are correct:
  \[
  \Phi(\rho;c)=
  -2(1-\rho)\log(1-\rho)+2\rho+2\rho\log c
  -4\rho\log(4\rho),
  \]
  \[
  c(1-\rho)=16\rho^2,
  \qquad
  f(c)=2\rho_*(c)-2\log(1-\rho_*(c)).
  \]
- For `omega=s=delta=1`, I independently reproduced
  \[
  q_*=0.766248608162,\quad c_*=1.140929199308,\quad
  \gamma_*^{(\infty)}=0.065411106579.
  \]

This part of the entry is solid.

### 2. `calculation-1` is now best read as a diagnostic failure mode

The `calculation-1` code appears to construct the even-sector Hamiltonian
correctly. The matrix elements for `(a+a^\dagger)^4` are right:

- diagonal: `6m^2 + 6m + 3`,
- `Delta n = 2`: `(4m+6) sqrt((m+1)(m+2))`,
- `Delta n = 4`: `sqrt((m+1)(m+2)(m+3)(m+4))`.

The first-crossing data are useful because they expose the zeros of `A_1`, but
they are not a robust measure of the connected/disconnected transition. The
entry now says this clearly. I agree with Opus's diagnosis that the small-`K`
first crossings should not be extrapolated.

One minor reporting issue remains: the table still gives
`gamma_*^{exact, first}` to many digits even though this column is now known to
be a fragile diagnostic. Since the refined section supersedes it, this is not a
scientific error, but rounding that column more aggressively or labeling it
"fragile first crossing" would better match the interpretation.

### 3. The refined `calculation-2` is a meaningful improvement, but its precision is overstated

The logic of `calculation-2` is sound at the level of physical diagnosis:

- `B(gamma)=<0|exp(-sH)|4K>` need not be sign definite.
- Zeros of `B` produce artificial divergences in
  `log(A_2 exp(-K delta)/A_1)`.
- Double precision is not reliable when `B` falls to the `1e-21` scale or
  below.
- A zero-insensitive envelope of `A_1=|B|^2` is a reasonable diagnostic for the
  smooth competition of the two weights.

I verified from the saved `diagnostic_K15.csv` that the reported envelope
crossing is obtained by the stated procedure: the local maxima of `log A_1` are
interpolated, and the envelope crosses `log(A_2 e^{-K delta})` at
`gamma = 0.03847643076885681` for `K=15`.

However, the numerical precision in the refined section should be read with
caution:

- The extended-precision scan uses `GAMMA_PTS = 90` on a geometric grid from
  `0.004` to `0.13`. Near `gamma ~ 0.04`, adjacent grid points are separated by
  about four percent. Linear interpolation can print many digits, but the
  envelope location is not determined to those digits.
- The envelope is built from grid-detected local maxima. Unless those maxima are
  refined, the envelope itself inherits grid-resolution error.
- The cutoff-convergence claim is demonstrated only at `K=14`. This is a useful
  spot check, but it does not establish cutoff convergence for all
  `K=10,...,18`, especially since `nmax=90` is relatively close to the initial
  state `|4K>` at `K=18`.
- The entry says the resulting crossover is cutoff-converged, then cites the
  `K=14` check. I would phrase this as "a cutoff-convergence spot check at
  `K=14` is stable" rather than as a global convergence statement.

The robust conclusion should be qualitative/semiquantitative:

> With a zero-insensitive envelope diagnostic, the crossover stays nonzero and
> lies near `gamma_* ~ 0.04` for `K=10,...,18`.

I would not treat the table entries such as `0.039340186` as high-precision
numbers.

### 4. The `s`-average cross-check is weaker than the entry suggests

The refined section says the `s`-average gives the same qualitative picture. In
the saved `crossover_summary.csv`, this is mostly true for `K=10,...,17`, but
the `K=18` `s_average_gamma_star` is

```text
0.008423238284818177
```

which is in the same collapsed band as the unreliable first-crossing diagnostic,
not in the stable nonzero plateau band. Opus notes this caveat in its addendum,
but the entry itself does not. The entry would be more accurate if it described
the `s`-average as a partial/qualitative cross-check with a visible large-`K`
failure, rather than as fully corroborating the plateau.

### 5. The Opus report is high quality and mostly still current

I agree with Opus's original checks of the analytic calculation and with the
addendum's diagnosis of the first-crossing problem. The addendum's main
conclusion -- that the apparent downward drift was an artifact and that the
refined diagnostic supports a nonzero plateau -- is consistent with the saved
outputs.

My main addition to the Opus report is a numerical-analysis caveat: the envelope
diagnostic is good enough to establish the existence of a nonzero plateau in
this finite-`K` range, but not good enough to justify the quoted precision of the
individual crossover values.

## Section-by-Section Notes

### Holographic Motivation

The physical motivation is qualitative, and I found no internal inconsistency.
The analogy between disconnected saddles and matter self-annihilation is clearly
marked as heuristic. The scaling discussion, with `O(1/G)` matter and `O(G)`
quartic coupling, is consistent with later setting `g = gamma/K`.

One provenance/publishing nit: the byline still says `Reports from ...` even
though reports now exist. Before publication, this should be replaced with links
to the report files.

### Oscillator Model

The definitions of `A_1` and `A_2` are clean. The Cauchy-Schwarz claim is
correct: for
\[
|\psi\rangle=e^{-sH}(a^\dagger)^{4K}|0\rangle,
\]
one has `A_2=<psi|psi>` and `A_1=|<0|psi>|^2 <= A_2`.

The interpretation of `A_1` as disconnected and `A_2` as connected is coherent
within the toy model.

### Perturbative Calculation

No algebraic issues found. The caveat to keep visible is that the calculation
does not resum the number-conserving and `Delta n = +/-2` pieces of
`(a+a^\dagger)^4`; later numerical structure is therefore not surprising.

### Numerical Calculation

The old first-crossing calculation is now correctly demoted to a warning sign.
The `K=10` options plot is a useful addition because it shows directly how
zeros of `A_1` create multiple crossings.

The refined numerical study is a good next step, but I recommend softening
precision claims and adding one sentence about the `K=18` `s`-average outlier.

### Outlook

The entry currently has only a header for Outlook. The natural takeaway to put
there, after the refined numerics, would be:

- the leading large-`K` perturbative calculation predicts a transition at
  nonzero `gamma`;
- finite-`K` exact first crossings are contaminated by zeros of `A_1`;
- a zero-insensitive envelope diagnostic supports a nonzero crossover plateau
  near `gamma ~ 0.04` for the accessible range;
- the precise asymptotic value remains method-dependent and needs a more
  analytic understanding of the oscillatory transition amplitude.

## What I Ran

- Imported both numerical scripts successfully in the current scientific venv.
- Independently recomputed the large-`K` constants and finite-`K` perturbative
  roots for `K=1,4,10,16`; all matched the saved values.
- Verified from `calculation-2/outputs/diagnostic_K15.csv` that the saved
  envelope curve gives the reported `K=15` envelope crossing.
- Started a full `calculation-2` rerun in `/tmp` and a separate
  extended-precision spot check. Both were interrupted after several minutes
  because they were no longer lightweight checks. I therefore did not complete a
  full independent rerun of `calculation-2`.

## Recommended Follow-Ups

1. Increase `GAMMA_PTS` or refine local maxima of `|B|` with a local optimizer
   before quoting envelope crossing values.
2. Add cutoff checks for at least `K=10,16,18`, not only `K=14`.
3. Report envelope crossovers to two or three significant figures unless a
   resolution/convergence study supports more.
4. Mention the `K=18` `s`-average outlier in the entry or omit the `s`-average
   from the main narrative.
5. Consider an analytic explanation of the zeros of `B(gamma)`, perhaps by
   separating the monotone `Delta n=-4` path from loops involving `Delta n=0`
   and `Delta n=+/-2` processes.
