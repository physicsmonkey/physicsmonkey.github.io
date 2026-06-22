# Report on `cvsd.md` — "A Toy Model of Wormholes"

> Reviewer: Claude Code (Opus 4.8)<br>
> Date: June 22, 2026<br>
> Scope: correctness of the entry and `calculation-1`, per the guidelines in `journal.md`<br>
> Previous reports: none found in `reports/`.

## Summary verdict

The entry is internally consistent and, to the extent I could check, correct.
I re-derived every analytic step of the perturbative and large-$K$ sections by
hand and they all reproduce the formulas in the entry exactly. I re-ran
`calculation-1` in an isolated copy (no entry files modified) and reproduced the
tabulated numbers to all printed digits, and I independently re-derived the
$K=1$ crossing from a from-scratch operator construction (explicit $a,a^\dagger$
matrices + `scipy.linalg.expm`, full Hilbert space rather than the even-sector
analytic matrix elements) — it agrees to 8 digits.

The one substantive issue is **not** an error but a question of how much the
finite-$K$ "exact first crossing" data can bear: at $K\gtrsim 13$ the tabulated
exact crossings are **not converged in the cutoff** $n_{\max}$, with up to ~21%
spread, and at $K=16$ the value actually printed in the table (the
$n_{\max}=220$ entry) is the cutoff outlier. The entry already warns that this
column is contaminated by near-zeros of $A_1$, so the conclusions are not
overstated, but the table presents these unconverged numbers to 8 significant
figures and the Outlook's tentative "drifting downward" reading rests partly on
them. Details and a suggested more robust diagnostic are below.

---

## Section-by-section assessment

### Holographic motivation

This is framing/motivation, not a calculation, so there is nothing to verify
numerically. The physical analogy is stated carefully and with appropriate
hedging ("toy model", confidence Medium). The mapping it sets up — disconnected
$\leftrightarrow$ matter self-annihilating within each side, connected
$\leftrightarrow$ worldlines threading the wormhole — is exactly what the
oscillator model later encodes via the non-conservation of $a^\dagger a$ under
the quartic term. No issues. The $O(1/G)$ matter / $O(G)$ coupling scaling that
motivates the later $g=\gamma/K$, "amount of matter" $\sim K$ choice is
consistent.

### Oscillator model

Definitions of $A_1$, $A_2$, $|0\rangle$, and the disconnected/connected
identifications are clean. Two checkable claims, both correct:

- **Cauchy–Schwarz $A_1(s)\le A_2(s)$.** With
  $|\psi\rangle=e^{-sH}(a^\dagger)^{4K}|0\rangle$ one has
  $A_2=\langle\psi|\psi\rangle$ and $A_1=|\langle 0|\psi\rangle|^2$, so
  $A_1\le\langle 0|0\rangle\,\langle\psi|\psi\rangle=A_2$ since
  $\langle 0|0\rangle=1$. ✓
- **Endpoint values.** At $s=0$, $A_2=\langle 0|a^{4K}(a^\dagger)^{4K}|0\rangle
  =(4K)!$ and $A_1=|\langle 0|(a^\dagger)^{4K}|0\rangle|^2=0$ (and likewise
  $A_1=0$ at $\gamma=0$). ✓

The $\Delta I = K\delta$ weighting and the resulting competition
$A_1$ vs. $A_2 e^{-\Delta I}$ are well posed: since $A_2\ge A_1$, a positive
$\delta$ is exactly what is needed to make the comparison nontrivial.

### Perturbative calculation

I reproduced each formula in imaginary-time (interaction-picture) perturbation
theory.

- **$s\to\infty$ limits** of $A_1$ and $A_2$ via
  $e^{-sH}\to e^{-sE_{\rm gs}}|{\rm gs}\rangle\langle{\rm gs}|$, and the ratio
  $A_1/A_2\to|\langle 0|{\rm gs}\rangle|^2\to 1$ at small $g$. ✓
- **$A_1$.** The leading process needs $K$ insertions of the $a^4$ piece of
  $(a+a^\dagger)^4$; with $I(s)=\int_0^s e^{-4\omega\tau}d\tau$ and
  $\langle 0|a^{4K}(a^\dagger)^{4K}|0\rangle=(4K)!$ one gets
  $A_1=\big[(4K)!/K!\big]^2 [gI]^{2K}$. ✓
- **$A_2$ per-sector amplitudes $c_j$.** At order $g^j$ the only path to the
  $n-4j$ sector is $j$ factors of $a^4$ (any $j$ vertices summing to $\Delta
  n=-4j$ must all be the $-4$ piece), and the $a^4$ coefficient in
  $(a+a^\dagger)^4$ is $1$, giving the stated $c_j$. ✓ The text is commendably
  explicit that $A_2\simeq\sum_j|c_j|^2$ keeps the *leading process per sector*
  rather than a fixed order in $g$.
- **Ratio and rate function.** Setting $\ell=K-j$ gives the stated
  $A_2/A_1=\sum_\ell \frac{1}{(4\ell)!}\big[\frac{K!}{(K-\ell)!}\big]^2 y^{2\ell}$
  with $y=e^{-4\omega s}/(gI)$. Writing $y=Kc$, $\ell=\rho K$ and applying
  Stirling, the $\log K$ terms cancel and I recover
  $\Phi(\rho;c)=-2(1-\rho)\log(1-\rho)+2\rho+2\rho\log c-4\rho\log(4\rho)$
  **exactly**. ✓
- **Saddle.** $\Phi'(\rho)=0\Rightarrow c(1-\rho)=16\rho^2$, and back-substitution
  collapses to $f(c)=2\rho_*-2\log(1-\rho_*)$. ✓
- **Transition.** $f(c_*)=\delta$ with $q=1-\rho_*$ gives
  $\delta/2=1-q-\log q$, $c_*=16(1-q)^2/q$, and inverting
  $c(s)=c_*$ yields $s_*=\frac{1}{4\omega}\log(1+\frac{4\omega}{\gamma c_*})$ and
  the dual $\gamma_*$. ✓ The small-$\delta$ asymptotic $c_*\simeq\delta^2$ also
  checks: $q=1-\epsilon\Rightarrow\delta/2\approx2\epsilon$, so
  $c_*\approx16\epsilon^2\approx\delta^2$. ✓

This section is solid. One presentation nit: the claim that $A_2\le A_1$-type
bounds, the per-sector truncation, and the omission of the number-conserving and
$\Delta n=\pm2$ pieces of $(a+a^\dagger)^4$ are all controlled only at leading
order in $g=\gamma/K$; the entry says as much, but a one-line statement of the
*expected* size of the neglected pieces (they shift level energies at $O(g)=
O(1/K)$, i.e. an $O(1)$ shift of the $K$-extensive exponent) would make the
later "finite-$K$ has more structure" observation feel less like a surprise.

### Numerical calculation

**Reproduction.** Running `finite_k_numerics.py` in an isolated `/tmp` copy
reproduced the large-$K$ values $q_*=0.7662486082$, $c_*=1.140929199$,
$\gamma_*^{(\infty)}=0.06541110658$ and the entire $K=1\ldots16$ table to all
printed digits (e.g. $K=16$: exact $0.006113557$, pert $0.05618363$). ✓ The
even-sector Hamiltonian matrix elements ($6m^2+6m+3$ diagonal,
$(4m+6)\sqrt{(m+1)(m+2)}$ at $\Delta n=2$, $\sqrt{(m+1)\cdots(m+4)}$ at
$\Delta n=4$) are correct, and my independent full-space brute force confirms
the $K=1$ crossing ($0.107654492$). ✓

**Issue — cutoff convergence at large $K$.** The tabulated "exact" column uses
$n_{\max}=220$ and is reported to 8 sig figs, but it is not converged in
$n_{\max}$ for $K\gtrsim 13$. Spread across $n_{\max}\in\{80,120,160,220\}$:

| $K$ | $n_{80}$ | $n_{120}$ | $n_{160}$ | $n_{220}$ | spread |
|----|---------|----------|----------|----------|-------|
| ≤12 | — | — | — | — | <0.1% |
| 13 | 0.0049633 | 0.0049502 | 0.0048988 | 0.0049624 | 1.3% |
| 14 | 0.0052401 | 0.0047921 | 0.0053898 | 0.0054159 | 12.0% |
| 15 | 0.0064700 | 0.0064728 | 0.0064769 | 0.0063391 | 2.1% |
| 16 | 0.0076205 | 0.0076156 | 0.0076103 | **0.0061136** | 20.8% |

The non-monotonic, non-converging behavior (and, at $K=16$, the fact that the
*tabulated* $n_{\max}=220$ number is the lone outlier while
$n_{\max}=80/160/160$ cluster near $0.00762$) is a direct symptom of the
mechanism the entry names: the "first crossing" is being set by a near-zero of
$A_1=|B|^2$, where $B(\gamma)=\langle 0|e^{-sH}(a^\dagger)^{4K}|0\rangle$ changes
sign. Near such a zero $\log(A_2 e^{-K\delta}/A_1)\to+\infty$ and the crossing
location depends on the precise sub-leading spectrum, which is exactly what
$n_{\max}$ perturbs. So these are not well-defined crossover couplings at the
quoted precision.

This does not invalidate the entry's stated conclusion — it explicitly says the
exact first-crossing column "should not be read as a clean approach to the
large-$K$ perturbative saddle" and is "a useful warning." But two refinements
would make the section more honest and more informative:

1. **Report uncertainty, not 8 digits.** For $K\ge13$ the values are good to
   ~1–2 digits at best; quoting them to 8 figures (and drawing them as sharp
   points in `finite_k_crossings.svg`) overstates precision. Error bars from the
   $n_{\max}$ spread, or rounding to the converged digits, would be more
   faithful.

2. **Use a crossing diagnostic that is insensitive to $A_1$ zeros.** Because the
   physical question is whether the *connected* weight $A_2 e^{-K\delta}$ exceeds
   the *disconnected* weight $A_1$ in the regime where both are smooth, the
   "first" crossing as $\gamma$ increases from $0$ is the fragile choice: at
   small $\gamma$, $A_1\sim[\gamma I/K]^{2K}$ is tiny and noisy and its zeros
   dominate. Candidates that would sharpen the Outlook: track the *last* crossing
   (large-$\gamma$ side, where $A_1$ is generically nonzero); compare
   $\log A_2$ vs $\log A_1 + K\delta$ on a grid and locate the crossing of the
   *smoothed/enveloped* curves; or compare against $|B|$ rather than $|B|^2$ and
   handle the sign explicitly. The `k10_options_vs_gamma.svg` panel already shows
   the two weights are smooth and well-separated for $K=10$ — extending exactly
   that smoothed comparison (rather than the first-crossing finder) to
   $K=13\ldots16$ would tell you directly whether a genuine crossover survives.

**Minor.** The Outlook's suggestion that the crossover "potentially even
drift[s] downward as $K$ is increased" is, on the present data, partly an
artifact of (1): the apparent drop at $K\ge13$ coincides exactly with the onset
of the $A_1$-near-zero regime and the loss of cutoff convergence, not obviously
with a physical trend. I would treat the monotone, well-converged $K\le12$ band
(crossings settling around $\gamma_*\sim0.03$, a factor ~2 below
$\gamma_*^{(\infty)}=0.065$) as the trustworthy statement, and label the
$K\ge13$ points as cutoff/zero-limited.

### Outlook

The "mixed" summary is fair and appropriately cautious. The open question it
poses — does a nonzero-$\gamma$ crossover persist as $K\to\infty$? — is the right
one, but as noted, answering it convincingly needs the more robust diagnostic
above, because the current first-crossing data degrade precisely in the
large-$K$ regime that matters. Nothing here is wrong; it is the natural place
where the suggested numerical refinement would pay off.

---

## What I ran

- Re-derived by hand: Cauchy–Schwarz, $s\to\infty$ limits, $A_1$, $c_j$, $A_2$,
  the $A_2/A_1$ sum, $\Phi(\rho;c)$, the saddle $c(1-\rho)=16\rho^2$, $f(c)$, the
  transition relations, and $c_*\simeq\delta^2$. All match the entry.
- Ran `calculation-1/finite_k_numerics.py` in `/tmp/cvsd_check` (entry files
  untouched): reproduced the full table, the large-$K$ constants, and produced
  the convergence spread table above. Runtime ~5 s with the venv at
  `~/.venvs/sci`.
- Independent from-scratch operator construction (explicit $a,a^\dagger$,
  `expm`, full — not even-sector — space): $K=1$ first crossing
  $0.107654492$ vs. entry $0.10765449$. ✓

## Things I did not / could not fully check

- I did not push $n_{\max}$ high enough to *converge* the $K\ge13$ exact crossings
  (this is the point: they may not converge to a single value in the
  first-crossing definition). A definitive resolution needs the alternative
  diagnostic, which I did not implement so as not to alter the calculation.
- I did not verify the large-$K$ saddle is the global (vs. local) maximum of
  $\Phi$ over $\rho\in[0,1]$ analytically for all $c$; at the transition
  ($c_*\approx1.14$, $\rho_*\approx0.234$) it is a clean interior maximum, which
  is the regime that matters.

---

## Addendum (same day): robust crossover and the source of the "drift"

After discussion with Brian, I built `calculation-2` to follow up on the two
refinements suggested above. It resolves what was going on at $K\ge13$ and, I
think, changes the headline of the Outlook. (This addendum is part of the
report; per the journal guidelines I did not edit the entry — these are findings
to fold in.)

### Two distinct problems with the first crossing, both now pinned down

1. **$A_1$ has genuine zeros.** $B(\gamma)=\langle0|e^{-sH}(a^\dagger)^{4K}|0\rangle$
   is genuinely not sign-definite. Verified at **60-digit precision** (mpmath):
   $B$ changes sign at a discrete set of finite $\gamma$ (competing multi-step
   paths cancel), including at small $\gamma$, well beyond the leading monomial
   $B\sim(-1)^K\gamma^K$. So $\log A_1\to-\infty$ on a real, discrete set and the
   "first crossing" is set by whichever zero comes first — not the physics.

2. **Double precision cannot reach the crossover region for $K\ge12$.** This is
   the part I had not appreciated in the main report. In the relevant region
   $|B|$ falls below $\sim10^{-21}$, and a *double-precision* eigendecomposition
   then returns deterministic round-off **noise** for $B$ — which is why it
   *reproduces across $n_{\max}$* (identical banded round-off) yet is wrong in
   magnitude and sometimes in sign. Side-by-side at $K=13$:

   | $\gamma$ | double precision $B$ | 60-digit $B$ | |
   |---------|---------------------|--------------|---|
   | $10^{-4}$ | $-1.3\times10^{-46}$ | $-4.2\times10^{-51}$ | wrong magnitude |
   | $2\times10^{-3}$ | $-3.0\times10^{-27}$ | $+3.4\times10^{-34}$ | wrong sign |
   | $10^{-2}$ | $+7.3\times10^{-23}$ | $-1.3\times10^{-25}$ | wrong sign |
   | $3.4\times10^{-2}$ | $+4.54\times10^{-21}$ | $+4.54\times10^{-21}$ | exact |
   | $5\times10^{-2}$ | $-4.4557\times10^{-18}$ | $-4.4557\times10^{-18}$ | exact |

   So raising $n_{\max}$ never helped because the obstruction is the
   floating-point floor, not the Hilbert-space truncation. This means the
   $K\ge13$ rows of the entry's exact column are **not just imprecise but, in the
   small-$\gamma$ region, not computed at all** by the double-precision code.

### The fix and the result

`calculation-2/robust_crossover.py` computes $A_1$ and $A_2$ at $s=1$ in
**extended precision** (mpmath `eigsy`, 40 digits) and defines the crossover from
the **smooth envelope** of $A_1$: the local maxima of $|B|$ trace the envelope
$E$ (zeros are where the oscillation vanishes), so interpolating $\log A_1$
through its peaks and crossing it with the smooth $\log(A_2e^{-K\delta})$ gives a
single, zero-insensitive $\gamma_*(K)$. The construction is illustrated in
`calculation-2/figures/diagnostic_K15.svg`.

Result, $K=10\ldots18$ (full data in `calculation-2/outputs/crossover_summary.csv`):

| $K$ | envelope (ext. prec., $s{=}1$) | $\gamma$-window | $s$-average | first crossing (calc-1) |
|----|------|------|------|------|
| 10 | 0.0399 | 0.030 | 0.079 | 0.0338 |
| 11 | 0.0412 | 0.026 | 0.066 | 0.0306 |
| 12 | 0.0392 | 0.038 | 0.062 | 0.0305 |
| 13 | 0.0399 | 0.033 | 0.078 | **0.0050** |
| 14 | 0.0393 | 0.029 | 0.068 | **0.0048** |
| 15 | 0.0385 | 0.026 | 0.091 | **0.0064** |
| 16 | 0.0393 | 0.035 | 0.078 | **0.0076** |
| 17 | 0.0386 | 0.031 | 0.072 | **0.0071** |
| 18 | 0.0376 | 0.028 | 0.0084 | **0.0082** |

The envelope crossover is **cutoff-stable** — unchanged at $K=14$ over
$n_{\max}=80,120,160$, and at $K=16,18$ up to $n_{\max}=140$
(`convergence_check.csv`) — and it sits on a **flat, nonzero plateau
$\gamma_*\approx0.038$–$0.041$ across $K=10\ldots18$, with no downward drift.**
The apparent collapse of the first-crossing column at $K\ge13$ (bold) is exactly
the artifact of problems 1–2: the finder latches onto the small-$\gamma$ zeros
and the noise floor below them. See `figures/robust_crossover_vs_K.svg`.

*Precision note (added after the Codex (GPT-5) report
`report-2026-06-22-codex-gpt5.md`).* The many-digit agreement across $n_{\max}$
is reproducibility, not accuracy: it measures cutoff stability at a fixed
$\gamma$-grid. The dominant uncertainty is the grid resolution of the envelope
construction — refining the $90$-point grid moves the $K=14$ crossover over
$0.0389$–$0.0394$ (about $1\%$). So the individual $\gamma_*(K)$ are reliable
only to about two significant figures ($\approx0.039$), and earlier digit-heavy
quotes (including the $0.039340186$ above) overstate the precision. The entry's
refined-study table and wording were rounded/softened accordingly; the
qualitative conclusion (nonzero plateau, no drift) is unaffected.

### Bearing on the entry's Outlook

- The Outlook's tentative reading that the crossover "potentially even
  drift[s] downward as $K$ is increased" looks like it should be **retracted**:
  with a zero-insensitive, precision-correct diagnostic the crossover does not
  drift down; it plateaus at $\gamma_*\approx0.039$.
- That plateau is **below** the large-$K$ perturbative value
  $\gamma_*^{(\infty)}=0.0654$ — a genuine finite-$K$ suppression of about $40\%$
  — but it is clearly heading to a nonzero limit, not to zero. So the entry's
  central physical claim is *strengthened*: there is a finite-coupling window in
  which the connected (wormhole) weight dominates, and it survives to the largest
  $K$ I can resolve.
- The absolute crossover value is method-dependent at the $\sim30\%$ level
  (envelope $\approx0.039$, $\gamma$-window $\approx0.030$, calc-1 first crossing
  at small $K\approx0.033$), reflecting genuine ambiguity in "where the
  oscillating finite-$K$ amplitude crosses." The robust, method-independent
  statement is the *existence of a nonzero plateau with no downward drift*, not a
  precise number. The envelope is the most principled of the three.

### Caveats on the cross-checks

- The **$s$-average** confirms the no-drift trend but reports a higher band
  ($\sim0.06$–$0.09$) because $A_1(s),A_2(s)$ grow as $s$ decreases, so the
  integral over $s\in[0.7,1.3]$ is weighted toward the small-$s$ end rather than
  the nominal $s=1$. Its $K=18$ entry ($0.0084$) is a clear outlier: that
  cross-check runs in double precision, which — consistent with problem 2 —
  itself fails at $K=18$. So it should be read only as qualitative support.
- Everything here is at the leading model definition ($\omega=s=\delta=1$). I did
  not re-examine whether the plateau value moves with $\delta$ or $s$, which
  would be the natural next question if one wants the $K\to\infty$ limit of
  $\gamma_*$ quantitatively.

### What I ran (addendum)

- 60-digit mpmath cross-check of $B(\gamma)$ confirming the zeros are real and
  that double precision fails below $|B|\sim10^{-21}$.
- `calculation-2/robust_crossover.py`: extended-precision envelope crossover for
  $K=10\ldots18$ at $n_{\max}=90$, cutoff check at $K=14$
  ($n_{\max}=80/120/160$), plus the $\gamma$-window, $s$-average, and first-cross
  cross-checks. Runtime a few minutes in `~/.venvs/sci` (mpmath added there).
