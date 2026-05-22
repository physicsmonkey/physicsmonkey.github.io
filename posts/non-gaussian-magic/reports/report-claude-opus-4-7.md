# Report on `ngm.md`

> Reviewer: Claude Opus 4.7 (`claude-opus-4-7`), 2026-05-22
> Scope: section-by-section review of the entry, audit of `calculation-1`, and
> commentary on the two prior reports (`codex-gpt5-report.md` and
> `report-claude-sonnet-4-6.md`).
> No edits were made to the entry or supporting calculations.

## Executive assessment

The entry is in good shape. The framing is clear, the definition is natural,
and both tests are now substantively supported. Since the earlier Codex
report was written the entry has been revised in two important ways:

1. The Test 1 displayed sum (line 67) now correctly carries primed indices
   on $J$, which is the typo that Codex flagged.
2. The fixed-$O$ tail estimate in Test 2 has been rewritten. The current
   version derives the tail from the Hilbert–Schmidt identity
   $\sum_a \langle \psi | \mu(a) |\psi\rangle^2 = D$ plus a single-observable
   Haar concentration bound, rather than asserting a hypercontractive tail.
   This is exactly the gap Codex flagged, and the rewrite closes it cleanly
   (in my view more cleanly than the hypercontractive route the Sonnet
   report had also accepted).

The remaining issues I found are minor and editorial:

- A handful of small language typos (carried over from earlier drafts) are
  still present.
- A few claims would benefit from the qualifications already suggested by
  the prior reports (the mixed-Gaussian case, the hardness of alternating
  $U_G/U_C$ circuits, and the leap from finite-$N$ Test 1 results to a
  thermodynamic statement).
- One small piece of bookkeeping in the union bound at the end of Test 2.

I independently re-ran `calculation-1` with a fresh seed (`--seed 4242`,
trials = 3, sweeps = 4) on $N=8, 12$ and obtained min ratios in
$0.36$–$0.63$ and max ratios in $3.0$–$4.7$, consistent with the reported
ranges in the entry and with the two prior reports' independent reruns.

---

## Header and setup

The Hermiticity phase $i^{W(W-1)/2}$ is the correct one for the ordered
product $\psi_1^{a_1}\cdots\psi_N^{a_N}$, and the anticommutation
$\{\psi_i,\psi_j\}=2\delta_{ij}$ is the standard convention. The text says
"$N$ is taken to be even", so the implicit even-$N$ assumption that both
prior reports asked to be made explicit *is* already stated at the end of
the setup paragraph — only the wording is unobtrusive; promoting it to a
short stand-alone sentence would help.

Carry-over language nits (still present):

- Line 21: "thanks for the fermionic analog" → "thanks to the fermionic
  analog".
- Line 27: "I only know of a one paper" → "I only know of one paper".
- Line 38: "in marginally increases" → "it marginally increases".
- Line 44: "insenstive" → "insensitive".

These were flagged by the Codex report and have not yet been picked up.

## Cases 1 and 2

Both cases are stated correctly and at the right level of detail. The Codex
and Sonnet reports' qualification on "superpositions of low-weight strings"
(efficiency needs polynomially many terms) is reasonable; the sentence as
written is fine if read in context, so this is a stylistic call.

## Combining the cases

The argument is correct. A Gaussian unitary $U_G$ conjugates a weight-$W$
Majorana string into a linear combination of at most $\binom{N}{W}$
weight-$W$ strings (the weight is preserved because the antisymmetric
$W$-th exterior power of $O$ acts on the string), and each resulting string
can be evaluated in the stabilizer state $U_C|0\rangle^{\otimes N/2}$ in
polynomial time. The reversed order $U_C U_G$ works because
$|0\rangle^{\otimes N/2}$ is both Gaussian and stabilizer and because a
Clifford conjugation maps a single Majorana string to a single (possibly
high-weight) Majorana string, which the determinant method can handle.

The Sonnet report rightly flagged that the closing sentence — "What we
cannot generically do is combine more than one instance of $U_G$ and
$U_C$" — is stated without a reference. As written it is most naturally
read as a (very plausible) folk claim about the hardness of alternating
fermionic Clifford / Gaussian circuits; I agree it would be worth either a
hedge ("we conjecture") or a citation. The recent paper referenced at line
27 (2505.06336) is the natural place to look.

## Why?

The motivation is well-told. The statistical $O(N)$ symmetry of the SYK
ensemble follows from the rotational invariance of the Gaussian
distribution on $J_{ijkl}$, since the action $\psi_i \to \sum_j O_{ij}
\psi_j$ corresponds to $J \to O^{\otimes 4} J$, which preserves the
distribution.

## Non-Gaussian magic

The definition $M_{NG}(\rho) = \min_{U_G} M(U_G \rho U_G^\dagger)$ is
natural and the pure-Gaussian and ansatz cases both reduce to zero by the
arguments given.

The point that the Sonnet report raised about *mixed* Gaussian states is
worth restating: vanishing of $M_{NG}$ on mixed Gaussian states requires
the underlying magic measure to vanish on all Gaussian states, not just
stabilizer states. For pure Gaussian states this is immediate (a
Bogoliubov transformation brings them to $|0\rangle^{\otimes N/2}$). For
the robustness of magic and the stabilizer Rényi entropy applied to
thermal free-fermion states, the vanishing statement is plausible but the
entry does not justify it. Since the application is SYK (a non-Gaussian
state), this does not affect the downstream argument; a one-line
qualification would tighten the statement, however.

The parameter-counting argument is correctly labelled as motivation. I
agree with Codex that it is suggestive rather than decisive — a structured
objective can in principle be moved a lot by a lower-dimensional group
action.

## Test 1

**Mathematical content.** The reduction from the SRE at $\alpha=2$ in the
thermal state to the quartic sum $Q_4(J) = \sum_{i<j<k<l} J_{ijkl}^4$
follows the large-$N$ relation $\langle \mu_{W=4}\rangle \propto J_{ijkl}$
to leading order. The $\alpha=1$ case is the squared Frobenius norm, which
is genuinely $O(N)$-invariant, and the $\alpha=2$ case is the $\ell^4$
norm on the orbit of the symmetric group / orthogonal group — non-trivial
under rotations. Constant-factor stability of $Q_4$ implies only an
additive constant shift in $\log Q_4$ and hence a vanishing change in the
SRE *density*, which is the relevant quantity in the thermodynamic limit.

**Code audit.** I read `rotate_l4.py` carefully:

- The antisymmetric 4-tensor is materialized over all permutations with
  the correct sign in `random_antisymmetric_4form`, and the independent
  components are sampled as standard Gaussians (the SYK normalization
  cancels in the ratios reported).
- `rotate_one_axis` correctly grabs `old_p`/`old_q` from the *unmodified*
  input tensor before writing to the copy `out`, so there is no aliasing
  bug. Sequential application along all four axes implements
  $J \to O^{\otimes 4} J$ correctly because rotations on disjoint indices
  commute.
- The equal-component lower bound $Q_2^2/\binom{N}{4}$ is the Cauchy–Schwarz
  / Jensen lower bound for the $\ell^4$ norm at fixed $\ell^2$ norm, which
  is itself $O(N)$-invariant.

**Numerical match.** The means in `outputs/summary.json` reproduce the
tables in the entry to the quoted precision:

| $N$ | reported min / init | summary mean | reported max / init | summary mean |
|-----|---------------------|--------------|---------------------|--------------|
| 8   | 0.50                | 0.5005       | 5.33                | 5.332        |
| 10  | 0.52                | 0.5155       | 4.76                | 4.760        |
| 12  | 0.57                | 0.5653       | 4.67                | 4.668        |
| 14  | 0.60                | 0.6017       | 4.34                | 4.340        |

The entry's "about $0.33$" for the lower-bound ratio at the larger sizes
matches the summary ($0.3346$ at $N=14$).

**Independent rerun.** With `--seed 4242 --trials 3 --sweeps 4` on
$N\in\{8,12\}$ I obtained per-trial min ratios in $0.36$–$0.63$ and max
ratios in $3.0$–$4.7$. These are consistent with the reported means given
the small trial counts in the rerun. (The Sonnet rerun reported a similar
qualitative agreement at $N=8,10$.)

**Interpretation.** I echo Codex's caution: the Givens search is a local
coordinate optimizer and is not a certificate of the global $O(N)$
extremum, so the conclusion is best read as "Gaussian rotations within
reach of the local search modify the $\alpha=2$ weight-4 sum by only an
order-one factor in this range of $N$." The slow drift of the min ratio
from $0.50$ at $N=8$ to $0.60$ at $N=14$ leaves open whether it converges
to a constant strictly below 1 or drifts to 1; with only four $N$ values
this is not distinguishable from the data. The outlook's claim of
"evidence that the $\alpha=2$ SRE per fermion is invariant in the
thermodynamic limit" is appropriately hedged ("we obtained evidence"), but
the qualification in the body of Test 1 — "preliminary evidence" — is
stronger and should be preferred.

## Test 2

The current proof is substantively different from the version reviewed in
the Codex report. Codex's main objection was that the stretched-exponential
tail $\Pr[F_2 > t] \leq \exp(-c t^{1/4})$ was asserted from the polynomial
degree alone, without a derivation. The current entry no longer uses that
route; it derives the fixed-$O$ tail directly. I verified the new argument
step by step.

**Haar invariance.** $U_O$ is a unitary on the fermion Hilbert space, so
$U_O|\psi\rangle$ is Haar-random whenever $|\psi\rangle$ is; this reduces
the analysis to $O=I$ for each fixed $O$. ✓

**Identity and mean.** The identity string contributes $1$, and for each
non-identity string the fourth moment of $\langle\psi|\mu(a)|\psi\rangle$
is $O(D^{-2})$ by the Weingarten formula, giving $\mathbb{E}F_2 = O(1)$. ✓

**The Hilbert–Schmidt identity.** The expansion of $|\psi\rangle\langle\psi|$
in the orthogonal basis $\{\mu(a)\}$ with $\mathrm{tr}(\mu(a)\mu(b)) = D
\delta_{ab}$ gives
$$
1 = \|\,|\psi\rangle\langle\psi|\,\|_F^2 = \frac{1}{D}\sum_a \langle\psi|\mu(a)|\psi\rangle^2,
$$
so $\sum_a X_a^2 = D$ with $X_a = \langle\psi|U_O^\dagger\mu(a)U_O|\psi\rangle$.
The entry's identity is correct. ✓

**From sum-of-fourth-powers to a max.** Using $F_2 = \sum_a X_a^4$ and
$\sum_a X_a^2 = D$ with $X_{\mathrm{id}} = 1$,
$$
F_2 - 1 = \sum_{a\neq 0} X_a^4 \leq \max_{a\neq 0} X_a^2 \cdot \sum_{a\neq 0} X_a^2 \leq \max_{a\neq 0} X_a^2 \cdot (D-1).
$$
So $F_2 > t \Rightarrow \max_{a\neq 0} X_a^2 > (t-1)/(D-1)$, matching the
entry. ✓

**Single-observable concentration.** For each non-identity string,
$A = U_O^\dagger \mu(a) U_O$ is a traceless Hermitian involution. The
function $\psi \mapsto \langle\psi|A|\psi\rangle$ is $2$-Lipschitz on the
unit complex sphere (since $\|A\|=1$), so Levy's lemma gives
$$
\Pr_\psi\!\left[|\langle\psi|A|\psi\rangle| > s\right] \leq 2 e^{-c D s^2}
$$
for a universal constant $c$. The entry's Beta-distribution remark is an
equivalent route. ✓

Combining with $s^2 \approx t/D$ and union-bounding over the $D^2-1$
non-identity strings gives $\Pr[F_2 > t] \leq C D^2 e^{-c t}$ for the
relevant range of $t$. ✓

**Deterministic bound for $t > D$.** As pointed out in the entry,
$|X_a|\leq 1$ and $\sum_a X_a^2 = D$ together imply $F_2 \leq D$ almost
surely, so the event $\{F_2 > t\}$ is empty when $t > D$. ✓

**Lipschitz bound and net argument.** The crude
$|F_2(O) - F_2(O')| \leq C N 2^N \|O - O'\|$ follows from
$\sum_a 1 = 2^N$, $W\leq N$, and the operator-norm change of a single
weight-$W$ string under $O \to O'$ being $O(W\|O-O'\|)$. The
$\epsilon$-net of $O(N)$ at scale $\epsilon$ has
$|\mathcal{N}_\epsilon| \leq (C/\epsilon)^{N(N-1)/2}$. Setting
$t = e^{\delta N}$, $\epsilon \sim t/(N 2^N)$, the log of the net size is
$O(N^3)$, while the per-point tail contributes $-c\, e^{\delta N}$ in the
exponent. The doubly-exponential gap is the right side of the comparison,
so $\sup_O F_2 \leq e^{\delta N}$ with high probability for every fixed
$\delta > 0$. ✓

**A small bookkeeping comment.** The entry's summary line of the union
bound reads
$$
\exp\!\Big[\tfrac{N(N-1)}{2}\log\!\big(\tfrac{CN2^N}{e^{\delta N}}\big) - c\, e^{\delta N}\Big]
$$
"up to harmless changes of constants in the exponent and the extra
fixed-$O$ prefactor $D^2 = e^{O(N)}$." The $D^2$ prefactor and the
overall logarithm of the net size are both polynomial in $N$ (or at most
linear in $N$ via $\log 2^N$), so they are indeed dominated by
$c\,e^{\delta N}$; the parenthetical disclaimer is fine. If the
argument were to be tightened for a paper, separating out the prefactor
explicitly would make the bookkeeping more transparent.

Subject to the standard caveats about extending from Haar-random states to
SYK eigenstates (the entry is explicit about this), the Test 2 argument is
sound.

## Outlook

The outlook is calibrated, but I agree with Codex that "we gave two
calculations that suggest that SYK-4 states have high non-Gaussian magic"
is slightly stronger than what was shown for Test 1 alone — Test 1 gives
finite-$N$ evidence consistent with constant-factor stability rather than
demonstrating thermodynamic-limit invariance. The "evidence that the
$\alpha=2$ SRE per fermion is invariant in the thermodynamic limit"
phrasing is the right one and could be repeated in the summary sentence.

## Comments on the prior reports

- The **Codex report** correctly identified the Test 1 primed-index typo
  (now fixed) and the proof gap in the original Test 2 (also addressed by
  the rewrite). Its suggested revisions (1)–(4) are well chosen; items
  (1) and (3) appear to have been acted on, while items (2) and (4) (the
  even-$N$ note and the wording in the outlook) are still pending. The
  language nits are also still present.

- The **Sonnet report** accepted the *original* Test 2 proof via the
  hypercontractive route (tail $\exp(-c t^{1/4})$, valid degree-$8$
  polynomial concentration on the complex sphere). That argument is
  correct as far as it goes; the current rewrite achieves the same end via
  a more elementary route and is self-contained. The Sonnet report's main
  substantive flag — that the alternating-circuit hardness claim is stated
  without proof — has not been addressed and remains the most natural
  thing to revise next.

- Both prior reports raised the mixed-Gaussian-state subtlety in the
  definition of $M_{NG}$; the entry has not yet incorporated that
  qualification.

## Suggested revisions (cumulative across the three reports)

1. Sweep up the language nits ("thanks to", "one paper", "it marginally",
   "insensitive").
2. Hedge or cite the claim that alternating $U_G/U_C$ circuits cannot
   generically be classically simulated.
3. Add one sentence clarifying that the $M_{NG} = 0$ statement for
   Gaussian states holds rigorously for the pure case and depends on a
   per-measure check for mixtures.
4. Optionally, tighten the outlook wording to match the more careful
   "preliminary evidence" phrasing used in the body of Test 1.

None of these are essential to the correctness of the entry; the
mathematical core is in good shape.
