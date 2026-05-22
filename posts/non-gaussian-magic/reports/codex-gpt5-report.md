# Report on `ngm.md`

> Codex with GPT-5, May 22, 2026
> Scope: section-by-section assessment of `ngm.md` and spot check of `calculation-1`.
> Status: critical review only; no edits were made to the entry or supporting calculation.

## Executive assessment

The entry is coherent and the main physical idea is plausible: quotienting ordinary stabilizer magic by Gaussian rotations is a natural way to separate genuinely non-Gaussian/non-Clifford structure from basis-dependent Gaussian structure. The numerical results quoted in Test 1 match the stored output files in `calculation-1`, and a small smoke rerun of the script succeeded.

I found two issues worth addressing before treating the note as final:

1. A notational typo in Test 1: the rotated tensor formula uses the unprimed tensor indices inside the primed sum.
2. A proof gap in Test 2: the stated stretched-exponential tail bound for the whole random variable `F_2` is plausible but not derived from the preceding facts as written. Supplying a reference or a short moment/hypercontractive estimate would strengthen the argument.

## Header and setup

The definition of Majorana strings is standard and the Hermiticity phase is correct for ordered products. The text implicitly assumes even `N` later when writing the Hilbert space dimension as `2^{N/2}` and the product state `|0>^{\otimes N/2}`. It would be helpful to say near the start that `N` is even when a Hilbert-space representation is needed.

Minor language notes:

- "thanks for the fermionic analog" should be "thanks to the fermionic analog".
- "I only know of a one paper" should be "I only know of one paper".

## Case 1

The Gaussian-state tractability claim is sound for low-weight Majorana strings via Wick/Pfaffian evaluation. The determinant-method sentence is plausible as a pointer to stronger single-string/high-weight evaluation, though I did not independently verify the cited paper.

One possible precision improvement: for an arbitrary superposition, efficiency requires either low weight with polynomially many terms or another compact representation. The current phrase "superpositions of low-weight strings" is fine if read as polynomial-size superpositions.

## Case 2

The stabilizer-state claim is consistent with the fermionic Gottesman-Knill analogy. The statement about mixtures is appropriately qualified: efficiency follows if the mixture has an efficient sampler or other efficient decomposition.

## Combining the cases

The argument for `U_G U_C |0>` is correct for constant-weight observables. A Gaussian unitary maps a weight-`W` Majorana monomial into a linear combination of weight-`W` monomials, with at most `binom(N,W)` independent components, and the stabilizer expectation values can then be evaluated term by term.

The switched-order argument is also basically correct: Clifford conjugation maps a single Majorana string to a single Majorana string, possibly high weight, and the remaining Gaussian expectation can be evaluated efficiently. The caveat about not generically combining multiple alternating layers is reasonable.

## Why?

The motivation is clear. The statistical `O(N)` symmetry of the SYK ensemble is the right symmetry to compare with Gaussian rotations of Majorana modes.

Minor language notes:

- "in marginally increases" should be "it marginally increases".
- "insenstive" should be "insensitive".

## Non-Gaussian magic

The definition

```tex
M_{NG}(\rho)=\min_{U_G} M(U_G\rho U_G^\dagger)
```

is natural. If `M` is a stabilizer magic monotone that vanishes on convex mixtures of fermionic stabilizer states, then Gaussian states should indeed have zero non-Gaussian magic: a Gaussian unitary brings a general Gaussian state to a product thermal/Fock-diagonal normal form, which is a mixture of occupation-basis stabilizer states. For pure Gaussian states this is immediate.

The parameter-counting argument is suggestive rather than decisive. It is fine as motivation, but a very structured objective can sometimes be strongly affected by a lower-dimensional group action. The later tests help, but the entry should keep the distinction between heuristic parameter counting and proof.

## Test 1

The numerical calculation is internally consistent. I checked that:

- `rotate_l4.py` samples independent antisymmetric 4-form components as standard Gaussians.
- The rotation is applied along each tensor index, preserving the antisymmetric 4-form structure.
- The reported means in `ngm.md` agree with `outputs/summary.json`.
- The equal-component lower bound is `Q_2^2 / binom(N,4)`, as stated.
- A smoke rerun with `N=8`, two trials, and two sweeps completed successfully and reproduced the expected qualitative behavior.

There is one formula typo in the entry. In the displayed equation around lines 66-68, the inner tensor should be `J_{i'j'k'l'}`, not `J_{ijkl}`:

```tex
\sum_{ijkl}\left|\sum_{i'j'k'l'}
O_{ii'}O_{jj'}O_{kk'}O_{ll'}J_{i'j'k'l'}\right|^{2\alpha}.
```

The interpretation of the calculation should remain cautious. The Givens search is a coordinate/local search and not a certificate of the global optimum. The entry says this explicitly. The claim that the increase is "no worse than constant factor" should be read as an observation for `N <= 14`, not as an asymptotic conclusion.

## Test 2

The concentration-plus-net strategy is the right shape for the claim, and the net-size/Lipschitz bookkeeping is plausible. The final implication, "for every fixed `delta > 0` with high probability", is the standard formulation of `log sup_O F_2 = o(N)` in probability.

The main proof gap is the fixed-`O` tail estimate. The text states that because `F_2` is a degree-8 polynomial, standard hypercontractive concentration gives

```tex
Pr[F_2>t] <= exp(-c t^{1/4}).
```

As written, this does not follow just from `E F_2 = O(1)` and degree. One needs either:

- a cited hypercontractive inequality on the complex sphere plus a bound such as `||F_2||_2 = O(1)`, or
- a direct high-moment estimate for `F_2`, strong enough to union-bound over the `O(N)` net.

This is likely fixable, but it should be made explicit because the entire asymptotic proof rests on getting a stretched-exponential tail at `t = exp(delta N)`.

The Lipschitz bound

```tex
|F_2(O)-F_2(O')| <= C N 2^N ||O-O'||
```

is crude but sufficient for the net argument if the tail bound is valid. A fully rigorous version should track whether the constant is uniform over all weights and over the chosen normalization of the Majorana-string basis, but no scaling obstruction is apparent.

## Outlook

The outlook accurately summarizes the intended lesson, but it somewhat overstates the status of Test 1. Test 1 gives finite-size numerical evidence for an order-one change in the weight-4 `alpha=2` sum; it does not by itself show thermodynamic invariance of the SRE density for SYK thermal states. The random-state Test 2 is closer to a proof, subject to the concentration gap above.

## Suggested revisions

1. Fix the primed-index typo in the Test 1 rotation formula.
2. Add the even-`N` convention near the setup.
3. Add a citation or derivation for the fixed-`O` tail bound in Test 2, including the needed norm estimate for `F_2`.
4. Soften the wording in the outlook from "we showed" to "the tests support" for the SYK thermal-state thermodynamic claim, unless an additional analytic argument is added.
