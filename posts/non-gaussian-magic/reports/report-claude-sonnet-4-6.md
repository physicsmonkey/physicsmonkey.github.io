# Report on "Non-Gaussian Magic"

> Reviewer: Claude Sonnet 4.6 (claude-sonnet-4-6), May 22, 2026
> Scope: entry ngm.md and calculation-1

---

## Summary

The entry introduces a notion of non-Gaussian magic, motivates it via the SYK model's statistical $O(N)$ symmetry, and presents two supporting tests. The theoretical framing is sound, the numerical calculation is correctly implemented and internally consistent, and the proof in Test 2 is logically complete. I flag several clarifications and one substantive gap (the claim about simulation hardness for alternating circuits), but no mathematical errors.

I independently reran calculation-1 with different seeds (`--n-values 8 10 --trials 4 --sweeps 4 --seed 9999`) and obtained min ratios of $0.53$–$0.57$ and max ratios of $3.8$–$6.5$ for $N=8,10$, consistent with the reported results.

---

## Section: Introduction and Cases 1 and 2

**Algebraic setup.** The Majorana anticommutation relation $\{\psi_i,\psi_j\}=2\delta_{ij}$ and the Majorana string definition with phase $i^{W(W-1)/2}$ (ensuring Hermiticity) are standard and correct.

**Case 1.** The appeal to Wick's theorem for Gaussian states and the reference to determinant-based evaluation of single-string expectation values are correct. One small clarification: the statement "superpositions of low-weight strings" should specify that the number of terms in the superposition is polynomial for the efficiency claim to hold.

**Case 2.** The fermionic Gottesman-Knill theorem and the sampling-based efficiency for mixtures of stabilizer states are standard results, correctly stated.

---

## Section: Combining the Cases

**Order $U_G U_C$.** Applying $U_G^\dagger$ to a weight-$W$ string gives a linear combination of at most $\binom{N}{W}$ weight-$W$ strings. Each coefficient is a degree-$W$ polynomial in the entries of $O$, and there are exactly $\binom{N}{W}$ distinct weight-$W$ Majorana strings. The "at most" bound is correct (some coefficients can vanish). Each such string can then be evaluated exactly in the stabilizer state $U_C|0\rangle$ in polynomial time. This argument is correct.

**Order $U_C U_G$.** The entry argues that $|0\rangle^{\otimes N/2}$ is simultaneously Gaussian and stabilizer, so $U_G|0\rangle$ is Gaussian, and a single low-weight string after applying $U_C^\dagger$ (giving a single, possibly heavy Majorana string) can be evaluated in the resulting Gaussian state by determinant methods. This is correct.

**Notation.** The state $|0\rangle^{\otimes N/2}$ implicitly requires $N$ even. This should be stated explicitly.

**Hardness of alternating circuits.** The sentence "What we cannot generically do is combine more than one instance of $U_G$ and $U_C$" is stated as a fact without proof or reference. This appears to be a conjecture or a known folk theorem. If it is known, a reference would be helpful; if it is a conjecture, it should be labeled as such.

---

## Section: Why?

The motivation is clear and well-written. The SYK-4 Hamiltonian and its disorder distribution are standard. The robustness of magic is a well-defined measure. The observation that SYK has a statistical $O(N)$ symmetry that is a Gaussian unitary is correct: the SYK ensemble is invariant under $\psi_i \to \sum_j O_{ij} \psi_j$ for $O\in O(N)$, since this maps $J_{ijkl} \to (O^T)^{\otimes 4} J$, which has the same distribution by rotational invariance of the Gaussian ensemble.

---

## Section: Non-Gaussian Magic

**Definition.** $M_{\rm NG}(\rho)=\min_{U_G} M(U_G \rho U_G^\dagger)$ is natural and well-posed.

**Pure Gaussian states have $M_{\rm NG}=0$.** For any pure Gaussian state $|\psi\rangle$ (a Slater determinant), there exists a Bogoliubov transformation $U_G$ mapping it to the vacuum $|0\rangle^{\otimes N/2}$, which is a stabilizer state with $M=0$. So $M_{\rm NG}(|\psi\rangle)=0$. ✓

**Mixed Gaussian states.** The entry asserts that Gaussian states have zero non-Gaussian magic. For mixed Gaussian states (e.g., free-fermion thermal states), this is less immediate: $U_G \rho U_G^\dagger$ is again Gaussian for any $U_G$, but achieving $M=0$ requires mapping $\rho$ to a state on which the magic measure vanishes. This holds if the magic measure $M$ vanishes on all Gaussian (not just stabilizer) states. Whether robustness of magic or the stabilizer Renyi entropy vanishes for all mixed Gaussian states is not discussed. Since the main application is to SYK thermal states (which are non-Gaussian), this gap does not affect the core claims, but the definition statement should be qualified.

**Ansatz $|\psi\rangle=U_G U_C|0\rangle$.** Applying $U_G^\dagger$ gives $U_C|0\rangle$, a stabilizer state. So $M_{\rm NG}(|\psi\rangle)=0$ by definition. ✓

**Parameter counting.** $O(N)$ has dimension $N(N-1)/2$, while the independent components of $J_{ijkl}$ number $\binom{N}{4}\sim N^4/24$. So the orbit under $O(N)$ is far lower-dimensional than the full space of coupling tensors at large $N$. This is a plausible but informal argument; the numerical tests give it quantitative content.

---

## Section: Test 1

**Connection to SRE.** The leading-$N$ approximation $\langle \mu_{W=4}\rangle \propto J_{i_1 i_2 i_3 i_4}$ in the SYK thermal state rests on large-$N$ factorization (see the cited reference 2104.03336). Taking this on faith, the sum $\sum_{i<j<k<l} J_{ijkl}^4$ is indeed the $\alpha=2$ SRE quantity up to normalization, and its behavior under $O(N)$ rotations speaks to the stability of the SRE per fermion. The leap from "the quartic sum changes by at most a constant factor" to "the SRE per fermion is invariant in the thermodynamic limit" is correct, since $\text{SRE} \sim \log Q_4$ and a constant multiplicative change in $Q_4$ produces only a constant additive change in the SRE, leaving the density unchanged.

**Numerical implementation.** The code in `rotate_l4.py` is correct:

- The antisymmetric 4-tensor is stored with all permutations, and the Givens rotation is applied axis by axis. Sequential application over all 4 axes correctly computes $(O^{\otimes 4} T)$, since rotations on different indices commute.
- The `rotate_one_axis` function stores `old_p = tensor[p_tuple]` and `old_q = tensor[q_tuple]` as views into the original (un-modified) array before writing to the copy `out`, so there is no aliasing error.
- The equal-component lower bound $Q_2^2/\binom{N}{4}$ follows from Jensen's inequality ($\sum x_i^4 \geq (\sum x_i^2)^2/n$) with $Q_2$ fixed by $O(N)$ invariance of the Frobenius norm. ✓
- The reported ratios are verified to match the CSV output and summary JSON (checked manually for $N=8$, seed 9234).

**Reported numbers.** Spot-checks confirm the tables in the entry match `summary.json`:

| $N$ | min/init (entry) | min/init (JSON mean) | max/init (entry) | max/init (JSON mean) |
|-----|------------------|----------------------|------------------|----------------------|
| 8   | 0.50             | 0.5005               | 5.33             | 5.332                |
| 10  | 0.52             | 0.5155               | 4.76             | 4.760                |
| 12  | 0.57             | 0.5653               | 4.67             | 4.668                |
| 14  | 0.60             | 0.6017               | 4.34             | 4.340                |

**Independent rerun (different seeds).** Running with `--seed 9999 --trials 4 --sweeps 4` on $N=8,10$ gives min ratios of $0.53$–$0.59$ and max ratios of $3.8$–$6.5$, consistent with the reported range.

**Caveat on convergence.** With 8 sweeps over $N(N-1)/2$ pairs and a 25-point grid, the search is a heuristic. The entry correctly notes it is not a certificate of the global optimum. The slow upward trend in the min ratio with $N$ (0.50 to 0.60 over $N=8$ to $14$) could indicate a slow approach to 1 rather than convergence to a constant below 1; this cannot be distinguished from the current system sizes. The entry's characterization as "preliminary evidence" is appropriately cautious.

---

## Section: Test 2

The proof claims: for Haar-random $|\psi\rangle$, $\log \sup_{O\in O(N)} F_2(|\psi\rangle;O) = o(N)$ with high probability.

**Step 1: Fixed-$O$ distribution.** By Haar invariance, $U_O |\psi\rangle$ is also Haar-random, so $F_2(|\psi\rangle;O)$ has the same distribution for every fixed $O$. ✓

**Step 2: Expected value $O(1)$.** The identity string contributes 1. For a non-identity Majorana string $\mu(a)$ (which is traceless and satisfies $\mu(a)^2 = I$, so $\|\mu(a)\|_F^2 = D$), the standard Weingarten formula gives $\mathbb{E}_\psi[\langle\psi|\mu(a)|\psi\rangle^4] = O(D^{-2})$. With $2^N - 1 = D^2 - 1$ non-identity strings, the sum contributes $O(1)$. So $\mathbb{E}_\psi F_2 = O(1)$. ✓

**Step 3: Concentration.** $F_2$ is a degree-8 polynomial in the $2D$ real parameters of $|\psi\rangle$. Standard hypercontractive concentration for degree-$d$ polynomials on the real sphere gives tails decaying as $\exp(-c\, t^{2/d})$; for $d=8$ this is $\exp(-c\,t^{1/4})$. The entry's claim is consistent with this. ✓

**Step 4: Covering number.** $O(N)$ has dimension $N(N-1)/2$ and the bound $|\mathcal{N}_\epsilon| \leq (C/\epsilon)^{N(N-1)/2}$ is standard for compact Lie groups. ✓

**Step 5: Lipschitz constant.** The entry claims $|F_2(\psi;O) - F_2(\psi;O')| \leq CN2^N \|O-O'\|$. For a weight-$W$ string, the rotated expectation value changes by at most $O(W\|O-O'\|)$ in operator norm, so the $\ell^4$ term changes by $O(W\|O-O'\|)$. Summing over $2^N$ strings with $W\leq N$ gives the stated bound. This Lipschitz constant is exponential in $N$, but this is acceptable since the concentration tail is doubly exponential. ✓

**Step 6: Union bound.** Setting $t = e^{\delta N}$ and $\epsilon \sim e^{\delta N}/(N 2^N)$:

$$\log |\mathcal{N}_\epsilon| = O(N^2 \log(N 2^N / e^{\delta N})) = O(N^3)$$

which is polynomial in $N$. The concentration bound contributes $-c e^{\delta N/4}$, which is super-exponential. So the probability of $\sup_O F_2 > e^{\delta N}$ tends to zero. ✓

**Scope.** This proof applies to Haar-random pure states, not to SYK thermal states. The entry is explicit about this ("a simple version of the corresponding statement for random pure states"). The connection to SYK would require additional structure, e.g., showing that SYK states behave like pseudorandom states in the relevant sense. This is an interesting open direction left implicit.

---

## Section: Outlook

The summary correctly represents what was shown: two tests supporting high non-Gaussian magic in SYK-4. The connection to the holographic interpretation is motivated but is the part of the argument that rests most on physics intuition rather than rigorous proof.

---

## Overall Assessment

The entry is scientifically sound. The main theoretical claims are either correct or properly labeled as preliminary. No mathematical errors were found. The key limitations are:

1. The hardness of alternating $U_G$/$U_C$ circuits is asserted without proof or reference — this should be flagged as a claim.
2. The $M_{\rm NG}=0$ statement for mixed Gaussian states is stated more broadly than the pure-state argument supports.
3. The $N$-range of the numerics is modest; the "constant factor" characterization of the min/max variability is plausible but not established.
4. Test 2 is for random pure states; the SYK application is by analogy.

None of these gaps invalidate the entry's conclusions. Items 1 and 2 are the most worth addressing in a future revision.
