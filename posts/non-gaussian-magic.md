---
title: "Non-Gaussian Magic"
date: 2026-05-22
tags: [SYK, magic, fermions]
description: "A notion of magic invariant under Gaussian rotations, with two tests suggesting SYK-4 has high non-Gaussian magic."
---

# Non-Gaussian Magic

> Brian Swingle + Claude Code with Opus 4.7 + Codex with GPT 5.5 
> May 22, 2026
> Confidence: High
> Reports from [Sonnet 4.6](non-gaussian-magic/reports/report-claude-sonnet-4-6.html), [Opus 4.7](non-gaussian-magic/reports/report-claude-opus-4-7.html), [GPT 5.5](non-gaussian-magic/reports/codex-gpt5-report.html)
> Acknowledgements: Helpful discussions with Val Bettaque, Juan Maldacena, and Anna Biggs

For systems of many fermions, there are two limits where the description becomes classically tractable. For the purpose of this discussion, I'll take classically tractable to mean that expectation values of low-weight fermion operators can be computed with time and memory that scale polynomially with the number of fermions.

Let's label our Majorana fermion modes $\psi_i$ with the index $i$ running from $1$ to $N$. Each of these mode operators is Hermitian and squares to the identity, and they collectively obey the algebraic constraint
$$ \{ \psi_i , \psi_j \} = 2 \delta_{ij}.$$
A complete basis of Hermitian operators is provided by the Majorana strings,
$$ \mu(a)=i^{W(W-1)/2} \psi_1^{a_1} \cdots \psi_N^{a_N},$$
where the binary $N$-component vector $a$ specifies the string and the phase factor ensures $\mu$ is Hermitian. $W = \sum_i a_i$ is the weight of the string. Throughout $N$ is taken to be even.

## Case 1
If the fermions are in a Gaussian state, then thanks to Wick's theorem we can evaluate the expectation of any low-weight Majorana string in terms of the $W=2$ object, $\langle \psi_i \psi_j \rangle$. The same is true for superpositions of low-weight strings.

There is also an efficient [determinant-based method](https://arxiv.org/abs/2412.05367) for evaluating such expectation values for a single string of any weight, or superpositions of such strings with polynomially many terms.

## Case 2
If the fermions are in a pure fermionic stabilizer state, then thanks to the fermionic analog of the Gottesman-Knill theorem, we can also efficiently compute any Majorana string expectation value.

For general mixtures of pure stabilizer states, the computation is efficient if one can efficiently sample pure stabilizer states from the mixture.

## Combining the cases

It is interesting to ask if we can combine these two tools in some way. I only know of one [paper](https://arxiv.org/abs/2505.06336) on this.

In my own thinking, I considered a simple composition of the two cases. Start with a fixed pure stabilizer state, apply a fermionic Clifford unitary $U_C$, and then apply a fermionic Gaussian unitary $U_G$,
$$ |\psi \rangle = U_G U_C |0\rangle^{\otimes N/2}.$$

Given a weight $W$ string $\mu$, the expectation value in $|\psi\rangle$ can be evaluated by using the Heisenberg picture to express $U_G^\dagger \mu U_G$ as a sum of at most $\binom{N}{W}$ (the number of weight-$W$ strings, since $U_G$ doesn't change the weight) strings. Each string expectation value can then be efficiently evaluated using Clifford/stabilizer technology, so the overall time and memory cost are poly($N$) provided $W$ is constant.

What if we switched the order? Then one can still simulate expectation values of low-weight strings provided $|0\rangle^{\otimes N/2}$ is both Gaussian and stabilizer. Briefly, the Clifford unitary turns a single low-weight string into an arbitrary weight string, but a single such string can be efficiently evaluated using the determinant methods mentioned above. What we cannot generically do is combine more than one instance of $U_G$ and $U_C$.

## Why?

This kind of construction might be useful in that it marginally increases our ability to classically simulate some many-fermion states. I was motivated to think about it because of a recent [paper](https://arxiv.org/abs/2602.12339) with my PhD student Valerie Bettaque in which we studied measures of non-stabilizerness in fermionic thermal states of the Sachdev-Ye-Kitaev model at large $N$ (see also this [line of work](https://arxiv.org/abs/2509.17417)). Non-stabilizerness is also known as quantum [magic](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.71.022316). 

The SYK model with $p=4$-body interactions is 
$$ H = \frac{1}{4!} \sum_{ijkl} J_{ijkl} \psi_i \psi_j \psi_k \psi_l $$
where $J_{ijkl}$ is a Gaussian random coupling with zero mean and variance $3! J^2/N^3$. As such, the SYK model really refers to an ensemble of Hamiltonians, with the example above being one instance of the ensemble. We wanted to study non-stabilizerness/magic in this model to both get some analytical handle on magic at large $N$ in a strongly interacting system and to better understand the manifestation of magic in holographic duality. We primarily considered a measure of magic known as [robustness of magic](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.090501), and I think we achieved both goals.

However, there is one aspect of our analysis that I found unsatisfying. The SYK model has a lovely statistical $O(N)$ symmetry that is essentially a Gaussian unitary transformation, so one expects any self-averaging physical property to be insensitive to such rotations. However, the measures of magic we considered are all built on a fixed stabilizer basis, essentially a fixed choice of Majorana strings (as opposed to Gaussian rotations of them). So I started thinking about how to formally adapt the magic measures to account for this.

## Non-Gaussian magic

Given any magic measure $M(\rho)$, we can define the associated non-Gaussian magic measure as 
$$ M_{NG}(\rho) = \min_{U_G} M(U_G \rho U_G^\dagger).$$
In essence, we consider all possible Gaussian rotations of the state and take the minimum magic among those.

This definition has the virtue that Gaussian states, which typically have lots of magic (as we showed), have zero non-Gaussian magic. More generally, the ansatz $|\psi\rangle = U_G U_C |0\rangle^{\otimes N/2}$ can have lots of magic but its non-Gaussian magic is zero. This explains our choice of ansatz above.

Gaussian unitaries implement $O(N)$ rotations of the fermion modes,
$$ U_G \psi_i U_G^\dagger = \sum_j O_{ij} \psi_j,$$
and they have $N(N-1)/2$ continuous parameters. Just from parameter counting, we therefore expect that such Gaussian unitaries cannot substantially reduce the magic in thermal states of the 4-body (or higher) SYK model. If true, this would establish that there is no strong dependence on the precise choice of fermion modes when defining magic. That would be very satisfying if we want to interpret magic in terms of coarse-grained properties, since we likewise don't expect these to depend strongly on precisely how we define the modes or the precise SYK instance.

## Test 1

One simple test of this hypothesis is provided by looking at SYK thermal states. Many measures of magic are related to sums of powers of string expectation values,
$$ \sum_a \langle \mu(a) \rangle^{2\alpha}.$$
For example, consider the case where $\mu$ is weight $4$. In the SYK thermal state at large $N$, the statistical $O(N)$ symmetry requires, to leading order in $1/N$,
$$ \langle \mu_{W=4} \rangle \propto J_{i_1 i_2 i_3 i_4},$$
where the indices $i_1\cdots i_4$ are those in $\mu$ (see also this [paper](https://arxiv.org/abs/2104.03336)).

So for this set of strings, the dependence on Gaussian unitaries reduces to asking how
$$ \sum_{ijkl} \left|\sum_{i'j'k'l'} O_{ii'} O_{jj'} O_{kk'} O_{ll'} J_{i'j'k'l'}\right|^{2\alpha}$$
depends on $O_{ii'}$.

For $\alpha=1$, the sum is an $O(N)$ invariant, so it does not depend on the Gaussian rotation. However, for $\alpha=2$, the Gaussian rotation can in principle modify the sum.

To get a first numerical sense of the size of this effect, I implemented a Jacobi-style search over Givens rotations in [calculation-1](non-gaussian-magic/calculation-1/readme.html). The calculation studies the independent antisymmetric components and tracks
$$ Q_4(J) = \sum_{i<j<k<l} J_{ijkl}^4 $$
under $O(N)$ rotations. The overall SYK normalization is irrelevant for the ratios, so the independent tensor entries were sampled as standard Gaussians. The search is not a proof of the global optimum, but it does preserve the $O(N)$ constraint exactly and gives a useful stress test.

With $8$ disorder realizations and $8$ sweeps of the Givens search, the optimized minimum divided by the initial value was approximately
$$
\begin{array}{c|cccc}
N & 8 & 10 & 12 & 14 \\
\hline
Q_4^{\rm min}/Q_4^{\rm initial} & 0.50 & 0.52 & 0.57 & 0.60
\end{array}
$$
while the optimized maximum divided by the initial value was approximately
$$
\begin{array}{c|cccc}
N & 8 & 10 & 12 & 14 \\
\hline
Q_4^{\rm max}/Q_4^{\rm initial} & 5.33 & 4.76 & 4.67 & 4.34 .
\end{array}
$$
For comparison, the equal-component lower bound $Q_2^2/\binom{N}{4}$ was about $0.33$ of the initial value at the larger sizes, so the search did not come close to completely flattening the tensor components. This preliminary evidence suggests that the $\alpha=2$ quantity is not invariant, but the ability of Gaussian rotations to reduce it is only an order-one effect in these system sizes, and may weaken slowly with $N$. The ability to increase $Q_4$ is larger than the ability to decrease it (although still no worse than constant factor), reflecting an inherent asymmetry of the $l^4$ norm under the $O(N)$ orbit.

In particular, if the ability of $O(N)$ rotations to modify these sums-of-strings is limited to constant factors, this would only change measures like the [stabilizer Renyi entropy](https://arxiv.org/abs/2106.12587) (SRE) by a constant shift. So the Renyi entropy per fermion would be unchanged in the thermodynamic limit.

Agent note: Codex (GPT-5) implemented `calculation-1` and edited this section of the entry to report the numerical results from that calculation on May 22, 2026.

## Test 2 

Here is a simple version of the corresponding statement for random pure states. Let $D=2^{N/2}$ be the Hilbert space dimension, and define
$$
F_2(|\psi\rangle;O)=\sum_a \left|\langle \psi|U_O^\dagger \mu(a) U_O|\psi\rangle\right|^4,
$$
where $O\in O(N)$ is the Gaussian rotation implemented by $U_O$. Up to an additive normalization independent of $O$, the stabilizer Renyi entropy at $\alpha=2$ is $-\log F_2$. Hence the non-Gaussian version minimizes the entropy by maximizing $F_2$ over $O(N)$.

The claim is that for a Haar-random pure state,
$$
\log \sup_{O\in O(N)} F_2(|\psi\rangle;O)=o(N)
$$
with high probability. Since $F_2(|\psi\rangle;I)\geq 1$ because of the identity string, this implies
$$
0 \leq S_2(|\psi\rangle)-S_{2,NG}(|\psi\rangle)
= \log \frac{\sup_{O\in O(N)}F_2(|\psi\rangle;O)}{F_2(|\psi\rangle;I)}
= o(N)
$$
with high probability. Thus the Gaussian optimization cannot change the stabilizer Renyi entropy density of a random pure state.

The proof is a concentration-plus-net argument. First fix a Gaussian rotation $O$. Since Haar measure is unitary invariant, $F_2(|\psi\rangle;O)$ has the same distribution for every fixed $O$. For the identity string the contribution is $1$. For each non-identity Majorana string $\mu(a)$, the random variable $X_a=\langle \psi|\mu(a)|\psi\rangle$ has mean zero and fourth moment of order $D^{-2}$. Since there are $D^2-1$ non-identity strings,
$$
\mathbb{E}_\psi F_2(|\psi\rangle;O)=O(1).
$$
We can get a useful tail bound without invoking a general polynomial concentration theorem. The rotated operators
$U_O^\dagger \mu(a) U_O$ form another complete orthogonal Majorana-string basis. Hence, for every pure state,
$$
\sum_a \left|\langle \psi|U_O^\dagger \mu(a) U_O|\psi\rangle\right|^2=D,
$$
because this is just the Hilbert-Schmidt identity for the expansion of
$|\psi\rangle\langle\psi|$ in an orthogonal operator basis with
$\mathrm{tr}(\mu(a)\mu(b))=D\delta_{ab}$. Therefore, for $t>2$,
$$
F_2(|\psi\rangle;O)>t
\quad\Longrightarrow\quad
\max_{a\neq 0}
\left|\langle \psi|U_O^\dagger \mu(a) U_O|\psi\rangle\right|^2
\gt \frac{t-1}{D-1}.
$$
For each non-identity string, $A=U_O^\dagger \mu(a)U_O$ is a traceless Hermitian unitary. If $|\psi\rangle$ is Haar random, then
$\langle\psi|A|\psi\rangle=2P-1$, where $P$ is the total Haar weight in the positive eigenspace of $A$. Since the positive and negative eigenspaces both have dimension $D/2$, $P$ has the beta distribution $\mathrm{Beta}(D/2,D/2)$, and a standard Chernoff or Levy concentration bound gives
$$
\Pr_\psi\!\left[
\left|\langle\psi|A|\psi\rangle\right|>s
\right]
\leq 2e^{-cDs^2}
$$
for a universal constant $c>0$. Union bounding over the $D^2-1$ non-identity strings gives, for every fixed $O$ and $C\leq t\leq D$,
$$
\Pr_\psi\!\left[F_2(|\psi\rangle;O)>t\right]
\leq C D^2 e^{-c t}.
$$
For $t>D$, the deterministic bound $F_2\leq D$ is already enough. Thus the fixed-$O$ tail is exponentially small in $t$ throughout the only nontrivial regime needed below.

Now cover $O(N)$ by an $\epsilon$-net $\mathcal{N}_\epsilon$ in operator norm. One may take
$$
|\mathcal{N}_\epsilon|\leq \left(\frac{C}{\epsilon}\right)^{N(N-1)/2}.
$$
The function $F_2$ is Lipschitz as a function of $O$. Indeed, a weight-$W$ string transforms in the $W$th exterior-power representation, so changing $O$ by $\epsilon$ changes the rotated string by at most $O(W\epsilon)$ in operator norm. Since $W\leq N$ and there are $2^N$ strings,
$$
|F_2(|\psi\rangle;O)-F_2(|\psi\rangle;O')|
\leq C N 2^N \|O-O'\|
$$
uniformly in $|\psi\rangle$. This bound is crude but sufficient.

Fix any $\delta>0$ and set $t=e^{\delta N}$. Choose
$$
\epsilon \sim \frac{t}{N2^N}.
$$
Then if $\sup_O F_2(|\psi\rangle;O)>t$, some point of the net has $F_2$ larger than a constant multiple of $t$. By the fixed-$O$ tail bound and the union bound,
$$
\Pr_\psi\!\left[\sup_{O\in O(N)}F_2(|\psi\rangle;O)>e^{\delta N}\right]
\leq
\exp\!\left[
\frac{N(N-1)}{2}\log\!\left(\frac{C N2^N}{e^{\delta N}}\right)
-c e^{\delta N}
\right],
$$
up to harmless changes of constants in the exponent and the extra fixed-$O$ prefactor $D^2=e^{O(N)}$. If $e^{\delta N}>D$, the desired bound is instead deterministic. Otherwise, the logarithm of the net size is only polynomial in $N$, while the concentration term is exponentially large in $N$. Therefore the probability tends to zero for every fixed $\delta>0$.

Equivalently, with high probability,
$$
\sup_{O\in O(N)}F_2(|\psi\rangle;O)\leq e^{\delta N}
$$
for every fixed $\delta>0$, which is the desired statement that $\log\sup_O F_2=o(N)$. This proves that the non-Gaussian stabilizer Renyi entropy differs from the ordinary stabilizer Renyi entropy by a subextensive amount for Haar-random pure fermion states.

Agent note: Codex (GPT-5) replaced the Test 2 TODO with this proof argument on May 22, 2026, and later revised the fixed-$O$ tail estimate using the Hilbert-Schmidt identity and a single-observable Haar concentration bound.

## Outlook

Where are we? We defined a notion of non-Gaussian magic and discussed its application to thermal states of the SYK model. Focusing on the 4-body case, we gave two calculations that suggest that SYK-4 states have high non-Gaussian magic. In the first calculation, we explicitly studied how much expectation values could be modified by Gaussian rotations. Based on a small range of $N$ values, we obtained evidence that the $\alpha=2$ SRE per fermion is invariant in the thermodynamic limit. In the second test, we sketched an argument that the SRE per fermion for random states is also invariant in the thermodynamic limit. While SYK energy eigenstates are not literally Haar random, they are plausibly highly random (in the spirit of eigenstate thermalization). In my view, this makes the [connection](https://arxiv.org/abs/2403.07056) between magic and gravitational back-reaction more sensible, since the measures are insensitive to microscopic details that, arguably, are not captured in [coarse-grained gravity calculations](https://arxiv.org/abs/2602.12339).
