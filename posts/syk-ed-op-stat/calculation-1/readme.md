# Calculation 1: ED versus the two-replica saddle

This calculation tests the proposed relation
$$
\frac{1}{N}\sum_i\mathbb E_J[
\zeta_{Ai}(\tau)\zeta_{Ai}(\tau')]
\mathrel{\stackrel{?}{\propto}}G_{12}(\tau,\tau')
$$
for the $q=4$ SYK model. The main comparison uses $N=20$,
$A=\{1,2,3\}$, $\beta J=0.5$, 24 disorder realizations, nine imaginary-time
points, and 80 time slices for the saddle. The coupling convention is
$\mathbb E[J_I^2]=3!J^2/N^3$ for the Hamiltonian and Majorana normalization in
the entry.

`compare_ed_path_integral.py` contains both parts of the calculation:

- Exact diagonalization uses Jordan--Wigner Majoranas and diagonalizes the two
  fermion-parity blocks of dimension 512. It evaluates the trace defining
  $\zeta_{Ai}$ directly, including the factor $\psi_i=\gamma_i/\sqrt{2}$.
- The path-integral part is a NumPy port of the $R=2$, $q=4$
  `WeightedReplicas` iteration in
  [SYKRE.jl](https://github.com/vbettaque/SYKRE.jl), inspected at commit
  `fafda41aa385a89a820b00b18b48a60e25eab453`. It solves
  $\Sigma=J^2G^3$ and uses the weighted periodic/antiperiodic propagator at
  $w=|A|/N=0.15$.

Run the checks and main calculation from this directory with

```bash
python3 compare_ed_path_integral.py --self-test
python3 compare_ed_path_integral.py --outdir outputs
```

The self-test checks the Clifford algebra, replica block signs, finiteness of a
small saddle, and the parity-block trace against a brute-force dense trace.

## Result

The main output is `outputs/comparison.png`; numerical arrays (including all
$\zeta$ samples) are in `outputs/comparison_data.npz`, and scalar diagnostics
are in `outputs/summary.json`. The time-lattice refinement results are recorded
in `outputs/path_convergence.json`.

The result has two parts:

1. The literal raw proportionality is not supported. A one-parameter fit
   through the origin gives cosine similarity $-0.0229$ and relative residual
   $0.9997$. The saddle has a large time-independent component
   ($G_{12}\simeq0.88$), whereas the ED surface changes sign and has nearly zero
   time average.
2. The time-dependent components agree well. After subtracting the time
   average of each surface and fitting one amplitude, the shape correlation is
   $0.9901$ and the relative residual is $0.1407$. The aggregate ED
   signal-to-noise ratio is $11.4$.

The centered saddle shape is numerically stable: using 40, 80, and 120 time
slices gives normalized-shape overlaps of 0.99920 (40 versus 80), 0.99975 (80
versus 120), and 0.99905 (40 versus 120). Thus the raw mismatch is not a
time-discretization effect. It points instead to a missing prescription for
the periodic zero mode or a subtraction/normalization in the proposed
identification. The centered agreement is evidence for the claimed
time-dependent structure, but should not be described as verification of the
unmodified equation.

The follow-up [periodic versus antiperiodic sector
note](sector-decomposition-note.md) localizes the mismatch to the periodic
sector and motivates calculation 2.

Status: implemented and run by Codex (GPT-5) on July 29, 2026. The entry was
also edited with these results.
