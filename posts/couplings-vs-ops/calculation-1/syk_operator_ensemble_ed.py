#!/usr/bin/env python3
"""Exact diagonalization tests for coupling vs operator ensembles in SYK.

The script builds dense Majorana matrices for small even N, samples p=4 SYK
Hamiltonians, computes thermal one-point functions for all strings of a fixed
weight W, and compares the fixed-H operator averages M1 and M2 to the disorder
ensemble expectations described in the journal entry.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


PAULI_I = np.array([[1, 0], [0, 1]], dtype=complex)
PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)


@dataclass
class SummaryRow:
    n_majorana: int
    weight: int
    beta: float
    n_samples: int
    n_strings: int
    mean_m1: float
    std_m1: float
    predicted_std_m1: float
    std_m1_over_prediction: float
    mean_m2: float
    sigma2_mean_over_strings: float
    rel_std_m2: float
    predicted_rel_std_m2: float
    rel_std_m2_over_prediction: float
    mean_abs_offdiag_corr: float
    max_abs_offdiag_corr: float
    mean_excess_kurtosis: float


def kron_all(factors: list[np.ndarray]) -> np.ndarray:
    out = factors[0]
    for factor in factors[1:]:
        out = np.kron(out, factor)
    return out


def majoranas(n_majorana: int) -> list[np.ndarray]:
    """Jordan-Wigner Majoranas with {psi_i, psi_j}=2 delta_ij."""
    if n_majorana % 2:
        raise ValueError("n_majorana must be even")
    n_qubits = n_majorana // 2
    gammas: list[np.ndarray] = []
    for site in range(n_qubits):
        z_string = [PAULI_Z] * site
        tail = [PAULI_I] * (n_qubits - site - 1)
        gammas.append(kron_all(z_string + [PAULI_X] + tail))
        gammas.append(kron_all(z_string + [PAULI_Y] + tail))
    return gammas


def string_operator(gammas: list[np.ndarray], indices: tuple[int, ...]) -> np.ndarray:
    dim = gammas[0].shape[0]
    op = np.eye(dim, dtype=complex)
    for index in indices:
        op = op @ gammas[index]
    phase = (1j) ** (len(indices) * (len(indices) - 1) // 2)
    op = phase * op
    return 0.5 * (op + op.conj().T)


def apply_majorana_to_basis_state(
    state: int, gamma_index: int, n_qubits: int
) -> tuple[int, complex]:
    """Apply one Jordan-Wigner Majorana to a computational basis state.

    The dense construction above uses kron factors ordered by site, so site 0
    is the most significant bit of the integer basis label.
    """
    site = gamma_index // 2
    bit_shift = n_qubits - 1 - site
    bit = (state >> bit_shift) & 1
    prefix_parity = 0
    for q in range(site):
        prefix_parity ^= (state >> (n_qubits - 1 - q)) & 1

    phase = -1.0 + 0.0j if prefix_parity else 1.0 + 0.0j
    if gamma_index % 2:
        phase *= -1.0j if bit else 1.0j
    return state ^ (1 << bit_shift), phase


def string_action(
    n_majorana: int, indices: tuple[int, ...]
) -> tuple[np.ndarray, np.ndarray]:
    """Return perm, phase for mu(indices)|b> = phase[b]|perm[b]>."""
    n_qubits = n_majorana // 2
    dim = 2**n_qubits
    perm = np.empty(dim, dtype=np.int64)
    phases = np.empty(dim, dtype=complex)
    hermitian_phase = (1j) ** (len(indices) * (len(indices) - 1) // 2)

    for basis_state in range(dim):
        state = basis_state
        phase = 1.0 + 0.0j
        # Matrix product gamma_i gamma_j acts on kets from right to left.
        for gamma_index in reversed(indices):
            state, extra_phase = apply_majorana_to_basis_state(
                state, gamma_index, n_qubits
            )
            phase *= extra_phase
        perm[basis_state] = state
        phases[basis_state] = hermitian_phase * phase
    return perm, phases


def action_to_dense(perm: np.ndarray, phases: np.ndarray) -> np.ndarray:
    dim = len(perm)
    op = np.zeros((dim, dim), dtype=complex)
    op[perm, np.arange(dim)] = phases
    return op


def build_string_actions(
    n_majorana: int, weight: int
) -> tuple[list[tuple[int, ...]], list[tuple[np.ndarray, np.ndarray]]]:
    strings = list(itertools.combinations(range(n_majorana), weight))
    actions = [string_action(n_majorana, inds) for inds in strings]
    return strings, actions


def sample_hamiltonian(
    rng: np.random.Generator,
    h_actions: list[tuple[np.ndarray, np.ndarray]],
    n_majorana: int,
    coupling_scale: float,
) -> np.ndarray:
    # From the entry: E[J_I^2] = (p-1)! / N^(p-1) / 2^p for p=4.
    coupling_std = coupling_scale * math.sqrt(math.factorial(3) / n_majorana**3 / 2**4)
    coeffs = rng.normal(0.0, coupling_std, size=len(h_actions))
    dim = len(h_actions[0][0])
    cols = np.arange(dim)
    hamiltonian = np.zeros((dim, dim), dtype=complex)
    for coeff, (perm, phases) in zip(coeffs, h_actions):
        hamiltonian[perm, cols] += coeff * phases
    return 0.5 * (hamiltonian + hamiltonian.conj().T)


def thermal_one_points(
    hamiltonian: np.ndarray,
    string_actions: list[tuple[np.ndarray, np.ndarray]],
    beta: float,
) -> np.ndarray:
    energies, vectors = np.linalg.eigh(hamiltonian)
    shifted = energies - energies.min()
    boltz = np.exp(-beta * shifted)
    rho = (vectors * boltz) @ vectors.conj().T / boltz.sum()
    values = np.array(
        [np.sum(rho[np.arange(len(perm)), perm] * phases).real for perm, phases in string_actions]
    )
    return values


def offdiag_corr_stats(values_by_sample: np.ndarray) -> tuple[float, float]:
    if values_by_sample.shape[0] < 3 or values_by_sample.shape[1] < 2:
        return float("nan"), float("nan")
    corr = np.corrcoef(values_by_sample, rowvar=False)
    mask = ~np.eye(corr.shape[0], dtype=bool)
    offdiag = np.abs(corr[mask])
    return float(np.mean(offdiag)), float(np.max(offdiag))


def excess_kurtosis(values_by_sample: np.ndarray) -> float:
    centered = values_by_sample - values_by_sample.mean(axis=0, keepdims=True)
    var = np.mean(centered**2, axis=0)
    good = var > 1e-20
    if not np.any(good):
        return float("nan")
    fourth = np.mean(centered[:, good] ** 4, axis=0)
    return float(np.mean(fourth / var[good] ** 2 - 3.0))


def run_case(
    rng: np.random.Generator,
    n_majorana: int,
    weight: int,
    beta: float,
    n_samples: int,
    coupling_scale: float,
) -> tuple[SummaryRow, np.ndarray]:
    _, h_actions = build_string_actions(n_majorana, 4)
    _, measured_actions = build_string_actions(n_majorana, weight)

    values = []
    for _ in range(n_samples):
        hamiltonian = sample_hamiltonian(rng, h_actions, n_majorana, coupling_scale)
        values.append(thermal_one_points(hamiltonian, measured_actions, beta))
    values_by_sample = np.vstack(values)

    m1 = values_by_sample.mean(axis=1)
    m2 = np.mean(values_by_sample**2, axis=1)
    sigma2_by_string = np.mean(values_by_sample**2, axis=0)
    sigma2_mean = float(np.mean(sigma2_by_string))
    n_strings = values_by_sample.shape[1]

    std_m1 = float(np.std(m1, ddof=1))
    predicted_std_m1 = math.sqrt(sigma2_mean / n_strings)
    std_m1_ratio = std_m1 / predicted_std_m1 if predicted_std_m1 > 0 else float("nan")

    mean_m2 = float(np.mean(m2))
    rel_std_m2 = float(np.std(m2, ddof=1) / mean_m2) if mean_m2 > 0 else float("nan")
    predicted_rel_std_m2 = math.sqrt(2.0 / n_strings)
    rel_std_m2_ratio = rel_std_m2 / predicted_rel_std_m2

    mean_abs_corr, max_abs_corr = offdiag_corr_stats(values_by_sample)
    kurt = excess_kurtosis(values_by_sample)

    row = SummaryRow(
        n_majorana=n_majorana,
        weight=weight,
        beta=beta,
        n_samples=n_samples,
        n_strings=n_strings,
        mean_m1=float(np.mean(m1)),
        std_m1=std_m1,
        predicted_std_m1=predicted_std_m1,
        std_m1_over_prediction=std_m1_ratio,
        mean_m2=mean_m2,
        sigma2_mean_over_strings=sigma2_mean,
        rel_std_m2=rel_std_m2,
        predicted_rel_std_m2=predicted_rel_std_m2,
        rel_std_m2_over_prediction=rel_std_m2_ratio,
        mean_abs_offdiag_corr=mean_abs_corr,
        max_abs_offdiag_corr=max_abs_corr,
        mean_excess_kurtosis=kurt,
    )
    return row, values_by_sample


def write_csv(path: Path, rows: list[SummaryRow]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-values", type=int, nargs="+", default=[8, 10, 12])
    parser.add_argument("--weights", type=int, nargs="+", default=[4])
    parser.add_argument("--beta", type=float, default=4.0)
    parser.add_argument("--samples", type=int, default=80)
    parser.add_argument("--coupling-scale", type=float, default=4.0)
    parser.add_argument("--seed", type=int, default=20260606)
    parser.add_argument("--outdir", type=Path, default=Path("outputs"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows: list[SummaryRow] = []
    for n_majorana in args.n_values:
        for weight in args.weights:
            row, values = run_case(
                rng=rng,
                n_majorana=n_majorana,
                weight=weight,
                beta=args.beta,
                n_samples=args.samples,
                coupling_scale=args.coupling_scale,
            )
            rows.append(row)
            np.save(
                args.outdir / f"xi_N{n_majorana}_W{weight}_beta{args.beta:g}.npy",
                values,
            )
            print(json.dumps(asdict(row), indent=2))

    write_csv(args.outdir / "summary.csv", rows)
    with (args.outdir / "run_config.json").open("w") as handle:
        json.dump(vars(args) | {"outdir": str(args.outdir)}, handle, indent=2)


if __name__ == "__main__":
    main()
