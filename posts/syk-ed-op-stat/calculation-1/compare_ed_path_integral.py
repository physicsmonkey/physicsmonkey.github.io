#!/usr/bin/env python3
"""Compare an SYK replica saddle with a finite-N exact-diagonalization probe.

The conventions match ``../syk-ed-op-stat.md``:

    {psi_i, psi_j} = delta_ij,
    H = sum_I i^(q/2) J_I prod_{i in I} psi_i,
    E[J_I^2] = (q-1)! J^2 / N^(q-1),                         (q = 4)

Internally gamma_i = sqrt(2) psi_i, so a normalized Majorana string is
mu_A = i^(W(W-1)/2) prod_{i in A} gamma_i.

The path-integral solver is a NumPy port of the R=2 ``WeightedReplicas``
Schwinger-Dyson iteration in https://github.com/vbettaque/SYKRE.jl.  The ED
calculation uses fermion-parity blocks and evaluates

    zeta_Ai(tau) = Tr[e^(-beta H) psi_i(tau) mu_A] / Z

without constructing dense operator products.  Since the journal claim only
fixes an overall proportionality constant, the two bilocal surfaces are
compared after a least-squares amplitude fit.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


@dataclass
class SaddleResult:
    G_blocks: np.ndarray
    iterations: int
    relative_update: float
    converged: bool


@dataclass
class ComparisonSummary:
    n_majorana: int
    weight: int
    relative_weight: float
    beta_j: float
    disorder_samples: int
    time_points: int
    path_slices: int
    path_iterations: int
    path_relative_update: float
    path_converged: bool
    raw_fitted_amplitude: float
    raw_cosine_similarity: float
    raw_relative_residual: float
    centered_fitted_amplitude: float
    centered_shape_correlation: float
    centered_relative_residual: float
    ed_imaginary_fraction: float
    ed_signal_to_noise: float
    seed: int


def apply_majorana(
    states: np.ndarray, gamma_index: int, n_qubits: int
) -> tuple[np.ndarray, np.ndarray]:
    """Apply a Jordan-Wigner gamma to computational-basis states."""
    site = gamma_index // 2
    bit_shift = n_qubits - 1 - site
    bits = (states >> bit_shift) & 1
    prefix = states >> (n_qubits - site)
    parity = np.fromiter(
        (int(value).bit_count() & 1 for value in prefix),
        dtype=np.int8,
        count=len(states),
    )
    phases = np.where(parity == 0, 1.0 + 0.0j, -1.0 + 0.0j)
    if gamma_index % 2:
        phases *= np.where(bits == 0, 1.0j, -1.0j)
    return states ^ (1 << bit_shift), phases


def string_action(
    states: np.ndarray, indices: tuple[int, ...], n_majorana: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return full-basis target states and phases for mu_indices |state>."""
    targets = states.copy()
    phases = np.ones(len(states), dtype=complex)
    # Operators in a matrix product act on kets from right to left.
    for gamma_index in reversed(indices):
        targets, extra = apply_majorana(targets, gamma_index, n_majorana // 2)
        phases *= extra
    phases *= (1j) ** (len(indices) * (len(indices) - 1) // 2)
    return targets, phases


class ParityBasis:
    """Even/odd computational bases and sparse Majorana-string actions."""

    def __init__(self, n_majorana: int):
        if n_majorana % 2:
            raise ValueError("N must be even")
        self.n_majorana = n_majorana
        self.dim = 2 ** (n_majorana // 2)
        states = np.arange(self.dim, dtype=np.int64)
        parity = np.fromiter(
            (int(state).bit_count() & 1 for state in states),
            dtype=np.int8,
            count=self.dim,
        )
        self.even = states[parity == 0]
        self.odd = states[parity == 1]
        self.local_even = np.full(self.dim, -1, dtype=np.int64)
        self.local_odd = np.full(self.dim, -1, dtype=np.int64)
        self.local_even[self.even] = np.arange(len(self.even))
        self.local_odd[self.odd] = np.arange(len(self.odd))

    def action_between(
        self, indices: tuple[int, ...], source_parity: int
    ) -> tuple[np.ndarray, np.ndarray]:
        sources = self.even if source_parity == 0 else self.odd
        targets, phases = string_action(sources, indices, self.n_majorana)
        target_local = (
            self.local_odd[targets]
            if (source_parity + len(indices)) % 2
            else self.local_even[targets]
        )
        if np.any(target_local < 0):
            raise RuntimeError("Majorana string did not have the expected parity")
        return target_local, phases


def build_hamiltonian_actions(
    basis: ParityBasis,
) -> tuple[list[tuple[int, ...]], list[tuple[np.ndarray, np.ndarray]], list[tuple[np.ndarray, np.ndarray]]]:
    strings = list(itertools.combinations(range(basis.n_majorana), 4))
    even_actions = [basis.action_between(indices, 0) for indices in strings]
    odd_actions = [basis.action_between(indices, 1) for indices in strings]
    return strings, even_actions, odd_actions


def sample_hamiltonian_blocks(
    rng: np.random.Generator,
    basis: ParityBasis,
    even_actions: list[tuple[np.ndarray, np.ndarray]],
    odd_actions: list[tuple[np.ndarray, np.ndarray]],
    coupling_j: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample q=4 SYK and return its even and odd parity blocks."""
    # H = sum J_I mu_I / 4 for q=4 in the gamma normalization.
    coefficient_std = (
        coupling_j
        * math.sqrt(math.factorial(3) / basis.n_majorana**3)
        / 4.0
    )
    coefficients = rng.normal(0.0, coefficient_std, len(even_actions))
    block_dim = len(basis.even)
    columns = np.arange(block_dim)
    h_even = np.zeros((block_dim, block_dim), dtype=complex)
    h_odd = np.zeros_like(h_even)
    for coefficient, (even_action, odd_action) in zip(
        coefficients, zip(even_actions, odd_actions)
    ):
        even_rows, even_phases = even_action
        odd_rows, odd_phases = odd_action
        h_even[even_rows, columns] += coefficient * even_phases
        h_odd[odd_rows, columns] += coefficient * odd_phases
    # This also removes harmless roundoff asymmetry.
    h_even = 0.5 * (h_even + h_even.conj().T)
    h_odd = 0.5 * (h_odd + h_odd.conj().T)
    return h_even, h_odd


def thermal_kernels(
    energies: np.ndarray,
    vectors: np.ndarray,
    times: np.ndarray,
) -> list[np.ndarray]:
    """Return exp(-t H) at each requested t, with shifted energies."""
    return [
        (vectors * np.exp(-time * energies)) @ vectors.conj().T
        for time in times
    ]


def zeta_grid(
    h_even: np.ndarray,
    h_odd: np.ndarray,
    basis: ParityBasis,
    measured_indices: tuple[int, ...],
    beta: float,
    times: np.ndarray,
) -> np.ndarray:
    """Compute zeta_Ai(tau) for every i and every tau in ``times``."""
    energies_even, vectors_even = np.linalg.eigh(h_even)
    energies_odd, vectors_odd = np.linalg.eigh(h_odd)
    energy_shift = min(energies_even[0], energies_odd[0])
    energies_even -= energy_shift
    energies_odd -= energy_shift

    # A grid invariant under tau -> beta-tau lets the same kernels serve both
    # factors in the trace.  Include beta, which is needed when tau=0.
    kernel_times = np.unique(np.concatenate([times, beta - times]))
    even_kernels = thermal_kernels(energies_even, vectors_even, kernel_times)
    odd_kernels = thermal_kernels(energies_odd, vectors_odd, kernel_times)
    even_at = {float(time): kernel for time, kernel in zip(kernel_times, even_kernels)}
    odd_at = {float(time): kernel for time, kernel in zip(kernel_times, odd_kernels)}

    partition_function = float(
        np.exp(-beta * energies_even).sum()
        + np.exp(-beta * energies_odd).sum()
    )
    mu_even_rows, mu_even_phases = basis.action_between(measured_indices, 0)
    mu_odd_rows, mu_odd_phases = basis.action_between(measured_indices, 1)

    values = np.empty((basis.n_majorana, len(times)), dtype=complex)
    for gamma_index in range(basis.n_majorana):
        gamma_even_rows, gamma_even_phases = basis.action_between(
            (gamma_index,), 0
        )
        gamma_odd_rows, gamma_odd_phases = basis.action_between(
            (gamma_index,), 1
        )
        gamma_even_phases = gamma_even_phases / math.sqrt(2.0)
        gamma_odd_phases = gamma_odd_phases / math.sqrt(2.0)

        for time_index, tau in enumerate(times):
            # Trace beginning in the even sector:
            # sum_{a even,c odd} K_e(beta-tau)[a,p_i(c)] p_i(c)
            #                      K_o(tau)[c,p_A(a)] p_A(a).
            even_left = even_at[float(beta - tau)][:, gamma_odd_rows]
            odd_right = odd_at[float(tau)][:, mu_even_rows].T
            trace_even = np.sum(
                even_left
                * odd_right
                * mu_even_phases[:, None]
                * gamma_odd_phases[None, :]
            )

            odd_left = odd_at[float(beta - tau)][:, gamma_even_rows]
            even_right = even_at[float(tau)][:, mu_odd_rows].T
            trace_odd = np.sum(
                odd_left
                * even_right
                * mu_odd_phases[:, None]
                * gamma_even_phases[None, :]
            )
            values[gamma_index, time_index] = (
                trace_even + trace_odd
            ) / partition_function
    return values


def replica_full_matrix(blocks: np.ndarray) -> np.ndarray:
    """Convert R=2 replica-circulant blocks to [[A,-B],[B,A]]."""
    diagonal, off_diagonal = blocks
    return np.block(
        [[diagonal, -off_diagonal], [off_diagonal, diagonal]]
    )


def replica_blocks(matrix: np.ndarray) -> np.ndarray:
    length = matrix.shape[0] // 2
    return np.stack([matrix[:length, :length], matrix[length:, :length]])


def discrete_derivative(length: int, periodic: bool) -> np.ndarray:
    derivative = np.eye(length)
    derivative[np.arange(1, length), np.arange(length - 1)] = -1.0
    derivative[0, -1] = -1.0 if periodic else 1.0
    zeros = np.zeros_like(derivative)
    return np.block([[derivative, zeros], [zeros, derivative]])


def initial_replica_green(length: int) -> np.ndarray:
    diagonal = np.empty((length, length))
    rows, columns = np.indices((length, length))
    diagonal[:] = 0.5 * np.sign(rows - columns)
    off_diagonal = np.full((length, length), 0.5)
    return np.stack([diagonal, off_diagonal])


def solve_weighted_replica_saddle(
    beta: float,
    relative_weight: float,
    length: int,
    coupling_j: float = 1.0,
    mixing: float = 0.01,
    tolerance: float = 1e-8,
    max_iterations: int = 4000,
) -> SaddleResult:
    """Port of SYKRE.jl's q=4, R=2 WeightedReplicas iteration."""
    if not 0.0 <= relative_weight <= 1.0:
        raise ValueError("relative weight must lie in [0,1]")
    green = initial_replica_green(length)
    derivative_minus = discrete_derivative(length, periodic=False)
    derivative_plus = discrete_derivative(length, periodic=True)
    delta_tau_squared = (beta / length) ** 2
    relative_update = math.inf

    for iteration in range(1, max_iterations + 1):
        sigma = coupling_j**2 * green**3
        sigma_full = replica_full_matrix(sigma)
        propagator_minus = derivative_minus - delta_tau_squared * sigma_full
        propagator_plus = derivative_plus - delta_tau_squared * sigma_full
        new_full = -(
            relative_weight * np.linalg.inv(propagator_plus)
            + (1.0 - relative_weight) * np.linalg.inv(propagator_minus)
        ).T
        new_green = replica_blocks(new_full)
        mixed = mixing * new_green + (1.0 - mixing) * green
        relative_update = float(
            np.linalg.norm(mixed - green) / np.linalg.norm(green)
        )
        green = mixed
        if relative_update < tolerance:
            return SaddleResult(green, iteration, relative_update, True)
    return SaddleResult(
        green, max_iterations, relative_update, False
    )


def fit_surfaces(
    ed_surface: np.ndarray, path_surface: np.ndarray
) -> tuple[float, float, float]:
    ed = np.asarray(ed_surface.real).ravel()
    path = np.asarray(path_surface).ravel()
    amplitude = float(np.dot(ed, path) / np.dot(path, path))
    residual = ed - amplitude * path
    relative_residual = float(np.linalg.norm(residual) / np.linalg.norm(ed))
    cosine = float(np.dot(ed, path) / (np.linalg.norm(ed) * np.linalg.norm(path)))
    return amplitude, cosine, relative_residual


def save_figure(
    output_path: Path,
    fractions: np.ndarray,
    path_surface: np.ndarray,
    ed_surface: np.ndarray,
    ed_standard_error: np.ndarray,
    centered_amplitude: float,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    extent = [
        float(fractions[0]),
        float(fractions[-1]),
        float(fractions[0]),
        float(fractions[-1]),
    ]
    centered_path = path_surface - path_surface.mean()
    centered_ed = ed_surface.real - ed_surface.real.mean()
    fitted_path = centered_amplitude * centered_path
    difference = centered_ed - fitted_path
    scale = max(np.max(np.abs(ed_surface.real)), np.max(np.abs(fitted_path)))

    figure, axes = plt.subplots(1, 5, figsize=(16.0, 3.1), constrained_layout=True)
    panels = [
        (path_surface, r"raw saddle $G_{12}$", "viridis", None, None),
        (fitted_path, r"fitted $a(G_{12}-\overline{G}_{12})$", "RdBu_r", -scale, scale),
        (ed_surface.real, r"ED $\mathbb{E}[\sum_i\zeta_i\zeta_i']/N$", "RdBu_r", -scale, scale),
        (difference, "centered ED minus fit", "RdBu_r", -scale, scale),
        (ed_standard_error, "ED standard error", "viridis", 0.0, None),
    ]
    for axis, (data, title, cmap, vmin, vmax) in zip(axes, panels):
        image = axis.imshow(
            data,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        axis.set_title(title, fontsize=10)
        axis.set_xlabel(r"$\tau'/\beta$")
        axis.set_ylabel(r"$\tau/\beta$")
        figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def run(args: argparse.Namespace) -> ComparisonSummary:
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    basis = ParityBasis(args.n_majorana)
    if args.weight % 2 != 1:
        raise ValueError("The measured string weight must be odd")
    if args.weight > args.n_majorana:
        raise ValueError("The measured string weight cannot exceed N")
    measured_indices = tuple(range(args.weight))
    _, even_actions, odd_actions = build_hamiltonian_actions(basis)

    fractions = np.arange(args.time_points, dtype=float) / args.time_points
    times = args.beta * fractions
    samples = np.empty(
        (args.samples, args.n_majorana, args.time_points), dtype=complex
    )
    for sample_index in range(args.samples):
        h_even, h_odd = sample_hamiltonian_blocks(
            rng, basis, even_actions, odd_actions, args.coupling_j
        )
        samples[sample_index] = zeta_grid(
            h_even,
            h_odd,
            basis,
            measured_indices,
            args.beta,
            times,
        )
        print(
            f"ED sample {sample_index + 1}/{args.samples}",
            flush=True,
        )

    sample_surfaces = np.einsum(
        "sit,siu->stu", samples, samples, optimize=True
    ) / args.n_majorana
    ed_surface = sample_surfaces.mean(axis=0)
    ed_standard_error = (
        sample_surfaces.real.std(axis=0, ddof=1) / math.sqrt(args.samples)
        if args.samples > 1
        else np.full(ed_surface.shape, np.nan)
    )

    saddle = solve_weighted_replica_saddle(
        beta=args.beta,
        relative_weight=args.weight / args.n_majorana,
        length=args.path_slices,
        coupling_j=args.coupling_j,
        mixing=args.path_mixing,
        tolerance=args.path_tolerance,
        max_iterations=args.path_max_iterations,
    )
    path_indices = np.rint(fractions * args.path_slices).astype(int)
    # SYKRE.jl calls blocks[:,:,2] G12 in its R=2 saddle notebook.
    path_surface_lattice = saddle.G_blocks[1][np.ix_(path_indices, path_indices)]
    # At finite lattice spacing the forward-difference discretization violates
    # G12(t,t') = G12(t',t) by O(Delta tau).  Restore the exact replica symmetry
    # before comparing with the manifestly symmetric ED product.
    path_surface = 0.5 * (path_surface_lattice + path_surface_lattice.T)
    raw_amplitude, raw_cosine, raw_residual = fit_surfaces(
        ed_surface, path_surface
    )
    centered_ed = ed_surface.real - ed_surface.real.mean()
    centered_path = path_surface - path_surface.mean()
    centered_amplitude, centered_correlation, centered_residual = fit_surfaces(
        centered_ed, centered_path
    )
    imaginary_fraction = float(
        np.linalg.norm(ed_surface.imag) / np.linalg.norm(ed_surface.real)
    )
    signal_to_noise = float(
        np.linalg.norm(ed_surface.real)
        / np.linalg.norm(ed_standard_error)
    )

    summary = ComparisonSummary(
        n_majorana=args.n_majorana,
        weight=args.weight,
        relative_weight=args.weight / args.n_majorana,
        beta_j=args.beta * args.coupling_j,
        disorder_samples=args.samples,
        time_points=args.time_points,
        path_slices=args.path_slices,
        path_iterations=saddle.iterations,
        path_relative_update=saddle.relative_update,
        path_converged=saddle.converged,
        raw_fitted_amplitude=raw_amplitude,
        raw_cosine_similarity=raw_cosine,
        raw_relative_residual=raw_residual,
        centered_fitted_amplitude=centered_amplitude,
        centered_shape_correlation=centered_correlation,
        centered_relative_residual=centered_residual,
        ed_imaginary_fraction=imaginary_fraction,
        ed_signal_to_noise=signal_to_noise,
        seed=args.seed,
    )

    np.savez_compressed(
        args.outdir / "comparison_data.npz",
        fractions=fractions,
        times=times,
        zeta_samples=samples,
        ed_surface=ed_surface,
        ed_standard_error=ed_standard_error,
        path_surface=path_surface,
        path_surface_lattice=path_surface_lattice,
        path_G_blocks=saddle.G_blocks,
    )
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(asdict(summary), handle, indent=2)
    save_figure(
        args.outdir / "comparison.png",
        fractions,
        path_surface,
        ed_surface,
        ed_standard_error,
        centered_amplitude,
    )
    print(json.dumps(asdict(summary), indent=2))
    return summary


def self_test() -> None:
    basis = ParityBasis(6)
    states = np.arange(basis.dim)
    for gamma_index in range(6):
        twice, phase_2 = string_action(states, (gamma_index, gamma_index), 6)
        if not np.array_equal(twice, states) or not np.allclose(phase_2, 1j):
            # string_action includes the normalized-string phase i for W=2.
            raise AssertionError("unexpected two-gamma string convention")
        target, phase = apply_majorana(states, gamma_index, 3)
        back, back_phase = apply_majorana(target, gamma_index, 3)
        if not np.array_equal(back, states) or not np.allclose(
            phase * back_phase, 1.0
        ):
            raise AssertionError("gamma_i^2 != 1")

    saddle = solve_weighted_replica_saddle(
        beta=0.5,
        relative_weight=0.15,
        length=12,
        mixing=0.02,
        tolerance=1e-6,
        max_iterations=1500,
    )
    full_initial = replica_full_matrix(initial_replica_green(12))
    if not np.allclose(full_initial + full_initial.T, 0.0):
        raise AssertionError("replica block conversion has the wrong signs")
    if not np.all(np.isfinite(saddle.G_blocks)):
        raise AssertionError("replica iteration produced non-finite values")

    # Compare the parity-block trace contraction with a direct dense trace.
    beta = 0.5
    times = np.array([0.0, 0.17, 0.4])
    _, even_actions, odd_actions = build_hamiltonian_actions(basis)
    h_even, h_odd = sample_hamiltonian_blocks(
        np.random.default_rng(2), basis, even_actions, odd_actions, 1.0
    )
    zeta = zeta_grid(h_even, h_odd, basis, (0, 1, 2), beta, times)
    h_full = np.zeros((basis.dim, basis.dim), dtype=complex)
    h_full[np.ix_(basis.even, basis.even)] = h_even
    h_full[np.ix_(basis.odd, basis.odd)] = h_odd
    energies, vectors = np.linalg.eigh(h_full)
    energies -= energies[0]

    def dense_string(indices: tuple[int, ...]) -> np.ndarray:
        source = np.arange(basis.dim)
        target, phase = string_action(source, indices, basis.n_majorana)
        operator = np.zeros((basis.dim, basis.dim), dtype=complex)
        operator[target, source] = phase
        return operator

    kernels = {
        float(time): (vectors * np.exp(-time * energies)) @ vectors.conj().T
        for time in np.unique(np.concatenate([times, beta - times]))
    }
    partition_function = np.exp(-beta * energies).sum()
    mu = dense_string((0, 1, 2))
    for gamma_index in range(basis.n_majorana):
        psi = dense_string((gamma_index,)) / math.sqrt(2.0)
        for time_index, tau in enumerate(times):
            direct = (
                np.trace(kernels[float(beta - tau)] @ psi @ kernels[float(tau)] @ mu)
                / partition_function
            )
            if not np.allclose(zeta[gamma_index, time_index], direct, atol=2e-14):
                raise AssertionError("parity-block zeta does not match dense trace")
    print("self-test passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-majorana", type=int, default=20)
    parser.add_argument("--weight", type=int, default=3)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--coupling-j", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=24)
    parser.add_argument("--time-points", type=int, default=9)
    parser.add_argument("--path-slices", type=int, default=80)
    parser.add_argument("--path-mixing", type=float, default=0.01)
    parser.add_argument("--path-tolerance", type=float, default=1e-8)
    parser.add_argument("--path-max-iterations", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--outdir", type=Path, default=Path("outputs"))
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
    else:
        run(args)


if __name__ == "__main__":
    main()
