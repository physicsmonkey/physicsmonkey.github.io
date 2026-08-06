#!/usr/bin/env python3
"""Direct, normalization-fixed tests of the conditional G_11 propagators.

We use an even parent string A={1,2,3,4}. Two identical probe fermions are
inserted on replica 1, so every trace is fermion-parity even. Probes outside A
are compared with the AP conditional G_11; probes inside A are compared with
the P conditional G_11.

The ED observable is the annealed ratio

  E[X_A;i^(11)(tau,tau') X_A] / E[X_A^2],

which is the direct derivative of the two-replica moment and fixes the absolute
normalization. No amplitude is fitted.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CALCULATION_1 = HERE.parent / "calculation-1"


def load_calculation_1_module():
    module_path = CALCULATION_1 / "compare_ed_path_integral.py"
    spec = importlib.util.spec_from_file_location("compare_ed_path_integral", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def dense_action(
    rows: np.ndarray, phases: np.ndarray
) -> np.ndarray:
    matrix = np.zeros((len(rows), len(rows)), dtype=complex)
    matrix[rows, np.arange(len(rows))] = phases
    return matrix


def kernel_from_eigensystem(
    energies: np.ndarray, vectors: np.ndarray, time: float
) -> np.ndarray:
    return (vectors * np.exp(-time * energies)) @ vectors.conj().T


def conjugated_kernel_average(
    kernel_opposite: np.ndarray,
    actions_from_parity: list[tuple[np.ndarray, np.ndarray]],
    actions_from_opposite: list[tuple[np.ndarray, np.ndarray]],
    indices: range,
) -> np.ndarray:
    """Average psi_i K psi_i in one computational-basis parity block."""
    result = np.zeros_like(kernel_opposite)
    columns = np.arange(kernel_opposite.shape[0])
    for index in indices:
        opposite_rows, phases_from_parity = actions_from_parity[index]
        parity_rows, phases_from_opposite = actions_from_opposite[index]
        # Right multiplication by psi maps parity columns to opposite rows;
        # left multiplication maps opposite rows back to parity rows.
        temp = (
            kernel_opposite[:, opposite_rows]
            * phases_from_parity[None, :]
            / math.sqrt(2.0)
        )
        contribution = np.empty_like(temp)
        contribution[parity_rows, :] = (
            phases_from_opposite[:, None] * temp / math.sqrt(2.0)
        )
        result += contribution
    return result / len(indices)


def trace_from_eigenbasis(
    energies: np.ndarray,
    q_eigen: np.ndarray,
    mu_eigen: np.ndarray,
    beta_minus_tau: float,
    tau_prime: float,
) -> complex:
    left = np.exp(-beta_minus_tau * energies)
    right = np.exp(-tau_prime * energies)
    return np.sum(
        left[:, None] * q_eigen * right[None, :] * mu_eigen.T
    )


def same_side_sample(
    calc1,
    h_even: np.ndarray,
    h_odd: np.ndarray,
    basis,
    measured_indices: tuple[int, ...],
    beta: float,
    times: np.ndarray,
) -> tuple[complex, np.ndarray, np.ndarray]:
    """Return X_A and the inside/outside time-ordered trace surfaces."""
    energies_even, vectors_even = np.linalg.eigh(h_even)
    energies_odd, vectors_odd = np.linalg.eigh(h_odd)

    mu_even_rows, mu_even_phases = basis.action_between(measured_indices, 0)
    mu_odd_rows, mu_odd_phases = basis.action_between(measured_indices, 1)
    mu_even = dense_action(mu_even_rows, mu_even_phases)
    mu_odd = dense_action(mu_odd_rows, mu_odd_phases)
    mu_even_eigen = vectors_even.conj().T @ mu_even @ vectors_even
    mu_odd_eigen = vectors_odd.conj().T @ mu_odd @ vectors_odd

    x_a = (
        np.sum(np.exp(-beta * energies_even) * np.diag(mu_even_eigen))
        + np.sum(np.exp(-beta * energies_odd) * np.diag(mu_odd_eigen))
    )

    actions_even = [
        basis.action_between((index,), 0)
        for index in range(basis.n_majorana)
    ]
    actions_odd = [
        basis.action_between((index,), 1)
        for index in range(basis.n_majorana)
    ]
    inside_indices = range(len(measured_indices))
    outside_indices = range(len(measured_indices), basis.n_majorana)
    inside = np.zeros((len(times), len(times)), dtype=complex)
    outside = np.zeros_like(inside)

    # The code's forward-difference convention takes the negative one-sided
    # contact value on the matrix diagonal.
    np.fill_diagonal(inside, -0.5 * x_a)
    np.fill_diagonal(outside, -0.5 * x_a)

    spacings = np.diff(times)
    if len(spacings) and not np.allclose(spacings, spacings[0]):
        raise ValueError("same_side_sample currently requires a uniform time grid")
    delta_values = {
        delta_index: times[delta_index] - times[0]
        for delta_index in range(1, len(times))
    }
    q_eigen: dict[tuple[str, int, int], np.ndarray] = {}
    for delta_index, delta in delta_values.items():
        kernel_even = kernel_from_eigensystem(
            energies_even, vectors_even, delta
        )
        kernel_odd = kernel_from_eigensystem(
            energies_odd, vectors_odd, delta
        )
        for label, indices in (
            ("inside", inside_indices),
            ("outside", outside_indices),
        ):
            # Result in the even block has an odd intermediate state.
            q_even = conjugated_kernel_average(
                kernel_odd, actions_even, actions_odd, indices
            )
            q_odd = conjugated_kernel_average(
                kernel_even, actions_odd, actions_even, indices
            )
            q_eigen[(label, 0, delta_index)] = (
                vectors_even.conj().T @ q_even @ vectors_even
            )
            q_eigen[(label, 1, delta_index)] = (
                vectors_odd.conj().T @ q_odd @ vectors_odd
            )

    for time_index, tau in enumerate(times):
        for prime_index in range(time_index):
            tau_prime = times[prime_index]
            delta_index = time_index - prime_index
            for label, surface in (("inside", inside), ("outside", outside)):
                value = trace_from_eigenbasis(
                    energies_even,
                    q_eigen[(label, 0, delta_index)],
                    mu_even_eigen,
                    beta - tau,
                    tau_prime,
                )
                value += trace_from_eigenbasis(
                    energies_odd,
                    q_eigen[(label, 1, delta_index)],
                    mu_odd_eigen,
                    beta - tau,
                    tau_prime,
                )
                surface[time_index, prime_index] = value
                surface[prime_index, time_index] = -value
    return x_a, inside, outside


def ratio_and_jackknife(
    x_a_samples: np.ndarray,
    trace_samples: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    denominator_samples = np.real(x_a_samples**2)
    numerator_samples = np.real(
        trace_samples * x_a_samples[:, None, None]
    )
    ratio = numerator_samples.mean(axis=0) / denominator_samples.mean()
    leave_one_out = []
    for omitted in range(len(x_a_samples)):
        numerator = np.delete(numerator_samples, omitted, axis=0).mean(axis=0)
        denominator = np.delete(
            denominator_samples, omitted, axis=0
        ).mean()
        leave_one_out.append(numerator / denominator)
    leave_one_out_array = np.asarray(leave_one_out)
    mean = leave_one_out_array.mean(axis=0)
    standard_error = np.sqrt(
        (len(leave_one_out_array) - 1)
        * np.mean((leave_one_out_array - mean) ** 2, axis=0)
    )
    return ratio, standard_error


def comparison_metrics(
    ed: np.ndarray,
    path: np.ndarray,
    standard_error: np.ndarray,
) -> dict[str, float]:
    mask = ~np.eye(ed.shape[0], dtype=bool)
    difference = ed - path
    return {
        "ed_imaginary_fraction": 0.0,
        "relative_residual_all_entries": float(
            np.linalg.norm(difference) / np.linalg.norm(ed)
        ),
        "relative_residual_off_diagonal": float(
            np.linalg.norm(difference[mask]) / np.linalg.norm(ed[mask])
        ),
        "difference_to_standard_error_all_entries": float(
            np.linalg.norm(difference) / np.linalg.norm(standard_error)
        ),
        "off_diagonal_cosine_similarity": float(
            np.sum(ed[mask] * path[mask])
            / (np.linalg.norm(ed[mask]) * np.linalg.norm(path[mask]))
        ),
        "ed_contact_mean": float(np.diag(ed).mean()),
        "path_contact_mean": float(np.diag(path).mean()),
        "max_absolute_difference": float(np.max(np.abs(difference))),
    }


def conditional_propagators(
    calc1,
    green_blocks: np.ndarray,
    beta: float,
    coupling_j: float,
) -> tuple[np.ndarray, np.ndarray]:
    length = green_blocks.shape[1]
    sigma_full = calc1.replica_full_matrix(coupling_j**2 * green_blocks**3)
    delta_tau_squared = (beta / length) ** 2
    periodic = calc1.replica_blocks(
        -np.linalg.inv(
            calc1.discrete_derivative(length, periodic=True)
            - delta_tau_squared * sigma_full
        ).T
    )
    antiperiodic = calc1.replica_blocks(
        -np.linalg.inv(
            calc1.discrete_derivative(length, periodic=False)
            - delta_tau_squared * sigma_full
        ).T
    )
    return periodic, antiperiodic


def save_figure(
    path: Path,
    fractions: np.ndarray,
    cases: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]],
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
    figure, axes = plt.subplots(
        len(cases), 4, figsize=(12.8, 3.1 * len(cases)), constrained_layout=True
    )
    if len(cases) == 1:
        axes = axes[None, :]
    for row, (label, path_surface, ed_surface, error) in enumerate(cases):
        difference = ed_surface - path_surface
        scale = max(np.max(np.abs(path_surface)), np.max(np.abs(ed_surface)))
        panels = [
            (path_surface, f"saddle {label}", "RdBu_r", -scale, scale),
            (ed_surface, f"ED {label}", "RdBu_r", -scale, scale),
            (difference, "ED minus saddle", "RdBu_r", -scale, scale),
            (error, "ED jackknife error", "viridis", 0.0, float(error.max())),
        ]
        for axis, (data, title, cmap, vmin, vmax) in zip(axes[row], panels):
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
    figure.savefig(path, dpi=180)
    plt.close(figure)


def run(args: argparse.Namespace) -> dict:
    calc1 = load_calculation_1_module()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.samples < 2:
        raise ValueError("at least two disorder samples are needed for jackknife errors")
    if args.time_points < 2:
        raise ValueError("at least two time points are needed")
    if args.weight % 2:
        raise ValueError("the parent string weight must be even")
    if not 0 < args.weight < args.n_majorana:
        raise ValueError("weight must lie strictly between zero and N")
    basis = calc1.ParityBasis(args.n_majorana)
    _, even_actions, odd_actions = calc1.build_hamiltonian_actions(basis)
    measured_indices = tuple(range(args.weight))
    fractions = np.arange(args.time_points, dtype=float) / args.time_points
    times = args.beta * fractions
    rng = np.random.default_rng(args.seed)

    x_a_samples = np.empty(args.samples, dtype=complex)
    inside_trace_samples = np.empty(
        (args.samples, args.time_points, args.time_points), dtype=complex
    )
    outside_trace_samples = np.empty_like(inside_trace_samples)
    for sample_index in range(args.samples):
        h_even, h_odd = calc1.sample_hamiltonian_blocks(
            rng, basis, even_actions, odd_actions, args.coupling_j
        )
        (
            x_a_samples[sample_index],
            inside_trace_samples[sample_index],
            outside_trace_samples[sample_index],
        ) = same_side_sample(
            calc1,
            h_even,
            h_odd,
            basis,
            measured_indices,
            args.beta,
            times,
        )
        print(f"ED sample {sample_index + 1}/{args.samples}", flush=True)

    inside_ratio, inside_error = ratio_and_jackknife(
        x_a_samples, inside_trace_samples
    )
    outside_ratio, outside_error = ratio_and_jackknife(
        x_a_samples, outside_trace_samples
    )

    relative_weight = args.weight / args.n_majorana
    if args.path_slices % args.time_points:
        raise ValueError(
            "path-slices must be divisible by time-points so the grids align"
        )
    saddle_coarse = calc1.solve_weighted_replica_saddle(
        beta=args.beta,
        relative_weight=relative_weight,
        length=args.path_slices,
        coupling_j=args.coupling_j,
        mixing=args.path_mixing,
        tolerance=args.path_tolerance,
        max_iterations=args.path_max_iterations,
    )
    saddle_fine = calc1.solve_weighted_replica_saddle(
        beta=args.beta,
        relative_weight=relative_weight,
        length=2 * args.path_slices,
        coupling_j=args.coupling_j,
        mixing=args.path_mixing,
        tolerance=args.path_tolerance,
        max_iterations=args.path_max_iterations,
    )
    periodic_coarse, antiperiodic_coarse = conditional_propagators(
        calc1, saddle_coarse.G_blocks, args.beta, args.coupling_j
    )
    periodic_fine, antiperiodic_fine = conditional_propagators(
        calc1, saddle_fine.G_blocks, args.beta, args.coupling_j
    )
    path_indices_coarse = (
        np.arange(args.time_points) * args.path_slices // args.time_points
    )
    path_indices_fine = 2 * path_indices_coarse
    p_11_coarse = periodic_coarse[0][
        np.ix_(path_indices_coarse, path_indices_coarse)
    ]
    ap_11_coarse = antiperiodic_coarse[0][
        np.ix_(path_indices_coarse, path_indices_coarse)
    ]
    p_11_fine = periodic_fine[0][
        np.ix_(path_indices_fine, path_indices_fine)
    ]
    ap_11_fine = antiperiodic_fine[0][
        np.ix_(path_indices_fine, path_indices_fine)
    ]
    # The forward-difference derivative is first-order accurate. The aligned
    # grids therefore give a parameter-free Richardson continuum estimate.
    p_11 = 2.0 * p_11_fine - p_11_coarse
    ap_11 = 2.0 * ap_11_fine - ap_11_coarse

    summary = {
        "n_majorana": args.n_majorana,
        "weight": args.weight,
        "relative_weight": relative_weight,
        "beta_j": args.beta * args.coupling_j,
        "disorder_samples": args.samples,
        "time_points": args.time_points,
        "path_slices_coarse": args.path_slices,
        "path_slices_fine": 2 * args.path_slices,
        "continuum_extrapolation": "2 G_(2L) - G_L",
        "seed": args.seed,
        "path_coarse_iterations": saddle_coarse.iterations,
        "path_coarse_relative_update": saddle_coarse.relative_update,
        "path_coarse_converged": saddle_coarse.converged,
        "path_fine_iterations": saddle_fine.iterations,
        "path_fine_relative_update": saddle_fine.relative_update,
        "path_fine_converged": saddle_fine.converged,
        "denominator_mean_xa_squared": float(np.mean(np.real(x_a_samples**2))),
        "inside_p": comparison_metrics(inside_ratio, p_11, inside_error),
        "inside_p_coarse": comparison_metrics(
            inside_ratio, p_11_coarse, inside_error
        ),
        "inside_p_fine": comparison_metrics(inside_ratio, p_11_fine, inside_error),
        "outside_ap": comparison_metrics(outside_ratio, ap_11, outside_error),
        "outside_ap_coarse": comparison_metrics(
            outside_ratio, ap_11_coarse, outside_error
        ),
        "outside_ap_fine": comparison_metrics(
            outside_ratio, ap_11_fine, outside_error
        ),
    }
    np.savez_compressed(
        args.outdir / "same_side_data.npz",
        fractions=fractions,
        times=times,
        x_a_samples=x_a_samples,
        inside_trace_samples=inside_trace_samples,
        outside_trace_samples=outside_trace_samples,
        inside_ratio=inside_ratio,
        inside_standard_error=inside_error,
        outside_ratio=outside_ratio,
        outside_standard_error=outside_error,
        periodic_G11=p_11,
        antiperiodic_G11=ap_11,
        periodic_G11_coarse=p_11_coarse,
        antiperiodic_G11_coarse=ap_11_coarse,
        periodic_G11_fine=p_11_fine,
        antiperiodic_G11_fine=ap_11_fine,
        saddle_G_blocks_coarse=saddle_coarse.G_blocks,
        saddle_G_blocks_fine=saddle_fine.G_blocks,
        periodic_G_blocks_coarse=periodic_coarse,
        antiperiodic_G_blocks_coarse=antiperiodic_coarse,
        periodic_G_blocks_fine=periodic_fine,
        antiperiodic_G_blocks_fine=antiperiodic_fine,
    )
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    save_figure(
        args.outdir / "same_side_comparison.png",
        fractions,
        [
            (r"$G^{\rm P}_{11}$, $i\in A$", p_11, inside_ratio, inside_error),
            (
                r"$G^{\rm AP}_{11}$, $i\notin A$",
                ap_11,
                outside_ratio,
                outside_error,
            ),
        ],
    )
    print(json.dumps(summary, indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-majorana", type=int, default=20)
    parser.add_argument("--weight", type=int, default=4)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--coupling-j", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=24)
    parser.add_argument("--time-points", type=int, default=9)
    parser.add_argument("--path-slices", type=int, default=90)
    parser.add_argument("--path-mixing", type=float, default=0.01)
    parser.add_argument("--path-tolerance", type=float, default=1e-8)
    parser.add_argument("--path-max-iterations", type=int, default=6000)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--outdir", type=Path, default=Path("outputs"))
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
