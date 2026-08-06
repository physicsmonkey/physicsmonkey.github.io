#!/usr/bin/env python3
"""Shared machinery for normalization-fixed cross-replica comparisons.

The background string A is odd, so Tr(exp(-beta H) mu_A) vanishes.  A probe
collision instead produces the even string B=A symmetric-difference {i}:

    psi_i mu_A = c_i mu_B.

This module evaluates the annealed ratio

    E[Y_Ai(tau) Y_Ai(tau')] / E[X_B^2]

from unnormalized traces and compares it, without fitting, with

    c_i^2 G_12(tau,tau') / G_12(0,0).
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
CALCULATION_1 = HERE / "calculation-1"


def load_calculation_1_module():
    module_path = CALCULATION_1 / "compare_ed_path_integral.py"
    spec = importlib.util.spec_from_file_location(
        "compare_ed_path_integral_moment", module_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def thermal_partition_and_even_string_traces(
    h_even: np.ndarray,
    h_odd: np.ndarray,
    basis,
    even_strings: list[tuple[int, ...]],
    beta: float,
) -> tuple[float, np.ndarray]:
    """Return Z and Tr(exp(-beta H) mu_B) for a list of even B."""
    energies_even, vectors_even = np.linalg.eigh(h_even)
    energies_odd, vectors_odd = np.linalg.eigh(h_odd)
    weights_even = np.exp(-beta * energies_even)
    weights_odd = np.exp(-beta * energies_odd)
    partition = float(weights_even.sum() + weights_odd.sum())
    kernel_even = (vectors_even * weights_even) @ vectors_even.conj().T
    kernel_odd = (vectors_odd * weights_odd) @ vectors_odd.conj().T
    columns = np.arange(len(basis.even))
    traces = np.empty(len(even_strings), dtype=complex)
    for string_index, indices in enumerate(even_strings):
        even_rows, even_phases = basis.action_between(indices, 0)
        odd_rows, odd_phases = basis.action_between(indices, 1)
        traces[string_index] = (
            np.sum(kernel_even[columns, even_rows] * even_phases)
            + np.sum(kernel_odd[columns, odd_rows] * odd_phases)
        )
    return partition, traces


def collision_coefficients(
    basis,
    measured_indices: tuple[int, ...],
) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    """Find c_i in psi_i mu_A = c_i mu_(A symmetric-difference {i})."""
    if len(measured_indices) % 2 != 1:
        raise ValueError("the cross-replica background A must have odd weight")
    measured = set(measured_indices)
    mu_rows, mu_phases = basis.action_between(measured_indices, 0)
    coefficients = np.empty(basis.n_majorana, dtype=complex)
    even_strings = []
    for index in range(basis.n_majorana):
        string = tuple(sorted(measured.symmetric_difference({index})))
        even_strings.append(string)
        psi_rows, psi_phases = basis.action_between((index,), 1)
        combined_rows = psi_rows[mu_rows]
        combined_phases = (
            mu_phases * psi_phases[mu_rows] / math.sqrt(2.0)
        )
        string_rows, string_phases = basis.action_between(string, 0)
        if not np.array_equal(combined_rows, string_rows):
            raise RuntimeError("collision string action has inconsistent rows")
        ratios = combined_phases / string_phases
        if not np.allclose(ratios, ratios[0], atol=2e-14):
            raise RuntimeError("collision coefficient is not state independent")
        coefficients[index] = ratios[0]
    return coefficients, even_strings


def ratio_and_jackknife(
    numerator_samples: np.ndarray,
    denominator_samples: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if len(numerator_samples) < 2:
        raise ValueError("at least two samples are needed for jackknife errors")
    ratio = numerator_samples.mean(axis=0) / denominator_samples.mean()
    leave_one_out = []
    for omitted in range(len(numerator_samples)):
        keep = np.arange(len(numerator_samples)) != omitted
        leave_one_out.append(
            numerator_samples[keep].mean(axis=0)
            / denominator_samples[keep].mean()
        )
    estimates = np.asarray(leave_one_out)
    mean = estimates.mean(axis=0)
    error = np.sqrt(
        (len(estimates) - 1) * np.mean((estimates - mean) ** 2, axis=0)
    )
    return ratio, error


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


def symmetric_sample(
    blocks: np.ndarray,
    indices: np.ndarray,
) -> np.ndarray:
    lattice = blocks[1][np.ix_(indices, indices)]
    return 0.5 * (lattice + lattice.T)


def comparison_metrics(
    ed: np.ndarray,
    prediction: np.ndarray,
    error: np.ndarray,
) -> dict[str, float]:
    difference = ed - prediction
    mask = np.ones(ed.shape, dtype=bool)
    mask[0, 0] = False
    return {
        "relative_residual_all_entries": float(
            np.linalg.norm(difference) / np.linalg.norm(ed)
        ),
        "relative_residual_away_from_collision": float(
            np.linalg.norm(difference[mask]) / np.linalg.norm(ed[mask])
        ),
        "cosine_similarity": float(
            np.sum(ed * prediction)
            / (np.linalg.norm(ed) * np.linalg.norm(prediction))
        ),
        "difference_to_jackknife_error_norm": float(
            np.linalg.norm(difference) / np.linalg.norm(error)
        ),
        "ed_collision_value": float(ed[0, 0]),
        "predicted_collision_value": float(prediction[0, 0]),
        "max_absolute_difference": float(np.max(np.abs(difference))),
    }


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
    for row, (label, prediction, ed, error) in enumerate(cases):
        difference = ed - prediction
        scale = max(np.max(np.abs(prediction)), np.max(np.abs(ed)))
        difference_scale = float(np.max(np.abs(difference)))
        panels = [
            (prediction, f"normalized saddle {label}", "RdBu_r", -scale, scale),
            (ed, f"moment-normalized ED {label}", "RdBu_r", -scale, scale),
            (
                difference,
                "ED minus prediction",
                "RdBu_r",
                -difference_scale,
                difference_scale,
            ),
            (
                error,
                "ED jackknife error",
                "viridis",
                0.0,
                float(np.max(error)),
            ),
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
    source = np.load(args.input)
    zeta_samples = source["zeta_samples"]
    fractions = source["fractions"]
    n_samples, n_majorana, time_points = zeta_samples.shape
    if args.samples is not None:
        n_samples = min(args.samples, n_samples)
        zeta_samples = zeta_samples[:n_samples]
    if args.weight % 2 != 1:
        raise ValueError("A must have odd weight")
    if n_majorana != args.n_majorana:
        raise ValueError("input N does not match --n-majorana")
    if args.path_slices % time_points:
        raise ValueError("path-slices must be divisible by the number of ED times")

    basis = calc1.ParityBasis(n_majorana)
    _, even_actions, odd_actions = calc1.build_hamiltonian_actions(basis)
    measured_indices = tuple(range(args.weight))
    coefficients, even_strings = collision_coefficients(
        basis, measured_indices
    )
    rng = np.random.default_rng(args.seed)
    partitions = np.empty(n_samples)
    collision_traces = np.empty((n_samples, n_majorana), dtype=complex)
    for sample_index in range(n_samples):
        h_even, h_odd = calc1.sample_hamiltonian_blocks(
            rng, basis, even_actions, odd_actions, args.coupling_j
        )
        partitions[sample_index], collision_traces[sample_index] = (
            thermal_partition_and_even_string_traces(
                h_even,
                h_odd,
                basis,
                even_strings,
                args.beta,
            )
        )
        print(f"moment sample {sample_index + 1}/{n_samples}", flush=True)

    raw_probe_traces = zeta_samples * partitions[:, None, None]
    expected_collisions = (
        collision_traces * coefficients[None, :]
    )
    collision_replay_error = float(
        np.linalg.norm(raw_probe_traces[:, :, 0] - expected_collisions)
        / np.linalg.norm(expected_collisions)
    )
    if collision_replay_error > 2e-10:
        raise RuntimeError(
            "saved zeta ensemble does not match the regenerated Hamiltonians; "
            f"relative collision error is {collision_replay_error:.3e}"
        )

    if args.sector == "outside":
        group_indices = np.arange(args.weight, n_majorana)
        sector_label = "AP"
    else:
        group_indices = np.arange(args.weight)
        sector_label = "P"
    groups: dict[str, np.ndarray] = {"sector_average": group_indices}
    if args.probe is not None:
        probe_index = args.probe - 1
        if probe_index not in group_indices:
            raise ValueError("the requested fixed probe is outside the sector")
        groups[f"probe_{args.probe}"] = np.asarray([probe_index])

    group_data: dict[str, dict[str, np.ndarray | float]] = {}
    for label, indices in groups.items():
        selected = raw_probe_traces[:, indices, :]
        numerator_samples = np.einsum(
            "sit,siu->stu", selected, selected, optimize=True
        ).real / len(indices)
        denominator_samples = np.mean(
            np.real(collision_traces[:, indices] ** 2), axis=1
        )
        ratio, error = ratio_and_jackknife(
            numerator_samples, denominator_samples
        )
        collision_squared = float(
            np.real(np.mean(coefficients[indices] ** 2))
        )
        group_data[label] = {
            "numerator_samples": numerator_samples,
            "denominator_samples": denominator_samples,
            "ratio": ratio,
            "error": error,
            "collision_squared": collision_squared,
        }

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary: dict = {
        "n_majorana": n_majorana,
        "weight": args.weight,
        "sector": args.sector,
        "conditional_boundary_condition": sector_label,
        "beta_j": args.beta * args.coupling_j,
        "disorder_samples": n_samples,
        "time_points": time_points,
        "path_slices_coarse": args.path_slices,
        "path_slices_fine": 2 * args.path_slices,
        "continuum_extrapolation": "2 G_(2L) - G_L",
        "seed": args.seed,
        "collision_replay_relative_error": collision_replay_error,
        "groups": {},
        "path_cases": {},
    }
    arrays: dict[str, np.ndarray] = {
        "fractions": fractions,
        "partitions": partitions,
        "collision_coefficients": coefficients,
        "collision_traces": collision_traces,
        "raw_probe_traces": raw_probe_traces,
    }
    for label, data in group_data.items():
        summary["groups"][label] = {
            "indices_one_based": [int(index + 1) for index in groups[label]],
            "collision_coefficient_squared": data["collision_squared"],
            "denominator_mean": float(
                np.mean(data["denominator_samples"])
            ),
        }
        arrays[f"{label}_numerator_samples"] = data["numerator_samples"]
        arrays[f"{label}_denominator_samples"] = data["denominator_samples"]
        arrays[f"{label}_ratio"] = data["ratio"]
        arrays[f"{label}_jackknife_error"] = data["error"]

    indices_coarse = (
        np.arange(time_points) * args.path_slices // time_points
    )
    indices_fine = 2 * indices_coarse
    for relative_weight in args.path_weights:
        saddles = []
        conditionals = []
        for length in (args.path_slices, 2 * args.path_slices):
            saddle = calc1.solve_weighted_replica_saddle(
                beta=args.beta,
                relative_weight=relative_weight,
                length=length,
                coupling_j=args.coupling_j,
                mixing=args.path_mixing,
                tolerance=args.path_tolerance,
                max_iterations=args.path_max_iterations,
            )
            periodic, antiperiodic = conditional_propagators(
                calc1, saddle.G_blocks, args.beta, args.coupling_j
            )
            saddles.append(saddle)
            conditionals.append(
                antiperiodic if args.sector == "outside" else periodic
            )
        coarse = symmetric_sample(conditionals[0], indices_coarse)
        fine = symmetric_sample(conditionals[1], indices_fine)
        continuum = 2.0 * fine - coarse
        key = f"{relative_weight:.3f}".rstrip("0").rstrip(".").replace(".", "p")
        case_summary = {
            "relative_weight": relative_weight,
            "coarse_iterations": saddles[0].iterations,
            "coarse_converged": saddles[0].converged,
            "fine_iterations": saddles[1].iterations,
            "fine_converged": saddles[1].converged,
            "groups": {},
        }
        figure_cases = []
        for label, data in group_data.items():
            prediction = (
                data["collision_squared"] * continuum / continuum[0, 0]
            )
            metrics = comparison_metrics(
                data["ratio"], prediction, data["error"]
            )
            case_summary["groups"][label] = metrics
            arrays[f"prediction_w{key}_{label}"] = prediction
            figure_cases.append(
                (
                    label.replace("_", " "),
                    prediction,
                    data["ratio"],
                    data["error"],
                )
            )
        summary["path_cases"][key] = case_summary
        arrays[f"conditional_G12_coarse_w{key}"] = coarse
        arrays[f"conditional_G12_fine_w{key}"] = fine
        arrays[f"conditional_G12_continuum_w{key}"] = continuum
        save_figure(
            args.outdir / f"moment_normalized_w{key}.png",
            fractions,
            figure_cases,
        )

    np.savez_compressed(
        args.outdir / "moment_normalized_data.npz", **arrays
    )
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))
    return summary


def parser_with_defaults(
    *,
    default_input: Path,
    default_weight: int,
    default_sector: str,
    default_path_weights: list[float],
    default_probe: int | None,
    default_outdir: Path,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--n-majorana", type=int, default=20)
    parser.add_argument("--weight", type=int, default=default_weight)
    parser.add_argument(
        "--sector",
        choices=("inside", "outside"),
        default=default_sector,
    )
    parser.add_argument("--probe", type=int, default=default_probe)
    parser.add_argument(
        "--path-weights",
        type=float,
        nargs="+",
        default=default_path_weights,
    )
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--coupling-j", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=None)
    parser.add_argument("--path-slices", type=int, default=90)
    parser.add_argument("--path-mixing", type=float, default=0.01)
    parser.add_argument("--path-tolerance", type=float, default=1e-8)
    parser.add_argument("--path-max-iterations", type=int, default=6000)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--outdir", type=Path, default=default_outdir)
    return parser
