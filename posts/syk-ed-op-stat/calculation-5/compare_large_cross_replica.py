#!/usr/bin/env python3
"""Normalization-fixed N=24 tests of the conditional G_12 propagators.

The even W=8 string used for the G_11 test cannot be used directly: a
single probe on either replica would leave an odd trace.  We instead use the
two neighboring odd backgrounds

    A_- = {1,...,7},   i outside A_-  (AP),
    A_+ = {1,...,9},   i inside A_+   (P),

whose probe collisions both leave weight-eight strings.  The ED ratios are
compared without a fitted scale to c_i^2 G_12/G_12(0,0).

The implementation reuses each Hamiltonian eigensystem for both backgrounds
and all requested temperatures.  Probe operators are transformed to the
energy basis once per flavor, avoiding the dense thermal-kernel list used by
the original N=20 calculation.
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_modules():
    large = load_module("large_same_side_cross", HERE / "compare_large_same_side.py")
    calc1, calc4 = large.load_calculation_modules()
    moment = load_module(
        "cross_replica_moment_large",
        PARENT / "calculation-2" / "cross_replica_moment_normalization.py",
    )
    return large, calc1, calc4, moment


def cross_action_in_eigenbasis(
    rows: np.ndarray,
    phases: np.ndarray,
    source_vectors: np.ndarray,
    target_vectors: np.ndarray,
) -> np.ndarray:
    """Transform an operator mapping source parity to target parity."""
    acted = np.empty_like(source_vectors)
    acted[rows, :] = phases[:, None] * source_vectors
    return target_vectors.conj().T @ acted


def odd_background_matrices(
    basis,
    measured_indices: tuple[int, ...],
    vectors_even: np.ndarray,
    vectors_odd: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return mu_A in the odd<-even and even<-odd energy bases."""
    rows, phases = basis.action_between(measured_indices, 0)
    odd_even = cross_action_in_eigenbasis(
        rows, phases, vectors_even, vectors_odd
    )
    # mu_A is Hermitian in the convention used throughout the calculation.
    return odd_even, odd_even.conj().T


def probe_matrix(
    basis,
    index: int,
    vectors_even: np.ndarray,
    vectors_odd: np.ndarray,
) -> np.ndarray:
    """Return psi_i in the even<-odd energy basis."""
    rows, phases = basis.action_between((index,), 1)
    return cross_action_in_eigenbasis(
        rows, phases / math.sqrt(2.0), vectors_odd, vectors_even
    )


def traces_from_eigensystem(
    energies_even: np.ndarray,
    vectors_even: np.ndarray,
    energies_odd: np.ndarray,
    vectors_odd: np.ndarray,
    basis,
    backgrounds: dict[str, tuple[int, ...]],
    betas: np.ndarray,
    fractions: np.ndarray,
) -> dict[str, np.ndarray]:
    """Evaluate unnormalized Y_Ai(tau) for all backgrounds, i, and beta."""
    background_matrices = {
        label: odd_background_matrices(
            basis, indices, vectors_even, vectors_odd
        )
        for label, indices in backgrounds.items()
    }
    result = {
        label: np.empty(
            (len(betas), basis.n_majorana, len(fractions)),
            dtype=np.complex128,
        )
        for label in backgrounds
    }
    thermal_weights = []
    for beta in betas:
        times = beta * fractions
        thermal_weights.append(
            (
                np.exp(-(beta - times[:, None]) * energies_even[None, :]),
                np.exp(-times[:, None] * energies_even[None, :]),
                np.exp(-(beta - times[:, None]) * energies_odd[None, :]),
                np.exp(-times[:, None] * energies_odd[None, :]),
            )
        )

    for index in range(basis.n_majorana):
        even_odd = probe_matrix(
            basis, index, vectors_even, vectors_odd
        )
        odd_even = even_odd.conj().T
        for label, (mu_odd_even, mu_even_odd) in background_matrices.items():
            contraction_even = even_odd * mu_odd_even.T
            contraction_odd = odd_even * mu_even_odd.T
            for beta_index, (left_e, right_e, left_o, right_o) in enumerate(
                thermal_weights
            ):
                # einsum evaluates l(t)^T Q r(t) for every requested time.
                result[label][beta_index, index] = (
                    np.einsum(
                        "ta,ab,tb->t",
                        left_e,
                        contraction_even,
                        right_o,
                        optimize=True,
                    )
                    + np.einsum(
                        "ta,ab,tb->t",
                        left_o,
                        contraction_odd,
                        right_e,
                        optimize=True,
                    )
                )
        del even_odd, odd_even
    return result


def ratio_and_jackknife(
    numerator_samples: np.ndarray,
    denominator_samples: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ratio = numerator_samples.mean(axis=0) / denominator_samples.mean()
    if len(numerator_samples) < 2:
        return ratio, np.full(ratio.shape, np.nan)
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


def ed_ratio(
    trace_samples: np.ndarray,
    group_indices: np.ndarray,
    coefficients: np.ndarray,
) -> dict[str, np.ndarray | float]:
    selected = trace_samples[:, group_indices, :]
    numerator_samples = np.einsum(
        "sit,siu->stu", selected, selected, optimize=True
    ).real / len(group_indices)
    collision_traces = selected[:, :, 0] / coefficients[group_indices][None, :]
    denominator_samples = np.mean(
        np.real(collision_traces**2), axis=1
    )
    ratio, error = ratio_and_jackknife(
        numerator_samples, denominator_samples
    )
    collision_squared = float(np.real(np.mean(coefficients[group_indices] ** 2)))
    denominator_ess = float(
        denominator_samples.sum() ** 2 / np.sum(denominator_samples**2)
    )
    return {
        "numerator_samples": numerator_samples,
        "denominator_samples": denominator_samples,
        "ratio": ratio,
        "error": error,
        "collision_squared": collision_squared,
        "denominator_effective_sample_size": denominator_ess,
    }


def validate_evaluator(large, calc1) -> float:
    """Compare the energy-basis evaluator with calculation 1 at small N."""
    n_majorana = 8
    beta = 0.7
    fractions = np.arange(5, dtype=float) / 5
    basis = calc1.ParityBasis(n_majorana)
    _, even_actions, odd_actions = calc1.build_hamiltonian_actions(basis)
    h_even, h_odd = calc1.sample_hamiltonian_blocks(
        np.random.default_rng(77123), basis, even_actions, odd_actions, 1.0
    )
    energies_even, vectors_even = np.linalg.eigh(h_even)
    energies_odd, vectors_odd = np.linalg.eigh(h_odd)
    background = tuple(range(3))
    candidate = traces_from_eigensystem(
        energies_even,
        vectors_even,
        energies_odd,
        vectors_odd,
        basis,
        {"test": background},
        np.asarray([beta]),
        fractions,
    )["test"][0]
    reference_normalized = calc1.zeta_grid(
        h_even,
        h_odd,
        basis,
        background,
        beta,
        beta * fractions,
    )
    partition = float(
        np.exp(-beta * energies_even).sum()
        + np.exp(-beta * energies_odd).sum()
    )
    reference = reference_normalized * partition
    return float(np.max(np.abs(candidate - reference)))


def format_number(value: float) -> str:
    return f"{value:.6g}".replace("-", "m").replace(".", "p")


def run(args: argparse.Namespace) -> dict:
    large, calc1, calc4, moment = load_modules()
    if args.validate:
        validation_error = validate_evaluator(large, calc1)
        print(
            f"cross-replica evaluator validation max error {validation_error:.3e}",
            flush=True,
        )
        if validation_error > 2e-12:
            raise RuntimeError("cross-replica evaluator validation failed")
    if args.outside_weight % 2 != 1 or args.inside_weight % 2 != 1:
        raise ValueError("both cross-replica backgrounds must have odd weight")
    if not args.outside_weight < args.target_weight < args.inside_weight:
        raise ValueError("weights must obey outside < target < inside")
    if args.path_slices % args.time_points:
        raise ValueError("path-slices must be divisible by time-points")

    args.outdir.mkdir(parents=True, exist_ok=True)
    betas = np.asarray(args.betas, dtype=float)
    fractions = np.arange(args.time_points, dtype=float) / args.time_points
    basis = calc1.ParityBasis(args.n_majorana)
    actions, action_source = large.compact_actions(
        calc1, basis, args.action_cache
    )
    backgrounds = {
        "outside_ap": tuple(range(args.outside_weight)),
        "inside_p": tuple(range(args.inside_weight)),
    }
    coefficients = {}
    for label, background in backgrounds.items():
        coefficients[label], _ = moment.collision_coefficients(
            basis, background
        )

    traces = {
        label: np.empty(
            (
                len(betas),
                args.samples,
                args.n_majorana,
                args.time_points,
            ),
            dtype=np.complex128,
        )
        for label in backgrounds
    }
    rng = np.random.default_rng(args.seed)
    sample_seconds = []
    ed_started = time.perf_counter()
    for sample_index in range(args.samples):
        started = time.perf_counter()
        h_even, h_odd = large.sample_hamiltonian_blocks_compact(
            rng, basis, actions, args.coupling_j
        )
        energies_even, vectors_even = np.linalg.eigh(h_even)
        energies_odd, vectors_odd = np.linalg.eigh(h_odd)
        del h_even, h_odd
        sample_traces = traces_from_eigensystem(
            energies_even,
            vectors_even,
            energies_odd,
            vectors_odd,
            basis,
            backgrounds,
            betas,
            fractions,
        )
        for label in backgrounds:
            traces[label][:, sample_index] = sample_traces[label]
        del energies_even, energies_odd, vectors_even, vectors_odd, sample_traces
        gc.collect()
        sample_seconds.append(time.perf_counter() - started)
        print(
            f"ED sample {sample_index + 1}/{args.samples}: "
            f"{sample_seconds[-1]:.2f} s, peak RSS {large.peak_rss_mib():.0f} MiB",
            flush=True,
        )

    groups = {
        "outside_ap": np.arange(args.outside_weight, args.n_majorana),
        "inside_p": np.arange(args.inside_weight),
    }
    ed_data: dict[tuple[int, str], dict[str, np.ndarray | float]] = {}
    for beta_index in range(len(betas)):
        for label in backgrounds:
            ed_data[(beta_index, label)] = ed_ratio(
                traces[label][beta_index],
                groups[label],
                coefficients[label],
            )

    unique_weights = sorted(
        {
            args.outside_weight / args.n_majorana,
            args.target_weight / args.n_majorana,
            args.inside_weight / args.n_majorana,
        }
    )
    summary: dict = {
        "n_majorana": args.n_majorana,
        "outside_weight": args.outside_weight,
        "target_weight": args.target_weight,
        "inside_weight": args.inside_weight,
        "disorder_samples": args.samples,
        "time_points": args.time_points,
        "betas": betas.tolist(),
        "coupling_j": args.coupling_j,
        "path_slices_coarse": args.path_slices,
        "path_slices_fine": 2 * args.path_slices,
        "continuum_extrapolation": "2 G_(2L) - G_L",
        "seed": args.seed,
        "compact_action_source": action_source,
        "compact_action_mib": actions.nbytes / 2**20,
        "ed_total_seconds": time.perf_counter() - ed_started,
        "mean_seconds_per_sample": float(np.mean(sample_seconds)),
        "sample_seconds": sample_seconds,
        "peak_rss_mib": large.peak_rss_mib(),
        "temperatures": {},
    }
    arrays: dict[str, np.ndarray] = {
        "betas": betas,
        "fractions": fractions,
        "outside_collision_coefficients": coefficients["outside_ap"],
        "inside_collision_coefficients": coefficients["inside_p"],
        "outside_probe_traces": traces["outside_ap"],
        "inside_probe_traces": traces["inside_p"],
    }

    indices_coarse = np.arange(args.time_points) * args.path_slices // args.time_points
    indices_fine = 2 * indices_coarse
    for beta_index, beta in enumerate(betas):
        beta_key = format_number(float(beta * args.coupling_j))
        temperature_summary = {"beta_j": float(beta * args.coupling_j), "sectors": {}}
        conditional_by_weight = {}
        for relative_weight in unique_weights:
            saddles = []
            conditionals = []
            for length in (args.path_slices, 2 * args.path_slices):
                saddle = calc1.solve_weighted_replica_saddle(
                    beta=float(beta),
                    relative_weight=relative_weight,
                    length=length,
                    coupling_j=args.coupling_j,
                    mixing=args.path_mixing,
                    tolerance=args.path_tolerance,
                    max_iterations=args.path_max_iterations,
                )
                periodic, antiperiodic = calc4.conditional_propagators(
                    calc1, saddle.G_blocks, float(beta), args.coupling_j
                )
                saddles.append(saddle)
                conditionals.append((periodic, antiperiodic))
            conditional_by_weight[relative_weight] = (saddles, conditionals)

        figure_cases = []
        for label in ("inside_p", "outside_ap"):
            literal_weight = (
                args.inside_weight if label == "inside_p" else args.outside_weight
            ) / args.n_majorana
            tested_weights = [literal_weight]
            target_relative_weight = args.target_weight / args.n_majorana
            if target_relative_weight != literal_weight:
                tested_weights.append(target_relative_weight)
            data = ed_data[(beta_index, label)]
            sector_summary = {
                "background_weight": len(backgrounds[label]),
                "probe_count": len(groups[label]),
                "probe_indices_one_based": [int(index + 1) for index in groups[label]],
                "collision_coefficient_squared": data["collision_squared"],
                "denominator_effective_sample_size": data[
                    "denominator_effective_sample_size"
                ],
                "path_cases": {},
            }
            arrays[f"ratio_beta{beta_key}_{label}"] = data["ratio"]
            arrays[f"jackknife_error_beta{beta_key}_{label}"] = data["error"]
            arrays[f"numerator_samples_beta{beta_key}_{label}"] = data[
                "numerator_samples"
            ]
            arrays[f"denominator_samples_beta{beta_key}_{label}"] = data[
                "denominator_samples"
            ]
            for relative_weight in tested_weights:
                saddles, conditionals = conditional_by_weight[relative_weight]
                boundary_index = 0 if label == "inside_p" else 1
                coarse = moment.symmetric_sample(
                    conditionals[0][boundary_index], indices_coarse
                )
                fine = moment.symmetric_sample(
                    conditionals[1][boundary_index], indices_fine
                )
                continuum = 2.0 * fine - coarse
                prediction = data["collision_squared"] * continuum / continuum[0, 0]
                metrics = moment.comparison_metrics(
                    data["ratio"], prediction, data["error"]
                )
                weight_key = format_number(relative_weight)
                sector_summary["path_cases"][weight_key] = {
                    "relative_weight": relative_weight,
                    "coarse_iterations": saddles[0].iterations,
                    "coarse_converged": saddles[0].converged,
                    "fine_iterations": saddles[1].iterations,
                    "fine_converged": saddles[1].converged,
                    **metrics,
                }
                arrays[
                    f"conditional_G12_beta{beta_key}_{label}_w{weight_key}"
                ] = continuum
                arrays[f"prediction_beta{beta_key}_{label}_w{weight_key}"] = prediction
                if relative_weight == target_relative_weight:
                    figure_cases.append(
                        (
                            r"$G^{\rm P}_{12}$, $W=9$" if label == "inside_p" else r"$G^{\rm AP}_{12}$, $W=7$",
                            prediction,
                            data["ratio"],
                            data["error"],
                        )
                    )
            temperature_summary["sectors"][label] = sector_summary
        summary["temperatures"][beta_key] = temperature_summary
        moment.save_figure(
            args.outdir / f"cross_replica_comparison_beta{beta_key}.png",
            fractions,
            figure_cases,
        )

    np.savez_compressed(args.outdir / "large_cross_replica_data.npz", **arrays)
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-majorana", type=int, default=24)
    parser.add_argument("--outside-weight", type=int, default=7)
    parser.add_argument("--target-weight", type=int, default=8)
    parser.add_argument("--inside-weight", type=int, default=9)
    parser.add_argument("--betas", type=float, nargs="+", default=[0.5, 1.0])
    parser.add_argument("--coupling-j", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--time-points", type=int, default=9)
    parser.add_argument("--path-slices", type=int, default=90)
    parser.add_argument("--path-mixing", type=float, default=0.01)
    parser.add_argument("--path-tolerance", type=float, default=1e-8)
    parser.add_argument("--path-max-iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument(
        "--action-cache",
        type=Path,
        default=Path("/tmp/syk_ed_op_stat_compact_actions_N24.npz"),
    )
    parser.add_argument("--validate", action="store_true")
    parser.add_argument(
        "--outdir", type=Path, default=HERE / "cross_replica_outputs_4"
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
