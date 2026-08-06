#!/usr/bin/env python3
"""Large-N pilot for normalization-fixed same-side correlators.

This extends calculation 4 to N=24 and W=8.  Two changes keep the exact
calculation practical:

1. q=4 Hamiltonian actions are stored as uint16 permutations and int8 phase
   codes rather than Python lists of int64/complex128 arrays.
2. Same-side trace intermediates are consumed one time separation at a time
   rather than retained for the entire time grid.
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import itertools
import json
import math
import os
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CALCULATION_1 = HERE.parent / "calculation-1"
CALCULATION_4 = HERE.parent / "calculation-4"
PHASE_VALUES = np.asarray([1.0, -1.0, 1.0j, -1.0j], dtype=np.complex128)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_calculation_modules():
    calc1 = load_module(
        "compare_ed_path_integral_large",
        CALCULATION_1 / "compare_ed_path_integral.py",
    )
    calc4 = load_module(
        "compare_same_side_large",
        CALCULATION_4 / "compare_same_side.py",
    )
    return calc1, calc4


def peak_rss_mib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def encode_phases(phases: np.ndarray) -> np.ndarray:
    codes = np.full(phases.shape, -1, dtype=np.int8)
    for code, value in enumerate(PHASE_VALUES):
        codes[np.isclose(phases, value, atol=1e-13)] = code
    if np.any(codes < 0):
        raise RuntimeError("Majorana action contained a noncanonical phase")
    return codes


@dataclass
class CompactActions:
    n_majorana: int
    rows_even: np.ndarray
    rows_odd: np.ndarray
    phases_even: np.ndarray
    phases_odd: np.ndarray

    @property
    def n_terms(self) -> int:
        return self.rows_even.shape[0]

    @property
    def nbytes(self) -> int:
        return sum(
            array.nbytes
            for array in (
                self.rows_even,
                self.rows_odd,
                self.phases_even,
                self.phases_odd,
            )
        )


def build_compact_actions(calc1, basis) -> CompactActions:
    block_dim = len(basis.even)
    if block_dim > np.iinfo(np.uint16).max:
        raise ValueError("uint16 action rows require parity blocks below 65536")
    n_terms = math.comb(basis.n_majorana, 4)
    rows_even = np.empty((n_terms, block_dim), dtype=np.uint16)
    rows_odd = np.empty_like(rows_even)
    phases_even = np.empty((n_terms, block_dim), dtype=np.int8)
    phases_odd = np.empty_like(phases_even)
    for term_index, indices in enumerate(
        itertools.combinations(range(basis.n_majorana), 4)
    ):
        even_rows, even_phases = basis.action_between(indices, 0)
        odd_rows, odd_phases = basis.action_between(indices, 1)
        rows_even[term_index] = even_rows
        rows_odd[term_index] = odd_rows
        phases_even[term_index] = encode_phases(even_phases)
        phases_odd[term_index] = encode_phases(odd_phases)
    return CompactActions(
        basis.n_majorana,
        rows_even,
        rows_odd,
        phases_even,
        phases_odd,
    )


def save_action_cache(path: Path, actions: CompactActions) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        n_majorana=np.asarray(actions.n_majorana),
        rows_even=actions.rows_even,
        rows_odd=actions.rows_odd,
        phases_even=actions.phases_even,
        phases_odd=actions.phases_odd,
    )


def load_action_cache(path: Path, n_majorana: int) -> CompactActions:
    source = np.load(path)
    saved_n = int(source["n_majorana"])
    if saved_n != n_majorana:
        raise ValueError(
            f"action cache has N={saved_n}, requested N={n_majorana}"
        )
    return CompactActions(
        saved_n,
        source["rows_even"],
        source["rows_odd"],
        source["phases_even"],
        source["phases_odd"],
    )


def compact_actions(calc1, basis, cache: Path | None) -> tuple[CompactActions, str]:
    if cache is not None and cache.exists():
        return load_action_cache(cache, basis.n_majorana), "loaded"
    actions = build_compact_actions(calc1, basis)
    if cache is not None:
        save_action_cache(cache, actions)
    return actions, "built"


def sample_hamiltonian_blocks_compact(
    rng: np.random.Generator,
    basis,
    actions: CompactActions,
    coupling_j: float,
) -> tuple[np.ndarray, np.ndarray]:
    coefficient_std = (
        coupling_j
        * math.sqrt(math.factorial(3) / basis.n_majorana**3)
        / 4.0
    )
    coefficients = rng.normal(0.0, coefficient_std, actions.n_terms)
    block_dim = len(basis.even)
    columns = np.arange(block_dim)
    h_even = np.zeros((block_dim, block_dim), dtype=np.complex128)
    h_odd = np.zeros_like(h_even)
    for term_index, coefficient in enumerate(coefficients):
        h_even[actions.rows_even[term_index], columns] += (
            coefficient * PHASE_VALUES[actions.phases_even[term_index]]
        )
        h_odd[actions.rows_odd[term_index], columns] += (
            coefficient * PHASE_VALUES[actions.phases_odd[term_index]]
        )
    h_even = 0.5 * (h_even + h_even.conj().T)
    h_odd = 0.5 * (h_odd + h_odd.conj().T)
    return h_even, h_odd


def action_in_eigenbasis(
    rows: np.ndarray,
    phases: np.ndarray,
    vectors: np.ndarray,
) -> np.ndarray:
    acted = np.empty_like(vectors)
    acted[rows, :] = phases[:, None] * vectors
    return vectors.conj().T @ acted


def kernel_from_eigensystem(
    energies: np.ndarray,
    vectors: np.ndarray,
    interval: float,
) -> np.ndarray:
    return (vectors * np.exp(-interval * energies)) @ vectors.conj().T


def conjugated_kernel_average(
    kernel_opposite: np.ndarray,
    actions_from_parity: list[tuple[np.ndarray, np.ndarray]],
    actions_from_opposite: list[tuple[np.ndarray, np.ndarray]],
    indices: range,
) -> np.ndarray:
    result = np.zeros_like(kernel_opposite)
    for index in indices:
        opposite_rows, phases_from_parity = actions_from_parity[index]
        parity_rows, phases_from_opposite = actions_from_opposite[index]
        temp = kernel_opposite[:, opposite_rows].copy()
        temp *= phases_from_parity[None, :] / math.sqrt(2.0)
        temp *= phases_from_opposite[:, None] / math.sqrt(2.0)
        result[parity_rows, :] += temp
        del temp
    result /= len(indices)
    return result


def trace_from_eigenbasis(
    energies: np.ndarray,
    q_eigen: np.ndarray,
    mu_eigen: np.ndarray,
    beta_minus_tau: float,
    tau_prime: float,
) -> complex:
    left = np.exp(-beta_minus_tau * energies)
    right = np.exp(-tau_prime * energies)
    return np.sum(left[:, None] * q_eigen * right[None, :] * mu_eigen.T)


def transformed_conjugated_kernel(
    kernel_opposite: np.ndarray,
    actions_from_parity: list[tuple[np.ndarray, np.ndarray]],
    actions_from_opposite: list[tuple[np.ndarray, np.ndarray]],
    indices: range,
    vectors: np.ndarray,
) -> np.ndarray:
    q_computational = conjugated_kernel_average(
        kernel_opposite,
        actions_from_parity,
        actions_from_opposite,
        indices,
    )
    temp = q_computational @ vectors
    del q_computational
    result = vectors.conj().T @ temp
    del temp
    return result


def same_side_sample_streaming(
    h_even: np.ndarray,
    h_odd: np.ndarray,
    basis,
    measured_indices: tuple[int, ...],
    beta: float,
    times: np.ndarray,
) -> tuple[complex, np.ndarray, np.ndarray, dict[str, float]]:
    timings = {
        "diagonalization_seconds": 0.0,
        "background_transform_seconds": 0.0,
        "kernel_seconds": 0.0,
        "probe_transform_seconds": 0.0,
        "trace_contraction_seconds": 0.0,
    }
    started = time.perf_counter()
    energies_even, vectors_even = np.linalg.eigh(h_even)
    energies_odd, vectors_odd = np.linalg.eigh(h_odd)
    timings["diagonalization_seconds"] = time.perf_counter() - started

    started = time.perf_counter()
    mu_even_rows, mu_even_phases = basis.action_between(measured_indices, 0)
    mu_odd_rows, mu_odd_phases = basis.action_between(measured_indices, 1)
    mu_even_eigen = action_in_eigenbasis(
        mu_even_rows, mu_even_phases, vectors_even
    )
    mu_odd_eigen = action_in_eigenbasis(
        mu_odd_rows, mu_odd_phases, vectors_odd
    )
    x_a = (
        np.sum(np.exp(-beta * energies_even) * np.diag(mu_even_eigen))
        + np.sum(np.exp(-beta * energies_odd) * np.diag(mu_odd_eigen))
    )
    timings["background_transform_seconds"] = time.perf_counter() - started

    actions_even = [
        basis.action_between((index,), 0)
        for index in range(basis.n_majorana)
    ]
    actions_odd = [
        basis.action_between((index,), 1)
        for index in range(basis.n_majorana)
    ]
    groups = (
        ("inside", range(len(measured_indices))),
        (
            "outside",
            range(len(measured_indices), basis.n_majorana),
        ),
    )
    surfaces = {
        "inside": np.zeros((len(times), len(times)), dtype=np.complex128),
        "outside": np.zeros((len(times), len(times)), dtype=np.complex128),
    }
    for surface in surfaces.values():
        np.fill_diagonal(surface, -0.5 * x_a)

    spacings = np.diff(times)
    if len(spacings) and not np.allclose(spacings, spacings[0]):
        raise ValueError("the streaming evaluator requires a uniform time grid")
    for delta_index in range(1, len(times)):
        delta = times[delta_index] - times[0]
        started = time.perf_counter()
        kernel_even = kernel_from_eigensystem(
            energies_even, vectors_even, delta
        )
        kernel_odd = kernel_from_eigensystem(
            energies_odd, vectors_odd, delta
        )
        timings["kernel_seconds"] += time.perf_counter() - started

        pairs = [
            (time_index, time_index - delta_index)
            for time_index in range(delta_index, len(times))
        ]
        for label, indices in groups:
            started = time.perf_counter()
            q_even_eigen = transformed_conjugated_kernel(
                kernel_odd,
                actions_even,
                actions_odd,
                indices,
                vectors_even,
            )
            timings["probe_transform_seconds"] += (
                time.perf_counter() - started
            )
            started = time.perf_counter()
            values = []
            for time_index, prime_index in pairs:
                values.append(
                    trace_from_eigenbasis(
                        energies_even,
                        q_even_eigen,
                        mu_even_eigen,
                        beta - times[time_index],
                        times[prime_index],
                    )
                )
            timings["trace_contraction_seconds"] += (
                time.perf_counter() - started
            )
            del q_even_eigen

            started = time.perf_counter()
            q_odd_eigen = transformed_conjugated_kernel(
                kernel_even,
                actions_odd,
                actions_even,
                indices,
                vectors_odd,
            )
            timings["probe_transform_seconds"] += (
                time.perf_counter() - started
            )
            started = time.perf_counter()
            for pair_index, (time_index, prime_index) in enumerate(pairs):
                values[pair_index] += trace_from_eigenbasis(
                    energies_odd,
                    q_odd_eigen,
                    mu_odd_eigen,
                    beta - times[time_index],
                    times[prime_index],
                )
                surfaces[label][time_index, prime_index] = values[pair_index]
                surfaces[label][prime_index, time_index] = -values[pair_index]
            timings["trace_contraction_seconds"] += (
                time.perf_counter() - started
            )
            del q_odd_eigen
        del kernel_even, kernel_odd
        gc.collect()
    return x_a, surfaces["inside"], surfaces["outside"], timings


def ratio_and_jackknife(
    x_a_samples: np.ndarray,
    trace_samples: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    denominator_samples = np.real(x_a_samples**2)
    numerator_samples = np.real(
        trace_samples * x_a_samples[:, None, None]
    )
    ratio = numerator_samples.mean(axis=0) / denominator_samples.mean()
    if len(x_a_samples) < 2:
        return ratio, np.full(ratio.shape, np.nan)
    leave_one_out = []
    for omitted in range(len(x_a_samples)):
        keep = np.arange(len(x_a_samples)) != omitted
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


def comparison_metrics(
    ed: np.ndarray,
    path: np.ndarray,
    error: np.ndarray,
) -> dict[str, float]:
    mask = ~np.eye(ed.shape[0], dtype=bool)
    difference = ed - path
    metrics = {
        "relative_residual_all_entries": float(
            np.linalg.norm(difference) / np.linalg.norm(ed)
        ),
        "relative_residual_off_diagonal": float(
            np.linalg.norm(difference[mask]) / np.linalg.norm(ed[mask])
        ),
        "off_diagonal_cosine_similarity": float(
            np.sum(ed[mask] * path[mask])
            / (np.linalg.norm(ed[mask]) * np.linalg.norm(path[mask]))
        ),
        "ed_contact_mean": float(np.diag(ed).mean()),
        "path_contact_mean": float(np.diag(path).mean()),
        "max_absolute_difference": float(np.max(np.abs(difference))),
    }
    if np.all(np.isfinite(error)) and np.linalg.norm(error):
        metrics["difference_to_jackknife_error_norm"] = float(
            np.linalg.norm(difference) / np.linalg.norm(error)
        )
    return metrics


def validate_compact_hamiltonian(calc1, n_majorana: int = 10) -> float:
    basis = calc1.ParityBasis(n_majorana)
    _, even_actions, odd_actions = calc1.build_hamiltonian_actions(basis)
    compact = build_compact_actions(calc1, basis)
    rng_reference = np.random.default_rng(918273)
    rng_compact = np.random.default_rng(918273)
    reference = calc1.sample_hamiltonian_blocks(
        rng_reference, basis, even_actions, odd_actions, 1.0
    )
    candidate = sample_hamiltonian_blocks_compact(
        rng_compact, basis, compact, 1.0
    )
    return max(
        float(np.max(np.abs(left - right)))
        for left, right in zip(reference, candidate)
    )


def run(args: argparse.Namespace) -> dict:
    calc1, calc4 = load_calculation_modules()
    if args.validate_compact:
        validation_error = validate_compact_hamiltonian(calc1)
        print(
            f"compact Hamiltonian validation max error {validation_error:.3e}",
            flush=True,
        )
        if validation_error > 2e-14:
            raise RuntimeError("compact Hamiltonian validation failed")
    if args.weight % 2:
        raise ValueError("the parent string weight must be even")
    if not 0 < args.weight < args.n_majorana:
        raise ValueError("weight must lie strictly between zero and N")
    if args.time_points < 2:
        raise ValueError("at least two time points are required")

    args.outdir.mkdir(parents=True, exist_ok=True)
    basis = calc1.ParityBasis(args.n_majorana)
    started = time.perf_counter()
    actions, action_source = compact_actions(
        calc1, basis, args.action_cache
    )
    action_seconds = time.perf_counter() - started
    print(
        f"compact actions {action_source}: {actions.n_terms} terms, "
        f"{actions.nbytes / 2**20:.1f} MiB, {action_seconds:.2f} s",
        flush=True,
    )

    measured_indices = tuple(range(args.weight))
    fractions = np.arange(args.time_points, dtype=float) / args.time_points
    times = args.beta * fractions
    rng = np.random.default_rng(args.seed)
    x_a_samples = np.empty(args.samples, dtype=np.complex128)
    inside_trace_samples = np.empty(
        (args.samples, args.time_points, args.time_points),
        dtype=np.complex128,
    )
    outside_trace_samples = np.empty_like(inside_trace_samples)
    sample_timings = []
    pilot_started = time.perf_counter()
    for sample_index in range(args.samples):
        sample_started = time.perf_counter()
        h_started = time.perf_counter()
        h_even, h_odd = sample_hamiltonian_blocks_compact(
            rng, basis, actions, args.coupling_j
        )
        h_seconds = time.perf_counter() - h_started
        (
            x_a_samples[sample_index],
            inside_trace_samples[sample_index],
            outside_trace_samples[sample_index],
            stages,
        ) = same_side_sample_streaming(
            h_even,
            h_odd,
            basis,
            measured_indices,
            args.beta,
            times,
        )
        del h_even, h_odd
        gc.collect()
        stages["hamiltonian_seconds"] = h_seconds
        stages["total_seconds"] = time.perf_counter() - sample_started
        stages["peak_rss_mib"] = peak_rss_mib()
        sample_timings.append(stages)
        print(
            f"ED sample {sample_index + 1}/{args.samples}: "
            f"{stages['total_seconds']:.2f} s, "
            f"peak RSS {stages['peak_rss_mib']:.0f} MiB",
            flush=True,
        )
    ed_seconds = time.perf_counter() - pilot_started

    inside_ratio, inside_error = ratio_and_jackknife(
        x_a_samples, inside_trace_samples
    )
    outside_ratio, outside_error = ratio_and_jackknife(
        x_a_samples, outside_trace_samples
    )
    denominator_samples = np.real(x_a_samples**2)
    denominator_ess = float(
        denominator_samples.sum() ** 2
        / np.sum(denominator_samples**2)
    )
    summary: dict = {
        "n_majorana": args.n_majorana,
        "weight": args.weight,
        "relative_weight": args.weight / args.n_majorana,
        "beta_j": args.beta * args.coupling_j,
        "disorder_samples": args.samples,
        "time_points": args.time_points,
        "seed": args.seed,
        "compact_action_source": action_source,
        "compact_action_mib": actions.nbytes / 2**20,
        "compact_action_seconds": action_seconds,
        "ed_total_seconds": ed_seconds,
        "mean_seconds_per_sample": float(
            np.mean([item["total_seconds"] for item in sample_timings])
        ),
        "sample_timings": sample_timings,
        "peak_rss_mib": peak_rss_mib(),
        "denominator_mean_xa_squared": float(denominator_samples.mean()),
        "denominator_effective_sample_size": denominator_ess,
    }
    arrays: dict[str, np.ndarray] = {
        "fractions": fractions,
        "times": times,
        "x_a_samples": x_a_samples,
        "inside_trace_samples": inside_trace_samples,
        "outside_trace_samples": outside_trace_samples,
        "inside_ratio": inside_ratio,
        "inside_standard_error": inside_error,
        "outside_ratio": outside_ratio,
        "outside_standard_error": outside_error,
    }

    if not args.skip_saddle:
        if args.path_slices % args.time_points:
            raise ValueError(
                "path-slices must be divisible by time-points"
            )
        relative_weight = args.weight / args.n_majorana
        saddle_started = time.perf_counter()
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
        periodic_coarse, antiperiodic_coarse = calc4.conditional_propagators(
            calc1, saddle_coarse.G_blocks, args.beta, args.coupling_j
        )
        periodic_fine, antiperiodic_fine = calc4.conditional_propagators(
            calc1, saddle_fine.G_blocks, args.beta, args.coupling_j
        )
        indices_coarse = (
            np.arange(args.time_points)
            * args.path_slices
            // args.time_points
        )
        indices_fine = 2 * indices_coarse
        p_coarse = periodic_coarse[0][
            np.ix_(indices_coarse, indices_coarse)
        ]
        ap_coarse = antiperiodic_coarse[0][
            np.ix_(indices_coarse, indices_coarse)
        ]
        p_fine = periodic_fine[0][np.ix_(indices_fine, indices_fine)]
        ap_fine = antiperiodic_fine[0][np.ix_(indices_fine, indices_fine)]
        p_continuum = 2.0 * p_fine - p_coarse
        ap_continuum = 2.0 * ap_fine - ap_coarse
        summary.update(
            {
                "path_slices_coarse": args.path_slices,
                "path_slices_fine": 2 * args.path_slices,
                "saddle_seconds": time.perf_counter() - saddle_started,
                "path_coarse_converged": saddle_coarse.converged,
                "path_fine_converged": saddle_fine.converged,
                "inside_p": comparison_metrics(
                    inside_ratio, p_continuum, inside_error
                ),
                "outside_ap": comparison_metrics(
                    outside_ratio, ap_continuum, outside_error
                ),
            }
        )
        arrays.update(
            {
                "periodic_G11": p_continuum,
                "antiperiodic_G11": ap_continuum,
                "periodic_G11_coarse": p_coarse,
                "antiperiodic_G11_coarse": ap_coarse,
                "periodic_G11_fine": p_fine,
                "antiperiodic_G11_fine": ap_fine,
                "saddle_G_blocks_coarse": saddle_coarse.G_blocks,
                "saddle_G_blocks_fine": saddle_fine.G_blocks,
            }
        )
        calc4.save_figure(
            args.outdir / "large_same_side_comparison.png",
            fractions,
            [
                (
                    r"$G^{\rm P}_{11}$, $i\in A$",
                    p_continuum,
                    inside_ratio,
                    inside_error,
                ),
                (
                    r"$G^{\rm AP}_{11}$, $i\notin A$",
                    ap_continuum,
                    outside_ratio,
                    outside_error,
                ),
            ],
        )

    np.savez_compressed(
        args.outdir / "large_same_side_data.npz", **arrays
    )
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-majorana", type=int, default=24)
    parser.add_argument("--weight", type=int, default=8)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--coupling-j", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--time-points", type=int, default=5)
    parser.add_argument("--path-slices", type=int, default=60)
    parser.add_argument("--path-mixing", type=float, default=0.01)
    parser.add_argument("--path-tolerance", type=float, default=1e-8)
    parser.add_argument("--path-max-iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument(
        "--action-cache",
        type=Path,
        default=Path("/tmp/syk_ed_op_stat_compact_actions_N24.npz"),
    )
    parser.add_argument("--skip-saddle", action="store_true")
    parser.add_argument("--validate-compact", action="store_true")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=HERE / "pilot_outputs",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
