#!/usr/bin/env python3
"""Compare probes with i in A against the conditional periodic propagator.

The principal case is A={1,2,3,4,5} and i=5. At coincident time the probe
removes gamma_5 from mu_A and leaves a weight-four string, which couples to the
q=4 Hamiltonian at first order in beta J. We compare both the fixed i=5 probe
and the average over all five i in A with connected P-sector propagators at
w=5/20 and w=4/20.
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


def fit_metrics(
    ed_surface: np.ndarray,
    path_surface: np.ndarray,
    ed_standard_error: np.ndarray,
) -> dict[str, float]:
    ed = np.asarray(ed_surface).real
    path = np.asarray(path_surface).real
    amplitude = float(np.sum(ed * path) / np.sum(path * path))
    residual_surface = ed - amplitude * path
    raw_cosine = float(
        np.sum(ed * path) / (np.linalg.norm(ed) * np.linalg.norm(path))
    )

    ed_normalized = ed / ed[0, 0]
    path_normalized = path / path[0, 0]
    origin_cosine = float(
        np.sum(ed_normalized * path_normalized)
        / (np.linalg.norm(ed_normalized) * np.linalg.norm(path_normalized))
    )

    ed_centered = ed - ed.mean()
    path_centered = path - path.mean()
    centered_amplitude = float(
        np.sum(ed_centered * path_centered) / np.sum(path_centered * path_centered)
    )
    centered_residual_surface = ed_centered - centered_amplitude * path_centered
    centered_correlation = float(
        np.sum(ed_centered * path_centered)
        / (np.linalg.norm(ed_centered) * np.linalg.norm(path_centered))
    )

    return {
        "ed_time_average": float(ed.mean()),
        "ed_centered_rms": float(np.sqrt(np.mean(ed_centered**2))),
        "path_time_average": float(path.mean()),
        "path_centered_rms": float(np.sqrt(np.mean(path_centered**2))),
        "raw_fitted_amplitude": amplitude,
        "raw_cosine_similarity": raw_cosine,
        "raw_relative_residual": float(
            np.linalg.norm(residual_surface) / np.linalg.norm(ed)
        ),
        "raw_residual_to_ed_standard_error": float(
            np.linalg.norm(residual_surface) / np.linalg.norm(ed_standard_error)
        ),
        "origin_normalized_cosine_similarity": origin_cosine,
        "origin_normalized_relative_residual": float(
            np.linalg.norm(ed_normalized - path_normalized)
            / np.linalg.norm(ed_normalized)
        ),
        "centered_fitted_amplitude": centered_amplitude,
        "centered_shape_correlation": centered_correlation,
        "centered_relative_residual": float(
            np.linalg.norm(centered_residual_surface) / np.linalg.norm(ed_centered)
        ),
        "ed_signal_to_noise": float(
            np.linalg.norm(ed) / np.linalg.norm(ed_standard_error)
        ),
    }


def jackknife_shape_errors(
    sample_surfaces: np.ndarray, path_surface: np.ndarray
) -> dict[str, float]:
    """Jackknife errors for the origin-normalized and centered comparisons."""
    estimates = []
    for omitted in range(len(sample_surfaces)):
        ed = np.delete(sample_surfaces, omitted, axis=0).mean(axis=0)
        ed_normalized = ed / ed[0, 0]
        path_normalized = path_surface / path_surface[0, 0]
        origin_residual = float(
            np.linalg.norm(ed_normalized - path_normalized)
            / np.linalg.norm(ed_normalized)
        )
        ed_centered = ed - ed.mean()
        path_centered = path_surface - path_surface.mean()
        amplitude = float(
            np.sum(ed_centered * path_centered)
            / np.sum(path_centered * path_centered)
        )
        centered_residual = float(
            np.linalg.norm(ed_centered - amplitude * path_centered)
            / np.linalg.norm(ed_centered)
        )
        centered_correlation = float(
            np.sum(ed_centered * path_centered)
            / (np.linalg.norm(ed_centered) * np.linalg.norm(path_centered))
        )
        estimates.append(
            [origin_residual, centered_residual, centered_correlation]
        )
    estimates_array = np.asarray(estimates)
    jackknife_mean = estimates_array.mean(axis=0)
    standard_errors = np.sqrt(
        (len(estimates_array) - 1)
        * np.mean((estimates_array - jackknife_mean) ** 2, axis=0)
    )
    return {
        "origin_normalized_residual_jackknife_se": float(standard_errors[0]),
        "centered_residual_jackknife_se": float(standard_errors[1]),
        "centered_correlation_jackknife_se": float(standard_errors[2]),
    }


def surface_statistics(
    zeta_samples: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-sample products, their mean, and standard error."""
    sample_surfaces = np.einsum(
        "sit,siu->stu", zeta_samples, zeta_samples, optimize=True
    ).real / zeta_samples.shape[1]
    mean = sample_surfaces.mean(axis=0)
    standard_error = (
        sample_surfaces.std(axis=0, ddof=1) / math.sqrt(len(sample_surfaces))
        if len(sample_surfaces) > 1
        else np.full(mean.shape, np.nan)
    )
    return sample_surfaces, mean, standard_error


def conditional_periodic(
    calc1,
    green_blocks: np.ndarray,
    beta: float,
    coupling_j: float,
) -> tuple[np.ndarray, np.ndarray, float]:
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
    return periodic, antiperiodic, delta_tau_squared


def save_comparison_figure(
    path: Path,
    fractions: np.ndarray,
    p_surface: np.ndarray,
    cases: list[
        tuple[str, np.ndarray, np.ndarray, dict[str, float]]
    ],
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
        len(cases), 4, figsize=(12.8, 3.05 * len(cases)), constrained_layout=True
    )
    if len(cases) == 1:
        axes = axes[None, :]

    for row, (label, ed, standard_error, metrics) in enumerate(cases):
        fitted = metrics["raw_fitted_amplitude"] * p_surface
        difference = ed - fitted
        scale = max(np.max(np.abs(fitted)), np.max(np.abs(ed)))
        panels = [
            (fitted, rf"fitted $aG^{{\rm P}}_{{12}}$ ({label})", "RdBu_r", -scale, scale),
            (ed, rf"ED {label}", "RdBu_r", -scale, scale),
            (difference, "ED minus P fit", "RdBu_r", -scale, scale),
            (
                standard_error,
                "ED standard error",
                "viridis",
                0.0,
                float(np.max(standard_error)),
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


def save_normalized_figure(
    path: Path,
    fractions: np.ndarray,
    p_surface: np.ndarray,
    cases: list[tuple[str, np.ndarray]],
) -> None:
    """Show small deviations from the value at the origin."""
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p_deviation = p_surface / p_surface[0, 0] - 1.0
    extent = [
        float(fractions[0]),
        float(fractions[-1]),
        float(fractions[0]),
        float(fractions[-1]),
    ]
    figure, axes = plt.subplots(
        len(cases), 3, figsize=(9.8, 3.05 * len(cases)), constrained_layout=True
    )
    if len(cases) == 1:
        axes = axes[None, :]
    for row, (label, ed) in enumerate(cases):
        ed_deviation = ed / ed[0, 0] - 1.0
        difference = ed_deviation - p_deviation
        scale = max(np.max(np.abs(p_deviation)), np.max(np.abs(ed_deviation)))
        panels = [
            (p_deviation, r"$G^{\rm P}_{12}/G^{\rm P}_{12}(0,0)-1$"),
            (ed_deviation, rf"ED {label}, normalized minus 1"),
            (difference, "normalized ED minus P"),
        ]
        for axis, (data, title) in zip(axes[row], panels):
            image = axis.imshow(
                data,
                origin="lower",
                extent=extent,
                aspect="equal",
                cmap="RdBu_r",
                vmin=-scale,
                vmax=scale,
            )
            axis.set_title(title, fontsize=10)
            axis.set_xlabel(r"$\tau'/\beta$")
            axis.set_ylabel(r"$\tau/\beta$")
            figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def format_weight(weight: float) -> str:
    return f"{weight:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def run(args: argparse.Namespace) -> dict:
    calc1 = load_calculation_1_module()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.weight % 2 != 1:
        raise ValueError("A must have odd weight")
    if args.probe < 1 or args.probe > args.weight:
        raise ValueError("the one-based probe index must belong to A")

    basis = calc1.ParityBasis(args.n_majorana)
    _, even_actions, odd_actions = calc1.build_hamiltonian_actions(basis)
    rng = np.random.default_rng(args.seed)
    fractions = np.arange(args.time_points, dtype=float) / args.time_points
    times = args.beta * fractions
    measured_indices = tuple(range(args.weight))
    zeta_samples = np.empty(
        (args.samples, args.n_majorana, args.time_points), dtype=complex
    )

    for sample_index in range(args.samples):
        h_even, h_odd = calc1.sample_hamiltonian_blocks(
            rng, basis, even_actions, odd_actions, args.coupling_j
        )
        zeta_samples[sample_index] = calc1.zeta_grid(
            h_even,
            h_odd,
            basis,
            measured_indices,
            args.beta,
            times,
        )
        print(f"ED sample {sample_index + 1}/{args.samples}", flush=True)

    fixed_samples = zeta_samples[:, args.probe - 1 : args.probe, :]
    inside_samples = zeta_samples[:, : args.weight, :]
    fixed_sample_surfaces, fixed_surface, fixed_error = surface_statistics(
        fixed_samples
    )
    inside_sample_surfaces, inside_surface, inside_error = surface_statistics(
        inside_samples
    )

    summary: dict = {
        "n_majorana": args.n_majorana,
        "weight": args.weight,
        "probe_one_based": args.probe,
        "beta_j": args.beta * args.coupling_j,
        "disorder_samples": args.samples,
        "time_points": args.time_points,
        "path_slices": args.path_slices,
        "seed": args.seed,
        "path_cases": {},
    }
    arrays: dict[str, np.ndarray] = {
        "fractions": fractions,
        "times": times,
        "zeta_samples": zeta_samples,
        "fixed_probe_surface": fixed_surface,
        "fixed_probe_standard_error": fixed_error,
        "fixed_probe_sample_surfaces": fixed_sample_surfaces,
        "inside_average_surface": inside_surface,
        "inside_average_standard_error": inside_error,
        "inside_average_sample_surfaces": inside_sample_surfaces,
    }

    for relative_weight in args.path_weights:
        saddle = calc1.solve_weighted_replica_saddle(
            beta=args.beta,
            relative_weight=relative_weight,
            length=args.path_slices,
            coupling_j=args.coupling_j,
            mixing=args.path_mixing,
            tolerance=args.path_tolerance,
            max_iterations=args.path_max_iterations,
        )
        periodic, antiperiodic, _ = conditional_periodic(
            calc1, saddle.G_blocks, args.beta, args.coupling_j
        )
        reconstructed = (
            relative_weight * periodic + (1.0 - relative_weight) * antiperiodic
        )
        reconstruction_error = float(
            np.linalg.norm(reconstructed - saddle.G_blocks)
            / np.linalg.norm(saddle.G_blocks)
        )
        path_indices = np.rint(fractions * args.path_slices).astype(int)
        p_lattice = periodic[1][np.ix_(path_indices, path_indices)]
        p_surface = 0.5 * (p_lattice + p_lattice.T)
        key = format_weight(relative_weight)
        fixed_metrics = fit_metrics(fixed_surface, p_surface, fixed_error)
        inside_metrics = fit_metrics(inside_surface, p_surface, inside_error)
        fixed_metrics.update(
            jackknife_shape_errors(fixed_sample_surfaces, p_surface)
        )
        inside_metrics.update(
            jackknife_shape_errors(inside_sample_surfaces, p_surface)
        )
        summary["path_cases"][key] = {
            "relative_weight": relative_weight,
            "path_iterations": saddle.iterations,
            "path_relative_update": saddle.relative_update,
            "path_converged": saddle.converged,
            "saddle_reconstruction_relative_error": reconstruction_error,
            "fixed_probe": fixed_metrics,
            "inside_average": inside_metrics,
        }
        arrays[f"periodic_surface_w{key}"] = p_surface
        arrays[f"periodic_surface_lattice_w{key}"] = p_lattice
        arrays[f"periodic_G_blocks_w{key}"] = periodic
        arrays[f"antiperiodic_G_blocks_w{key}"] = antiperiodic
        arrays[f"saddle_G_blocks_w{key}"] = saddle.G_blocks

        figure_cases = [
            (rf"$C_{{i={args.probe}}}$", fixed_surface, fixed_error, fixed_metrics),
            (r"$C_{\mathrm{in}}$", inside_surface, inside_error, inside_metrics),
        ]
        save_comparison_figure(
            args.outdir / f"inside_p_comparison_w{key}.png",
            fractions,
            p_surface,
            figure_cases,
        )
        save_normalized_figure(
            args.outdir / f"inside_p_normalized_w{key}.png",
            fractions,
            p_surface,
            [
                (rf"$C_{{i={args.probe}}}$", fixed_surface),
                (r"$C_{\mathrm{in}}$", inside_surface),
            ],
        )

    np.savez_compressed(args.outdir / "inside_p_data.npz", **arrays)
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-majorana", type=int, default=20)
    parser.add_argument("--weight", type=int, default=5)
    parser.add_argument("--probe", type=int, default=5)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--coupling-j", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=24)
    parser.add_argument("--time-points", type=int, default=9)
    parser.add_argument("--path-slices", type=int, default=80)
    parser.add_argument("--path-weights", type=float, nargs="+", default=[0.25, 0.20])
    parser.add_argument("--path-mixing", type=float, default=0.01)
    parser.add_argument("--path-tolerance", type=float, default=1e-8)
    parser.add_argument("--path-max-iterations", type=int, default=6000)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--outdir", type=Path, default=Path("outputs"))
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
