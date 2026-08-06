#!/usr/bin/env python3
"""Compare the i-not-in-A ED correlator with the conditional AP propagator.

This is a reproducible post-processing calculation based on the raw zeta
samples and converged connected saddle saved by calculation 1. It separates
the antiperiodic propagator before the weighted P/AP average and compares its
replica-off-diagonal block directly with

    C_out(t,t') = average_{i not in A} E[zeta_Ai(t) zeta_Ai(t')].
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
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


@dataclass
class SectorSummary:
    n_majorana: int
    weight: int
    outside_fermions: int
    beta_j: float
    disorder_samples: int
    time_points: int
    path_slices: int
    path_relative_weight: float
    path_iterations: int | None
    path_relative_update: float | None
    path_converged: bool
    saddle_reconstruction_relative_error: float
    ap_time_average: float
    ap_centered_rms: float
    ed_outside_time_average: float
    ed_outside_centered_rms: float
    raw_fitted_amplitude: float
    raw_cosine_similarity: float
    raw_relative_residual: float
    raw_residual_to_ed_standard_error: float
    origin_normalized_cosine_similarity: float
    origin_normalized_relative_residual: float
    centered_fitted_amplitude: float
    centered_shape_correlation: float
    centered_relative_residual: float
    origin_normalized_residual_jackknife_se: float
    centered_residual_jackknife_se: float
    centered_correlation_jackknife_se: float
    ed_outside_signal_to_noise: float


def fit_through_origin(
    target: np.ndarray, predictor: np.ndarray
) -> tuple[float, float, float]:
    target_flat = np.asarray(target).ravel()
    predictor_flat = np.asarray(predictor).ravel()
    amplitude = float(
        np.dot(target_flat, predictor_flat) / np.dot(predictor_flat, predictor_flat)
    )
    cosine = float(
        np.dot(target_flat, predictor_flat)
        / (np.linalg.norm(target_flat) * np.linalg.norm(predictor_flat))
    )
    residual = float(
        np.linalg.norm(target_flat - amplitude * predictor_flat)
        / np.linalg.norm(target_flat)
    )
    return amplitude, cosine, residual


def conditional_propagators(
    green_blocks: np.ndarray,
    beta: float,
    coupling_j: float,
):
    calc1 = load_calculation_1_module()
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


def jackknife_shape_errors(
    sample_surfaces: np.ndarray, path_surface: np.ndarray
) -> tuple[float, float, float]:
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
    means = estimates_array.mean(axis=0)
    errors = np.sqrt(
        (len(estimates_array) - 1)
        * np.mean((estimates_array - means) ** 2, axis=0)
    )
    return float(errors[0]), float(errors[1]), float(errors[2])


def save_figure(
    path: Path,
    fractions: np.ndarray,
    ap_surface: np.ndarray,
    ed_surface: np.ndarray,
    ed_standard_error: np.ndarray,
    amplitude: float,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fitted = amplitude * ap_surface
    difference = ed_surface - fitted
    scale = max(np.max(np.abs(fitted)), np.max(np.abs(ed_surface)))
    error_scale = np.max(ed_standard_error)
    extent = [
        float(fractions[0]),
        float(fractions[-1]),
        float(fractions[0]),
        float(fractions[-1]),
    ]

    figure, axes = plt.subplots(1, 4, figsize=(12.8, 3.15), constrained_layout=True)
    panels = [
        (fitted, r"fitted $aG^{\mathrm{AP}}_{12}$", "RdBu_r", -scale, scale),
        (
            ed_surface,
            r"ED $C_{\mathrm{out}}$",
            "RdBu_r",
            -scale,
            scale,
        ),
        (difference, "ED minus AP fit", "RdBu_r", -scale, scale),
        (
            ed_standard_error,
            "ED standard error",
            "viridis",
            0.0,
            error_scale,
        ),
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
    figure.savefig(path, dpi=180)
    plt.close(figure)


def run(args: argparse.Namespace) -> SectorSummary:
    args.outdir.mkdir(parents=True, exist_ok=True)
    source = np.load(args.input)
    zeta_samples = source["zeta_samples"]
    fractions = source["fractions"]
    saved_green_blocks = source["path_G_blocks"]
    n_samples, n_majorana, time_points = zeta_samples.shape
    if not 0 < args.weight < n_majorana:
        raise ValueError("weight must lie strictly between zero and N")

    relative_weight = (
        args.path_weight
        if args.path_weight is not None
        else args.weight / n_majorana
    )
    if args.path_weight is None:
        green_blocks = saved_green_blocks
        path_iterations = None
        path_relative_update = None
        path_converged = True
    else:
        calc1 = load_calculation_1_module()
        saddle = calc1.solve_weighted_replica_saddle(
            beta=args.beta,
            relative_weight=relative_weight,
            length=saved_green_blocks.shape[1],
            coupling_j=args.coupling_j,
            mixing=args.path_mixing,
            tolerance=args.path_tolerance,
            max_iterations=args.path_max_iterations,
        )
        green_blocks = saddle.G_blocks
        path_iterations = saddle.iterations
        path_relative_update = saddle.relative_update
        path_converged = saddle.converged
    periodic, antiperiodic = conditional_propagators(
        green_blocks, args.beta, args.coupling_j
    )
    reconstructed = relative_weight * periodic + (1.0 - relative_weight) * antiperiodic
    reconstruction_error = float(
        np.linalg.norm(reconstructed - green_blocks) / np.linalg.norm(green_blocks)
    )

    path_slices = green_blocks.shape[1]
    path_indices = np.rint(fractions * path_slices).astype(int)
    ap_lattice = antiperiodic[1][np.ix_(path_indices, path_indices)]
    ap_surface = 0.5 * (ap_lattice + ap_lattice.T)

    outside = zeta_samples[:, args.weight :, :]
    outside_count = outside.shape[1]
    sample_surfaces = np.einsum(
        "sit,siu->stu", outside, outside, optimize=True
    ).real / outside_count
    ed_surface = sample_surfaces.mean(axis=0)
    ed_standard_error = (
        sample_surfaces.std(axis=0, ddof=1) / math.sqrt(n_samples)
        if n_samples > 1
        else np.full(ed_surface.shape, np.nan)
    )

    raw_amplitude, raw_cosine, raw_residual = fit_through_origin(
        ed_surface, ap_surface
    )
    normalized_ed = ed_surface / ed_surface[0, 0]
    normalized_ap = ap_surface / ap_surface[0, 0]
    normalized_cosine = float(
        np.dot(normalized_ed.ravel(), normalized_ap.ravel())
        / (np.linalg.norm(normalized_ed) * np.linalg.norm(normalized_ap))
    )
    normalized_residual = float(
        np.linalg.norm(normalized_ed - normalized_ap)
        / np.linalg.norm(normalized_ed)
    )
    centered_ed = ed_surface - ed_surface.mean()
    centered_ap = ap_surface - ap_surface.mean()
    centered_amplitude, centered_correlation, centered_residual = fit_through_origin(
        centered_ed, centered_ap
    )
    (
        origin_residual_error,
        centered_residual_error,
        centered_correlation_error,
    ) = jackknife_shape_errors(sample_surfaces, ap_surface)

    summary = SectorSummary(
        n_majorana=n_majorana,
        weight=args.weight,
        outside_fermions=outside_count,
        beta_j=args.beta * args.coupling_j,
        disorder_samples=n_samples,
        time_points=time_points,
        path_slices=path_slices,
        path_relative_weight=relative_weight,
        path_iterations=path_iterations,
        path_relative_update=path_relative_update,
        path_converged=path_converged,
        saddle_reconstruction_relative_error=reconstruction_error,
        ap_time_average=float(ap_surface.mean()),
        ap_centered_rms=float(np.sqrt(np.mean(centered_ap**2))),
        ed_outside_time_average=float(ed_surface.mean()),
        ed_outside_centered_rms=float(np.sqrt(np.mean(centered_ed**2))),
        raw_fitted_amplitude=raw_amplitude,
        raw_cosine_similarity=raw_cosine,
        raw_relative_residual=raw_residual,
        raw_residual_to_ed_standard_error=float(
            np.linalg.norm(ed_surface - raw_amplitude * ap_surface)
            / np.linalg.norm(ed_standard_error)
        ),
        origin_normalized_cosine_similarity=normalized_cosine,
        origin_normalized_relative_residual=normalized_residual,
        centered_fitted_amplitude=centered_amplitude,
        centered_shape_correlation=centered_correlation,
        centered_relative_residual=centered_residual,
        origin_normalized_residual_jackknife_se=origin_residual_error,
        centered_residual_jackknife_se=centered_residual_error,
        centered_correlation_jackknife_se=centered_correlation_error,
        ed_outside_signal_to_noise=float(
            np.linalg.norm(ed_surface) / np.linalg.norm(ed_standard_error)
        ),
    )

    np.savez_compressed(
        args.outdir / "outside_ap_data.npz",
        fractions=fractions,
        ap_surface=ap_surface,
        ap_surface_lattice=ap_lattice,
        ed_outside_surface=ed_surface,
        ed_outside_standard_error=ed_standard_error,
        ed_outside_sample_surfaces=sample_surfaces,
        periodic_G_blocks=periodic,
        antiperiodic_G_blocks=antiperiodic,
    )
    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(asdict(summary), handle, indent=2)
    save_figure(
        args.outdir / "outside_ap_comparison.png",
        fractions,
        ap_surface,
        ed_surface,
        ed_standard_error,
        raw_amplitude,
    )
    print(json.dumps(asdict(summary), indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=CALCULATION_1 / "outputs" / "comparison_data.npz",
    )
    parser.add_argument("--weight", type=int, default=3)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--coupling-j", type=float, default=1.0)
    parser.add_argument(
        "--path-weight",
        type=float,
        default=None,
        help="solve a connected saddle at this w instead of using the saved saddle",
    )
    parser.add_argument("--path-mixing", type=float, default=0.01)
    parser.add_argument("--path-tolerance", type=float, default=1e-8)
    parser.add_argument("--path-max-iterations", type=int, default=6000)
    parser.add_argument("--outdir", type=Path, default=Path("outputs"))
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
