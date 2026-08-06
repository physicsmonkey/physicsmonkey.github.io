#!/usr/bin/env python3
"""Centered-shape diagnostics for the saved N=24 cross-replica surfaces.

The periodic G_12 comparison is dominated by its collision-normalized constant.
This postprocessor removes the surface mean, fits one shape amplitude, and uses
delete-one-disorder-sample jackknife errors for the resulting global metrics.
The scalar jackknife retains all correlations among time points without trying
to invert a singular 81-by-81 covariance matrix from only 16 samples.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def centered_metrics(ed: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    ed_centered = ed - ed.mean()
    prediction_centered = prediction - prediction.mean()
    fitted_amplitude = float(
        np.sum(ed_centered * prediction_centered)
        / np.sum(prediction_centered**2)
    )
    fitted_difference = ed_centered - fitted_amplitude * prediction_centered
    return {
        "ed_mean": float(ed.mean()),
        "prediction_mean": float(prediction.mean()),
        "ed_centered_rms": float(np.sqrt(np.mean(ed_centered**2))),
        "prediction_centered_rms": float(
            np.sqrt(np.mean(prediction_centered**2))
        ),
        "centered_fitted_amplitude": fitted_amplitude,
        "centered_correlation": float(
            np.sum(ed_centered * prediction_centered)
            / (
                np.linalg.norm(ed_centered)
                * np.linalg.norm(prediction_centered)
            )
        ),
        "centered_fitted_relative_residual": float(
            np.linalg.norm(fitted_difference) / np.linalg.norm(ed_centered)
        ),
        "centered_no_fit_relative_residual": float(
            np.linalg.norm(ed_centered - prediction_centered)
            / np.linalg.norm(ed_centered)
        ),
    }


def jackknife_metrics(
    numerator_samples: np.ndarray,
    denominator_samples: np.ndarray,
    prediction: np.ndarray,
) -> dict[str, float]:
    estimates = []
    n_samples = len(denominator_samples)
    for omitted in range(n_samples):
        keep = np.arange(n_samples) != omitted
        ed = (
            numerator_samples[keep].mean(axis=0)
            / denominator_samples[keep].mean()
        )
        metrics = centered_metrics(ed, prediction)
        estimates.append(
            [
                metrics["centered_fitted_amplitude"],
                metrics["centered_correlation"],
                metrics["centered_fitted_relative_residual"],
                metrics["centered_no_fit_relative_residual"],
            ]
        )
    estimates = np.asarray(estimates)
    jackknife_mean = estimates.mean(axis=0)
    standard_errors = np.sqrt(
        (n_samples - 1)
        * np.mean((estimates - jackknife_mean) ** 2, axis=0)
    )
    names = (
        "centered_fitted_amplitude_jackknife_se",
        "centered_correlation_jackknife_se",
        "centered_fitted_relative_residual_jackknife_se",
        "centered_no_fit_relative_residual_jackknife_se",
    )
    return {name: float(value) for name, value in zip(names, standard_errors)}


def run(data_path: Path, summary_path: Path, output_path: Path) -> dict:
    summary = json.loads(summary_path.read_text())
    output = {
        "disorder_samples": summary["disorder_samples"],
        "method": (
            "global centered metrics with delete-one-sample jackknife; "
            "the resampling preserves correlations among all time points"
        ),
        "covariance_note": (
            "A full 81-entry covariance has rank at most 15 with 16 samples, "
            "so no pseudoinverse chi-square or p-value is reported."
        ),
        "temperatures": {},
    }
    with np.load(data_path) as data:
        for beta_key, temperature in summary["temperatures"].items():
            temperature_output = {"beta_j": temperature["beta_j"], "sectors": {}}
            for label, sector in temperature["sectors"].items():
                ratio = data[f"ratio_beta{beta_key}_{label}"]
                numerator_samples = data[f"numerator_samples_beta{beta_key}_{label}"]
                denominator_samples = data[
                    f"denominator_samples_beta{beta_key}_{label}"
                ]
                path_output = {}
                for weight_key, path_case in sector["path_cases"].items():
                    prediction = data[
                        f"prediction_beta{beta_key}_{label}_w{weight_key}"
                    ]
                    path_output[weight_key] = {
                        "relative_weight": path_case["relative_weight"],
                        **centered_metrics(ratio, prediction),
                        **jackknife_metrics(
                            numerator_samples,
                            denominator_samples,
                            prediction,
                        ),
                    }
                temperature_output["sectors"][label] = {
                    "background_weight": sector["background_weight"],
                    "path_cases": path_output,
                }
            output["temperatures"][beta_key] = temperature_output

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    source = HERE / "cross_replica_outputs_16"
    parser.add_argument(
        "--data",
        type=Path,
        default=source / "large_cross_replica_data.npz",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=source / "summary.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=source / "shape_diagnostics.json",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.data, args.summary, args.output)
