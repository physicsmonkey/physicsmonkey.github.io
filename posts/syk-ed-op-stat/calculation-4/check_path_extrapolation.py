#!/usr/bin/env python3
"""Reproduce the L_tau=180,360 same-side P extrapolation check."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run() -> dict:
    calc1 = load_module(
        "path_check_calc1",
        HERE.parent / "calculation-1" / "compare_ed_path_integral.py",
    )
    calc4 = load_module("path_check_calc4", HERE / "compare_same_side.py")
    beta = 0.5
    relative_weight = 0.2
    with np.load(HERE / "outputs" / "same_side_data.npz") as data:
        ed = data["inside_ratio"].copy()
        error = data["inside_standard_error"].copy()
        periodic_180 = data["periodic_G11_fine"].copy()

    saddle_360 = calc1.solve_weighted_replica_saddle(
        beta=beta,
        relative_weight=relative_weight,
        length=360,
        coupling_j=1.0,
        mixing=0.01,
        tolerance=1e-8,
        max_iterations=6000,
    )
    if not saddle_360.converged:
        raise RuntimeError("L_tau=360 saddle did not converge")
    periodic_360, _ = calc4.conditional_propagators(
        calc1,
        saddle_360.G_blocks,
        beta,
        1.0,
    )
    indices = np.arange(ed.shape[0]) * 360 // ed.shape[0]
    periodic_360_sampled = periodic_360[0][np.ix_(indices, indices)]
    extrapolated = 2.0 * periodic_360_sampled - periodic_180
    metrics = calc4.comparison_metrics(ed, extrapolated, error)
    output = {
        "beta_j": beta,
        "relative_weight": relative_weight,
        "coarse_length": 180,
        "fine_length": 360,
        "continuum_extrapolation": "2 G_360 - G_180",
        "fine_iterations": saddle_360.iterations,
        "fine_relative_update": saddle_360.relative_update,
        "fine_converged": saddle_360.converged,
        "inside_p": metrics,
    }
    output_path = HERE / "outputs" / "path_extrapolation_180_360.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))
    return output


if __name__ == "__main__":
    run()
