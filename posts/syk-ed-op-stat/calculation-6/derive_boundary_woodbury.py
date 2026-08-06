#!/usr/bin/env python3
"""Check the exact finite-lattice relation between conditional P/AP kernels.

For two replicas the forward-difference matrices obey

    A_AP = A_P + 2 E F^T,

where A_s = D_s - dt^2 Sigma and E/F select the first/last time site on
each replica.  With G_s = -A_s^{-T}, Woodbury gives the exact identity

    C_AP = (I - 2 K_P)^{-1} C_P,

where C_P = E^T G_P E and K_P = E^T G_P F.  This script checks the identity
against saved saddles and records the approach to the continuum collision
product G^P_12(0,0) G^AP_12(0,0) = -1/4.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def replica_parameters(matrix: np.ndarray) -> tuple[float, float, float]:
    """Return diagonal, antisymmetric off-diagonal, and structure error."""
    diagonal = 0.5 * (matrix[0, 0] + matrix[1, 1])
    off_diagonal = 0.5 * (matrix[1, 0] - matrix[0, 1])
    reconstructed = np.array(
        [[diagonal, -off_diagonal], [off_diagonal, diagonal]]
    )
    error = float(np.max(np.abs(matrix - reconstructed)))
    return float(diagonal), float(off_diagonal), error


def check_grid(calc1, periodic: np.ndarray, antiperiodic: np.ndarray) -> dict:
    length = periodic.shape[1]
    periodic_full = calc1.replica_full_matrix(periodic)
    antiperiodic_full = calc1.replica_full_matrix(antiperiodic)
    first = np.asarray([0, length])
    last = np.asarray([length - 1, 2 * length - 1])

    collision_p = periodic_full[np.ix_(first, first)]
    corner_p = periodic_full[np.ix_(first, last)]
    collision_ap = antiperiodic_full[np.ix_(first, first)]
    woodbury = np.linalg.solve(
        np.eye(2) - 2.0 * corner_p,
        collision_p,
    )

    k, m_p, collision_structure_error = replica_parameters(collision_p)
    a, b, corner_structure_error = replica_parameters(corner_p)
    k_ap, m_ap, ap_structure_error = replica_parameters(collision_ap)
    denominator = (1.0 - 2.0 * a) ** 2 + 4.0 * b**2
    explicit_m_ap = (
        (1.0 - 2.0 * a) * m_p + 2.0 * b * k
    ) / denominator

    return {
        "length": length,
        "periodic_contact_k": k,
        "periodic_collision_m": m_p,
        "periodic_corner_a": a,
        "periodic_corner_b": b,
        "antiperiodic_contact_k": k_ap,
        "antiperiodic_collision_m": m_ap,
        "collision_product": m_p * m_ap,
        "collision_product_plus_quarter": m_p * m_ap + 0.25,
        "explicit_formula_m_ap": explicit_m_ap,
        "explicit_formula_error": abs(explicit_m_ap - m_ap),
        "woodbury_max_error": float(np.max(np.abs(woodbury - collision_ap))),
        "periodic_collision_structure_error": collision_structure_error,
        "periodic_corner_structure_error": corner_structure_error,
        "antiperiodic_collision_structure_error": ap_structure_error,
    }


def richardson_summary(coarse: dict, fine: dict) -> dict:
    keys = (
        "periodic_contact_k",
        "periodic_collision_m",
        "periodic_corner_a",
        "periodic_corner_b",
        "antiperiodic_contact_k",
        "antiperiodic_collision_m",
    )
    extrapolated = {
        key: 2.0 * fine[key] - coarse[key]
        for key in keys
    }
    m_p = extrapolated["periodic_collision_m"]
    m_ap = extrapolated["antiperiodic_collision_m"]
    extrapolated.update(
        {
            "collision_product": m_p * m_ap,
            "collision_product_plus_quarter": m_p * m_ap + 0.25,
            "contact_k_plus_half": extrapolated["periodic_contact_k"] + 0.5,
            "corner_a_minus_half": extrapolated["periodic_corner_a"] - 0.5,
            "corner_b_minus_collision_m": (
                extrapolated["periodic_corner_b"] - m_p
            ),
            "continuum_formula_m_ap": -1.0 / (4.0 * m_p),
            "continuum_formula_error": abs(m_ap + 1.0 / (4.0 * m_p)),
        }
    )
    return extrapolated


def conditionals_from_saved_saddle(calc1, calc4, path: Path, beta: float):
    with np.load(path) as data:
        result = []
        for label in ("coarse", "fine"):
            periodic, antiperiodic = calc4.conditional_propagators(
                calc1,
                data[f"saddle_G_blocks_{label}"],
                beta,
                1.0,
            )
            result.append((periodic, antiperiodic))
    return result


def run(outdir: Path) -> dict:
    calc1 = load_module(
        "woodbury_calc1",
        ENTRY / "calculation-1" / "compare_ed_path_integral.py",
    )
    calc4 = load_module(
        "woodbury_calc4",
        ENTRY / "calculation-4" / "compare_same_side.py",
    )

    cases = []
    with np.load(
        ENTRY / "calculation-4" / "outputs" / "same_side_data.npz"
    ) as data:
        cases.append(
            (
                "N20_beta0p5_w0p2",
                0.5,
                0.2,
                [
                    (
                        data["periodic_G_blocks_coarse"],
                        data["antiperiodic_G_blocks_coarse"],
                    ),
                    (
                        data["periodic_G_blocks_fine"],
                        data["antiperiodic_G_blocks_fine"],
                    ),
                ],
            )
        )

    for beta_key, beta in (("0p5", 0.5), ("1", 1.0)):
        source = (
            ENTRY
            / "calculation-5"
            / f"outputs_16_beta{beta_key}"
            / "large_same_side_data.npz"
        )
        cases.append(
            (
                f"N24_beta{beta_key}_w1over3",
                beta,
                1.0 / 3.0,
                conditionals_from_saved_saddle(calc1, calc4, source, beta),
            )
        )

    output = {"identity": "C_AP = (I - 2 K_P)^(-1) C_P", "cases": {}}
    for label, beta, weight, grids in cases:
        grid_checks = [
            check_grid(calc1, periodic, antiperiodic)
            for periodic, antiperiodic in grids
        ]
        output["cases"][label] = {
            "beta_j": beta,
            "relative_weight": weight,
            "grids": grid_checks,
            "richardson_2fine_minus_coarse": richardson_summary(*grid_checks),
        }

    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / "woodbury_check.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=HERE / "outputs")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args().outdir)
