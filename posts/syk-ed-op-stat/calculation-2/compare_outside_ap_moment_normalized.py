#!/usr/bin/env python3
"""Moment-normalized version of the calculation-2 AP comparison."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE / "cross_replica_moment_normalization.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "cross_replica_moment_normalization_calc2", MODULE_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


if __name__ == "__main__":
    moment = load_module()
    parser = moment.parser_with_defaults(
        default_input=HERE / "ed_inputs_100" / "comparison_data.npz",
        default_weight=3,
        default_sector="outside",
        default_path_weights=[0.15, 0.20],
        default_probe=None,
        default_outdir=HERE / "outputs_100_moment_normalized",
    )
    moment.run(parser.parse_args())
