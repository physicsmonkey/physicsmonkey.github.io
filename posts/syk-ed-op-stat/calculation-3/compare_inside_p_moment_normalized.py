#!/usr/bin/env python3
"""Moment-normalized version of the calculation-3 P comparison."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
MODULE_PATH = (
    HERE.parent / "calculation-2" / "cross_replica_moment_normalization.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "cross_replica_moment_normalization_calc3", MODULE_PATH
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
        default_input=HERE / "outputs_100" / "inside_p_data.npz",
        default_weight=5,
        default_sector="inside",
        default_path_weights=[0.25, 0.20],
        default_probe=5,
        default_outdir=HERE / "outputs_100_moment_normalized",
    )
    moment.run(parser.parse_args())
