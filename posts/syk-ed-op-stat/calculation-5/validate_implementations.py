#!/usr/bin/env python3
"""Archive the compact, streaming, and cross-replica implementation checks."""

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


def streaming_error(large, calc1, calc4) -> float:
    n_majorana = 10
    beta = 0.7
    times = beta * np.arange(5) / 5
    measured_indices = (0, 1, 2, 3)
    basis = calc1.ParityBasis(n_majorana)
    _, even_actions, odd_actions = calc1.build_hamiltonian_actions(basis)
    h_even, h_odd = calc1.sample_hamiltonian_blocks(
        np.random.default_rng(271828),
        basis,
        even_actions,
        odd_actions,
        1.0,
    )
    reference = calc4.same_side_sample(
        calc1,
        h_even,
        h_odd,
        basis,
        measured_indices,
        beta,
        times,
    )
    candidate = large.same_side_sample_streaming(
        h_even,
        h_odd,
        basis,
        measured_indices,
        beta,
        times,
    )[:3]
    return max(
        float(np.max(np.abs(left - right)))
        for left, right in zip(candidate, reference)
    )


def run() -> dict:
    large = load_module(
        "validate_large_same_side",
        HERE / "compare_large_same_side.py",
    )
    calc1, calc4 = large.load_calculation_modules()
    cross = load_module(
        "validate_large_cross_replica",
        HERE / "compare_large_cross_replica.py",
    )
    output = {
        "compact_hamiltonian_max_error": large.validate_compact_hamiltonian(calc1),
        "streaming_same_side_max_error": streaming_error(large, calc1, calc4),
        "cross_replica_trace_max_error": cross.validate_evaluator(large, calc1),
    }
    output_path = HERE / "validation.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))
    return output


if __name__ == "__main__":
    run()
