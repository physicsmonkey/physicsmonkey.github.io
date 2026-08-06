#!/usr/bin/env python3
"""Validate the spectral same-side evaluator against direct dense traces."""

from __future__ import annotations

import importlib.util
import json
import math
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


def dense_string(calc1, basis, indices: tuple[int, ...]) -> np.ndarray:
    states = np.arange(basis.dim)
    rows, phases = calc1.string_action(states, indices, basis.n_majorana)
    matrix = np.zeros((basis.dim, basis.dim), dtype=complex)
    matrix[rows, states] = phases
    return matrix


def direct_surfaces(calc1, h_full, basis, measured_indices, beta, times):
    energies, vectors = np.linalg.eigh(h_full)

    def kernel(time: float) -> np.ndarray:
        return (vectors * np.exp(-time * energies)) @ vectors.conj().T

    mu = dense_string(calc1, basis, measured_indices)
    x_a = np.trace(kernel(beta) @ mu)
    surfaces = []
    for group in (
        range(len(measured_indices)),
        range(len(measured_indices), basis.n_majorana),
    ):
        surface = np.zeros((len(times), len(times)), dtype=complex)
        np.fill_diagonal(surface, -0.5 * x_a)
        for time_index, tau in enumerate(times):
            for prime_index in range(time_index):
                tau_prime = times[prime_index]
                values = []
                for index in group:
                    psi = dense_string(calc1, basis, (index,)) / math.sqrt(2.0)
                    values.append(
                        np.trace(
                            kernel(beta - tau)
                            @ psi
                            @ kernel(tau - tau_prime)
                            @ psi
                            @ kernel(tau_prime)
                            @ mu
                        )
                    )
                value = np.mean(values)
                surface[time_index, prime_index] = value
                surface[prime_index, time_index] = -value
        surfaces.append(surface)
    return x_a, surfaces[0], surfaces[1]


def run() -> dict:
    calc1 = load_module(
        "validate_same_side_calc1",
        HERE.parent / "calculation-1" / "compare_ed_path_integral.py",
    )
    calc4 = load_module(
        "validate_same_side_calc4",
        HERE / "compare_same_side.py",
    )
    n_majorana = 6
    beta = 0.7
    times = beta * np.arange(4) / 4
    measured_indices = (0, 1, 2, 3)
    basis = calc1.ParityBasis(n_majorana)
    _, even_actions, odd_actions = calc1.build_hamiltonian_actions(basis)
    h_even, h_odd = calc1.sample_hamiltonian_blocks(
        np.random.default_rng(314159),
        basis,
        even_actions,
        odd_actions,
        1.0,
    )
    candidate = calc4.same_side_sample(
        calc1,
        h_even,
        h_odd,
        basis,
        measured_indices,
        beta,
        times,
    )
    h_full = np.zeros((basis.dim, basis.dim), dtype=complex)
    h_full[np.ix_(basis.even, basis.even)] = h_even
    h_full[np.ix_(basis.odd, basis.odd)] = h_odd
    reference = direct_surfaces(
        calc1,
        h_full,
        basis,
        measured_indices,
        beta,
        times,
    )
    component_errors = [
        float(np.max(np.abs(left - right)))
        for left, right in zip(candidate, reference)
    ]
    output = {
        "n_majorana": n_majorana,
        "beta_j": beta,
        "time_points": len(times),
        "x_a_max_error": component_errors[0],
        "inside_surface_max_error": component_errors[1],
        "outside_surface_max_error": component_errors[2],
        "overall_max_error": max(component_errors),
    }
    output_path = HERE / "validation.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))
    return output


if __name__ == "__main__":
    run()
