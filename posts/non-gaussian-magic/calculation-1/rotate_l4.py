#!/usr/bin/env python3
"""Numerical test of the O(N) variability of the quartic 4-form norm.

The entry asks how much

    sum_{i<j<k<l} (O.J)_{ijkl}^4

can be changed by an orthogonal rotation O, for a random SYK-like antisymmetric
4-tensor J.  The overall SYK variance is irrelevant for the ratios reported
here, so the independent components are sampled as standard Gaussians.

The optimizer is a Jacobi-style coordinate search over Givens rotations.  It is
not a certificate of the global optimum, but it is reproducible and gives a
useful stress test of how much the alpha=2 quantity moves under accessible
Gaussian rotations.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


@dataclass
class TrialResult:
    n: int
    seed: int
    m_components: int
    initial_l2: float
    initial_l4: float
    equal_component_lower_bound: float
    best_min_l4: float
    best_max_l4: float
    min_ratio: float
    max_ratio: float
    lower_bound_ratio: float
    min_to_lower_bound: float
    max_to_lower_bound: float


def independent_combos(n: int) -> list[tuple[int, int, int, int]]:
    return list(itertools.combinations(range(n), 4))


def permutation_sign(items: tuple[int, ...]) -> int:
    inversions = 0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            inversions += items[i] > items[j]
    return -1 if inversions % 2 else 1


def random_antisymmetric_4form(n: int, rng: np.random.Generator) -> np.ndarray:
    tensor = np.zeros((n, n, n, n), dtype=np.float64)
    for combo in independent_combos(n):
        value = rng.normal()
        for perm in itertools.permutations(combo):
            tensor[perm] = permutation_sign(perm) * value
    return tensor


def independent_values(tensor: np.ndarray, combos: list[tuple[int, int, int, int]]) -> np.ndarray:
    return np.array([tensor[c] for c in combos], dtype=np.float64)


def moments(tensor: np.ndarray, combos: list[tuple[int, int, int, int]]) -> tuple[float, float]:
    values = independent_values(tensor, combos)
    return float(np.sum(values * values)), float(np.sum(values**4))


def rotate_one_axis(tensor: np.ndarray, p: int, q: int, theta: float, axis: int) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    out = np.array(tensor, copy=True)
    p_slice = [slice(None)] * tensor.ndim
    q_slice = [slice(None)] * tensor.ndim
    p_slice[axis] = p
    q_slice[axis] = q
    p_tuple = tuple(p_slice)
    q_tuple = tuple(q_slice)
    old_p = tensor[p_tuple]
    old_q = tensor[q_tuple]
    out[p_tuple] = c * old_p + s * old_q
    out[q_tuple] = -s * old_p + c * old_q
    return out


def rotate_pair(tensor: np.ndarray, p: int, q: int, theta: float) -> np.ndarray:
    rotated = tensor
    for axis in range(4):
        rotated = rotate_one_axis(rotated, p, q, theta, axis)
    return rotated


def optimize_by_givens(
    tensor: np.ndarray,
    combos: list[tuple[int, int, int, int]],
    mode: str,
    sweeps: int,
    rng: np.random.Generator,
    tol: float = 1e-10,
) -> tuple[np.ndarray, float]:
    if mode not in {"min", "max"}:
        raise ValueError("mode must be 'min' or 'max'")

    pairs = list(itertools.combinations(range(tensor.shape[0]), 2))
    current = np.array(tensor, copy=True)
    _, current_l4 = moments(current, combos)

    for _ in range(sweeps):
        before_sweep = current_l4
        rng.shuffle(pairs)
        for p, q in pairs:
            def signed_objective(theta: float) -> float:
                _, l4 = moments(rotate_pair(current, p, q, theta), combos)
                return l4 if mode == "min" else -l4

            # A coarse scan avoids getting trapped by the local interval around
            # theta=0, then bounded scalar minimization refines the best angle.
            grid = np.linspace(-0.5 * math.pi, 0.5 * math.pi, 25, endpoint=False)
            signed_values = np.array([signed_objective(theta) for theta in grid])
            best_index = int(np.argmin(signed_values))
            center = float(grid[best_index])
            step = math.pi / 25.0
            result = minimize_scalar(
                signed_objective,
                bounds=(center - step, center + step),
                method="bounded",
                options={"xatol": 1e-5},
            )
            candidate_l4 = result.fun if mode == "min" else -result.fun
            improved = candidate_l4 < current_l4 - tol if mode == "min" else candidate_l4 > current_l4 + tol
            if improved:
                current = rotate_pair(current, p, q, float(result.x))
                current_l4 = float(candidate_l4)

        relative_change = abs(current_l4 - before_sweep) / max(abs(before_sweep), 1.0)
        if relative_change < 1e-6:
            break

    return current, current_l4


def run_trial(n: int, seed: int, sweeps: int) -> TrialResult:
    combos = independent_combos(n)
    rng = np.random.default_rng(seed)
    tensor = random_antisymmetric_4form(n, rng)
    initial_l2, initial_l4 = moments(tensor, combos)
    lower_bound = initial_l2 * initial_l2 / len(combos)

    # Use separate optimizer RNGs so min/max are reproducible and independent.
    _, min_l4 = optimize_by_givens(tensor, combos, "min", sweeps, np.random.default_rng(seed + 10_000))
    _, max_l4 = optimize_by_givens(tensor, combos, "max", sweeps, np.random.default_rng(seed + 20_000))

    return TrialResult(
        n=n,
        seed=seed,
        m_components=len(combos),
        initial_l2=initial_l2,
        initial_l4=initial_l4,
        equal_component_lower_bound=lower_bound,
        best_min_l4=min_l4,
        best_max_l4=max_l4,
        min_ratio=min_l4 / initial_l4,
        max_ratio=max_l4 / initial_l4,
        lower_bound_ratio=lower_bound / initial_l4,
        min_to_lower_bound=min_l4 / lower_bound,
        max_to_lower_bound=max_l4 / lower_bound,
    )


def write_csv(path: Path, results: list[TrialResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(results[0]).keys()))
        writer.writeheader()
        for result in results:
            writer.writerow(asdict(result))


def write_summary(path: Path, results: list[TrialResult], args: argparse.Namespace) -> None:
    grouped: dict[int, list[TrialResult]] = {}
    for result in results:
        grouped.setdefault(result.n, []).append(result)

    summary = {
        "description": "Givens-rotation coordinate search for the O(N) variability of sum J_ijkl^4.",
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "by_n": {},
    }
    for n, n_results in grouped.items():
        summary["by_n"][str(n)] = {}
        for key in ["min_ratio", "max_ratio", "lower_bound_ratio", "min_to_lower_bound", "max_to_lower_bound"]:
            values = np.array([getattr(result, key) for result in n_results], dtype=float)
            summary["by_n"][str(n)][key] = {
                "mean": float(np.mean(values)),
                "std": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }

    path.write_text(json.dumps(summary, indent=2) + "\n")


def write_plot(path: Path, results: list[TrialResult]) -> None:
    grouped: dict[int, list[TrialResult]] = {}
    for result in results:
        grouped.setdefault(result.n, []).append(result)
    ns = sorted(grouped)

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for key, label, marker in [
        ("min_ratio", "optimized minimum / initial", "o"),
        ("max_ratio", "optimized maximum / initial", "s"),
        ("lower_bound_ratio", "equal-component lower bound / initial", "^"),
    ]:
        means = []
        errors = []
        for n in ns:
            values = np.array([getattr(result, key) for result in grouped[n]], dtype=float)
            means.append(float(np.mean(values)))
            errors.append(float(np.std(values, ddof=1)) if len(values) > 1 else 0.0)
        ax.errorbar(ns, means, yerr=errors, marker=marker, capsize=3, label=label)

    ax.axhline(1.0, color="black", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("N Majorana modes")
    ax.set_ylabel(r"ratio for $\sum_{i<j<k<l} J_{ijkl}^4$")
    ax.set_title(r"$O(N)$ variability of the $\alpha=2$ weight-4 sum")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-values", nargs="+", type=int, default=[8, 10, 12, 14])
    parser.add_argument("--trials", type=int, default=8)
    parser.add_argument("--sweeps", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    results: list[TrialResult] = []
    for n in args.n_values:
        if n < 4:
            raise ValueError("all N values must be at least 4")
        for trial in range(args.trials):
            seed = args.seed + 1000 * n + trial
            result = run_trial(n, seed, args.sweeps)
            results.append(result)
            print(
                f"N={n:2d} seed={seed}: "
                f"min/initial={result.min_ratio:.4f}, "
                f"max/initial={result.max_ratio:.4f}, "
                f"lower/initial={result.lower_bound_ratio:.4f}",
                flush=True,
            )

    write_csv(args.out_dir / "givens_search_results.csv", results)
    write_summary(args.out_dir / "summary.json", results, args)
    write_plot(args.out_dir / "ratio_plot.png", results)


if __name__ == "__main__":
    main()
