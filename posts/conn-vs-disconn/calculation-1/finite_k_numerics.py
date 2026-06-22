#!/usr/bin/env python3
r"""Finite-K oscillator numerics for the connected/disconnected toy model.

The calculation uses the even Fock-space truncation of

    H = omega a^\dagger a + (gamma / K) (a + a^\dagger)^4

and compares the exact cutoff crossover with the finite-K perturbative formula
and the large-K prediction discussed in the entry.  The script assumes NumPy,
SciPy, and Matplotlib are available.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
from typing import Callable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq


OMEGA = 1.0
S = 1.0
DELTA = 1.0
K_VALUES = list(range(1, 17))
NMAX_VALUES = [80, 120, 160, 220]
FINAL_NMAX = NMAX_VALUES[-1]
GAMMA_MIN = 1.0e-4
GAMMA_MAX = 3.0
SCAN_POINTS = 121
OPTION_PLOT_K = 10


def build_hamiltonian(k: int, gamma: float, nmax: int) -> np.ndarray:
    """Build the even-sector cutoff Hamiltonian in the |0>, |2>, ... basis."""
    if nmax < 4 * k:
        raise ValueError(f"nmax={nmax} is too small for K={k}")

    m = np.arange(0, nmax + 1, 2, dtype=float)
    dim = len(m)
    g = gamma / k

    diag = OMEGA * m + g * (6.0 * m * m + 6.0 * m + 3.0)
    h = np.diag(diag)

    m1 = m[:-1]
    off1 = g * (4.0 * m1 + 6.0) * np.sqrt((m1 + 1.0) * (m1 + 2.0))
    idx = np.arange(dim - 1)
    h[idx, idx + 1] = off1
    h[idx + 1, idx] = off1

    m2 = m[:-2]
    off2 = g * np.sqrt((m2 + 1.0) * (m2 + 2.0) * (m2 + 3.0) * (m2 + 4.0))
    idx = np.arange(dim - 2)
    h[idx, idx + 2] = off2
    h[idx + 2, idx] = off2

    return h


def exact_log_weighted_ratio(k: int, gamma: float, nmax: int) -> float:
    """Return log(A2 exp(-K delta) / A1) in the cutoff Hilbert space."""
    log_a1, log_a2_weighted = exact_log_options(k, gamma, nmax)
    return log_a2_weighted - log_a1


def exact_log_options(k: int, gamma: float, nmax: int) -> tuple[float, float]:
    r"""Return common-factor-stripped log(A1), log(A2 exp(-K delta)).

    The omitted common factor is (4K)! because
    (a^\dagger)^(4K)|0> = sqrt((4K)!) |4K>.
    """
    h = build_hamiltonian(k, gamma, nmax)
    evals, evecs = eigh(h, check_finite=False)
    n_index = 2 * k  # n=4K, while the even basis is 0,2,4,...

    weights = np.exp(-S * evals)
    v0 = evecs[0, :]
    vn = evecs[n_index, :]
    b = float(np.dot(v0 * vn, weights))
    m = float(np.dot(vn * vn, weights * weights))

    if b == 0.0 or m <= 0.0:
        return -math.inf, -math.inf
    return 2.0 * math.log(abs(b)), math.log(m) - k * DELTA


def perturbative_log_weighted_ratio(k: int, gamma: float) -> float:
    y = 4.0 * OMEGA * math.exp(-4.0 * OMEGA * S)
    y /= (gamma / k) * (1.0 - math.exp(-4.0 * OMEGA * S))

    terms = []
    log_y = math.log(y)
    log_falling = 0.0
    for ell in range(k + 1):
        if ell > 0:
            log_falling += math.log(k - ell + 1)
        terms.append(2.0 * log_falling + 2.0 * ell * log_y - math.lgamma(4 * ell + 1))

    max_term = max(terms)
    log_total = max_term + math.log(sum(math.exp(t - max_term) for t in terms))
    return log_total - k * DELTA


def large_k_gamma_star() -> tuple[float, float, float]:
    q = brentq(lambda x: 1.0 - x - math.log(x) - DELTA / 2.0, 1.0e-15, 1.0)
    c_star = 16.0 * (1.0 - q) ** 2 / q
    gamma_star = 4.0 * OMEGA * math.exp(-4.0 * OMEGA * S)
    gamma_star /= c_star * (1.0 - math.exp(-4.0 * OMEGA * S))
    return q, c_star, gamma_star


def first_crossing(func: Callable[[float], float]) -> float | None:
    """Find the first positive-to-negative crossing as gamma increases."""
    grid = np.geomspace(GAMMA_MIN, GAMMA_MAX, SCAN_POINTS)
    previous_gamma = float(grid[0])
    previous_value = func(previous_gamma)

    for gamma_np in grid[1:]:
        gamma = float(gamma_np)
        value = func(gamma)
        if previous_value > 0.0 and value <= 0.0:
            return float(brentq(func, previous_gamma, gamma, xtol=1.0e-12, rtol=1.0e-12))
        previous_gamma = gamma
        previous_value = value
    return None


def write_csv(path: Path, rows: list[dict[str, float | int | None]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_crossing_figure(
    summary_rows: list[dict[str, float | int | None]],
    gamma_large: float,
    out_path: Path,
) -> None:
    final_rows = [r for r in summary_rows if r["nmax"] == FINAL_NMAX]
    k = np.array([float(r["K"]) for r in final_rows])
    exact = np.array([float(r["exact_gamma_crossing"]) for r in final_rows])
    pert = np.array([float(r["perturbative_gamma_crossing"]) for r in final_rows])

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.plot(k, exact, "o-", color="#4477aa", label=fr"exact cutoff, $n_{{\max}}={FINAL_NMAX}$")
    ax.plot(k, pert, "s-", color="#cc6677", label="finite-$K$ perturbative")
    ax.axhline(gamma_large, color="black", linestyle="--", linewidth=1.4, label="large-$K$")
    ax.set_xlabel("$K$")
    ax.set_ylabel(r"crossover $\gamma_*$")
    ax.set_xlim(min(K_VALUES) - 0.4, max(K_VALUES) + 0.4)
    ax.set_ylim(0.0, max(0.12, 1.08 * float(np.nanmax(exact))))
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(out_path)
    plt.close(fig)


def write_convergence_figure(
    summary_rows: list[dict[str, float | int | None]],
    out_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    for nmax in NMAX_VALUES:
        rows = [r for r in summary_rows if r["nmax"] == nmax]
        k = [float(r["K"]) for r in rows]
        exact = [float(r["exact_gamma_crossing"]) for r in rows]
        ax.plot(k, exact, "o-", linewidth=1.2, markersize=3.5, label=fr"$n_{{\max}}={nmax}$")
    ax.set_xlabel("$K$")
    ax.set_ylabel(r"exact cutoff $\gamma_*$")
    ax.set_xlim(min(K_VALUES) - 0.4, max(K_VALUES) + 0.4)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.savefig(out_path)
    plt.close(fig)


def write_k10_options(
    rows: list[dict[str, float | int | None]],
    out_path: Path,
) -> None:
    gamma = np.array([float(r["gamma"]) for r in rows])
    log_a1 = np.array([float(r["log_A1_over_factorial"]) for r in rows])
    log_a2 = np.array([float(r["log_A2_weighted_over_factorial"]) for r in rows])
    ratio = log_a2 - log_a1

    fig, (ax, ax_ratio) = plt.subplots(
        2,
        1,
        figsize=(7.2, 5.4),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )
    ax.plot(gamma, log_a1 / OPTION_PLOT_K, color="#4477aa", label=r"$A_1$")
    ax.plot(gamma, log_a2 / OPTION_PLOT_K, color="#cc6677", label=r"$A_2 e^{-K\delta}$")
    ax.set_xscale("log")
    ax.set_ylabel(r"$K^{-1}\log(\mathrm{weight}/(4K)!)$")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.set_title(fr"Exact cutoff weights, $K={OPTION_PLOT_K}$, $n_{{\max}}={FINAL_NMAX}$")

    ax_ratio.axhline(0.0, color="black", linewidth=1.0, linestyle="--")
    ax_ratio.plot(gamma, ratio / OPTION_PLOT_K, color="#228833")
    ax_ratio.set_xscale("log")
    ax_ratio.set_xlabel(r"$\gamma$")
    ax_ratio.set_ylabel(r"$K^{-1}\log\frac{A_2e^{-K\delta}}{A_1}$")
    ax_ratio.grid(True, alpha=0.25)
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parent
    out_dir = here / "outputs"
    fig_dir = here / "figures"
    out_dir.mkdir(exist_ok=True)
    fig_dir.mkdir(exist_ok=True)

    q_star, c_star, gamma_large = large_k_gamma_star()

    summary_rows: list[dict[str, float | int | None]] = []
    for k in K_VALUES:
        perturbative = first_crossing(lambda gamma, kk=k: perturbative_log_weighted_ratio(kk, gamma))
        for nmax in NMAX_VALUES:
            exact = first_crossing(lambda gamma, kk=k, nn=nmax: exact_log_weighted_ratio(kk, gamma, nn))
            summary_rows.append(
                {
                    "K": k,
                    "nmax": nmax,
                    "exact_gamma_crossing": exact,
                    "perturbative_gamma_crossing": perturbative,
                    "large_K_gamma_crossing": gamma_large,
                }
            )
            print(f"K={k:2d}, nmax={nmax:3d}: exact={exact:.10g}, pert={perturbative:.10g}")

    write_csv(out_dir / "summary.csv", summary_rows)

    curve_rows: list[dict[str, float | int | None]] = []
    curve_grid = np.geomspace(GAMMA_MIN, 1.5, 100)
    for k in K_VALUES:
        for gamma in curve_grid:
            curve_rows.append(
                {
                    "K": k,
                    "gamma": float(gamma),
                    "exact_log_ratio": exact_log_weighted_ratio(k, float(gamma), FINAL_NMAX),
                    "perturbative_log_ratio": perturbative_log_weighted_ratio(k, float(gamma)),
                }
            )
    write_csv(out_dir / "curves.csv", curve_rows)

    option_rows: list[dict[str, float | int | None]] = []
    option_grid = np.geomspace(1.0e-3, 0.3, 240)
    for gamma in option_grid:
        log_a1, log_a2_weighted = exact_log_options(
            OPTION_PLOT_K, float(gamma), FINAL_NMAX
        )
        option_rows.append(
            {
                "K": OPTION_PLOT_K,
                "nmax": FINAL_NMAX,
                "gamma": float(gamma),
                "log_A1_over_factorial": log_a1,
                "log_A2_weighted_over_factorial": log_a2_weighted,
                "log_weighted_ratio": log_a2_weighted - log_a1,
            }
        )
    write_csv(out_dir / "k10_options.csv", option_rows)

    write_csv(
        out_dir / "large_k_prediction.csv",
        [{"q_star": q_star, "c_star": c_star, "gamma_star": gamma_large}],
    )

    write_crossing_figure(summary_rows, gamma_large, fig_dir / "finite_k_crossings.svg")
    write_convergence_figure(summary_rows, fig_dir / "cutoff_convergence.svg")
    write_k10_options(option_rows, fig_dir / "k10_options_vs_gamma.svg")


if __name__ == "__main__":
    main()
