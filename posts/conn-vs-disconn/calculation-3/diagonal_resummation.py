#!/usr/bin/env python3
r"""Diagonal-resummed perturbation theory for the toy oscillator.

This calculation treats the diagonal part of ``(a+a^\dagger)^4`` exactly and
keeps the leading monotone ``a^4`` lowering path into each ``4``-spaced Fock
sector.  It is a simple analytic improvement over the bare perturbative formula
in ``cvsd.md`` because the diagonal quartic self-energy is ``O(K)`` when
``n=O(K)`` and ``g=gamma/K``.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from mpmath import mp
from scipy.optimize import brentq


OMEGA = 1.0
S = 1.0
DELTA = 1.0
MP_DPS = 80

K_VALUES = list(range(1, 41)) + [50, 60, 80]
GAMMA_MIN = 1.0e-4
GAMMA_MAX = 1.0


def diagonal_energies(k: int, gamma: float) -> list[mp.mpf]:
    """Diagonal energies of |4q>, q=0,...,K, with g=gamma/K."""
    g = mp.mpf(gamma) / k
    return [
        4 * OMEGA * q + g * (96 * q * q + 24 * q + 3)
        for q in range(k + 1)
    ]


def consecutive_divided_differences(energies: list[mp.mpf]) -> list[list[mp.mpf]]:
    """Divided-difference table for f(E)=exp(-s E) on consecutive intervals."""
    cols: list[list[mp.mpf]] = [[mp.e ** (-S * e) for e in energies]]
    k = len(energies) - 1
    for order in range(1, k + 1):
        prev = cols[-1]
        col = []
        for start in range(k + 1 - order):
            denom = energies[start + order] - energies[start]
            col.append((prev[start + 1] - prev[start]) / denom)
        cols.append(col)
    return cols


def log_sector_amplitudes(k: int, gamma: float) -> list[float]:
    r"""Return logs of the common-factor-stripped sector amplitudes.

    The initial state ``(a^\dagger)^(4K)|0>`` contributes a common factor
    ``sqrt((4K)!)`` to every sector amplitude.  This common factor cancels in
    ``A2/A1`` and is omitted here.

    The amplitude into ``|4(K-j)>`` is approximated by

        (gamma/K)^j sqrt((4K)!/(4K-4j)!) [E_{K-j},...,E_K] exp(-s x),

    where the bracket is a divided difference.
    """
    mp.dps = MP_DPS
    energies = diagonal_energies(k, gamma)
    dd = consecutive_divided_differences(energies)
    g = mp.mpf(gamma) / k

    logs = []
    for j in range(k + 1):
        value = dd[j][k - j]
        lowering_log = mp.mpf("0.5") * (
            mp.loggamma(4 * k + 1) - mp.loggamma(4 * (k - j) + 1)
        )
        if j:
            lowering_log += j * mp.log(abs(g))
        logs.append(float(lowering_log + mp.log(abs(value))))
    return logs


def logsumexp(values: Iterable[float]) -> float:
    values = list(values)
    max_value = max(values)
    return max_value + math.log(sum(math.exp(v - max_value) for v in values))


def diagonal_resummed_log_ratio(k: int, gamma: float) -> float:
    """Return log(A2 exp(-K delta)/A1) in the diagonal-resummed approximation."""
    logs = log_sector_amplitudes(k, gamma)
    log_a2 = logsumexp(2.0 * x for x in logs) - k * DELTA
    log_a1 = 2.0 * logs[-1]
    return log_a2 - log_a1


def diagonal_resummed_options(k: int, gamma: float) -> tuple[float, float]:
    """Return common-factor-stripped log(A1), log(A2 exp(-K delta))."""
    logs = log_sector_amplitudes(k, gamma)
    return 2.0 * logs[-1], logsumexp(2.0 * x for x in logs) - k * DELTA


def diagonal_resummed_crossing(k: int) -> float | None:
    grid = np.geomspace(GAMMA_MIN, GAMMA_MAX, 80)
    previous_gamma = float(grid[0])
    previous_value = diagonal_resummed_log_ratio(k, previous_gamma)
    for gamma_np in grid[1:]:
        gamma = float(gamma_np)
        value = diagonal_resummed_log_ratio(k, gamma)
        if previous_value > 0.0 and value <= 0.0:
            return float(
                brentq(
                    lambda x, kk=k: diagonal_resummed_log_ratio(kk, x),
                    previous_gamma,
                    gamma,
                    xtol=1.0e-11,
                    rtol=1.0e-11,
                )
            )
        previous_gamma = gamma
        previous_value = value
    return None


def bare_perturbative_crossing(k: int) -> float:
    def log_ratio(gamma: float) -> float:
        y = 4.0 * OMEGA * math.exp(-4.0 * OMEGA * S)
        y /= (gamma / k) * (1.0 - math.exp(-4.0 * OMEGA * S))
        terms = []
        log_falling = 0.0
        log_y = math.log(y)
        for ell in range(k + 1):
            if ell > 0:
                log_falling += math.log(k - ell + 1)
            terms.append(2.0 * log_falling + 2.0 * ell * log_y - math.lgamma(4 * ell + 1))
        return logsumexp(terms) - k * DELTA

    return float(brentq(log_ratio, GAMMA_MIN, GAMMA_MAX))


def large_k_gamma_star() -> float:
    q = brentq(lambda x: 1.0 - x - math.log(x) - DELTA / 2.0, 1.0e-15, 1.0)
    c_star = 16.0 * (1.0 - q) ** 2 / q
    return 4.0 * OMEGA * math.exp(-4.0 * OMEGA * S) / (
        c_star * (1.0 - math.exp(-4.0 * OMEGA * S))
    )


def read_refined_summary(path: Path) -> dict[int, dict[str, float]]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as f:
        rows = csv.DictReader(f)
        return {
            int(row["K"]): {
                key: float(value)
                for key, value in row.items()
                if key != "K" and value not in ("", "None")
            }
            for row in rows
        }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_figures(out_dir: Path, summary: list[dict], comparison: list[dict]) -> None:
    k = np.array([r["K"] for r in summary], dtype=float)
    diag = np.array([r["diagonal_resummed_gamma_star"] for r in summary], dtype=float)
    bare = np.array([r["bare_perturbative_gamma_star"] for r in summary], dtype=float)
    gamma_inf = summary[0]["large_K_gamma_star"]

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.plot(k, diag, "o-", color="#228833", label="diagonal-resummed")
    ax.plot(k, bare, "s--", color="#cc6677", markersize=3, label="bare finite-$K$")
    ax.axhline(gamma_inf, color="black", linestyle="--", linewidth=1.2, label="bare large-$K$")
    if comparison:
        kc = np.array([r["K"] for r in comparison], dtype=float)
        env = np.array([r["refined_envelope_gamma_star"] for r in comparison], dtype=float)
        ax.plot(kc, env, "^-", color="#4477aa", label="refined envelope")
    ax.set_xlabel("$K$")
    ax.set_ylabel(r"crossover $\gamma_*$")
    ax.set_ylim(0.0, 0.075)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(out_dir / "diagonal_resummed_vs_K.svg")
    plt.close(fig)

    k0 = 15
    grid = np.geomspace(0.006, 0.09, 120)
    rows = []
    for gamma in grid:
        log_a1, log_a2w = diagonal_resummed_options(k0, float(gamma))
        rows.append((float(gamma), log_a1 / k0, log_a2w / k0))

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.plot([r[0] for r in rows], [r[1] for r in rows], color="#228833", label=r"$A_1$ diag-resummed")
    ax.plot([r[0] for r in rows], [r[2] for r in rows], color="#cc6677", label=r"$A_2 e^{-K\delta}$ diag-resummed")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$K^{-1}\log(\mathrm{weight}/(4K)!)$")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.set_title(r"Diagonal-resummed weights at $K=15$")
    fig.savefig(out_dir / "diagonal_resummed_K15_rates.svg")
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parent
    output_dir = here / "outputs"
    figure_dir = here / "figures"
    output_dir.mkdir(exist_ok=True)
    figure_dir.mkdir(exist_ok=True)

    refined = read_refined_summary(here.parent / "calculation-2" / "outputs" / "crossover_summary.csv")
    gamma_inf = large_k_gamma_star()

    summary = []
    for k in K_VALUES:
        diag = diagonal_resummed_crossing(k)
        bare = bare_perturbative_crossing(k)
        summary.append(
            {
                "K": k,
                "diagonal_resummed_gamma_star": diag,
                "bare_perturbative_gamma_star": bare,
                "large_K_gamma_star": gamma_inf,
            }
        )
        print(f"K={k:2d}: diag={diag:.10g}, bare={bare:.10g}")

    comparison = []
    for row in summary:
        k = row["K"]
        if k in refined:
            comparison.append(
                {
                    "K": k,
                    "diagonal_resummed_gamma_star": row["diagonal_resummed_gamma_star"],
                    "bare_perturbative_gamma_star": row["bare_perturbative_gamma_star"],
                    "refined_envelope_gamma_star": refined[k]["envelope_gamma_star"],
                    "refined_window_gamma_star": refined[k]["window_gamma_star"],
                    "large_K_gamma_star": gamma_inf,
                }
            )

    write_csv(output_dir / "diagonal_resummed_summary.csv", summary)
    if comparison:
        write_csv(output_dir / "comparison_to_refined.csv", comparison)
    make_figures(figure_dir, summary, comparison)


if __name__ == "__main__":
    main()
