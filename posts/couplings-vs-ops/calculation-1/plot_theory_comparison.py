#!/usr/bin/env python3
"""Plot ED/theory comparison for the operator-ensemble calculation."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-couplings-vs-ops")

import matplotlib.pyplot as plt
import numpy as np


def load_rows(path: Path) -> list[dict[str, float]]:
    with path.open() as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            rows.append(
                {
                    "N": float(row["n_majorana"]),
                    "samples": float(row["n_samples"]),
                    "m1_ratio": float(row["std_m1_over_prediction"]),
                    "m2_ratio": float(row["rel_std_m2_over_prediction"]),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("outputs_beta1_highstats/summary.csv"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("figures"))
    args = parser.parse_args()

    rows = sorted(load_rows(args.summary), key=lambda row: row["N"])
    n = np.array([row["N"] for row in rows])
    inv_n = 1.0 / n
    samples = np.array([row["samples"] for row in rows])
    m1 = np.array([row["m1_ratio"] for row in rows])
    m2 = np.array([row["m2_ratio"] for row in rows])
    # Delta-method uncertainty for a sample standard deviation ratio.
    m1_err = m1 / np.sqrt(2.0 * (samples - 1.0))
    m2_err = m2 / np.sqrt(2.0 * (samples - 1.0))

    args.outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axhline(1.0, color="0.25", linewidth=1.2, linestyle="--", label="ideal joint Gaussian")
    ax.errorbar(
        inv_n,
        m1,
        yerr=m1_err,
        marker="o",
        linewidth=1.5,
        capsize=3,
        label=r"$\mathrm{std}(M_1)$ / theory",
    )
    ax.errorbar(
        inv_n,
        m2,
        yerr=m2_err,
        marker="s",
        linewidth=1.5,
        capsize=3,
        label=r"$\mathrm{relstd}(M_2)$ / theory",
    )
    ax.set_xlabel(r"$1/N$")
    ax.set_ylabel("measured / theory")
    ax.set_title(r"Operator-ensemble self-averaging, $p=4$, $W=4$, $\beta=1$")
    ax.set_xticks(inv_n)
    ax.set_xticklabels([f"{x:.3f}" for x in inv_n])
    ax.set_ylim(0.85, max(1.7, float(np.max(m2 + m2_err)) + 0.05))
    ax.legend(frameon=False)
    ax.grid(True, linewidth=0.5, alpha=0.35)
    fig.tight_layout()

    png = args.outdir / "theory_comparison_vs_invN.png"
    svg = args.outdir / "theory_comparison_vs_invN.svg"
    fig.savefig(png, dpi=220)
    fig.savefig(svg)
    print(png)
    print(svg)


if __name__ == "__main__":
    main()
