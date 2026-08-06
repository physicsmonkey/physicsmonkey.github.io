#!/usr/bin/env python3
"""Make direct visual diagnostics for the N=24 G_11 and G_12 tests.

Each row shows the same comparison in three complementary ways: an ED-versus-
saddle parity plot, a signal-scaled difference surface, and the difference in
units of the pointwise jackknife error.  The latter is descriptive rather than
a chi-squared statistic because the surface entries are strongly correlated.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


@dataclass
class Case:
    label: str
    prediction: np.ndarray
    ed: np.ndarray
    error: np.ndarray
    mask: np.ndarray


def metrics(case: Case) -> tuple[float, float, float]:
    difference = case.ed - case.prediction
    relative_residual = float(
        np.linalg.norm(difference[case.mask])
        / np.linalg.norm(case.ed[case.mask])
    )
    error_norm = float(np.linalg.norm(case.error))
    difference_to_error = (
        float(np.linalg.norm(difference) / error_norm)
        if error_norm > 0
        else float("nan")
    )
    max_signal_scaled_difference = float(
        np.max(np.abs(difference)) / np.max(np.abs(case.prediction))
    )
    return relative_residual, difference_to_error, max_signal_scaled_difference


def make_figure(path: Path, fractions: np.ndarray, cases: list[Case], beta_j: float) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(
        len(cases),
        3,
        figsize=(11.8, 2.85 * len(cases)),
        constrained_layout=True,
    )
    extent = [
        float(fractions[0]),
        float(fractions[-1]),
        float(fractions[0]),
        float(fractions[-1]),
    ]
    for row, case in enumerate(cases):
        residual, difference_to_error, max_scaled_difference = metrics(case)
        prediction = case.prediction[case.mask]
        ed = case.ed[case.mask]
        error = case.error[case.mask]
        amplitude = float(
            max(np.max(np.abs(case.prediction)), np.max(np.abs(case.ed)))
        )

        parity_axis = axes[row, 0]
        parity_axis.fill_between(
            [-amplitude, amplitude],
            [-amplitude - 0.01 * amplitude, amplitude - 0.01 * amplitude],
            [-amplitude + 0.01 * amplitude, amplitude + 0.01 * amplitude],
            color="0.90",
            label=r"$\pm1\%$ of peak",
        )
        parity_axis.errorbar(
            prediction,
            ed,
            yerr=error,
            fmt="o",
            markersize=2.6,
            elinewidth=0.55,
            capsize=0,
            alpha=0.68,
            color="#1769aa",
        )
        parity_axis.plot(
            [-amplitude, amplitude],
            [-amplitude, amplitude],
            color="black",
            linewidth=0.9,
        )
        combined = np.concatenate([prediction, ed])
        padding = 0.05 * max(float(np.ptp(combined)), amplitude)
        limits = [float(combined.min() - padding), float(combined.max() + padding)]
        parity_axis.set_xlim(limits)
        parity_axis.set_ylim(limits)
        parity_axis.set_aspect("equal", adjustable="box")
        parity_axis.set_xlabel("saddle prediction")
        parity_axis.set_ylabel("ED")
        parity_axis.set_title(case.label, fontsize=10)
        parity_axis.text(
            0.03,
            0.97,
            rf"residual $={100 * residual:.3g}\%$" + "\n" + rf"$\|\Delta\|/\|\mathrm{{SE}}\|={difference_to_error:.2g}$",
            transform=parity_axis.transAxes,
            va="top",
            fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82},
        )

        difference = case.ed - case.prediction
        scaled_difference = 100.0 * difference / amplitude
        difference_limit = max(
            float(np.max(np.abs(scaled_difference))), 1e-12
        )
        image = axes[row, 1].imshow(
            scaled_difference,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap="RdBu_r",
            vmin=-difference_limit,
            vmax=difference_limit,
        )
        axes[row, 1].set_title(
            rf"$(\mathrm{{ED}}-\mathrm{{saddle}})/\max|\mathrm{{saddle}}|$ (\%)"
            + "\n"
            + rf"max $={100 * max_scaled_difference:.3g}\%$",
            fontsize=9,
        )
        figure.colorbar(image, ax=axes[row, 1], fraction=0.046, pad=0.04)

        pulls = np.full(difference.shape, np.nan)
        reliable = case.error > 10 * np.finfo(float).eps * amplitude
        pulls[reliable] = difference[reliable] / case.error[reliable]
        image = axes[row, 2].imshow(
            pulls,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap="RdBu_r",
            vmin=-5.0,
            vmax=5.0,
        )
        axes[row, 2].set_title(
            r"$(\mathrm{ED}-\mathrm{saddle})/\mathrm{SE}$ (clipped at $\pm5$)",
            fontsize=9,
        )
        figure.colorbar(image, ax=axes[row, 2], fraction=0.046, pad=0.04)
        for axis in axes[row, 1:]:
            axis.set_xlabel(r"$\tau'/\beta$")
            axis.set_ylabel(r"$\tau/\beta$")

    figure.suptitle(
        rf"$N=24$, $W/N=1/3$, $\beta J={beta_j:g}$: parameter-free component tests",
        fontsize=13,
    )
    figure.savefig(path, dpi=190)
    plt.close(figure)


def run(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    with np.load(args.cross_data) as cross:
        fractions = cross["fractions"].copy()
        for beta, same_path in (
            (0.5, args.same_beta_0p5),
            (1.0, args.same_beta_1),
        ):
            beta_key = "0p5" if beta == 0.5 else "1"
            with np.load(same_path) as same:
                off_diagonal = ~np.eye(len(fractions), dtype=bool)
                away_from_collision = np.ones(
                    (len(fractions), len(fractions)), dtype=bool
                )
                away_from_collision[0, 0] = False
                cases = [
                    Case(
                        r"$G^{\rm P}_{11}$ (inside $W=8$)",
                        same["periodic_G11"],
                        same["inside_ratio"],
                        same["inside_standard_error"],
                        off_diagonal,
                    ),
                    Case(
                        r"$G^{\rm AP}_{11}$ (outside $W=8$)",
                        same["antiperiodic_G11"],
                        same["outside_ratio"],
                        same["outside_standard_error"],
                        off_diagonal,
                    ),
                    Case(
                        r"normalized $G^{\rm P}_{12}$ (inside $W=9$)",
                        cross[f"prediction_beta{beta_key}_inside_p_w0p333333"],
                        cross[f"ratio_beta{beta_key}_inside_p"],
                        cross[f"jackknife_error_beta{beta_key}_inside_p"],
                        away_from_collision,
                    ),
                    Case(
                        r"normalized $G^{\rm AP}_{12}$ (outside $W=7$)",
                        cross[f"prediction_beta{beta_key}_outside_ap_w0p333333"],
                        cross[f"ratio_beta{beta_key}_outside_ap"],
                        cross[f"jackknife_error_beta{beta_key}_outside_ap"],
                        away_from_collision,
                    ),
                ]
                make_figure(
                    args.outdir / f"component_agreement_beta{beta_key}.png",
                    fractions,
                    cases,
                    beta,
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--same-beta-0p5",
        type=Path,
        default=HERE / "outputs_16_beta0p5" / "large_same_side_data.npz",
    )
    parser.add_argument(
        "--same-beta-1",
        type=Path,
        default=HERE / "outputs_16_beta1" / "large_same_side_data.npz",
    )
    parser.add_argument(
        "--cross-data",
        type=Path,
        default=HERE / "cross_replica_outputs_16" / "large_cross_replica_data.npz",
    )
    parser.add_argument(
        "--outdir", type=Path, default=HERE / "component_diagnostics"
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
