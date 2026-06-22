#!/usr/bin/env python3
r"""Robust (zero-insensitive) crossover for the connected/disconnected toy model.

Background
----------
`calculation-1` locates the connected/disconnected crossover by the *first*
crossing of ``A2 exp(-K delta) / A1`` as ``gamma`` increases, at fixed ``s=1``.
For ``K >= 13`` that observable is not well defined:

  * The disconnected amplitude ``B(gamma) = <0|e^{-sH}|4K>`` (so that
    ``A1 = |B|^2``) is genuinely *not* sign definite.  Confirmed at 60-digit
    precision, ``B`` has real zeros at finite ``gamma`` (competing multi-step
    paths cancel), so ``log A1`` genuinely spikes to ``-inf`` on a discrete set
    and the "first crossing" is set by whichever zero comes first, not by the
    smooth competition of the two weights.
  * Worse, in *double* precision ``B`` cannot even be computed in the crossover
    region for ``K >= 12``: there ``|B|`` falls below ~1e-21 and the value
    returned is deterministic round-off noise (it reproduces across ``nmax`` yet
    is wrong in magnitude and sometimes in sign).  Raising ``nmax`` does not
    help because it does not raise the floating-point floor.

This script fixes both problems:

  * It computes ``A1`` and ``A2`` at ``s=1`` in *extended precision* (mpmath),
    so the small-``gamma`` weights are the true values rather than noise.
  * It replaces the fragile "first crossing" with a crossover defined from the
    *smooth envelope* of ``A1``.  Writing ``B = E(gamma) * oscillation``, the
    local maxima of ``|B|`` trace the smooth envelope ``E`` (the zeros are where
    the oscillation vanishes).  Interpolating ``log A1`` through its peaks and
    crossing it with the smooth ``log(A2 e^{-K delta})`` gives a single,
    zero-insensitive crossover ``gamma_*(K)``.

Cross-checks (same headline conclusion):
  * ``window``  -- a multiplicative running average of ``A1`` in ``gamma``
    (linear average; dominated by the peaks, hence the envelope).
  * ``s_avg``   -- average of ``A1(s)``, ``A2(s)`` over a window of ``s`` around
    the nominal ``s=1`` (double precision; the zeros of ``B(gamma,s)`` are
    isolated in ``s`` so the ``s``-integral is smooth).  NOTE: because
    ``A1(s), A2(s)`` grow as ``s`` decreases, this integral is weighted toward
    the small-``s`` end of the window and therefore reports a crossover at a
    somewhat *larger* ``gamma`` than the fixed ``s=1`` value; it is included
    only to confirm the qualitative trend (no downward drift), not the absolute
    location.
  * ``first``   -- the original `calculation-1` double-precision first crossing,
    reproduced here purely for contrast.

All runs use omega = s = delta = 1.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
from mpmath import mp, mpf, matrix, eigsy, sqrt, exp, log as mlog
from scipy.linalg import eigh
from scipy.integrate import simpson
from scipy.optimize import brentq

OMEGA = 1.0
S = 1.0
DELTA = 1.0
MP_DPS = 40

K_VALUES = list(range(10, 19))
NMAX_MP = 90                       # converged for the extended-precision s=1 scan
# cutoff-convergence spot checks: (K, [nmax values]).  K=16,18 added to confirm
# that nmax=90 is adequate even when |4K> is close to the cutoff.
CONV_SPEC = [(14, [80, 120, 160]), (16, [90, 140]), (18, [90, 140])]
DIAG_K = 15                        # K used for the method-illustration figure

GAMMA_LO = 4.0e-3
GAMMA_HI = 0.13
GAMMA_PTS = 90

WINDOW_FACTOR = 1.18               # multiplicative half-width for the gamma window
S_LO, S_HI, S_PTS = 0.7, 1.3, 41   # s-average window


# --------------------------------------------------------------------------
# Extended-precision weights at s = 1
# --------------------------------------------------------------------------
def _mp_eigsystem(k: int, gamma: float, nmax: int):
    """Even-sector H = omega a^dag a + (gamma/K)(a+a^dag)^4; return (E, Q, dim)."""
    mp.dps = MP_DPS
    g = mpf(gamma) / k
    ms = list(range(0, nmax + 1, 2))
    dim = len(ms)
    h = matrix(dim, dim)
    for i, mi in enumerate(ms):
        m = mpf(mi)
        h[i, i] = OMEGA * m + g * (6 * m * m + 6 * m + 3)
    for i in range(dim - 1):
        m = mpf(ms[i])
        v = g * (4 * m + 6) * sqrt((m + 1) * (m + 2))
        h[i, i + 1] = v
        h[i + 1, i] = v
    for i in range(dim - 2):
        m = mpf(ms[i])
        v = g * sqrt((m + 1) * (m + 2) * (m + 3) * (m + 4))
        h[i, i + 2] = v
        h[i + 2, i] = v
    E, Q = eigsy(h)
    return E, Q, dim


def mp_log_weights(k: int, gamma: float, nmax: int) -> tuple[float, float]:
    r"""Return (log A1, log(A2 e^{-K delta})) at s=1, common factor (4K)! stripped.

    B  = <0|e^{-sH}|4K>,  A1 = B^2;  A2 = <4K|e^{-2sH}|4K>.
    """
    E, Q, dim = _mp_eigsystem(k, gamma, nmax)
    n = 2 * k  # even-basis index of |4K>
    b = mp.fsum(Q[0, a] * Q[n, a] * exp(-S * E[a]) for a in range(dim))
    a2 = mp.fsum(Q[n, a] ** 2 * exp(-2 * S * E[a]) for a in range(dim))
    log_a1 = 2.0 * float(mlog(abs(b))) if b != 0 else -1.0e18
    log_a2w = float(mlog(a2)) - k * DELTA
    return log_a1, log_a2w


def scan_grid(k: int, nmax: int):
    """Return (gamma, log_gamma, logA1, logA2w) arrays on the working grid."""
    g = np.geomspace(GAMMA_LO, GAMMA_HI, GAMMA_PTS)
    lg = np.log(g)
    l1 = np.empty(GAMMA_PTS)
    l2 = np.empty(GAMMA_PTS)
    for i, gg in enumerate(g):
        l1[i], l2[i] = mp_log_weights(k, float(gg), nmax)
    return g, lg, l1, l2


# --------------------------------------------------------------------------
# Crossover definitions
# --------------------------------------------------------------------------
def _crossings_from_below(lg, diff):
    """gamma values where ``diff`` (= disconnected - connected) goes - -> +."""
    out = []
    for i in range(len(lg) - 1):
        if diff[i] < 0.0 and diff[i + 1] >= 0.0:
            t = -diff[i] / (diff[i + 1] - diff[i])
            out.append(math.exp(lg[i] + t * (lg[i + 1] - lg[i])))
    return out


def envelope_curve(lg, l1):
    """Upper envelope of log A1: interpolate through interior local maxima."""
    pk = [i for i in range(1, len(l1) - 1) if l1[i] >= l1[i - 1] and l1[i] >= l1[i + 1]]
    if len(pk) >= 2:
        return np.interp(lg, lg[pk], l1[pk]), pk
    return l1.copy(), pk


def envelope_crossover(lg, l1, l2):
    env, pk = envelope_curve(lg, l1)
    cr = _crossings_from_below(lg, env - l2)
    return (cr[-1] if cr else None), env, pk


def window_crossover(g, lg, l1, l2):
    """Crossover from a multiplicative running average of A1 (linear space)."""
    a1 = np.exp(l1 - l1.max())  # shift for numerical safety; constant cancels in log
    avg = np.empty_like(a1)
    for i in range(len(g)):
        sel = (g >= g[i] / WINDOW_FACTOR) & (g <= g[i] * WINDOW_FACTOR)
        avg[i] = a1[sel].mean()
    log_avg = np.log(avg) + l1.max()
    cr = _crossings_from_below(lg, log_avg - l2)
    return cr[-1] if cr else None


# --------------------------------------------------------------------------
# Double-precision cross-checks
# --------------------------------------------------------------------------
def _np_eigsystem(k: int, gamma: float, nmax: int):
    m = np.arange(0, nmax + 1, 2, dtype=float)
    g = gamma / k
    h = np.diag(OMEGA * m + g * (6 * m * m + 6 * m + 3))
    m1 = m[:-1]
    o1 = g * (4 * m1 + 6) * np.sqrt((m1 + 1) * (m1 + 2))
    idx = np.arange(len(m) - 1)
    h[idx, idx + 1] = o1
    h[idx + 1, idx] = o1
    m2 = m[:-2]
    o2 = g * np.sqrt((m2 + 1) * (m2 + 2) * (m2 + 3) * (m2 + 4))
    idx = np.arange(len(m) - 2)
    h[idx, idx + 2] = o2
    h[idx + 2, idx] = o2
    ev, U = eigh(h, check_finite=False)
    return ev, U[0, :], U[2 * k, :]


def s_average_crossover(k: int, nmax: int) -> float | None:
    sgrid = np.linspace(S_LO, S_HI, S_PTS)

    def diff(gamma: float) -> float:
        ev, v0, vn = _np_eigsystem(k, gamma, nmax)
        a1 = np.array([np.dot(v0 * vn, np.exp(-s * ev)) ** 2 for s in sgrid])
        a2 = np.array([np.dot(vn * vn, np.exp(-2 * s * ev)) for s in sgrid])
        return (math.log(simpson(a2, x=sgrid)) - k * DELTA) - math.log(simpson(a1, x=sgrid))

    lo, hi = 5.0e-3, 0.3
    if diff(lo) <= 0 or diff(hi) >= 0:
        return None
    return float(brentq(diff, lo, hi, xtol=1e-10, rtol=1e-10))


def first_crossover(k: int, nmax: int) -> float | None:
    """Original calculation-1 double-precision first crossing, for contrast."""
    def ratio(gamma: float) -> float:
        ev, v0, vn = _np_eigsystem(k, gamma, nmax)
        w = np.exp(-S * ev)
        b = float(np.dot(v0 * vn, w))
        a2 = float(np.dot(vn * vn, w * w))
        if b == 0.0 or a2 <= 0.0:
            return -math.inf
        return (math.log(a2) - k * DELTA) - 2.0 * math.log(abs(b))

    grid = np.geomspace(1.0e-4, 3.0, 121)
    prev_g, prev_v = float(grid[0]), ratio(float(grid[0]))
    for gn in grid[1:]:
        g = float(gn)
        v = ratio(g)
        if prev_v > 0.0 and v <= 0.0:
            return float(brentq(ratio, prev_g, g, xtol=1e-12, rtol=1e-12))
        prev_g, prev_v = g, v
    return None


def large_k_gamma_star() -> float:
    q = brentq(lambda x: 1.0 - x - math.log(x) - DELTA / 2.0, 1e-15, 1.0)
    c_star = 16.0 * (1.0 - q) ** 2 / q
    return 4.0 * OMEGA * math.exp(-4.0 * OMEGA * S) / (c_star * (1.0 - math.exp(-4.0 * OMEGA * S)))


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    here = Path(__file__).resolve().parent
    out = here / "outputs"
    fig = here / "figures"
    out.mkdir(exist_ok=True)
    fig.mkdir(exist_ok=True)

    gamma_inf = large_k_gamma_star()
    print(f"large-K gamma_* = {gamma_inf:.6f}\n")

    summary = []
    diag = None
    for k in K_VALUES:
        g, lg, l1, l2 = scan_grid(k, NMAX_MP)
        env_star, env, pk = envelope_crossover(lg, l1, l2)
        win_star = window_crossover(g, lg, l1, l2)
        s_star = s_average_crossover(k, 240)
        first_star = first_crossover(k, 240)
        summary.append(
            {
                "K": k,
                "envelope_gamma_star": env_star,
                "window_gamma_star": win_star,
                "s_average_gamma_star": s_star,
                "first_crossing_gamma_star": first_star,
                "large_K_gamma_star": gamma_inf,
            }
        )
        print(
            f"K={k:2d}: envelope={env_star!s:>9.9} window={win_star!s:>9.9} "
            f"s_avg={s_star!s:>9.9} first={first_star!s:>9.9}"
        )
        if k == DIAG_K:
            diag = (g, l1, l2, env, pk, env_star)

    write_csv(out / "crossover_summary.csv", summary)

    # cutoff-convergence spot checks (extended precision)
    conv = []
    for conv_k, nmaxes in CONV_SPEC:
        for nmax in nmaxes:
            g, lg, l1, l2 = scan_grid(conv_k, nmax)
            env_star, _, _ = envelope_crossover(lg, l1, l2)
            conv.append({"K": conv_k, "nmax": nmax, "envelope_gamma_star": env_star})
            print(f"[conv] K={conv_k} nmax={nmax}: envelope={env_star}")
    write_csv(out / "convergence_check.csv", conv)

    # diagnostic data for the illustration figure
    g, l1, l2, env, pk, env_star = diag
    write_csv(
        out / f"diagnostic_K{DIAG_K}.csv",
        [
            {"gamma": float(gg), "logA1": float(a), "logA2w": float(b), "envelope": float(e)}
            for gg, a, b, e in zip(g, l1, l2, env)
        ],
    )

    _figures(fig, summary, gamma_inf, (g, l1, l2, env, pk, env_star))
    print("\ndone")


def _figures(fig: Path, summary: list[dict], gamma_inf: float, diag) -> None:
    import matplotlib.pyplot as plt

    ks = np.array([r["K"] for r in summary], dtype=float)

    def col(name):
        return np.array([np.nan if r[name] is None else r[name] for r in summary])

    f, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.plot(ks, col("envelope_gamma_star"), "o-", color="#228833",
            label="envelope (extended precision, $s{=}1$)")
    ax.plot(ks, col("window_gamma_star"), "^--", color="#66ccee", lw=1.1,
            label=r"$\gamma$-window average")
    ax.plot(ks, col("s_average_gamma_star"), "v:", color="#aa3377", lw=1.1,
            label="$s$-average (weights $s{<}1$)")
    ax.plot(ks, col("first_crossing_gamma_star"), "s-", color="#cc6677", lw=1.0,
            label="first crossing (calc-1, double prec.)")
    ax.axhline(gamma_inf, color="black", ls="--", lw=1.3, label="large-$K$")
    ax.set_xlabel("$K$")
    ax.set_ylabel(r"crossover $\gamma_*$")
    ax.set_ylim(0.0, 0.12)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    f.savefig(fig / "robust_crossover_vs_K.svg")
    plt.close(f)

    g, l1, l2, env, pk, env_star = diag
    f, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    k = DIAG_K
    ax.plot(g, l1 / k, color="#4477aa", lw=1.0, alpha=0.7, label=r"$K^{-1}\log A_1$ (exact)")
    ax.plot(g, env / k, color="#228833", lw=1.8, label=r"$K^{-1}\log A_1$ envelope")
    ax.plot(g, l2 / k, color="#cc6677", lw=1.8, label=r"$K^{-1}\log(A_2 e^{-K\delta})$")
    if len(pk):
        ax.plot(g[pk], l1[pk] / k, "o", ms=3, color="#228833")
    if env_star is not None:
        ax.axvline(env_star, color="black", ls=":", lw=1.2)
        ax.text(env_star, ax.get_ylim()[0], fr"  $\gamma_*\approx{env_star:.3f}$", va="bottom")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel("rate")
    ax.set_title(fr"Envelope crossover, $K={DIAG_K}$, $n_{{\max}}={NMAX_MP}$ (extended precision)")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    f.savefig(fig / f"diagnostic_K{DIAG_K}.svg")
    plt.close(f)


if __name__ == "__main__":
    main()
