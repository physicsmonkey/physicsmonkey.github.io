#!/usr/bin/env python3
"""Scan the P/AP collision product over beta*J and w on the connected branch.

Tests:
  1. product = G^P_12(0,0) * G^AP_12(0,0) vs -1/4 (Richardson 2G_2L - G_L)
  2. zero-mode prediction  G^P_12(0,0) ?= 1 / (integral of Sigma_12)
  3. first-order prediction G^AP_12(0,0) ?= [G^AP_d Sigma_12 G^AP_d](0,0)
     where G^AP_d = -(D_AP - dt^2 Sigma_11)^{-T} is the AP diagonal-only
     propagator.
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np

CALC1 = (
    Path(__file__).resolve().parents[1]
    / "calculation-1"
    / "compare_ed_path_integral.py"
)
spec = importlib.util.spec_from_file_location("calc1", CALC1)
calc1 = importlib.util.module_from_spec(spec)
sys.modules["calc1"] = calc1
spec.loader.exec_module(calc1)


def conditional_collisions(beta, w, length, J=1.0, mixing=0.01, tol=1e-8, itmax=20000):
    saddle = calc1.solve_weighted_replica_saddle(
        beta=beta, relative_weight=w, length=length, coupling_j=J,
        mixing=mixing, tolerance=tol, max_iterations=itmax)
    if not saddle.converged:
        raise RuntimeError(f"no convergence beta={beta} w={w} L={length}")
    G = saddle.G_blocks
    sigma = J**2 * G**3
    sigma_full = calc1.replica_full_matrix(sigma)
    dt2 = (beta / length) ** 2
    P = calc1.replica_blocks(-np.linalg.inv(
        calc1.discrete_derivative(length, periodic=True) - dt2 * sigma_full).T)
    AP = calc1.replica_blocks(-np.linalg.inv(
        calc1.discrete_derivative(length, periodic=False) - dt2 * sigma_full).T)
    gp = 0.5 * (P[1] + P[1].T)[0, 0]
    gap = 0.5 * (AP[1] + AP[1].T)[0, 0]

    # zero-mode prediction for the P off-diagonal collision value
    dt = beta / length
    sig12_int = dt * dt * sigma[1].sum()

    # first-order-in-Sigma_12 prediction for the AP off-diagonal
    D_ap_single = np.eye(length)
    D_ap_single[np.arange(1, length), np.arange(length - 1)] = -1.0
    D_ap_single[0, -1] = 1.0
    Gd = -np.linalg.inv(D_ap_single - dt * dt * sigma[0]).T
    ap_first = (dt * dt) * (Gd @ sigma[1] @ Gd)
    ap_first_col = 0.5 * (ap_first + ap_first.T)[0, 0]
    return gp, gap, sig12_int, ap_first_col, saddle.iterations


rows = []
for betaJ in (0.25, 0.5, 1.0, 2.0):
    for w in (0.1, 0.2, 0.3):
        try:
            vals = {}
            for L in (90, 180):
                vals[L] = conditional_collisions(betaJ, w, L)
            gp = 2 * vals[180][0] - vals[90][0]
            gap = 2 * vals[180][1] - vals[90][1]
            s12 = 2 * vals[180][2] - vals[90][2]
            apf = 2 * vals[180][3] - vals[90][3]
            rows.append({
                "betaJ": betaJ, "w": w,
                "GP00": gp, "GAP00": gap,
                "product": gp * gap,
                "product_minus_quarter": gp * gap + 0.25,
                "zero_mode_pred_GP00": 1.0 / s12 if s12 else float("nan"),
                "first_order_pred_GAP00": apf,
                "iters": vals[180][4],
            })
            r = rows[-1]
            print(f"betaJ={betaJ:5.2f} w={w:4.2f}  GP00={gp:+.6f}  GAP00={gap:+.6f}  "
                  f"prod={gp*gap:+.6f}  dev={gp*gap+0.25:+.2e}  "
                  f"ZM_pred={r['zero_mode_pred_GP00']:+.6f}  AP1_pred={apf:+.6f}",
                  flush=True)
        except RuntimeError as err:
            print("SKIP", err, flush=True)

out = Path(__file__).parent / "outputs" / "quarter_product_scan.json"
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(rows, indent=2))
print("saved", out)
