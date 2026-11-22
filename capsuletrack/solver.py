"""6-DOF capsule pose estimation.

1. Global seed: evaluate the field model on a workspace grid; at each grid point the
   best orientation is solved in closed form (weighted orthogonal Procrustes / Kabsch),
   giving a residual map. The best few grid points seed step 2.
2. Local refinement: Levenberg-Marquardt on noise-whitened residuals over
   (position, rotation vector).
3. Cramer-Rao lower bound from the whitened Jacobian at the true pose, used to check
   estimator efficiency.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

from .system import Patch


def procrustes(Bw, bw):
    """R minimising ||Bw R - bw||_F over SO(3). Bw, bw: (N,3) weighted rows."""
    U, _, Vt = np.linalg.svd(Bw.T @ bw)
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    return U @ D @ Vt


def seed(patch: Patch, b, sigma, grid, model="loop", k=3):
    """Rotation-invariant seeding: match log field magnitudes |b_i| (one per coil) on a
    grid, then fit orientation in closed form at the best grid points. The log of the
    magnitude is far smoother near the patch than the raw vector residual."""
    Bg = patch.fields(grid, model=model)            # (G,N,3)
    lm = np.log(np.linalg.norm(b, axis=1))
    lg = np.log(np.linalg.norm(Bg, axis=2))
    costs = np.sum((lg - lm[None, :]) ** 2, axis=1)
    idx = np.argsort(costs)[:k]
    w = 1 / sigma
    out = []
    for i in idx:
        p0 = grid[i]
        # refine position on log-magnitudes before fitting orientation
        f = lambda x: np.log(np.linalg.norm(patch.fields(x, model=model)[0], axis=1)) - lm
        p0 = least_squares(f, p0, method="lm", x_scale=0.01, max_nfev=200).x
        R0 = procrustes(patch.fields(p0, model=model)[0] * w, b * w)
        out.append((p0, R0))
    return out


def _residual(x, patch, b, sigma, model):
    R = Rotation.from_rotvec(x[3:]).as_matrix()
    B = patch.fields(x[:3], model=model)[0]
    return ((B @ R - b) / sigma).ravel()


def solve(patch: Patch, b, sigma, grid, model="loop", k=3, grid_fine=None, chi2_gate=1000.0):
    """Seeded LM. If the normalised chi-square of the best fit exceeds chi2_gate
    (gross misfit -> wrong local minimum), re-seed from grid_fine with 4x more seeds."""
    p, R, best = _solve(patch, b, sigma, grid, model, k)
    dof = max(b.size - 6, 1)
    best.retried = False
    if grid_fine is not None and 2 * best.cost / dof > chi2_gate:
        p2, R2, best2 = _solve(patch, b, sigma, grid_fine, model, 4 * k)
        best2.retried = True
        if best2.cost < best.cost:
            return p2, R2, best2
    return p, R, best


def _solve(patch, b, sigma, grid, model, k):
    best = None
    for p0, R0 in seed(patch, b, sigma, grid, model=model, k=k):
        x0 = np.concatenate([p0, Rotation.from_matrix(R0).as_rotvec()])
        r = least_squares(_residual, x0, args=(patch, b, sigma, model), method="lm",
                          x_scale=np.array([0.01] * 3 + [0.1] * 3), max_nfev=400)
        if best is None or r.cost < best.cost:
            best = r
    p = best.x[:3]
    R = Rotation.from_rotvec(best.x[3:]).as_matrix()
    return p, R, best


def crlb(patch: Patch, p, R, sigma, h=1e-5):
    """Cramer-Rao bound: (position RMS std [m], orientation RMS std [rad])."""
    x = np.concatenate([p, Rotation.from_matrix(R).as_rotvec()])
    b0 = patch.fields(p)[0] @ R
    f0 = _residual(x, patch, b0, sigma, "loop")
    J = np.empty((f0.size, 6))
    for j in range(6):
        dx = np.zeros(6)
        dx[j] = h
        J[:, j] = (_residual(x + dx, patch, b0, sigma, "loop") - f0) / h
    cov = np.linalg.pinv(J.T @ J)
    return float(np.sqrt(np.trace(cov[:3, :3]))), float(np.sqrt(np.trace(cov[3:, 3:])))


def orientation_error_deg(R_est, R_true):
    return float(np.degrees(Rotation.from_matrix(R_est @ R_true.T).magnitude()))
