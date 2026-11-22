"""Quasi-static magnetic field models.

* Exact field of a circular current loop (complete elliptic integrals; Simpson et al.,
  NASA/TM-2001-210854 formulation) - used as ground truth and as the solver model.
* Point magnetic dipole - used to quantify model-mismatch error.

A planar multi-layer spiral transmitter coil is represented as a set of concentric
coaxial loops (one per turn), see coils.py.
"""
from __future__ import annotations

import numpy as np
from scipy.special import ellipe, ellipk

MU0 = 4e-7 * np.pi


def _orthonormal_frame(n):
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(n, a)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    return u, v, n


def loop_field_local(rho, z, a, I=1.0):
    """B_rho, B_z (T) of a loop radius a in its own frame. Vectorised over rho, z, a."""
    rho = np.asarray(rho, float)
    z = np.asarray(z, float)
    a = np.asarray(a, float)
    r2 = rho ** 2 + z ** 2
    alpha2 = a ** 2 + r2 - 2 * a * rho
    beta2 = a ** 2 + r2 + 2 * a * rho
    beta = np.sqrt(beta2)
    k2 = 1.0 - alpha2 / beta2
    K = ellipk(k2)
    E = ellipe(k2)
    C = MU0 * I / np.pi
    Bz = C / (2 * alpha2 * beta) * ((a ** 2 - r2) * E + alpha2 * K)
    with np.errstate(invalid="ignore", divide="ignore"):
        Brho = C * z / (2 * alpha2 * beta * rho) * ((a ** 2 + r2) * E - alpha2 * K)
    Brho = np.where(rho < 1e-12, 0.0, Brho)
    return Brho, Bz


def coil_field(points, center, normal, radii, offsets, I=1.0, weights=None):
    """Field (T) of a coil = coaxial loops with given radii and axial offsets.

    points: (P,3). radii, offsets: (L,) arrays (offset along normal). Returns (P,3).
    """
    points = np.atleast_2d(points)
    u, v, n = _orthonormal_frame(normal)
    d = points - np.asarray(center)
    x = d @ u
    y = d @ v
    zc = d @ n
    rho = np.hypot(x, y)
    # broadcast (P, L)
    Brho, Bz = loop_field_local(rho[:, None], zc[:, None] - offsets[None, :], radii[None, :], I)
    if weights is not None:
        Brho = Brho * weights[None, :]
        Bz = Bz * weights[None, :]
    Brho = Brho.sum(axis=1)
    Bz = Bz.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx = np.where(rho > 1e-12, x / rho, 0.0)
        cy = np.where(rho > 1e-12, y / rho, 0.0)
    return (Brho * cx)[:, None] * u + (Brho * cy)[:, None] * v + Bz[:, None] * n


def dipole_field(points, center, moment):
    """Point-dipole field (T). moment: vector (A m^2)."""
    points = np.atleast_2d(points)
    r = points - np.asarray(center)
    rn = np.linalg.norm(r, axis=1)[:, None]
    rh = r / rn
    m = np.asarray(moment)
    return MU0 / (4 * np.pi) * (3 * rh * (rh @ m)[:, None] - m) / rn ** 3
