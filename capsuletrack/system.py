"""Tracking system: transmitter patch + capsule receiver + measurement model.

Architecture (patient-worn transmitter, capsule-borne receiver):
  * N planar coils on the patch, each driven with its own sinusoidal tone
    (frequency-division multiplexing) in the 1-50 kHz band.
  * The capsule carries a 3-axis ferrite-cored receive coil. Lock-in demodulation of
    each tone gives the 3-component field of every transmitter in the capsule frame:
        b_i = R^T B_i(p) + noise,  i = 1..N   (3N real measurements)
  * Unknowns: capsule position p (3) and orientation R (3) -> 6 DOF.

Receiver noise is converted from a voltage noise density to an equivalent field noise:
    sigma_B(f) = e_n sqrt(BW) / (N_r A_r mu_app 2 pi f)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .coils import SpiralCoil
from .field import coil_field, dipole_field


@dataclass
class Receiver:
    n_turns: int = 400
    radius: float = 1.25e-3      # m, winding radius on ferrite core
    mu_app: float = 15.0         # apparent permeability of the short ferrite rod
    e_n: float = 5e-9            # V/sqrt(Hz), input-referred (coil thermal + LNA)
    bandwidth: float = 10.0      # Hz, lock-in equivalent noise bandwidth per tone

    def sigma_B(self, f):
        area = np.pi * self.radius ** 2
        v_n = self.e_n * np.sqrt(self.bandwidth)
        return v_n / (self.n_turns * area * self.mu_app * 2 * np.pi * np.asarray(f))


@dataclass
class Patch:
    centers: np.ndarray
    normals: np.ndarray
    freqs: np.ndarray
    current_pk: float = 0.08     # A, per coil
    coil: SpiralCoil = field(default_factory=SpiralCoil)

    @property
    def n(self):
        return len(self.centers)

    def fields(self, points, model="loop"):
        """(P, N, 3) field of every coil at every point, for the drive current."""
        points = np.atleast_2d(points)
        out = np.empty((len(points), self.n, 3))
        if model == "loop":
            r, z = self.coil.loops()
            w = self.coil.weights
            for i in range(self.n):
                out[:, i, :] = coil_field(points, self.centers[i], self.normals[i], r, z,
                                          I=self.current_pk, weights=w)
        else:
            m = self.coil.moment_per_amp() * self.current_pk
            r, z = self.coil.loops()
            zc = float(np.average(z, weights=self.coil.weights))
            for i in range(self.n):
                c = self.centers[i] + zc * self.normals[i]
                out[:, i, :] = dipole_field(points, c, m * self.normals[i])
        return out


def default_freqs(n, f0=14e3, df=2.5e3):
    return f0 + df * np.arange(n)


@dataclass
class ErrorModel:
    """Measurement and model error budget (all 1-sigma)."""
    ambient_T: float = 1e-9      # ambient interference floor in the tone bands (T)
    gain_rel: float = 2e-3       # per-measurement gain/phase-drift error (relative)
    coil_pos_m: float = 0.5e-3   # patch coil-centre placement uncertainty (m, per axis)
    coil_tilt_deg: float = 1.0   # patch coil-normal tilt uncertainty (deg)

    def sigma(self, rx: Receiver, freqs, b_abs=None):
        s2 = rx.sigma_B(freqs) ** 2 + self.ambient_T ** 2
        s2 = np.broadcast_to(s2[:, None], (len(freqs), 3)).copy()
        if b_abs is not None:
            s2 = s2 + (self.gain_rel * b_abs) ** 2
        return np.sqrt(s2)


def perturb_patch(patch: Patch, em: ErrorModel, rng) -> Patch:
    """True patch geometry differs from the nominal model by placement errors."""
    C = patch.centers + rng.normal(scale=em.coil_pos_m, size=patch.centers.shape)
    N = patch.normals + rng.normal(scale=np.radians(em.coil_tilt_deg), size=patch.normals.shape)
    N = N / np.linalg.norm(N, axis=1, keepdims=True)
    return Patch(C, N, patch.freqs, patch.current_pk, patch.coil)


def measure(patch: Patch, p, R, rx: Receiver, rng=None, em: ErrorModel | None = None):
    """Capsule-frame measurements (N, 3) from the TRUE patch, with noise if rng given."""
    B = patch.fields(p)[0]               # (N,3) world frame
    b = B @ R                            # R^T B for each row
    if rng is not None:
        em = em or ErrorModel(0.0, 0.0, 0.0, 0.0)
        s = np.sqrt(rx.sigma_B(patch.freqs) ** 2 + em.ambient_T ** 2)[:, None]
        b = b * (1 + em.gain_rel * rng.normal(size=b.shape))
        b = b + rng.normal(size=b.shape) * s
    return b
