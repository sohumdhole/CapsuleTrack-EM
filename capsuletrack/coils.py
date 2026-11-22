"""Planar multi-layer spiral transmitter coil for a flexible patch.

Electrical parameters are computed from geometry, not assumed:
  * inductance  = sum of turn self-inductances + all mutual inductances (Maxwell's
    coaxial-loop formula, complete elliptic integrals); the rectangular trace is
    represented by its geometric-mean-distance radius 0.2235 (w + t)
  * DC resistance from copper resistivity; AC skin depth of Cu is checked
  * magnetic moment per amp = sum of turn areas
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.special import ellipe, ellipk

MU0 = 4e-7 * np.pi
RHO_CU = 1.72e-8  # ohm m at 20 C


def mutual_coaxial(a, b, d):
    """Mutual inductance (H) of two coaxial circular loops radii a, b, separation d."""
    k2 = 4 * a * b / ((a + b) ** 2 + d ** 2)
    k = np.sqrt(k2)
    return MU0 * np.sqrt(a * b) * ((2 / k - k) * ellipk(k2) - 2 / k * ellipe(k2))


@dataclass
class SpiralCoil:
    r_out: float = 20e-3          # outer radius of outer turn centreline (m)
    turns_per_layer: tuple = (16.25, 16.25, 16.25, 16.125)
    pitch: float = 0.45e-3         # trace pitch (m)
    trace_w: float = 0.25e-3       # trace width (m)
    cu_t: float = 35e-6            # copper thickness, 1 oz (m)
    layer_z: tuple = (0.0, -0.2e-3, -0.4e-3, -0.6e-3)  # copper layer offsets (m)
    _cache: dict = field(default_factory=dict, repr=False)

    def loops(self):
        """Equivalent coaxial loops: (radii, axial offsets)."""
        if "loops" not in self._cache:
            radii, offs = [], []
            for n_t, z in zip(self.turns_per_layer, self.layer_z):
                n_full = int(np.floor(n_t))
                rs = self.r_out - self.pitch * np.arange(n_full)
                frac = n_t - n_full
                radii.extend(rs)
                offs.extend([z] * n_full)
                if frac > 0:  # partial turn: weight approximated by area fraction later
                    radii.append(self.r_out - self.pitch * n_full)
                    offs.append(z)
            self._cache["loops"] = (np.array(radii), np.array(offs))
            # weights for partial turns
            w = []
            for n_t in self.turns_per_layer:
                n_full = int(np.floor(n_t))
                w.extend([1.0] * n_full)
                if n_t - n_full > 0:
                    w.append(n_t - n_full)
            self._cache["w"] = np.array(w)
        return self._cache["loops"]

    @property
    def weights(self):
        self.loops()
        return self._cache["w"]

    @property
    def n_turns(self):
        return float(sum(self.turns_per_layer))

    @property
    def r_in(self):
        return self.r_out - self.pitch * max(self.turns_per_layer)

    def inductance(self):
        r, z = self.loops()
        w = self.weights
        gmd = 0.2235 * (self.trace_w + self.cu_t)
        L = 0.0
        for i in range(len(r)):
            L += w[i] ** 2 * MU0 * r[i] * (np.log(8 * r[i] / gmd) - 2.0)
            for j in range(i + 1, len(r)):
                L += 2 * w[i] * w[j] * mutual_coaxial(r[i], r[j], abs(z[i] - z[j]))
        return float(L)

    def length(self):
        r, _ = self.loops()
        return float(np.sum(2 * np.pi * r * self.weights))

    def resistance_dc(self):
        return RHO_CU * self.length() / (self.trace_w * self.cu_t)

    def moment_per_amp(self):
        r, _ = self.loops()
        return float(np.sum(np.pi * r ** 2 * self.weights))

    def summary(self):
        return dict(r_out_mm=self.r_out * 1e3, r_in_mm=self.r_in * 1e3,
                    turns=self.n_turns, layers=len(self.turns_per_layer),
                    length_m=self.length(), R_dc_ohm=self.resistance_dc(),
                    L_uH=self.inductance() * 1e6,
                    moment_Am2_per_A=self.moment_per_amp())


def cu_skin_depth(f):
    return np.sqrt(2 * RHO_CU / (2 * np.pi * f * MU0))
