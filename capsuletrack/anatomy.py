"""Subject-specific abdomen geometry from ANSUR II and patch placement.

Coordinate frame (metres): origin on the body axis at omphalion (navel) height,
x = left-right, y = anterior (+) / posterior (-), z = superior (+).

Abdomen cross-section: ellipse with semi-axes a = waistbreadth/2 (x) and
b = waistdepth/2 (y), extruded along z. The flexible patch conforms to the anterior
surface, so each subject gets a different coil geometry (positions + normals).

Capsule workspace (GI tract proxy): the region inside the ellipse inset by an
abdominal-wall thickness, anterior of y = -0.3 b (excludes spine/retroperitoneum),
and within +/- 0.09 * stature of the navel. This is a geometric proxy, not an
organ segmentation; ANSUR II carries no internal anatomy.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANSUR = ROOT / "data" / "raw" / "ansur"
COLS = ["subjectid", "Gender", "waistbreadth", "waistdepth", "waistcircumference",
        "stature", "weightkg", "abdominalextensiondepthsitting"]


def load_ansur() -> pd.DataFrame:
    frames = []
    for fn in ["ANSUR II MALE Public.csv", "ANSUR II FEMALE Public.csv"]:
        d = pd.read_csv(ANSUR / fn, encoding="latin1").rename(columns={"SubjectId": "subjectid"})
        frames.append(d[COLS])
    d = pd.concat(frames, ignore_index=True)
    d["mass_kg"] = d["weightkg"] / 10.0
    d["bmi"] = d["mass_kg"] / (d["stature"] / 1000) ** 2
    return d


@dataclass
class Body:
    a: float          # semi-axis x (m)
    b: float          # semi-axis y (m)
    stature: float    # m
    wall: float = 0.025  # abdominal wall thickness (m), fixed assumption

    @classmethod
    def from_row(cls, row, wall=0.025):
        return cls(a=row.waistbreadth / 2000, b=row.waistdepth / 2000,
                   stature=row.stature / 1000, wall=wall)

    # --- anterior surface parametrised by arc length s from the midline ---
    def _arc_table(self, n=2000):
        t = np.linspace(-np.pi / 2, np.pi / 2, n)   # t=0 -> anterior midline
        x = self.a * np.sin(t)
        y = self.b * np.cos(t)
        ds = np.hypot(np.gradient(x), np.gradient(y))
        s = np.cumsum(ds)
        s -= np.interp(0.0, t, s)
        return t, s

    def surface_point(self, s_arc):
        t_tab, s_tab = self._arc_table()
        t = np.interp(s_arc, s_tab, t_tab)
        p = np.array([self.a * np.sin(t), self.b * np.cos(t)])
        n = np.array([np.sin(t) / self.a, np.cos(t) / self.b])  # outward normal
        n /= np.linalg.norm(n)
        return p, n

    def patch_coils(self, n_cols=4, n_rows=2, pitch=0.060, standoff=0.001):
        """Coil centres and outward normals for a conformal patch centred on the navel."""
        cols = (np.arange(n_cols) - (n_cols - 1) / 2) * pitch
        rows = (np.arange(n_rows) - (n_rows - 1) / 2) * pitch
        C, N = [], []
        for zr in rows:
            for sc in cols:
                p2, n2 = self.surface_point(sc)
                C.append([p2[0] + standoff * n2[0], p2[1] + standoff * n2[1], zr])
                N.append([n2[0], n2[1], 0.0])
        return np.array(C), np.array(N)

    def sample_workspace(self, n, rng):
        ai, bi = self.a - self.wall, self.b - self.wall
        H = 0.09 * self.stature
        out = []
        while len(out) < n:
            x = rng.uniform(-ai, ai)
            y = rng.uniform(-0.3 * self.b, bi)
            if (x / ai) ** 2 + (y / bi) ** 2 <= 1.0:
                out.append([x, y, rng.uniform(-H, H)])
        return np.array(out)

    def workspace_grid(self, n=7):
        ai, bi = self.a - self.wall, self.b - self.wall
        H = 0.09 * self.stature
        g = []
        for x in np.linspace(-ai, ai, n):
            for y in np.linspace(-0.3 * self.b, bi, n):
                if (x / ai) ** 2 + (y / bi) ** 2 <= 1.05:
                    for z in np.linspace(-H, H, n):
                        g.append([x, y, z])
        return np.array(g)
