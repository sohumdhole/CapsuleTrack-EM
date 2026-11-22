"""Tissue electromagnetics at 1-50 kHz from the Gabriel (1996) 4-Cole-Cole model.

Quantities computed per tissue and frequency:
  * conductivity sigma (S/m) and relative permittivity eps_r'
  * loss ratio sigma / (w eps0 eps_r')  (>> 1 -> conduction dominates)
  * electromagnetic skin depth delta (general lossy-medium expression)
  * induction number N = (L / delta)^2 for a path length L through tissue, the
    first-order relative size of the eddy-current (secondary) magnetic field.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

EPS0 = 8.8541878128e-12
MU0 = 4e-7 * np.pi
DATA = Path(__file__).resolve().parents[1] / "data" / "raw" / "tissue" / "gabriel1996_colecole.csv"

# Optional: map tissue -> LF conductivity (S/m) to override Gabriel values, e.g. with
# the IT'IS V4.2 low-frequency conductivity table. Empty by default.
lf_conductivity_override: dict[str, float] = {}

ABDOMINAL = ["Skin (Wet)", "Fat (Not Infiltrated)", "Muscle", "Small Intestine",
             "Colon", "Stomach", "Liver", "Blood", "Kidney", "Spleen", "Bladder"]


def load_params(path: Path = DATA) -> pd.DataFrame:
    return pd.read_csv(path).set_index("tissue")


def complex_permittivity(f, p: pd.Series) -> np.ndarray:
    """Relative complex permittivity eps' - j eps'' (Gabriel et al. 1996, eq. 3)."""
    f = np.asarray(f, dtype=float)
    w = 2 * np.pi * f
    taus = [p.tau1_ps * 1e-12, p.tau2_ns * 1e-9, p.tau3_us * 1e-6, p.tau4_ms * 1e-3]
    dels = [p.del1, p.del2, p.del3, p.del4]
    alfs = [p.alf1, p.alf2, p.alf3, p.alf4]
    eps = np.full(f.shape, p.ef, dtype=complex)
    for d, t, a in zip(dels, taus, alfs):
        eps = eps + d / (1 + (1j * w * t) ** (1 - a))
    eps = eps + p.sig / (1j * w * EPS0)
    return eps


def conductivity(f, p: pd.Series) -> np.ndarray:
    name = p.name
    if name in lf_conductivity_override:
        return np.full(np.shape(f), lf_conductivity_override[name], dtype=float)
    eps = complex_permittivity(f, p)
    return -2 * np.pi * np.asarray(f) * EPS0 * eps.imag


def rel_permittivity(f, p: pd.Series) -> np.ndarray:
    return complex_permittivity(f, p).real


def skin_depth(f, sigma, eps_r) -> np.ndarray:
    """Penetration depth 1/alpha for a lossy dielectric (mu = mu0)."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    eps = EPS0 * np.asarray(eps_r)
    loss = sigma / (w * eps)
    alpha = w * np.sqrt(MU0 * eps / 2 * (np.sqrt(1 + loss ** 2) - 1))
    return 1.0 / alpha


def induction_number(f, sigma, eps_r, L) -> np.ndarray:
    return (np.asarray(L) / skin_depth(f, sigma, eps_r)) ** 2


def table(freqs, tissues=None, L=0.30) -> pd.DataFrame:
    """Long-format table of tissue EM quantities over frequency."""
    P = load_params()
    tissues = tissues or ABDOMINAL
    rows = []
    for t in tissues:
        p = P.loc[t]
        s = conductivity(freqs, p)
        e = rel_permittivity(freqs, p)
        d = skin_depth(freqs, s, e)
        for fi, si, ei, di in zip(freqs, s, e, d):
            rows.append(dict(tissue=t, f_hz=fi, sigma_S_per_m=si, eps_r=ei,
                             loss_ratio=si / (2 * np.pi * fi * EPS0 * ei),
                             skin_depth_m=di, induction_number=(L / di) ** 2))
    return pd.DataFrame(rows)
