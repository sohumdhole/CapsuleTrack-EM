"""Monte Carlo accuracy study across ANSUR II body sizes."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from .anatomy import Body
from .system import ErrorModel, Patch, Receiver, default_freqs, measure, perturb_patch
from . import solver as S

LAYOUTS = {4: (2, 2), 6: (3, 2), 8: (4, 2), 12: (4, 3)}


@dataclass
class Config:
    n_coils: int = 8
    current_pk: float = 0.08
    f0: float = 14e3
    df: float = 2.5e3
    pitch: float = 0.060
    model: str = "loop"            # solver model: "loop" (exact) or "dipole"
    nominal_body: bool = False     # solver assumes a fixed 50th-percentile body shape
    em: ErrorModel = field(default_factory=ErrorModel)
    rx: Receiver = field(default_factory=Receiver)
    poses: int = 8
    grid_n: int = 6


def make_patch(body: Body, cfg: Config) -> Patch:
    cols, rows = LAYOUTS[cfg.n_coils]
    C, N = body.patch_coils(n_cols=cols, n_rows=rows, pitch=cfg.pitch)
    return Patch(C, N, default_freqs(cfg.n_coils, cfg.f0, cfg.df), cfg.current_pk)


NOMINAL = dict(waistbreadth=316.0, waistdepth=227.0, stature=1715.0)


def run_subject(args):
    row, cfg, seed = args
    row = pd.Series(row)
    rng = np.random.default_rng(seed)
    body = Body.from_row(row)
    patch_nom = make_patch(body, cfg)
    if cfg.nominal_body:
        nb = Body(a=NOMINAL["waistbreadth"] / 2000, b=NOMINAL["waistdepth"] / 2000,
                  stature=NOMINAL["stature"] / 1000)
        # solver believes the patch sits on the nominal body, re-centred to the true navel
        pn = make_patch(nb, cfg)
        shift = np.array([0.0, body.b - nb.b, 0.0])
        patch_model = Patch(pn.centers + shift, pn.normals, pn.freqs, pn.current_pk)
    else:
        patch_model = patch_nom
    patch_true = perturb_patch(patch_nom, cfg.em, rng)
    grid = body.workspace_grid(cfg.grid_n)
    grid_fine = body.workspace_grid(10)
    out = []
    for k in range(cfg.poses):
        p = body.sample_workspace(1, rng)[0]
        R = Rotation.random(random_state=rng.integers(1 << 31)).as_matrix()
        b = measure(patch_true, p, R, cfg.rx, rng, cfg.em)
        sig = cfg.em.sigma(cfg.rx, patch_true.freqs, np.abs(b))
        pe, Re, res = S.solve(patch_model, b, sig, grid, model=cfg.model, grid_fine=grid_fine)
        sig_true = cfg.em.sigma(cfg.rx, patch_true.freqs, np.abs(patch_true.fields(p)[0] @ R))
        crlb_p, _ = S.crlb(patch_nom, p, R, sig_true)
        anterior_depth = body.b - p[1]
        dmin = np.min(np.linalg.norm(patch_true.centers - p, axis=1))
        out.append(dict(subjectid=row.subjectid, gender=row.Gender,
                        waistbreadth=row.waistbreadth, waistdepth=row.waistdepth,
                        waistcircumference=row.waistcircumference, bmi=row.bmi,
                        px=p[0], py=p[1], pz=p[2], depth_mm=anterior_depth * 1e3,
                        dist_nearest_coil_mm=dmin * 1e3,
                        pos_err_mm=np.linalg.norm(pe - p) * 1e3,
                        ori_err_deg=S.orientation_error_deg(Re, R),
                        crlb_mm=crlb_p * 1e3, nfev=res.nfev,
                        chi2_dof=2 * res.cost / max(b.size - 6, 1), retried=res.retried))
    return out


def run(df: pd.DataFrame, cfg: Config, seed=0, processes=2) -> pd.DataFrame:
    args = [(row, cfg, seed * 100003 + i) for i, row in enumerate(df.to_dict("records"))]
    with Pool(processes) as pool:
        chunks = pool.map(run_subject, args, chunksize=4)
    return pd.DataFrame([r for c in chunks for r in c])


def sample_subjects(d: pd.DataFrame, n_per_sex: int, seed=0) -> pd.DataFrame:
    """Random sample per sex (male/female analysed as separate strata)."""
    parts = [g.sample(n=min(n_per_sex, len(g)), random_state=seed)
             for _, g in d.groupby("Gender")]
    return pd.concat(parts, ignore_index=True)


def summarize(r: pd.DataFrame, by=None) -> pd.DataFrame:
    def f(g):
        return pd.Series(dict(n=len(g), pos_median_mm=g.pos_err_mm.median(),
                              pos_rmse_mm=np.sqrt(np.mean(g.pos_err_mm ** 2)),
                              pos_p95_mm=g.pos_err_mm.quantile(0.95),
                              frac_le_1mm=(g.pos_err_mm <= 1.0).mean(),
                              gross_fail_rate_gt10mm=(g.pos_err_mm > 10.0).mean(),
                              ori_median_deg=g.ori_err_deg.median(),
                              crlb_median_mm=g.crlb_mm.median()))
    if by is None:
        return f(r).to_frame().T
    return r.groupby(by, observed=True).apply(f, include_groups=False).reset_index()
