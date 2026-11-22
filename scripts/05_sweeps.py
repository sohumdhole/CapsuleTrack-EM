"""Design-space sweeps on a fixed ANSUR II subset (30 subjects per sex x 5 poses).

Answers: how many coils, how much noise/drift, how much patch-placement error, how
much drive current, does a dipole model suffice, and does per-patient geometry matter?
Outputs: results/tables/sweeps.csv, results/figures/sweep_*.png
"""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from capsuletrack.anatomy import load_ansur  # noqa: E402
from capsuletrack import montecarlo as M  # noqa: E402
from capsuletrack.system import ErrorModel  # noqa: E402
from capsuletrack.plotstyle import apply, SERIES  # noqa: E402

TAB, FIG = ROOT / "results" / "tables", ROOT / "results" / "figures"

if __name__ == "__main__":
    plt = apply()
    d = load_ansur()
    subj = M.sample_subjects(d, 30, seed=7)
    base = M.Config(poses=5)
    em0 = ErrorModel()
    runs = []

    def go(study, label, value, cfg):
        r = M.run(subj, cfg, seed=3)
        s = M.summarize(r).assign(study=study, label=label, value=value)
        runs.append(s)
        print(study, label, s[["pos_median_mm", "pos_p95_mm", "frac_le_1mm", "crlb_median_mm"]].round(3).values)

    for n in [4, 6, 8, 12]:
        go("n_coils", f"{n} coils", n, replace(base, n_coils=n))
    for a in [0.1e-9, 0.3e-9, 1e-9, 3e-9, 10e-9]:
        go("ambient_nT", f"{a*1e9:g} nT", a * 1e9, replace(base, em=replace(em0, ambient_T=a)))
    for g in [0.0, 0.001, 0.002, 0.005, 0.01]:
        go("gain_rel_pct", f"{g*100:g} %", g * 100, replace(base, em=replace(em0, gain_rel=g)))
    for pos, tilt in [(0.0, 0.0), (0.25e-3, 0.5), (0.5e-3, 1.0), (1e-3, 2.0), (2e-3, 4.0)]:
        go("patch_geometry", f"{pos*1e3:g} mm / {tilt:g} deg", pos * 1e3,
           replace(base, em=replace(em0, coil_pos_m=pos, coil_tilt_deg=tilt)))
    for i in [0.02, 0.04, 0.08, 0.16]:
        go("drive_mA", f"{i*1e3:g} mA", i * 1e3, replace(base, current_pk=i))
    go("model", "exact loop model", 0, base)
    go("model", "point-dipole model", 1, replace(base, model="dipole"))
    go("body_shape", "per-patient ellipse (ANSUR)", 0, base)
    go("body_shape", "fixed nominal body", 1, replace(base, nominal_body=True))
    go("noise_only", "thermal + ambient only", 0,
       replace(base, em=replace(em0, gain_rel=0.0, coil_pos_m=0.0, coil_tilt_deg=0.0)))

    sw = pd.concat(runs, ignore_index=True)
    sw.to_csv(TAB / "sweeps.csv", index=False)

    def panel(study, xlabel, fname, title, logx=False):
        s = sw[sw.study == study]
        fig, ax = plt.subplots(figsize=(5.6, 3.4))
        ax.plot(s.value, s.pos_median_mm, "-o", color=SERIES[0], label="Median")
        ax.plot(s.value, s.pos_p95_mm, "-s", color=SERIES[1], label="95th percentile")
        ax.axhline(1.0, color="#6b6a63", lw=1, ls="--")
        if logx:
            ax.set_xscale("symlog", linthresh=0.1)
        ax.set_yscale("log"); ax.set_xlabel(xlabel); ax.set_ylabel("Position error (mm)")
        ax.set_title(title); ax.legend(fontsize=8)
        fig.savefig(FIG / fname); plt.close(fig)

    panel("n_coils", "Number of patch coils", "sweep_n_coils.png", "More coils, lower error")
    panel("ambient_nT", "Ambient interference (nT)", "sweep_ambient.png", "Sensitivity to ambient interference", True)
    panel("gain_rel_pct", "Gain / drift error (%)", "sweep_gain.png", "Sensitivity to gain drift")
    panel("patch_geometry", "Coil placement error (mm, with 2x deg tilt)", "sweep_geometry.png",
          "Patch geometry knowledge dominates")
    panel("drive_mA", "Drive current per coil (mA peak)", "sweep_current.png", "Drive current vs accuracy")

    cmp = sw[sw.study.isin(["model", "body_shape", "noise_only"])]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    y = np.arange(len(cmp))
    ax.barh(y, cmp.pos_median_mm, color=SERIES[0], height=0.6)
    ax.set_yticks(y, cmp.label, fontsize=8); ax.set_xscale("log")
    ax.axvline(1.0, color="#6b6a63", lw=1, ls="--")
    for yi, v in zip(y, cmp.pos_median_mm):
        ax.text(v * 1.08, yi, f"{v:.2f} mm", va="center", fontsize=8, color="#1f1f1e")
    ax.set_xlabel("Median position error (mm)"); ax.set_title("Model choices")
    fig.savefig(FIG / "sweep_model_choices.png"); plt.close(fig)
