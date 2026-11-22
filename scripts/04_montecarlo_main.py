"""Deliverable 6 - Monte Carlo accuracy study across ANSUR II body sizes (main design).

300 randomly sampled subjects per sex x 8 random capsule poses each, nominal error
budget (ErrorModel defaults). Outputs:
  results/tables/mc_main_raw.csv, mc_main_summary.csv,
  mc_by_waistdepth_quartile.csv, mc_by_depth_band.csv
  results/figures/mc_error_vs_depth.png, mc_by_body_size.png, mc_error_cdf.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from capsuletrack.anatomy import load_ansur  # noqa: E402
from capsuletrack import montecarlo as M  # noqa: E402
from capsuletrack.plotstyle import apply, SERIES  # noqa: E402

TAB, FIG = ROOT / "results" / "tables", ROOT / "results" / "figures"
N_PER_SEX = int(sys.argv[1]) if len(sys.argv) > 1 else 300

if __name__ == "__main__":
    plt = apply()
    d = load_ansur()
    subj = M.sample_subjects(d, N_PER_SEX, seed=42)
    r = M.run(subj, M.Config(), seed=1)
    r.to_csv(TAB / "mc_main_raw.csv", index=False)

    overall = M.summarize(r).assign(group="All")
    by_sex = M.summarize(r, "gender").rename(columns={"gender": "group"})
    pd.concat([overall, by_sex]).to_csv(TAB / "mc_main_summary.csv", index=False)

    r["waistdepth_quartile"] = pd.qcut(r.waistdepth, 4, labels=["Q1 (slimmest)", "Q2", "Q3", "Q4 (largest)"])
    q = M.summarize(r, "waistdepth_quartile")
    edges = r.groupby("waistdepth_quartile", observed=True).waistdepth.agg(["min", "max"]).reset_index()
    q = q.merge(edges, on="waistdepth_quartile")
    q.to_csv(TAB / "mc_by_waistdepth_quartile.csv", index=False)

    r["depth_band"] = pd.cut(r.depth_mm, [0, 50, 100, 150, 200, 400],
                             labels=["<50", "50-100", "100-150", "150-200", ">200"])
    db = M.summarize(r, "depth_band")
    db.to_csv(TAB / "mc_by_depth_band.csv", index=False)

    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    for g, c in zip(["Female", "Male"], SERIES):
        s = r[r.gender == g]
        ax.scatter(s.depth_mm, s.pos_err_mm, s=8, alpha=0.35, color=c, label=g, edgecolor="none")
    ax.axhline(1.0, color="#6b6a63", lw=1, ls="--")
    ax.text(5, 1.08, "1 mm target", fontsize=8, color="#6b6a63")
    ax.set_yscale("log"); ax.set_xlabel("Capsule depth below anterior skin (mm)")
    ax.set_ylabel("Position error (mm)")
    ax.set_title("Error grows with depth below the patch")
    ax.legend(fontsize=8, markerscale=2)
    fig.savefig(FIG / "mc_error_vs_depth.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    x = np.arange(len(q))
    ax.bar(x - 0.2, q.pos_median_mm, 0.38, color=SERIES[0], label="Median")
    ax.bar(x + 0.2, q.pos_p95_mm, 0.38, color=SERIES[1], label="95th percentile")
    ax.set_xticks(x, [f"{a}\n{int(lo)}-{int(hi)} mm" for a, lo, hi in zip(q.waistdepth_quartile, q["min"], q["max"])], fontsize=8)
    ax.axhline(1.0, color="#6b6a63", lw=1, ls="--")
    ax.set_ylabel("Position error (mm)"); ax.set_title("Accuracy by ANSUR II waist-depth quartile")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "mc_by_body_size.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for g, c in zip(["Female", "Male"], SERIES):
        e = np.sort(r[r.gender == g].pos_err_mm)
        ax.plot(e, np.arange(1, len(e) + 1) / len(e), color=c, label=g)
    ax.axvline(1.0, color="#6b6a63", lw=1, ls="--")
    ax.set_xscale("log"); ax.set_xlabel("Position error (mm)"); ax.set_ylabel("Fraction of poses")
    ax.set_title("Cumulative position error"); ax.legend(fontsize=8)
    fig.savefig(FIG / "mc_error_cdf.png"); plt.close(fig)

    print(pd.concat([overall, by_sex]).round(3).to_string(index=False))
    print(q.round(3).to_string(index=False))
    print(db.round(3).to_string(index=False))
