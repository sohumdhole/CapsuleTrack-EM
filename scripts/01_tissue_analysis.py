"""Deliverable 1 - tissue analysis at 1-50 kHz.

Outputs
  results/tables/tissue_em_properties.csv   sigma, eps_r, loss ratio, skin depth per tissue/f
  results/tables/eddy_bias_by_frequency.csv worst-case eddy-current field error and
                                            equivalent position bias vs frequency and depth
  results/figures/tissue_conductivity.png, tissue_skin_depth.png, eddy_bias.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from capsuletrack import tissue as T  # noqa: E402
from capsuletrack.plotstyle import apply, SERIES  # noqa: E402

plt = apply()
TAB, FIG = ROOT / "results" / "tables", ROOT / "results" / "figures"

freqs = np.geomspace(1e3, 50e3, 60)
tab = T.table(freqs)
tab.to_csv(TAB / "tissue_em_properties.csv", index=False)

# Headline table at design-relevant frequencies
key = T.table(np.array([1e3, 10e3, 20e3, 30e3, 50e3]))
key.to_csv(TAB / "tissue_em_key_frequencies.csv", index=False)

show = ["Fat (Not Infiltrated)", "Colon", "Muscle", "Small Intestine", "Blood"]
fig, ax = plt.subplots(figsize=(6.4, 3.8))
for t, c in zip(show, SERIES):
    d = tab[tab.tissue == t]
    ax.plot(d.f_hz / 1e3, d.sigma_S_per_m, color=c, label=t)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Frequency (kHz)"); ax.set_ylabel("Conductivity (S/m)")
ax.set_title("Tissue conductivity is almost flat across 1-50 kHz")
ax.legend(fontsize=8, ncol=2)
fig.savefig(FIG / "tissue_conductivity.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.4, 3.8))
for t, c in zip(show, SERIES):
    d = tab[tab.tissue == t]
    ax.plot(d.f_hz / 1e3, d.skin_depth_m, color=c, label=t)
ax.axhline(0.30, color="#6b6a63", lw=1, ls="--")
ax.text(1.1, 0.33, "abdomen depth ~0.3 m", fontsize=8, color="#6b6a63")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Frequency (kHz)"); ax.set_ylabel("Skin depth (m)")
ax.set_title("Skin depth is metres: the abdomen is quasi-static at 1-50 kHz")
ax.legend(fontsize=8, ncol=2)
fig.savefig(FIG / "tissue_skin_depth.png"); plt.close(fig)

# Worst-case eddy-current error: whole path filled with the most conductive bulk
# abdominal tissue (small intestine; blood shown as an upper bound).
P = T.load_params()
rows = []
for t in ["Small Intestine", "Blood"]:
    s = T.conductivity(freqs, P.loc[t]); e = T.rel_permittivity(freqs, P.loc[t])
    d = T.skin_depth(freqs, s, e)
    for L in [0.10, 0.20, 0.30]:
        eps = (L / d) ** 2                   # relative secondary field (order of magnitude)
        bias = L * eps / 3                   # |B| ~ r^-3  ->  dr/r = dB/B / 3
        for fi, ei, bi in zip(freqs, eps, bias):
            rows.append(dict(tissue=t, depth_m=L, f_hz=fi, rel_field_error=ei,
                             pos_bias_mm_inphase_worstcase=bi * 1e3))
eb = pd.DataFrame(rows)
eb.to_csv(TAB / "eddy_bias_by_frequency.csv", index=False)

fig, ax = plt.subplots(figsize=(6.4, 3.8))
for L, c in zip([0.10, 0.20, 0.30], SERIES):
    d = eb[(eb.tissue == "Small Intestine") & (eb.depth_m == L)]
    ax.plot(d.f_hz / 1e3, d.pos_bias_mm_inphase_worstcase, color=c, label=f"{int(L*100)} cm path")
ax.axvspan(14, 31.5, color="#2a78d6", alpha=0.08)
ax.text(15, 2e-3, "chosen band\n14-31.5 kHz", fontsize=8, color="#1c5cab")
ax.axhline(1.0, color="#6b6a63", lw=1, ls="--")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Frequency (kHz)"); ax.set_ylabel("Worst-case position bias (mm)")
ax.set_title("Eddy-current bias through small-intestine-filled path")
ax.legend(fontsize=8, loc="upper left")
fig.savefig(FIG / "eddy_bias.png"); plt.close(fig)

k = key[key.tissue.isin(["Small Intestine", "Muscle", "Fat (Not Infiltrated)"])]
print(k[["tissue", "f_hz", "sigma_S_per_m", "loss_ratio", "skin_depth_m"]].round(3).to_string(index=False))
w = eb[(eb.tissue == "Small Intestine") & (eb.depth_m == 0.2)]
for f in [14e3, 31.5e3, 50e3]:
    i = np.argmin(abs(w.f_hz - f))
    print(f"bias @ {f/1e3:.1f} kHz, 20 cm: {w.pos_bias_mm_inphase_worstcase.iloc[i]:.3f} mm")
