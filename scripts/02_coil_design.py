"""Coil electrical design, power budget and ICNIRP 2010 exposure check.

Outputs
  results/tables/coil_design.csv        geometry-derived L, R, moment
  results/tables/drive_per_channel.csv  per-tone resonant C, Q, power
  results/tables/exposure_check.csv     B at skin and induced E vs ICNIRP 2010 limits
  results/figures/field_vs_depth.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from capsuletrack.coils import SpiralCoil, cu_skin_depth  # noqa: E402
from capsuletrack.field import coil_field  # noqa: E402
from capsuletrack.system import Receiver, ErrorModel, default_freqs  # noqa: E402
from capsuletrack.plotstyle import apply, SERIES  # noqa: E402

plt = apply()
TAB, FIG = ROOT / "results" / "tables", ROOT / "results" / "figures"

coil = SpiralCoil()
I_PK = 0.08
R_ON = 0.28          # H-bridge high+low side on-resistance (ohm), DRV8837-class
R_SENSE = 0.5        # current-sense resistor (ohm)
V_RAIL = 3.3        # H-bridge supply rail (V)
RECORD_H = 12.0      # design recording duration for battery sizing (h)

s = coil.summary()
s["cu_skin_depth_um_at_31.5kHz"] = cu_skin_depth(31.5e3) * 1e6
s["cu_thickness_um"] = coil.cu_t * 1e6
s["drive_current_pk_A"] = I_PK
s["moment_pk_Am2"] = coil.moment_per_amp() * I_PK
pd.DataFrame([s]).to_csv(TAB / "coil_design.csv", index=False)

L, R = coil.inductance(), coil.resistance_dc()
freqs = default_freqs(8)
rows = []
for i, f in enumerate(freqs):
    w = 2 * np.pi * f
    C = 1 / (w ** 2 * L)
    R_tot = R + R_ON + R_SENSE
    p_coil = 0.5 * I_PK ** 2 * R
    p_loss = 0.5 * I_PK ** 2 * (R_ON + R_SENSE)
    # Phase-shifted full bridge from a 3.3 V rail: three-level output with pulse width
    # phi gives fundamental amplitude V1 = (4 Vrail / pi) sin(phi / 2). Firmware sets phi
    # from the sense-resistor current (closed loop).
    v1 = I_PK * R_tot
    phi = 2 * np.degrees(np.arcsin(v1 * np.pi / (4 * V_RAIL)))
    rows.append(dict(channel=i + 1, f_hz=f, L_uH=L * 1e6, C_res_nF=C * 1e9,
                     X_L_ohm=w * L, R_total_ohm=R_tot, Q=w * L / R_tot,
                     V1_fundamental_V=v1, phase_shift_deg=phi, P_coil_mW=p_coil * 1e3, P_driver_mW=p_loss * 1e3))
drv = pd.DataFrame(rows)
drv.to_csv(TAB / "drive_per_channel.csv", index=False)
P_tot = drv.P_coil_mW.sum() + drv.P_driver_mW.sum()
P_sys = P_tot / 0.85 + 30          # 85 % converter efficiency + 30 mW MCU/radio budget
E_wh = P_sys * RECORD_H / 1e3
budget = dict(P_coils_and_drivers_mW=P_tot, P_system_mW=P_sys, record_hours=RECORD_H,
              energy_Wh=E_wh, battery_mAh_at_3V7=E_wh / 3.7 * 1e3)
pd.DataFrame([budget]).to_csv(TAB / "power_budget.csv", index=False)

# --- ICNIRP 2010 general-public check (1 Hz - 100 kHz guidelines) ---
# Basic restriction, 3 kHz-10 MHz, all tissues of head and body: E_int = 1.35e-4 f V/m (rms)
# Reference level, 3 kHz-10 MHz: B = 2.7e-5 T (rms)
r, z = coil.loops()
wts = coil.weights
exp_rows = []
for gap in [0.5e-3, 2e-3, 5e-3]:
    plane_z = min(coil.layer_z) - gap   # below the bottom copper layer, skin side
    rr = np.linspace(0, 40e-3, 400)
    pts = np.c_[rr, np.zeros_like(rr), np.full_like(rr, plane_z)]
    B = coil_field(pts, [0, 0, 0], [0, 0, 1], r, z, I=I_PK, weights=wts)
    Bmag = np.linalg.norm(B, axis=1)
    # flux through disc radius rho in that plane, Faraday: E(rho) = w Phi / (2 pi rho)
    Bz = B[:, 2]
    Phi = np.concatenate([[0], np.cumsum(2 * np.pi * rr[1:] * Bz[1:] * np.diff(rr))])
    for f in [freqs.min(), freqs.max()]:
        w = 2 * np.pi * f
        with np.errstate(divide="ignore", invalid="ignore"):
            E = np.where(rr > 0, w * np.abs(Phi) / (2 * np.pi * rr), 0)
        exp_rows.append(dict(gap_below_patch_mm=gap * 1e3, f_hz=f,
                             B_max_rms_uT=Bmag.max() / np.sqrt(2) * 1e6,
                             B_ref_level_uT=27.0,
                             E_induced_max_rms_V_per_m=E.max() / np.sqrt(2),
                             E_basic_restriction_V_per_m=1.35e-4 * f,
                             E_margin=1.35e-4 * f / (E.max() / np.sqrt(2))))
exp = pd.DataFrame(exp_rows)
exp.to_csv(TAB / "exposure_check.csv", index=False)

# --- field magnitude vs depth on coil axis, against the noise floor ---
depth = np.linspace(0.02, 0.30, 200)
pts = np.c_[np.zeros_like(depth), np.zeros_like(depth), -depth]
B = np.linalg.norm(coil_field(pts, [0, 0, 0], [0, 0, 1], r, z, I=I_PK, weights=wts), axis=1)
rx, em = Receiver(), ErrorModel()
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.plot(depth * 1e3, B * 1e9, color=SERIES[0], label="One patch coil, 80 mA pk, on axis")
ax.axhline(rx.sigma_B(freqs.min()) * 1e9, color=SERIES[1], lw=1.5, label="Receiver thermal noise (10 Hz BW)")
ax.axhline(em.ambient_T * 1e9, color=SERIES[2], lw=1.5, label="Assumed ambient interference, 1 nT")
ax.set_yscale("log"); ax.set_xlabel("Depth below patch (mm)"); ax.set_ylabel("|B| (nT)")
ax.set_title("Signal vs noise floor along the coil axis")
ax.legend(fontsize=8)
fig.savefig(FIG / "field_vs_depth.png"); plt.close(fig)

print(pd.DataFrame([s]).T.round(4).to_string())
print(drv.round(2).to_string(index=False))
print(budget)
print(exp.round(3).to_string(index=False))
