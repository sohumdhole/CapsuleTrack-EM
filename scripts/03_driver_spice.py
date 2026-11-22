"""Deliverable 2 - transmitter coil driver simulation (SPICE).

Circuit (one channel): phase-shifted full-bridge from a 3.3 V rail -> series resonant
capacitor -> patch coil (L, R from coils.py) -> 0.5 ohm current-sense resistor.
The bridge legs are modelled as ideal square-wave sources with the H-bridge on-
resistance split across the two legs (DRV8837-class, 0.28 ohm HS+LS).

Generates LTspice-compatible files (hardware/ltspice/*.cir, *.asc) and runs the same
netlists in ngspice for: transient (current amplitude, rms, harmonics, power), AC
resonance sweep, and a component-tolerance Monte Carlo.
Outputs: results/tables/spice_*.csv, results/figures/spice_*.png
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from capsuletrack.coils import SpiralCoil  # noqa: E402
from capsuletrack.plotstyle import apply, SERIES  # noqa: E402

plt = apply()
HW = ROOT / "hardware" / "ltspice"
TAB, FIG = ROOT / "results" / "tables", ROOT / "results" / "figures"
WORK = Path("/tmp/ct_spice"); WORK.mkdir(exist_ok=True)

coil = SpiralCoil()
L = coil.inductance(); RC = coil.resistance_dc()
RON, RSENSE, VRAIL, IPK = 0.28, 0.5, 3.3, 0.08


def design(f):
    C = 1 / ((2 * np.pi * f) ** 2 * L)
    v1 = IPK * (RC + RON + RSENSE)
    phi = 2 * np.degrees(np.arcsin(v1 * np.pi / (4 * VRAIL)))
    return C, phi


def tran_netlist(f, C, phi, Lv=L, Rv=RC, control=True):
    T = 1 / f
    s = f"""* CapsuleTrack-EM patch coil driver - one channel, f = {f/1e3:g} kHz
* Phase-shifted full bridge (3-level output), series-resonant coil, current sense.
.param F={f:.6g} VRAIL={VRAIL} PHI={phi:.4g} L={Lv:.6g} RCOIL={Rv:.6g} C={C:.6g}
.param RON={RON} RSENSE={RSENSE}
VA a 0 PULSE(0 {{VRAIL}} 0 20n 20n {{0.5/F-20n}} {{1/F}})
VB b 0 PULSE(0 {{VRAIL}} {{PHI/360/F}} 20n 20n {{0.5/F-20n}} {{1/F}})
RONA a n1 {{RON/2}}
C1 n1 n2 {{C}}
L1 n2 n3 {{L}}
RCOIL n3 n4 {{RCOIL}}
VIM n4 n5 0
RSENSE n5 nb {{RSENSE}}
RONB nb b {{RON/2}}
.tran 0.05u {60*T:.6g} {40*T:.6g} 0.05u
.meas tran Ipk MAX i(VIM) FROM={40*T:.6g} TO={60*T:.6g}
.meas tran Irms RMS i(VIM) FROM={40*T:.6g} TO={60*T:.6g}
.four {{F}} i(VIM)
"""
    if control:
        s += f"""
.control
run
meas tran ipk MAX i(VIM) FROM={40*T:.6g} TO={60*T:.6g}
meas tran irms RMS i(VIM) FROM={40*T:.6g} TO={60*T:.6g}
fourier {f:.6g} i(VIM)
wrdata {WORK}/tran_{int(f)}.txt i(VIM) v(a,b)
.endc
"""
    return s + ".end\n"


def ac_netlist(f, C, Lv=L, Rv=RC, tag="nom"):
    return f"""* AC sweep - series RLC channel at {f/1e3:g} kHz
VS a 0 AC 1
RON a n1 {RON}
C1 n1 n2 {C:.6g}
L1 n2 n3 {Lv:.6g}
RCOIL n3 n4 {Rv:.6g}
VIM n4 n5 0
RSENSE n5 0 {RSENSE}
.control
ac dec 400 1k 100k
wrdata {WORK}/ac_{tag}.txt mag(i(VIM))
.endc
.end
"""


def run(netlist, name):
    p = WORK / f"{name}.cir"
    p.write_text(netlist)
    out = subprocess.run(["ngspice", "-b", str(p)], capture_output=True, text=True, timeout=300)
    return out.stdout + out.stderr


def asc(f, C, phi):
    T = 1 / f
    return f"""Version 4
SHEET 1 1120 520
FLAG 96 176 a
FLAG 96 256 0
FLAG 256 176 a
FLAG 256 256 n1
FLAG 400 176 n1
FLAG 400 240 n2
FLAG 544 176 n2
FLAG 544 256 n3
FLAG 688 176 n3
FLAG 688 256 nb
FLAG 832 176 nb
FLAG 832 256 b
FLAG 960 176 b
FLAG 960 256 0
SYMBOL voltage 96 160 R0
SYMATTR InstName VA
SYMATTR Value PULSE(0 {{VRAIL}} 0 20n 20n {{0.5/F-20n}} {{1/F}})
SYMBOL res 240 160 R0
SYMATTR InstName RONA
SYMATTR Value {{RON/2}}
SYMBOL cap 384 176 R0
SYMATTR InstName C1
SYMATTR Value {{C}}
SYMBOL ind 528 160 R0
SYMATTR InstName L1
SYMATTR Value {{L}}
SYMATTR SpiceLine Rser={{RCOIL}}
SYMBOL res 672 160 R0
SYMATTR InstName RSENSE
SYMATTR Value {{RSENSE}}
SYMBOL res 816 160 R0
SYMATTR InstName RONB
SYMATTR Value {{RON/2}}
SYMBOL voltage 960 160 R0
SYMATTR InstName VB
SYMATTR Value PULSE(0 {{VRAIL}} {{PHI/360/F}} 20n 20n {{0.5/F-20n}} {{1/F}})
TEXT 80 336 Left 2 !.param F={f:.6g} VRAIL={VRAIL} PHI={phi:.4g} L={L:.6g} RCOIL={RC:.6g} C={C:.6g} RON={RON} RSENSE={RSENSE}
TEXT 80 368 Left 2 !.tran 0 {60*T:.6g} {40*T:.6g} 0.05u
TEXT 80 400 Left 2 !.meas tran Ipk MAX I(L1) FROM={40*T:.6g} TO={60*T:.6g}
TEXT 80 432 Left 2 !.four {{F}} I(L1)
TEXT 80 472 Left 2 ;CapsuleTrack-EM: phase-shifted full bridge (VA-VB) -> series-resonant patch coil, {f/1e3:g} kHz
"""


rows, harm_rows = [], []
for f in [14e3, 31.5e3]:
    C, phi = design(f)
    net = tran_netlist(f, C, phi)
    (HW / f"coil_driver_{int(f/1e3*10)/10:g}kHz.cir").write_text(tran_netlist(f, C, phi, control=False))
    (HW / f"coil_driver_{int(f/1e3*10)/10:g}kHz.asc").write_text(asc(f, C, phi))
    out = run(net, f"tran_{int(f)}")
    get = lambda k: float(re.search(rf"^{k}\s*=\s*([-+0-9.eE]+)", out, re.M | re.I).group(1))
    ipk, irms = get("ipk"), get("irms")
    pcoil = irms ** 2 * RC
    # harmonics table from ngspice fourier output
    harm = re.findall(r"^\s*(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+[-0-9.eE+]+\s+[-0-9.eE+]+", out, re.M)
    thd = re.search(r"THD:\s*([0-9.eE+-]+)", out)
    for h in harm:
        harm_rows.append(dict(f_hz=f, harmonic=int(h[0]), freq_hz=float(h[1]), magnitude_A=float(h[2])))
    rows.append(dict(f_hz=f, C_nF=C * 1e9, phase_shift_deg=phi, I_pk_A=ipk, I_rms_A=irms,
                     P_coil_mW=pcoil * 1e3, THD_pct=float(thd.group(1)) if thd else np.nan))
tran = pd.DataFrame(rows)
tran.to_csv(TAB / "spice_transient.csv", index=False)
pd.DataFrame(harm_rows).to_csv(TAB / "spice_harmonics.csv", index=False)

# waveform figure (14 kHz)
w = np.loadtxt(WORK / "tran_14000.txt")
t, i_l, v_ab = w[:, 0], w[:, 1], w[:, 3]
sel = t >= t[-1] - 3 / 14e3
fig, ax = plt.subplots(2, 1, figsize=(6.4, 4.2), sharex=True)
ax[0].plot((t[sel] - t[sel][0]) * 1e6, v_ab[sel], color=SERIES[1], lw=1.5)
ax[0].set_ylabel("V bridge (V)"); ax[0].set_title("14 kHz channel: 3-level bridge voltage and coil current")
ax[1].plot((t[sel] - t[sel][0]) * 1e6, i_l[sel] * 1e3, color=SERIES[0], lw=1.5)
ax[1].set_ylabel("Coil current (mA)"); ax[1].set_xlabel("Time (us)")
fig.savefig(FIG / "spice_waveform_14kHz.png"); plt.close(fig)

# AC resonance + tolerance Monte Carlo (C +/-5 %, L +/-5 %, coil R +0..20 % for heating)
f0 = 14e3
C0, _ = design(f0)
run(ac_netlist(f0, C0, tag="nom"), "ac_nom")
nom = np.loadtxt(WORK / "ac_nom.txt")
rng = np.random.default_rng(0)
tol = []
fig, ax = plt.subplots(figsize=(6.4, 3.8))
for k in range(40):
    Ck = C0 * (1 + rng.uniform(-0.05, 0.05)); Lk = L * (1 + rng.uniform(-0.05, 0.05))
    Rk = RC * (1 + rng.uniform(0, 0.20))
    run(ac_netlist(f0, Ck, Lk, Rk, tag=f"t{k}"), f"ac_t{k}")
    a = np.loadtxt(WORK / f"ac_t{k}.txt")
    i_f0 = np.interp(f0, a[:, 0], a[:, 1]); i_nom = np.interp(f0, nom[:, 0], nom[:, 1])
    tol.append(dict(run=k, C_nF=Ck * 1e9, L_uH=Lk * 1e6, R_ohm=Rk, I_rel_at_f0=i_f0 / i_nom))
    ax.plot(a[:, 0] / 1e3, a[:, 1], color="#b7d3f6", lw=0.8)
ax.plot(nom[:, 0] / 1e3, nom[:, 1], color=SERIES[0], lw=2, label="Nominal")
ax.axvline(14, color="#6b6a63", lw=1, ls="--")
ax.set_xscale("log"); ax.set_xlabel("Frequency (kHz)"); ax.set_ylabel("|I| per 1 V drive (A)")
ax.set_title("Low-Q (~1.3) resonance: tolerances shift current, not frequency")
ax.legend(fontsize=8)
fig.savefig(FIG / "spice_ac_tolerance.png"); plt.close(fig)
tol = pd.DataFrame(tol)
tol.to_csv(TAB / "spice_tolerance_mc.csv", index=False)

print(tran.round(4).to_string(index=False))
print("tolerance: current at f0 relative to nominal  min %.3f  max %.3f" % (tol.I_rel_at_f0.min(), tol.I_rel_at_f0.max()))
