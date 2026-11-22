"""Deliverable 3 - KiCad 7 design files, generated from the same coil parameters.

  hardware/kicad/patch/  4-layer flexible transmitter patch: 8 planar spiral coils
                         (geometry from capsuletrack.coils.SpiralCoil), inter-layer vias,
                         per-coil 2-pad lead connector, board outline.
  hardware/kicad/driver/ block-level driver schematic: MCU, 8 H-bridge channels with
                         series-resonant capacitor, current sensing, power and debug.
  hardware/kicad/bom.csv bill of materials.

Design notes
  * Each coil is one continuous conductor F.Cu (in) -> In1 (out) -> In2 (in) -> B.Cu (out);
    all four layers circulate in the same sense, so the turns add.
  * The two lead pads of a coil therefore share one net (COILn) - KiCad has no way to
    represent the series inductance inside copper without a net-tie.
  * Odd turn counts per layer (n + 1/4, n + 1/8) put every via at a different angle, so
    through-vias never land on another layer's trace.
"""
import sys
import uuid
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from capsuletrack.coils import SpiralCoil  # noqa: E402

KI = ROOT / "hardware" / "kicad"
coil = SpiralCoil()
U = lambda: str(uuid.uuid4())

# ----------------------------------------------------------------------------- PCB
PITCH, NCOL, NROW = 60.0, 4, 2
W, H = 260.0, 150.0
X0, Y0 = 20.0, 20.0                      # board origin on the sheet
CX = X0 + W / 2 - (NCOL - 1) / 2 * PITCH - 12
CY = Y0 + H / 2 - (NROW - 1) / 2 * PITCH
TW = coil.trace_w * 1e3
P = coil.pitch * 1e3
ROUT = coil.r_out * 1e3
LAYERS = ["F.Cu", "In1.Cu", "In2.Cu", "B.Cu"]
VIA = (0.6, 0.3)
STUB = 0.8


def pt(cx, cy, r, th):
    return cx + r * np.cos(th), cy - r * np.sin(th)


segs, vias, fps = [], [], []
nets = ["COIL%d" % (i + 1) for i in range(NCOL * NROW)]


def seg(a, b, layer, net):
    segs.append(f'  (segment (start {a[0]:.4f} {a[1]:.4f}) (end {b[0]:.4f} {b[1]:.4f}) '
                f'(width {TW}) (layer "{layer}") (net {net}) (uuid "{U()}"))')


def via(a, net):
    vias.append(f'  (via (at {a[0]:.4f} {a[1]:.4f}) (size {VIA[0]}) (drill {VIA[1]}) '
                f'(layers "F.Cu" "B.Cu") (net {net}) (uuid "{U()}"))')


def spiral(cx, cy, r0, th0, turns, inward, layer, net, step=np.radians(6)):
    n = max(2, int(np.ceil(2 * np.pi * turns / step)))
    th = np.linspace(th0, th0 + 2 * np.pi * turns, n + 1)
    dr = -P if inward else P
    r = r0 + dr * (th - th0) / (2 * np.pi)
    xy = [pt(cx, cy, ri, ti) for ri, ti in zip(r, th)]
    for a, b in zip(xy[:-1], xy[1:]):
        seg(a, b, layer, net)
    return r[-1], th[-1]


for k in range(NCOL * NROW):
    row, col = divmod(k, NCOL)
    cx, cy = CX + col * PITCH, CY + row * PITCH
    net = k + 1
    r, th = ROUT, 0.0
    ends = []
    for li, (layer, turns) in enumerate(zip(LAYERS, coil.turns_per_layer)):
        inward = li % 2 == 0
        r_start, th_start = r, th
        r, th = spiral(cx, cy, r, th, turns, inward, layer, net)
        # stub to the transition via (inside the coil after an inward layer, outside after outward)
        rv = r - STUB if inward else r + STUB
        a, b = pt(cx, cy, r, th), pt(cx, cy, rv, th)
        seg(a, b, layer, net)
        if li < 3:
            via(b, net)
            seg(b, a, LAYERS[li + 1], net)   # next layer starts from the same via
        ends.append((rv, th))
    # terminal A: F.Cu start (r_out, 0) -> pad 1 ; terminal B: final via (B.Cu) -> F.Cu -> pad 2
    rv, th = ends[-1]
    vb = pt(cx, cy, rv, th)
    via(vb, net)
    jx, jy = cx + ROUT + 4.6, cy
    p1, p2 = (jx, jy - 1.25), (jx, jy + 1.25)
    seg(pt(cx, cy, ROUT, 0.0), (cx + ROUT + 1.2, cy), "F.Cu", net)
    seg((cx + ROUT + 1.2, cy), p1, "F.Cu", net)
    seg(vb, p2, "F.Cu", net)
    fps.append(f'''  (footprint "CapsuleTrack:CoilLead_2Pad" (layer "F.Cu") (uuid "{U()}") (at {jx:.4f} {jy:.4f})
    (attr smd)
    (fp_text reference "J{k+1}" (at 2.6 0 90) (layer "F.SilkS") (uuid "{U()}") (effects (font (size 0.8 0.8) (thickness 0.12))))
    (fp_text value "COIL{k+1}_LEADS" (at -2.4 0 90) (layer "F.Fab") (uuid "{U()}") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (fp_rect (start -0.9 -2.3) (end 0.9 2.3) (stroke (width 0.05) (type default)) (fill none) (layer "F.CrtYd") (uuid "{U()}"))
    (pad "1" smd rect (at 0 -1.25) (size 1.2 1.5) (layers "F.Cu" "F.Paste" "F.Mask") (net {net} "{nets[k]}") (uuid "{U()}"))
    (pad "2" smd rect (at 0 1.25) (size 1.2 1.5) (layers "F.Cu" "F.Paste" "F.Mask") (net {net} "{nets[k]}") (uuid "{U()}"))
  )''')

texts = [f'  (gr_text "COIL {k+1}  {14 + 2.5*k:g} kHz" (at {CX + (k % NCOL)*PITCH:.2f} {CY + (k//NCOL)*PITCH + ROUT + 4:.2f}) '
         f'(layer "F.SilkS") (uuid "{U()}") (effects (font (size 1.2 1.2) (thickness 0.18))))'
         for k in range(NCOL * NROW)]
texts.append(f'  (gr_text "CapsuleTrack-EM transmitter patch  |  4-layer flex  |  8 x 65-turn planar coils" '
             f'(at {X0 + W/2:.1f} {Y0 + H - 6:.1f}) (layer "F.SilkS") (uuid "{U()}") (effects (font (size 1.5 1.5) (thickness 0.2))))')

net_decl = '  (net 0 "")\n' + "\n".join(f'  (net {i+1} "{n}")' for i, n in enumerate(nets))
pcb = f'''(kicad_pcb (version 20221018) (generator pcbnew)
  (general (thickness 0.31))
  (paper "A3")
  (title_block (title "CapsuleTrack-EM transmitter patch") (rev "A") (company "S. Dhole - research prototype, not fabricated"))
  (layers
    (0 "F.Cu" signal) (1 "In1.Cu" signal) (2 "In2.Cu" signal) (31 "B.Cu" signal)
    (32 "B.Adhes" user "B.Adhesive") (33 "F.Adhes" user "F.Adhesive")
    (34 "B.Paste" user) (35 "F.Paste" user)
    (36 "B.SilkS" user "B.Silkscreen") (37 "F.SilkS" user "F.Silkscreen")
    (38 "B.Mask" user) (39 "F.Mask" user)
    (40 "Dwgs.User" user "User.Drawings") (41 "Cmts.User" user "User.Comments")
    (44 "Edge.Cuts" user) (45 "Margin" user)
    (46 "B.CrtYd" user "B.Courtyard") (47 "F.CrtYd" user "F.Courtyard")
    (48 "B.Fab" user) (49 "F.Fab" user)
  )
  (setup
    (stackup
      (layer "F.SilkS" (type "Top Silk Screen"))
      (layer "F.Mask" (type "Top Solder Mask") (thickness 0.025))
      (layer "F.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 1" (type "core") (thickness 0.05) (material "Polyimide") (epsilon_r 3.3) (loss_tangent 0.002))
      (layer "In1.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 2" (type "prepreg") (thickness 0.05) (material "Polyimide") (epsilon_r 3.3) (loss_tangent 0.002))
      (layer "In2.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 3" (type "core") (thickness 0.05) (material "Polyimide") (epsilon_r 3.3) (loss_tangent 0.002))
      (layer "B.Cu" (type "copper") (thickness 0.035))
      (layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.025))
      (layer "B.SilkS" (type "Bottom Silk Screen"))
      (copper_finish "ENIG")
      (dielectric_constraints no)
    )
    (pad_to_mask_clearance 0)
  )
{net_decl}
  (gr_rect (start {X0} {Y0}) (end {X0 + W} {Y0 + H}) (stroke (width 0.15) (type default)) (fill none) (layer "Edge.Cuts") (uuid "{U()}"))
{chr(10).join(texts)}
{chr(10).join(fps)}
{chr(10).join(segs)}
{chr(10).join(vias)}
)
'''
(KI / "patch" / "capsuletrack_patch.kicad_pcb").write_text(pcb.replace("(uuid ", "(tstamp "))
(KI / "patch" / "capsuletrack_patch.kicad_pro").write_text('{"board": {"design_settings": {"rules": {"min_clearance": 0.15, "min_track_width": 0.15, "min_via_diameter": 0.5, "min_through_hole_diameter": 0.25}}}, "meta": {"filename": "capsuletrack_patch.kicad_pro", "version": 1}}\n')
pretty = KI / "patch" / "CapsuleTrack.pretty"
pretty.mkdir(exist_ok=True)
(pretty / "CoilLead_2Pad.kicad_mod").write_text("""(footprint "CoilLead_2Pad" (version 20221018) (generator pcbnew)
  (layer "F.Cu")
  (attr smd)
  (fp_text reference "REF**" (at 2.6 0 90) (layer "F.SilkS") (effects (font (size 0.8 0.8) (thickness 0.12))))
  (fp_text value "CoilLead_2Pad" (at -2.4 0 90) (layer "F.Fab") (effects (font (size 0.6 0.6) (thickness 0.1))))
  (fp_rect (start -0.9 -2.3) (end 0.9 2.3) (stroke (width 0.05) (type default)) (fill none) (layer "F.CrtYd"))
  (pad "1" smd rect (at 0 -1.25) (size 1.2 1.5) (layers "F.Cu" "F.Paste" "F.Mask"))
  (pad "2" smd rect (at 0 1.25) (size 1.2 1.5) (layers "F.Cu" "F.Paste" "F.Mask"))
)
""")
(KI / "patch" / "fp-lib-table").write_text('(fp_lib_table\n  (lib (name "CapsuleTrack")(type "KiCad")(uri "${KIPRJMOD}/CapsuleTrack.pretty")(options "")(descr "CapsuleTrack-EM footprints"))\n)\n')
print("PCB:", len(segs), "segments,", len(vias), "vias")

# ----------------------------------------------------------------------- schematic
LIB = {}


def lib_symbol(name, pins, w=10.16, ref="U", body=True):
    """pins: list of (number, name, side) with side in L/R. Returns lib symbol + pin map."""
    left = [p for p in pins if p[2] == "L"]
    right = [p for p in pins if p[2] == "R"]
    n = max(len(left), len(right), 1)
    h = 2.54 * (n + 1)
    pin_xy, s = {}, []
    for i, (num, nm, _) in enumerate(left):
        y = h / 2 - 2.54 * (i + 1)
        pin_xy[num] = (-w / 2 - 2.54, y)
        s.append(f'(pin passive line (at {-w/2 - 2.54:.2f} {y:.2f} 0) (length 2.54) (name "{nm}" (effects (font (size 1.0 1.0)))) (number "{num}" (effects (font (size 1.0 1.0)))))')
    for i, (num, nm, _) in enumerate(right):
        y = h / 2 - 2.54 * (i + 1)
        pin_xy[num] = (w / 2 + 2.54, y)
        s.append(f'(pin passive line (at {w/2 + 2.54:.2f} {y:.2f} 180) (length 2.54) (name "{nm}" (effects (font (size 1.0 1.0)))) (number "{num}" (effects (font (size 1.0 1.0)))))')
    txt = f'''(symbol "CT:{name}" (in_bom yes) (on_board yes)
      (property "Reference" "{ref}" (at 0 {h/2 + 1.5:.2f} 0) (effects (font (size 1.27 1.27))))
      (property "Value" "{name}" (at 0 {-h/2 - 1.5:.2f} 0) (effects (font (size 1.27 1.27))))
      (property "Footprint" "" (at 0 0 0) (effects (font (size 1.27 1.27)) hide))
      (property "Datasheet" "" (at 0 0 0) (effects (font (size 1.27 1.27)) hide))
      (symbol "{name}_0_1" (rectangle (start {-w/2:.2f} {h/2:.2f}) (end {w/2:.2f} {-h/2:.2f}) (stroke (width 0.254) (type default)) (fill (type background))))
      (symbol "{name}_1_1" {" ".join(s)}))'''
    LIB[name] = (txt, pin_xy)


lib_symbol("R", [("1", "~", "L"), ("2", "~", "R")], w=5.08, ref="R")
lib_symbol("C", [("1", "~", "L"), ("2", "~", "R")], w=5.08, ref="C")
lib_symbol("Conn_01x02", [("1", "Pin_1", "R"), ("2", "Pin_2", "R")], w=5.08, ref="J")
lib_symbol("Conn_SWD_5", [("1", "VTREF", "R"), ("2", "SWDIO", "R"), ("3", "SWCLK", "R"), ("4", "NRST", "R"), ("5", "GND", "R")], w=10.16, ref="J")
lib_symbol("DRV8837", [("8", "VCC", "L"), ("6", "IN1", "L"), ("5", "IN2", "L"), ("7", "nSLEEP", "L"),
                       ("1", "VM", "R"), ("2", "OUT1", "R"), ("3", "OUT2", "R"), ("4", "GND", "R")], w=15.24)
lib_symbol("CurrentSenseAmp", [("1", "IN+", "L"), ("2", "IN-", "L"), ("3", "REF", "L"),
                               ("4", "VS", "R"), ("5", "OUT", "R"), ("6", "GND", "R")], w=15.24)
lib_symbol("Buck_Reg", [("1", "VIN", "L"), ("2", "EN", "L"), ("3", "VOUT", "R"), ("4", "GND", "R")], w=10.16)
lib_symbol("LDO_3V3", [("1", "VIN", "L"), ("2", "VOUT", "R"), ("3", "GND", "R")], w=10.16)
mcu_pins = [(str(i + 1), f"PWM{i+1}A", "R") for i in range(8)] + [(str(i + 9), f"PWM{i+1}B", "R") for i in range(8)]
mcu_pins += [(str(i + 17), f"ISNS{i+1}", "L") for i in range(8)]
mcu_pins += [("25", "DRV_EN", "L"), ("26", "VDD", "L"), ("27", "VSS", "L"), ("28", "SWDIO", "L"),
             ("29", "SWCLK", "L"), ("30", "NRST", "L"), ("31", "BUCK_EN", "L"), ("32", "VREF_IN", "L")]
lib_symbol("MCU_STM32G4_class", mcu_pins, w=20.32)

SHEET = U()
placed, labels = [], []
refcount = {}


def place(sym, x, y, value, conn, ref_prefix=None, footprint=""):
    """conn: dict pin_number -> net label name."""
    txt, pin_xy = LIB[sym]
    x, y = round((x + 30) / 2.54) * 2.54, round((y + 25) / 2.54) * 2.54
    small = len(pin_xy) <= 2
    dy = 5.08 if small else 8
    ref_prefix = ref_prefix or {"R": "R", "C": "C", "Conn_01x02": "J", "Conn_SWD_5": "J"}.get(sym, "U")
    refcount[ref_prefix] = refcount.get(ref_prefix, 0) + 1
    ref = f"{ref_prefix}{refcount[ref_prefix]}"
    pins = " ".join(f'(pin "{n}" (uuid "{U()}"))' for n in pin_xy)
    placed.append(f'''  (symbol (lib_id "CT:{sym}") (at {x:.2f} {y:.2f} 0) (unit 1) (in_bom yes) (on_board yes) (dnp no) (uuid "{U()}")
    (property "Reference" "{ref}" (at {x:.2f} {y - dy:.2f} 0) (effects (font (size 1.27 1.27))))
    (property "Value" "{value}" (at {x:.2f} {y + dy:.2f} 0) (effects (font (size 1.0 1.0))))
    (property "Footprint" "{footprint}" (at {x:.2f} {y:.2f} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (at {x:.2f} {y:.2f} 0) (effects (font (size 1.27 1.27)) hide))
    {pins}
    (instances (project "capsuletrack_driver" (path "/{SHEET}" (reference "{ref}") (unit 1)))))''')
    for num, net in conn.items():
        px, py = pin_xy[num]
        lx, ly = x + px, y - py
        ang = 180 if px < 0 else 0
        just = "right" if px < 0 else "left"
        labels.append(f'  (label "{net}" (at {lx:.2f} {ly:.2f} {ang}) (fields_autoplaced) (effects (font (size 1.0 1.0)) (justify {just})) (uuid "{U()}"))')
    return ref


bom = []
# MCU + support
mcu_conn = {**{str(i + 1): f"PWM{i+1}_A" for i in range(8)}, **{str(i + 9): f"PWM{i+1}_B" for i in range(8)},
            **{str(i + 17): f"ISNS{i+1}" for i in range(8)},
            "25": "DRV_EN", "26": "+3V3", "27": "GND", "28": "SWDIO", "29": "SWCLK", "30": "NRST",
            "31": "BUCK_EN", "32": "VREF"}
place("MCU_STM32G4_class", 70, 130, "STM32G474RE (example)", mcu_conn)
bom.append(("U1", "STM32G474RET6 (example MCU: HRTIM + TIM1/TIM8 for >=16 PWM outputs, 5 ADCs)", 1, "LQFP-64", "STMicroelectronics"))
place("LDO_3V3", 40, 30, "3.3 V LDO", {"1": "VBAT", "2": "+3V3", "3": "GND"})
bom.append(("U2", "3.3 V low-dropout regulator, >=150 mA", 1, "SOT-23-5", "e.g. TI TLV75533"))
place("Buck_Reg", 40, 60, "3.3 V buck, driver rail", {"1": "VBAT", "2": "BUCK_EN", "3": "VDRV", "4": "GND"})
bom.append(("U3", "Step-down converter for H-bridge rail VDRV (3.3 V, >=600 mA)", 1, "SOT-563", "e.g. TI TPS62A01"))
place("Conn_01x02", 5, 30, "Battery 3.7 V Li-Po", {"1": "VBAT", "2": "GND"})
place("Conn_SWD_5", 15, 130, "SWD debug", {"1": "+3V3", "2": "SWDIO", "3": "SWCLK", "4": "NRST", "5": "GND"})
for j, net in enumerate(["VBAT", "+3V3", "VDRV", "VDRV", "+3V3"]):
    place("C", 75, 20 + 14 * j, f"10uF ({net})", {"1": net, "2": "GND"})
    bom.append((f"C{refcount['C']}", f"10 uF X7R decoupling on {net}", 1, "0603", "generic"))
place("R", 105, 20, "10k", {"1": "+3V3", "2": "VREF"}); place("R", 105, 34, "10k", {"1": "VREF", "2": "GND"})
place("R", 105, 48, "10k", {"1": "DRV_EN", "2": "GND"})
bom.append(("R1-R3", "10 kOhm 1% (VREF divider, DRV_EN pull-down)", 3, "0402", "generic"))

# 8 driver channels
caps = [636.25, 458.05, 345.44, 269.78, 216.50, 177.58, 148.28, 125.68]
for i in range(8):
    ch = i + 1
    x0 = 150 + (i % 4) * 85
    y0 = 40 + (i // 4) * 100
    place("DRV8837", x0, y0, "DRV8837", {"8": "+3V3", "6": f"PWM{ch}_A", "5": f"PWM{ch}_B", "7": "DRV_EN",
                                          "1": "VDRV", "2": f"CH{ch}_OUT1", "3": f"CH{ch}_OUT2", "4": "GND"})
    place("C", x0 + 22, y0 - 14, f"{caps[i]:.0f}nF C0G", {"1": f"CH{ch}_OUT1", "2": f"CH{ch}_C"})
    place("Conn_01x02", x0 + 44, y0 - 14, f"Coil {ch} ({14 + 2.5*i:g} kHz)", {"1": f"CH{ch}_C", "2": f"CH{ch}_RET"})
    place("R", x0 + 22, y0 + 12, "0.5R 1%", {"1": f"CH{ch}_RET", "2": f"CH{ch}_OUT2"})
    place("CurrentSenseAmp", x0 + 15, y0 + 30, "INA240-class", {"1": f"CH{ch}_RET", "2": f"CH{ch}_OUT2", "3": "VREF",
                                                                "4": "+3V3", "5": f"ISNS{ch}", "6": "GND"})
    place("C", x0 + 40, y0 + 30, "100nF", {"1": "VDRV", "2": "GND"})
bom += [("U4-U11", "DRV8837 H-bridge (one per coil channel)", 8, "WSON-8 2x2", "Texas Instruments"),
        ("Cres1-8", "Series-resonant capacitor, C0G/NP0, 50 V, values 636/458/345/270/217/178/148/126 nF (trim per coil)", 8, "1206", "e.g. Murata GRM C0G"),
        ("Rs1-8", "0.5 Ohm 1% current-sense resistor", 8, "0805", "generic"),
        ("U12-U19", "Bidirectional current-sense amplifier, PWM-rejecting (INA240-class)", 8, "TSSOP-8", "Texas Instruments"),
        ("C 100n x8", "100 nF X7R bulk at each H-bridge VM pin", 8, "0402", "generic"),
        ("Jcoil1-8", "2-pin lead connector to patch coil (twisted pair harness)", 8, "JST SH 1.0 mm", "JST"),
        ("PATCH", "4-layer polyimide flex patch, 8 x 65-turn planar coils (hardware/kicad/patch)", 1, "260 x 150 mm", "flex PCB fab"),
        ("BAT", "3.7 V Li-Po pouch cell, >=1500 mAh (12 h at nominal 80 mA drive)", 1, "", "generic"),
        ("RX", "Capsule receiver: 3-axis ferrite-cored coil, 400 turns, 2.5 mm winding diameter (concept)", 1, "", "custom")]

root = U()
sch = f'''(kicad_sch (version 20230121) (generator eeschema)
  (uuid "{SHEET}")
  (paper "A2")
  (title_block (title "CapsuleTrack-EM patch driver - block-level schematic") (rev "A")
    (comment 1 "Design intent: verify IC pinouts vs datasheets before layout")
    (comment 2 "Channel: H-bridge -> C0G series C -> coil -> 0.5R sense"))
  (lib_symbols
    {chr(10).join(v[0] for v in LIB.values())}
  )
{chr(10).join(labels)}
{chr(10).join(placed)}
  (sheet_instances (path "/" (page "1")))
)
'''
(KI / "driver" / "capsuletrack_driver.kicad_sch").write_text(sch)
(KI / "driver" / "capsuletrack_driver.kicad_pro").write_text('{"meta": {"filename": "capsuletrack_driver.kicad_pro", "version": 1}}\n')
(KI / "driver" / "sym-lib-table").write_text("(sym_lib_table\n)\n")

import csv  # noqa: E402
with open(KI / "bom.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Reference", "Description", "Qty", "Package", "Manufacturer / example part"])
    w.writerows(bom)
print("schematic:", len(placed), "symbols,", len(labels), "labels")
