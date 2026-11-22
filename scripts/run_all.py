"""Reproduce every table, figure and hardware file in order."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for s in ["01_tissue_analysis.py", "02_coil_design.py", "03_driver_spice.py",
          "04_montecarlo_main.py", "05_sweeps.py", "06_kicad.py"]:
    print(f"\n=== {s} ===", flush=True)
    subprocess.run([sys.executable, str(ROOT / "scripts" / s)], check=True)
