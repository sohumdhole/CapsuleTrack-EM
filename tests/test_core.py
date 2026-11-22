import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from capsuletrack import tissue as T  # noqa: E402
from capsuletrack.anatomy import Body, load_ansur  # noqa: E402
from capsuletrack.coils import SpiralCoil, mutual_coaxial  # noqa: E402
from capsuletrack.field import coil_field, dipole_field, MU0  # noqa: E402
from capsuletrack.system import Patch, Receiver, default_freqs, measure  # noqa: E402
from capsuletrack import solver as S  # noqa: E402


def test_loop_on_axis_matches_analytic():
    a, z = 0.02, 0.07
    B = coil_field(np.array([[0, 0, z]]), [0, 0, 0], [0, 0, 1], np.array([a]), np.array([0.0]))
    assert np.isclose(B[0, 2], MU0 * a ** 2 / (2 * (a ** 2 + z ** 2) ** 1.5), rtol=1e-9)
    assert abs(B[0, 0]) < 1e-15 and abs(B[0, 1]) < 1e-15


def test_loop_far_field_matches_dipole():
    p = np.array([[0.3, -0.2, 0.4]])
    B = coil_field(p, [0, 0, 0], [0, 0, 1], np.array([0.01]), np.array([0.0]))
    D = dipole_field(p, [0, 0, 0], [0, 0, np.pi * 0.01 ** 2])
    assert np.allclose(B, D, rtol=1e-3)


def test_mutual_inductance_symmetry_and_limit():
    assert np.isclose(mutual_coaxial(0.02, 0.015, 0.003), mutual_coaxial(0.015, 0.02, 0.003))
    # far apart -> dipole-dipole coupling mu0 pi a^2 b^2 / (2 d^3)
    a, b, d = 0.01, 0.01, 1.0
    assert np.isclose(mutual_coaxial(a, b, d), MU0 * np.pi * a ** 2 * b ** 2 / (2 * d ** 3), rtol=1e-3)


def test_muscle_conductivity_10khz():
    p = T.load_params().loc["Muscle"]
    assert 0.30 < T.conductivity(np.array([1e4]), p)[0] < 0.38


def test_quasistatic_regime():
    p = T.load_params().loc["Small Intestine"]
    f = np.array([50e3])
    d = T.skin_depth(f, T.conductivity(f, p), T.rel_permittivity(f, p))
    assert d[0] > 2.0          # metres >> abdomen


def test_coil_parameters_reasonable():
    s = SpiralCoil().summary()
    assert 100 < s["L_uH"] < 400 and 5 < s["R_dc_ohm"] < 30


def test_ansur_loaded():
    d = load_ansur()
    assert (d.Gender == "Male").sum() == 4082 and (d.Gender == "Female").sum() == 1986


def test_solver_recovers_noiseless_pose():
    d = load_ansur()
    body = Body.from_row(d.iloc[10])
    C, N = body.patch_coils()
    patch = Patch(C, N, default_freqs(len(C)))
    rng = np.random.default_rng(1)
    p = body.sample_workspace(1, rng)[0]
    R = Rotation.random(random_state=3).as_matrix()
    b = measure(patch, p, R, Receiver())
    sig = np.full(b.shape, 1e-9)
    pe, Re, _ = S.solve(patch, b, sig, body.workspace_grid(6))
    assert np.linalg.norm(pe - p) < 1e-5
    assert S.orientation_error_deg(Re, R) < 1e-3
