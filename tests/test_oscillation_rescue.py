"""Tests for the period-2 oscillation rescue in the minimum-pressure iteration.

The problematic/control profiles are real adjacent-day soundings (1995) copied
verbatim from notebooks/illustrate_numerical_instability.ipynb. Before the rescue,
the problematic profile locked into a period-2 pressure cycle (~961.0 <-> ~961.7
hPa, amplitude above the 0.5 hPa tolerance) and returned NaN with IFL=2; the
control profile converged normally. The rescue collapses a detected 2-cycle to its
midpoint so the loop exits through the normal convergence path, while leaving all
previously-converging profiles bit-identical.
"""

import numpy as np

import tcpyPI

PI_OPTIONS = dict(CKCD=0.9, ascent_flag=0, diss_flag=1, ptop=50, miss_handle=1)

# fmt: off
PROFILE_PROBLEMATIC = np.array([
    [1000.0, 20.266, 9.65], [975.0, 18.11, 9.517], [950.0, 15.978, 9.459],
    [925.0, 13.94, 9.205], [900.0, 12.576, 7.784], [875.0, 11.213, 6.561],
    [850.0, 9.669, 5.952], [825.0, 8.084, 5.517], [800.0, 6.516, 4.976],
    [775.0, 4.993, 4.332], [750.0, 4.374, 3.102], [700.0, 4.93, 1.207],
    [650.0, 1.852, 1.521], [600.0, -2.047, 1.263], [550.0, -6.917, 0.857],
    [500.0, -13.015, 0.556], [450.0, -19.698, 0.33], [400.0, -26.881, 0.211],
    [350.0, -34.178, 0.104], [300.0, -42.822, 0.073], [250.0, -51.573, 0.037],
    [225.0, -54.68, 0.027], [200.0, -57.207, 0.018], [175.0, -58.943, 0.01],
    [150.0, -60.751, 0.005], [125.0, -63.255, 0.004], [100.0, -63.277, 0.003],
    [70.0, -61.5, 0.003], [50.0, -58.453, 0.003], [30.0, -53.327, 0.003],
    [20.0, -50.156, 0.003], [10.0, -46.602, 0.003], [7.0, -42.283, 0.003],
    [5.0, -35.389, 0.003], [3.0, -22.349, 0.003], [2.0, -14.944, 0.004],
    [1.0, -11.886, 0.004],
])
SST_PROBLEMATIC, MSL_PROBLEMATIC = 25.267, 1012.9

PROFILE_CONTROL = np.array([
    [1000.0, 21.326, 11.215], [975.0, 19.192, 11.011], [950.0, 17.315, 10.556],
    [925.0, 15.676, 9.931], [900.0, 14.286, 8.639], [875.0, 12.957, 7.58],
    [850.0, 11.495, 6.97], [825.0, 10.019, 6.439], [800.0, 8.575, 5.734],
    [775.0, 6.993, 5.198], [750.0, 5.373, 4.752], [700.0, 3.13, 3.122],
    [650.0, 0.534, 2.562], [600.0, -2.573, 1.442], [550.0, -6.616, 0.865],
    [500.0, -12.233, 0.61], [450.0, -18.954, 0.42], [400.0, -26.098, 0.292],
    [350.0, -33.337, 0.184], [300.0, -40.639, 0.112], [250.0, -49.373, 0.057],
    [225.0, -52.236, 0.03], [200.0, -55.706, 0.017], [175.0, -57.323, 0.01],
    [150.0, -60.084, 0.005], [125.0, -65.224, 0.004], [100.0, -65.651, 0.003],
    [70.0, -63.457, 0.003], [50.0, -59.522, 0.003], [30.0, -53.504, 0.003],
    [20.0, -49.926, 0.003], [10.0, -45.925, 0.003], [7.0, -41.476, 0.003],
    [5.0, -34.68, 0.003], [3.0, -23.495, 0.003], [2.0, -16.575, 0.004],
    [1.0, -12.898, 0.004],
])
SST_CONTROL, MSL_CONTROL = 25.292, 1011.4
# fmt: on


def _run(profile, sst, msl):
    P, TC, R = profile[:, 0], profile[:, 1], profile[:, 2]
    return tcpyPI.pi(sst, msl, P, TC, R, **PI_OPTIONS)


def test_problematic_1995_profile_rescued():
    # Previously: period-2 oscillation -> NaN with IFL=2. The rescue collapses the
    # cycle to its midpoint and the solver converges through the normal path.
    vmax, pmin, ifl, to, otl = _run(PROFILE_PROBLEMATIC, SST_PROBLEMATIC, MSL_PROBLEMATIC)
    assert ifl == 1
    assert np.isfinite(vmax) and 60.0 < vmax < 90.0
    assert np.isfinite(pmin) and 925.0 < pmin < 950.0
    assert np.isfinite(to) and np.isfinite(otl)


def test_control_1995_profile_bit_identical():
    # The adjacent-day control profile converged before the rescue existed; its
    # outputs are pinned here (values captured pre-rescue) to guard bit-identity
    # for converging profiles.
    vmax, pmin, ifl, to, otl = _run(PROFILE_CONTROL, SST_CONTROL, MSL_CONTROL)
    assert ifl == 1
    np.testing.assert_allclose(vmax, 68.98260512391965, rtol=1e-13)
    np.testing.assert_allclose(pmin, 943.0252390748838, rtol=1e-13)
    np.testing.assert_allclose(to, 214.15121722308447, rtol=1e-13)
    np.testing.assert_allclose(otl, 159.82630589536842, rtol=1e-13)


def test_fuzz_marginal_profiles():
    # Stress test around the marginal regime: perturb the problematic sounding and
    # require self-consistent outputs — valid flags, and never NaN-with-success.
    rng = np.random.default_rng(0)
    P = PROFILE_PROBLEMATIC[:, 0]
    for _ in range(300):
        tc = PROFILE_PROBLEMATIC[:, 1] + rng.uniform(-0.5, 0.5, P.size)
        r = PROFILE_PROBLEMATIC[:, 2] * rng.uniform(0.8, 1.2, P.size)
        sst = SST_PROBLEMATIC + rng.uniform(-0.5, 0.5)
        msl = MSL_PROBLEMATIC + rng.uniform(-5.0, 5.0)
        vmax, pmin, ifl, to, otl = tcpyPI.pi(sst, msl, P, tc, r, **PI_OPTIONS)
        assert ifl in (0, 1, 2, 3)
        if ifl == 1:
            assert np.isfinite(vmax) and np.isfinite(pmin)
            assert np.isfinite(to) and np.isfinite(otl)
        else:
            assert np.isnan(vmax)
