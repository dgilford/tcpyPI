"""Convergence tests for marginal (near-neutral) soundings — GitHub issue #77.

Under pcmin.m's legacy CAPE rule (integrate signed buoyancy to the highest
positive grid level, clamp at zero) these soundings made the minimum-pressure
map discontinuous, and the iteration locked into a period-2 cycle (IFL=2 after
the iteration cap). With CAPE defined as the maximum running buoyancy work
(max-W; see cape()), the map is continuous and all of them converge through
the normal path. Full analysis: discontinuity_analysis/.

Profiles: the issue-77 sounding (Hurricane Francine's pre-landfall environment
in a -1.6 K counterfactual-SST sensitivity member, from PR #77's failing test)
and two real adjacent-day 1995 soundings copied verbatim from
notebooks/illustrate_numerical_instability.ipynb (the marginal one previously
cycled at ~961.0 <-> ~961.7 hPa; the control always converged).
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


def test_problematic_1995_profile_converges():
    # Previously: period-2 oscillation -> NaN with IFL=2 (then midpoint-rescued).
    # Under max-W CAPE the pressure map is continuous and the solver converges
    # through the normal path. Values pinned as a regression guard.
    vmax, pmin, ifl, to, otl = _run(PROFILE_PROBLEMATIC, SST_PROBLEMATIC, MSL_PROBLEMATIC)
    assert ifl == 1
    np.testing.assert_allclose(vmax, 74.37826294760555, rtol=1e-12)
    np.testing.assert_allclose(pmin, 936.9133348865204, rtol=1e-12)
    np.testing.assert_allclose(to, 212.78200681976475, rtol=1e-12)
    np.testing.assert_allclose(otl, 155.29600137949106, rtol=1e-12)


def test_control_1995_profile_converges():
    # The adjacent-day control profile always converged. max-W CAPE moves its
    # fixed point by ~0.004 hPa (legacy pins: VMAX=68.98260512, PMIN=943.02523907);
    # current values pinned as a regression guard.
    vmax, pmin, ifl, to, otl = _run(PROFILE_CONTROL, SST_CONTROL, MSL_CONTROL)
    assert ifl == 1
    np.testing.assert_allclose(vmax, 68.98207820269134, rtol=1e-12)
    np.testing.assert_allclose(pmin, 943.0291251377307, rtol=1e-12)
    np.testing.assert_allclose(to, 214.1513826661174, rtol=1e-12)
    np.testing.assert_allclose(otl, 159.8278039307987, rtol=1e-12)


# fmt: off
PROFILE_FRANCINE = np.array([
    [1000.0, 25.260956, 1.0783693e01], [975.0, 23.078949, 1.0704287e01],
    [950.0, 20.881561, 1.0680210e01], [925.0, 18.652649, 1.0617845e01],
    [900.0, 16.615143, 1.0320683e01], [875.0, 14.63147, 9.8811483e00],
    [850.0, 12.794586, 9.1884289e00], [825.0, 11.79306, 7.1884680e00],
    [800.0, 11.01236, 5.6963191e00], [775.0, 10.847565, 3.5568204e00],
    [750.0, 10.354492, 1.5912720e00], [700.0, 7.807007, 1.0433695e00],
    [650.0, 5.473297, 5.9723043e-01], [600.0, 2.5278625, 4.3974420e-01],
    [550.0, -1.532135, 4.8722979e-01], [500.0, -7.3142395, 5.8590513e-01],
    [450.0, -13.635345, 4.5599860e-01], [400.0, -20.613602, 3.1293562e-01],
    [350.0, -28.588928, 1.9222400e-01], [300.0, -37.270096, 9.7611703e-02],
    [250.0, -45.40825, 3.3851895e-02], [225.0, -50.047455, 2.4188591e-02],
    [200.0, -54.84575, 1.9636340e-02], [175.0, -59.173737, 1.3214327e-02],
    [150.0, -62.878662, 7.2453087e-03], [125.0, -65.78009, 4.3027173e-03],
    [100.0, -69.30669, 3.7014042e-03], [70.0, -65.512024, 2.9414182e-03],
    [50.0, -60.76361, 2.7806845e-03], [30.0, -55.100464, 2.8306348e-03],
    [20.0, -51.38333, 2.9053832e-03], [10.0, -43.95462, 2.9848625e-03],
    [7.0, -39.709717, 3.0779461e-03], [5.0, -33.193268, 3.1315640e-03],
    [3.0, -21.576843, 3.2939769e-03], [2.0, -15.402496, 3.3872949e-03],
    [1.0, -13.186676, 3.7360021e-03],
])
SST_FRANCINE, MSL_FRANCINE = 28.20263671875, 1014.9654541015625
# fmt: on


def test_francine_issue77_profile_converges():
    # The exact sounding of GitHub issue #77 / PR #77: Hurricane Francine's
    # pre-landfall environment in a -1.6 K counterfactual-SST sensitivity
    # member. Under the legacy CAPE rule the iteration locked into a bit-exact
    # period-2 cycle (950.653 <-> 951.279 hPa) and returned IFL=2.
    vmax, pmin, ifl, to, otl = _run(PROFILE_FRANCINE, SST_FRANCINE, MSL_FRANCINE)
    assert ifl == 1
    np.testing.assert_allclose(vmax, 83.8498407491696, rtol=1e-12)
    np.testing.assert_allclose(pmin, 920.8856926042054, rtol=1e-12)
    np.testing.assert_allclose(to, 206.09980227350385, rtol=1e-12)
    np.testing.assert_allclose(otl, 115.99623060103113, rtol=1e-12)


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
        elif ifl == 2:
            # non-convergence always returns missing values
            assert np.isnan(vmax)
        # ifl in (0, 3): a flag tripped by the environmental-CAPE call persists
        # (issue #78) and may accompany finite outputs; no NaN requirement.
