import numpy as np
import pytest

import tcpyPI


def _example_profile():
    # Same profile as the example in `tcpyPI.pi` docstring.
    SSTC = 30
    MSL = 1010
    level_data = np.array(
        [
            [1000, 28, 18],
            [975, 25, 18],
            [950, 24, 16],
            [925, 23, 13],
            [900, 22, 12],
            [875, 20, 11],
            [850, 19, 10],
            [825, 18, 10],
            [800, 16, 9],
            [775, 15, 8],
            [750, 13, 7],
            [700, 11, 4],
            [650, 8, 3],
            [600, 5, 1.7],
            [550, 2, 1.2],
            [500, -2, 1.7],
            [450, -6, 0.7],
            [400, -11, 0.2],
            [350, -18, 0.15],
            [300, -27, 0.10],
            [250, -37, 0.11],
            [225, -43, 0.08],
            [200, -49, 0.05],
            [175, -57, 0.03],
            [150, -65, 0.014],
            [125, -73, 0.005],
            [100, -79, 0.003],
            [70, -73, 0.002],
            [50, -64, 0.002],
        ]
    )
    P = level_data[:, 0]
    TC = level_data[:, 1]
    R = level_data[:, 2]
    return SSTC, MSL, P, TC, R


def test_pi_outflow_source_cape_star_and_cape_env_smoke():
    SSTC, MSL, P, TC, R = _example_profile()

    vmax1, pmin1, ifl1, t01, otl1 = tcpyPI.pi(
        SSTC, MSL, P, TC, R, outflow_source="cape_star"
    )
    vmax2, pmin2, ifl2, t02, otl2 = tcpyPI.pi(
        SSTC, MSL, P, TC, R, outflow_source="cape_env"
    )

    assert ifl1 == 1
    assert ifl2 == 1
    assert np.isfinite(vmax1)
    assert np.isfinite(vmax2)
    assert np.isfinite(pmin1)
    assert np.isfinite(pmin2)
    assert np.isfinite(t01)
    assert np.isfinite(t02)
    assert np.isfinite(otl1)
    assert np.isfinite(otl2)


def _no_env_lnb_profile():
    # Column from data/sample_data.nc (MERRA-2 2004, September, 45N 30E), rounded.
    # The environmental parcel has no level of neutral buoyancy at convergence,
    # while the saturated (CAPE*) parcel does (GitHub issue #90).
    SSTC = 20.51
    MSL = 1018.81
    level_data = np.array(
        [
            [1000.0, 17.68, 8.8073], [975.0, 16.44, 7.8708], [950.0, 15.02, 7.4699],
            [925.0, 13.64, 6.985], [900.0, 12.41, 6.455], [875.0, 11.09, 6.0506],
            [850.0, 9.83, 5.5673], [825.0, 8.57, 4.9546], [800.0, 7.39, 4.3265],
            [775.0, 6.06, 3.8513], [750.0, 4.65, 3.4322], [725.0, 3.2, 2.9884],
            [700.0, 1.71, 2.5357], [650.0, -1.51, 1.9899], [600.0, -5.3, 1.5573],
            [550.0, -9.65, 1.1423], [500.0, -14.49, 0.8354], [450.0, -19.92, 0.6208],
            [400.0, -26.43, 0.4179], [350.0, -33.69, 0.2529], [300.0, -41.73, 0.1358],
            [250.0, -49.59, 0.065], [200.0, -54.46, 0.0267], [150.0, -55.76, 0.0072],
            [100.0, -58.47, 0.0031], [70.0, -58.34, 0.0027], [50.0, -57.57, 0.0028],
            [40.0, -56.86, 0.0029], [30.0, -54.41, 0.003], [20.0, -50.38, 0.0031],
            [10.0, -43.0, 0.0034],
        ]
    )
    return SSTC, MSL, level_data[:, 0], level_data[:, 1], level_data[:, 2]


def test_cape_env_without_env_lnb_returns_missing_ifl4():
    SSTC, MSL, P, TC, R = _no_env_lnb_profile()

    vmax, pmin, ifl, t0, otl = tcpyPI.pi(SSTC, MSL, P, TC, R, outflow_source="cape_env")
    assert ifl == 4
    assert np.isnan(vmax) and np.isnan(pmin) and np.isnan(t0) and np.isnan(otl)

    # the same column is well-defined under the default CAPE* outflow
    vmax, pmin, ifl, t0, otl = tcpyPI.pi(SSTC, MSL, P, TC, R, outflow_source="cape_star")
    assert ifl == 1
    assert np.all(np.isfinite([vmax, pmin, t0, otl]))
    assert t0 < 250.0 and otl < 400.0  # an upper-tropospheric outflow


def test_cape_env_without_env_lnb_keeps_pi_when_outflow_not_used():
    SSTC, MSL, P, TC, R = _no_env_lnb_profile()

    # diss_flag=0: RAT=1, so VMAX/PMIN do not depend on TO and match cape_star
    env = tcpyPI.pi(SSTC, MSL, P, TC, R, outflow_source="cape_env", diss_flag=0)
    star = tcpyPI.pi(SSTC, MSL, P, TC, R, outflow_source="cape_star", diss_flag=0)
    assert env[2] == 1
    assert env[:2] == star[:2]
    assert np.isnan(env[3]) and np.isnan(env[4])

    # CAPE* = 0 (cold SST): VMAX = 0 and PMIN = MSL regardless of the outflow
    vmax, pmin, ifl, t0, otl = tcpyPI.pi(14.0, MSL, P, TC, R, outflow_source="cape_env")
    assert ifl == 1
    assert vmax == 0.0
    assert pmin == pytest.approx(MSL)
    assert np.isnan(t0) and np.isnan(otl)


def test_cape_star_without_lnb_reports_missing_outflow():
    # SST too cold for the saturated parcel to be buoyant: CAPE* = 0
    SSTC, MSL, P, TC, R = _example_profile()
    vmax, pmin, ifl, t0, otl = tcpyPI.pi(15.0, MSL, P, TC, R, outflow_source="cape_star")
    assert ifl == 1
    assert vmax == 0.0
    assert np.isfinite(pmin)
    # TO=T[0] and OTL=0 are sentinels from cape(), not an outflow
    assert np.isnan(t0) and np.isnan(otl)


def test_pi_outflow_source_invalid_value_raises():
    SSTC, MSL, P, TC, R = _example_profile()
    with pytest.raises(ValueError):
        tcpyPI.pi(SSTC, MSL, P, TC, R, outflow_source="not_a_real_mode")

