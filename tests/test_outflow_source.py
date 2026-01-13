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


def test_pi_outflow_source_invalid_value_raises():
    SSTC, MSL, P, TC, R = _example_profile()
    with pytest.raises(ValueError):
        tcpyPI.pi(SSTC, MSL, P, TC, R, outflow_source="not_a_real_mode")

