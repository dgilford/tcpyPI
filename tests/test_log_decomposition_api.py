import numpy as np

import tcpyPI
from tcpyPI import utilities


def test_log_decompose_pi_matches_utilities_scalar():
    out = tcpyPI.log_decompose_pi(70, 300, 200, CKCD=0.9, sst_units="K")
    expected = utilities.decompose_pi(70, 300, 200, CKCD=0.9)
    np.testing.assert_allclose(out, expected, rtol=0, atol=1e-12)


def test_log_decompose_pi_sst_units_celsius():
    out_k = tcpyPI.log_decompose_pi(70, 300, 200, CKCD=0.9, sst_units="K")
    out_c = tcpyPI.log_decompose_pi(70, 26.85, 200, CKCD=0.9, sst_units="C")
    np.testing.assert_allclose(out_k, out_c, rtol=0, atol=1e-12)


def test_log_decompose_pi_vectorized_inputs():
    pi = np.array([70.0, 0.0])
    lnpi, lneff, lndiseq, lnCKCD = tcpyPI.log_decompose_pi(
        pi, 300.0, 200.0, CKCD=0.9, sst_units="K"
    )

    assert np.isfinite(lnpi[0])
    assert np.isfinite(lneff[0])
    assert np.isfinite(lndiseq[0])

    assert np.isnan(lnpi[1])
    assert np.isfinite(lneff[1])
    assert np.isnan(lndiseq[1])

    assert np.isfinite(lnCKCD)


def test_pi_log_decomposition_smoke_and_consistency():
    # Example from `tcpyPI.pi` docstring
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

    out = tcpyPI.pi_log_decomposition(SSTC, MSL, P, TC, R, CKCD=0.9)

    for key in ["vmax", "pmin", "ifl", "t0", "otl", "lnpi", "lneff", "lndiseq", "lnCKCD"]:
        assert key in out

    vmax, pmin, ifl, t0, otl = tcpyPI.pi(SSTC, MSL, P, TC, R, CKCD=0.9)
    assert out["vmax"] == vmax
    assert out["pmin"] == pmin
    assert out["ifl"] == ifl
    assert out["t0"] == t0
    assert out["otl"] == otl

    expected_lnpi, expected_lneff, expected_lndiseq, expected_lnCKCD = tcpyPI.log_decompose_pi(
        vmax, SSTC, t0, CKCD=0.9, sst_units="C"
    )
    np.testing.assert_allclose(
        [out["lnpi"], out["lneff"], out["lndiseq"], out["lnCKCD"]],
        [expected_lnpi, expected_lneff, expected_lndiseq, expected_lnCKCD],
        rtol=0,
        atol=1e-12,
    )
