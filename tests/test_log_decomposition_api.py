import numpy as np

import tcpyPI
from tcpyPI import utilities


def test_pi_log_decomposition_matches_utilities_scalar():
    out = tcpyPI.pi_log_decomposition(70, 300, 200, CKCD=0.9, sst_units="K")
    expected = utilities.decompose_pi(70, 300, 200, CKCD=0.9)
    np.testing.assert_allclose(out, expected, rtol=0, atol=1e-12)


def test_pi_log_decomposition_sst_units_celsius():
    out_k = tcpyPI.pi_log_decomposition(70, 300, 200, CKCD=0.9, sst_units="K")
    out_c = tcpyPI.pi_log_decomposition(70, 26.85, 200, CKCD=0.9, sst_units="C")
    np.testing.assert_allclose(out_k, out_c, rtol=0, atol=1e-12)


def test_pi_log_decomposition_vectorized_inputs():
    pi = np.array([70.0, 0.0])
    lnpi, lneff, lndiseq, lnCKCD = tcpyPI.pi_log_decomposition(
        pi, 300.0, 200.0, CKCD=0.9, sst_units="K"
    )

    assert np.isfinite(lnpi[0])
    assert np.isfinite(lneff[0])
    assert np.isfinite(lndiseq[0])

    assert np.isnan(lnpi[1])       # pi <= 0
    assert np.isfinite(lneff[1])
    assert np.isnan(lndiseq[1])

    assert np.isfinite(lnCKCD)


def test_pi_log_decomposition_from_pi_outputs_closes():
    # Decompose the outputs of pi() directly (SST in Celsius); the terms must
    # reconstruct 2*ln(vmax). This replaces the former run-and-decompose wrapper.
    SSTC = 30
    MSL = 1010
    level_data = np.array([
        [1000, 28, 18], [975, 25, 18], [950, 24, 16], [925, 23, 13], [900, 22, 12],
        [875, 20, 11], [850, 19, 10], [825, 18, 10], [800, 16, 9], [775, 15, 8],
        [750, 13, 7], [700, 11, 4], [650, 8, 3], [600, 5, 1.7], [550, 2, 1.2],
        [500, -2, 1.7], [450, -6, 0.7], [400, -11, 0.2], [350, -18, 0.15],
        [300, -27, 0.10], [250, -37, 0.11], [225, -43, 0.08], [200, -49, 0.05],
        [175, -57, 0.03], [150, -65, 0.014], [125, -73, 0.005], [100, -79, 0.003],
        [70, -73, 0.002], [50, -64, 0.002],
    ])
    P, TC, R = level_data[:, 0], level_data[:, 1], level_data[:, 2]
    vmax, pmin, ifl, t0, otl = tcpyPI.pi(SSTC, MSL, P, TC, R, CKCD=0.9)

    lnpi, lneff, lndiseq, lnCKCD = tcpyPI.pi_log_decomposition(
        vmax, SSTC, t0, CKCD=0.9, sst_units="C"
    )
    # decomposition closes: lneff + lndiseq + lnCKCD == lnpi == 2*ln(vmax)
    np.testing.assert_allclose(lneff + lndiseq + lnCKCD, lnpi, rtol=0, atol=1e-9)
    np.testing.assert_allclose(lnpi, 2.0 * np.log(vmax), rtol=0, atol=1e-9)
