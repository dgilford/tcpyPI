import numpy as np

import tcpyPI


def test_gpi_en04_scalar_matches_hand_calc():
    abs_vort = 2.0e-5  # s^-1
    rh_mid = 60.0      # %
    v_shear = 10.0     # m/s
    v_pot = 60.0       # m/s

    eta = 1.0e5 * abs_vort
    rh = rh_mid / 50.0
    vp = v_pot / 70.0
    shear_term = (1.0 + 0.1 * v_shear) ** (-2.0)
    expected = (eta**3.0) * (rh**3.0) * (vp**3.0) * shear_term

    out = tcpyPI.genesis_potential_index(abs_vort, rh_mid, v_shear, v_pot, formulation="en04")
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_gpi_c07_pi_thresholding_reduces_value_when_vp_low():
    abs_vort = 2.0e-5
    rh_mid = 60.0
    v_shear = 10.0
    v_pot = 40.0  # close to the 35 m/s threshold in the c07 variant

    en04 = tcpyPI.genesis_potential_index(abs_vort, rh_mid, v_shear, v_pot, formulation="en04")
    c07 = tcpyPI.genesis_potential_index(abs_vort, rh_mid, v_shear, v_pot, formulation="c07")
    assert c07 < en04


def test_gpi_broadcasting_and_invalid_inputs():
    out = tcpyPI.genesis_potential_index(
        abs_vort=np.array([2.0e-5, 2.0e-5]),
        rh_mid=60.0,
        v_shear=np.array([5.0, 10.0]),
        v_pot=60.0,
        formulation="en04",
    )
    assert out.shape == (2,)

    out2 = tcpyPI.genesis_potential_index(
        abs_vort=np.array([-1.0e-5, 2.0e-5]),
        rh_mid=60.0,
        v_shear=10.0,
        v_pot=60.0,
        formulation="en04",
    )
    assert np.isnan(out2[0])
    assert np.isfinite(out2[1])

    with np.testing.assert_raises(ValueError):
        tcpyPI.genesis_potential_index(2.0e-5, 60.0, 10.0, 60.0, formulation="nope")
