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
    expected = (eta**1.5) * (rh**3.0) * (vp**3.0) * shear_term

    out = tcpyPI.genesis_potential_index(abs_vort, rh_mid, v_shear, v_pot, formulation="en04")
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_gpi_e10_scalar_matches_hand_calc():
    # Emanuel (2010, JAMES, Eq. 11): |eta|^3 chi^-4/3 max(Vp-35,0)^2 (25+Vshear)^-4
    abs_vort = 5.0e-5  # s^-1
    chi = 1.2          # nondimensional midlevel entropy deficit
    v_shear = 10.0     # m/s
    v_pot = 60.0       # m/s

    expected = (
        (abs(abs_vort) ** 3.0)
        * (chi ** (-4.0 / 3.0))
        * (max(v_pot - 35.0, 0.0) ** 2.0)
        * ((25.0 + v_shear) ** (-4.0))
    )
    out = tcpyPI.genesis_potential_index(
        abs_vort, v_shear=v_shear, v_pot=v_pot, formulation="e10", chi=chi
    )
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_gpi_e10_threshold_invalid_chi_and_missing_chi():
    # V_pot below the 35 m/s threshold -> exactly zero.
    z = tcpyPI.genesis_potential_index(
        5.0e-5, v_shear=10.0, v_pot=30.0, formulation="e10", chi=1.2
    )
    assert z == 0.0

    # Non-positive chi (entropy deficit) -> NaN (chi^-4/3 undefined).
    assert np.isnan(
        tcpyPI.genesis_potential_index(
            5.0e-5, v_shear=10.0, v_pot=60.0, formulation="e10", chi=0.0
        )
    )

    # e10 without chi raises.
    with np.testing.assert_raises(ValueError):
        tcpyPI.genesis_potential_index(
            5.0e-5, v_shear=10.0, v_pot=60.0, formulation="e10"
        )


def test_gpi_broadcasting_and_invalid_inputs():
    out = tcpyPI.genesis_potential_index(
        abs_vort=np.array([2.0e-5, 2.0e-5]),
        rh_mid=60.0,
        v_shear=np.array([5.0, 10.0]),
        v_pot=60.0,
        formulation="en04",
    )
    assert out.shape == (2,)

    # Signed absolute vorticity uses its magnitude internally (published |η|), so a
    # negative input yields the same finite value as its positive magnitude.
    out2 = tcpyPI.genesis_potential_index(
        abs_vort=np.array([-1.0e-5, 2.0e-5]),
        rh_mid=60.0,
        v_shear=10.0,
        v_pot=60.0,
        formulation="en04",
    )
    pos = tcpyPI.genesis_potential_index(1.0e-5, 60.0, 10.0, 60.0, formulation="en04")
    assert np.isfinite(out2[0])
    assert np.isfinite(out2[1])
    np.testing.assert_allclose(out2[0], pos, rtol=0, atol=0)

    # NaN vorticity still yields NaN.
    assert np.isnan(
        tcpyPI.genesis_potential_index(np.nan, 60.0, 10.0, 60.0, formulation="en04")
    )

    with np.testing.assert_raises(ValueError):
        tcpyPI.genesis_potential_index(2.0e-5, 60.0, 10.0, 60.0, formulation="nope")


def test_gpi_log_decomposition_en04_closure():
    a, rh, vs, vp = 2.0e-5, 60.0, 10.0, 60.0
    d = tcpyPI.gpi_log_decomposition(a, rh_mid=rh, v_shear=vs, v_pot=vp, formulation="en04")
    assert np.isclose(d["lngpi"], d["vorticity"] + d["humidity"] + d["pi"] + d["shear"])
    assert np.isclose(
        d["lngpi"], np.log(tcpyPI.genesis_potential_index(a, rh, vs, vp, formulation="en04"))
    )


def test_gpi_log_decomposition_e10_closure_and_threshold():
    a, ch, vs, vp = 5.0e-5, 1.2, 10.0, 60.0
    d = tcpyPI.gpi_log_decomposition(a, v_shear=vs, v_pot=vp, formulation="e10", chi=ch)
    assert np.isclose(d["lngpi"], d["vorticity"] + d["chi"] + d["pi"] + d["shear"])
    assert np.isclose(
        d["lngpi"],
        np.log(tcpyPI.genesis_potential_index(a, v_shear=vs, v_pot=vp, formulation="e10", chi=ch)),
    )
    # V_PI below the 35 m/s threshold -> undefined (NaN)
    dth = tcpyPI.gpi_log_decomposition(a, v_shear=vs, v_pot=30.0, formulation="e10", chi=ch)
    assert np.isnan(dth["lngpi"])
    # e10 without chi raises
    with np.testing.assert_raises(ValueError):
        tcpyPI.gpi_log_decomposition(a, v_shear=vs, v_pot=vp, formulation="e10")
