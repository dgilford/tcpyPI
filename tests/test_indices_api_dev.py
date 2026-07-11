import numpy as np

import tcpyPI
import tcpyPI.utilities as util


def _pi_example_profile():
    # Same profile as the example in `tcpyPI.pi` docstring / other tests.
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


def test_vi_te12_scalar_matches_hand_calc():
    v_shear = 10.0
    v_pot = 50.0
    chi = 0.5
    expected = (v_shear / v_pot) * chi
    out = tcpyPI.ventilation_index(v_shear, v_pot, chi, formulation="te12")
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_chi_m_te12_matches_hand_calc():
    s_m_star = 350.0
    s_m = 340.0
    s_sst_star = 360.0
    s_b = 345.0
    expected = (s_m_star - s_m) / (s_sst_star - s_b)
    out = tcpyPI.entropy_deficit_te12(s_m_star, s_m, s_sst_star, s_b)
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_vi_te12_can_compute_chi_m_from_entropies():
    v_shear = 10.0
    v_pot = 50.0
    s_m_star = 350.0
    s_m = 340.0
    s_sst_star = 360.0
    s_b = 345.0
    chi_m = (s_m_star - s_m) / (s_sst_star - s_b)
    expected = (v_shear / v_pot) * chi_m
    out = tcpyPI.ventilation_index(
        v_shear,
        v_pot,
        formulation="te12",
        s_m_star=s_m_star,
        s_m=s_m,
        s_sst_star=s_sst_star,
        s_b=s_b,
    )
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_chi_m_te12_from_profile_matches_entropy_algebra():
    SSTC, MSL, P, TC, R = _pi_example_profile()

    # Midlevel is 600 hPa in this example profile.
    mid_i = int(np.where(P == 600)[0][0])
    p_mid = 600.0

    T_mid_C = float(TC[mid_i])
    T_mid_K = T_mid_C + 273.15
    R_mid = float(R[mid_i]) * 0.001  # g/kg -> g/g

    # sm: environmental entropy at midlevels
    s_m = util.entropy_S(T_mid_K, R_mid, p_mid)

    # s*m: saturation entropy at midlevels (using environmental temperature, RH=100%)
    es_mid = util.es_cc(T_mid_C)
    R_mid_sat = util.rv(es_mid, p_mid)
    s_m_star = util.entropy_S(T_mid_K, R_mid_sat, p_mid)

    # s*SST: saturation entropy at SST at surface pressure
    SST_K = SSTC + 273.15
    es_sst = util.es_cc(SSTC)
    R_sst_sat = util.rv(es_sst, float(MSL))
    s_sst_star = util.entropy_S(SST_K, R_sst_sat, float(MSL))

    # sb: boundary-layer entropy; default uses 2m (here proxied by lowest level T/q)
    T2m_C = float(TC[0])
    T2m_K = T2m_C + 273.15
    R2m = float(R[0]) * 0.001
    s_b = util.entropy_S(T2m_K, R2m, float(MSL))

    expected = (s_m_star - s_m) / (s_sst_star - s_b)

    out = tcpyPI.entropy_deficit_te12_from_profile(
        P,
        TC,
        R,
        T_units="C",
        q_units="g/kg",
        SST=SSTC,
        SST_units="C",
        psfc_hpa=MSL,
        T2m=T2m_C,
        q2m=float(R[0]),
        entropy_method="emanuel94",
    )

    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_chi_m_te12_from_profile_bryan2008_runs():
    SSTC, MSL, P, TC, R = _pi_example_profile()
    out = tcpyPI.entropy_deficit_te12_from_profile(
        P,
        TC,
        R,
        T_units="C",
        q_units="g/kg",
        SST=SSTC,
        SST_units="C",
        psfc_hpa=MSL,
        T2m=float(TC[0]),
        q2m=float(R[0]),
        entropy_method="bryan2008",
    )
    assert np.isfinite(out)


def test_chi_m_te12_from_profile_moist_adiabat_from_sst_runs_emanuel94():
    SSTC, MSL, P, TC, R = _pi_example_profile()
    out = tcpyPI.entropy_deficit_te12_from_profile(
        P,
        TC,
        R,
        T_units="C",
        q_units="g/kg",
        SST=SSTC,
        SST_units="C",
        psfc_hpa=MSL,
        T2m=float(TC[0]),
        q2m=float(R[0]),
        entropy_method="emanuel94",
        s_m_star_source="moist_adiabat_from_sst",
    )
    assert np.isfinite(out)


def test_chi_m_te12_from_profile_moist_adiabat_from_sst_runs_bryan2008():
    SSTC, MSL, P, TC, R = _pi_example_profile()
    out = tcpyPI.entropy_deficit_te12_from_profile(
        P,
        TC,
        R,
        T_units="C",
        q_units="g/kg",
        SST=SSTC,
        SST_units="C",
        psfc_hpa=MSL,
        T2m=float(TC[0]),
        q2m=float(R[0]),
        entropy_method="bryan2008",
        s_m_star_source="moist_adiabat_from_sst",
    )
    assert np.isfinite(out)


def test_gpi_en04_scalar_matches_hand_calc():
    # Example values in typical ranges.
    abs_vort = 2.0e-5  # s^-1
    rh_mid = 60.0      # %
    v_shear = 10.0     # m/s
    v_pot = 60.0       # m/s

    eta = 1.0e5 * abs_vort
    rh = rh_mid / 50.0
    vp = v_pot / 70.0
    shear_term = (1.0 + 0.1 * v_shear) ** (-2.0)
    expected = (eta**1.5) * (rh**3.0) * (vp**3.0) * shear_term

    out = tcpyPI.genesis_potential_index(
        abs_vort, rh_mid, v_shear, v_pot, formulation="en04"
    )
    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_gpi_e10_scalar_matches_hand_calc():
    # Emanuel (2010, JAMES, Eq. 11): |eta|^3 chi^-4/3 max(Vp-35,0)^2 (25+Vshear)^-4
    abs_vort = 5.0e-5
    chi = 1.2
    v_shear = 10.0
    v_pot = 60.0

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


def test_pdi_e05_si_and_scaled():
    vmax = np.array([10.0, 20.0, 30.0])  # m/s
    dt = 3600.0  # seconds
    core = float(np.sum((vmax**3.0) * dt))

    out_si = tcpyPI.power_dissipation_index(vmax, dt, formulation="e05_si", wind_units="m/s", dt_units="s")
    out_scaled = tcpyPI.power_dissipation_index(vmax, dt, formulation="e05_1e11", wind_units="m/s", dt_units="s")

    np.testing.assert_allclose(out_si, core, rtol=0, atol=0)
    np.testing.assert_allclose(out_scaled, core / 1.0e11, rtol=0, atol=0)
