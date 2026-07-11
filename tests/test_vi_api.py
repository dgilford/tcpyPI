import numpy as np

import tcpyPI
import tcpyPI.utilities as util


def _pi_example_profile():
    # Same profile as the example in `tcpyPI.pi` docstring / other tests.
    sst_c = 30
    msl_hpa = 1010
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
    p_hpa = level_data[:, 0]
    t_c = level_data[:, 1]
    r_gkg = level_data[:, 2]
    return sst_c, msl_hpa, p_hpa, t_c, r_gkg


def test_vi_te12_scalar_matches_hand_calc():
    v_shear = 10.0
    v_pot = 50.0
    chi_m = 0.5
    expected = (v_shear / v_pot) * chi_m
    out = tcpyPI.ventilation_index(v_shear, v_pot, chi_m, formulation="te12")
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


def test_chi_m_te12_from_profile_env_T_at_pmid_matches_entropy_algebra():
    sst_c, msl_hpa, p, tc, r = _pi_example_profile()

    p_mid_hpa = 600.0
    mid_i = int(np.where(p == p_mid_hpa)[0][0])

    t_mid_c = float(tc[mid_i])
    t_mid_k = t_mid_c + 273.15
    r_mid = float(r[mid_i]) * 0.001  # g/kg -> g/g

    s_m = util.entropy_S(t_mid_k, r_mid, p_mid_hpa)

    es_mid = util.es_cc(t_mid_c)
    r_mid_sat = util.rv(es_mid, p_mid_hpa)
    s_m_star = util.entropy_S(t_mid_k, r_mid_sat, p_mid_hpa)

    sst_k = sst_c + 273.15
    es_sst = util.es_cc(sst_c)
    r_sst_sat = util.rv(es_sst, float(msl_hpa))
    s_sst_star = util.entropy_S(sst_k, r_sst_sat, float(msl_hpa))

    t2m_c = float(tc[0])
    t2m_k = t2m_c + 273.15
    r2m = float(r[0]) * 0.001
    s_b = util.entropy_S(t2m_k, r2m, float(msl_hpa))

    expected = (s_m_star - s_m) / (s_sst_star - s_b)

    out = tcpyPI.entropy_deficit_te12_from_profile(
        p,
        tc,
        r,
        T_units="C",
        q_units="g/kg",
        SST=sst_c,
        SST_units="C",
        psfc_hpa=msl_hpa,
        T2m=t2m_c,
        q2m=float(r[0]),
        entropy_method="emanuel94",
        s_m_star_source="env_T_at_pmid",
    )

    np.testing.assert_allclose(out, expected, rtol=0, atol=0)


def test_chi_m_te12_from_profile_moist_adiabat_from_sst_runs():
    sst_c, msl_hpa, p, tc, r = _pi_example_profile()
    out = tcpyPI.entropy_deficit_te12_from_profile(
        p,
        tc,
        r,
        T_units="C",
        q_units="g/kg",
        SST=sst_c,
        SST_units="C",
        psfc_hpa=msl_hpa,
        T2m=float(tc[0]),
        q2m=float(r[0]),
        entropy_method="emanuel94",
        s_m_star_source="moist_adiabat_from_sst",
    )
    assert np.isfinite(out)


def test_vi_log_decomposition_closure_and_handcalc():
    us, vp, ch = 12.0, 75.0, 0.5
    d = tcpyPI.vi_log_decomposition(us, vp, ch)
    assert np.isclose(d["shear"], np.log(us))
    assert np.isclose(d["pi"], -np.log(vp))
    assert np.isclose(d["chi"], np.log(ch))
    # contributions sum to lnvi, and lnvi == ln(ventilation_index)
    assert np.isclose(d["lnvi"], d["shear"] + d["pi"] + d["chi"])
    assert np.isclose(d["lnvi"], np.log(tcpyPI.ventilation_index(us, vp, ch)))
    # non-positive input -> nan
    dn = tcpyPI.vi_log_decomposition(12.0, 75.0, 0.0)
    assert np.isnan(dn["lnvi"]) and np.isnan(dn["chi"])
