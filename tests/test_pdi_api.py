import numpy as np

import tcpyPI


def test_pdi_e05_si_and_scaled():
    vmax = np.array([10.0, 20.0, 30.0])  # m/s
    dt = 3600.0  # seconds
    core = float(np.sum((vmax**3.0) * dt))

    out_si = tcpyPI.power_dissipation_index(
        vmax, dt, formulation="e05_si", wind_units="m/s", dt_units="s"
    )
    out_scaled = tcpyPI.power_dissipation_index(
        vmax, dt, formulation="e05_1e11", wind_units="m/s", dt_units="s"
    )

    np.testing.assert_allclose(out_si, core, rtol=0, atol=0)
    np.testing.assert_allclose(out_scaled, core / 1.0e11, rtol=0, atol=0)


def test_pdi_knots_hours_matches_ms_seconds():
    vmax_kt = np.array([20.0, 30.0, 40.0])
    dt_h = 1.0

    out_kt_h = tcpyPI.power_dissipation_index(vmax_kt, dt_h, wind_units="kt", dt_units="h")

    vmax_ms = vmax_kt * 0.5144444444444445
    out_ms_s = tcpyPI.power_dissipation_index(
        vmax_ms, dt_h * 3600.0, wind_units="m/s", dt_units="s"
    )
    np.testing.assert_allclose(out_kt_h, out_ms_s, rtol=0, atol=0)


def test_pdi_nan_policy_omit_drops_missing_steps():
    vmax = np.array([10.0, np.nan, 30.0])
    dt = np.array([3600.0, 3600.0, 3600.0])

    out_omit = tcpyPI.power_dissipation_index(vmax, dt, nan_policy="omit")
    expected = float(np.nansum((vmax**3.0) * dt))
    np.testing.assert_allclose(out_omit, expected, rtol=0, atol=0)


def test_pdi_nan_policy_propagate_returns_nan():
    vmax = np.array([10.0, np.nan, 30.0])
    dt = np.array([3600.0, 3600.0, 3600.0])
    out = tcpyPI.power_dissipation_index(vmax, dt, nan_policy="propagate")
    assert np.isnan(out)
