"""Tests for the MCP server tool functions (tcpyPI.mcp_server).

The tool bodies are plain functions, so most tests need neither fastmcp nor a
running server; the grid tool needs xarray and the registration smoke test
needs fastmcp (both guarded with importorskip so basic-tests can collect this
file in minimal environments).
"""

import math
import os

import numpy as np
import pytest

import tcpyPI
from tcpyPI import mcp_server

# Compact ocean sounding (from the pi() docstring example, thinned) — same
# fixture as tests/test_pi_field.py.
P = [1000.0, 975, 950, 925, 900, 850, 800, 750, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50]
TC = [28.0, 25, 24, 23, 22, 19, 16, 13, 11, 5, -2, -11, -27, -37, -49, -65, -79, -73, -64]
R = [18.0, 18, 16, 13, 12, 10, 9, 7, 4, 1.7, 1.7, 0.2, 0.1, 0.05, 0.014, 0.003, 0.002, 0.002, 0.002]

DATA_FILE = os.path.join("data", "sample_data.nc")
EXPECTED_RAW_OUTPUT = os.path.join("data", "raw_sample_output.nc")


def test_compute_pi_matches_direct_call():
    out = mcp_server.compute_pi(28.5, 1005.0, P, TC, R)
    assert out["status"] == "success"

    vmax, pmin, ifl, t0, otl = tcpyPI.pi(28.5, 1005.0, np.array(P), np.array(TC), np.array(R))
    assert out["vmax_ms"] == pytest.approx(float(vmax), rel=0, abs=0)
    assert out["pmin_hpa"] == pytest.approx(float(pmin), rel=0, abs=0)
    assert out["ifl"] == int(ifl)
    assert out["t0_k"] == pytest.approx(float(t0), rel=0, abs=0)
    assert out["otl_hpa"] == pytest.approx(float(otl), rel=0, abs=0)


def test_compute_pi_contract_fields():
    out = mcp_server.compute_pi(28.5, 1005.0, P, TC, R)
    assert out["trustworthy"] is (out["ifl"] == 1)
    assert "converged" in out["ifl_meaning"]

    prov = out["provenance"]
    assert prov["doi"] == "10.5281/zenodo.3756005"
    assert prov["tcpypi_version"]
    assert isinstance(prov["numba"], bool)
    assert prov["settings"]["CKCD"] == 0.9
    assert prov["settings"]["V_reduc"] == 0.8
    assert prov["settings"]["outflow_source"] == "cape_star"
    assert len(prov["citations"]) == 2


def test_compute_pi_nan_sst_flags_untrustworthy_and_is_json_safe():
    out = mcp_server.compute_pi(float("nan"), 1005.0, P, TC, R)
    assert out["status"] == "success"
    assert out["ifl"] != 1
    assert out["trustworthy"] is False
    # NaN outputs must be serialized as None, never NaN
    for key in ("vmax_ms", "pmin_hpa", "t0_k", "otl_hpa"):
        assert out[key] is None or math.isfinite(out[key])


def test_compute_pi_accepts_json_null_for_missing_values():
    """JSON cannot carry NaN: null (None) must map to NaN, not be rejected."""
    out = mcp_server.compute_pi(None, 1005.0, P, TC, R)
    assert out["status"] == "success"
    assert out["trustworthy"] is False

    r_with_gap = list(R)
    r_with_gap[3] = None
    out2 = mcp_server.compute_pi(28.5, 1005.0, P, TC, r_with_gap)
    assert out2["status"] == "success"
    assert out2["ifl"] in _valid_ifl_values()


def _valid_ifl_values():
    return set(mcp_server._IFL_MEANING)


def test_compute_pi_length_mismatch_is_structured_error():
    out = mcp_server.compute_pi(28.5, 1005.0, P, TC[:-1], R)
    assert out["status"] == "error"
    assert "equal lengths" in out["message"]


def test_compute_pi_oversize_profile_redirects_to_grid_tool():
    n = mcp_server._MAX_INLINE_LEVELS + 1
    out = mcp_server.compute_pi(28.5, 1005.0, [1000.0] * n, [20.0] * n, [10.0] * n)
    assert out["status"] == "error"
    assert "compute_pi_grid" in out["message"]


def test_decompose_pi_terms_sum_in_log_space():
    pi_out = mcp_server.compute_pi(28.5, 1005.0, P, TC, R)
    out = mcp_server.decompose_pi(pi_out["vmax_ms"], 28.5, pi_out["t0_k"])
    assert out["status"] == "success"
    assert out["lnpi"] == pytest.approx(out["lneff"] + out["lndiseq"] + out["lnCKCD"], rel=1e-12)
    assert out["lnpi"] == pytest.approx(2.0 * math.log(pi_out["vmax_ms"]), rel=1e-12)


def test_ventilation_index_matches_direct_call():
    out = mcp_server.ventilation_index(10.0, 70.0, 0.8)
    assert out["status"] == "success"
    expected = float(tcpyPI.ventilation_index(10.0, 70.0, 0.8))
    assert out["ventilation_index"] == pytest.approx(expected, rel=0, abs=0)


def test_ventilation_index_accepts_lists():
    out = mcp_server.ventilation_index([10.0, 20.0], [70.0, 70.0], [0.8, 0.8])
    assert out["status"] == "success"
    assert len(out["ventilation_index"]) == 2
    assert out["ventilation_index"][1] == pytest.approx(2 * out["ventilation_index"][0])


def test_gpi_en04_matches_direct_call():
    out = mcp_server.genesis_potential_index(5e-5, rh_mid_pct=70.0, v_shear_ms=10.0, v_pot_ms=70.0)
    assert out["status"] == "success"
    expected = float(tcpyPI.genesis_potential_index(5e-5, 70.0, 10.0, 70.0))
    assert out["gpi"] == pytest.approx(expected, rel=0, abs=0)


def test_gpi_e10_requires_chi():
    out = mcp_server.genesis_potential_index(
        5e-5, v_shear_ms=10.0, v_pot_ms=70.0, formulation="e10"
    )
    assert out["status"] == "error"
    assert out["parameter"] == "chi"


def test_gpi_en04_requires_rh():
    out = mcp_server.genesis_potential_index(5e-5, v_shear_ms=10.0, v_pot_ms=70.0)
    assert out["status"] == "error"
    assert out["parameter"] == "rh_mid_pct"


# ---------------------------------------------------------------------------
# Grid tool (needs xarray)
# ---------------------------------------------------------------------------


def test_compute_pi_grid_round_trip(tmp_path):
    xr = pytest.importorskip("xarray")
    out_file = str(tmp_path / "pi_out.nc")

    out = mcp_server.compute_pi_grid(DATA_FILE, out_file)
    assert out["status"] == "success", out
    assert out["output_file"] == out_file

    # written PI fields must match the tracked regression reference
    got = xr.open_dataset(out_file)
    expected = xr.open_dataset(EXPECTED_RAW_OUTPUT)
    for var in ("vmax", "pmin", "ifl", "t0", "otl"):
        xr.testing.assert_allclose(got[var], expected[var], rtol=1e-13, atol=1e-9)

    # summary contract
    assert out["n_columns"] == int(expected["ifl"].size)
    assert sum(out["ifl_histogram"].values()) == out["n_columns"]
    assert out["n_trustworthy"] == int((np.asarray(expected["ifl"]) == 1).sum())
    stats = out["vmax_ms_over_trustworthy"]
    assert stats["min"] <= stats["mean"] <= stats["max"]
    assert out["provenance"]["settings"]["CKCD"] == 0.9


def test_compute_pi_grid_rejects_wrong_units(tmp_path):
    xr = pytest.importorskip("xarray")
    ds = xr.open_dataset(DATA_FILE)
    ds["sst"].attrs["units"] = "K"  # wrong on purpose; values unchanged
    bad_file = str(tmp_path / "bad_units.nc")
    ds.to_netcdf(bad_file)

    out = mcp_server.compute_pi_grid(bad_file, str(tmp_path / "out.nc"))
    assert out["status"] == "error"
    assert out["variable"] == "sst"
    assert "never converts" in out["message"]


def test_compute_pi_grid_rejects_missing_units(tmp_path):
    xr = pytest.importorskip("xarray")
    ds = xr.open_dataset(DATA_FILE)
    del ds["r"].attrs["units"]
    bad_file = str(tmp_path / "no_units.nc")
    ds.to_netcdf(bad_file)

    out = mcp_server.compute_pi_grid(bad_file, str(tmp_path / "out.nc"))
    assert out["status"] == "error"
    assert out["variable"] == "r"
    assert "no `units` attribute" in out["message"]


def test_compute_pi_grid_missing_input_is_structured_error(tmp_path):
    pytest.importorskip("xarray")
    out = mcp_server.compute_pi_grid(str(tmp_path / "nope.nc"), str(tmp_path / "out.nc"))
    assert out["status"] == "error"
    assert out["parameter"] == "input_nc"


def test_compute_pi_grid_refuses_overwriting_input(tmp_path):
    pytest.importorskip("xarray")
    out = mcp_server.compute_pi_grid(DATA_FILE, DATA_FILE)
    assert out["status"] == "error"
    assert out["parameter"] == "output_nc"


# ---------------------------------------------------------------------------
# Server registration (needs fastmcp)
# ---------------------------------------------------------------------------


def test_create_server_registers_all_five_tools():
    pytest.importorskip("fastmcp")
    import asyncio

    server = mcp_server.create_server()
    tools = asyncio.run(server.list_tools())
    assert {t.name for t in tools} == {
        "compute_pi",
        "compute_pi_grid",
        "decompose_pi",
        "ventilation_index",
        "genesis_potential_index",
    }
