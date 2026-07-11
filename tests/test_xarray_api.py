"""Tests for the packaged whole-field xarray entry point (tcpyPI.xarray.pi_dataset)."""

import os

import pytest

xr = pytest.importorskip("xarray")

from tcpyPI.xarray import pi_dataset  # noqa: E402

DATA_FILE = os.path.join("data", "sample_data.nc")
EXPECTED_RAW_OUTPUT = os.path.join("data", "raw_sample_output.nc")

PI_VARS = ["vmax", "pmin", "ifl", "t0", "otl"]


def test_pi_dataset_matches_reference_pins():
    """pi_dataset reproduces the tracked run_sample regression output exactly."""
    ds = xr.open_dataset(DATA_FILE)
    out = pi_dataset(ds, dim="p", CKCD=0.9)

    expected = xr.open_dataset(EXPECTED_RAW_OUTPUT)
    for var in PI_VARS:
        xr.testing.assert_allclose(out[var], expected[var], rtol=1e-13, atol=1e-9)


def test_pi_dataset_output_metadata():
    """Outputs carry standard_name/units attrs and only the PI variables."""
    ds = xr.open_dataset(DATA_FILE)
    out = pi_dataset(ds)

    assert sorted(out.data_vars) == sorted(PI_VARS)
    assert out.vmax.attrs["units"] == "m/s"
    assert out.pmin.attrs["units"] == "hPa"
    assert out.t0.attrs["units"] == "K"
    assert out.otl.attrs["units"] == "hPa"
    assert out.ifl.attrs["standard_name"] == "pyPI Flag"


def test_pi_dataset_missing_variable_errors():
    """A dataset missing required variables raises a ValueError naming them."""
    ds = xr.open_dataset(DATA_FILE).drop_vars("sst")
    with pytest.raises(ValueError, match="sst"):
        pi_dataset(ds)


def test_pi_dataset_missing_dim_errors():
    """A wrong vertical-coordinate name raises a ValueError naming it."""
    ds = xr.open_dataset(DATA_FILE)
    with pytest.raises(ValueError, match="plev"):
        pi_dataset(ds, dim="plev")
