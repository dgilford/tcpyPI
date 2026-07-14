import os
from importlib import util

import numpy as np
import xarray as xr

# Import run_sample.py directly from the file
spec = util.spec_from_file_location("run_sample", "run_sample.py")
if spec is None:
    raise ImportError("Failed to import run_sample.py")
run_sample = util.module_from_spec(spec)
if spec.loader is None:
    raise ImportError("Failed to get loader for run_sample.py")
spec.loader.exec_module(run_sample)

# Now we can use run_sample's functions
run_sample_dataset = run_sample.run_sample_dataset
run_sample_analyses = run_sample.run_sample_analyses

# Define paths to input and reference output files
DATA_FILE = os.path.join("data", "sample_data.nc")
EXPECTED_RAW_OUTPUT = os.path.join("data", "raw_sample_output.nc")
EXPECTED_FULL_OUTPUT = os.path.join("data", "full_sample_output.nc")

# The pinned .nc outputs predate the max-running-CAPE LNB selection (issue
# #77).  Two kinds of differences are expected until the pins are
# regenerated:
#   * iteration-path noise: the pressure loop stops at a slightly different
#     iterate when CAPE changes in the last bits (different summation
#     order), moving outputs by up to ~0.1 m/s / ~0.1 hPa;
#   * intentional changes: columns whose energetically optimal outflow
#     level differs from the topmost positive-buoyancy level (measured:
#     ~1-2.5% of columns on this sample, depending on the variable).
# So compare in three parts: the NaN/failure pattern must match exactly,
# differences beyond the per-variable noise scale must stay a small
# fraction, and the intentional changes must stay within physical bounds.
#
# Per variable: (noise_atol, max_changed_fraction, tail_bound, nan_frac).
# nan_frac > 0 only for log-derived analysis variables: columns whose vmax
# moves between exactly 0 and a small positive value legitimately flip
# their logs between non-finite and finite (measured 0.05% here).  The raw
# outputs' NaN patterns (the failure flags) must remain identical.
THRESHOLDS = {
    # raw PI outputs
    "vmax": (0.1, 0.03, 20.0, 0.0),  # m/s
    "pmin": (0.1, 0.03, 5.0, 0.0),  # hPa
    "ifl": (0.0, 0.0, 0.0, 0.0),  # flags must be identical
    "t0": (0.1, 0.03, 80.0, 0.0),  # K
    "otl": (1.0, 0.03, 800.0, 0.0),  # hPa
    # analysis variables
    "eff": (1e-3, 0.03, np.inf, 0.005),
    "diseq": (10.0, 0.03, np.inf, 0.005),  # J/kg; diverges as eff -> 0
    "lnpi": (1e-3, 0.03, np.inf, 0.005),
    "lneff": (1e-3, 0.03, np.inf, 0.005),
    "lndiseq": (1e-2, 0.03, np.inf, 0.005),
    "lnCKCD": (0.0, 0.0, 0.0, 0.0),  # a constant
}


def _assert_matches_pin(actual, expected, name):
    noise_atol, max_frac, tail_bound, nan_frac = THRESHOLDS[name]
    a = np.asarray(actual, dtype=float)
    b = np.asarray(expected, dtype=float)
    # failure/NaN pattern must match (exactly, unless log-derived)
    pattern_mismatch = np.isfinite(a) != np.isfinite(b)
    assert pattern_mismatch.mean() <= nan_frac, (
        f"{name}: NaN pattern changed for {100 * pattern_mismatch.mean():.3f}% "
        f"of points (allowed {100 * nan_frac:.1f}%)"
    )
    finite = np.isfinite(a) & np.isfinite(b)
    d = np.abs(a - b)[finite]
    changed = d > noise_atol
    frac = changed.mean() if d.size else 0.0
    assert frac <= max_frac, (
        f"{name}: {100 * frac:.2f}% of points differ from the pinned output "
        f"beyond the {noise_atol:g} noise scale "
        f"(allowed {100 * max_frac:.1f}% for the LNB-selection change)"
    )
    if changed.any():
        assert d[changed].max() <= tail_bound, (
            f"{name}: max change {d[changed].max():.4g} exceeds bound {tail_bound}"
        )


def test_run_sample_dataset():
    """Test run_sample_dataset matches expected output."""
    ds_out = run_sample_dataset(DATA_FILE, dim="p", CKCD=0.9)

    # Load expected results and compare
    expected_ds = xr.open_dataset(EXPECTED_RAW_OUTPUT)
    for var in ["vmax", "pmin", "ifl", "t0", "otl"]:
        _assert_matches_pin(ds_out[var].values, expected_ds[var].values, var)


def test_run_sample_analyses():
    """Test run_sample_analyses matches expected output."""
    # Get input dataset from first step
    ds = run_sample_dataset(DATA_FILE, dim="p", CKCD=0.9)

    # Run analyses
    ds_out = run_sample_analyses(ds, CKCD=0.9)

    # Load expected results and compare
    expected_ds = xr.open_dataset(EXPECTED_FULL_OUTPUT)
    # Only compare the analysis variables (not the input variables)
    analysis_vars = ["eff", "diseq", "lnpi", "lneff", "lndiseq", "lnCKCD"]
    for var in analysis_vars:
        _assert_matches_pin(ds_out[var].values, expected_ds[var].values, var)
