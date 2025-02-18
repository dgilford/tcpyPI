import os
from importlib import util

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


def test_run_sample_dataset():
    """Test run_sample_dataset matches expected output."""
    ds_out = run_sample_dataset(DATA_FILE, dim="p", CKCD=0.9)

    # Load expected results and compare
    expected_ds = xr.open_dataset(EXPECTED_RAW_OUTPUT)
    xr.testing.assert_allclose(ds_out, expected_ds, rtol=1e-13, atol=1e-9)


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
        xr.testing.assert_allclose(
            ds_out[var], expected_ds[var], rtol=1e-13, atol=1e-9
        )
