import os
from importlib import util

import numpy as np
import pytest
import xarray as xr

from tcpyPI.pi import pi

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
        xr.testing.assert_allclose(ds_out[var], expected_ds[var], rtol=1e-13, atol=1e-9)


@pytest.mark.xfail(reason="TODO: fix this test")
def test_for_lack_of_pi_convergence():
    """Test case where pi does not converge.

    The value of PM settles to alternating between 950.6533454438319 and
    951.2790079838533, until the iteration limit is reached.
    """
    params = {
        "SSTC": 28.20263671875,
        "MSL": 1014.9654541015625,
        "P": np.array(
            [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750]
            + [700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 225]
            + [200, 175, 150, 125, 100, 70, 50, 30, 20, 10, 7]
            + [5, 3, 2, 1]
        ),
        "TC": np.array(
            [25.260956, 23.078949, 20.881561, 18.652649, 16.615143]
            + [14.63147, 12.794586, 11.79306, 11.01236, 10.847565]
            + [10.354492, 7.807007, 5.473297, 2.5278625, -1.532135]
            + [-7.3142395, -13.635345, -20.613602, -28.588928, -37.270096]
            + [-45.40825, -50.047455, -54.84575, -59.173737, -62.878662]
            + [-65.78009, -69.30669, -65.512024, -60.76361, -55.100464]
            + [-51.38333, -43.95462, -39.709717, -33.193268, -21.576843]
            + [-15.402496, -13.186676]
        ),
        "R": np.array(
            [1.0783693e01, 1.0704287e01, 1.0680210e01, 1.0617845e01]
            + [1.0320683e01, 9.8811483e00, 9.1884289e00, 7.1884680e00]
            + [5.6963191e00, 3.5568204e00, 1.5912720e00, 1.0433695e00]
            + [5.9723043e-01, 4.3974420e-01, 4.8722979e-01, 5.8590513e-01]
            + [4.5599860e-01, 3.1293562e-01, 1.9222400e-01, 9.7611703e-02]
            + [3.3851895e-02, 2.4188591e-02, 1.9636340e-02, 1.3214327e-02]
            + [7.2453087e-03, 4.3027173e-03, 3.7014042e-03, 2.9414182e-03]
            + [2.7806845e-03, 2.8306348e-03, 2.9053832e-03, 2.9848625e-03]
            + [3.0779461e-03, 3.1315640e-03, 3.2939769e-03, 3.3872949e-03]
            + [3.7360021e-03]
        ),
    }
    VMAX, PMIN, IFL, TO, OTL = pi(**params)
    assert IFL == 1
