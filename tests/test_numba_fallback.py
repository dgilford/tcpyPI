"""tcpyPI must remain importable and correct when Numba is unavailable.

Runs a fresh interpreter with the `numba` import blocked and checks that
importing tcpyPI warns (RuntimeWarning), and that pi() and pi_field() still
produce correct (pure-Python) results.
"""

import subprocess
import sys
from pathlib import Path

_SRC = str(Path(__file__).resolve().parents[1] / "src")

_SCRIPT = r"""
import sys

class _BlockNumba:
    def find_spec(self, name, path=None, target=None):
        if name == "numba" or name.startswith("numba."):
            raise ImportError("numba import blocked for fallback test")
        return None

sys.meta_path.insert(0, _BlockNumba())
sys.path.insert(0, {src!r})

import warnings
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    import tcpyPI
assert any(
    issubclass(w.category, RuntimeWarning) and "pure-Python" in str(w.message)
    for w in caught
), "expected a pure-Python-mode RuntimeWarning on import without numba"

import numpy as np
P = np.array([1000.0, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50])
TC = np.array([28.0, 22, 16, 11, 5, -2, -11, -27, -49, -79, -64])
R = np.array([18.0, 12, 9, 4, 1.7, 1.7, 0.2, 0.1, 0.05, 0.003, 0.002])

vmax, pmin, ifl, to, otl = tcpyPI.pi(30.0, 1010.0, P, TC, R)
assert ifl == 1 and np.isfinite(vmax), (vmax, ifl)

# pi_field falls back to a per-column loop and must match pi() exactly.
fv, fp, fi, ft, fo = tcpyPI.pi_field(np.array([[30.0, np.nan]]), 1010.0, P, TC, R)
assert fi[0, 0] == 1 and fi[0, 1] == 0
assert fv[0, 0] == vmax
print("fallback OK")
"""


def test_import_and_compute_without_numba():
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT.format(src=_SRC)],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    assert "fallback OK" in proc.stdout
