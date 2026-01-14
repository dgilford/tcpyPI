"""Relative intensity utilities.

Relative intensity is commonly defined as the ratio of a storm's realized
intensity to its local potential intensity (PI), often denoted by Greek nu:

    ν = V / V_PI

This is a lightweight convenience function to keep unit-agnostic ratio logic
and NaN handling consistent across analyses.
"""

from typing import Any, Optional

import numpy as np


def relative_intensity(
    v: Any,
    v_pot: Any,
    *,
    nan_policy: str = "propagate",
    clip: Optional[tuple[float, float]] = None,
):
    """Compute relative intensity ν = V / V_PI.

    Parameters
    ----------
    v : float or array-like
        Realized intensity time series (e.g., best-track Vmax).
    v_pot : float or array-like
        Potential intensity time series (same wind units as `v`).
    nan_policy : {"propagate", "omit"}, default="propagate"
        - ``"propagate"``: RI is NaN where inputs are NaN or where V_PI is invalid.
        - ``"omit"``: returns the ratio where possible, otherwise NaN (same as propagate
          for elementwise ratios; kept for API symmetry with other utilities).
    clip : (low, high), optional
        If provided, clip RI to the interval [low, high].

    Returns
    -------
    float or numpy.ndarray
        Relative intensity ν (unitless).
    """
    v_arr, v_pot_arr = np.broadcast_arrays(
        np.asarray(v, dtype=float), np.asarray(v_pot, dtype=float)
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        ri = v_arr / v_pot_arr

    # Mask nonphysical / undefined ratios.
    bad = (~np.isfinite(v_arr)) | (~np.isfinite(v_pot_arr)) | (v_pot_arr <= 0.0)
    ri = np.where(bad, np.nan, ri)

    if nan_policy not in {"propagate", "omit"}:
        raise ValueError("Invalid nan_policy; expected 'propagate' or 'omit'.")

    if clip is not None:
        lo, hi = clip
        ri = np.clip(ri, lo, hi)

    return float(ri) if ri.ndim == 0 else ri
