"""Power Dissipation Index (PDI) utilities.

PDI is commonly defined as the time-integral of the cube of maximum sustained
surface wind speed over a tropical cyclone's lifetime:

    PDI = ∫ Vmax(t)^3 dt

This module provides simple hard-coded variants for output scaling and unit
conventions. The formulation labels use ``e05`` to reflect Emanuel (2005)-style
usage and plotting conventions.
"""

from typing import Any

import numpy as np

from .numba import njit


@njit()
def _pdi_numba(vmax_cubed_sum_dt, formulation_flag=0):
    """Numba-compiled PDI scaling implementation.

    Parameters
    ----------
    vmax_cubed_sum_dt : float
        The core integral/sum in SI-like units before reporting/scaling.
    formulation_flag : int, default=0
        - 0: "e05_si" (return raw ∑ V^3 dt)
        - 1: "e05_1e11" (return (∑ V^3 dt) / 1e11)
    """
    if formulation_flag == 0:
        return vmax_cubed_sum_dt
    if formulation_flag == 1:
        return vmax_cubed_sum_dt / 1.0e11
    raise ValueError("Invalid formulation_flag for PDI.")


def power_dissipation_index(
    vmax: Any,
    dt: Any,
    *,
    formulation: str = "e05_si",
    wind_units: str = "m/s",
    dt_units: str = "s",
    nan_policy: str = "propagate",
):
    """Compute the Power Dissipation Index (PDI) from a time series.

    Parameters
    ----------
    vmax : array-like
        Maximum sustained surface wind speed time series.
    dt : float or array-like
        Time step(s) corresponding to `vmax`. May be a scalar or an array broadcastable
        to the same shape as `vmax`.
    formulation : {"e05_si", "e05_1e11"}, default="e05_si"
        Output scaling convention:

        - ``"e05_si"``: return raw ``sum(Vmax^3 * dt)`` in (m^3/s^2) if inputs are m/s and seconds.
        - ``"e05_1e11"``: return ``sum(Vmax^3 * dt) / 1e11`` (common plotting scale).
    wind_units : {"m/s", "kt"}, default="m/s"
        Units of `vmax`. Converted to m/s internally.
    dt_units : {"s", "h"}, default="s"
        Units of `dt`. Converted to seconds internally.
    nan_policy : {"propagate", "omit"}, default="propagate"
        How to handle NaNs in the time series:

        - ``"propagate"``: any NaN in ``V^3 * dt`` yields NaN PDI (default).
        - ``"omit"``: drop NaN contributions via ``np.nansum``.

    Returns
    -------
    float
        Power Dissipation Index (PDI) according to the selected formulation.
    """
    if formulation == "e05_si":
        formulation_flag = 0
    elif formulation == "e05_1e11":
        formulation_flag = 1
    else:
        raise ValueError(
            f"Invalid formulation={formulation!r}; expected one of: 'e05_si', 'e05_1e11'."
        )

    vmax_arr = np.asarray(vmax, dtype=float)
    dt_arr = np.asarray(dt, dtype=float)

    if wind_units == "m/s":
        vmax_ms = vmax_arr
    elif wind_units == "kt":
        vmax_ms = vmax_arr * 0.5144444444444445
    else:
        raise ValueError(f"Invalid wind_units={wind_units!r}; expected 'm/s' or 'kt'.")

    if dt_units == "s":
        dt_s = dt_arr
    elif dt_units == "h":
        dt_s = dt_arr * 3600.0
    else:
        raise ValueError(f"Invalid dt_units={dt_units!r}; expected 's' or 'h'.")

    if np.any(dt_s < 0):
        raise ValueError("dt must be non-negative.")

    vmax_cubed = vmax_ms**3.0
    vmax_cubed, dt_s = np.broadcast_arrays(vmax_cubed, dt_s)
    contrib = vmax_cubed * dt_s

    if nan_policy == "propagate":
        core = float(np.sum(contrib))
    elif nan_policy == "omit":
        core = float(np.nansum(contrib))
    else:
        raise ValueError("Invalid nan_policy; expected 'propagate' or 'omit'.")

    return float(_pdi_numba(core, formulation_flag=formulation_flag))


# Short alias matching the module name.
pdi = power_dissipation_index

