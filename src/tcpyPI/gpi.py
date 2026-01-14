"""Genesis Potential Index (GPI) diagnostics.

This module provides hard-coded GPI formulations behind a simple switch.

GPI is an environment-only genesis metric typically combining:
- low-level absolute vorticity,
- midlevel humidity,
- vertical wind shear,
- potential intensity.

Different publications use different exponents/normalizations, so we expose
`formulation` as an explicit selector.
"""

from typing import Any

import numpy as np

from .numba import njit


@njit()
def _gpi_numba(abs_vort, rh_mid, v_shear, v_pot, formulation_flag=0):
    """Numba-compiled GPI implementation.

    Parameters
    ----------
    formulation_flag : int, default=0
        - 0: "en04" (Emanuel & Nolan 2004-style)
        - 1: "c07" (Camargo et al. 2007-style; PI thresholded)
    """
    if (abs_vort < 0.0) or np.isnan(abs_vort):
        # Absolute vorticity magnitude is expected non-negative.
        return np.nan

    # Common normalizations used in GPI formulations:
    # - absolute vorticity scaled by 1e5 s^-1
    # - midlevel RH scaled by 50%
    # - potential intensity scaled by 70 m/s
    # - shear damped by (1 + 0.1*Vshear)^-2
    eta = 1.0e5 * abs_vort  # scale to O(1) magnitude for typical tropical values
    rh = rh_mid / 50.0      # normalize humidity
    shear_term = (1.0 + 0.1 * v_shear) ** (-2.0)  # stronger shear reduces genesis likelihood

    if formulation_flag == 0:
        # Emanuel & Nolan (2004)-style: (Vp/70)^3 with no explicit PI threshold.
        vp = v_pot / 70.0
        return (eta**3.0) * (rh**3.0) * (vp**3.0) * shear_term

    if formulation_flag == 1:
        # Camargo et al. (2007)-style variant: uses a PI threshold (Vp - 35 m/s).
        # The subtraction penalizes environments with low PI that are unlikely to support genesis.
        vp_eff = max(v_pot - 35.0, 0.0) / 70.0
        return (eta**3.0) * (rh**3.0) * (vp_eff**3.0) * shear_term

    raise ValueError("Invalid formulation_flag for GPI.")


def genesis_potential_index(
    abs_vort: Any,
    rh_mid: Any,
    v_shear: Any,
    v_pot: Any,
    *,
    formulation: str = "en04",
):
    """Compute the Genesis Potential Index (GPI).

    Parameters
    ----------
    abs_vort : float or array-like
        Absolute vorticity magnitude at low levels (commonly 850 hPa) in s^-1.
        Many published GPIs use ``|η|`` (magnitude).
    rh_mid : float or array-like
        Midlevel relative humidity (%), commonly near 700 hPa.
    v_shear : float or array-like
        Deep-layer vertical wind shear magnitude (m/s), commonly 850–200 hPa.
    v_pot : float or array-like
        Potential intensity wind speed (m/s), typically from :func:`tcpyPI.pi`.
    formulation : {"en04", "c07"}, default="en04"
        Which hard-coded formulation to use:

        - ``"en04"``: Emanuel & Nolan (2004)-style
        - ``"c07"``: Camargo et al. (2007)-style variant with PI thresholding

    Returns
    -------
    float or numpy.ndarray
        GPI (unitless, relative index).

    Notes
    -----
    GPI is not a probability; it is an empirical index intended for spatial/seasonal
    comparisons. Ensure you use inputs at the levels and with the units expected
    by the chosen formulation.
    """
    if formulation == "en04":
        formulation_flag = 0
    elif formulation == "c07":
        formulation_flag = 1
    else:
        raise ValueError(
            f"Invalid formulation={formulation!r}; expected one of: 'en04', 'c07'."
        )

    abs_vort_arr, rh_mid_arr, v_shear_arr, v_pot_arr = np.broadcast_arrays(
        np.asarray(abs_vort, dtype=float),
        np.asarray(rh_mid, dtype=float),
        np.asarray(v_shear, dtype=float),
        np.asarray(v_pot, dtype=float),
    )

    out = np.empty_like(abs_vort_arr, dtype=float)
    it = np.nditer(
        [abs_vort_arr, rh_mid_arr, v_shear_arr, v_pot_arr, out],
        op_flags=[["readonly"]] * 4 + [["writeonly"]],
    )
    for eta, rh, vs, vp, o in it:
        o[...] = _gpi_numba(
            float(eta),
            float(rh),
            float(vs),
            float(vp),
            formulation_flag=formulation_flag,
        )

    return float(out) if out.ndim == 0 else out


# Short alias matching the module name.
gpi = genesis_potential_index

