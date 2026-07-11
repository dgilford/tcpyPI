"""Genesis Potential Index (GPI) diagnostics.

This module provides hard-coded GPI formulations behind a simple switch.

GPI is an environment-only genesis metric typically combining:
- low-level absolute vorticity,
- midlevel humidity (or midlevel moist-entropy deficit),
- vertical wind shear,
- potential intensity.

Different publications use different exponents/normalizations and even different
thermodynamic ingredients, so we expose `formulation` as an explicit selector:

- ``"en04"``: Emanuel & Nolan (2004) / Camargo et al. (2007), using midlevel RH.
- ``"e10"``: Emanuel (2010, JAMES, Eq. 11), using the nondimensional midlevel
  moist-entropy deficit ``chi_m`` (Tang & Emanuel 2012, Eq. 2) instead of RH.

GPI is a relative index (defined up to a constant of proportionality); interpret
magnitudes comparatively (spatially/seasonally), not absolutely.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .numba import njit


@njit()
def _gpi_en04_numba(abs_vort, rh_mid, v_shear, v_pot):
    """Emanuel & Nolan (2004) / Camargo et al. (2007) GPI.

    GPI = |1e5 eta|^{3/2} (RH/50)^3 (Vp/70)^3 (1 + 0.1 Vshear)^{-2}
    """
    if np.isnan(abs_vort) or np.isnan(rh_mid) or np.isnan(v_shear) or np.isnan(v_pot):
        return np.nan

    # Published GPIs use the magnitude |eta| (e.g. |1e5 eta|^{3/2}), so take abs()
    # to accept signed absolute vorticity (f+zeta is negative in the S. Hemisphere).
    eta = 1.0e5 * abs(abs_vort)  # scaled to O(1) for typical tropical values
    rh = rh_mid / 50.0  # normalize humidity
    vp = v_pot / 70.0  # normalize potential intensity
    shear_term = (1.0 + 0.1 * v_shear) ** (-2.0)  # stronger shear reduces genesis
    return (eta**1.5) * (rh**3.0) * (vp**3.0) * shear_term


@njit()
def _gpi_e10_numba(abs_vort, chi, v_shear, v_pot):
    """Emanuel (2010, JAMES, Eq. 11) GPI.

    GPI = |eta|^3 * chi_m^{-4/3} * max(Vp - 35, 0)^2 * (25 + Vshear)^{-4}

    where ``chi_m`` is the nondimensional midlevel (600 hPa) moist-entropy deficit
    of Tang & Emanuel (2012, Eq. 2): chi_m = (s*_m - s_m) / (s*_SST - s_b).
    The shear is the 850-250 hPa magnitude and ``eta`` is the 850 hPa absolute
    vorticity (s^-1). This form replaces the EN04 relative-humidity factor with a
    thermodynamically motivated entropy deficit and does not use the 1e5/70/50
    normalizations, so its absolute magnitude differs from ``en04`` (both are
    relative indices).
    """
    if np.isnan(abs_vort) or np.isnan(chi) or np.isnan(v_shear) or np.isnan(v_pot):
        return np.nan
    if chi <= 0.0:
        # chi^{-4/3} is undefined for a non-positive entropy deficit.
        return np.nan

    eta = abs(abs_vort)
    chi_term = chi ** (-4.0 / 3.0)
    vp_term = max(v_pot - 35.0, 0.0) ** 2.0
    shear_term = (25.0 + v_shear) ** (-4.0)
    return (eta**3.0) * chi_term * vp_term * shear_term


def _broadcast_apply(kernel, a, b, c, d):
    """Broadcast four inputs and apply a scalar GPI ``kernel`` elementwise."""
    a_arr, b_arr, c_arr, d_arr = np.broadcast_arrays(
        np.asarray(a, dtype=float),
        np.asarray(b, dtype=float),
        np.asarray(c, dtype=float),
        np.asarray(d, dtype=float),
    )
    out = np.empty_like(a_arr, dtype=float)
    it = np.nditer(
        [a_arr, b_arr, c_arr, d_arr, out],
        op_flags=[["readonly"]] * 4 + [["writeonly"]],
    )
    for av, bv, cv, dv, o in it:
        o[...] = kernel(float(av), float(bv), float(cv), float(dv))
    return float(out) if out.ndim == 0 else out


def genesis_potential_index(
    abs_vort: Any,
    rh_mid: Any = None,
    v_shear: Any = None,
    v_pot: Any = None,
    *,
    formulation: str = "en04",
    chi: Any | None = None,
):
    """Compute the Genesis Potential Index (GPI).

    Parameters
    ----------
    abs_vort : float or array-like
        Absolute vorticity at low levels (commonly 850 hPa) in s^-1. Published GPIs
        use the magnitude ``|eta|``, so the sign is ignored (``|abs_vort|`` is used
        internally) and signed input such as ``f + zeta`` works in both hemispheres.
    rh_mid : float or array-like, optional
        Midlevel relative humidity (%), commonly near 700 hPa. Required for
        ``formulation="en04"``; unused by ``"e10"``.
    v_shear : float or array-like
        Deep-layer vertical wind shear magnitude (m/s), commonly 850-200 hPa
        (``en04``) or 850-250 hPa (``e10``).
    v_pot : float or array-like
        Potential intensity wind speed (m/s), typically from :func:`tcpyPI.pi`.
    formulation : {"en04", "e10"}, default="en04"
        Which formulation to use:

        - ``"en04"``: Emanuel & Nolan (2004) / Camargo et al. (2007), RH-based.
        - ``"e10"``: Emanuel (2010, JAMES, Eq. 11), entropy-deficit-based
          (requires ``chi``).
    chi : float or array-like, optional
        Nondimensional midlevel (600 hPa) moist-entropy deficit ``chi_m`` of Tang &
        Emanuel (2012, Eq. 2). Required for ``formulation="e10"``; unused by
        ``"en04"``. Compute e.g. with
        :func:`tcpyPI.entropy_deficit_te12_from_profile`.

    Returns
    -------
    float or numpy.ndarray
        GPI (unitless, relative index).

    Notes
    -----
    GPI is not a probability; it is an empirical index intended for spatial/seasonal
    comparisons and is defined up to a constant of proportionality (so ``en04`` and
    ``e10`` magnitudes are not directly comparable to each other). Ensure you use
    inputs at the levels and with the units expected by the chosen formulation.

    Examples
    --------
        >>> import tcpyPI
        >>> round(float(tcpyPI.genesis_potential_index(
        ...     2.0e-5, rh_mid=60.0, v_shear=10.0, v_pot=60.0, formulation="en04")), 6)
        0.769464
        >>> tcpyPI.genesis_potential_index(
        ...     2.0e-5, v_shear=10.0, v_pot=30.0, formulation="e10", chi=1.0)
        0.0
    """
    if formulation == "en04":
        if rh_mid is None:
            raise ValueError("formulation='en04' requires rh_mid (midlevel RH, %).")
        if v_shear is None or v_pot is None:
            raise ValueError("v_shear and v_pot are required.")
        return _broadcast_apply(_gpi_en04_numba, abs_vort, rh_mid, v_shear, v_pot)

    if formulation == "e10":
        if chi is None:
            raise ValueError(
                "formulation='e10' requires chi (nondimensional midlevel moist-entropy "
                "deficit chi_m; e.g. from tcpyPI.entropy_deficit_te12_from_profile)."
            )
        if v_shear is None or v_pot is None:
            raise ValueError("v_shear and v_pot are required.")
        return _broadcast_apply(_gpi_e10_numba, abs_vort, chi, v_shear, v_pot)

    raise ValueError(f"Invalid formulation={formulation!r}; expected one of: 'en04', 'e10'.")


# Short alias matching the module name.
gpi = genesis_potential_index


def gpi_log_decomposition(
    abs_vort: Any,
    rh_mid: Any = None,
    v_shear: Any = None,
    v_pot: Any = None,
    *,
    formulation: str = "en04",
    chi: Any | None = None,
):
    """Additive log-space decomposition of the Genesis Potential Index.

    GPI is multiplicative, so it separates additively in log space (cf. the PI and
    ventilation-index log-decompositions).

    For ``"en04"``:

    .. math::
        \\ln\\mathrm{GPI} = \\tfrac{3}{2}\\ln(10^5|\\eta|) + 3\\ln\\tfrac{RH}{50}
        + 3\\ln\\tfrac{V_{PI}}{70} - 2\\ln(1+0.1\\,V_{shear})

    For ``"e10"`` (Emanuel 2010):

    .. math::
        \\ln\\mathrm{GPI} = 3\\ln|\\eta| - \\tfrac{4}{3}\\ln\\chi_m
        + 2\\ln\\max(V_{PI}-35,0) - 4\\ln(25+V_{shear})

    Parameters match :func:`genesis_potential_index` (``chi`` required for ``e10``).

    Returns
    -------
    dict
        The signed additive contributions that **sum to** ``lngpi``. Keys:

        - ``en04``: ``vorticity``, ``humidity``, ``pi``, ``shear``, ``lngpi``
        - ``e10`` : ``vorticity``, ``chi``, ``pi``, ``shear``, ``lngpi``

        Values are ``nan`` where a term is undefined (non-positive argument to a
        log; for ``e10`` this includes ``V_PI <= 35``).

    Examples
    --------
        >>> import tcpyPI
        >>> d = tcpyPI.gpi_log_decomposition(2.0e-5, rh_mid=60.0, v_shear=10.0, v_pot=60.0)
        >>> round(d["lngpi"] - (d["vorticity"] + d["humidity"] + d["pi"] + d["shear"]), 12)
        0.0
    """

    def _s(a):
        return float(a) if a.ndim == 0 else a

    if formulation == "en04":
        if rh_mid is None:
            raise ValueError("formulation='en04' requires rh_mid (midlevel RH, %).")
        if v_shear is None or v_pot is None:
            raise ValueError("v_shear and v_pot are required.")
        eta, rh, vsh, vp = np.broadcast_arrays(
            np.abs(np.asarray(abs_vort, dtype=float)),
            np.asarray(rh_mid, dtype=float),
            np.asarray(v_shear, dtype=float),
            np.asarray(v_pot, dtype=float),
        )
        bad = (eta <= 0) | (rh <= 0) | (vp <= 0) | (1.0 + 0.1 * vsh <= 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            vort = np.where(bad, np.nan, 1.5 * np.log(1.0e5 * eta))
            hum = np.where(bad, np.nan, 3.0 * np.log(rh / 50.0))
            pi_t = np.where(bad, np.nan, 3.0 * np.log(vp / 70.0))
            sh = np.where(bad, np.nan, -2.0 * np.log(1.0 + 0.1 * vsh))
        lngpi = vort + hum + pi_t + sh
        return {
            "lngpi": _s(lngpi),
            "vorticity": _s(vort),
            "humidity": _s(hum),
            "pi": _s(pi_t),
            "shear": _s(sh),
        }

    if formulation == "e10":
        if chi is None:
            raise ValueError(
                "formulation='e10' requires chi (nondimensional midlevel moist-entropy "
                "deficit chi_m; e.g. from tcpyPI.entropy_deficit_te12_from_profile)."
            )
        if v_shear is None or v_pot is None:
            raise ValueError("v_shear and v_pot are required.")
        eta, ch, vsh, vp = np.broadcast_arrays(
            np.abs(np.asarray(abs_vort, dtype=float)),
            np.asarray(chi, dtype=float),
            np.asarray(v_shear, dtype=float),
            np.asarray(v_pot, dtype=float),
        )
        vp_excess = vp - 35.0
        bad = (eta <= 0) | (ch <= 0) | (vp_excess <= 0) | (25.0 + vsh <= 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            vort = np.where(bad, np.nan, 3.0 * np.log(eta))
            chi_t = np.where(bad, np.nan, (-4.0 / 3.0) * np.log(ch))
            pi_t = np.where(bad, np.nan, 2.0 * np.log(vp_excess))
            sh = np.where(bad, np.nan, -4.0 * np.log(25.0 + vsh))
        lngpi = vort + chi_t + pi_t + sh
        return {
            "lngpi": _s(lngpi),
            "vorticity": _s(vort),
            "chi": _s(chi_t),
            "pi": _s(pi_t),
            "shear": _s(sh),
        }

    raise ValueError(f"Invalid formulation={formulation!r}; expected one of: 'en04', 'e10'.")
