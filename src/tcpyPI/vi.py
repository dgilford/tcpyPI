"""Ventilation index (VI) diagnostics.

This module provides a hard-coded "ventilation index" formulation commonly used
in tropical cyclone genesis/intensity diagnostics.

Design
------
- Public API accepts a string `formulation` so additional published variants can
  be added later.
- Currently, only the Tang & Emanuel (2012) ventilation index is implemented.

Tang & Emanuel (2012)
---------------------
Tang & Emanuel (2012) define a ventilation index ``Λ`` and midlevel entropy deficit
``χm`` (Eqs. 1–2).

The core definitions are:

- Ventilation index:
    ``Λ = (u_shear / u_PI) * χm``

- Entropy deficit:
    ``χm = (s*m - sm) / (s*SST - sb)``

where ``u_shear`` is bulk (850–200 hPa) shear magnitude, ``u_PI`` is potential
intensity, ``sm`` is environmental entropy at a midlevel (commonly 600 hPa),
``s*m`` is saturation entropy at that midlevel in the storm core, ``s*SST`` is
saturation entropy at the sea surface temperature, and ``sb`` is boundary-layer
entropy.
"""

from typing import Any

import numpy as np

from . import constants, utilities
from .numba import njit


@njit()
def _vi_numba(v_shear, v_pot, chi_m, formulation_flag=0):
    """Numba-compiled VI implementation.

    Parameters
    ----------
    formulation_flag : int, default=0
        - 0: "te12" (Λ = (u_shear / u_PI) * χm)
    """
    if formulation_flag != 0:
        # Numba-friendly error (no f-strings).
        raise ValueError("Invalid formulation_flag for VI.")

    # Scientific rationale:
    # - Stronger environmental shear ventilates the core, hindering development.
    # - Larger potential intensity Vp implies a more favorable thermodynamic ceiling.
    # - The entropy deficit term χ penalizes dry/low-entropy midlevels.
    #
    # Tang & Emanuel (2012): Λ = (u_shear / u_PI) * χm
    if (v_pot <= 0.0) or np.isnan(v_pot):
        return np.nan
    return (v_shear / v_pot) * chi_m


@njit()
def _chi_m_te12_numba(s_m_star, s_m, s_sst_star, s_b):
    """Numba-compiled entropy deficit χm (Tang & Emanuel 2012).

    Implements:
        χm = (s*m - sm) / (s*SST - sb)
    """
    denom = s_sst_star - s_b
    if denom == 0.0 or np.isnan(denom):
        return np.nan
    return (s_m_star - s_m) / denom


def entropy_deficit_te12(
    s_m_star: Any,
    s_m: Any,
    s_sst_star: Any,
    s_b: Any,
):
    """Compute the midlevel entropy deficit χm (Tang & Emanuel 2012).

    Parameters
    ----------
    s_m_star : float or array-like
        Saturation entropy at the midlevel of interest in the storm inner core
        (often 600 hPa).
    s_m : float or array-like
        Environmental entropy at the same midlevel.
    s_sst_star : float or array-like
        Saturation entropy at the sea surface temperature.
    s_b : float or array-like
        Boundary-layer entropy.

    Returns
    -------
    float or numpy.ndarray
        χm (unitless).
    """
    s_m_star_arr, s_m_arr, s_sst_star_arr, s_b_arr = np.broadcast_arrays(
        np.asarray(s_m_star, dtype=float),
        np.asarray(s_m, dtype=float),
        np.asarray(s_sst_star, dtype=float),
        np.asarray(s_b, dtype=float),
    )

    out = np.empty_like(s_m_star_arr, dtype=float)
    it = np.nditer(
        [s_m_star_arr, s_m_arr, s_sst_star_arr, s_b_arr, out],
        op_flags=[["readonly"]] * 4 + [["writeonly"]],  # type: ignore[arg-type]  # per-operand flag lists are valid at runtime
    )
    for sm_star, sm, sst_star, sb, o in it:
        o[...] = _chi_m_te12_numba(float(sm_star), float(sm), float(sst_star), float(sb))
    return float(out) if out.ndim == 0 else out


@njit()
def _interp_linear_numba(x, y, x_target):
    """1D linear interpolation for monotonic x (Numba-friendly).

    This assumes `x` is monotonically decreasing (e.g., pressure levels from surface
    upward). If `x_target` is outside the bounds, returns NaN.
    """
    n = len(x)
    if n < 2:
        return np.nan

    x0 = x[0]
    x1 = x[n - 1]
    # For decreasing x: x0 >= ... >= x1
    if (x_target > x0) or (x_target < x1):
        return np.nan

    for i in range(n - 1):
        xa = x[i]
        xb = x[i + 1]
        if (xa >= x_target) and (x_target >= xb):
            ya = y[i]
            yb = y[i + 1]
            if xa == xb:
                return ya
            w = (x_target - xa) / (xb - xa)
            return ya + w * (yb - ya)

    return np.nan


@njit()
def _sat_mixing_ratio_numba(Tk, P_hpa):
    """Saturation mixing ratio (g/g) at (Tk, P_hpa) using Clausius-Clapeyron."""
    Tc = utilities.T_ktoC(Tk)
    es = utilities.es_cc(Tc)
    return utilities.rv(es, P_hpa)


@njit()
def _sat_entropy_numba(Tk, P_hpa, entropy_method_flag):
    """Saturation entropy (RH=100%) at (Tk, P_hpa) for the selected entropy method."""
    R_sat = _sat_mixing_ratio_numba(Tk, P_hpa)
    if entropy_method_flag == 0:
        return _entropy_method_emanuel94_numba(Tk, R_sat, P_hpa)
    return _entropy_method_bryan2008_numba(Tk, R_sat, P_hpa, 1.0)


@njit()
def _moist_adiabat_T_at_p_numba(
    s_sat_surface,
    p_target_hpa,
    entropy_method_flag,
    T_low_k=150.0,
    T_high_k=350.0,
    tol_s=1.0e-6,
    max_iter=80,
):
    """Invert saturation entropy to find the moist-adiabatic temperature at p_target.

    Finds T at `p_target_hpa` such that:
        s_sat(T, p_target) = s_sat_surface

    This treats a saturated moist adiabat as a curve of constant saturation entropy,
    matching the entropy method selected by `entropy_method_flag`.
    """
    s_low = _sat_entropy_numba(T_low_k, p_target_hpa, entropy_method_flag)
    if np.isnan(s_low):
        return np.nan

    # Some entropy formulations become undefined if saturation vapor pressure exceeds
    # ambient pressure; search downward for a valid high bracket.
    hi = T_high_k
    s_high = _sat_entropy_numba(hi, p_target_hpa, entropy_method_flag)
    for _ in range(60):
        if not np.isnan(s_high):
            break
        hi -= 5.0
        if hi <= T_low_k + 1.0:
            return np.nan
        s_high = _sat_entropy_numba(hi, p_target_hpa, entropy_method_flag)
    if np.isnan(s_high):
        return np.nan
    if (s_sat_surface < s_low) or (s_sat_surface > s_high):
        return np.nan

    lo = T_low_k
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        s_mid = _sat_entropy_numba(mid, p_target_hpa, entropy_method_flag)
        if np.isnan(s_mid):
            return np.nan
        if np.abs(s_mid - s_sat_surface) <= tol_s:
            return mid
        if s_mid < s_sat_surface:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@njit()
def _entropy_method_emanuel94_numba(Tk, R_gg, P_hpa):
    """Moist entropy per unit mass dry air as used elsewhere in tcpyPI (`entropy_S`)."""
    return utilities.entropy_S(Tk, R_gg, P_hpa)


@njit()
def _entropy_method_bryan2008_numba(Tk, R_gg, P_hpa, RH):
    """Pseudo-adiabatic entropy used by Tang & Emanuel (2012) (their Eq. 3).

    This follows the variable definitions given in TE12 immediately after Eq. (3):
    - cp: specific heat at constant pressure for dry air
    - T: temperature (K)
    - Rd: gas constant for dry air
    - pd: partial pressure of dry air
    - Lvo: latent heat of vaporization (set to 2.555e6 J/kg in TE12)
    - rv: water vapor mixing ratio (g/g)
    - Rv: gas constant for water vapor
    - H: relative humidity (0–1)

    Notes
    -----
    TE12 note: Lvo is set to 2.555e6 J/kg to compensate for neglecting the entropy
    of water vapor.
    """
    # Clamp RH to physically meaningful bounds.
    if RH <= 0.0 or np.isnan(RH):
        return np.nan
    if RH > 1.0:
        RH = 1.0

    # Vapor pressure from mixing ratio, using TE12 variable naming.
    e = utilities.ev(R_gg, P_hpa)
    pd = P_hpa - e
    if pd <= 0.0 or np.isnan(pd):
        return np.nan

    Lvo = 2.555e6
    return (
        constants.CPD * np.log(Tk)
        - constants.RD * np.log(pd)
        + (Lvo * R_gg) / Tk
        - R_gg * constants.RV * np.log(RH)
    )


@njit()
def _relative_humidity_numba(Tk, R_gg, P_hpa):
    """Compute RH from (T, r, p) using the same CC saturation formula as tcpyPI."""
    e = utilities.ev(R_gg, P_hpa)
    es = utilities.es_cc(utilities.T_ktoC(Tk))
    if es <= 0.0 or np.isnan(es):
        return np.nan
    rh = e / es
    if rh > 1.0:
        rh = 1.0
    if rh < 0.0:
        rh = 0.0
    return rh


@njit()
def _chi_m_te12_from_profile_numba(
    P_hpa,
    T_k,
    R_gg,
    SST_k,
    psfc_hpa,
    T2m_k,
    R2m_gg,
    entropy_method_flag=0,
    sb_source_flag=0,
    s_m_star_source_flag=0,
    p_mid_hpa=600.0,
):
    """Compute χm from raw profile inputs (Numba core).

    Parameters
    ----------
    sb_source_flag : int, default=0
        - 0: "t2m" (preferred; uses T2m/R2m at psfc_hpa; expects provided)
        - 1: "layer_1000_900" (mean entropy over 1000–900 hPa; fallback to lowest level)
        - 2: "lowest_level" (uses first profile level)
    s_m_star_source_flag : int, default=0
        - 0: "env_T_at_pmid" (default; uses environmental T at p_mid with RH=100%)
        - 1: "moist_adiabat_from_sst" (lift saturated parcel from SST at psfc to p_mid)
    """
    # 1) Midlevel environmental entropy sm at p_mid.
    T_mid = _interp_linear_numba(P_hpa, T_k, p_mid_hpa)
    R_mid = _interp_linear_numba(P_hpa, R_gg, p_mid_hpa)
    if np.isnan(T_mid) or np.isnan(R_mid):
        return np.nan

    if entropy_method_flag == 0:
        s_m = _entropy_method_emanuel94_numba(T_mid, R_mid, p_mid_hpa)
    else:
        RH_mid = _relative_humidity_numba(T_mid, R_mid, p_mid_hpa)
        s_m = _entropy_method_bryan2008_numba(T_mid, R_mid, p_mid_hpa, RH_mid)

    # 2) Saturation entropy at SST at the surface pressure s*SST.
    # This is the surface saturated entropy used by the "moist_adiabat_from_sst"
    # pathway (constant saturation entropy ascent).
    s_sst_star = _sat_entropy_numba(SST_k, psfc_hpa, entropy_method_flag)
    if np.isnan(s_sst_star):
        return np.nan

    # 3) Saturation entropy at midlevel s*m.
    if s_m_star_source_flag == 0:
        # "env_T_at_pmid": use environmental T at p_mid, set RH=100%.
        s_m_star = _sat_entropy_numba(T_mid, p_mid_hpa, entropy_method_flag)
    else:
        # "moist_adiabat_from_sst": lift a saturated parcel from SST at psfc to p_mid.
        T_mid_star = _moist_adiabat_T_at_p_numba(
            s_sst_star,
            p_mid_hpa,
            entropy_method_flag,
            T_low_k=150.0,
            T_high_k=400.0,
        )
        if np.isnan(T_mid_star):
            return np.nan
        s_m_star = _sat_entropy_numba(T_mid_star, p_mid_hpa, entropy_method_flag)
        if np.isnan(s_m_star):
            return np.nan

    # 4) Boundary-layer entropy sb.
    if sb_source_flag == 0:
        # "t2m": use 2m air properties (T2m, q2m) at surface pressure.
        if np.isnan(T2m_k) or np.isnan(R2m_gg):
            return np.nan
        if entropy_method_flag == 0:
            s_b = _entropy_method_emanuel94_numba(T2m_k, R2m_gg, psfc_hpa)
        else:
            RH_b = _relative_humidity_numba(T2m_k, R2m_gg, psfc_hpa)
            s_b = _entropy_method_bryan2008_numba(T2m_k, R2m_gg, psfc_hpa, RH_b)
    elif sb_source_flag == 2:
        # "lowest_level": use first profile level.
        if entropy_method_flag == 0:
            s_b = _entropy_method_emanuel94_numba(T_k[0], R_gg[0], P_hpa[0])
        else:
            RH0 = _relative_humidity_numba(T_k[0], R_gg[0], P_hpa[0])
            s_b = _entropy_method_bryan2008_numba(T_k[0], R_gg[0], P_hpa[0], RH0)
    else:
        # "layer_1000_900": average entropy between 1000 and 900 hPa (inclusive).
        s_sum = 0.0
        n = 0
        for i in range(len(P_hpa)):
            p = P_hpa[i]
            if (p <= 1000.0) and (p >= 900.0):
                if entropy_method_flag == 0:
                    s_sum += _entropy_method_emanuel94_numba(T_k[i], R_gg[i], p)
                else:
                    RHi = _relative_humidity_numba(T_k[i], R_gg[i], p)
                    s_sum += _entropy_method_bryan2008_numba(T_k[i], R_gg[i], p, RHi)
                n += 1
        if n == 0:
            # fallback: lowest level
            if entropy_method_flag == 0:
                s_b = _entropy_method_emanuel94_numba(T_k[0], R_gg[0], P_hpa[0])
            else:
                RH0 = _relative_humidity_numba(T_k[0], R_gg[0], P_hpa[0])
                s_b = _entropy_method_bryan2008_numba(T_k[0], R_gg[0], P_hpa[0], RH0)
        else:
            s_b = s_sum / n

    return _chi_m_te12_numba(s_m_star, s_m, s_sst_star, s_b)


def entropy_deficit_te12_from_profile(
    P_hpa: Any,
    T: Any,
    q: Any,
    *,
    T_units: str = "K",
    q_units: str = "g/kg",
    SST: Any | None = None,
    SST_units: str = "K",
    psfc_hpa: Any | None = None,
    T2m: Any | None = None,
    q2m: Any | None = None,
    sb_source: str = "t2m",
    sb_fallback: tuple[str, ...] = ("layer_1000_900", "lowest_level"),
    p_mid_hpa: float = 600.0,
    s_m_star_source: str = "env_T_at_pmid",
    entropy_method: str = "emanuel94",
):
    """Compute χm from raw profile inputs using TE12 Eq. (2).

    This computes the entropy deficit term:
        χm = (s*m − sm) / (s*SST − sb)

    Parameters
    ----------
    P_hpa : array-like
        Pressure levels (hPa), ordered from high pressure (near-surface) to low pressure.
    T : array-like
        Temperature profile. Units given by `T_units`.
    q : array-like
        Specific humidity or mixing ratio profile. Units given by `q_units`.
        For consistency with `tcpyPI.pi`, defaults assume `q` is mixing ratio in g/kg.
    T_units : {"K", "C"}, default="K"
        Units of `T`.
    q_units : {"g/kg", "g/g"}, default="g/kg"
        Units of `q`.
    SST : float, optional
        Sea surface temperature. If omitted, uses `T2m` if provided.
    SST_units : {"K", "C"}, default="K"
        Units of `SST`.
    psfc_hpa : float, optional
        Surface pressure (hPa). If omitted, uses the first element of `P_hpa`.
    T2m : float, optional
        2m air temperature.
    q2m : float, optional
        2m humidity in units given by `q_units`.
    sb_source : {"t2m", "layer_1000_900", "lowest_level"}, default="t2m"
        Primary source for boundary-layer entropy `sb`.
    sb_fallback : tuple of str, default=("layer_1000_900","lowest_level")
        Fallbacks if `sb_source="t2m"` but `T2m/q2m` are not provided.
    p_mid_hpa : float, default=600.0
        Midlevel pressure used in TE12 (hPa).
    s_m_star_source : {"env_T_at_pmid", "moist_adiabat_from_sst"}, default="env_T_at_pmid"
        How to compute the inner-core saturation entropy at midlevels `s*m`.

        - ``"env_T_at_pmid"``: use environmental `T(p_mid)` with RH=100% (default).
        - ``"moist_adiabat_from_sst"``: lift a saturated parcel from ``SST`` at ``psfc_hpa`` to
          ``p_mid_hpa`` along a moist adiabat (treated as constant saturation entropy), then
          evaluate saturation entropy at that lifted temperature.
    entropy_method : {"emanuel94", "bryan2008"}, default="emanuel94"
        Which entropy formulation to use inside TE12 Eq. (2).

        TE12 state they use the pseudo-adiabatic entropy from Bryan (2008) (their Eq. 3).
        Use ``entropy_method="bryan2008"`` to match TE12; the default remains
        ``"emanuel94"`` for consistency with the rest of tcpyPI.

    Returns
    -------
    float
        χm (unitless).
    """
    P_arr = np.asarray(P_hpa, dtype=float)
    T_arr = np.asarray(T, dtype=float)
    q_arr = np.asarray(q, dtype=float)

    if T_units == "K":
        T_k = T_arr
    elif T_units == "C":
        T_k = T_arr + 273.15
    else:
        raise ValueError(f"Invalid T_units={T_units!r}; expected 'K' or 'C'.")

    if q_units == "g/kg":
        R_gg = q_arr * 0.001
    elif q_units == "g/g":
        R_gg = q_arr
    else:
        raise ValueError(f"Invalid q_units={q_units!r}; expected 'g/kg' or 'g/g'.")

    if psfc_hpa is None:
        psfc = float(P_arr[0])
    else:
        psfc = float(psfc_hpa)

    if SST is None:
        if T2m is None:
            raise ValueError("SST must be provided unless T2m is provided as a fallback.")
        SST_val = T2m
        SST_units_eff = T_units
    else:
        SST_val = SST
        SST_units_eff = SST_units

    SST_val_arr = np.asarray(SST_val, dtype=float)
    if SST_units_eff == "K":
        SST_k = float(SST_val_arr)
    elif SST_units_eff == "C":
        SST_k = float(SST_val_arr + 273.15)
    else:
        raise ValueError(f"Invalid SST_units={SST_units_eff!r}; expected 'K' or 'C'.")

    T2m_k = np.nan
    R2m_gg = np.nan
    if T2m is not None and q2m is not None:
        T2m_arr = np.asarray(T2m, dtype=float)
        q2m_arr = np.asarray(q2m, dtype=float)
        T2m_k = float(T2m_arr if T_units == "K" else (T2m_arr + 273.15))
        R2m_gg = float(q2m_arr * 0.001 if q_units == "g/kg" else q2m_arr)

    # Map sb_source to flag, with fallbacks only for the t2m pathway.
    def _sb_flag(name: str) -> int:
        if name == "t2m":
            return 0
        if name == "layer_1000_900":
            return 1
        if name == "lowest_level":
            return 2
        raise ValueError(
            f"Invalid sb_source={name!r}; expected 't2m', 'layer_1000_900', or 'lowest_level'."
        )

    if sb_source == "t2m" and (T2m is None or q2m is None):
        # Try fallbacks in order.
        used = None
        for cand in sb_fallback:
            used = cand
            sb_flag = _sb_flag(cand)
            break
        if used is None:
            raise ValueError("sb_source='t2m' requires T2m/q2m or a non-empty sb_fallback.")
    else:
        sb_flag = _sb_flag(sb_source)

    if s_m_star_source == "env_T_at_pmid":
        s_m_star_flag = 0
    elif s_m_star_source == "moist_adiabat_from_sst":
        s_m_star_flag = 1
    else:
        raise ValueError(
            "Invalid s_m_star_source; expected 'env_T_at_pmid' or 'moist_adiabat_from_sst'."
        )

    if entropy_method == "emanuel94":
        entropy_method_flag = 0
    elif entropy_method == "bryan2008":
        entropy_method_flag = 1
    else:
        raise ValueError("Invalid entropy_method; expected 'emanuel94' or 'bryan2008'.")

    return float(
        _chi_m_te12_from_profile_numba(
            P_arr,
            T_k,
            R_gg,
            float(SST_k),
            float(psfc),
            float(T2m_k),
            float(R2m_gg),
            entropy_method_flag=entropy_method_flag,
            sb_source_flag=sb_flag,
            s_m_star_source_flag=s_m_star_flag,
            p_mid_hpa=float(p_mid_hpa),
        )
    )


def ventilation_index(
    v_shear: Any,
    v_pot: Any,
    chi_m: Any = None,
    *,
    s_m_star: Any = None,
    s_m: Any = None,
    s_sst_star: Any = None,
    s_b: Any = None,
    formulation: str = "te12",
):
    """Compute the ventilation index (VI).

    Parameters
    ----------
    v_shear : float or array-like
        Deep-layer vertical wind shear magnitude (m/s).
    v_pot : float or array-like
        Potential intensity wind speed (m/s), typically from :func:`tcpyPI.pi`.
    chi_m : float or array-like, optional
        Midlevel entropy deficit term (unitless). If omitted, you must provide
        `s_m_star`, `s_m`, `s_sst_star`, and `s_b` to compute χm.
    s_m_star, s_m, s_sst_star, s_b : float or array-like, optional
        Inputs to compute χm (see :func:`entropy_deficit_te12`).
    formulation : {"te12"}, default="te12"
        Which hard-coded VI formulation to use.

        - ``"te12"``: Λ = (u_shear / u_PI) * χm

    Returns
    -------
    float or numpy.ndarray
        Ventilation index (unitless).

    Notes
    -----
    Citation naming:
    - ``"te12"`` is intended to represent the Tang & Emanuel ventilation framework.
      If you provide the entropy components, χm is computed using the TE12 definition.
    """
    if formulation == "te12":
        formulation_flag = 0
    else:
        raise ValueError(f"Invalid formulation={formulation!r}; expected one of: 'te12'.")

    if chi_m is None:
        missing = [
            name
            for name, value in [
                ("s_m_star", s_m_star),
                ("s_m", s_m),
                ("s_sst_star", s_sst_star),
                ("s_b", s_b),
            ]
            if value is None
        ]
        if missing:
            raise ValueError(
                "chi_m is required unless all of "
                f"s_m_star/s_m/s_sst_star/s_b are provided; missing: {', '.join(missing)}"
            )
        chi_m = entropy_deficit_te12(s_m_star, s_m, s_sst_star, s_b)

    v_shear_arr, v_pot_arr, chi_arr = np.broadcast_arrays(
        np.asarray(v_shear, dtype=float),
        np.asarray(v_pot, dtype=float),
        np.asarray(chi_m, dtype=float),
    )

    # Vectorize by calling the compiled scalar function elementwise.
    # This keeps behavior consistent across scalar and array inputs.
    out = np.empty_like(v_shear_arr, dtype=float)
    it = np.nditer(
        [v_shear_arr, v_pot_arr, chi_arr, out],
        op_flags=[["readonly"]] * 3 + [["writeonly"]],  # type: ignore[arg-type]  # per-operand flag lists are valid at runtime
    )
    for vs, vp, ch, o in it:
        o[...] = _vi_numba(float(vs), float(vp), float(ch), formulation_flag=formulation_flag)
    return float(out) if out.ndim == 0 else out


# Short alias matching the module name.
vi = ventilation_index


def vi_log_decomposition(u_shear, v_pot, chi_m):
    """Additive log-space decomposition of the TE12 ventilation index.

    Because :math:`\\Lambda = (u_{shear}/V_{PI})\\,\\chi_m` is multiplicative, it
    separates additively in log space (cf. the PI log-decomposition):

    .. math:: \\ln\\Lambda = \\ln u_{shear} - \\ln V_{PI} + \\ln\\chi_m

    Parameters
    ----------
    u_shear : float or array-like
        Deep-layer (e.g. 850-200 hPa) vertical wind shear magnitude (m/s).
    v_pot : float or array-like
        Potential intensity wind speed (m/s), e.g. from :func:`tcpyPI.pi`.
    chi_m : float or array-like
        Nondimensional midlevel moist-entropy deficit (TE12), e.g. from
        :func:`entropy_deficit_te12_from_profile`.

    Returns
    -------
    dict
        The signed additive contributions that **sum to** ``lnvi``:

        - ``shear``  = ``+ln(u_shear)``
        - ``pi``     = ``-ln(V_PI)``
        - ``chi``    = ``+ln(chi_m)``
        - ``lnvi``   = ``ln(Lambda)`` (their sum)

        Values are ``nan`` where an input is non-positive or non-finite (the
        decomposition is only defined for strictly positive terms).

    Examples
    --------
        >>> import tcpyPI
        >>> d = tcpyPI.vi_log_decomposition(10.0, 80.0, 0.4)
        >>> round(d["lnvi"] - (d["shear"] + d["pi"] + d["chi"]), 12)
        0.0
    """
    us, vp, ch = np.broadcast_arrays(
        np.asarray(u_shear, dtype=float),
        np.asarray(v_pot, dtype=float),
        np.asarray(chi_m, dtype=float),
    )
    bad = (us <= 0) | (vp <= 0) | (ch <= 0) | ~np.isfinite(us) | ~np.isfinite(vp) | ~np.isfinite(ch)
    with np.errstate(divide="ignore", invalid="ignore"):
        shear = np.where(bad, np.nan, np.log(us))
        pi_t = np.where(bad, np.nan, -np.log(vp))
        chi_t = np.where(bad, np.nan, np.log(ch))
    lnvi = shear + pi_t + chi_t

    def _s(a):
        return float(a) if a.ndim == 0 else a

    return {"lnvi": _s(lnvi), "shear": _s(shear), "pi": _s(pi_t), "chi": _s(chi_t)}
