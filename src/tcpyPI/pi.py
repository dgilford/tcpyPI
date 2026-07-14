"""pyPI: Potential Intensity Calculations in Python

Adapted from pcmin.m by Kerry Emanuel (ftp://texmex.mit.edu/pub/emanuel/TCMAX)
Originally Updated by Daniel Gilford for
Gilford et al. (2017) -- https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-16-0827.1
Gilford et al. (2019) -- https://journals.ametsoc.org/doi/10.1175/MWR-D-19-0021.1

Adapted for Python (pyPI/tcpypi) by Daniel Gilford, PhD (Climate Central, dgilford@climatecentral.org)
Full pyPI documentation, module validation, and sample code provided at:
********************** https://github.com/dgilford/tcpyPI ***************************

Download with the python Package index from the command line with:
   > pip install tcpypi

Last updated 7/11/2026, v1.4

Revision History:
  Revised on 9/24/2005 by K. Emanuel to fix convergence problems at high pressure
    --Converted to MATLAB  5/17/2008
  Revised 7/20/15 by D. Gilford to output the LNB
  Revised 8/4/16 by D. Gilford to include lack of convergence if SST < 5C for TO/LNB
  Revised 8/5/16 by D. Gilford to fix the "cape()" function output and include LNB
  Revised 10/3/16 by D. Gilford to set LNB to the pressure-weighted crossing of buoyancy from negative to positive (the zero-line)
  Revised 1/7/2017 by Luke Davis to fix and improper definition of JMIN, which should be the lowest profile level at or above the parcel level in the calculation of CAPE. Kerry is grateful to Luke and to Tim Merlis for pointing out this error. In practice, the error was usually less than 1 hPa in central pressure.
    --Converted to Python  04/2020
  Revised 4/10/2020 by D. Rothenberg (daniel@danielrothenberg.com) for Numba optimization
  Revised 4/13/2020 by D. Gilford to add new handling of missing profile data
  Revised 6/17/2020 by D. Gilford for auxilary files
  Revised 8/5/2020 by D. Gilford for auxilary files
  Revised 8/14/2020 by D. Gilford for python packaging
  Revised 10/15/2020 by D. Gilford to add missing SST-->IFL=0 condition/check
  Revised 2/1/2021 by D. Gilford to validate units of input SSTs/T profile (should be Celsius)
  Modernized 2/20/2025 by B. Mares and D. Gilford
  Revised 10/7/2025 by D. Gilford updating non-convergence flag==>2
  Revised 1/13/2026 by D. Gilford:
    - Added `outflow_source` option ("cape_star" vs "cape_env") to control how outflow temperature/level are defined.
    - Added log decomposition API (`log_decompose_pi`, `pi_log_decomposition`) and updated sample workflow/docs to use it.
    - Standardized main-module docstrings (NumPy style) and expanded README guidance on sensitivity/configuration options.
  Revised 7/11/2026 by D. Gilford (pi()/cape() input handling + log-decomposition consolidation):
    - `pi()` now sorts the input sounding to decreasing pressure, so profiles supplied in any vertical order are handled correctly.
    - `pi()` returns missing (IFL=3) when the pressure array or MSL contains NaNs, instead of an unconverged value incorrectly flagged IFL=1.
    - `cape()` now retains all levels with pressure > ptop (previously the level nearest ptop could be dropped when no level fell exactly at ptop).
    - Under miss_handle=0, the parcel is now lifted from the lowest valid level (previously a NaN in the lowest level(s) produced a NaN parcel and unconverged output).
    - Moved the PI log-decomposition into this module as `pi_log_decomposition` (from the former analyses.py) and removed the redundant run-and-decompose wrapper, matching the `vi_log_decomposition`/`gpi_log_decomposition` APIs.
    - Initialized the output flag before the environmental-CAPE call (as in pcmin.m) so a flag tripped by that call is no longer overwritten (GitHub issue #78).
    - Added a period-2 oscillation rescue to the minimum-pressure iteration: marginal near-neutral soundings whose fixed-point map locks into a 2-cycle (previously non-convergent, IFL=2) are collapsed to the cycle midpoint and converge (IFL=1). Previously-converging profiles are unchanged bit-for-bit.
  Revised 7/13/2026 by B. Mares (CAPE as maximum buoyancy work; GitHub issue #77):
    - `cape()` now defines CAPE as the maximum of the running signed buoyancy-work
      integral over candidate outflow levels (max-W), and the LNB as the level
      attaining it, instead of pcmin.m's integrate-to-the-highest-positive-level
      rule. The legacy rule is a discontinuous function of the inputs (a marginal
      buoyancy sign flicker aloft toggles the intervening negative area, jumping
      CAPE by O(100 J/kg)); the maximum is continuous, coincides with the legacy
      value on single-crossing profiles, and guarantees the pressure iteration a
      fixed point. Deliberately NOT bit-identical to pcmin.m: results change
      materially only on marginal multi-layer profiles (1.9% of the 2024 sample,
      no flag changes). Full analysis and validation: discontinuity_analysis/.
    - Removed the period-2 oscillation rescue: with a continuous pressure map the
      2-cycles it patched over no longer occur (both historical failing soundings
      now converge in a handful of iterations through the normal path).
"""
# doctest: +ELLIPSIS

# import required packages
import numpy as np

from . import constants, utilities
from .numba import guvectorize, njit


# define the function to calculate CAPE
@njit()
def cape(TP, RP, PP, T, R, P, ascent_flag=0, ptop=50, miss_handle=1):
    """Calculate the CAPE of a parcel given parcel and environmental conditions.

    function [CAPED,TOB,LNB,IFLAG]= cape(TP,RP,PP,T,R,P,ascent_flag=0,ptop=50,miss_handle=1)

    This function calculates the CAPE of a parcel given parcel pressure PP (hPa),
    temperature TP (K) and mixing ratio RP (gram/gram) and given a sounding
    of temperature (T in K) and mixing ratio (R in gram/gram) as a function
    of pressure (P in hPa). CAPED is the calculated value of CAPE following
    Emanuel 1994 (E94) Equation 6.3.6 and TOB is the temperature at the
    level of neutral buoyancy ("LNB") for the displaced parcel.

    Parameters
    ----------
    TP : float
        Parcel temperature (K).
    RP : float
        Parcel mixing ratio (gram/gram).
    PP : float
        Parcel pressure (hPa).
    T : array
        Environmental temperature profile (K).
    R : array
        Environmental mixing ratio profile (gram/gram).
    P : array
        Environmental pressure profile (hPa).

        The arrays MUST be arranged so that the lowest index corresponds to the
        lowest model level, with increasing index corresponding to decreasing
        pressure.
    ascent_flag : int, default=0
        Adjustable constant fraction for buoyancy of displaced parcels.

        - ``0``: reversible ascent
        - ``1``: pseudo-adiabatic ascent
    ptop : float, default=50
        Pressure below which sounding is ignored (hPa).
    miss_handle : int, default=1
        Flag for handling missing values.

        - ``0``: ignore NaN (BE02 default); NaNs in the profile are ignored and
          CAPE may still be calculated.
        - ``1``: return missing values (pyPI default); any NaNs set outputs to
          missing with ``IFLAG=3``.

        NOTE: If any missing values are between the lowest valid level and `ptop`,
        outputs are set to missing with ``IFLAG=3`` regardless of `miss_handle`.

    Returns
    -------
    tuple
        ``(CAPED, TOB, LNB, IFLAG)`` where:

        - CAPED : float
            Convective available potential energy (J/kg): the maximum of the
            running signed buoyancy-work integral over candidate outflow
            levels (grid levels and interpolated buoyancy zero-crossings).
            Non-negative by construction (the launch level, with zero work,
            is always a candidate).
        - TOB : float
            Temperature at level of neutral buoyancy (K).
        - LNB : float
            Level of neutral buoyancy pressure (hPa): the level attaining the
            work maximum (an interpolated buoyancy zero-crossing, or a grid
            level at the retained-profile top). ``0`` when no positive work
            is available.
        - IFLAG : int
            Status flag:

            - ``1``: success
            - ``0``: improper sounding/parcel
            - ``2``: did not converge
            - ``3``: missing values in input profile
    """
    #
    #   ***  Handle missing values   ***
    #

    # find if any values are missing in the temperature or mixing ratio array
    valid_i = ~np.isnan(T)
    first_valid = np.where(valid_i)[0][0]
    # Are there missing values? If so, assess according to flag
    if np.sum(valid_i) != len(P):
        # if not allowed, set IFLAG=3 and return missing CAPE
        if miss_handle != 0:
            CAPED = np.nan
            TOB = np.nan
            LNB = np.nan
            IFLAG = 3
            # Return the unsuitable values
            return (CAPED, TOB, LNB, IFLAG)
        else:
            # if allowed, but there are missing values between the lowest existing level
            # and ptop, then set IFLAG=3 and return missing CAPE
            if np.sum(np.isnan(T[first_valid : len(P)]) > 0):
                CAPED = np.nan
                TOB = np.nan
                LNB = np.nan
                IFLAG = 3
                # Return the unsuitable values
                return (CAPED, TOB, LNB, IFLAG)
            else:
                first_lvl = first_valid
    else:
        first_lvl = 0

    # Populate new environmental profiles removing values above ptop and
    # find new number, N, of profile levels with which to calculate CAPE.
    # Keep every level with P > ptop (matching pcmin); with the profile ordered
    # by decreasing pressure these are the leading N entries. Using argmin(|P-ptop|)
    # would drop the topmost retained level whenever no level sits exactly at ptop.
    N = int(np.count_nonzero(P > ptop))

    P = P[first_lvl:N]
    T = T[first_lvl:N]
    R = R[first_lvl:N]
    nlvl = len(P)
    TVRDIF = np.zeros((nlvl,))

    #
    #   ***  Run checks   ***
    #

    # CHECK: is the input profile ordered with increasing pressure? If not, return missing CAPE
    if P[2] - P[1] > 0:
        CAPED = 0
        TOB = np.nan
        LNB = np.nan
        IFLAG = 0
        # Return the unsuitable values
        return (CAPED, TOB, LNB, IFLAG)

    # CHECK: Is the input parcel suitable? If not, return missing CAPE
    if (RP < 1e-6) or (TP < 200):
        CAPED = 0
        TOB = np.nan
        LNB = np.nan
        IFLAG = 0
        # Return the unsuitable values
        return (CAPED, TOB, LNB, IFLAG)

    #
    #  ***  Define various parcel quantities, including reversible   ***
    #  ***                       entropy, S                          ***
    #
    TPC = utilities.T_ktoC(TP)  # Parcel temperature in Celsius
    ESP = utilities.es_cc(TPC)  # Parcel's saturated vapor pressure
    EVP = utilities.ev(RP, PP)  # Parcel's partial vapor pressure
    RH = EVP / ESP  # Parcel's relative humidity
    RH = min(RH, 1.0)  # ensure that the relatively humidity does not exceed 1.0
    # calculate reversible total specific entropy per unit mass of dry air (E94, EQN. 4.5.9)
    S = utilities.entropy_S(TP, RP, PP)

    #
    #   ***  Estimate lifted condensation level pressure, PLCL   ***
    #     Based on E94 "calcsound.f" code at http://texmex.mit.edu/pub/emanuel/BOOK/
    #     see also https://psl.noaa.gov/data/composites/day/calculation.html
    #
    #   NOTE: Modern PLCL calculations are made following the exact expressions of Romps (2017),
    #   see https://journals.ametsoc.org/doi/pdf/10.1175/JAS-D-17-0102.1
    #   and Python PLCL code at http://romps.berkeley.edu/papers/pubdata/2016/lcl/lcl.py
    #
    PLCL = utilities.e_pLCL(TP, RH, PP)

    # Initial default values before loop
    CAPED = 0
    TOB = T[0]
    IFLAG = 1
    # Values to help loop
    jmin = int(1e6)

    #
    #   ***  Begin updraft loop   ***
    #

    # loop over each level in the profile
    for j in range(nlvl):
        # jmin is the index of the lowest pressure level evaluated in the loop
        jmin = int(min(jmin, j))

        #
        #   *** Calculate Parcel quantities BELOW lifted condensation level   ***
        #
        if P[j] >= PLCL:
            # Parcel temperature at this pressure
            TG = TP * (P[j] / PP) ** (constants.RD / constants.CPD)
            # Parcel Mixing ratio
            RG = RP
            # Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
            TLVR = utilities.Trho(TG, RG, RG)
            TVENV = utilities.Trho(T[j], R[j], R[j])
            # Bouyancy of the parcel in the environment (Proxy of E94, EQN. 6.1.5)
            TVRDIF[j,] = TLVR - TVENV

        #
        #   *** Calculate Parcel quantities ABOVE lifted condensation level   ***
        #
        else:
            TG, RG, IFLAG = solve_temperature_from_entropy(S=S, P=P[j], RP=RP, T_initial=T[j])
            if IFLAG == 2:  # Did not converge
                CAPED = 0
                TOB = T[0]
                LNB = P[0]
                # Return the uncoverged values
                return (CAPED, TOB, LNB, IFLAG)

            #
            #   *** Calculate buoyancy   ***
            #
            # Parcel total mixing ratio: either reversible (ascent_flag=0) or pseudo-adiabatic (ascent_flag=1)
            RMEAN = ascent_flag * RG + (1 - ascent_flag) * RP
            # Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
            TLVR = utilities.Trho(TG, RMEAN, RG)
            TENV = utilities.Trho(T[j], R[j], R[j])
            # Bouyancy of the parcel in the environment (Proxy of E94, EQN. 6.1.5)
            TVRDIF[j,] = TLVR - TENV

    #
    #  ***  CAPE as maximum buoyancy work (max-W)  ***
    #
    # Accumulate the running signed work integral W(p_t) of the
    # piecewise-linear buoyancy profile (the same trapezoid increments and
    # interpolated zero-crossing partial areas as pcmin.m), and take CAPE as
    # its maximum over candidate outflow levels p_t: every grid level above
    # the parcel, plus every interpolated sign crossing of the buoyancy.
    #
    # This replaces pcmin.m's terminal-point rule (integrate to the HIGHEST
    # grid level with positive buoyancy, clamp the signed result at zero),
    # which is a discontinuous function of the inputs: an arbitrarily small
    # buoyancy sign change at a marginal level aloft drags the intervening
    # negative area in or out of the integral all at once, producing finite
    # CAPE jumps, a discontinuous pressure map, and the period-2 solver
    # failures of GitHub issue #77. The maximum is continuous in the inputs
    # (Berge's maximum theorem) even when the maximizing level changes, and
    # is identical to the pcmin.m value for single-crossing buoyancy
    # profiles -- the overwhelmingly common case. See
    # discontinuity_analysis/ for the full derivation and validation.
    #
    # The launch level (zero accumulated work) is always a candidate, so
    # CAPED >= 0 with no final clamp.

    # signed work from the parcel pressure PP to the first retained level
    W = constants.RD * (PP - P[jmin]) / (PP + P[jmin]) * TVRDIF[jmin]
    BESTW = 0.0
    BESTP = 0.0
    BESTT = T[0]
    if W > BESTW:
        BESTW = W
        BESTP = P[jmin]
        BESTT = T[jmin]
    for j in range(jmin + 1, nlvl, 1):
        B0 = TVRDIF[j - 1]
        B1 = TVRDIF[j]
        if B0 * B1 < 0.0:
            # candidate: interpolated buoyancy zero-crossing between the
            # levels (the same PINB/partial-area formulas as pcmin.m)
            PC = (P[j] * B0 - P[j - 1] * B1) / (B0 - B1)
            WC = W + constants.RD * B0 * (P[j - 1] - PC) / (P[j - 1] + PC)
            if WC > BESTW:
                BESTW = WC
                BESTP = PC
                BESTT = (T[j - 1] * (PC - P[j]) + T[j] * (P[j - 1] - PC)) / (P[j - 1] - P[j])
        # candidate: the grid level itself (relevant at the retained-profile
        # top, where buoyancy may still be positive)
        W = W + constants.RD * (B1 + B0) * (P[j - 1] - P[j]) / (P[j] + P[j - 1])
        if W > BESTW:
            BESTW = W
            BESTP = P[j]
            BESTT = T[j]

    # CHECK: Is any positive work available? If not, return zero CAPE
    # (mirrors the historical INB==0 return: LNB=0, TOB=T[0])
    if BESTW <= 0.0:
        CAPED = 0
        TOB = T[0]
        LNB = 0
        return (CAPED, TOB, LNB, IFLAG)

    CAPED = BESTW
    TOB = BESTT
    LNB = BESTP
    # set the flag to OK if procedure reached this point
    IFLAG = 1
    # Return the calculated outputs to the above program level
    return (CAPED, TOB, LNB, IFLAG)


@njit()
def solve_temperature_from_entropy(S, P, RP, T_initial):
    """Compute the temperature corresponding to a given value for saturated entropy.

    For given pressure P and dry mixing ratio RP, saturated entropy is a function
    of temperature SG(T). This function uses Newton-Raphson iteration to invert
    saturated entropy to find the temperature T that satisfies SG(T) = S for
    a target value of saturated entropy S.

    Parameters
    ----------
    S : float
        Target saturated entropy (J/kg/K).
    P : float
        Ambient pressure (hPa).
    RP : float
        Parcel mixing ratio for the dry component (g/g).
    T_initial : float
        Initial guess for temperature (K).

    Returns
    -------
    tuple
        ``(TG, RG, IFLAG)`` where:

        - TG : float
            Temperature (K) satisfying ``SG(T) = S``.
        - RG : float
            Computed saturated mixing ratio (g/g).
        - IFLAG : int
            Status flag:

            - ``1``: success
            - ``2``: did not converge

    Examples
    --------
        >>> S = 4000
        >>> P = 1000
        >>> RP = 0.01
        >>> TG, RG, IFLAG = solve_temperature_from_entropy(S, P, RP, T_initial=300)
        >>> print(f"TG: {TG}, RG: {RG}, IFLAG: {IFLAG}")
        TG: 292.676..., RG: 0.0144..., IFLAG: 1
    """
    # Initial default values before loop
    TGNEW = T_initial
    TJC = utilities.T_ktoC(T_initial)
    ES = utilities.es_cc(TJC)
    RG = utilities.rv(ES, P)

    #
    #   ***  Iteratively calculate lifted parcel temperature and mixing   ***
    #   ***                ratio for reversible ascent                    ***
    #

    # set loop counter and initial condition
    NC = 0
    TG = 0

    # loop until loop converges or bails out
    while (np.abs(TGNEW - TG)) > 0.001:
        # Parcel temperature and mixing ratio during this iteration
        TG = TGNEW
        TC = utilities.T_ktoC(TG)
        ENEW = utilities.es_cc(TC)
        RG = utilities.rv(ENEW, P)

        # increase iteration count in the loop
        NC += 1

        #
        #   ***  Calculate estimates of the rates of change of the entropy    ***
        #   ***           with temperature at constant pressure               ***
        #

        ALV = utilities.Lv(TC)
        # calculate the rate of change of entropy with temperature, s_ell
        SL = (constants.CPD + RP * constants.CL + ALV * ALV * RG / (constants.RV * TG * TG)) / TG
        EM = utilities.ev(RG, P)
        # calculate the saturated entropy, s_k, noting r_T=RP and
        # the last term vanishes with saturation, i.e. RH=1
        SG = (
            (constants.CPD + RP * constants.CL) * np.log(TG)
            - constants.RD * np.log(P - EM)
            + ALV * RG / TG
        )
        # convergence speed (AP, step in entropy fraction) varies as a function of
        # number of iterations
        if NC < 3:
            # converge slowly with a smaller step
            AP = 0.3
        else:
            # speed the process with a larger step when nearing convergence
            AP = 1.0
        # find the new temperature in the iteration
        TGNEW = TG + AP * (S - SG) / SL

        #
        #   ***   If the routine does not converge, set IFLAG=2 and bail out   ***
        #
        if (NC > 500) or (ENEW > (P - 1)):
            IFLAG = 2  # Did not converge
            return TG, RG, IFLAG

    IFLAG = 1  # Success
    return TG, RG, IFLAG


# define the function to calculate PI
@njit()
def _pi_numba(
    SSTC,
    MSL,
    P,
    TC,
    R,
    CKCD=0.9,
    ascent_flag=0,
    diss_flag=1,
    V_reduc=0.8,
    ptop=50,
    miss_handle=1,
    outflow_source_flag=0,
):
    """Internal Numba-compiled implementation for :func:`tcpyPI.pi`.

    Notes
    -----
    `outflow_source_flag` selects how outflow (`TO`, `OTL`) is defined:
    - 0: ``"cape_star"`` (default)
    - 1: ``"cape_env"``
    """

    # convert units
    SSTK = utilities.T_Ctok(SSTC)  # SST in kelvin
    T = utilities.T_Ctok(TC)  # Temperature profile in kelvin
    R = R * 0.001  # Mixing ratio profile in g/g

    # CHECK 0: is MSL missing? If so, set IFL=3 and return missing PI.
    # (Mirrors the pi() wrapper guard with the same precedence; needed for the
    # whole-field pi_field path, where per-cell MSL NaNs reach this compiled core
    # directly. Untaken for valid MSL, so existing results are unchanged.)
    if np.isnan(MSL):
        VMAX = np.nan
        PMIN = np.nan
        IFL = 3
        TO = np.nan
        OTL = np.nan
        return (VMAX, PMIN, IFL, TO, OTL)

    # CHECK 1a: do SSTs exceed 5C?
    # CHECK 1b: are SSTs less than 100C (if not, indicative of input in kelvin)
    # If not, set IFL=0 and return missing PI
    if SSTC <= 5.0 or SSTC > 100:
        VMAX = np.nan
        PMIN = np.nan
        IFL = 0
        TO = np.nan
        OTL = np.nan
        return (VMAX, PMIN, IFL, TO, OTL)
    # CHECK 1c: are SSTs missing? If so, set IFL=0 and return missing PI
    elif np.isnan(SSTC):
        VMAX = np.nan
        PMIN = np.nan
        IFL = 0
        TO = np.nan
        OTL = np.nan
        return (VMAX, PMIN, IFL, TO, OTL)

    # CHECK 2a: do Temperature profiles exceed 100K?
    # CHECK 2b: are Temperatures in Celsius less than 100C (if not, indicative of input in kelvin)
    # If not, set IFL=0 and return missing PI
    if np.min(T) <= 100 or np.max(TC) > 100:
        VMAX = np.nan
        PMIN = np.nan
        IFL = 0
        TO = np.nan
        OTL = np.nan
        return (VMAX, PMIN, IFL, TO, OTL)

    # Set Missing mixing ratios to zero g/g, following Kerry's BE02 algorithm
    R[np.isnan(R)] = 0.0

    # Saturated water vapor pressure
    # from Clausius-Clapeyron relation/August-Roche-Magnus formula
    ES0 = utilities.es_cc(SSTC)

    # define the level from which parcels are lifted. Normally the lowest level
    # (index 0). Under miss_handle=0 (ignore-NaN), skip any leading missing
    # temperature levels (e.g. sub-surface levels over topography) so the parcel is
    # lifted from the lowest *valid* level rather than a NaN.
    NK = 0
    if miss_handle == 0:
        valid_T = ~np.isnan(T)
        if np.any(valid_T):
            NK = int(np.where(valid_T)[0][0])

    #
    #   ***   Find environmental CAPE ***
    #
    # Default flag for the CAPE calculations. Initialized BEFORE the environmental
    # CAPE call (as in pcmin.m) so a flag tripped by CAPEA is not clobbered by a
    # later reset (see GitHub issue #78).
    IFL = 1
    TP = T[NK]
    RP = R[NK]
    PP = P[NK]
    result = cape(TP, RP, PP, T, R, P, ascent_flag, ptop, miss_handle)
    CAPEA = result[0]
    IFLAG = result[3]
    # if the CAPE function tripped a flag, set the output IFL to it
    if IFLAG != 1:
        IFL = int(IFLAG)

    #
    #   ***   Begin iteration to find mimimum pressure   ***
    #

    # set loop counter and initial condition
    NP = 0  # loop counter
    PM = 970.0
    PMOLD = PM  # initial condition from minimum pressure
    PNEW = 0.0  # initial condition from minimum pressure

    TO_ENV = np.nan
    OTL_ENV = np.nan

    # loop until convergence or bail out
    while np.abs(PNEW - PMOLD) > 0.5:
        #
        #   ***  Find CAPE at radius of maximum winds   ***
        #
        TP = T[NK]
        PP = min(PM, 1000.0)
        # find the mixing ratio with the average of the lowest level pressure and MSL
        RP = constants.EPS * R[NK] * MSL / (PP * (constants.EPS + R[NK]) - R[NK] * MSL)
        CAPEM, TOB_ENV, LNB_ENV, IFLAG = cape(TP, RP, PP, T, R, P, ascent_flag, ptop, miss_handle)
        # if the CAPE function tripped a different flag, set the output IFL to it
        if IFLAG != 1:
            IFL = int(IFLAG)
        TO_ENV = TOB_ENV
        OTL_ENV = LNB_ENV

        #
        #  ***  Find saturation CAPE at radius of maximum winds    ***
        #  *** Note that TO and OTL are found with this assumption ***
        #
        TP = SSTK
        PP = min(PM, 1000.0)
        RP = utilities.rv(ES0, PP)
        result = cape(TP, RP, PP, T, R, P, ascent_flag, ptop, miss_handle)
        CAPEMS, TOMS, LNBS, IFLAG = result
        # if the CAPE function tripped a flag, set the output IFL to it
        if IFLAG != 1:
            IFL = int(IFLAG)
        # Store the outflow temperature and level of neutral bouyancy at the outflow level (OTL)
        TO = TOMS
        OTL = LNBS
        # Optionally use CAPEenv (final-iteration) outflow definition.
        if outflow_source_flag == 1:
            TO = TO_ENV
            OTL = OTL_ENV
        # Calculate the proxy for TC efficiency (BE02, EQN. 1-3)
        RAT = SSTK / TO
        # If dissipative heating is "off", TC efficiency proxy is set to 1.0 (BE02, pg. 3)
        if diss_flag == 0:
            RAT = 1.0

        #
        #  ***  Initial estimate of pressure at the radius of maximum winds  ***
        #
        RS0 = RP
        # Lowest level and Sea-surface Density Temperature (E94, EQN. 4.3.1 and 6.3.7)
        TV0 = utilities.Trho(T[NK], R[NK], R[NK])
        TVSST = utilities.Trho(SSTK, RS0, RS0)
        # Average Surface Density Temperature, e.g. 1/2*[Tv(Tsfc)+Tv(sst)]
        TVAV = 0.5 * (TV0 + TVSST)
        # Converge toward CAPE*-CAPEM (BE02, EQN 3-4)
        CAT = (CAPEM - CAPEA) + 0.5 * CKCD * RAT * (CAPEMS - CAPEM)
        CAT = max(CAT, 0.0)
        # Iterate on pressure
        PNEW = MSL * np.exp(-CAT / (constants.RD * TVAV))

        #
        #   ***  Test for convergence (setup for possible next while iteration)  ***
        #
        # store the previous step's pressure
        PMOLD = PM
        # store the current step's pressure
        PM = PNEW

        # increase iteration count in the loop
        NP += 1
        #
        #   ***   If the routine does not converge, set IFL=2 and return missing PI   ***
        #
        if (NP > 200) or (PM < 400):
            VMAX = np.nan
            PMIN = np.nan
            IFL = 2
            TO = np.nan
            OTL = np.nan
            return (VMAX, PMIN, IFL, TO, OTL)

    # Once converged, set potential intensity at the radius of maximum winds
    CATFAC = 0.5 * (1.0 + 1 / constants.b)
    CAT = (CAPEM - CAPEA) + CKCD * RAT * CATFAC * (CAPEMS - CAPEM)
    CAT = max(CAT, 0.0)

    # Calculate the minimum pressure at the eye of the storm
    # BE02 EQN. 4
    PMIN = MSL * np.exp(-CAT / (constants.RD * TVAV))

    # Calculate the potential intensity at the radius of maximum winds
    # BE02 EQN. 3, reduced by some fraction (default 20%) to account for the reduction
    # of 10-m winds from gradient wind speeds (Emanuel 2000, Powell 1980)
    FAC = max(0.0, (CAPEMS - CAPEM))
    VMAX = V_reduc * np.sqrt(CKCD * RAT * FAC)

    # Return the calculated outputs to the above program level
    return (VMAX, PMIN, IFL, TO, OTL)


def _outflow_source_flag(outflow_source):
    """Map the ``outflow_source`` keyword to the compiled core's integer flag.

    Shared by :func:`pi` and :func:`pi_field` so the accepted values and error
    message cannot drift between the two entry points.
    """
    if outflow_source == "cape_star":
        return 0
    if outflow_source == "cape_env":
        return 1
    raise ValueError(
        f"Invalid outflow_source={outflow_source!r}; expected 'cape_star' or 'cape_env'."
    )


def _sort_profile_descending(P, TC, R):
    """Sort profile arrays to decreasing pressure along the trailing level axis.

    Shared by :func:`pi` (1-D profiles) and :func:`pi_field` (N-D fields with the
    level dimension last); already-descending input is returned untouched.
    """
    order = np.argsort(P)[::-1]
    if np.array_equal(order, np.arange(P.size)):
        return P, TC, R
    return P[order], np.take(TC, order, axis=-1), np.take(R, order, axis=-1)


def pi(
    SSTC,
    MSL,
    P,
    TC,
    R,
    CKCD=0.9,
    ascent_flag=0,
    diss_flag=1,
    V_reduc=0.8,
    ptop=50,
    miss_handle=1,
    outflow_source="cape_star",
):
    r"""Calculate maximum potential intensity given environmental conditions.

    function [VMAX,PMIN,IFL,TO,OTL] = pi(SSTC,MSL,P,TC,R,CKCD=0.9,ascent_flag=0,diss_flag=1,V_reduc=0.8,ptop=50,miss_handle=0,outflow_source="cape_star")

    This function calculates the maximum wind speed and minimum central pressure
    achievable in tropical cyclones, given a sounding and a sea surface temperature.

    Thermodynamic and dynamic technical backgrounds (and calculations) are found in
    Bister and Emanuel (2002; BE02) and Emanuel's "Atmospheric Convection" (E94; 1994; ISBN: 978-0195066302).

    Parameters
    ----------
    SSTC : float
        Sea surface temperature (C).
    MSL : float
        Mean sea level pressure (hPa).
    P : array
        Pressure levels (hPa).
    TC : array
        Temperature profile (C).
    R : array
        Mixing ratio profile (g/kg).

        The arrays MUST be arranged so that the lowest index corresponds to the
        lowest model level, with increasing index corresponding to decreasing
        pressure. The temperature sounding should extend to at least the tropopause
        and preferably to the lower stratosphere; mixing ratios are not important
        above the boundary layer. Missing mixing ratios can be replaced by zeros.
    CKCD : float, default=0.9
        Ratio of C_k to C_D (unitless).

        The ratio of the exchange coefficients of enthalpy and momentum flux
        (e.g. see Bister and Emanuel 1998, EQN. 17-18). More discussion on Ck/Cd is
        found in Emanuel (2003). Default is 0.9 based on e.g. Wing et al. (2015).
    ascent_flag : int, default=0
        Adjustable constant fraction for buoyancy (unitless).

        - ``0``: reversible ascent
        - ``1``: pseudo-adiabatic ascent
    diss_flag : int, default=1
        Switch for dissipative heating.

        - ``1``: allowed
        - ``0``: disallowed

        See Bister and Emanuel (1998) for inclusion of dissipative heating.
    V_reduc : float, default=0.8
        Reduction factor from gradient winds to 10 m winds (unitless).
        See Emanuel (2000) and Powell (1980).
    ptop : float, default=50
        Pressure below which sounding is ignored (hPa).
    miss_handle : int, default=1
        Flag for handling missing values.

        - ``0``: ignore NaN (BE02 default); NaN values in the profile are ignored
          and PI may still be calculated.
        - ``1``: return missing values (pyPI default); any NaNs set outputs to
          missing with ``IFL=3``.

        NOTE: If any missing values are between the lowest valid level and `ptop`,
        outputs are set to missing with ``IFL=3`` regardless of `miss_handle`.
    outflow_source : {"cape_star", "cape_env"}, default="cape_star"
        Which CAPE calculation supplies the outflow temperature and level:

        - ``"cape_star"``: use the saturated CAPE* calculation (parcel saturated at
          SST) (default; matches BE02/pcmin behavior).
        - ``"cape_env"``: use the environmental CAPE calculation (CAPEenv) from the
          final convergence iteration (as discussed in Gilford et al. 2021).

    Returns
    -------
    tuple
        ``(VMAX, PMIN, IFL, TO, OTL)`` where:

        - VMAX : float
            Maximum surface wind speed (m/s), reduced by `V_reduc`.
        - PMIN : float
            Minimum central pressure (hPa).
        - IFL : int
            Status flag:

            - ``1``: success
            - ``0``: invalid input (e.g. SST <= 5 C, Kelvin-like units, improper sounding)
            - ``2``: routine did not converge (pressure iteration or CAPE solver)
            - ``3``: missing data in the input profile

            Notes on flag semantics:

            - With CAPE defined as the maximum buoyancy work (see `cape`), the
              pressure iteration's map is continuous and marginal near-neutral
              soundings that previously locked into period-2 oscillations
              (returning ``IFL=2``; GitHub issue #77) now converge normally
              (``IFL=1``). This is a documented, deliberate improvement over
              pcmin.m, which returns missing for such soundings.
            - A flag tripped by the *environmental*-CAPE calculation persists
              (see GitHub issue #78), so ``IFL=0`` or ``IFL=3`` can accompany
              finite ``VMAX``/``PMIN`` on an otherwise-converged column (e.g. a
              degenerate near-dry surface parcel). Treat ``IFL==1`` as the
              "fully trustworthy output" gate. (pcmin.m collapses these cases to
              its generic non-convergence flag instead.)
        - TO : float
            Outflow temperature (K).
        - OTL : float
            Outflow temperature level (hPa): the level of neutral buoyancy of
            the saturated parcel, defined as the level attaining the maximum
            buoyancy work (see `cape`).

    Examples
    --------
        >>> SSTC = 30
        >>> MSL = 1010
        >>> level_data = np.array(
        ...     [
        ...         [1000, 28, 18],
        ...         [975, 25, 18],
        ...         [950, 24, 16],
        ...         [925, 23, 13],
        ...         [900, 22, 12],
        ...         [875, 20, 11],
        ...         [850, 19, 10],
        ...         [825, 18, 10],
        ...         [800, 16, 9],
        ...         [775, 15, 8],
        ...         [750, 13, 7],
        ...         [700, 11, 4],
        ...         [650, 8, 3],
        ...         [600, 5, 1.7],
        ...         [550, 2, 1.2],
        ...         [500, -2, 1.7],
        ...         [450, -6, 0.7],
        ...         [400, -11, 0.2],
        ...         [350, -18, 0.15],
        ...         [300, -27, 0.10],
        ...         [250, -37, 0.11],
        ...         [225, -43, 0.08],
        ...         [200, -49, 0.05],
        ...         [175, -57, 0.03],
        ...         [150, -65, 0.014],
        ...         [125, -73, 0.005],
        ...         [100, -79, 0.003],
        ...         [70, -73, 0.002],
        ...         [50, -64, 0.002],
        ...     ]
        ... )
        >>> P = level_data[:, 0]
        >>> TC = level_data[:, 1]
        >>> R = level_data[:, 2]
        >>> VMAX, PMIN, IFL, TO, OTL = pi(SSTC, MSL, P, TC, R)
        >>> print(f"VMAX: {VMAX}\nPMIN: {PMIN}\nIFL: {IFL}\nTO: {TO}\nOTL: {OTL}")
        VMAX: 82.4845...
        PMIN: 900.2039...
        IFL: 1
        TO: 197.1666...
        OTL: 84.9169...

    Exceptional cases:
        - SST is missing, e.g. over land
            >>> pi(np.nan, MSL, P, TC, R)
            (nan, nan, 0, nan, nan)

        - SST is given in Kelvin
            >>> pi(300, MSL, P, TC, R)
            (nan, nan, 0, nan, nan)

        - Temperatures contain a nan
            >>> nan_in_levels = np.zeros_like(P)
            >>> nan_in_levels[5] = np.nan
            >>> pi(SSTC, MSL, P, TC + nan_in_levels, R)
            (nan, nan, 3, nan, nan)

        - Mixing ratios contain a nan
            >>> pi(SSTC, MSL, P, TC, R + nan_in_levels)
            (82.4845..., 900.2039..., 1, 197.1666..., 84.9169...)
    """
    outflow_source_flag = _outflow_source_flag(outflow_source)

    P = np.asarray(P)
    TC = np.asarray(TC)
    R = np.asarray(R)

    # Guard against NaNs in the shared pressure array, which poison the CAPE
    # integration (the loop can otherwise exit early and return numerically
    # plausible values with IFL=1). Fail fast to IFL=3 (missing data). NaN MSL is
    # handled by CHECK 0 inside the compiled core, the single owner of that policy
    # for both pi() and pi_field().
    if np.any(np.isnan(P)):
        return (np.nan, np.nan, 3, np.nan, np.nan)

    # Order-agnostic: the algorithm requires lowest index = lowest level
    # (decreasing pressure). Sort the profile here so any input ordering is
    # accepted; already-descending input is left untouched.
    P, TC, R = _sort_profile_descending(P, TC, R)

    return _pi_numba(
        SSTC,
        MSL,
        P,
        TC,
        R,
        CKCD=CKCD,
        ascent_flag=ascent_flag,
        diss_flag=diss_flag,
        V_reduc=V_reduc,
        ptop=ptop,
        miss_handle=miss_handle,
        outflow_source_flag=outflow_source_flag,
    )


if guvectorize is not None:
    # Whole-field generalized ufunc around the compiled per-column core. Processing
    # the entire field inside compiled, multithreaded code removes the per-column
    # Python-loop overhead of xarray.apply_ufunc(..., vectorize=True) (the speedup
    # highlighted by GitHub issue #86). Same per-column math as pi(): results are
    # bit-identical for identical (float64) inputs.
    @guvectorize(
        [
            "void(float64, float64, float64[:], float64[:], float64[:], "
            "float64, int64, int64, float64, float64, int64, int64, "
            "float64[:], float64[:], int64[:], float64[:], float64[:])"
        ],
        "(),(),(n),(n),(n),(),(),(),(),(),(),()->(),(),(),(),()",
        target="parallel",
        cache=True,
    )
    def _pi_field_gufunc(
        sstc,
        msl,
        p,
        tc,
        r,
        ckcd,
        ascent_flag,
        diss_flag,
        v_reduc,
        ptop,
        miss_handle,
        outflow_source_flag,
        vmax,
        pmin,
        ifl,
        to,
        otl,
    ):
        res = _pi_numba(
            sstc,
            msl,
            p,
            tc,
            r,
            ckcd,
            ascent_flag,
            diss_flag,
            v_reduc,
            ptop,
            miss_handle,
            outflow_source_flag,
        )
        vmax[0] = res[0]
        pmin[0] = res[1]
        ifl[0] = res[2]
        to[0] = res[3]
        otl[0] = res[4]
else:  # pragma: no cover - exercised only with TCPYPI_DISABLE_NUMBA=1
    _pi_field_gufunc = None


def pi_field(
    SSTC,
    MSL,
    P,
    TC,
    R,
    *,
    CKCD=0.9,
    ascent_flag=0,
    diss_flag=1,
    V_reduc=0.8,
    ptop=50,
    miss_handle=1,
    outflow_source="cape_star",
):
    r"""Calculate potential intensity over a whole gridded field at once.

    A fast, compiled, multithreaded path over many columns: the entire field is
    processed inside a Numba generalized ufunc (``target='parallel'``), removing
    the per-column Python-loop overhead of ``xarray.apply_ufunc(pi, ...,
    vectorize=True)``. The per-column math is identical to :func:`pi` (the same
    compiled core), so for identical float64 inputs the results are bit-identical.
    Motivated by the whole-field benchmark in GitHub issue #86 — implemented
    natively so only tcpyPI code executes.

    Parameters
    ----------
    SSTC : array-like
        Sea surface temperatures (C); any shape (e.g. ``(month, lat, lon)``).
        NaN cells (e.g. land) return missing with ``IFL=0``.
    MSL : array-like
        Mean sea level pressures (hPa), broadcastable against `SSTC`.
        NaN cells return missing with ``IFL=3``.
    P : array
        1-D pressure levels (hPa), shared by all columns. Any vertical order is
        accepted (sorted internally to decreasing pressure).
    TC, R : array-like
        Temperature (C) and mixing ratio (g/kg) profiles. The **level dimension
        must be the last (trailing) axis**; leading axes broadcast against `SSTC`.
        (With xarray, ``xr.apply_ufunc(pi_field, ..., input_core_dims=[[], [],
        ['p'], ['p'], ['p']], ...)`` arranges this automatically — no
        ``vectorize=True`` needed.)
    CKCD, ascent_flag, diss_flag, V_reduc, ptop, miss_handle, outflow_source
        As in :func:`pi`.

    Returns
    -------
    tuple of numpy.ndarray
        ``(VMAX, PMIN, IFL, TO, OTL)`` with the broadcast field shape; `IFL` is an
        integer array. See :func:`pi` for definitions, units, and flag meanings.

    Notes
    -----
    Inputs are converted to float64 (use :func:`pi` directly if you need another
    precision path). If Numba is disabled/unavailable, falls back to a per-column
    loop over :func:`pi` (slow but correct).

    Examples
    --------
        >>> field = pi_field(np.array([[30.0, np.nan]]), 1010.0,
        ...                  np.array([1000., 900, 800, 700, 600, 500, 400, 300, 200, 100, 50]),
        ...                  np.array([28., 22, 16, 11, 5, -2, -11, -27, -49, -79, -64]),
        ...                  np.array([18., 12, 9, 4, 1.7, 1.7, 0.2, 0.1, 0.05, 0.003, 0.002]))
        >>> field[2].tolist()  # IFL: ocean cell converges, NaN-SST cell is missing
        [[1, 0]]
    """
    outflow_source_flag = _outflow_source_flag(outflow_source)

    SSTC = np.asarray(SSTC, dtype=np.float64)
    MSL = np.asarray(MSL, dtype=np.float64)
    # ascontiguousarray for the profile arrays: xarray.apply_ufunc hands us
    # transposed (strided) views; a contiguous layout lets the parallel gufunc
    # stream each column from adjacent memory. Memory layout only — values (and
    # therefore results) are unchanged.
    P = np.ascontiguousarray(P, dtype=np.float64)
    TC = np.ascontiguousarray(TC, dtype=np.float64)
    R = np.ascontiguousarray(R, dtype=np.float64)

    if P.ndim != 1:
        raise ValueError(
            "pi_field expects a single 1-D pressure-level array shared by all "
            "columns; for per-column pressure grids call pi() per column."
        )
    if TC.shape[-1] != P.size or R.shape[-1] != P.size:
        raise ValueError(
            "The level dimension of TC and R must be the last axis and match P "
            f"(got P: {P.size}, TC: {TC.shape[-1]}, R: {R.shape[-1]})."
        )

    field_shape = np.broadcast_shapes(SSTC.shape, MSL.shape, TC.shape[:-1], R.shape[:-1])

    # Global guard (mirrors pi()): a NaN anywhere in the shared pressure array
    # poisons every column, so return an all-missing field with IFL=3.
    if np.any(np.isnan(P)):
        nanfield = np.full(field_shape, np.nan)
        return (
            nanfield,
            nanfield.copy(),
            np.full(field_shape, 3, dtype=np.int64),
            nanfield.copy(),
            nanfield.copy(),
        )

    # Order-agnostic: sort the shared profile to decreasing pressure once.
    P, TC, R = _sort_profile_descending(P, TC, R)

    if _pi_field_gufunc is None:
        # Pure-Python fallback (TCPYPI_DISABLE_NUMBA=1 or Numba unavailable).
        b_sst, b_msl = np.broadcast_arrays(
            np.broadcast_to(SSTC, field_shape), np.broadcast_to(MSL, field_shape)
        )
        b_tc = np.broadcast_to(TC, field_shape + (P.size,))
        b_r = np.broadcast_to(R, field_shape + (P.size,))
        vmax = np.full(field_shape, np.nan)
        pmin = np.full(field_shape, np.nan)
        ifl = np.zeros(field_shape, dtype=np.int64)
        to = np.full(field_shape, np.nan)
        otl = np.full(field_shape, np.nan)
        for idx in np.ndindex(field_shape):
            out = pi(
                b_sst[idx],
                b_msl[idx],
                P,
                b_tc[idx],
                b_r[idx],
                CKCD=CKCD,
                ascent_flag=ascent_flag,
                diss_flag=diss_flag,
                V_reduc=V_reduc,
                ptop=ptop,
                miss_handle=miss_handle,
                outflow_source=outflow_source,
            )
            vmax[idx], pmin[idx], ifl[idx], to[idx], otl[idx] = out
        return (vmax, pmin, ifl, to, otl)

    return _pi_field_gufunc(
        SSTC,
        MSL,
        P,
        TC,
        R,
        float(CKCD),
        int(ascent_flag),
        int(diss_flag),
        float(V_reduc),
        float(ptop),
        int(miss_handle),
        int(outflow_source_flag),
    )


def pi_log_decomposition(pi, sst, t0, CKCD=0.9, *, sst_units="K"):
    r"""Log-decompose potential intensity into efficiency, disequilibrium, and Ck/Cd.

    Following Wing et al. (2015, Eq. 2), PI separates additively in log space:

    .. math:: \ln(V^2) = \ln(C_k/C_D) + \ln\!\frac{T_s - T_0}{T_0} + \ln(\text{disequilibrium})

    This is the PI analogue of :func:`tcpyPI.vi_log_decomposition` and
    :func:`tcpyPI.gpi_log_decomposition`. Inputs may be scalars or arrays (NumPy
    broadcasting); invalid physical states yield ``nan``.

    Parameters
    ----------
    pi : float or array-like
        Potential intensity wind speed (m/s), e.g. the first output of :func:`pi`.
    sst : float or array-like
        Sea surface temperature, in units given by ``sst_units``.
    t0 : float or array-like
        Outflow temperature (K).
    CKCD : float, default=0.9
        Ratio of exchange coefficients (Ck/Cd).
    sst_units : {"K", "C"}, default="K"
        Units for ``sst``. :func:`pi` takes/returns SST-related quantities in
        Celsius, so pass ``sst_units="C"`` when decomposing its outputs directly.

    Returns
    -------
    tuple
        ``(lnpi, lneff, lndiseq, lnCKCD)`` where ``lnpi = ln(pi**2) = 2*ln(pi)``.
        A 4-tuple (rather than a dict like the VI/GPI decompositions) so it composes
        with :func:`xarray.apply_ufunc` for gridded PI fields.

    Notes
    -----
    Matches the scalar :func:`tcpyPI.utilities.decompose_pi`, including edge cases:
    if efficiency <= 0, all terms are ``nan`` except ``lnCKCD``; if efficiency > 0
    but pi <= 0, ``lneff`` is returned while ``lnpi`` and ``lndiseq`` are ``nan``.

    Examples
    --------
        >>> pi_log_decomposition(70, 300, 200, CKCD=0.9) == utilities.decompose_pi(70, 300, 200, 0.9)
        True
    """
    units = str(sst_units).upper()
    if units not in ("K", "C"):
        raise ValueError(f"Unsupported sst_units={sst_units!r}; expected 'C' or 'K'.")

    pi_arr = np.asarray(pi, dtype=float)
    t0_arr = np.asarray(t0, dtype=float)
    sst_arr = np.asarray(sst, dtype=float)
    # Convert SST to kelvin if given in Celsius (matches utilities.T_Ctok).
    sst_k = sst_arr + 273.15 if units == "C" else sst_arr

    # Fast path: preserve the exact scalar behavior of utilities.decompose_pi.
    if pi_arr.ndim == 0 and t0_arr.ndim == 0 and sst_k.ndim == 0:
        return utilities.decompose_pi(float(pi_arr), float(sst_k), float(t0_arr), CKCD=CKCD)

    pi_arr, sst_k, t0_arr = np.broadcast_arrays(pi_arr, sst_k, t0_arr)
    lnCKCD = float(np.log(CKCD))
    efficiency = (sst_k - t0_arr) / t0_arr
    valid_eff = efficiency > 0
    valid = valid_eff & (pi_arr > 0)

    lneff = np.full(efficiency.shape, np.nan)
    lneff[valid_eff] = np.log(efficiency[valid_eff])
    # `lnpi` requires efficiency > 0 as well as pi > 0 so the array path matches the
    # scalar decompose_pi (all terms nan except lnCKCD when efficiency <= 0).
    lnpi = np.full(pi_arr.shape, np.nan)
    lnpi[valid] = 2.0 * np.log(pi_arr[valid])
    lndiseq = np.full(pi_arr.shape, np.nan)
    lndiseq[valid] = lnpi[valid] - lneff[valid] - lnCKCD

    def _s(a):
        return float(a) if a.ndim == 0 else a

    return _s(lnpi), _s(lneff), _s(lndiseq), lnCKCD
