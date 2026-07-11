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
"""
# doctest: +ELLIPSIS

# import required packages
import numpy as np
from .numba import njit
from . import constants
from . import utilities

# define the function to calculate CAPE
@njit()
def cape(TP,RP,PP,T,R,P,ascent_flag=0,ptop=50,miss_handle=1):
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
            Convective available potential energy (J/kg).
        - TOB : float
            Temperature at level of neutral buoyancy (K).
        - LNB : float
            Level of neutral buoyancy pressure (hPa).
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
    valid_i=~np.isnan(T)
    first_valid=np.where(valid_i)[0][0]
    # Are there missing values? If so, assess according to flag
    if (np.sum(valid_i) != len(P)):
        # if not allowed, set IFLAG=3 and return missing CAPE
        if (miss_handle != 0):
            CAPED=np.nan
            TOB=np.nan
            LNB=np.nan
            IFLAG=3
            # Return the unsuitable values
            return(CAPED,TOB,LNB,IFLAG)
        else:
            # if allowed, but there are missing values between the lowest existing level
            # and ptop, then set IFLAG=3 and return missing CAPE
            if np.sum(np.isnan(T[first_valid:len(P)])>0):
                CAPED=np.nan
                TOB=np.nan
                LNB=np.nan
                IFLAG=3
                # Return the unsuitable values
                return(CAPED,TOB,LNB,IFLAG)
            else:
                first_lvl=first_valid
    else:
        first_lvl=0

    # Populate new environmental profiles removing values above ptop and
    # find new number, N, of profile levels with which to calculate CAPE.
    # Keep every level with P > ptop (matching pcmin); with the profile ordered
    # by decreasing pressure these are the leading N entries. Using argmin(|P-ptop|)
    # would drop the topmost retained level whenever no level sits exactly at ptop.
    N=int(np.count_nonzero(P > ptop))

    P=P[first_lvl:N]
    T=T[first_lvl:N]
    R=R[first_lvl:N]
    nlvl=len(P)
    TVRDIF = np.zeros((nlvl,))
    
    #
    #   ***  Run checks   ***
    #
    
    # CHECK: is the input profile ordered with increasing pressure? If not, return missing CAPE
    if (P[2]-P[1] > 0):
        CAPED=0
        TOB=np.nan
        LNB=np.nan
        IFLAG=0
        # Return the unsuitable values
        return(CAPED,TOB,LNB,IFLAG)

    # CHECK: Is the input parcel suitable? If not, return missing CAPE
    if ((RP < 1e-6) or (TP < 200)):
        CAPED=0
        TOB=np.nan
        LNB=np.nan
        IFLAG=0
        # Return the unsuitable values
        return(CAPED,TOB,LNB,IFLAG)
    
    #
    #  ***  Define various parcel quantities, including reversible   ***
    #  ***                       entropy, S                          ***
    #                         
    TPC=utilities.T_ktoC(TP)                 # Parcel temperature in Celsius
    ESP=utilities.es_cc(TPC)                # Parcel's saturated vapor pressure
    EVP=utilities.ev(RP,PP)                 # Parcel's partial vapor pressure
    RH=EVP/ESP                              # Parcel's relative humidity
    RH=min([RH,1.0])                        # ensure that the relatively humidity does not exceed 1.0
    # calculate reversible total specific entropy per unit mass of dry air (E94, EQN. 4.5.9)
    S=utilities.entropy_S(TP,RP,PP)
    
    
    #
    #   ***  Estimate lifted condensation level pressure, PLCL   ***
    #     Based on E94 "calcsound.f" code at http://texmex.mit.edu/pub/emanuel/BOOK/
    #     see also https://psl.noaa.gov/data/composites/day/calculation.html
    #
    #   NOTE: Modern PLCL calculations are made following the exact expressions of Romps (2017),
    #   see https://journals.ametsoc.org/doi/pdf/10.1175/JAS-D-17-0102.1
    #   and Python PLCL code at http://romps.berkeley.edu/papers/pubdata/2016/lcl/lcl.py
    #
    PLCL=utilities.e_pLCL(TP,RH,PP)
    
    # Initial default values before loop
    CAPED=0
    TOB=T[0]
    IFLAG=1
    # Values to help loop
    NCMAX=0
    jmin=int(1e6)
    
    #
    #   ***  Begin updraft loop   ***
    #

    # loop over each level in the profile
    for j in range(nlvl):
        
        # jmin is the index of the lowest pressure level evaluated in the loop
        jmin=int(min([jmin,j]))
    
        #
        #   *** Calculate Parcel quantities BELOW lifted condensation level   ***
        #
        if (P[j] >= PLCL):
            # Parcel temperature at this pressure
            TG=TP*(P[j]/PP)**(constants.RD/constants.CPD)
            # Parcel Mixing ratio
            RG=RP
            # Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
            TLVR=utilities.Trho(TG,RG,RG)
            TVENV=utilities.Trho(T[j],R[j],R[j])
            # Bouyancy of the parcel in the environment (Proxy of E94, EQN. 6.1.5)
            TVRDIF[j,]=TLVR-TVENV
            
        #
        #   *** Calculate Parcel quantities ABOVE lifted condensation level   ***
        # 
        else:
            TG, RG, IFLAG = solve_temperature_from_entropy(S=S, P=P[j], RP=RP, T_initial=T[j])
            if IFLAG == 2:  # Did not converge
                CAPED=0
                TOB=T[0]
                LNB=P[0]
                # Return the uncoverged values
                return(CAPED,TOB,LNB,IFLAG)

            #
            #   *** Calculate buoyancy   ***
            #
            # Parcel total mixing ratio: either reversible (ascent_flag=0) or pseudo-adiabatic (ascent_flag=1)
            RMEAN=ascent_flag*RG+(1-ascent_flag)*RP
            # Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
            TLVR=utilities.Trho(TG,RMEAN,RG)
            TENV=utilities.Trho(T[j],R[j],R[j])
            # Bouyancy of the parcel in the environment (Proxy of E94, EQN. 6.1.5)
            TVRDIF[j,]=TLVR-TENV
            

    #
    #  ***  Begin loop to find Positive areas (PA) and Negative areas (NA) ***
    #                  ***  and CAPE from reversible ascent ***
    NA=0.0
    PA=0.0
    
    #
    #   ***  Find maximum level of positive buoyancy, INB    ***
    #
    INB=0
    for j in range(nlvl-1, jmin, -1):
        if (TVRDIF[j] > 0):
            INB=max([INB,j])
            
    # CHECK: Is the LNB higher than the surface? If not, return zero CAPE  
    if (INB==0):
        CAPED=0
        TOB=T[0]
        LNB=P[INB]
#         TOB=np.nan
        LNB=0
        # Return the unconverged values
        return(CAPED,TOB,LNB,IFLAG)
    
    # if check is passed, continue with the CAPE calculation
    else:
    
    #
    #   ***  Find positive and negative areas and CAPE  ***
    #                  via E94, EQN. 6.3.6)
    #
        for j in range(jmin+1, INB+1, 1):
            PFAC=constants.RD*(TVRDIF[j]+TVRDIF[j-1])*(P[j-1]-P[j])/(P[j]+P[j-1])
            PA=PA+max([PFAC,0.0])
            NA=NA-min([PFAC,0.0])

    #
    #   ***   Find area between parcel pressure and first level above it ***
    #
        PMA=(PP+P[jmin])
        PFAC=constants.RD*(PP-P[jmin])/PMA
        PA=PA+PFAC*max([TVRDIF[jmin],0.0])
        NA=NA-PFAC*min([TVRDIF[jmin],0.0])
        
    #
    #   ***   Find residual positive area above INB and TO  ***
    #         and finalize estimate of LNB and its temperature
    #
        PAT=0.0
        TOB=T[INB]
        LNB=P[INB]
        if (INB < nlvl-1):
            PINB=(P[INB+1]*TVRDIF[INB]-P[INB]*TVRDIF[INB+1])/(TVRDIF[INB]-TVRDIF[INB+1])
            LNB=PINB
            PAT=constants.RD*TVRDIF[INB]*(P[INB]-PINB)/(P[INB]+PINB)
            TOB=(T[INB]*(PINB-P[INB+1])+T[INB+1]*(P[INB]-PINB))/(P[INB]-P[INB+1])
    
    #
    #   ***   Find CAPE  ***
    #            
        CAPED=PA+PAT-NA
        CAPED=max([CAPED,0.0])
        # set the flag to OK if procedure reached this point
        IFLAG=1
        # Return the calculated outputs to the above program level 
        return(CAPED,TOB,LNB,IFLAG)


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
    TGNEW=T_initial
    TJC=utilities.T_ktoC(T_initial)
    ES=utilities.es_cc(TJC)
    RG=utilities.rv(ES,P)
    
    #
    #   ***  Iteratively calculate lifted parcel temperature and mixing   ***
    #   ***                ratio for reversible ascent                    ***
    #
    
    # set loop counter and initial condition
    NC=0
    TG=0

    # loop until loop converges or bails out
    while ((np.abs(TGNEW-TG)) > 0.001):
    
        # Parcel temperature and mixing ratio during this iteration
        TG=TGNEW
        TC=utilities.T_ktoC(TG)
        ENEW=utilities.es_cc(TC)
        RG=utilities.rv(ENEW,P)
        
        # increase iteration count in the loop
        NC += 1
        
        #
        #   ***  Calculate estimates of the rates of change of the entropy    ***
        #   ***           with temperature at constant pressure               ***
        #

        ALV=utilities.Lv(TC)
        # calculate the rate of change of entropy with temperature, s_ell
        SL=(constants.CPD+RP*constants.CL+ALV*ALV*RG/(constants.RV*TG*TG))/TG
        EM=utilities.ev(RG,P)
        # calculate the saturated entropy, s_k, noting r_T=RP and
        # the last term vanishes with saturation, i.e. RH=1
        SG=(constants.CPD+RP*constants.CL)*np.log(TG)-constants.RD*np.log(P-EM)+ALV*RG/TG
        # convergence speed (AP, step in entropy fraction) varies as a function of 
        # number of iterations
        if (NC < 3):
            # converge slowly with a smaller step
            AP=0.3
        else:
            # speed the process with a larger step when nearing convergence
            AP=1.0
        # find the new temperature in the iteration
        TGNEW=TG+AP*(S-SG)/SL
        
        #
        #   ***   If the routine does not converge, set IFLAG=2 and bail out   ***
        #
        if (NC > 500) or (ENEW > (P-1)):
            IFLAG = 2  # Did not converge
            return TG, RG, IFLAG
        
        # store the number of iterations
        NCMAX=NC

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
    SSTK=utilities.T_Ctok(SSTC) # SST in kelvin
    T=utilities.T_Ctok(TC)      # Temperature profile in kelvin
    R=R*0.001                   # Mixing ratio profile in g/g

    # CHECK 1a: do SSTs exceed 5C?
    # CHECK 1b: are SSTs less than 100C (if not, indicative of input in kelvin)
    # If not, set IFL=0 and return missing PI
    if (SSTC <= 5.0 or SSTC>100):
        VMAX=np.nan
        PMIN=np.nan
        IFL=0
        TO=np.nan
        OTL=np.nan
        return(VMAX,PMIN,IFL,TO,OTL)
    # CHECK 1c: are SSTs missing? If so, set IFL=0 and return missing PI
    elif (np.isnan(SSTC)==True):
        VMAX=np.nan
        PMIN=np.nan
        IFL=0
        TO=np.nan
        OTL=np.nan
        return(VMAX,PMIN,IFL,TO,OTL)

    # CHECK 2a: do Temperature profiles exceed 100K?
    # CHECK 2b: are Temperatures in Celsius less than 100C (if not, indicative of input in kelvin)
    # If not, set IFL=0 and return missing PI
    if (np.min(T) <= 100 or np.max(TC)>100):
        VMAX=np.nan
        PMIN=np.nan
        IFL=0
        TO=np.nan
        OTL=np.nan
        return(VMAX,PMIN,IFL,TO,OTL)
    
    # Set Missing mixing ratios to zero g/g, following Kerry's BE02 algorithm
    R[np.isnan(R)]=0.
    
    # Saturated water vapor pressure
    # from Clausius-Clapeyron relation/August-Roche-Magnus formula
    ES0=utilities.es_cc(SSTC)

    # define the level from which parcels are lifted. Normally the lowest level
    # (index 0). Under miss_handle=0 (ignore-NaN), skip any leading missing
    # temperature levels (e.g. sub-surface levels over topography) so the parcel is
    # lifted from the lowest *valid* level rather than a NaN.
    NK=0
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
    IFL=int(1)
    TP=T[NK]
    RP=R[NK]
    PP=P[NK]
    result = cape(TP,RP,PP,T,R,P,ascent_flag,ptop,miss_handle)
    CAPEA = result[0]
    IFLAG = result[3]
    # if the CAPE function tripped a flag, set the output IFL to it
    if (IFLAG != 1):
        IFL=int(IFLAG)

    #
    #   ***   Begin iteration to find mimimum pressure   ***
    #

    # set loop counter and initial condition
    NP=0         # loop counter
    PM=970.0
    PMOLD=PM     # initial condition from minimum pressure
    PNEW=0.0     # initial condition from minimum pressure

    TO_ENV=np.nan
    OTL_ENV=np.nan

    # loop until convergence or bail out
    while (np.abs(PNEW-PMOLD) > 0.5):
        
        #
        #   ***  Find CAPE at radius of maximum winds   ***
        #
        TP=T[NK]
        PP=min([PM,1000.0])
        # find the mixing ratio with the average of the lowest level pressure and MSL
        RP=constants.EPS*R[NK]*MSL/(PP*(constants.EPS+R[NK])-R[NK]*MSL)
        CAPEM, TOB_ENV, LNB_ENV, IFLAG = cape(
            TP,RP,PP,T,R,P,ascent_flag,ptop,miss_handle
        )
        # if the CAPE function tripped a different flag, set the output IFL to it
        if (IFLAG != 1):
            IFL=int(IFLAG)
        TO_ENV=TOB_ENV
        OTL_ENV=LNB_ENV
        
        #
        #  ***  Find saturation CAPE at radius of maximum winds    ***
        #  *** Note that TO and OTL are found with this assumption ***
        #
        TP=SSTK
        PP=min([PM,1000.0])
        RP=utilities.rv(ES0,PP)
        result = cape(TP,RP,PP,T,R,P,ascent_flag,ptop,miss_handle)
        CAPEMS, TOMS, LNBS, IFLAG = result
        # if the CAPE function tripped a flag, set the output IFL to it
        if (IFLAG != 1):
            IFL=int(IFLAG)
        # Store the outflow temperature and level of neutral bouyancy at the outflow level (OTL)
        TO=TOMS   
        OTL=LNBS
        # Optionally use CAPEenv (final-iteration) outflow definition.
        if outflow_source_flag == 1:
            TO=TO_ENV
            OTL=OTL_ENV
        # Calculate the proxy for TC efficiency (BE02, EQN. 1-3)
        RAT=SSTK/TO
        # If dissipative heating is "off", TC efficiency proxy is set to 1.0 (BE02, pg. 3)
        if (diss_flag == 0):
            RAT=1.0
        
        #
        #  ***  Initial estimate of pressure at the radius of maximum winds  ***
        #
        RS0=RP
        # Lowest level and Sea-surface Density Temperature (E94, EQN. 4.3.1 and 6.3.7)
        TV0=utilities.Trho(T[NK],R[NK],R[NK])
        TVSST=utilities.Trho(SSTK,RS0,RS0)
        # Average Surface Density Temperature, e.g. 1/2*[Tv(Tsfc)+Tv(sst)]
        TVAV=0.5*(TV0+TVSST)
        # Converge toward CAPE*-CAPEM (BE02, EQN 3-4)
        CAT=(CAPEM-CAPEA)+0.5*CKCD*RAT*(CAPEMS-CAPEM)
        CAT=max([CAT,0.0])
        # Iterate on pressure
        PNEW=MSL*np.exp(-CAT/(constants.RD*TVAV))
        
        #
        #   ***  Test for convergence (setup for possible next while iteration)  ***
        #
        # store the previous step's pressure       
        PMOLD=PM
        # store the current step's pressure
        PM=PNEW
        
        # increase iteration count in the loop
        NP += 1
        #
        #   ***   If the routine does not converge, set IFL=2 and return missing PI   ***
        #
        if (NP > 200)  or (PM < 400):
            VMAX=np.nan
            PMIN=np.nan
            IFL=2
            TO=np.nan
            OTL=np.nan
            return(VMAX,PMIN,IFL,TO,OTL)
    
    # Once converged, set potential intensity at the radius of maximum winds
    CATFAC=0.5*(1.+1/constants.b)
    CAT=(CAPEM-CAPEA)+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
    CAT=max([CAT,0.0])
    
    # Calculate the minimum pressure at the eye of the storm
    # BE02 EQN. 4
    PMIN=MSL*np.exp(-CAT/(constants.RD*TVAV))
                 
    # Calculate the potential intensity at the radius of maximum winds
    # BE02 EQN. 3, reduced by some fraction (default 20%) to account for the reduction 
    # of 10-m winds from gradient wind speeds (Emanuel 2000, Powell 1980)
    FAC=max([0.0,(CAPEMS-CAPEM)])
    VMAX=V_reduc*np.sqrt(CKCD*RAT*FAC)
        
    # Return the calculated outputs to the above program level
    return(VMAX,PMIN,IFL,TO,OTL)


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
            - ``0``: no convergence
            - ``2``: CAPE routine failed to converge
            - ``3``: CAPE routine failed due to missing data
        - TO : float
            Outflow temperature (K).
        - OTL : float
            Outflow temperature level (hPa), defined as the level of neutral buoyancy.

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
    if outflow_source == "cape_star":
        outflow_source_flag = 0
    elif outflow_source == "cape_env":
        outflow_source_flag = 1
    else:
        raise ValueError(
            f"Invalid outflow_source={outflow_source!r}; expected 'cape_star' or 'cape_env'."
        )

    P = np.asarray(P)
    TC = np.asarray(TC)
    R = np.asarray(R)

    # Guard against NaNs the Numba core cannot flag safely. A NaN in the pressure
    # array or in MSL poisons the CAPE integration and the pressure-convergence
    # loop, which can otherwise exit early and return numerically plausible values
    # with IFL=1. Fail fast to IFL=3 (missing data), matching NaN handling of the
    # temperature profile.
    if np.any(np.isnan(P)) or np.any(np.isnan(np.asarray(MSL, dtype=float))):
        return (np.nan, np.nan, 3, np.nan, np.nan)

    # Order-agnostic: the algorithm requires lowest index = lowest level
    # (decreasing pressure). Sort the profile here so any input ordering is
    # accepted; already-descending input is left untouched.
    order = np.argsort(P)[::-1]
    if not np.array_equal(order, np.arange(P.size)):
        P = P[order]
        TC = TC[order]
        R = R[order]

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
