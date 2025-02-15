"""Helper functions used throughout tcpyPI"""
# doctest: +ELLIPSIS

import numpy as np

from . import constants
from .numba import njit

# ---------------------- Longitude conversion ---------------------- %

def convert_lon_to180(lon360):
    """Convert longitudes from 0-360E to 180W-180E.
    
    Args:
        lon360 (float): Longitude in 0-360E format
        
    Returns:
        float: Longitude in 180W-180E format
        
    Examples:
        >>> convert_lon_to180(0)
        0
        >>> convert_lon_to180(180)
        180
        >>> convert_lon_to180(181)
        -179
        >>> convert_lon_to180(360)
        0
    """
    lon180 = -((-lon360 + 180) % 360 - 180)
    return(lon180)

def convert_lon_to360(lon180):
    """Convert longitudes from 180W-180E to 0-360E.
    
    Args:
        lon180 (float): Longitude in 180W-180E format
        
    Returns:
        float: Longitude in 0-360E format
        
    Examples:
        >>> convert_lon_to360(-1)
        359
        >>> convert_lon_to360(0)
        0
        >>> convert_lon_to360(1)
        1
        >>> convert_lon_to360(-180)
        180
        >>> convert_lon_to360(180)
        180
    """
    lon360 = lon180 % 360
    return(lon360)

# ---------------------- Unit Conversions ---------------------- %

@njit()
def T_ktoC(Tk):
    """Convert kelvin to degrees Celsius.
    
    Args:
        Tk (float): Temperature in kelvin
        
    Returns:
        float: Temperature in degrees Celsius
        
    Examples:
        >>> T_ktoC(273.15)
        0.0
        >>> T_ktoC(283.15)
        10.0
    """
    return (Tk-273.15)

@njit()
def T_Ctok(TC):
    """Convert degrees Celsius to kelvin.
    
    Args:
        TC (float): Temperature in degrees Celsius
        
    Returns:
        float: Temperature in kelvin
        
    Examples:
        >>> T_Ctok(0)
        273.15
        >>> T_Ctok(10)
        283.15
    """
    return (TC+273.15)

# ---------------------- Thermodynamic Calculations ---------------------- %

@njit()
def es_cc(TC):
    """Calculate saturated water vapor pressure from Clausius-Clapeyron relation/August-Roche-Magnus formula.
    
    Args:
        TC (float): Temperature in Celsius
        
    Returns:
        float: Saturation vapor pressure in hPa
        
    Examples:
        >>> es_cc(20)
        23.369...
        >>> es_cc(0)
        6.112...
    """
    return 6.112*np.exp(17.67*TC/(243.5+TC))

@njit()
def Lv(TC):
    """Calculate latent heat of vaporization as a function of temperature.
    
    Args:
        TC (float): Temperature in Celsius
        
    Returns:
        float: Latent heat of vaporization in J/kg
        
    Examples:
        >>> Lv(0)
        2501000.0
        >>> Lv(20)
        2488400.0
    """
    return constants.ALV0+constants.CPVMCL*TC

@njit()
def ev(R,P):
    """Calculate parcel vapor pressure.
    
    Args:
        R (float): Mixing ratio in gram/gram
        P (float): Pressure in hPa
        
    Returns:
        float: Vapor pressure in hPa
        
    Examples:
        >>> ev(0.01, 1000)
        15.823...
        >>> ev(0.02, 1000)
        31.154...
    """
    return R*P/(constants.EPS+R)

@njit()
def rv(E,P):
    """Calculate parcel mixing ratio.
    
    Args:
        E (float): Vapor pressure in hPa
        P (float): Pressure in hPa
        
    Returns:
        float: Mixing ratio in gram/gram
        
    Examples:
        >>> rv(15.942, 1000)
        0.010076...
        >>> rv(31.250, 1000)
        0.020063...
    """
    return constants.EPS*E/(P-E)

@njit()
def entropy_S(T,R,P):
    """Calculate total specific entropy per unit mass of dry air (E94, EQN. 4.5.9).
    
    Args:
        T (float): Temperature in kelvin
        R (float): Mixing ratio in gram/gram
        P (float): Pressure in hPa
        
    Returns:
        float: Specific entropy
        
    Examples:
        >>> entropy_S(300, 0.01, 1000)
        3987.17...
        >>> entropy_S(290, 0.005, 900)
        3868.01...
    """
    EV=ev(R,P)
    ES=es_cc(T-273.15)
    RH=min([EV/ES,1.0])
    ALV=Lv(T-273.15)
    S=(constants.CPD+R*constants.CL)*np.log(T)-constants.RD*np.log(P-EV)+ALV*R/T-R*constants.RV*np.log(RH)
    return(S)

@njit()
def Trho(T,RT,R):
    """Calculate density temperature in K.
    
    Args:
        T (float): Temperature in kelvin
        RT (float): Total water content mixing ratio in gram/gram
        R (float): Parcel water vapor mixing ratio in gram/gram
        
    Returns:
        float: Density temperature in kelvin
        
    Examples:
        >>> Trho(300, 0.02, 0.01)
        298.84...
        >>> Trho(290, 0.01, 0.005)
        289.43...
    """
    return T*(1.+R/constants.EPS)/(1.+RT)

@njit()
def e_pLCL(TP,RH,PP):
    """Calculate empirical lifting condensation level pressure (pLCL).
    
    Args:
        TP (float): Parcel temperature in kelvin
        RH (float): Relative humidity (0-1)
        PP (float): Parcel pressure in hPa
        
    Returns:
        float: LCL pressure in hPa
        
    Examples:
        >>> e_pLCL(300, 0.8, 1000)
        948.70...
        >>> e_pLCL(290, 0.7, 900)
        830.83...
    """
    return PP*(RH**(TP/(constants.A-constants.B*RH-TP)))

# ---------------------- Analyses ---------------------- %

@njit()
def pi_efficiency(sstk,t0):
    """Calculate TC efficiency.
    
    Args:
        sstk (float): Sea surface temperature in kelvin
        t0 (float): Outflow temperature in kelvin
        
    Returns:
        float: TC efficiency
        
    Examples:
        >>> pi_efficiency(300, 200)
        0.5
        >>> pi_efficiency(295, 210)
        0.404...
    """
    efficiency=(sstk-t0)/t0
    return(efficiency)

@njit()
def pi_diseq_resid(pi,sstk,t0,CKCD=0.9):
    """Calculate TC air-sea thermodynamic disequilibrium as a residual from Bister and Emanuel (1998; BE98) EQN. 21.
    
    Args:
        pi (float): Potential intensity in m/s
        sstk (float): Sea surface temperature in kelvin
        t0 (float): Outflow temperature in kelvin
        CKCD (float, optional): Ratio of exchange coefficients. Defaults to 0.9.
        
    Returns:
        float: Disequilibrium
        
    Examples:
        >>> pi_diseq_resid(70, 300, 200, 0.9)
        10888.88...
        >>> pi_diseq_resid(50, 295, 210, 0.9)
        6862.74...
    """
    efficiency=(sstk-t0)/t0
    diseq=pi**2/(CKCD*efficiency)
    return(diseq)

@njit()
def decompose_pi(pi,sstk,t0,CKCD=0.9):
    """Perform decomposition of TC PI terms from Wing et al. (2015), EQN. 2.
    
    Args:
        pi (float): Potential intensity in m/s
        sstk (float): Sea surface temperature in kelvin
        t0 (float): Outflow temperature in kelvin
        CKCD (float, optional): Ratio of exchange coefficients. Defaults to 0.9.
        
    Returns:
        tuple: (lnpi, lneff, lndiseq, lnCKCD)
        
    Examples:
        >>> result = decompose_pi(70, 300, 200, 0.9)
        >>> [x for x in result]
        [8.496..., -0.693..., 9.295..., -0.105...]
        >>> result = decompose_pi(50, 295, 210, 0.9)
        >>> [x for x in result]
        [7.824..., -0.904..., 8.833..., -0.105...]

    Exceptional cases:
        - Efficiency is non-positive

            >>> decompose_pi(70, 300, 300, 0.9)
            (nan, nan, nan, -0.105...)

        - Potential intensity is non-positive

            >>> decompose_pi(0, 300, 200, 0.9)
            (nan, -0.693..., nan, -0.105...)
    """
    # the natural log of Ck/CD is a constant
    lnCKCD=np.log(CKCD)
    # calculate efficiency with SST (K) and Outflow temperature (K)
    efficiency=(sstk-t0)/t0
    # Is efficiency greater than 0? If not, return missing
    if (efficiency>0):
        # calculate disequilibrium with the BE98 equality
        diseq=pi**2/(CKCD*efficiency)
        lneff=np.log(efficiency)
        # Is potential intensity greater than 0? If not, return missing
        if pi>0:
            lnpi=2*np.log(pi)
            # the log of disequilibrium is calculated as a residual
            lndiseq=lnpi-lneff-lnCKCD
            return(lnpi,lneff,lndiseq,lnCKCD)
        else:
            return(np.nan,lneff,np.nan,lnCKCD)
    else:
        return(np.nan,np.nan,np.nan,lnCKCD)