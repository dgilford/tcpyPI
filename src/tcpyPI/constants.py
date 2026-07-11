"""Meteorological and algorithmic constants for tcpyPI.

Values follow Emanuel's *Atmospheric Convection* (E94) and the Bister & Emanuel
(2002) MATLAB implementation (``pcmin.m``). These are module-level attributes
(``from tcpyPI import constants``); each is also annotated inline below with units.

Constants
---------
CPD : float
    Specific heat of dry air at constant pressure [J/kg/K].
CPV : float
    Specific heat of water vapor at constant pressure [J/kg/K].
CL : float
    Specific heat of liquid water [J/kg/K]. Set to Emanuel's **modified** value of
    2500 (rather than the physical ~4190) for internal consistency with the pcmin
    entropy formulation; this is why :func:`tcpyPI.utilities.Lv` is not the standard
    Kirchhoff latent heat.
CPVMCL : float
    ``CPV - CL`` [J/kg/K]; sets the temperature dependence of the latent heat.
RV : float
    Gas constant for water vapor [J/kg/K].
RD : float
    Gas constant for dry air [J/kg/K].
EPS : float
    ``RD / RV``, the ratio of gas constants (~0.622) [unitless].
ALV0 : float
    Latent heat of vaporization at 0 degrees C [J/kg].
A, B : float
    Empirical parameters in the lifted-condensation-level pressure estimate
    (E94 4.6.24; see :func:`tcpyPI.utilities.e_pLCL`) [unitless].
b : float
    Exponent for the assumed azimuthal velocity profile in the eye,
    ``V = V_m (r/r_m)**b`` (Emanuel 1995, Eq. 25); used in the PI pressure factor
    [unitless].
ptop : float
    Default pressure below which the sounding is ignored in the CAPE integration
    [hPa].
"""

#
#   ***  Define constants   ***
#

# Thermodynamic Constants
CPD = 1005.7        # [J/kg/K] Specific heat of dry air at constant pressure
CPV = 1870.0        # [J/kg/K] Specific heat of water vapor at constant pressure
# CL = 4190.0       # [J/kg/K] Specific heat of liquid water (physical value; unused)
CL = 2500.0         # [J/kg/K] Modified specific heat of liquid water (Emanuel/pcmin value)
CPVMCL = CPV - CL   # [J/kg/K] CPV - CL
RV = 461.5          # [J/kg/K] gas constant of water vapor
RD = 287.04         # [J/kg/K] gas constant of dry air
EPS = RD / RV       # [unitless] epsilon, the ratio of gas constants (~0.622)
ALV0 = 2.501e6      # [J/kg] latent heat of vaporization at 0 degrees C

# pLCL Empirical Parameters (E94 4.6.24; see utilities.e_pLCL)
A = 1669.0          # [unitless] empirical LCL-pressure parameter
B = 122.0           # [unitless] empirical LCL-pressure parameter

# PI Auxiliaries
b = 2.0             # [unitless] azimuthal-velocity exponent, V=V_m(r/r_m)**b (Emanuel 1995, Eq. 25)
ptop = 50           # [hPa] default pressure below which the sounding is ignored
