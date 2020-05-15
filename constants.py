# Constants file supporting pyPI: Potential Intensity Calculation in Python
#
# Last updated 5/14/20
# 

#
#   ***  Define constants   ***
#

# Thermodynamic Constants
CPD=1005.7       # [J/kg.K] Specific heat of dry air at constant pressure
CPV=1870.0       # [J/kg.K] Specific heat of water vapor at constant pressure
#CL=4190.0       # [J/kg.K] Specific heat of liquid water
CL=2500.0        # [J/kg.K] Modified specific heat of liquid water
CPVMCL=CPV-CL
RV=461.5         # [J/kg.K] gas constant of water vapor
RD=287.04        # [J/kg.K] gas constant of dry air
EPS=RD/RV        # [unitless] epsilon, the ratio of gas constants
ALV0=2.501e6     # [J/kg] latent heat of vaporization at 0 degrees C

# pLCL Empirical Parameters
A=1669.0
B=122.0

# PI Auxiliaries
b=2.0        # Exponent for estimating azimuthal velocity in the eye, V=V_m(r/r_m)**b (Emanuel 1995, EQN. 25)
ptop=50      # Pressure below which sounding is ignored (hPa)