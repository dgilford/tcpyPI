# This module stores functions used throughout pyPI
# Last updated 4/16/2020
#

# import numpy to use in functions
import numpy as np

# ---------------------- Longitude conversion ---------------------- %

# define function to convert longitudes from 0-360E to 180W-180E, and vice-versa
def convert_lon_to180(lon360):
    lon180 = (lon360-1e-6 + 180) % 360 - 180
    return(lon180)

def convert_lon_to360(lon180):
    lon360 = lon180 % 360
    return(lon360)

# ---------------------- Analyses ---------------------- %

# define function to calculate TC efficiency
def pi_effiency(sstk,t0):
    # calculate efficiency with SST (K) and Outflow temperature (K)
    efficiency=(sstk-t0)/t0
    return(efficiency)

# define function to calculate TC air-sea thermodynamic disequilibrium
# as a residual from Bister and Emanuel (1998; BE98) EQN. 21
def pi_diseq_resid(pi,sstk,t0,CKCD=0.9):
    # calculate efficiency with SST (K) and Outflow temperature (K)
    efficiency=(sstk-t0)/t0
    # calculate disequilibrium with the BE98 equality
    diseq=pi**2/(CKCD*efficiency)
    return(diseq)

# define function to perform decomposition of TC PI terms
# from Wing et al. (2015), EQN. 2
def decompose_pi(pi,sstk,t0,CKCD=0.9):
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