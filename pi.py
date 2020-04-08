# Potential Intensity Calculation for Python
# -----------------------------------------------------------------------------------
# Adapted from pcmin.m by Kerry Emanuel (ftp://texmex.mit.edu/pub/emanuel/TCMAX)
# Originally Updated by Daniel Gilford for
# Gilford et al. (2017) -- https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-16-0827.1
# Gilford et al. (2019) -- https://journals.ametsoc.org/doi/10.1175/MWR-D-19-0021.1
# 
# Adapted for Python (NAME OF REPOSITORY) by Daniel Gilford, PhD (Rutgers U., daniel.gilford@rutgers.edu)
# Last updated 4/1/2020
# -----------------------------------------------------------------------------------
#
# Revision History:
#   Revised on 9/24/2005 by K. Emanuel to fix convergence problems at high pressure
#     Converted to MATLAB  5/17/2008
#
#   Revised 7/20/15 by D. Gilford to output the LNB
#   Revised 8/4/16 by D. Gilford to include lack of convergence if SST < 5C for TO/LNB
#   Revised 8/5/16 by D. Gilford to fix the "cape()" function output and include LNB
#   Revised 10/3/16 by D. Gilford to set LNB to the pressure-weighted crossing of buoyancy from negative to positive (the zero-line)
#     Converted to Python  4/1/2020
# -----------------------------------------------------------------------------------
# 
#   ***    This function calculates the maximum wind speed         ***
#   ***             and mimimum central pressure                   ***
#   ***    achievable in tropical cyclones, given a sounding       ***
#   ***             and a sea surface temperature.                 ***
#
#   Most thermodynamic and dynamic technical background/calculations are found in Bister 
#   and Emanuel (2002; BE02) and Emanuel's "Atmospheric Convection" (E94; 1994; ISBN: 978-0195066302)
#
#  INPUT:   SSTC: Sea surface temperature (C)
#
#           MSL: Mean Sea level pressure (hPa)
#
#           P,T,R: One-dimensional arrays 
#             containing pressure (hPa), temperature (C),
#             and mixing ratio (g/kg). The arrays MUST be
#             arranged so that the lowest index corresponds
#             to the lowest model level, with increasing index
#             corresponding to decreasing pressure. The temperature
#             sounding should extend to at least the tropopause and 
#             preferably to the lower stratosphere, however the
#             mixing ratios are not important above the boundary
#             layer. Missing mixing ratios can be replaced by zeros.
#
#           CKCD: Ratio of C_k to C_D (unitless number)
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!! STOPPED HERE, FILL OUT REST OF INPUTS and documentation
#
#  OUTPUT:  PMIN is the minimum central pressure, in mb
#
#           VMAX is the maximum surface wind speed, in m/s
#                  (reduced to reflect surface drag)
#
#           TO is the outflow temperature (K)
#
#           LNB is the level of neutral bouyancy where the outflow temperature
#               is found (hPa), i.e. where buoyancy is actually equal to zero under the 
#               condition of an air parcel that is saturated at sea level pressure
#
#           IFL is a flag: A value of 1 means OK; a value of 0
#              indicates no convergence; a value of 2
#              means that the CAPE routine failed to converge
#
# -----------------------------------------------------------------------------------
#
#

# import packages we need
import numpy as np

# # define the function to calculate PI
# def pi(SSTC,MSL,P,T,R,CKCD=0.9,ascent_flag=0,diss_flag=1,V_reduc=0.8):
    
#     # convert units
#     SSTK=SSTC+273.15 # SST in kelvin
#     TK=T+273.15      # Temperature profile in kelvin
#     Rgg=R*0.001   # Mixing ratio profile in gm/gm

#     # CHECK 1: do SSTs exceed 5C? If not, return missing PI
#     if (SSTC <= 5.0):
#         VMAX=np.nan
#         PMIN=np.nan
#         IFL=0
#         TO=np.nan
#         LNB=np.nan
#         return(VMAX,PMIN,IFL,TO,LNB)
#     end

#     # CHECK 2: do Temperature profiles exceed 100K? If not, return missing PI
#     if (np.min(T) <= 100):
#         VMAX=np.nan
#         PMIN=np.nan
#         IFL=0
#         TO=np.nan
#         LNB=np.nan
#         return(VMAX,PMIN,IFL,TO,LNB)
#     end

#     # Constants
#     NK=0  # level from which parcels lifted
#     b=2.0 # Exponent in assumed profile of azimuthal velocity in eye, V=V_m(r/r_m)**b

#     # Initial saturated water vapor pressure, from Clausius-Clapeyron relation/August-Roche-Magnus formula
#     ES0=6.112*np.exp(17.67*SSTC/(243.5+SSTC))






#     # Return the calculated outputs to the above program level
#     return(VMAX,PMIN,IFL,TO,LNB)





# define the function to calculate CAPE
def cape(TP,RP,PP,T,R,P,ascent_flag=0,ptop=50):

#     function [CAPED,TOB,LNB,IFLAG]= cape(TP,RP,PP,T,R,P,ascent_flag=0,ptop=50)
#
#       This function calculates the CAPE of a parcel given parcel pressure PP (hPa), 
#       temperature TP (K) and mixing ratio RP (gram/gram) and given a sounding
#       of temperature (T in K) and mixing ratio (R in gram/gram) as a function
#       of pressure (P in hPa). CAPED is the calculated value of CAPE following
#       Emanuel 1994 (E94) Equation 6.3.6 and TOB is the temperature at the
#       level of neutral buoyancy ("LNB") for the displaced parcel. IFLAG is a flag
#       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
#       not run owing to improper sounding (e.g. no water vapor at parcel level).
#       IFLAG=2 indicates that the routine did not converge.           
#
#  INPUT:   TP,RP,PP: floating point numbers of Parcel pressure (hPa), 
#             temperature (K), and mixing ratio (gram/gram)
#
#           T,R,P: One-dimensional arrays 
#             containing environmental pressure (hPa), temperature (K),
#             and mixing ratio (gram/gram) profiles.
#
#           ascent_flag: Adjustable constant integer for buoyancy of displaced  
#           parcels, where 0=Reversible ascent;  1=Pseudo-adiabatic ascent
#
#           ptop: Pressure below which sounding is ignored (hPa)
#
#
#  OUTPUT:  CAPED (J/kg) is Convective Available Potential Energy of an air parcel
#           consistent with its parcel and environmental properties.
#
#           TOB is the Temperature (K) at the level of neutral bouyancy 
#           for the displaced air parcel
#
#           LNB is the pressure level of neutral bouyancy (hPa) for the 
#           displaced air parcel
#
#           IFLAG is a flag where the value of 1 means OK; a value of 0
#              indicates an improper sounding or parcel; a value of 2
#              means that the routine failed to converge
#

    # Populate new environmental profiles removing values above ptop and
    # find new number, N, of profile levels with which to calculate CAPE
    try:
        N=np.where((P-ptop) <= 0)[0][0]
    except:
        N=len(P)
    P=P[:N]
    T=T[:N]
    R=R[:N]
    TVRDIF=np.zeros((N,),dtype='float64')

    # CHECK: Is the input parcel suitable? If not, return missing PI
    if ((RP < 1e-6) or (TP < 200)):
        CAPED=0
        TOB=np.nan
        LNB=np.nan
        IFLAG=0
        # Return the unsuitable values
        return(CAPED,TOB,LNB,IFLAG)

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
    
    #
    #  ***  Define various parcel quantities, including reversible   ***
    #  ***                       entropy, S                          ***
    #                         
    TPC=TP-273.15    # Parcel temperature in Celsius
    ESP=6.112*np.exp(17.67*TPC/(243.5+TPC)) # Parcel's saturated vapor pressure
    EVP=RP*PP/(EPS+RP) # Parcel's partial vapor pressure
    RH=EVP/ESP # Parcel's relative humidity
    RH=np.min([float(RH),1.0]) # ensure that the relatively humidity does not exceed 1
    ALV=ALV0+CPVMCL*TPC # calculate the latent heat of vaporization, dependant on temperature
    # calculate reversible total specific entropy per unit mass of dry air (E94, EQN. 4.5.9)
    S=(CPD+RP*CL)*np.log(TP)-RD*np.log(PP-EVP)+ALV*RP/TP-RP*RV*np.log(RH)
    
    #
    #   ***  Estimate lifted condensation level pressure, PLCL   ***
    #     Based on E94 "calcsound.f" code at http://texmex.mit.edu/pub/emanuel/BOOK/
    #     see also https://psl.noaa.gov/data/composites/day/calculation.html
    #
    #   NOTE: Modern PLCL calculations are made following the exact expressions of Romps (2017),
    #   see https://journals.ametsoc.org/doi/pdf/10.1175/JAS-D-17-0102.1
    #   and Python PLCL code at http://romps.berkeley.edu/papers/pubdata/2016/lcl/lcl.py
    #
    CHI=TP/(1669.0-122.0*RH-TP)
    PLCL=PP*(RH**CHI)
    
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
    for j in range(N):
        
        jmin=int(np.min([jmin,j]))
    
        #
        #   *** Calculate Parcel quantities BELOW lifted condensation level   ***
        #
        if (P[j] >= PLCL):
            # Parcel temperature at this pressure
            TG=TP*(P[j]/PP)**(RD/CPD)
            # Parcel Mixing ratio
            RG=RP
            # Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
            TLVR=TG*(1.+RG/EPS)/(1.+RG)
            TVENV=T[j]*(1.+R[j]/EPS)/(1+R[j])
            # Bouyancy of the parcel in the environment (Proxy of E94, EQN. 6.1.5)
            TVRDIF[j,]=TLVR-TVENV
            
        #
        #   *** Calculate Parcel quantities ABOVE lifted condensation level   ***
        # 
        else:
            TGNEW=T[j]
            TJC=T[j]-273.15
            ES=6.112*np.exp(17.67*TJC/(243.5+TJC))
            RG=EPS*ES/(P[j]-ES)
            
            #
            #   ***  Iteratively calculate lifted parcel temperature and mixing   ***
            #   ***                ratio for reversible ascent                    ***
            #
            
            # set loop counter and initial condition
            NC=0
            TG=0

            # loop until we converge or we break out
            while ((np.abs(TGNEW-TG)) > 0.001):
            
                # Parcel temperature and mixing ratio during this iteration
                TG=TGNEW
                TC=TG-273.15
                ENEW=6.112*np.exp(17.67*TC/(243.5+TC))
                RG=EPS*ENEW/(P[j]-ENEW)
                
                # increase iteration count in the loop
                NC=NC+1
                
                #
                #   ***  Calculate estimates of the rates of change of the entropy    ***
                #   ***           with temperature at constant pressure               ***
                #

                ALV=ALV0+CPVMCL*TC
                SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
                EM=RG*P[j]/(EPS+RG)
                # (last term vanishes with saturation, i.e. RH=1)
                SG=(CPD+RP*CL)*np.log(TG)-RD*np.log(P[j]-EM)+ALV*RG/TG
                # convergence speed (AP, step in entropy fraction) varies as a function of 
                # number of iterations
                if (NC < 3):
                    AP=0.3
                else:
                    AP=1.0
                # find the new temperature in the iteration
                TGNEW=TG+AP*(S-SG)/SL
                
                #
                #   ***   If the routine does not converge, set IFLAG=2 and bail out   ***
                #
                if (NC > 500) or (ENEW > (P[j]-1)):
                    CAPED=0
                    TOB=T[0]
                    LNB=P[0]
                    IFLAG=2
                    # Return the uncoverged values
                    return(CAPED,TOB,LNB,IFLAG)
                
                # store the number of iterations
                NCMAX=NC
                
            #
            #   *** Calculate buoyancy   ***
            #
            # Parcel Mixing ratio: either reversible (ascent_flag=0) or pseudo-adiabatic (ascent_flag=1)
            RMEAN=ascent_flag*RG+(1-ascent_flag)*RP
            # Parcel and Environmental Density Temperatures at this pressure (E94, EQN. 4.3.1 and 6.3.7)
            TLVR=TG*(1.+RG/EPS)/(1.+RMEAN)
            TENV=T[j]*(1.+R[j]/EPS)/(1.+R[j])
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
    for j in range(N-1, jmin, -1):
        if (TVRDIF[j] > 0):
            INB=np.max([INB,j])
            
    # CHECK: Is the LNB higher than the surface? If not, return zero CAPE  
    if (INB==0):
        CAPED=0
        TOB=T[0]
        LNB=P[0]
        IFLAG=2;
        # Return the uncoverged values
        return(CAPED,TOB,LNB,IFLAG)
    
    # if we pass the check, continue with the CAPE calculation
    else:
    
    #
    #   ***  Find positive and negative areas and CAPE  ***
    #                  via E94, EQN. 6.3.6)
    #
        for j in range(jmin+1, INB+1, 1):
            PFAC=RD*(TVRDIF[j]+TVRDIF[j-1])*(P[j-1]-P[j])/(P[j]+P[j-1])
            PA=PA+np.max([PFAC,0.0])
            NA=NA-np.min([PFAC,0.0])

    #
    #   ***   Find area between parcel pressure and first level above it ***
    #
        PMA=(PP+P[jmin])
        PFAC=RD*(PP-P[jmin])/PMA
        PA=PA+PFAC*np.max([TVRDIF[jmin],0.0])
        NA=NA-PFAC*np.min([TVRDIF[jmin],0.0])
        
    #
    #   ***   Find residual positive area above INB and TO  ***
    #         and finalize estimate of LNB and its temperature
    #
        PAT=0.0
        TOB=T[INB]
        LNB=P[INB]
        if (INB < N):
            PINB=(P[INB+1]*TVRDIF[INB]-P[INB]*TVRDIF[INB+1])/(TVRDIF[INB]-TVRDIF[INB+1])
            LNB=PINB;
            PAT=RD*TVRDIF[INB]*(P[INB]-PINB)/(P[INB]+PINB);
            TOB=(T[INB]*(PINB-P[INB+1])+T[INB+1]*(P[INB]-PINB))/(P[INB]-P[INB+1]);
    
    #
    #   ***   Find CAPE  ***
    #            
        CAPED=PA+PAT-NA;
        CAPED=np.max([CAPED,0.0]);
        # set the flag to OK if we reached this far
        IFLAG=1
        # Return the calculated outputs to the above program level 
        return(CAPED,TOB,LNB,IFLAG)

    
    