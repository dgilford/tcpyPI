%
          function [PMIN,VMAX,TO,LNB,IFL]= pc_min(SST,PSL,P,T,R)
%
%   Revised on 9/24/2005 to fix convergence problems at high pressure
%     Converted to MATLAB  5/17/2008
%
%   Revised 7/20/15 by D. Gilford to output the LNB
%   Revised 8/4/16 by D. Gilford to include lack of convergence if SST < 5C for TO/LNB
%   Revised 8/5/16 by D. Gilford to fix the "cape()" function output and include LNB
%   Revised 10/3/16 by D. Gilford to set the LNB to the pressure-weighted crossing of the buoyancy from negative to positive (the zero-line)
%
%   ***    This function calculates the maximum wind speed         ***
%   ***             and mimimum central pressure                   ***
%   ***    achievable in tropical cyclones, given a sounding       ***
%   ***             and a sea surface temperature.                 ***
%
%  INPUT:   SST: Sea surface temperature in C
%
%           PSL: Sea level pressure (mb)
%
%           P,T,R: One-dimensional arrays 
%             containing pressure (mb), temperature (C),
%             and mixing ratio (g/kg). The arrays MUST be
%             arranged so that the lowest index corresponds
%             to the lowest model level, with increasing index
%             corresponding to decreasing pressure. The temperature
%             sounding should extend to at least the tropopause and 
%             preferably to the lower stratosphere, however the
%             mixing ratios are not important above the boundary
%             layer. Missing mixing ratios can be replaced by zeros.
%
%
%  OUTPUT:  PMIN is the minimum central pressure, in mb
%
%           VMAX is the maximum surface wind speed, in m/s
%                  (reduced to reflect surface drag)
%
%           TO is the outflow temperature (K)
%
%	    LNB is the level of neutral bouyancy where the outflow temperature
%		is found (hPa), i.e. where buoyancy is actually equal to zero under the 
%		condition of an air parcel that is saturated at sea level pressure
%
%           IFL is a flag: A value of 1 means OK; a value of 0
%              indicates no convergence (hypercane); a value of 2
%              means that the CAPE routine failed.
%
%-----------------------------------------------------------------------------
%
%   ***   Adjustable constant: Ratio of C_k to C_D    ***
%
	CKCD=0.9;
%
%   ***   Adjustable constant for buoyancy of displaced parcels:  ***
%   ***    0=Reversible ascent;  1=Pseudo-adiabatic ascent        ***
%
    SIG=0.0;
%
%   ***  Adjustable switch: if IDISS = 0, no dissipative heating is   ***
%   ***     allowed; otherwise, it is                                 ***
%
	IDISS=1;
%
%   ***  Exponent, b, in assumed profile of azimuthal velocity in eye,   ***
%   ***   V=V_m(r/r_m)^b. Used only in calculation of central pressure   ***
%
	b=2.0;
%
%   *** Set level from which parcels lifted   ***
%
	NK=1;
%
%   *** Factor to reduce gradient wind to 10 m wind   ***
%
	VREDUC=0.8;
%
%--------------------------------------------------------------------------
%
	SSTK=SST+273.15;
	TOMS=230.0;
%
    if SST <= 5.0
	 VMAX=0.0;
	 PMIN=0.0;
	 IFL=0;
	 TO=0;
	 LNB=0;
	 return
    end
%
	ES0=6.112.*exp(17.67.*SST./(243.5+SST));
%    
	 R=R.*0.001;
	 T=T+273.15;
%
    if min(T) <= 100.0
	  VMAX=0.0;
	  PMIN=0.0;
	  IFL=0;
	  return
    end
%
%   ***   Default value   ***
%
	IFL=1;
%
	NP=0;
	PM=970.0;
	PMOLD=PM;
    PNEW=0.0;
%
%   ***   Find environmental CAPE *** 
%
      TP=T(NK);
      RP=R(NK);
      PP=P(NK);
      [CAPEA, ~, ~, IFLAG]= cape(TP,RP,PP,T,R,P,SIG);
      if IFLAG ~= 1
          IFL=2;
      end    
%
%   ***   Begin iteration to find mimimum pressure   ***
%
      while (abs(PNEW-PMOLD)) > 0.5
%
%   ***  Find CAPE at radius of maximum winds   ***
%
      TP=T(NK);
      PP=min(PM,1000.0);
      RP=0.622.*R(NK).*PSL./(PP.*(0.622+R(NK))-R(NK).*PSL);
      [CAPEM, TOM, ~, IFLAG]=cape(TP,RP,PP,T,R,P,SIG);
      if IFLAG ~= 1
          IFL=2;
      end    
%
%  ***  Find saturation CAPE at radius of maximum winds   ***
%
      TP=SSTK;
      PP=min(PM,1000.0);
      RP=0.622.*ES0./(PP-ES0);
      [CAPEMS, TOMS, LNB, IFLAG]=cape(TP,RP,PP,T,R,P,SIG);
      TO=TOMS;
      if IFLAG ~= 1
          IFL=2;
      end    
      RAT=SSTK/TOMS;
      if IDISS == 0
          RAT=1.0;
      end    
      
%
%  ***  Initial estimate of minimum pressure   ***
%
    RS0=RP;
    TV1=T(1).*(1.+R(1)/0.622)./(1.+R(1));
	TVAV=0.5.*(TV1+SSTK.*(1.+RS0./0.622)/(1.+RS0));
	CAT=CAPEM-CAPEA+0.5.*CKCD.*RAT.*(CAPEMS-CAPEM);
	CAT=max(CAT,0.0);
	PNEW=PSL.*exp(-CAT./(287.04.*TVAV));
%
%   ***  Test for convergence   ***
%
	PMOLD=PM;
	PM=PNEW;
	NP=NP+1;
    if NP > 200  || PM < 400
	  PMIN=PSL;
	  VMAX=0;
	  IFL=0;
	  return
    end
%
      end
%      
	 CATFAC=0.5.*(1.+1./b);
	 CAT=CAPEM-CAPEA+CKCD.*RAT.*CATFAC.*(CAPEMS-CAPEM);
	 CAT=max(CAT,0.0);
	 PMIN=PSL.*exp(-CAT./(287.04.*TVAV));
%
     FAC=max(0.0,(CAPEMS-CAPEM));
	 VMAX=VREDUC.*sqrt(CKCD.*RAT.*FAC);
%
%  
     function [CAPED,TOB,LNB,IFLAG]= cape(TP,RP,PP,T,R,P,SIG)
%
%     This function calculates the CAPE of a parcel with pressure PP (mb), 
%       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
%       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
%       of pressure (P in mb). CAPED is
%       the calculated value of CAPE and TOB is the temperature at the
%       level of neutral buoyancy.  IFLAG is a flag
%       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
%       not run owing to improper sounding (e.g.no water vapor at parcel level).
%       IFLAG=2 indicates that routine did not converge.                 
%-------------------------------------------------------------------------
      ptop=50;   %  Pressure below which sounding is ignored
%------------------------------------------------------------------------      
      Nold=max(size(P));
      N=1;
      for i=Nold:-1:1,
          if P(i) > ptop
              N=max(N,i);
              break
          end    
      end   
      if N < Nold
          P(N+1:Nold)=[];
          T(N+1:Nold)=[];
          R(N+1:Nold)=[];
      end    
      TVRDIF=zeros(1,N);
%
%   ***   Default values   ***
%      
      CAPED=0.0;
      TOB=T(1);
      IFLAG=1;
%
%   ***   Check that sounding is suitable    ***
%
      if RP < 1e-6 || TP < 200
       IFLAG=0;
       return
      end            
%
%   ***   Assign values of thermodynamic constants     ***
%
      CPD=1005.7;
      CPV=1870.0;
%     CL=4190.0;
      CL=2500.0;
      CPVMCL=CPV-CL;
      RV=461.5;
      RD=287.04;
      EPS=RD./RV;
      ALV0=2.501e6;
%
%   ***  Define various parcel quantities, including reversible   ***
%   ***                       entropy, S.                         ***
%                           
      TPC=TP-273.15;
      ESP=6.112*exp(17.67.*TPC/(243.5+TPC));
      EVP=RP*PP/(EPS+RP);
      RH=EVP/ESP;
      RH=min(RH,1.0);
      ALV=ALV0+CPVMCL*TPC;
      S=(CPD+RP*CL)*log(TP)-RD*log(PP-EVP)+...
         ALV*RP./TP-RP*RV*log(RH);            
%
%   ***  Find lifted condensation pressure, PLCL   ***
%     
	CHI=TP/(1669.0-122.0*RH-TP);
	PLCL=PP*(RH^CHI);
%
%   ***  Begin updraft loop   ***
%
	NCMAX=0;
%
	JMIN=1e6;
%    
    for J=1:N,
%
	JMIN=min(JMIN,J);
%
%    ***  Parcel quantities below lifted condensation level   ***
%	 
     if P(J) >= PLCL
	  TG=TP*(P(J)./PP)^(RD/CPD);
	  RG=RP;
%
%   ***   Calculate buoyancy   ***
%  
	  TLVR=TG*(1.+RG/EPS)./(1.+RG);
	  TVRDIF(J)=TLVR-T(J).*(1.+R(J)/EPS)/(1+R(J));
     else
%
%   ***  Parcel quantities above lifted condensation level  ***
%	 
	  TGNEW=T(J);          
	  TJC=T(J)-273.15; 
	  ES=6.112*exp(17.67*TJC/(243.5+TJC));
	  RG=EPS*ES/(P(J)-ES);
%
%   ***  Iteratively calculate lifted parcel temperature and mixing   ***
%   ***                ratio for reversible ascent                    ***
%
	  NC=0;
      TG=0.0;
%    
      while (abs(TGNEW-TG)) > 0.001
%
	   TG=TGNEW;
	   TC=TG-273.15;
	   ENEW=6.112*exp(17.67*TC./(243.5+TC));
       RG=EPS*ENEW/(P(J)-ENEW);   
%          
	  NC=NC+1;
%
%   ***  Calculate estimates of the rates of change of the entropy    ***
%   ***           with temperature at constant pressure               ***
%  
	  ALV=ALV0+CPVMCL*(TG-273.15);
	  SL=(CPD+RP*CL+ALV*ALV*RG./(RV*TG*TG))/TG;
	  EM=RG*P(J)/(EPS+RG);
	  SG=(CPD+RP*CL)*log(TG)-RD*log(P(J)-EM)+ ...
          ALV*RG/TG;
          if NC < 3
	   AP=0.3;
          else
	   AP=1.0;
          end
	  TGNEW=TG+AP*(S-SG)/SL;  
%
%   ***   Bail out if things get out of hand   ***
%
       if NC > 500 || ENEW > (P(J)-1)
            IFLAG=2;
            return
       end
%       
       end
%       
       NCMAX=max(NC,NCMAX);
%
%   *** Calculate buoyancy   ***
%
      RMEAN=SIG*RG+(1-SIG)*RP;
	  TLVR=TG*(1.+RG/EPS)/(1.+RMEAN);
	  TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J));
     end
    end
%
%  ***  Begin loop to find NA, PA, and CAPE from reversible ascent ***
%
	NA=0.0;
	PA=0.0;
%
%   ***  Find maximum level of positive buoyancy, INB    ***
%
	INB=1;
    for J=N:-1:JMIN;
        if TVRDIF(J) > 0
            INB=max(INB,J);
        end
    end   
    if INB == 1
	LNB=0;
        return
    end    
%
%   ***  Find positive and negative areas and CAPE  ***
%
    if INB > 1
        for J=(JMIN+1):INB
            PFAC=RD*(TVRDIF(J)+TVRDIF(J-1))*(P(J-1)-P(J))/(P(J)+P(J-1));
            PA=PA+max(PFAC,0.0);
            NA=NA-min(PFAC,0.0);
        end    
  
%   ***   Find area between parcel pressure and first level above it ***
%
        PMA=(PP+P(JMIN)) ;
        PFAC=RD*(PP-P(JMIN))/PMA;
        PA=PA+PFAC*max(TVRDIF(JMIN),0.0);
        NA=NA-PFAC*min(TVRDIF(JMIN),0.0);
%
%   ***   Find residual positive area above INB and TO  ***
%
       PAT=0.0;
       TOB=T(INB);
       LNB=P(INB);
       if INB < N
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/ ...
         (TVRDIF(INB)-TVRDIF(INB+1));
		LNB=PINB;
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB);
	    TOB=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/ ...
         (P(INB)-P(INB+1));
       end

%
%   ***   Find CAPE  ***
%            
	 CAPED=PA+PAT-NA;
	 CAPED=max(CAPED,0.0);
    end
%
	return
