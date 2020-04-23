%% Description and Setup
% This MATLAB script calculates Potential Intensity using the BE02 MATLAB
% algorithm and outputs the data into "sample_data.nc" for comparison with
% the Python PI algorithm (pi.py). The sample data are average monthly values
% from MERRA2 in 2004.

% setup
close all
clear
% add paths to run the code
path=pwd;
addpath(strcat(path,'/..'))
addpath(strcat(path,'/../data/'))
addpath(strcat(path,'/../data/gse17_mat_files/'))

%% Load Data and Prep for PI calculation

% define the data filenames
slp_file='merra2_slp_2004_monthly_2.5x2.5.mat';
T_file='merra2_T_2004_monthly_2.5x2.5.mat';
TS_file='merra2_TS_2004_monthly_2.5x2.5.mat';
q_file='merra2_q_2004_monthly_2.5x2.5.mat';

% load the land-sea mask
lsm_shift=ncread('lsmask_2.5x2.5.nc','lsm');
lsm_lon=ncread('lsmask_2.5x2.5.nc','longitude');
% shift the land-sea mask
lsm_lon=lsm_lon-180;
minlon=find(lsm_lon==min(abs(lsm_lon)));
lsm=circshift(lsm_shift,[-minlon+1 0]);
clear minlon lsm_lon

% load data
slp_dat=load(slp_file);     % in Pa
T_dat=load(T_file);         % in K
TS_dat=load(TS_file);       % in K
q_dat=load(q_file);         % in kg/kg
clear slp_file T_file TS_file q_file

% define grids (pressure in hPa)
lat=slp_dat.latgrid; lon=slp_dat.longrid; p=T_dat.lvlgrid;
nlon=length(lon); nlat=length(lat);

% convert to the appropriate units and store
q=q_dat.q.*1e3;                     % kg/kg	-->	g/kg	(approximating q == r/(1+r) ~~ r)
sst=squeeze(TS_dat.TS)-273.15;		% K	-->	C	(skin temperature -->	sst)
T=T_dat.T-273.15;                   % K	-->	C
msl=slp_dat.slp./100;               % Pa -->	hPa

% choose a month to calculate TC PI, subset the data
% mon=9;
% q=squeeze(q(:,:,:,mon));
% sst=squeeze(sst(:,:,:,mon));
% T=squeeze(T(:,:,:,mon));
% msl=squeeze(msl(:,:,mon));

%% Calculate Potential Intensity with the BE02 algorithm (pc_min.m)

% inputs:
% sst 
% msl
% T
% q
% 
% outputs:
% VMAX
% PMIN
% TO
% LNB
% IFL

% create the output arrays
TO=nan(nlon,nlat,12);
LNB=nan(nlon,nlat,12);
VMAX=nan(nlon,nlat,12);
IFL=nan(nlon,nlat,12);
PMIN=nan(nlon,nlat,12);

% calculate TC PI for the average month
for x=1:nlon
	for y=1:nlat
        for m=1:12
            if (squeeze(sst(x,y,m)) > 0)
                [PMIN(x,y,m),VMAX(x,y,m),TO(x,y,m),LNB(x,y,m),IFL(x,y,m)]= ...
                    pc_min(squeeze(sst(x,y,m)),squeeze(msl(x,y,m)), ...
                    squeeze(p),squeeze(T(x,y,:,m)),squeeze(q(x,y,:,m)));
            end
		end
    end
end

%% Test speed of PI calculation
xi=1;
yi=38;
mi=1;
for k=1:100
    t = cputime;
    [check,check2,check3,check4]=pc_min(squeeze(sst(xi,yi,mi)),squeeze(msl(xi,yi,mi)), ...
                    squeeze(p),squeeze(T(xi,yi,:,mi)),squeeze(q(xi,yi,:,mi)));
    e = cputime-t;
    time_elapsed(k,1)=e;
end

nanmean(time_elapsed)
hist(time_elapsed)

%% Apply Land-sea mask

% figure(2)
% contourf(lon,lat,lsm')
% colorbar

for x=1:nlon
	for y=1:nlat
        if (lsm(x,y)==1)
            % for output variables
            PMIN(x,y,:)=NaN;
            VMAX(x,y,:)=NaN;
            TO(x,y,:)=NaN;
            LNB(x,y,:)=NaN;
            IFL(x,y,:)=NaN;
            % for input variables
            sst(x,y,:)=NaN;
            msl(x,y,:)=NaN;
            T(x,y,:,:)=NaN;
            q(x,y,:,:)=NaN;
        end
    end
end

%% Plot the VMAX_average

figure(1)
contourf(lon,lat,nanmean(VMAX,3)')
hold on
    colorbar
    
%% save out the variables for use in pyPI

% choose name to save out, delete if it already exists
if isfile('../data/sample_data.nc')
     % File exists.
     delete '../data/sample_data.nc'
     nc_savepath='../data/sample_data.nc';
end

% Save the Data

% dimensions
nccreate(nc_savepath,'p',...
          'Dimensions',{'p',length(p)},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'p', 'standard_name', 'Atmospheric Pressure');
ncwriteatt(nc_savepath, 'p', 'units', 'hPa');
ncwrite(nc_savepath,'p',p)

nccreate(nc_savepath,'lat',...
          'Dimensions',{'lat',length(lat)},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'lat', 'standard_name', 'Latitude');
ncwriteatt(nc_savepath, 'lat', 'units', 'degrees');
ncwrite(nc_savepath,'lat',lat)

nccreate(nc_savepath,'lon',...
          'Dimensions',{'lon',length(lon)},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'lon', 'standard_name', 'Longitude');
ncwriteatt(nc_savepath, 'lon', 'units', 'degrees');
ncwrite(nc_savepath,'lon',lon)

month=1:12;
nccreate(nc_savepath,'month',...
          'Dimensions',{'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'month', 'standard_name', 'Month');
ncwriteatt(nc_savepath, 'month', 'units', 'Month Number');
ncwrite(nc_savepath,'month',month)

% inputs
nccreate(nc_savepath,'lsm',...
          'Dimensions',{'lon',length(lon),'lat',length(lat)},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'lsm', 'standard_name', 'ERA-I Land-sea Mask');
ncwriteatt(nc_savepath, 'lsm', 'units', '0=Ocean, 1=Land');
ncwrite(nc_savepath,'lsm',lsm)

nccreate(nc_savepath,'sst',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'sst', 'standard_name', 'Sea Surface Temperature');
ncwriteatt(nc_savepath, 'sst', 'units', 'degrees C');
ncwrite(nc_savepath,'sst',sst)

nccreate(nc_savepath,'msl',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'msl', 'standard_name', 'Mean Sea Level Pressure');
ncwriteatt(nc_savepath, 'msl', 'units', 'hPa');
ncwrite(nc_savepath,'msl',msl)

nccreate(nc_savepath,'t',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'p',length(p),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 't', 'standard_name', 'Atmospheric Temperature');
ncwriteatt(nc_savepath, 't', 'units', 'degrees C');
ncwrite(nc_savepath,'t',T)

nccreate(nc_savepath,'q',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'p',length(p),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'q', 'standard_name', 'Specific Humidity');
ncwriteatt(nc_savepath, 'q', 'units', 'g/kg');
ncwrite(nc_savepath,'q',q)

% outputs
nccreate(nc_savepath,'Vmax',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'Vmax', 'standard_name', 'Maximum Potential Intensity');
ncwriteatt(nc_savepath, 'Vmax', 'units', 'm/s');
ncwrite(nc_savepath,'Vmax',VMAX)

nccreate(nc_savepath,'To',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'To', 'standard_name', 'Outflow Temperature');
ncwriteatt(nc_savepath, 'To', 'units', 'kelvin');
ncwrite(nc_savepath,'To',TO)

nccreate(nc_savepath,'Pmin',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'Pmin', 'standard_name', 'Minimum Central Pressure');
ncwriteatt(nc_savepath, 'Pmin', 'units', 'hPa');
ncwrite(nc_savepath,'Pmin',PMIN)

nccreate(nc_savepath,'LNB',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'LNB', 'standard_name', 'Level of Neutral Bouyancy for a Parcel Reversibly Lifted from Sea Level');
ncwriteatt(nc_savepath, 'LNB', 'units', 'hPa');
ncwrite(nc_savepath,'LNB',LNB)

nccreate(nc_savepath,'PI_flag',...
          'Dimensions',{'lon',length(lon),'lat',length(lat),'month',12},...
          'Format','netcdf4_classic')
ncwriteatt(nc_savepath, 'PI_flag', 'standard_name', 'Flag for BE02 algorithm');
ncwrite(nc_savepath,'PI_flag',IFL)

% show the file that was written
ncdisp(nc_savepath)
