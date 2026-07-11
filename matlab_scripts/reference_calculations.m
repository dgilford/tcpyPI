%% Description and Setup
% Generates BE02 MATLAB reference potential-intensity outputs for validating pyPI.
%
% Reads the pyPI 2024 sample (data/sample_data.nc), runs the BE02 MATLAB algorithm
% (pc_min.m) over every ocean column, and writes the reference outputs to
% data/matlab_pi_reference_2024.nc. The Python notebook verify_pi.ipynb merges this
% file with sample_data.nc and compares pyPI vs MATLAB.
%
% Run this script from the matlab_scripts/ directory (it needs pc_min.m alongside).
%
% NOTE on dimension order: MATLAB's ncread returns netCDF variables with their
% dimension order reversed, so sample_data.nc variable t (month,p,lat,lon) is read
% here as [nlon,nlat,np,nmon]. We write the outputs with dims {lon,lat,month} so
% that xarray reads them back as (month,lat,lon), matching sample_data.nc.

close all
clear

% paths
addpath(pwd)                                    % pc_min.m lives alongside this script

% Locate the input sample. Works whether run from matlab_scripts/ (repo layout) or
% from a flat folder with all files uploaded together (e.g. MATLAB Online).
in_file = fullfile('..', 'data', 'sample_data.nc');
if ~isfile(in_file)
    in_file = 'sample_data.nc';
end
out_dir = fullfile('..', 'data');
if ~isfolder(out_dir)
    out_dir = '.';
end
out_file = fullfile(out_dir, 'matlab_pi_reference_2024.nc');

%% Load the pyPI sample (already in pyPI units: SST/T in C, MSL in hPa, r in g/kg)

p     = ncread(in_file, 'p');       % [np]                 hPa
lat   = ncread(in_file, 'lat');     % [nlat]               degrees north
lon   = ncread(in_file, 'lon');     % [nlon]               degrees east (0-360)
month = ncread(in_file, 'month');   % [nmon]               month number (1-12)
sst   = ncread(in_file, 'sst');     % [nlon,nlat,nmon]     C   (NaN over land)
msl   = ncread(in_file, 'msl');     % [nlon,nlat,nmon]     hPa
T     = ncread(in_file, 't');       % [nlon,nlat,np,nmon]  C
r     = ncread(in_file, 'r');       % [nlon,nlat,np,nmon]  g/kg (water-vapor mixing ratio)

nlon = length(lon); nlat = length(lat); nmon = length(month);

%% Calculate Potential Intensity with the BE02 algorithm (pc_min.m)
% pc_min signature: [PMIN,VMAX,TO,LNB,IFL] = pc_min(sst, msl, p, T_column, q_column)

VMAX = nan(nlon, nlat, nmon);
PMIN = nan(nlon, nlat, nmon);
TO   = nan(nlon, nlat, nmon);
LNB  = nan(nlon, nlat, nmon);
IFL  = nan(nlon, nlat, nmon);

for x = 1:nlon
    for y = 1:nlat
        for m = 1:nmon
            % SST is NaN over land (masked in the sample); skip those columns
            if (~isnan(sst(x, y, m)) && sst(x, y, m) > 0)
                [PMIN(x, y, m), VMAX(x, y, m), TO(x, y, m), LNB(x, y, m), IFL(x, y, m)] = ...
                    pc_min(squeeze(sst(x, y, m)), squeeze(msl(x, y, m)), ...
                           squeeze(p), squeeze(T(x, y, :, m)), squeeze(r(x, y, :, m)));
            end
        end
    end
end

%% Write the reference outputs (dims {lon,lat,month} -> xarray reads month,lat,lon)

if isfile(out_file)
    delete(out_file)
end

% coordinates (written so the file aligns with sample_data.nc on merge)
nccreate(out_file, 'p',     'Dimensions', {'p', np_len(p)},       'Format', 'netcdf4_classic')
ncwrite(out_file, 'p', p);   ncwriteatt(out_file, 'p', 'units', 'hPa');
nccreate(out_file, 'lat',   'Dimensions', {'lat', nlat},          'Format', 'netcdf4_classic')
ncwrite(out_file, 'lat', lat);   ncwriteatt(out_file, 'lat', 'units', 'degrees_north');
nccreate(out_file, 'lon',   'Dimensions', {'lon', nlon},          'Format', 'netcdf4_classic')
ncwrite(out_file, 'lon', lon);   ncwriteatt(out_file, 'lon', 'units', 'degrees_east');
nccreate(out_file, 'month', 'Dimensions', {'month', nmon},        'Format', 'netcdf4_classic')
ncwrite(out_file, 'month', month);

% helper to create+write a 3-D (lon,lat,month) output variable with attributes
dims3 = {'lon', nlon, 'lat', nlat, 'month', nmon};
write3(out_file, 'Vmax', VMAX, dims3, 'Maximum Potential Intensity (BE02 MATLAB)', 'm/s');
write3(out_file, 'Pmin', PMIN, dims3, 'Minimum Central Pressure (BE02 MATLAB)', 'hPa');
write3(out_file, 'To',   TO,   dims3, 'Outflow Temperature (BE02 MATLAB)', 'K');
write3(out_file, 'LNB',  LNB,  dims3, 'Level of Neutral Buoyancy (BE02 MATLAB)', 'hPa');
write3(out_file, 'PI_flag', IFL, dims3, 'BE02 MATLAB status flag', '');

ncdisp(out_file)

%% local helper functions
function n = np_len(p)
    n = length(p);
end

function write3(file, name, data, dims, standard_name, units)
    nccreate(file, name, 'Dimensions', dims, 'Format', 'netcdf4_classic')
    ncwrite(file, name, data);
    ncwriteatt(file, name, 'standard_name', standard_name);
    if ~isempty(units)
        ncwriteatt(file, name, 'units', units);
    end
end
