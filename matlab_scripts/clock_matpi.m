% This code is designed to mimic the behavior of the clock_pypi.ipynb
% notebook, in order to compare the timings between the MATLAB PI calculations
% and the Python calculations performed by pyPI

% It should be noted that the MATLAB PI function (pc_min.m) is not
% vectorized. This results in the need for "for loops" which calculate PI
% over each sample profile, greatly increasing the time of execution. In 
% contrast, pyPI may be used in combination with xarray.apply_ufunc to
% calculte PI in a vectorized manner. However, this does not necessarily
% imply the pyPI calculation is faster.

%% Set up
close all
clear

%% Load the Data

dat_loc='../data/sample_data.nc';

% load the grids
p=ncread(dat_loc,'p');
lat=ncread(dat_loc,'lat');
lon=ncread(dat_loc,'lon');
month=ncread(dat_loc,'month');
% get the number of locations
nlat=length(lat);
nlon=length(lon);
nmon=length(month);

% load the base data we will randomly sample from
t=ncread(dat_loc,'t'); % in C
q=ncread(dat_loc,'q'); % in g/kg
msl=ncread(dat_loc,'msl'); % in hPa
sst=ncread(dat_loc,'sst'); % in C

%% Develop Sample Parameters

rng(2) % fix the random seed
% copy the sample sizes from the Python code
nsamps=[     1,    216,   1331,   4096,   9261,  17576,  29791,  46656, ...
        68921,  97336, 132651, 175616, 226981, 287496, 357911, 438976, ...
       531441, 636056, 753571, 884736];

%% Time the Runs

% create the array to store the time elapsed
nruns=10;
time_elapsed=nan(nruns,length(nsamps));

% get the starting time
tstart=tic;

% loop over the array of number of samples
for j=1:length(nsamps)
    
    % get random samples
    rlat=randsample(nlat,nsamps(j),true);
    rlon=randsample(nlon,nsamps(j),true);
    rmon=randsample(nmon,nsamps(j),true);
    
    % loop over the number of runs
    for i=1:nruns
        % allocate the data storage arrays
        PMIN=nan(nsamps(j),1);
        VMAX=nan(nsamps(j),1);
        TO=nan(nsamps(j),1);
        LNB=nan(nsamps(j),1);
        IFL=nan(nsamps(j),1);
        
        % start the timer
        tic
        % calculate potential intensity
        for s=1:nsamps(j)
            [PMIN(s,1),VMAX(s,1),TO(s,1),LNB(s,1),IFL(s,1)] = ...
                pc_min(squeeze(sst(rlon(s),rlat(s),rmon(s))), ...
                                squeeze(msl(rlon(s),rlat(s),rmon(s))), ...
                                p, ...
                                squeeze(t(rlon(s),rlat(s),:,rmon(s))), ...
                                squeeze(q(rlon(s),rlat(s),:,rmon(s))));
        end
        clear PMIN VMAX TO LNB IFL
                            
        % store how much time has elapsed
        time_elapsed(i,j)=toc;
    
    end
    clear rlat rlon rmon
    % print to screen how many iterations we have done
    fprintf(strcat(int2str(j),' Time Elapsed=',num2str(toc(tstart)),'\n'))
end

%% Plot the results

% fit a curve and calculate the correlation
z = polyfit(nsamps,mean(time_elapsed,1), 1);
x = linspace(0, 1e6, 1000);
y = polyval(z,x);
Rsquared=corrcoef(nsamps,mean(time_elapsed,1)).^2;
% print the correlation to screen
fprintf(strcat('R-squared=',num2str(Rsquared(2,1)),'\n'))

% plot the result
h1=figure(1);
hold on
    errorbar(nsamps,mean(time_elapsed,1),2.*std(time_elapsed,1),'k')
    plot(nsamps,mean(time_elapsed,1),'r','LineWidth',4)
    plot(x,y,'--k')
    ylabel('Number of Samples')
    xlabel('Time per BE02 MATLAB run (s)')
    grid
    xticks(0:1e5:1e6)
    xticklabels({'0','100k','200k','300k','400k','500k','600k','700k','800k','900k','1Mil.'})
    legend('2-sigma','Mean Runtime',strcat('Fit=',num2str(round(z(1)*100000,2)),'s per 100k runs'), ...
        'Location','Northwest')
hold off
