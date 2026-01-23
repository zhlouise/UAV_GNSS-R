% ------------------------------------------------------------------------
% Code to simulate and plot Fresnel zones from ephemeris data, given the
% desired location and time of observation. 
% ------------------------------------------------------------------------

clc;
clear;
close all;

%% Prompt user for inputs

% Ephemeris data
[eph_filename, eph_path] = uigetfile('.nav','Select the ephemeris file.');
nav = gt.Gnav([eph_path, eph_filename]);

% Desired location
location_prompt = {'Enter location latitude (decimal degrees): ', ...
    'Enter location longitude (decimal degrees): ', ...
    'Enter location ellipsoidal height (m): ', ...
    'Enter desired flight altitude above ground (m): '};
location_input = inputdlg(location_prompt, 'Input Location');
pos_llh = str2double(location_input(1:3))';
gpos = gt.Gpos(pos_llh, 'llh');
altitude_ag = str2double(location_input{4});

% Desired simulation time in UTC (HKT-8)
time_prompt = {'Year (YYYY): ', 'Month (1-12): ', 'Day: ', 'Hour (0-23): ', ...
    'Minute (0-59): ', 'Second (0-59): '};
time_input = inputdlg(time_prompt, 'Input Simulation Time in UTC');
sim_time = gt.Gtime(str2double(time_input)', 1);

%% Obtain satellite elevation and azimuth angles

sats = 1:143; 
gsat = gt.Gsat(sim_time, sats, nav, gt.C.EPHOPT_BRDC);

% Calculate azimuth and elevation angles
[~, ex, ey, ez] = rtklib.geodist(gsat.x, gsat.y, gsat.z, gpos.xyz);
[az, el] = rtklib.satazel(gpos.llh, ex, ey, ez);

% Filter out satellites with negative elevation angle (below the horizon)
visibleIdx = el > 0;
ele_visible = el(visibleIdx);
azi_visible = az(visibleIdx);

%% Obtain simulated Fresnel zones

[c_lat, c_lon, a, b] = cal_fresnel_zones(pos_llh, altitude_ag, ele_visible, azi_visible);

% Plot Fresnel zones
% fresnel_zone_plot handles the georeferenced plotting on the map
fresnel_zone_plot(c_lat, c_lon, a, b, azi_visible, ele_visible);