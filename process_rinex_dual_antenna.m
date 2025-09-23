% ------------------------------------------------------------------------
% Code to process measurement-level data from rinex obs files. Requires
% dual antenna data input.
% ------------------------------------------------------------------------

clc;
clear;
close all;

%% Prompt user for input data files

% Zenith/upper antenna data
[zenith_filename, zenith_path] = uigetfile('.obs','Select the zenith antenna obs file.');
obs_zenith = gt.Gobs([zenith_path, zenith_filename]);

% Nadir/lower antenna data
[nadir_filename, nadir_path] = uigetfile('.obs','Select the nadir antenna obs file.');
obs_nadir = gt.Gobs([nadir_path, nadir_filename]);

% Ephemeris data
[eph_filename, eph_path] = uigetfile('.nav','Select the ephemeris file.');
nav = gt.Gnav([eph_path, eph_filename]);

%% Select constelations
obs_zenith = obs_zenith.selectSat(obs_zenith.sys==gt.C.SYS_GPS | obs_zenith.sys==gt.C.SYS_QZS | ...
    obs_zenith.sys==gt.C.SYS_CMP | obs_zenith.sys==gt.C.SYS_GAL | obs_zenith.sys==gt.C.SYS_GLO);
obs_nadir = obs_nadir.selectSat(obs_nadir.sys==gt.C.SYS_GPS | obs_nadir.sys==gt.C.SYS_QZS | ...
    obs_nadir.sys==gt.C.SYS_CMP | obs_nadir.sys==gt.C.SYS_GAL | obs_nadir.sys==gt.C.SYS_GLO);

%% Pre-processing
% Rough Hong Kong coordinates for initial guess & correction modelling
pos_ini = gt.Gpos(rtklib.llh2xyz([22.3193,114.1694,0]),"xyz"); 

sat_zenith = gt.Gsat(obs_zenith, nav);
sat_zenith.setRcvPos(pos_ini);

sat_nadir = gt.Gsat(obs_nadir, nav);
sat_nadir.setRcvPos(pos_ini);

%% Observation processing - zenith antenna
zenith_data = cell(obs_zenith.n, 11);

% Data processing loop
for idt = 1:obs_zenith.n

    % Initial positioning (uses raw pseudorange without corrections)
    rawP = obs_zenith.L1.P(idt,:);
    sv_zenith = [obs_zenith.sat', rawP', sat_zenith.x(idt,:)', sat_zenith.y(idt,:)', sat_zenith.z(idt,:)'];
    % Some SV not avaliable in Eph/pseudorange, excluded
    exclusion_bool = any([isnan(sv_zenith(:,2)), isnan(sv_zenith(:,3))],2);
    sv_zenith(exclusion_bool,:) = [];
    [pos_est_zenith, pr_resi_zenith] = least_squares_positioning(sv_zenith, pos_ini.xyz);

    % Populate zenith_data
    % Col-1: GPS Time
    zenith_data{idt,1} = obs_zenith.time.tow(idt);
    % Col-2: Positioning Solution (Ordinary Least Squares)
    zenith_data{idt,2} = rtklib.xyz2llh(pos_est_zenith);
    % Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
    zenith_data{idt,3} = sv_zenith(:,1);
    % Col-4: Satellite Elevation Angle
    [~, ex, ey, ez] = rtklib.geodist(sv_zenith(:,3),sv_zenith(:,4),sv_zenith(:,5), pos_est_zenith);
    [az, el] = rtklib.satazel(zenith_data{idt,2}, ex, ey, ez);
    zenith_data{idt,4} = el;
    % Col-5: Satellite Azimuth Angle (NED, North as zero)
    zenith_data{idt,5} = az; 
    % Col-6,7,8: Full-ephemeris PRN, El, Az (Pending to be completed)
    zenith_data{idt,6} = nan;
    zenith_data{idt,7} = nan;
    zenith_data{idt,8} = nan;
    % Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
    zenith_data{idt,9} = [obs_zenith.L1.S(idt,~exclusion_bool)',pr_resi_zenith,...
        obs_zenith.L1.L(idt,~exclusion_bool)'];
    % Col-10: L1 PRN, Raw Pseudorange, SV-XYZ-Position in ECEF for positioning
    zenith_data{idt,10} = sv_zenith;
    % Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
    zenith_data{idt,11} = [sv_zenith(:,1), obs_zenith.L1.D(idt,~exclusion_bool)',sat_zenith.vx(idt,~exclusion_bool)'...
        ,sat_zenith.vy(idt,~exclusion_bool)',sat_zenith.vz(idt,~exclusion_bool)'];

    fprintf("Processing zenith data => %d/%d\n",idt,obs_zenith.n);
end

% Plot on map
figure();
temp_pos_llh = cell2mat(zenith_data(:,2));
geoplot(temp_pos_llh(:,1), temp_pos_llh(:,2), 'b.', 'MarkerSize', 5);