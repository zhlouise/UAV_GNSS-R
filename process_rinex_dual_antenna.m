% ------------------------------------------------------------------------
% Code to process measurement-level data from rinex obs files. Requires
% dual antenna data input.
%
% Data cells in the following structure:
% Col-1: GPS Time
% Col-2: Positioning Solution (Ordinary Least Squares - corrected by rtklib)
% Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
% Col-4: Satellite Elevation Angle
% Col-5: Satellite Azimuth Angle (NED, North as zero)
% Col-6,7,8: Full-ephemeris PRN, El, Az
% Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
% Col-10: L1 PRN, Pseudorange with receiver_bias correction, SV-XYZ-Position in ECEF for positioning
% Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
% ------------------------------------------------------------------------

clc;
clear;
close all;
c = gt.C.CLIGHT; % Speed of light

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

% Load rtklib config file for spp
opt = gt.Gopt("MatRTKLIB-main/examples/data/kinematic/spp.conf");

% Full SV data from Ephemeris
obsf_GE = gt.Gobs('MatRTKLIB-main/examples/data/pseudo_GEfull.obs');
obsf_BJ = gt.Gobs('MatRTKLIB-main/examples/data/pseudo_BJfull.obs');

%% Observation processing - zenith antenna
zenith_data = cell(obs_zenith.n, 11);
% Calculate position and pre-process residuals
[sol_zenith, stat_zenith] = gt.Gfun.pntpos(obs_zenith, nav, opt); % Positioning calculation using rtklib
sat_zenith.setRcvPos(sol_zenith.pos); % Reset the estimated position to be the solved position (for more accurate trop+ion correction)
obs_zenith = obs_zenith.residuals(sat_zenith); % Compute the pre-fit residuals (with sat_bias+trop+ion removed)
% Retrieve time group delay (TGD)
TGD = nav.getTGD(sat_zenith.sat);
% Get the abosolute receiver clock bias instead of the reciver clock bias
% relative to GPS
dtr = [sol_zenith.dtr(:,1),sol_zenith.dtr(:,2:end)+sol_zenith.dtr(:,1)];

% Data processing loop
for idt = 1:obs_zenith.n

    % Extraction of raw pseudorange (without corrections)
    rawP = obs_zenith.L1.P(idt,:);
    sv_zenith = [obs_zenith.sat', rawP', sat_zenith.x(idt,:)', sat_zenith.y(idt,:)', sat_zenith.z(idt,:)'];
    % Some SV not avaliable in Eph/pseudorange, excluded
    exclusion_bool = any([isnan(sv_zenith(:,2)), isnan(sv_zenith(:,3))],2);
    sv_zenith(exclusion_bool,:) = [];

    % Populate zenith_data
    % Col-1: GPS Time
    zenith_data{idt,1} = obs_zenith.time.tow(idt);

    % Col-2: Positioning Solution (Ordinary Least Squares)
    pos_est_zenith = sol_zenith.pos.llh(idt,:);
    zenith_data{idt,2} = pos_est_zenith;
    
    % Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
    zenith_data{idt,3} = sv_zenith(:,1);

    % Col-4: Satellite Elevation Angle
    [~, ex, ey, ez] = rtklib.geodist(sv_zenith(:,3),sv_zenith(:,4),sv_zenith(:,5), sol_zenith.pos.xyz(idt,:));
    [az, el] = rtklib.satazel(pos_est_zenith, ex, ey, ez);
    zenith_data{idt,4} = el;

    % Col-5: Satellite Azimuth Angle (NED, North as zero)
    zenith_data{idt,5} = az; 

    % Col-6,7,8: Full-ephemeris PRN, El, Az
    obsc = obs_zenith.selectTime(idt); % Current observation
    obsf_GE.time = obsc.time;
    obsf_BJ.time = obsc.time;
    satf_GE = gt.Gsat(obsf_GE, nav);
    satf_BJ = gt.Gsat(obsf_BJ, nav);
    % Combine GE and BJ data
    satf.x = [satf_GE.x,satf_BJ.x];
    satf.y = [satf_GE.y,satf_BJ.y];
    satf.z = [satf_GE.z,satf_BJ.z];
    satf.sat = [satf_GE.sat,satf_BJ.sat];
    [~, ex_f, ey_f, ez_f] = rtklib.geodist(satf.x', satf.y', ...
        satf.z', sol_zenith.pos.xyz(idt,:));
    [az_f, el_f] = rtklib.satazel(pos_est_zenith, ex_f, ey_f, ez_f);
    temp_sv_full = [satf.sat',el_f,az_f];
    temp_sv_full = temp_sv_full(temp_sv_full(:,2)>0,:); % Filter out satellites below the horizon
    % Save into data cell
    zenith_data{idt,6} = temp_sv_full(:,1);
    zenith_data{idt,7} = temp_sv_full(:,2);
    zenith_data{idt,8} = temp_sv_full(:,3);

    % Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
    temp_sys = log2(obs_zenith.sys(~exclusion_bool));
    temp_sys(temp_sys==0) = 1; % Map the constellations to receiver bias columns (1=GPS,2=GLONASS,3=Galileo,4=QZSS,5=Beidou)
    pr_resi_zenith = obs_zenith.L1.resPc(idt,~exclusion_bool) - TGD(~exclusion_bool) - ...
        dtr(idt,temp_sys)*c; % Pre-fit residual - TGD - receiver_bias
    zenith_data{idt,9} = [obs_zenith.L1.S(idt,~exclusion_bool)',pr_resi_zenith',...
        obs_zenith.L1.L(idt,~exclusion_bool)'];

    % Col-10: L1 PRN, Pseudorange with receiver_bias correction, SV-XYZ-Position in ECEF for positioning
    zenith_data{idt,10} = [sv_zenith(:,1),sv_zenith(:,2)-dtr(idt,temp_sys)*c,sv_zenith(:,3:5)]; % Include receiver bias correction for pseudorange.
    
    % Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
    zenith_data{idt,11} = [sv_zenith(:,1), obs_zenith.L1.D(idt,~exclusion_bool)',sat_zenith.vx(idt,~exclusion_bool)'...
        ,sat_zenith.vy(idt,~exclusion_bool)',sat_zenith.vz(idt,~exclusion_bool)'];

    fprintf("Processing zenith data => %d/%d\n",idt,obs_zenith.n);
end

% Plot on map
figure();
temp_pos_llh = cell2mat(zenith_data(:,2));
geoplot(temp_pos_llh(:,1), temp_pos_llh(:,2), 'b.', 'MarkerSize', 5);

%% Observation processing - nadir antenna
nadir_data = cell(obs_nadir.n, 11);
% Calculate position and pre-process residuals
[sol_nadir, stat_nadir] = gt.Gfun.pntpos(obs_nadir, nav, opt); % Positioning calculation using rtklib
sat_nadir.setRcvPos(sol_nadir.pos); % Reset the estimated position to be the solved position (for more accurate trop+ion correction)
obs_nadir = obs_nadir.residuals(sat_nadir); % Compute the pre-fit residuals (with sat_bias+trop+ion removed)
% Retrieve time group delay (TGD)
TGD = nav.getTGD(sat_nadir.sat);
% Get the abosolute receiver clock bias instead of the reciver clock bias
% relative to GPS
dtr = [sol_nadir.dtr(:,1),sol_nadir.dtr(:,2:end)+sol_nadir.dtr(:,1)];

% Data processing loop
for idt = 1:obs_nadir.n

    % Extraction of raw pseudorange (without corrections)
    rawP = obs_nadir.L1.P(idt,:);
    sv_nadir = [obs_nadir.sat', rawP', sat_nadir.x(idt,:)', sat_nadir.y(idt,:)', sat_nadir.z(idt,:)'];
    % Some SV not avaliable in Eph/pseudorange, excluded
    exclusion_bool = any([isnan(sv_nadir(:,2)), isnan(sv_nadir(:,3))],2);
    sv_nadir(exclusion_bool,:) = [];

    % Populate nadir_data
    % Col-1: GPS Time
    nadir_data{idt,1} = obs_nadir.time.tow(idt);

    % Col-2: Positioning Solution (Ordinary Least Squares)
    pos_est_nadir = sol_nadir.pos.llh(idt,:);
    nadir_data{idt,2} = pos_est_nadir;
    
    % Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
    nadir_data{idt,3} = sv_nadir(:,1);

    % Col-4: Satellite Elevation Angle
    [~, ex, ey, ez] = rtklib.geodist(sv_nadir(:,3),sv_nadir(:,4),sv_nadir(:,5), sol_nadir.pos.xyz(idt,:));
    [az, el] = rtklib.satazel(pos_est_nadir, ex, ey, ez);
    nadir_data{idt,4} = el;

    % Col-5: Satellite Azimuth Angle (NED, North as zero)
    nadir_data{idt,5} = az; 

    % Col-6,7,8: Full-ephemeris PRN, El, Az
    obsc = obs_nadir.selectTime(idt); % Current observation
    obsf_GE.time = obsc.time;
    obsf_BJ.time = obsc.time;
    satf_GE = gt.Gsat(obsf_GE, nav);
    satf_BJ = gt.Gsat(obsf_BJ, nav);
    % Combine GE and BJ data
    satf.x = [satf_GE.x,satf_BJ.x];
    satf.y = [satf_GE.y,satf_BJ.y];
    satf.z = [satf_GE.z,satf_BJ.z];
    satf.sat = [satf_GE.sat,satf_BJ.sat];
    [~, ex_f, ey_f, ez_f] = rtklib.geodist(satf.x', satf.y', ...
        satf.z', sol_nadir.pos.xyz(idt,:));
    [az_f, el_f] = rtklib.satazel(pos_est_nadir, ex_f, ey_f, ez_f);
    temp_sv_full = [satf.sat',el_f,az_f];
    temp_sv_full = temp_sv_full(temp_sv_full(:,2)>0,:); % Filter out satellites below the horizon
    % Save into data cell
    nadir_data{idt,6} = temp_sv_full(:,1);
    nadir_data{idt,7} = temp_sv_full(:,2);
    nadir_data{idt,8} = temp_sv_full(:,3);

    % Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
    temp_sys = log2(obs_nadir.sys(~exclusion_bool));
    temp_sys(temp_sys==0) = 1; % Map the constellations to receiver bias columns (1=GPS,2=GLONASS,3=Galileo,4=QZSS,5=Beidou)
    pr_resi_nadir = obs_nadir.L1.resPc(idt,~exclusion_bool) - TGD(~exclusion_bool) - ...
        dtr(idt,temp_sys)*c; % Pre-fit residual - TGD - receiver_bias
    nadir_data{idt,9} = [obs_nadir.L1.S(idt,~exclusion_bool)',pr_resi_nadir',...
        obs_nadir.L1.L(idt,~exclusion_bool)'];

    % Col-10: L1 PRN, Pseudorange with receiver_bias correction, SV-XYZ-Position in ECEF for positioning
    nadir_data{idt,10} = [sv_nadir(:,1),sv_nadir(:,2)-dtr(idt,temp_sys)*c,sv_nadir(:,3:5)]; % Include receiver bias correction for pseudorange.;
    
    % Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
    nadir_data{idt,11} = [sv_nadir(:,1), obs_nadir.L1.D(idt,~exclusion_bool)',sat_nadir.vx(idt,~exclusion_bool)'...
        ,sat_nadir.vy(idt,~exclusion_bool)',sat_nadir.vz(idt,~exclusion_bool)'];

    fprintf("Processing nadir data => %d/%d\n",idt,obs_nadir.n);
end

% Plot on map
figure();
temp_pos_llh = cell2mat(nadir_data(:,2));
geoplot(temp_pos_llh(:,1), temp_pos_llh(:,2), 'b.', 'MarkerSize', 5);

%% Save results
dual_antenna_data = {zenith_data, nadir_data};

[processed_filename, processed_pathname] = uiputfile('.mat', 'Save processed dual-antenna data.');

if isequal(processed_filename, 0) || isequal(processed_pathname, 0)
    disp('Save operation cancelled.');
else
    % Create the full file path
    full_path = fullfile(processed_pathname, processed_filename);
    
    % Perform the save operation
    save(full_path, 'dual_antenna_data');
    disp(['Data saved to: ', full_path]);
end