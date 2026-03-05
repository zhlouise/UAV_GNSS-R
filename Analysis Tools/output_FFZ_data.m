% ------------------------------------------------------------------------
% Code to output data relating to a First Fresnel zone as a table:
% Col 1: GPS Time of the Week (aligned to the flight time from flight log)
% Col 2: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
% Col 3: Centroid Latitude
% Col 4: Centroid Longitude
% Col 5: Ellipse Azimuth (in degrees)
% Col 6: Ellipse Semi-major Axis (in meters)
% Col 7: Ellipse Semi-minor Axis (in meters)
% ------------------------------------------------------------------------

clc;
clear;
close all;

% Prompt for input of processed dual antenna data files
[filename, path] = uigetfile('.mat','Select the processed dual antenna .mat file.');
load([path, filename]);
zenith_data = dual_antenna_data{1,1};
nadir_data = dual_antenna_data{1,2};

% Prompt for input of Pixhawk ulog file
[filename, path] = uigetfile('.ulg','Select the .ulg file recorded on Pixhawk flight controller.');
ulogOBJ = ulogreader([path, filename]);
flight_log = readTopicMsgs(ulogOBJ);

% Sync zenith_data and nadir_data with flight_log timestamps
[zenith_data_new,nadir_data_new] = sync_dual_antenna_data(zenith_data, nadir_data, flight_log);

% Extract measurement data
[gps_time, prn_zenith, ele_angle_zenith, azi_angle_zenith, CNR_zenith, pr_resi_zenith, ...
    carrier_ph_zenith, pseudorange_zenith, dopp_zenith] = extract_measurements(zenith_data_new);
[~, prn_nadir, ele_angle_nadir, azi_angle_nadir, CNR_nadir, pr_resi_nadir, ...
    carrier_ph_nadir, pseudorange_nadir, dopp_nadir] = extract_measurements(nadir_data_new);

% Extract ground truth
ground_truth = extract_ground_truth(gps_time, flight_log);

% Calculate the above ground altitude
altitude_ag = cal_altitude_above_ground_HK(ground_truth,'../example_data/Whole_HK_DTM_5m.asc');

% Find the indexes of prns that exist in both zenith_data and nadir_data
[~,id_zenith, id_nadir] = intersect(prn_zenith, prn_nadir);

% Obtain the First Fresnel Zones of the reflected signals
[centroid_lat, centroid_lon, a, b] = cal_fresnel_zones(ground_truth, altitude_ag, ...
    ele_angle_zenith(:,id_zenith), azi_angle_zenith(:,id_zenith));

% Prepare table data (all matrices: time x PRN)
rep_gps_time = repmat(gps_time, 1, height(prn_zenith(id_zenith)));
rep_prn = repmat(prn_zenith(id_zenith).', height(gps_time), 1);
azi_list = azi_angle_zenith(:, id_zenith);

% Prepare table data (all matrices: time x PRN)
T = table(rep_gps_time(:), rep_prn(:), centroid_lat(:), centroid_lon(:), azi_list(:), a(:), b(:), ...
    'VariableNames', {'GPSTimeOfWeek', 'SatellitePRN', 'CentroidLat', ...
    'CentroidLon', 'EllipseAzimuth_deg', 'SemiMajorAxis_m', 'SemiMinorAxis_m'});
T = rmmissing(T);  % Removes rows with any NaN

% Prompt for saving of T
[filename, path] = uiputfile({'*.csv','CSV files (*.csv)'; ...
                               '*.xlsx','Excel files (*.xlsx)'; ...
                               '*.*','All files (*.*)'}, ...
                               'Save Fresnel Zone Table As...');
writetable(T, fullfile(path, filename));