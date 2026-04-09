% ------------------------------------------------------------------------
% Main code to analyse measurement-level data processed from rinex obs files.
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

% Prompt for input of processed dual antenna data files
[filename, path] = uigetfile('.mat','Select the processed dual antenna .mat file.');
load([path, filename]);
zenith_data = dual_antenna_data{1,1};
nadir_data = dual_antenna_data{1,2};

% Prompt for input of Pixhawk ulog file
[filename, path] = uigetfile('.ulg','Select the .ulg file recorded on Pixhawk flight controller.');
ulogOBJ = ulogreader([path, filename]);
flight_log = readTopicMsgs(ulogOBJ);

%% Extract and process dual antenna information

% Sync zenith_data and nadir_data with flight_log timestamps
[zenith_data_new,nadir_data_new] = sync_dual_antenna_data(zenith_data, nadir_data, flight_log);

% Extract measurement data
[gps_time, prn_zenith, ele_angle_zenith, azi_angle_zenith, CNR_zenith, pr_resi_zenith, ...
    carrier_ph_zenith, pseudorange_zenith, dopp_zenith] = extract_measurements(zenith_data_new);
[~, prn_nadir, ele_angle_nadir, azi_angle_nadir, CNR_nadir, pr_resi_nadir, ...
    carrier_ph_nadir, pseudorange_nadir, dopp_nadir] = extract_measurements(nadir_data_new);

% Extract ground truth
ground_truth = extract_ground_truth(gps_time, flight_log);

% % Calculate pseudorange errors
% zenith_pr_error = cal_SD_pr_error(zenith_data_new, ground_truth);
% nadir_pr_error = cal_SD_pr_error(nadir_data_new, ground_truth);

% Calculate the above ground altitude
altitude_ag = cal_altitude_above_ground_HK(ground_truth,'example_data/Whole_HK_DTM_5m.asc');

% Find the indexes of prns that exist in both zenith_data and nadir_data
[~,id_zenith, id_nadir] = intersect(prn_zenith, prn_nadir);

% Obtain the First Fresnel Zones of the reflected signals
[centroid_lat, centroid_lon, a, b] = cal_fresnel_zones(ground_truth, altitude_ag, ...
    ele_angle_zenith(:,id_zenith), azi_angle_zenith(:,id_zenith));

%% Satellite sky plots

validIdx = ~isnan(azi_angle_zenith) & ~isnan(ele_angle_zenith) & ele_angle_zenith > 0;
figure();
skyplot(azi_angle_zenith(validIdx), ele_angle_zenith(validIdx));
title('Satellites Received by Zenith Receiver');
data_skyplot(ele_angle_zenith, azi_angle_zenith, CNR_zenith);
title('RHCP CNR');

validIdx = ~isnan(azi_angle_nadir) & ~isnan(ele_angle_nadir) & ele_angle_nadir > 0;
figure();
skyplot(azi_angle_nadir(validIdx), ele_angle_nadir(validIdx));
title('Satellites Received by Nadir Receiver');
data_skyplot(ele_angle_nadir, azi_angle_nadir, CNR_nadir);
title('LHCP CNR');

%% Percentage of SV received

% Obtain the SV ratios for both zenith and nadir data
sv_ratio_zenith = NaN(length(gps_time),1);
sv_ratio_nadir = NaN(length(gps_time),1);
sv_ratio_n2z = NaN(length(gps_time),1); % Ratio of SV received by nadir antenna over SV received by zenith antenna
for epoch = 1:length(gps_time)
    sv_ratio_zenith(epoch) = length(zenith_data_new{epoch,3})/length(zenith_data_new{epoch,6});
    sv_ratio_nadir(epoch) = length(nadir_data_new{epoch,3})/length(nadir_data_new{epoch,6});
    sv_ratio_n2z(epoch) = length(nadir_data_new{epoch,3})/length(zenith_data_new{epoch,3});
end

% Histogram plot
figure();
grid on;
hold on;
histogram(sv_ratio_zenith,'Normalization','pdf','BinWidth',0.01);
histogram(sv_ratio_nadir,'Normalization','pdf','BinWidth',0.01);
hold off;
xlim([0.3 0.9]);
xlabel('Satellite Reception Ratio');
ylabel('Normalized Probability Distribution');
legend('Zenith RHCP Antenna', 'Nadir LHCP Antenna');

% Plot ratio of SV received by nadir antenna over SV received by zenith
% antenna on satellite imagery map
sat_img_plot(ground_truth, sv_ratio_n2z); % Use zenith SPP result as georeference coordinates
clim([0.5 2.5]);
colormap('turbo');
title('Nadir-Received SV/Zenith-Received SV');

%% CNR ratio

CNR_ratio = CNR_nadir(:,id_nadir)./CNR_zenith(:,id_zenith);
CNR_ratio_clean = CNR_ratio;
CNR_ratio_clean(altitude_ag<35,:) = NaN;
fresnel_zone_heatmap(centroid_lat,centroid_lon, a, b, azi_angle_zenith(:,id_zenith), CNR_ratio_clean);

%% CNR difference

CNR_diff = CNR_zenith(:,id_zenith) - CNR_nadir(:,id_nadir);
CNR_diff_clean = CNR_diff;
CNR_diff_clean(altitude_ag<35,:) = NaN;
fresnel_zone_heatmap(centroid_lat,centroid_lon, a, b, azi_angle_zenith(:,id_zenith), CNR_diff_clean);

%% Pseudorange error

figure();
hold on;
grid on;
histogram(zenith_pr_error(zenith_pr_error<100 & zenith_pr_error>-100 & zenith_pr_error~=0), ...
    'Normalization','pdf','BinWidth',1);
histogram(nadir_pr_error(nadir_pr_error<100 & nadir_pr_error>-100 & nadir_pr_error~=0), ...
    'Normalization','pdf','BinWidth',1);
hold off;
xlabel('Single Difference Pseudorange Error (m)');
ylabel('Normalized Probability Distribution');
legend('Zenith RHCP Antenna', 'Nadir LHCP Antenna');

%% Vegetation delay

reflection_delay = estimate_delay(altitude_ag, ele_angle_nadir(:,id_nadir));
reflection_delay_clean = reflection_delay;
reflection_delay_clean(abs(reflection_delay_clean)>150) = NaN;
pr_diff = pseudorange_nadir(:,id_nadir)-pseudorange_zenith(:,id_zenith);
pr_diff_clean = pr_diff;
pr_diff_clean(abs(pr_diff_clean)>300) = NaN;
clock_and_delay_est = pr_diff_clean-reflection_delay_clean;
est_clock_delay = mean(clock_and_delay_est, 2, 'omitnan'); % Estimated receiver clock delay for each epoch
veg_delay = clock_and_delay_est - est_clock_delay;

% Estimated vegetation delay on FFZ
veg_delay_clean = veg_delay;
veg_delay_clean(altitude_ag<40,:) = NaN;
fresnel_zone_heatmap(centroid_lat,centroid_lon, a, b, azi_angle_zenith(:,id_zenith), veg_delay_clean);
% gridded_heatmap(centroid_lat, centroid_lon, veg_delay_clean);