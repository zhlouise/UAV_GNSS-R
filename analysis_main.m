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

% Prompt for input of processed dual antenna files
[filename, path] = uigetfile('.mat','Select the processed dual antenna .mat file.');
load([path, filename]);
zenith_data = dual_antenna_data{1,1};
nadir_data = dual_antenna_data{1,2};

% Sync zenith_data and nadir_data
[zenith_data_new,nadir_data_new] = sync_dual_antenna_data(zenith_data, nadir_data);

% Extract measurement data
[gps_time, prn_zenith, ele_angle_zenith, azi_angle_zenith, CNR_zenith, pr_resi_zenith, ...
    carrier_ph_zenith, pseudorange_zenith, dopp_zenith] = extract_measurements(zenith_data_new);
[~, prn_nadir, ele_angle_nadir, azi_angle_nadir, CNR_nadir, pr_resi_nadir, ...
    carrier_ph_nadir, pseudorange_nadir, dopp_nadir] = extract_measurements(nadir_data_new);

%% Percentage of SV received

sv_ratio_z = NaN(height(zenith_data),1); % Zenith data
for epoch = 1:height(zenith_data)
    sv_ratio_z(epoch) = length(zenith_data{epoch,3})/length(zenith_data{epoch,6});
end

sv_ratio_n = NaN(height(nadir_data),1); % Nadir data
for epoch = 1:height(nadir_data)
    sv_ratio_n(epoch) = length(nadir_data{epoch,3})/length(nadir_data{epoch,6});
end

% Plot SV reception ratio distributions
figure();
grid on;
hold on;
histogram(sv_ratio_z,'Normalization','pdf','BinWidth',0.01);
histogram(sv_ratio_n,'Normalization','pdf','BinWidth',0.01);
hold off;
xlabel('Satellite Reception Ratio');
ylabel('Probability Distribution');
legend('Zenith Data', 'Nadir Data');

%% Pseudorange difference between zenith and nadir (nadir - zenith)

pr_difference = cell(height(zenith_data),1);

% Use the timestamps on zenith_data as a standard
for epoch = 1:height(zenith_data)
    
    temp_time = round(zenith_data{epoch,1});
    % Align nadir_data timestamp to the zenith_data timestamp
    idn = find(round(cell2mat(nadir_data(:,1)))==temp_time);
    if isempty(idn)
        % No nadir data in this timestamp
        continue;
    end

    % Only use satellites that are received by both the zenith and nadir
    % antenna
    temp_prn_list = intersect(zenith_data{epoch,3}, nadir_data{idn,3});
    
    % Column 1: prn
    % Column 2: pseudorange difference
    % Column 3: elevation angle (zenith as standard)
    % Column 4: azimuth angle (zenith as standard)
    % Column 5: CNR ratio (nadir/zenith)
    pr_difference{epoch} = [temp_prn_list, ...
        nadir_data{idn,10}(ismember(nadir_data{idn,3},temp_prn_list),2)-zenith_data{epoch,10}(ismember(zenith_data{epoch,3},temp_prn_list),2), ...
        zenith_data{epoch,4}(ismember(zenith_data{epoch,3},temp_prn_list)), ...
        zenith_data{epoch,5}(ismember(zenith_data{epoch,3},temp_prn_list)), ...
        nadir_data{idn,9}(ismember(nadir_data{idn,3},temp_prn_list),1)./zenith_data{epoch,9}(ismember(zenith_data{epoch,3},temp_prn_list),1)];

end

% Compile pr_difference into matrix form (remove epoch dimension)
pr_difference_mat = vertcat(pr_difference{:,1});

% Plot pseudorange difference vs elevation angle
figure();
grid on;
scatter(pr_difference_mat(:,3), pr_difference_mat(:,2),'filled');
xlabel('Elevation Angle (From Zenith Data)');
ylabel('Difference in Nadir and Zenith Pseudorange (m)');

% Plot pseudorange difference vs azimuth (on polar coordinates)
figure();
grid on;
polarscatter(deg2rad(pr_difference_mat(:,4)), pr_difference_mat(:,2),'filled');
hold on;
% Configure the polar axis
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.RAxisLocation = 0;
% Configure figure and axis titles
thetaticks([0 30 60 90 120 150 180 210 240 270 300 330]);
thetaticklabels({'North','30','60','East','120','150','South','210','240','West','300','330'});
title('Difference in Nadir and Zenith Pseudorange (m)');
hold off;

% Plot pseudorange difference vs CNR ratio
figure();
grid on;
scatter(pr_difference_mat(:,5), pr_difference_mat(:,2),'filled');
xlabel('Difference in Zenith and Nadir CNR (dB-Hz)');
ylabel('Difference in Nadir and Zenith Pseudorange (m)');

% Plot CNR ratio distribution
figure();
grid on;
histogram(pr_difference_mat(:,5),'BinWidth',0.01)
xlim([0, 2.5])
xlabel('CNR ratio between Nadir and Zenith');
ylabel('Probability Distribution');

