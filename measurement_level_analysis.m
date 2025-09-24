% ------------------------------------------------------------------------
% Code to analyse measurement-level data processed from rinex obs files.
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

    pr_difference{epoch} = [temp_prn_list, nadir_data{idn,10}(ismember(nadir_data{idn,3},temp_prn_list),2) - ...
        zenith_data{epoch,10}(ismember(zenith_data{epoch,3},temp_prn_list),2)];

end

