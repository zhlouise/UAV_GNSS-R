function [gps_time, prn, ele_angle, azi_angle, CNR, pr_resi, carrier_ph, pseudorange, dopp] = extract_measurements(data_all)
% -------------------------------------------------------------------------
% Extract individual GNSS measurements from the processed rinex data cell.
% Input: 
%   data_all (n*11 cell, n = number of epochs)
%       Col-1: GPS Time
%       Col-2: Positioning Solution (Ordinary Least Squares - corrected by rtklib)
%       Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
%       Col-4: Satellite Elevation Angle
%       Col-5: Satellite Azimuth Angle (NED, North as zero)
%       Col-6,7,8: Full-ephemeris PRN, El, Az
%       Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
%       Col-10: L1 PRN, Pseudorange with receiver_bias correction, SV-XYZ-Position in ECEF for positioning
%       Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
% Output: 
%   (All measurement outputs are n*m matrices, where n = number of epochs 
%   and m = total number of satellites visible throughout the observation)
%   gps_time: timestamps (n)
%   prn: all of the unique prns throughout the observation (m)
%   ele_angle: elevation angle
%   azi_angle: azimuth angle
%   CNR: carrier to noise ratio (C/N0)
%   pr_resi: pseudorange residual
%   carrier_ph: carrier phase
%   pseudorange: raw pseudorange with only receiver bias correction
%   dopp: doppler measurement in Hz
% -------------------------------------------------------------------------

gps_time = cell2mat(data_all(:,1));
prn = unique(vertcat(data_all{:,3}));

% Pre-allocate ele_angle, azi_angle, CNR, pr_resi, carrier_ph, pseudorange,
% and dopp matrices
ele_angle = NaN(height(data_all),length(prn));
azi_angle = NaN(height(data_all),length(prn));
CNR = NaN(height(data_all),length(prn));
pr_resi = NaN(height(data_all),length(prn));
carrier_ph = NaN(height(data_all),length(prn));
pseudorange = NaN(height(data_all),length(prn));
dopp = NaN(height(data_all),length(prn));

% For each timestamp, find the index of received satellites in prn and
% populate the measurement matrices
for epoch = 1:height(data_all)

    % Find the index of each received satellite in prn
    [~,idx] = ismember(data_all{epoch,3},prn);

    % If no data, proceed to next loop
    if isempty(idx)
        continue;
    end

    % Populate measurement matrices
    ele_angle(epoch,idx) = data_all{epoch,4};
    azi_angle(epoch,idx) = data_all{epoch,5};
    CNR(epoch,idx) = data_all{epoch,9}(:,1);
    pr_resi(epoch,idx) = data_all{epoch,9}(:,2);
    carrier_ph(epoch,idx) = data_all{epoch,9}(:,3);
    pseudorange(epoch,idx) = data_all{epoch,10}(:,2);
    dopp(epoch,idx) = data_all{epoch,11}(:,2);

end