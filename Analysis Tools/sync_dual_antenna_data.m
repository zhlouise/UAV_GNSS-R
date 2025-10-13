function [zenith_data_new,nadir_data_new] = sync_dual_antenna_data(zenith_data, nadir_data)
% -------------------------------------------------------------------------
% Sync both the zenith_data and nadir_data by the gps timestamp, such that
% the output data cells have the same dimensions.
% Input: 
%   zenith_data & nadir_data (n*11 cell, n = number of epochs)
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
%   zenith_data_new & nadir_data_new
%       The new data cells have the same dimensions. The first timestamp in
%       both cells is the most recent timestamp of the starting timestamps
%       in zenith_data and nadir_data while the last timestamp in both
%       cells is the least recent timestamp of the ending timestamps in
%       zenith_data and nadir_data.
% -------------------------------------------------------------------------

timestamp_start = max(round(zenith_data{1,1}),round(nadir_data{1,1}));
timestamp_end = min(round(zenith_data{end,1}),round(nadir_data{end,1}));

% Crop existing zenith_data
zenith_data_new = zenith_data(find(round([zenith_data{:,1}])==timestamp_start):...
find(round([zenith_data{:,1}])==timestamp_end),:);

% Align nadir_data with the timestamp of zenith_data_new
nadir_data_new = cell(height(zenith_data_new), 11);
for epoch = 1:height(zenith_data_new)
    nadir_data_new(epoch,:) = nadir_data(round([nadir_data{:,1}])==round(zenith_data{epoch,1}), :);
end

end