function [zenith_data_new,nadir_data_new] = sync_dual_antenna_data(zenith_data, nadir_data, flight_log)
% -------------------------------------------------------------------------
% Sync both the zenith_data and nadir_data by the gps timestamp of the 
% flight log, such that the output data cells have the same dimensions and 
% is of the same period as when the UAV was armed.
% Input: 
%   zenith_data & nadir_data (n*11 cell, n = number of epochs)
%       Col-1: GPS Time
%       Col-2: Positioning Solution (Ordinary Least Squares - corrected by rtklib)
%       Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
%       Col-4: Satellite Elevation Angles
%       Col-5: Satellite Azimuth Angle (NED, North as zero)
%       Col-6,7,8: Full-ephemeris PRN, El, Az
%       Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
%       Col-10: L1 PRN, Pseudorange with receiver_bias correction, SV-XYZ-Position in ECEF for positioning
%       Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
%   flight_log: table of topic messages from .ulg file logged by Pixhawk
% Output: 
%   zenith_data_new & nadir_data_new
%       The new data cells have the same dimensions. The first timestamp in
%       both cells is the first timestamp within the flight_log and the
%       last timestamp in both cells is the last timestamp within the
%       flight_log.
% -------------------------------------------------------------------------

% Retrieve timestamps in utc seconds
timestamp_start_utc = flight_log.TopicMessages{flight_log.TopicNames=='vehicle_gps_position', 1}{1,'time_utc_usec'};
timestamp_end_utc = flight_log.TopicMessages{flight_log.TopicNames=='vehicle_gps_position', 1}{end,'time_utc_usec'};

% Convert timestamps to GPS TOW
timestamp_start_gps = utcsec2gpstow(timestamp_start_utc/(1e+6)); % ulog utc time is in microseconds
timestamp_end_gps = utcsec2gpstow(timestamp_end_utc/(1e+6));

% Crop existing zenith_data
zenith_data_new = zenith_data(find(round([zenith_data{:,1}])==timestamp_start_gps):...
find(round([zenith_data{:,1}])==timestamp_end_gps),:);

% Align nadir_data with the timestamp of zenith_data_new
nadir_data_new = cell(height(zenith_data_new), 11);
for epoch = 1:height(zenith_data_new)
    nadir_data_new(epoch,:) = nadir_data(round([nadir_data{:,1}])==round(zenith_data{epoch,1}), :);
end

end