function ground_truth = extract_ground_truth(gps_time, flight_log)
% -------------------------------------------------------------------------
% Extract the UAV ground truth in LLH coordinates from the Pixhawk log. 
% Input: 
%   gps_time: time of the week (n*1 array)
%   flight_log: table of topic messages from .ulg file logged by Pixhawk
% Output: 
%   ground_truth: UAV position in LLH coordinates (n*3 array)
% -------------------------------------------------------------------------

% Pre-allocate ground_truth array
ground_truth = NaN(height(gps_time), 3);

% Convert all timestamps in flight_log's gps output from utc to gps tow
utc_timestamp = flight_log.TopicMessages{flight_log.TopicNames=='vehicle_gps_position', 1}{:,'time_utc_usec'};
gps_timestamp = utcsec2gpstow(utc_timestamp./(1e+6)); % ulog utc time is in microseconds

% Sample the GPS solutions in 1 Hz (according to the input gps_time from
% the dual antenna setup) to output the Pixhawk system time of each sample
[~, sampled_idx] = ismember(round(gps_time), gps_timestamp);
sampled_system_timestamp = flight_log.TopicMessages{flight_log.TopicNames=='vehicle_gps_position', 1}.timestamp(sampled_idx);

% Extract ground_truth from flight_log's vehicle_global_position (filtered
% positioning solution) at the sampled system timestamps
for epoch = 1:height(gps_time)
    % Find the timestamp in vehicle_global_position that is the closest to
    % the sampled timestamp from vehicle_gps_position
    [~,idx] = min(abs(flight_log.TopicMessages{flight_log.TopicNames=='vehicle_global_position', 1}.timestamp - ...
        sampled_system_timestamp(epoch))); 
    ground_truth(epoch,:) = flight_log.TopicMessages{flight_log.TopicNames=='vehicle_global_position', 1}{idx, 2:4};
end

end

