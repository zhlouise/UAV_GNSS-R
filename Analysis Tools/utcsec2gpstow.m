function gps_tow = utcsec2gpstow(utc_sec)
% -------------------------------------------------------------------------
% Converts UTC seconds into gps time of the week (in seconds). 
% Input: 
%   utc_sec: UTC second (scalar)
% Output:
%   gps_tow: GPS time of the week in seconds (scalar)
% -------------------------------------------------------------------------

t = datetime(utc_sec, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
utc_time = [year(t), month(t), day(t), hour(t), minute(t), second(t)];

[gps_tow, ~] = rtklib.epoch2tow(utc_time,1);

end

