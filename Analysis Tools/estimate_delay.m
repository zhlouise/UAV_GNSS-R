function reflection_delay = estimate_delay(altitude_ag, ele_angle)
% -------------------------------------------------------------------------
% Estimates the delay of the reflected signals received by the nadir
% antenna through calculating 2*H*sin(theta), where H is the above ground
% height of the receiver and theta is the elevation angle. 
% (All n*m matrices indicate: n = number of epochs and m = total number of 
% satellites visible throughout the observation)
% Input: 
%   ele_angle: satellite elevation angle (n*m matrix) 
%   azi_angle: satellite azimuth angle (n*m matrix)
%   
% Output:
%   reflection_delay (n*m matrix) 
% -------------------------------------------------------------------------

% Replicate altitude_ag for all the satellite instances for matrix 
% operations
altitude_ag_rep = repmat(altitude_ag, 1, size(ele_angle,2));

reflection_delay = 2*altitude_ag_rep./sind(ele_angle);

end