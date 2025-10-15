function [centroid, a, b] = cal_fresnel_zones(pos_llh, altitude_ag, ele_angle, azi_angle)
% -------------------------------------------------------------------------
% Calculates the center, semi-major axis, and semi-minor axis of the
% First Fresnel Zones (FFZs) for all satellites during the observation.
% (All n*m matrices indicate: n = number of epochs and m = total number of 
% satellites visible throughout the observation)
% Input: 
%   pos_llh: coordinates in latitude, longitude, and height (n*3 matrix)
%   altitude_ag: above-ground altitude (n*1 matrix)
%   ele_angle: satellite elevation angle (n*m matrix) 
%   azi_angle: satellite azimuth angle (n*m matrix)
%   
% Output:
%   centroid: centroid of the FFZ ellipse in LLH (n*m matrix) 
%   a: semi-major axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
%   b: semi-minor axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
% -------------------------------------------------------------------------

f = gt.C.FREQ1; % L1/E1/B1C frequency
c = gt.C.CLIGHT; % Speed of light
lambda = c/f; % Signal wavelength

% Replicate altitude_ag for all the satellite instances for matrix
% operation
altitude_ag_rep = repmat(altitude_ag,1,size(ele_angle,2));


% Calculation of FFZ center
% x_c is the horizontal distance between the receiver and the centroid
% along the projection of the vector between the receiver and the satellite
% on to the ground tangent plane
x_c = (lambda/2+altitude_ag_rep.*sind(ele_angle)).*cosd(ele_angle)./(sind(ele_angle)).^2;


end

