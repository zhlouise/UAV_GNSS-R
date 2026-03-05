function [centroid_lat, centroid_lon, a, b] = cal_fresnel_zones(pos_llh, altitude_ag, ele_angle, azi_angle)
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
%   centroid_lat: centroid latitude of the FFZ ellipse (n*m matrix) 
%   centroid_lon: centroid longitude of the FFZ ellipse (n*m matrix) 
%   a: semi-major axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
%   b: semi-minor axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
% -------------------------------------------------------------------------

f = gt.C.FREQ1; % L1/E1/B1C frequency
c = gt.C.CLIGHT; % Speed of light
lambda = c/f; % Signal wavelength

% Filter out invalid ele_angles
ele_angle(ele_angle<0) = NaN;

% Replicate altitude_ag and latitude, longitude and height of pos_llh 
% for all the satellite instances for matrix operations
altitude_ag_rep = repmat(altitude_ag, 1, size(ele_angle,2));
ref_lat_rep = repmat(pos_llh(:,1), 1, size(ele_angle,2));
ref_lon_rep = repmat(pos_llh(:,2), 1, size(ele_angle,2));
ref_height_rep = repmat(pos_llh(:,3), 1, size(ele_angle,2));

% Calculation of FFZ center
% x_c is the horizontal distance between the receiver and the centroid
% along the projection of the vector between the receiver and the satellite
% on to the ground tangent plane
x_c = (lambda/2+altitude_ag_rep.*sind(ele_angle)).*cosd(ele_angle)./(sind(ele_angle)).^2;
% Centroid coordinates in ENU
centroid_east = x_c.*sind(azi_angle);
centroid_north = x_c.*cosd(azi_angle);
% Convert ENU to LLH
[centroid_lat, centroid_lon, ~] = enu2geodetic(centroid_east, centroid_north, ...
    zeros(size(centroid_east)), ref_lat_rep, ref_lon_rep, ref_height_rep, wgs84Ellipsoid);

% Calculate semi-major axis
a = sqrt(lambda.*altitude_ag_rep.*sind(ele_angle)+lambda^2/4)./(sind(ele_angle)).^2;

% Calculate semi-minor axis
b = sqrt(lambda.*altitude_ag_rep.*sind(ele_angle)+lambda^2/4)./sind(ele_angle);

end