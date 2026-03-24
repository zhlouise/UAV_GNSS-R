function fresnel_zone_plot(centroid_lat,centroid_lon, a, b, azi_angle, data)
% -------------------------------------------------------------------------
% Visualizes fresnel zones on satellite imagery. This function also 
% overlays measurement values onto a satellite map, where each value is 
% represented by color at its corresponding FFZ.
% Input: 
%   centroid_lat: centroid latitude of the FFZ ellipse (n*m matrix) 
%   centroid_lon: centroid longitude of the FFZ ellipse (n*m matrix) 
%   a: semi-major axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
%   b: semi-minor axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
%   azi_angle: satellite azimuth angle (n*m matrix)
%   data: measurement values to be visualized (n*m matrix)
% -------------------------------------------------------------------------

% Flatten inputs to vectors
centroid_lat_v = centroid_lat(:);
centroid_lon_v = centroid_lon(:);
a_v = a(:);
b_v = b(:);
az_v = azi_angle(:);
data_v = data(:);

% Filter out invalid data
valid_idx = ~isnan(centroid_lat_v) & ~isnan(centroid_lon_v) & ~isnan(a_v) & ~isnan(b_v) & ~isnan(az_v) & ~isnan(data_v) & (a_v>0) & (b_v>0);
centroid_lat_v = centroid_lat_v(valid_idx);
centroid_lon_v = centroid_lon_v(valid_idx);
a_v = a_v(valid_idx);
b_v = b_v(valid_idx);
az_v = az_v(valid_idx);
data_v = data_v(valid_idx);
num_ellipses = numel(centroid_lat_v);

% Plot centroids using sat_img_plot
pos_llh = [centroid_lat_v, centroid_lon_v, zeros(num_ellipses,1)];
sat_img_plot(pos_llh, data_v);
hold on;

% Setup for ellipse plotting
num_points = 30;
t = linspace(0, 2*pi, num_points);
cmap = turbo(256);
data_min = min(data_v);
data_max = max(data_v);
% For cases when all the data is just one value to prevent errors
if data_min == data_max
    data_max = data_min + 1;
end

% Map data values to RGB colors
mapValueToRGB = @(val) interp1(linspace(data_min,data_max,size(cmap,1))', cmap, val, 'linear', 'extrap');

% Plot each ellipse
for id = 1:num_ellipses
    
    % Convert azimuth (clockwise from North) to rotation angle (counterclockwise from East)
    theta = 90 - az_v(id);
    % Create ellipse in ENU coordinates (parametric expression of an ellipse)
    x = a_v(id)*cos(t); % East
    y = b_v(id)*sin(t); % North
    % Rotate ellipse
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % 2D rotation matrix
    ellipse_enu = (R*[x;y]).';
    % Convert from ENU to llh
    [ellipse_llh, ~] = rtklib.enu2llh([ellipse_enu, zeros(height(ellipse_enu),1)], [centroid_lat_v(id),centroid_lon_v(id),0]);
    % Get color for this ellipse
    rgb = mapValueToRGB(data_v(id));
    % Plot ellipse boundary
    geoplot(ellipse_llh(:,1), ellipse_llh(:,2), 'Color', rgb, 'LineWidth', 1);

end

hold off;

end

