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
latv = centroid_lat(:);
lonv = centroid_lon(:);
av = a(:);
bv = b(:);
azv = azi_angle(:);
dv = data(:);
nEllipses = numel(latv);

% Filter valid data
validIdx = ~isnan(latv) & ~isnan(lonv) & ~isnan(av) & ~isnan(bv) & ~isnan(azv) & ~isnan(dv) & (av>0) & (bv>0);

% Plot centroids using sat_img_plot
pos_llh = [latv, lonv, zeros(nEllipses,1)];
sat_img_plot(pos_llh(validIdx,:), dv(validIdx));
hold on;

% Setup for ellipse plotting
nPoints = 96;
t = linspace(0, 2*pi, nPoints);
cmap = jet(256);
dataMin = min(dv(validIdx));
dataMax = max(dv(validIdx));
if dataMin == dataMax
    dataMax = dataMin + 1;
end

% Map data values to RGB colors
mapValueToRGB = @(val) interp1(linspace(dataMin,dataMax,size(cmap,1))', cmap, val, 'linear', 'extrap');

% Plot each ellipse
for k = 1:nEllipses
    if ~validIdx(k), continue; end
    
    a_k = av(k);
    b_k = bv(k);
    az_k = azv(k);
    lat0 = latv(k);
    lon0 = lonv(k);
    
    % Rotation angle: convert azimuth (CW from North) to theta (CCW from East)
    theta_k = 90 - az_k;
    
    % Create ellipse in ENU coordinates (meters)
    x = a_k * cos(t);  % east
    y = b_k * sin(t);  % north
    
    % Rotate ellipse
    R = [cosd(theta_k) -sind(theta_k); sind(theta_k) cosd(theta_k)];
    pts = R * [x; y];
    east = pts(1,:);
    north = pts(2,:);
    
    % Convert meters to degrees (approximate)
    dlat = north / 111320;
    dlon = east / (111320 * cosd(lat0));
    latPts = lat0 + dlat;
    lonPts = lon0 + dlon;
    
    % Get color for this ellipse
    rgb = mapValueToRGB(dv(k));
    
    % Plot ellipse boundary
    geoplot(latPts, lonPts, 'Color', rgb, 'LineWidth', 1.5);
end

end

