function [lat_grid, lon_grid, data_grid, lat_lim, lon_lim] = fresnel_zone_heatmap(centroid_lat, centroid_lon, a, b, azi_angle, data)
% -------------------------------------------------------------------------
% Rasterizes First Fresnel Zone (FFZ) ellipses using precomputed geometry
% to a georeferenced heatmap. Where ellipses overlap, each cell value is
% the mean of the provided per-ellipse values. This function does so by
% rasterizing the FFZ data onto a georeferenced grid. 
%
% Inputs:
%   centroid_lat: centroid latitude of the FFZ ellipse (n*m matrix) 
%   centroid_lon: centroid longitude of the FFZ ellipse (n*m matrix) 
%   a: semi-major axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
%   b: semi-minor axis of the FFZ in meters (ENU coordinate) (n*m matrix) 
%   azi_angle: satellite azimuth angle (n*m matrix)
%   data: measurement values to be visualized (n*m matrix)
% Outputs:
%   lat_grid: latitude of all grid cells to be plotted (vector)
%   lon_grid: longitude of all grid cells to be plotted (vector)
%   data_grid: data value of all grid cells to be plotted (vector)
%   lat_lim: latitude limits of the grid [lat_min lat_max]
%   lon_lim: longitude limits of the grid [lon_min lon_max]
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
num_data = numel(centroid_lat_v); % Number of data points

% Determine grid bounds (with added padding based on maximum ellipse size,
% which is an approximate conversion from distance to degrees in llh)
lat_min = min(centroid_lat_v) - max(a_v)/111320;
lat_max = max(centroid_lat_v) + max(a_v)/111320;
lon_min = min(centroid_lon_v) - max(a_v)/(111320*cosd(mean(centroid_lat_v))); % Using mean lat for approximation
lon_max = max(centroid_lon_v) + max(a_v)/(111320*cosd(mean(centroid_lat_v)));

% Grid resolution setting
resolution = min(b_v)/2; % In meters
resolution_deg = resolution/111320; % Approximated in degrees
fprintf('Grid resolution: %.6f degrees\n', resolution_deg);

% Allocate grid
grid_lat = lat_min:resolution_deg:lat_max;
grid_lon = lon_min:resolution_deg:lon_max;

% Initialize accumulation cell
value_list = cell(length(grid_lat), length(grid_lon)); 

% Loop for each FFZ to calculate data value at each grid point within the FFZ
for id = 1:num_data
    
    % Define bounds to minimize points to search in the grid
    lat_lowerbound = centroid_lat_v(id) - a_v(id)/111320;
    lat_upperbound = centroid_lat_v(id) + a_v(id)/111320;
    lon_lowerbound = centroid_lon_v(id) - a_v(id)/(111320*cosd(centroid_lat_v(id)));
    lon_upperbound = centroid_lon_v(id) + a_v(id)/(111320*cosd(centroid_lat_v(id)));

    % Find grid points within bound
    [searchgrid_lat, searchgrid_lon] = meshgrid(grid_lat(grid_lat>=lat_lowerbound & grid_lat<=lat_upperbound), ...
        grid_lon(grid_lon>=lon_lowerbound & grid_lon<=lon_upperbound));
    searchgrid_lat_v = searchgrid_lat(:);
    searchgrid_lon_v = searchgrid_lon(:);
    
    % Convert search grid points to ENU coordinates
    [searchgrid_enu, ~] = rtklib.llh2enu([searchgrid_lat_v,searchgrid_lon_v,zeros(length(searchgrid_lat_v),1)],...
        [centroid_lat_v(id),centroid_lon_v(id),0]);
    east = searchgrid_enu(:,1);
    north = searchgrid_enu(:,2);

    % Convert azimuth (clockwise from North) to rotation angle (counterclockwise from East)
    theta = 90 - az_v(id);
    
    % Check in ENU coordinates which points are inside the (rotated) ellipse
    inside = ((east.*cosd(theta) + north.*sind(theta))./a_v(id)).^2 + ...
        ((east.*sind(theta) - north.*cosd(theta))./b_v(id)).^2 <= 1; % Equation of a rotated ellipse
    inside_lat = searchgrid_lat_v(inside);
    inside_lon = searchgrid_lon_v(inside);

    % Find the indices of the points inside the ellipse in the accumulation arrays
    lat_idx = 1+round((inside_lat-lat_min)./resolution_deg);
    lon_idx = 1+round((inside_lon-lon_min)./resolution_deg);

    % Update the accumulation array
    idx = sub2ind(size(value_list), lat_idx, lon_idx);
    for k = 1:numel(idx)
        value_list{idx(k)}(end+1) = data_v(id); % Append in the grid data
    end
end

% Calculate the mean for overlapping grid points
heatmap_data = cellfun(@mean, value_list);

% % Calculate the standard deviation for overlapping grid points
% heatmap_data = cellfun(@std, value_list);
% heatmap_data(heatmap_data==0) = NaN; % Remove grids with only 1 ellipse point

% Create vectors for all valid grid points with their values for output
[lon, lat] = meshgrid(grid_lon, grid_lat);
lat_grid = lat(~isnan(heatmap_data));
lon_grid = lon(~isnan(heatmap_data));
data_grid = heatmap_data(~isnan(heatmap_data));
lat_lim = [lat_min lat_max];
lon_lim = [lon_min lon_max];

% Plot heatmap on satellite imagery
figure();
geolimits(lat_lim, lon_lim);
geobasemap('satellite');
hold on;
% Determine marker size based on grid resolution
marker_size = max(1, 100*resolution_deg);
geoscatter(lat_grid, lon_grid, marker_size, data_grid, 'filled');
colormap(turbo);
colorbar;
hold off;

end