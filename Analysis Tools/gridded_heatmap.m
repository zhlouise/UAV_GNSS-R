function gridded_heatmap(centroid_lat, centroid_lon, data)
% -------------------------------------------------------------------------
% From the location information of the First Fresnel Zone centroid, divide
% the map region into 5m*5m grids, where each grid contains the data of FFZ
% centroids that falls onto the grid. A georeferenced heatmap of the grids
% is plotted, where the grid color indicates the mean/std of the data. 
% 
% Inputs:
%   centroid_lat: centroid latitude of the FFZ ellipse (n*m matrix) 
%   centroid_lon: centroid longitude of the FFZ ellipse (n*m matrix) 
%   data: measurement values to be visualized (n*m matrix)
% -------------------------------------------------------------------------

% Flatten inputs to vectors
centroid_lat_v = centroid_lat(:);
centroid_lon_v = centroid_lon(:);
data_v = data(:);

% Filter out invalid data
valid_idx = ~isnan(centroid_lat_v) & ~isnan(centroid_lon_v) & ~isnan(data_v);
centroid_lat_v = centroid_lat_v(valid_idx);
centroid_lon_v = centroid_lon_v(valid_idx);
data_v = data_v(valid_idx);
num_data = numel(data_v); % Number of data points

% Grid resolution setting
resolution = 5; % In meters
resolution_lat_deg = resolution/111320; % Approximated in degrees
resolution_lon_deg = resolution/(111320*cosd(mean(centroid_lat_v))); % Using mean lat for approximation

% Determine grid bounds (with added padding of one cell)
lat_min = min(centroid_lat_v) - resolution_lat_deg;
lat_max = max(centroid_lat_v) + resolution_lat_deg;
lon_min = min(centroid_lon_v) - resolution_lon_deg;
lon_max = max(centroid_lon_v) + resolution_lon_deg;

% Allocate grid
grid_lat = lat_min:resolution_lat_deg:lat_max;
grid_lon = lon_min:resolution_lon_deg:lon_max;

% Initialize accumulation cell
value_list = cell(length(grid_lat), length(grid_lon));

% Map each centroid to a grid cell index
lat_idx = 1+floor((centroid_lat_v-lat_min)./resolution_lat_deg);
lon_idx = 1+floor((centroid_lon_v-lon_min)./resolution_lon_deg);
idx = sub2ind(size(value_list), lat_idx, lon_idx);

% Loop for each FFZ to assign it into a cell
for id = 1:num_data
    value_list{idx(id)}(end+1) = data_v(id); % Append in the grid data
end

% Calculate the mean for overlapping grid points
heatmap_data = cellfun(@mean, value_list);

% % Calculate the standard deviation for overlapping grid points
% heatmap_data = cellfun(@std, value_list);
% heatmap_data(heatmap_data==0) = NaN; % Remove grids with only 1 ellipse point

% Create vectors for plotting
[lon, lat] = meshgrid(grid_lon, grid_lat);
lat_grid = lat(~isnan(heatmap_data));
lon_grid = lon(~isnan(heatmap_data));
data_grid = heatmap_data(~isnan(heatmap_data));
lat_lim = [lat_min lat_max];
lon_lim = [lon_min lon_max];

% Plot heatmap on satellite imagery
figure();
geolimits(lat_lim, lon_lim);
% addCustomBasemap('GoogleSat', 'https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}')
geobasemap('satellite');
% geobasemap('GoogleSat');
hold on;
geoscatter(lat_grid+resolution_lat_deg/2, lon_grid+resolution_lon_deg/2, 50, data_grid,'square','filled','MarkerFaceAlpha', 0.75);
colormap(turbo);
colorbar;
hold off;

end