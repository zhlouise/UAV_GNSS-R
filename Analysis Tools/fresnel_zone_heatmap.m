function fresnel_zone_heatmap(centroid_lat, centroid_lon, a, b, azi_angle, data)
% -------------------------------------------------------------------------
% Rasterizes First Fresnel Zone (FFZ) ellipses using precomputed geometry
% to a georeferenced heatmap. Where ellipses overlap, each cell value is
% the mean of the provided per-ellipse values.
% Inputs:
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

% Filter valid data
validIdx = ~isnan(latv) & ~isnan(lonv) & ~isnan(av) & ~isnan(bv) & ~isnan(azv) & ~isnan(dv) & (av>0) & (bv>0);

if ~any(validIdx)
    error('No valid data points found');
end

% Extract valid data only
validLat = latv(validIdx);
validLon = lonv(validIdx);
validA = av(validIdx);
validB = bv(validIdx);
validAz = azv(validIdx);
validData = dv(validIdx);
nValid = numel(validLat);

fprintf('Processing %d ellipses...\n', nValid);

% Determine grid bounds with some padding
latMin = min(validLat);
latMax = max(validLat);
lonMin = min(validLon);
lonMax = max(validLon);

% Add padding based on maximum ellipse size
maxSemiMajor = max(validA);
paddingLat = maxSemiMajor / 111320;
paddingLon = maxSemiMajor / (111320 * cosd(mean(validLat)));

latMin = latMin - paddingLat;
latMax = latMax + paddingLat;
lonMin = lonMin - paddingLon;
lonMax = lonMax + paddingLon;

% OPTIMIZATION 1: Use adaptive resolution
meanSemiMinor = mean(validB);
resolution = meanSemiMinor / 8;  % More reasonable resolution
resolutionDeg = resolution / 111320;

fprintf('Grid resolution: %.6f degrees\n', resolutionDeg);

% Create grid
latGrid = latMin:resolutionDeg:latMax;
lonGrid = lonMin:resolutionDeg:lonMax;
[nRows, nCols] = deal(length(latGrid), length(lonGrid));

fprintf('Grid size: %d x %d = %d total points\n', nRows, nCols, nRows*nCols);

% Initialize accumulation arrays
valueSum = zeros(nRows, nCols);
weightSum = zeros(nRows, nCols);

% OPTIMIZATION 2: Precompute constants
metersPerDegLat = 111320;
R_earth = 6371000; % Earth radius in meters

% Process each valid ellipse with progress indicator
fprintf('Processing ellipses: ');
progressInterval = max(1, floor(nValid/10));

for k = 1:nValid
    if mod(k, progressInterval) == 0
        fprintf('%d... ', k);
    end
    
    a_k = validA(k);
    b_k = validB(k);
    az_k = validAz(k);
    lat0 = validLat(k);
    lon0 = validLon(k);
    dataVal = validData(k);
    
    % Convert azimuth to rotation angle (CCW from East)
    theta_k = 90 - az_k;
    
    % OPTIMIZATION 3: More efficient bounding box calculation
    metersPerDegLon = 111320 * cosd(lat0);
    halfHeightDeg = a_k / metersPerDegLat;  % Use semi-major for conservative bounds
    halfWidthDeg = a_k / metersPerDegLon;
    
    bboxLat = [lat0 - halfHeightDeg, lat0 + halfHeightDeg];
    bboxLon = [lon0 - halfWidthDeg, lon0 + halfWidthDeg];
    
    % Find grid indices within bounding box
    latIdx = find(latGrid >= bboxLat(1) & latGrid <= bboxLat(2));
    lonIdx = find(lonGrid >= bboxLon(1) & lonGrid <= bboxLon(2));
    
    if isempty(latIdx) || isempty(lonIdx)
        continue;
    end
    
    % OPTIMIZATION 4: Vectorize the inner loop
    [I, J] = meshgrid(latIdx, lonIdx);
    I = I(:);
    J = J(:);
    
    gridLats = latGrid(I);
    gridLons = lonGrid(J);
    
    % Convert grid points to ENU coordinates
    north = (gridLats - lat0) * metersPerDegLat;
    east = (gridLons - lon0) * metersPerDegLon;
    
    % Rotate points to ellipse-aligned coordinate system
    cosTheta = cosd(theta_k);
    sinTheta = sind(theta_k);
    
    x_rot = east * cosTheta + north * sinTheta;
    y_rot = -east * sinTheta + north * cosTheta;
    
    % Check which points are inside ellipse
    inside = (x_rot./a_k).^2 + (y_rot./b_k).^2 <= 1;
    
    % Update accumulation arrays for points inside ellipse
    validIndices = find(inside);
    for idx = 1:length(validIndices)
        pos = validIndices(idx);
        valueSum(I(pos), J(pos)) = valueSum(I(pos), J(pos)) + dataVal;
        weightSum(I(pos), J(pos)) = weightSum(I(pos), J(pos)) + 1;
    end
end

fprintf('Done!\n');

% Compute mean values where ellipses overlap
heatmapData = valueSum ./ weightSum;
heatmapData(weightSum == 0) = NaN;

% Create vectors for all valid grid points with their values
[LON, LAT] = meshgrid(lonGrid, latGrid);
validHeatmap = ~isnan(heatmapData);

% Plot satellite imagery background
figure();
geolimits([latMin latMax], [lonMin lonMax]);
geobasemap('satellite');
hold on;

% Use geoscatter for compatible geographic plotting
latPoints = LAT(validHeatmap);
lonPoints = LON(validHeatmap);
dataPoints = heatmapData(validHeatmap);

% Determine marker size based on grid resolution
markerSize = max(1, 80 * resolutionDeg);

h = geoscatter(latPoints, lonPoints, markerSize, dataPoints, 'filled');
h.MarkerEdgeAlpha = 0.3;
h.MarkerFaceAlpha = 0.7;

% Add colorbar and styling
colormap(hsv);
colorbar;

hold off;

end