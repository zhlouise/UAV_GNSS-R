function altitude_ag = cal_altitude_above_ground_HK(pos_llh,dtm_path)
% -------------------------------------------------------------------------
% Calculates the above ground altitude from the WGS84 LLH coordinates and
% the local digital terrain model (CODE SPECIFIED FOR HK DTM USAGE). 
% Input: 
%   pos_llh: coordinates in latitude, longitude, and height (n*3 matrix)
%   dtm_path: filepath for the digital terrain model (asc format)
%   
% Output:
%   altitude_ag: above-ground altitude (n*1 matrix)
% -------------------------------------------------------------------------

fid = fopen(dtm_path);
header = textscan(fid, '%s %f', 6);
fclose(fid);

ncols = header{2}(1);
nrows = header{2}(2);
xllcorner = header{2}(3);
yllcorner = header{2}(4);
cellsize = header{2}(5);
NODATA_value = header{2}(6);

% Read elevation grid
dtm = readmatrix(dtm_path, 'NumHeaderLines', 6);

x_grid = xllcorner + cellsize * (0:ncols-1);         % HK80 Easting
y_grid = yllcorner + cellsize * (nrows-1:-1:0);      % HK80 Northing

% Convert LLH (WGS84) to HK80 grid using Mapping Toolbox
hk80 = projcrs(2326); % HK80 Grid EPSG:2326

% Apply the additive constants in WGS84 to HK80 geographic coordinates conversion
lat_HK80 = pos_llh(:,1) + 5.5/3600;  % +5.5 seconds to latitude
lon_HK80 = pos_llh(:,2) - 8.8/3600;  % -8.8 seconds to longitude

% Convert HK80 latitude and longidude to HK80 easting and northing
[x_hk80, y_hk80] = projfwd(hk80, lat_HK80, lon_HK80);

% Interpolate ground elevation from dtm
ground_elevation = interp2(x_grid, y_grid, dtm, x_hk80, y_hk80, 'linear', NODATA_value);
ground_elevation(ground_elevation == NODATA_value) = NaN; % Handle nodata

% Get geoid undulation (distance between MSL and WGS84 ellipsoid)
N = egm96geoid(pos_llh(:,1), pos_llh(:,2)); % [lat], [lon] in WGS84

% Orthometric height (approximated height above MSL) = ellipsoidal height
% (height obtained from GNSS) - geoid undulation
MSL_to_HKPD_offset = 1.30; % MSL is about 1.30m above HKPD zero
height_HKPD = pos_llh(:,3) - N + MSL_to_HKPD_offset;

% Calculate above-ground altitude
altitude_ag = height_HKPD - ground_elevation;

end