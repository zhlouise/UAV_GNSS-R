% -------------------------------------------------------------------------
% A converter tool that converts downloaded aerial imagery into 
% georeferenced TIFF file. Requires the input of the longitude and latitude
% of the four corners of the aerial imagery. 
% This tool assumes that the aerial imagery is already rotated with north
% being the up side of the image.
% -------------------------------------------------------------------------

clc;
clear;
close all;

% Prompt for input of the aerial imagery
[filename, path] = uigetfile('.jpg','Select the aerial imagery.');
img = imread([path, filename]);
[rows, cols, ~] = size(img);

% Prompt for the geographic coordinates of the 4 corners of the image
geo_coordinates_prompt = {'Enter corner 1 longitude (decimal degrees): ', ...
    'Enter corner 1 latitude (decimal degrees): ', ...
    'Enter corner 2 longitude (decimal degrees): ', ...
    'Enter corner 2 latitude (decimal degrees): ', ...
    'Enter corner 3 longitude (decimal degrees): ', ...
    'Enter corner 3 latitude (decimal degrees): ', ...
    'Enter corner 4 longitude (decimal degrees): ', ...
    'Enter corner 4 latitude (decimal degrees): '};
geo_coordinates_input = inputdlg(geo_coordinates_prompt, 'Input Corner Coordinates');
% Order: top left; top right; bottom right; bottom left
geo_coordinates = [
    str2double(geo_coordinates_input{1}), str2double(geo_coordinates_input{2});
    str2double(geo_coordinates_input{3}), str2double(geo_coordinates_input{4});
    str2double(geo_coordinates_input{5}), str2double(geo_coordinates_input{6});
    str2double(geo_coordinates_input{7}), str2double(geo_coordinates_input{8});
];

lon_lim = [min(geo_coordinates(:,1)) max(geo_coordinates(:,1))];
lat_lim = [min(geo_coordinates(:,2)) max(geo_coordinates(:,2))];

% Create the geographic referencing object
R = georefcells(lat_lim, lon_lim, [rows cols], 'ColumnsStartFrom', 'north');

% Save as a GeoTIFF
[output_filename, output_pathname] = uiputfile('.tif', 'Save as GeoTIFF.');
geotiffwrite([output_pathname, output_filename], img, R);