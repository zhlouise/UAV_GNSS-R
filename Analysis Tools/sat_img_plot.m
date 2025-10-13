function sat_img_plot(pos_llh, data)
% -------------------------------------------------------------------------
% Visualizes georeferenced measurement data on satellite imagery. This 
% function overlays measurement values onto a satellite map, where each
% value is represented by color at its corresponding location.
% 
% Input: 
%   pos_llh: coordinates in latitude, longitude, and height (n*3 matrix)
%   data: measurement values to be visualized (n*3 matrix)
% -------------------------------------------------------------------------

figure();

geoplot(pos_llh(:,1), pos_llh(:,2), data, 'MarkerSize', 10);
geobasemap('satellite');

end

