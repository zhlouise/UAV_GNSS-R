function data_skyplot(ele_angle, azi_angle, data)
% -------------------------------------------------------------------------
% Visualizes satellite measurement data on a fisheye skyplot. This function
% overlays measurement values onto a skyplot, where each value is
% represented by color.
% Input: 
%   ele_angle: satellite elevation angle (n*m matrix)
%   azi_angle: azimuth angle (n*m matrix)
%   data: measurement values to be visualized (n*m matrix)
% -------------------------------------------------------------------------

figure();

for i = 1:size(ele_angle, 2)
    polarscatter(deg2rad(azi_angle(:,i)), ele_angle(:,i), 75, data(:,i),'filled');
    hold on;
end

colorbar;
clim([0 55]);

% Configure the polar axis
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.RAxisLocation = 270;
pax.RLim = [0, 90];
pax.RDir = 'reverse';

% Configure figure and axis titles
thetaticks([0 30 60 90 120 150 180 210 240 270 300 330]);
thetaticklabels({'N','30','60','E','120','150','S','210','240','W','300','330'});
set(gca,'color','none');

hold off;

end

