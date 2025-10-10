% ------------------------------------------------------------------------
% Code to analyse measurement-level data processed from rinex obs files.
%
% Data cells in the following structure:
% Col-1: GPS Time
% Col-2: Positioning Solution (Ordinary Least Squares - corrected by rtklib)
% Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
% Col-4: Satellite Elevation Angle
% Col-5: Satellite Azimuth Angle (NED, North as zero)
% Col-6,7,8: Full-ephemeris PRN, El, Az
% Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
% Col-10: L1 PRN, Pseudorange with receiver_bias correction, SV-XYZ-Position in ECEF for positioning
% Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
% ------------------------------------------------------------------------

clc;
clear;
close all;

% Prompt for input of processed dual antenna files
[filename, path] = uigetfile('.mat','Select the processed dual antenna .mat file.');
load([path, filename]);
zenith_data = dual_antenna_data{1,1};
nadir_data = dual_antenna_data{1,2};

%% Percentage of SV received

sv_ratio_z = NaN(height(zenith_data),1); % Zenith data
for epoch = 1:height(zenith_data)
    sv_ratio_z(epoch) = length(zenith_data{epoch,3})/length(zenith_data{epoch,6});
end

sv_ratio_n = NaN(height(nadir_data),1); % Nadir data
for epoch = 1:height(nadir_data)
    sv_ratio_n(epoch) = length(nadir_data{epoch,3})/length(nadir_data{epoch,6});
end

% Plot SV reception ratio distributions
figure();
grid on;
hold on;
histogram(sv_ratio_z,'Normalization','pdf','BinWidth',0.01);
histogram(sv_ratio_n,'Normalization','pdf','BinWidth',0.01);
hold off;
xlabel('Satellite Reception Ratio');
ylabel('Probability Distribution');
legend('Zenith Data', 'Nadir Data');

%% Pseudorange difference between zenith and nadir (nadir - zenith)

pr_difference = cell(height(zenith_data),1);

% Use the timestamps on zenith_data as a standard
for epoch = 1:height(zenith_data)
    
    temp_time = round(zenith_data{epoch,1});
    % Align nadir_data timestamp to the zenith_data timestamp
    idn = find(round(cell2mat(nadir_data(:,1)))==temp_time);
    if isempty(idn)
        % No nadir data in this timestamp
        continue;
    end

    % Only use satellites that are received by both the zenith and nadir
    % antenna
    temp_prn_list = intersect(zenith_data{epoch,3}, nadir_data{idn,3});
    
    % Column 1: prn
    % Column 2: pseudorange difference
    % Column 3: elevation angle (zenith as standard)
    % Column 4: azimuth angle (zenith as standard)
    % Column 5: CNR ratio (nadir/zenith)
    pr_difference{epoch} = [temp_prn_list, ...
        nadir_data{idn,10}(ismember(nadir_data{idn,3},temp_prn_list),2)-zenith_data{epoch,10}(ismember(zenith_data{epoch,3},temp_prn_list),2), ...
        zenith_data{epoch,4}(ismember(zenith_data{epoch,3},temp_prn_list)), ...
        zenith_data{epoch,5}(ismember(zenith_data{epoch,3},temp_prn_list)), ...
        nadir_data{idn,9}(ismember(nadir_data{idn,3},temp_prn_list),1)./zenith_data{epoch,9}(ismember(zenith_data{epoch,3},temp_prn_list),1)];

end

% Compile pr_difference into matrix form (remove epoch dimension)
pr_difference_mat = vertcat(pr_difference{:,1});

% Plot pseudorange difference vs elevation angle
figure();
grid on;
scatter(pr_difference_mat(:,3), pr_difference_mat(:,2),'filled');
xlabel('Elevation Angle (From Zenith Data)');
ylabel('Difference in Nadir and Zenith Pseudorange (m)');

% Plot pseudorange difference vs azimuth (on polar coordinates)
figure();
grid on;
polarscatter(deg2rad(pr_difference_mat(:,4)), pr_difference_mat(:,2),'filled');
hold on;
% Configure the polar axis
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.RAxisLocation = 0;
% Configure figure and axis titles
thetaticks([0 30 60 90 120 150 180 210 240 270 300 330]);
thetaticklabels({'North','30','60','East','120','150','South','210','240','West','300','330'});
title('Difference in Nadir and Zenith Pseudorange (m)');
hold off;

% Plot pseudorange difference vs CNR ratio
figure();
grid on;
scatter(pr_difference_mat(:,5), pr_difference_mat(:,2),'filled');
xlabel('Difference in Zenith and Nadir CNR (dB-Hz)');
ylabel('Difference in Nadir and Zenith Pseudorange (m)');

% Plot CNR ratio distribution
figure();
grid on;
histogram(pr_difference_mat(:,5),'BinWidth',0.01)
xlim([0, 2.5])
xlabel('CNR ratio between Nadir and Zenith');
ylabel('Probability Distribution');

%% Estimate dielectric constants

% Constants
c = 299792458; % speed of light (m/s)
f_l1 = 1575.42e6; % GPS L1 frequency (Hz)
lambda = c/f_l1; % wavelength (m)

% Assume calibration constant C = 1 (no calibration) or calibrate using water
% For now, we'll set C=1 and you can add water calibration later
C = 1;

% Initialize storage for dielectric constant results
dielectric_results = cell(height(zenith_data), 1);

% Use the timestamps on zenith_data as a standard
for epoch = 1:height(zenith_data)
    
    temp_time = round(zenith_data{epoch,1});
    % Align nadir_data timestamp to the zenith_data timestamp
    idn = find(round(cell2mat(nadir_data(:,1))) == temp_time);
    if isempty(idn)
        continue;
    end

    % Get common satellites between zenith and nadir
    temp_prn_list = intersect(zenith_data{epoch,3}, nadir_data{idn,3});
    if isempty(temp_prn_list)
        continue;
    end
    
    % Get indices for common PRNs
    [~, idx_z] = ismember(temp_prn_list, zenith_data{epoch,3});
    [~, idx_n] = ismember(temp_prn_list, nadir_data{idn,3});
    
    % Initialize results for this epoch
    epoch_results = [];
    
    for i = 1:length(temp_prn_list)
        prn = temp_prn_list(i);
        
        if prn <= 32 % GPS satellite only
            % Extract data for this satellite
            elev_angle = zenith_data{epoch,4}(idx_z(i)); % elevation angle from zenith
            theta = 90 - elev_angle; % incidence angle (degrees)
            theta_rad = deg2rad(theta); % convert to radians
            
            % Use CNR as proxy for SNR (convert from dB to linear)
            snr_zenith_linear = 10^(zenith_data{epoch,9}(idx_z(i),1)/10); % direct RHCP signal
            snr_nadir_linear = 10^(nadir_data{idn,9}(idx_n(i),1)/10); % reflected LHCP signal
            
            % Calculate power ratio (Eq. 16 in paper)
            power_ratio = snr_nadir_linear / snr_zenith_linear;
            
            % Get satellite and receiver positions for geometric calculations
            % Assuming you have satellite positions in ECEF from col 10
            sv_pos = zenith_data{epoch,10}(idx_z(i), 3:5); % satellite ECEF position
            rx_pos = zenith_data{epoch,2}(1:3); % receiver ECEF position from OLS solution
            
            % Calculate distances
            R3 = norm(sv_pos - rx_pos); % direct path distance (satellite to receiver)
            
            R1_plus_R2_approx = nadir_data{epoch,10}(idx_n(i),2); % approximate reflected path length is pseudorange
            
            % Calculate the geometric factor from Eq. 16
            geometric_factor = (R3^2) / (4 * R1_plus_R2_approx^2);
            
            % Rearrange Eq. 16 to solve for (R_vv - R_hh)^2
            fresnel_diff_squared = power_ratio / (geometric_factor * C);
            
            % Now we need to solve for dielectric constant using Fresnel equations
            % This requires numerical solution since it's transcendental
            
            % Define the function to minimize (difference between calculated and measured)
            fresnel_error = @(eps_r) calculate_fresnel_error(eps_r, theta_rad, fresnel_diff_squared);
            
            % Set bounds for dielectric constant (typical range for soils: 3-80)
            eps_r_min = 3;
            eps_r_max = 80;
            
            % Use fminbnd to find the dielectric constant
            options = optimset('Display', 'off', 'TolX', 1e-3);
            try
                [eps_r_est, fval] = fminbnd(fresnel_error, eps_r_min, eps_r_max, options);
                
                % Only accept reasonable solutions
                if fval < 0.1 && eps_r_est >= eps_r_min && eps_r_est <= eps_r_max
                    epoch_results = [epoch_results; prn, theta, eps_r_est, fval, power_ratio];
                end
            catch
                % Skip if optimization fails
                continue;
            end
        end
    end
    
    dielectric_results{epoch} = epoch_results;
end

% Compile all results
all_dielectric_results = [];
for epoch = 1:length(dielectric_results)
    if ~isempty(dielectric_results{epoch})
        all_dielectric_results = [all_dielectric_results; ...
            repmat(zenith_data{epoch,1}, size(dielectric_results{epoch}, 1), 1), ...
            dielectric_results{epoch}];
    end
end

% Plot results
if ~isempty(all_dielectric_results)
    figure;
    histogram(all_dielectric_results(:,3), 50);
    xlabel('Dielectric Constant');
    ylabel('Count');
    grid on;
end

%% Helper function to calculate Fresnel coefficients error
function error = calculate_fresnel_error(eps_r, theta, target_value)
    % Calculate Fresnel coefficients for horizontal and vertical polarization
    % Eqs. 7 and 8 from the paper
    
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    
    % Horizontal polarization (HH)
    R_hh = (cos_theta - sqrt(eps_r - sin_theta^2)) / ...
           (cos_theta + sqrt(eps_r - sin_theta^2));
    
    % Vertical polarization (VV)  
    R_vv = (eps_r * cos_theta - sqrt(eps_r - sin_theta^2)) / ...
           (eps_r * cos_theta + sqrt(eps_r - sin_theta^2));
    
    % Calculate the squared difference (Eq. 16)
    calculated_value = abs(R_vv - R_hh)^2;
    
    % Return the squared error
    error = (calculated_value - target_value)^2;
end

%% Optional: Convert dielectric constant to soil moisture using Topp's equation
if ~isempty(all_dielectric_results)
    % Topp's equation (empirical relationship)
    % θ_v = -5.3e-2 + 2.92e-2*ε_r - 5.5e-4*ε_r^2 + 4.3e-6*ε_r^3
    soil_moisture = -5.3e-2 + 2.92e-2 * all_dielectric_results(:,3) - ...
                    5.5e-4 * all_dielectric_results(:,3).^2 + ...
                    4.3e-6 * all_dielectric_results(:,3).^3;
    
    % Ensure soil moisture is between 0 and 1 (0-100%)
    soil_moisture = max(0, min(1, soil_moisture));
    
    figure();
    scatter(all_dielectric_results(:,1), soil_moisture * 100, 20, 'filled');
    xlabel('Time');
    ylabel('Soil Moisture Content (%)');
    title('Estimated Soil Moisture Content vs Time');
    grid on;
    
    fprintf('\nSoil Moisture Statistics:\n');
    fprintf('Mean: %.1f%%\n', mean(soil_moisture)*100);
    fprintf('Range: %.1f%% to %.1f%%\n', min(soil_moisture)*100, max(soil_moisture)*100);
end