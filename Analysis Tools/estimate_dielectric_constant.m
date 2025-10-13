function [outputArg1,outputArg2] = estimate_dielectric_constant(inputArg1,inputArg2)
% NOT WORKED ON YET, DO NOT USE
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
end

