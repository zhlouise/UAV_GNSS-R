function [pos_est, pr_resi] = least_squares_positioning(sv_data, pos_ini)

% ------------------------------------------------------------------------
% Function to perform simple least square positioning.
% Input: 
% 1. sv_data from each epoch
%        Column 1: Satellite ID
%        Column 2: Pseudorange
%        Column 3-5: Satellite position in ECEF
% 2. pos_ini: initial positioning guess in ECEF
% Output: 
% 1. pos_est - estimated position in ECEF
% 2. pr_resi - pseudorange residual
% ------------------------------------------------------------------------

satID = sv_data(:,1);
numSat = length(satID);
pseudoranges = sv_data(:,2);
satPos = sv_data(:,3:5);

constellation = zeros(numSat, 1);
H_t = zeros(numSat, 5); % [GPS, GLO, GAL, QZSS, BDS]
for i = 1:numSat
    if satID(i) <= 32        % GPS
        constellation(i) = 1;
        H_t(i,:) = [1,0,0,0,0];
    elseif satID(i) <= 59    % GLONASS
        constellation(i) = 2;
        H_t(i,:) = [1,1,0,0,0];
    elseif satID(i) <= 95    % Galileo
        constellation(i) = 3;
        H_t(i,:) = [1,0,1,0,0];
    elseif satID(i) <= 105   % QZSS
        constellation(i) = 4;
        H_t(i,:) = [1,0,0,1,0];
    elseif satID(i) <= 150   % BeiDou
        constellation(i) = 5;
        H_t(i,:) = [1,0,0,0,1];
    end
end

% Keep only active constellations
active_constellations = any(H_t, 1);
H_t = H_t(:, active_constellations);
num_const = sum(active_constellations);

if numSat < (3 + num_const)
    disp("Available measurements less than 4, returning NaN for positioning solution.");
    pos_est = [NaN, NaN, NaN];
    pr_resi = NaN(numSat,1);
    return;
end

% Pre-allocate geometry matrix and pseudorange residuals
state_size = 3 + num_const;
geometry_matrix = zeros(numSat,state_size);
pr_resi = zeros(numSat,1);

% Accuracy requirement
accuracy_req = 1e-5;

% Initial estimates [estimated receiver location in ECEF, clockbias]
pos_estimate = [pos_ini.'; zeros(num_const,1)];

% LS iterative process
num_iter = 0;
while true
    
    % For each satellite, calculate its predicted pseudorange
    for i = 1:numSat
        pr_predicted = norm(satPos(i,:)-pos_estimate(1:3)');
        pr_measured = pseudoranges(i);
        const_bias = H_t(i,:)*pos_estimate(4:end);
        pr_resi(i) = pr_measured - (pr_predicted + const_bias);
        % Calculate geometry matrix
        geometry_matrix(i,:) = [(pos_estimate(1:3)'-satPos(i,:))/pr_predicted, H_t(i,:)];
    end
    pos_predicted = pos_estimate + (geometry_matrix'*eye(numSat)*geometry_matrix)\...
        geometry_matrix'*eye(numSat)*pr_resi;

    pos_residual = norm(pos_predicted - pos_estimate);

    if pos_residual < accuracy_req
        pos_est = pos_predicted(1:3)';
        break;
    else
        pos_estimate = pos_predicted;
    end
    
    num_iter = num_iter+1;
    % Break if cannot converge
    if num_iter > 20
        pos_est = [NaN, NaN, NaN];
        break;
    end

end

end