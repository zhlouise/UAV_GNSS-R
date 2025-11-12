function pr_error = cal_SD_pr_error(data_all, ground_truth)
% -------------------------------------------------------------------------
% Calculates the pseudorange error from single difference of pseudorange
% Input: 
%   data_all (n*11 cell, n = number of epochs)
%       Col-1: GPS Time
%       Col-2: Positioning Solution (Ordinary Least Squares - corrected by rtklib)
%       Col-3: Satellite PRN (1-32 GPS, 33-59 GLONASS, 60-95 Galileo, 96-105 QZSS, 106-150 Beidou)
%       Col-4: Satellite Elevation Angle
%       Col-5: Satellite Azimuth Angle (NED, North as zero)
%       Col-6,7,8: Full-ephemeris PRN, El, Az
%       Col-9: L1 C/N0, Pseudorange Residual, Carrier Phase [***some SV not available in Eph, excluded C/N0]
%       Col-10: L1 PRN, Pseudorange with receiver_bias correction, SV-XYZ-Position in ECEF for positioning
%       Col-11: L1 PRN, Doppler, SV-XYZ-Velocity in ECEF for velocity estimation
%   ground_truth: in LLH coordinates (n*3 matrix)
% Output:
%   pr_error: n*m matrix
% -------------------------------------------------------------------------

% Prompt for observation and ephemeris data
[obs_filename, obs_path] = uigetfile('.obs','Select the observation file for calculating pseudorange error.');
obs = gt.Gobs([obs_path, obs_filename]);
[eph_filename, eph_path] = uigetfile('.nav','Select the ephemeris file for calculating pseudorange error.');
nav = gt.Gnav([eph_path, eph_filename]);

% Concatenate current obs to match the time period represented in data_all
[~, idx] = ismember(round([data_all{:,1}]), round(obs.time.tow));
obs = obs.selectTime(idx);

% Select constellations
obs = obs.selectSat(obs.sys==gt.C.SYS_GPS | obs.sys==gt.C.SYS_QZS | ...
    obs.sys==gt.C.SYS_CMP | obs.sys==gt.C.SYS_GAL | obs.sys==gt.C.SYS_GLO);

% Calculate satellite positions and residuals
gt_xyz = rtklib.llh2xyz(ground_truth); % Convert ground_truth from LLH to ECEF
pos_ini = gt.Gpos(gt_xyz,"xyz"); 
sat = gt.Gsat(obs, nav);
sat.setRcvPos(pos_ini);
obs = obs.residuals(sat); % Compute residuals

corrected_pseudorange = cell(length(data_all),1);
for idt = 1:obs.n
    % Pseudorange correction
    resP = obs.L1.resPc(idt,:) - nav.getTGD(sat.sat);
    [temp_range, ~, ~, ~] = rtklib.geodist(sat.x(idt,:), sat.y(idt,:), sat.z(idt,:), pos_ini.xyz(idt,:));
    temp_PrC = resP + temp_range;
    sv_data = [obs.sat', temp_PrC', temp_range'];
    % Some SV not avaliable in Eph/pseudorange, excluded
    exclusion_bool = any([isnan(sv_data(:,1)), isnan(sv_data(:,2))],2);
    sv_data(exclusion_bool,:) = [];
    % Populate corrected_pseudorange
    corrected_pseudorange{idt,1} = sv_data;

    % Processing progress bar
    fprintf('Processed corrected pseudorange for epoch ... %d/%d\n', idt, obs.n);
end

% Separate constellations by prn
prns = unique(vertcat(data_all{:,3}));
prn_GPS = prns(prns>=1 & prns<=32);
prn_GLONASS = prns(prns>=33 & prns<=59);
prn_Galileo = prns(prns>=60 & prns<=95);
prn_QZSS = prns(prns>=96 & prns<=105);
prn_Beidou = prns(prns>=106 & prns<=150);

% Initialize pr_error matrix
pr_error = NaN(length(data_all),length(prns));

for idt = 1:length(data_all)

    % Separate data by constellation in the following structure: 
    % [prn, elevation angle, corrected pseudorange, range]

    % GPS
    [~, idx_GPS] = ismember(cell2mat(data_all(idt,3)),prn_GPS);
    GPS_data = [data_all{idt,3}(idx_GPS>0),data_all{idt,4}(idx_GPS>0),corrected_pseudorange{idt,1}(idx_GPS>0,2),...
        corrected_pseudorange{idt,1}(idx_GPS>0,3)];
    % GLONASS
    [~, idx_GLONASS] = ismember(cell2mat(data_all(idt,3)),prn_GLONASS);
    GLONASS_data = [data_all{idt,3}(idx_GLONASS>0),data_all{idt,4}(idx_GLONASS>0),corrected_pseudorange{idt,1}(idx_GLONASS>0,2),...
        corrected_pseudorange{idt,1}(idx_GLONASS>0,3)];
    % Galileo
    [~, idx_Galileo] = ismember(cell2mat(data_all(idt,3)),prn_Galileo);
    Galileo_data = [data_all{idt,3}(idx_Galileo>0),data_all{idt,4}(idx_Galileo>0),corrected_pseudorange{idt,1}(idx_Galileo>0,2),...
        corrected_pseudorange{idt,1}(idx_Galileo>0,3)];
    % QZSS
    [~, idx_QZSS] = ismember(cell2mat(data_all(idt,3)),prn_QZSS);
    QZSS_data = [data_all{idt,3}(idx_QZSS>0),data_all{idt,4}(idx_QZSS>0),corrected_pseudorange{idt,1}(idx_QZSS>0,2),...
        corrected_pseudorange{idt,1}(idx_QZSS>0,3)];
    % Beidou
    [~, idx_Beidou] = ismember(cell2mat(data_all(idt,3)),prn_Beidou);
    Beidou_data = [data_all{idt,3}(idx_Beidou>0),data_all{idt,4}(idx_Beidou>0),corrected_pseudorange{idt,1}(idx_Beidou>0,2),...
        corrected_pseudorange{idt,1}(idx_Beidou>0,3)];

    % Calculate pseudorange error
    try % Supress any error -> pr_error entry becomes NaN
        % Calculate SD pseudorange error (Master satellite from highest
        % elevation angle) 
        pr_error_GPS = (GPS_data(:,3)-GPS_data(GPS_data(:,2)==max(GPS_data(:,2)),3))...
            -(GPS_data(:,4)-GPS_data(GPS_data(:,2)==max(GPS_data(:,2)),4));
        pr_error_GLONASS = GLONASS_data(:,3) - GLONASS_data(GLONASS_data(:,2)==max(GLONASS_data(:,2)),3)...
            -(GLONASS_data(:,4)-GLONASS_data(GLONASS_data(:,2)==max(GLONASS_data(:,2)),4));
        pr_error_Galileo = Galileo_data(:,3) - Galileo_data(Galileo_data(:,2)==max(Galileo_data(:,2)),3)...
            -(Galileo_data(:,4)-Galileo_data(Galileo_data(:,2)==max(Galileo_data(:,2)),4));
        pr_error_QZSS = QZSS_data(:,3) - QZSS_data(QZSS_data(:,2)==max(QZSS_data(:,2)),3)...
            -(QZSS_data(:,4)-QZSS_data(QZSS_data(:,2)==max(QZSS_data(:,2)),4));
        pr_error_Beidou = Beidou_data(:,3) - Beidou_data(Beidou_data(:,2)==max(Beidou_data(:,2)),3)...
            -(Beidou_data(:,4)-Beidou_data(Beidou_data(:,2)==max(Beidou_data(:,2)),4));

        % Populate pr_error
        pr_error(idt,ismember(prns,data_all{idt,3}(idx_GPS>0))) = pr_error_GPS(:,1);
        pr_error(idt,ismember(prns,data_all{idt,3}(idx_GLONASS>0))) = pr_error_GLONASS(:,1);
        pr_error(idt,ismember(prns,data_all{idt,3}(idx_Galileo>0))) = pr_error_Galileo(:,1);
        pr_error(idt,ismember(prns,data_all{idt,3}(idx_QZSS>0))) = pr_error_QZSS(:,1);
        pr_error(idt,ismember(prns,data_all{idt,3}(idx_Beidou>0))) = pr_error_Beidou(:,1);
    catch
        continue;
    end

end

end
