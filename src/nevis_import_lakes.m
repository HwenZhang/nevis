%% Script to import lake data into NEVIS model
% This script imports lake data from Mechanistic Lake-Drainage Catalogue for 2022 (Stevens et al., In review)
% and integrates it into the NEVIS ice sheet hydrology model
% input: ./data/mechanistic_drainage_catalogue_2022-main/catalogues/environs_lakes_2022B_250505_archive.mat
% output: an array of lake objects that contain lake position, volume, drainage time, etc.
% Author: Hanwen Zhang  
% Date: 2026-01-20

function lakes = nevis_import_lakes(yr)

    clc
    if nargin < 1
        yr = 2022; % default year
    end
    yr_str = num2str(yr);
    % struct to hold lake data, each array corresponds to a feature
    lakes = struct();
    %% setup paths and directories
    
    oo.root = './';                                % filename root
    oo.data = '../nevis/data';                     % code directory 
    oo.dataset = ['mechanistic_drainage_catalogue_', yr_str, '-main/catalogues/environs_lakes_', yr_str, 'B_250505_archive.mat'];
                                                 % dataset name  
    oo.dn = [oo.root, 'data/datasets/nevis_2022_140km/', oo.dataset];   % catalogue under the reference dataset package
    load(oo.dn);                                   % load lake
    
    % filter out lakes which drain through moulins and hydrofractures, 1=HF, 2=moulin, 3=overflow
    overflow_inds = find(environs_lakes.drainage_type_num == 1 | environs_lakes.drainage_type_num == 2); % environs_lakes.drainage_type_num == 2

    % assign lake drainage type number
    lakes.id = environs_lakes.laketypeing_dates(overflow_inds,1);
    lakes.x_m = 1e3*environs_lakes.X_km(overflow_inds)';
    lakes.y_m = 1e3*environs_lakes.Y_km(overflow_inds)';
    lakes.drainage_type_num = environs_lakes.drainage_type_num(overflow_inds)';
    lakes.drainage_time_start = environs_lakes.laketypeing_dates(overflow_inds,3);
    lakes.drainage_time_end = environs_lakes.laketypeing_dates(overflow_inds,4);

    % [n_lakes_total, n_times] = size(environs_lakes.lake_volume_interp);
    % n_lakes_total = length(overflow_inds);
    n_lakes = length(overflow_inds);
    lakes.volume_pre_drainage_m3 = zeros(n_lakes, 1);
    lakes.volume_post_drainage_m3 = zeros(n_lakes, 1);
    lakes.max_surface_area_m2 = environs_lakes.max_surface_area(overflow_inds)';  % in m^2

    for i = 1:n_lakes
        row = overflow_inds(i);  % lake index in the original dataset
        col_start = lakes.drainage_time_start(i);  % time for drainage start, corresponds to the third column of laketypeing_dates
        col_end = lakes.drainage_time_end(i);  % time for drainage end, corresponds to the fourth column of laketypeing_dates
        % disp(['Processing lake ID: ', num2str(lakes.id(i)), ', Index: ', num2str(i), '/', num2str(n_lakes)]);

        % disp(['  Drainage time: ', num2str(col_start), ' to ', num2str(col_end)]);
        lakes.volume_pre_drainage_m3(i) = environs_lakes.lake_volume_interp(row, col_start);
        % disp(['  Pre-drainage volume (1e6 m^3): ', num2str(1e-6*lakes.volume_pre_drainage_m3(i))]);
        lakes.volume_post_drainage_m3(i) = environs_lakes.lake_volume_interp(row, col_end);
        % disp(['  Post-drainage volume (1e6 m^3): ', num2str(1e-6*lakes.volume_post_drainage_m3(i))]);
    end
        
    lakes.drainage_duration_days = lakes.drainage_time_end - lakes.drainage_time_start;
    % Override HF drainage duration to 2 hours
    hf_mask = lakes.drainage_type_num == 1;
    lakes.drainage_duration_days(hf_mask) = 2/24;
    lakes.volume_m3 = lakes.volume_pre_drainage_m3;
    % lakes.volume_m3 = lakes.volume_pre_drainage_m3 - lakes.volume_post_drainage_m3;
    % disp(lakes.volume_m3);
end