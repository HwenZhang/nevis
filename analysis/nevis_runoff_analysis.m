%% Calculate total RACMO surface runoff and lake drainage volumes
clear; clc; close all;

%% Parameter settings
oo.root = './';
oo.code = './src';                             % code directory   
oo.dataset = 'nevis_regional';
addpath(oo.code);
oo.dn = [oo.root, 'data/', oo.dataset, '/'];

%% Load RACMO runoff data
% load([oo.dn '/' 'runoff_2010_nevis140.mat']);
load([oo.dn '/' 'runoff_2008_nevis140.mat']);
rr = runoff_2008_nevis140;

%% Load grid data
load([oo.dn '/' 'geometry']); % load BedMachine geometry (previously built)
dd = geometry; dd.skip = 6;

%% Calculate grid area
dx = abs(dd.X_m(2,1) - dd.X_m(1,1)) * dd.skip;  % grid spacing (m)
dy = abs(dd.Y_m(1,2) - dd.Y_m(1,1)) * dd.skip;  % grid spacing (m)
cell_area = dx * dy;  % single cell area (m^2)

X_skip = dd.X_m(1:dd.skip:end, 1:dd.skip:end);
Y_skip = dd.Y_m(1:dd.skip:end, 1:dd.skip:end);
S_skip = dd.S_m(1:dd.skip:end, 1:dd.skip:end);
B_skip = dd.B_m(1:dd.skip:end, 1:dd.skip:end);
H_skip = S_skip - B_skip;

% Valid ice region (non-NaN and ice thickness > 0)
valid_ice = ~isnan(H_skip) & H_skip > 0;
n_valid_cells = sum(valid_ice(:));
total_ice_area = n_valid_cells * cell_area;  % total ice area (m^2)

fprintf('========== Grid Information ==========\n');
fprintf('Grid spacing: dx=%.1f m, dy=%.1f m\n', dx, dy);
fprintf('Single cell area: %.2e m^2 (%.2f km^2)\n', cell_area, cell_area/1e6);
fprintf('Number of valid ice cells: %d\n', n_valid_cells);
fprintf('Total ice area: %.2e m^2 (%.2f km^2)\n', total_ice_area, total_ice_area/1e6);

%% Calculate annual total RACMO surface runoff
fprintf('\n========== RACMO Surface Runoff ==========\n');
fprintf('runoff_2008_nevis140 size: %s\n', mat2str(size(rr)));
% Distributed runoff [time x space] or [space x time]
[n1, n2] = size(rr);
fprintf('Distributed runoff data size: %d x %d\n', n1, n2);
n_days = n1;
seconds_per_day = 86400;

% The unit of the runoff data is mm w.e./day
daily_runoff = 1e-3*rr;  % m/day
total_runoff_per_cell = sum(daily_runoff, 1, 'omitnan');  % annual total per cell (m)
total_runoff_volume = sum(total_runoff_per_cell(:)) * cell_area;

fprintf('Annual total runoff volume: %.2e m^3 (%.2f km^3)\n', total_runoff_volume, total_runoff_volume/1e9);


%% Calculate total Lake Drainage volume
fprintf('\n========== Lake Drainage ==========\n');

lakes = nevis_import_lakes();

n_lakes = length(lakes.volume_m3);
total_lake_volume = sum(lakes.volume_m3);

fprintf('Number of lakes: %d\n', n_lakes);
fprintf('Individual lake volumes (m^3):\n');
for i = 1:n_lakes
    fprintf('  Lake %d: %.2e m^3 (%.4f km^3)\n', i, lakes.volume_m3(i), lakes.volume_m3(i)/1e9);
end
fprintf('Total lake volume: %.2e m^3 (%.4f km^3)\n', total_lake_volume, total_lake_volume/1e9);

%% Comparison
fprintf('\n========== Summary ==========\n');
if exist('total_runoff_volume', 'var')
    fprintf('RACMO Surface Runoff annual total: %.2e m^3 (%.4f km^3)\n', total_runoff_volume, total_runoff_volume/1e9);
end
fprintf('Lake Drainage total: %.2e m^3 (%.4f km^3)\n', total_lake_volume, total_lake_volume/1e9);

if exist('total_runoff_volume', 'var')
    ratio = total_lake_volume / total_runoff_volume * 100;
    fprintf('Lake Drainage / Surface Runoff = %.2f%%\n', ratio);
end