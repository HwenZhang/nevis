%% Plot the RACMO runoff over Nevis region
clear; clc; close all;

%% Parameter settings
oo.root = './';
oo.code = './src';                             % code directory   
oo.dataset = 'nevis_regional';
addpath(oo.code);
oo.dn = [oo.root, 'data/', oo.dataset, '/'];

%% Load RACMO runoff data
% load([oo.dn '/' 'runoff_2010_nevis140.mat']);
load([oo.dn '/' 'runoff_2022_nevis_noSK_300m.mat']);
load([oo.dn '/' 'runoff_2023_nevis_noSK_300m.mat']);

% data structure: runoff_2022_nevis, runoff_2023_nevis
runoff_years = {runoff_2022_nevis, runoff_2023_nevis};
year_labels = [2022, 2023];

X_m = runoff_2023_nevis.X_m;
nx = size(X_m, 1);
ny = size(X_m, 2);
Y_m = runoff_2023_nevis.Y_m;
xx = X_m(:,1); % x coordinates
yy = Y_m(1,:); % y coordinates

nd = size(runoff_2023_nevis.runoff, 1);

%% Load bed topography data for contours
load([oo.dn '/' 'morlighem_for_nevis_140km']); % load Morlighem bedmap (previously collated)
dd = morlighem_for_nevis_140km; dd.skip = 6;

dx = abs(dd.X_m(2,1) - dd.X_m(1,1));  % grid spacing (m)
dy = abs(dd.Y_m(1,2) - dd.Y_m(1,1));  % grid spacing (m)
cell_area = dx * dy;  % single cell area (m^2)

X_skip = dd.X_m(1:dd.skip:end, 1:dd.skip:end);
Y_skip = dd.Y_m(1:dd.skip:end, 1:dd.skip:end);
S_skip = dd.S_m(1:dd.skip:end, 1:dd.skip:end);
B_skip = dd.B_m(1:dd.skip:end, 1:dd.skip:end);
H_skip = S_skip - B_skip;

X_km_to_m = dd.X_m;  % convert to meters
Y_km_to_m = dd.Y_m;  % convert to meters
X_skip_m = X_skip;  % convert skip grid to meters
Y_skip_m = Y_skip;  % convert skip grid to meters

%% Process each year
for yi = 1:length(year_labels)
    year = year_labels(yi);
    runoff_data_year = runoff_years{yi};
    nd_year = size(runoff_data_year.runoff, 1);
    
    runoff_nevis140 = zeros(nd_year, size(X_skip_m,1)*size(X_skip_m,2));
    total_runoff_volume_day = zeros(nd_year,1);

    %% load the lake drainage data
    lakes = nevis_import_lakes(year);                     % import lake data
    % multiple lakes from the catalogue
    V_l = lakes.volume_m3;          % volume of lakes, set to zero for spinup       
    t_drainage = 0.5*(lakes.drainage_time_start + lakes.drainage_time_end);  % time of lake drainages (assumed to be the middle time of the Gaussian)
    t_duration = 0.5*(lakes.drainage_duration_days);  % duration of lake drainages
    lake_drainage_volume = zeros(nd_year,1);  % accumulated lake drainage volume
    accumulated_runoff = zeros(nd_year,1);    % accumulated surface runoff volume

    for d=1:nd_year
        % load runoff data for day d
        runoff_data = reshape(runoff_data_year.runoff(d,:), [nx, ny]); % mm w
        % calculate total runoff volume for this day
        total_runoff_volume_day(d) = sum(runoff_data(:), 'omitnan') * 1e-3 * cell_area; % m^3
        
        % accumulated runoff
        if d > 1
            accumulated_runoff(d) = accumulated_runoff(d-1) + total_runoff_volume_day(d);
        else
            accumulated_runoff(d) = total_runoff_volume_day(d);
        end
        
        % accumulated lake drainage volume
        if d > 1
            lake_drainage_volume(d) = lake_drainage_volume(d-1);
        end
        for li = 1:length(V_l)
            t0 = t_drainage(li);
            if t0 >= d && t0 < d+1  % drainage at d
                lake_drainage_volume(d) = lake_drainage_volume(d) + V_l(li);
            end
        end
        
        % Interpolate runoff data onto full bed topography grid
        runoff_on_topo_grid = interp2(xx, yy, runoff_data', X_km_to_m, Y_km_to_m, 'linear');
        runoff_data_interpolated = interp2(xx, yy, runoff_data', X_skip_m, Y_skip_m, 'linear');
        % make nan points zero
        nan_points = isnan(runoff_data_interpolated);
        runoff_data_interpolated(nan_points) = 0;

        runoff_nevis140(d,:)= runoff_data_interpolated(:)';
    end
    % save runoff data for this year
    varname = ['runoff_' num2str(year) '_nevis140'];
    eval([varname ' = runoff_nevis140;']); % rename runoff_nevis140 to include year
    save([oo.dn 'runoff_' num2str(year) '_nevis140.mat'], varname);

    fprintf('Year %d: Total runoff volume: %.4f km^3\n', year, accumulated_runoff(end)/1e9);
    fprintf('Year %d: Total lake drainage volume: %.4f km^3\n', year, lake_drainage_volume(end)/1e9);
    fprintf('The ratio of lake drainage to surface runoff is %.2f%%\n', 100*lake_drainage_volume(end)/accumulated_runoff(end)); 

    %% Plotting for this year
    fig = figure('Units', 'pixels', 'Position', [100, 100, 800, 600]);
    
    % 手动设置 subplot 位置，留出 title 空间
    ax1 = axes('Position', [0.1, 0.55, 0.85, 0.38]);
    plot(1:nd_year, accumulated_runoff, 'b-', 'LineWidth', 2);
    hold on
    plot(1:nd_year, lake_drainage_volume, 'r-', 'LineWidth', 2);
    h=legend('RACMO Runoff Volume', 'Lake Drainage Volume');
    h.Location = 'northwest';
    xlabel('Day of Year');
    ylabel('Total Runoff Volume (m^3)');
    title(['Total RACMO Surface Runoff Volume over Nevis Region (' num2str(year) ')']);
    xlim([150 300]);
    grid on;

    ax2 = axes('Position', [0.1, 0.08, 0.85, 0.38]);
    plot(1:nd_year, total_runoff_volume_day, 'b-', 'LineWidth', 2);
    hold on;
    lake_drainage_rate = [lake_drainage_volume(1); diff(lake_drainage_volume)];
    plot(1:nd_year, lake_drainage_rate+eps, 'r-', 'LineWidth', 2);
    h=legend('RACMO Runoff Rate', 'Lake Drainage Rate');
    h.Location = 'northwest';
    xlabel('Day of Year');
    ylabel('Runoff Rate (m^3/day)');
    title(['Daily RACMO Surface Runoff Rate over Nevis Region (' num2str(year) ')']);
    xlim([150 300]);
    grid on;

    % save as jpeg
    drawnow;
    frame = getframe(fig);
    imwrite(frame.cdata, ['nevis_runoff_' num2str(year) '.jpg'], 'Quality', 95);
end

%% Plotting
% figure;
% pcolor(X_m, Y_m, runoff_data);
% hold on
% shading interp;
% colorbar;
% contour(X_m, Y_m, runoff_data, 'k');
% xlabel('X (m)');
% ylabel('Y (m)');
% title('RACMO Surface Runoff over Nevis Region (2023)');

% plotting bed topography contours
% figure;
% for d = 1:nd
%     runoff_data = reshape(runoff_2023_nevis140(d,:), [size(X_skip_m,1), size(X_skip_m,2)]); % mm w
%     % find all nan points and put a scatter point there for better visualization
%     na_points = isnan(runoff_data);
%     nan_count = sum(na_points(:));
%     total_nan_days = 0;
%     if nan_count > 0
%         total_nan_days = total_nan_days + 1;
%         fprintf('Day %3d: found %d NaN points.\n', d, nan_count);
%     end
%     hold on;
%     pcolor(X_skip_m, Y_skip_m, runoff_data);
%     shading interp;
%     scatter(X_skip_m(na_points), Y_skip_m(na_points), 10, 'r.');
%     colorbar;
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     title(['RACMO Surface Runoff over Nevis Region (Day ' num2str(d) ')']);
%     pause(0.1);
%     clf;
% end