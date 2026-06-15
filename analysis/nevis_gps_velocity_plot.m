%% GPS velocity comparison: model vs observations per cluster
%  Each cluster gets two panels: vertical speed w and horizontal speed U
clc; clear; close all;

%% Settings
casename = 'nevis_regional_rebuilt_legacy_forward';
% casename = 'n2d_regional_racmo_eps1e_02_kappa5e_11_mu5e0_partition5e_01_k01e_01_spinup';
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
load(fullfile(repo_root, 'results', casename, casename))
if exist('region', 'var') ~= 1
    region = struct;
end
oo.fn = ['/',casename];
oo.rn = [oo.root,oo.results,oo.fn];
oo.code = fullfile(repo_root, 'src');
addpath(oo.code);

%% Cluster selection
clusters_to_plot = [1 2 3 4 5];

%% Time range
tmin_yr = 0.4;
tmax_yr = tmin_yr + 0.4;
tmin = tmin_yr * 365;
tmax = tmax_yr * 365;

%% Read time series
t = (ps.t/pd.td) * [tt.t];  % days

%% GPS data at pts_ni
pts_hs  = ps.h  * [tt.pts_hs];
pts_hb  = ps.hb * [tt.pts_hb];
pts_u   = ps.u_b * pd.ty * [tt.pts_u];   % m/yr
pts_phi = ps.phi * [tt.pts_phi];
pts_N   = ps.phi*aa.phi_0(oo.pts_ni) - pts_phi;

%% Map GPS stations to pts_ni rows
[~, gps_rows] = ismember(pp.ni_gps, oo.pts_ni);
gps_rows = gps_rows(gps_rows > 0);

%% Cluster GPS stations
cluster_centers = [0, 0; 20, -15; 40, -25];
gps_x_km = (ps.x/1e3) * gg.nx(pp.ni_gps);
gps_y_km = (ps.x/1e3) * gg.ny(pp.ni_gps);
gps_cluster = zeros(1, length(pp.ni_gps));
for ig = 1:length(pp.ni_gps)
    dists = sum((cluster_centers - [gps_x_km(ig), gps_y_km(ig)]).^2, 2);
    [~, gps_cluster(ig)] = min(dists);
end

cluster_colors = [0.85 0.33 0.10; ...   % vermillion
                  0.13 0.55 0.13; ...   % forest green
                  0.12 0.47 0.71; ...   % steel blue
                  0.58 0.00 0.83; ...   % purple
                  0.72 0.53 0.04];      % dark gold
cluster_labels = {'C1','C2','C3','C4','C5'};

% two exceptions
dist_to_c2 = sqrt((gps_x_km - 17).^2 + (gps_y_km + 25).^2);
[~, idx_c2] = min(dist_to_c2);
gps_cluster(idx_c2) = 4;
dist_to_c1 = sqrt((gps_x_km - 5).^2 + (gps_y_km + 10).^2);
[~, idx_c1] = min(dist_to_c1);
gps_cluster(idx_c1) = 5;

%% Load observations from the dataset package recorded in the result
stations = nevis_load_result_stations(region, oo);
n_stations = length(stations);

%% Compute model vertical speed d(h_b+h_s)/dt
dt_days = diff(t);
hb_hs = pts_hb + pts_hs;
dhdt_model = diff(hb_hs, 1, 2) ./ dt_days * 365;  % m/yr
t_mid = 0.5*(t(1:end-1) + t(2:end));

%% Create figure: 1 row per cluster, 2 columns (U left, w right)
n_clust = length(clusters_to_plot);
f = figure('Position', [100, 100, 1200, 250*n_clust]);
tl = tiledlayout(f, n_clust, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, ['GPS velocity comparison: model vs observations (' casename ')'], 'FontSize', 14, 'Interpreter', 'none');

for ic = 1:n_clust
    kc = clusters_to_plot(ic);
    col     = cluster_colors(kc,:);
    col_dark = col * 0.4;
    members = gps_rows(gps_cluster == kc);
    members_obs = find(gps_cluster == kc);
    label = cluster_labels{kc};

    %% --- Left panel: horizontal speed U ---
    ax_u = nexttile(tl, (ic-1)*2 + 1);
    hold on;

    % shade drainage events as light blue bands
    for i_drain = 1:length(pp.t_drainage)
        td_d = pp.t_drainage(i_drain) * ps.t / pd.td;  % center time in days
        ts_d = pp.t_duration(i_drain) * ps.t / pd.td;   % duration in days
        patch(ax_u, [td_d-0.5*ts_d, td_d+0.5*ts_d, td_d+0.5*ts_d, td_d-0.5*ts_d], ...
            [1e-4 1e-4 1e6 1e6], [0.7 0.85 1], ...
            'EdgeColor','none', 'FaceAlpha',0.4, 'HandleVisibility','off');
    end

    % plot(ax_u, NaN, NaN, '--', 'Color', col, 'LineWidth', 1.2, 'DisplayName', 'Obs');
    % plot(ax_u, NaN, NaN, '-', 'Color', col_dark, 'LineWidth', 1.2, 'DisplayName', 'Model');

    % individual station scatter
    for is = 1:n_stations
        if gps_cluster(is) ~= kc, continue; end
        scatter(ax_u, stations(is).t_doy, stations(is).u_s, 1, ...
            col, 'filled', 'MarkerFaceAlpha', 0.02, 'HandleVisibility','off');
    end

    % model speed
    if ~isempty(members)
        plot(ax_u, t, mean(pts_u(members,:),1), '-', ...
            'Color', col_dark, 'LineWidth', 1.2, 'DisplayName', 'Model');
    end

    % cluster-mean observed speed
    if ~isempty(members_obs)
        t_all = vertcat(stations(members_obs).t_doy);
        u_all = vertcat(stations(members_obs).u_s);
        valid = ~isnan(t_all) & ~isnan(u_all);
        t_all = t_all(valid); u_all = u_all(valid);
        [t_sorted, si] = sort(t_all);
        u_sorted = u_all(si);
        u_smooth = movmean(u_sorted, max(1, min(6, length(u_sorted))), 'omitnan');
        plot(ax_u, t_sorted, u_smooth, '-', 'Color', col, 'LineWidth', 1.2, 'DisplayName', 'Obs');
    end

    ylabel('U [m/yr]');
    if ic == n_clust, xlabel('t [d]'); end
    title([label ': horizontal speed']);
    xlim([tmin tmax]); ylim([0 300]);
    grid on; grid minor;
    legend(ax_u, 'Location','northeast', 'FontSize', 8);

    %% --- Right panel: vertical speed w ---
    ax_w = nexttile(tl, (ic-1)*2 + 2);
    hold on;

    % shade drainage events as light blue bands
    for i_drain = 1:length(pp.t_drainage)
        td_d = pp.t_drainage(i_drain) * ps.t / pd.td;  % center time in days
        ts_d = pp.t_duration(i_drain) * ps.t / pd.td;   % duration in days
        patch(ax_w, [td_d-0.5*ts_d, td_d+0.5*ts_d, td_d+0.5*ts_d, td_d-0.5*ts_d], ...
            [-1e6 -1e6 1e6 1e6], [0.7 0.85 1], ...
            'EdgeColor','none', 'FaceAlpha',0.4, 'HandleVisibility','off');
    end

    % plot(ax_w, NaN, NaN, '--', 'Color', col, 'LineWidth', 1.2, 'DisplayName', 'Obs');
    % plot(ax_w, NaN, NaN, '-', 'Color', col_dark, 'LineWidth', 1.2, 'DisplayName', 'Model');
    
    % individual station scatter
    for is = 1:n_stations
        if gps_cluster(is) ~= kc, continue; end
        scatter(ax_w, stations(is).t_doy, stations(is).w_s, 1, ...
            col, 'filled', 'MarkerFaceAlpha', 0.02, 'HandleVisibility','off');
    end
    
    % model d(h_b+h_s)/dt
    if ~isempty(members)
        plot(ax_w, t_mid, mean(dhdt_model(members,:),1), '-', ...
            'Color', col_dark, 'LineWidth', 1.2, 'DisplayName', 'Model');
    end

    % cluster-mean observed w_s
    t_all = []; w_all = [];
    for im = members_obs
        if isfield(stations(im), 'w_s') && ~isempty(stations(im).w_s)
            t_all = [t_all; stations(im).t_doy(:)];
            w_all = [w_all; stations(im).w_s(:)];
        end
    end
    if ~isempty(t_all)
        valid = ~isnan(t_all) & ~isnan(w_all);
        t_all = t_all(valid); w_all = w_all(valid);
        [t_sorted, si] = sort(t_all);
        w_sorted = w_all(si);
        w_smooth = movmean(w_sorted, max(1, round(length(w_sorted)/100)), 'omitnan');
        plot(ax_w, t_sorted, w_smooth, '-', 'Color', col, 'LineWidth', 1.2, 'DisplayName', 'Obs');
    end

    ylabel('w [m/yr]');
    if ic == n_clust, xlabel('t [d]'); end
    title([label ': vertical speed']);
    xlim([tmin tmax]); ylim([-100 100]);
    grid on; grid minor;
    legend(ax_w, 'Location','northeast', 'FontSize', 8);
end

saveas(f, fullfile(repo_root, [casename '_gps_comparison.png']));
