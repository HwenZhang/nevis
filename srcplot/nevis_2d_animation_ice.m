%% Animation script for NEVIS 2D ice-hydrology coupled model
% Produces a video with time series (left) and spatial fields (right):
%   Left:  Q fluxes, effective pressure, sheet/channel, ice speed
%   Right: net flux, cavity thickness, effective pressure, 
%          channel cross-section, ice speed, ice speed change
clc; clear; close all;

%% Settings
casename = 'nevis_2022_140km_botright30_drainage';
% casename = 'n2d_regional_racmo_eps1e_02_kappa5e_11_mu5e0_partition5e_01_k01e_01_spinup';
script_dir = fileparts(mfilename('fullpath'));
repo_root = script_dir;
result_dir = fullfile(repo_root, 'results', casename);
result_file = fullfile(result_dir, [casename '.mat']);
if exist(result_file, 'file') ~= 2
    error('nevis_2d_animation_ice:MissingResultFile', ...
        'Missing result file: %s', result_file);
end
load(result_file)
if exist('region', 'var') ~= 1
    region = struct;
end
oo.fn = ['/',casename];
oo.rn = result_dir;
oo.code = fullfile(repo_root, 'src');
filepath = [oo.rn filesep];
addpath(oo.code);
formatSpec = '%04d';

%% Colormap
n = 256;
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% Cluster selection (which clusters to plot; e.g. [4 5] or 1:5 for all)
clusters_to_plot = [5];

%% Time range
t = (ps.t/(24*60*60))*[tt.t];               % dimensional time series (days)
tspan = (ps.t/pd.td)*oo.t_span;
tmin = 0.4*365*pd.td/ps.t;
tmax = 1.0*365*pd.td/ps.t;
tmin_d = tmin*ps.t/pd.td;
tmax_d = tmax*ps.t/pd.td;                   % time range for the plot

% Controllable time step: pick frames on a regular dt_plot grid and snap each
% to the nearest saved timestep. Set dt_plot smaller for a slower animation.
dt_plot = 1.0; % use 1.0 for fast animation, 0.1 for slow ones
t_plot = tmin_d:dt_plot:tmax_d;

frame_indices = zeros(size(t_plot));
for k = 1:length(t_plot)
    [~, idx] = min(abs(tspan - t_plot(k)));
    frame_indices(k) = idx;
end
frame_indices = unique(frame_indices);

%% Read time series data
t = (ps.t/pd.td) * [tt.t];  % dimensional time (days)

Q_b_in  = pd.Q_0 * [tt.Qb_in];
Q_b_dec = ps.h*ps.x^2/ps.t * [tt.Qb_dec];
Q_in    = ps.Q * [tt.Q_in];
Q_out   = ps.Q * ([tt.Q_out]-pp.c6*[tt.E]);
Q_out_Q = ps.Q * [tt.Q_outQ];
Q_out_q = ps.Q * [tt.Q_outq];
Q_out_b2 = ps.Q0 * [tt.Q_outb];
E       = (ps.m*ps.x^2) * [tt.E];
N_ts    = (ps.phi) * [tt.N];  % effective pressure time series
hs_ts   = ps.x^2*ps.h * [tt.hs];
he_ts   = ps.x^2*ps.h * [tt.he];
S_ts    = ps.x*ps.S * [tt.S];
A_total = ps.x^2 * sum(gg.Dx.*gg.Dy);
U       = ps.u_b * [tt.U];

h_b = ps.hb*[tt.pts_hb];       %
p_b = ps.phi*[tt.pts_pb];      %
V_b = ps.x^2*ps.hb*[tt.Vb];

%% Read initial frame
nframe = frame_indices(1);
vva = load([filepath num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ);

xx = (ps.x/1e3)*gg.nx;  % x grid in km
yy = (ps.x/1e3)*gg.ny;

%% Compute initial ice velocity
uxn = gg.nmeanx2(:,gg.es2)*vva.u(gg.es2);
vyn = gg.nmeany2(:,gg.fs2)*vva.v(gg.fs2);
Un  = sqrt(uxn.^2 + vyn.^2);
Un(gg.nout) = NaN;

vux0 = ps.u_b*pd.ty * reshape(uxn,gg.nI,gg.nJ);
vuy0 = ps.u_b*pd.ty * reshape(vyn,gg.nI,gg.nJ);
vUU0 = sqrt(vux0.^2 + vuy0.^2);

%% Create figure layout
f = figure('Position', [50, 50, 1800, 900]);
mainLayout = tiledlayout(f, 1, 5);
mainLayout.TileSpacing = 'compact';
mainLayout.Padding = 'compact';

% Left panel: time series (7 rows)
leftLayout = tiledlayout(mainLayout, 7, 1);
leftLayout.Layout.Tile = 1;
leftLayout.Layout.TileSpan = [1, 2];
leftLayout.TileSpacing = 'compact';
leftLayout.Padding = 'compact';

% Right panel: spatial fields (4 rows x 2 cols)
rightLayout = tiledlayout(mainLayout, 4, 2);
rightLayout.Layout.Tile = 3;
rightLayout.Layout.TileSpan = [1, 3];
rightLayout.TileSpacing = 'compact';
rightLayout.Padding = 'compact';

%% ===== Left layout: time series =====
tframe_d = vva.t*ps.t/pd.td;

% (a) Fluxes
ax_a = nexttile(leftLayout);
hold on;
% shade drainage events as light blue bands
for i_drain = 1:length(pp.t_drainage)
    td_d = pp.t_drainage(i_drain) * ps.t / pd.td;  % center time in days
    ts_d = pp.t_duration(i_drain) * ps.t / pd.td;   % duration in days
    patch(ax_a, [td_d-0.5*ts_d, td_d+0.5*ts_d, td_d+0.5*ts_d, td_d-0.5*ts_d], ...
          [1e-4 1e-4 1e6 1e6], [0.7 0.85 1], ...
          'EdgeColor','none', 'FaceAlpha',0.4, 'HandleVisibility','off');
end
plot(ax_a, t, Q_out_Q + Q_out_q,'-', 'Color',[0,0.25,0], 'LineWidth',2.5);
plot(ax_a, t, Q_out_b2,'b--', 'LineWidth',1.5);
plot(ax_a, t, E,'k-.', 'LineWidth',1.5);
x1 = xline(tframe_d,'k--','LineWidth',1.5);
xlabel('t [d]'); ylabel('Q [m^3/s]');
legend('Q_{out}','Q_{out,b}','Q_{in}','NumColumns',2,'Location','northeast');
text(0.025,0.85,'(a) flux','Units','normalized','FontSize',12);
xlim([tmin_d tmax_d]); set(gca,'YScale','log'); ylim([1e-2 1e4]);
grid on; grid minor;

% (b) Effective pressure
ax_b = nexttile(leftLayout);
plot(ax_b, t, N_ts/1e6,'-', 'LineWidth',1.5);
hold on;
x2 = xline(tframe_d,'--k','LineWidth',1.5);
xlabel('t [d]'); ylabel('N [MPa]');
legend('averaged N','Location','northeast');
text(0.025,0.85,'(b) effective pressure','Units','normalized','FontSize',12);
xlim([tmin_d tmax_d]); 
ylim([0.5 2.5]);
grid on; grid minor;

% (c) Sheet thickness & channel
ax_c = nexttile(leftLayout);
yyaxis(ax_c,'left');
plot(ax_c, t, hs_ts./A_total,'b-', t, he_ts./A_total,'k-', 'LineWidth',1.5);
ylabel('Average h [m]');
yyaxis(ax_c,'right');
plot(ax_c, t, S_ts/A_total,'r-', 'LineWidth',1.5);
hold on; x3 = xline(tframe_d,'--k','LineWidth',1.5);
ylabel('S/A [m]');
legend('h_{cav}','h_e','S','NumColumns',3,'Location','northeast');
text(0.025,0.85,'(c) averaged h and S','Units','normalized','FontSize',12);
xlim([tmin_d tmax_d]); grid on; grid minor;

% (d) Ice speed time series (if available)
ax_d = nexttile(leftLayout);
yyaxis left
plot(ax_d, t, U*pd.ty, 'b-', 'LineWidth',1.5);
x4 = xline(tframe_d,'--k','LineWidth',1.5);
ylabel('U [m/yr]');
text(0.025,0.85,'(d) blister volume and ice speed','Units','normalized','FontSize',12);
xlim([tmin_d tmax_d]); 
ylim([50 150]); 
grid on; grid minor;

yyaxis right
plot(ax_d, t, V_b, 'r-', 'LineWidth',1.5);
ylabel('Vb [m^3]');

%% GPS cluster setup (shared by panels e and f)
pts_hs = ps.h * [tt.pts_hs];  % dimensional hs at all pts_ni
pts_hb = ps.hb * [tt.pts_hb];  % dimensional hb at all pts_ni
pts_u = ps.u_b * pd.ty * [tt.pts_u];  % dimensional speed at pts_ni (m/yr)
pts_phi = ps.phi * [tt.pts_phi];  % dimensional phi at pts_ni (Pa)
pts_N = (ps.phi*aa.phi_0(oo.pts_ni) - pts_phi);  % effective pressure at pts_ni (MPa)

% find which rows of pts_hb correspond to GPS stations
[~, gps_rows] = ismember(pp.ni_gps, oo.pts_ni);
gps_rows = gps_rows(gps_rows > 0);  % keep only matched indices
% cluster GPS stations by proximity to 3 centers (in km)
cluster_centers = [0, 0; 20, -15; 40, -25];  % [x, y] in km
gps_x_km = (ps.x/1e3) * gg.nx(pp.ni_gps);
gps_y_km = (ps.x/1e3) * gg.ny(pp.ni_gps);
gps_cluster = zeros(1, length(pp.ni_gps));
for ig = 1:length(pp.ni_gps)
    dists = sum((cluster_centers - [gps_x_km(ig), gps_y_km(ig)]).^2, 2);
    [~, gps_cluster(ig)] = min(dists);
end

% TODO: if some clusters have no members, reassign the nearest GPS to that cluster
cluster_colors = [0.85 0.33 0.10; ...   % vermillion
                  0.13 0.55 0.13; ...   % forest green
                  0.12 0.47 0.71; ...   % steel blue
                  0.58 0.00 0.83; ...   % purple
                  0.72 0.53 0.04];      % dark gold
cluster_labels = {'C1','C2','C3','C4','C5'};
% two exceptions:
% find the nearest to (17, -25)
dist_to_c2 = sqrt((gps_x_km - 17).^2 + (gps_y_km + 25).^2);
[~, idx_c2] = min(dist_to_c2);
gps_cluster(idx_c2) = 4;  % assign to cluster 4
% find the nearest to (5, -10)
dist_to_c1 = sqrt((gps_x_km - 5).^2 + (gps_y_km + 10).^2);
[~, idx_c1] = min(dist_to_c1);
gps_cluster(idx_c1) = 5;  % assign to cluster 5

% Load observations from the dataset package recorded in the result.
stations = nevis_load_result_stations(region, oo);
n_stations = length(stations);

%% (e) GPS station blister thicknesses h_b + h_s
ax_e = nexttile(leftLayout);
hold on;
yyaxis left
for kc = clusters_to_plot
    members = gps_rows(gps_cluster == kc);
    if isempty(members), continue; end
    plot(ax_e, t, mean(pts_hb(members,:)+pts_hs(members,:),1), '-', ...
        'Color', cluster_colors(kc,:), 'LineWidth', 1.5, ...
        'DisplayName', [cluster_labels{kc}]);
end
ylabel('h_b + h_s [m]');
x5 = xline(tframe_d,'--k','LineWidth',1.5);
text(0.025,0.85,'(e) h_b+h_s & N at GPS','Units','normalized','FontSize',12);
xlim([tmin_d tmax_d]);

yyaxis right
for kc = clusters_to_plot
    members = gps_rows(gps_cluster == kc);
    if isempty(members), continue; end
    plot(ax_e, t, mean(pts_N(members,:),1)/1e6, '--', ...
        'Color', cluster_colors(kc,:), 'LineWidth', 1.5, ...
        'DisplayName', [cluster_labels{kc} ' N']);
end
ylabel('N [MPa]');

grid on; grid minor;
legend(ax_e, 'Location', 'southwest', 'NumColumns', 2, 'FontSize', 7);

% %% (f) GPS station surface vertical speed w (model + obs)
% ax_f = nexttile(leftLayout);
% hold on;
% % compute model d(h_b+h_s)/dt in m/yr via central differences
% dt_days = diff(t);  % time steps in days
% hb_hs = pts_hb + pts_hs;  % (n_gps x n_t) dimensional h_b + h_s
% dhdt_model = diff(hb_hs, 1, 2) ./ dt_days * 365;  % m/yr
% t_mid = 0.5*(t(1:end-1) + t(2:end));  % midpoint times
% % dummy handles for legend grouping
% plot(ax_f, NaN, NaN, '-', 'Color', [0.6 0.6 0.6], ...
%     'LineWidth', 0.8, 'DisplayName', 'Obs (light)');
% plot(ax_f, NaN, NaN, '-', 'Color', [0.2 0.2 0.2], ...
%     'LineWidth', 1.2, 'DisplayName', 'Model (dark)');
% % observation scatter w_s
% % for is = 1:n_stations
% %     kc = gps_cluster(is);
% %     if ~ismember(kc, clusters_to_plot), continue; end
% %     if isfield(stations(is), 'w_s') && ~isempty(stations(is).w_s)
% %         scatter(ax_f, stations(is).t_doy, stations(is).w_s, 8, ...
% %             cluster_colors(kc,:), 'filled', 'MarkerFaceAlpha', 0.3, ...
% %             'HandleVisibility', 'off');
% %     end
% % end
% % overlay cluster-mean observed w_s
% for kc = clusters_to_plot
%     members_obs = find(gps_cluster == kc);
%     if isempty(members_obs), continue; end
%     t_all = []; w_all = [];
%     for im = members_obs
%         if isfield(stations(im), 'w_s') && ~isempty(stations(im).w_s)
%             t_all = [t_all; stations(im).t_doy(:)];
%             w_all = [w_all; stations(im).w_s(:)];
%         end
%     end
%     if isempty(t_all), continue; end
%     [t_sorted, si] = sort(t_all);
%     w_sorted = w_all(si);
%     if length(t_sorted) > 1
%         w_smooth = movmean(w_sorted, 1);
%     else
%         w_smooth = w_sorted;
%     end
%     plot(ax_f, t_sorted, w_smooth, '-', 'Color', cluster_colors(kc,:), ...
%         'LineWidth', 0.8, ...
%         'DisplayName', [cluster_labels{kc} ' obs']);
% end
% % model d(h_b+h_s)/dt (darkened cluster color, solid thick)
% for kc = clusters_to_plot
%     members = gps_rows(gps_cluster == kc);
%     if isempty(members), continue; end
%     dark_col = cluster_colors(kc,:) * 0.55;
%     plot(ax_f, t_mid, mean(dhdt_model(members,:),1), '-', ...
%         'Color', dark_col, 'LineWidth', 1.0, ...
%         'DisplayName', [cluster_labels{kc} ' model']);
% end
% ylabel('w [m/yr]');
% x6 = xline(tframe_d,'--k','LineWidth',1.5);
% text(0.025,0.85,'(f) vertical speed at GPS','Units','normalized','FontSize',12);
% xlim([tmin_d tmax_d]);
% ylim([-100 100]);
% grid on; grid minor;
% legend(ax_f, 'Location', 'southwest', 'NumColumns', 2, 'FontSize', 7);

% %% (g) GPS station ice speed U (model + obs)
% ax_g = nexttile(leftLayout);
% hold on;
% % dummy handles for legend grouping
% h_obs_dummy = plot(ax_g, NaN, NaN, '-', 'Color', [0.6 0.6 0.6], ...
%     'LineWidth', 0.8, 'DisplayName', 'Obs (light)');
% h_mod_dummy = plot(ax_g, NaN, NaN, '-', 'Color', [0.2 0.2 0.2], ...
%     'LineWidth', 1.2, 'DisplayName', 'Model (dark)');
% % observation scatter
% for is = 1:n_stations
%     kc = gps_cluster(is);
%     if ~ismember(kc, clusters_to_plot), continue; end
%     scatter(ax_g, stations(is).t_doy, stations(is).u_s, 8, ...
%         cluster_colors(kc,:), 'filled', 'MarkerFaceAlpha', 0.3, ...
%         'HandleVisibility', 'off');
% end
% % overlay cluster-mean observed speed
% for kc = clusters_to_plot
%     members_obs = find(gps_cluster == kc);
%     if isempty(members_obs), continue; end
%     t_all = vertcat(stations(members_obs).t_doy);
%     u_all = vertcat(stations(members_obs).u_s);
%     [t_sorted, si] = sort(t_all);
%     u_sorted = u_all(si);
%     if length(t_sorted) > 6
%         u_smooth = movmean(u_sorted, 6);
%     else
%         u_smooth = u_sorted;
%     end
%     plot(ax_g, t_sorted, u_smooth, '-', 'Color', cluster_colors(kc,:), ...
%         'LineWidth', 0.8, ...
%         'DisplayName', [cluster_labels{kc} ' obs']);
% end
% % model speed (darkened cluster color, solid thick, plotted last = on top)
% for kc = clusters_to_plot
%     members = gps_rows(gps_cluster == kc);
%     if isempty(members), continue; end
%     dark_col = cluster_colors(kc,:) * 0.55;  % darken for contrast
%     plot(ax_g, t, mean(pts_u(members,:),1), '-', ...
%         'Color', dark_col, 'LineWidth', 0.8, ...
%         'DisplayName', [cluster_labels{kc} ' model']);
% end
% ylabel('U [m/yr]');
% x7 = xline(tframe_d,'--k','LineWidth',1.5);
% xlabel('t [d]');
% text(0.025,0.85,'(g) U at GPS','Units','normalized','FontSize',12);
% xlim([tmin_d tmax_d]);
% ylim([0 300]);
% grid on; grid minor;
% legend(ax_g, 'Location', 'southwest', 'NumColumns', 2, 'FontSize', 7);

xlines = [x1, x2, x3, x4, x5];

%% ===== Right layout: spatial fields =====

% (1) Blister sheet thickness hb
ax = nexttile(rightLayout);
zhe = (ps.hb)*reshape(vva.hb,gg.nI,gg.nJ); 
pblister = pcolor(ax,xx,yy,zhe); 
hold on
% plot all GPS as small magenta circles
% plot(ax,xx(pp.ni_gps),yy(pp.ni_gps),'o', 'MarkerSize',3,'MarkerFaceColor',"m", 'MarkerEdgeColor', 'k'); % mark all GPS
% highlight GPS stations colored by cluster
for kc = clusters_to_plot
    members_kc = find(gps_cluster == kc);
    for km = 1:length(members_kc)
        plot(ax, xx(pp.ni_gps(members_kc(km))), yy(pp.ni_gps(members_kc(km))), ...
            's', 'MarkerSize', 4, 'MarkerFaceColor', cluster_colors(kc,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
    end
end
% add an annotation for the lake
% text(ax,xx(oo.pts_ni(end))+1,yy(oo.pts_ni(end))-1,'GPS','FontSize',10,'Color','k');

set(pblister,'linestyle','none'); % shading interp
cx = colorbar();
colormap(ax,cmap)
cx.Label.String = 'h_b [ m ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([-1.0 1.0]);
hold on

% lake input marker
time = vva.t;
xl = pp.x_l*(ps.x/10^3);
yl = pp.y_l*(ps.x/10^3);
Q_peak = max((ps.Q0)*pp.V_l ./ (sqrt(2*pi)*pp.t_drainage));  % peak flux for normalization
markersize_l = lake_input(time,pp,ps);
norm_size = markersize_l / (Q_peak + eps);          % normalise to [0,~1]
sz = 5 + 295 * norm_size.^0.3;                     % power-law: small inputs stay visible
hscatter = scatter(ax, xl, yl, sz, 'ko', 'filled', 'MarkerFaceColor',[0.2 0.5 0.9], 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', 'k'); % mark lake location

zpb = (ps.phi)*reshape(vva.pb,gg.nI,gg.nJ); 
[Cb,pblister_contour] = contour(ax,xx,yy,zpb,'linecolor','k','linewidth',0.5);

title('blister sheet and pb contour');
ylabel('y (km)')
axis equal
set(gca, 'XLimMode','manual', 'YLimMode','manual');

% (2) Cavity sheet thickness
zphi = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ); 
ax2 = nexttile(rightLayout);
vva.hs(gg.nout) = NaN;
zhs = ps.hs * reshape(vva.hs,gg.nI,gg.nJ);
phs = pcolor(ax2, xx, yy, zhs);
set(phs,'linestyle','none');
cx = colorbar();
cx.Label.String = 'h_s [m]';
clim([0 0.15]);
hold on;
[~, phs_contour] = contour(ax2, xx, yy, zphi, 'linecolor','k','linewidth',0.5);
title('(g) cavity sheet'); ylabel('y (km)');
axis equal; axis tight;
set(gca, 'XLimMode','manual', 'YLimMode','manual');

% (3) Effective pressure
ax3 = nexttile(rightLayout);
zN = (ps.phi/1e6) * reshape(aa.phi_0 - vva.phi, gg.nI, gg.nJ);
peff = pcolor(ax3, xx, yy, zN);
set(peff,'linestyle','none');
cx = colorbar(); colormap(ax3, cmap);
cx.Label.String = 'N [MPa]';
clim([-3 3]);
title('(h) effective pressure'); ylabel('y (km)');
axis equal; axis tight;
set(gca, 'XLimMode','manual', 'YLimMode','manual');

% (4) Channel cross-section area
ax4 = nexttile(rightLayout);
zS = ps.S * reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr), gg.nI, gg.nJ);
pS = pcolor(ax4, xx, yy, zS);
set(pS,'linestyle','none');
ax4.ColorScale = 'log';
cx = colorbar();
cx.Label.String = 'S [m^2]';
clim([1e-4 1e1]);
cx.Ticks = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
title('(i) channel cross section'); ylabel('y (km)');
axis equal; axis tight;
set(gca, 'XLimMode','manual', 'YLimMode','manual');

% (5) Ice speed
ax5 = nexttile(rightLayout);
pvel = pcolor(ax5, xx, yy, vUU0);
set(pvel,'linestyle','none');
cx = colorbar();
cx.Label.String = 'U [m/yr]';
clim([0 200]);
hold on;
skip = 5;  % quiver skip for clarity
vel_field = quiver(ax5, xx(1:skip:end,1:skip:end), yy(1:skip:end,1:skip:end), ...
    vux0(1:skip:end,1:skip:end), vuy0(1:skip:end,1:skip:end), 0.3, 'k');
title('(j) ice speed'); ylabel('y (km)');
axis equal; axis tight;
set(gca, 'XLimMode','manual', 'YLimMode','manual');

% (6) Ice speed change
ax6 = nexttile(rightLayout);
pvel_delta = pcolor(ax6, xx, yy, zeros(gg.nI,gg.nJ));
set(pvel_delta,'linestyle','none');
cx = colorbar(); colormap(ax6, cmap);
% cx.Label.String = '\Delta U [m/yr]';
cx.Label.String = '%';
clim([-200 200]);
ttext = text(ax6, 0.5, 0.9, ['t = ' num2str(tframe_d,'%.1f') ' d'], ...
    'Units','normalized','FontSize',14,'FontWeight','bold');
title('(k) ice speed change %'); ylabel('y (km)'); xlabel('x (km)');
axis equal; axis tight;
set(gca, 'XLimMode','manual', 'YLimMode','manual');

% (7) Plot principal stress magnitudes and directions as crosses
ax7 = nexttile(rightLayout);
stress_scale = ps.eta*ps.u/ps.x / 1e3; % convert to kPa
[tauxx_s,tauyy_s,tauxy_s,~,~,~,~,~,~,sigma1_s,sigma2_s,t1_s,~] = nevis_stresses(aa.H,vva.u,vva.v,aa.phi_0-vva.phi,aa,pp,gg,oo);
sigma1_kPa = sigma1_s * stress_scale;
sigma1_kPa(gg.nout) = NaN;
p_tau = pcolor(ax7, xx, yy, reshape(sigma1_kPa, gg.nI, gg.nJ));
shading interp;
set(p_tau,'linestyle','none'); colormap(ax7, jet);
cx = colorbar(); cx.Label.String = '\sigma_1 [kPa]';
clim([0 200]);
hold on;
skip7 = 8;
ii7 = 1:skip7:gg.nI; jj7 = 1:skip7:gg.nJ;
xq7 = xx(ii7,jj7); yq7 = yy(ii7,jj7);
s1_s = stress_scale * reshape(sigma1_s, gg.nI, gg.nJ); s1_q = s1_s(ii7,jj7);
s2_s = stress_scale * reshape(sigma2_s, gg.nI, gg.nJ); s2_q = s2_s(ii7,jj7);
dx1_s = reshape(cos(t1_s), gg.nI, gg.nJ); dx1_q = dx1_s(ii7,jj7);
dy1_s = reshape(sin(t1_s), gg.nI, gg.nJ); dy1_q = dy1_s(ii7,jj7);
dx2_s = reshape(cos(t1_s+pi/2), gg.nI, gg.nJ); dx2_q = dx2_s(ii7,jj7);
dy2_s = reshape(sin(t1_s+pi/2), gg.nI, gg.nJ); dy2_q = dy2_s(ii7,jj7);
q_s1p = quiver(ax7, xq7, yq7,  s1_q.*dx1_q,  s1_q.*dy1_q, 0.5, 'r', 'LineWidth', 1);
q_s1n = quiver(ax7, xq7, yq7, -s1_q.*dx1_q, -s1_q.*dy1_q, 0.5, 'r', 'LineWidth', 1, 'ShowArrowHead', 'off');
q_s2p = quiver(ax7, xq7, yq7,  s2_q.*dx2_q,  s2_q.*dy2_q, 0.5, 'b', 'LineWidth', 1);
q_s2n = quiver(ax7, xq7, yq7, -s2_q.*dx2_q, -s2_q.*dy2_q, 0.5, 'b', 'LineWidth', 1, 'ShowArrowHead', 'off');
title('(l) \sigma_1 + principal stresses'); ylabel('y (km)'); xlabel('x (km)');
axis equal; axis tight;
set(gca, 'XLimMode','manual', 'YLimMode','manual');

% Save initial stress state for computing change
sigma1_0 = sigma1_s;
sigma2_0 = sigma2_s;
t1_0 = t1_s;

% (8) Plot principal stress change
ax8 = nexttile(rightLayout);
dsigma1_kPa = zeros(size(sigma1_s)) * stress_scale;
dsigma1_kPa(gg.nout) = NaN;
p_dtau = pcolor(ax8, xx, yy, reshape(dsigma1_kPa, gg.nI, gg.nJ));
shading interp;
set(p_dtau,'linestyle','none'); colormap(ax8, cmap);
cx = colorbar(); cx.Label.String = '\Delta\sigma_1 [kPa]';
clim([-50 50]);
hold on;
skip8 = skip7;
ii8 = ii7; jj8 = jj7;
xq8 = xx(ii8,jj8); yq8 = yy(ii8,jj8);
ds1_q = zeros(size(s1_q)); ds2_q = zeros(size(s2_q));
dx1_q0 = dx1_q; dy1_q0 = dy1_q;
dx2_q0 = dx2_q; dy2_q0 = dy2_q;
% q_ds1p = quiver(ax8, xq8, yq8,  ds1_q.*dx1_q0,  ds1_q.*dy1_q0, 0.5, 'r', 'LineWidth', 1);
% q_ds1n = quiver(ax8, xq8, yq8, -ds1_q.*dx1_q0, -ds1_q.*dy1_q0, 0.5, 'r', 'LineWidth', 1, 'ShowArrowHead', 'off');
% q_ds2p = quiver(ax8, xq8, yq8,  ds2_q.*dx2_q0,  ds2_q.*dy2_q0, 0.5, 'b', 'LineWidth', 1);
% q_ds2n = quiver(ax8, xq8, yq8, -ds2_q.*dx2_q0, -ds2_q.*dy2_q0, 0.5, 'b', 'LineWidth', 1, 'ShowArrowHead', 'off');
title('(m) \Delta\sigma_1 + principal stress change'); ylabel('y (km)'); xlabel('x (km)');
axis equal; axis tight;
set(gca, 'XLimMode','manual', 'YLimMode','manual');

%% ===== Animation loop =====
video_dir = fullfile(repo_root, 'results', 'videos');
if exist(video_dir, 'dir') ~= 7
    mkdir(video_dir);
end
v = VideoWriter(fullfile(video_dir, casename), 'MPEG-4');
v.FrameRate = 3;
open(v);

for i_idx = 1:length(frame_indices)
    if ~isvalid(f)
        disp('Figure closed. Animation stopped.');
        break
    end
    i_t = frame_indices(i_idx);
    fprintf('Frame %d / %d ...\n', i_idx, length(frame_indices));

    % Load timestep
    vva = load([filepath num2str(i_t,formatSpec)], 'vv');
    vva = vva.vv;
    aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
    [vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);
    qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ);

    % Ice velocity
    uxn = gg.nmeanx2(:,gg.es2)*vva.u(gg.es2);
    vyn = gg.nmeany2(:,gg.fs2)*vva.v(gg.fs2);
    Un  = sqrt(uxn.^2 + vyn.^2);
    Un(gg.nout) = NaN;
    vux = ps.u_b*pd.ty * reshape(uxn,gg.nI,gg.nJ);
    vuy = ps.u_b*pd.ty * reshape(vyn,gg.nI,gg.nJ);
    vUU = sqrt(vux.^2 + vuy.^2);

    % Compute stresses for this frame
    [~,~,~,~,~,~,~,~,~,sigma1_s,sigma2_s,t1_s,~] = nevis_stresses(aa.H,vva.u,vva.v,aa.phi_0-vva.phi,aa,pp,gg,oo);
    sigma1_kPa = sigma1_s * stress_scale;
    sigma1_kPa(gg.nout) = NaN;
    p_tau.CData = reshape(sigma1_kPa, gg.nI, gg.nJ);

    % Update principal stress quivers
    s1_s = stress_scale * reshape(sigma1_s, gg.nI, gg.nJ); s1_q = s1_s(ii7,jj7);
    s2_s = stress_scale * reshape(sigma2_s, gg.nI, gg.nJ); s2_q = s2_s(ii7,jj7);
    dx1_s = reshape(cos(t1_s), gg.nI, gg.nJ); dx1_q = dx1_s(ii7,jj7);
    dy1_s = reshape(sin(t1_s), gg.nI, gg.nJ); dy1_q = dy1_s(ii7,jj7);
    dx2_s = reshape(cos(t1_s+pi/2), gg.nI, gg.nJ); dx2_q = dx2_s(ii7,jj7);
    dy2_s = reshape(sin(t1_s+pi/2), gg.nI, gg.nJ); dy2_q = dy2_s(ii7,jj7);
    q_s1p.UData =  s1_q.*dx1_q; q_s1p.VData =  s1_q.*dy1_q;
    q_s1n.UData = -s1_q.*dx1_q; q_s1n.VData = -s1_q.*dy1_q;
    q_s2p.UData =  s2_q.*dx2_q; q_s2p.VData =  s2_q.*dy2_q;
    q_s2n.UData = -s2_q.*dx2_q; q_s2n.VData = -s2_q.*dy2_q;

    % Update stress change panel (ax8)
    dsigma1_kPa = (sigma1_s - sigma1_0) * stress_scale;
    dsigma1_kPa(gg.nout) = NaN;
    p_dtau.CData = reshape(dsigma1_kPa, gg.nI, gg.nJ);
    ds1_s = stress_scale * reshape(sigma1_s - sigma1_0, gg.nI, gg.nJ); ds1_q = ds1_s(ii8,jj8);
    ds2_s = stress_scale * reshape(sigma2_s - sigma2_0, gg.nI, gg.nJ); ds2_q = ds2_s(ii8,jj8);
    ddx1_s = reshape(cos(t1_s), gg.nI, gg.nJ); ddx1_q = ddx1_s(ii8,jj8);
    ddy1_s = reshape(sin(t1_s), gg.nI, gg.nJ); ddy1_q = ddy1_s(ii8,jj8);
    ddx2_s = reshape(cos(t1_s+pi/2), gg.nI, gg.nJ); ddx2_q = ddx2_s(ii8,jj8);
    ddy2_s = reshape(sin(t1_s+pi/2), gg.nI, gg.nJ); ddy2_q = ddy2_s(ii8,jj8);
    % q_ds1p.UData =  ds1_q.*ddx1_q; q_ds1p.VData =  ds1_q.*ddy1_q;
    % q_ds1n.UData = -ds1_q.*ddx1_q; q_ds1n.VData = -ds1_q.*ddy1_q;
    % q_ds2p.UData =  ds2_q.*ddx2_q; q_ds2p.VData =  ds2_q.*ddy2_q;
    % q_ds2n.UData = -ds2_q.*ddx2_q; q_ds2n.VData = -ds2_q.*ddy2_q;

    % Update right panels
    pblister.CData = (ps.hb)*reshape(vva.hb,gg.nI,gg.nJ); 
    pblister_contour.ZData = (ps.phi)*reshape(vva.pb,gg.nI,gg.nJ);

    phs.CData = ps.hs * reshape(vva.hs,gg.nI,gg.nJ);
    phs_contour.ZData = ps.phi * reshape(vva.phi,gg.nI,gg.nJ);

    peff.CData = (ps.phi/1e6) * reshape(aa.phi_0 - vva.phi, gg.nI, gg.nJ);
    pS.CData = ps.S * reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + ...
               gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr), gg.nI, gg.nJ);

    % update scatter marker size
    time = vva.t;
    markersize_l = lake_input(time,pp,ps);
    norm_size = markersize_l / (Q_peak + eps);
    hscatter.SizeData = 5 + 295 * norm_size.^0.3;

    pvel.CData = vUU;
    vel_field.UData = vux(1:skip:end,1:skip:end);
    vel_field.VData = vuy(1:skip:end,1:skip:end);
    pvel_delta.CData = 100*(vUU - vUU0)./vUU0;

    % Update time markers
    tframe_d = vva.t*ps.t/pd.td;
    for xl = xlines
        xl.Value = tframe_d;
    end
    ttext.String = ['t = ' num2str(tframe_d,'%.1f') ' d'];

    drawnow;
    frame = getframe(f);
    writeVideo(v, frame);
end

try close(v); catch; end
fprintf('Video saved: %s\n', v.Filename);

function out = lake_input(t,pp,ps)
    out = (ps.Q0)*pp.V_l./(sqrt(2*pi)*pp.t_drainage).*...
        exp(-0.5*((t-pp.t_drainage)./pp.t_duration).^2);
end
