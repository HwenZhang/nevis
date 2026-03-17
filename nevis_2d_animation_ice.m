%% Animation script for NEVIS 2D ice-hydrology coupled model
% Produces a video with time series (left) and spatial fields (right):
%   Left:  Q fluxes, effective pressure, sheet/channel, ice speed
%   Right: net flux, cavity thickness, effective pressure, 
%          channel cross-section, ice speed, ice speed change
clc; clear; close all;

%% Settings
casename = 'n2d_region_ice_Cinv_test2_drainage';
% casename = 'n2d_region_ice_Cinv_test2';
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];
oo.rn = [oo.root,oo.results,oo.fn];
oo.code = '../nevis/src';
filepath = [oo.rn,'/'];
addpath(oo.code);
formatSpec = '%04d';

%% Colormap
n = 256;
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% Time range
tmin_yr = 0.0;           % start time in years
tmax_yr = tmin_yr + 0.45;  % end time in years
tmin = tmin_yr * 365;  % in days
tmax = tmax_yr * 365;

tspan_d = (ps.t/pd.td) * oo.t_span;  % t_span in days
index = find(tspan_d >= tmin & tspan_d <= tmax);
if isempty(index)
    error('No frames found in the specified time range [%.1f, %.1f] days', tmin, tmax);
end
fprintf('Found %d frames in time range [%.1f, %.1f] days\n', length(index), tmin, tmax);

%% Read time series data
t = (ps.t/pd.td) * [tt.t];  % dimensional time (days)

Q_b_in  = pd.Q_0 * [tt.Qb_in];
Q_b_dec = ps.h*ps.x^2/ps.t * [tt.Qb_dec];
Q_in    = ps.Q * [tt.Q_in];
Q_out   = ps.Q * [tt.Q_out];
Q_out_Q = ps.Q * [tt.Q_outQ];
Q_out_q = ps.Q * [tt.Q_outq];
Q_out_b = ps.Q * [tt.Q_outb];
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
nframe = index(1);
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

% Left panel: time series (4 rows)
leftLayout = tiledlayout(mainLayout, 4, 1);
leftLayout.Layout.Tile = 1;
leftLayout.Layout.TileSpan = [1, 2];
leftLayout.TileSpacing = 'compact';
leftLayout.Padding = 'compact';

% Right panel: spatial fields (3 rows x 2 cols)
rightLayout = tiledlayout(mainLayout, 3, 2);
rightLayout.Layout.Tile = 3;
rightLayout.Layout.TileSpan = [1, 3];
rightLayout.TileSpacing = 'compact';
rightLayout.Padding = 'compact';

%% ===== Left layout: time series =====
tframe_d = vva.t*ps.t/pd.td;

% (a) Fluxes
ax_a = nexttile(leftLayout);
plot(ax_a, t, Q_b_dec,'r-', 'LineWidth',1.5);
% plot(ax_a, t, Q_b_in,'b-', t, Q_b_dec,'r-', 'LineWidth',1.5);
hold on;
plot(ax_a, t, Q_out_Q + Q_out_q + Q_out_b,'-', 'Color',[0,0.25,0], 'LineWidth',2.5);
% plot(ax_a, t, Q_out_q,'b--', 'LineWidth',1.5);
plot(ax_a, t, E,'k-.', 'LineWidth',1.5);
x1 = xline(tframe_d,'--k','LineWidth',1.5);
xlabel('t [d]'); ylabel('Q [m^3/s]');
legend('Q_{b,relax}','Q_{out}','Q_{in}','NumColumns',2,'Location','southeast');
text(0.025,0.85,'(a) flux','Units','normalized','FontSize',12);
xlim([tmin tmax]); set(gca,'YScale','log'); ylim([1e-2 1e4]);
grid on; grid minor;

% (b) Effective pressure
ax_b = nexttile(leftLayout);
plot(ax_b, t, N_ts/1e6,'-', 'LineWidth',1.5);
hold on;
x2 = xline(tframe_d,'--k','LineWidth',1.5);
xlabel('t [d]'); ylabel('N [MPa]');
legend('averaged N','Location','southeast');
text(0.025,0.85,'(b) effective pressure','Units','normalized','FontSize',12);
xlim([tmin tmax]); 
ylim([0 1.1]);
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
legend('h_{cav}','h_e','S','NumColumns',3,'Location','southeast');
text(0.025,0.85,'(c) averaged h and S','Units','normalized','FontSize',12);
xlim([tmin tmax]); grid on; grid minor;

% (d) Ice speed time series (if available)
ax_d = nexttile(leftLayout);
yyaxis left
plot(ax_d, t, U*pd.ty, 'b-', 'LineWidth',1.5);
x4 = xline(tframe_d,'--k','LineWidth',1.5);
xlabel('t [d]'); ylabel('U [m/yr]');
text(0.025,0.85,'(d) blister volume and ice speed','Units','normalized','FontSize',12);
xlim([tmin tmax]); 
ylim([0 200]); 
grid on; grid minor;

yyaxis right
plot(ax_d, t, V_b, 'r-', 'LineWidth',1.5);
ylabel('Vb [m^3]');

xlines = [x1, x2, x3, x4];

%% ===== Right layout: spatial fields =====

% (1) Blister sheet thickness hb
ax1 = nexttile(rightLayout);
zhb = ps.hb * reshape(vva.hb,gg.nI,gg.nJ);
zphi = ps.phi * reshape(vva.phi,gg.nI,gg.nJ);
phb = pcolor(ax1, xx, yy, zhb);
set(phb,'linestyle','none');
cx = colorbar(); colormap(ax1, cmap);
clim([-1 1]);
cx.Label.String = 'h_b [m]';
hold on;
[~, phb_contour] = contour(ax1, xx, yy, zphi, 'linecolor','k','linewidth',0.5);
title('(e) blister sheet h_b'); ylabel('y (km)');
axis equal; axis tight;

% (2) Cavity sheet thickness
ax2 = nexttile(rightLayout);
vva.hs(gg.nout) = NaN;
zhs = ps.hs * reshape(vva.hs,gg.nI,gg.nJ);
phs = pcolor(ax2, xx, yy, zhs);
set(phs,'linestyle','none');
cx = colorbar();colormap(ax2, cmap);
cx.Label.String = 'h_s [m]';
clim([0 0.1]);
hold on;
[~, phs_contour] = contour(ax2, xx, yy, zphi, 'linecolor','k','linewidth',0.5);
title('(f) cavity sheet'); ylabel('y (km)');
axis equal; axis tight;

% (3) Effective pressure
ax3 = nexttile(rightLayout);
zN = (ps.phi/1e6) * reshape(aa.phi_0 - vva.phi, gg.nI, gg.nJ);
peff = pcolor(ax3, xx, yy, zN);
set(peff,'linestyle','none');
cx = colorbar(); colormap(ax3, cmap);
cx.Label.String = 'N [MPa]';
clim([-3 3]);
title('(g) effective pressure'); ylabel('y (km)');
axis equal; axis tight;

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
title('(h) channel cross section'); ylabel('y (km)');
axis equal; axis tight;

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
title('(i) ice speed'); ylabel('y (km)');
axis equal; axis tight;

% (6) Ice speed change
ax6 = nexttile(rightLayout);
pvel_delta = pcolor(ax6, xx, yy, zeros(gg.nI,gg.nJ));
set(pvel_delta,'linestyle','none');
cx = colorbar(); colormap(ax6, cmap);
cx.Label.String = '\Delta U [m/yr]';
clim([-100 100]);
ttext = text(ax6, 0.5, 0.9, ['t = ' num2str(tframe_d,'%.1f') ' d'], ...
    'Units','normalized','FontSize',14,'FontWeight','bold');
title('(j) ice speed change'); ylabel('y (km)'); xlabel('x (km)');
axis equal; axis tight;

%% ===== Animation loop =====
v = VideoWriter(['./results/videos/' casename],'MPEG-4');
v.FrameRate = 3;
open(v);

for i_idx = 1:length(index)
    if ~isvalid(f)
        disp('Figure closed. Animation stopped.');
        break
    end
    i_t = index(i_idx);
    fprintf('Frame %d / %d ...\n', i_idx, length(index));

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

    % Update right panels
    phb.CData = ps.hb * reshape(vva.hb,gg.nI,gg.nJ);
    phb_contour.ZData = ps.phi * reshape(vva.phi,gg.nI,gg.nJ);

    phs.CData = ps.hs * reshape(vva.hs,gg.nI,gg.nJ);
    phs_contour.ZData = ps.phi * reshape(vva.phi,gg.nI,gg.nJ);

    peff.CData = (ps.phi/1e6) * reshape(aa.phi_0 - vva.phi, gg.nI, gg.nJ);
    pS.CData = ps.S * reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + ...
               gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr), gg.nI, gg.nJ);

    pvel.CData = vUU;
    vel_field.UData = vux(1:skip:end,1:skip:end);
    vel_field.VData = vuy(1:skip:end,1:skip:end);
    pvel_delta.CData = vUU - vUU0;

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