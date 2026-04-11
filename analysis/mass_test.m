%% Animation script for NEVIS 2D ice-hydrology coupled model
% Produces a video with time series (left) and spatial fields (right):
%   Left:  Q fluxes, effective pressure, sheet/channel, ice speed
%   Right: net flux, cavity thickness, effective pressure, 
%          channel cross-section, ice speed, ice speed change
clc; clear; close all;

%% Settings
casename = 'n2d_regional_eps1e_02_kappa1e_10_mu2e1_partition5e_01_spinup';
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
tmax_yr = tmin_yr + 1;  % end time in years
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


%% ===== Left layout: time series =====
tframe_d = vva.t*ps.t/pd.td;
tc_s = 0.5*(t(1:end-1)+t(2:end))*pd.td;
dvdt = gradient(hs_ts+S_ts+V_b, t*pd.td);

% (a) Fluxes
figure();
hold on;
plot(t, Q_out_Q + Q_out_q,'-', 'Color',[0,0.25,0], 'LineWidth',2.5);
plot(t, Q_out_b2,'b--', 'LineWidth',1.5);
plot(t, E,'k-.', 'LineWidth',1.5);
plot(t, -dvdt);


xlabel('t [d]'); ylabel('Q [m^3/s]');
legend('Q_{b,relax}','Q_{out}','Q_{out,b}','Q_{in}','NumColumns',2,'Location','southeast');
text(0.025,0.85,'(a) flux','Units','normalized','FontSize',12);
xlim([tmin tmax]); set(gca,'YScale','log'); ylim([1e-2 1e4]);
grid on; grid minor;

plot(t, hs_ts)