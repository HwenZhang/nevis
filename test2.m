%% Import data
casename = 'n2d_100m3s_alpha1e_5_kappa1e_10_mu1e1_V1e7_highres';  % specify the case name
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                         % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code

t = (ps.t/(24*60*60))*[tt.t];               % dimensional time series (days)
tspan = (ps.t/pd.td)*oo.t_span;
tmin = 2.9*365*pd.td/ps.t;
tmax = 3.1*365*pd.td/ps.t;
tmin_d = tmin*ps.t/pd.td; 
tmax_d = tmax*ps.t/pd.td;                   % time range for the plot
[~,t_init] = min(abs(tspan-365*2.9));             % initial time step
[~,t_end] = min(abs(tspan-3.1*365));              % final time step
% t_end = 1200;

% extrpolate Q_out_Q0
% Q_out_Q0_interp = interp1(t0, Q_out_Q0, t, 'spline', 'extrap');
% Q_out_q0_interp = interp1(t0, Q_out_q0, t, 'spline', 'extrap');

Q_b_in = pd.Q_0*[tt.Qb_in];               % dimensional influx (m^3/s)
Q_b_dec = ps.h*ps.x^2/ps.t*[tt.Qb_dec];   % dimensional relaxation (m^3/s)

Q_in = ps.Q*[tt.Q_in];          % dimensional influx (m^3/s)
Q_out = ps.Q*[tt.Q_out];        % dimensional outflux (m^3/s)
Qb_out = ps.Q0*[tt.Qb_out];     % dimensional outflux from blister sheet mass conservation (m^3/s)

Q_out_b = ps.Q0*[tt.Q_outb];    % dimensional blister outflux (m^3/s)
Q_out_Q = ps.Q*[tt.Q_outQ];     % dimensional channel outflux (m^3/s)
Q_out_q = ps.Q*[tt.Q_outq];     % dimensional sheet outflux (m^3/s)

m = (ps.m*ps.x^2)*[tt.m];      % dimensional melting rate (m^3/s)
E = (ps.m*ps.x^2)*[tt.E];      % dimensional source terms  (m^3/s)
h_b = ps.hb*[tt.pts_hb];       %
p_b = ps.phi*[tt.pts_pb];      %
V_b = ps.x^2*ps.hb*[tt.Vb];
V_c = ps.x*ps.S*[tt.S];
V_s = ps.x^2*ps.h*[tt.hs];

% phi = (ps.phi)*[tt.phi];     % dimensional hydrulic potential (MPa)
N = (ps.phi)*[tt.N];           % dimensional effective stress (MPa)
hs = ps.x^2*ps.h*[tt.hs];      % integrated hs (m^3)
hs_b = ps.h*[tt.hs_b];         % integrated hs (m^3)
he = ps.x^2*ps.h*[tt.he];      % integrated he (m^3)
p_w = ps.phi*[tt.pwb];     % dimensional hydraulic potential at the lake (MPa)

Sx_b = ps.S*[tt.Sx_b];
Sy_b = ps.S*[tt.Sy_b];
Ss_b = ps.S*[tt.Ss_b];
Sr_b = ps.S*[tt.Sr_b];
S_ave = (0.25*(Sx_b.^2 + Sy_b.^2 + Ss_b.^2 + Sr_b.^2)).^0.5;


A = ps.x^2*sum(gg.Dx.*gg.Dy);

if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])    
    pts_phi = (ps.phi/10^6)*[tt.pts_phi];
    pts_hs = ps.hs*[tt.pts_hs];
    pts_N = (ps.phi/10^6)*(aa.phi_0(oo.pts_ni)*[tt.t].^0 - [tt.pts_phi]);
    pts_pw = (ps.phi/10^6)*([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_prat = ([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0)./...
               (aa.phi_0(oo.pts_ni)*[tt.t].^0-aa.phi_a(oo.pts_ni)*[tt.t].^0);
end

%% 创建主 figure
f = figure('Position', [100, 100, 700, 600]);   % wider to accommodate two columns
% 改成 3 行 2 列 布局
mainLayout = tiledlayout(f, 3, 2);
mainLayout.TileSpacing = 'compact';
mainLayout.Padding = 'compact';



%% panel (a) — tile 1
ax = nexttile(mainLayout, 1);
plot(ax, t, Q_b_in, 'b-', t, Q_b_dec, 'r-', 'LineWidth', 1.5); hold on;
plot(ax, t, Q_out_b, '--r', 'LineWidth', 1.5);
plot(ax, t, Q_out_Q, '--', 'Color', [0,0.5,0], 'LineWidth', 1.5);
plot(ax, t, Q_out_q, '--', 'Color', [0,0,1], 'LineWidth', 1.5);
plot(ax, t, E, '-.', 'LineWidth', 1.5);
tframe = 0;
xline(ax, tframe * ps.t/pd.td, '--k', 'LineWidth', 1.5);
xlim(ax, [tmin_d tmax_d]);
set(ax, 'YScale', 'log');
ylim(ax, [1e-2 1e4]);
xlabel(ax, 't [ d ]');
ylabel(ax, 'Q [ m^3/s ]');
text(ax, 0, 1.08, '(a) flux', 'Units', 'normalized', 'FontSize', 14, 'Clipping', 'off');

% legend below
h = legend(ax, 'Q_{b,in}', 'Q_{b,relax}', 'Q_{outb}', '\DeltaQ_{outQ}', '\DeltaQ_{outq}', 'Q_{in}');
h.Orientation = 'horizontal';
h.Location = 'southoutside';
h.Box = 'off';

grid(ax, 'on');
grid(ax, 'minor');

%% panel (b) — tile 2
ax = nexttile(mainLayout, 2);
if exist('pts_N','var')
    plot(ax, t, pts_N(1,:), '-', 'LineWidth', 1.5); hold on;
end
plot(ax, t, N/1e6, '-', 'LineWidth', 1.5);
xline(ax, tframe * ps.t/pd.td, '--k', 'LineWidth', 1.5);

xlabel(ax, 't [ d ]');
ylabel(ax, 'N [ MPa ]');
text(ax, 0, 1.08, '(b) effective pressure', 'Units', 'normalized', 'FontSize', 14, 'Clipping', 'off');

h = legend(ax, 'N at the blister', 'averaged N');
h.Orientation = 'horizontal';
h.Location = 'southoutside';
h.Box = 'off';

xlim(ax, [tmin_d tmax_d]);
grid(ax, 'on');
grid(ax, 'minor');

% %% panel (d) — tile 3
% ax = nexttile(mainLayout, 3);
% yyaxis(ax, 'left');
% plot(ax, t, hs_b, 'b-', 'LineWidth', 1.5);
% ylabel(ax, 'h at blister [ m ]');

% yyaxis(ax, 'right');
% plot(ax, t, 0.25*(Sx_b+Sy_b+Ss_b+Sr_b), 'r-', 'LineWidth', 1.5);
% ylabel(ax, 'S at blister [ m^2 ]');

% xline(ax, tframe * ps.t/pd.td, '--k', 'LineWidth', 1.5);
% xlim(ax, [tmin_d tmax_d]);
% xlabel(ax, 't [ d ]');
% text(ax, 0.025, 0.8, '(d) h and S at the lake', 'Units', 'normalized', 'FontSize', 14);

% h = legend(ax, 'h_{cav}', 'S');
% h.Orientation = 'horizontal';
% h.Location = 'southoutside';
% h.Box = 'off';

% grid(ax, 'on');
% grid(ax, 'minor');

%% panel (e) — tile 4
ax = nexttile(mainLayout, 3);
yyaxis(ax, 'left');
plot(ax, t, h_b(1,:), 'b-', 'LineWidth', 1.5);
ylabel(ax, 'h [ m ]');

yyaxis(ax, 'right');
plot(ax, t, p_b(1,:)/1e6, 'r-', 'LineWidth', 1.5);
hold on;
plot(ax, t, p_w/1e6, '-.', 'LineWidth', 1.5);
ylabel(ax, 'p [ MPa ]');
ylim([4 10]);

xline(ax, tframe * ps.t/pd.td, '--k', 'LineWidth', 1.5);
xlim(ax, [tmin_d tmax_d]);
xlabel(ax, 't [ d ]');
text(ax, 0, 1.08, '(e) h_b and p_b at the lake', 'Units', 'normalized', 'FontSize', 14, 'Clipping', 'off');

h = legend(ax, 'h_b', 'p_b', 'p_w');
h.Orientation = 'horizontal';
h.Location = 'southoutside';
h.Box = 'off';

grid(ax, 'on');
grid(ax, 'minor');

%% panel (f) — tile 5
ax = nexttile(mainLayout, 4);
% yyaxis(ax, 'left');
plot(ax, t, V_b, 'r-', 'LineWidth', 1.5);
hold on;
plot(ax, t, V_c, 'b-', 'LineWidth', 1.5);
plot(ax, t, V_s, '-', 'Color', [0, 0.5, 0], 'LineWidth', 1.5);

xline(ax, tframe * ps.t/pd.td, '--k', 'LineWidth', 1.5);
xlim(ax, [tmin_d tmax_d]);
ylim(ax, [1e5 1e8]);
xlabel(ax, 't [ d ]');
ylabel(ax, 'V [ m ]');
text(ax, 0, 1.08, '(f) volume', 'Units', 'normalized', 'FontSize', 14, 'Clipping', 'off');

h = legend(ax, 'V_b', 'V_c', 'V_s');
h.Orientation = 'horizontal';
h.Location = 'southoutside';
h.Box = 'off';

set(ax, 'YScale', 'log');
grid(ax, 'on');
grid(ax, 'minor');

%%
ax = nexttile(mainLayout, 5, [1 2]);  % spans both columns on bottom row