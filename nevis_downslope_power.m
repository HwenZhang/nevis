%% Import necessary libraries
casename = 'n1d_Vconst_cg0_00_a0_kh0_ks0_mu5e0_c1_V1e7';
load(['./results/' casename '/' casename])
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code
  
dt = oo.dt;
tmin = 0*pd.td/ps.t;
tmax = 400*pd.td/ps.t;
tmin_d = tmin*ps.t/pd.td; 
tmax_d = tmax*ps.t/pd.td;                        % time range for the plot

%% colormap
n = 256; % number of colors
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% read in the screenshot at the intial timestep
formatSpec = '%04d';
tframe = 300.0*pd.td/ps.t; % time frame for the screenshot
nframe = round(tframe/1.0); % frame number
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
pp.hb_reg1 = 0;
pp.hb_reg2 = 0;
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ);
qs = ps.qs*vv2.qs;
qb = ps.qs*vv2.qb;
qQ = ps.qs*vv2.qQ;
Q = ps.qs*vv2.Q;

xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% read in the time series
t = (ps.t/(24*60*60))*[tt.t];             % dimensional time series (days)
Q_b_in = pd.Q_0*[tt.Qb_in];               % dimensional influx (m^3/s)
Q_b_dec = ps.h*ps.x^2/ps.t*[tt.Qb_dec];   % dimensional relaxation (m^3/s)

Q_in = ps.Q*[tt.Q_in];          % dimensional influx (m^3/s)
Q_out = ps.Q*[tt.Q_out];        % dimensional outflux (m^3/s)
Qb_out = ps.Q0*[tt.Qb_out];     % dimensional outflux from blister sheet mass conservation (m^3/s)

Q_out_b = ps.Q0*[tt.Q_outb];    % dimensional sheet outflux (m^3/s)
Q_out_Q = ps.Q*[tt.Q_outQ];     % dimensional channel outflux (m^3/s)
Q_out_q = ps.Q*[tt.Q_outq];     % dimensional blister outflux (m^3/s)

m = (ps.m*ps.x^2)*[tt.m];      % dimensional melting rate (m^3/s)
E = (ps.m*ps.x^2)*[tt.E];      % dimensional source terms  (m^3/s)
h_b = ps.hb*[tt.pts_hb];       %
p_b = ps.phi*[tt.pts_pb];      %
V_b = ps.x^2*ps.hb*[tt.Vb];
R_b = ps.x*[tt.Rb];
R_b = smoothdata(ps.x*[tt.Rb],'gaussian',101); % radius of the blister (m)
hb_max = ps.hb*[tt.hb_max]; % maximum blister base thickness (m)

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

S = ps.x*ps.S*[tt.S];
A = ps.x^2*sum(gg.Dx.*gg.Dy);

if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])    
    pts_phi = (ps.phi/10^6)*[tt.pts_phi];
    pts_hs = ps.hs*[tt.pts_hs];
    pts_N = (ps.phi/10^6)*(aa.phi_0(oo.pts_ni)*[tt.t].^0 - [tt.pts_phi]);
    pts_pw = (ps.phi/10^6)*([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_prat = ([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0)./...
               (aa.phi_0(oo.pts_ni)*[tt.t].^0-aa.phi_a(oo.pts_ni)*[tt.t].^0);
end

%% log-log plot of the power law
inject_time = 51;
growth_phase = t > inject_time; % Define the growth phase for the power law fit
t0 = t(growth_phase); % Filter out times before 200 days
R_b0 = R_b(growth_phase); % Filter out corresponding R_b values
hb_max0 = hb_max(growth_phase); % Filter out corresponding hb_max values
V2d = V_b(growth_phase)/1e5; % Filter out corresponding V_b values
xf = (0.75)^(2/3)*(pd.rho_w*pd.g*0.01/pd.mu)^(1/3)*(V_b/1e5).^(2/3).*((t-inject_time)*pd.td).^(1/3);
Q0 = (ps.Q0*ps.t).*pp.V_l./(pp.t_duration*ps.t)/1e5;

% Rb_a = (8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((t0-inject_time)*pd.td).^(0.4);

figure("Position",[100,100,800,300]);
subplot(1,1,1)
loglog(t0-inject_time, R_b0, 'k-', 'LineWidth', 2);
hold on
loglog(t-inject_time, xf, 'b--', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Blister radius(front) (m)');
title('Blister radius(front) vs time');
grid on;
set(gca, 'FontSize', 14);

% yyaxis right
%         plot(t0-inject_time, V2d, 'g-', 'LineWidth', 2);
%         hold on
%         plot(t0-inject_time, Q0*((t0-inject_time)*pd.td), 'm--', 'LineWidth', 2);
%         ylabel('Blister volume (km^3)');
% legend('Blister radius(front)', 'Analytical solution x_f=(0.75V_0)^{2/3}(\frac{\rho_w g \theta}{\mu})t^{1/3}', 'Location', 'NorthWest','Interpreter','latex');
legend( ...
  'Blister radius (front)', ...
  '$x_f=(0.75V_0)^{2/3}\!\bigl(\frac{\rho_w g \theta}{\mu}\bigr)t^{1/3}$', ...
  'Location','NorthWest', ...
  'Interpreter','latex' ...
);


% Save the figure
saveas(gcf, ['./blister_power_law.png']);