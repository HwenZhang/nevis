%% Import necessary libraries
% casename = oo.casename;
casename = 'n2d_100m3s_alpha1e_5_kappa1e_10_mu1e1_V1e7_highres';  % specify the case name
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                         % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code
  
%% colormap
n = 256; % number of color
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% read in the screenshot at the intial timestep
formatSpec = '%04d';
tframe = 1.0
nframe = 1;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
pp.deltap_reg = 5.0e7/ps.phi; % regularization parameter pressure difference (Pa)
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ + 0*vv2.Q);

xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% read in the time series

% steady state outflux
casename = 'n2d_100m3s_alpha1e_5_kappa1e_10_mu1e1_V0e7_test';
tt0 = load(['./results/' casename '/' casename],'tt');
tt0 = tt0.tt; % load the initial time step data
Q_out_Q0 = ps.Q*[tt0.Q_outQ];     % dimensional channel outflux (m^3/s)
Q_out_q0 = ps.Q*[tt0.Q_outq];     % dimensional sheet outflux (m^3/s)
t0 = (ps.t/(24*60*60))*[tt0.t];

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
Q_out_Q0_interp = interp1(t0, Q_out_Q0, t, 'spline', 'extrap');
Q_out_q0_interp = interp1(t0, Q_out_q0, t, 'spline', 'extrap');

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
% R_b = ps.x*[tt.Rb];

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

%% 创建主 figure
f = figure('Position', [100, 100, 600, 600]);
% 创建主 tiledlayout：1 行 2 列（左侧一列、右侧一列）
mainLayout = tiledlayout(f, 5, 1);
mainLayout.TileSpacing = 'compact';
mainLayout.Padding = 'compact';

%% 左侧子 layout：6 行 1 列（垂直）
% panel (a)
ax = nexttile(mainLayout);
plot(ax,t,Q_b_in,'b-',t,Q_b_dec,'r-',LineWidth=1.5);
hold on;
% plot(ax,t,Q_out+Q_out_b,color=[0,0.5,0],LineStyle='-',LineWidth=1.5);

plot(ax,t,Q_out_b,color=[1,0,0],LineStyle='--',LineWidth=1.5);
plot(ax,t,Q_out_Q-Q_out_Q0_interp,color=[0,0.5,0],LineStyle='--',LineWidth=1.5);
plot(ax,t,Q_out_q-Q_out_q0_interp,color=[0,0,1],LineStyle='--',LineWidth=1.5);

plot(ax,t,E,color=[0,0,0],LineStyle='-.',LineWidth=1.5);

x1 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line

xlabel('t [ d ]');
ylabel('Q [ m^3/s ]');
h=legend('Q_{b,in}','Q_{b,relax}','Q_{outb}','\Delta Q_{outQ}','\Delta Q_{outq}','Q_{in}','NumColumns',2);
h.Location='southwest';
text(0.025,0.8,'(a) flux','Units','normalized','FontSize',14)

xlim([tmin_d tmax_d])
set(gca, 'YScale', 'log')
ylim([1e-2 1e4])
% yticks([1e-2,1e1,1e2,1e3,1e4])
grid on
grid minor

% panel (b)
ax = nexttile(mainLayout);
if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])  
    plot(ax,t,pts_N(1,:),'-',LineWidth=1.5);
end
hold on
plot(ax,t,N/1e6,'-',LineWidth=1.5);
x2 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
xlabel('t [ d ]');
ylabel('N [ MPa ]');
h = legend('N at the blister','averaged N');
text(0.025,0.8,'(b) effective pressure','Units','normalized','FontSize',14)
h.Location='southeast';
xlim([tmin_d tmax_d])
% ylim([0 2])
grid on    
grid minor

% panel (d)
ax = nexttile(mainLayout);
yyaxis left
    plot(ax,t,hs_b,'b-',LineWidth=1.5);    
    % plot(t,hs./A,'b-',t,he./A,'b--',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('h at blister [ m ]');
    % ylim([0.02 0.03])
    text(0.025,0.8,'(d) h and S at the lake','Units','normalized','FontSize',14)

yyaxis right
    plot(ax,t,0.25*(Sx_b+Sy_b+Ss_b+Sr_b),'r-',LineWidth=1.5);
    hold on
    x4 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
    ylabel('S at blister [ m^2 ]');
    xlim([tmin_d tmax_d])
    h=legend('h_{cav}','S','NumColumns',2);
    h.Location='southeast';
    % ylim([0 0.05])
    grid on
    grid minor

% panel (e)
ax = nexttile(mainLayout);
% hb_analytical = 3*V_b./(pi*R_b.^2);
yyaxis left
    plot(ax,t,h_b(1,:),'b-',LineWidth=1.5);
    hold on
    % plot(ax,t,hb_analytical,'b--',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('h [ m ]');
    text(0.025,0.8,'(e) h_b and p_b at the lake','Units','normalized','FontSize',14)

yyaxis right
    plot(ax,t,p_b(1,:)/1e6,'r-',LineWidth=1.5); % pressure at the blister
    hold on
    x5 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
    xlabel('t [ d ]');
    ylabel('p_b [ MPa ]');
    xlim([tmin_d tmax_d])
    grid on
    grid minor

% panel (f)
ax = nexttile(mainLayout);
yyaxis left
    plot(ax,t,p_b(1,:)/1e6,'b-',LineWidth=1.5); % blister pressure + \rho_i g h
    hold on
    % plot(ax,t,pb_analytical/1e6,'b--',LineWidth=1.5); % blister pressure + \rho_i g h
    plot(ax,t,p_w/1e6,'b-.',LineWidth=1.5); 

    xlim([tmin_d tmax_d])
    ylim([0 15])
    text(0.025,0.8,'(f) pressure','Units','normalized','FontSize',14)
    xlabel('t [ d ]');
    ylabel('N [ MPa ]');
    
    grid on
    grid minor  
yyaxis right
% Rb_analytical = (3/pi*V_b./h_b).^(1/2);
    plot(ax,t,V_b,'r-',LineWidth=1.5);
    hold on
    % plot(ax,t,Rb_analytical,'r--',LineWidth=1.5);
    % ylim([0 1.1e8])
    ylabel('V_b [ m ]');
    legend('p_b','p_w','V_b','NumColumns',2,location='southeast')

    x6 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line