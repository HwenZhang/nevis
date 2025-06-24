%% Import necessary libraries
% casename = oo.casename;
casename = 'n1d_Vconst_cg0_00_a0_kh0_ks0_mu5e0_c1_V1e7_DX62_5';
load(['./results/' casename '/' casename])
% oo.fn = ['/',casename];                         % filename (same as casename)
% oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
% oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/']; 
addpath(oo.code);                                % add path to code
  
dt = 0.5*pd.td;
tmin = 0*pd.td/ps.t;
tmax = 2000*pd.td/ps.t;
tmin_d = tmin*ps.t/pd.td; 
tmax_d = tmax*ps.t/pd.td;                        % time range for the plot
t_init = 200; t_end = 450;                       % time range for animation
oo.input_constant = 0;

xmax_km = 100; % maximum x value in km
xmin_km = 0;  % minimum x value in km

%% colormap
n = 256; % number of colors
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% read in the screenshot at the intial timestep
formatSpec = '%04d';
tframe = 1.0*pd.td/ps.t; % time frame for the screenshot
nframe = round(tframe/1.0); % frame number
nframe = 1;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge

qs = ps.qs*vv2.qs;
qe = ps.qs*vv2.qe;
qb = ps.qs*vv2.qb;
qQ = ps.qs*vv2.qQ;
Q = ps.qs*vv2.Q;
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qb + vv2.qQ);

xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% read in the time series
t = (ps.t/(24*60*60))*[tt.t];             % dimensional time series (days)
tc = 0.5*(t(1:end-1) + t(2:end));         % center time for plotting
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
% R_b = sgolayfilt(ps.x*[tt.Rb],3,11); % radius of the blister (m)
R_b = smoothdata(ps.x*[tt.Rb],'movmean',101); % radius of the blister (m)
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
f = figure('Position', [100, 100, 1600, 800]);
% 创建主 tiledlayout：1 行 2 列（左侧一列、右侧一列）
mainLayout = tiledlayout(f, 1, 5);
mainLayout.TileSpacing = 'compact';
mainLayout.Padding = 'compact';

leftLayout = tiledlayout(mainLayout, 6, 1);
leftLayout.Layout.Tile = 1;
leftLayout.Layout.TileSpan = [1, 2];

leftLayout.TileSpacing = 'compact';
leftLayout.Padding = 'compact';

rightLayout = tiledlayout(mainLayout, 3, 2);
rightLayout.Layout.Tile = 3;
rightLayout.Layout.TileSpan = [1, 3];

rightLayout.TileSpacing = 'compact';
rightLayout.Padding = 'compact';

%% LEFT layout
% panel (a)
ax = nexttile(leftLayout);
plot(ax,t,Q_b_in,'b-',LineWidth=1.5);
hold on;
plot(ax,t,Q_b_dec,'r-',LineWidth=1.5);
plot(ax,t,Q_out+Q_out_b,color=[0,0.5,0],LineStyle='-',LineWidth=1.5);
plot(ax,t,E,color=[0,0,0],LineStyle='-',LineWidth=1.5);

x1 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line

xlabel('t [ d ]');
ylabel('Q [ m^3/s ]');
h=legend('Q_{b,in}','Q_{b,relax}','Q_{out}','moulin','NumColumns',2);
h.Location='southwest';
text(0.025,0.8,'(a) flux','Units','normalized','FontSize',14)

xlim([tmin_d tmax_d])
set(gca, 'YScale', 'log')
ylim([1e0 1e4])
yticks([1e0,1e1,1e2,1e3,1e4])
grid on
grid minor

% panel (b)
ax = nexttile(leftLayout);
if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])  
    plot(ax,t,pts_N,'-',LineWidth=1.5);
end
hold on
plot(ax,t,N/1e6,'-',LineWidth=1.5);
x2 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
xlabel('t [ d ]');
ylabel('N [ MPa ]');
h = legend('N at the blister','averaged N');
text(0.025,0.8,'(b) effective pressure','Units','normalized','FontSize',14)
h.Location='southwest';
xlim([tmin_d tmax_d])
% ylim([0 2])
grid on    
grid minor

% panel (c)
ax = nexttile(leftLayout);
yyaxis left
    plot(ax,t,hs./A,'b-',t,he./A,'k-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average h [ m ]');
    ylim([0 0.3])
    text(0.025,0.8,'(c) averaged h and S','Units','normalized','FontSize',14)
yyaxis right
    plot(ax,t,S/A,'r-',LineWidth=1.5);
    hold on
    x3 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
    ylabel('Channel h [ m ]');
    h=legend('h_{cav}','h_{e}','S','NumColumns',2);
    h.Location='southwest';
    xlim([tmin_d tmax_d])
    % ylim([0 0.05])
    grid on
    grid minor

% panel (d)
ax = nexttile(leftLayout);
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
    h.Location='southwest';
    % ylim([0 0.05])
    grid on
    grid minor

% panel (e)
ax = nexttile(leftLayout);
% hb_analytical = 3*V_b./(pi*R_b.^2);
yyaxis left
    plot(ax,t,h_b,'b-',LineWidth=1.5);
    hold on
    % plot(ax,t,hb_analytical,'b--',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('h [ m ]');
    text(0.025,0.8,'(e) h_b and p_b at the lake','Units','normalized','FontSize',14)

yyaxis right
    plot(ax,t,R_b,'r-',LineWidth=1.5);
    hold on
    x5 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
    xlabel('t [ d ]');
    ylabel('R_b [ m ]');
    xlim([tmin_d tmax_d])
    grid on
    grid minor

% panel (f)
ax = nexttile(leftLayout); 
yyaxis left
    plot(ax,t,p_b/1e6,'b-',LineWidth=1.5); % blister pressure + \rho_i g h
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
    dRdt = gradient(R_b,t);
    plot(ax,t,dRdt,'r-',LineWidth=1.5);
    hold on
    % plot(ax,t,Rb_analytical,'r--',LineWidth=1.5);
    ylim([0 100])
    ylabel('dV_b/dt [ m/d ]');
    legend('p_b','p_w','dV_b/dt','NumColumns',2,location='northeast')

    x6 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line

%% RIGHT layout：3 行 2 列
% panel(a): total flux
ax = nexttile(rightLayout);
% yyaxis left
pqs = plot(ax,xx(gg.ns),qs(gg.ns),'b-','LineWidth',1.5);
hold on
pqb = plot(ax,xx(gg.ns),qb(gg.ns),'b--','LineWidth',1.5);
pqQ = plot(ax,xx(gg.ns),qQ(gg.ns),'b-.','LineWidth',1.5);
pflux = plot(ax, xx(gg.ns),qnet(gg.ns),'k-','LineWidth',1.5);

title('flux','FontSize',14);
ylim([0 2*ps.qs])
ylabel('q [m^3/s]')
xlim([xmin_km xmax_km])
xlabel('x (km)')
legend([pqs pqb pqQ pflux],'q_s','q_b','q_Q','q_{net}','NumColumns',2,'Location','northwest');
grid on
grid minor
if vva.t*ps.t/pd.td <= 700.0
    time = [num2str(700.0-vva.t*ps.t/pd.td,'%.1f') ' d' ' before drainage'];
else
    time = [num2str(vva.t*ps.t/pd.td-700.0,'%.1f') ' d' ' after drainage'];
end
ttext = text(0.1,8,time,"Position",[0.5 0.9],"Units","normalized");
ttext.FontSize=18;
text(0.025,0.1,'(1) flux','Units','normalized','FontSize',14)

% yyaxis right
%     pphi1 = plot(ax,xx,ps.phi/1e6*vva.phi,'r-','LineWidth',1.5);
%     ylim([0 10]);
%     xlim([xmin_km xmax_km])
%     xlabel('x (km)')
%     ylabel('\phi [ MPa ]')
%     grid on
%     grid minor

% panel(f): effective pressure
ax = nexttile(rightLayout);
L = 5e4; % length of the domain [m]
b = 4e3 - 0.01*xx*1e3;   % bed slope angle = 0.01
s = (1060)*(1-1e3*xx/L).^0.5; % ice surface topography 
s = b+1000*ones(size(s)); % constant ice surface

H = max(s-b,0); % ice thickness
plot(ax,xx(gg.nin),s(gg.nin),'b-','LineWidth',1.5);
hold on
plot(ax,xx(gg.nin),b(gg.nin),'k-','LineWidth',1.5);
xlabel('x (km)')
ylabel('z [m]')
title('ice surface and bed','FontSize',14);
% ylim([0 1100])
xlim([xmin_km xmax_km])
grid on
grid minor
text(0.025,0.1,'(2) ice surface and bed','Units','normalized','FontSize',14)

% panel(b): cavity sheet thickness
ax = nexttile(rightLayout);
yyaxis left
    phs = plot(ax,xx(gg.ns),ps.hs*vva.hs(gg.ns),'b-','LineWidth',1.5); 
    title('cavity sheet and \phi contour','FontSize',14);
    ylim([0 0.1]);
    ylabel('h_s [m]')
yyaxis right
    pphi2 = plot(ax,xx(gg.nin),ps.phi/1e6*vva.phi(gg.nin),'r-','LineWidth',1.5);
    hold on
    xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
    xlabel('x (km)')
    ylabel('\phi [ MPa ]')
    ylim([0 10]);
    xlim([xmin_km xmax_km])
    grid on
    grid minor
text(0.025,0.1,'(3) cavity sheet','Units','normalized','FontSize',14)

% panel(c): blister sheet thickness
ax = nexttile(rightLayout);
yyaxis left
    pblister = plot(ax,xx(gg.ns),ps.hb*vva.hb(gg.ns),'b-','LineWidth',1.5); 
    ylim([-0.1 0.1]);
    ylabel('h_b [m]')
    title('blister sheet and pb','FontSize',14);
yyaxis right
    pphi3 = plot(xx(gg.nin),ps.phi/1e6*(vva.pb(gg.nin)+aa.phi_a(gg.nin)),'r-','LineWidth',1.5);
    hold on
    xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
    xlabel('x (km)')
    ylabel('\phi [ MPa ]')
    % ylim([0 10]);
    xlim([xmin_km xmax_km])
    grid on
    grid minor
text(0.725,0.5,'(4) blister sheet','Units','normalized','FontSize',14)

% panel(d): cross-sectional area
ax = nexttile(rightLayout);
Smean = (0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr)); 
% yyaxis left
pS = plot(ax,xx(gg.nin),ps.S*Smean(gg.nin),'b-','LineWidth',1.5); 
ylim([0 8])
title('channel cross section','FontSize',14);
ylabel('y [km]')

% yyaxis right
%     pphi4 = plot(ax,xx(gg.nin),ps.phi/1e6*vva.pb(gg.nin),'r-','LineWidth',1.5);
%     hold on
%     xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
%     ylabel('\phi [ MPa ]')
%     ylim([0 10]);
xlabel('x (km)')
xlim([xmin_km xmax_km])
grid on
grid minor
text(0.025,0.1,'(5) channel cross section','Units','normalized','FontSize',14)

% panel(e): effective pressure
ax = nexttile(rightLayout);
N = ((aa.phi_0-vva.phi)); 
peff = plot(ax,xx(gg.nin),(ps.phi/10^6)*N(gg.nin),'k-','LineWidth',1.5);  
yline(0,'-k','LineWidth',1.5); % horizontal line at y=0
text(0.03,0.55, 'N=0','Units','normalized','FontSize',14,'Color','k');
% id_max_hb = find(vva.hb(gg.nin)==max(vva.hb(gg.nin)));
% xline(xx(id_max_hb),'--k','LineWidth',1.5); % dashed line

title('effective pressure','FontSize',14);
xlabel('x (km)')
ylabel('N [MPa]')
hold on
% add moulins
% if ~isfield(pp,'ni_m'), pp.ni_m = []; end
% x = (ps.x/10^3)*gg.nx(pp.ni_m);
% mscale = 100;
% for i_m = 1:length(pp.ni_m)
%     if aa.E(pp.ni_m(i_m))>0
%         plot(x(i_m),0,'ko','Markersize',8+aa.E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%     else
%         plot(x(i_m),0,'ko','Markersize',8,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%     end
% end
axis tight
ylim([-3 3]);
grid on
grid minor
text(0.025,0.1,'(6) effective pressure','Units','normalized','FontSize',14)

%% make video
v = VideoWriter(['./results/videos/' casename],'MPEG-4');
% v = VideoWriter(['./results/' oo.casename],'MPEG-4');
v.FrameRate = 1;
open(v)
for i_t = t_init:t_end
    disp(['Frame ',num2str(i_t-t_init),' / ',num2str(t_end-t_init),' ...']);
    %% load timestep
    vva = load([path num2str(i_t,formatSpec)], 'vv');
    vva = vva.vv;
    aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);

    N = (ps.phi/10^6)*(aa.phi_0-vva.phi);
    pb = (ps.phi/10^6)*(vva.pb-(aa.phi_0-aa.phi_a));
    [vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge

    qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qb + vv2.qQ);
    qs = ps.qs*vv2.qs;
    qb = ps.qs*vv2.qb;
    qQ = ps.qs*vv2.qQ;

    pqs.YData = qs(gg.ns);
    pqb.YData = qb(gg.ns);
    pqQ.YData = qQ(gg.ns);
    pflux.YData = qnet(gg.ns);

    phs.YData = (ps.hs)*vva.hs(gg.ns);
    pblister.YData = (ps.hb)*vva.hb(gg.ns); 
    Smean = (ps.S)*0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr);
    pS.YData = Smean(gg.nin);
    peff.YData = N(gg.nin);
    % ppb.YData = pb(gg.nin);

    % pphi1.YData = (ps.phi/1e6)*vva.phi(gg.nin);
    pphi2.YData = (ps.phi/1e6)*vva.phi(gg.nin);
    pphi3.YData = (ps.phi/1e6)*(vva.pb(gg.nin)+aa.phi_a(gg.nin));
    % pphi4.YData = (ps.phi/1e6)*vva.pb(gg.nin);

    if vva.t*ps.t/pd.td <= 700.0
        time = [num2str(700.0-vva.t*ps.t/pd.td,'%.1f') ' d' ' before drainage'];
    else
        time = [num2str(vva.t*ps.t/pd.td-700.0,'%.1f') ' d' ' after drainage'];
    end
    set(ttext,{'string'},{time})  %notice the column vector of new values

    x1.Value = vva.t*ps.t/pd.td;
    x2.Value = vva.t*ps.t/pd.td;
    x3.Value = vva.t*ps.t/pd.td;
    x4.Value = vva.t*ps.t/pd.td;
    x5.Value = vva.t*ps.t/pd.td;
    x6.Value = vva.t*ps.t/pd.td;
    
    disp(min(vva.hb));
    refreshdata
    drawnow
    % pause(0.2)
    frame = getframe(gcf);
    writeVideo(v,frame)
end
close(v)