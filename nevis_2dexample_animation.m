%% Import necessary libraries
% casename = oo.casename;
casename = 'n2d_30mm_cg0_00_kappa1e_5_kh0_ks0_mu1e6_c1_V1e7_test'; 

load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                         % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code
  
dt = 0.5*pd.td/ps.t;

% used for runs with initial condition
tmin = 0*pd.td/ps.t;
tmax = 2*365*pd.td/ps.t;
tmin_d = tmin*ps.t/pd.td; 
tmax_d = tmax*ps.t/pd.td;                       % time range for the plot
t_init = 365-40; t_end = 365+80;                      % time range for animation

%% colormap
n = 256; % number of colors
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% read in the screenshot at the intial timestep
formatSpec = '%04d';
tframe = 1.0;
nframe = 1;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ + 0*vv2.Q);

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

%% 左侧子 layout：6 行 1 列（垂直）
% panel (a)
ax = nexttile(leftLayout);
plot(ax,t,Q_b_in,'b-',t,Q_b_dec,'r-',LineWidth=1.5);
hold on;
plot(ax,t,Q_out+Q_out_b,color=[0,0.5,0],LineStyle='-',LineWidth=1.5);

plot(ax,t,Q_out_b,color=[1,0,0],LineStyle='--',LineWidth=1.5);
plot(ax,t,Q_out_Q,color=[0,1,0],LineStyle='--',LineWidth=1.5);
plot(ax,t,Q_out_q,color=[0,0,1],LineStyle='--',LineWidth=1.5);

plot(ax,t,E,color=[0,0,0],LineStyle='-',LineWidth=1.5);

x1 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line

xlabel('t [ d ]');
ylabel('Q [ m^3/s ]');
h=legend('Q_{b,in}','Q_{b,relax}','Q_{out}','Q_{outb}','Q_{outQ}','Q_{outq}','Q_{in}','NumColumns',2);
h.Location='southwest';
text(0.025,0.8,'(a) flux','Units','normalized','FontSize',14)

xlim([tmin_d tmax_d])
set(gca, 'YScale', 'log')
ylim([1e1 1e4])
% yticks([1e2,1e3,1e4])
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
    plot(ax,t,p_b/1e6,'r-',LineWidth=1.5);
    hold on
    x5 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line
    xlabel('t [ d ]');
    ylabel('p_b [ MPa ]');
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
    plot(ax,t,V_b,'r-',LineWidth=1.5);
    hold on
    % plot(ax,t,Rb_analytical,'r--',LineWidth=1.5);
    ylim([0 1.1e9])
    ylabel('V_b [ m ]');
    legend('p_b','p_w','V_b','NumColumns',2,location='southwest')

    x6 = xline(tframe*ps.t/pd.td,'--k','LineWidth',1.5); % dashed line

%% right sublayout：3 行 2 列

% total flux
ax = nexttile(rightLayout);
zq = reshape(qnet,gg.nI,gg.nJ);
zphi = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ); 
pqnet = pcolor(ax,xx,yy,zq);
set(pqnet,'linestyle','none'); % shading interp
cx = colorbar();
colormap(parula)
clim([0 5e-1])
cx.Label.String = 'q_{all} [ m^2 s^{-1} ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 

hold on
[C4,pqnet_contour] = contour(ax,xx,yy,zphi,'linecolor','k','linewidth',0.5);

title('net flux','FontSize',14);
ylabel('y (km)')
axis equal

%% cavity sheet thickness
ax = nexttile(rightLayout);
zhs = (ps.hs)*reshape(vva.hs,gg.nI,gg.nJ); 
phs = pcolor(ax,xx,yy,zhs); 
set(phs,'linestyle','none'); % shading interp
cx = colorbar();
cx.Label.String = 'h_s [ m ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.3 0.5]; 
clim([0 0.1]);

hold on
[Cphi,phs_contour] = contour(ax,xx,yy,zphi,'linecolor','k','linewidth',0.5);

title('cavity sheet and \phi contour','FontSize',14);
ylabel('y (km)')
axis equal

%% blister sheet thickness
ax = nexttile(rightLayout);
zhe = (ps.hb)*reshape(vva.hb,gg.nI,gg.nJ); 
pblister = pcolor(ax,xx,yy,zhe); 
set(pblister,'linestyle','none'); % shading interp
cx = colorbar();
cx.Label.String = 'h_b [ m ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([-0.1 0.1]);
hold on

zpb = (ps.phi)*reshape(vva.pb,gg.nI,gg.nJ); 
[Cb,pblister_contour] = contour(ax,xx,yy,zpb,'linecolor','k','linewidth',0.5);

title('blister sheet and pb contour','FontSize',14);
ylabel('y (km)')
axis equal

%% cross-sectional area
ax = nexttile(rightLayout);
zS = (ps.S)*reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr) ,gg.nI,gg.nJ); 
pS = pcolor(ax,xx,yy,zS); 
set(pS,'linestyle','none');
% shading interp
cx = colorbar();
cx.Label.String = 'S [ m^2 ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
% clim([0 11]); 

title('channel cross section','FontSize',14);
ylabel('y (km)')
xlim([0 50])
ylim([0 10])
axis equal
axis tight

%% effective pressure
ax = nexttile(rightLayout);
zN = (ps.phi/10^6)*reshape((aa.phi_0-vva.phi),gg.nI,gg.nJ); 
peff = pcolor(ax,xx,yy,zN); 
set(peff,'linestyle','none'); % shading interp
hold on

% xlabel('x (km)')
cx = colorbar();
colormap(ax,cmap)
cx.Label.String = 'N [ MPa ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([-3 3]); 

% add moulins
if ~isfield(pp,'ni_m'), pp.ni_m = []; end
x = (ps.x/10^3)*gg.nx(pp.ni_m);
y = (ps.x/10^3)*gg.ny(pp.ni_m);
mscale = 100;
for i_m = 1:length(pp.ni_m)
    if aa.E(pp.ni_m(i_m))>0
        plot(x(i_m),y(i_m),'ko','Markersize',4+aa.E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
    else
        plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
    end
end

title('effective pressure','FontSize',14);
ylabel('y (km)')
axis equal
axis tight

%% effective pressure
ax = nexttile(rightLayout);
zpb = (ps.phi/10^6)*reshape((vva.pb-(aa.phi_0-aa.phi_a)),gg.nI,gg.nJ); 
ppb = pcolor(xx,yy,zpb); 
set(ppb,'linestyle','none');
% shading interp
cx = colorbar();
colormap(ax,cmap)
cx.Label.String = 'p_b [ MPa ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
pbmax = max(abs(vva.pb-(aa.phi_0-aa.phi_a)));
clim([-1 1]); 

time = ['t=' num2str(vva.t*ps.t/pd.td,'%.1f') ' d'];
ttext = text(0.1,8,time,"Position",[0.5 0.8],"Units","normalized");
ttext.FontSize=16;

title('blister pressure','FontSize',14);
ylabel('y (km)')
xlabel('x (km)')
axis equal
axis tight

%% make video
v = VideoWriter(['./results/videos/' casename '_2'],'MPEG-4');
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
    qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ + 0*vv2.Q);

    pqnet.CData = reshape(qnet,gg.nI,gg.nJ); 
    pqnet_contour.ZData = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ);

    phs.CData = (ps.hs)*reshape(vva.hs,gg.nI,gg.nJ);
    phs_contour.ZData = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ);

    pblister.CData = (ps.hb)*reshape(vva.hb,gg.nI,gg.nJ); 
    pblister_contour.ZData = (ps.phi)*reshape(vva.pb,gg.nI,gg.nJ);

    pS.CData = (ps.S)*reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr),gg.nI,gg.nJ);
    peff.CData = reshape(N,gg.nI,gg.nJ); 
    ppb.CData = reshape(pb,gg.nI,gg.nJ);

    set(ttext,{'string'},{['t=' num2str(vva.t*ps.t/(24*60*60),'%.1f'), ' d']})  %notice the column vector of new values
    
    x1.Value = vva.t*ps.t/pd.td;
    x2.Value = vva.t*ps.t/pd.td;
    x3.Value = vva.t*ps.t/pd.td;
    x4.Value = vva.t*ps.t/pd.td;
    x5.Value = vva.t*ps.t/pd.td;
    x6.Value = vva.t*ps.t/pd.td;
    disp(max(vva.hb));
    disp(min(vva.hb));
    refreshdata
    drawnow
    % pause(0.2)
    frame = getframe(gcf);
    writeVideo(v,frame)
end

close(v)
% end

% img = getframe(gcf);
% imwrite(img.cdata, ['./results/figures/' oo.casename, '_2dplot.png']);