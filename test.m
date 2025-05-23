%% Import necessary libraries
casename = 'nevis_regional_RACMO_kappa1e_11_kls0_mu1e1_c01e0';
casename = 'nevis_regional_RACMO_kappa1e_11_kls0_mu1e1_c01e0';
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                         % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code

tmin = 0*ps.t/pd.td; 
tmax = 40*ps.t/pd.td;                           % time range for the plot
t_init = 1450; t_end = 1600;                    % time range for animation

%% colormap
n = 256; % number of colors
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% read in the screenshot at the intial timestep
formatSpec = '%04d';
tframe = 21.0;
nframe = floor(tframe/0.1)+1;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ + vv2.Q);

xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

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

%% 右侧子 layout：3 行 2 列

% total flux
ax = nexttile(rightLayout);
zq = reshape(qnet,gg.nI,gg.nJ);
zphi = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ); 
pqnet = pcolor(ax,xx,yy,zq);
set(pqnet,'linestyle','none'); % shading interp
cx = colorbar();
colormap(parula)
clim([0 1])
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
clim([0 1.5]);
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
clim([0 10]); 

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
ttext = text(0.1,8,time,"Position",[0.7 0.8],"Units","normalized");
ttext.FontSize=16;

title('blister pressure','FontSize',14);
ylabel('y (km)')
xlabel('x (km)')
axis equal
axis tight

% save as png
% exportgraphics(f,[oo.rn '/nevis_regional_test_2009_140km_mu1e1_kappa1e_11_Vl1e8_td200_t210d.png'],'Resolution',600)