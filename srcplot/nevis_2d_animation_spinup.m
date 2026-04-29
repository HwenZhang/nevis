%% Import necessary libraries
% casename = oo.casename;
casename = 'n2d_region_ice_Cinv_test';  % specify the case name
oo.plot_residual = 0;
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
xx0 = xx;
yy0 = yy;
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% read in the time series
t = (ps.t/(24*60*60))*[tt.t];               % dimensional time series (days)
tspan = (ps.t/pd.td)*oo.t_span;
tmin = 0.55*365*pd.td/ps.t;
tmax = 1.0*365*pd.td/ps.t;
tmin_d = tmin*ps.t/pd.td; 
tmax_d = tmax*ps.t/pd.td;                   % time range for the plot

dt_plot = 1.0;
t_plot = tmin_d:dt_plot:tmax_d;

frame_indices = zeros(size(t_plot));
for k = 1:length(t_plot)
    [~, idx] = min(abs(tspan - t_plot(k)));
    frame_indices(k) = idx;
end
frame_indices = unique(frame_indices); % 去重，避免重复帧

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
R_b = ps.x*[tt.Rb];

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

f = figure('Position', [100, 100, 1200, 900]);
rightLayout = tiledlayout(f, 3, 2);
rightLayout.TileSpacing = 'compact';
rightLayout.Padding = 'compact';

%% right sublayout

% total flux
ax = nexttile(rightLayout);
zq = reshape(qnet,gg.nI,gg.nJ);
zphi = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ); 
pqnet = pcolor(ax,xx,yy,zq);
set(pqnet,'linestyle','none'); % shading interp
ax.ColorScale = 'log';
cx = colorbar();
colormap(parula)
clim([1e-6 5e-1])
cx.Label.String = 'q_{all} [ m^2 s^{-1} ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 

hold on
[C4,pqnet_contour] = contour(ax,xx,yy,zphi,'linecolor','k','linewidth',0.5);
% quiver: water flows down the gradient of phi
xx_grid = reshape(xx0, gg.nI, gg.nJ);
yy_grid = reshape(yy0, gg.nI, gg.nJ);
nout_mask = reshape(false(size(xx0)), gg.nI, gg.nJ);
nout_mask(gg.nout) = true;  % exterior nodes
zphi_clean = zphi;
zphi_clean(nout_mask) = NaN;
% fill NaN so gradient doesn't propagate them to interior neighbors
zphi_fill = fillmissing(zphi_clean, 'nearest');
% ndgrid: dim1=x (rows), dim2=y (cols)
% gradient(F, h_col, h_row) returns [dF/dcol, dF/drow] = [dF/dy, dF/dx]
xvec = xx_grid(:,1);  % x varies along dim 1
yvec = yy_grid(1,:);  % y varies along dim 2
[dphidy, dphidx] = gradient(zphi_fill, yvec, xvec);
qskip = 6; % subsample for readability
idx_q = 1:qskip:gg.nI;
jdx_q = 1:qskip:gg.nJ;
qmag = sqrt(dphidx.^2 + dphidy.^2);
qmag(qmag < eps) = 1;
qu = -dphidx(idx_q,jdx_q)./qmag(idx_q,jdx_q);
qv = -dphidy(idx_q,jdx_q)./qmag(idx_q,jdx_q);
% hide arrows at exterior nodes
nout_sub = nout_mask(idx_q, jdx_q);
qu(nout_sub) = NaN;
qv(nout_sub) = NaN;
hquiver = quiver(ax, xx_grid(idx_q,jdx_q), yy_grid(idx_q,jdx_q), ...
    qu, qv, 0.6, 'w', 'LineWidth', 1.0);
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
shading interp
set(pblister,'linestyle','none'); % shading interp
cx = colorbar();
cx.Label.String = 'h_b [ m ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([-0.01 0.01]);
% highlight all the nodes where hb<-0.005
hold on
% [C_neg, h_neg] = contour(ax, xx, yy, zhe, [-0.005 -0.005], 'r-', 'LineWidth', 2);

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
ax.ColorScale = 'log';

cx = colorbar();
cx.Label.String = 'S [ m^2 ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([1e-4 1e1]);
cx.Ticks = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];

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
% if ~isfield(pp,'ni_m'), pp.ni_m = []; end
% x = (ps.x/10^3)*gg.nx(pp.ni_m);
% y = (ps.x/10^3)*gg.ny(pp.ni_m);
% mscale = 100;
% for i_m = 1:length(pp.ni_m)
%     if aa.E(pp.ni_m(i_m))>0
%         plot(x(i_m),y(i_m),'ko','Markersize',4+aa.E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
%     else
%         plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
%     end
% end

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
v = VideoWriter(['./results/videos/' casename],'MPEG-4');
% v = VideoWriter(['./results/' oo.casename],'MPEG-4');
v.FrameRate = 3;
open(v)

for i_idx = 1:length(frame_indices)
    if ~isvalid(f)
        disp('Figure closed by user. Animation stopped.');
        break
    end
    i_t = frame_indices(i_idx);
    disp(['Frame ',num2str(i_idx),' / ',num2str(length(frame_indices)),' ...']);
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

    % update quiver arrows
    zphi_new = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ);
    zphi_new(nout_mask) = NaN;
    zphi_fill = fillmissing(zphi_new, 'nearest');
    [dphidy, dphidx] = gradient(zphi_fill, yvec, xvec);
    qmag = sqrt(dphidx.^2 + dphidy.^2);
    qmag(qmag < eps) = 1;
    qu = -dphidx(idx_q,jdx_q)./qmag(idx_q,jdx_q);
    qv = -dphidy(idx_q,jdx_q)./qmag(idx_q,jdx_q);
    qu(nout_sub) = NaN;
    qv(nout_sub) = NaN;
    hquiver.UData = qu;
    hquiver.VData = qv;

    phs.CData = (ps.hs)*reshape(vva.hs,gg.nI,gg.nJ);
    phs_contour.ZData = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ);

    % h_neg.ZData = (ps.hb)*reshape(vva.hb,gg.nI,gg.nJ);
    pblister.CData = (ps.hb)*reshape(vva.hb,gg.nI,gg.nJ); 
    pblister_contour.ZData = (ps.phi)*reshape(vva.pb,gg.nI,gg.nJ);

    pS.CData = (ps.S)*reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr),gg.nI,gg.nJ);
    peff.CData = reshape(N,gg.nI,gg.nJ); 
    ppb.CData = reshape(pb,gg.nI,gg.nJ);

    set(ttext,{'string'},{['t=' num2str(vva.t*ps.t/(24*60*60),'%.1f'), ' d']})  %notice the column vector of new values

    disp(vva.t*ps.t/pd.td);
    refreshdata
    drawnow
    % pause(0.2)
    frame = getframe(gcf);
    writeVideo(v,frame)
end

% close VideoWriter if still open
try
    close(v);
catch
    % video writer already closed or not open
end

% img = getframe(gcf);
% imwrite(img.cdata, ['./results/figures/' oo.casename, '_2dplot.png']);