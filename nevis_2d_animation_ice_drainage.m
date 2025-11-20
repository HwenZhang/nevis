%% Import necessary libraries
casename = 'n2d_moulin0e00_eps1e_01_kappa1e_10_mu1e1_V1e7_drainage';
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                         % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code
plot_delta = 1;                                 % plot change from initial condition

L = 100e3;
W = 50e3;
%% colormap
n = 256; % number of color
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% read in the screenshot at the intial timestep
formatSpec = '%04d';
nframe = 1;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny; 
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% Calculate node means of ice velocities
uxn = gg.nmeanx(:,gg.es2)*vva.u(gg.es2); % x-component of ice velocity at nodes
vyn = gg.nmeany(:,gg.fs2)*vva.v(gg.fs2); % y-component of ice velocity at nodes
Un = ((uxn).^2+(vyn).^2).^(1/2); % ice speed at nodes
Un(gg.nout) = NaN;

vux = ps.u_b*pd.ty*reshape(uxn,gg.nI,gg.nJ); % x-component of ice velocity at nodes in m/yr
vuy = ps.u_b*pd.ty*reshape(vyn,gg.nI,gg.nJ); % y-component of ice velocity at nodes in m/yr
vUU = ps.u_b*pd.ty*reshape(Un,gg.nI,gg.nJ);   % ice speed at nodes in m/yr

vphi = reshape(vva.phi*ps.phi,gg.nI,gg.nJ);
vphi0 = reshape(aa.phi_0*ps.phi,gg.nI,gg.nJ);

vux_init = vux;
vuy_init = vuy;
vUU_init = vUU;
vphi_init = vphi;
vhs_init = reshape(vva.hs*ps.h,gg.nI,gg.nJ);

skip1 = 2;
skip2 = 2;

%% Plot the velocity field with time
f=figure(); clf;
f.Position = [100 100 1200 1200];
% velocity map
subplot(3,2,1);
pvel = pcolor(xx,yy,vUU); 
shading interp; colorbar; colormap(cmap);
caxis([-20 20]);
clabels = colorbar;
clabels.Label.String = 'Ice velocity (m/yr)';
hold on;
skip = skip1;
vel_field = quiver(xx(1:skip:end, 1:skip:end),...
                   yy(1:skip:end, 1:skip:end),...
                   vux(1:skip:end, 1:skip:end),...
                   vuy(1:skip:end, 1:skip:end), 0,'k');
xlabel('x (km)','FontSize',10);
ylabel('y (km)','FontSize',10);
title('Ice velocity (m/yr)','FontSize',10);
set(gca,'FontSize',10);
axis equal;

% velocity delta map
subplot(3,2,2);
pvel_delta = pcolor(xx,yy,vUU); 
shading interp; colorbar; colormap(cmap);
caxis([-20 20]);
clabels = colorbar;
clabels.Label.String = 'Ice velocity change (m/yr)';
hold on;
skip = skip2;
vel_field_delta = quiver(xx(1:skip:end, 1:skip:end),...
                         yy(1:skip:end, 1:skip:end),...
                         vux(1:skip:end, 1:skip:end),...
                         vuy(1:skip:end, 1:skip:end), 0.3,'k'); 

ttext = text(0.1*L/10^3,0.9*W/10^3,['Time = ' num2str(vva.t*ps.t/pd.td,'%.1f') ' days'],'FontSize',10,'Color','k');
xlabel('x (km)','FontSize',10);
ylabel('y (km)','FontSize',10);
title('Ice velocity change (m/yr)','FontSize',10);
set(gca,'FontSize',10);
axis equal;

% sheet thickness map
subplot(3,2,3);
phs = pcolor(xx,yy,reshape(vva.hs*ps.h,gg.nI,gg.nJ)); 
shading interp; 
colorbar;colormap(cmap);
caxis([-0.1 0.1]);
clabels = colorbar;
clabels.Label.String = 'Sheet thickness (m)';
xlabel('x (km)','FontSize',10);
ylabel('y (km)','FontSize',10);
title('Sheet thickness (m)','FontSize',10);
axis equal;

% sheet thickness delta map
subplot(3,2,4);
phs_delta = pcolor(xx,yy,reshape(vva.hs*ps.h,gg.nI,gg.nJ)); 
shading interp; 
colorbar;colormap(cmap);
caxis([-0.05 0.05]);
clabels = colorbar;
clabels.Label.String = 'Sheet thickness change (m)';
xlabel('x (km)','FontSize',10);
ylabel('y (km)','FontSize',10);
title('Sheet thickness change (m)','FontSize',10);
axis equal;

% effective pressure map
subplot(3,2,5);
pphi = pcolor(xx,yy,(vphi0-vphi)/1e6); 
shading interp; 
colorbar;colormap(cmap);
clim([-2 2]);
clabels = colorbar;
clabels.Label.String = 'Effective pressure (MPa)';
hold on;
[Cb,phb] = contour(xx,yy,reshape(vva.hb*ps.hb,gg.nI,gg.nJ),[0.01 0.01],'k');
clabel(Cb, phb);
xlabel('x (km)','FontSize',10);
ylabel('y (km)','FontSize',10);
title('Effective pressure (MPa)','FontSize',10);
axis equal;

subplot(3,2,6);
pphi_delta = pcolor(xx,yy,(vphi0-vphi)/1e6); 
pphi_delta = pcolor(xx,yy,ps.hb*reshape(vva.hb,gg.nI,gg.nJ)); 
shading interp; 
colorbar;colormap(cmap);
clim([-0.5 0.5]);
clabels = colorbar;
clabels.Label.String = 'Effective pressure change (MPa)';
hold on;
[Cb,phb] = contour(xx,yy,reshape(vva.hb*ps.hb,gg.nI,gg.nJ),[0.005 0.005],'k','LineWidth',1);
% clabel(Cb, phb);
xlabel('x (km)','FontSize',10);
ylabel('y (km)','FontSize',10);
title('Effective pressure change (MPa)','FontSize',10);
axis equal;

%% Animate over all timesteps
tmin = 10*365 + 0.19*365;
tdrain = 10*365 + 0.2*365;
tmax = tmin+0.15*365;
index = find(oo.t_span>=tmin*pd.td/ps.t & oo.t_span<=tmax*pd.td/ps.t);
% n0 = index; % starting frame number
v = VideoWriter(['./presentation/animation/' casename],'MPEG-4');
v.FrameRate = 1;
open(v);
for nframe = index(1:4:end)

    vva = load([path num2str(nframe,formatSpec)], 'vv');
    vva = vva.vv;
    uxn = gg.nmeanx(:,gg.es2)*vva.u(gg.es2); % x-component of ice velocity at nodes
    vyn = gg.nmeany(:,gg.fs2)*vva.v(gg.fs2); % y-component of ice velocity at nodes
    Un = ((uxn).^2+(vyn).^2).^(1/2); % ice speed at nodes
    Un(gg.nout) = NaN;  

    vux = ps.u_b*pd.ty*reshape(uxn,gg.nI,gg.nJ); % x-component of ice velocity at nodes in m/yr
    vuy = ps.u_b*pd.ty*reshape(vyn,gg.nI,gg.nJ); % y-component of ice velocity at nodes in m/yr

    vUU = (vux.^2+vuy.^2).^(1/2);   % ice speed at nodes in m/yr
    vphi = reshape(vva.phi*ps.phi/1e6,gg.nI,gg.nJ);
    vphi_0 = reshape(aa.phi_0*ps.phi/1e6,gg.nI,gg.nJ);
    vhs = reshape(vva.hs*ps.h,gg.nI,gg.nJ);
    vN = vphi_0-vphi;

    if nframe==index(1)
        vN0 = vphi_0-vphi;
        vux0 = vux;
        vuy0 = vuy;
        vUU0 = vUU;
        vhs0 = vhs;
    end
    
    phs.CData = vhs;
    pvel.CData = vUU;
    pphi.CData = vN;
    phs_delta.CData = vhs - vhs0;
    pvel_delta.CData = vUU - vUU0;
    % pphi_delta.CData = vN - vN0;
    pphi_delta.CData = ps.hb*reshape(vva.hb,gg.nI,gg.nJ);
    phb.ZData = reshape(vva.hb*ps.hb,gg.nI,gg.nJ);

    ttext.String = ['Time = ' num2str(vva.t*ps.t/pd.td-tdrain,'%.1f') ' days'];
    
    % vel_centerline.YData = vux(:,centerline_j);

    delete(vel_field);
    subplot(3,2,1);
    hold on;
    scale_factor = 0.5; % adjust the scale of the arrows
    skip = skip1;
    vel_field = quiver(xx(1:skip:end, 1:skip:end),...
                     yy(1:skip:end, 1:skip:end), ...
                     vux(1:skip:end, 1:skip:end)*scale_factor,...
                     vuy(1:skip:end, 1:skip:end)*scale_factor, 0, "Color", 'k');

    delete(vel_field_delta);
    subplot(3,2,2);
    hold on;
    skip = skip2;
    scale_factor = 0.5; % adjust the scale of the arrows
    vel_field_delta = quiver(xx(1:skip:end, 1:skip:end),...
                             yy(1:skip:end, 1:skip:end),...
                            (vux(1:skip:end, 1:skip:end) - vux0(1:skip:end, 1:skip:end))*scale_factor,...
                            (vuy(1:skip:end, 1:skip:end) - vuy0(1:skip:end, 1:skip:end))*scale_factor, 0, "Color", 'k');
    refreshdata
    drawnow
    % pause(0.2)
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

test;