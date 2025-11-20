%% Import necessary libraries
casename = 'n2d_moulin1e1_eps1e_02_kappa1e_10_mu1e1_spinup';
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                         % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code
oo.t_span = (0:1:365*5)*pd.td/ps.t; 
plot_delta = 1;                                 % plot change from initial condition

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

%% Plot the velocity field with time
figure(2); clf;
% zu = (ps.x/ps.t*pd.ty)*reshape(vva.u,gg.nI,gg.nJ);
% zv = (ps.x/ps.t*pd.ty)*reshape(vva.v,gg.nI,gg.nJ); 
subplot(3,1,1);
pvel = pcolor(xx,yy,vUU); 
shading interp; colorbar; colormap(cmap);
caxis([-30 30]);
hold on;
vel_field = quiver(xx,yy,vux,vuy, 0.3,'k'); 

ttext = text(0.1*L/10^3,0.9*W/10^3,['Time = ' num2str(vva.t*ps.t/pd.td,'%.1f') ' days'],'FontSize',14,'Color','k');

xlabel('x (km)','FontSize',14);
ylabel('y (km)','FontSize',14);
title('Ice velocity change (m/yr)','FontSize',14);
set(gca,'FontSize',12);

subplot(3,1,2);
phs = pcolor(xx,yy,reshape(vva.hs*ps.h,gg.nI,gg.nJ)); 
shading interp; 
colorbar;colormap(cmap);
caxis([-0.1 0.1]);
clabels = colorbar;
clabels.Label.String = 'Sheet thickness change (m)';
xlabel('x (km)','FontSize',14);
ylabel('y (km)','FontSize',14);
title('Sheet thickness change (m)','FontSize',14);

subplot(3,1,3);
pphi = pcolor(xx,yy,(vphi0-vphi)/1e6); 
shading interp; 
colorbar;colormap(cmap);
caxis([-1 1]);
clabels = colorbar;
clabels.Label.String = 'Effective pressure change (MPa)';
hold on;
[Cb,phb] = contour(xx,yy,reshape(vva.hb*ps.hb,gg.nI,gg.nJ),[0.01, 0.1],'k');

xlabel('x (km)','FontSize',14);
ylabel('y (km)','FontSize',14);
title('Effective pressure change (MPa)','FontSize',14);

%% Animate over all timesteps
tmin = 0*365;
tmax = tmin+1.1*365;
index = find(oo.t_span>=tmin*pd.td/ps.t & oo.t_span<=tmax*pd.td/ps.t);
% n0 = index; % starting frame number
v = VideoWriter(['./test'],'MPEG-4');
v.FrameRate = 1;
open(v);
for nframe = index

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
    vel_field.UData = vux;
    vel_field.VData = vuy;
    pphi.CData = vN;
    phb.ZData = reshape(vva.hb*ps.hb,gg.nI,gg.nJ);
    ttext.String = ['Time = ' num2str(vva.t*ps.t/pd.td,'%.1f') ' days'];
    
    refreshdata
    drawnow
    % pause(0.2)
    frame = getframe(gcf);
%     writeVideo(v,frame);
end
% close(v);