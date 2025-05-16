%% set up a figure
figure();

fs = 12;
set(gcf,'position',[0 0 800 800]);
casename = oo.casename;
oo.fn = ['/',oo.casename];                      % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
path = [oo.rn,'/'];

formatSpec = '%04d';
nframe = 751;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ + vv2.Q);
% effective pressure
N = (ps.phi/10^6)*(aa.phi_0-vva.phi);
pb = (ps.phi/10^6)*(vva.pb-(aa.phi_0-aa.phi_a));
pb_net = (ps.phi/10^6)*(vva.pb);
xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% total flux
ax(1) = subplot(6,1,1);
zq = reshape(qnet,gg.nI,gg.nJ);
zphi = (ps.phi)*reshape(vva.phi,gg.nI,gg.nJ); 
p4 = pcolor(xx,yy,zq);
set(p4,'linestyle','none'); %
% shading interp
hold on

ylabel('y (km)')
cx = colorbar();
clim([0 2])
cx.Label.String = 'q_{all} [ m^2 s^{-1} ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
[C4,p4_contour] = contour(xx,yy,zphi,'linecolor','r','linewidth',0.5);
% clim([0 0.1]);
axis equal

%% cavity sheet thickness
ax(2) = subplot(6,1,2);
zhs = (ps.hs)*reshape(vva.hs,gg.nI,gg.nJ); 
p1 = pcolor(xx,yy,zhs); 
set(p1,'linestyle','none'); % 
% shading interp
hold on
[Ce,pe_contour] = contour(xx,yy,zphi,'linecolor','r','linewidth',0.5);

ylabel('y (km)')
cx = colorbar();
cx.Label.String = 'h_s [ m ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([0 0.1]);
axis equal

%% elastic sheet thickness
ax(3) = subplot(6,1,3);
zhe = (ps.hb)*reshape(vva.hb,gg.nI,gg.nJ); 
pe = pcolor(xx,yy,zhe); 
set(pe,'linestyle','none'); % 
cx = colorbar();
cx.Label.String = 'h_b [ m ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
% clim([0 ps.hs*max([tt.pts_he])]);
clim([0 ps.hb*max(vva.hb)]);

zpb = (ps.phi)*reshape(vva.pb,gg.nI,gg.nJ); 
% shading interp
hold on
[C1,p1_contour] = contour(xx,yy,zpb,'linecolor','r','linewidth',0.5);
ylabel('y (km)')
axis equal

%% cross-sectional area
ax(4) = subplot(6,1,4);
zS = (ps.S)*reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr) ,gg.nI,gg.nJ); 
p2 = pcolor(xx,yy,zS); 
set(p2,'linestyle','none');
% shading interp

ylabel('y (km)')
cx = colorbar();
cx.Label.String = 'S [ m^2 ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([0 10]); 
xlim([0 50])
ylim([0 10])
axis equal
axis tight

%% effective pressure
ax(5) = subplot(6,1,5);
zN = reshape(N+1*pb,gg.nI,gg.nJ); 
p3 = pcolor(xx,yy,zN); set(p3,'linestyle','none');
% shading interp

ylabel('y (km)')
% xlabel('x (km)')
cx = colorbar();
cx.Label.String = 'N [ MPa ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([-3 3]); 

% add lakes
mscale = 100;  % amount to scale moulin input by to make size of moulin dot
if ~isfield(pp,'ni_l'), pp.ni_l = []; end
vbscale = 1e5;  % amount to scale moulin input by to make size of moulin dot
x = (ps.x/10^3)*gg.nx(pp.ni_l);
y = (ps.x/10^3)*gg.ny(pp.ni_l);
hold on
% sca = scatter(x,y,ps.h*vva.Vb(pp.ni_l)./mscale,'MarkerFaceColor',[1, 1, 1],'MarkerFaceAlpha',.3); % mark moulins

% add moulins
if ~isfield(pp,'ni_m'), pp.ni_m = []; end
x = (ps.x/10^3)*gg.nx(pp.ni_m);
y = (ps.x/10^3)*gg.ny(pp.ni_m);
for i_m = 1:length(pp.ni_m),
    if aa.E(pp.ni_m(i_m))>0,
        plot(x(i_m),y(i_m),'ko','Markersize',4+aa.E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
    else
        plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
    end
end
axis equal
axis tight

%% effective pressure
ax(6) = subplot(6,1,6);
zpb = reshape(pb,gg.nI,gg.nJ); 
p6 = pcolor(xx,yy,zpb); 
set(p6,'linestyle','none');
% shading interp

ylabel('y (km)')
xlabel('x (km)')
cx = colorbar();
cx.Label.String = 'p_b [ MPa ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([-max(pb) max(pb)]); 

time = ['t=' num2str(vva.t*ps.t/pd.td,'%.1f') ' d'];
ttext = text(0.1,8,time,"Position",[2 2]);
ttext.FontSize=16;

axis equal
axis tight

img = getframe(gcf);
imwrite(img.cdata, ['./results/figures/' oo.casename, '_2dplot.png']);