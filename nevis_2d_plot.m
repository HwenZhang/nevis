% function nevis_2d_plot(tt,vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]


%% set up a figure
figure();

fs = 12;
path = '/Users/eart0487/nevis/nevis_2d_example/';
formatSpec = '%04d';
nframe = 752;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
L = 5e4;                               % length of the domain [m]

xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  

%% cavity sheet thickness
ax(1) = subplot(3,1,1);
% gg = load([path num2str(nframe,formatSpec)], 'gg');
zz = (ps.hs)*reshape(vva.hs,gg.nI,gg.nJ); 
hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % 
shading interp;
ylabel('y (km)')
cx = colorbar();
cx.Label.String = 'h_s [ m ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
axis equal
% clim([0 10]);

ax(2) = subplot(3,1,2);
zz = (ps.S)*reshape(0.25*(gg.nmeanx*vva.Sx + gg.nmeany*vva.Sy + gg.nmeans*vva.Ss + gg.nmeanr*vva.Sr) ,gg.nI,gg.nJ); 
hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % 
shading interp;
ylabel('y (km)')

cx = colorbar();
cx.Label.String = 'S [ m^2 ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
axis equal
% clim([0 50]); 

if ~isfield(pp,'ni_m'), pp.ni_m = []; end
    mscale = 100;  % amount to scale moulin input by to make size of moulin dot
    x = (ps.x/10^3)*gg.nx(pp.ni_m);
    y = (ps.x/10^3)*gg.ny(pp.ni_m);
    hold on;
    for i_m = 1:length(pp.ni_m),
        if aa.E(pp.ni_m(i_m))>0,
            plot(x(i_m),y(i_m),'ko','Markersize',4+aa.E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
        else
            plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
        end
    end

ax(3) = subplot(3,1,3);
zz = (ps.phi/10^6)*reshape(aa.phi_0-vva.phi,gg.nI,gg.nJ); 
hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % 
shading interp;
ylabel('y (km)')
xlabel('x (km)')
cx = colorbar();
cx.Label.String = 'N [ MPa ]'; 
cx.Label.Units = 'normalized'; 
cx.Label.Position = [2.2 0.5]; 
clim([-3 3]); 
axis equal

time = ['t=' num2str(vva.t*ps.t/(24*60*60),'%.1f') ' d'];
ttext = text(0.1,8,time);
ttext.FontSize=12;

% ax(4) = subplot(4,1,4);
% total flux

% add lakes
if ~isfield(pp,'ni_l'), pp.ni_l = []; end
mscale = 1e5;  % amount to scale moulin input by to make size of moulin dot
x = (ps.x/10^3)*gg.nx(pp.ni_l);
y = (ps.x/10^3)*gg.ny(pp.ni_l);
hold on;
scatter(x,y,ps.V*vva.Vb(pp.ni_l)./mscale,'MarkerFaceColor',[1, 1, 1],'MarkerFaceAlpha',.3); % mark moulins

% add moulins
if ~isfield(pp,'ni_m'), pp.ni_m = []; end
mscale = 100;  % amount to scale moulin input by to make size of moulin dot
x = (ps.x/10^3)*gg.nx(pp.ni_m);
y = (ps.x/10^3)*gg.ny(pp.ni_m);
hold on;
for i_m = 1:length(pp.ni_m),
    if aa.E(pp.ni_m(i_m))>0,
        plot(x(i_m),y(i_m),'ko','Markersize',4+aa.E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
    else
        plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
    end
end
% end

img = getframe(gcf);
imwrite(img.cdata, ['./figures/' oo.casename, '_2dplot.png']);