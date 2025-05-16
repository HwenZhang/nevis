% plots discharge
% requires vv,aa,gg,pp,ps,oo
% 
% 21 August 21 : taken from nevis_plot [ similar to figure_arolla etc ]
% [ abandon ? can do same with nevis_plot (oo.topography) ]

% fn = '';
% load([fn,'/0000']); load(fn);
fn = [oo.root,oo.fn];

if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end
if isfield(oo,'halfcmap'), halfcmap = oo.halfcmap; else halfcmap = 0; end

%% contour levels
z_conts = -1000:50:5000;
b_conts = z_conts;
s_conts = z_conts;
p_conts = z_conts/100;

%% extract variables
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
oo.evaluate_variables = 1; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,gg,vv2);

%% axes limits
axx = (ps.x/10^3)*[gg.xl gg.xr gg.yb gg.yt];

%get rid of points outside domain
nx(gg.nout) = NaN;
ex(gg.eout) = NaN;
fx(gg.fout) = NaN;
cx(gg.cout) = NaN;
    
%boundary curve
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary

%% set up figure
figure(1); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 15 15]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
ax = axes('position',[0.1 0.1 0.8 0.8]);

%% adjust axes positions to account for aspect ratio of figure
aspect = (axx(4)-axx(3))/(axx(2)-axx(1)); % aspect of axes limits
for ai = 1:length(ax);
pos = get(ax(ai),'position'); tmp= get(gcf,'position'); sca = tmp(4)/tmp(3);
if aspect>(pos(4)/pos(3))*sca, 
    pos = pos + [pos(4)*sca*(1-1/aspect)/2 0 -pos(3)+pos(4)*sca*(1/aspect) 0]; % constrained by y
else
    pos = pos + [0 pos(3)/sca*(1-aspect)/2 0 -pos(4)+pos(3)/sca*aspect]; % constrained by x
end 
set(ax(ai),'Position',pos);
end

%% dischrage
axes(ax);
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny; 
    zz = ps.qs*reshape(qs+qe+qQ,gg.nI,gg.nJ);   
    cax = [10^(-5) 10^(0)]; 
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    tmp = get(cx(ai),'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
%     text(1.1,0,'m^2 s^{-1}','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','left');
    text(-0.05,0,'q [ m^2 s^{-1} ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');

    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
    % add moulins
    if ~isfield(pp,'ni_m'), pp.ni_m = []; end
    mscale = 100;  %amount to scale moulin input by to make size of moulin dot
    x = (ps.x/10^3)*nx(pp.ni_m);
    y = (ps.x/10^3)*ny(pp.ni_m);
    hold on;
    for i_m = 1:length(pp.ni_m),
        if E(pp.ni_m(i_m))>0,
            plot(x(i_m),y(i_m),'ko','Markersize',4+E(pp.ni_m(i_m))/mscale,'MarkerFaceColor',1*[1 1 1]); % mark moulins   
        else
            plot(x(i_m),y(i_m),'ko','Markersize',4,'MarkerFaceColor',0.8*[1 1 1]); % mark moulins  
        end
    end
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 

    if reversey, set(gca,'YDir','reverse'); end

%% colormap   
load cmapq2; colormap(cmap);

%% axes labels 
axes(ax); 
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
   
print(gcf,'-depsc2',[fn,'_fig3']);