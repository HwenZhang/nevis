% plots topography and discharge
% requires vv,aa,gg,pp,ps,oo
% 
% 31 August: taken from nevis_fig1, combining plots of topography and q

% fn = '';
% load([fn,'/0000']); load(fn);
fn = [oo.root,oo.fn];

if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end
if isfield(oo,'halfcmap'), halfcmap = oo.halfcmap; else halfcmap = 0; end
if strcmp(version,'8.4.0.150421 (R2014b)'), ver = 2; else ver = 1; end % matlab version for text input

%% elevation range
z_range = ps.z*[min(aa.b) max(aa.s)];
% z_range = [-1499 1500];
z_range = [-999 1000];
z_range = [-2999 3000];

%% contour levels
z_conts = -2000:50:5000;
z_conts = -3000:100:3000;
b_conts = z_conts;
s_conts = z_conts;
p_conts = z_conts/100;

%% axes limits
axx = (ps.x/10^3)*[gg.xl gg.xr gg.yb gg.yt];

%% extract variables
if ~isfield(vv,'t'), vv.t = 0; end
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,gg,vv2);

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
axes_size = [0.1 0.1 0.8 0.8;
            0.1 0.1 0.8 0.8];

%adjust axes positions to account for aspect ratio of figure
aspect = (axx(4)-axx(3))/(axx(2)-axx(1)); % aspect of axes limits
for ai = 1:size(axes_size,1);
pos = axes_size(ai,:); tmp= get(gcf,'position'); sca = tmp(4)/tmp(3);
if aspect>(pos(4)/pos(3))*sca, 
    pos = pos + [pos(4)*sca*(1-1/aspect)/2 0 -pos(3)+pos(4)*sca*(1/aspect) 0]; % constrained by y
else
    pos = pos + [0 pos(3)/sca*(1-aspect)/2 0 -pos(4)+pos(3)/sca*aspect]; % constrained by x
end 
axes_size(ai,:) = pos;
end

%% topography
ai = 1;
ax(ai) = axes('position',axes_size(ai,:));
    xx = (ps.x/10^3)*gg.nx; yy = (ps.x/10^3)*gg.ny;  % keep outside points here to prevent edge effects disappearing
    zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=z_range(1)) = NaN; % zz(isnan(zz)) = z_range(1);
    cax = z_range; if halfcmap, cax = [cax(1)-1.02*(cax(2)-cax(1)) cax(2)]; end
    sq = 0.02; cax = [cax(1)-sq*(cax(2)-cax(1)) cax(1)+(1+sq)*(cax(2)-cax(1))]; % expand cax if not going to use last entry
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
%     hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); shading interp;
    hand = imagesc(xx(:,1),yy(1,:),zz'); set(gca,'YDir','normal');
    caxis(cax); 
    axis image;
    axis(axx);
    
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    ccx(ai) = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     text(1.1,0,'m','units','normalized','Parent',cx(ai),'VerticalAlignment','bottom','HorizontalAlignment','left');
    if ver==1,
    text(-0.05,0,'z [ m ]','units','normalized','Parent',ccx(ai),'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    ccx(ai).Label.String = 'z [ m ]'; ccx(ai).Label.Units = 'normalized'; ccx(ai).Label.Position = [-0.05 0]; ccx(ai).Label.VerticalAlignment = 'bottom'; ccx(ai).Label.HorizontalAlignment = 'right';
    end
    
%     % add bed contours
%     zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=-400) = NaN;
%     fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
%     xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
%     hold on; contour(xx(xrg,yrg),yy(xrg,yrg),zzz,b_conts,'color',0.5*[1 1 1]);

    % add surface contours
    s(gg.nout) = NaN;
    zz = ps.z*reshape(s,gg.nI,gg.nJ);
    fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
    xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
    hold on; contour(xx(xrg,yrg),yy(xrg,yrg),zzz,s_conts,'color',0.5*[1 1 1],'linewidth',1);
%     contour(xx,yy,ps.z*reshape(s-b,gg.nI,gg.nJ),1*[1 1],'color','k','linewidth',1); % mark margin

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
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 

    if reversey, set(gca,'YDir','reverse'); end

%% discharge
ai = 2;
ax(ai) = axes('position',axes_size(ai,:));
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny; 
    zz = ps.qs*reshape(qs+qe+qQ,gg.nI,gg.nJ);   
    cax = [10^(-5) 10^(-2)]; 
    cax = [10^(-7) 10^(-4)]; 

    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    alpha(0.5); % transparency
    
%     % add pressure contours
%     hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
    if reversey, set(gca,'YDir','reverse'); end
    axis off;
    
%% colormaps    
load cmapland2; cmap([1 end],:) = 0.5*[1 1 1; 1 1 1]; 
cmaps(1).cmap = cmap;
load cmapq2; cmaps(2).cmap = cmap;
% load cmapgray; cmaps(2).cmap = cmap;
% load cmapred; cmaps(2).cmap = cmap;
[CLim] = splitcmapn([ax(1) ax(2)],cmaps);
if ver==1,
set(ccx(1),'XLim',CLim(1,:)); % set(cx(2),'XLim',CLim(2,:));
set(ccx(1),'XLim',z_range);
else
ccx(1).XLim = CLim(1,:); ccx(1).XLim = z_range;
end
    
%% axes labels 
axes(ax(1)); 
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    
axes(ax(2)); 
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    
print(gcf,'-depsc2',[fn,'_fig1a'],'-opengl','-loose');
 