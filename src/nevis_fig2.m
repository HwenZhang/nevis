% plots topography
% requires vv,aa,gg,pp,ps,oo
% 
% 21 August 21 : taken from nevis_plot [ similar to figure_arolla etc ]
% [ abandon ? can do same with nevis_plot (oo.topography) ]

% fn = '';
% load([fn,'/0000']); load(fn);
fn = [oo.root,oo.fn];

if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end
if isfield(oo,'halfcmap'), halfcmap = oo.halfcmap; else halfcmap = 0; end

%% elevation range
z_range = ps.z*[min(aa.b) max(aa.s)];
z_range = [-999 1000]; halfcmap = 0;

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

%% topography
axes(ax);
    xx = (ps.x/10^3)*gg.nx; yy = (ps.x/10^3)*gg.ny;  % keep outside points here to prevent edge effects disappearing
    zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=-400) = NaN;
    cax = z_range; if halfcmap, cax = [cax(1)-1.02*(cax(2)-cax(1)) cax(2)]; end
    sq = 0.02; cax = [cax(1)-sq*(cax(2)-cax(1)) cax(1)+(1+sq)*(cax(2)-cax(1))]; % expand cax if not going to use last entry
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
%     hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); shading interp;
    hand = imagesc(xx(:,1),yy(1,:),zz'); set(gca,'YDir','normal');
    caxis(cax); 
    axis image;
    axis(axx);
    
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     text(1.1,0,'m','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','left');
    text(-0.05,0,'z [ m ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');

%     % add bed contours
%     zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=-400) = NaN;
%     fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
%     xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
%     hold on; contour(xx(xrg,yrg),yy(xrg,yrg),zzz,b_conts,'color',0.5*[1 1 1]);

    % add surface contours
    zz = ps.z*reshape(s,gg.nI,gg.nJ);
    fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
    xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
    hold on; contour(xx(xrg,yrg),yy(xrg,yrg),zzz,s_conts,'color',0.5*[1 1 1],'linewidth',1);
%     contour(xx,yy,ps.z*reshape(s-b,gg.nI,gg.nJ),10*[1 1],'color','k','linewidth',1); % mark margin

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

    % add moulin labels
    if isfield(pp,'num_m'),
    for i_m = 1:length(pp.ni_m),
        label = num2str(pp.num_m(i_m));
        x_off = 0.01*(axx(2)-axx(1)); y_off = 0.01*(axx(4)-axx(3));
        text(x(i_m)+x_off,y(i_m)+y_off,label,'fontsize',6,'HorizontalAlignment','left','VerticalAlignment','bottom');
    end
    end
    
    % add outlets
    plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 

    if reversey, set(gca,'YDir','reverse'); end

%% colormap   
load cmapland2; cmap([1 end],:) = 0.5*[1 1 1; 1 1 1]; colormap(cmap); set(cx,'XLim',z_range);

%% axes labels 
axes(ax); 
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
   
print(gcf,'-depsc2',[fn,'_fig2']);