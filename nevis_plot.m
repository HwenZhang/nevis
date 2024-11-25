%function nevis_plot(vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]
% nevis_plot(vv,aa,pp,ps,gg,oo)
% comment first line to use as script taking inputs from current workspace
% plots discharge by default
% 
% 21 August 2014: taken from figure_tyne, similar to figure_arolla etc 
% 29 Jan 2015: edited labelling for new matlab version

%% options
fn = [oo.root,oo.fn];
if isfield(oo,'save_plot'), save_plot = oo.save_plot; else save_plot = 0; end
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end
% if strcmp(version,'8.4.0.150421 (R2014b)'), ver = 2; else ver = 1; end % matlab version for text input
ver = 2; % default to version 2 [2 Oct 2019]

%% figures to plot
if isfield(oo,'discharge'), discharge = oo.discharge; else discharge = 1; end
if isfield(oo,'topography'), topography = oo.topography; else topography = 0; end
if isfield(oo,'discharge_lines'), discharge_lines = oo.discharge_lines; else discharge_lines = 0; end
if isfield(oo,'velocity'), velocity = oo.velocity; else velocity = 0; end
if isfield(oo,'input'), input = oo.input; else input = 0; end
if isfield(oo,'dissipation'), dissipation = oo.dissipation; else dissipation = 0; end
if isfield(oo,'thickness'), thickness = oo.thickness; else thickness = 0; end
if isfield(oo,'sheet'), sheet = oo.sheet; else sheet = 0; end
if isfield(oo,'elastic'), elastic = oo.elastic; else elastic = 0; end
if isfield(oo,'area'), area = oo.area; else area = 0; end
if isfield(oo,'pressure'), pressure = oo.pressure; else pressure = 0; end

figis = [];
if discharge, figis = [figis 1]; end
if discharge_lines, figis = [figis 2]; end
if velocity, figis = [figis 3]; end
if pressure, figis = [figis 4]; end
if sheet, figis = [figis 5]; end
if elastic, figis = [figis 6]; end
if area, figis = [figis 7]; end
if topography, figis = [figis 8]; end
if dissipation, figis = [figis 9]; end
if thickness, figis = [figis 10]; end
if input, figis = [figis 11]; end
% figis',

%% elevation range
z_range = [-999 1000];
% z_range = [2000 3000];
z_range = [-1999 2000];

%% contour levels
z_conts = -1000:100:5000;
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
if ~isempty(gg.n1),
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary
else x_out = []; y_out = [];
end

for figi = figis,

%% set up figure
figure(figi); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5+5*figi 2 15 15]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.8 0.8];

%adjust axes positions to account for aspect ratio of figure
aspect = (axx(4)-axx(3))/(axx(2)-axx(1)); % aspect of axes limits
pos = axes_size; tmp = get(gcf,'position'); sca = tmp(4)/tmp(3);
if aspect>(pos(4)/pos(3))*sca, 
    pos = pos + [pos(4)*sca*(1-1/aspect)/2 0 -pos(3)+pos(4)*sca*(1/aspect) 0]; % constrained by y
else
    pos = pos + [0 pos(3)/sca*(1-aspect)/2 0 -pos(4)+pos(3)/sca*aspect]; % constrained by x
end 
axes_size = pos;

%% plot
if figi==1, %discharge
axes('position',axes_size);

    % discharge
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;
    zz = ps.qs*reshape(qs+qe+qQ,gg.nI,gg.nJ);   
    cax = [10^(-6) 10^(-3)]; 
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
    if ver==1
    text(-0.05,0,'q [ m^2 s^{-1} ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2
    cx.Label.String = 'q [ m^2 s^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
    % add pressure contours
%     hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
    % add smoothed pressure contours
    zz = (ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ);
    fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
    xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
    hold on; hand = contour((ps.x/10^3)*nx(xrg,yrg),(ps.x/10^3)*ny(xrg,yrg),zzz,p_conts,'linecolor','r','linewidth',1);
    
    % add moulins
    if ~isfield(pp,'ni_m'), pp.ni_m = []; end
    mscale = 1e3;  %amount to scale moulin input by to make size of moulin dot
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
%     plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
    
    load cmapq2; colormap(cmap);
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_discharge']); end
end

if figi==2, %discharge_lines
clf; axes('position',axes_size);
    hold on;
    
    threshs = [1 10 100]; % threshold values for plotting channels
    widths = [1 3 5]; % widths corresponding to thresholds
    
    disp('Plotting channels as lines (may take a while) ...');
    for iQ = 1:length(threshs),
    tmp = (ps.Q*abs(Qx)>threshs(iQ)); nis = gg.econnect(tmp,[1 2]);
    plot((ps.x/10^3)*gg.nx(nis)',(ps.x/10^3)*gg.ny(nis)','b','linewidth',widths(iQ));
    tmp = (ps.Q*abs(Qy)>threshs(iQ)); nis = gg.fconnect(tmp,[1 2]);
    plot((ps.x/10^3)*gg.nx(nis)',(ps.x/10^3)*gg.ny(nis)','b','linewidth',widths(iQ));
    tmp = (ps.Q*abs(Qs)>threshs(iQ)); nis = [gg.fconnect(gg.cconnect(tmp,1),1) gg.fconnect(gg.cconnect(tmp,2),2)];
    plot((ps.x/10^3)*gg.nx(nis)',(ps.x/10^3)*gg.ny(nis)','b','linewidth',widths(iQ));
    tmp = (ps.Q*abs(Qr)>threshs(iQ)); nis = [gg.fconnect(gg.cconnect(tmp,1),2) gg.fconnect(gg.cconnect(tmp,2),1)];
    plot((ps.x/10^3)*gg.nx(nis)',(ps.x/10^3)*gg.ny(nis)','b','linewidth',widths(iQ));
    end
    disp('Done');

%     % channels
%     widfac = .1; widmin = .2; widmax = 2;
%     wide = widfac*ps.Q*abs(Qx); wide = min(max(widmin,wide),widmax); nnz(wide>widmin),
%     for ei = gg.es',
%         wi = wide(ei); if wi>widmin, ni = find(gg.eddx(ei,:)); xx = (ps.x/10^3)*nx(ni)'; yy = (ps.x/10^3)*ny(ni)'; plot(xx,yy,'b','linewidth',wi); end
%     end
%     widf = widfac*ps.Q*abs(Qy); widf = min(max(widmin,widf),widmax); nnz(widf>widmin),
%     for fi = gg.fs',
%         wi = widf(fi); if wi>widmin, ni = find(gg.fddy(fi,:)); xx = (ps.x/10^3)*nx(ni)'; yy = (ps.x/10^3)*ny(ni)'; plot(xx,yy,'b','linewidth',wi); end
%     end
%     wids = widfac*ps.Q*abs(Qs); wids = min(max(widmin,wids),widmax); nnz(wids>widmin),
%     widr = widfac*ps.Q*abs(Qr); widr = min(max(widmin,widr),widmax); nnz(widr>widmin),
%     for ci = gg.cs',
%         wi = wids(ci); if wi>widmin, ni = find(gg.cdds(ci,:)); xx = (ps.x/10^3)*nx(ni)'; yy = (ps.x/10^3)*ny(ni)'; plot(xx,yy,'b','linewidth',wi); end
%         wi = widr(ci); if wi>widmin, ni = find(gg.cddr(ci,:)); xx = (ps.x/10^3)*nx(ni)'; yy = (ps.x/10^3)*ny(ni)'; plot(xx,yy,'b','linewidth',wi); end
%     end
    
    %     axis([11 127 30 129]);
    text(0,1,'Q','Units','normalized','VerticalAlignment','bottom');
    axis image;
    axis(axx);
    
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
    
    box on;
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_discharge_lines']); end
end

if figi==3, %velocity
clf; axes('position',axes_size);

    % velocity
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.u_b*pd.ty)*reshape(aa.Ub,gg.nI,gg.nJ); xx(zz==0) = NaN; xx(isnan(zz)) = NaN;
    cax = [1e1 1e4]; 
    
%     % pcolor version
%     zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
%     hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
%     caxis(cax); 
%     axis image;
%     axis(axx);
%     
%     % color bar
%     tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
%     cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     if ver==1,
%     text(-0.05,0,'U [ m y^{-1} ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
%     elseif ver==2,
%     cx.Label.String = 'U [ m y^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
%     end
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
    if ver==1,
    text(-0.05,0,'U [ m y^{-1} ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'U [ m y^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
%     % add boundary
%     hold on; plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
%     
    % add boundary (alternative method)
    ns = union(intersect([gg.bdy.nlbdy; gg.bdy.nrbdy; gg.bdy.ntopbdy; gg.bdy.nbotbdy],gg.ns),gg.n1); % pts on boundary
    x = gg.nx(ns); y = gg.ny(ns); is = nevis_orderboundary(x,y);
    hold on; plot((ps.x/10^3)*x(is),(ps.x/10^3)*y(is),'k-');
    
    load cmapu1; colormap(cmap);
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_velocity']); end
end

if figi==4, %pressure
clf; axes('position',axes_size);

    % effective pressure
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.phi/10^6)*reshape(phi_0-phi,gg.nI,gg.nJ); 
    cax = [-3 3];
    cax = [-2 2];
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    if ver==1,
        text(-0.05,0,'N [ MPa ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'N [ MPa ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
    % add boundary (alternative method)
    ns = union(intersect([gg.bdy.nlbdy; gg.bdy.nrbdy; gg.bdy.ntopbdy; gg.bdy.nbotbdy],gg.ns),gg.n1); % pts on boundary
    x = gg.nx(ns); y = gg.ny(ns); is = nevis_orderboundary(x,y);
    hold on; plot((ps.x/10^3)*x(is),(ps.x/10^3)*y(is),'k-');
    
    load cmapbluered; colormap(cmap(end:-1:1,:));
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_pressure']); end
end

if figi==5, %sheet
clf; axes('position',axes_size);

    % sheet thickness
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.h*100)*reshape(hs,gg.nI,gg.nJ); 
    cax = [0 50]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    if ver==1,
        text(-0.05,0,'h_s [ cm ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'h_s [ cm ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
    load cmapq2; colormap(cmap);
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_sheet']); end
end

if figi==6, %elastic
clf; axes('position',axes_size);

    % elastic sheet thickness
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.h*100)*reshape(he,gg.nI,gg.nJ); 
    cax = [0 500]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    if ver==1,
    text(-0.05,0,'h_e [ cm ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'h_e [ cm ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
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

    load cmapq2; colormap(cmap);
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_elastic']); end
end

if figi==7, %area
clf; axes('position',axes_size);

    % channel volume, converted to sheet thickness
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.h*100)*reshape(VS./(gg.Dx.*gg.Dy),gg.nI,gg.nJ); 
    cax = [0 .1]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    if ver==1,
    text(-0.05,0,'h_S [ cm ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'h_S [ cm ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
    % add boundary (alternative method)
    ns = union(intersect([gg.bdy.nlbdy; gg.bdy.nrbdy; gg.bdy.ntopbdy; gg.bdy.nbotbdy],gg.ns),gg.n1); % pts on boundary
    x = gg.nx(ns); y = gg.ny(ns); is = nevis_orderboundary(x,y);
    hold on; plot((ps.x/10^3)*x(is),(ps.x/10^3)*y(is),'k-');
    
    load cmapq2; colormap(cmap);
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_area']); end
end

if figi==8, %topography
clf; axes('position',axes_size);

    % bed
    xx = (ps.x/10^3)*gg.nx; yy = (ps.x/10^3)*gg.ny;  % keep outside points here to prevent edge effects disappearing
    zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=-9400) = NaN;
    cax = z_range; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
%     hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); shading interp;
    hand = imagesc(xx(:,1),yy(1,:),zz'); set(gca,'YDir','normal');
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    if ver==1,
    text(-0.05,0,'z [ m ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'z [ m ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
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
    
%     % add outlets
%     plot((ps.x/10^3)*nx(gg.nbdy),(ps.x/10^3)*ny(gg.nbdy),'ko','markersize',4,'markerfacecolor','y');
%     
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
%     
    % add boundary (alternative method)
    ns = union(intersect([gg.bdy.nlbdy; gg.bdy.nrbdy; gg.bdy.ntopbdy; gg.bdy.nbotbdy],gg.ns),gg.n1); % pts on boundary
    x = gg.nx(ns); y = gg.ny(ns); is = nevis_orderboundary(x,y);
    plot((ps.x/10^3)*x(is),(ps.x/10^3)*y(is),'k-');

    load cmapland2; cmap(end,:) = .4*[1 1 1]; colormap(cmap); 
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_topography']); end
end

if figi==9, %dissipation
clf; axes('position',axes_size);

    % dissipation
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.Xi/ps.x)*reshape(Xi,gg.nI,gg.nJ); 
    cax = [10^(-4) 10]; 
    
    % logarithmic pcolor
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis
    hand = pcolor(xx,yy,log10(zz)); set(hand,'linestyle','none'); % shading interp;
    caxis(log10(cax)); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
    if ver==1,
    text(-0.05,0,'\Xi [ J m^{-2} s^{-1} ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = '\Xi [ J m^{-2} s^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
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
    text(0,1,'\Xi','Units','normalized','VerticalAlignment','bottom');
    
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_dissipation']); end
end

if figi==10, %thickness
clf; axes('position',axes_size);
    
    % ice thickness
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = ps.z*reshape(H,gg.nI,gg.nJ); 
    cax = [0 1500]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    if ver==1,
        text(-0.05,0,'H [ m ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'H [ m ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_thickness']); end
end

if figi==11, %input
clf; axes('position',axes_size);

    % surface input
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.m*10^3*pd.td)*reshape(E,gg.nI,gg.nJ); 
    cax = [0 80]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
    if ver==1,
        text(-0.05,0,'m [ mm d^{-1} ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
    elseif ver==2,
    cx.Label.String = 'm [ mm d^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';
    end
    
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
    text(0,1,'\Xi','Units','normalized','VerticalAlignment','bottom');
    
    if reversey, set(gca,'YDir','reverse'); end
    text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
    text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');
    if save_plot, print(gcf,'-depsc2',[fn,'_input']); end
end

end