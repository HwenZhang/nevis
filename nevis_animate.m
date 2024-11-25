function nevis_animate(fn,tis,mode,make_movie)
% animates solution in directory fn
%
% 21 August 2014: taken from animate_arolla

%% mode
if nargin<3, mode = 1; end %
% 1 discharge
% 2 pressure
% 3 elastic sheet
% 4 cavity sheet
% 5 area

% could copy others from nevis_plot

if nargin<4
    make_movie = 0; 
end
% 0 don't make movie
% 1 make movie using avifile
% 2 make movie using VideoWriter and getframe
% 3 make figures

%% options
fps = 4;   

%% load initial timestep
if exist([fn,'/0000']), load([fn,'/0000']); else load(fn); end
load([fn,'/0001']);
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end

%% open movie file
if mode == 1, fnm = 'movie_q'; 
elseif mode == 2, fnm = 'movie_N';
end
if make_movie == 1,
    mov = avifile([fn,'_',fnm],'fps',fps); 
elseif make_movie == 2,
    obj = VideoWriter([fn,'_',fnm]); 
    obj.FrameRate = fps;
    open(obj);
end

%% contour levels
z_conts = -1000:50:5000;
% z_conts = -1000:10:5000;
z_conts = -3000:100:3000;
b_conts = z_conts;
s_conts = z_conts;
p_conts = z_conts/100;

%% axes limits
axx = (ps.x/10^3)*[gg.xl gg.xr gg.yb gg.yt];

%% extract initial variables
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,gg);

%get rid of points outside domain
nx(gg.nout) = NaN;
ex(gg.eout) = NaN;
fx(gg.fout) = NaN;
cx(gg.cout) = NaN;
    
%boundary curve
x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
if ~isempty(x_out), tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); end % reorder to follow boundary

%% set up figure
figure(1); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 15 15]);
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

%% plot frames
for i_t = 1:length(tis)
disp(['Frame ',num2str(i_t),' / ',num2str(length(tis)),' ...']);

%% load timestep
load([fn,'/',int2four(tis(i_t))]);

%% extract new variables
if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,vv2);
clear tt vv vv2

%% plot
clf; axes('position',axes_size);

% [ could copy other things here from nevis_plot ]
if mode == 1
% %%
    % discharge
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
    tmp = get(cx,'XTick'); labels = 10.^tmp; set(cx,'XTickLabel',labels);
%     text(-0.05,0,'|q| [ m^2 s^{-1} ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
%     cx.Label.String = 'q [ m^2 s^{-1} ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';

    load cmapq2; colormap(cmap);
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
elseif mode == 2
%%
    % effective pressure
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.phi/10^6)*reshape(phi_0-phi,gg.nI,gg.nJ); 
    cax = [-1 1]; % cax = [-.5 .5]; % [-3 3]; 
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     text(-0.05,0,'N [ MPa ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');
%     cx.Label.String = 'N [ MPa ]'; cx.Label.Units = 'normalized'; cx.Label.Position = [-0.05 0]; cx.Label.VerticalAlignment = 'bottom'; cx.Label.HorizontalAlignment = 'right';

    load cmapbluered; colormap(cmap(end:-1:1,:));
elseif mode == 3
%%
    % elastic sheet
    xx = (ps.x/10^3)*nx; yy = (ps.x/10^3)*ny;  
    zz = (ps.h*100)*reshape(he,gg.nI,gg.nJ); 
    cax = [0 50]; %[0 100]
    
    zz = max(cax(1),min(zz,cax(2))); % cap so as to avoid going off color axis     
    hand = pcolor(xx,yy,zz); set(hand,'linestyle','none'); % shading interp;
    caxis(cax); 
    axis image;
    axis(axx);
    
    % color bar
    tmp = get(gca,'position'); cbar_size = [tmp(3)*3/16 0.01  tmp(3)*5/8 0.02];
    cx = colorbar('peer',gca,'horizontal','position',[tmp(1) tmp(2)+tmp(4) 0 0]+cbar_size,'XAxislocation','top');
%     text(-0.05,0,'h_e [ cm ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');

    load cmapq2; colormap(cmap);
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
elseif mode == 4
%%
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
%     text(-0.05,0,'h_s [ cm ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');

    load cmapq2; colormap(cmap);
    
    % add pressure contours
    hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
   
elseif mode == 5    
%%
    % area
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
%     text(-0.05,0,'h_S [ cm ]','units','normalized','Parent',cx,'VerticalAlignment','bottom','HorizontalAlignment','right');

    load cmapq2; colormap(cmap);
    
    % add pressure contours
%     hold on; hand = contour((ps.x/10^3)*nx,(ps.x/10^3)*ny,(ps.phi/10^6)*reshape(phi,gg.nI,gg.nJ),p_conts,'linecolor','r','linewidth',0.5);
    
end

%     % add bed contours
%     zz = ps.z*reshape(b,gg.nI,gg.nJ); zz(zz<=-400) = NaN;
%     fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
%     xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
%     hold on; contour(xx(xrg,yrg),yy(xrg,yrg),zzz,b_conts,'color',0.5*[1 1 1]);
% 
%     % add surface contours
%     zz = ps.z*reshape(s,gg.nI,gg.nJ);
%     fil = fspecial('average'); zzz = conv2(zz,fil,'valid'); %zzz = zz;
%     xrg = (size(zz,1)-size(zzz,1))/2+(1:size(zzz,1)); yrg = (size(zz,2)-size(zzz,2))/2+(1:size(zzz,2));
%     hold on; contour(xx(xrg,yrg),yy(xrg,yrg),zzz,b_conts,'color',0.5*[1 1 1],'linewidth',1);
% %     contour(xx,yy,ps.z*reshape(s-b,gg.nI,gg.nJ),1*[1 1],'color','k','linewidth',1); % mark margin

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
    
    %% formatting   
    if reversey, set(gca,'YDir','reverse'); end
%     text(-0.08,0.5,'y [ km ]','HorizontalAlignment','center','units','normalized','rotation',90);
%     text(0.5,-0.08/aspect,'x [ km ]','HorizontalAlignment','center','VerticalAlignment','bottom','units','normalized');

    text(0.05,0.95,['t = ',num2str(round((ps.t/(24*60*60))*t*100)/100)],'units','normalized');
    set(gcf,'Color','w');
    drawnow; shg;
    

if make_movie == 1,
%% add frame to movie
disp('Saving ...');
pause(0.05); mov = addframe(mov,gcf); pause(0.05);
elseif make_movie == 2
disp('Saving ...');
pause(0.05); set(gcf,'units','pixels'); rect = get(gcf,'Position'); frame = getframe(gcf,[0 0 rect(3) rect(4)]); writeVideo(obj,frame); pause(0.05);
elseif make_movie == 3,
%% save figure
disp('Saving ...');
print(gcf,'-depsc2',[fn,'/',fnm,'_',num2str(tis(i_t))]);
end
    
end

%% close movie
if make_movie ==1, 
mov = close(mov);
elseif make_movie == 2
close(obj);
end

disp('Done');
end