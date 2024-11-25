function nevis_plot_topo(dd,oo)
% plot topography in struct dd

if nargin<2. oo = struct; end
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 0; end
if isfield(oo,'z_range'), z_range = oo.z_range; else z_range = [-999 1000]; end
if isfield(oo,'z_conts'), z_conts = oo.z_conts; else z_conts = -1000:50:5000; end
sc_x = 1000; % plot distances in km

s = dd.s; 
b = dd.b; 
x = (1/sc_x)*dd.x; 
y = (1/sc_x)*dd.y; 
[yy,xx] = meshgrid(y,x);
% s = s.*(s>b)+b.*(s<=b);  % ice depth derived from s and b cleared up
H = s-b;
  
%% plot 
figure; clf; set(gcf,'Units','centimeters','Position',[5 2 20 20]);  

    %% image map
    imagesc(x,y,H',z_range); set(gca,'YDir','normal'); %title('Ice depth');
%     imagesc(x,y,s',z_range); set(gca,'YDir','normal'); %title('Surface');
%     imagesc(x,y,b',z_range); set(gca,'YDir','normal'); %title('Bed'); 
    
    load cmapland2; colormap(cmap); caxis(z_range); colorbar; 
    xlabel('x [ km ]');
    ylabel('y [ km ]');
    axis image; 
    
    %% bed contours
    fil = fspecial('average'); z = conv2(b,fil,'valid'); 
    xrg = (size(b,1)-size(z,1))/2+(1:size(z,1)); yrg = (size(b,2)-size(z,2))/2+(1:size(z,2));
    hold on; contour(xx(xrg,yrg),yy(xrg,yrg),z,z_conts,'color',0.5*[1 1 1]);
  
    %% surface contours
    fil = fspecial('average'); z = conv2(s,fil,'valid');  % fil = fspecial('average',[9 9]);
    xrg = (size(s,1)-size(z,1))/2+(1:size(z,1)); yrg = (size(s,2)-size(z,2))/2+(1:size(z,2));
    hold on; contour(xx(xrg,yrg),yy(xrg,yrg),z,z_conts,'k','linewidth',1);
    
    %% moulins
    if isfield(dd,'x_m'),
        hold on; plot((1/sc_x)*dd.x_m,(1/sc_x)*dd.y_m,'ko','Markersize',4,'MarkerFaceColor','w'); 
    end
    
    %% outlets
    if isfield(dd,'x_d'),
        plot((1/sc_x)*dd.x_d,(1/sc_x)*dd.y_d,'ko','markerfacecolor','y'); 
    end
    
    if reversey, set(gca,'YDir','reverse'); end
end