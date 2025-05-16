% 22 Sept 2014: test routing algorithm

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
addpath(oo.code);       % add path to code

%% parameters
% [ put non-default parameters and options here ]
pd = struct;
[pd,oo] = nevis_defaults(pd,oo);
[ps,pp] = nevis_nondimension(pd);

% %% grid and geometry
% x = linspace(0,(10000/ps.x),50); y = linspace(0,(10000/ps.x),50);
% gg = nevis_grid(x,y,oo);
% b = (0/ps.z)*gg.nx;
% s = (1000/ps.z)*2*max(0.5^2-(gg.nx-0.5).^2-(gg.ny-0.5).^2,0).^(1/2);
% 
% %% mask grid
% gg = nevis_mask(gg,find(s-b<=0)); 
% gg.n1m = gg.n1;                 % margin boundary nodes
% gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes
% 
% %% initialize
% [aa,vv] = nevis_initialize(b,s,gg,pp,oo);
% vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a); 
% vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);
% 
% [Q,low,Qx,Qy,Qs,Qr,nconnect4] = nevis_route(aa.phi_a+(1000/910)*(aa.phi_0-aa.phi_a),ones(gg.nIJ,1),gg,oo);
% 
% %% plot
% figure(1); clf;
%     %upstream area
%     xx = (ps.x/10^3)*gg.nx; yy = (ps.x/10^3)*gg.ny;  
%     zz = ps.x.^2*reshape(Q,gg.nI,gg.nJ); zz(gg.nout) = NaN;
% %     cax = [0 2*mean(zz(~isnan(zz)))];
%     cax = [0 1e6];
%     imagesc(xx(:,1),yy(1,:),zz');
%     xlabel('x');
%     ylabel('y');
%     colorbar; caxis(cax);
% 
% return;

% load('../paakitsoq/topo2_100'); dd = tt;
% 
% %% grid and geometry
% gg = nevis_grid(length(dd.x),length(dd.y),1,dd.x(1)/ps.x,dd.x(end)/ps.x,dd.y(1)/ps.x,dd.y(end)/ps.x,oo);
% b = reshape(dd.b/ps.z,gg.nIJ,1);
% s = reshape(dd.s/ps.z,gg.nIJ,1);
% %% adjust bed near margin to have minimum thickness [ trust surface more than bed; surface is NaN for icefree areas ]
% Hmin = (1+1e-3)*15/ps.z; tmp = (b>s-Hmin); b(tmp) = s(tmp)-Hmin; 
% 
% %% mask with minimum thickness
% H = max(s-b,0);
% Hmin = 15/ps.z; 
% gg = nevis_mask(gg,find(H<Hmin)); 
% gg.n1m = gg.n1;
% 
% %% re-mask to area inside mask
% H(~dd.mask) = 0;
% nout = find(H<Hmin);
% gg = nevis_mask(gg,nout);
% gg.n1m = intersect(gg.n1m,gg.n1); % pts on boundary that are there because of minimum thickness rather than mask
% 
% %% label boundary nodes
% gg = nevis_label(gg,gg.n1m);
% oo.adjust_boundaries = 1;
% 
% %% initialize
% [aa,vv] = nevis_initialize(b,s,gg,pp,oo);
% vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a); 
% vv.hs = 0.1/ps.hs*ones(gg.nIJ,1);
% 
% %% calculate upstream areas
% [Q,low,~,~,~,~,nconnect4,nconnect5,drains] = nevis_route(aa.phi_a+(1000/910)*(aa.phi_0-aa.phi_a),ones(gg.nIJ,1),gg,oo);
% 
% %% plot
% figure(1); clf;
%     %upstream area
%     xx = (ps.x/10^3)*gg.nx; yy = (ps.x/10^3)*gg.ny;  
%     zz = ps.x.^2*reshape(Q,gg.nI,gg.nJ); zz(gg.nout) = NaN;
% %     cax = [0 2*mean(zz(~isnan(zz)))];
%     cax = [0 1e6];
%     imagesc(xx(:,1),yy(1,:),zz');
%     xlabel('x');
%     ylabel('y');
%     colorbar; caxis(cax);
%     load cmapq2; colormap(cmap);
%     
%     hold on;
%     % add lows
% 	ni = find(low);
%     plot(xx(ni),yy(ni),'ko','Markersize',8,'MarkerFaceColor',0.8*[1 1 1]); 
%  
%     %boundary curve
%     x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
%     tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary    
%     % add boundary
%     plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 
%     
%     %plot catchment boundaries
%     lows = find(low==1);
%     for il = 1:length(lows),
% %         x = gg.nx(drains(il,:)); y = gg.ny(drains(il,:));
% %         x = [x(1)-1e-3 x(1)+1e-3 x]; y = [y(1)-1e-3 y(1)-1e-3 y]; % add close point so as not to have too few points
% %         k = convhull(x,y);
% %         plot((ps.x/10^3)*x(k),(ps.x/10^3)*y(k),'k-');
%         bw = zeros(gg.nI,gg.nJ);
%         bw(drains(il,:)) = 1;
%         B = bwboundaries(bw);
%         tmp = B{1}(:,1)+gg.nI*(B{1}(:,2)-1);
%         x = gg.nx(tmp); y = gg.ny(tmp);
%         plot((ps.x/10^3)*x,(ps.x/10^3)*y,'r-');
%     end
%     
% 
% return;  
    

load('../gap/topo_140901'); % topography input, cropped

%% grid and geometry
gg = nevis_grid(dd.x(dd.xi_crop)/ps.x,dd.y(dd.yi_crop)/ps.x,oo); ps.x0 = dd.x0; ps.y0 = dd.y0;
b = reshape(dd.b(dd.xi_crop,dd.yi_crop)/ps.z,gg.nIJ,1);
s = reshape(dd.s(dd.xi_crop,dd.yi_crop)/ps.z,gg.nIJ,1);
h = fspecial('gaussian',3,1); s_smooth = filter2(h,reshape(s,gg.nI,gg.nJ)); % smooth surface for surface routing

%% mask with minimum thickness and mask
b(~dd.mask(dd.xi_crop,dd.yi_crop)) = NaN; s(~dd.mask(dd.xi_crop,dd.yi_crop)) = NaN; 
H = max(s-b,0);
H(isnan(H)) = 0;
Hmin = 15/ps.z; 
gg = nevis_mask(gg,find(H<Hmin)); 
gg.n1m = gg.n1(H(gg.n1)<100/ps.z); 

%% label boundary nodes
gg = nevis_label(gg,gg.n1m); oo.adjust_boundaries = 1;

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a); 
vv.hs = 0.1/ps.hs*ones(gg.nIJ,1);

%% calculate upstream areas
% [Q,low,~,~,~,~,nconnect4,nconnect5,drains] = nevis_route(aa.phi_a+(1000/910)*(aa.phi_0-aa.phi_a),ones(gg.nIJ,1),gg,oo);
[Q,low,~,~,~,~,nconnect4,nconnect5,drains] = nevis_route(reshape(s_smooth,gg.nIJ,1),ones(gg.nIJ,1),gg,oo);

%% plot
figure(1); clf;
    %upstream area
    xx = (ps.x/10^3)*gg.nx; yy = (ps.x/10^3)*gg.ny;  
    zz = ps.x.^2*reshape(Q,gg.nI,gg.nJ); zz(gg.nout) = NaN;
%     cax = [0 2*mean(zz(~isnan(zz)))];
    cax = [0 1e6];
    imagesc(xx(:,1),yy(1,:),zz');
    xlabel('x');
    ylabel('y');
    colorbar; caxis(cax);
    load cmapq2; colormap(cmap);
    
    hold on;
    % add lows
	ni = find(low);
    plot(xx(ni),yy(ni),'ko','Markersize',8,'MarkerFaceColor',0.8*[1 1 1]); 
 
    %boundary curve
    x_out = gg.nx(gg.n1); y_out = gg.ny(gg.n1);
    tmp = nevis_orderboundary(x_out,y_out); x_out = x_out(tmp); y_out = y_out(tmp); % reorder to follow boundary    
    % add boundary
    plot((ps.x/10^3)*x_out,(ps.x/10^3)*y_out,'k','linewidth',1.5); 

