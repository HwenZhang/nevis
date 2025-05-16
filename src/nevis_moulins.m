function [ni_m,sum_m] = nevis_moulins(x_m,y_m,gg,oo)
% find indices and locations of moulins nearest coordinates [x_m y_m]
%   also calculates matrix sum_m [n_m-by-nIJ] to describe Voronoi cell for each moulin [each row
%   corresponds to a moulin, having entries 1 for nodes in its cell, 0.5 for nodes on boundary, and 0
%   for nodes outside] [eg area_m = sum_m*(Dx.*Dy) gives the area of each moulin's cell]
% inputs  
%   x_m,y_m   coordinates of moulins  
%   gg     	  grid structre containing node coordinates nx and ny
%   oo        [optional] options
% outputs 
%   n_m       number of moulins [ie. length(x_m)]
%   x_m,y_m   [n_m-by-1] coordinates of moulin grid points
%   ni_m      [n_m-by-1] indices of moulin nodes
%   sum_m     [n_m-by-nIJ] matrix to sum quantities defined on nodes over each moulin's grid cell
%
% 21 August 2014 : taken from hydro_locate_m2

% determine the type of moulin prescribed
if nargin<4, oo = struct; end
if ~isfield(oo,'random_moulins'), oo.random_moulins = 0; end % pick random moulin coordinates [value should be desired number of moulins]
if ~isfield(oo,'keep_all_moulins'), oo.keep_all_moulins = 0; end % keep all moulins even if on same node
if ~isfield(oo,'move_moulins'), oo.move_moulins = 0; end % move moulins that are outside domain to nearest position inside domain
if ~isfield(oo,'prescribed_catchments'), oo.prescribed_catchments = 0; end % use prescribed shapefile S to give catchments of moulins

%% pick random moulin coordinates
if oo.random_moulins
    xl = min(min(gg.nx)); xr = max(max(gg.nx)); 
    yb = min(min(gg.ny)); yt = max(max(gg.ny)); 
    n_m = oo.random_moulins;          % number of moulins
    x_m = xl + (xr-xl)*rand(n_m,1);   % x locations
    y_m = yb + (yt-yb)*rand(n_m,1);   % y locations
%     % or
%     x_m = xl + (xr-xl)*rand(n_m,1).^(1/2);   %x locations
%     y_m = yb + (yt-yb)*rand(n_m,1);   %y locations
%     % or 
%     x_m = xl + (xr-xl)*rand(n_m,1).^(1/2);   %x locations
%     f = @(x) 0.0001*x+x.^5;         %x = 0:0.001:1; plot(x,0.5+0.5*f(x-0.5)/f(0.5));
%     y_m = yb + (yt-yb)*( 0.5+0.5*f(rand(n_m,1)-0.5)/f(0.5) );   %dimensionless y locations
end

%% find indices close to moulins
nx = gg.nx; 
ny = gg.ny;
if oo.move_moulins % restrict to active nodes
    nx = nx(gg.ns); 
    ny = ny(gg.ns); 
end 
n_m = length(x_m); % amount of moulins 
ni_m = zeros(n_m,1); % indices close to each moulin
for i_m = 1:n_m
% ni_m(i_m) = find(nx>=x_m(i_m) & ny>=y_m(i_m),1,'first');
% find the closest index 
    [~,tmp] = min((reshape(nx,[],1)-x_m(i_m)).^2+...
        (reshape(ny,[],1)-y_m(i_m)).^2); 
    if oo.move_moulins
        ni_m(i_m) = gg.ns(tmp); % don't understand this
    else 
        ni_m(i_m) = tmp; 
    end
end
% if ~oo.keep_all_moulins, [ni_m,~,tmp] = unique(ni_m); end; %ni_m = ni_m(tmp); end % [ 26 August 2014: this doesn't look to be actually getting rid of repeated moulins, line below would be better? ]
if ~oo.keep_all_moulins
    ni_m = unique(ni_m,'stable'); 
end % [ 26 August 2014: changed from above line, option 'stable' doesn't work on older versions of matlab ]
n_m = length(ni_m);
x_m = gg.nx(ni_m); 
y_m = gg.ny(ni_m);

%% set up catchments for each moulin
sum_m = zeros(length(ni_m),gg.nIJ);

% x,y node coord 
x = reshape(gg.nx,gg.nIJ,1); 
y = reshape(gg.ny,gg.nIJ,1);

% coordinates of corners that are far away [to prevent any voronoi cells that matter being unbounded]

% an extension to the computational domain
big = 0.5*(max(x)-min(x) + max(y)-min(y));
big = 0.2*(max(x)-min(x) + max(y)-min(y)); % change 7 Nov
big = 0.05*max(max(x)-min(x),max(y)-min(y)); % change 6 Feb 2016

% corners of the voronoi grid
x_far = [max(x)+big; max(x)+big; min(x)-big; min(x)-big];
y_far = [max(y)+big; min(y)-big; max(y)+big; min(y)-big];

% calculate voronoi cells [ last four correspond to far away corners and are irrelevant ]
[v,c] = voronoin([x_m y_m;x_far y_far]);

% find nodes within each moulin cell
for i = 1:length(ni_m)
    [in] = inpolygon(x,y,v(c{i},1),v(c{i},2));
    sum_m(i,:) = in';
end

%% prescribed catchments
if oo.prescribed_catchments
    sum_m = zeros(length(ni_m),gg.nIJ);
    for i=1:length(ni_m)
        S = oo.catchments;
        xx = S(i).x(1:end-1);
        yy = S(i).y(1:end-1);
        [in] = inpolygon(x,y,xx(~isnan(xx)),yy(~isnan(yy)));
        sum_m(i,:) = in';
    end
end

% split contibution of nodes that are on the boundary of cells
temp = sum(sum_m,1); % number of polygons each node is connected to (0-3) 
temp(temp==0) = inf;
sum_m = sum_m*sparse(1:length(temp),1:length(temp),temp.^(-1)); %

%% plot voronoi cells and nodes
% figure;
%     voronoi([x_m; x_far],[y_m; y_far]); % voronoi cell
%     hold on
%     plot(x_m,y_m,'ko',x_far,y_far,'ko','MarkerSize',16);
%     % distniguish between different types of nodes
%     sz = 12;
%     plot(x(temp==inf),y(temp==inf),'g.','markersize',sz); % external
%     % nodes feeding 1 moulin
%     plot(x(temp==1),y(temp==1),'r.','markersize',sz); 
%     % nodes feeding 2 moulins
%     plot(x(temp==2),y(temp==2),'b.','markersize',sz);
%     % nodes feeding 3 moulins
%     plot(x(temp==3),y(temp==3),'k.','markersize',sz); 
%     % axis range
%     % axis([min(x) max(x) min(y) max(y)]); 
% shg; 

end