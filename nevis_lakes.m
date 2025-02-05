function [ni_l,sum_l] = nevis_lakes(x_l,y_l,gg,oo)
% 27/11/2024 By Hanwen Zhang, modified from nevis_moulins.m

% find indices and locations of lake's nearest coordinates [x_m y_m]
% also calculates matrix sum_m [n_m-by-nIJ] to describe Voronoi cell for each lake [each row
% corresponds to a lake, having entries 1 for nodes in its cell, 0.5 for nodes on boundary, and 0
% for nodes outside] [eg area_m = sum_m*(Dx.*Dy) gives the area of each lake's cell]
% inputs  
%   x_m,y_m   coordinates of lakes 
%   gg     	  grid structre containing node coordinates nx and ny
%   oo        [optional] options
% outputs 
%   n_m       number of lakes [ie. length(x_m)]
%   x_m,y_m   [n_m-by-1] coordinates of lake grid points
%   ni_m      [n_m-by-1] indices of lake nodes
%   sum_m     [n_m-by-nIJ] matrix to sum quantities defined on nodes over each lake's grid cell
%
% 21 August 2014 : taken from hydro_locate_m2

% determine the type of lake prescribed
if nargin<4, oo = struct; end
% pick random lake coordinates [value should be desired number of lakes]
if ~isfield(oo,'random_lakes'), oo.random_lakes = 0; end 
% keep all lakes even if on same node
if ~isfield(oo,'keep_all_lakes'), oo.keep_all_lakes = 0; end 
% move lakes that are outside domain to nearest position inside domain
if ~isfield(oo,'move_lakes'), oo.move_lakes = 0; end 
% use prescribed shapefile S to give catchments of lakes
if ~isfield(oo,'prescribed_lake_catchments'), oo.prescribed_lake_catchments = 0; end 
% use lake coordinates from a specific data file
if ~isfield(oo,'lake_output_file'), oo.lake_output_file = 0; end 

%% pick random lake coordinates
if oo.random_lakes
    xl = min(min(gg.nx)); xr = max(max(gg.nx)); 
    yb = min(min(gg.ny)); yt = max(max(gg.ny)); 
    n_l = oo.random_lakes;          % number of lakes
    x_l = xl + (xr-xl)*rand(n_l,1);   % x locations
    y_l = yb + (yt-yb)*rand(n_l,1);   % y locations
end

%% Read in lake coordinates
% if oo.lake_output_file
%     path = '/Users/eart0487/nevis/nevis_1d_example/';
%     [xl,yl] = load(path, 'lake coordinates');
% end

%% find indices close to lakes
nx = gg.nx; 
ny = gg.ny;
if oo.move_lakes % restrict to active nodes
    nx = nx(gg.ns); 
    ny = ny(gg.ns); 
end 
n_l = length(x_l); % amount of lakes 
ni_l = zeros(n_l,1); % indices close to each lake
for i_l = 1:n_l
% find the closest index 
    [~,tmp] = min((reshape(nx,[],1)-x_l(i_l)).^2+...
        (reshape(ny,[],1)-y_l(i_l)).^2); 
    if oo.move_lakes
        ni_l(i_l) = gg.ns(tmp); % don't understand this
    else 
        ni_l(i_l) = tmp; 
    end
end

if ~oo.keep_all_lakes
    ni_l = unique(ni_l,'stable'); 
end

n_l = length(ni_l);
x_l = gg.nx(ni_l); 
y_l = gg.ny(ni_l);

%% set up catchments for each lake
sum_l = zeros(length(ni_l),gg.nIJ);

% x,y node coord 
x = reshape(gg.nx,gg.nIJ,1); 
y = reshape(gg.ny,gg.nIJ,1);

% coordinates of corners that are far away [to prevent any voronoi cells that matter being unbounded]

% an extension to the computational domain
big = 0.5*(max(x)-min(x) + max(y)-min(y));
% big = 0.2*(max(x)-min(x) + max(y)-min(y)); % change 7 Nov
% big = 0.05*max(max(x)-min(x),max(y)-min(y)); % change 6 Feb 2016

% corners of the voronoi grid
x_far = [max(x)+big; max(x)+big; min(x)-big; min(x)-big];
y_far = [max(y)+big; min(y)-big; max(y)+big; min(y)-big];

% calculate voronoi cells [ last four correspond to far away corners and are irrelevant ]
[v,c] = voronoin([x_l y_l;x_far y_far]);

% find nodes within each lake cell
for i = 1:length(ni_l)
    [in] = inpolygon(x,y,v(c{i},1),v(c{i},2));
    sum_l(i,:) = in';
end

%% prescribed catchments
if oo.prescribed_lake_catchments
    sum_l = zeros(length(ni_l),gg.nIJ);
    for i=1:length(ni_l)
        S = oo.lake_catchments;
        xx = S(i).x(1:end-1);
        yy = S(i).y(1:end-1);
        [in] = inpolygon(x,y,xx(~isnan(xx)),yy(~isnan(yy)));
        sum_l(i,:) = in';
    end
end

% split contibution of nodes that are on the boundary of cells
temp = sum(sum_l,1); % number of polygons each node is connected to (0-3) 
temp(temp==0) = inf;
sum_l = sum_l*sparse(1:length(temp),1:length(temp),temp.^(-1)); %

%% plot voronoi cells and nodes
% figure;
%     voronoi([x_m; x_far],[y_m; y_far]); % voronoi cell
%     hold on
%     plot(x_m,y_m,'ko',x_far,y_far,'ko','MarkerSize',16);
%     % distniguish between different types of nodes
%     sz = 12;
%     plot(x(temp==inf),y(temp==inf),'g.','markersize',sz); % external
%     % nodes feeding 1 lake
%     plot(x(temp==1),y(temp==1),'r.','markersize',sz); 
%     % nodes feeding 2 lakes
%     plot(x(temp==2),y(temp==2),'b.','markersize',sz);
%     % nodes feeding 3 lakes
%     plot(x(temp==3),y(temp==3),'k.','markersize',sz); 
%     % axis range
%     % axis([min(x) max(x) min(y) max(y)]); 
% shg; 

end