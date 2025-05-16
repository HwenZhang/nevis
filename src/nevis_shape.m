function [n1,n2,e1,e2,f1,f2] = nevis_shape(gg,ns,oo)
% identify boundary of region with nodes ns
% 
% inputs
%   gg grid structure [see below for contents]
%   ns 
%   oo [optional] option structure with method to use [see below for options]
% outputs:
%   n1  boundary nodes [inside]
%   n2  boundary nodes [outside]
%   e1  boundary edges
%   f1  boundary edges
%   e2  along boundary edges
%   f2  along boundary edges
% 
% 10 Dec 2014

if nargin<4, oo = struct; end
if any(ns==0), ns = find(ns); end % interpret ns as logical if any entries are zero

%% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*I ]
na = (1:gg.nIJ)';
ea = (1:gg.eIJ)';
fa = (1:gg.fIJ)';
ca = (1:gg.cIJ)';

%% identify boundary nodes and edges
temp = 0*ones(gg.nIJ,1); temp(ns) = 1;
e0 = find(gg.emean*temp==0); % edges connected to all outside nodes
e1 = find(gg.emean*temp==0.5); % edges connected to one outside node : boundary edges
f0 = find(gg.fmean*temp==0); % edges connected to all outside nodes
f1 = find(gg.fmean*temp==0.5); % edges connected to one outside node : boundary edges
c0 = find(gg.cmean*temp==0); % corners connected to all outside nodes
c1 = find(gg.cmean*temp==0.75); % corners connected to one outside node 
c2 = find(gg.cmean*temp==0.5); % corners connected to two outside nodes
c3 = find(gg.cmean*temp==0.25); % corners connected to three outside nodes

temp = zeros(gg.eIJ,1); temp(e1) = 1; n1x = find(gg.nmeanx*temp>0);  
temp = zeros(gg.fIJ,1); temp(f1) = 1; n1y = find(gg.nmeany*temp>0); 
n1 = intersect([n1x;n1y],ns); % nodes connected by an edge to an outside node : boundary nodes
n2 = setdiff([n1x;n1y],ns); % nodes connected by an edge to an inside node : boundary nodes

temp = zeros(gg.nIJ,1); temp(n1) = 1;
e2 = find(gg.emean*temp==1);  % edges connected to two edge nodes : along-boundary edges
f2 = find(gg.fmean*temp==1);  % edges connected to two edge nodes : along-boundary edges

end
