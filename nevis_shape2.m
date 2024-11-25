function [n1,n2,e1,e2,f1,f2,s1,s2,r1,r2] = nevis_shape2(gg,ns,oo)
% identify boundary of region with nodes ns
% 10 April 2016: nevis_shape2 is same as nevis_shape but included boundary connections
%   across diagonals
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

%% identify boundary nodes and edges
temp = 0*ones(gg.nIJ,1); temp(ns) = 1;
e0 = find(gg.emean*temp==0); % edges connected to all outside nodes
e1 = find(gg.emean*temp==0.5); % edges connected to one outside node : boundary edges
f0 = find(gg.fmean*temp==0); % edges connected to all outside nodes
f1 = find(gg.fmean*temp==0.5); % edges connected to one outside node : boundary edges
s0 = find(gg.cmeans*temp==0); % corner edges connected to all outside nodes
s1 = find(gg.cmeans*temp==0.5); % corner edges connected to one outside node : boundary edges
r0 = find(gg.cmeanr*temp==0); % corner edges connected to all outside nodes
r1 = find(gg.cmeanr*temp==0.5); % corner edges connected to one outside node : boundary edges

temp = zeros(gg.eIJ,1); temp(e1) = 1; n1x = find(gg.nmeanx*temp>0);  
temp = zeros(gg.fIJ,1); temp(f1) = 1; n1y = find(gg.nmeany*temp>0); 
temp = zeros(gg.cIJ,1); temp(s1) = 1; n1s = find(gg.nmeans*temp>0);  
temp = zeros(gg.cIJ,1); temp(r1) = 1; n1r = find(gg.nmeanr*temp>0);  
n1 = intersect([n1x;n1y;n1s;n1r],ns); % nodes connected by an edge to an outside node : boundary nodes
n2 = setdiff([n1x;n1y;n1s;n1r],ns); % nodes connected by an edge to an inside node : boundary nodes

temp = zeros(gg.nIJ,1); temp(n1) = 1;
e2 = find(gg.emean*temp==1);  % edges connected to two edge nodes : along-boundary edges
f2 = find(gg.fmean*temp==1);  % edges connected to two edge nodes : along-boundary edges
s2 = find(gg.cmeans*temp==1);  % corner edges connected to two edge nodes : along-boundary edges
r2 = find(gg.cmeanr*temp==1);  % corner edges connected to two edge nodes : along-boundary edges


end
