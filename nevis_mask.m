function [gg] = nevis_mask(gg,nout,oo)
% mask to identify boundaries based on outside node indices nout
% combine with label to label boundary nodes and edges
% 
% inputs
%   gg grid structure [see below for contents]
%   nout indices of exterior nodes
%   oo [optional] option structure with method to use [see below for options]
% outputs:
%   gg grid structure [added/modified labels]
% 
% 31 July 2014 : taken from hydro_mask2, but alterered a lot

if nargin<3, oo = struct; end

%% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*I ]
ns = (1:gg.nIJ)';
es = (1:gg.eIJ)';
fs = (1:gg.fIJ)';
cs = (1:gg.cIJ)';

%% identify boundary nodes and edges
temp = ones(gg.nIJ,1); temp(nout) = 0;
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
n1 = setdiff([n1x;n1y],nout); % nodes connected by an edge to an outside node : boundary nodes

temp = zeros(gg.nIJ,1); temp(n1) = 1;
e2 = find(gg.emean*temp==1);  % edges connected to two edge nodes : along-boundary edges
f2 = find(gg.fmean*temp==1);  % edges connected to two edge nodes : along-boundary edges

% [ should perhaps have different labels for diagonal channels, so can have ones in one direction and
% not the other ]

%% prescribe flux on all boundaries as default
nbdy = [];
ebdy = [e1]; eout = [e0];
fbdy = [f1]; fout = [f0];
cbdy = []; cout = unique([c0;c1;c2;c3]);

% nbdy = setdiff(nbdy,nout);   
% ebdy = setdiff(ebdy,eout);   
% fbdy = setdiff(fbdy,fout);   
% cbdy = setdiff(cbdy,cout); 

%% define active and interior labels
ns = setdiff(ns,nout);   
es = setdiff(es,eout);   
fs = setdiff(fs,fout);   
cs = setdiff(cs,cout);  

nin = setdiff(ns,nbdy);   
ein = setdiff(es,ebdy);  
fin = setdiff(fs,fbdy);  
cin = setdiff(cs,cbdy);  

gg.ns = ns;
gg.es = es;
gg.fs = fs;
gg.cs = cs;
gg.nbdy = nbdy;
gg.ebdy = ebdy;
gg.cbdy = cbdy;
gg.fbdy = fbdy;
gg.nin = nin;
gg.ein = ein;
gg.fin = fin;
gg.cin = cin;
gg.nout = nout;
gg.eout = eout;
gg.fout = fout;
gg.cout = cout;
gg.n1 = n1;
gg.e0 = e0;
gg.e1 = e1;
gg.e2 = e2;
gg.f0 = f0;
gg.f1 = f1;
gg.f2 = f2;
gg.c0 = c0;
gg.c1 = c1;
gg.c2 = c2;
gg.c3 = c3;

% %% boundary labels [ may not really be needed ? ]
% gg.bb.ns = n1;
% gg.bb.es = e1;
% gg.bb.fs = f1;
% gg.bb.es2 = e2;
% gg.bb.fs2 = f2;
% gg.bb.cs = unique([c1; c2; c3]);

end
