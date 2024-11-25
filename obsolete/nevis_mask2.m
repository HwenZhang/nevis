function [gg] = nevis_mask2(gg,d_inds,oo)
% mask2 makes the boundary labels and list of boundary indices
% 
% inputs
%   gg grid structure [see below for contents]
%   d_inds indices of nodes at which pressure prescribed
%   oo [optional] option structure with method to use [see below for options]
% outputs:
%   gg grid structure [added/modified labels]
% 
% 31 July 2014 : taken from hydro_mask2, but alterered a lot, split into
%   mask and mask2

if nargin<4, oo = struct; end

%% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*I ]
ns = (1:gg.nIJ)';
es = (1:gg.eIJ)';
fs = (1:gg.fIJ)';
cs = (1:gg.cIJ)';

%% labels of boundary nodes, edges, corners
nout = gg.nout;
n1 = gg.n1;
e0 = gg.e0;
e1 = gg.e1;
e2 = gg.e2;
f0 = gg.f0;
f1 = gg.f1;
f2 = gg.f2;
c0 = gg.c0;
c1 = gg.c1;
c2 = gg.c2;
c3 = gg.c3;

%% prescribe pressure at nodes specified by d_inds [ eliminate connected boundary edges ]
nbdy = intersect(n1,d_inds);    % only include d_inds that are on the boundary

ed1 = find(gg.emean*sparse(nbdy,1,1,gg.nIJ,1)==0.5); ed1 = intersect(ed1,e1); % connected boundary edges
ed2 = find(gg.emean*sparse(nbdy,1,1,gg.nIJ,1)==1); ed2 = intersect(ed2,e2); % connected along-boundary edges
ed = union(ed1,ed2);
fd1 = find(gg.fmean*sparse(nbdy,1,1,gg.nIJ,1)==0.5); fd1 = intersect(fd1,f1);
fd2 = find(gg.fmean*sparse(nbdy,1,1,gg.nIJ,1)==1); fd2 = intersect(fd2,f2); % connected along-boundary edges
fd = union(fd1,fd2);

ebdy = setdiff(e1,ed); eout = union(e0,ed); 
fbdy = setdiff(f1,fd); fout = union(f0,fd);
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

end
