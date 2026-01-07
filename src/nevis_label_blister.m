function [gg] = nevis_label_blister(gg,d_inds,oo)
% label makes the boundary labels and list of boundary indices
% 
% inputs
%   gg grid structure [see below for contents]
%   d_inds indices of nodes at which pressure prescribed
%   oo [optional] option structure with method to use [see below for options]
% outputs:
%   gg grid structure [added/modified labels]

if nargin<4, oo = struct; end

%% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*I ]
ns = (1:gg.nIJ)';
es = (1:gg.eIJ)';
fs = (1:gg.fIJ)';
cs = (1:gg.cIJ)';

%% labels of boundary nodes, edges, corners
nout = gg.nout_blister;
n1 = gg.n1_blister;
e0 = gg.e0_blister;
e1 = gg.e1_blister;
e2 = gg.e2_blister;
f0 = gg.f0_blister;
f1 = gg.f1_blister;
f2 = gg.f2_blister;
c0 = gg.c0_blister;
c1 = gg.c1_blister;
c2 = gg.c2_blister;
c3 = gg.c3_blister;
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
% 20 Nov 2025: The template below is for hydrology models
ns = setdiff(ns,nout);  %  
es = setdiff(es,eout);   
fs = setdiff(fs,fout);   
cs = setdiff(cs,cout);  

nin = setdiff(ns,nbdy);   
ein = setdiff(es,ebdy);  
fin = setdiff(fs,fbdy);  
cin = setdiff(cs,cbdy);  

% internal+bound elements
gg.ns_blister = ns;
gg.es_blister = es;
gg.fs_blister = fs;
gg.cs_blister = cs;

% bound elements
gg.nbdy_blister = nbdy;
gg.ebdy_blister = ebdy;
gg.cbdy_blister = cbdy;
gg.fbdy_blister = fbdy;

% internal elements
gg.nin_blister = nin;
gg.ein_blister = ein;
gg.fin_blister = fin;
gg.cin_blister = cin;

% external elements
gg.nout_blister = nout;
gg.eout_blister = eout;
gg.fout_blister = fout;
gg.cout_blister = cout;

end
