function gg = nevis_label_ice(gg,oo)
% update grid with boundary labels for solving for ice velocity
% 
% 10 Oct 2025 - edited from ssa/mask_and_label.m and nevis_label, to create grid structure
%   with boundary node/edge labels for ice velocity solve, which may be
%   different from those needed by the hydrology model

if nargin<2, oo = struct; end

%% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*I ]
ns = (1:gg.nIJ)';
es = (1:gg.eIJ)';
fs = (1:gg.fIJ)';
cs = (1:gg.cIJ)';

%% labels of boundary nodes, edges, corners (from nevis_mask)
nout = gg.nout;
n1 = gg.n1; % nodes connected by an edge to an outside node : boundary nodes
e0 = gg.e0; % edges connected to all outside nodes
e1 = gg.e1; % edges connected to one outside node : boundary edges
e2 = gg.e2; % edges connected to two edge nodes : along-boundary edges
f0 = gg.f0; % edges connected to all outside nodes
f1 = gg.f1; % edges connected to one outside node : boundary edges
f2 = gg.f2; % edges connected to two edge nodes : along-boundary edges
c0 = gg.c0; % corners connected to all outside nodes
c1 = gg.c1; % corners connected to one outside node 
c2 = gg.c2; % corners connected to two outside nodes
c3 = gg.c3; % corners connected to three outside nodes

temp = zeros(gg.eIJ,1); temp(e1) = 1; n1x = find(gg.nmeanx*temp>0); n1x = setdiff(n1x,nout); % nodes connected by an e edge to an outside node : x boundary nodes 
temp = zeros(gg.fIJ,1); temp(f1) = 1; n1y = find(gg.nmeany*temp>0); n1y = setdiff(n1y,nout); % nodes connected by an f edge to an outside node : y boundary nodes 

%% default boundary labels for whole grid
nlbdy = gg.bdy.nlbdy; nrbdy = gg.bdy.nrbdy; ntopbdy = gg.bdy.ntopbdy; nbotbdy = gg.bdy.nbotbdy;
elbdy = gg.bdy.elbdy; erbdy = gg.bdy.erbdy; etopbdy = gg.bdy.etopbdy; ebotbdy = gg.bdy.ebotbdy;
flbdy = gg.bdy.flbdy; frbdy = gg.bdy.frbdy; ftopbdy = gg.bdy.ftopbdy; fbotbdy = gg.bdy.fbotbdy;
clbdy = gg.bdy.clbdy; crbdy = gg.bdy.crbdy; ctopbdy = gg.bdy.ctopbdy; cbotbdy = gg.bdy.cbotbdy;

if isfield(oo,'boundary_method') && strcmp(oo.boundary_method,'per_tblr')
% % [ appropriate for periodic top and bottom boundaries, periodic left and right ]
% ebdy = unique([]);       
% fbdy = unique([]);         
% nbdyx = unique([]); 
% nbdyy = unique([]); 
% cbdy = unique([]); 
% eout = unique([]);
% fout = unique([]);
% cout = unique([]);

elseif isfield(oo,'boundary_method') && strcmp(oo.boundary_method,'vel_tblr')
% [ appropriate for prescribing velocity on all boundaries ]
ebdy = unique([etopbdy; ebotbdy; elbdy; erbdy]);       
fbdy = unique([ftopbdy; fbotbdy; flbdy; frbdy]);         
nbdyx = unique([]); 
nbdyy = unique([]); 
cbdy = unique([]); 
eout = unique([]);
fout = unique([]);
cout = unique([ctopbdy; cbotbdy; clbdy; crbdy]);

elseif isfield(oo,'boundary_method') && strcmp(oo.boundary_method,'stress_tblr')
% [ appropriate for prescribing normal and tangential stress on all boundaries]
ebdy = unique([]);
fbdy = unique([]);
nbdyx = unique([nlbdy; nrbdy]); 
nbdyy = unique([nbotbdy; ntopbdy]); 
cbdy = unique([cbotbdy; ctopbdy; clbdy; crbdy]);
eout = unique([ebotbdy; etopbdy; elbdy; erbdy]);
fout = unique([fbotbdy; ftopbdy; flbdy; frbdy]);
cout = unique([]);

else
% [ appropriate for prescribing velocity on left boundary, reflective or periodic on top and bottom, and normal and tangential stress on right]
ebdy = unique([elbdy]);
fbdy = unique([flbdy]);
nbdyx = unique([nrbdy]); 
nbdyy = unique([]); 
cbdy = unique([crbdy]);
eout = unique([erbdy]);
fout = unique([]);
cout = unique([clbdy]);

end
% [ could add more alternatives here ]

% % [ appropriate for prescribing velocity on left boundary, normal velocity and tangential stress on top and bottom, and normal and tangential stress on right]
% ebdy = unique([elbdy]);
% fbdy = unique([flbdy; ftopbdy; fbotbdy]);
% nbdyx = unique([nrbdy]); 
% nbdyy = unique([]); 
% cbdy = unique([ctopbdy; cbotbdy; crbdy]);
% eout = unique([erbdy]);
% fout = unique([]);
% cout = unique([clbdy]);

% % [ appropriate for prescribing velocity on left, top and bottom boundaries, and normal and tangential stress on right]
% ebdy = unique([elbdy; etopbdy; ebotbdy]);
% fbdy = unique([flbdy; ftopbdy; fbotbdy]);
% nbdyx = unique([nrbdy]); 
% nbdyy = unique([]); 
% cbdy = unique([crbdy]);
% eout = unique([erbdy]);
% fout = unique([]);
% cout = unique([clbdy]);


%% modification to boundary labels for masked region
% [ appropriate for prescribing velocity on all masked boundaries ]
nbdyx = [nbdyx];
nbdyy = [nbdyy];
ebdy = unique([ebdy; e1; e2]); 
fbdy = unique([fbdy; f1; f2]); 
cbdy = [cbdy]; 

nout = unique([nout]);
eout = unique([eout; e0]);
fout = unique([fout; f0]);
cout = unique([cout; c0; c3; c2]);

% [ appropriate for prescribing stress on labelled nodes ]
% [ not included yet ]

%% restrict labels to active area (i.e. remove points now assigned as out)
ns = setdiff(ns,nout);   
es = setdiff(es,eout);   
fs = setdiff(fs,fout);   
cs = setdiff(cs,cout);  

nbdyx = setdiff(nbdyx,nout); 
nbdyy = setdiff(nbdyy,nout); 
nbdy = unique([nbdyx; nbdyy]);   
ebdy = setdiff(ebdy,eout);   
fbdy = setdiff(fbdy,fout);   
cbdy = setdiff(cbdy,cout); 

nin = setdiff(ns,nbdy);   
ein = setdiff(es,ebdy);  
fin = setdiff(fs,fbdy);  
cin = setdiff(cs,cbdy);   

% [ possibly want to redefine mean operators here ? ]

%% assign labels to gg (ice labels with suffix 2 may be slightly different to hydrology labels)
gg.ns2 = ns;
gg.es2 = es;
gg.fs2 = fs;
gg.cs2 = cs;
gg.nbdyx2 = nbdyx;
gg.nbdyy2 = nbdyy;
gg.ebdy2 = ebdy;
gg.cbdy2 = cbdy;
gg.fbdy2 = fbdy;
gg.nin2 = nin;
gg.ein2 = ein;
gg.fin2 = fin;
gg.cin2 = cin;
gg.nout2 = nout;
gg.eout2 = eout;
gg.fout2 = fout;
gg.cout2 = cout;

end
