function [gg] = nevis_grid(nI,nJ,K,xl,xr,yb,yt,oo)
% setup grid and matrix operators
% inputs
%   nI,nJ,K size of grid [x,y,z]
%   xl,xr,yb,yt left and right, bottom and top boundary coordinates
%   oo [optional] option structure
% outputs
%   gg grid struct [see below for contents]
%
% 30 July 2014 : taken from hydro_setup
%   29 August 2014: edited to treat reflective conditions differently, using all
%       edges rather than pointing back to inner edges
%   14 Jan 2015: edited to allow one dimensional input for y, in which case
%       y should be width; also need oo.yperiodic in this case.
%   9 April 2016: note that bottom and top are confusingly named, and inconsistent 
%       with yb and yt.  yb is y coordinate of node with lowest y-index. 'top' edge
%       refers to the edge connecting to the node with a lower y-index (i.e. above in a matrix format) which
%       is actually closer to the yb boundary.
%   28/11/2024: edit the notation, now "top" refers to the top of the physical domain

%% alternative to use existing grid if fewer inputs provided [ unfinished? ]
if nargin<=3,
% alternative input
% [gg] = nevis_grid(x,y,oo)
x = nI; y = nJ; 
if nargin<3, oo = struct; else oo = K; end

nI = length(x); nJ = length(y); K = 1;
xl = min(x); xr = max(x); 
if nJ==1, yb = 0; yt = y; else yb = min(y); yt = max(y); end % use y as width if y has length 1
elseif nargin<8, oo = struct;
end

%% options
if ~isfield(oo,'yperiodic'), oo.yperiodic = 0; end
if ~isfield(oo,'xperiodic'), oo.xperiodic = 0; end

disp('nevis_grid: Setting up grid ...');
disp(['  ',num2str(nI),'-by-',num2str(nJ),' nodes']);
disp(['  ',num2str(K),' vertical layers']);
disp(['  ',num2str(xl),'<x<',num2str(xr),', ',num2str(yb),'<y<',num2str(yt)]);

%% setup rectangular grid
% number of nodes (cells) (n)
nIJ = nI*nJ;
% number of xedges (e)
eI = nI+1;
eJ = nJ;
eIJ = eI*eJ;
% number of yedges (f)
fI = nI;
fJ = nJ+1;
fIJ = fI*fJ;
% number of corners (c)
cI = nI+1;
cJ = nJ+1;
cIJ = cI*cJ;

% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*nI ]
ns = (1:nI*nJ)';
es = (1:eI*eJ)';
fs = (1:fI*fJ)';
cs = (1:cI*cJ)';

%% vertical grid
zz = (0:K-1)/(K-1);    % z = s-(s-b)*zz

%% connections [ these assume periodic boundary conditions ]
% nconnent -- node(cell) connection
% econnect -- x-edge connection
% fconnect -- y-edge connection
% cconnect -- corner connection

% nconnect : [left edge, right edge, bottom edge, top edge]
nconnect = zeros(nIJ,4);
for i = 1:nI
    for j = 1:nJ
        nconnect(i+(j-1)*nI,:) = [ i+(j-1)*eI i+1+(j-1)*eI i+(j-1)*fI i+j*fI ]; 
    end
end

% econnect : [left edge, right edge, bottom edge, top edge]
econnect = zeros(eIJ,4);
for i = 1:eI
    for j = 1:eJ
        if i==1,            econnect(i+(j-1)*eI,:) = [ nI+(j-1)*nI i+(j-1)*nI i+(j-1)*cI i+j*cI ];
        elseif i==eI,       econnect(i+(j-1)*eI,:) = [ i-1+(j-1)*nI 1+(j-1)*nI i+(j-1)*cI i+j*cI ];
        else                econnect(i+(j-1)*eI,:) = [ i-1+(j-1)*nI i+(j-1)*nI i+(j-1)*cI i+j*cI ]; 
        end
    end
end
% correct for reflective edges if the boundary is non-periodic
if ~oo.xperiodic
for j = 1:eJ
    for i = 1, econnect(i+(j-1)*eI,[1 2]) = [ (i+1)+(j-1)*nI i+(j-1)*nI ]; end
    for i = eI, econnect(i+(j-1)*eI,[1 2]) = [ i-1+(j-1)*nI (i-2)+(j-1)*nI ]; end
end
end


% fconnect : [bottom node, top node, left corner, right corner]
fconnect = zeros(fIJ,4);
for i = 1:fI
    for j = 1:fJ
        if j==1,            fconnect(i+(j-1)*fI,:) = [ i+(nJ-1)*nI i+(j-1)*nI i+(j-1)*cI i+1+(j-1)*cI ]; 
        elseif j==fJ,       fconnect(i+(j-1)*fI,:) = [ i+(j-2)*nI i+0*nI i+(j-1)*cI i+1+(j-1)*cI ]; 
        else                fconnect(i+(j-1)*fI,:) = [ i+(j-2)*nI i+(j-1)*nI i+(j-1)*cI i+1+(j-1)*cI ]; 
        end
    end
end
% correct for reflective edges if the boundary is non-periodic
if ~oo.yperiodic
for i = 1:fI
    for j = 1,           fconnect(i+(j-1)*fI,[1 2]) = [ i+(j)*nI i+(j-1)*nI ]; end
    for j = fJ,          fconnect(i+(j-1)*fI,[1 2]) = [ i+(j-2)*nI i+(j-3)*nI ]; end
end
end

% cconnect : [left edge, right edge, bottom edge, top edge]
cconnect = zeros(cIJ,4);
for i = 1:cI
    for j = 1:cJ
        if i==1, 
            if j==1,        cconnect(i+(j-1)*cI,:) = [ fI+(j-1)*fI i+(j-1)*fI i+(eJ-1)*eI i+(j-1)*eI ];
            elseif j==cJ,   cconnect(i+(j-1)*cI,:) = [ fI+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+0*eI ];
            else            cconnect(i+(j-1)*cI,:) = [ fI+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+(j-1)*eI ];
            end
        elseif i==cI, 
            if j==1,        cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI 1+(j-1)*fI i+(eJ-1)*eI i+(j-1)*eI ];
            elseif j==cJ,   cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI 1+(j-1)*fI i+(j-2)*eI i+0*eI ];
            else            cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI 1+(j-1)*fI i+(j-2)*eI i+(j-1)*eI ];
            end
        elseif j==1,        cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI i+(j-1)*fI i+(eJ-1)*eI i+(j-1)*eI ];
        elseif j==cJ,       cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+0*eI ];
        else                cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+(j-1)*eI ];
 
        end
    end
end
% correct for reflective edges if the boundary is non-periodic
if ~oo.yperiodic
for i = 1:cI,
    for j = 1,    cconnect(i+(j-1)*cI,[3 4]) = [ i+(j)*eI i+(j-1)*eI ]; end
    for j = cJ,   cconnect(i+(j-1)*cI,[3 4]) = [ i+(j-2)*eI i+(j-3)*eI ]; end
end
end
if ~oo.xperiodic
for j = 1:cJ,
    for i = 1,    cconnect(i+(j-1)*cI,[1 2]) = [ (i+1)+(j-1)*fI i+(j-1)*fI ]; end
    for i = cI,   cconnect(i+(j-1)*cI,[1 2]) = [ (i-1)+(j-1)*fI (i-2)+(j-1)*fI ]; end
end
end

%% boundaries
% boundary labels
etopbdy = (1:eI)'; ebotbdy = (1:eI)'+(eJ-1)*eI;           %side boundary xedges
elbdy = 1+((1:eJ)'-1)*eI; erbdy = eI+((1:eJ)'-1)*eI;      %end boundary xedges
ftopbdy = (1:fI)'; fbotbdy = (1:fI)'+(fJ-1)*fI;           %side boundary yedges
flbdy = 1+((1:fJ)'-1)*fI; frbdy = fI+((1:fJ)'-1)*fI;      %end boundary yedges
ntopbdy = (1:nI)'; nbotbdy = (1:nI)'+(nJ-1)*nI;           %side boundary nodes
nlbdy = 1+((1:nJ)'-1)*nI; nrbdy = nI+((1:nJ)'-1)*nI;      %end boundary nodes
ctopbdy = (1:cI)'; cbotbdy = (1:cI)'+(cJ-1)*cI;           %side boundary corners
clbdy = 1+((1:cJ)'-1)*cI; crbdy = cI+((1:cJ)'-1)*cI;      %end boundary corners

bdy.nlbdy = nlbdy; bdy.nrbdy = nrbdy; bdy.ntopbdy = ntopbdy; bdy.nbotbdy = nbotbdy;
bdy.elbdy = elbdy; bdy.erbdy = erbdy; bdy.etopbdy = etopbdy; bdy.ebotbdy = ebotbdy;
bdy.flbdy = flbdy; bdy.frbdy = frbdy; bdy.ftopbdy = ftopbdy; bdy.fbotbdy = fbotbdy;
bdy.clbdy = clbdy; bdy.crbdy = crbdy; bdy.ctopbdy = ctopbdy; bdy.cbotbdy = cbotbdy;

% margin nodes/edges/corners
n1 = unique([nlbdy; nrbdy; ntopbdy; nbotbdy]);
e0 = [];
e1 = unique([elbdy; erbdy]);
e2 = unique([etopbdy; ebotbdy]);
f0 = [];
f1 = unique([ftopbdy; fbotbdy]);
f2 = unique([flbdy; frbdy]);
c0 = [];
c1 = [];
c2 = setdiff(unique([clbdy; crbdy; ctopbdy; cbotbdy]),[1 cI 1+(cJ-1)*cI cI+(cJ-1)*cI]');
c3 = [1 cI 1+(cJ-1)*cI cI+(cJ-1)*cI]'; % "corner" corners

%% coordinates [ grid is xl<x<xr, yb<y<yt ]
% [nx,ny,ex,ey,fx,fy,cx,cy] = hydro_coords(nI,nJ,xl,xr,yb,yt,oo);    %[ see hydro_coords.m ]
% [ this will make nodes on boundaries ]
ex = (-1/2+(0:nI))/(nI-1);
fy = (-1/2+(0:nJ))/(nJ-1);
if oo.xperiodic, ex = (-1/2+(0:nI))/(nI); end % [ this will make nodes on boundaries at left and top and bottom, for periodic bcs in x ]
if oo.yperiodic, fy = (-1/2+(0:nJ))/(nJ); end % [ this will make nodes on boundaries at left and right and bottom, for periodic bcs in y ]
nx = (ex(1:end-1)+ex(2:end))/2;
ny = (fy(1:end-1)+fy(2:end))/2;
ey = ny;
fx = nx;
cx = ex;
cy = fy;

%     % [ this will make edges on boundaries ]
%     ex = (0:nI)/nI;
%     fy = (0:nJ)/nJ;
%     nx = (ex(1:end-1)+ex(2:end))/2;
%     ny = (fy(1:end-1)+fy(2:end))/2;
%     ey = ny;
%     fx = nx;
%     cx = ex;
%     cy = fy;

% scale to size
nx = xl+(xr-xl)*nx;
ny = yb+(yt-yb)*ny;
ex = xl+(xr-xl)*ex;
ey = yb+(yt-yb)*ey;
fx = xl+(xr-xl)*fx;
fy = yb+(yt-yb)*fy;
cx = xl+(xr-xl)*cx;
cy = yb+(yt-yb)*cy;

[ny,nx] = meshgrid(ny,nx); % node grid (nI*nJ)
[ey,ex] = meshgrid(ey,ex); % xedge grid (nI+1)*nJ
[fy,fx] = meshgrid(fy,fx); % yedge grid (nI)*(nJ+1)
[cy,cx] = meshgrid(cy,cx); % corner grid (nI+1)*(nJ+1)
    

%% grid spacing [ needed for hydrology ]
% Dx & Dy are defined on nodes
% Dx : x size of nodes [nIJ-by-1]
% Dy : y size of nodes [nIJ-by-1]
Dx = ex(nconnect(:,2))-ex(nconnect(:,1)); 
Dy = fy(nconnect(:,4))-fy(nconnect(:,3));

%% operators [ these work as matrix multipliers ; may have to be careful about boundary locations ]

% nmeanx : node mean over 2 adjacent xedges [nIJ-by-eIJ]  
% nmeany : node mean over 2 adjacent yedges [nIJ-by-fIJ]
% nmeanc : node mean over 4 adjacent corners [nIJ-by-cIJ]
nmeanx = sparse([ns; ns],[nconnect(:,1); nconnect(:,2)],...
    [1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,eIJ);
nmeany = sparse([ns; ns],[nconnect(:,3); nconnect(:,4)],...
    [1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,fIJ);
nmeanc = sparse([ns; ns; ns; ns],...
    [econnect(nconnect(:,1),3); econnect(nconnect(:,1),4);...
    econnect(nconnect(:,2),3); econnect(nconnect(:,2),4)],...
    [1/4*ones(nIJ,1); 1/4*ones(nIJ,1); 1/4*ones(nIJ,1); 1/4*ones(nIJ,1)],...
    nIJ,cIJ);

% emean : xedge mean over 2 adjacent nodes [eIJ-by-nIJ]   
% [ careful at left and right boundaries ]
emean = sparse([es; es],[econnect(:,1); econnect(:,2)],...
    [1/2*ones(eIJ,1); 1/2*ones(eIJ,1)],eIJ,nIJ);

% fmean : yedge mean over 2 adjacent nodes [fIJ-by-nIJ]  
% [ careful at top and bottom boundaries ]
fmean = sparse([fs; fs],[fconnect(:,1); fconnect(:,2)],...
    [1/2*ones(fIJ,1); 1/2*ones(fIJ,1)],fIJ,nIJ);

% cmean : corner mean over 4 adjacent nodes [cIJ-by-nIJ] 
% [ careful at boundaries ]
cmean = sparse([cs; cs; cs; cs],...
    [fconnect(cconnect(:,1),1); fconnect(cconnect(:,1),2);...
    fconnect(cconnect(:,2),1); fconnect(cconnect(:,2),2)],...
    [1/4*ones(cIJ,1); 1/4*ones(cIJ,1); 1/4*ones(cIJ,1); 1/4*ones(cIJ,1)],...
    cIJ,nIJ);

% cmeans : corner mean over 2 adjacent nodes [cIJ-by-nIJ] 
% cmeanr : corner mean over 2 adjacent nodes [cIJ-by-nIJ]
% [ careful at boundaries ]
cmeans = sparse([cs; cs],...
    [fconnect(cconnect(:,1),1); fconnect(cconnect(:,2),2)],...
    [1/2*ones(cIJ,1); 1/2*ones(cIJ,1)],cIJ,nIJ);
cmeanr = sparse([cs; cs],...
    [fconnect(cconnect(:,1),2); fconnect(cconnect(:,2),1)],...
    [1/2*ones(cIJ,1); 1/2*ones(cIJ,1)],cIJ,nIJ);

% nmeans : node mean over 2 adjacent corners [nIJ-by-cIJ] (diagonal)
% nmeanr : node mean over 2 adjacent corners [nIJ-by-cIJ]
nmeans = sparse([ns; ns],...
    [econnect(nconnect(:,1),3); econnect(nconnect(:,2),4)],...
    [1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,cIJ);
nmeanr = sparse([ns; ns],...
    [econnect(nconnect(:,1),4); econnect(nconnect(:,2),3)],...
    [1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,cIJ);

% nddx : x derivative over node [nIJ-by-eIJ]
% nddy : y derivative over node [nIJ-by-fIJ]
dx = Dx; 
nddx = sparse([ns; ns],...
    [nconnect(:,1); nconnect(:,2)],[-dx.^(-1); dx.^(-1)],nIJ,eIJ);
dy = Dy;
nddy = sparse([ns; ns],...
    [nconnect(:,3); nconnect(:,4)],[-dy.^(-1); dy.^(-1)],nIJ,fIJ);

% eddx : x derivative over xedge [eIJ-by-nIJ]      
% eddy : y derivative over xedge [eIJ-by-cIJ]    
dx = abs(nx(econnect(:,2))-nx(econnect(:,1))); 
dy = abs(cy(econnect(:,4))-cy(econnect(:,3)));
% derivatives on bounds [ correct left and right boundaries ]
if oo.xperiodic 
    dx(elbdy) = nx(econnect(elbdy,2))-nx(econnect(elbdy,1))+(xr-xl); 
    dx(erbdy) = nx(econnect(erbdy,2))-nx(econnect(erbdy,1))+(xr-xl); 
end %
eddx = sparse([es; es],...
    [econnect(:,1); econnect(:,2)],[-dx.^(-1); dx.^(-1)],eIJ,nIJ);
eddy = sparse([es; es],...
    [econnect(:,3); econnect(:,4)],[-dy.^(-1); dy.^(-1)],eIJ,cIJ);

% fddy : y derivative over yedge [fIJ-by-nIJ]       
% fddx : x derivative over yedge [fIJ-by-cIJ]       
dy = abs(ny(fconnect(:,2))-ny(fconnect(:,1)));
dx = abs(cx(fconnect(:,4))-cx(fconnect(:,3)));
if oo.yperiodic %[ correct top and bottom boundaries ]
    dy(ftopbdy) = ny(fconnect(ftopbdy,2))-ny(fconnect(ftopbdy,1))+(yt-yb); 
    dy(fbotbdy) = ny(fconnect(fbotbdy,2))-ny(fconnect(fbotbdy,1))+(yt-yb); 
end 
fddy = sparse([fs; fs],...
    [fconnect(:,1); fconnect(:,2)],[-dy.^(-1); dy.^(-1)],fIJ,nIJ);
fddx = sparse([fs; fs],...
    [fconnect(:,3); fconnect(:,4)],[-dx.^(-1); dx.^(-1)],fIJ,cIJ);

% cddx : x derivative over corner [cIJ-by-fIJ]      
% cddy : y derivative over corner [cIJ-by-eIJ]     
dx = abs(fx(cconnect(:,2))-fx(cconnect(:,1))); 
dy = abs(ey(cconnect(:,4))-ey(cconnect(:,3))); 
if oo.xperiodic %[ correct left and right boundaries ]
    dx(clbdy) = fx(cconnect(clbdy,2))-fx(cconnect(clbdy,1))+(xr-xl);  
    dx(crbdy) = fx(cconnect(crbdy,2))-fx(cconnect(crbdy,1))+(xr-xl); 
end        
if oo.yperiodic %[ correct top and bottom boundaries ]
    dy(ctopbdy) = ey(cconnect(ctopbdy,4))-ey(cconnect(ctopbdy,3))+(yt-yb); 
    dy(cbotbdy) = ey(cconnect(cbotbdy,4))-ey(cconnect(cbotbdy,3))+(yt-yb); 
end  
cddx = sparse([cs; cs],...
    [cconnect(:,1); cconnect(:,2)],[-dx.^(-1); dx.^(-1)],cIJ,fIJ);
cddy = sparse([cs; cs],...
    [cconnect(:,3); cconnect(:,4)],[-dy.^(-1); dy.^(-1)],cIJ,eIJ);

% cdds : s derivative over corner [cIJ-by-nIJ]      
% cddr : s derivative over corner [cIJ-by-nIJ] 
ds = (dx.^2+dy.^2).^(1/2); 
dr = (dx.^2+dy.^2).^(1/2);
cdds = sparse([cs; cs],...
    [fconnect(cconnect(:,1),1); fconnect(cconnect(:,2),2)],...
    [-ds.^(-1); ds.^(-1)],cIJ,nIJ);
cddr = sparse([cs; cs],...
    [fconnect(cconnect(:,1),2); fconnect(cconnect(:,2),1)],...
    [-dr.^(-1); dr.^(-1)],cIJ,nIJ);

% ndds : corner divergence over node [nIJ-by-cIJ]
% nddr : corner divergence over node [nIJ-by-cIJ]
ndds = sparse([ns; ns],...
    [econnect(nconnect(:,1),3); econnect(nconnect(:,2),4)],...
    [-Dx.^(-1).*Dy.^(-1); Dx.^(-1).*Dy.^(-1)],nIJ,cIJ);
nddr = sparse([ns; ns],...
    [econnect(nconnect(:,1),4); econnect(nconnect(:,2),3)],...
    [-Dx.^(-1).*Dy.^(-1); Dx.^(-1).*Dy.^(-1)],nIJ,cIJ);

%% default no boundaries 
nbdy = []; nout = [];
ebdy = []; eout = [];
fbdy = []; fout = [];
cbdy = []; cout = [];

%% define active and interior labels
ns = setdiff(ns,nout);   
es = setdiff(es,eout);   
fs = setdiff(fs,fout);   
cs = setdiff(cs,cout);  

nin = setdiff(ns,nbdy);   
ein = setdiff(es,ebdy);  
fin = setdiff(fs,fbdy);  
cin = setdiff(cs,cbdy);  

%% outputs
gg.xl = xl;
gg.xr = xr;
gg.yb = yb;
gg.yt = yt;
gg.nI = nI;
gg.nJ = nJ;
gg.nIJ = nIJ;
gg.eI = eI;
gg.eJ = eJ;
gg.eIJ = eIJ;
gg.fI = fI;
gg.fJ = fJ;
gg.fIJ = fIJ;
gg.cI = cI;
gg.cJ = cJ;
gg.cIJ = cIJ;

gg.K = K; % number of vertical layers
gg.bdy = bdy; % strcut for boundaries
gg.Dx = Dx; % dx in grid
gg.Dy = Dy; % dy in grid
gg.Ds = (Dx.^2+Dy.^2).^(1/2); % ds along channel
gg.Dr = gg.Ds;

% computational grid
gg.nx = nx; % node
gg.ny = ny;
gg.ex = ex; % x-edge
gg.ey = ey; 
gg.fx = fx; % y-edge
gg.fy = fy;
gg.cx = cx; % corner
gg.cy = cy;
gg.zz = zz; % vertical

% matrix operators
gg.nmeanx = nmeanx;
gg.nmeany = nmeany;
gg.nmeanc = nmeanc;
gg.emean = emean;
gg.fmean = fmean;
gg.cmean = cmean;
gg.nddx = nddx;
gg.nddy = nddy;
gg.eddx = eddx;
gg.eddy = eddy;
gg.fddx = fddx;
gg.fddy = fddy;
gg.cddx = cddx;
gg.cddy = cddy;
gg.nconnect = nconnect;  
gg.econnect = econnect;  
gg.fconnect = fconnect; 
gg.cconnect = cconnect;  
gg.cdds = cdds;
gg.cddr = cddr;
gg.ndds = ndds;
gg.nddr = nddr;
gg.nmeans = nmeans;
gg.nmeanr = nmeanr;
gg.cmeans = cmeans;
gg.cmeanr = cmeanr;

% element index
gg.ns = ns;
gg.es = es;
gg.fs = fs;
gg.cs = cs;

% boundary element index
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

% boundary element label
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

disp('nevis_grid: Done');
end
