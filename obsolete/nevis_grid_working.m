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
%   29 August, edited fconnect and cconnect for yreflective, though might not matter
%      also added option of xreflective
    

%% alternative to use existing grid if fewer inputs provided [ unfinished? ]
if nargin<=3,
% alternative input
% [gg] = nevis_grid(x,y,oo)
nI = length(x); nJ = length(y); K = 1;
xl = min(x); xr = max(x); yb = min(y); yt = max(y);
if nargin<3, oo = struct; end
elseif nargin<8, oo = struct;
end

if ~isfield(oo,'yreflective'), oo.yreflective = 0; end
if ~isfield(oo,'xreflective'), oo.xreflective = 0; end

disp('Setting up grid :');
disp(['  ',num2str(nI),'-by-',num2str(nJ),' nodes']);
disp(['  ',num2str(K),' vertical layers']);
disp(['  ',num2str(xl),'<x<',num2str(xr),', ',num2str(yb),'<y<',num2str(yt)]);
disp('...');

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

% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*I ]
ns = (1:nI*nJ)';
es = (1:eI*eJ)';
fs = (1:fI*fJ)';
cs = (1:cI*cJ)';

%% vertical grid
zz = (0:K-1)/(K-1);    % z = s-(s-b)*zz

%% connections [ these assume periodic boundary conditions ]
% nconnect : [left edge, right edge, top edge, bottom edge]
nconnect = zeros(nIJ,4);
for i = 1:nI
    for j = 1:nJ
        nconnect(i+(j-1)*nI,:) = [ i+(j-1)*eI i+1+(j-1)*eI i+(j-1)*fI i+j*fI ]; 
    end
end
if oo.yreflective, 
for i = 1:nI
        for j = 1, nconnect(i+(j-1)*nI,[3 4]) = [ i+(j)*fI i+j*fI ]; end
        for j = nJ, nconnect(i+(j-1)*nI,[3 4]) = [ i+(j-1)*fI i+(j-1)*fI ]; end
end
end
if oo.xreflective
for j = 1:nJ
        for i = 1, nconnect(i+(j-1)*nI,[1 2]) = [ (i+1)+(j-1)*eI i+1+(j-1)*eI ]; end
        for i = nI, nconnect(i+(j-1)*nI,[1 2]) = [ i+(j-1)*eI i+(j-1)*eI ]; end
end
end

% econnect : [left node, right node, top corner, bottom corner]
econnect = zeros(eIJ,4);
for i = 1:eI
    for j = 1:eJ
        if i==1,            econnect(i+(j-1)*eI,:) = [ nI+(j-1)*nI i+(j-1)*nI i+(j-1)*cI i+j*cI ];
        elseif i==eI,       econnect(i+(j-1)*eI,:) = [ i-1+(j-1)*nI 1+(j-1)*nI i+(j-1)*cI i+j*cI ];
        else                econnect(i+(j-1)*eI,:) = [ i-1+(j-1)*nI i+(j-1)*nI i+(j-1)*cI i+j*cI ]; 
        end
    end
end
if oo.yreflective
for i = 1:eI, 
    for j = 1, econnect(i+(j-1)*eI,[3 4]) = [ i+(j)*cI i+j*cI ]; end
    for j = eJ, econnect(i+(j-1)*eI,[3 4]) = [ i+(j-1)*cI i+(j-1)*cI  ]; end
end
end
if oo.xreflective
for j = 1:eJ, 
    for i = 1, econnect(i+(j-1)*eI,[1 2]) = [ (i+1)+(j-1)*nI i+(j-1)*nI ]; end
    for i = eI, econnect(i+(j-1)*eI,[1 2]) = [ i-1+(j-1)*nI (i-2)+(j-1)*nI ]; end
end
end


% fconnect : [top node, bottom node, left corner, right corner]
fconnect = zeros(fIJ,4);
for i = 1:fI
    for j = 1:fJ
        if j==1,            fconnect(i+(j-1)*fI,:) = [ i+(nJ-1)*nI i+(j-1)*nI i+(j-1)*cI i+1+(j-1)*cI ]; 
        elseif j==fJ,       fconnect(i+(j-1)*fI,:) = [ i+(j-2)*nI i+0*nI i+(j-1)*cI i+1+(j-1)*cI ]; 
        else                fconnect(i+(j-1)*fI,:) = [ i+(j-2)*nI i+(j-1)*nI i+(j-1)*cI i+1+(j-1)*cI ]; 
        end
    end
end
% if oo.yreflective, top and bottom edges should not be accessed so should not matter
if oo.yreflective
for i = 1:fI
    for j = 1,           fconnect(i+(j-1)*fI,[1 2]) = [ i+(j)*nI i+(j-1)*nI ]; end
    for j = fJ,          fconnect(i+(j-1)*fI,[1 2]) = [ i+(j-2)*nI i+(j-3)*nI ]; end
end
end
if oo.xreflective
for j = 1:fJ
    for i = 1,           fconnect(i+(j-1)*fI,[3 4]) = [ (i+1)+(j-1)*cI i+1+(j-1)*cI  ]; end
    for i = fI,          fconnect(i+(j-1)*fI,[3 4]) = [ i+(j-1)*cI (i)+(j-1)*cI  ]; end
end
end
% cconnect : [left edge, right edge, top edge, bottom edge]
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
% if oo.yreflective, top and bottom edges should not be accessed so should not matter
if oo.yreflective
for i = 1:cI,
    for j = 1,    cconnect(i+(j-1)*cI,[3 4]) = [ i+(j)*eI i+(j-1)*eI ]; end
    for j = cJ,   cconnect(i+(j-1)*cI,[3 4]) = [ i+(j-2)*eI i+(j-3)*eI ]; end
end
end
if oo.xreflective
for j = 1:cJ,
    for i = 1,    cconnect(i+(j-1)*cI,[1 2]) = [ (i+1)+(j-1)*fI i+(j-1)*fI ]; end
    for i = cI,   cconnect(i+(j-1)*cI,[1 2]) = [ (i-1)+(j-1)*fI (i-2)+(j-1)*fI ]; end
end
end

% % cconnectn : [left bottom node, right bottom node, left top node, right top node]
% cconnect = zeros(cIJ,4);
% for i = 1:cI
%     for j = 1:cJ
%         if i==1, 
%             if j==1,        cconnect(i+(j-1)*cI,:) = [ fI+(j-1)*fI i+(j-1)*fI i+(eJ-1)*eI i+(j-1)*eI ];
%             elseif j==cJ,   cconnect(i+(j-1)*cI,:) = [ fI+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+0*eI ];
%             else            cconnect(i+(j-1)*cI,:) = [ fI+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+(j-1)*eI ];
%             end
%         elseif i==cI, 
%             if j==1,        cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI 1+(j-1)*fI i+(eJ-1)*eI i+(j-1)*eI ];
%             elseif j==cJ,   cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI 1+(j-1)*fI i+(j-2)*eI i+0*eI ];
%             else            cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI 1+(j-1)*fI i+(j-2)*eI i+(j-1)*eI ];
%             end
%         elseif j==1,        cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI i+(j-1)*fI i+(eJ-1)*eI i+(j-1)*eI ];
%         elseif j==cJ,       cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+0*eI ];
%         else                cconnect(i+(j-1)*cI,:) = [ (i-1)+(j-1)*fI i+(j-1)*fI i+(j-2)*eI i+(j-1)*eI ];
%  
%         end
%     end
% end

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

%% coordinates [ grid is xl<x<xr, yb<y<yt ]
% [nx,ny,ex,ey,fx,fy,cx,cy] = hydro_coords(nI,nJ,xl,xr,yb,yt,oo);    %[ see hydro_coords.m ]
% [ this will make nodes on boundaries ]
ex = (-1/2+(0:nI))/(nI-1);
fy = (-1/2+(0:nJ))/(nJ-1);
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
%     % [ this will make nodes on boundaries at left and right and bottom, for periodic bcs in y ]
%     ex = (-1/2+(0:nI))/(nI-1);
%     fy = (-1/2+(0:nJ))/(nJ);
%     nx = (ex(1:end-1)+ex(2:end))/2;
%     ny = (fy(1:end-1)+fy(2:end))/2;
%     ey = ny;
%     fx = nx;
%     cx = ex;
%     cy = fy;
%     % [ this will make nodes on boundaries at left and bottom, for periodic bcs in x and y ]
%     ex = (-1/2+(0:nI))/(nI);
%     fy = (-1/2+(0:nJ))/(nJ);
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

[ny,nx] = meshgrid(ny,nx);
[ey,ex] = meshgrid(ey,ex);
[fy,fx] = meshgrid(fy,fx);
[cy,cx] = meshgrid(cy,cx);
    

%% grid spacing [ needed for hydrology ]
% Dx : x size of nodes [nIJ-by-1]
% Dy : y size of nodes [nIJ-by-1]
Dx = ex(nconnect(:,2))-ex(nconnect(:,1));
Dy = fy(nconnect(:,4))-fy(nconnect(:,3));
if oo.yreflective,
    Dy(ntopbdy) = fy(nconnect(ntopbdy,4))-fy(ftopbdy);
    Dy(nbotbdy) = fy(fbotbdy)-fy(nconnect(nbotbdy,3));
end
if oo.xreflective,
    Dx(nlbdy) = ex(nconnect(nlbdy,2))-ex(elbdy);
    Dx(nrbdy) = ex(erbdy)-ex(nconnect(nrbdy,1));
end

%% operators [ these work as matrix multipliers ; may have to be careful about boundary locations ]

% nmeanx : node mean over 2 adjacent xedges [nIJ-by-eIJ]  
% nmeany : node mean over 2 adjacent yedges [nIJ-by-fIJ]
% nmeanx : node mean over 4 adjacent corners [nIJ-by-cIJ]
nmeanx = sparse([ns; ns],[nconnect(:,1); nconnect(:,2)],[1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,eIJ);
nmeany = sparse([ns; ns],[nconnect(:,3); nconnect(:,4)],[1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,fIJ);
nmeanc = sparse([ns; ns; ns; ns],[econnect(nconnect(:,1),3); econnect(nconnect(:,1),4); econnect(nconnect(:,2),3); econnect(nconnect(:,2),4)],[1/4*ones(nIJ,1); 1/4*ones(nIJ,1); 1/4*ones(nIJ,1); 1/4*ones(nIJ,1)],nIJ,cIJ);
% emean : xedge mean over 2 adjacent nodes [eIJ-by-nIJ]   [ careful at left and right boundaries ]
emean = sparse([es; es],[econnect(:,1); econnect(:,2)],[1/2*ones(eIJ,1); 1/2*ones(eIJ,1)],eIJ,nIJ);
% fmean : yedge mean over 2 adjacent nodes [fIJ-by-nIJ]  [ careful at top and bottom boundaries ]
fmean = sparse([fs; fs],[fconnect(:,1); fconnect(:,2)],[1/2*ones(fIJ,1); 1/2*ones(fIJ,1)],fIJ,nIJ);
% cmean : corner mean over 4 adjacent nodes [cIJ-by-nIJ] [ careful at boundaries ]
cmean = sparse([cs; cs; cs; cs],[fconnect(cconnect(:,1),1); fconnect(cconnect(:,1),2); fconnect(cconnect(:,2),1); fconnect(cconnect(:,2),2)],[1/4*ones(cIJ,1); 1/4*ones(cIJ,1); 1/4*ones(cIJ,1); 1/4*ones(cIJ,1)],cIJ,nIJ);
% nmeans : node mean over 2 adjacent corners [nIJ-by-cIJ]
% nmeanr : node mean over 2 adjacent corners [nIJ-by-cIJ]
nmeans = sparse([ns; ns],[econnect(nconnect(:,1),3); econnect(nconnect(:,2),4)],[1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,cIJ);
nmeanr = sparse([ns; ns],[econnect(nconnect(:,1),4); econnect(nconnect(:,2),3)],[1/2*ones(nIJ,1); 1/2*ones(nIJ,1)],nIJ,cIJ);
% for reflective bcs, s and r must be switched on boundary : these definitions assume nmeans and nmeanr will always be used together
if oo.yreflective, 
for i = 1:length(ntopbdy),
    nmeans(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),3); econnect(nconnect(ntopbdy(i),2),4)]) = [0; 1];
    nmeanr(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),4); econnect(nconnect(ntopbdy(i),2),3)]) = [1; 0];
end
for i = 1:length(nbotbdy),
    nmeans(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),3); econnect(nconnect(nbotbdy(i),2),4)]) = [1; 0];
    nmeanr(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),4); econnect(nconnect(nbotbdy(i),2),3)]) = [0; 1];
end
end
if oo.xreflective, 
for i = 1:length(nlbdy),
    nmeans(nlbdy(i),[econnect(nconnect(nlbdy(i),1),3); econnect(nconnect(nlbdy(i),2),4)]) = [0; 1];
    nmeanr(nlbdy(i),[econnect(nconnect(nlbdy(i),1),4); econnect(nconnect(nlbdy(i),2),3)]) = [0; 1];
end
for i = 1:length(nrbdy),
    nmeans(nrbdy(i),[econnect(nconnect(nrbdy(i),1),3); econnect(nconnect(nrbdy(i),2),4)]) = [1; 0];
    nmeanr(nrbdy(i),[econnect(nconnect(nrbdy(i),1),4); econnect(nconnect(nrbdy(i),2),3)]) = [1; 0];
end
end
if oo.yreflective && oo.xreflective, 
% note the two entries given by econnect here should have the same index
% for the corner nodes, so only the second of the entries assigned here
% will be used
for i = 1,
    nmeans(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),3); econnect(nconnect(ntopbdy(i),2),4)]) = [0; 2];
    nmeanr(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),4); econnect(nconnect(ntopbdy(i),2),3)]) = [0; 0];
end
for i = length(ntopbdy),
    nmeans(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),3); econnect(nconnect(ntopbdy(i),2),4)]) = [0; 0];
    nmeanr(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),4); econnect(nconnect(ntopbdy(i),2),3)]) = [0; 2];
end
for i = 1,
    nmeans(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),3); econnect(nconnect(nbotbdy(i),2),4)]) = [0; 0];
    nmeanr(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),4); econnect(nconnect(nbotbdy(i),2),3)]) = [0; 2];
end
for i = length(nbotbdy),
    nmeans(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),3); econnect(nconnect(nbotbdy(i),2),4)]) = [0; 2];
    nmeanr(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),4); econnect(nconnect(nbotbdy(i),2),3)]) = [0; 0];
end
end

% nddx : x derivative over node [nIJ-by-eIJ]
% nddy : y derivative over node [nIJ-by-fIJ]
dx = Dx; 
nddx = sparse([ns; ns],[nconnect(:,1); nconnect(:,2)],[-dx.^(-1); dx.^(-1)],nIJ,eIJ);
dy = Dy;
nddy = sparse([ns; ns],[nconnect(:,3); nconnect(:,4)],[-dy.^(-1); dy.^(-1)],nIJ,fIJ);
if oo.yreflective, 
    nddy = nddy + sparse(ntopbdy,nconnect(ntopbdy,4),2*dy(ntopbdy).^(-1),nIJ,fIJ) ...
            + sparse(nbotbdy,nconnect(nbotbdy,3),-2*dy(nbotbdy).^(-1),nIJ,fIJ); 
end
if oo.xreflective, 
    nddx = nddx + sparse(nlbdy,nconnect(nlbdy,2),2*dx(nlbdy).^(-1),nIJ,eIJ) ...
            + sparse(nrbdy,nconnect(nrbdy,1),-2*dx(nrbdy).^(-1),nIJ,eIJ); 
end
% eddx : x derivative over xedge [eIJ-by-nIJ]      
% eddy : y derivative over xedge [eIJ-by-cIJ]       
dx = nx(econnect(:,2))-nx(econnect(:,1)); dx(elbdy) = dx(elbdy)+(xr-xl); dx(erbdy) = dx(erbdy)+(xr-xl); dx(dx==0) = inf;        %[ correct left and right boundaries ]
eddx = sparse([es; es],[econnect(:,1); econnect(:,2)],[-dx.^(-1); dx.^(-1)],eIJ,nIJ);
dy = cy(econnect(:,4))-cy(econnect(:,3));
if oo.yreflective, dy(etopbdy) = cy(econnect(etopbdy,4))-cy(ctopbdy);  dy(ebotbdy) = cy(cbotbdy)-cy(econnect(ebotbdy,3)); end  
eddy = sparse([es; es],[econnect(:,3); econnect(:,4)],[-dy.^(-1); dy.^(-1)],eIJ,cIJ);
if oo.yreflective, 
    eddy = eddy + sparse(etopbdy,econnect(etopbdy,4),2*dy(etopbdy).^(-1),eIJ,cIJ) ... 
            + sparse(ebotbdy,econnect(ebotbdy,3),-2*dy(ebotbdy).^(-1),eIJ,cIJ); 
end
% fddy : y derivative over yedge [fIJ-by-nIJ]       
% fddx : x derivative over yedge [fIJ-by-cIJ]       
dy = ny(fconnect(:,2))-ny(fconnect(:,1)); dy(ftopbdy) = dy(ftopbdy)+(yt-yb); dy(fbotbdy) = dy(fbotbdy)+(yt-yb); dy(dy==0) = inf; %[ correct top and bottom boundaries ]
fddy = sparse([fs; fs],[fconnect(:,1); fconnect(:,2)],[-dy.^(-1); dy.^(-1)],fIJ,nIJ);
dx = cx(fconnect(:,4))-cx(fconnect(:,3));
if oo.xreflective, dx(flbdy) = cx(fconnect(flbdy,4))-cx(clbdy);  dx(frbdy) = cx(crbdy)-cx(fconnect(frbdy,3)); end  
fddx = sparse([fs; fs],[fconnect(:,3); fconnect(:,4)],[-dx.^(-1); dx.^(-1)],fIJ,cIJ);
if oo.xreflective, 
    fddx = fddx + sparse(flbdy,fconnect(flbdy,4),2*dx(flbdy).^(-1),fIJ,cIJ) ... 
            + sparse(frbdy,fconnect(frbdy,3),-2*dx(frbdy).^(-1),fIJ,cIJ); 
end
% cddx : x derivative over corner [cIJ-by-fIJ]      
% cddy : y derivative over corner [cIJ-by-eIJ]     
dx = fx(cconnect(:,2))-fx(cconnect(:,1)); dx(clbdy) = dx(clbdy)+(xr-xl);  dx(crbdy) = dx(crbdy)+(xr-xl); dx(dx==0) = inf;         %[ correct left and right boundaries ]
cddx = sparse([cs; cs],[cconnect(:,1); cconnect(:,2)],[-dx.^(-1); dx.^(-1)],cIJ,fIJ);
dy = ey(cconnect(:,4))-ey(cconnect(:,3)); dy(ctopbdy) = dy(ctopbdy)+(yt-yb); dy(cbotbdy) = dy(cbotbdy)+(yt-yb); dy(dy==0) = inf;  %[ correct top and bottom boundaries ]
cddy = sparse([cs; cs],[cconnect(:,3); cconnect(:,4)],[-dy.^(-1); dy.^(-1)],cIJ,eIJ);
% cdds : s derivative over corner [cIJ-by-nIJ]      
% cddr : s derivative over corner [cIJ-by-nIJ] 
ds = (dx.^2+dy.^2).^(1/2); dr = (dx.^2+dy.^2).^(1/2);
cdds = sparse([cs; cs],[fconnect(cconnect(:,1),1); fconnect(cconnect(:,2),2)],[-ds.^(-1); ds.^(-1)],cIJ,nIJ);
cddr = sparse([cs; cs],[fconnect(cconnect(:,1),2); fconnect(cconnect(:,2),1)],[-dr.^(-1); dr.^(-1)],cIJ,nIJ);
% ndds : corner divergence over node [nIJ-by-cIJ]
% nddr : corner divergence over node [nIJ-by-cIJ]
ndds = sparse([ns; ns],[econnect(nconnect(:,1),3); econnect(nconnect(:,2),4)],[-Dx.^(-1).*Dy.^(-1); Dx.^(-1).*Dy.^(-1)],nIJ,cIJ);
nddr = sparse([ns; ns],[econnect(nconnect(:,1),4); econnect(nconnect(:,2),3)],[-Dx.^(-1).*Dy.^(-1); Dx.^(-1).*Dy.^(-1)],nIJ,cIJ);
% for reflective bcs, s and r must be switched on boundary : these definitions assume ndds and nddr will always be used together
if oo.yreflective, 
for i = 1:length(ntopbdy),
    ndds(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),3); econnect(nconnect(ntopbdy(i),2),4)]) = [0*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1); 2*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1)];
    nddr(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),4); econnect(nconnect(ntopbdy(i),2),3)]) = [-2*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1); 0*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1)];
end
for i = 1:length(nbotbdy),
    ndds(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),3); econnect(nconnect(nbotbdy(i),2),4)]) = [-2*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1); 0*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1)];
    nddr(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),4); econnect(nconnect(nbotbdy(i),2),3)]) = [0*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1); 2*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1)];
end
end
if oo.xreflective, 
for i = 1:length(nlbdy),
    ndds(nlbdy(i),[econnect(nconnect(nlbdy(i),1),3); econnect(nconnect(nlbdy(i),2),4)]) = [0*Dx(nlbdy(i)).^(-1).*Dy(nlbdy(i)).^(-1); 2*Dx(nlbdy(i)).^(-1).*Dy(nlbdy(i)).^(-1)];
    nddr(nlbdy(i),[econnect(nconnect(nlbdy(i),1),4); econnect(nconnect(nlbdy(i),2),3)]) = [0*Dx(nlbdy(i)).^(-1).*Dy(nlbdy(i)).^(-1); 2*Dx(nlbdy(i)).^(-1).*Dy(nlbdy(i)).^(-1)];
end
for i = 1:length(nrbdy),
    ndds(nrbdy(i),[econnect(nconnect(nrbdy(i),1),3); econnect(nconnect(nrbdy(i),2),4)]) = [-2*Dx(nrbdy(i)).^(-1).*Dy(nrbdy(i)).^(-1); 0*Dx(nrbdy(i)).^(-1).*Dy(nrbdy(i)).^(-1)];
    nddr(nrbdy(i),[econnect(nconnect(nrbdy(i),1),4); econnect(nconnect(nrbdy(i),2),3)]) = [-2*Dx(nrbdy(i)).^(-1).*Dy(nrbdy(i)).^(-1); 0*Dx(nrbdy(i)).^(-1).*Dy(nrbdy(i)).^(-1)];
end
end
if oo.yreflective && oo.xreflective, 
% note the two entries given by econnect here should have the same index
% for the corner nodes, so only the second of the entries assigned here
% will be used
for i = 1,
    ndds(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),3); econnect(nconnect(ntopbdy(i),2),4)]) = [0; 4*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1)];
    nddr(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),4); econnect(nconnect(ntopbdy(i),2),3)]) = [0; 0*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1)];
end
for i = length(ntopbdy),
    ndds(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),3); econnect(nconnect(ntopbdy(i),2),4)]) = [00*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1)];
    nddr(ntopbdy(i),[econnect(nconnect(ntopbdy(i),1),4); econnect(nconnect(ntopbdy(i),2),3)]) = [0; -4*Dx(ntopbdy(i)).^(-1).*Dy(ntopbdy(i)).^(-1);]; 
end
for i = 1,
    ndds(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),3); econnect(nconnect(nbotbdy(i),2),4)]) = [0; 0*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1)];
    nddr(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),4); econnect(nconnect(nbotbdy(i),2),3)]) = [0; 4*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1)];
end
for i = length(nbotbdy),
    ndds(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),3); econnect(nconnect(nbotbdy(i),2),4)]) = [0; -4*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1)];
    nddr(nbotbdy(i),[econnect(nconnect(nbotbdy(i),1),4); econnect(nconnect(nbotbdy(i),2),3)]) = [0; 0*Dx(nbotbdy(i)).^(-1).*Dy(nbotbdy(i)).^(-1)];
end
end


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
gg.K = K;
gg.bdy = bdy;
gg.Dx = Dx;
gg.Dy = Dy;
gg.Ds = (Dx.^2+Dy.^2).^(1/2);
gg.Dr = gg.Ds;
gg.nx = nx;
gg.ny = ny;
gg.ex = ex;
gg.ey = ey;
gg.fx = fx;
gg.fy = fy;
gg.cx = cx;
gg.cy = cy;
gg.zz = zz;

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

disp('Done');
end
