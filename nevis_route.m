function [Q,low,Qx,Qy,Qs,Qr,nconnect4,nconnect5,drains] = nevis_route(phi,m,gg,oo)
% steady state routing down prescried potential phi with source m
% Inputs
%   phi     potential on nodes
%   m       source on nodes [set to ones(size(phi)) for upstream area]
%   gg      grid structure
%   oo      options structure [optional]
% Outputs
%   Q       cumulative discharge on each node from upstream
%   low     indicator function for local minima
%   Qx,Qy,Qs,Qr     cumulative discharge on each edge/corner from upstream
%   nconnect4       list of most downstream nodes for each node
%   nconnect5       list of lows for each node
%   drains       matrix assigning each node to a low (rows correspond to nodes indicated in low)
%
% 22 Sept 2014: (D8 algotithm?) 
%   26 Oct: added nconnect5 and drains
%   18 Nov: edit to allow flow over lows

%% options
if nargin<4, oo = struct; end
if ~isfield(oo,'include_drains'), oo.include_drains = 0; end % calculate matrix assigning each node to a low
if ~isfield(oo,'prescribed_drains'), oo.prescribed_drains = 0; end % use prescribed sink locations rather than claculating lows [ include list of lows in oo ] [ this method doesn't work ]
if ~isfield(oo,'lows'), oo.lows = []; end

if oo.prescribed_drains, lows = oo.lows; end

disp('nevis_route: Routing ...');

%neigbhouring 4 nodes [left, top, right, bottom]
nconnect1 = [gg.econnect(gg.nconnect(:,1),1) gg.fconnect(gg.nconnect(:,3),1) gg.econnect(gg.nconnect(:,2),2) gg.fconnect(gg.nconnect(:,4),2)];

%neigbhouring 8 nodes [left, top left, top, top right, right, bottom right, bottom, bottom left]
nconnect2 =  [nconnect1(:,1) nconnect1(nconnect1(:,1),2) nconnect1(:,2) nconnect1(nconnect1(:,2),3) nconnect1(:,3) nconnect1(nconnect1(:,3),4) nconnect1(:,4) nconnect1(nconnect1(:,4),1)];
clear nconnect1;

%neighbouring edges and corners [left edge, top left corner, top edge, top right corner, right edge, bottom right corner, bottom edge, bottom left
%corner]
nconnect3 =  [gg.nconnect(:,1) gg.econnect(gg.nconnect(:,1),3) gg.nconnect(:,3) gg.econnect(gg.nconnect(:,2),3) gg.nconnect(:,2) gg.econnect(gg.nconnect(:,2),4) gg.nconnect(:,4) gg.econnect(gg.nconnect(:,1),4)];

%gradients in each of 8 compass points
Psix = gg.eddx*phi;
Psiy = gg.fddy*phi;
Psis = gg.cdds*phi;
Psir = gg.cddr*phi;
Psi = [-Psix(nconnect3(:,1)) -Psis(nconnect3(:,2)) -Psiy(nconnect3(:,3)) Psir(nconnect3(:,4)) Psix(nconnect3(:,5)) Psis(nconnect3(:,6)) Psiy(nconnect3(:,7)) -Psir(nconnect3(:,8))];
clear Psix Psiy Psis Psir

%downstream nodes [ may be overridden below ]
[~,down] = min(Psi,[],2);
nconnect4 = zeros(gg.nIJ,1);
for ni = 1:gg.nIJ
nconnect4(ni) = nconnect2(ni,down(ni));
end
% clear nconnect2 down;

%work through nodes and find final downstream node
nconnect5 = nconnect4;
drains = logical(speye(gg.nIJ,gg.nIJ));
low = false(gg.nIJ,1);
if oo.prescribed_drains, low(lows) = 1; end
[~,nis] = sort(phi(gg.ns),1,'descend'); nis2 = nis;
done = false(gg.nIJ,1);
i = 1;
while i<=length(nis),
%     ni1 = gg.ns(nis(i)); ni2 = nconnect4(ni1); % downstream node [ may be an already passed node, replace with next 4 lines 18 Nov? ]
    ni1 = gg.ns(nis(i)); 
    if ~done(ni1),
%     Psi(ni1,done(nconnect2(ni1,:))) = inf; % effecitvely elimiate already passed nodes from consideration
    [~,tmp] = min(Psi(ni1,:));
    ni2 = nconnect2(ni1,tmp); % most downstream node
%     [down(ni1) tmp],
    down(ni1) = tmp; nconnect4(ni1) = ni2; % downstream node not including already labelled upstream node 
   
    
%     nevis_plot_grid(gg,[ni1,ni2]); shg;
%     phi([ni1,ni2]),
%     [ni1,ni2],
    
    %assign upstream nodes and label lows % [ done test in line below may not be needed any more with new defn of ni2 18 Nov ]
%     if (phi(ni2)<=phi(ni1) && ~done(ni2)) || (oo.prescribed_drains && ~low(ni1)) % if prescribed_drains, drain into downstream node even if it's not actually downstream [ may be problem if phi(ni2) and phi(ni1) are equal as water may get stuck ]
    if ~isnan(phi(ni2))
    if ( (phi(ni2)<=phi(ni1)) || (oo.prescribed_drains && ~low(ni1)) ) 
    elseif (oo.prescribed_drains && ~low(ni1))
        phi(ni1) = phi(ni2); % raise the level to that of the lowest neighour, and transfer
    else
        low(ni1) = 1; % label node as a local minimum
    end
    end
    if oo.include_drains,
        drains(ni2,:) = drains(ni2,:)+drains(ni1,:); drains(ni1,:) = 0; % assign node indices of upstream nodes into downstream node
        drains(ni2,ni1) = 1; % add current node to downstream nodes
        nconnect5(logical(drains(ni2,:))) = ni2; % assign downstream node to all upstream nodes
    end

%     [ni1 ni2 done(ni1) done(ni2)],
    if done(ni2) && ~done(ni1), done(ni2)=0; i = find(gg.ns(nis)==ni2); else i = i+1; end % go back to ni2 if it's already been done
    done(ni1) = 1;
    else
        i = i+1;
    end
    nis2(i) = ni1;
end

%work through nodes summing source
Q = m.*gg.Dx.*gg.Dy;
Qx = zeros(gg.eIJ,1);
Qy = zeros(gg.fIJ,1);
Qs = zeros(gg.cIJ,1);
Qr = zeros(gg.cIJ,1);
[~,nis] = sort(phi(gg.ns),1,'descend');
done = false(gg.nIJ,1);
i = 1;
while i<=length(nis)
    ni1 = gg.ns(nis(i)); 
    ni1 = nis2(i); % nis2 is new ordering as determined from above
%     if ~done(ni1),
    ni2 = nconnect4(ni1); %downstream node
%     nevis_plot_grid(gg,[ni1,ni2]); shg;
%     phi([ni1,ni2]),
%     [ni1,ni2],

    %discharge on nodes
    if ~low(ni1)
        Q(ni2) = Q(ni2) + Q(ni1); % add discharge to downstream node
    end
    
    if ~low(ni1)
    %discharge on edges/corners
    if down(ni1)==1 || down(ni1)==5, Qx(nconnect3(ni1,down(ni1))) = Q(ni1);
    elseif down(ni1)==3 || down(ni1)==7, Qy(nconnect3(ni1,down(ni1))) = Q(ni1);
    elseif down(ni1)==2 || down(ni1)==6, Qs(nconnect3(ni1,down(ni1))) = Q(ni1);
    elseif down(ni1)==4 || down(ni1)==8, Qr(nconnect3(ni1,down(ni1))) = Q(ni1);
    end
    end
    
%     if done(ni2) && ~done(ni1), done(ni2)=0; i = find(gg.ns(nis)==ni2); else i = i+1; end % go back to ni2 if it's already been done
%     done(ni1) = 1;
%     else
%         i = i+1;
%     end
end
drains = drains(low,:);

disp('nevis_route: Done');

end


    

