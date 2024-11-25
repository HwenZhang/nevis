function nevis_route()
% nevis_route(elev,quant,gg,oo)
% 9 April 2016 - new routing algorithm

%% setup grid
% oo.xperiodic = 1; oo.yperiodic = 1;
[gg] = nevis_grid(100,150,1,-1,1,-2,1);

% elevation to use for routing
nr = (gg.nx.^2+gg.ny.^2).^(1/2);
elev = nr.*exp(-nr.^2);
elev2 = reshape(elev,gg.nIJ,[]);

% quantity to route
quant = ones(gg.nI,gg.nJ);
quant2 = reshape(quant,gg.nIJ,[]);

% label edge nodes
edge = false(gg.nIJ,1);
edge(gg.n1) = 1;

%% test
pd = struct;
[pd,oo] = nevis_defaults(pd);
[ps,pp] = nevis_nondimension(pd);
load('../gap/topo_140901'); % topography input, cropped
gg = nevis_grid(dd.x(dd.xi_crop)/ps.x,dd.y(dd.yi_crop)/ps.x,oo); ps.x0 = dd.x0; ps.y0 = dd.y0;
b = reshape(dd.b(dd.xi_crop,dd.yi_crop)/ps.z,gg.nIJ,1);
s = reshape(dd.s(dd.xi_crop,dd.yi_crop)/ps.z,gg.nIJ,1);
h = fspecial('gaussian',3,1); s_smooth = filter2(h,reshape(s,gg.nI,gg.nJ)); % smooth surface for surface routing
% mask with minimum thickness and mask
b(~dd.mask(dd.xi_crop,dd.yi_crop)) = NaN; s(~dd.mask(dd.xi_crop,dd.yi_crop)) = NaN; 
H = max(s-b,0);
H(isnan(H)) = 0;
Hmin = 15/ps.z; 
gg = nevis_mask(gg,find(H<Hmin)); 
gg.n1m = gg.n1(H(gg.n1)<100/ps.z); 
% label boundary nodes
gg = nevis_label(gg,gg.n1m); oo.adjust_boundaries = 1;

elev2 = reshape(s_smooth,gg.nIJ,[]);
quant2 = reshape(ones(gg.nIJ,1),gg.nIJ,[]);
quant2(gg.nout) = 0;

% label edge nodes
edge = false(gg.nIJ,1);
[n1,~,~,~,~,~,~,~,~,~] = nevis_shape2(gg,gg.ns); % nodes bordering active nodes
edge(union(n1,gg.nout)) = 1;

%%
% figure(1); clf; imagesc(gg.nx(:,1),gg.ny(1,:),reshape(elev2,gg.nI,gg.nJ)'); colorbar; drawnow; shg; % test plot

% node indices
nis = 1:gg.nIJ; 
ni_edge = nis(edge);

% neigbhouring 4 nodes [left, top, right, bottom]
nconnect1 = [gg.econnect(gg.nconnect(:,1),1) gg.fconnect(gg.nconnect(:,3),1) gg.econnect(gg.nconnect(:,2),2) gg.fconnect(gg.nconnect(:,4),2)];

% neigbhouring 8 nodes [left, top left, top, top right, right, bottom right, bottom, bottom left]
nconnect2 =  [nconnect1(:,1) nconnect1(nconnect1(:,1),2) nconnect1(:,2) nconnect1(nconnect1(:,2),3) nconnect1(:,3) nconnect1(nconnect1(:,3),4) nconnect1(:,4) nconnect1(nconnect1(:,4),1)];

% difference operator in each of 8 compass points [left, top left, top, top right, right, bottom right, bottom, bottom left]
diffs_op = sparse(kron([(0:7) (0:7)]',gg.nIJ*ones(gg.nIJ,1))'+repmat(1:gg.nIJ,1,16),[repmat(1:gg.nIJ,1,8) nconnect2(:)'],[-ones(1,8*gg.nIJ) ones(1,8*gg.nIJ)],8*gg.nIJ,gg.nIJ,16*gg.nIJ);


% neighbouring edges and corners [left edge, top left corner, top edge, top right corner, right edge, bottom right corner, bottom edge, bottom left
%corner]
nconnect3 =  [gg.nconnect(:,1) gg.econnect(gg.nconnect(:,1),3) gg.nconnect(:,3) gg.econnect(gg.nconnect(:,2),3) gg.nconnect(:,2) gg.econnect(gg.nconnect(:,2),4) gg.nconnect(:,4) gg.econnect(gg.nconnect(:,1),4)];

% normal derivative operator in each of 8 compass points [left, top left, top, top right, right, bottom right, bottom, bottom left]
grads_op = [-gg.eddx(nconnect3(:,1),:); -gg.cdds(nconnect3(:,2),:); -gg.fddy(nconnect3(:,3),:); gg.cddr(nconnect3(:,4),:); gg.eddx(nconnect3(:,5),:); gg.cdds(nconnect3(:,6),:); gg.fddy(nconnect3(:,7),:); -gg.cddr(nconnect3(:,8),:)];


% calculate differences
diffs = reshape(diffs_op*elev2,[],8);

% calculate normal derivatives
grads = reshape(grads_op*elev2,[],8);

% % use D4 (only 4 compass points rather than 8)
% nconnect2 = nconnect2(:,1:2:7);
% diffs = diffs(:,1:2:7); 
% grads = grads(:,1:2:7); 

% immediate downstream node
% [diff,down] = min(diffs,[],2); % use differences
[diff,down] = min(grads,[],2); % use gradients
nconnect4 = nconnect2((1:gg.nIJ)'+(down-1)*gg.nIJ);

% label local minima nodes [ where all normal derivatives are positive ]
localmin = (diff>=0) & ~edge;
ni_localmin = nis(localmin);
% display(['Local minima: ',num2str(nis(localmin))]);
display(['Number of local minima: ',num2str(length(ni_localmin))]);

% % label sink nodes
% sink = ones(gg.nIJ,1);
% sink(localmin) = 0;
% sink(edge) = 0;

% connect sinks to themselves
nconnect4(localmin) = nis(localmin);
nconnect4(edge) = nis(edge);

% matrix of immediate downstream connection
A = sparse(nis,nconnect4,ones(1,gg.nIJ));

% % matrix of all downstream connections
% B = A;
% for m = 2:gg.nIJ,
%     C = A^m; [m nnz(C*sink)],
%     if nnz(C*sink)==0, break; 
%     else B = max(B,C);
%     end
% end
% % spy(B) 
% % n_up = sum(B,1); % number of upstream cells
% % n_down = sum(B,2); % number of downstream cells
% % plot(nis(localmin),n_up(localmin),'o',nis(edge),n_up(edge),'o');

% matrix of final downstream connection
C = A^gg.nIJ;

% final downstream connection
nconnect5 = NaN*ones(gg.nIJ,1);
[i,j,~] = find(C);
nconnect5(i) = j;

% downstream node including diversion of local minima [ run through local minima
% and divert to lowest node on border of catchment area ]
nconnect6 = nconnect4; % downstream connections as before
% now edit entries for local minima
disp(['Diverting local minima ...']); tic;
clf;
for i = 1:length(ni_localmin), 
    ni = ni_localmin(i);
    % disp(['Diverting local minimum at node ',num2str(ni),' ...']);
    ns = find(C(:,ni)); % nodes in catchment area of local minimum
%     [~,n2,~,~,~,~] = nevis_shape(gg,ns); % nodes bordering catchment area
    [~,n2,~,~,~,~,~,~,~,~] = nevis_shape2(gg,ns); % nodes bordering catchment area
    [~,tmp] = min(elev2(n2)); ni2 = n2(tmp); % find lowest node on border of catchment area
    nconnect6(ni) = ni2; % divert local minimum to the new node
    C(:,nconnect5(ni2)) = max(C(:,nconnect5(ni2)),C(:,ni)); % add catchment area to that of new node
    hold on; plot(gg.nx(ns),gg.ny(ns),'.',gg.nx(ni),gg.ny(ni),'mo',gg.nx(ni2),gg.ny(ni2),'mx'); drawnow; shg; % test plot
    % disp('Done');
end
disp(['Done [',num2str(toc),' seconds ]']);

% re-label sink nodes (not including local minima)
sink = ones(gg.nIJ,1);
sink(edge) = 0;
sink_mat = sparse(1:gg.nIJ,1:gg.nIJ,sink);

% matrix of immediate downstream connection including diversion
A = sparse(nis,nconnect6,ones(1,gg.nIJ));

% matrix of all downstream connections inclusion diversion
disp('Calculating all connections ...'); tic;
B = A; % each row to indicate every node reachable downstream 
C = A; % each ro to indicate the node reached in one downstream step
for m = 1:log2(gg.nIJ)+1; %2:gg.nIJ,
%     C = A^m; % slower
%     C = C*A; %[m nnz(C*sink)], % increase number of downstream steps in C by one
    C = C*C; [m nnz(C*sink)], % increase number of downstream steps in C exponentially
    if nnz(C*sink)==0, break; 
    else
%         B = max(B,C); % either move according to C, or do any of moves possible in one fewer steps, already included in B
        B = (speye(gg.nIJ)+C)*B;  %either stay put or move according to C, then do any of moves possible in one fewer steps
    end
    if m>=log2(gg.nIJ), ind = (C*sink)>0; plot(gg.nx(ind),gg.ny(ind),'x'); end
end
disp(['Done [',num2str(toc),' seconds ]']);
B = logical(B);
% spy(B) 
% n_up = sum(B,1); % number of upstream cells
% n_down = sum(B,2); % number of downstream cells
% plot(nis(localmin),n_up(localmin),'o',nis(edge),n_up(edge),'o');

% matrix of final downstream connection including diversion
disp('Calculating final connections ...'); tic;
C = A^gg.nIJ;
disp(['Done [',num2str(toc),' seconds ]']);

% final downstream connection including diversion
nconnect7 = NaN*ones(gg.nIJ,1);
[i,j,~] = find(C);
nconnect7(i) = j;

% maxtrix of upstream nodes
upstream_mat = B';

% matrix of downtream nodes
downstream_mat = B;

% upstream area
area = gg.Dx.*gg.Dy;
upstream_area = upstream_mat*area;

%routed quantity
routed_quant = upstream_mat*quant2;

%% plot
figure(2); clf;
%     zz = log10(upstream_area); cax = [-3 1];
    zz = routed_quant; cax = [0 100];
    imagesc(gg.nx(:,1),gg.ny(1,:),reshape(zz,gg.nI,gg.nJ)'); 
    set(gca,'YDir','normal');
    xlabel('x');
    ylabel('y');
    colorbar; caxis(cax);
    load cmapq2; colormap(cmap);

    ni = ni_localmin; ni2 = nconnect6(ni);     
    hold on; plot(gg.nx(ni),gg.ny(ni),'ro',gg.nx(ni2),gg.ny(ni2),'go'); % add local minima nodes

    ni = ni_edge;
    hold on; plot(gg.nx(ni),gg.ny(ni),'mo');     % add edge nodes
    
    % add some random routing paths  
    for i = 1:20,
        ni = ceil(rand(1)*gg.nIJ);
        ns = find(downstream_mat(ni,:));
        hold all; plot(gg.nx(ns),gg.ny(ns),'.');
    end

gg.nIJ,
sum(routed_quant(ni_edge)),
nnz(~edge)+length(n1), 
sum(routed_quant(n1)),

end