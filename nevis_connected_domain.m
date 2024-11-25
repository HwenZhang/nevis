function ns = nevis_connected_domain(gg,x,y)
% find the connected domain closes to the point x,y
% from the points in gg.ns
% Input:
%   gg  structure with current nodes gg.ns
%   x,y coordinates of point in desired domain
% Output:
%   ns  indices of connected nodes
%
% 13 October 2015

[~,tmp] = min((x-gg.nx(gg.ns)).^2+(y-gg.ny(gg.ns)).^2); tmp = gg.ns(tmp); % find starting node
in = false(gg.nIJ,1); in(tmp) = 1; % indicator matrix of connected nodes
con = gg.nmeanx(:,gg.ein)*gg.emean(gg.ein,:) + gg.nmeany(:,gg.fin)*gg.fmean(gg.fin,:); % average across nodes

num_pts = NaN*ones(10^5,1);
num_pts(1) = 1;
for i = 2:length(num_pts),
    tmp = con*in;
    in = tmp>0;
    num_pts(i) = nnz(in);
    disp(['Number of connected points: ',num2str(num_pts(i))]);
    if num_pts(i)<=num_pts(i-1), break; end
end

ns = find(in);

end