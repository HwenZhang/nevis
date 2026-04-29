function [z_pts,ws_pts,ts] = nevis_surface_motion(x_pts,y_pts,tis,fn)
% calculate vertical surface motion as measured by monitoring stations at 
% points x_pts,y_pts using ouput from nevis model at time indices tis.
% inputs
%   x_pts,y_pts     dimensionless coordinates of monitoring station
%   tis             indices of timesteps
%   fn              filename of saved solution
% outputs
%   z_pts           vertical elevation of monitoring station, relative to initial
%   ws_pts          vertical velocity at monitoring station
%   ts              timesteps
% z_pts is scaled with [u]*[z]/[x]*[t], ws_pts is scaled with [u]*[z]/[x]
%   
% This does an explicit integration of the vertical surface velocity ws_pts
% calculated using nevis_surface_velocity() using data at saved timesteps. 
% The calculation could be made more accurate by using all timesteps (of 
% which there are typically more than those that are saved), eg by doing 
% this calculation on the fly and storing ws_pts.
%
% Note this does not include additional uplift and subsidence due to 
% hydrological processes at the bed
%
% 7 Apr 2026

% load parameter and grid data
load(fn,'aa','gg','pp','oo','tt'); 

% if y_pts is empty, interpret input as node indices, otherwise find node
% indices closes to the given points
if isempty(y_pts)
    ni_pts = x_pts;
else
    nx = gg.nx; ny = gg.ny;
    ni_pts = NaN*x_pts;
    for i = 1:length(ni_pts)
        [~,tmp] = min((reshape(nx,[],1)-x_pts(i)).^2+(reshape(ny,[],1)-y_pts(i)).^2);
        ni_pts(i) = tmp;
    end
end

% if no timesteps are given, use length of tt
if isempty(tis)
    tis = 1:length(tt);
end

% load saved solutions and calculate surface velocity
ts = zeros(1,length(tis));
ws_pts = zeros(length(x_pts),1)*ts.^0; % initialise
for i_t = 1:length(tis)
    if ~exist([fn,'/',num2str(tis(i_t),'%04.0f'),'.mat']) 
    % stop if no file exists 
        ws_pts(:,i_t) = NaN; ts(i_t) = NaN; 
        break
    end
    load([fn,'/',num2str(tis(i_t),'%04.0f')]); % load solution
    t = vv.t; u = vv.u; v = vv.v;
    ws = nevis_surface_velocity(aa.H,u,v,aa,pp,gg,oo); % calculate surface vertical velocity
    ws_pts(:,i_t) = ws(ni_pts);
    ts(i_t) = t;
end

% calculate z_pts
z_pts = zeros(length(x_pts),1)*ts.^0; % initialise
for i_t = 2:length(ts)
    dt = ts(i_t)-ts(i_t-1);
    z_pts(:,i_t) = z_pts(:,i_t-1)+ws_pts(i_t).*dt;
end

end


