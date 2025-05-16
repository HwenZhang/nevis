function [x_tran,y_tran,s_tran,phi_tran] = nevis_phi_transect(x,y,N,phi,gg)
% [x_tran,y_tran,phi_tran] = nevis_phi_transect(x1,y1,x2,y2,N,phi,gg)
%
% Interpolate potential along a straight line segment between coordinates x
% and y, with N points
% Inputs:
%   x,y           coordinates of end points of transect (2-by-1 vectors, eg
%   obtained from nevis_get_transect)
%   N             number of points to include on transect
%   phi           pressure on nodes
% Outputs:
%   x_tran,y_tran coordinates of transect
%   s_tran        arc length along transect
%   phi_tran      pressure on transect (interpolated)
%   
% 13 March 2015 IJH

% reshape potential onto a matrix
phi = reshape(phi,gg.nI,gg.nJ);  
xx = reshape(gg.nx,gg.nI,gg.nJ);
yy = reshape(gg.ny,gg.nI,gg.nJ);

% intepolate
x_tran = linspace(x(1),x(2),N);
y_tran = linspace(y(1),y(2),N);
dx = x_tran(2:end)-x_tran(1:end-1);
dy = y_tran(2:end)-y_tran(1:end-1);
s_tran = [0 cumsum((dx.^2+dy.^2).^(1/2))];
phi_tran = interp2(xx',yy',phi',x_tran,y_tran);

end