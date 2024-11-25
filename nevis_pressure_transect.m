function [x_tran,y_tran,s_tran,phi_tran,pw_tran,pi_tran,k_tran] = nevis_pressure_transect(x,y,N,phi,phi_0,phi_a,gg)
% [x_tran,y_tran,s_tran,phi_tran,pw_tran,pi_tran,k_tran] = nevis_pressure_transect(x,y,N,phi,phi_0,phi_a,gg)
%
% Interpolate pressure along a straight line segment between coordinates x
% and y, with N points
% Inputs:
%   x,y           coordinates of end points of transect (2-by-1 vectors, eg
%   obtained from nevis_get_transect)
%   N             number of points to include on transect
%   phi           potential on nodes
%   
% Outputs:
%   x_tran,y_tran coordinates of transect
%   s_tran        arc length along transect
%   phi_tran      potential on transect (interpolated)
%   pw_tran       water pressure on transect (interploated)
%   pi_tran       ice pressure on transect (interploated)
%   k_tran        pressure ratio pw_tran/pi_tran 
%   
% 1 May 2015 IJH - edited from nevis_phi_transect

% reshape pressure onto a matrix
pi = reshape(phi_0-phi_a,gg.nI,gg.nJ);  
pw = reshape(phi-phi_a,gg.nI,gg.nJ);  
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
pw_tran = interp2(xx',yy',pw',x_tran,y_tran);
pi_tran = interp2(xx',yy',pi',x_tran,y_tran);
k_tran = pw_tran./pi_tran;

end