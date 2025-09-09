
load([oo.dn '/' 'morlighem_for_nevis_140km']); % load Morlighem bedmap (previously collated)
dd = morlighem_for_nevis_140km; dd.skip = 6;
gg = nevis_grid(dd.X_km(1:dd.skip:end,1)/ps.x,dd.Y_km(1,1:dd.skip:end)/ps.x,oo); 
b = reshape(dd.B_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
s = reshape(dd.S_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);

x = dd.X_km(1:dd.skip:end,1)/ps.x;
y = dd.Y_km(1,1:dd.skip:end)/ps.x;

pcolor(xx,yy,reshape(ps.z*s,gg.nI,gg.nJ));
colorbar();
hold on
% Overlay a contour plot on top of the pcolor plot
% the contours are 500 1000 and 1500
contour(xx, yy, reshape(ps.z*s, gg.nI, gg.nJ), [500 1000 1500], 'LineColor', 'k');