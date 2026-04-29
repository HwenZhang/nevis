function ws = nevis_surface_velocity(H,u,v,aa,pp,gg,oo)
% to calculate vertical surface velocity from nevis SSA solution found
% using nevis_velocity
% inputs
%   H       ice thickness on nodes [nIJ-by-1]
%   u       x component of velocity on edges [eIJ-by-1]
%   v       y component of velocity on edges [fIJ-by-1]
%   aa      prescribed fields and boundary terms [struct]
%   pp      parameters [struct]
%   gg      grid and discretisation operators [struct]
%   oo      options [struct]
% outputs 
%   ws       z component of surface velocity on nodes [nIJ-by-1]
%
% ws is scaled with [u][z]/[x]
% the time integral of this quantity, plus the blister thickness h_b, 
% should give vertical position measured by a GNSS station fixed to the ice
% surface
% 7 Apr 2026

% operators
nmeanx = gg.nmeanx;
nmeany = gg.nmeany;
nddx = gg.nddx;
nddy = gg.nddy;
eddx = gg.eddx;
fddy = gg.fddy;

% topography
b = aa.b;       % bed elevation on nodes [nIJ-by-1]
s = b+H;        % surface elevation [nIJ-by-1]

% vertical velocity at bed (ignores melting and vertical motion associated
% with drainage system expanding/contracting)
wb = nmeanx*(u.*(eddx*b))+nmeany*(v.*(fddy*b));

% vertical velocity at surface
ws = wb + H.*(nddx*u + nddy*v);

end
