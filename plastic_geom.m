function zs = plastic_geom(x,zb,tauc,x_offset)
%   Finds surface elevation zs for bed elevation zb and bed yield stress tauc
%   x scaled with x0, zb,zs scaled with z0, tauc scaled with rho_i*g*z0^2/x0
%   optional input offset is a small distance to offset the margin from the
%   end of vector x

    xx = x(end)-x(end:-1:1);
    zb = zb(end:-1:1);
    
    if nargin > 3 && x_offset>0,
        x1 = x_offset;
    else
        x1 = eps;
    end
    zs1 = zb(1)+(2*tauc*x1)^(1/2);
    [xx,y] = ode15s(@y_xx,xx,zs1);
    
    zs = [y(end:-1:1)];

    function s_xx = y_xx(t,y)
        b = interp1(xx,zb,t);
        s_xx = tauc.*(y - b).^(-1);
    end
end
