function out = k_b(hs,pp,opts)
    % if constant k, k=k0
    % else k = C h_s^3
    % dimensionless permeability k as a function of dimensionless sheet thickness h
    temp = pp.c49*(hs).^3; % dimensional value/scale 
    out = opts.constant_k + (1-opts.constant_k)*temp;
    % out = ones(size(hs));
end