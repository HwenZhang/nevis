function out = Vb_reg(Vb,pp,opts)   
    % if opts.hyperbolic_reg, out = tanh(pp.V_b_sigma*Vb); end
    % out = 1./(exp(pp.V_b_sigma.*(-Vb+pp.V_b_reg))+1);
    out = tanh(pp.V_b_sigma*(Vb-pp.V_b_reg));
end