function [C1_hat_dim, C2_hat_dim] = nevis_inv_partition(C_total, partition_ratio, u_obs_noisy, v_obs_noisy, N_current, aa, pp, gg, oo, ps)
    % Partition the total slipperiness coefficient C_total into two components C1 and C2 based on the specified partition ratio.
    % C_total: the total slipperiness coefficient (nondimensional)
    % partition_ratio: the ratio of C1 to C_total (between 0 and 1), e.g., 0.25 means that 25% of the basal shear stress is from C1 and can be reduced by lowering N, while 75% is from C2 and is unaffected by N.
    
    % Ensure that the partition ratio is between 0 and 1
    if partition_ratio < 0 || partition_ratio > 1
        error('Partition ratio must be between 0 and 1');
    end
    
    aa.C = C_total;
    [u_inv, v_inv] = nevis_velocity(aa.H, u_obs_noisy, v_obs_noisy, N_current, aa, pp, gg, oo);
    U_slide = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);
    % 2) compute basal shear stress from velocity and current N
    tau_b = slide_fun_local(U_slide, N_current, C_total, pp.mu*ones(size(aa.b)),pp).*U_slide;  % basal shear stress magnitude from sliding law
    % 3) partition basal shear stress into C and C2 contributions (25% for C and 75% for C2 for simplicity)
    % define a ratio of mu N to tau_b
    % ratio = pp.mu * N_current ./ (tau_b + pp.N_slide_reg);  % add reg to avoid div by zero
    % taub_c = tau_b.*(ratio - 1)./ratio;
    % taub_c2 = tau_b - taub_c;
    taub_c = partition_ratio * tau_b;
    taub_c2 = tau_b - taub_c;

    C1_hat = (U_slide./(taub_c.^pp.n_slide) - U_slide./(pp.mu.*N_current).^pp.n_slide).^(-1/pp.n_slide);
    C1_hat = max(C1_hat, 1e-16);
    C2_hat = taub_c2./(U_slide.^(1/pp.n_slide));

    fprintf('Partitioning C into C1 and C2 with partition ratio = %.2f\n', partition_ratio);
    fprintf(' Max C contribution = %.3e, max C2 contribution = %.3e\n', max(C1_hat), max(C2_hat));
    fprintf(' Min C contribution = %.3e, min C2 contribution = %.3e\n', min(C1_hat), min(C2_hat));

    %% Dimensionalize C1 and C2 for better interpretability (units of Pa s/m)
    C1_hat_dim = C1_hat * (ps.tau / ps.u_b^(1/pp.n_slide));
    C2_hat_dim = C2_hat * (ps.tau / ps.u_b^(1/pp.n_slide));
    fprintf('  Dimensional C1 max = %.3e Pa s/m, C2 max = %.3e Pa s/m\n', max(C1_hat_dim), max(C2_hat_dim));
    fprintf('  Dimensional C1 min = %.3e Pa s/m, C2 min = %.3e Pa s/m\n', min(C1_hat_dim), min(C2_hat_dim));

end

%% ============================================================
%  Local sliding law (same as in nevis_velocity)
%  ============================================================
function tau_b_over_Ub = slide_fun_local(Ub, N, C, mu, pp)
    % MUST match nevis_velocity slide_fun exactly (same regularisation)
    n = pp.n_slide;
    Nr = max(N, pp.N_slide_reg);  % max-based (matches nevis_velocity)
    Ur = max(Ub, pp.Ub_reg);      % max-based (matches nevis_velocity)
    tau_b_over_Ub = Nr .* Ur.^(1/n-1) .* ...
        (mu.^(-n).*Ur + C.^(-n).*Nr.^n).^(-1/n);
end
