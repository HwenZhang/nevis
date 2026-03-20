%% NEVIS C-field inversion script
% Invert for the spatial basal friction coefficient C(x,y) under a
% pure-Weertman sliding law (mu = Inf, C2 = 0):
%
%   tau_b = C * Ub^(1/n)
%
% using adjoint-based L-BFGS minimisation with Tikhonov continuation.
% So this is essentially an inversion for C2 in the sliding law, but we call it C-field inversion for simplicity.
%
% Author: Hanwen Zhang
% Date: Feb 2026

clc; clear; close all; format compact;

%% ============================================================
%  1. Load spinup results and set up NEVIS framework
%  ============================================================
oo.root    = './';
oo.code    = './src';
oo.results = 'results';
oo.dataset = 'nevis_regional';
oo.casename = 'n2d_region_ice_inversion_test';
oo.fn = ['/', oo.casename];
oo.rn = [oo.root, oo.results, oo.fn];
oo.in = [oo.root, 'data/', oo.casename];
oo.dn = [oo.root, 'data/', oo.dataset, '/'];
addpath(oo.code);

load([oo.in, oo.fn], 'pp','pd','ps','gg','aa','oo');
[pd, ps, pp, oo] = nevis_update_parameters_ice(pd, ps, pp, oo);

% Pure Weertman: no Coulomb cap, no C2 term
pd.mu_s = Inf;
pp.mu   = Inf;
pp.C2   = 0;

% Picard iterator
oo.iter_max             = 100;
oo.tol_vel              = 1e-6;
oo.verb                 = 0;
oo.display_norms        = 0;
oo.boundary_method      = 'stress_l_vel_tbl';
oo.mask_boundary_method = 'stress_free';

% Regularisation
pp.eps_reg     = 1e-1;
pp.Ub_reg      = 1e-16;
pp.N_slide_reg = 1e-16;
pp.taud_reg    = 1e-16;

%% ============================================================
%  2. Load observed velocity and effective pressure
%  ============================================================
gg = nevis_label_ice_test(gg, oo);
load([oo.dn 'measures_for_nevis_140km.mat']);
dd    = measures_for_nevis_140km;
u_obs = dd.u_obs_dim / (ps.u_b * pd.ty);
v_obs = dd.v_obs_dim / (ps.u_b * pd.ty);
u_obs = gg.emean2 * u_obs(:);
v_obs = gg.fmean2 * v_obs(:);
u_obs(gg.eout2) = NaN;
v_obs(gg.fout2) = NaN;

% Observation mask: interior edges only, exclude Dirichlet boundaries
obs_mask_e = zeros(gg.eIJ, 1);  obs_mask_e(gg.ein2) = 1;  obs_mask_e(gg.ebdy2) = 0;
obs_mask_f = zeros(gg.fIJ, 1);  obs_mask_f(gg.fin2) = 1;  obs_mask_f(gg.fbdy2) = 0;

% Effective pressure from spinup (year-end snapshot)
init_cond = load(['./data/' oo.casename '/' '0365.mat']);
vv = init_cond.vv;
N = aa.phi_0 - vv.phi;
N(isnan(N)) = 0;

%% ============================================================
%  3. Inversion settings
%  ============================================================
opts_inv.u0_reg = 1e-1;

% Continuation schedule (coarse → fine)
opts_inv.alpha_schedule    = [1e-3, 1e-4, 1e-9, 1e-12];
opts_inv.gamma_schedule    = [1e-6, 1e-7, 1e-8, 1e-9];
opts_inv.max_iter_schedule = [25, 25, 25, 25];
opts_inv.max_iter_total    = 200;
opts_inv.J_tol             = 1e-3;
opts_inv.dJ_tol            = 1e-6;
opts_inv.alpha             = opts_inv.alpha_schedule(1);
opts_inv.gamma             = opts_inv.gamma_schedule(1);
opts_inv.verbose           = false;
opts_inv.speed_misfit      = true;

history.iter  = [];
history.J     = [];
history.normg = [];

% Uniform prior: c = log(C) = 0 → C = 1
c_prior = zeros(gg.nIJ, 1);

% Adaptive Tikhonov operator: down-weight smooth regularisation in shear zones
dudx_e    = gg.eddx * (gg.nmeanx2(:,gg.es2) * u_obs(gg.es2));
dvdy_f    = gg.fddy * (gg.nmeany2(:,gg.fs2) * v_obs(gg.fs2));
grad_mag_e = abs(dudx_e(gg.ein2));
grad_mag_f = abs(dvdy_f(gg.fin2));
beta_adaptive = 0.1;
w_e = 1 ./ (1 + beta_adaptive * grad_mag_e / mean(grad_mag_e));
w_f = 1 ./ (1 + beta_adaptive * grad_mag_f / mean(grad_mag_f));
Lx = spdiags(w_e, 0, length(gg.ein2), length(gg.ein2)) * gg.eddx(gg.ein2, :);
Ly = spdiags(w_f, 0, length(gg.fin2), length(gg.fin2)) * gg.fddy(gg.fin2, :);
L  = [Lx; Ly];

%% ============================================================
%  4. Run inversion (L-BFGS with Tikhonov continuation)
%  ============================================================
c_stage    = c_prior;
iter_total = 0;

for k = 1:min(numel(opts_inv.alpha_schedule), numel(opts_inv.gamma_schedule))
    opts_inv.alpha = opts_inv.alpha_schedule(k);
    opts_inv.gamma = opts_inv.gamma_schedule(k);
    max_iter_k     = opts_inv.max_iter_schedule(min(k, end));

    fprintf('\n--- Stage %d/%d: alpha=%.2e, gamma=%.2e, max_iter=%d ---\n', k, ...
        numel(opts_inv.alpha_schedule), opts_inv.alpha, opts_inv.gamma, max_iter_k);

    obj_fun = @(c) objective_and_grad(c, u_obs, v_obs, ...
        obs_mask_e, obs_mask_f, N, aa, pp, ps, pd, gg, oo, L, c_prior, opts_inv);

    opts_lbfgs = optimoptions('fminunc', ...
        'Algorithm',                'quasi-newton', ...
        'HessianApproximation',     {'lbfgs', 20}, ...
        'SpecifyObjectiveGradient', true, ...
        'Display',                  'iter', ...
        'MaxIterations',            max_iter_k, ...
        'MaxFunctionEvaluations',   max_iter_k * 20, ...
        'OptimalityTolerance',      1e-6, ...
        'StepTolerance',            1e-12, ...
        'ObjectiveLimit',           1e-4);

    [c_stage, J_stage, ~, output_stage] = fminunc(obj_fun, c_stage, opts_lbfgs);

    iter_total = iter_total + output_stage.iterations;
    history.iter(end+1,1)  = iter_total;
    history.J(end+1,1)     = J_stage;
    history.normg(end+1,1) = output_stage.firstorderopt;

    fprintf('Stage %d: J=%.4e, grad=%.2e, iters=%d (total=%d)\n', ...
        k, J_stage, output_stage.firstorderopt, output_stage.iterations, iter_total);

    if J_stage < opts_inv.J_tol,                      fprintf('Early stop: J < J_tol\n');           break; end
    if output_stage.firstorderopt < opts_inv.dJ_tol,  fprintf('Early stop: grad < dJ_tol\n');       break; end
    if iter_total >= opts_inv.max_iter_total,          fprintf('Early stop: max_iter_total\n');      break; end
end

c_hat     = c_stage;
C_hat     = exp(c_hat);
C_hat_dim = C_hat * (ps.tau / ps.u_b^(1/pp.n_slide));

save('./data/C_inversion_runoff_30_C2_results.mat', ...
     'C_hat_dim', 'c_hat', 'history', 'opts_inv', 'J_stage');
fprintf('\nSaved results. Final J = %.4e\n', J_stage);

%% ============================================================
%  5. Plot results
%  ============================================================
aa_inv   = aa;  aa_inv.C = C_hat;
u_inv    = u_obs;  v_inv = v_obs;
u_inv(isnan(u_inv) & ~ismember((1:gg.eIJ)', gg.eout2)) = 0;
v_inv(isnan(v_inv) & ~ismember((1:gg.fIJ)', gg.fout2)) = 0;
[u_inv, v_inv] = nevis_velocity(aa.H, u_inv, v_inv, N, aa_inv, pp, gg, oo);

Un_obs     = sqrt((gg.nmeanx2(:,gg.es2)*u_obs(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_obs(gg.fs2)).^2);
Un_inv     = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);
Un_obs_dim = Un_obs * ps.u_b * pd.ty;
Un_inv_dim = Un_inv * ps.u_b * pd.ty;

figure('Position', [100 100 1400 380]);
subplot(1,4,1);
z = reshape(log10(C_hat_dim), gg.nI, gg.nJ);  z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('Inverted log_{10}(C)'); axis equal tight;

subplot(1,4,2);
z = reshape(Un_obs_dim, gg.nI, gg.nJ);  z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar; clim([0 1000]);
title('Observed speed [m/yr]'); axis equal tight;

subplot(1,4,3);
z = reshape(Un_inv_dim, gg.nI, gg.nJ);  z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar; clim([0 1000]);
title('Modeled speed [m/yr]'); axis equal tight;

subplot(1,4,4);
z = reshape(Un_inv_dim - Un_obs_dim, gg.nI, gg.nJ);  z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar; clim([-100 100]);
colormap(gca, interp1([0 0.5 1],[0 0 1; 1 1 1; 1 0 0],linspace(0,1,256)));
title('Speed misfit [m/yr]'); axis equal tight;

figure('Position', [100 540 800 330]);
subplot(1,2,1);
semilogy(history.iter, history.J, 'o-', 'LineWidth', 1.5);
xlabel('Total iterations'); ylabel('J'); grid on; title('Cost function');
subplot(1,2,2);
semilogy(history.iter, history.normg, 's-', 'LineWidth', 1.5);
xlabel('Total iterations'); ylabel('firstorderopt'); grid on; title('Optimality');
drawnow;

nin2 = gg.nin2;
RMSE = sqrt(mean((Un_inv_dim(nin2) - Un_obs_dim(nin2)).^2, 'omitnan'));
fprintf('\n=== Inversion summary ===\n');
fprintf('  Final J        : %.4e\n', J_stage);
fprintf('  RMSE speed     : %.1f m/yr\n', RMSE);
fprintf('  c = log(C) range: [%.2f, %.2f]\n', min(c_hat(nin2)), max(c_hat(nin2)));
fprintf('  C_hat_dim range : [%.2e, %.2e]\n', min(C_hat_dim(nin2)), max(C_hat_dim(nin2)));

%% ============================================================
%  Objective function and adjoint gradient
%  ============================================================
function [J, g] = objective_and_grad(c, u_obs, v_obs, ...
    obs_mask_e, obs_mask_f, N, aa, pp, ps, pd, gg, oo, L, c_prior, opts_inv)
    % Evaluate cost function J(c) and gradient dJ/dc using adjoint method
    % c = log(C) on all nIJ nodes
    assert(all(isfinite(u_obs(gg.ebdy2))), 'u_obs has NaN at Dirichlet boundaries');

    C = exp(c);
    aa_loc = aa;
    aa_loc.C = C;
    
    % --- 1. Forward solve ---
    % Use observed velocity as initial guess (warm start).
    % Starting from u=v=0 causes Picard to hit max iterations.
    % Observed velocity is close to the solution, so convergence is fast.
    u = u_obs;  % includes NaN at eout2, values at ebdy2
    v = v_obs;  % includes NaN at fout2, values at fbdy2
    % Ensure interior edges have finite initial guess
    u(isnan(u) & ~ismember((1:gg.eIJ)', gg.eout2)) = 0;
    v(isnan(v) & ~ismember((1:gg.fIJ)', gg.fout2)) = 0;

    [u, v] = nevis_velocity(aa.H, u, v, N, aa_loc, pp, gg, oo);

    % --- Check: boundary values must match observations exactly ---
    % nevis_velocity never modifies ebdy2/fbdy2 edges (they are treated as
    % prescribed BCs via the rhs), so u(ebdy2) should still equal u_obs(ebdy2).
    % A mismatch here means the initial guess was wrong or there is a bug.
    ebdy_err_u = max(abs(u(gg.ebdy2) - u_obs(gg.ebdy2)));
    fbdy_err_v = max(abs(v(gg.fbdy2) - v_obs(gg.fbdy2)));
    if ebdy_err_u > 1e-10 || fbdy_err_v > 1e-10
        warning('objective_and_grad: boundary mismatch! max|u(ebdy)-u_obs(ebdy)|=%.2e, max|v(fbdy)-v_obs(fbdy)|=%.2e', ...
            ebdy_err_u, fbdy_err_v);
        % Force consistency: overwrite with observed BCs
        u(gg.ebdy2) = u_obs(gg.ebdy2);
        v(gg.fbdy2) = v_obs(gg.fbdy2);
    end

    % --- 2. Velocity misfit ---
    ein = gg.ein2; 
    fin = gg.fin2;
    u0 = opts_inv.u0_reg;
    
    if opts_inv.speed_misfit
        % --- Speed-based misfit (on nodes) ---
        % U_mod and U_obs defined on nodes via renormalized mean operators
        nmx = gg.nmeanx2; 
        nmy = gg.nmeany2;
        ux_mod = nmx(:,gg.es2)*u(gg.es2);   
        vy_mod = nmy(:,gg.fs2)*v(gg.fs2);
        ux_obs = nmx(:,gg.es2)*u_obs(gg.es2); 
        vy_obs = nmy(:,gg.fs2)*v_obs(gg.fs2);
        U_mod = sqrt(ux_mod.^2 + vy_mod.^2);       % speed at nodes
        U_obs = sqrt(ux_obs.^2 + vy_obs.^2);       % observed speed at nodes
        % node-level mask: active if ANY adjacent edge is observed
        obs_mask_n = double(nmx(:,ein)*obs_mask_e(ein) + nmy(:,fin)*obs_mask_f(fin) > 0);
        dU = obs_mask_n .* (U_mod - U_obs) ./ (U_obs + u0);
        Jmis = 0.5 * (dU' * dU);
    else
        % --- Component-wise misfit (on edges) ---
        du = obs_mask_e(ein) .* (u(ein) - u_obs(ein)) ./ (abs(u_obs(ein)) + u0);
        dv = obs_mask_f(fin) .* (v(fin) - v_obs(fin)) ./ (abs(v_obs(fin)) + u0);
        Jmis = 0.5 * (du'*du + dv'*dv);
    end
    
    % Diagnostic: print speed statistics (only when verbose)
    if opts_inv.verbose
        fprintf('Jmis = %.2e, Jreg = %.2e, Jdamp = %.2e\n', Jmis, ...
            0.5*opts_inv.alpha*(c'*(L'*L)*c), 0.5*opts_inv.gamma*sum((c-c_prior).^2));
    end

    % --- 3. Regularisation ---
    Jreg = 0.5 * opts_inv.alpha * (c'*(L'*L)*c); % Tikhonov regularisation on spatial gradients
    Jdamp = 0.5 * opts_inv.gamma * sum((c - c_prior).^2); % damping to prior
    
    J = Jmis + Jreg + Jdamp;
    
    if nargout < 2, return; end  % only J requested
    
    % --- 4. Adjoint RHS: dJ/d[u,v] ---
    dJdu = zeros(gg.eIJ,1);
    dJdv = zeros(gg.fIJ,1);
    
    if opts_inv.speed_misfit
        % Speed misfit: dJ/dU * dU/d[u,v]
        % dJ/dU_j = (U_mod_j - U_obs_j) / (U_obs_j + u0)^2 * mask_j
        dJdU = obs_mask_n .* (U_mod - U_obs) ./ (U_obs + u0).^2;  % [nIJ x 1]
        % dU/du_k = nmeanx(j,k) * ux_mod(j) / max(U_mod(j), eps)
        U_mod_safe = max(U_mod, 1e-40);
        w_u = dJdU .* ux_mod ./ U_mod_safe;  % [nIJ x 1]
        w_v = dJdU .* vy_mod ./ U_mod_safe;  % [nIJ x 1]
        dJdu(ein) = nmx(gg.ns2,ein)' * w_u(gg.ns2);  % chain rule: nodes -> x-edges
        dJdv(fin) = nmy(gg.ns2,fin)' * w_v(gg.ns2);  % chain rule: nodes -> y-edges
    else
        % Component misfit: straightforward
        dJdu(ein) = obs_mask_e(ein) .* (u(ein) - u_obs(ein)) ./ (abs(u_obs(ein)) + u0).^2;
        dJdv(fin) = obs_mask_f(fin) .* (v(fin) - v_obs(fin)) ./ (abs(v_obs(fin)) + u0).^2;
    end
    
    % --- 5. Assemble Jacobian with Newton correction for exact adjoint ---
    % The Picard Jacobian A freezes slide_fun(U) w.r.t. U.
    % The EXACT Jacobian includes: d(slide_fun)/dU * dU/d[u,v]
    % Without this Newton correction, the adjoint gradient is approximate
    % and the Taylor test will NOT show 2nd-order convergence.
    ns = gg.ns2;
    nIJ_loc = gg.nIJ;
    ne = length(ein); nf = length(fin);
    
    [A, ~, ops] = assemble_velocity_jacobian(u, v, N, aa_loc, pp, gg, oo);
    % Use renormalized operators from Jacobian assembly — identical to those
    % used inside nevis_velocity/assemble_velocity_jacobian to compute slide_fun.
    % Using raw gg.emean2/nmeanx2 would give wrong dF/dC and wrong U.
    nmeanx    = ops.nmeanx;  % renormalized node←x-edge average [nIJ x eIJ]
    nmeany    = ops.nmeany;  % renormalized node←y-edge average [nIJ x fIJ]
    emean_loc = ops.emean;   % renormalized edge←node average   [eIJ x nIJ]
    fmean_loc = ops.fmean;   % renormalized edge←node average   [fIJ x nIJ]
    
    % Speed at nodes (uses same renormalized operators as assemble_velocity_jacobian)
    U = sqrt((nmeanx(:,gg.es2)*u(gg.es2)).^2 + (nmeany(:,gg.fs2)*v(gg.fs2)).^2);
    
    % Newton correction: d(slide_fun)/dU * dU/d[u,v]
    % F_fric_x(i) = -c61 * [emean(i,:) * S(U)] * u(i)
    % dF_fric_x(i)/du(k) = ... (Picard, in A)
    %   + -c61 * u(i) * sum_j emean(i,j) * dS/dU(j) * dU(j)/du(k)   (Newton)
    % where dU(j)/du(k) = nmeanx(j,k) * ux(j) / U(j)
    dslide_dU = compute_dslide_dU(U, N, C, pp.mu, pp);  % [nIJ x 1]
    ux_n = nmeanx(:, gg.es2) * u(gg.es2);  % x-velocity at nodes
    vy_n = nmeany(:, gg.fs2) * v(gg.fs2);  % y-velocity at nodes
    U_safe = max(U, pp.Ub_reg);  % max-based (matches slide_fun_local)
    
    % Node-level weights for Newton correction
    w_ux = dslide_dU .* ux_n ./ U_safe;  % [nIJ x 1]
    w_vy = dslide_dU .* vy_n ./ U_safe;  % [nIJ x 1]
    
    % Newton correction blocks: dA_{ab} = -c61 * diag(a) * [e/f]mean * diag(w_b) * nmean_b
    dA_uu = -pp.c61 * spdiags(u(ein),0,ne,ne) * emean_loc(ein,:) * spdiags(w_ux,0,nIJ_loc,nIJ_loc) * nmeanx(:,ein);
    dA_uv = -pp.c61 * spdiags(u(ein),0,ne,ne) * emean_loc(ein,:) * spdiags(w_vy,0,nIJ_loc,nIJ_loc) * nmeany(:,fin);
    dA_vu = -pp.c61 * spdiags(v(fin),0,nf,nf) * fmean_loc(fin,:) * spdiags(w_ux,0,nIJ_loc,nIJ_loc) * nmeanx(:,ein);
    dA_vv = -pp.c61 * spdiags(v(fin),0,nf,nf) * fmean_loc(fin,:) * spdiags(w_vy,0,nIJ_loc,nIJ_loc) * nmeany(:,fin);
    
    A_full = A + [dA_uu dA_uv; dA_vu dA_vv];
    
    % Adjoint solve: (dF/d[u,v])^T * lambda = -(dJ/d[u,v])^T
    rhs_adj = -[dJdu(ein); dJdv(fin)];
    lambda = A_full' \ rhs_adj;
    lambda_u = lambda(1:ne);
    lambda_v = lambda(ne+1:end);
    
    % --- 6. Gradient dJ/dC from friction term ---
    % dF/dC_j = -c61 * [e/f]mean(:,j) * dslide_dC(j) * [u/v]
    dslide_dC = compute_dslide_dC(U, N, C, pp.mu, pp);  % [nIJ x 1]
    
    dFu_dC = -pp.c61 * spdiags(u(ein), 0, ne, ne) * ...
             emean_loc(ein,:) * spdiags(dslide_dC, 0, nIJ_loc, nIJ_loc);
    dFv_dC = -pp.c61 * spdiags(v(fin), 0, nf, nf) * ...
             fmean_loc(fin,:) * spdiags(dslide_dC, 0, nIJ_loc, nIJ_loc);
    
    dJdC = (lambda_u' * dFu_dC + lambda_v' * dFv_dC)';  % [nIJ x 1]
    
    % --- 7. Chain rule: dJ/dc = C .* dJ/dC + reg gradient ---
    g = C .* dJdC;
    g = g + opts_inv.alpha * (L'*L) * c;
    g = g + opts_inv.gamma * (c - c_prior);
end

%% ============================================================
%  Assemble velocity Jacobian (extracted from nevis_velocity logic)
%  ============================================================
function [A, rhs, ops] = assemble_velocity_jacobian(u, v, N, aa, pp, gg, oo)
    % Assemble the linear system matrix A and rhs for the SSA velocity solve
    % at the current viscosity (i.e. the last Picard iteration).
    % This is the same system solved in nevis_velocity.m
    % ops: struct containing renormalized mean operators for use in adjoint
    
    nIJ = gg.nIJ; eIJ = gg.eIJ; fIJ = gg.fIJ; cIJ = gg.cIJ;
    ns = gg.ns2; es = gg.es2; fs = gg.fs2; cs = gg.cs2;
    nin = gg.nin2; ein = gg.ein2; fin = gg.fin2; cin = gg.cin2;
    nbdyx = gg.nbdyx2; nbdyy = gg.nbdyy2;
    ebdy = gg.ebdy2; fbdy = gg.fbdy2; cbdy = gg.cbdy2;
    
    % redefine mean operators
    nmeanx = gg.nmeanx; nmeany = gg.nmeany; nmeanc = gg.nmeanc;
    emean = gg.emean; fmean = gg.fmean; cmean = gg.cmean;
    nddx = gg.nddx; nddy = gg.nddy;
    eddx = gg.eddx; eddy = gg.eddy;
    fddx = gg.fddx; fddy = gg.fddy;
    cddx = gg.cddx; cddy = gg.cddy;
    
    temp = nmeanx(:,es)*ones(length(es),1); temp(temp==0) = inf;
    nmeanx = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeanx;
    temp = nmeany(:,fs)*ones(length(fs),1); temp(temp==0) = inf;
    nmeany = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeany;
    temp = nmeanc(:,cs)*ones(length(cs),1); temp(temp==0) = inf;
    nmeanc = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeanc;
    temp = emean(:,ns)*ones(length(ns),1); temp(temp==0) = inf;
    emean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*emean;
    temp = fmean(:,ns)*ones(length(ns),1); temp(temp==0) = inf;
    fmean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*fmean;
    temp = cmean(:,ns)*ones(length(ns),1); temp(temp==0) = inf;
    cmean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*cmean;
    gg.nmeanc = nmeanc;
    
    % fields
    H = aa.H; b = aa.b; s = aa.s;
    C = aa.C; mu = pp.mu*ones(size(b));
    
    % viscosity at current u,v
    etabar = calculate_viscosity_local(u, v, pp, gg);
    
    % speed
    U = sqrt((nmeanx(:,es)*u(es)).^2 + (nmeany(:,fs)*v(fs)).^2);
    
    % sliding function
    tau_b_over_Ub = slide_fun_local(U(ns), N(ns), C(ns), mu(ns), pp);
    
    % stress BCs
    Txx = zeros(nIJ,1); Tyy = zeros(nIJ,1); Txy = zeros(cIJ,1);
    
    % assembly
    H_eta = sparse(nin,nin,H(nin).*etabar(nin),nIJ,nIJ);
    H_eta_c = sparse(cin,cin,cmean(cin,ns)*(H(ns).*etabar(ns)),cIJ,cIJ);
    
    % x-momentum
    Fx_res = -pp.c60*(emean(ein,ns)*H(ns)).*(eddx(ein,ns)*s(ns));
    if ~isempty(ebdy), Fx_res = Fx_res + pp.c62*(eddx(ein,nin)*(4*H_eta(nin,nin)*nddx(nin,ebdy)) + eddy(ein,cin)*(H_eta_c(cin,cin)*cddy(cin,ebdy)))*u(ebdy); end
    if ~isempty(fbdy), Fx_res = Fx_res + pp.c62*(eddx(ein,nin)*(2*H_eta(nin,nin)*nddy(nin,fbdy)) + eddy(ein,cin)*(H_eta_c(cin,cin)*cddx(cin,fbdy)))*v(fbdy); end
    if ~isempty(nbdyx), Fx_res = Fx_res + pp.c62*(eddx(ein,nbdyx)*(H(nbdyx).*Txx(nbdyx))); end
    if ~isempty(cbdy), Fx_res = Fx_res + pp.c62*(eddy(ein,cbdy)*((cmean(cbdy,ns)*H(ns)).*Txy(cbdy))); end
    
    Fx_ub = pp.c61*sparse(1:length(ein),1:length(ein),-emean(ein,ns)*tau_b_over_Ub,length(ein),length(ein)) ...
        + pp.c62*(eddx(ein,nin)*(4*H_eta(nin,nin)*nddx(nin,ein)) + eddy(ein,cin)*(H_eta_c(cin,cin)*cddy(cin,ein)));
    Fx_vb = pp.c62*(eddx(ein,nin)*(2*H_eta(nin,nin)*nddy(nin,fin)) + eddy(ein,cin)*(H_eta_c(cin,cin)*cddx(cin,fin)));
    
    % y-momentum
    Fy_res = -pp.c60*(fmean(fin,ns)*H(ns)).*(fddy(fin,ns)*s(ns));
    if ~isempty(ebdy), Fy_res = Fy_res + pp.c62*(fddy(fin,nin)*(2*H_eta(nin,nin)*nddx(nin,ebdy)) + fddx(fin,cin)*(H_eta_c(cin,cin)*cddy(cin,ebdy)))*u(ebdy); end
    if ~isempty(fbdy), Fy_res = Fy_res + pp.c62*(fddy(fin,nin)*(4*H_eta(nin,nin)*nddy(nin,fbdy)) + fddx(fin,cin)*(H_eta_c(cin,cin)*cddx(cin,fbdy)))*v(fbdy); end
    if ~isempty(nbdyy), Fy_res = Fy_res + pp.c62*(fddy(fin,nbdyy)*(H(nbdyy).*Tyy(nbdyy))); end
    if ~isempty(cbdy), Fy_res = Fy_res + pp.c62*(fddx(fin,cbdy)*((cmean(cbdy,ns)*H(ns)).*Txy(cbdy))); end
    
    Fy_ub = pp.c62*(fddy(fin,nin)*(2*H_eta(nin,nin)*nddx(nin,ein)) + fddx(fin,cin)*(H_eta_c(cin,cin)*cddy(cin,ein)));
    Fy_vb = pp.c61*sparse(1:length(fin),1:length(fin),-fmean(fin,ns)*tau_b_over_Ub,length(fin),length(fin)) ...
        + pp.c62*(fddy(fin,nin)*(4*H_eta(nin,nin)*nddy(nin,fin)) + fddx(fin,cin)*(H_eta_c(cin,cin)*cddx(cin,fin)));
    
    A = [Fx_ub Fx_vb; Fy_ub Fy_vb];
    rhs = -[Fx_res; Fy_res];

    % ===== Viscosity Newton correction (non-zero only when n_glen ~= 1) =====
    % When n > 1, eta = A^{-1/n} * I2^{(1/n-1)/2} depends on strain rates,
    % so d(eta)/d[u,v] ≠ 0.  The Picard Jacobian above freezes eta;
    % this block adds the missing  dF/d(eta) * d(eta)/d[u,v]  terms.
    if pp.n_glen ~= 1
        ne = length(ein); nf = length(fin);

        % --- Strain rates (same definitions as calculate_viscosity_local) ---
        epsxx   = nddx(:,es)*u(es);                                 % [nIJ x 1]
        epsyy   = nddy(:,fs)*v(fs);                                 % [nIJ x 1]
        epsxy_c = 0.5*(cddy(:,es)*u(es) + cddx(:,fs)*v(fs));       % corners [cIJ x 1]
        epsxy_n = nmeanc(:,cs)*epsxy_c(cs);                         % averaged to nodes [nIJ x 1]

        % --- Second invariant and derivative factor ---
        I2 = pp.eps_reg^2 + epsxx.^2 + epsyy.^2 + epsxx.*epsyy + epsxy_n.^2;
        p_exp = (1/pp.n_glen - 1) / 2;          % exponent factor
        beta  = p_exp * etabar ./ I2;            % d(eta)/d(I2) at each node [nIJ x 1]

        % --- d(eta)/du  [nIJ x ne] ---
        % d(I2)/du_k = (2*epsxx + epsyy)*nddx(:,k) + epsxy_n*[nmeanc*cddy](:,k)
        Dxy_u = nmeanc(:,cs) * cddy(cs,ein);    % maps u(ein) -> epsxy contrib at nodes
        deta_du = spdiags(beta, 0, nIJ, nIJ) * ...
            ( spdiags(2*epsxx + epsyy, 0, nIJ, nIJ) * nddx(:,ein) + ...
              spdiags(epsxy_n,         0, nIJ, nIJ) * Dxy_u );

        % --- d(eta)/dv  [nIJ x nf] ---
        Dxy_v = nmeanc(:,cs) * cddx(cs,fin);    % maps v(fin) -> epsxy contrib at nodes
        deta_dv = spdiags(beta, 0, nIJ, nIJ) * ...
            ( spdiags(epsxx + 2*epsyy, 0, nIJ, nIJ) * nddy(:,fin) + ...
              spdiags(epsxy_n,         0, nIJ, nIJ) * Dxy_v );

        % --- Stress coefficients at current u,v ---
        tau_xx   = 4*epsxx + 2*epsyy;                           % node: x-mom [nIJ x 1]
        tau_yy   = 2*epsxx + 4*epsyy;                           % node: y-mom [nIJ x 1]
        tau_xy_c = cddy(:,es)*u(es) + cddx(:,fs)*v(fs);        % corner: 2*epsxy [cIJ x 1]

        % --- d(H*eta at corners)/d[u,v] via cmean ---
        H_diag     = spdiags(H, 0, nIJ, nIJ);
        dHeta_c_du = cmean(cin,:) * H_diag * deta_du;           % [length(cin) x ne]
        dHeta_c_dv = cmean(cin,:) * H_diag * deta_dv;           % [length(cin) x nf]

        % --- Diagonal weight matrices ---
        Hn_tau_xx  = spdiags(H(nin) .* tau_xx(nin), 0, length(nin), length(nin));
        Hn_tau_yy  = spdiags(H(nin) .* tau_yy(nin), 0, length(nin), length(nin));
        tau_xy_diag = spdiags(tau_xy_c(cin),         0, length(cin), length(cin));

        % --- x-momentum viscosity Newton correction ---
        dVisc_Fx_du = pp.c62 * ( eddx(ein,nin) * Hn_tau_xx  * deta_du(nin,:) ...
                               + eddy(ein,cin) * tau_xy_diag * dHeta_c_du );
        dVisc_Fx_dv = pp.c62 * ( eddx(ein,nin) * Hn_tau_xx  * deta_dv(nin,:) ...
                               + eddy(ein,cin) * tau_xy_diag * dHeta_c_dv );

        % --- y-momentum viscosity Newton correction ---
        dVisc_Fy_du = pp.c62 * ( fddy(fin,nin) * Hn_tau_yy  * deta_du(nin,:) ...
                               + fddx(fin,cin) * tau_xy_diag * dHeta_c_du );
        dVisc_Fy_dv = pp.c62 * ( fddy(fin,nin) * Hn_tau_yy  * deta_dv(nin,:) ...
                               + fddx(fin,cin) * tau_xy_diag * dHeta_c_dv );

        A = A + [dVisc_Fx_du dVisc_Fx_dv; dVisc_Fy_du dVisc_Fy_dv];
    end

    % Export renormalized operators for use in adjoint gradient computation
    ops.emean  = emean;   % renormalized edge←node average [eIJ x nIJ]
    ops.fmean  = fmean;   % renormalized edge←node average [fIJ x nIJ]
    ops.nmeanx = nmeanx;  % renormalized node←x-edge average [nIJ x eIJ]
    ops.nmeany = nmeany;  % renormalized node←y-edge average [nIJ x fIJ]
    ops.nmeanc = nmeanc;  % renormalized node←corner average [nIJ x cIJ]
    ops.cmean  = cmean;   % renormalized corner←node average [cIJ x nIJ]
end

%% ============================================================
%  d(slide_fun)/dC  (analytical derivative of sliding law w.r.t. C)
%  ============================================================
function dslide_dC = compute_dslide_dC(U, N, C, mu, pp)
    % d(tau_b/Ub)/dC — max-based regularisation (matches nevis_velocity slide_fun)
    n = pp.n_slide;
    Ur = max(U, pp.Ub_reg);          % max-based (matches slide_fun_local)
    Nr = max(N, pp.N_slide_reg);     % max-based (matches slide_fun_local)
    
    bracket = mu.^(-n) .* Ur + C.^(-n) .* Nr.^n;
    
    dslide_dC = Nr .* Ur.^(1/n - 1) .* C.^(-n-1) .* Nr.^n .* bracket.^(-1/n - 1);
    % Note: C2 term has no C dependence, so no contribution
end

%% ============================================================
%  d(slide_fun)/dU  (derivative of sliding law w.r.t. speed Ub)
%  ============================================================
function dslide_dU = compute_dslide_dU(U, N, C, mu, pp)
    % d(tau_b/Ub)/dU — max-based regularisation (matches nevis_velocity slide_fun)
    % Ur = max(U, Ub_reg), Nr = max(N, N_slide_reg)
    % Chain rule: dslide/dU = dslide/dUr * dUr/dU, where dUr/dU = 1(U >= Ub_reg)
    %
    % Note: dslide_dU < 0 for n > 1 (shear thinning); = 0 for n = 1
    
    n = pp.n_slide;
    Ur = max(U, pp.Ub_reg);       % max-based (matches slide_fun_local)
    Nr = max(N, pp.N_slide_reg);  % max-based (matches slide_fun_local)
    
    bracket = mu.^(-n) .* Ur + C.^(-n) .* Nr.^n;
    
    % d(slide)/dUr:
    % Term 1: Nr * (1/n-1) * Ur^(1/n-2) * bracket^(-1/n)
    term1 = (1/n - 1) .* Nr .* Ur.^(1/n - 2) .* bracket.^(-1/n);
    
    % Term 2: Nr * Ur^(1/n-1) * (-1/n) * bracket^(-1/n-1) * mu^(-n)
    term2 = -(1/n) .* mu.^(-n) .* Nr .* Ur.^(1/n - 1) .* bracket.^(-1/n - 1);
    
    % Chain rule: dUr/dU = 1 if U >= Ub_reg, else 0 (subgradient of max)
    dslide_dU = (term1 + term2) .* double(U >= pp.Ub_reg);
end

%% ============================================================
%  Local viscosity calculation (same as in nevis_velocity)
%  ============================================================
function eta = calculate_viscosity_local(u, v, pp, gg)
    epsxx = gg.nddx(:,gg.es2)*u(gg.es2);
    epsyy = gg.nddy(:,gg.fs2)*v(gg.fs2);
    epsxy = 0.5*(gg.cddy(:,gg.es2)*u(gg.es2) + gg.cddx(:,gg.fs2)*v(gg.fs2));
    eta = pp.A_glen.^(-1/pp.n_glen) .* ...
        (pp.eps_reg^2 + epsxx.^2 + epsyy.^2 + epsxx.*epsyy + ...
        (gg.nmeanc(:,gg.cs2)*epsxy(gg.cs2)).^2).^((1/pp.n_glen-1)/2);
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
