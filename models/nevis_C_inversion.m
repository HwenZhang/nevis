%% NEVIS C-field inversion script
% Invert for the spatial basal friction coefficient C(x,y) using adjoint method
% Given: observed surface velocity (u_obs, v_obs), ice geometry (H, s, b),
%        effective pressure N, Coulomb coefficient mu
% Control variable: c = log(C) on nodes
%
% The sliding law is:
%   tau_b/Ub = N * Ub^(1/n-1) * (mu^(-n)*Ub + C^(-n)*N^n)^(-1/n) + C2*Ub^(1/n-1)
%
% Author: Hanwen Zhang
% Date: Feb 2026

clc; clear; close all; format compact;

%% ============================================================
%  1. Load spinup results and set up NEVIS framework
%  ============================================================
oo.root = './';
oo.code = './src';
oo.results = 'results';
oo.dataset = 'nevis_regional';
oo.casename = 'n2d_region_ice_inversion_test';
oo.fn = ['/',oo.casename];
oo.rn = [oo.root,oo.results,oo.fn];
oo.dn = [oo.root, 'data/', oo.dataset, '/'];
addpath(oo.code);

% load saved spinup state
load([oo.rn, oo.fn], 'pp','pd','ps','gg','aa','oo');

pd.u_b = 100/pd.ty;
pd.mu_s = 0.5;
pd.n_glen = 1;
pd.n_slide = 1;
eps = 0.05; % ratio of membrane stress to driving stress, used to set the ice softness to ensure 
% that membrane stresses are small in this test case, eps rises, A decreases, and viscosity increases
pd.A_glen = 1/2/((eps)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make membrane stress terms of dimensionless size eps in momentum equation
% Force-update ps.eta from the new pd.A_glen (loaded ps may have eta from a spinup with different eps).
% nevis_add_parameters_ice only sets ps.eta if it is absent, so we must clear it first.
if isfield(ps,'eta'), ps = rmfield(ps,'eta'); end
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo);
% NOTE: nevis_add_parameters_ice sets pp.A_glen = 1 (hardcoded), ignoring pd.A_glen.
% The eps-based pd.A_glen only affects ps.eta → pp.c62 (membrane stress coefficient).

% Picard iteration settings for velocity solver
oo.iter_max = 100;      % increase from default 10 — essential for convergence
oo.tol_vel = 1e-6;      % tighter convergence tolerance
oo.display_norms = 0;   % set to 1 for debugging
oo.verb = 0;            % suppress Picard iteration messages from nevis_velocity

fprintf('\n=== Key dimensionless parameters ===\n');
fprintf('pp.n_glen  = %d\n', pp.n_glen);
fprintf('pp.n_slide = %d\n', pp.n_slide);
fprintf('pp.A_glen  = %.4e (always 1 from nevis_add_parameters_ice)\n', pp.A_glen);
fprintf('pp.mu      = %.4e\n', pp.mu);
fprintf('pp.c60     = %.4e (driving stress)\n', pp.c60);
fprintf('pp.c61     = %.4e (friction)\n', pp.c61);
fprintf('pp.c62     = %.4e (membrane stress, controlled by eps=%.2f)\n', pp.c62, eps);
fprintf('pp.C (default) = %.4e\n', pp.C);
fprintf('pp.C2      = %.4e\n', pp.C2);
fprintf('================================\n\n');

% add regularisation parameters to pp if not already present (for backward compatibility)
if ~isfield(pp,'eps_reg'), pp.eps_reg = 1e-16; end % regularisation on strain rates
if ~isfield(pp,'Ub_reg'), pp.Ub_reg = 1e-16; end % regularisation on sliding speed (max-based, matches nevis_velocity)
if ~isfield(pp,'N_slide_reg'), pp.N_slide_reg = 1e-16; end % regularisation on effective pressure (max-based, matches nevis_velocity)
if ~isfield(pp,'taud_reg'), pp.taud_reg = 1e-16; end % regularisation on basal shear stress [ may not be needed ? ]
if ~isfield(pp,'C2'), pp.C2 = 0; end % added power-law coefficient in sliding law

%% ============================================================
%  2. Load or generate observed velocity
%  ============================================================
% Option A: load observed velocity from data file
load([oo.dn 'measures_for_nevis_140km.mat']);  % [m/yr]
dd = measures_for_nevis_140km;
u_obs = dd.u_obs_dim / (ps.u_b * pd.ty);  % non-dimensionalise
v_obs = dd.v_obs_dim / (ps.u_b * pd.ty);
u_obs = gg.emean2*u_obs(:);
v_obs = gg.fmean2*v_obs(:);
u_obs(gg.eout2) = NaN;
v_obs(gg.fout2) = NaN;

% observation mask: only fit where we have valid observations
% Exclude Dirichlet boundary edges — their velocity is prescribed, not controlled by C
obs_mask_e = zeros(gg.eIJ,1);  obs_mask_e(gg.ein2) = 1;
obs_mask_f = zeros(gg.fIJ,1);  obs_mask_f(gg.fin2) = 1;
obs_mask_e(gg.ebdy2) = 0;  % exclude Dirichlet x-boundary edges
obs_mask_f(gg.fbdy2) = 0;  % exclude Dirichlet y-boundary edges
fprintf('Excluded %d ebdy + %d fbdy Dirichlet edges from misfit\n', ...
    length(gg.ebdy2), length(gg.fbdy2));

% effective pressure for forward model (use spinup or assume N=1)
oo.initname = 'n2d_region_melt_meanperms1_Hreg1000_kappa1e_10_ks1e_03_mu1e1_spinup';
oo.initname = 'n2d_region_ice_inversion_test';
init_cond = load(['./results/' oo.initname '/' '0036.mat']); % load initial condition
vv = init_cond.vv;                                % load state variables from the initial 
% N_obs = 0.1*ones(gg.nIJ,1);
N_obs = (aa.phi_0 - vv.phi);

% no noise for real observations
u_obs_noisy = u_obs;
v_obs_noisy = v_obs;

% Option B: use forward model with a known C field as synthetic truth
% fprintf('Generating synthetic observations...\n');
% C_true = aa.C;  % true C field from spinup
% N_obs = ones(gg.nIJ,1);  % effective pressure (dimensionless)
% 
% u0 = zeros(gg.eIJ,1);
% v0 = zeros(gg.fIJ,1);
% u0(gg.eout2) = NaN;
% v0(gg.fout2) = NaN;
% [u_obs, v_obs] = nevis_velocity(aa.H, u0, v0, N_obs, aa, pp, gg, oo);
% 
% % add noise (optional)
% noise_level = 0.0;  % fraction of signal
% u_obs_noisy = u_obs;
% v_obs_noisy = v_obs;
% u_obs_noisy(gg.ein2) = u_obs(gg.ein2) .* (1 + noise_level*randn(length(gg.ein2),1));
% v_obs_noisy(gg.fin2) = v_obs(gg.fin2) .* (1 + noise_level*randn(length(gg.fin2),1));
% 
% % observation mask: only fit where we have valid observations
% obs_mask_e = zeros(gg.eIJ,1);  obs_mask_e(gg.ein2) = 1;
% obs_mask_f = zeros(gg.fIJ,1);  obs_mask_f(gg.fin2) = 1;

%% ============================================================
%  3. Inversion settings
%  ============================================================
opts_inv.u0_reg = 2e-1;      % velocity scale for relative misfit (dimensionless)

% Continuation schedule (coarse -> fine): start strongly regularized, then relax
opts_inv.alpha_schedule = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-7];
opts_inv.gamma_schedule = [1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-10];
opts_inv.alpha_schedule = [1e-3, 3e-4, 1e-4, 3e-5, 1e-6];
opts_inv.gamma_schedule = [1e-6, 3e-7, 1e-7, 3e-8, 1e-9];
opts_inv.alpha_schedule = [1e-3, 1e-4, 1e-5, 1e-6];
opts_inv.gamma_schedule = [1e-6, 1e-7, 1e-8, 1e-9];
opts_inv.max_iter_schedule = [30, 25, 5, 1];  % more iters when reg is strong
% Iteration controls (per stage)
opts_inv.max_iter_stage = 25;    % max iterations per continuation stage
opts_inv.max_iter_total = 100;    % safety cap across all stages

% Convergence thresholds (per stage)
opts_inv.J_tol = 1e-3;
opts_inv.dJ_tol = 1e-6;

% Outer loop: iterative C-N coupling
opts_inv.max_outer_iter = 5;      % max C-N coupling iterations
opts_inv.C_tol = 1e-3;             % relative change in C for outer convergence
opts_inv.N_tol = 1e-3;             % relative change in N for outer convergence

% Initialize alpha/gamma with first schedule values (needed for initial test call)
opts_inv.alpha = opts_inv.alpha_schedule(1);
opts_inv.gamma = opts_inv.gamma_schedule(1);
opts_inv.verbose = false;
opts_inv.check_grad = false;  % set true to run Taylor test & per-component FD check

% History
history.stage = [];
history.iter  = [];
history.J     = [];
history.normg = [];

% Outer loop history
outer_history.outer_iter = [];
outer_history.J_final = [];
outer_history.C_change = [];
outer_history.N_change = [];

% initial guess for C (perturbed from truth for synthetic test)
% load(['./data/C_inversion_results.mat'], 'C_hat_dim');
% C_init = 0.01*C_hat_dim * (ps.u_b^(1/pp.n_slide) / ps.tau);

C_init = ones(gg.nIJ,1);  % uniform prior
c_prior = log(C_init);

% U_obs_nodes = sqrt((gg.nmeanx2(:,gg.es2)*u_obs(gg.es2)).^2 + ...
%     (gg.nmeany2(:,gg.fs2)*v_obs(gg.fs2)).^2);
% % Surface gradient: nodes -> edges (eddx, fddy), then edges -> nodes (nmeanx2, nmeany2)
% dsdx_e = gg.eddx * aa.s;   % ds/dx at x-edges [eIJ x 1]
% dsdy_f = gg.fddy * aa.s;   % ds/dy at y-edges [fIJ x 1]
% dsdx_n = gg.nmeanx2(:,gg.es2) * dsdx_e(gg.es2);  % average to nodes
% dsdy_n = gg.nmeany2(:,gg.fs2) * dsdy_f(gg.fs2);  % average to nodes
% tau_d = pp.c60 * aa.H .* sqrt(dsdx_n.^2 + dsdy_n.^2);
% C_init = tau_d ./ (max(U_obs_nodes, 0.01).^(1/pp.n_slide));
% C_init = max(C_init, 0.01); 
% C_init = min(C_init, 100); 
% c_prior = log(C_init);

% plot the prior C field
% C_plot = reshape(C_init, gg.nI, gg.nJ);
% xx = reshape(gg.nx, gg.nI, gg.nJ);
% yy = reshape(gg.ny, gg.nI, gg.nJ);
% figure; pcolor(xx, yy, reshape(C_init, gg.nI, gg.nJ)); shading flat; colorbar;
% return;

% build gradient operator L for regularisation (using grid finite differences)
% L * c gives spatial gradient of c (defined on nodes)
% eddx: [eIJ x nIJ] - x derivative at x-edges
% fddy: [fIJ x nIJ] - y derivative at y-edges
Lx = gg.eddx(gg.ein2, :);  % x-gradient at interior x-edges [length(ein2) x nIJ]
Ly = gg.fddy(gg.fin2, :);  % y-gradient at interior y-edges [length(fin2) x nIJ]
L = [Lx; Ly];

%% ============================================================
%  4. Run iterative C-N coupled inversion
%  ============================================================
fprintf('\n========================================\n');
fprintf('Starting iterative C-N coupled inversion\n');
fprintf('Max outer iterations: %d\n', opts_inv.max_outer_iter);
fprintf('========================================\n\n');

c0 = c_prior;
C_prev = exp(c0);
N_current = N_obs;  % Initial N from input
vv_hydro = [];       % hydrology state (empty = fresh start on first iteration)

for outer_iter = 1:opts_inv.max_outer_iter
    fprintf('\n╔════════════════════════════════════════╗\n');
    fprintf('║  OUTER ITERATION %d/%d                  ║\n', outer_iter, opts_inv.max_outer_iter);
    fprintf('╚════════════════════════════════════════╝\n\n');

    % DEBUG: test objective at initial point (only first outer iteration)
    if outer_iter == 1
        fprintf('Testing objective at initial point...\n');
        fprintf('NaN in u_obs(ein2): %d / %d\n', sum(isnan(u_obs_noisy(gg.ein2))), length(gg.ein2));
        fprintf('NaN in v_obs(fin2): %d / %d\n', sum(isnan(v_obs_noisy(gg.fin2))), length(gg.fin2));
        fprintf('NaN in N_current: %d / %d\n', sum(isnan(N_current)), length(N_current));
        fprintf('NaN in c0: %d / %d\n', sum(isnan(c0)), length(c0));

        % Mask out NaN observations
        obs_mask_e(isnan(u_obs_noisy)) = 0;
        obs_mask_f(isnan(v_obs_noisy)) = 0;
        % Replace NaN in obs with 0 (masked anyway)
        u_obs_noisy(isnan(u_obs_noisy)) = 0;
        v_obs_noisy(isnan(v_obs_noisy)) = 0;
    end

    % Replace NaN in N_current (softplus Reg_Ni handles N<0 gracefully in sliding law)
    N_current(isnan(N_current)) = 0;

    if outer_iter == 1
        fprintf('  After masking: obs_mask_e sum = %d, obs_mask_f sum = %d\n', ...
            sum(obs_mask_e), sum(obs_mask_f));
    end

    % --- Apply adaptive regularization AFTER NaN cleaning ---
    % Compute velocity gradient magnitude at edges using cleaned observations
    dudx_e = gg.eddx * (gg.nmeanx2(:,gg.es2) * u_obs_noisy(gg.es2));  % du/dx at x-edges
    dvdy_f = gg.fddy * (gg.nmeany2(:,gg.fs2) * v_obs_noisy(gg.fs2));  % dv/dy at y-edges
    grad_mag_e = abs(dudx_e(gg.ein2));  % velocity gradient at x-edges
    grad_mag_f = abs(dvdy_f(gg.fin2));  % velocity gradient at y-edges

    % Adaptive weight: w = 1 / (1 + beta * |grad(u)|/mean(|grad(u)|))
    % High gradient -> small w -> less regularization
    beta_adaptive = 0.0;  % tuning parameter (increase to reduce reg more in shear zones)
    w_e = 1 ./ (1 + beta_adaptive * grad_mag_e / mean(grad_mag_e));
    w_f = 1 ./ (1 + beta_adaptive * grad_mag_f / mean(grad_mag_f));

    % Apply weights to regularization operator
    Lx_base = gg.eddx(gg.ein2, :);  % unweighted x-gradient operator
    Ly_base = gg.fddy(gg.fin2, :);  % unweighted y-gradient operator
    Lx = spdiags(w_e, 0, length(gg.ein2), length(gg.ein2)) * Lx_base;
    Ly = spdiags(w_f, 0, length(gg.fin2), length(gg.fin2)) * Ly_base;
    L = [Lx; Ly];  % weighted regularization operator

    fprintf('Adaptive regularization: weight range [%.3f, %.3f]\n', min([w_e; w_f]), max([w_e; w_f]));

    % Recreate objective with cleaned data and current N
    obj_fun = @(c) objective_and_grad(c, u_obs_noisy, v_obs_noisy, ...
        obs_mask_e, obs_mask_f, N_current, aa, pp, ps, pd, gg, oo, L, c_prior, opts_inv);

    if outer_iter == 1
        [J0, g0] = obj_fun(c0);
        fprintf('  J(c0) = %.6e, norm(g) = %.6e\n', J0, norm(g0));
        if isnan(J0) || isinf(J0)
            error('Objective still undefined. Check forward model.');
        end

      if opts_inv.check_grad
        % === Gradient verification (finite-difference Taylor test) ===
        fprintf('\n--- Gradient verification (Taylor test) ---\n');
        % Pick a random perturbation direction
        rng(42);
        dc = randn(size(c0));
        dc = dc / norm(dc);  % unit direction
        
        % Directional derivative from adjoint
        dJdc_adj = g0' * dc;
        
        % Taylor test: J(c0 + h*dc) = J(c0) + h*dJ/dc*dc + O(h^2)
        % First-order:  |J(c0+h*dc) - J0|           should be O(h)
        % With gradient: |J(c0+h*dc) - J0 - h*g'*dc| should be O(h^2)
        h_vals = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
        fprintf('  %12s  %14s  %14s  %10s  %10s\n', ...
            'h', '|J(c+h*dc)-J0|', '|...-h*g''dc|', 'ratio_1st', 'ratio_2nd');
        err1_prev = NaN; err2_prev = NaN;
        r2_vals = [];
        for ih = 1:length(h_vals)
            h = h_vals(ih);
            J_pert = obj_fun(c0 + h * dc);
            err1 = abs(J_pert - J0);
            err2 = abs(J_pert - J0 - h * dJdc_adj);
            if ih > 1
                r1 = err1_prev / err1;
                r2 = err2_prev / err2;
                r2_vals(end+1) = r2; %#ok<AGROW>
            else
                r1 = NaN; r2 = NaN;
            end
            fprintf('  %12.1e  %14.6e  %14.6e  %10.2f  %10.2f\n', ...
                h, err1, err2, r1, r2);
            err1_prev = err1; err2_prev = err2;
        end
        fprintf('  Expected ratios: ~10 (1st order), ~100 (2nd order)\n');
        fprintf('  If ratio_2nd ≈ 100, adjoint gradient is correct.\n');
        fprintf('  If ratio_2nd ≈ 10 or erratic, gradient has a bug.\n');
        
        % Estimate gradient relative error from Taylor test
        % At smallest h: err2 ≈ h*|ε|, err1 ≈ h*|g·dc|
        % So relative error ≈ err2/err1
        grad_rel_err = err2_prev / max(err1_prev, 1e-30);
        fprintf('\n  Estimated gradient relative error: %.2f%%\n', 100*grad_rel_err);
        
        if grad_rel_err < 0.05
            fprintf('  ✅ Gradient accurate to ~%.1f%% (acceptable for L-BFGS)\n', 100*grad_rel_err);
            if grad_rel_err > 0.005
                fprintf('     Note: small residual error likely from Picard solver tolerance.\n');
                fprintf('     To reduce: tighten oo.tol_vel (currently %.1e).\n\n', oo.tol_vel);
            else
                fprintf('\n');
            end
        elseif grad_rel_err < 0.20
            fprintf('  ⚠️  Gradient has ~%.0f%% error. Usable but suboptimal.\n', 100*grad_rel_err);
            fprintf('     Consider tightening oo.tol_vel or checking adjoint.\n\n');
        else
            fprintf('  ❌ GRADIENT ERROR TOO LARGE: %.0f%%\n', 100*grad_rel_err);
            error('Gradient relative error > 20%%. Fix the adjoint or tighten Picard tolerance before proceeding.');
        end
        
        % Also check a few individual components via central differences
        fprintf('--- Per-component gradient check (central FD, 20 random nodes) ---\n');
        n_check = 20;
        idx_check = randperm(length(c0), n_check);
        h_fd = 1e-5;
        fprintf('  %6s  %14s  %14s  %14s  %10s\n', ...
            'node', 'g_adjoint', 'g_FD_central', 'abs_diff', 'rel_diff');
        g_match = true;
        for ic = 1:n_check
            idx = idx_check(ic);
            c_plus = c0; c_plus(idx) = c_plus(idx) + h_fd;
            c_minus = c0; c_minus(idx) = c_minus(idx) - h_fd;
            J_plus = obj_fun(c_plus);
            J_minus = obj_fun(c_minus);
            g_fd = (J_plus - J_minus) / (2 * h_fd);
            g_adj = g0(idx);
            abs_diff = abs(g_adj - g_fd);
            rel_diff = abs_diff / max(abs(g_adj), abs(g_fd) + 1e-30);
            fprintf('  %6d  %14.6e  %14.6e  %14.6e  %10.4f\n', ...
                idx, g_adj, g_fd, abs_diff, rel_diff);
            if rel_diff > 0.1 && abs_diff > 1e-8
                g_match = false;
            end
        end
        if g_match
            fprintf('  ✅ Gradient matches finite differences (all rel_diff < 10%%)\n\n');
        else
            fprintf('  ❌ GRADIENT MISMATCH DETECTED\n');
            error('Gradient verification failed (per-component FD check): adjoint gradient does not match central finite differences (rel_diff > 10%%). Fix the adjoint before proceeding.');
        end
      end % check_grad
    end

    fprintf('Starting C-field inversion (continuation) with current N...\n');
    
    c_stage = c0;
    iter_total = 0;
    J_hat = NaN; exitflag = NaN; output = struct();

    for k = 1:min(numel(opts_inv.alpha_schedule), numel(opts_inv.gamma_schedule))
        opts_inv.alpha = opts_inv.alpha_schedule(k);
        opts_inv.gamma = opts_inv.gamma_schedule(k);
        max_iter_k = opts_inv.max_iter_schedule(min(k, end));  % per-stage limit

        fprintf('\n--- Stage %d/%d: alpha=%.3e, gamma=%.3e, max_iter=%d ---\n', k, ...
            min(numel(opts_inv.alpha_schedule), numel(opts_inv.gamma_schedule)), opts_inv.alpha, opts_inv.gamma, max_iter_k);

        % Recreate objective with current regularization weights
        obj_fun = @(c) objective_and_grad(c, u_obs_noisy, v_obs_noisy, ...
            obs_mask_e, obs_mask_f, N_current, aa, pp, ps, pd, gg, oo, L, c_prior, opts_inv);
        opts = optimoptions('fminunc', ...
            'Algorithm',                'quasi-newton', ...
            'HessianApproximation',     {'lbfgs', 20}, ...
            'SpecifyObjectiveGradient', true, ...
            'Display',                  'iter', ...
            'MaxIterations',            max_iter_k, ...
            'MaxFunctionEvaluations',   max_iter_k * 20, ...
            'OptimalityTolerance',      1e-6, ...
            'StepTolerance',            1e-12, ...
            'ObjectiveLimit',           1e-4);

        [c_stage, J_stage, exitflag_stage, output_stage] = fminunc(obj_fun, c_stage, opts);

        % Bookkeeping
        iter_total = iter_total + output_stage.iterations;
        J_hat = J_stage;
        exitflag = exitflag_stage;
        output = output_stage;

        % Record history (lightweight: no extra forward solve)
        history.stage(end+1,1) = k;
        history.iter(end+1,1)  = iter_total;
        history.J(end+1,1)     = J_stage;
        history.normg(end+1,1) = output_stage.firstorderopt;

        fprintf('Stage %d done: J=%.6e, firstorderopt=%.3e, iters=%d (total=%d)\n', ...
            k, J_stage, output_stage.firstorderopt, output_stage.iterations, iter_total);

        % Early stop conditions
        if J_stage < opts_inv.J_tol
            fprintf('*** Early stop: J < J_tol (%.2e) ***\n', opts_inv.J_tol);
            break;
        end
        if output_stage.firstorderopt < opts_inv.dJ_tol
            fprintf('*** Early stop: firstorderopt < dJ_tol (%.2e) ***\n', opts_inv.dJ_tol);
            break;
        end
        if iter_total >= opts_inv.max_iter_total
            fprintf('*** Early stop: reached max_iter_total=%d ***\n', opts_inv.max_iter_total);
            break;
        end
    end

    c_hat = c_stage;
    C_hat = exp(c_hat);

    % --- Compute change in C from previous iteration ---
    C_change = norm(C_hat - C_prev) / norm(C_prev);
    fprintf('\nC relative change: %.4e\n', C_change);

    % Record outer iteration history
    outer_history.outer_iter(end+1,1) = outer_iter;
    outer_history.J_final(end+1,1) = J_hat;
    outer_history.C_change(end+1,1) = C_change;
    
    C_hat_dim = C_hat * (ps.tau / ps.u_b^(1/pp.n_slide));
    save(sprintf('./data/C_inversion_outer%d.mat', outer_iter), ...
        'C_hat', 'C_hat_dim', 'N_current', 'outer_history', 'J_hat');
    fprintf('  Saved intermediate results to C_inversion_outer%d.mat\n', outer_iter);
    
    % --- Update N using forward model with inverted C ---
    fprintf('Running forward model to update N...\n');
    N_old = N_current;  % save old N for convergence check
    
    aa.C = C_hat;
    [u_inv, v_inv] = nevis_velocity(aa.H, u_obs_noisy, v_obs_noisy, N_current, aa, pp, gg, oo);
    C_dim = C_hat * (ps.tau / ps.u_b^(1/pp.n_slide));
    vv.hydro.u = u_inv;
    vv.hydro.v = v_inv;
    [N_new, vv_hydro] = nevis_run_fwd_hydrology(C_dim, C_hat, vv_hydro, pd.mu_s);
    % N_new = N_old; 

    N_new(isnan(N_new)) = 0;
    % Compute change in N
    N_change = norm(N_new - N_old) / norm(N_old);
    fprintf('  N relative change: %.4e\n', N_change);
    outer_history.N_change(end+1,1) = N_change;
    
    % --- Check outer loop convergence ---
    fprintf('\n--- Outer iteration %d summary ---\n', outer_iter);
    fprintf('  J = %.6e\n', J_hat);
    fprintf('  C change = %.4e (tol = %.4e)\n', C_change, opts_inv.C_tol);
    fprintf('  N change = %.4e (tol = %.4e)\n', N_change, opts_inv.N_tol);
    
    if C_change < opts_inv.C_tol && N_change < opts_inv.N_tol
        fprintf('\n*** Outer loop converged! ***\n');
        break;
    end
    
    if outer_iter >= opts_inv.max_outer_iter
        fprintf('\n*** Reached max outer iterations ***\n');
        break;
    end
    
    % --- Prepare for next outer iteration ---
    C_prev = C_hat;
    N_current = N_new;
    c0 = c_hat;  % Warm start from previous solution
    
    fprintf('\nPreparing for next outer iteration...\n');
end

fprintf('\n========================================\n');
fprintf('Iterative C-N inversion complete!\n');
fprintf('Total outer iterations: %d\n', outer_iter);
fprintf('========================================\n\n');

%% Compute final C and print
fprintf('Inversion complete. Final J = %.6e, exit flag = %d\n', J_hat, exitflag);

% Dimensionalise C for interpretation, save results
C_hat_dim = C_hat * (ps.tau / ps.u_b^(1/pp.n_slide));
save('./data/C_inversion_results.mat', 'C_hat_dim', 'c_hat', 'history', 'opts_inv', 'J_hat', 'exitflag');

figure('Name','Convergence','Position',[100 100 1200 400]);
subplot(1,3,1);
semilogy(history.iter, history.J, 'o-','LineWidth',1.5);
xlabel('Total iterations'); ylabel('J'); grid on;
title('Cost by continuation stage');

subplot(1,3,2);
semilogy(history.iter, history.normg, 's-','LineWidth',1.5);
xlabel('Total iterations'); ylabel('firstorderopt'); grid on;
title('Optimality measure');

subplot(1,3,3);
yyaxis left;
semilogy(outer_history.outer_iter, outer_history.C_change, 'o-', 'LineWidth', 1.5);
ylabel('Relative C change'); grid on;
yyaxis right;
semilogy(outer_history.outer_iter, outer_history.N_change, 's-', 'LineWidth', 1.5);
ylabel('Relative N change');
xlabel('Outer iteration');
title('C-N coupling convergence');

%% ============================================================
%  5. Plot results
%  ============================================================
C_prior = C_init;
% C_prior = exp(C_prior);
% compute final velocity (warm start from obs)
aa_inv = aa;
aa_inv.C = C_hat;
u_inv = u_obs_noisy; v_inv = v_obs_noisy;  % warm start — use NaN-cleaned obs
% Ensure interior edges that are still NaN get a finite initial guess
u_inv(isnan(u_inv) & ~ismember((1:gg.eIJ)', gg.eout2)) = 0;
v_inv(isnan(v_inv) & ~ismember((1:gg.fIJ)', gg.fout2)) = 0;
[u_inv, v_inv] = nevis_velocity(aa.H, u_inv, v_inv, N_current, aa_inv, pp, gg, oo);

figure('Position',[100 100 1200 500]);
C_scale = ps.tau / ps.u^(1/pp.n_slide);  % [Pa s^(1/n) m^(-1/n)]

% (a) Prior C
subplot(1,3,1);
zC_prior = reshape(log10(C_prior * (ps.tau / ps.u_b^(1/pp.n_slide))), gg.nI, gg.nJ);
zC_prior(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, zC_prior); shading flat; colorbar;
title('(a) Prior log_{10}(C)'); xlabel('x'); ylabel('y'); axis equal tight;

% (b) Inverted C
subplot(1,3,2);
zC_hat = reshape(log10(C_hat_dim), gg.nI, gg.nJ);
zC_hat(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, zC_hat); shading flat; colorbar;
title('(b) Inverted log_{10}(C)'); xlabel('x'); ylabel('y'); axis equal tight;

% (c) Speed comparison (obs vs inverted)
subplot(1,3,3);
Un_obs = sqrt((gg.nmeanx2(:,gg.es2)*u_obs(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_obs(gg.fs2)).^2);
Un_inv = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);
Un_misfit = Un_inv - Un_obs;
Un_misfit(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, ps.u_b*pd.ty*reshape(Un_misfit,gg.nI,gg.nJ)); shading flat; colorbar;
title('(c) Speed misfit (m/yr)'); xlabel('x'); ylabel('y'); axis equal tight;
caxis([-50 50]);  % symmetric colorbar
drawnow;

%% ============================================================
%  5b. Spatial error diagnostics
%  ============================================================
% Node-based speed
Un_obs = sqrt((gg.nmeanx2(:,gg.es2)*u_obs(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_obs(gg.fs2)).^2);
Un_inv = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);

% Dimensional speeds [m/yr]
Un_obs_dim = Un_obs * ps.u_b * pd.ty;
Un_inv_dim = Un_inv * ps.u_b * pd.ty;
Un_misfit_dim = Un_inv_dim - Un_obs_dim;

% Relative error per node
rel_err = abs(Un_misfit_dim) ./ max(Un_obs_dim, 1);  % floor at 1 m/yr

% Find worst nodes (speed too low)
too_low = find(Un_misfit_dim < -20 & Un_obs_dim > 50);  % model > 20 m/yr too slow, obs > 50 m/yr
fprintf('\n=== Spatial Error Diagnostics ===\n');
fprintf('Nodes where model is >20 m/yr too slow (obs > 50 m/yr): %d\n', length(too_low));

if ~isempty(too_low)
    fprintf('\n  Top 20 worst under-predictions:\n');
    [~, idx_sort] = sort(Un_misfit_dim(too_low), 'ascend');
    n_show = min(20, length(idx_sort));
    for ii = 1:n_show
        jj = too_low(idx_sort(ii));
        fprintf('    node %5d: obs=%.1f, mod=%.1f, misfit=%.1f m/yr, C=%.3e, N=%.3e\n', ...
            jj, Un_obs_dim(jj), Un_inv_dim(jj), Un_misfit_dim(jj), ...
            C_hat(jj), N_current(jj));
    end
    
    % Statistics in the slow-bias region
    fprintf('\n  Slow-bias region statistics:\n');
    fprintf('    C range: [%.3e, %.3e]\n', min(C_hat(too_low)), max(C_hat(too_low)));
    fprintf('    N range: [%.3e, %.3e]\n', min(N_current(too_low)), max(N_current(too_low)));
    fprintf('    C_prior range: [%.3e, %.3e]\n', min(C_init(too_low)), max(C_init(too_low)));
end

% Find nodes where speed is too high
too_high = find(Un_misfit_dim > 20 & Un_obs_dim > 50);
fprintf('Nodes where model is >20 m/yr too fast (obs > 50 m/yr): %d\n', length(too_high));

figure('Name', 'Error Diagnostics', 'Position', [100 600 1600 500]);

% (a) Observed speed
subplot(2,4,1);
z = reshape(Un_obs_dim, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('(a) Observed speed [m/yr]'); axis equal tight;
clim([0 1000]);

% (b) Modeled speed
subplot(2,4,2);
z = reshape(Un_inv_dim, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('(b) Modeled speed [m/yr]'); axis equal tight;
clim([0 1000]);

% (c) Misfit (mod - obs)
subplot(2,4,3);
z = reshape(Un_misfit_dim, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('(c) Misfit [m/yr]'); axis equal tight;
caxis([-100 100]); colormap(gca, interp1([0 0.5 1],[0 0 1; 1 1 1; 1 0 0],linspace(0,1,256)));

% (d) Relative error
subplot(2,4,4);
z = reshape(100*rel_err, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('(d) Relative error [%]'); axis equal tight;
caxis([0 50]);

% (e) Inverted C (log10)
subplot(2,4,5);
z = reshape(log10(C_hat_dim), gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('(e) Inverted log_{10}(C_{dim})'); axis equal tight;

% (f) N field
subplot(2,4,6);
z = reshape(N_current * ps.phi, gg.nI, gg.nJ) / 1e6; z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('(f) N [MPa]'); axis equal tight;

% (g) Scatter: obs speed vs C
subplot(2,4,7);
nin = gg.nin2;
scatter(Un_obs_dim(nin), log10(C_hat_dim(nin)), 5, Un_misfit_dim(nin), 'filled');
xlabel('Observed speed [m/yr]'); ylabel('log_{10}(C_{dim})');
colorbar; caxis([-50 50]); colormap(gca, interp1([0 0.5 1],[0 0 1; 1 1 1; 1 0 0],linspace(0,1,256)));
title('(g) C vs speed (color=misfit)');

% (h) Scatter: N vs misfit
subplot(2,4,8);
scatter(N_current(nin)*ps.phi/1e6, Un_misfit_dim(nin), 5, log10(C_hat_dim(nin)), 'filled');
xlabel('N [MPa]'); ylabel('Speed misfit [m/yr]');
colorbar; title('(h) N vs misfit (color=log_{10}C)');

drawnow;

% --- Print summary of possible causes ---
fprintf('\n=== Possible Causes of Speed Under-prediction ===\n');

% Check 1: C hitting bounds
c_hat_range = [min(c_hat(nin)), max(c_hat(nin))];
fprintf('1. c = log(C) range: [%.2f, %.2f]\n', c_hat_range(1), c_hat_range(2));
if c_hat_range(2) > 10
    fprintf('   ⚠️  C is very large (>%.0e) — may need upper bound or stronger regularization\n', exp(10));
end
if c_hat_range(1) < -10
    fprintf('   ⚠️  C is very small (<%.0e) — may need lower bound\n', exp(-10));
end

% Check 2: N too high in slow regions
if ~isempty(too_low)
    N_too_low = N_current(too_low) * ps.phi;  % dimensional [Pa]
    fprintf('2. N in slow-bias region: mean=%.2f MPa, max=%.2f MPa\n', ...
        mean(N_too_low)/1e6, max(N_too_low)/1e6);
    if mean(N_too_low) > 1e6
        fprintf('   ⚠️  N > 1 MPa in under-predicted region → sliding law suppresses speed\n');
        fprintf('   Fix: check if N_obs is realistic, or reduce N\n');
    end
end

% Check 3: Regularization too strong
% Compute J components from final solution (use cleaned obs to avoid NaN)
Lc = L * c_hat;
Jmis_final = 0.5 * (sum((obs_mask_e(gg.ein2) .* (u_inv(gg.ein2) - u_obs_noisy(gg.ein2)) ./ (abs(u_obs_noisy(gg.ein2)) + opts_inv.u0_reg)).^2) + ...
                     sum((obs_mask_f(gg.fin2) .* (v_inv(gg.fin2) - v_obs_noisy(gg.fin2)) ./ (abs(v_obs_noisy(gg.fin2)) + opts_inv.u0_reg)).^2));
Jreg_final = 0.5 * opts_inv.alpha * (Lc' * Lc);
Jdamp_final = 0.5 * opts_inv.gamma * sum((c_hat - c_prior).^2);
fprintf('3. Final alpha=%.2e, gamma=%.2e\n', opts_inv.alpha, opts_inv.gamma);
fprintf('   Jmis=%.2e, Jreg=%.2e, Jdamp=%.2e\n', ...
    Jmis_final, Jreg_final, Jdamp_final);
if Jreg_final > 0.5 * Jmis_final
    fprintf('   ⚠️  Jreg > 0.5*Jmis — regularization may be too strong\n');
end

% Check 5: Sliding regime diagnostic
% Coulomb-Weertman law:
%   C large → Coulomb limit:  tau_b ≈ mu*N (MAXIMUM friction for given N)
%   C small → Weertman limit: tau_b ≈ C*Ub^(1/n) (friction decreases with C)
%   d(tau_b/Ub)/dC > 0 ALWAYS: increasing C always increases friction
%   To reduce friction → decrease C
fprintf('\n5. Sliding regime diagnostic (mu=%.2f, n_slide=%d):\n', pp.mu, pp.n_slide);
n_sl = pp.n_slide;
U_nodes = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);
Ub_n = max(U_nodes, pp.Ub_reg);      % max-based (matches nevis_velocity)
Nr_n = max(N_current, pp.N_slide_reg); % max-based (matches nevis_velocity)

% Regime ratio R = C^{-n}*N^n / (mu^{-n}*Ub)
%   R >> 1 → Weertman limit (C small, tau_b ≈ C*Ub^{1/n}, friction controllable)
%   R << 1 → Coulomb limit  (C large, tau_b ≈ mu*N, friction at maximum)
regime_ratio = (C_hat.^(-n_sl) .* Nr_n.^n_sl) ./ (pp.mu.^(-n_sl) .* Ub_n);
fprintf('   Regime ratio R = C^{-n}N^n / (mu^{-n}Ub):  R>>1=Weertman, R<<1=Coulomb\n');
fprintf('     All interior: median=%.2e, mean=%.2e\n', ...
    median(regime_ratio(gg.nin2)), mean(regime_ratio(gg.nin2)));

if ~isempty(too_low)
    fprintf('     Slow-bias region: median=%.2e, mean=%.2e\n', ...
        median(regime_ratio(too_low)), mean(regime_ratio(too_low)));
    
    n_coulomb = sum(regime_ratio(too_low) < 0.1);  % C large → max friction
    n_weertman = sum(regime_ratio(too_low) > 10);   % C small → friction = C*Ub^{1/n}
    n_transition = length(too_low) - n_coulomb - n_weertman;
    
    fprintf('     Coulomb limit  (R<0.1, C large): %d / %d (%.0f%%)\n', ...
        n_coulomb, length(too_low), 100*n_coulomb/length(too_low));
    fprintf('     Weertman limit (R>10,  C small): %d / %d (%.0f%%)\n', ...
        n_weertman, length(too_low), 100*n_weertman/length(too_low));
    fprintf('     Transition     (0.1<R<10):       %d / %d (%.0f%%)\n', ...
        n_transition, length(too_low), 100*n_transition/length(too_low));
    
    fprintf('\n     C at slow-bias nodes: median=%.3e, range=[%.3e, %.3e]\n', ...
        median(C_hat(too_low)), min(C_hat(too_low)), max(C_hat(too_low)));
    
    if n_coulomb > 0.3 * length(too_low)
        fprintf('   ⚠️  Many slow nodes have LARGE C → Coulomb limit (tau_b ≈ mu*N)\n');
        fprintf('      Large C = MAXIMUM friction. To speed up, C should DECREASE.\n');
        fprintf('      The optimizer drove C the WRONG direction!\n');
        fprintf('      → Likely cause: adjoint gradient has wrong sign. Run Taylor test.\n');
    end
    
    if n_weertman > 0.3 * length(too_low)
        fprintf('   ℹ️  Many slow nodes have SMALL C → Weertman limit (tau_b ≈ C*Ub^{1/n})\n');
        fprintf('      Friction is already low (controlled by C), but speed still too slow.\n');
        % Compute actual friction at Weertman-limit nodes
        weertman_idx = too_low(regime_ratio(too_low) > 10);
        tau_b_weertman = C_hat(weertman_idx) .* Ub_n(weertman_idx).^(1/n_sl) * ps.phi;
        fprintf('      tau_b in Weertman region: mean=%.0f kPa, max=%.0f kPa\n', ...
            mean(tau_b_weertman)/1e3, max(tau_b_weertman)/1e3);
        fprintf('      → Speed limited by something other than friction:\n');
        fprintf('        (a) Membrane stress too weak (eps=%.4f too small)\n', eps);
        fprintf('        (b) Driving stress locally low\n');
        fprintf('        (c) Boundary conditions constraining flow\n');
    end
    
    % Driving stress for context
    dsdx_e = gg.eddx * aa.s;
    dsdy_f = gg.fddy * aa.s;
    dsdx_n = gg.nmeanx2(:,gg.es2) * dsdx_e(gg.es2);
    dsdy_n = gg.nmeany2(:,gg.fs2) * dsdy_f(gg.fs2);
    tau_d = pp.c60 * aa.H .* sqrt(dsdx_n.^2 + dsdy_n.^2) * ps.phi;
    tau_coulomb_cap = pp.mu * Nr_n * ps.phi;
    
    fprintf('\n     Force comparison at slow-bias nodes:\n');
    fprintf('       Driving stress:  mean=%.0f kPa, max=%.0f kPa\n', ...
        mean(tau_d(too_low))/1e3, max(tau_d(too_low))/1e3);
    fprintf('       Coulomb cap mu*N: mean=%.0f kPa, max=%.0f kPa\n', ...
        mean(tau_coulomb_cap(too_low))/1e3, max(tau_coulomb_cap(too_low))/1e3);
    if mean(tau_coulomb_cap(too_low)) < mean(tau_d(too_low))
        fprintf('       mu*N < tau_d → Coulomb cap below driving stress → fast flow IS possible\n');
        fprintf('       If speed still too low with large C, gradient sign is likely wrong.\n');
    else
        fprintf('       mu*N > tau_d → even Coulomb cap exceeds driving stress\n');
        fprintf('       Reduce N or mu to allow fast flow.\n');
    end
end

% Check 6: Membrane stress diagnostic
fprintf('\n6. Membrane stress check:\n');
fprintf('   eps (membrane/driving ratio) = %.4f\n', eps);
fprintf('   c60 (driving) = %.4e\n', pp.c60);
fprintf('   c61 (friction) = %.4e\n', pp.c61);
fprintf('   c62 (membrane) = %.4e\n', pp.c62);
fprintf('   c62/c60 = %.4e (effective membrane importance)\n', pp.c62/pp.c60);
fprintf('   c62/c61 = %.4e (membrane vs friction ratio)\n', pp.c62/pp.c61);

% Estimate force balance at worst nodes
if ~isempty(too_low)
    % Driving stress
    dsdx_e_diag = gg.eddx * aa.s;
    dsdy_f_diag = gg.fddy * aa.s;
    dsdx_n_diag = gg.nmeanx2(:,gg.es2) * dsdx_e_diag(gg.es2);
    dsdy_n_diag = gg.nmeany2(:,gg.fs2) * dsdy_f_diag(gg.fs2);
    tau_d_nondim = pp.c60 * aa.H .* sqrt(dsdx_n_diag.^2 + dsdy_n_diag.^2);
    
    % Friction stress (Coulomb limit)
    Nr_diag = max(N_current, pp.N_slide_reg);  % max-based (matches slide_fun_local)
    tau_f_nondim = pp.c61 * pp.mu * Nr_diag;
    
    % Membrane stress (rough estimate: c62 * H * eta * d2u/dx2)
    % Just estimate order of magnitude
    epsxx = gg.nddx(:,gg.es2) * u_inv(gg.es2);
    epsyy = gg.nddy(:,gg.fs2) * v_inv(gg.fs2);
    eta = pp.A_glen^(-1/pp.n_glen) * (pp.eps_reg^2 + epsxx.^2 + epsyy.^2 + epsxx.*epsyy).^((1/pp.n_glen-1)/2);
    tau_m_nondim = pp.c62 * aa.H .* eta .* sqrt(epsxx.^2 + epsyy.^2);
    
    fprintf('   At slow-bias nodes (nondimensional):\n');
    fprintf('     tau_driving:  mean=%.3e, max=%.3e\n', mean(tau_d_nondim(too_low)), max(tau_d_nondim(too_low)));
    fprintf('     tau_friction: mean=%.3e, max=%.3e (Coulomb cap)\n', mean(tau_f_nondim(too_low)), max(tau_f_nondim(too_low)));
    fprintf('     tau_membrane: mean=%.3e, max=%.3e\n', mean(tau_m_nondim(too_low)), max(tau_m_nondim(too_low)));
    fprintf('     tau_d - tau_f (imbalance): mean=%.3e\n', mean(tau_d_nondim(too_low) - tau_f_nondim(too_low)));
    
    if mean(tau_m_nondim(too_low)) < 0.1 * mean(tau_d_nondim(too_low) - tau_f_nondim(too_low))
        fprintf('   ⚠️  Membrane stress << driving-friction imbalance\n');
        fprintf('      The model cannot transfer momentum laterally.\n');
        fprintf('      eps=%.4f is too small. Try eps=0.1 or higher.\n', eps);
    end
end
% Check 4: Boundary effects
fprintf('4. ebdy2: %d edges, fbdy2: %d edges\n', length(gg.ebdy2), length(gg.fbdy2));
if ~isempty(too_low)
    % Check if slow nodes are near boundary
    x_too_low = gg.nx(too_low);
    y_too_low = gg.ny(too_low);
    x_range = [min(gg.nx(nin)), max(gg.nx(nin))];
    y_range = [min(gg.ny(nin)), max(gg.ny(nin))];
    near_bdy = (x_too_low < x_range(1) + 0.1*(x_range(2)-x_range(1))) | ...
               (x_too_low > x_range(2) - 0.1*(x_range(2)-x_range(1))) | ...
               (y_too_low < y_range(1) + 0.1*(y_range(2)-y_range(1))) | ...
               (y_too_low > y_range(2) - 0.1*(y_range(2)-y_range(1)));
    fprintf('   Slow-bias nodes near boundary: %d / %d (%.0f%%)\n', ...
        sum(near_bdy), length(too_low), 100*sum(near_bdy)/length(too_low));
    if sum(near_bdy) > 0.5 * length(too_low)
        fprintf('   ⚠️  Most under-predictions near boundary → BC problem\n');
    end
end

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

    % Implement Dirichlet boundary condition here
    % u(gg.ebdy2) = u_obs(gg.ebdy2);
    % v(gg.fbdy2) = v_obs(gg.fbdy2);

    % --- 2. Velocity misfit ---
    ein = gg.ein2; 
    fin = gg.fin2;
    u0 = opts_inv.u0_reg;
    
    du = obs_mask_e(ein) .* (u(ein) - u_obs(ein)) ./ (abs(u_obs(ein)) + u0);
    dv = obs_mask_f(fin) .* (v(fin) - v_obs(fin)) ./ (abs(v_obs(fin)) + u0);
    
    Jmis = 0.5 * (du'*du + dv'*dv); % residual misfit
    
    % Diagnostic: print speed statistics (only when verbose)
    if opts_inv.verbose
        fast_mask = abs(u_obs(gg.ein2)) > 10*opts_inv.u0_reg;  
        slow_mask = abs(u_obs(gg.ein2)) < opts_inv.u0_reg;
        fprintf('Fast region: %d edges, slow region: %d edges\n', sum(fast_mask), sum(slow_mask));
        U_obs_e = abs(u_obs(ein)) * ps.u_b * pd.ty;  % [m/yr]
        U_mod_e = abs(u(ein)) * ps.u_b * pd.ty;
        fprintf('Speed [m/yr]: obs [%.1f, %.1f], mod [%.1f, %.1f]\n', ...
            min(U_obs_e), max(U_obs_e), min(U_mod_e), max(U_mod_e));
        fprintf('Relative error: mean = %.1f%%, median = %.1f%%\n', ...
            100*mean(abs(du)), 100*median(abs(du)));
        fprintf('Jmis = %.2e, Jreg = %.2e, Jdamp = %.2e\n', Jmis, ...
            0.5*opts_inv.alpha*(c'*(L'*L)*c), 0.5*opts_inv.gamma*sum((c-c_prior).^2));
    end

    % --- 3. Regularisation ---
    Jreg = 0.5 * opts_inv.alpha * (c'*(L'*L)*c); % Tikhonov regularisation on spatial gradients
    Jdamp = 0.5 * opts_inv.gamma * sum((c - c_prior).^2); % damping to prior
    
    J = Jmis + Jreg + Jdamp;
    
    if nargout < 2, return; end  % only J requested
    
    % --- 4. Adjoint RHS: dJ/d[u,v] ---
    % dJ/du_i = (u_i - u_obs_i) / (|u_obs_i| + u0)^2 * mask_i
    dJdu = zeros(gg.eIJ,1);
    dJdv = zeros(gg.fIJ,1);
    dJdu(ein) = obs_mask_e(ein) .* (u(ein) - u_obs(ein)) ./ (abs(u_obs(ein)) + u0).^2;
    dJdv(fin) = obs_mask_f(fin) .* (v(fin) - v_obs(fin)) ./ (abs(v_obs(fin)) + u0).^2;
    
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
    % Export renormalized operators for use in adjoint gradient computation
    ops.emean  = emean;   % renormalized edge←node average [eIJ x nIJ]
    ops.fmean  = fmean;   % renormalized edge←node average [fIJ x nIJ]
    ops.nmeanx = nmeanx;  % renormalized node←x-edge average [nIJ x eIJ]
    ops.nmeany = nmeany;  % renormalized node←y-edge average [nIJ x fIJ]
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
    
    % Term 3: C2 * (1/n-1) * Ur^(1/n-2)
    term3 = (1/n - 1) .* pp.C2 .* Ur.^(1/n - 2);
    
    % Chain rule: dUr/dU = 1 if U >= Ub_reg, else 0 (subgradient of max)
    dslide_dU = (term1 + term2 + term3) .* double(U >= pp.Ub_reg);
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
        (mu.^(-n).*Ur + C.^(-n).*Nr.^n).^(-1/n) + ...
        pp.C2 * Ur.^(1/n-1);
end
