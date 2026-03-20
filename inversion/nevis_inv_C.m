%% NEVIS C-field inversion script
% Invert for the spatial basal friction coefficient C(x,y) using adjoint method
% Given: observed surface velocity (u_obs, v_obs), 
%        ice geometry (H, s, b),
%        effective pressure N, 
%        Coulomb coefficient mu
% Control variable: c = log(C) on nodes

% The sliding law is:
%   tau_b/Ub = N * Ub^(1/n-1) * (mu^(-n)*Ub + C^(-n)*N^n)^(-1/n) + C2*Ub^(1/n-1)
%   C2 is assumed to be 0 for the moment
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
% note: this is a standard case used to initialise the 
oo.casename = 'n2d_region_ice_inversion_test'; 
oo.fn = ['/',oo.casename];
oo.rn = [oo.root,oo.results,oo.fn];
oo.in = [oo.root, 'data/', oo.casename];
oo.dn = [oo.root, 'data/', oo.dataset, '/'];
addpath(oo.code);

% set the surface runoff options
racmo = 0;  % use RACMO runoff data (if 0, use simple sinusoidal function)
distributed = 1;  % if 1, distribute runoff across the domain; if 0, input to discrete moulins

% load saved spinup state
load([oo.in, oo.fn], 'pp','pd','ps','gg','aa','oo');
[pd,ps,pp,oo] = nevis_update_parameters_ice(pd,ps,pp,oo);

% Picard iteration settings for velocity solver
oo.iter_max = 100;      % increase from default 10 — essential for convergence
oo.tol_vel = 1e-6;      % convergence tolerance
oo.display_norms = 0;   % set to 1 for debugging
oo.verb = 0;            % suppress Picard iteration messages from nevis_velocity
oo.display_norms = 0;
oo.boundary_method = 'stress_l_vel_tbl';
oo.mask_boundary_method = 'stress_free';  % allow outlet glaciers to flow freely at mask boundary

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
pp.eps_reg = 1e-1;      % regularisation on strain rates
pp.Ub_reg = 1e-16;      % regularisation on sliding speed (max-based, matches nevis_velocity)
pp.N_slide_reg = 1e-16; % regularisation on effective pressure (max-based, matches nevis_velocity)
pp.taud_reg = 1e-16;    % regularisation on basal shear stress [ may not be needed ? ]
pp.C2 = 0;              % added power-law coefficient in sliding law

%% ============================================================
%  2. Load or generate observed velocity
%  ============================================================
% Option A: load observed velocity from data file
gg = nevis_label_ice_test(gg, oo);
load([oo.dn 'measures_for_nevis_140km.mat']);  % [m/yr]
dd = measures_for_nevis_140km;
u_obs = dd.u_obs_dim / (ps.u_b * pd.ty);  % non-dimensionalise
v_obs = dd.v_obs_dim / (ps.u_b * pd.ty);
u_obs = gg.emean2*u_obs(:);               % project edge vel onto nodes
v_obs = gg.fmean2*v_obs(:);
u_obs(gg.eout2) = NaN;
v_obs(gg.fout2) = NaN;

% Observation mask: only fit where we have valid observations
% Exclude Dirichlet boundary edges — their velocity is prescribed, not controlled by C
obs_mask_e = zeros(gg.eIJ,1);  
obs_mask_e(gg.ein2) = 1;
obs_mask_f = zeros(gg.fIJ,1); 
obs_mask_f(gg.fin2) = 1;
obs_mask_e(gg.ebdy2) = 0;  % exclude Dirichlet x-boundary edges
obs_mask_f(gg.fbdy2) = 0;  % exclude Dirichlet y-boundary edges
fprintf('Excluded %d ebdy + %d fbdy Dirichlet edges from misfit\n', ...
    length(gg.ebdy2), length(gg.fbdy2));

% effective pressure for forward model (use spinup or assume N=1)
oo.initname = 'n2d_region_ice_inversion_test';
init_cond = load(['./data/' oo.initname '/' '0365.mat']); 
                           % load initial condition
vv = init_cond.vv;         % load state variables from the initial 
% N_obs = max(aa.phi_0 - vv.phi, pp.N_slide_reg);
N_obs = aa.phi_0 - vv.phi;
% no noise for real observations
u_obs_noisy = u_obs;
v_obs_noisy = v_obs;

%% ============================================================
%  3. Inversion settings
%  ============================================================
opts_inv.u0_reg = 1e-1;      % velocity scale for relative misfit (dimensionless)

% Continuation schedule (coarse -> fine): start strongly regularized, then relax
% opts_inv.alpha_schedule = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-7];
% opts_inv.gamma_schedule = [1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-10];
% opts_inv.alpha_schedule = [1e-3, 3e-4, 1e-4, 3e-5, 1e-6];
% opts_inv.gamma_schedule = [1e-6, 3e-7, 1e-7, 3e-8, 1e-9];
opts_inv.alpha_schedule = [1e-3, 1e-4, 1e-9, 1e-12];
opts_inv.gamma_schedule = [1e-6, 1e-7, 1e-8, 1e-9];
opts_inv.max_iter_schedule = [25, 5, 5, 5];  % more iters when reg is strong
% Iteration controls (per stage)
opts_inv.max_iter_stage = 50;    % max iterations per continuation stage
opts_inv.max_iter_total = 200;    % safety cap across all stages

% Convergence thresholds (per stage)
opts_inv.J_tol = 1e-3;
opts_inv.dJ_tol = 1e-6;

% Initialize alpha/gamma with first schedule values (needed for initial test call)
opts_inv.alpha = opts_inv.alpha_schedule(1);
opts_inv.gamma = opts_inv.gamma_schedule(1);
opts_inv.verbose = false;
opts_inv.check_grad = false;  % set true to run Taylor test & per-component FD check
opts_inv.speed_misfit = true; % true: misfit on speed |u|; false: misfit on components (u,v)

% History
history.stage = [];
history.iter  = [];
history.J     = [];
history.normg = [];

% initial guess for C (perturbed from truth for synthetic test)
% load(['./data/C_inversion_results.mat'], 'C_hat_dim');
% C_init = 0.01*C_hat_dim * (ps.u_b^(1/pp.n_slide) / ps.tau);
C_init = ones(gg.nIJ,1);  % uniform prior
c_prior = log(C_init);

% build gradient operator L for regularisation (using grid finite differences)
% L * c gives spatial gradient of c (defined on nodes)
Lx = gg.eddx(gg.ein2, :);  % x-gradient at interior x-edges [length(ein2) x nIJ]
Ly = gg.fddy(gg.fin2, :);  % y-gradient at interior y-edges [length(fin2) x nIJ]
L = [Lx; Ly];

%% ============================================================
%  4. Run iterative C-N coupled inversion
%  ============================================================
fprintf('\n========================================\n');
fprintf('Starting iterative C-N coupled inversion\n');
fprintf('========================================\n\n');

c0 = c_prior;
C_prev = exp(c0);
N_current = N_obs;   % Initial N from input
vv_hydro = [];       % hydrology state (empty = fresh start on first iteration)

% DEBUG: test objective at initial point
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


% Replace NaN in N_current (softplus Reg_Ni handles N<0 gracefully in sliding law)
N_current(isnan(N_current)) = 0;

fprintf('  After masking: obs_mask_e sum = %d, obs_mask_f sum = %d\n', ...
    sum(obs_mask_e), sum(obs_mask_f));

% --- Apply adaptive regularization AFTER NaN cleaning ---
% Compute velocity gradient magnitude at edges using cleaned observations
dudx_e = gg.eddx * (gg.nmeanx2(:,gg.es2) * u_obs_noisy(gg.es2));  % du/dx at x-edges
dvdy_f = gg.fddy * (gg.nmeany2(:,gg.fs2) * v_obs_noisy(gg.fs2));  % dv/dy at y-edges
grad_mag_e = abs(dudx_e(gg.ein2));  % velocity gradient at x-edges
grad_mag_f = abs(dvdy_f(gg.fin2));  % velocity gradient at y-edges

% Adaptive weight: w = 1 / (1 + beta * |grad(u)|/mean(|grad(u)|))
% High gradient -> small w -> less regularization
beta_adaptive = 0.1;  % tuning parameter (increase to reduce reg more in shear zones)
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

[J0, g0] = obj_fun(c0);
fprintf('  J(c0) = %.6e, norm(g) = %.6e\n', J0, norm(g0));
if isnan(J0) || isinf(J0)
    error('Objective still undefined. Check forward model.');
end

%% check the gradient calculator
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
    err1_prev = NaN; 
    err2_prev = NaN;
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
C_hat_dim = C_hat * (ps.u_b^(1/pp.n_slide) / ps.tau);  % dimensional C [Pa s/m]
fprintf('\nInversion complete: final J=%.6e, exitflag=%d\n', J_hat, exitflag);

%% save the relelvant fields so that we can reconstruct the C field with arbitrary C1/C2 partitioning in the post-processing step without rerunning the forward model
save('./data/C_inversion_results.mat', 'C_hat', 'u_obs_noisy', 'v_obs_noisy', 'N_current', 'aa', 'pp', 'gg', 'oo', 'ps', 'exitflag');
% save('./data/C_inversion_C2_results.mat', 'c_hat', 'C1_hat', 'C2_hat', 'C1_hat_dim', 'C2_hat_dim', 'history', 'opts_inv', 'J_hat', 'exitflag');

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

% Node-based speed
Un_obs = sqrt((gg.nmeanx2(:,gg.es2)*u_obs(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_obs(gg.fs2)).^2);
Un_inv = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);
Un_misfit = Un_inv - Un_obs;
Un_misfit(gg.nout) = NaN;

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

% Merge section 5 and 5b into one compact diagnostic figure
cmap_div = interp1([0 0.5 1], [0 0 1; 1 1 1; 1 0 0], linspace(0,1,256));
fig = figure('Name', 'Inversion Diagnostics', 'Position', [60 60 1700 880], 'Color', 'w');
tlo = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% (1) Prior C
ax = nexttile(tlo, 1);
zC_prior = reshape(log10(C_prior * (ps.tau / ps.u_b^(1/pp.n_slide))), gg.nI, gg.nJ);
zC_prior(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, zC_prior); shading(ax, 'flat'); colorbar(ax);
title(ax, 'Prior log_{10}(C)'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');

% (2) Inverted C
ax = nexttile(tlo, 2);
z = reshape(log10(C_hat_dim), gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z); shading(ax, 'flat'); colorbar(ax);
title(ax, 'Inverted log_{10}(C)'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');

% (3) Observed speed
ax = nexttile(tlo, 3);
z = reshape(Un_obs_dim, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z); shading(ax, 'flat'); colorbar(ax);
title(ax, 'Observed speed [m/yr]'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');
clim(ax, [0 300]);

% (4) Modeled speed
ax = nexttile(tlo, 4);
z = reshape(Un_inv_dim, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z); shading(ax, 'flat'); colorbar(ax);
title(ax, 'Modeled speed [m/yr]'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');
clim(ax, [0 300]);

% (5) Misfit
ax = nexttile(tlo, 5);
z = reshape(Un_misfit_dim, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z); shading(ax, 'flat'); colorbar(ax);
title(ax, 'Misfit [m/yr]'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');
caxis(ax, [-100 100]); colormap(ax, cmap_div);

% (6) Relative error
ax = nexttile(tlo, 6);
z = reshape(100*rel_err, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z); shading(ax, 'flat'); colorbar(ax);
title(ax, 'Relative error [%]'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');
caxis(ax, [0 50]);

% (7) N field
ax = nexttile(tlo, 7);
z = reshape(N_current * ps.phi, gg.nI, gg.nJ) / 1e6; z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z); shading(ax, 'flat'); colorbar(ax);
title(ax, 'N [MPa]'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');

% (8) Speed misfit highlight (same variable, tighter color range)
ax = nexttile(tlo, 8);
z = reshape(Un_misfit_dim, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z); shading(ax, 'flat'); colorbar(ax);
title(ax, 'Misfit (focus) [m/yr]'); xlabel(ax, 'x'); ylabel(ax, 'y'); axis(ax, 'equal'); axis(ax, 'tight');
caxis(ax, [-50 50]); colormap(ax, cmap_div);

% (9) Scatter: obs speed vs C
ax = nexttile(tlo, 9);
nin = gg.nin2;
scatter(ax, Un_obs_dim(nin), log10(C_hat_dim(nin)), 6, Un_misfit_dim(nin), 'filled');
xlabel(ax, 'Observed speed [m/yr]'); ylabel(ax, 'log_{10}(C_{dim})');
title(ax, 'C vs speed (color=misfit)'); colorbar(ax); caxis(ax, [-50 50]); colormap(ax, cmap_div);
grid(ax, 'on'); box(ax, 'on');

% (10) Scatter: N vs misfit
ax = nexttile(tlo, 10);
scatter(ax, N_current(nin)*ps.phi/1e6, Un_misfit_dim(nin), 6, log10(C_hat_dim(nin)), 'filled');
xlabel(ax, 'N [MPa]'); ylabel(ax, 'Speed misfit [m/yr]');
title(ax, 'N vs misfit (color=log_{10}C)'); colorbar(ax);
grid(ax, 'on'); box(ax, 'on');

% (11) and (12): summary text panel spanning two tiles
ax = nexttile(tlo, 11, [1 2]);
axis(ax, 'off');
rmse_speed = sqrt(mean((Un_inv_dim(nin) - Un_obs_dim(nin)).^2, 'omitnan'));
txt = sprintf(['Diagnostics summary\n', ...
    'RMSE speed: %.2f m/yr\n', ...
    'C_{dim} range: [%.2e, %.2e]\n', ...
    'N range: [%.2f, %.2f] MPa\n', ...
    'Too slow nodes (<-20 m/yr): %d\n', ...
    'Too fast nodes (>20 m/yr): %d'], ...
    rmse_speed, min(C_hat_dim(nin)), max(C_hat_dim(nin)), ...
    min(N_current(nin)*ps.phi/1e6), max(N_current(nin)*ps.phi/1e6), ...
    numel(too_low), numel(too_high));
text(ax, 0.02, 0.9, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 11, 'FontName', 'Helvetica');

set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 10);
drawnow;

%% ============================================================
%  5c. Basal friction distribution
%  ============================================================
% % Compute tau_b on nodes using the same sliding law as nevis_velocity (slide_fun)
% U_nodes_fric = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + ...
%                     (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);  % [nIJ x 1]
% Ub_fric  = max(U_nodes_fric, pp.Ub_reg);      % regularised speed
% Nr_fric  = max(N_current,    pp.N_slide_reg);  % regularised N
% C_fric   = C_hat;                             % inverted C
% mu_fric  = pp.mu * ones(gg.nIJ, 1);           % Coulomb coefficient
% 
% % tau_b / Ub  (dimensionless)
% taub_over_Ub = Nr_fric .* Ub_fric.^(1/pp.n_slide-1) .* ...
%     (mu_fric.^(-pp.n_slide) .* Ub_fric + C_fric.^(-pp.n_slide) .* Nr_fric.^pp.n_slide).^(-1/pp.n_slide) ...
%     + pp.C2 * Ub_fric.^(1/pp.n_slide-1);
% 
% % tau_b  [dimensional: Pa]
% tau_b_dim = pp.c61 * taub_over_Ub .* U_nodes_fric * ps.phi;   % [Pa]
% tau_b_kPa = tau_b_dim / 1e3;                                   % [kPa]
% 
% % Driving stress [kPa] (for comparison)
% dsdx_fric = gg.nmeanx2(:,gg.es2) * (gg.eddx(gg.es2,:) * aa.s);
% dsdy_fric = gg.nmeany2(:,gg.fs2) * (gg.fddy(gg.fs2,:) * aa.s);
% tau_d_kPa = pp.c60 * aa.H .* sqrt(dsdx_fric.^2 + dsdy_fric.^2) * ps.phi / 1e3;  % [kPa]
% 
% figure('Name','Basal Friction Distribution','Position',[100 100 1600 400]);
% 
% % (a) Basal friction tau_b [kPa]
% subplot(1,4,1);
% z = reshape(tau_b_kPa, gg.nI, gg.nJ); z(gg.nout) = NaN;
% pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
% clim([0 200]);
% title('\tau_b [kPa]'); xlabel('x'); ylabel('y'); axis equal tight;
% 
% % (b) Driving stress tau_d [kPa]
% subplot(1,4,2);
% z = reshape(tau_d_kPa, gg.nI, gg.nJ); z(gg.nout) = NaN;
% pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
% clim([0 200]);
% title('\tau_d [kPa]'); xlabel('x'); ylabel('y'); axis equal tight;
% 
% % (c) Friction fraction tau_b / tau_d
% subplot(1,4,3);
% fric_frac = tau_b_kPa ./ max(tau_d_kPa, 0.1);  % avoid div by zero
% fric_frac(gg.nout) = NaN;
% z = reshape(fric_frac, gg.nI, gg.nJ); z(gg.nout) = NaN;
% pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
% clim([0 2]);
% title('\tau_b / \tau_d'); xlabel('x'); ylabel('y'); axis equal tight;
% 
% % (d) Coulomb cap mu*N [kPa]  (maximum possible friction)
% subplot(1,4,4);
% coulomb_cap_kPa = pp.mu * Nr_fric * ps.phi / 1e3;  % [kPa]
% z = reshape(coulomb_cap_kPa, gg.nI, gg.nJ); z(gg.nout) = NaN;
% pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
% clim([0 200]);
% title('\mu N [kPa] (Coulomb cap)'); xlabel('x'); ylabel('y'); axis equal tight;
% 
% colormap(jet);
% drawnow;

fprintf('\n=== Basal Friction Summary ===\n');
nin2 = gg.nin2;
fprintf('  tau_b [kPa]:  mean=%.1f, median=%.1f, max=%.1f\n', ...
    mean(tau_b_kPa(nin2)), median(tau_b_kPa(nin2)), max(tau_b_kPa(nin2)));
fprintf('  tau_d [kPa]:  mean=%.1f, median=%.1f, max=%.1f\n', ...
    mean(tau_d_kPa(nin2)), median(tau_d_kPa(nin2)), max(tau_d_kPa(nin2)));
fprintf('  tau_b/tau_d:  mean=%.2f, median=%.2f\n', ...
    mean(fric_frac(nin2)), median(fric_frac(nin2)));
fprintf('  mu*N [kPa]:   mean=%.1f, median=%.1f\n', ...
    mean(coulomb_cap_kPa(nin2)), median(coulomb_cap_kPa(nin2)));
fprintf('  Nodes at Coulomb limit (tau_b > 0.9*mu*N): %d / %d\n', ...
    sum(tau_b_kPa(nin2) > 0.9*coulomb_cap_kPa(nin2)), length(nin2));

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
