function [inv, vv_hydro, summary] = nevis_inv_C(config_file)
%NEVIS_INV_C Invert a regional NEVIS C field from an inversion config.
%
% Usage:
%   [inv, vv_hydro, summary] = nevis_inv_C('./cases/templates/nevis_inversion.m');
% Prefer nevis_run_inversion(config_file) for user-facing workflows where
% cfg.sliding_law selects the inversion backend.
%
% The config declares all data inputs and output locations. Missing files are
% errors; this function does not fall back to legacy ./data files.
%
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

if nargin < 1 || isempty(config_file)
    error('nevis_inv_C:MissingConfig', ...
        'Call nevis_inv_C with an inversion config file.');
end

format compact;
cfg = load_inversion_config(config_file);
root = cfg.root;
code = cfg.code;
addpath(code);
addpath(fullfile(root, 'inversion'));

%% ============================================================
%  1. Load spinup results and set up NEVIS framework
%  ============================================================
state_file = resolve_project_path(cfg.source.state_file, cfg);
if exist(state_file, 'file') ~= 2
    error('nevis_inv_C:MissingSourceState', ...
        'Configured inversion source state not found: %s', state_file);
end
state = load(state_file, 'pp', 'pd', 'ps', 'gg', 'aa', 'oo');
required_state = {'pp', 'pd', 'ps', 'gg', 'aa', 'oo'};
for k_req = 1:numel(required_state)
    if ~isfield(state, required_state{k_req})
        error('nevis_inv_C:MissingSourceVariable', ...
            'Source state "%s" must contain variable %s.', ...
            state_file, required_state{k_req});
    end
end
pp = state.pp;
pd = state.pd;
ps = state.ps;
gg = state.gg;
aa = state.aa;
oo = state.oo;

oo.root = root;
oo.code = code;
oo.results = cfg.results;
oo.dataset = cfg.dataset.name;
oo.dataset_root = cfg.dataset.root;
oo.dn = cfg.dataset.root;
oo.casename = cfg.source.casename;
oo.fn = ['/', oo.casename];
oo.rn = fullfile(oo.root, oo.results, oo.casename);

% set the surface runoff options
racmo = cfg.forward_hydrology.racmo_runoff;
distributed = cfg.forward_hydrology.distributed_input;
partition_ratio = cfg.partition_ratio;
[pd,ps,pp,oo] = nevis_update_parameters_ice(pd,ps,pp,oo);

% Picard iteration settings for velocity solver
oo.iter_max = cfg.solver.iter_max;
oo.tol_vel = cfg.solver.tol_vel;
oo.display_norms = cfg.solver.display_norms;
oo.verb = cfg.solver.verb;
oo.boundary_method = cfg.solver.boundary_method;
oo.mask_boundary_method = cfg.solver.mask_boundary_method;
oo.partition_ratio = partition_ratio;

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

pp.eps_reg = cfg.regularization.eps_reg;
pp.Ub_reg = cfg.regularization.Ub_reg;
pp.N_slide_reg = cfg.regularization.N_slide_reg;
pp.taud_reg = cfg.regularization.taud_reg;
pp.C2 = cfg.regularization.C2;

%% ============================================================
%  2. Load or generate observed velocity
%  ============================================================
% Option A: load observed velocity from data file
gg = nevis_label_ice_test(gg, oo);
velocity_file = resolve_dataset_path(cfg.velocity.file, cfg);
if exist(velocity_file, 'file') ~= 2
    error('nevis_inv_C:MissingVelocityFile', ...
        'Configured velocity file not found: %s', velocity_file);
end
velocity_data = load(velocity_file);
velocity_var = cfg.velocity.variable;
if ~isfield(velocity_data, velocity_var)
    error('nevis_inv_C:MissingVelocityVariable', ...
        'Velocity file "%s" does not contain variable "%s".', ...
        velocity_file, velocity_var);
end
dd = velocity_data.(velocity_var);
u_field = cfg.velocity.u_field;
v_field = cfg.velocity.v_field;
if ~isfield(dd, u_field) || ~isfield(dd, v_field)
    error('nevis_inv_C:MissingVelocityFields', ...
        'Velocity variable must contain fields "%s" and "%s".', ...
        u_field, v_field);
end
u_obs = dd.(u_field) / (ps.u_b * pd.ty);  % non-dimensionalise
v_obs = dd.(v_field) / (ps.u_b * pd.ty);
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

% Initial effective pressure for the inversion.
[N_obs, initial_hydrology_source] = nevis_inversion_initial_effective_pressure( ...
    cfg.initial_hydrology, aa, gg, cfg);
% no noise for real observations
u_obs_noisy = u_obs;
v_obs_noisy = v_obs;

%% ============================================================
%  3. Inversion settings
%  ============================================================
opts_inv = cfg.opts_inv;
opts_inv.alpha = opts_inv.alpha_schedule(1);
opts_inv.gamma = opts_inv.gamma_schedule(1);

% History
history.stage = [];
history.iter  = [];
history.J     = [];
history.normg = [];

C_init = build_initial_C(cfg.prior, gg);
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
beta_adaptive = cfg.regularization.adaptive_beta;
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
C_hat_dim = C_hat * (ps.tau / ps.u_b^(1/pp.n_slide));  % dimensional C [Pa s/m]
fprintf('\nInversion complete: final J=%.6e, exitflag=%d\n', J_hat, exitflag);

%% save the relelvant fields so that we can reconstruct the C field with arbitrary C1/C2 partitioning in the post-processing step without rerunning the forward model
inversion_file = resolve_dataset_path(cfg.output.inversion_file, cfg);
hydrology_file = resolve_dataset_path(cfg.output.initial_hydrology_file, cfg);
ensure_parent_dir(inversion_file);
ensure_parent_dir(hydrology_file);
save(inversion_file, 'C_hat', 'C_hat_dim', 'u_obs_noisy', 'v_obs_noisy', ...
    'N_current', 'aa', 'pp', 'gg', 'oo', 'ps', 'exitflag', ...
    'partition_ratio', 'history', 'opts_inv', 'J_hat', 'cfg', ...
    'initial_hydrology_source');
% save('./data/C_inversion_C2_results.mat', 'c_hat', 'C1_hat', 'C2_hat', 'C1_hat_dim', 'C2_hat_dim', 'history', 'opts_inv', 'J_hat', 'exitflag');
fprintf('Saved inversion result to %s\n', inversion_file);

%% run a forward case with the inverted C so it can be reused later
% Inversion uses mu=inf (pure Weertman); forward uses mu=0.25 (Coulomb+Weertman)
% with C partitioned into C1 (Coulomb) and C2 (power-law) via partition_ratio
mu_s_fwd = cfg.forward_hydrology.mu_s;
fprintf('Preparing reusable hydrology state (mu_s=%.2f, partition_ratio=%.2f).\n', mu_s_fwd, partition_ratio);
aa_fwd = aa;
aa_fwd.C = C_hat;
u_fwd = u_obs_noisy;
v_fwd = v_obs_noisy;
u_fwd(isnan(u_fwd) & ~ismember((1:gg.eIJ)', gg.eout2)) = 0;
v_fwd(isnan(v_fwd) & ~ismember((1:gg.fIJ)', gg.fout2)) = 0;
[u_fwd, v_fwd] = nevis_velocity(aa.H, u_fwd, v_fwd, N_current, aa_fwd, pp, gg, oo);

inv = struct('C_hat', C_hat, ...
             'C_hat_dim', C_hat_dim, ...
             'u_obs_noisy', u_obs_noisy, ...
             'v_obs_noisy', v_obs_noisy, ...
             'N_current', N_current, ...
             'aa', aa, ...
             'pp', pp, ...
             'gg', gg, ...
             'oo', oo, ...
             'ps', ps, ...
             'partition_ratio', partition_ratio);

vv_hydro = struct('u', u_fwd, 'v', v_fwd, 'N', N_current);
if cfg.forward_hydrology.enabled
    output_casename = cfg.forward_hydrology.output_casename;
    [N_forward, vv_hydro] = nevis_run_fwd_hydrology(inv, vv_hydro, ...
        mu_s_fwd, racmo, distributed, output_casename);
    vv_hydro.N = N_forward;
end
save(hydrology_file, 'vv_hydro', 'C_hat', 'C_hat_dim', 'inv', 'cfg');
fprintf('Saved reusable initial hydrology to %s\n', hydrology_file);

summary = struct;
summary.config_file = cfg.config_file;
summary.inversion_file = inversion_file;
summary.initial_hydrology_file = hydrology_file;
summary.final_J = J_hat;
summary.exitflag = exitflag;
summary.partition_ratio = partition_ratio;

%% ============================================================
%  5. Plot results
%  ============================================================
if cfg.plot.enabled
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

function cfg = load_inversion_config(config_file)
    if isstruct(config_file)
        cfg = config_file;
        if ~isfield(cfg, 'config_file')
            cfg.config_file = '<struct>';
        end
        validate_inversion_config(cfg);
        return
    end

    if exist(config_file, 'file') ~= 2
        error('nevis_inv_C:MissingConfig', ...
            'Inversion config not found: %s', config_file);
    end

    [~, ~, ext] = fileparts(config_file);
    switch lower(ext)
        case '.m'
            run(config_file);
            if exist('cfg', 'var') ~= 1
                error('nevis_inv_C:MissingCfg', ...
                    'M-file inversion config must assign variable cfg.');
            end
        case '.mat'
            data = load(config_file);
            if ~isfield(data, 'cfg')
                error('nevis_inv_C:MissingCfg', ...
                    'MAT inversion config must contain variable cfg.');
            end
            cfg = data.cfg;
        case '.json'
            cfg = jsondecode(fileread(config_file));
        otherwise
            error('nevis_inv_C:UnsupportedConfig', ...
                'Unsupported inversion config extension: %s', ext);
    end

    if ~isfield(cfg, 'config_file') || isempty(cfg.config_file)
        cfg.config_file = config_file;
    end
    validate_inversion_config(cfg);
end

function validate_inversion_config(cfg)
    required = {'name', 'root', 'code', 'results', 'dataset', 'source', ...
        'velocity', 'initial_hydrology', 'output', 'solver', ...
        'regularization', 'opts_inv', 'prior', 'forward_hydrology', ...
        'plot', 'partition_ratio'};
    for i = 1:numel(required)
        if ~isfield(cfg, required{i}) || isempty(cfg.(required{i}))
            error('nevis_inv_C:MissingConfigField', ...
                'Missing inversion config field cfg.%s.', required{i});
        end
    end
    require_fields(cfg.dataset, {'name', 'root', 'manifest'}, 'cfg.dataset');
    require_fields(cfg.source, {'casename', 'state_file'}, 'cfg.source');
    require_fields(cfg.velocity, {'file', 'variable', 'u_field', 'v_field'}, 'cfg.velocity');
    require_initial_hydrology_fields(cfg.initial_hydrology);
    require_fields(cfg.output, {'inversion_file', 'initial_hydrology_file'}, 'cfg.output');
    require_fields(cfg.solver, {'iter_max', 'tol_vel', 'display_norms', ...
        'verb', 'boundary_method', 'mask_boundary_method'}, 'cfg.solver');
    require_fields(cfg.regularization, {'eps_reg', 'Ub_reg', ...
        'N_slide_reg', 'taud_reg', 'C2', 'adaptive_beta'}, ...
        'cfg.regularization');
    require_fields(cfg.opts_inv, {'u0_reg', 'alpha_schedule', ...
        'gamma_schedule', 'max_iter_schedule', 'max_iter_stage', ...
        'max_iter_total', 'J_tol', 'dJ_tol', 'verbose', 'check_grad', ...
        'speed_misfit'}, 'cfg.opts_inv');
    require_fields(cfg.prior, {'mode'}, 'cfg.prior');
    require_fields(cfg.forward_hydrology, {'enabled', 'mu_s', ...
        'racmo_runoff', 'distributed_input', 'output_casename'}, ...
        'cfg.forward_hydrology');
    require_fields(cfg.plot, {'enabled'}, 'cfg.plot');

    if numel(cfg.opts_inv.alpha_schedule) ~= numel(cfg.opts_inv.gamma_schedule)
        error('nevis_inv_C:ScheduleSizeMismatch', ...
            'cfg.opts_inv.alpha_schedule and gamma_schedule must have the same length.');
    end
    if isempty(cfg.opts_inv.max_iter_schedule)
        error('nevis_inv_C:MissingConfigField', ...
            'cfg.opts_inv.max_iter_schedule must not be empty.');
    end
end

function require_initial_hydrology_fields(initial_hydrology)
    require_fields(initial_hydrology, {'mode'}, 'cfg.initial_hydrology');
    mode = lower(initial_hydrology.mode);
    switch mode
        case 'file'
            require_fields(initial_hydrology, {'file', 'variable'}, ...
                'cfg.initial_hydrology');
        case 'k_factor'
            require_fields(initial_hydrology, {'k_factor'}, ...
                'cfg.initial_hydrology');
        case {'constant_n', 'uniform_n'}
            require_fields(initial_hydrology, {'N'}, 'cfg.initial_hydrology');
        otherwise
            error('nevis_inv_C:UnsupportedInitialHydrologyMode', ...
                ['Unsupported cfg.initial_hydrology.mode "%s". Use "file", ' ...
                 '"k_factor", or "constant_N".'], initial_hydrology.mode);
    end
end

function require_fields(s, names, label)
    for i = 1:numel(names)
        if ~isfield(s, names{i}) || isempty(s.(names{i}))
            error('nevis_inv_C:MissingConfigField', ...
                'Missing inversion config field %s.%s.', label, names{i});
        end
    end
end

function C_init = build_initial_C(prior, gg)
    switch lower(prior.mode)
        case 'uniform'
            require_fields(prior, {'C'}, 'cfg.prior');
            C_init = prior.C * ones(gg.nIJ, 1);
        otherwise
            error('nevis_inv_C:UnsupportedPrior', ...
                'Unsupported cfg.prior.mode: %s', prior.mode);
    end
end

function path_out = resolve_dataset_path(path_in, cfg)
    if isempty(path_in)
        path_out = '';
        return
    end
    if is_absolute_path(path_in)
        path_out = path_in;
    else
        path_out = fullfile(cfg.dataset.root, path_in);
    end
    root_canon = canonical_path(cfg.dataset.root);
    path_canon = canonical_path(path_out);
    if ~startsWith([path_canon filesep], [root_canon filesep])
        error('nevis_inv_C:PathOutsideDataset', ...
            'Configured dataset path is outside cfg.dataset.root: %s', path_out);
    end
end

function path_out = resolve_project_path(path_in, cfg)
    if isempty(path_in)
        path_out = '';
        return
    end
    if is_absolute_path(path_in)
        path_out = path_in;
    else
        path_out = fullfile(cfg.root, path_in);
    end
end

function ensure_parent_dir(path_in)
    parent = fileparts(path_in);
    if ~isempty(parent) && exist(parent, 'dir') ~= 7
        mkdir(parent);
    end
end

function tf = is_absolute_path(path_in)
    tf = startsWith(path_in, filesep) || ...
        ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]', 'once'));
end

function path_out = canonical_path(path_in)
    if exist(path_in, 'file') || exist(path_in, 'dir')
        path_out = char(java.io.File(path_in).getCanonicalPath());
    else
        [parent, name, ext] = fileparts(path_in);
        if isempty(parent)
            parent = '.';
        end
        path_out = fullfile(char(java.io.File(parent).getCanonicalPath()), [name ext]);
    end
end
