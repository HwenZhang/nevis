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
oo.casename = 'n2d_region_ice_meanperms1_Hreg1000_kappa1e_10_mu2e1_spinup';
oo.fn = ['/',oo.casename];
oo.rn = [oo.root,oo.results,oo.fn];
oo.dn = [oo.root, 'data/', oo.dataset, '/'];
addpath(oo.code);

% load saved spinup state
load([oo.rn, oo.fn], 'pp','pd','ps','gg','aa','oo');

pd.mu_s = 0.5;
pd.n_glen = 1;
eps = 0.01; % ratio of membrane stress to driving stress, used to set the ice softness to ensure that membrane stresses are small in this test case, eps rises, A decreases, and viscosity increases
pd.A_glen = 1/2/((eps)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make membrane stress terms of dimensionless size eps in momentum equation
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo);

% pp.mu = 0.01;
% pp.n_slide = 3;
% add regularisation parameters to pp if not already present (for backward compatibility)
if ~isfield(pp,'eps_reg'), pp.eps_reg = 1e-16; end % regularisation on strain rates
if ~isfield(pp,'Ub_reg'), pp.Ub_reg = 1e-16; end % regularisation on sliding speed
if ~isfield(pp,'N_slide_reg'), pp.N_slide_reg = 1e-16; end % regularisation on effective pressure in sliding law
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
oo.initname = 'n2d_region_ice_meanperms1_Hreg1000_kappa1e_10_mu2e1_bdtest';
init_cond = load(['./results/' oo.initname '/' '0365.mat']); % load initial condition
vv = init_cond.vv;                                % load state variables from the initial 
% N_obs = ones(gg.nIJ,1);
N_obs = aa.phi_0 - vv.phi;

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
opts_inv.u0_reg = 1e-2;      % velocity scale for relative misfit (dimensionless)

% Continuation schedule (coarse -> fine): start strongly regularized, then relax
opts_inv.alpha_schedule = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6];
opts_inv.gamma_schedule = [1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9];
opts_inv.alpha_schedule = [1e-3, 3e-4, 1e-4, 3e-5, 1e-5];
opts_inv.gamma_schedule = [1e-6, 3e-7, 1e-7, 3e-8, 1e-8];
% Iteration controls (per stage)
opts_inv.max_iter_stage = 50;     % max iterations per continuation stage
opts_inv.max_iter_total = 300;    % safety cap across all stages

% Convergence thresholds (per stage)
opts_inv.J_tol = 1e-3;
opts_inv.dJ_tol = 1e-6;

% Initialize alpha/gamma with first schedule values (needed for initial test call)
opts_inv.alpha = opts_inv.alpha_schedule(1);
opts_inv.gamma = opts_inv.gamma_schedule(1);

% History
history.stage = [];
history.iter  = [];
history.J     = [];
history.normg = [];

% initial guess for C (perturbed from truth for synthetic test)
C_init = 1 * ones(gg.nIJ,1);  % uniform prior
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
L = [Lx; Ly];  % size: [length(ein2)+length(fin2), nIJ]

% nevis_plot_grid_ice(gg); return; 

%% ============================================================
%  4. Run inversion with fminunc
%  ============================================================
c0 = c_prior;

obj_fun = @(c) objective_and_grad(c, u_obs_noisy, v_obs_noisy, ...
    obs_mask_e, obs_mask_f, N_obs, aa, pp, gg, oo, L, c_prior, opts_inv);

% DEBUG: test objective at initial point
fprintf('Testing objective at initial point...\n');
fprintf('  NaN in u_obs(ein2): %d / %d\n', sum(isnan(u_obs_noisy(gg.ein2))), length(gg.ein2));
fprintf('  NaN in v_obs(fin2): %d / %d\n', sum(isnan(v_obs_noisy(gg.fin2))), length(gg.fin2));
fprintf('  NaN in N_obs: %d / %d\n', sum(isnan(N_obs)), length(N_obs));
fprintf('  NaN in c0: %d / %d\n', sum(isnan(c0)), length(c0));

% Mask out NaN observations
obs_mask_e(isnan(u_obs_noisy)) = 0;
obs_mask_f(isnan(v_obs_noisy)) = 0;
% Replace NaN in obs with 0 (masked anyway)
u_obs_noisy(isnan(u_obs_noisy)) = 0;
v_obs_noisy(isnan(v_obs_noisy)) = 0;
% Replace NaN/negative N_obs with small positive value (N must be > 0 for sliding law)
N_obs(isnan(N_obs) | N_obs < 0) = pp.N_slide_reg;

fprintf('  After masking: obs_mask_e sum = %d, obs_mask_f sum = %d\n', ...
    sum(obs_mask_e), sum(obs_mask_f));

% Recreate objective with cleaned data
obj_fun = @(c) objective_and_grad(c, u_obs_noisy, v_obs_noisy, ...
    obs_mask_e, obs_mask_f, N_obs, aa, pp, gg, oo, L, c_prior, opts_inv);

[J0, g0] = obj_fun(c0);
fprintf('  J(c0) = %.6e, norm(g) = %.6e\n', J0, norm(g0));
if isnan(J0) || isinf(J0)
    error('Objective still undefined. Check forward model.');
end

fprintf('Starting C-field inversion (continuation)...\n');

c_stage = c0;
iter_total = 0;
J_hat = NaN; exitflag = NaN; output = struct();

for k = 1:min(numel(opts_inv.alpha_schedule), numel(opts_inv.gamma_schedule))
    opts_inv.alpha = opts_inv.alpha_schedule(k);
    opts_inv.gamma = opts_inv.gamma_schedule(k);

    fprintf('\n--- Stage %d/%d: alpha=%.3e, gamma=%.3e ---\n', k, ...
        min(numel(opts_inv.alpha_schedule), numel(opts_inv.gamma_schedule)), opts_inv.alpha, opts_inv.gamma);

    % Recreate objective with current regularization weights
    obj_fun = @(c) objective_and_grad(c, u_obs_noisy, v_obs_noisy, ...
        obs_mask_e, obs_mask_f, N_obs, aa, pp, gg, oo, L, c_prior, opts_inv);
    opts = optimoptions('fminunc', ...
        'Algorithm', 'quasi-newton', ...
        'HessianApproximation', 'lbfgs', ...
        'SpecifyObjectiveGradient', true, ...
        'Display', 'iter', ...
        'MaxIterations', opts_inv.max_iter_stage, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-10);

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

% Compute final C and print
C_hat = exp(c_hat);
fprintf('Inversion complete. Final J = %.6e, exit flag = %d\n', J_hat, exitflag);

% Dimensionalise C for interpretation, save results
C_hat_dim = C_hat * (ps.tau / ps.u_b^(1/pp.n_slide));
save('./data/C_inversion_results.mat', 'C_hat_dim', 'c_hat', 'history', 'opts_inv', 'J_hat', 'exitflag');

% Plot lightweight convergence diagnostics
figure('Name','Convergence (lightweight)','Position',[100 100 900 350]);
subplot(1,2,1);
semilogy(history.iter, history.J, 'o-','LineWidth',1.5);
xlabel('Total iterations'); ylabel('J'); grid on;
title('Cost by continuation stage');

subplot(1,2,2);
semilogy(history.iter, history.normg, 's-','LineWidth',1.5);
xlabel('Total iterations'); ylabel('firstorderopt'); grid on;
title('Optimality measure');

%% ============================================================
%  5. Plot results
%  ============================================================
C_prior = C_init;
% C_prior = exp(C_prior);
% compute final velocity
aa_inv = aa;
aa_inv.C = C_hat;
u_inv = zeros(gg.eIJ,1); v_inv = zeros(gg.fIJ,1);
u_inv(gg.eout2) = NaN; v_inv(gg.fout2) = NaN;
[u_inv, v_inv] = nevis_velocity(aa.H, u_inv, v_inv, N_obs, aa_inv, pp, gg, oo);

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
%  Objective function and adjoint gradient
%  ============================================================
function [J, g] = objective_and_grad(c, u_obs, v_obs, ...
    obs_mask_e, obs_mask_f, N, aa, pp, gg, oo, L, c_prior, opts_inv)
    % Evaluate cost function J(c) and gradient dJ/dc using adjoint method
    % c = log(C) on all nIJ nodes
    
    C = exp(c);
    aa_loc = aa;
    aa_loc.C = C;
    
    % --- 1. Forward solve ---
    u = zeros(gg.eIJ,1); 
    v = zeros(gg.fIJ,1);
    u(gg.eout2) = NaN; 
    v(gg.fout2) = NaN;
    [u, v] = nevis_velocity(aa.H, u, v, N, aa_loc, pp, gg, oo);
    
    % --- 2. Velocity misfit ---
    ein = gg.ein2; 
    fin = gg.fin2;
    u0 = opts_inv.u0_reg;
    
    du = obs_mask_e(ein) .* (u(ein) - u_obs(ein)) ./ (abs(u_obs(ein)) + u0);
    dv = obs_mask_f(fin) .* (v(fin) - v_obs(fin)) ./ (abs(v_obs(fin)) + u0);
    
    Jmis = 0.5 * (du'*du + dv'*dv); % residual misfit
    
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
    
    % --- 5. Assemble Jacobian (same matrix as in nevis_velocity's linear system) ---
    % We need the final-iteration Jacobian: A * [u;v] = rhs
    % The adjoint equation is: A^T * lambda = -dJ/d[u,v]
    [A, ~] = assemble_velocity_jacobian(u, v, N, aa_loc, pp, gg, oo);
    
    rhs_adj = -[dJdu(ein); dJdv(fin)];
    lambda = A' \ rhs_adj;
    lambda_u = lambda(1:length(ein));
    lambda_v = lambda(length(ein)+1:end);
    
    % --- 6. Gradient dJ/dC from friction term ---
    % The friction contribution to the residual at edge ein is:
    %   F_friction_u = -c61 * emean(ein,:) * (tau_b_over_Ub .* N_nodes) * u(ein)
    % We need d(slide_fun)/dC
    ns = gg.ns2;
    nmeanx = gg.nmeanx2; nmeany = gg.nmeany2;
    emean_loc = gg.emean2; fmean_loc = gg.fmean2;
    
    U = sqrt((nmeanx(:,gg.es2)*u(gg.es2)).^2 + (nmeany(:,gg.fs2)*v(gg.fs2)).^2);
    
    dslide_dC = compute_dslide_dC(U, N, C, pp.mu, pp);  % [nIJ x 1]
    
    % dF_u/dC_j = -c61 * emean(ein, j) * dslide_dC(j) * u(ein)  (sum over j via emean)
    % dJ/dC_j = lambda_u^T * dF_u/dC_j + lambda_v^T * dF_v/dC_j
    
    % For u-equation contribution
    dFu_dC = -pp.c61 * spdiags(u(ein), 0, length(ein), length(ein)) * ...
             emean_loc(ein,:) * spdiags(dslide_dC, 0, gg.nIJ, gg.nIJ);
    % For v-equation contribution  
    dFv_dC = -pp.c61 * spdiags(v(fin), 0, length(fin), length(fin)) * ...
             fmean_loc(fin,:) * spdiags(dslide_dC, 0, gg.nIJ, gg.nIJ);
    
    dJdC = (lambda_u' * dFu_dC + lambda_v' * dFv_dC)';  % [nIJ x 1]
    
    % --- 7. Chain rule: dJ/dc = C .* dJ/dC + reg gradient ---
    g = C .* dJdC;
    g = g + opts_inv.alpha * (L'*L) * c;
    g = g + opts_inv.gamma * (c - c_prior);
end

%% ============================================================
%  Assemble velocity Jacobian (extracted from nevis_velocity logic)
%  ============================================================
function [A, rhs] = assemble_velocity_jacobian(u, v, N, aa, pp, gg, oo)
    % Assemble the linear system matrix A and rhs for the SSA velocity solve
    % at the current viscosity (i.e. the last Picard iteration).
    % This is the same system solved in nevis_velocity.m
    
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
end

%% ============================================================
%  d(slide_fun)/dC  (analytical derivative of sliding law w.r.t. C)
%  ============================================================
function dslide_dC = compute_dslide_dC(U, N, C, mu, pp)
    % tau_b/Ub = N * Ub^(1/n-1) * (mu^(-n)*Ub + C^(-n)*N^n)^(-1/n) + C2*Ub^(1/n-1)
    % d(tau_b/Ub)/dC = N * Ub^(1/n-1) * C^(-n-1) * N^n * (mu^(-n)*Ub + C^(-n)*N^n)^(-1/n-1)
    
    n = pp.n_slide;
    Ub = max(U, pp.Ub_reg);
    Np = max(N, pp.N_slide_reg);
    
    bracket = mu.^(-n) .* Ub + C.^(-n) .* Np.^n;
    
    dslide_dC = Np .* Ub.^(1/n - 1) .* C.^(-n-1) .* Np.^n .* bracket.^(-1/n - 1);
    % Note: C2 term has no C dependence, so no contribution
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
    n = pp.n_slide;
    Ubr = max(Ub, pp.Ub_reg);
    Nr = max(N, pp.N_slide_reg);
    tau_b_over_Ub = Nr .* Ubr.^(1/n-1) .* ...
        (mu.^(-n).*Ubr + C.^(-n).*Nr.^n).^(-1/n) + ...
        pp.C2 * Ubr.^(1/n-1);
end
