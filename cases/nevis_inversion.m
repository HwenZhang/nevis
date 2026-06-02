% Template regional NEVIS inversion configuration.
%
% Run with:
%   [inv, vv_hydro, summary] = nevis_run_inversion('./cases/nevis_inversion.m');
%
% The inversion reads reusable regional observations from cfg.dataset.root and
% writes derived inversion products back into the same dataset package.

cfg = struct;
cfg.name = 'weertman_ns3_default';
cfg.root = './';
cfg.code = './src';
cfg.results = 'results';

%% Sliding law
% Available choices:
%   'weertman'              fixed-N inversion with the Weertman-style C law
%   'regularized_coulomb'   iterative C-N inversion using regularized Coulomb sliding
% For parameter sweeps, keep cfg.name and cfg.output.* unique for each run.
cfg.sliding_law = 'weertman';

%% Dataset binding
cfg.dataset.name = 'nevis_2022_140km';
cfg.dataset.root = './data/datasets/nevis_2022_140km';
cfg.dataset.manifest = 'dataset_manifest.m';

%% Source model state used to initialize grid, parameters, and geometry
% This file must contain pp, pd, ps, gg, aa, and oo.
cfg.source.casename = 'n2d_region_ice_inversion_ns3';
cfg.source.state_file = './data/redundant/n2d_region_ice_inversion_ns3/n2d_region_ice_inversion_ns3.mat';

%% Observed velocity from the dataset package
cfg.velocity.file = 'velocity.mat';
cfg.velocity.variable = 'velocity';
cfg.velocity.u_field = 'u_obs_dim';
cfg.velocity.v_field = 'v_obs_dim';

%% Initial effective pressure source
% Available choices:
%   'k_factor'   no external file; phi = phi_a + k_factor*(phi_0 - phi_a)
%   'constant_N' no external file; uniform dimensionless N
%   'file'       load variable from MAT file; variable must contain N or phi
cfg.initial_hydrology.mode = 'k_factor';
cfg.initial_hydrology.k_factor = 0.9;
% cfg.initial_hydrology.mode = 'constant_N';
% cfg.initial_hydrology.N = 1;
% cfg.initial_hydrology.mode = 'file';
% cfg.initial_hydrology.file = './data/redundant/n2d_region_ice_inversion_ns3/0365.mat';
% cfg.initial_hydrology.variable = 'vv';

%% Output products, stored inside cfg.dataset.root
cfg.output.inversion_file = 'inversions/weertman_ns3_default/C_inversion_results.mat';
cfg.output.initial_hydrology_file = 'inversions/weertman_ns3_default/velocity_inverted.mat';
cfg.output.save_intermediate = true;
cfg.output.intermediate_prefix = 'inversions/weertman_ns3_default/C_inversion_outer';
cfg.output.intermediate_hydrology_file = 'inversions/weertman_ns3_default/velocity_inverted_outer.mat';

%% Inversion and forward initialization options
% nevis_inv_C does not provide hidden defaults for these fields. Every solver,
% regularization, prior, optimizer, plotting, and forward-hydrology option
% used by the inversion should be declared here.
cfg.partition_ratio = 0.5;

cfg.solver.iter_max = 100;
cfg.solver.tol_vel = 1e-6;
cfg.solver.display_norms = 0;
cfg.solver.verb = 0;
cfg.solver.boundary_method = 'stress_l_vel_tbl';
cfg.solver.mask_boundary_method = 'stress_free';

cfg.regularization.eps_reg = 1e-1;
cfg.regularization.Ub_reg = 1e-16;
cfg.regularization.N_slide_reg = 1e-16;
cfg.regularization.taud_reg = 1e-16;
cfg.regularization.C2 = 0;
cfg.regularization.adaptive_beta = 0.1;

cfg.prior.mode = 'uniform';
cfg.prior.C = 1;

cfg.opts_inv.u0_reg = 1e-1;
cfg.opts_inv.alpha_schedule = [1e-3, 1e-4, 1e-9, 1e-12];
cfg.opts_inv.gamma_schedule = [1e-6, 1e-7, 1e-8, 1e-9];
cfg.opts_inv.max_iter_schedule = [60, 50, 30, 30];
cfg.opts_inv.max_iter_stage = 50;
cfg.opts_inv.max_iter_total = 200;
cfg.opts_inv.J_tol = 1e-3;
cfg.opts_inv.dJ_tol = 1e-6;
cfg.opts_inv.verbose = false;
cfg.opts_inv.check_grad = false;
cfg.opts_inv.speed_misfit = true;

% Used by cfg.sliding_law = 'regularized_coulomb'.
cfg.opts_inv.max_outer_iter = 3;
cfg.opts_inv.C_tol = 1e-3;
cfg.opts_inv.N_tol = 1e-3;

% false saves a reusable vv_hydro from the inverted velocity and initial N.
% true additionally runs nevis_run_fwd_hydrology, which is slower. For
% regularized_coulomb this controls whether each outer loop updates N.
cfg.forward_hydrology.enabled = false;
cfg.forward_hydrology.mu_s = 0.25;
cfg.forward_hydrology.racmo_runoff = false;
cfg.forward_hydrology.distributed_input = true;
cfg.forward_hydrology.output_casename = 'weertman_ns3_default_forward_hydrology';

cfg.plot.enabled = true;
