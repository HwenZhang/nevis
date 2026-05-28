% Template idealized NEVIS case configuration.
%
% The idealised case has no external topography, velocity, inversion, or runoff data.
% The geometry, moulins, lakes, and forcing are generated analytically.
% 
% Copy this file to ./cases/<new_casename>.m, set cfg.casename, and call
% n2d_idealized_template('<new_casename>').
%
% Data isolation rule:
% - The config file is the only case-specific input.
% - No missing analytical option falls back to another case.

cfg = struct;

%% Case identity
cfg.casename = 'nevis_idealized';
cfg.case_root = './cases';

%% Model and solver options: cfg.oo
cfg.oo.evaluate_variables = 1;      % expand derived variables after the run
cfg.oo.input_gaussian = 1;          % use Gaussian point-input spreading
cfg.oo.relaxation_term = 1;         % 0: alpha*hb leakage; 1: kappa*deltap*hb leakage
cfg.oo.initial_condition = 1;       % 1: analytical steady pressure/sheet state
cfg.oo.mean_perms = 1;              % use mean permeability averaging
cfg.oo.include_blister = 1;         % include blister component
cfg.oo.include_ice = 0;             % no coupled ice velocity solve
cfg.oo.display_residual = 0;        % nonlinear residual text output
cfg.oo.adjust_boundaries = 1;       % allow pressure boundary updates
cfg.oo.yperiodic = 1;               % periodic transverse direction
cfg.oo.xperiodic = 0;               % non-periodic downslope direction
cfg.oo.keep_all_moulins = 0;        % merge repeated moulin grid nodes
cfg.oo.random_moulins = 0;          % 0 disables random moulin placement
cfg.oo.random_lakes = 0;            % 0 disables random lake placement
cfg.oo.surface_runoff = 0;          % analytical moulin input only
cfg.oo.RACMO_runoff = 0;            % no gridded RACMO runoff file
cfg.oo.distributed_input = 0;       % no distributed runoff forcing
cfg.oo.save_timesteps = 1;          % save timestep state files
cfg.oo.save_pts_all = 1;            % save time series at lake/moulin nodes

%% Dimensional physical parameters: cfg.pd
% Function handles are evaluated after scalar fields have been applied.
cfg.pd.alpha_b = 0;                 % blister relaxation rate for relaxation_term = 0
cfg.pd.kappa_b = 1e-10;             % pressure relaxation coefficient for relaxation_term = 1
cfg.pd.c0 = 1;                      % leakage switch used by nondimensionalization
cfg.pd.mu = 1000.0;                 % water viscosity, Pa s
cfg.pd.Ye = 8.8e9;                  % Young's modulus, Pa
cfg.pd.B = @(pd) pd.Ye*(1e3)^3/(12*(1-0.33^2)); % bending stiffness, Pa m^3
cfg.pd.E_lapse = @(pd) 30/1000/pd.td/10^3;      % lapse-rate melt scale
cfg.pd.hb_reg1 = 5e-3;              % primary hb regularization, m
cfg.pd.hb_reg2 = 1e-3;              % secondary hb regularization, m
cfg.pd.N_reg1 = 1e4;                % effective pressure regularization, Pa
cfg.pd.deltap_reg = 1e4;            % pressure-difference regularization, Pa
cfg.pd.B_reg = 0;                   % bending regularization
cfg.pd.G = 0.01;                    % geothermal heat flux, W/m^2
cfg.pd.melt = @(pd) pd.G/pd.rho_w/pd.L;         % geothermal melt, m/s

%% Nondimensional overrides: cfg.pp
% Empty by default; keep this struct so parameter overrides have the same
% entry point as regional cases.
cfg.pp = struct;

%% Analytical geometry: cfg.geometry
% mode choices:
%   'hewitt_2013_sqrt' uses flat bed and s = H0*(1 - x/L)^surface_exponent.
%   'function' uses cfg.geometry.bed_function and surface_function.
%
% For mode = 'function', the functions return dimensional meters and are
% called as:
%   bed_m = bed_function(X_m, Y_m, L, W, pd, ps, pp, gg, oo, cfg)
%   surf_m = surface_function(X_m, Y_m, L, W, pd, ps, pp, gg, oo, cfg)
% where X_m and Y_m are grid-sized coordinate arrays in meters. The return
% value can be scalar, grid-sized, or one value per node.
cfg.geometry.mode = 'hewitt_2013_sqrt';
cfg.geometry.length_m = 5e4;             % domain length, m
cfg.geometry.width_ratio = 0.4;          % width = width_ratio * length
cfg.geometry.nx = 201;                   % number of x nodes
cfg.geometry.ny = 80;                    % number of y nodes
cfg.geometry.bed_elevation_m = 0;        % flat bed elevation, m
cfg.geometry.surface_head_m = 1060;      % upstream surface elevation scale, m
cfg.geometry.surface_exponent = 0.5;     % square-root profile exponent
cfg.geometry.min_thickness_m = 0;        % mask nodes at or below this thickness

% Example custom analytical geometry:
% cfg.geometry.mode = 'function';
% cfg.geometry.bed_function = @(X,Y,L,W,pd,ps,pp,gg,oo,cfg) ...
%     -50 + 10*sin(2*pi*X/L).*cos(2*pi*Y/W);
% cfg.geometry.surface_function = @(X,Y,L,W,pd,ps,pp,gg,oo,cfg) ...
%     900 - 500*X/L + 20*cos(2*pi*Y/W);

%% Initial hydrology state: cfg.initial_hydrology
% mode choices:
%   'k_factor_cavity_equilibrium' prescribes phi = phi_a + k_f*(phi_0-phi_a)
%   and initializes sheet thickness from the cavity equilibrium expression.
%   'function' calls cfg.initial_hydrology.function to update vv and aa.
%
% For mode = 'function', the function must return [vv, aa] and is called as:
%   [vv, aa] = function(aa, vv, pd, ps, pp, gg, oo, cfg)
% This lets a case prescribe arbitrary initial phi, hs, hb, etc. The values
% should use NEVIS nondimensional units, matching aa/vv fields.
cfg.initial_hydrology.mode = 'k_factor_cavity_equilibrium';
cfg.initial_hydrology.k_factor = 1.00;   % fraction of overburden pressure

% Example custom initial condition:
% cfg.initial_hydrology.mode = 'function';
% cfg.initial_hydrology.function = @(aa,vv,pd,ps,pp,gg,oo,cfg) ...
%     deal(setfield(vv, 'phi', aa.phi_a + 0.8*(aa.phi_0-aa.phi_a)), aa);

%% Moulins: cfg.moulins
% mode choices:
%   'coordinates_fractional' places moulins at fractions of L and W.
%   'coordinates_m' places moulins from cfg.moulins.x_m/y_m in meters.
%   'function' calls cfg.moulins.position_function and expects x_m/y_m
%   coordinate arrays in meters.
cfg.moulins.mode = 'coordinates_fractional';
cfg.moulins.x_fraction = [0.25];
cfg.moulins.y_fraction = [0.50];

% Example positions in meters:
% cfg.moulins.mode = 'coordinates_m';
% cfg.moulins.x_m = [12e3 25e3 35e3];
% cfg.moulins.y_m = [5e3 10e3 12e3];
%
% Example analytical position function:
% cfg.moulins.mode = 'function';
% cfg.moulins.position_function = @(L,W,pd,ps,pp,gg,oo,cfg) ...
%     deal([0.2; 0.5; 0.8]*L, [0.5; 0.5; 0.5]*W);

% Prescribed direct moulin input, in m^3/s. Use a scalar for all moulins or
% one value per moulin after repeated/merged moulin nodes are removed.
cfg.moulins.input_m3_per_s = 0;

% Ramp timescale for cfg.moulins.input_m3_per_s.
cfg.moulins.input_ramp_days = 30;

% For fully custom input, set cfg.moulins.input_function. It returns m^3/s
% and is called as:
%   q_m3s = input_function(t, pd, ps, pp, gg, oo, cfg)
% where t is nondimensional model time. Return a scalar or one value per
% active moulin.
% cfg.moulins.input_function = @(t,pd,ps,pp,gg,oo,cfg) ...
%     [0; 5; 10] .* (1-exp(-t/(20*pd.td/ps.t)));

%% Supraglacial lakes: cfg.lakes
% mode choices:
%   'coordinates_fractional' places lakes at fractions of L and W.
cfg.lakes.mode = 'coordinates_fractional';
cfg.lakes.x_fraction = 0.50;
cfg.lakes.y_fraction = 0.50;
cfg.lakes.volume_m3 = 0e7;               % lake drainage volume, m^3
cfg.lakes.drainage_day = 0.5*365;        % center drainage day
cfg.lakes.duration_days = 0.025;         % drainage duration, days

%% Analytical input forcing: cfg.input
% runoff_max_mm_per_day:
%   Maximum prescribed moulin-equivalent melt rate for pp.meltE.
cfg.input.runoff_max_mm_per_day = 30;

% moulin_input_m3_per_s:
%   Legacy fallback for direct moulin input if cfg.moulins.input_m3_per_s is
%   not set. Prefer cfg.moulins.input_m3_per_s for new cases.
cfg.input.moulin_input_m3_per_s = 0;

% ramp_days:
%   E-folding time for runoff and moulin input ramp-up.
cfg.input.ramp_days = 30;

%% Run control: cfg.run
cfg.run.dt_days = 1/24;                   % initial/suggested timestep, days
cfg.run.t_span_days = 1:1:2*365;          % output times, days
