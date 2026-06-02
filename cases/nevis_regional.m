% Template regional NEVIS case configuration.
%
% Copy this file to ./cases/<new_casename>.m, set cfg.casename below, then
% bind the case to a shared cfg.dataset package. The standardized runner
% n2d_regional_template(casename) reads only this config file for
% case-specific parameters and reads heavy data through the dataset manifest.
%
% Data isolation rule:
% - cfg.case_root contains case-specific config files.
% - Every relative data file path is resolved under cfg.dataset.root.
% - Missing files, missing variables, or mismatched array sizes are errors.
% - The workflow never falls back to another case or to global default data.

cfg = struct;

%% Case identity
% cfg.casename:
%   Must exactly match the folder name passed to n2d_regional_template.
%   Example call: n2d_regional_template('nevis_regional')
cfg.casename = 'nevis_regional';

% cfg.case_root:
%   Directory containing case config files.
cfg.case_root = './cases';

%% Shared dataset binding
% cfg.dataset:
%   Heavy, reusable input data live in a shared immutable dataset package.
%   Parameter-sweep cases should point to the same dataset.root and change
%   only cfg.oo/cfg.pd/cfg.pp/cfg.run values.
cfg.dataset.name = 'nevis_2022_140km';
cfg.dataset.root = './data/datasets/nevis_2022_140km';
cfg.dataset.manifest = 'dataset_manifest.m';

%% Model and solver options: cfg.oo
% These fields are applied after nevis_defaults and before the case setup.
% Use numeric 0/1 values for legacy logical options unless otherwise noted.

% evaluate_variables:
%   1 expands derived variables after the final timestep using
%   nevis_backbone; 0 skips expansion.
cfg.oo.evaluate_variables = 1;

% include_ice:
%   1 couples hydrology to the ice velocity solve; 0 disables ice coupling
%   where supported by the core model.
cfg.oo.include_ice = 1;

% adjust_boundaries:
%   1 lets nevis_timesteps update pressure boundary nodes; 0 keeps the
%   initial boundary set fixed.
cfg.oo.adjust_boundaries = 1;

% input_gaussian:
%   1 spreads point input with the model's Gaussian input treatment; 0 uses
%   the legacy unsmoothed point treatment.
cfg.oo.input_gaussian = 1;

% relaxation_term:
%   0 uses an alpha*hb relaxation form; 1 uses an alpha*deltap*hb form.
cfg.oo.relaxation_term = 1;

% initial_condition:
%   0 uses a legacy saved initial condition when available in old workflows;
%   1 initializes a steady-state drainage system from configured data.
cfg.oo.initial_condition = 1;

% mean_perms:
%   1 uses mean permeability averaging; 0 disables that averaging.
cfg.oo.mean_perms = 1;

% modified_mean_perms:
%   1 uses the modified mean permeability rule; 0 uses the standard rule.
cfg.oo.modified_mean_perms = 0;

% display_residual:
%   1 prints nonlinear residual diagnostics; 0 keeps output quieter.
cfg.oo.display_residual = 0;

% visualize_vel:
%   1 draws velocity diagnostics during timestepping; 0 disables plotting.
cfg.oo.visualize_vel = 0;

% N_coupling:
%   1 couples effective pressure into sliding/ice dynamics; 0 removes that
%   coupling where supported.
cfg.oo.N_coupling = 1;

% U_coupling:
%   1 couples velocity into hydrology opening/melt terms; 0 removes that
%   coupling where supported.
cfg.oo.U_coupling = 1;

% boundary_method:
%   Ice-flow boundary treatment used by nevis_label_ice_test.
%   Available choices in the current code include:
%     'default'
%     'per_tblr'
%     'vel_tblr'
%     'stress_tblr'
%     'vel_r_stress_tbl'
%     'stress_l_vel_tbl'
cfg.oo.boundary_method = 'stress_l_vel_tbl';

% mask_boundary_method:
%   How to treat irregular mask boundaries for the ice solve.
%   Available choices in the current code include:
%     'velocity'
%     'stress_free'
cfg.oo.mask_boundary_method = 'stress_free';

% plot_residual:
%   1 plots nonlinear residual diagnostics; 0 disables those plots.
cfg.oo.plot_residual = 0;

% max_iter_new:
%   Maximum Newton iterations per nonlinear solve.
cfg.oo.max_iter_new = 50;

% partition_ratio:
%   Fraction of the inverted sliding coefficient assigned to C1 in
%   nevis_inv_partition. 1.0 means all inverted C goes into C1 and C2 is zero.
cfg.oo.partition_ratio = 1.0;

% iterative_solver:
%   1 uses the configured iterative linear solver path where supported;
%   0 uses the direct/default linear solver path.
cfg.oo.iterative_solver = 0;

%% Dimensional physical parameters: cfg.pd
% Units are dimensional SI unless stated otherwise. Function handles are
% evaluated after scalar cfg.pd fields have been applied, so derived
% parameters can depend on fields already present in pd.

% alpha_b:
%   Blister relaxation rate coefficient. Unit depends on relaxation_term;
%   0 disables this relaxation contribution.
cfg.pd.alpha_b = 0;

% kappa_b:
%   Blister/hydraulic relaxation coefficient used in the hydrology model.
cfg.pd.kappa_b = 5e-11;

% mu:
%   Effective viscosity or frictional resistance parameter used by the
%   coupled model.
cfg.pd.mu = 5.0;

% c_e_reg2:
%   Secondary regularization floor for the englacial/channel exchange term.
cfg.pd.c_e_reg2 = 0.00/1e3/9.81;

% u_b:
%   Characteristic basal sliding speed in m/s. This template sets 100 m/yr.
cfg.pd.u_b = @(pd) 100/pd.ty;

% sigma:
%   Conductivity/permeability exponent regularization parameter.
cfg.pd.sigma = 1e-3;

% h_r:
%   Sheet thickness scale used by regularization terms, in meters.
cfg.pd.h_r = 0.1;

% l_r:
%   Regularization length scale, in meters.
cfg.pd.l_r = 10;

% l_c:
%   Channel/cavity length scale, in meters.
cfg.pd.l_c = 1000;

% k_s:
%   Sheet conductivity or permeability coefficient.
cfg.pd.k_s = 0.1;

% tau_b:
%   Characteristic basal shear stress, in Pa.
cfg.pd.tau_b = 60e3;

% melt:
%   Basal melt rate, in m/s. This template combines geothermal and frictional
%   heating using pd.G, pd.u_b, pd.tau_b, pd.rho_w, and pd.L.
cfg.pd.melt = @(pd) (pd.G + pd.u_b*pd.tau_b)/pd.rho_w/pd.L;

% meltinterior:
%   Interior melt contribution. 0 disables the extra contribution.
cfg.pd.meltinterior = 0;

% hb_reg1:
%   Primary blister/cavity opening regularization thickness, in meters.
cfg.pd.hb_reg1 = 5e-3;

% hb_reg2:
%   Secondary blister/cavity opening regularization thickness, in meters.
cfg.pd.hb_reg2 = 1e-3;

% N_reg1:
%   Effective pressure regularization scale, in Pa.
cfg.pd.N_reg1 = 1e4;

% deltap_reg:
%   Pressure-difference regularization scale, in Pa.
cfg.pd.deltap_reg = 1e4;

% B_reg:
%   Elastic bending regularization/stiffness scale. This template uses a
%   plate stiffness expression with a 1000 m thickness and Poisson ratio 0.33.
cfg.pd.B_reg = @(pd) pd.Ye*(1000)^3/(12*(1-0.33)^2);

%% Nondimensional overrides: cfg.pp
% These are applied after nevis_nondimension and after ice parameter updates.

% eps_reg:
%   Nondimensional ice-flow regularization/viscosity floor.
cfg.pp.eps_reg = 0.01;

% Ub_reg:
%   Nondimensional sliding speed regularization floor.
cfg.pp.Ub_reg = 1e-16;

% N_slide_reg:
%   Nondimensional effective-pressure sliding regularization floor.
cfg.pp.N_slide_reg = 1e-16;

% taud_reg:
%   Nondimensional driving-stress regularization floor.
cfg.pp.taud_reg = 1e-16;

% C2:
%   Nondimensional pressure-insensitive sliding coefficient. Usually set by
%   partitioning; 0 disables the C2 contribution.
cfg.pp.C2 = 0;

%% Geometry data: cfg.geometry
% mode:
%   Geometry is currently read from a MATLAB struct variable in a .mat file.
% Required fields are configured below. Relative file paths are resolved under
% cfg.dataset.root and checked against cfg.dataset.manifest.
cfg.geometry.file = 'geometry.mat';
cfg.geometry.variable = 'geometry';
cfg.geometry.x_field = 'X_m';
cfg.geometry.y_field = 'Y_m';
cfg.geometry.bed_field = 'B_m';
cfg.geometry.surface_field = 'S_m';

% skip:
%   Downsampling stride used when building the model grid. Larger values
%   produce a coarser and cheaper model.
cfg.geometry.skip = 6;

% min_thickness_m:
%   Nodes with ice thickness <= this value are masked out.
cfg.geometry.min_thickness_m = 0;

%% Velocity observations: cfg.velocity
% mode choices:
%   'file' reads observed velocity from cfg.velocity.file.
% units choices:
%   'm_per_year' is the currently supported unit in the importer.
cfg.velocity.file = 'velocity.mat';
cfg.velocity.variable = 'velocity';
cfg.velocity.u_field = 'u_obs_dim';
cfg.velocity.v_field = 'v_obs_dim';
cfg.velocity.units = 'm_per_year';
cfg.velocity.mode = 'file';

%% Inversion data: cfg.inversion
% mode choices:
%   'file' reads a prior inversion result containing fields required by
%   nevis_inv_partition, including C_hat and observed velocities.
%   'uniform' uses spatially uniform C=1 and C2=0 without reading a file.
% file:
%   Dataset-relative path. For alternate inversion packages, use paths like
%   'inversions/<name>/C_inversion_results.mat'.
cfg.inversion.file = 'C_inversion_results.mat';
cfg.inversion.mode = 'file';

%% Initial hydrology state: cfg.initial_hydrology
% mode choices:
%   'file' reads cfg.initial_hydrology.variable from the configured .mat file.
% Required variable fields:
%   N is required; u and v are optional velocity fields.
% file:
%   Dataset-relative path. For alternate inversion packages, use paths like
%   'inversions/<name>/velocity_inverted.mat'.
cfg.initial_hydrology.file = 'velocity_inverted.mat';
cfg.initial_hydrology.variable = 'vv_hydro';
cfg.initial_hydrology.mode = 'file';

%% Moulin data: cfg.moulins
% mode choices:
%   'node_indices' reads precomputed pp.ni_m and pp.sum_m from a saved struct.
%   'coordinates' maps raw moulin x/y coordinates to the current model grid.
cfg.moulins.file = 'nevis_170207a.mat';
cfg.moulins.variable = 'pp';
cfg.moulins.mode = 'node_indices';

% ni_field, sum_field:
%   Field names used when mode = 'node_indices'.
cfg.moulins.ni_field = 'ni_m';
cfg.moulins.sum_field = 'sum_m';

% x_field, y_field:
%   Uncomment and set these when mode = 'coordinates'.
% cfg.moulins.x_field = 'x_m';
% cfg.moulins.y_field = 'y_m';

% keep_all:
%   false lets the importer drop or move moulins outside the active domain;
%   true requests that all configured moulins be retained where possible.
cfg.moulins.keep_all = false;

% move_to_domain:
%   true moves out-of-domain moulins to nearest valid locations where the
%   legacy moulin mapper supports it.
cfg.moulins.move_to_domain = true;

%% Lake data: cfg.lakes
% mode choices:
%   'environs_lakes_catalogue' reads the current environs_lakes catalogue
%   format used by this Nevis template.
%   'table' reads a generic table/struct with configured coordinate, timing,
%   and volume fields.
cfg.lakes.file = 'mechanistic_drainage_catalogue_2022-main/catalogues/environs_lakes_2022B_250505_archive.mat';
cfg.lakes.variable = 'environs_lakes';
cfg.lakes.mode = 'environs_lakes_catalogue';

% spinup_volume_scale:
%   Multiplier applied to lake volumes for spinup. 0 removes lake input during
%   spinup; 1 uses catalogue volumes.
cfg.lakes.spinup_volume_scale = 0;

% time_offset_days:
%   Day offset added to all configured lake drainage times.
cfg.lakes.time_offset_days = 0;

% For mode = 'table', configure fields such as:
% cfg.lakes.x_field = 'x_m';
% cfg.lakes.y_field = 'y_m';
% cfg.lakes.volume_field = 'volume_m3';
% cfg.lakes.drainage_start_field = 'drainage_time_start';
% cfg.lakes.drainage_end_field = 'drainage_time_end';
% cfg.lakes.duration_field = 'drainage_duration_days';

%% Runoff forcing: cfg.runoff
% mode choices:
%   'file' reads a time-by-node runoff matrix from cfg.runoff.file.
% units choices:
%   'mm_we_per_day' is the expected unit for the current importer.
cfg.runoff.file = 'racmo_runoff_2022.mat';
cfg.runoff.variable = 'racmo_runoff_2022';
cfg.runoff.units = 'mm_we_per_day';
cfg.runoff.mode = 'file';

% smooth_days:
%   Moving-average window in days used to construct region.runoff_smoothed.
%   Set to 0 or 1 to disable smoothing.
cfg.runoff.smooth_days = 5;

% use_smoothed:
%   true uses the smoothed runoff for the simulation; false stores smoothed
%   runoff for diagnostics but uses the raw configured runoff.
cfg.runoff.use_smoothed = false;

% subtract_lake_volume:
%   true scales runoff down by total lake volume / total runoff volume to
%   avoid double-counting lake water; false leaves runoff unchanged.
cfg.runoff.subtract_lake_volume = true;

% distributed_input:
%   true applies runoff as distributed surface input; false uses only point
%   input pathways where configured.
cfg.runoff.distributed_input = true;

% surface_runoff:
%   true enables surface runoff forcing in nevis_inputs.
cfg.runoff.surface_runoff = true;

% RACMO_runoff:
%   Legacy switch for RACMO-specific runoff pathways. This template keeps it
%   false because runoff is already preprocessed into a node matrix.
cfg.runoff.RACMO_runoff = false;

% meltE_mm_per_day:
%   Uniform extra melt forcing in mm water equivalent per day.
cfg.runoff.meltE_mm_per_day = 0;

% meltE_ramp_days:
%   E-folding ramp timescale for the uniform melt forcing.
cfg.runoff.meltE_ramp_days = 20;

%% Station/GPS output points: cfg.stations
% mode choices:
%   'file' reads station coordinates from cfg.stations.file.
cfg.stations.file = 'station_timeseries_2022.mat';
cfg.stations.variable = 'station_data';
cfg.stations.mode = 'file';
cfg.stations.x_field = 'x_m';
cfg.stations.y_field = 'y_m';

%% Run control: cfg.run
% dt_days:
%   Initial/suggested timestep in days. The adaptive timestepper may change
%   the actual internal step within its configured limits.
cfg.run.dt_days = 1/24;

% t_span_days:
%   Output/save times in day units.
cfg.run.t_span_days = 1:365;

% save_timesteps:
%   true writes timestep state files under ./results/<casename>/.
cfg.run.save_timesteps = true;

% save_pts_all:
%   true saves time series at lakes, moulins, and configured stations.
cfg.run.save_pts_all = true;
