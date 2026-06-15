% End-to-end NEVIS user workflow:
% dataset package -> inversion with automatic hydrology-only spinup -> forward model -> drainage case.
%
% Edit only the USER SETTINGS section for a normal run.
%
% WHERE TO TUNE PARAMETERS
% ------------------------
% Each case owns a subdirectory under ./cases/<case_name>/ holding four
% fully self-contained config files:
%
%   cases/<case_name>/<case_name>_inversion.m         (inversion config)
%   cases/<case_name>/<case_name>_inversion_spinup.m  (hydrology spinup case)
%   cases/<case_name>/<case_name>_forward.m           (forward regional case)
%   cases/<case_name>/<case_name>_drainage.m          (drainage case)
%
% On the first run for a new case_name the workflow materializes these
% files by inlining the templates and appending the workflow-specific
% bindings (dataset/output/source paths). They are SELF-CONTAINED: each
% file has every cfg.* assignment it needs, with no `run(template)` line
% and no runtime fallback to the template. From that point on the
% workflow leaves them alone, so:
%
%   - To change physics/numerics for THIS case: edit the per-case file.
%   - To start fresh from the current template: set reinit_case = true
%     below (overwrites the per-case files) or delete cases/<case_name>/.
%   - The original templates (cases/templates/nevis_inversion.m,
%     cases/templates/nevis_regional.m,
%     cases/templates/nevis_regional_inversion_spinup.m) are only used at
%     first-materialization or when reinit_case = true.
%
% Field index (per file):
%   <case>_inversion.m:        cfg.solver.* cfg.regularization.* cfg.prior.*
%                              cfg.opts_inv.* cfg.partition_ratio
%                              cfg.forward_hydrology.* cfg.plot.*
%   <case>_inversion_spinup.m: cfg.oo.* cfg.pd.* cfg.pp.* cfg.run.t_span_days
%   <case>_forward.m:          cfg.oo.* cfg.pd.* cfg.pp.* cfg.run.*
%   <case>_drainage.m:         cfg.oo.* cfg.pd.* cfg.pp.* cfg.run.*
%                              cfg.initial_hydrology.result_case
%
% See also cases/README.md for the case/dataset contract, and
% inversion/INVERSION.md for the inversion contract.

addpath('src');
addpath('src/workflow');
addpath('inversion');

%% USER SETTINGS
case_name = 'nevis_2022_140km_botright30';   % dataset/config prefix for this run
bbox_km = [-80.10, 90.15, -80.07, 30.18];                             
                                          % [] = default [-80.10, 60.15, -50.07, 30.18] km: [xmin, xmax, ymin, ymax]

run_dataset_build = true;                 % rebuild geometry/velocity/runoff/stations
run_inversion = true;                     % run inversion; auto-runs spinup if needed
run_forward = true;                       % run forward regional model (lakes off; spinup-like)
run_drainage = true;                      % run drainage variant (lakes on); needs forward's 0365.mat

reinit_case = false;                      % true: regenerate cases/<case_name>/ from templates
                                          %        (overwrites local edits in that dir)

% First-run-only seeds (baked into the per-case files; ignored on later runs
% unless reinit_case = true). To change them for an existing case, edit the
% per-case file directly.
sliding_law = 'weertman';                 % 'weertman' or 'regularized_coulomb'
spinup_days = 1:365;                      % hydrology-only spinup output days
spinup_timestep_file = '0365.mat';        % spinup timestep used to initialize inversion

%% PREPARE WORKFLOW
opts = nevis_demo_options(case_name, bbox_km);
opts.inversion.sliding_law = sliding_law;
opts.spinup.t_span_days = spinup_days;
opts.spinup.timestep_file = spinup_timestep_file;
opts.reinit_case = reinit_case;

workflow = nevis_prepare_demo_workflow(opts);

%% BUILD DATASET
if run_dataset_build
    dataset_report = nevis_build_demo_dataset(workflow);
end

%% RUN INVERSION
% Reads cases/<case_name>/<case_name>_inversion.m. That file inlines the
% inversion template at first materialization; edit it to tune solver,
% regularization, prior, opts_inv, partition_ratio, forward_hydrology, plot.
% The hydrology-only spinup is driven by
% cases/<case_name>/<case_name>_inversion_spinup.m (cfg.oo/cfg.pd/cfg.pp).
if run_inversion
    [inv, vv_hydro, summary] = nevis_run_inversion( ...
        workflow.inversion.config_file);
end

%% RUN FORWARD MODEL
% Reads cases/<case_name>/<case_name>_forward.m. Edit it to tune cfg.oo,
% cfg.pd, cfg.pp, cfg.run.*. Lakes are off, so this is a spinup-like year.
if run_forward
    [tt, vv, info, state] = nevis_run_regional(workflow.forward.case);
end

%% RUN DRAINAGE CASE
% Reads cases/<case_name>/<case_name>_drainage.m. This workflow-managed
% case starts from the forward run's day-365 timestep using
% cfg.initial_hydrology.mode='result_timestep' and turns the case-local
% 2022 lake drainage catalogue on.
if run_drainage
    [tt_d, vv_d, info_d, state_d] = nevis_run_regional( ...
        workflow.drainage.case);
end
