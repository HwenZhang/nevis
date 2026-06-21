% End-to-end driver for the migrated nevis_2022_140km regional case.
%
% This replaces the old top-level direct-run scripts:
%   n2d_regional_analytic_eps1e_02_kappa5e_11_mu5e0_partition5e_01_k01e_01_spinup.m
%   n2d_regional_analytic_V2e8_eps1e_02_kappa5e_11_mu5e0_partition5e_01_k01e_01_drainage_highelev.m
%
% Case-specific parameters now live in:
%   cases/nevis_2022_140km/nevis_2022_140km_inversion_spinup.m
%   cases/nevis_2022_140km/nevis_2022_140km_inversion.m
%   cases/nevis_2022_140km/nevis_2022_140km_forward.m
%   cases/nevis_2022_140km/nevis_2022_140km_drainage.m
%
% Edit only the USER SETTINGS section for normal runs.

addpath('src');
addpath('src/workflow');
addpath('inversion');

%% USER SETTINGS
case_name = 'nevis_2022_140km';

run_dataset_build = true;  % false preserves the existing reference package/manifest
run_inversion = true;       % auto-runs nevis_2022_140km_inversion_spinup if needed
run_forward = true;         % migrated spinup-like annual run, lakes off
run_drainage = true;        % migrated high-elevation manual drainage case

%% RESOLVE CASE FILES
case_dir = fullfile('cases', case_name);
inversion_config = fullfile(case_dir, [case_name, '_inversion.m']);
forward_case = [case_name, '_forward'];
drainage_case = [case_name, '_drainage'];

if exist(inversion_config, 'file') ~= 2
    error('run_nevis_2022_140km:MissingInversionConfig', ...
        'Missing inversion config: %s', inversion_config);
end

%% OPTIONAL DATASET BUILD
if run_dataset_build
    addpath('data/pipelines');
    dataset_root = fullfile('data', 'datasets', case_name);
    raw_root = fullfile('data', 'raw');
    bbox_km = [];
    dataset_report = build_nevis_dataset(dataset_root, raw_root, bbox_km);
end

%% RUN INVERSION
if run_inversion
    [inv, vv_hydro, summary] = nevis_run_inversion(inversion_config);
end

%% RUN FORWARD MODEL
if run_forward
    [tt, vv, info, state] = nevis_run_regional(forward_case);
end

%% RUN DRAINAGE CASE
if run_drainage
    [tt_d, vv_d, info_d, state_d] = nevis_run_regional(drainage_case);
end
