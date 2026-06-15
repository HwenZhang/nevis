function [inv, vv_hydro, summary] = nevis_run_inversion(config_file)
%NEVIS_RUN_INVERSION Dispatch an inversion config to the selected sliding law.
%
% Usage:
%   [inv, vv_hydro, summary] = nevis_run_inversion('./cases/templates/nevis_inversion.m');
%
% The config must define cfg.sliding_law. If cfg.initial_hydrology.mode is
% "file", this entry point also ensures the configured source state and
% hydrology file exist. Missing files are produced from the configured
% hydrology-only spinup case; no default or legacy case is used.

if nargin < 1 || isempty(config_file)
    error('nevis_run_inversion:MissingConfig', ...
        'Call nevis_run_inversion with an inversion config file.');
end

cfg = load_dispatch_config(config_file);
if ~isfield(cfg, 'root') || isempty(cfg.root)
    cfg.root = './';
end
if ~isfield(cfg, 'results') || isempty(cfg.results)
    cfg.results = 'results';
end
if ~isfield(cfg, 'sliding_law') || isempty(cfg.sliding_law)
    error('nevis_run_inversion:MissingSlidingLaw', ...
        'Inversion config must define cfg.sliding_law.');
end

if isfield(cfg, 'code') && ~isempty(cfg.code)
    addpath(cfg.code);
end
addpath(fullfile(cfg.root, 'inversion'));

cfg = ensure_file_mode_inputs(cfg);

sliding_law = lower(strrep(cfg.sliding_law, '-', '_'));
switch sliding_law
    case {'weertman', 'regularized_weertman', 'regularised_weertman'}
        [inv, vv_hydro, summary] = nevis_inv_C(cfg);
    case {'regularized_coulomb', 'regularised_coulomb', 'coulomb'}
        [inv, vv_hydro, summary] = nevis_inv_C1(cfg);
    otherwise
        error('nevis_run_inversion:UnsupportedSlidingLaw', ...
            ['Unsupported cfg.sliding_law "%s". Use "weertman" or ' ...
             '"regularized_coulomb".'], cfg.sliding_law);
end
summary.sliding_law = sliding_law;
end

function cfg = load_dispatch_config(config_file)
    if isstruct(config_file)
        cfg = config_file;
        if ~isfield(cfg, 'config_file')
            cfg.config_file = '<struct>';
        end
        return
    end

    if exist(config_file, 'file') ~= 2
        error('nevis_run_inversion:MissingConfig', ...
            'Inversion config not found: %s', config_file);
    end

    [~, ~, ext] = fileparts(config_file);
    switch lower(ext)
        case '.m'
            run(config_file);
            if exist('cfg', 'var') ~= 1
                error('nevis_run_inversion:MissingCfg', ...
                    'M-file inversion config must assign variable cfg.');
            end
        case '.mat'
            data = load(config_file);
            if ~isfield(data, 'cfg')
                error('nevis_run_inversion:MissingCfg', ...
                    'MAT inversion config must contain variable cfg.');
            end
            cfg = data.cfg;
        case '.json'
            cfg = jsondecode(fileread(config_file));
        otherwise
            error('nevis_run_inversion:UnsupportedConfig', ...
                'Unsupported inversion config extension: %s', ext);
    end
    cfg.config_file = config_file;
end

function cfg = ensure_file_mode_inputs(cfg)
    if ~isfield(cfg, 'root') || isempty(cfg.root)
        cfg.root = './';
    end
    if ~isfield(cfg, 'initial_hydrology') || ...
            ~isfield(cfg.initial_hydrology, 'mode') || ...
            ~strcmpi(cfg.initial_hydrology.mode, 'file')
        return
    end

    require_fields(cfg, {'source', 'spinup', 'dataset'}, 'cfg');
    require_fields(cfg.source, {'state_file'}, 'cfg.source');
    require_fields(cfg.initial_hydrology, {'file', 'variable'}, ...
        'cfg.initial_hydrology');
    require_fields(cfg.spinup, {'casename', 'timestep_file'}, 'cfg.spinup');

    source_file = resolve_dataset_or_project_path(cfg.source.state_file, cfg);
    hydrology_file = resolve_dataset_or_project_path( ...
        cfg.initial_hydrology.file, cfg);

    missing_source = exist(source_file, 'file') ~= 2;
    missing_hydrology = exist(hydrology_file, 'file') ~= 2;
    if ~(missing_source || missing_hydrology)
        cfg.source.state_file = source_file;
        cfg.initial_hydrology.file = hydrology_file;
        return
    end

    spinup_case = cfg.spinup.casename;
    assert_hydrology_only_spinup_config(cfg, spinup_case);

    if spinup_outputs_missing(cfg, spinup_case)
        fprintf(['Missing inversion file-mode input(s); running ' ...
            'hydrology-only spinup "%s".\n'], spinup_case);
        nevis_run_regional(spinup_case);
    else
        fprintf(['Missing inversion file-mode input(s); using existing ' ...
            'spinup output "%s" for export.\n'], spinup_case);
    end

    nevis_export_spinup_for_inversion(cfg);

    if exist(source_file, 'file') ~= 2
        error('nevis_run_inversion:MissingGeneratedSourceState', ...
            'Spinup export did not create cfg.source.state_file: %s', ...
            source_file);
    end
    if exist(hydrology_file, 'file') ~= 2
        error('nevis_run_inversion:MissingGeneratedInitialHydrology', ...
            ['Spinup export did not create cfg.initial_hydrology.file: ' ...
             '%s'], hydrology_file);
    end

    cfg.source.state_file = source_file;
    cfg.initial_hydrology.file = hydrology_file;
end

function tf = spinup_outputs_missing(cfg, spinup_case)
    spinup_file = resolve_spinup_run_file(cfg, spinup_case);
    timestep_file = resolve_spinup_timestep_file(cfg, spinup_case);
    tf = exist(spinup_file, 'file') ~= 2 || exist(timestep_file, 'file') ~= 2;
end

function assert_hydrology_only_spinup_config(cfg, spinup_case)
    spinup_config = nevis_resolve_case_config(spinup_case, cfg.root);
    spinup_cfg = nevis_load_region_config(spinup_config, cfg.root);
    if ~isfield(spinup_cfg, 'oo') || ...
            ~isfield(spinup_cfg.oo, 'include_ice') || ...
            spinup_cfg.oo.include_ice ~= 0
        error('nevis_run_inversion:SpinupMustBeHydrologyOnly', ...
            ['Configured inversion spinup "%s" must set ' ...
             'cfg.oo.include_ice = 0.'], spinup_case);
    end
    if ~strcmp(spinup_cfg.casename, spinup_case)
        error('nevis_run_inversion:SpinupCaseNameMismatch', ...
            'Spinup config casename "%s" does not match "%s".', ...
            spinup_cfg.casename, spinup_case);
    end
end

function path_out = resolve_spinup_run_file(cfg, spinup_case)
    if isfield(cfg.source, 'result_file') && ~isempty(cfg.source.result_file)
        path_out = resolve_project_path(cfg.source.result_file, cfg);
    else
        path_out = fullfile(cfg.root, cfg.results, spinup_case, ...
            [spinup_case, '.mat']);
    end
end

function path_out = resolve_spinup_timestep_file(cfg, spinup_case)
    path_in = cfg.spinup.timestep_file;
    if is_absolute_path(path_in) || startsWith(path_in, './') || ...
            startsWith(path_in, ['..' filesep])
        path_out = resolve_project_path(path_in, cfg);
    else
        path_out = fullfile(cfg.root, cfg.results, spinup_case, path_in);
    end
end

function path_out = resolve_dataset_or_project_path(path_in, cfg)
    if is_absolute_path(path_in)
        path_out = path_in;
    elseif startsWith(path_in, './') || startsWith(path_in, ['..' filesep])
        path_out = fullfile(cfg.root, path_in);
    else
        path_out = fullfile(cfg.dataset.root, path_in);
    end
end

function path_out = resolve_project_path(path_in, cfg)
    if is_absolute_path(path_in)
        path_out = path_in;
    else
        path_out = fullfile(cfg.root, path_in);
    end
end

function tf = is_absolute_path(path_in)
    tf = startsWith(path_in, filesep) || ...
        ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]', 'once'));
end

function require_fields(s, names, label)
    for i = 1:numel(names)
        if ~isfield(s, names{i}) || isempty(s.(names{i}))
            error('nevis_run_inversion:MissingConfigField', ...
                'Missing inversion config field %s.%s.', label, names{i});
        end
    end
end
