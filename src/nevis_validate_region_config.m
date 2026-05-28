function cfg = nevis_validate_region_config(cfg)
%NEVIS_VALIDATE_REGION_CONFIG Validate the strict per-case data contract.

required_top = {'casename', 'case_root', 'dataset', 'oo', 'pd', 'pp', ...
    'geometry', 'velocity', 'inversion', 'initial_hydrology', 'moulins', ...
    'lakes', 'runoff', 'stations', 'run'};
require_fields(cfg, required_top, 'cfg');

if exist(cfg.case_root, 'dir') ~= 7
    error('nevis_validate_region_config:MissingCaseRoot', ...
        'Case root does not exist: %s', cfg.case_root);
end
require_fields(cfg.dataset, {'name', 'root', 'manifest'}, 'cfg.dataset');
if exist(cfg.dataset.root, 'dir') ~= 7
    error('nevis_validate_region_config:MissingDatasetRoot', ...
        'Dataset root does not exist: %s', cfg.dataset.root);
end

manifest_file = nevis_region_resolve_path(cfg.dataset.manifest, cfg);
dataset = nevis_load_dataset_manifest(manifest_file);
if ~isfield(dataset, 'name') || ~strcmp(dataset.name, cfg.dataset.name)
    error('nevis_validate_region_config:DatasetNameMismatch', ...
        'Config dataset "%s" does not match manifest dataset "%s".', ...
        cfg.dataset.name, get_opt(dataset, 'name', '<missing>'));
end
cfg.dataset.manifest_file = manifest_file;
cfg.dataset.manifest_data = dataset;

require_fields(cfg.geometry, {'file', 'variable', 'x_field', 'y_field', ...
    'bed_field', 'surface_field', 'skip'}, 'cfg.geometry');
require_fields(cfg.velocity, {'file', 'variable', 'u_field', 'v_field', ...
    'units', 'mode'}, 'cfg.velocity');
require_fields(cfg.inversion, {'mode'}, 'cfg.inversion');
if ~strcmpi(cfg.inversion.mode, 'uniform')
    require_fields(cfg.inversion, {'file'}, 'cfg.inversion');
end
require_fields(cfg.initial_hydrology, {'mode'}, 'cfg.initial_hydrology');
require_fields(cfg.moulins, {'file', 'variable', 'mode'}, 'cfg.moulins');
require_fields(cfg.lakes, {'mode'}, 'cfg.lakes');
require_fields(cfg.runoff, {'file', 'variable', 'mode'}, 'cfg.runoff');
require_fields(cfg.stations, {'mode'}, 'cfg.stations');
require_fields(cfg.run, {'dt_days', 'save_timesteps', 'save_pts_all'}, ...
    'cfg.run');
if ~isfield(cfg.run, 't_span_days') && ~isfield(cfg.run, 't_span_after_start_days')
    error('nevis_validate_region_config:MissingField', ...
        'Missing required config field: cfg.run.t_span_days or cfg.run.t_span_after_start_days');
end

require_fields(cfg.oo, {'evaluate_variables', 'include_ice', ...
    'adjust_boundaries', 'input_gaussian', 'relaxation_term', ...
    'initial_condition', 'mean_perms', 'modified_mean_perms', ...
    'display_residual', 'visualize_vel', 'N_coupling', ...
    'U_coupling', 'boundary_method', 'mask_boundary_method', ...
    'plot_residual', 'max_iter_new', 'partition_ratio', ...
    'iterative_solver'}, 'cfg.oo');
if cfg.oo.partition_ratio < 0 || cfg.oo.partition_ratio > 1
    error('nevis_validate_region_config:InvalidPartitionRatio', ...
        'cfg.oo.partition_ratio must be between 0 and 1.');
end
require_fields(cfg.pd, {'alpha_b', 'kappa_b', 'mu', 'c_e_reg2', 'u_b', ...
    'sigma', 'h_r', 'l_r', 'l_c', 'k_s', 'tau_b', 'melt', ...
    'meltinterior', 'hb_reg1', 'hb_reg2', 'N_reg1', ...
    'deltap_reg', 'B_reg'}, 'cfg.pd');
require_fields(cfg.pp, {'eps_reg', 'Ub_reg', 'N_slide_reg', ...
    'taud_reg', 'C2'}, 'cfg.pp');

% Eagerly check all configured data files. Any missing configured file is an
% error, even if a legacy/default path exists elsewhere.
sections = {'geometry', 'velocity', 'inversion', 'initial_hydrology', ...
    'moulins', 'lakes', 'runoff', 'stations'};
for i = 1:numel(sections)
    section = sections{i};
    scfg = cfg.(section);
    mode = get_opt(scfg, 'mode', 'file');
    if any(strcmpi(mode, {'none', 'manual', 'coordinates', ...
            'result_timestep', 'k_factor', 'uniform'}))
        continue
    end
    require_fields(scfg, {'file'}, ['cfg.' section]);
    path = nevis_region_resolve_path(scfg.file, cfg);
    if exist(path, 'file') ~= 2
        error('nevis_validate_region_config:MissingDataFile', ...
            'Configured %s file does not exist for case "%s": %s', ...
            section, cfg.casename, path);
    end
    assert_manifest_file(dataset, section, scfg.file, cfg.casename);
end

switch lower(cfg.initial_hydrology.mode)
    case {'file', 'k_factor'}
        if strcmpi(cfg.initial_hydrology.mode, 'file')
            require_fields(cfg.initial_hydrology, {'file', 'variable'}, ...
                'cfg.initial_hydrology');
        end
    case 'result_timestep'
        require_fields(cfg.initial_hydrology, {'result_case', ...
            'timestep_file'}, 'cfg.initial_hydrology');
    otherwise
        error('nevis_validate_region_config:UnsupportedInitialHydrology', ...
            'Unsupported cfg.initial_hydrology.mode: %s', ...
            cfg.initial_hydrology.mode);
end

switch lower(cfg.lakes.mode)
    case {'environs_lakes_catalogue', 'table'}
        require_fields(cfg.lakes, {'file', 'variable'}, 'cfg.lakes');
    case 'manual'
        require_fields(cfg.lakes, {'x_m', 'y_m', 'volume_m3', ...
            'duration_days'}, 'cfg.lakes');
        if ~isfield(cfg.lakes, 'drainage_day') && ...
                ~isfield(cfg.lakes, 'drainage_after_start_days')
            error('nevis_validate_region_config:MissingField', ...
                'Missing required config field: cfg.lakes.drainage_day or cfg.lakes.drainage_after_start_days');
        end
    otherwise
        error('nevis_validate_region_config:UnsupportedLakes', ...
            'Unsupported cfg.lakes.mode: %s', cfg.lakes.mode);
end

switch lower(cfg.stations.mode)
    case 'file'
        require_fields(cfg.stations, {'file', 'variable'}, 'cfg.stations');
    case 'coordinates'
        require_fields(cfg.stations, {'x_m', 'y_m'}, 'cfg.stations');
    case 'none'
    otherwise
        error('nevis_validate_region_config:UnsupportedStations', ...
            'Unsupported cfg.stations.mode: %s', cfg.stations.mode);
end

end

function require_fields(s, names, label)
for i = 1:numel(names)
    name = names{i};
    if ~isfield(s, name) || isempty(s.(name))
        error('nevis_validate_region_config:MissingField', ...
            'Missing required config field: %s.%s', label, name);
    end
end
end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end

function assert_manifest_file(dataset, section, file, casename)
if ~isfield(dataset, section) || ~isstruct(dataset.(section)) || ...
        ~isfield(dataset.(section), 'file') || isempty(dataset.(section).file)
    error('nevis_validate_region_config:DatasetManifestMissingSection', ...
        'Dataset manifest for case "%s" does not declare section "%s".', ...
        casename, section);
end
declared = dataset.(section).file;
normalized_file = normalize_path_text(file);
if strcmp(normalized_file, normalize_path_text(declared))
    return
end
if isfield(dataset.(section), 'files')
    declared_files = dataset.(section).files;
    for i = 1:numel(declared_files)
        if strcmp(normalized_file, normalize_path_text(declared_files{i}))
            return
        end
    end
end
if isfield(dataset.(section), 'directory') && ...
        startsWith(normalized_file, [normalize_path_text(dataset.(section).directory) '/'])
    return
end
if ~strcmp(normalized_file, normalize_path_text(declared))
    error('nevis_validate_region_config:DatasetManifestMismatch', ...
        ['Config %s file "%s" does not match dataset manifest file "%s" ' ...
         'for case "%s".'], section, file, declared, casename);
end
end

function value = normalize_path_text(value)
value = strrep(value, '\', '/');
end
