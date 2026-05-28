function cfg = nevis_validate_idealized_config(cfg)
%NEVIS_VALIDATE_IDEALIZED_CONFIG Validate the analytical idealized case cfg.

required_top = {'casename', 'case_root', 'oo', 'pd', 'pp', 'geometry', ...
    'initial_hydrology', 'moulins', 'lakes', 'input', 'run'};
require_fields(cfg, required_top, 'cfg');

if exist(cfg.case_root, 'dir') ~= 7
    error('nevis_validate_idealized_config:MissingCaseRoot', ...
        'Case root does not exist: %s', cfg.case_root);
end

require_fields(cfg.oo, {'evaluate_variables', 'input_gaussian', ...
    'relaxation_term', 'initial_condition', 'mean_perms', ...
    'include_blister', 'include_ice', 'display_residual', ...
    'adjust_boundaries', 'yperiodic', 'xperiodic', ...
    'keep_all_moulins', 'random_moulins', 'random_lakes', ...
    'surface_runoff', 'RACMO_runoff', 'distributed_input', ...
    'save_timesteps', 'save_pts_all'}, 'cfg.oo');
require_fields(cfg.pd, {'alpha_b', 'kappa_b', 'c0', 'mu', 'Ye', 'B', ...
    'E_lapse', 'hb_reg1', 'hb_reg2', 'N_reg1', 'deltap_reg', ...
    'B_reg', 'G', 'melt'}, 'cfg.pd');
require_fields(cfg.geometry, {'mode', 'length_m', 'nx', 'ny', ...
    'min_thickness_m'}, 'cfg.geometry');
if ~isfield(cfg.geometry, 'width_m') && ~isfield(cfg.geometry, 'width_ratio')
    error('nevis_validate_idealized_config:MissingField', ...
        'Missing required config field: cfg.geometry.width_m or cfg.geometry.width_ratio');
end
require_fields(cfg.initial_hydrology, {'mode'}, 'cfg.initial_hydrology');
require_fields(cfg.moulins, {'mode'}, 'cfg.moulins');
require_fields(cfg.lakes, {'mode', 'x_fraction', 'y_fraction', ...
    'volume_m3', 'drainage_day', 'duration_days'}, 'cfg.lakes');
require_fields(cfg.input, {'runoff_max_mm_per_day', ...
    'moulin_input_m3_per_s', 'ramp_days'}, 'cfg.input');
require_fields(cfg.run, {'dt_days'}, 'cfg.run');
if ~isfield(cfg.run, 't_span_days') && ~isfield(cfg.run, 't_span_after_start_days')
    error('nevis_validate_idealized_config:MissingField', ...
        'Missing required config field: cfg.run.t_span_days or cfg.run.t_span_after_start_days');
end

switch lower(cfg.geometry.mode)
    case 'hewitt_2013_sqrt'
        require_fields(cfg.geometry, {'bed_elevation_m', ...
            'surface_head_m', 'surface_exponent'}, 'cfg.geometry');
    case 'function'
        require_fields(cfg.geometry, {'bed_function', 'surface_function'}, ...
            'cfg.geometry');
    otherwise
        error('nevis_validate_idealized_config:UnsupportedGeometry', ...
            'Unsupported cfg.geometry.mode: %s', cfg.geometry.mode);
end

switch lower(cfg.initial_hydrology.mode)
    case 'k_factor_cavity_equilibrium'
        require_fields(cfg.initial_hydrology, {'k_factor'}, ...
            'cfg.initial_hydrology');
    case 'function'
        require_fields(cfg.initial_hydrology, {'function'}, ...
            'cfg.initial_hydrology');
    case 'result_timestep'
        require_fields(cfg.initial_hydrology, {'result_case', ...
            'timestep_file'}, 'cfg.initial_hydrology');
    otherwise
        error('nevis_validate_idealized_config:UnsupportedInitialHydrology', ...
            'Unsupported cfg.initial_hydrology.mode: %s', ...
            cfg.initial_hydrology.mode);
end

switch lower(cfg.moulins.mode)
    case 'coordinates_fractional'
        require_fields(cfg.moulins, {'x_fraction', 'y_fraction'}, ...
            'cfg.moulins');
    case 'coordinates_m'
        require_fields(cfg.moulins, {'x_m', 'y_m'}, 'cfg.moulins');
    case 'function'
        require_fields(cfg.moulins, {'position_function'}, 'cfg.moulins');
    otherwise
        error('nevis_validate_idealized_config:UnsupportedMoulins', ...
            'Unsupported cfg.moulins.mode: %s', cfg.moulins.mode);
end

if ~strcmpi(cfg.lakes.mode, 'coordinates_fractional')
    error('nevis_validate_idealized_config:UnsupportedLakes', ...
        'Unsupported cfg.lakes.mode: %s', cfg.lakes.mode);
end

if isfield(cfg, 'stations') && ~isempty(cfg.stations)
    require_fields(cfg.stations, {'mode'}, 'cfg.stations');
    if strcmpi(cfg.stations.mode, 'downstream_lake')
        require_fields(cfg.stations, {'spacing_m'}, 'cfg.stations');
    else
        error('nevis_validate_idealized_config:UnsupportedStations', ...
            'Unsupported cfg.stations.mode: %s', cfg.stations.mode);
    end
end

end

function require_fields(s, names, label)
for i = 1:numel(names)
    name = names{i};
    if ~isfield(s, name) || isempty(s.(name))
        error('nevis_validate_idealized_config:MissingField', ...
            'Missing required config field: %s.%s', label, name);
    end
end
end
