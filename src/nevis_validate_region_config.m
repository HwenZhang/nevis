function cfg = nevis_validate_region_config(cfg)
%NEVIS_VALIDATE_REGION_CONFIG Validate the strict per-case data contract.

required_top = {'casename', 'data_root', 'oo', 'pd', 'pp', 'geometry', 'velocity', ...
    'inversion', 'initial_hydrology', 'moulins', 'lakes', 'runoff', ...
    'stations', 'run'};
require_fields(cfg, required_top, 'cfg');

if exist(cfg.data_root, 'dir') ~= 7
    error('nevis_validate_region_config:MissingDataRoot', ...
        'Case data_root does not exist: %s', cfg.data_root);
end

require_fields(cfg.geometry, {'file', 'variable', 'x_field', 'y_field', ...
    'bed_field', 'surface_field', 'skip'}, 'cfg.geometry');
require_fields(cfg.velocity, {'file', 'variable', 'u_field', 'v_field', ...
    'units', 'mode'}, 'cfg.velocity');
require_fields(cfg.inversion, {'file', 'mode'}, 'cfg.inversion');
require_fields(cfg.initial_hydrology, {'file', 'variable', 'mode'}, ...
    'cfg.initial_hydrology');
require_fields(cfg.moulins, {'file', 'variable', 'mode'}, 'cfg.moulins');
require_fields(cfg.lakes, {'file', 'variable', 'mode'}, 'cfg.lakes');
require_fields(cfg.runoff, {'file', 'variable', 'mode'}, 'cfg.runoff');
require_fields(cfg.stations, {'file', 'variable', 'mode'}, 'cfg.stations');
require_fields(cfg.run, {'dt_days', 't_span_days', 'save_timesteps', ...
    'save_pts_all'}, 'cfg.run');

require_fields(cfg.oo, {'evaluate_variables', 'include_ice', ...
    'adjust_boundaries', 'input_gaussian', 'relaxation_term', ...
    'initial_condition', 'mean_perms', 'modified_mean_perms', ...
    'display_residual', 'visualize_vel', 'N_coupling', ...
    'U_coupling', 'boundary_method', 'mask_boundary_method', ...
    'plot_residual', 'max_iter_new', 'partition_ratio', ...
    'iterative_solver'}, 'cfg.oo');
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
    if strcmpi(mode, 'none')
        continue
    end
    path = nevis_region_resolve_path(scfg.file, cfg);
    if exist(path, 'file') ~= 2
        error('nevis_validate_region_config:MissingDataFile', ...
            'Configured %s file does not exist for case "%s": %s', ...
            section, cfg.casename, path);
    end
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
