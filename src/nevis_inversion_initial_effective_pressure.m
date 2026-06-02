function [N_obs, source] = nevis_inversion_initial_effective_pressure(initial_hydrology, aa, gg, cfg)
%NEVIS_INVERSION_INITIAL_EFFECTIVE_PRESSURE Build initial N for inversion.
%
% Supported modes:
%   file       - load variable from a MAT file; accepts .N or .phi
%   k_factor   - set phi = phi_a + k_factor * (phi_0 - phi_a)
%   constant_N - set a uniform dimensionless N value

mode = lower(initial_hydrology.mode);
switch mode
    case 'file'
        if ~isfield(initial_hydrology, 'file') || isempty(initial_hydrology.file)
            error('nevis_inversion_initial_effective_pressure:MissingFile', ...
                'cfg.initial_hydrology.file is required when mode="file".');
        end
        if ~isfield(initial_hydrology, 'variable') || isempty(initial_hydrology.variable)
            error('nevis_inversion_initial_effective_pressure:MissingVariable', ...
                'cfg.initial_hydrology.variable is required when mode="file".');
        end

        init_file = resolve_project_path(initial_hydrology.file, cfg);
        if exist(init_file, 'file') ~= 2
            error('nevis_inversion_initial_effective_pressure:MissingFile', ...
                'Configured initial hydrology file not found: %s', init_file);
        end

        init_var = initial_hydrology.variable;
        data = load(init_file, init_var);
        if ~isfield(data, init_var)
            error('nevis_inversion_initial_effective_pressure:MissingVariable', ...
                'Initial hydrology file "%s" does not contain variable "%s".', ...
                init_file, init_var);
        end
        vv = data.(init_var);
        if isstruct(vv) && isfield(vv, 'N')
            N_obs = vv.N;
        elseif isstruct(vv) && isfield(vv, 'phi')
            N_obs = aa.phi_0 - vv.phi;
        else
            error('nevis_inversion_initial_effective_pressure:UnsupportedFileVariable', ...
                ['Initial hydrology variable "%s" must be a struct with field ' ...
                 'N or phi.'], init_var);
        end
        source = struct('mode', mode, 'file', init_file, 'variable', init_var);

    case 'k_factor'
        if ~isfield(initial_hydrology, 'k_factor') || isempty(initial_hydrology.k_factor)
            error('nevis_inversion_initial_effective_pressure:MissingKFactor', ...
                'cfg.initial_hydrology.k_factor is required when mode="k_factor".');
        end
        k_factor = initial_hydrology.k_factor;
        if k_factor < 0 || k_factor > 1
            error('nevis_inversion_initial_effective_pressure:InvalidKFactor', ...
                'cfg.initial_hydrology.k_factor must be between 0 and 1.');
        end
        phi = aa.phi_a + k_factor * (aa.phi_0 - aa.phi_a);
        N_obs = aa.phi_0 - phi;
        source = struct('mode', mode, 'k_factor', k_factor);

    case {'constant_n', 'uniform_n'}
        if ~isfield(initial_hydrology, 'N') || isempty(initial_hydrology.N)
            error('nevis_inversion_initial_effective_pressure:MissingN', ...
                'cfg.initial_hydrology.N is required when mode="constant_N".');
        end
        N_obs = initial_hydrology.N * ones(gg.nIJ, 1);
        source = struct('mode', mode, 'N', initial_hydrology.N);

    otherwise
        error('nevis_inversion_initial_effective_pressure:UnsupportedMode', ...
            ['Unsupported cfg.initial_hydrology.mode "%s". Use "file", ' ...
             '"k_factor", or "constant_N".'], initial_hydrology.mode);
end

N_obs = N_obs(:);
if numel(N_obs) ~= gg.nIJ
    error('nevis_inversion_initial_effective_pressure:SizeMismatch', ...
        'Initial effective pressure has %d values, but grid has %d nodes.', ...
        numel(N_obs), gg.nIJ);
end
N_obs(~isfinite(N_obs)) = 0;
end

function path_out = resolve_project_path(path_in, cfg)
    if startsWith(path_in, filesep) || ...
            ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]', 'once'))
        path_out = path_in;
    else
        path_out = fullfile(cfg.root, path_in);
    end
end
