function [vv, aa, region] = nevis_import_region_ice_velocity(cfg, vv, aa, gg, region)
%NEVIS_IMPORT_REGION_ICE_VELOCITY Optionally override ice velocity state.
%
% This is applied after observed velocity import and after hydrology-state
% restore. It is intended for forward runs that should keep hydrology from a
% spinup timestep but initialize/prescribe SSA Dirichlet velocity from an
% inverted velocity product.

if ~isfield(cfg, 'initial_ice_velocity') || ~isstruct(cfg.initial_ice_velocity)
    return
end

vcfg = cfg.initial_ice_velocity;
mode = get_opt(vcfg, 'mode', 'none');

switch lower(mode)
    case {'none', 'keep'}
        return
    case 'file'
        velocity_file = nevis_region_resolve_path(vcfg.file, cfg);
        if exist(velocity_file, 'file') ~= 2
            error('nevis_import_region_ice_velocity:MissingFile', ...
                'Initial ice velocity file not found: %s', velocity_file);
        end
        velocity_state = nevis_region_load_var(velocity_file, ...
            get_opt(vcfg, 'variable', ''));
    otherwise
        error('nevis_import_region_ice_velocity:UnsupportedMode', ...
            'Unsupported cfg.initial_ice_velocity.mode: %s', mode);
end

if ~isfield(velocity_state, 'u') || ~isfield(velocity_state, 'v')
    error('nevis_import_region_ice_velocity:MissingVelocity', ...
        'Initial ice velocity state must contain fields u and v.');
end
if numel(velocity_state.u) ~= gg.eIJ
    error('nevis_import_region_ice_velocity:USizeMismatch', ...
        'Initial u size (%d) does not match model eIJ (%d).', ...
        numel(velocity_state.u), gg.eIJ);
end
if numel(velocity_state.v) ~= gg.fIJ
    error('nevis_import_region_ice_velocity:VSizeMismatch', ...
        'Initial v size (%d) does not match model fIJ (%d).', ...
        numel(velocity_state.v), gg.fIJ);
end

vv.u = velocity_state.u(:);
vv.v = velocity_state.v(:);
vv.u(gg.eout2) = NaN;
vv.v(gg.fout2) = NaN;

% nevis_backbone uses vv.u/v on gg.ebdy2/fbdy2 for Dirichlet velocity
% boundaries. Keep aa.u_obs/v_obs aligned for diagnostics and downstream
% helpers that inspect the prescribed velocity field.
aa.u_obs = vv.u;
aa.v_obs = vv.v;

region.initial_ice_velocity = velocity_state;
region.initial_ice_velocity_file = velocity_file;
end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
