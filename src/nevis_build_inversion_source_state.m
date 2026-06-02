function [state_file, state] = nevis_build_inversion_source_state(config, state_file)
%NEVIS_BUILD_INVERSION_SOURCE_STATE Build a static source state for inversion.
%
% The inversion only needs pp, pd, ps, gg, aa, and oo on the same grid as the
% dataset velocity/geometry. This builder reconstructs that static state from
% the dataset geometry without running hydrology timesteps.

cfg = load_config(config);

if nargin < 2 || isempty(state_file)
    state_file = default_state_file(cfg);
end

root = cfg.root;
code = cfg.code;
addpath(code);

oo = struct;
oo.root = root;
oo.code = code;
oo.results = cfg.results;
oo.casename = cfg.source.casename;
oo.fn = ['/', oo.casename];
oo.rn = fullfile(oo.root, oo.results, oo.casename);
oo.dataset = cfg.dataset.name;
oo.dataset_root = cfg.dataset.root;
oo.dn = cfg.dataset.root;

source_cfg = build_source_cfg(cfg);

[pd, oo] = nevis_defaults([], oo);
[pd, oo] = nevis_apply_config_parameters(source_cfg, pd, oo);
ps = struct;
[ps, pp] = nevis_nondimension(pd, ps, oo);
[~, ~, pp] = nevis_apply_config_parameters(source_cfg, pd, oo, pp);

[gg, b, s, region, oo] = nevis_import_region_geometry(source_cfg, ps, oo);
[pd, ps, pp, oo] = nevis_update_parameters_ice(pd, ps, pp, oo);
[pd, oo, pp] = nevis_apply_config_parameters(source_cfg, pd, oo, pp);
gg = nevis_label_ice_test(gg, oo);
[aa, vv] = nevis_initialize(b, s, gg, pp, oo);
aa.phi_b = max(aa.phi_0, aa.phi_a);

state = struct;
state.cfg = cfg;
state.source_cfg = source_cfg;
state.pd = pd;
state.ps = ps;
state.pp = pp;
state.gg = gg;
state.aa = aa;
state.vv = vv;
state.oo = oo;
state.region = region;
state.state_file = state_file;

ensure_parent_dir(state_file);
save(state_file, 'pp', 'pd', 'ps', 'gg', 'aa', 'oo', 'vv', 'region', 'cfg', 'source_cfg');
end

function source_cfg = build_source_cfg(cfg)
source_cfg = struct;
source_cfg.casename = cfg.source.casename;
source_cfg.case_root = get_field(cfg, 'case_root', './cases');
source_cfg.dataset = cfg.dataset;
source_cfg.geometry = geometry_cfg(cfg);

source_cfg.oo.boundary_method = cfg.solver.boundary_method;
source_cfg.oo.mask_boundary_method = cfg.solver.mask_boundary_method;
source_cfg.oo.display_norms = cfg.solver.display_norms;
source_cfg.oo.verb = cfg.solver.verb;
source_cfg.oo.partition_ratio = cfg.partition_ratio;
source_cfg.oo.evaluate_variables = 1;
source_cfg.oo.include_ice = 1;
source_cfg.oo.adjust_boundaries = 1;
source_cfg.oo.input_gaussian = 1;
source_cfg.oo.relaxation_term = 1;
source_cfg.oo.initial_condition = 1;
source_cfg.oo.mean_perms = 1;
source_cfg.oo.modified_mean_perms = 0;
source_cfg.oo.display_residual = 0;
source_cfg.oo.N_coupling = 1;
source_cfg.oo.U_coupling = 1;

if isfield(cfg.source, 'oo')
    source_cfg.oo = merge_structs(source_cfg.oo, cfg.source.oo);
end
if isfield(cfg.source, 'pd')
    source_cfg.pd = cfg.source.pd;
end
if isfield(cfg.source, 'pp')
    source_cfg.pp = cfg.source.pp;
end
end

function gcfg = geometry_cfg(cfg)
if isfield(cfg.source, 'geometry') && isstruct(cfg.source.geometry)
    gcfg = cfg.source.geometry;
else
    manifest_file = fullfile(cfg.dataset.root, cfg.dataset.manifest);
    dataset = nevis_load_dataset_manifest(manifest_file);
    if ~isfield(dataset, 'geometry') || ~isstruct(dataset.geometry)
        error('nevis_build_inversion_source_state:MissingGeometryManifest', ...
            'Dataset manifest must define dataset.geometry.');
    end
    gcfg = dataset.geometry;
end

if ~isfield(gcfg, 'x_field') || isempty(gcfg.x_field)
    gcfg.x_field = 'X_m';
end
if ~isfield(gcfg, 'y_field') || isempty(gcfg.y_field)
    gcfg.y_field = 'Y_m';
end
if ~isfield(gcfg, 'bed_field') || isempty(gcfg.bed_field)
    gcfg.bed_field = 'B_m';
end
if ~isfield(gcfg, 'surface_field') || isempty(gcfg.surface_field)
    gcfg.surface_field = 'S_m';
end
if ~isfield(gcfg, 'skip') || isempty(gcfg.skip)
    gcfg.skip = 6;
end
if ~isfield(gcfg, 'min_thickness_m') || isempty(gcfg.min_thickness_m)
    gcfg.min_thickness_m = 0;
end
end

function cfg = load_config(config)
if isstruct(config)
    cfg = config;
    return
end
if exist(config, 'file') ~= 2
    error('nevis_build_inversion_source_state:MissingConfig', ...
        'Config not found: %s', config);
end
[~, ~, ext] = fileparts(config);
switch lower(ext)
    case '.m'
        run(config);
        if exist('cfg', 'var') ~= 1
            error('nevis_build_inversion_source_state:MissingCfg', ...
                'M-file config must assign variable cfg.');
        end
    case '.mat'
        data = load(config, 'cfg');
        if ~isfield(data, 'cfg')
            error('nevis_build_inversion_source_state:MissingCfg', ...
                'MAT config must contain variable cfg.');
        end
        cfg = data.cfg;
    case '.json'
        cfg = jsondecode(fileread(config));
    otherwise
        error('nevis_build_inversion_source_state:UnsupportedConfig', ...
            'Unsupported config extension: %s', ext);
end
end

function value = get_field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end

function merged = merge_structs(base, override)
merged = base;
names = fieldnames(override);
for i = 1:numel(names)
    merged.(names{i}) = override.(names{i});
end
end

function state_file = default_state_file(cfg)
state_file = fullfile(cfg.dataset.root, 'source_states', ...
    [cfg.source.casename, '.mat']);
end

function ensure_parent_dir(path_in)
parent = fileparts(path_in);
if ~isempty(parent) && exist(parent, 'dir') ~= 7
    mkdir(parent);
end
end
