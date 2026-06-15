function [hydrology_file, vv, source] = nevis_build_inversion_initial_hydrology(config, hydrology_file)
%NEVIS_BUILD_INVERSION_INITIAL_HYDROLOGY Write an inversion initial N file.
%
% The generated file contains vv.N and vv.phi on the inversion source grid.
% It can be consumed by cfg.initial_hydrology.mode='file' with variable='vv',
% or by a regional forward case after an inversion writes vv_hydro.

cfg = load_config(config);

if nargin < 2 || isempty(hydrology_file)
    hydrology_file = default_hydrology_file(cfg);
end

state_file = resolve_project_path(cfg.source.state_file, cfg);
if exist(state_file, 'file') ~= 2
    error('nevis_build_inversion_initial_hydrology:MissingSourceState', ...
        'Source state file not found: %s', state_file);
end

state = load(state_file, 'aa', 'gg');
if ~isfield(state, 'aa') || ~isfield(state, 'gg')
    error('nevis_build_inversion_initial_hydrology:InvalidSourceState', ...
        'Source state must contain aa and gg: %s', state_file);
end

generator = get_initial_hydrology_generator(cfg);
[N, source] = nevis_inversion_initial_effective_pressure( ...
    generator, state.aa, state.gg, cfg);

vv = struct;
vv.N = N(:);
vv.phi = state.aa.phi_0(:) - vv.N;
vv.source = source;
vv.source_state_file = state_file;

ensure_parent_dir(hydrology_file);
save(hydrology_file, 'vv', 'cfg', 'source');
end

function generator = get_initial_hydrology_generator(cfg)
if isfield(cfg.initial_hydrology, 'generator') && ...
        isstruct(cfg.initial_hydrology.generator)
    generator = cfg.initial_hydrology.generator;
else
    generator = cfg.initial_hydrology;
end
end

function cfg = load_config(config)
if isstruct(config)
    cfg = config;
    return
end
if exist(config, 'file') ~= 2
    error('nevis_build_inversion_initial_hydrology:MissingConfig', ...
        'Config not found: %s', config);
end
[~, ~, ext] = fileparts(config);
switch lower(ext)
    case '.m'
        run(config);
        if exist('cfg', 'var') ~= 1
            error('nevis_build_inversion_initial_hydrology:MissingCfg', ...
                'M-file config must assign variable cfg.');
        end
    case '.mat'
        data = load(config, 'cfg');
        if ~isfield(data, 'cfg')
            error('nevis_build_inversion_initial_hydrology:MissingCfg', ...
                'MAT config must contain variable cfg.');
        end
        cfg = data.cfg;
    case '.json'
        cfg = jsondecode(fileread(config));
    otherwise
        error('nevis_build_inversion_initial_hydrology:UnsupportedConfig', ...
            'Unsupported config extension: %s', ext);
end
end

function hydrology_file = default_hydrology_file(cfg)
hydrology_file = fullfile(cfg.dataset.root, 'initial_hydrology', ...
    [cfg.source.casename, '_initial_hydrology.mat']);
end

function path_out = resolve_project_path(path_in, cfg)
if startsWith(path_in, filesep) || ...
        ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]', 'once'))
    path_out = path_in;
else
    path_out = fullfile(cfg.root, path_in);
end
end

function ensure_parent_dir(path_in)
parent = fileparts(path_in);
if ~isempty(parent) && exist(parent, 'dir') ~= 7
    mkdir(parent);
end
end
