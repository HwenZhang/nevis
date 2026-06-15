function outputs = nevis_export_spinup_for_inversion(cfg)
%NEVIS_EXPORT_SPINUP_FOR_INVERSION Export inversion inputs from a spinup run.
%
% This function does NOT run the spinup. nevis_run_inversion calls it when
% file-mode inversion inputs are missing. It only reads a spinup that has
% already been run with:
%   nevis_run_regional('<spinup_case>')
%
% The existing spinup output provides two different inversion inputs:
%   1. source state: pp, pd, ps, gg, aa, oo on the exact inversion grid
%   2. initial hydrology: vv.N/vv.phi from a selected spinup timestep
%
% If the spinup output/timestep is missing, this function errors instead of
% falling back to a legacy or default file.

require_spinup_config(cfg);

source_cfg = cfg.source;
hydro_cfg = cfg.initial_hydrology;
spinup_cfg = get_field(cfg, 'spinup', struct);
spinup_case = get_field(spinup_cfg, 'casename', '');
if isempty(spinup_case)
    error('nevis_export_spinup_for_inversion:MissingSpinupCase', ...
        'cfg.spinup.casename is required.');
end

spinup_run_file = resolve_spinup_run_file(source_cfg, cfg, spinup_case);
spinup = load(spinup_run_file);
required_state = {'pp', 'pd', 'ps', 'gg', 'aa', 'oo'};
for i = 1:numel(required_state)
    if ~isfield(spinup, required_state{i})
        error('nevis_export_spinup_for_inversion:MissingSpinupVariable', ...
            'Spinup output "%s" must contain variable %s.', ...
            spinup_run_file, required_state{i});
    end
end
if ~isfield(spinup, 'vv')
    spinup.vv = struct;
end
if ~isfield(spinup, 'region')
    spinup.region = struct;
end

state_file = resolve_dataset_or_project_path(source_cfg.state_file, cfg);
ensure_parent_dir(state_file);
pp = spinup.pp;
pd = spinup.pd;
ps = spinup.ps;
gg = spinup.gg;
aa = spinup.aa;
oo = spinup.oo;
vv = spinup.vv;
region = spinup.region;
source_metadata = struct('mode', 'spinup_export', ...
    'spinup_case', spinup_case, 'run_file', spinup_run_file);
save(state_file, 'pp', 'pd', 'ps', 'gg', 'aa', 'oo', 'vv', ...
    'region', 'source_metadata');

hydrology_file = '';
if isfield(hydro_cfg, 'file') && ~isempty(hydro_cfg.file)
    hydrology_file = resolve_dataset_or_project_path(hydro_cfg.file, cfg);
    timestep_file = resolve_spinup_timestep_file(spinup_cfg, cfg, spinup_case);
    timestep = load(timestep_file, 'vv');
    if ~isfield(timestep, 'vv')
        error('nevis_export_spinup_for_inversion:MissingTimestepVV', ...
            'Spinup timestep "%s" must contain variable vv.', timestep_file);
    end

    vv = build_initial_hydrology(timestep.vv, aa, gg);
    vv.source = struct('mode', 'spinup_export', ...
        'spinup_case', spinup_case, 'timestep_file', timestep_file, ...
        'source_state_file', state_file);
    vv_hydro = vv;
    ensure_parent_dir(hydrology_file);
    save(hydrology_file, 'vv', 'vv_hydro', 'cfg');
end

outputs = struct;
outputs.spinup_case = spinup_case;
outputs.spinup_run_file = spinup_run_file;
outputs.source_state_file = state_file;
outputs.initial_hydrology_file = hydrology_file;
end

function require_spinup_config(cfg)
if ~isfield(cfg, 'spinup') || ~isstruct(cfg.spinup)
    error('nevis_export_spinup_for_inversion:MissingSpinup', ...
        'cfg.spinup.casename and cfg.spinup.timestep_file are required.');
end
if ~isfield(cfg, 'source') || ~isstruct(cfg.source)
    error('nevis_export_spinup_for_inversion:MissingSource', ...
        'cfg.source is required.');
end
if ~isfield(cfg.source, 'state_file') || isempty(cfg.source.state_file)
    error('nevis_export_spinup_for_inversion:MissingStateFile', ...
        'cfg.source.state_file is required.');
end
if ~isfield(cfg, 'initial_hydrology') || ~isstruct(cfg.initial_hydrology)
    error('nevis_export_spinup_for_inversion:MissingInitialHydrology', ...
        'cfg.initial_hydrology is required.');
end
if ~isfield(cfg.spinup, 'casename') || isempty(cfg.spinup.casename)
    error('nevis_export_spinup_for_inversion:MissingSpinupCase', ...
        'cfg.spinup.casename is required.');
end
if ~isfield(cfg.spinup, 'timestep_file') || isempty(cfg.spinup.timestep_file)
    error('nevis_export_spinup_for_inversion:MissingTimestepFile', ...
        'cfg.spinup.timestep_file is required.');
end
end

function vv_out = build_initial_hydrology(vv_in, aa, gg)
vv_out = struct;
if isfield(vv_in, 'N')
    N = vv_in.N(:);
elseif isfield(vv_in, 'phi')
    phi = vv_in.phi(:);
    if numel(phi) ~= gg.nIJ
        error('nevis_export_spinup_for_inversion:PhiSizeMismatch', ...
            'Spinup vv.phi has %d values; expected %d.', numel(phi), gg.nIJ);
    end
    N = aa.phi_0(:) - phi;
    vv_out.phi = phi;
else
    error('nevis_export_spinup_for_inversion:MissingHydrologyState', ...
        'Spinup timestep vv must contain N or phi.');
end

if numel(N) ~= gg.nIJ
    error('nevis_export_spinup_for_inversion:NSizeMismatch', ...
        'Spinup effective pressure has %d values; expected %d.', ...
        numel(N), gg.nIJ);
end
N(~isfinite(N)) = 0;
vv_out.N = N;
if ~isfield(vv_out, 'phi')
    vv_out.phi = aa.phi_0(:) - N;
end
if isfield(vv_in, 'u') && numel(vv_in.u) == gg.eIJ
    vv_out.u = vv_in.u;
end
if isfield(vv_in, 'v') && numel(vv_in.v) == gg.fIJ
    vv_out.v = vv_in.v;
end
end

function path_out = resolve_spinup_run_file(source_cfg, cfg, spinup_case)
if isfield(source_cfg, 'result_file') && ~isempty(source_cfg.result_file)
    path_out = resolve_project_path(source_cfg.result_file, cfg);
else
    path_out = fullfile(cfg.root, cfg.results, spinup_case, ...
        [spinup_case, '.mat']);
end
if exist(path_out, 'file') ~= 2
    error('nevis_export_spinup_for_inversion:MissingSpinupResult', ...
        ['Spinup output not found: %s\nRun first: ' ...
         'nevis_run_regional(''%s'')'], path_out, spinup_case);
end
end

function path_out = resolve_spinup_timestep_file(spinup_cfg, cfg, spinup_case)
path_in = spinup_cfg.timestep_file;
if is_absolute_path(path_in) || startsWith(path_in, './') || ...
        startsWith(path_in, ['..' filesep])
    path_out = resolve_project_path(path_in, cfg);
else
    path_out = fullfile(cfg.root, cfg.results, spinup_case, path_in);
end
if exist(path_out, 'file') ~= 2
    error('nevis_export_spinup_for_inversion:MissingSpinupTimestep', ...
        ['Spinup timestep not found: %s\nRun first: ' ...
         'nevis_run_regional(''%s'')'], path_out, spinup_case);
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

function ensure_parent_dir(path_in)
parent = fileparts(path_in);
if ~isempty(parent) && exist(parent, 'dir') ~= 7
    mkdir(parent);
end
end

function value = get_field(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
