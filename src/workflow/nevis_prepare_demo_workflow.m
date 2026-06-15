function workflow = nevis_prepare_demo_workflow(opts)
%NEVIS_PREPARE_DEMO_WORKFLOW Write configs for the demo workflow.
%
% The returned workflow struct contains the generated case names and config
% file paths needed by test.m. This function does not build data, run spinup,
% run inversion, or run the forward model.

opts = fill_missing_options(opts);
assert_valid_case_name(opts.case_name);

workflow = struct;
workflow.root = opts.root;
workflow.case_name = opts.case_name;
workflow.bbox_km = opts.bbox_km;
workflow.raw = opts.raw;

workflow.dataset.name = opts.dataset.name;
workflow.dataset.root = opts.dataset.root;
workflow.dataset.manifest = opts.dataset.manifest;

workflow.spinup.case = [opts.case_name, '_inversion_spinup'];
workflow.spinup.template_case = opts.spinup.template_case;
workflow.spinup.timestep_file = opts.spinup.timestep_file;
workflow.spinup.t_span_days = opts.spinup.t_span_days;

workflow.inversion.name = [opts.inversion.base_name, '_', ...
    opts.case_name, '_spinup_initialized'];
workflow.inversion.config_file = fullfile(opts.paths.case_dir, ...
    [opts.case_name, '_inversion.m']);
workflow.inversion.sliding_law = opts.inversion.sliding_law;
workflow.inversion.source_state_file = fullfile('source_states', ...
    [workflow.spinup.case, '_source.mat']);
workflow.inversion.initial_hydrology_file = fullfile('initial_hydrology', ...
    [workflow.spinup.case, '_', timestep_label(opts.spinup.timestep_file), ...
     '.mat']);
workflow.inversion.result_file = fullfile('inversions', ...
    workflow.inversion.name, 'C_inversion_results.mat');
workflow.inversion.hydrology_file = fullfile('inversions', ...
    workflow.inversion.name, 'velocity_inverted.mat');
workflow.inversion.intermediate_prefix = fullfile('inversions', ...
    workflow.inversion.name, 'C_inversion_outer');
workflow.inversion.intermediate_hydrology_file = fullfile('inversions', ...
    workflow.inversion.name, 'velocity_inverted_outer.mat');

workflow.forward.case = [opts.case_name, '_forward'];

workflow.drainage.case = [opts.case_name, opts.drainage.suffix];
workflow.drainage.template_case = opts.drainage.template_case;
workflow.drainage.initial_result_case = workflow.forward.case;
workflow.drainage.initial_timestep_file = opts.drainage.initial_timestep_file;
workflow.drainage.lake_catalogue_file = opts.drainage.lake_catalogue_file;
workflow.drainage.t_span_days = opts.drainage.t_span_days;

if isfield(opts, 'reinit_case') && ~isempty(opts.reinit_case)
    workflow.overwrite = opts.reinit_case;
else
    workflow.overwrite = false;
end

ensure_dir(opts.paths.cases_dir);
ensure_dir(opts.paths.case_dir);
ensure_dir(workflow.dataset.root);
nevis_write_regional_spinup_case(workflow, opts.paths.case_dir, ...
    opts.spinup.template_case);
nevis_write_regional_forward_case(workflow, opts.paths.case_dir, ...
    opts.forward.template_case);
nevis_write_regional_drainage_case(workflow, opts.paths.case_dir, ...
    opts.drainage.template_case);
nevis_write_inversion_config(workflow, opts.inversion.template_config);
end

function opts = fill_missing_options(opts)
defaults = nevis_demo_options;
opts = merge_struct(defaults, opts);
end

function out = merge_struct(defaults, values)
out = defaults;
names = fieldnames(values);
for i = 1:numel(names)
    name = names{i};
    if isstruct(values.(name)) && isfield(defaults, name) && ...
            isstruct(defaults.(name))
        out.(name) = merge_struct(defaults.(name), values.(name));
    else
        out.(name) = values.(name);
    end
end
end

function assert_valid_case_name(case_name)
if isempty(regexp(case_name, '^[A-Za-z]\w*$', 'once'))
    error('nevis_prepare_demo_workflow:InvalidCaseName', ...
        ['case_name must be a valid MATLAB-style identifier containing ' ...
         'letters, numbers, and underscores: %s'], case_name);
end
end

function label = timestep_label(timestep_file)
[~, label, ext] = fileparts(timestep_file);
if isempty(label) && ~isempty(ext)
    label = erase(ext, '.');
end
end

function ensure_dir(path_in)
if exist(path_in, 'dir') ~= 7
    mkdir(path_in);
end
end
