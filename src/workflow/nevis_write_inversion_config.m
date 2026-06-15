function config_file = nevis_write_inversion_config(workflow, template_config)
%NEVIS_WRITE_INVERSION_CONFIG Materialize a per-case inversion config as .m.
%
% Writes workflow.inversion.config_file (an .m file under the case_dir)
% containing the template content inlined verbatim followed by this
% workflow's overrides. Self-contained: editing it changes what the
% inversion runs, with no further dependency on the template.
%
% Idempotent: existing files are left alone so user edits survive.

config_file = workflow.inversion.config_file;
if exist(config_file, 'file') == 2 && ~workflow_force_overwrite(workflow)
    return
end

if exist(template_config, 'file') ~= 2
    error('nevis_write_inversion_config:MissingTemplate', ...
        'Inversion template config not found: %s', template_config);
end
template_text = fileread(template_config);

ensure_parent_dir(config_file);
fid = fopen(config_file, 'w');
if fid == -1
    error('nevis_write_inversion_config:FileWrite', ...
        'Could not write inversion config: %s', config_file);
end
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, '%% Per-case inversion config.\n');
fprintf(fid, '%% Generated from template: %s\n', template_config);
fprintf(fid, '%% Self-contained: edit freely. Subsequent workflow runs do not\n');
fprintf(fid, '%% overwrite this file. Delete it to regenerate from the template.\n\n');
fprintf(fid, '%% ====== inlined template ======\n\n');
fprintf(fid, '%s', template_text);
if ~endsWith(template_text, newline)
    fprintf(fid, '\n');
end
fprintf(fid, '\n%% ====== workflow overrides ======\n');
fprintf(fid, 'cfg.name = ''%s'';\n', workflow.inversion.name);
fprintf(fid, 'cfg.sliding_law = ''%s'';\n', workflow.inversion.sliding_law);
fprintf(fid, 'cfg.dataset.name = ''%s'';\n', workflow.dataset.name);
fprintf(fid, 'cfg.dataset.root = ''%s'';\n', workflow.dataset.root);
fprintf(fid, 'cfg.dataset.manifest = ''%s'';\n', workflow.dataset.manifest);
fprintf(fid, 'cfg.spinup.casename = ''%s'';\n', workflow.spinup.case);
fprintf(fid, 'cfg.spinup.timestep_file = ''%s'';\n', workflow.spinup.timestep_file);
fprintf(fid, 'cfg.source.casename = ''%s_source_from_spinup'';\n', workflow.dataset.name);
fprintf(fid, 'cfg.source.state_file = ''%s'';\n', workflow.inversion.source_state_file);
fprintf(fid, 'cfg.initial_hydrology.mode = ''file'';\n');
fprintf(fid, 'cfg.initial_hydrology.file = ''%s'';\n', workflow.inversion.initial_hydrology_file);
fprintf(fid, 'cfg.initial_hydrology.variable = ''vv'';\n');
fprintf(fid, 'cfg.output.inversion_file = ''%s'';\n', workflow.inversion.result_file);
fprintf(fid, 'cfg.output.initial_hydrology_file = ''%s'';\n', workflow.inversion.hydrology_file);
fprintf(fid, 'cfg.output.intermediate_prefix = ''%s'';\n', workflow.inversion.intermediate_prefix);
fprintf(fid, 'cfg.output.intermediate_hydrology_file = ''%s'';\n', workflow.inversion.intermediate_hydrology_file);
end

function tf = workflow_force_overwrite(workflow)
tf = isfield(workflow, 'overwrite') && isscalar(workflow.overwrite) && ...
    workflow.overwrite;
end

function ensure_parent_dir(path_in)
parent = fileparts(path_in);
if ~isempty(parent) && exist(parent, 'dir') ~= 7
    mkdir(parent);
end
end
