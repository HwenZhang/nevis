function case_file = nevis_write_regional_forward_case(workflow, case_dir, template_case)
%NEVIS_WRITE_REGIONAL_FORWARD_CASE Materialize a per-case forward config.
%
% Writes case_dir/<workflow.forward.case>.m containing the template content
% inlined verbatim followed by this workflow's overrides. The result is a
% fully self-contained case config: editing it changes what the forward
% model runs, with no further dependency on the template.
%
% Idempotent: existing files are left alone so user edits survive.

case_file = fullfile(case_dir, [workflow.forward.case, '.m']);
if exist(case_file, 'file') == 2 && ~workflow_force_overwrite(workflow)
    return
end

template_file = resolve_template_path(workflow, template_case);
template_text = fileread(template_file);

ensure_dir(case_dir);
fid = fopen(case_file, 'w');
if fid == -1
    error('nevis_write_regional_forward_case:FileWrite', ...
        'Could not write forward case config: %s', case_file);
end
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, '%% Per-case forward regional config.\n');
fprintf(fid, '%% Generated from template: %s\n', template_file);
fprintf(fid, '%% Self-contained: edit freely. Subsequent workflow runs do not\n');
fprintf(fid, '%% overwrite this file. Delete it to regenerate from the template.\n\n');
fprintf(fid, '%% ====== inlined template (%s) ======\n\n', template_case);
fprintf(fid, '%s', template_text);
if ~endsWith(template_text, newline)
    fprintf(fid, '\n');
end
fprintf(fid, '\n%% ====== workflow overrides ======\n');
fprintf(fid, 'cfg.casename = ''%s'';\n', workflow.forward.case);
fprintf(fid, 'cfg.dataset.name = ''%s'';\n', workflow.dataset.name);
fprintf(fid, 'cfg.dataset.root = ''%s'';\n', workflow.dataset.root);
fprintf(fid, 'cfg.dataset.manifest = ''%s'';\n', workflow.dataset.manifest);
fprintf(fid, 'cfg.moulins.file = ''moulins_2022_coordinates.mat'';\n');
fprintf(fid, 'cfg.moulins.variable = ''moulins'';\n');
fprintf(fid, 'cfg.moulins.mode = ''coordinates'';\n');
fprintf(fid, 'cfg.moulins.x_field = ''x_m'';\n');
fprintf(fid, 'cfg.moulins.y_field = ''y_m'';\n');
fprintf(fid, 'cfg.lakes.file = ''environs_lakes_2022B_250416.mat'';\n');
fprintf(fid, 'cfg.inversion.mode = ''file'';\n');
fprintf(fid, 'cfg.inversion.file = ''%s'';\n', workflow.inversion.result_file);
fprintf(fid, '\n');
fprintf(fid, '%% Forward starts from the hydrology-only spinup timestep, not\n');
fprintf(fid, '%% from the inversion helper file. This restores the saved hydrology\n');
fprintf(fid, '%% state (phi/hs) from the spinup result.\n');
fprintf(fid, 'cfg.initial_hydrology.mode = ''result_timestep'';\n');
fprintf(fid, 'cfg.initial_hydrology.file = '''';\n');
fprintf(fid, 'cfg.initial_hydrology.variable = '''';\n');
fprintf(fid, 'cfg.initial_hydrology.result_case = ''%s'';\n', ...
    workflow.spinup.case);
fprintf(fid, 'cfg.initial_hydrology.timestep_file = ''%s'';\n', ...
    workflow.spinup.timestep_file);
fprintf(fid, '\n');
fprintf(fid, '%% The coupled ice velocity initial state and Dirichlet velocity\n');
fprintf(fid, '%% boundaries come from the inversion output.\n');
fprintf(fid, 'cfg.initial_ice_velocity.mode = ''file'';\n');
fprintf(fid, 'cfg.initial_ice_velocity.file = ''%s'';\n', ...
    workflow.inversion.hydrology_file);
fprintf(fid, 'cfg.initial_ice_velocity.variable = ''vv_hydro'';\n');
end

function tf = workflow_force_overwrite(workflow)
tf = isfield(workflow, 'overwrite') && isscalar(workflow.overwrite) && ...
    workflow.overwrite;
end

function path = resolve_template_path(workflow, template_case)
root = './';
if isfield(workflow, 'root') && ~isempty(workflow.root)
    root = workflow.root;
end
path = fullfile(root, 'cases', 'templates', [template_case, '.m']);
if exist(path, 'file') ~= 2
    error('nevis_write_regional_forward_case:MissingTemplate', ...
        'Forward template case not found: %s', path);
end
end

function ensure_dir(path_in)
if exist(path_in, 'dir') ~= 7
    mkdir(path_in);
end
end
