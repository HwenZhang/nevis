function case_file = nevis_write_regional_drainage_case(workflow, case_dir, template_case)
%NEVIS_WRITE_REGIONAL_DRAINAGE_CASE Materialize a per-case drainage config.
%
% The drainage case starts hydrology from the forward run's saved day-365
% timestep using cfg.initial_hydrology.mode='result_timestep', initializes
% ice velocity from the inversion output, then turns the case-local 2022
% lake drainage catalogue on.
%
% Idempotent: existing files are left alone so user edits survive.

case_file = fullfile(case_dir, [workflow.drainage.case, '.m']);
if exist(case_file, 'file') == 2 && ~workflow_force_overwrite(workflow)
    return
end

template_file = resolve_template_path(workflow, template_case);
template_text = fileread(template_file);

ensure_dir(case_dir);
fid = fopen(case_file, 'w');
if fid == -1
    error('nevis_write_regional_drainage_case:FileWrite', ...
        'Could not write drainage case config: %s', case_file);
end
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, '%% Per-case drainage regional config.\n');
fprintf(fid, '%% Generated from template: %s\n', template_file);
fprintf(fid, '%% Self-contained: edit freely. Subsequent workflow runs do not\n');
fprintf(fid, '%% overwrite this file. Delete it to regenerate from the template.\n\n');
fprintf(fid, '%% ====== inlined template (%s) ======\n\n', template_case);
fprintf(fid, '%s', template_text);
if ~endsWith(template_text, newline)
    fprintf(fid, '\n');
end

fprintf(fid, '\n%% ====== workflow drainage overrides ======\n');
fprintf(fid, '%% Starts from the forward run timestep and turns the 2022 lake catalogue on.\n');
fprintf(fid, 'cfg.casename = ''%s'';\n', workflow.drainage.case);
fprintf(fid, 'cfg.dataset.name = ''%s'';\n', workflow.dataset.name);
fprintf(fid, 'cfg.dataset.root = ''%s'';\n', workflow.dataset.root);
fprintf(fid, 'cfg.dataset.manifest = ''%s'';\n', workflow.dataset.manifest);
fprintf(fid, 'cfg.moulins.file = ''moulins_2022_coordinates.mat'';\n');
fprintf(fid, 'cfg.moulins.variable = ''moulins'';\n');
fprintf(fid, 'cfg.moulins.mode = ''coordinates'';\n');
fprintf(fid, 'cfg.moulins.x_field = ''x_m'';\n');
fprintf(fid, 'cfg.moulins.y_field = ''y_m'';\n');
fprintf(fid, 'cfg.lakes.file = ''%s'';\n', workflow.drainage.lake_catalogue_file);
fprintf(fid, 'cfg.lakes.variable = ''environs_lakes'';\n');
fprintf(fid, 'cfg.lakes.mode = ''environs_lakes_catalogue'';\n');
fprintf(fid, 'cfg.lakes.spinup_volume_scale = 1;\n');
fprintf(fid, 'cfg.inversion.mode = ''file'';\n');
fprintf(fid, 'cfg.inversion.file = ''%s'';\n', workflow.inversion.result_file);
fprintf(fid, 'cfg.initial_hydrology.mode = ''result_timestep'';\n');
fprintf(fid, 'cfg.initial_hydrology.file = '''';\n');
fprintf(fid, 'cfg.initial_hydrology.variable = '''';\n');
fprintf(fid, 'cfg.initial_hydrology.result_case = ''%s'';\n', ...
    workflow.drainage.initial_result_case);
fprintf(fid, 'cfg.initial_hydrology.timestep_file = ''%s'';\n', ...
    workflow.drainage.initial_timestep_file);
fprintf(fid, 'cfg.initial_ice_velocity.mode = ''file'';\n');
fprintf(fid, 'cfg.initial_ice_velocity.file = ''%s'';\n', ...
    workflow.inversion.hydrology_file);
fprintf(fid, 'cfg.initial_ice_velocity.variable = ''vv_hydro'';\n');
fprintf(fid, 'cfg.run.t_span_days = %s;\n', ...
    numeric_vector_literal(workflow.drainage.t_span_days));
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
    error('nevis_write_regional_drainage_case:MissingTemplate', ...
        'Drainage template case not found: %s', path);
end
end

function ensure_dir(path_in)
if exist(path_in, 'dir') ~= 7
    mkdir(path_in);
end
end

function txt = numeric_vector_literal(values)
values = values(:).';
if isempty(values)
    txt = '[]';
else
    pieces = {};
    start_idx = 1;
    while start_idx <= numel(values)
        if start_idx == numel(values)
            stop_idx = start_idx;
        else
            step = values(start_idx + 1) - values(start_idx);
            stop_idx = start_idx + 1;
            while stop_idx < numel(values) && ...
                    abs((values(stop_idx + 1) - values(stop_idx)) - step) < 1e-10
                stop_idx = stop_idx + 1;
            end
        end
        pieces{end+1} = vector_piece_literal(values(start_idx:stop_idx)); %#ok<AGROW>
        start_idx = stop_idx + 1;
    end
    if numel(pieces) == 1
        txt = pieces{1};
    else
        txt = ['[', strjoin(pieces, ', '), ']'];
    end
end
end

function txt = vector_piece_literal(values)
if numel(values) == 1
    txt = number_literal(values(1));
    return
end
step = values(2) - values(1);
if abs(step - 1) < 1e-10
    txt = sprintf('%s:%s', number_literal(values(1)), ...
        number_literal(values(end)));
else
    txt = sprintf('%s:%s:%s', number_literal(values(1)), ...
        number_literal(step), number_literal(values(end)));
end
end

function txt = number_literal(value)
value = round(value * 1e12) / 1e12;
txt = sprintf('%.12g', value);
end
