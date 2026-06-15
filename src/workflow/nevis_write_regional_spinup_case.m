function case_file = nevis_write_regional_spinup_case(workflow, case_dir, template_case)
%NEVIS_WRITE_REGIONAL_SPINUP_CASE Materialize a per-case spinup config.
%
% Writes case_dir/<workflow.spinup.case>.m containing the template content
% inlined verbatim followed by this workflow's overrides. The template
% file path is recorded at the top so the lineage is auditable, but the
% generated file is fully self-contained: editing it changes what the
% spinup actually runs.
%
% Idempotent: if the target file already exists, returns the path
% unchanged so user edits survive subsequent workflow invocations. To
% regenerate from the current template, delete the file first or pass
% overwrite=true through the workflow.

case_file = fullfile(case_dir, [workflow.spinup.case, '.m']);
if exist(case_file, 'file') == 2 && ~workflow_force_overwrite(workflow)
    return
end

template_file = resolve_template_path(workflow, template_case);
template_text = fileread(template_file);
template_text = strip_template_run_block(template_text);

ensure_dir(case_dir);
fid = fopen(case_file, 'w');
if fid == -1
    error('nevis_write_regional_spinup_case:FileWrite', ...
        'Could not write spinup case config: %s', case_file);
end
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, '%% Per-case hydrology-only spinup config.\n');
fprintf(fid, '%% Generated from template: %s\n', template_file);
fprintf(fid, '%% Self-contained: edit freely. Subsequent workflow runs do not\n');
fprintf(fid, '%% overwrite this file. Delete it to regenerate from the template.\n\n');
fprintf(fid, '%% ====== inlined template (%s) ======\n\n', template_case);
fprintf(fid, '%s', template_text);
if ~endsWith(template_text, newline)
    fprintf(fid, '\n');
end
fprintf(fid, '\n%% ====== workflow overrides ======\n');
fprintf(fid, 'cfg.casename = ''%s'';\n', workflow.spinup.case);
fprintf(fid, 'cfg.dataset.name = ''%s'';\n', workflow.dataset.name);
fprintf(fid, 'cfg.dataset.root = ''%s'';\n', workflow.dataset.root);
fprintf(fid, 'cfg.dataset.manifest = ''%s'';\n', workflow.dataset.manifest);
fprintf(fid, 'cfg.moulins.file = ''moulins_2022_coordinates.mat'';\n');
fprintf(fid, 'cfg.moulins.variable = ''moulins'';\n');
fprintf(fid, 'cfg.moulins.mode = ''coordinates'';\n');
fprintf(fid, 'cfg.moulins.x_field = ''x_m'';\n');
fprintf(fid, 'cfg.moulins.y_field = ''y_m'';\n');
fprintf(fid, 'cfg.lakes.file = ''environs_lakes_2022B_250416.mat'';\n');
fprintf(fid, 'cfg.run.t_span_days = %s;\n', ...
    numeric_vector_literal(workflow.spinup.t_span_days));
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
    error('nevis_write_regional_spinup_case:MissingTemplate', ...
        'Spinup template case not found: %s', path);
end
end

function text = strip_template_run_block(text)
% Spinup template starts with `case_dir = ...` and `run(fullfile(case_dir, ...))`
% to chain into nevis_regional.m. Replace that block with an inlined copy
% of nevis_regional.m so the generated per-case file is fully self-contained.
lines = regexp(text, '\r?\n', 'split');
out = {};
i = 1;
nested_loaded = false;
while i <= numel(lines)
    line = lines{i};
    if ~nested_loaded && contains(line, 'case_dir = fileparts(mfilename')
        % Skip this line and the following run(fullfile(case_dir, '<file>.m')); line.
        if i + 1 <= numel(lines)
            run_line = strtrim(lines{i + 1});
            tok = regexp(run_line, "run\(fullfile\(case_dir, '([^']+)'\)\);?", ...
                'tokens', 'once');
            if ~isempty(tok)
                nested_path = fullfile('./cases', 'templates', tok{1});
                nested_text = fileread(nested_path);
                nested_text = strip_template_run_block(nested_text);
                out{end+1} = sprintf('%% ---- inlined %s ----', tok{1}); %#ok<AGROW>
                out{end+1} = nested_text; %#ok<AGROW>
                out{end+1} = sprintf('%% ---- end %s ----', tok{1}); %#ok<AGROW>
                i = i + 2;
                nested_loaded = true;
                continue
            end
        end
    end
    out{end+1} = line; %#ok<AGROW>
    i = i + 1;
end
text = strjoin(out, newline);
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
elseif numel(values) > 1 && all(abs(diff(values) - diff(values(1:2))) < eps)
    step = values(2) - values(1);
    if step == 1
        txt = sprintf('%g:%g', values(1), values(end));
    else
        txt = sprintf('%g:%g:%g', values(1), step, values(end));
    end
else
    txt = mat2str(values);
end
end
