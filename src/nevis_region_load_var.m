function value = nevis_region_load_var(filename, variable)
%NEVIS_REGION_LOAD_VAR Load a named variable, or the only variable in a MAT-file.

data = load(filename);

if nargin >= 2 && ~isempty(variable)
    if ~isfield(data, variable)
        error('nevis_region_load_var:MissingVariable', ...
            'Variable "%s" not found in "%s".', variable, filename);
    end
    value = data.(variable);
    return
end

names = fieldnames(data);
if numel(names) ~= 1
    error('nevis_region_load_var:AmbiguousFile', ...
        'File "%s" has %d variables; specify variable in config.', filename, numel(names));
end

value = data.(names{1});

end
