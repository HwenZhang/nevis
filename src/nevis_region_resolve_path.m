function path_out = nevis_region_resolve_path(path_in, cfg)
%NEVIS_REGION_RESOLVE_PATH Resolve and isolate paths under cfg.data_root.

if nargin < 2
    cfg = struct;
end

if isempty(path_in)
    path_out = '';
    return
end

if ~isfield(cfg, 'data_root') || isempty(cfg.data_root)
    error('nevis_region_resolve_path:MissingDataRoot', ...
        'cfg.data_root is required for strict data isolation.');
end

if is_absolute_path(path_in)
    path_out = path_in;
else
    path_out = fullfile(cfg.data_root, path_in);
end

root_canon = canonical_path(cfg.data_root);
path_canon = canonical_path(path_out);
if ~startsWith([path_canon filesep], [root_canon filesep])
    error('nevis_region_resolve_path:PathOutsideDataRoot', ...
        'Configured file is outside data_root. file="%s", data_root="%s"', ...
        path_out, cfg.data_root);
end

end

function tf = is_absolute_path(p)
tf = startsWith(p, filesep) || ...
    (~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once')));
end

function p = canonical_path(path_in)
if exist(path_in, 'file') || exist(path_in, 'dir')
    p = char(java.io.File(path_in).getCanonicalPath());
else
    [parent, name, ext] = fileparts(path_in);
    if isempty(parent)
        parent = '.';
    end
    p = fullfile(char(java.io.File(parent).getCanonicalPath()), [name ext]);
end
end
