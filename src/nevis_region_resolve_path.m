function path_out = nevis_region_resolve_path(path_in, cfg)
%NEVIS_REGION_RESOLVE_PATH Resolve and isolate regional dataset paths.

if nargin < 2
    cfg = struct;
end

if isempty(path_in)
    path_out = '';
    return
end

root = dataset_root(cfg);
if isempty(root)
    error('nevis_region_resolve_path:MissingDataRoot', ...
        'cfg.dataset.root, cfg.dataset_root, or cfg.case_root is required for strict data isolation.');
end

if is_absolute_path(path_in)
    path_out = path_in;
else
    path_out = fullfile(root, path_in);
end

root_canon = canonical_path(root);
path_canon = canonical_path(path_out);
if ~startsWith([path_canon filesep], [root_canon filesep])
    error('nevis_region_resolve_path:PathOutsideDataRoot', ...
        'Configured file is outside dataset root. file="%s", dataset_root="%s"', ...
        path_out, root);
end

end

function root = dataset_root(cfg)
root = '';
if isfield(cfg, 'dataset') && isstruct(cfg.dataset) && ...
        isfield(cfg.dataset, 'root') && ~isempty(cfg.dataset.root)
    root = cfg.dataset.root;
elseif isfield(cfg, 'dataset_root') && ~isempty(cfg.dataset_root)
    root = cfg.dataset_root;
elseif isfield(cfg, 'case_root') && ~isempty(cfg.case_root)
    root = cfg.case_root;
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
