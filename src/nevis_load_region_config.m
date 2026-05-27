function cfg = nevis_load_region_config(config_file, root)
%NEVIS_LOAD_REGION_CONFIG Load a strict per-case region config.
%
% The config is mandatory. Missing config fields and missing data files are
% errors; this function never falls back to legacy/default data.

if exist(config_file, 'file') ~= 2
    error('nevis_load_region_config:MissingConfig', ...
        'Region config not found: %s', config_file);
end

[~, ~, ext] = fileparts(config_file);
switch lower(ext)
    case '.mat'
        data = load(config_file);
        if isfield(data, 'cfg')
            cfg = data.cfg;
        else
            error('nevis_load_region_config:MissingCfg', ...
                'MAT config must contain a variable named cfg.');
        end
    case '.json'
        txt = fileread(config_file);
        cfg = jsondecode(txt);
    case '.m'
        cfg = [];
        run(config_file);
        if isempty(cfg)
            error('nevis_load_region_config:MissingCfg', ...
                'M-file config must assign a variable named cfg.');
        end
    otherwise
        error('nevis_load_region_config:UnsupportedConfig', ...
            'Unsupported region config extension: %s', ext);
end

if nargin >= 2 && ~isempty(root) && ~isfield(cfg, 'root')
    cfg.root = root;
end

cfg.config_file = config_file;
cfg = nevis_validate_region_config(cfg);

end
