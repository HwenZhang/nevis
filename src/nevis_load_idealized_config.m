function cfg = nevis_load_idealized_config(config_file, root)
%NEVIS_LOAD_IDEALIZED_CONFIG Load a strict analytical idealized case config.

if exist(config_file, 'file') ~= 2
    error('nevis_load_idealized_config:MissingConfig', ...
        'Idealized config not found: %s', config_file);
end

[~, ~, ext] = fileparts(config_file);
switch lower(ext)
    case '.m'
        cfg = [];
        run(config_file);
        if isempty(cfg)
            error('nevis_load_idealized_config:MissingCfg', ...
                'M-file config must assign a variable named cfg.');
        end
    case '.mat'
        data = load(config_file);
        if ~isfield(data, 'cfg')
            error('nevis_load_idealized_config:MissingCfg', ...
                'MAT config must contain a variable named cfg.');
        end
        cfg = data.cfg;
    case '.json'
        cfg = jsondecode(fileread(config_file));
    otherwise
        error('nevis_load_idealized_config:UnsupportedConfig', ...
            'Unsupported idealized config extension: %s', ext);
end

if nargin >= 2 && ~isempty(root) && ~isfield(cfg, 'root')
    cfg.root = root;
end

cfg.config_file = config_file;
cfg = nevis_validate_idealized_config(cfg);

end
