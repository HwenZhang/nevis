function [inv, vv_hydro, summary] = nevis_run_inversion(config_file)
%NEVIS_RUN_INVERSION Dispatch an inversion config to the selected sliding law.
%
% Usage:
%   [inv, vv_hydro, summary] = nevis_run_inversion('./cases/nevis_inversion.m');
%
% The config must define cfg.sliding_law. This keeps the user-facing
% inversion workflow fixed while allowing the sliding law to be case-specific.

if nargin < 1 || isempty(config_file)
    error('nevis_run_inversion:MissingConfig', ...
        'Call nevis_run_inversion with an inversion config file.');
end

cfg = load_dispatch_config(config_file);
if ~isfield(cfg, 'sliding_law') || isempty(cfg.sliding_law)
    error('nevis_run_inversion:MissingSlidingLaw', ...
        'Inversion config must define cfg.sliding_law.');
end

sliding_law = lower(strrep(cfg.sliding_law, '-', '_'));
switch sliding_law
    case {'weertman', 'regularized_weertman', 'regularised_weertman'}
        [inv, vv_hydro, summary] = nevis_inv_C(config_file);
    case {'regularized_coulomb', 'regularised_coulomb', 'coulomb'}
        [inv, vv_hydro, summary] = nevis_inv_C1(config_file);
    otherwise
        error('nevis_run_inversion:UnsupportedSlidingLaw', ...
            ['Unsupported cfg.sliding_law "%s". Use "weertman" or ' ...
             '"regularized_coulomb".'], cfg.sliding_law);
end
summary.sliding_law = sliding_law;
end

function cfg = load_dispatch_config(config_file)
    if exist(config_file, 'file') ~= 2
        error('nevis_run_inversion:MissingConfig', ...
            'Inversion config not found: %s', config_file);
    end

    [~, ~, ext] = fileparts(config_file);
    switch lower(ext)
        case '.m'
            run(config_file);
            if exist('cfg', 'var') ~= 1
                error('nevis_run_inversion:MissingCfg', ...
                    'M-file inversion config must assign variable cfg.');
            end
        case '.mat'
            data = load(config_file);
            if ~isfield(data, 'cfg')
                error('nevis_run_inversion:MissingCfg', ...
                    'MAT inversion config must contain variable cfg.');
            end
            cfg = data.cfg;
        case '.json'
            cfg = jsondecode(fileread(config_file));
        otherwise
            error('nevis_run_inversion:UnsupportedConfig', ...
                'Unsupported inversion config extension: %s', ext);
    end
end
