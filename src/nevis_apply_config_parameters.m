function [pd, oo, pp] = nevis_apply_config_parameters(cfg, pd, oo, pp)
%NEVIS_APPLY_CONFIG_PARAMETERS Apply case-bound parameters from cfg.
%
% cfg.oo, cfg.pd, and cfg.pp are the only supported case-specific parameter
% entry points for the standardized regional workflow. Function-handle values
% are evaluated against the current target struct, which allows derived
% physical parameters such as pd.melt to live in the config file.

if nargin < 4
    pp = struct;
end

if isfield(cfg, 'oo')
    oo = apply_fields(oo, cfg.oo);
end
if isfield(cfg, 'pd')
    pd = apply_fields(pd, cfg.pd);
end
if isfield(cfg, 'pp')
    pp = apply_fields(pp, cfg.pp);
end

end

function target = apply_fields(target, values)
names = fieldnames(values);
for i = 1:numel(names)
    name = names{i};
    value = values.(name);
    if isa(value, 'function_handle')
        target.(name) = value(target);
    else
        target.(name) = value;
    end
end
end
