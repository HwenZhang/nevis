function [tt, vv, info, state] = n2d_idealized_template(casename)
%N2D_IDEALIZED_TEMPLATE Run the standardized NEVIS idealized workflow.
%
% All case-specific analytical setup and parameters must live in
% ./cases/<casename>.m. This runner contains only the fixed workflow.

if nargin < 1 || isempty(casename)
    error('n2d_idealized_template:MissingCaseName', ...
        'Call n2d_idealized_template with a casename.');
end

format compact

root = './';
code = './src';
results = 'results';
addpath(code);

config_file = fullfile(root, 'cases', [casename, '.m']);
cfg = nevis_load_idealized_config(config_file, root);
if ~strcmp(cfg.casename, casename)
    error('n2d_idealized_template:CaseNameMismatch', ...
        'Config casename "%s" does not match requested case "%s".', ...
        cfg.casename, casename);
end

oo = struct;
oo.root = root;
oo.code = code;
oo.results = results;
oo.casename = cfg.casename;
oo.fn = ['/', oo.casename];
oo.rn = [oo.root, oo.results, oo.fn];
if exist(oo.rn, 'dir') ~= 7
    mkdir(oo.rn);
end
savefile = fullfile(oo.rn, oo.casename);

[pd, oo] = nevis_defaults([], oo);
[pd, oo] = nevis_apply_config_parameters(cfg, pd, oo);

ps = struct;
[ps, pp] = nevis_nondimension(pd, ps, oo);
[~, ~, pp] = nevis_apply_config_parameters(cfg, pd, oo, pp);

[pp, pd, ps, gg, aa, vv, oo, idealized] = ...
    nevis_setup_idealized_case(cfg, pd, ps, pp, oo);

save(savefile, 'pp', 'pd', 'ps', 'gg', 'aa', 'vv', 'oo', 'idealized');
[tt, vv, info] = nevis_timesteps(oo.t_span, vv, aa, pp, gg, oo);

aa = nevis_inputs(vv.t, aa, vv, pp, gg, oo);
oo.evaluate_variables = cfg.oo.evaluate_variables;
vv2 = nevis_backbone(inf, vv, vv, aa, pp, gg, oo);
vv2 = nevis_nodedischarge(vv2, aa, pp, gg, oo);

save(savefile, 'pp', 'pd', 'ps', 'gg', 'aa', 'vv', 'vv2', 'oo', ...
    'tt', 'info', 'idealized');

state = struct;
state.cfg = cfg;
state.pd = pd;
state.ps = ps;
state.pp = pp;
state.gg = gg;
state.aa = aa;
state.vv2 = vv2;
state.oo = oo;
state.idealized = idealized;
state.savefile = savefile;

end
