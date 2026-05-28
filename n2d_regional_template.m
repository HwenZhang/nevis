function [tt, vv, info, state] = n2d_regional_template(casename)
%N2D_REGIONAL_TEMPLATE Run the standardized NEVIS regional workflow.
%
% All case-specific parameters must live in ./cases/<casename>.m. Heavy data
% are referenced through the config's shared dataset manifest.

if nargin < 1 || isempty(casename)
    error('n2d_regional_template:MissingCaseName', ...
        'Call n2d_regional_template with a casename.');
end

format compact

% Set up the root directory
root = './';
code = './src';
results = 'results';
addpath(code);

% Read in the config file from ./cases/<casename>.m.
% Note: the user should guarantee the integrity of the data files, including the config
config_file = fullfile(root, 'cases', [casename, '.m']);
cfg = nevis_load_region_config(config_file, root);
if ~strcmp(cfg.casename, casename)
    error('n2d_regional_template:CaseNameMismatch', ...
        'Config casename "%s" does not match requested case "%s".', ...
        cfg.casename, casename);
end

% Set up the option struct
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

% import parameters from the config
[pd, oo] = nevis_defaults([], oo);
[pd, oo] = nevis_apply_config_parameters(cfg, pd, oo);
ps = struct;
[ps, pp] = nevis_nondimension(pd, ps, oo);
[~, ~, pp] = nevis_apply_config_parameters(cfg, pd, oo, pp);

% set up the regional case
[pp, pd, ps, gg, aa, vv, oo, region] = nevis_setup_region_case(cfg, pd, ps, pp, oo);
save(savefile, 'pp', 'pd', 'ps', 'gg', 'aa', 'vv', 'oo', 'region');

% simulation run
[tt, vv, info] = nevis_timesteps(oo.t_span, vv, aa, pp, gg, oo);

% postprocessing
aa = nevis_inputs(vv.t, aa, vv, pp, gg, oo);
oo.evaluate_variables = cfg.oo.evaluate_variables;
vv2 = nevis_backbone(inf, vv, vv, aa, pp, gg, oo);
vv2 = nevis_nodedischarge(vv2, aa, pp, gg, oo);

save(savefile, 'pp', 'pd', 'ps', 'gg', 'aa', 'vv', 'vv2', 'oo', ...
    'tt', 'info', 'region');

state = struct;
state.cfg = cfg;
state.pd = pd;
state.ps = ps;
state.pp = pp;
state.gg = gg;
state.aa = aa;
state.vv2 = vv2;
state.oo = oo;
state.region = region;
state.savefile = savefile;

end
