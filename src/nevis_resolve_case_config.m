function config_file = nevis_resolve_case_config(casename, root)
%NEVIS_RESOLVE_CASE_CONFIG Find a case config file under cases/.
%
% Per-case subdir layouts take precedence over the legacy flat layout, so
% stale flat files left behind by older workflow versions do not shadow
% the per-case files written by nevis_prepare_demo_workflow.
%
% Lookup order:
%   1. <root>/cases/<casename>/<casename>.m   (per-case subdir, direct hit)
%   2. <root>/cases/*/<casename>.m            (sibling subdir; e.g. a
%      spinup config bundled with its forward case)
%   3. <root>/cases/<casename>.m              (legacy flat layout;
%      retained so the original templates still resolve)
%
% Errors if not found, or if step 2 matches more than one location.

if nargin < 2 || isempty(root)
    root = './';
end

self_subdir = fullfile(root, 'cases', casename, [casename, '.m']);
if exist(self_subdir, 'file') == 2
    config_file = self_subdir;
    return
end

candidates = dir(fullfile(root, 'cases', '*', [casename, '.m']));
candidates = candidates(~strcmp({candidates.name}, ''));
if isscalar(candidates)
    config_file = fullfile(candidates(1).folder, candidates(1).name);
    return
elseif numel(candidates) > 1
    locs = strjoin(arrayfun(@(c) fullfile(c.folder, c.name), candidates, ...
        'UniformOutput', false), ', ');
    error('nevis_resolve_case_config:Ambiguous', ...
        'Casename "%s" found in multiple cases/ subdirectories: %s', ...
        casename, locs);
end

flat = fullfile(root, 'cases', [casename, '.m']);
if exist(flat, 'file') == 2
    config_file = flat;
    return
end

error('nevis_resolve_case_config:NotFound', ...
    ['Case config not found for "%s". Looked in cases/%s/%s.m, ' ...
     'cases/*/%s.m, and cases/%s.m.'], ...
    casename, casename, casename, casename, casename);
end
