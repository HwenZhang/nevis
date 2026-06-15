function out = nevis_channel_timeseries_plot(casename, varargin)
%NEVIS_CHANNEL_TIMESERIES_PLOT Plot channel area and sheet thickness at one point.
%
% out = nevis_channel_timeseries_plot(casename)
% out = nevis_channel_timeseries_plot(casename, 'PointXY', [x y], 'SaveFigure', true)
% out = nevis_channel_timeseries_plot(casename, 'TimeRange', [10 20])
%
% By default, the point is chosen near the middle of the model domain:
% among nodes in a central search radius, the node with the largest channel
% cross section over the saved frames is used. Coordinates supplied through
% PointXY are nondimensional model coordinates, matching gg.nx and gg.ny.

if nargin < 1 || isempty(casename)
    error('nevis_channel_timeseries_plot:MissingCaseName', ...
        'Call nevis_channel_timeseries_plot with a casename.');
end

opts = local_defaults();
opts = local_parse_options(opts, varargin{:});

case_file = fullfile(opts.Root, 'results', casename, [casename '.mat']);
if exist(case_file, 'file') ~= 2
    error('nevis_channel_timeseries_plot:MissingCaseFile', ...
        'Could not find case file: %s', case_file);
end

data = load(case_file, 'pd', 'ps', 'gg', 'oo');
pd = data.pd;
ps = data.ps;
gg = data.gg;
oo = data.oo;

if isfield(oo, 'code') && exist(oo.code, 'dir') == 7
    addpath(oo.code);
elseif exist(fullfile(opts.Root, 'src'), 'dir') == 7
    addpath(fullfile(opts.Root, 'src'));
end

result_dir = fullfile(opts.Root, 'results', casename);
frames = opts.Frames;
if isempty(frames)
    frames = local_saved_frames(result_dir);
end
if isempty(frames)
    error('nevis_channel_timeseries_plot:MissingTimesteps', ...
        'No saved timestep files were found in %s.', result_dir);
end

requested_point_node = local_requested_point_node(opts, gg);
candidate_nodes = local_candidate_nodes(opts, gg, requested_point_node);

n_frames = numel(frames);
n_candidates = numel(candidate_nodes);
t = NaN(n_frames, 1);
S_candidates = NaN(n_frames, n_candidates);
hs_candidates = NaN(n_frames, n_candidates);

for k = 1:n_frames
    frame_file = fullfile(result_dir, [int2four(frames(k)) '.mat']);
    if exist(frame_file, 'file') ~= 2
        error('nevis_channel_timeseries_plot:MissingFrame', ...
            'Could not find saved timestep file: %s', frame_file);
    end

    dframe = load(frame_file, 'vv');
    vv = dframe.vv;

    t(k) = local_time_days(vv.t, ps, pd);
    S_node = ps.S * 0.25 * (gg.nmeanx * vv.Sx + gg.nmeany * vv.Sy + ...
        gg.nmeans * vv.Ss + gg.nmeanr * vv.Sr);
    hs_node = ps.hs * vv.hs(:);

    S_candidates(k, :) = S_node(candidate_nodes);
    hs_candidates(k, :) = hs_node(candidate_nodes);
end

if ~isempty(opts.TimeRange)
    keep = t >= opts.TimeRange(1) & t <= opts.TimeRange(2);
    if ~any(keep)
        error('nevis_channel_timeseries_plot:EmptyTimeRange', ...
            'No saved timestep falls in TimeRange [%g %g] d.', ...
            opts.TimeRange(1), opts.TimeRange(2));
    end
    t = t(keep);
    frames = frames(keep);
    S_candidates = S_candidates(keep, :);
    hs_candidates = hs_candidates(keep, :);
end

if isempty(requested_point_node)
    [~, point_idx] = max(max(S_candidates, [], 1));
    point_node = candidate_nodes(point_idx);
else
    point_idx = 1;
    point_node = requested_point_node;
end

S_point = S_candidates(:, point_idx);
hs_point = hs_candidates(:, point_idx);
point_xy = [gg.nx(point_node), gg.ny(point_node)];
point_xy_km = (ps.x / 1e3) * point_xy;

fig = figure('Color', 'w', 'Position', [100 100 900 550]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl);
plot(ax1, t, S_point, 'k-', 'LineWidth', 1.5);
grid(ax1, 'on');
xlabel(ax1, 't [d]');
ylabel(ax1, 'S [m^2]');
title(ax1, sprintf('Channel cross section, x = %.2f km, y = %.2f km', ...
    point_xy_km(1), point_xy_km(2)));

ax2 = nexttile(tl);
plot(ax2, t, hs_point, 'k-', 'LineWidth', 1.5);
grid(ax2, 'on');
xlabel(ax2, 't [d]');
ylabel(ax2, 'h_s [m]');
title(ax2, 'Cavity sheet thickness at the same near-channel point');
linkaxes([ax1 ax2], 'x');

if opts.SaveFigure
    if isempty(opts.FigurePath)
        figure_dir = fullfile(opts.Root, 'figures');
        if exist(figure_dir, 'dir') ~= 7
            mkdir(figure_dir);
        end
        opts.FigurePath = fullfile(figure_dir, [casename '_channel_sheet_timeseries.png']);
    end
    exportgraphics(fig, opts.FigurePath, 'Resolution', 300);
end

out = struct;
out.casename = casename;
out.figure = fig;
out.t = t;
out.frames = frames(:);
out.point_node = point_node;
out.point_xy = point_xy;
out.point_xy_km = point_xy_km;
out.S_point = S_point;
out.hs_point = hs_point;
out.candidate_nodes = candidate_nodes(:);
out.figure_path = opts.FigurePath;

end

function opts = local_defaults()
opts = struct;
opts.Root = '.';
opts.Frames = [];
opts.TimeRange = [];
opts.PointNode = [];
opts.PointXY = [];
opts.SearchRadius = [];
opts.SaveFigure = false;
opts.FigurePath = '';
end

function opts = local_parse_options(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('nevis_channel_timeseries_plot:BadOptions', ...
        'Options must be supplied as name-value pairs.');
end

for i = 1:2:numel(varargin)
    name = varargin{i};
    value = varargin{i + 1};
    if ~ischar(name) && ~isstring(name)
        error('nevis_channel_timeseries_plot:BadOptionName', ...
            'Option names must be character vectors or strings.');
    end
    name = lower(char(name));
    switch name
        case 'root'
            opts.Root = char(value);
        case 'frames'
            opts.Frames = value;
        case 'timerange'
            opts.TimeRange = local_validate_time_range(value);
        case 'pointnode'
            opts.PointNode = value;
        case 'pointxy'
            opts.PointXY = value;
        case 'searchradius'
            opts.SearchRadius = value;
        case 'savefigure'
            opts.SaveFigure = logical(value);
        case 'figurepath'
            opts.FigurePath = char(value);
        otherwise
            error('nevis_channel_timeseries_plot:UnknownOption', ...
                'Unknown option "%s".', name);
    end
end
end

function time_range = local_validate_time_range(value)
time_range = value;
if ~isnumeric(time_range) || numel(time_range) ~= 2 || ...
        any(~isfinite(time_range)) || time_range(1) > time_range(2)
    error('nevis_channel_timeseries_plot:BadTimeRange', ...
        'TimeRange must be a finite numeric vector [tmin tmax] in days.');
end
time_range = time_range(:).';
end

function frames = local_saved_frames(result_dir)
files = dir(fullfile(result_dir, '*.mat'));
frames = [];
for i = 1:numel(files)
    token = regexp(files(i).name, '^(\d{4})\.mat$', 'tokens', 'once');
    if ~isempty(token)
        frames(end + 1) = str2double(token{1}); %#ok<AGROW>
    end
end
frames = sort(frames);
end

function candidate_nodes = local_candidate_nodes(opts, gg, requested_point_node)
if ~isempty(requested_point_node)
    candidate_nodes = requested_point_node;
    return
end

valid = true(gg.nIJ, 1);
valid(gg.nout) = false;

nx = gg.nx(:);
ny = gg.ny(:);
xc = 0.5 * (min(nx(valid)) + max(nx(valid)));
yc = 0.5 * (min(ny(valid)) + max(ny(valid)));

dx = nx - xc;
dy = ny - yc;
dist = (dx.^2 + dy.^2).^0.5;

if isempty(opts.SearchRadius)
    width = max(nx(valid)) - min(nx(valid));
    height = max(ny(valid)) - min(ny(valid));
    search_radius = 0.15 * (width^2 + height^2)^0.5;
else
    search_radius = opts.SearchRadius;
end

candidate_nodes = find(valid & dist <= search_radius);
if isempty(candidate_nodes)
    dist(~valid) = inf;
    [~, idx] = min(dist);
    candidate_nodes = idx;
end
end

function point_node = local_requested_point_node(opts, gg)
valid = true(gg.nIJ, 1);
valid(gg.nout) = false;

if ~isempty(opts.PointNode)
    point_node = opts.PointNode;
    if ~isscalar(point_node) || point_node ~= round(point_node) || ...
            point_node < 1 || point_node > gg.nIJ || ~valid(point_node)
        error('nevis_channel_timeseries_plot:BadPointNode', ...
            'PointNode must be one in-domain node index.');
    end
    return
end

if ~isempty(opts.PointXY)
    if numel(opts.PointXY) ~= 2
        error('nevis_channel_timeseries_plot:BadPointXY', ...
            'PointXY must be [x y] in model coordinates.');
    end
    dx = gg.nx(:) - opts.PointXY(1);
    dy = gg.ny(:) - opts.PointXY(2);
    dist2 = dx.^2 + dy.^2;
    dist2(~valid) = inf;
    [~, point_node] = min(dist2);
    return
end

point_node = [];
end

function t_days = local_time_days(t_model, ps, pd)
if isfield(pd, 'td') && ~isempty(pd.td)
    t_days = t_model * ps.t / pd.td;
else
    t_days = t_model * ps.t / (24 * 60 * 60);
end
end
