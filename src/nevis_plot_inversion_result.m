function summary = nevis_plot_inversion_result(inversion_file, output_dir)
%NEVIS_PLOT_INVERSION_RESULT Plot saved NEVIS inversion diagnostics.
%
% summary = nevis_plot_inversion_result(inversion_file)
% summary = nevis_plot_inversion_result(inversion_file, output_dir)

if nargin < 1 || isempty(inversion_file)
    error('nevis_plot_inversion_result:MissingInput', ...
        'Call with a saved C_inversion_results.mat file.');
end
if exist(inversion_file, 'file') ~= 2
    error('nevis_plot_inversion_result:MissingFile', ...
        'Inversion result not found: %s', inversion_file);
end
if nargin < 2 || isempty(output_dir)
    output_dir = fileparts(inversion_file);
end
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

data = load(inversion_file);
required = {'C_hat', 'C_hat_dim', 'u_obs_noisy', 'v_obs_noisy', ...
    'N_current', 'aa', 'pp', 'gg', 'oo', 'ps'};
for i = 1:numel(required)
    if ~isfield(data, required{i})
        error('nevis_plot_inversion_result:MissingVariable', ...
            'Inversion file must contain %s.', required{i});
    end
end

C_hat = data.C_hat;
C_hat_dim = data.C_hat_dim;
u_obs = data.u_obs_noisy;
v_obs = data.v_obs_noisy;
N_current = data.N_current;
aa = data.aa;
pp = data.pp;
gg = data.gg;
oo = data.oo;
ps = data.ps;
pd = load_source_pd(data, inversion_file);

aa_inv = aa;
aa_inv.C = C_hat;
u_inv = u_obs;
v_inv = v_obs;
u_inv(isnan(u_inv) & ~ismember((1:gg.eIJ)', gg.eout2)) = 0;
v_inv(isnan(v_inv) & ~ismember((1:gg.fIJ)', gg.fout2)) = 0;
[u_inv, v_inv] = nevis_velocity(aa.H, u_inv, v_inv, ...
    N_current, aa_inv, pp, gg, oo);

[speed_obs, speed_inv, speed_misfit] = node_speed_fields( ...
    u_obs, v_obs, u_inv, v_inv, gg, ps, pd);
rel_err = abs(speed_misfit) ./ max(speed_obs, 1);

nin = gg.nin2;
rmse_speed = sqrt(mean((speed_inv(nin) - speed_obs(nin)).^2, 'omitnan'));
too_low = find(speed_misfit < -20 & speed_obs > 50);
too_high = find(speed_misfit > 20 & speed_obs > 50);

summary = struct;
summary.inversion_file = inversion_file;
summary.output_dir = output_dir;
summary.rmse_speed_m_per_yr = rmse_speed;
summary.final_J = get_field(data, 'J_hat', NaN);
summary.exitflag = get_field(data, 'exitflag', NaN);
summary.too_slow_nodes = numel(too_low);
summary.too_fast_nodes = numel(too_high);
summary.diagnostics_png = fullfile(output_dir, 'inversion_diagnostics.png');
summary.convergence_png = fullfile(output_dir, 'inversion_convergence.png');

plot_diagnostics(summary.diagnostics_png, C_hat_dim, N_current, ...
    speed_obs, speed_inv, speed_misfit, rel_err, gg, ps, ...
    rmse_speed, too_low, too_high);

if isfield(data, 'history') && isstruct(data.history)
    plot_convergence(summary.convergence_png, data.history);
else
    summary.convergence_png = '';
end
end

function pd = load_source_pd(data, inversion_file)
if isfield(data, 'pd')
    pd = data.pd;
    return
end

if isfield(data, 'cfg') && isfield(data.cfg, 'source') && ...
        isfield(data.cfg.source, 'state_file')
    state_file = resolve_project_path(data.cfg.source.state_file, data.cfg);
    if exist(state_file, 'file') == 2
        state = load(state_file, 'pd');
        if isfield(state, 'pd')
            pd = state.pd;
            return
        end
    end
end

warning('nevis_plot_inversion_result:MissingPd', ...
    ['Could not load pd from source state for %s. Using 365-day year ' ...
     'only for velocity unit conversion.'], inversion_file);
pd = struct('ty', 365 * 24 * 3600);
end

function [speed_obs, speed_inv, speed_misfit] = node_speed_fields( ...
    u_obs, v_obs, u_inv, v_inv, gg, ps, pd)
Un_obs = sqrt((gg.nmeanx2(:, gg.es2) * u_obs(gg.es2)).^2 + ...
    (gg.nmeany2(:, gg.fs2) * v_obs(gg.fs2)).^2);
Un_inv = sqrt((gg.nmeanx2(:, gg.es2) * u_inv(gg.es2)).^2 + ...
    (gg.nmeany2(:, gg.fs2) * v_inv(gg.fs2)).^2);
speed_obs = Un_obs * ps.u_b * pd.ty;
speed_inv = Un_inv * ps.u_b * pd.ty;
speed_misfit = speed_inv - speed_obs;
speed_obs(gg.nout) = NaN;
speed_inv(gg.nout) = NaN;
speed_misfit(gg.nout) = NaN;
end

function plot_diagnostics(out_png, C_hat_dim, N_current, speed_obs, ...
    speed_inv, speed_misfit, rel_err, gg, ps, rmse_speed, too_low, too_high)
cmap_div = interp1([0 0.5 1], [0 0 1; 1 1 1; 1 0 0], ...
    linspace(0, 1, 256));

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [80 80 1600 900]);
tlo = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

plot_field(nexttile(tlo), gg, log10(C_hat_dim), ...
    'log10(C dim)', 'parula', []);
plot_field(nexttile(tlo), gg, N_current * ps.phi / 1e6, ...
    'N [MPa]', 'parula', []);
plot_field(nexttile(tlo), gg, speed_obs, ...
    'Observed speed [m/yr]', 'parula', [0 300]);
plot_field(nexttile(tlo), gg, speed_inv, ...
    'Inverted speed [m/yr]', 'parula', [0 300]);
plot_field(nexttile(tlo), gg, speed_misfit, ...
    'Speed misfit [m/yr]', cmap_div, [-100 100]);
plot_field(nexttile(tlo), gg, 100 * rel_err, ...
    'Relative speed error [%]', 'parula', [0 100]);

title(tlo, sprintf(['Inversion diagnostics: RMSE %.2f m/yr, ' ...
    'too slow %d, too fast %d'], rmse_speed, numel(too_low), ...
    numel(too_high)));

exportgraphics(fig, out_png, 'Resolution', 200);
close(fig);
end

function plot_field(ax, gg, values, plot_title, cmap, climits)
z = reshape(values, gg.nI, gg.nJ);
z(gg.nout) = NaN;
pcolor(ax, gg.nx, gg.ny, z);
shading(ax, 'flat');
axis(ax, 'equal');
axis(ax, 'tight');
title(ax, plot_title);
xlabel(ax, 'x');
ylabel(ax, 'y');
colorbar(ax);
if ischar(cmap)
    colormap(ax, cmap);
else
    colormap(ax, cmap);
end
if ~isempty(climits)
    clim(ax, climits);
end
end

function plot_convergence(out_png, history)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1000 420]);
tlo = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax = nexttile(tlo);
if isfield(history, 'iter') && isfield(history, 'J')
    semilogy(ax, history.iter, history.J, 'o-', 'LineWidth', 1.2);
else
    semilogy(ax, history.J, 'o-', 'LineWidth', 1.2);
end
xlabel(ax, 'Iteration');
ylabel(ax, 'J');
title(ax, 'Objective');
grid(ax, 'on');

ax = nexttile(tlo);
if isfield(history, 'normg')
    if isfield(history, 'iter')
        semilogy(ax, history.iter, history.normg, 'o-', 'LineWidth', 1.2);
    else
        semilogy(ax, history.normg, 'o-', 'LineWidth', 1.2);
    end
    ylabel(ax, 'First-order optimality');
else
    axis(ax, 'off');
end
xlabel(ax, 'Iteration');
title(ax, 'Gradient/optimality');
grid(ax, 'on');

exportgraphics(fig, out_png, 'Resolution', 200);
close(fig);
end

function value = get_field(s, name, default)
if isfield(s, name)
    value = s.(name);
else
    value = default;
end
end

function path_out = resolve_project_path(path_in, cfg)
if startsWith(path_in, filesep) || ...
        ~isempty(regexp(path_in, '^[A-Za-z]:[\\/]', 'once'))
    path_out = path_in;
else
    path_out = fullfile(cfg.root, path_in);
end
end
