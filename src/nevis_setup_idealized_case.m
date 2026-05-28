function [pp, pd, ps, gg, aa, vv, oo, idealized] = nevis_setup_idealized_case(cfg, pd, ps, pp, oo)
%NEVIS_SETUP_IDEALIZED_CASE Build an analytical idealized NEVIS case.

if ~isfield(cfg, 'casename') || isempty(cfg.casename)
    error('nevis_setup_idealized_case:MissingCaseName', ...
        'cfg.casename is required.');
end
if ~strcmp(cfg.casename, oo.casename)
    error('nevis_setup_idealized_case:CaseNameMismatch', ...
        'Config casename "%s" does not match script casename "%s".', ...
        cfg.casename, oo.casename);
end

casename = oo.casename;
rn = oo.rn;
fn = oo.fn;
oo.dataset = cfg.casename;
oo.case_root = cfg.case_root;
oo.dn = cfg.case_root;

gcfg = cfg.geometry;
L = gcfg.length_m;
W = get_width(gcfg, L);
x = linspace(0, L / ps.x, gcfg.nx);
y = linspace(0, W / ps.x, gcfg.ny);
gg = nevis_grid(x, y, oo);

[b, s] = idealized_geometry(gcfg, L, W, pd, ps, pp, gg, oo, cfg);
gg.bed = b;
gg.surface = s;

H = max(s - b, 0);
Hmin = gcfg.min_thickness_m / ps.z;
nout = find(H <= Hmin);
gg = nevis_mask(gg, nout);
gg = nevis_mask_blister(gg, nout, oo);
gg.n1m = gg.n1;
gg.n1_blister = gg.n1m;
gg = nevis_label(gg, gg.n1m);
gg = nevis_label_blister(gg, gg.n1_blister, oo);

[aa, vv] = nevis_initialize(b, s, gg, pp, oo);

[vv, aa, pd] = idealized_initial_hydrology(cfg.initial_hydrology, ...
    aa, vv, pd, ps, pp, gg, oo, cfg);
aa.phi_b = max(aa.phi_0, aa.phi_a);

mcfg = cfg.moulins;
[mx_m, my_m] = idealized_moulin_coordinates(mcfg, L, W, pd, ps, pp, gg, oo, cfg);
mx = mx_m / ps.x;
my = my_m / ps.x;
[pp.ni_m, pp.sum_m] = nevis_moulins(mx, my, gg, oo);

lcfg = cfg.lakes;
pp.x_l = lcfg.x_fraction * L / ps.x;
pp.y_l = lcfg.y_fraction * W / ps.x;
pp.V_l = lcfg.volume_m3 / (ps.Q0 * ps.t);
pp.t_drainage = lcfg.drainage_day * pd.td / ps.t;
if isfield(lcfg, 'drainage_after_start_days') && ~isempty(lcfg.drainage_after_start_days)
    pp.t_drainage = vv.t + lcfg.drainage_after_start_days * pd.td / ps.t;
end
pp.t_duration = lcfg.duration_days * pd.td / ps.t;
[pp.ni_l, pp.sum_l] = nevis_lakes(pp.x_l, pp.y_l, gg, oo);

pp = idealized_moulin_input(cfg, pp, pd, ps, gg, oo);
pp = idealized_station_points(cfg, pp, ps, gg, oo, L);

oo.dt = cfg.run.dt_days * pd.td / ps.t;
if isfield(cfg.run, 't_span_after_start_days') && ~isempty(cfg.run.t_span_after_start_days)
    oo.t_span = vv.t + cfg.run.t_span_after_start_days * pd.td / ps.t;
else
    oo.t_span = cfg.run.t_span_days * pd.td / ps.t;
end
ni_gps = get_opt(pp, 'ni_gps', []);
oo.pts_ni = unique([pp.ni_l(:); pp.ni_m(:); ni_gps(:)]', 'stable');

oo.casename = casename;
oo.dataset = cfg.casename;
oo.case_root = cfg.case_root;
oo.dn = cfg.case_root;
oo.rn = rn;
oo.fn = fn;

idealized = struct;
idealized.length_m = L;
idealized.width_m = W;
idealized.x = x;
idealized.y = y;
idealized.nout = nout;
idealized.moulin_xy_m = [mx_m(:), my_m(:)];
idealized.lake_xy = [pp.x_l, pp.y_l];

end

function W = get_width(gcfg, L)
if isfield(gcfg, 'width_m') && ~isempty(gcfg.width_m)
    W = gcfg.width_m;
else
    W = gcfg.width_ratio * L;
end
end

function [b, s] = idealized_geometry(gcfg, L, W, pd, ps, pp, gg, oo, cfg)
mode = lower(gcfg.mode);
X_m = ps.x * gg.nx;
Y_m = ps.x * gg.ny;
switch mode
    case 'hewitt_2013_sqrt'
        b = (gcfg.bed_elevation_m / ps.z) * ones(gg.nIJ, 1);
        s = (gcfg.surface_head_m / ps.z) * ...
            (1 - reshape(X_m, gg.nIJ, 1) / L) .^ gcfg.surface_exponent;
    case 'function'
        b_m = call_config_function(gcfg.bed_function, ...
            {X_m, Y_m, L, W, pd, ps, pp, gg, oo, cfg});
        s_m = call_config_function(gcfg.surface_function, ...
            {X_m, Y_m, L, W, pd, ps, pp, gg, oo, cfg});
        b = as_node_vector(b_m, gg, 'cfg.geometry.bed_function') / ps.z;
        s = as_node_vector(s_m, gg, 'cfg.geometry.surface_function') / ps.z;
    otherwise
        error('nevis_setup_idealized_case:UnsupportedGeometry', ...
            'Unsupported cfg.geometry.mode: %s', gcfg.mode);
end
end

function [vv, aa, pd] = idealized_initial_hydrology(icfg, aa, vv, pd, ps, pp, gg, oo, cfg)
mode = lower(icfg.mode);
switch mode
    case 'k_factor_cavity_equilibrium'
        pd.k_f = icfg.k_factor;
        vv.phi = aa.phi_a + pd.k_f * (aa.phi_0 - aa.phi_a);
        N = aa.phi_0 - vv.phi;
        vv.hs = ((((pd.u_b * pd.h_r / pd.l_r) ./ ...
            ((pd.u_b / pd.l_r) + (pd.K_c .* ((ps.phi * N) .^ 3))))) ./ ps.h);
    case 'function'
        [vv, aa] = icfg.function(aa, vv, pd, ps, pp, gg, oo, cfg);
    case 'result_timestep'
        init_file = fullfile(oo.root, oo.results, icfg.result_case, ...
            icfg.timestep_file);
        if exist(init_file, 'file') ~= 2
            error('nevis_setup_idealized_case:MissingInitialTimestep', ...
                'Initial timestep file not found: %s', init_file);
        end
        data = load(init_file, 'vv');
        if ~isfield(data, 'vv')
            error('nevis_setup_idealized_case:MissingInitialVV', ...
                'Initial timestep file must contain variable vv: %s', init_file);
        end
        vv = data.vv;
    otherwise
        error('nevis_setup_idealized_case:UnsupportedInitialHydrology', ...
            'Unsupported cfg.initial_hydrology.mode: %s', icfg.mode);
end
end

function [x_m, y_m] = idealized_moulin_coordinates(mcfg, L, W, pd, ps, pp, gg, oo, cfg)
switch lower(mcfg.mode)
    case 'coordinates_fractional'
        x_m = mcfg.x_fraction(:) * L;
        y_m = mcfg.y_fraction(:) * W;
    case 'coordinates_m'
        x_m = mcfg.x_m(:);
        y_m = mcfg.y_m(:);
    case 'function'
        [x_m, y_m] = call_position_function(mcfg.position_function, ...
            {L, W, pd, ps, pp, gg, oo, cfg});
        x_m = x_m(:);
        y_m = y_m(:);
    otherwise
        error('nevis_setup_idealized_case:UnsupportedMoulins', ...
            'Unsupported cfg.moulins.mode: %s', mcfg.mode);
end
if numel(x_m) ~= numel(y_m)
    error('nevis_setup_idealized_case:MoulinCoordinateSizeMismatch', ...
        'Moulin x and y coordinate arrays must have the same number of elements.');
end
end

function [x_m, y_m] = call_position_function(fn, args)
n = nargin(fn);
if n < 0 || n >= numel(args)
    [x_m, y_m] = fn(args{:});
else
    [x_m, y_m] = fn(args{1:n});
end
end

function pp = idealized_moulin_input(cfg, pp, pd, ps, gg, oo)
input_cfg = get_input_cfg(cfg);
pp.meltE = @(t) (input_cfg.runoff_max_mm_per_day / 1000 / pd.td / ps.m) .* ...
    (1 - exp(-t / (input_cfg.ramp_days * pd.td / ps.t)));

if isfield(cfg.moulins, 'input_function') && ~isempty(cfg.moulins.input_function)
    fn = cfg.moulins.input_function;
    pp.input_function = @(t) normalize_moulin_input( ...
        call_config_function(fn, {t, pd, ps, pp, gg, oo, cfg}), pp) ./ ...
        (ps.m * ps.x ^ 2);
else
    q = get_opt(cfg.moulins, 'input_m3_per_s', input_cfg.moulin_input_m3_per_s);
    ramp_days = get_opt(cfg.moulins, 'input_ramp_days', input_cfg.ramp_days);
    pp.input_function = @(t) normalize_moulin_input(q(:) .* ...
        (1 - exp(-t / (ramp_days * pd.td / ps.t))), pp) ./ ...
        (ps.m * ps.x ^ 2);
end
end

function pp = idealized_station_points(cfg, pp, ps, gg, oo, L)
if ~isfield(cfg, 'stations') || isempty(cfg.stations)
    pp.ni_gps = [];
    return
end
switch lower(cfg.stations.mode)
    case 'downstream_lake'
        if isempty(pp.x_l)
            pp.ni_gps = [];
            return
        end
        x = pp.x_l(1):cfg.stations.spacing_m / ps.x:L / ps.x;
        y = pp.y_l(1) * ones(size(x));
        pp.ni_gps = nevis_gps_array(x, y, gg, oo);
    otherwise
        error('nevis_setup_idealized_case:UnsupportedStations', ...
            'Unsupported cfg.stations.mode: %s', cfg.stations.mode);
end
end

function input_cfg = get_input_cfg(cfg)
if isfield(cfg, 'input') && ~isempty(cfg.input)
    input_cfg = cfg.input;
else
    input_cfg = struct;
end
input_cfg.runoff_max_mm_per_day = get_opt(input_cfg, 'runoff_max_mm_per_day', 0);
input_cfg.moulin_input_m3_per_s = get_opt(input_cfg, 'moulin_input_m3_per_s', 0);
input_cfg.ramp_days = get_opt(input_cfg, 'ramp_days', 30);
end

function q = normalize_moulin_input(q, pp)
q = q(:);
nm = numel(pp.ni_m);
if numel(q) == 1 && nm > 1
    q = repmat(q, nm, 1);
end
if numel(q) ~= nm
    error('nevis_setup_idealized_case:MoulinInputSizeMismatch', ...
        'Moulin input must be scalar or have one value per active moulin.');
end
end

function out = call_config_function(fn, args)
n = nargin(fn);
if n < 0 || n >= numel(args)
    out = fn(args{:});
else
    out = fn(args{1:n});
end
end

function values = as_node_vector(values, gg, label)
if isscalar(values)
    values = values * ones(gg.nIJ, 1);
elseif isequal(size(values), size(gg.nx))
    values = reshape(values, gg.nIJ, 1);
else
    values = values(:);
end
if numel(values) ~= gg.nIJ
    error('nevis_setup_idealized_case:AnalyticalFieldSizeMismatch', ...
        '%s must return a scalar, a grid-sized array, or one value per node.', ...
        label);
end
end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
