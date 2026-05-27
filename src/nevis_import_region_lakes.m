function [pp, lakes, region] = nevis_import_region_lakes(cfg, pp, pd, ps, gg, oo, region)
%NEVIS_IMPORT_REGION_LAKES Load lake catalogue/table and map lakes to nodes.

lcfg = cfg.lakes;
oo.random_lakes = 0;
mode = get_opt(lcfg, 'mode', 'table');
lake_file = nevis_region_resolve_path(get_opt(lcfg, 'file', ''), cfg);

if isempty(lake_file) || exist(lake_file, 'file') ~= 2
    error('nevis_import_region_lakes:MissingFile', ...
        'Lake file not found: %s', lake_file);
end

data = nevis_region_load_var(lake_file, get_opt(lcfg, 'variable', ''));

switch lower(mode)
    case 'environs_lakes_catalogue'
        lakes = import_environs_lakes(data);
    case 'table'
        lakes = import_lake_table(data, lcfg);
    otherwise
        error('nevis_import_region_lakes:UnsupportedMode', ...
            'Unsupported lake mode: %s', mode);
end

scale = get_opt(lcfg, 'spinup_volume_scale', 1);
time_offset = get_opt(lcfg, 'time_offset_days', 0);

pp.x_l = lakes.x_m(:) / ps.x;
pp.y_l = lakes.y_m(:) / ps.x;
pp.V_l = scale * lakes.volume_m3(:) / (ps.Q0 * ps.t);
pp.t_drainage = 0.5 * pd.td / ps.t * ...
    (lakes.drainage_time_start(:) + lakes.drainage_time_end(:)) + ...
    time_offset * pd.td / ps.t;
pp.t_duration = 0.5 * pd.td / ps.t * lakes.drainage_duration_days(:);

if isempty(pp.x_l)
    pp.ni_l = [];
    pp.sum_l = [];
else
    [pp.ni_l, pp.sum_l] = nevis_lakes(pp.x_l, pp.y_l, gg, oo);
end

region.lakes = lakes;
region.lake_file = lake_file;

end

function lakes = import_environs_lakes(environs_lakes)
inds = find(environs_lakes.drainage_type_num == 1 | ...
            environs_lakes.drainage_type_num == 2);
lakes = struct;
lakes.id = environs_lakes.laketypeing_dates(inds, 1);
lakes.x_m = 1e3 * environs_lakes.X_km(inds)';
lakes.y_m = 1e3 * environs_lakes.Y_km(inds)';
lakes.drainage_type_num = environs_lakes.drainage_type_num(inds)';
lakes.drainage_time_start = environs_lakes.laketypeing_dates(inds, 3);
lakes.drainage_time_end = environs_lakes.laketypeing_dates(inds, 4);
lakes.max_surface_area_m2 = environs_lakes.max_surface_area(inds)';

n_lakes = numel(inds);
lakes.volume_pre_drainage_m3 = zeros(n_lakes, 1);
lakes.volume_post_drainage_m3 = zeros(n_lakes, 1);
for i = 1:n_lakes
    row = inds(i);
    col_start = lakes.drainage_time_start(i);
    col_end = lakes.drainage_time_end(i);
    lakes.volume_pre_drainage_m3(i) = environs_lakes.lake_volume_interp(row, col_start);
    lakes.volume_post_drainage_m3(i) = environs_lakes.lake_volume_interp(row, col_end);
end

lakes.drainage_duration_days = lakes.drainage_time_end - lakes.drainage_time_start;
hf_mask = lakes.drainage_type_num == 1;
lakes.drainage_duration_days(hf_mask) = 2 / 24;
lakes.volume_m3 = lakes.volume_pre_drainage_m3;
end

function lakes = import_lake_table(data, lcfg)
lakes = struct;
lakes.x_m = get_table_field(data, get_opt(lcfg, 'x_field', 'x_m'));
lakes.y_m = get_table_field(data, get_opt(lcfg, 'y_field', 'y_m'));
lakes.volume_m3 = get_table_field(data, get_opt(lcfg, 'volume_field', 'volume_m3'));
lakes.drainage_time_start = get_table_field(data, get_opt(lcfg, 'start_field', 'drainage_time_start'));
lakes.drainage_time_end = get_table_field(data, get_opt(lcfg, 'end_field', 'drainage_time_end'));
lakes.drainage_duration_days = lakes.drainage_time_end - lakes.drainage_time_start;
end

function values = get_table_field(data, name)
if istable(data)
    values = data.(name);
elseif isstruct(data)
    values = data.(name);
else
    error('nevis_import_region_lakes:UnsupportedTable', ...
        'Lake table data must be a table or struct.');
end
values = values(:);
end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
