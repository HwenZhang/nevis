function [pp, oo, region] = nevis_import_region_stations(cfg, pp, gg, oo, ps, region)
%NEVIS_IMPORT_REGION_STATIONS Load station/GPS points and assemble output nodes.

scfg = cfg.stations;
station_file = nevis_region_resolve_path(get_opt(scfg, 'file', ''), cfg);
pp.ni_gps = [];

if ~isempty(station_file) && exist(station_file, 'file') == 2
    stations = nevis_region_load_var(station_file, get_opt(scfg, 'variable', ''));
    x = get_values(stations, get_opt(scfg, 'x_field', 'x_m')) / ps.x;
    y = get_values(stations, get_opt(scfg, 'y_field', 'y_m')) / ps.x;
    pp.x_gps = x;
    pp.y_gps = y;
    pp.ni_gps = nevis_gps_array(x, y, gg, oo);
    region.stations = stations;
    region.station_file = station_file;
else
    error('nevis_import_region_stations:MissingFile', ...
        'Station file not found: %s', station_file);
end

pts = [];
if isfield(pp, 'ni_l') && ~isempty(pp.ni_l), pts = [pts pp.ni_l']; end
if isfield(pp, 'ni_m') && ~isempty(pp.ni_m), pts = [pts pp.ni_m']; end
if isfield(pp, 'ni_gps') && ~isempty(pp.ni_gps), pts = [pts pp.ni_gps]; end
oo.pts_ni = pts;

end

function values = get_values(data, field)
if isstruct(data)
    if numel(data) > 1
        values = [data.(field)];
    else
        values = data.(field);
    end
elseif istable(data)
    values = data.(field);
else
    error('nevis_import_region_stations:UnsupportedData', ...
        'Station data must be a struct or table.');
end
values = values(:)';
end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
