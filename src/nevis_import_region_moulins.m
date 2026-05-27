function [pp, region] = nevis_import_region_moulins(cfg, pp, gg, oo, ps, region)
%NEVIS_IMPORT_REGION_MOULINS Load moulins from node indices or coordinates.

mcfg = cfg.moulins;
oo.keep_all_moulins = get_opt(mcfg, 'keep_all', false);
oo.move_moulins = get_opt(mcfg, 'move_to_domain', true);

mode = get_opt(mcfg, 'mode', 'coordinates');
moulin_file = nevis_region_resolve_path(get_opt(mcfg, 'file', ''), cfg);

if isempty(moulin_file) || exist(moulin_file, 'file') ~= 2
    error('nevis_import_region_moulins:MissingFile', ...
        'Moulin file not found: %s', moulin_file);
end

data = nevis_region_load_var(moulin_file, get_opt(mcfg, 'variable', ''));

switch lower(mode)
    case 'node_indices'
        ni_field = get_opt(mcfg, 'ni_field', 'ni_m');
        sum_field = get_opt(mcfg, 'sum_field', 'sum_m');
        if ~isfield(data, ni_field) || ~isfield(data, sum_field)
            error('nevis_import_region_moulins:MissingField', ...
                'Moulin file must contain fields "%s" and "%s".', ni_field, sum_field);
        end
        pp.ni_m = data.(ni_field);
        pp.sum_m = data.(sum_field);
    case 'coordinates'
        x = data.(get_opt(mcfg, 'x_field', 'x_m')) / ps.x;
        y = data.(get_opt(mcfg, 'y_field', 'y_m')) / ps.x;
        [pp.ni_m, pp.sum_m] = nevis_moulins(x, y, gg, oo);
    otherwise
        error('nevis_import_region_moulins:UnsupportedMode', ...
            'Unsupported moulin mode: %s', mode);
end

if ~isempty(pp.sum_m) && size(pp.sum_m, 2) ~= gg.nIJ
    error('nevis_import_region_moulins:SizeMismatch', ...
        'Moulin catchment matrix has %d columns; expected %d.', size(pp.sum_m, 2), gg.nIJ);
end

region.moulins = data;
region.moulin_file = moulin_file;

end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
