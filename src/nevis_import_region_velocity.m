function [vv, aa, region] = nevis_import_region_velocity(cfg, vv, aa, pd, ps, gg, region)
%NEVIS_IMPORT_REGION_VELOCITY Load observed velocity and map it to model edges.

vcfg = cfg.velocity;
mode = get_opt(vcfg, 'mode', 'file');
if strcmpi(mode, 'none')
    return
end

velocity_file = nevis_region_resolve_path(vcfg.file, cfg);
if exist(velocity_file, 'file') ~= 2
    error('nevis_import_region_velocity:MissingFile', ...
        'Velocity file not found: %s', velocity_file);
end

vel = nevis_region_load_var(velocity_file, get_opt(vcfg, 'variable', ''));
u = get_required_field(vel, vcfg.u_field, velocity_file);
v = get_required_field(vel, vcfg.v_field, velocity_file);

switch lower(get_opt(vcfg, 'units', 'm_per_year'))
    case {'m_per_year', 'm/yr', 'm_year'}
        un = u / (pd.ty * ps.u);
        vn = v / (pd.ty * ps.u);
    case {'m_per_second', 'm/s'}
        un = u / ps.u;
        vn = v / ps.u;
    case {'dimensionless', 'nondimensional'}
        un = u;
        vn = v;
    otherwise
        error('nevis_import_region_velocity:UnsupportedUnits', ...
            'Unsupported velocity units: %s', vcfg.units);
end

if numel(un) ~= gg.nIJ || numel(vn) ~= gg.nIJ
    error('nevis_import_region_velocity:SizeMismatch', ...
        'Velocity fields must have %d node values after preprocessing.', gg.nIJ);
end

un_filled = fill_nans_nearest(un(:), gg);
vn_filled = fill_nans_nearest(vn(:), gg);

vv.u = gg.emean2 * un_filled;
vv.v = gg.fmean2 * vn_filled;
vv.u(gg.eout2) = NaN;
vv.v(gg.fout2) = NaN;
aa.u_obs = vv.u;
aa.v_obs = vv.v;

region.velocity = vel;
region.velocity_file = velocity_file;

end

function values = fill_nans_nearest(values, gg)
nan_mask = isnan(values);
if any(nan_mask)
    valid_idx = find(~nan_mask);
    nan_idx = find(nan_mask);
    values(nan_idx) = griddata(gg.nx(valid_idx), gg.ny(valid_idx), ...
        values(valid_idx), gg.nx(nan_idx), gg.ny(nan_idx), 'nearest');
end
end

function value = get_required_field(s, name, filename)
if ~isfield(s, name)
    error('nevis_import_region_velocity:MissingField', ...
        'Velocity file "%s" is missing field "%s".', filename, name);
end
value = s.(name);
end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
