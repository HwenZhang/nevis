function [gg, b, s, region, oo] = nevis_import_region_geometry(cfg, ps, oo)
%NEVIS_IMPORT_REGION_GEOMETRY Load geometry, build grid, mask, and labels.

gcfg = cfg.geometry;
geometry_file = nevis_region_resolve_path(gcfg.file, cfg);
dd = nevis_region_load_var(geometry_file, get_opt(gcfg, 'variable', ''));

required = {gcfg.x_field, gcfg.y_field, gcfg.bed_field, gcfg.surface_field};
for i = 1:numel(required)
    if ~isfield(dd, required{i})
        error('nevis_import_region_geometry:MissingField', ...
            'Geometry file "%s" is missing field "%s".', geometry_file, required{i});
    end
end

skip = get_opt(gcfg, 'skip', 1);
dd.skip = skip;

X = dd.(gcfg.x_field);
Y = dd.(gcfg.y_field);
B = dd.(gcfg.bed_field);
S = dd.(gcfg.surface_field);

if ~isequal(size(X), size(Y), size(B), size(S))
    error('nevis_import_region_geometry:SizeMismatch', ...
        'Geometry X, Y, bed, and surface fields must have the same size.');
end

x = X(1:skip:end, 1) / ps.x;
y = Y(1, 1:skip:end) / ps.x;
gg = nevis_grid(x, y, oo);

b = reshape(B(1:skip:end, 1:skip:end) / ps.z, gg.nIJ, 1);
s = reshape(S(1:skip:end, 1:skip:end) / ps.z, gg.nIJ, 1);
gg.bed = b;
gg.surface = s;

H = max(s - b, 0);
Hmin = get_opt(gcfg, 'min_thickness_m', 0) / ps.z;
nout = find(H <= Hmin);
if isempty(nout)
    x_max = max(gg.nx);
    nout = find(abs(gg.nx - x_max) < 1e-10);
end

noutb = nout;
if isempty(noutb)
    x_max = max(max(gg.nx));
    y_max = max(max(gg.ny'));
    noutb = union(find(abs(gg.ny - y_max) < 1e-10), ...
                  find(abs(gg.nx - x_max) < 1e-10));
end

gg = nevis_mask(gg, nout);
gg = nevis_mask_blister(gg, noutb);
gg.n1m = gg.n1;
gg.n1_blister = gg.n1m;

gg = nevis_label(gg, gg.n1m);
gg = nevis_label_blister(gg, gg.n1_blister, oo);

region = struct;
region.geometry = dd;
region.geometry_file = geometry_file;
region.x = x;
region.y = y;

end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
