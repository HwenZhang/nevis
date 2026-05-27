function [pp, pd, ps, gg, aa, vv, oo, region] = nevis_setup_region_case(cfg, pd, ps, pp, oo)
%NEVIS_SETUP_REGION_CASE Standardized data import and model setup.

if ~isfield(cfg, 'casename') || isempty(cfg.casename)
    error('nevis_setup_region_case:MissingCaseName', ...
        'cfg.casename is required.');
end
if ~strcmp(cfg.casename, oo.casename)
    error('nevis_setup_region_case:CaseNameMismatch', ...
        'Config casename "%s" does not match script casename "%s".', ...
        cfg.casename, oo.casename);
end

casename = oo.casename;
rn = oo.rn;
fn = oo.fn;
oo.dataset = cfg.casename;
oo.dn = cfg.data_root;

[gg, b, s, region, oo] = nevis_import_region_geometry(cfg, ps, oo);

[pd, ps, pp, oo] = nevis_update_parameters_ice(pd, ps, pp, oo);
[pd, oo, pp] = nevis_apply_config_parameters(cfg, pd, oo, pp);
gg = nevis_label_ice_test(gg, oo);

[aa, vv] = nevis_initialize(b, s, gg, pp, oo);
[aa, region] = nevis_import_region_inversion(cfg, aa, pp, ps, oo, region);
[vv, aa, region] = nevis_import_region_initial_hydrology(cfg, vv, aa, pd, ps, pp, gg, region);
[vv, aa, region] = nevis_import_region_velocity(cfg, vv, aa, pd, ps, gg, region);

if isfield(region, 'initial_hydrology') && isstruct(region.initial_hydrology)
    vv_hydro = region.initial_hydrology;
    if isfield(vv_hydro, 'u') && numel(vv_hydro.u) == gg.eIJ
        vv.u = vv_hydro.u;
        vv.u(gg.eout2) = NaN;
    end
    if isfield(vv_hydro, 'v') && numel(vv_hydro.v) == gg.fIJ
        vv.v = vv_hydro.v;
        vv.v(gg.fout2) = NaN;
    end
    if isfield(vv, 'u'), aa.u_obs = vv.u; end
    if isfield(vv, 'v'), aa.v_obs = vv.v; end
end

aa.phi_b = max(aa.phi_0, aa.phi_a);

[pp, region] = nevis_import_region_moulins(cfg, pp, gg, oo, ps, region);
[pp, lakes, region] = nevis_import_region_lakes(cfg, pp, pd, ps, gg, oo, region);
[pp, oo, region] = nevis_import_region_runoff(cfg, pp, pd, ps, gg, lakes, oo, region);
[pp, oo, region] = nevis_import_region_stations(cfg, pp, gg, oo, ps, region);

oo.dt = get_nested(cfg, {'run', 'dt_days'}, 1/24) * pd.td / ps.t;
oo.save_timesteps = get_nested(cfg, {'run', 'save_timesteps'}, true);
oo.save_pts_all = get_nested(cfg, {'run', 'save_pts_all'}, true);
oo.t_span = get_nested(cfg, {'run', 't_span_days'}, 1:365) * pd.td / ps.t;
oo.casename = casename;
oo.dataset = cfg.casename;
oo.dn = cfg.data_root;
oo.rn = rn;
oo.fn = fn;

end

function value = get_nested(s, path, default)
value = default;
cur = s;
for i = 1:numel(path)
    name = path{i};
    if ~isstruct(cur) || ~isfield(cur, name) || isempty(cur.(name))
        return
    end
    cur = cur.(name);
end
value = cur;
end
