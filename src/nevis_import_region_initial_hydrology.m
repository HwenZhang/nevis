function [vv, aa, region] = nevis_import_region_initial_hydrology(cfg, vv, aa, pd, ps, pp, gg, region)
%NEVIS_IMPORT_REGION_INITIAL_HYDROLOGY Load hydrology/velocity initial state.

hcfg = cfg.initial_hydrology;
mode = get_opt(hcfg, 'mode', 'file');

if strcmpi(mode, 'k_factor')
    k_factor = get_opt(hcfg, 'k_factor', 0.9);
    N = (1 - k_factor) * (aa.phi_0 - aa.phi_a);
    vv.phi = aa.phi_0 - N;
    vv.hs = cavity_sheet_depth(N, pd, ps, pp);
    region.initial_hydrology = struct('mode', 'k_factor');
    return
end

if strcmpi(mode, 'result_timestep')
    init_file = fullfile('./results', hcfg.result_case, hcfg.timestep_file);
    if exist(init_file, 'file') ~= 2
        error('nevis_import_region_initial_hydrology:MissingInitialTimestep', ...
            'Initial timestep file not found: %s', init_file);
    end
    data = load(init_file, 'vv');
    if ~isfield(data, 'vv')
        error('nevis_import_region_initial_hydrology:MissingInitialVV', ...
            'Initial timestep file must contain variable vv: %s', init_file);
    end
    vv = data.vv;
    region.initial_hydrology = struct('mode', 'result_timestep', ...
        'result_case', hcfg.result_case, 'timestep_file', hcfg.timestep_file);
    region.initial_hydrology_file = init_file;
    return
end

hydro_file = nevis_region_resolve_path(hcfg.file, cfg);
if exist(hydro_file, 'file') ~= 2
    error('nevis_import_region_initial_hydrology:MissingFile', ...
        'Initial hydrology file not found: %s', hydro_file);
end

vv_hydro = nevis_region_load_var(hydro_file, get_opt(hcfg, 'variable', ''));
if isfield(vv_hydro, 'N')
    N = vv_hydro.N(:);
else
    error('nevis_import_region_initial_hydrology:MissingN', ...
        'Initial hydrology state must contain field N.');
end

if numel(N) ~= numel(aa.b)
    error('nevis_import_region_initial_hydrology:SizeMismatch', ...
        'Initial N size (%d) does not match model nodes (%d).', numel(N), numel(aa.b));
end

vv.phi = aa.phi_0 - N;
vv.hs = cavity_sheet_depth(N, pd, ps, pp);

if isfield(vv_hydro, 'u') && numel(vv_hydro.u) == gg.eIJ
    vv.u = vv_hydro.u;
    vv.u(gg.eout2) = NaN;
end
if isfield(vv_hydro, 'v') && numel(vv_hydro.v) == gg.fIJ
    vv.v = vv_hydro.v;
    vv.v(gg.fout2) = NaN;
end

region.initial_hydrology = vv_hydro;
region.initial_hydrology_file = hydro_file;

end

function hs = cavity_sheet_depth(N, pd, ps, pp)
hs = ((((pd.u_b * pd.h_r / pd.l_r) ./ ...
    ((pd.u_b / pd.l_r) + (pd.K_c .* ((ps.phi * N).^3))))) ./ ps.h);
hs = max(hs, 0);
if isfield(pp, 'hs_reg')
    hs = max(hs, pp.hs_reg);
end
end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
