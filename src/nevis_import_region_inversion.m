function [aa, region] = nevis_import_region_inversion(cfg, aa, pp, ps, oo, region)
%NEVIS_IMPORT_REGION_INVERSION Load or initialize sliding coefficients.

icfg = cfg.inversion;
mode = get_opt(icfg, 'mode', 'file');

if strcmpi(mode, 'uniform') || isempty(get_opt(icfg, 'file', ''))
    aa.C = ones(size(aa.b));
    aa.C2 = zeros(size(aa.b));
    region.inversion = struct('mode', 'uniform');
    return
end

inversion_file = nevis_region_resolve_path(icfg.file, cfg);
if exist(inversion_file, 'file') ~= 2
    error('nevis_import_region_inversion:MissingFile', ...
        'Inversion file not found: %s', inversion_file);
end

inv = load(inversion_file);
if isfield(inv, 'C_hat')
    partition_ratio = cfg.oo.partition_ratio;
    [C1_hat_dim, C2_hat_dim] = nevis_inv_partition( ...
        inv.C_hat, partition_ratio, inv.u_obs_noisy, inv.v_obs_noisy, ...
        inv.N_current, inv.aa, inv.pp, inv.gg, inv.oo, inv.ps);
else
    error('nevis_import_region_inversion:MissingC', ...
        ['Inversion file "%s" must contain C_hat so aa.C and aa.C2 ' ...
         'can be partitioned with cfg.oo.partition_ratio.'], inversion_file);
end

if numel(C1_hat_dim) ~= numel(aa.b)
    error('nevis_import_region_inversion:SizeMismatch', ...
        'Inversion C size (%d) does not match model nodes (%d).', ...
        numel(C1_hat_dim), numel(aa.b));
end

aa.C = C1_hat_dim(:) * (ps.u_b^(1 / pp.n_slide) / ps.tau);
aa.C2 = C2_hat_dim(:) * (ps.u_b^(1 / pp.n_slide) / ps.tau);
region.inversion = inv;
region.inversion.partition_ratio_applied = partition_ratio;
region.inversion.C1_hat_dim = C1_hat_dim(:);
region.inversion.C2_hat_dim = C2_hat_dim(:);
region.inversion_file = inversion_file;

end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
