function [pp, oo, region] = nevis_import_region_runoff(cfg, pp, pd, ps, gg, lakes, oo, region)
%NEVIS_IMPORT_REGION_RUNOFF Load runoff forcing and attach runoff/input functions.

rcfg = cfg.runoff;
oo.surface_runoff = get_opt(rcfg, 'surface_runoff', true);
oo.RACMO_runoff = get_opt(rcfg, 'RACMO_runoff', false);
oo.distributed_input = get_opt(rcfg, 'distributed_input', true);

pp.meltE = @(t) (get_opt(rcfg, 'meltE_mm_per_day', 0) / 1000 / pd.td / ps.m) .* ...
    (1 - exp(-t / (get_opt(rcfg, 'meltE_ramp_days', 20) * pd.td / ps.t)));
pp.input_function = @(t) 0;

runoff_file = nevis_region_resolve_path(get_opt(rcfg, 'file', ''), cfg);
if isempty(runoff_file) || exist(runoff_file, 'file') ~= 2
    error('nevis_import_region_runoff:MissingFile', ...
        'Runoff file not found: %s', runoff_file);
end

runoff_data = nevis_region_load_var(runoff_file, get_opt(rcfg, 'variable', ''));
runoff_data = double(runoff_data);

if size(runoff_data, 2) ~= gg.nIJ
    error('nevis_import_region_runoff:SizeMismatch', ...
        'Runoff must have %d columns, one for each model node.', gg.nIJ);
end

runoff_used = runoff_data;
smooth_days = get_opt(rcfg, 'smooth_days', 0);
if smooth_days > 1
    runoff_smoothed = smoothdata(runoff_data, 1, 'movmean', smooth_days);
    sum_raw = sum(runoff_data, 1);
    sum_smooth = sum(runoff_smoothed, 1);
    scale_col = ones(size(sum_raw));
    ii = sum_smooth > 0;
    scale_col(ii) = sum_raw(ii) ./ sum_smooth(ii);
    runoff_smoothed = runoff_smoothed .* scale_col;
    if get_opt(rcfg, 'use_smoothed', false)
        runoff_used = runoff_smoothed;
    end
    region.runoff_smoothed = runoff_smoothed;
end

lake_runoff_ratio = 0;
if get_opt(rcfg, 'subtract_lake_volume', false) && ~isempty(lakes.volume_m3)
    dx_m = median(abs(diff(unique(gg.nx(:))))) * ps.x;
    dy_m = median(abs(diff(unique(gg.ny(:))))) * ps.x;
    cell_area_m2 = dx_m * dy_m;
    total_lake_volume_m3 = sum(lakes.volume_m3);
    total_runoff_volume_m3 = sum(runoff_data(:)) / 1000 * cell_area_m2;
    if total_runoff_volume_m3 > 0
        lake_runoff_ratio = total_lake_volume_m3 / total_runoff_volume_m3;
    end
    disp(['Total lake volume: ' num2str(total_lake_volume_m3, '%.3e') ' m^3']);
    disp(['Total runoff volume: ' num2str(total_runoff_volume_m3, '%.3e') ' m^3']);
    disp(['Lake/Runoff ratio: ' num2str(lake_runoff_ratio, '%.4f')]);
end

runoff_used = max(0, (1 - lake_runoff_ratio) * runoff_used);
pp.runoff_function = @(t) runoff(((t * ps.t) / pd.td), runoff_used) ./ ps.m;

region.runoff = runoff_data;
region.runoff_used = runoff_used;
region.runoff_file = runoff_file;
region.lake_runoff_ratio = lake_runoff_ratio;

end

function value = get_opt(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
