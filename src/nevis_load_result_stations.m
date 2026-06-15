function stations = nevis_load_result_stations(region, oo)
%NEVIS_LOAD_RESULT_STATIONS Load station observations bound to a result.
%
% Prefer station data saved in region by the standardized regional workflow.
% For older result files, fall back to the dataset manifest declared by
% oo.dataset_root or oo.dn. This function never falls back to the old ./data
% root, so missing dataset bindings fail explicitly.

if nargin >= 1 && isstruct(region) && isfield(region, 'stations') && ...
        ~isempty(region.stations)
    stations = region.stations;
    return
end

if nargin >= 1 && isstruct(region) && isfield(region, 'station_file') && ...
        ~isempty(region.station_file)
    stations = load_station_file(region.station_file, '');
    return
end

dataset_root = '';
if nargin >= 2 && isstruct(oo)
    if isfield(oo, 'dataset_root') && ~isempty(oo.dataset_root)
        dataset_root = oo.dataset_root;
    elseif isfield(oo, 'dn') && ~isempty(oo.dn)
        dataset_root = oo.dn;
    end
end
if isempty(dataset_root)
    error('nevis_load_result_stations:MissingDatasetRoot', ...
        'Result does not contain region.stations, region.station_file, or oo.dataset_root.');
end

manifest_file = fullfile(dataset_root, 'dataset_manifest.m');
if exist(manifest_file, 'file') ~= 2
    error('nevis_load_result_stations:MissingDatasetManifest', ...
        'Dataset manifest not found: %s', manifest_file);
end
dataset = nevis_load_dataset_manifest(manifest_file);
if ~isfield(dataset, 'stations') || ~isfield(dataset.stations, 'file')
    error('nevis_load_result_stations:MissingStationsManifest', ...
        'Dataset manifest must declare dataset.stations.file.');
end

variable = '';
if isfield(dataset.stations, 'variable')
    variable = dataset.stations.variable;
end
stations = load_station_file(fullfile(dataset_root, dataset.stations.file), variable);
end

function stations = load_station_file(station_file, variable)
if exist(station_file, 'file') ~= 2
    error('nevis_load_result_stations:MissingStationFile', ...
        'Station file not found: %s', station_file);
end

data = load(station_file);
if ~isempty(variable)
    if ~isfield(data, variable)
        error('nevis_load_result_stations:MissingStationVariable', ...
            'Station file "%s" does not contain variable "%s".', ...
            station_file, variable);
    end
    stations = data.(variable);
    return
end

if isfield(data, 'station_data')
    stations = data.station_data;
    return
end

names = fieldnames(data);
if numel(names) ~= 1
    error('nevis_load_result_stations:AmbiguousStationFile', ...
        'Station file "%s" has multiple variables; declare dataset.stations.variable.', ...
        station_file);
end
stations = data.(names{1});
end
