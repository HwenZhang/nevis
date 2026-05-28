function dataset = nevis_load_dataset_manifest(manifest_file)
%NEVIS_LOAD_DATASET_MANIFEST Load a shared regional dataset manifest.

if exist(manifest_file, 'file') ~= 2
    error('nevis_load_dataset_manifest:MissingManifest', ...
        'Dataset manifest not found: %s', manifest_file);
end

[~, ~, ext] = fileparts(manifest_file);
switch lower(ext)
    case '.m'
        dataset = [];
        run(manifest_file);
        if isempty(dataset)
            error('nevis_load_dataset_manifest:MissingDataset', ...
                'Dataset manifest must assign a variable named dataset.');
        end
    case '.mat'
        data = load(manifest_file);
        if ~isfield(data, 'dataset')
            error('nevis_load_dataset_manifest:MissingDataset', ...
                'MAT dataset manifest must contain a variable named dataset.');
        end
        dataset = data.dataset;
    case '.json'
        dataset = jsondecode(fileread(manifest_file));
    otherwise
        error('nevis_load_dataset_manifest:UnsupportedManifest', ...
            'Unsupported dataset manifest extension: %s', ext);
end

dataset.manifest_file = manifest_file;

end
