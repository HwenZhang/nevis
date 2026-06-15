function lake_file = nevis_write_demo_lake_catalogue(workflow)
%NEVIS_WRITE_DEMO_LAKE_CATALOGUE Copy raw lake catalogue into case dataset.
%
% The workflow keeps runtime data isolated by placing the lake catalogue
% inside data/datasets/<case_name>/ rather than reading it from another
% dataset package.

source_file = workflow.raw.lake_catalogue;
if exist(source_file, 'file') ~= 2
    error('nevis_write_demo_lake_catalogue:MissingRawData', ...
        'Raw lake catalogue is missing: %s', source_file);
end

lake_file = fullfile(workflow.dataset.root, ...
    'environs_lakes_2022B_250416.mat');
ensure_parent_dir(lake_file);
copyfile(source_file, lake_file);
end

function ensure_parent_dir(path_in)
parent = fileparts(path_in);
if ~isempty(parent) && exist(parent, 'dir') ~= 7
    mkdir(parent);
end
end
