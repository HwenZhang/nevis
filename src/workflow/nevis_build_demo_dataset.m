function report = nevis_build_demo_dataset(workflow)
%NEVIS_BUILD_DEMO_DATASET Build a dataset package for the demo workflow.
%
% This wraps the data pipeline and copies the auxiliary files needed by the
% regional forward/inversion configs. It does not run model simulations.

addpath(fullfile(workflow.root, 'data', 'pipelines'));

report = build_nevis_dataset(workflow.dataset.root, workflow.raw.root, ...
    workflow.bbox_km);
nevis_write_demo_moulin_coordinates(workflow);
nevis_write_demo_lake_catalogue(workflow);
nevis_write_dataset_manifest(workflow);
end
