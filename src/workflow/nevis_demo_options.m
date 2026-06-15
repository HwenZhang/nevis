function opts = nevis_demo_options(case_name, bbox_km)
%NEVIS_DEMO_OPTIONS User-facing defaults for the NEVIS demo workflow.
%
% case_name is used as the dataset package name and as the prefix for
% generated spinup/forward case configs.

if nargin < 1 || isempty(case_name)
    case_name = 'nevis_2022_140km_rebuilt';
end
if nargin < 2
    bbox_km = [];
end

opts = struct;
opts.root = './';
opts.case_name = case_name;
opts.bbox_km = bbox_km;

opts.paths.cases_dir = fullfile(opts.root, 'cases');
opts.paths.case_dir = fullfile(opts.paths.cases_dir, case_name);
opts.paths.data_root = fullfile(opts.root, 'data', 'datasets');

opts.raw.root = fullfile(opts.root, 'data', 'raw');
opts.raw.moulin_shapefile = fullfile(opts.raw.root, 'gnss_2022', ...
    'NatalieTurner_Moulins_2022', 'moulins.shp');
opts.raw.lake_catalogue = fullfile(opts.raw.root, 'gnss_2022', ...
    'environs_lakes_2022B_250416.mat');
opts.raw.m1_lake_file = fullfile(opts.raw.root, 'bedmachine', 'lake.mat');

opts.dataset.name = case_name;
opts.dataset.root = fullfile(opts.paths.data_root, opts.dataset.name);
opts.dataset.manifest = 'dataset_manifest.m';

opts.spinup.template_case = 'nevis_regional_inversion_spinup';
opts.spinup.t_span_days = 1:365;
opts.spinup.timestep_file = '0365.mat';

opts.paths.templates_dir = fullfile(opts.paths.cases_dir, 'templates');
opts.inversion.template_config = fullfile(opts.paths.templates_dir, ...
    'nevis_inversion.m');
opts.inversion.base_name = 'weertman_ns3_default';
opts.inversion.sliding_law = 'weertman';

opts.forward.template_case = 'nevis_regional';

opts.drainage.template_case = 'nevis_regional';
opts.drainage.suffix = '_drainage';
opts.drainage.initial_timestep_file = '0365.mat';
opts.drainage.lake_catalogue_file = 'environs_lakes_2022B_250416.mat';
opts.drainage.t_span_days = [1:180, 180.01:0.01:220, 221:365];
end
