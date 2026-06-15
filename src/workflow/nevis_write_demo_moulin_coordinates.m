function moulin_file = nevis_write_demo_moulin_coordinates(workflow)
%NEVIS_WRITE_DEMO_MOULIN_COORDINATES Write moulin coordinates for a dataset.
%
% Reads the raw 2022 moulin shapefile and writes coordinates into the
% current case dataset package. Regional configs then use
% cfg.moulins.mode='coordinates' so catchments are rebuilt on the current
% domain. The model run never reads another dataset package.

source_file = workflow.raw.moulin_shapefile;
if exist(source_file, 'file') ~= 2
    error('nevis_write_demo_moulin_coordinates:MissingRawData', ...
        'Raw moulin shapefile is missing: %s', source_file);
end
if exist(workflow.raw.m1_lake_file, 'file') ~= 2
    error('nevis_write_demo_moulin_coordinates:MissingOriginData', ...
        'Raw M1 origin lake file is missing: %s', workflow.raw.m1_lake_file);
end

nt = shaperead(source_file);
[origin_x, origin_y] = m1_origin(workflow.raw.m1_lake_file);
radius = 6378137.0;
eccen = 0.08181919;
lat_true = 70;
lon_posy = -45;
[x, y] = polarstereo_fwd([nt.Y], [nt.X], radius, eccen, ...
    lat_true, lon_posy);

moulins = struct;
moulins.x_m = x(:) - origin_x;
moulins.y_m = y(:) - origin_y;
moulins.source_file = source_file;

moulin_file = fullfile(workflow.dataset.root, ...
    'moulins_2022_coordinates.mat');
ensure_parent_dir(moulin_file);
save(moulin_file, 'moulins');
end

function [ox, oy] = m1_origin(lake_file)
source = load(lake_file, 'lake');
moulin = source.lake(226, 4:5);
radius = 6378137.0;
eccen = 0.08181919;
lat_true = 70;
lon_posy = -45;
[ox, oy] = polarstereo_fwd(moulin(2), moulin(1), radius, eccen, ...
    lat_true, lon_posy);
end

function ensure_parent_dir(path_in)
parent = fileparts(path_in);
if ~isempty(parent) && exist(parent, 'dir') ~= 7
    mkdir(parent);
end
end
