% This script extends the field for specific time slices
% and saves the results back to the MATLAB mat file.

% load the filenames from the text file line by line
clear; clc;
addpath(genpath('./src'));  % add the nevis_python path
python_data = load('extend_field.mat');  % load the extension parameters
casename = python_data.casename;  % get the case name
load(['./results/' casename '/' casename]);

% load the data for each filename and extend the field
for i = 1:length(python_data.filenames)
    filename = python_data.filenames{i};
    vva = load(filename,'vv');
    vv = vva.vv;
    aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
    % extend the field for the current filename
    [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo);     % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
    save(filename,'vv2','-append');  % save the extended field back to the file
end