% This script extends the field for specific time slices
% and saves the results back to the MATLAB mat file.

script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);  % Go up one level to main directory
cd(parent_dir);  % Change to main directory

% load the filenames from the text file line by line
clear variables -except parent_dir script_dir; 
clc;

addpath(genpath('./src'));  % add the nevis_python path
python_data = load('extend_field.mat');  % load the extension parameters
casename = python_data.casename;  % get the case name
load(['./results/' casename '/' casename]);

% Older saved result files may not contain the blister-specific grid labels
% needed by nevis_backbone when expanding diagnostic flux fields.
if ~isfield(gg,'ebdy_blister')
    noutb = gg.nout;
    if isempty(noutb)
        y_max = max(gg.ny(:));
        x_max = max(gg.nx(:));
        noutb = union(find(abs(gg.ny - y_max) < 1e-10), find(abs(gg.nx - x_max) < 1e-10));
    end
    gg = nevis_mask_blister(gg,noutb,oo);
    if isfield(gg,'n1m')
        gg.n1_blister = gg.n1m;
    end
    gg = nevis_label_blister(gg,gg.n1_blister,oo);
end

if ~isfield(gg,'es2')
    gg.ns2 = gg.ns;
    gg.es2 = gg.es;
    gg.fs2 = gg.fs;
    gg.cs2 = gg.cs;
end

oo_diag = oo;
oo_diag.evaluate_residual = 0;
oo_diag.evaluate_jacobian = 0;

% load the data for each filename and extend the field
for i = 1:length(python_data.filenames)
    filename = python_data.filenames{i};
    vva = load(filename,'vv');
    vv = vva.vv;
    if ~isfield(vv,'u'), vv.u = zeros(gg.eIJ,1); end
    if ~isfield(vv,'v'), vv.v = zeros(gg.fIJ,1); end
    aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
    % extend the field for the current filename
    if ~isfield(pp,'B_reg'), pp.B_reg = 0; end
        
    [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo_diag); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
    save(filename,'vv2','-append');  % save the extended field back to the file
end
