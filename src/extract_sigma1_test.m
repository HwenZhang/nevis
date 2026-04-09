script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
cd(parent_dir);
casename = 'n2d_regional_racmo_eps1e_02_kappa5e_11_mu5e0_partition5e_01_k01e_01_spinup';
result_dir = fullfile('results', casename);
file1 = fullfile(parent_dir, result_dir, '0001.mat');
t1 = load(file1);
disp(fieldnames(t1));
if isfield(t1, 'ice_fields')
    disp(fieldnames(t1.ice_fields));
end
