
%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS 2-dimensional model for a regional domain of the Greenland ice sheet, specifically the NEVIS region in western Greenland. It is used to test the capability of the ice sheet module to simulate the steady-state ice velocity.
clc,clear
format compact

%% read in the initial condition
casename = 'n2d_regional_analytic_V2e8_eps1e_02_kappa5e_11_mu2e1_partition8e_01_k01e_01_drainage_highelev'; % drainage system filename
initname = 'n2d_regional_analytic_eps1e_02_kappa5e_11_mu2e1_partition8e_01_k01e_01_spinup';

%% parameters
data = load(['./results/' initname '/' initname]);
pd = data.pd;                                % load parameters from the initial condition
ps = data.ps;                                % load state variables from the initial condition
pp = data.pp;                                % load scaled parameters from the initial condition
aa = data.aa;                                % load state variables from the initial condition
oo = data.oo;                                % load options from the initial condition

oo.casename = casename;                      % drainage system filename
oo.initname = initname;                      % initial condition filename, for spinup
oo.display_residual = 0;
oo.visualize_vel = 0;
oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
addpath(oo.code);                              % add path to code
mkdir(oo.rn);                                  % create directory for results    
oo.iterative_solver = 0;

%% grid and geometry
load([oo.dn '/' 'morlighem_for_nevis_140km']); % load Morlighem bedmap (previously collated)
dd = morlighem_for_nevis_140km; dd.skip = 6;
gg = nevis_grid(dd.X_m(1:dd.skip:end,1)/ps.x,dd.Y_m(1,1:dd.skip:end)/ps.x,oo); 
b = reshape(dd.B_m(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
s = reshape(dd.S_m(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
% save a copy of the bed and surface
gg.bed = b;
gg.surface = s;

x = dd.X_m(1:dd.skip:end,1)/ps.x;
y = dd.Y_m(1,1:dd.skip:end)/ps.x;

%% mask with minimum ice thickness
H = max(s-b,0);
Hmin = 0/ps.z; 
% y_max = max(max(gg.ny'));
nout = find(H<=Hmin);
if isempty(nout)
    x_max = max(gg.nx);
    y_max = max(gg.ny);
    nout = find(abs(gg.nx - x_max) < 1e-10);
end
noutb = nout;
if isempty(noutb)
    x_max = max(max(gg.nx));
    y_max = max(max(gg.ny'));
    noutb = union(find(abs(gg.ny - y_max) < 1e-10), find(abs(gg.nx - x_max) < 1e-10));
end
gg = nevis_mask(gg,nout); 
gg = nevis_mask_blister(gg,noutb);
gg.n1m = gg.n1;                                   % label all edge nodes as boundary nodes for pressure

%% label boundary nodes
gg = nevis_label(gg,gg.n1m);
gg = nevis_label_blister(gg,gg.n1_blister,oo);    % label blister boundary nodes
oo.adjust_boundaries = 1;                         % enable option of changing conditions

gg = nevis_label_ice_test(gg, oo); % add boundary labels needed for ice velocity

%% plot grid
% nevis_plot_grid(gg); return;                    % check to see what grid looks like

%% initialize variables
init_cond = load(['./results/' oo.initname '/' '0365.mat']); % load initial condition
vv = init_cond.vv;                                % load state variables from the initial 

%% boundary conditions
aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric

%% moulins 
oo.density_moulins = 1;
oo.keep_all_moulins = 0;
oo.random_lakes = 0;
% load moulin locations
tmp = load([oo.dn '/' '/nevis_170207a.mat'],'pp');
pp.ni_m = tmp.pp.ni_m; pp.sum_m = tmp.pp.sum_m; clear tmp;

%% supraglacial lakes
lakes = nevis_import_lakes(2022);                     % import lake data
% multiple lakes from the catalogue
pp.x_l = [40e3]/ps.x;                          % x-coord of lakes
pp.y_l = [-15e3]/ps.x;                          % y-coord of lakes
pp.V_l = [200000000.00000003]/(ps.Q0*ps.t);                    % volume of lakes         
pp.t_drainage = pd.td/ps.t*[10];  % time of lake drainages (assumed to be the middle time of the Gaussian)
pp.t_duration = pd.td/ps.t*[0.25]; % duration of lake drainages
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments

%% surface runoff
load([oo.dn '/' 'runoff_2022_nevis140.mat']);                   % load data for year of interest (previously collated)
% smooth the runoff data and scale it to match the total runoff volume in the original data, to avoid spurious spikes in the input that can cause numerical issues. The smoothing window can be adjusted with oo.runoff_smooth_days (default is 5 days).
if ~isfield(oo,'runoff_smooth_days'), oo.runoff_smooth_days = 5; end
runoff_smoothed = smoothdata(runoff_2022_nevis140, 1, 'movmean', oo.runoff_smooth_days);
sum_raw = sum(runoff_2022_nevis140, 1);
sum_smooth = sum(runoff_smoothed, 1);
scale_col = ones(size(sum_raw));
ii = sum_smooth > 0;
scale_col(ii) = sum_raw(ii) ./ sum_smooth(ii);
runoff_2022_nevis140_smooth = runoff_smoothed .* scale_col;

% RACMO distributed input
oo.surface_runoff = 1;                           % If set to 1 turns on surface runoff input
                                                                % 0 is prescribe moulin input with a function
oo.RACMO_runoff = 0;                               % If set to 1 turns on RACMO surface runoff data
                                                                % 0 is a simple sinusoidal function
oo.distributed_input = 1;                     % If set to 1 distributes input across the domain
                                                                % 0 inputs to discrete moulins

pp.meltE = @(t) (0/1000/pd.td/ps.m)*(1-0*exp(-t/(30*pd.td/ps.t))); 
dd = morlighem_for_nevis_140km;
% Calculate ratio of total lake volume to total annual runoff volume
dx_m = abs(dd.X_m(2,1) - dd.X_m(1,1)) * dd.skip;     % grid spacing in x (m)
dy_m = abs(dd.Y_m(1,2) - dd.Y_m(1,1)) * dd.skip;     % grid spacing in y (m)
cell_area_m2 = dx_m * dy_m;                            % grid cell area (m^2)
total_lake_volume_m3 = sum(lakes.volume_m3);            % total lake volume (m^3)
total_runoff_volume_m3 = sum(runoff_2022_nevis140(:)) / 1000 * cell_area_m2;  % mm w.e. -> m, * area -> m^3
lake_runoff_ratio = total_lake_volume_m3 / total_runoff_volume_m3;
disp(['Total lake volume: ' num2str(total_lake_volume_m3, '%.3e') ' m^3']);
disp(['Total runoff volume: ' num2str(total_runoff_volume_m3, '%.3e') ' m^3']);
disp(['Lake/Runoff ratio: ' num2str(lake_runoff_ratio, '%.4f')]);

pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td), (1-lake_runoff_ratio)*runoff_2022_nevis140)./ps.m;  % distributed input (m/sec)
pp.input_function = @(t) 0;

%% timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 

% Add GPS station points downstream of the moulin every 5km
pp.ni_gps = nevis_gps_array([40e3,40e3,20e3,0,-40e3]/ps.x, [-15e3,-5e3,-15e3,-25e3,-30e3]/ps.x, gg, oo); % GPS station points
oo.pts_ni = [pp.ni_l' pp.ni_m' pp.ni_gps];    
% oo.t_span = [(1:1:365*1)*pd.td/ps.t];
oo.t_span = [(1:1:9)*pd.td/ps.t (9.5:0.1:10.5)*pd.td/ps.t (11:1:140)*pd.td/ps.t];

%% save initial parameters
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','vv','oo');
[tt,vv,info] = nevis_timesteps(oo.t_span,vv,aa,pp,gg,oo);

%% expand/update variables
aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); 
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); 
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','oo','tt');

%% Simple animate
% nevis_regional_animation

