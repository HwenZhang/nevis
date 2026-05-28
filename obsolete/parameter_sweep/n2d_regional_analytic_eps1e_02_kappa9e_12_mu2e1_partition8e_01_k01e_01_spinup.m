
%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS 2-dimensional model for a regional domain of the Greenland ice sheet, specifically the NEVIS region in western Greenland. It is used to test the capability of the ice sheet module to simulate the steady-state ice velocity.
clc,clear
format compact

%% read in the initial condition
oo.root = './';                                % filename root
oo.code = './src';                             % code directory   
oo.results = 'results';                        % path to the results folders
oo.dataset = 'nevis_regional';                 % dataset name
oo.casename = 'n2d_regional_analytic_eps1e_02_kappa9e_12_mu2e1_partition8e_01_k01e_01_spinup';                     % casename
oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
oo.dn = [oo.root, 'data/', oo.dataset, '/'];   % path to the data
addpath(oo.code);                              % add path to code
mkdir(oo.rn);                                  % create directory for results 

%% parameters
[pd,oo] = nevis_defaults([],oo);  

oo.evaluate_variables = 1;
oo.input_gaussian = 1;
oo.relaxation_term = 1;                         % 0 is alpha hb, 1 is alpha deltap hb
oo.initial_condition = 1;                       % 1 is default condition from 0365.mat, 0 is using steady-state drainage system, winter or summertime
oo.mean_perms = 1;
oo.modified_mean_perms = 0;
oo.display_residual = 0;
oo.visualize_vel = 0;
oo.N_coupling = 1; % turn on effective pressure coupling                                   
oo.U_coupling = 1; % turn on basal sliding coupling
oo.boundary_method = 'stress_l_vel_tbl';
% oo.boundary_method = 'stress_tblr';
oo.mask_boundary_method = 'stress_free';
oo.plot_residual = 0;
oo.max_iter_new = 50;
% ratio of reducable basal shear stress(C1) to total basal shear stress (C_total) for the inversion test, e.g., 0.25 means that 25% of the basal shear stress is from C1 and can be reduced by lowering N, while 75% is from C2 and is unaffected by N.
oo.partition_ratio = 0.8;
oo.iterative_solver = 0;

pd.alpha_b = 0;                                 % relaxation rate (s^-1)
pd.kappa_b = 9e-12;                           % relaxation coeff

% alter default parmaeters 
pd.mu = 20.0;                                    % water viscosity (Pa s)
pd.c_e_reg2 = 0.00/1e3/9.81;                    % elastic sheet thickness [m/Pa]
pd.u_b = 100/pd.ty;                             % sliding speed [m/s]
pd.sigma = 1e-3;                                % englacial void fraction
pd.h_r = 0.1;                                   % roughness height [m]
pd.l_r = 10;                                    % roughness length [m]
pd.l_c = 1000;                                  % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 0.1;                                  % sheet permeability constant
pd.tau_b = 60e3;                                % driving stress [Pa]
pd.melt = pd.G/pd.rho_w/pd.L;                   % geothermal heat derived basal melt [m/s]
pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  
                                                % geothermal heat + frictional heating derived basal melt [m/s]
pd.meltinterior = 0*((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; 
                                                % flux of basal melt up to the ~icedivide (200 km) [m2/s]

pd.hb_reg1 = 5e-3;                              % Regularisation parameter for hb in blister propagation
pd.hb_reg2 = 1e-3;                              % Regularisation parameter for hb in blister leakage
pd.N_reg1 = 1e4;                                % Regularisation parameter for N in blister leakage
pd.deltap_reg = 1e4;                            % Regularisation parameter for deltap
pd.B_reg = pd.Ye*(1000)^3/(12*(1-0.33)^2);      % Reg parameter for bending stiffness

% non-dimensionalise
ps = struct;
[ps,pp] = nevis_nondimension(pd,ps,oo);

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
gg.n1_blister = gg.n1m;

%% label boundary nodes
gg = nevis_label(gg,gg.n1m);
gg = nevis_label_blister(gg,gg.n1_blister,oo);    % label blister boundary nodes
oo.adjust_boundaries = 1;                         % enable option of changing conditions
% oo.boundary_method = 'vel_tblr';

%% add parameters and boundary labels for ice velocity
[pd,ps,pp,oo] = nevis_update_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity
gg = nevis_label_ice_test(gg, oo); % add boundary labels needed for ice velocity

if ~isfield(pp,'eps_reg'), pp.eps_reg = 0.01; end % regularisation on strain rates
if ~isfield(pp,'Ub_reg'), pp.Ub_reg = 1e-16; end % regularisation on sliding speed (max-based, matches nevis_velocity)
if ~isfield(pp,'N_slide_reg'), pp.N_slide_reg = 1e-16; end % regularisation on effective pressure (max-based, matches nevis_velocity)
if ~isfield(pp,'taud_reg'), pp.taud_reg = 1e-16; end % regularisation on basal shear stress [ may not be needed ? ]
if ~isfield(pp,'C2'), pp.C2 = 0; end % added power-law coefficient in sliding law

%% load the slipperiness field for the inversion test
inv = load(['./data/C_inversion_results.mat']);
% partition the total slipperiness coefficient C_total into two components C1 and C2 based on the specified partition ratio, dimensionalize C1 and C2 with the inversion scales. The dimensional C1 and C2 will be nondimensionalized again with the forward model scales below
[C1_hat_dim, C2_hat_dim] = nevis_inv_partition(inv.C_hat, oo.partition_ratio, inv.u_obs_noisy, inv.v_obs_noisy, inv.N_current, inv.aa, inv.pp, inv.gg, inv.oo, inv.ps);

%% initialize variables
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);         % default initialisation
aa.C = C1_hat_dim * (ps.u_b^(1/pp.n_slide) / ps.tau);
aa.C2 = C2_hat_dim * (ps.u_b^(1/pp.n_slide) / ps.tau); % added power-law coefficient in sliding law

% plot the C field
% figure;
% pcolor(gg.nx, gg.ny, reshape(aa.C, gg.nI, gg.nJ)); shading flat; colorbar();
% title('Nondimensional slipperiness field C');

load(['./data/velocity_inverted.mat'], 'vv_hydro');
% pd.k_f = 0.9;                                     % percent overburden (k-factor) 
% vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
% N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 

N = vv_hydro.N;     
vv.phi = aa.phi_0 - N;
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

% Use observed velocity
load([oo.dn '/' 'measures_for_nevis_140km']); 
dd = measures_for_nevis_140km;
un = dd.u_obs_dim/(pd.ty*ps.u);
vn = dd.v_obs_dim/(pd.ty*ps.u);

un_filled = un(:);
vn_filled = vn(:);

% find nans on es2,fs2 and interpolate with the nearest value
nan_mask_u = isnan(un_filled);
if any(nan_mask_u)
    valid_idx = find(~nan_mask_u);
    nan_idx = find(nan_mask_u);
    un_filled(nan_idx) = griddata(gg.nx(valid_idx), gg.ny(valid_idx), ...
                                   un_filled(valid_idx), gg.nx(nan_idx), gg.ny(nan_idx), 'nearest');
end
nan_mask_v = isnan(vn_filled);
if any(nan_mask_v)
    valid_idx = find(~nan_mask_v);
    nan_idx = find(nan_mask_v);
    vn_filled(nan_idx) = griddata(gg.nx(valid_idx), gg.ny(valid_idx), ...
                                   vn_filled(valid_idx), gg.nx(nan_idx), gg.ny(nan_idx), 'nearest');
end

vv.u = gg.emean2*un_filled;
vv.v = gg.fmean2*vn_filled;
vv.u(gg.eout2) = NaN;
vv.v(gg.fout2) = NaN;
aa.u_obs = vv.u;
aa.v_obs = vv.v;
% vv.u(gg.ebdy2) = 0*aa.u_obs(gg.ebdy2);
% vv.v(gg.fbdy2) = 0*aa.v_obs(gg.fbdy2);

% Option 3: used inverted velocity
vv.u = vv_hydro.u;
vv.v = vv_hydro.v;
vv.u(gg.eout2) = NaN;
vv.v(gg.fout2) = NaN;
aa.u_obs = vv.u;
aa.v_obs = vv.v;

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
pp.x_l = lakes.x_m/ps.x;                          % x-coord of lakes
pp.y_l = lakes.y_m/ps.x;                          % y-coord of lakes
pp.V_l = 0*lakes.volume_m3/(ps.Q0*ps.t);          % volume of lakes, set to zero for spinup       
pp.t_drainage = 0.5*pd.td/ps.t*(lakes.drainage_time_start + lakes.drainage_time_end);  % time of lake drainages (assumed to be the middle time of the Gaussian)
pp.t_duration = 0.5*pd.td/ps.t*(lakes.drainage_duration_days);  % duration of lake drainages
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments

%% surface runoff
load([oo.dn '/' 'runoff_2022_nevis140.mat']);                   % load data for year of interest (previously collated)
% Smooth runoff variability along time (rows = days) while preserving
% annual volume at each spatial point (columns).
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

pp.meltE = @(t) (0/1000/pd.td/ps.m)*(1-exp(-t/(20*pd.td/ps.t))); 

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
oo.t_span = [(1:1:365*1)*pd.td/ps.t];
% oo.t_span = [(1:1:59)*pd.td/ps.t (59.5:0.001:60.5)*pd.td/ps.t (61:1:120)*pd.td/ps.t];

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

