
%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS 1-dimensional model for an idealised ice sheet
% The surface runoff is set up by a prescribed function
% The surface and bed profiles are the same as in Hewitt 2013

% Author: Auto-generated script  
% Date: 2025-09
% Parameters: oo.mean_perms=1, pd.H_reg_value=pd.Ye*(1000)^3/(12*(1-0.33)^2)
format compact

%% read in the initial condition
oo.root = './';                                % filename root
oo.code = '../nevis/src';                      % code directory   
oo.results = 'results';                        % path to the results folders
oo.dataset = 'nevis_regional';                 % dataset name
oo.casename = 'n2d_region_meanperms1_Hreg1000_kappa1e_09_mu1e3_spinup';    % casename
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
oo.display_residual =0;
pd.alpha_b = 0;                                 % relaxation rate (s^-1)
pd.kappa_b = 1e-09;                             % relaxation coeff

% alter default parmaeters 
pd.mu = 1000.0;                                    % water viscosity (Pa s)
pd.c_e_reg2 = 0.00/1e3/9.81;                    % elastic sheet thickness [m/Pa]
pd.u_b = 100/pd.ty;                             % sliding speed [m/s]
pd.sigma = 1e-3;                                % englacial void fraction
pd.h_r = 0.1;                                   % roughness height [m]
pd.l_r = 10;                                    % roughness length [m]
pd.l_c = 1000;                                  % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 1e-3;                                  % sheet permeability constant
pd.tau_b = 60e3;                                % driving stress [Pa]
pd.melt = pd.G/pd.rho_w/pd.L;                   % geothermal heat derived basal melt [m/s]
pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  
                                                % geothermal heat + frictional heating derived basal melt [m/s]
pd.meltinterior = ((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; 
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
gg = nevis_grid(dd.X_km(1:dd.skip:end,1)/ps.x,dd.Y_km(1,1:dd.skip:end)/ps.x,oo); 
b = reshape(dd.B_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
s = reshape(dd.S_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);

x = dd.X_km(1:dd.skip:end,1)/ps.x;
y = dd.Y_km(1,1:dd.skip:end)/ps.x;
 
%% mask with minimum ice thickness
H = max(s-b,0);
Hmin = 0/ps.z; 
nout = find(H<=Hmin);
gg = nevis_mask(gg,nout); 
gg.n1m = gg.n1;                                   % label all edge nodes as boundary nodes for pressure

%% label boundary nodes
gg = nevis_label(gg,gg.n1m);
oo.adjust_boundaries = 1;                         % enable option of changing conditions

%% plot grid
% nevis_plot_grid(gg); return;                    % check to see what grid looks like

%% initialize variables
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);         % default initialisation
pd.k_f = 0.9;                                     % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

%% boundary conditions
aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric

%% moulins 
oo.density_moulins = 1;
oo.keep_all_moulins = 0;
% load moulin locations
tmp = load([oo.dn '/' '/nevis_170207a.mat'],'pp');
pp.ni_m = tmp.pp.ni_m; pp.sum_m = tmp.pp.sum_m; clear tmp;

%% supraglacial lakes
 pp.x_l = [0.5*(x(1)+x(end))];                                   % x-coord of lakes
 pp.y_l = [0.5*(y(1)+y(end))];                                   % y-coord of lakes
 pp.V_l = [0e8/(ps.Q0*ps.t)];                                    % volume of lakes         
 pp.t_drainage = [300*pd.td/ps.t];                               % time of lake drainages (assumed to be the middle time of the Gaussian)
 pp.t_duration = [0.25*pd.td/ps.t];                              % duration of lake drainages, 6hr
 [pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments  

%% surface runoff
load([oo.dn '/' 'runoff_2009_nevis140.mat']);               % load data for year of interest (previously collated)
% RACMO distributed input
oo.surface_runoff = 1;                                           % If set to 1 turns on surface runoff input
                                                                 % 0 is prescribe moulin input with a function
oo.RACMO_runoff = 1;                                             % If set to 1 turns on RACMO surface runoff data
                                                                 % 0 is a simple sinusoidal function
oo.distributed_input = 0;                                        % If set to 1 distributes input across the domain
                                                                 % 0 inputs to discrete moulins

% pp.meltE = @(t) (60/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 
pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2009_nevis140)./ps.m;  % distributed input (m/sec)
oo.pts_ni = [pp.ni_m; pp.ni_l];                                   % points to save moulin and lake pressures
% pp.input_function = @(t) runoff_moulins(((t*ps.t)/pd.td),runoff_2009_nevis140,pp.sum_m,gg.Dx(1))./ps.m; % RACMO moulin input (m3/sec)

%% timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.t_span = (1:1:365)*pd.td/ps.t; % 1 year of timesteps

% Add GPS station points downstream of the moulin every 5km
% pp.ni_gps = nevis_gps_array((pp.x_l:5e3/ps.x:L/ps.x), pp.y_l, gg, oo); % GPS station points
% oo.pts_ni = [pp.ni_l pp.ni_m pp.ni_gps];    
% oo.t_span = vv.t + [(1:1:19)*pd.td/ps.t (19.9:0.001:20.1)*pd.td/ps.t (21:1:1.5*365)*pd.td/ps.t];

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
