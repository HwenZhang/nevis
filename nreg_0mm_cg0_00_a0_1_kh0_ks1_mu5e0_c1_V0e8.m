%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS regional model for the outlet glacier in Stevens et al. 2018
% The surface runoff data can be either from RACMO or a prescribed function
% The surface and bed profiles are the same as in Stevens et al. 2018
% Initial condition can be either from the default condition (2008 DOY 365) or
% from a steady-state drainage system, either summertime or wintertime

% Author: Hanwen Zhang 
% Date: 2025-05
format compact

%% setup paths and directories
clc,clear
oo.root = './';                                % filename root
oo.code = '../nevis/src';                      % code directory  
oo.results = 'results';                        % path to the results folders
oo.dataset = 'nevis_regional';                 % dataset name     
oo.casename = 'nreg_0mm_cg0_00_a0_1_kh0_ks1_mu5e0_c1_V0e8';           
                                               % casename
oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
oo.dn = [oo.root, 'data/', oo.dataset, '/'];   % path to the data
addpath(oo.code);                              % add path to code
mkdir(oo.rn);                                  % create directory for results    

%% parameters
% default parameters 
[pd,oo] = nevis_defaults([],oo);  

oo.evaluate_variables = 1;
oo.use_modified_N = 0;
oo.input_gaussian = 1;
oo.relaxation_term = 0;                         % 0 is alpha hb, 1 is alpha deltap hb
oo.initial_condition = 1;                       % 0 is default condition from 0365.mat, 1 is using steady-state drainage system, wither summertime or wintertime

% leakage term
if oo.relaxation_term == 0                      % 0: exponential decay: -\alpha_0(1+h/hc+S/Sc) h_b         
    pd.alpha_b = 1.0/(10*pd.td);                % relaxation rate (s^-1)
    pd.kappa_b = 0;                             % relaxation coeff 
    pd.m_l=0;
elseif oo.relaxation_term == 1                  % 1: proportional to pressure diff and thickness: -\kappa/\mu(p_b-p_w)h_b
    pd.alpha_b = 0;                             % relaxation rate (s^-1)
    pd.kappa_b = 1e-11;                         % relaxation coeff 
    pd.m_l=1;
elseif oo.relaxation_term == 2                  % 2: channel control, enhanced at channels: -\alpha_0 (\tanh(S/S_c))
    pd.alpha_b = 1.0/(10*pd.td);                % relaxation rate (s^-1)
    pd.S_crit = 0.1;                            % critical cross section (m^2), below which there is no leakage to the drainage system
    pd.m_l=0;
end

% alter default parmaeters 
pd.mu = 5.0e0;                                  % water viscosity (Pa s)
pd.Ye = 8.8e9;                                  % Young's modulus (Pa)
pd.B = pd.Ye*(1e3)^3/(12*(1-0.33^2));           % bending stiffness (Pa m^3)
pd.c_e_reg2 = 0.00/1e3/9.81;                    % elastic sheet thickness [m/Pa]
pd.N_reg2 = 1e4; % 1e3                          % regularisation pressure for elastic sheet thickness 
pd.u_b = 100/pd.ty;                             % sliding speed [m/s]
pd.sigma = 1e-3;                                % englacial void fraction
pd.h_r = 0.1;                                   % roughness height [m]
pd.l_r = 10;                                    % roughness length [m]
pd.l_c = 1000;                                  % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 1e-3;                                  % sheet permeability constant
pd.tau_b = 60e3; %60 kPa                        % driving stress [Pa]
pd.melt = pd.G/pd.rho_w/pd.L;                   % geothermal heat derived basal melt [m/s]
pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  % geothermal heat + frictional heating derived basal melt [m/s]
pd.meltinterior = ((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; % flux of basal melt up to the ~icedivide (200 km) [m2/s]
pd.E_lapse = 40/1000/pd.td/10^3;
pd.kl_s = 1.0;                                  % leakage dependence on S
pd.kl_h = 1.0;                                  % leakage dependence on h
pd.c0 = 1.0;                                    % constant for leakage dependence on S and h, default is 1.0

% non-dimensionalise
ps = struct;
[ps,pp] = nevis_nondimension(pd,ps,oo);   

%% grid and geometry
load([oo.dn, 'morlighem_for_nevis_140km']);       % load Morlighem bedmap (previously collated)
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
% aa.phi = aa.phi_a(gg.nbdy)+k_factor*(aa.phi_0(gg.nbdy)-aa.phi_a(gg.nbdy));    % prescribed boundary pressure
% aa.phi_b = aa.phi_0;                            % prescribed boundary pressure at overburden
aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric (LAS 18 Nov. 2015)

%% moulins 
oo.density_moulins = 1;
oo.keep_all_moulins = 0;

% create new moulin locations:
% [pp.ni_m,pp.sum_m] = nevis_moulins_density_joughin([],[],gg,oo,aa,ps);

% load consistent location moulins for all runs:
temp = load([oo.dn, 'nevis_170207a.mat'], 'pp');
pp_temp = temp.pp;
pp.ni_m = pp_temp.ni_m;   % consistent locations
pp.sum_m = pp_temp.sum_m; % consistent locations

%% supraglacial lakes
 pp.x_l = [1e4/ps.x];                                            % x-coord of lakes
 pp.y_l = [-1e4/ps.x];                                           % y-coord of lakes
 pp.V_l = [0e8/(ps.Q0*ps.t)];                                    % volume of lakes         
 pp.t_drainage = [500*pd.td/ps.t];                               % time of lake drainages (assumed to be the middle time of the Gaussian)
 pp.t_duration = [0.25*pd.td/ps.t];                              % duration of lake drainages, 6hr
 [pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments 

%% parameters for timesteping
 % hourly timesteps, save timesteps, save moulin pressures
 oo.dt = 0.01/24*pd.td/ps.t; % hourly timesteps
 oo.save_timesteps = 1;      % save timesteps
 oo.save_pts_all = 1; 
 oo.pts_ni = pp.ni_l;   

%% surface input
% to include distributed runoff from e.g. RACMO, define function
% runoff(t,gg) to return runoff (m/s) at time t (s), at each point on the
% grid (ie runoff(t,gg) should return a vector of size gg.nIJ-by-1), then
% include as:
% pp.runoff_function = @(t) runoff(ps.t*t,gg)/ps.m; oo.runoff_function = 1;

% RACMO distributed input
oo.surface_runoff = 1;                          % If set to 1 turns on RACMO moulin input (m3/sec)
                                                % If set to 0, moulin input is collected from surface runoff
oo.RACMO_runoff = 0;                            % If set to 1, use RACMO input as below; 
                                                % If set to 0, use prescribed runoff function 
oo.distributed_input = 0;                       % If set to 1 turns on distributed input
                                                % If set to 0 , input is collected into moulins 

load([oo.dn, 'runoff_2009_nevis140.mat']);      % load data for year of interest (previously collated)
pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2009_nevis140)./ps.m;  % RACMO input (m/sec)
pp.meltE = @(t) (0/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));
pp.input_function = @(t) runoff_moulins(((t*ps.t)/pd.td),runoff_2009_nevis140,pp.sum_m,gg.Dx(1))./ps.m; % RACMO moulin input (m3/sec)

%% Timestep 
if oo.initial_condition
    % load initial condition from 2008 DOY 365 (default)
    load(['./initial_condition/', 'nevis_init_default0365.mat'],'vv')
else
    % use steady-state drainage system, either summertime or wintertime
    load(['./initial_condition/', 'nevis_init_condition_summer_alpha0_01.mat'],'vv')
end
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.pts_ni = pp.ni_l;                            % save lake pressures
t_span = (1:1:1000)*pd.td/ps.t;               % time span for simulation (in ps.t)

%% save initial parameters
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','vv','oo');
[tt,vv,info] = nevis_timesteps(t_span,vv,aa,pp,gg,oo);

%% expand/update variables
aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','oo','tt');

%% Simple animate
nevis_regional_animation