%% Run file for Hewitt (2013) nevis model as used in Stevens et al., JGR, 2018.
% 140km-wide region with RACMO runoff input starting on DOY 1 2009.
% Starting condition is saved values from 2008 DOY 365 (2008 = a spin up year). 
format compact
clc,clear
oo.root = './';                                % filename root
oo.code = '../nevis';                          % code directory  
oo.results = 'results';                        % path to the results folders
% oo.casename = 'nevis_2009_140km_alpha0_2';         % casename
oo.casename = 'test_2009_140km_mu1e_3_kappa0';     % casename

                                               % casename
oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
addpath(oo.code);                              % add path to code
mkdir(oo.rn);                                  % create directory for results    

%% parameters
% default parameters
[pd,oo] = nevis_defaults([],oo);  

oo.evaluate_variables = 1;
oo.use_modified_N = 0;
oo.input_gaussian = 1;
oo.input_function = 0;                         % If set to 1 turns on RACMO moulin input (m3/sec)
oo.relaxation_term = 1;                        % 0 is alpha hb, 1 is alpha deltap hb

pd.alpha_b = 1.0/(5*pd.td);                    % relaxation rate (s^-1)
pd.mu = 1.0e-3;                                % water viscosity (Pa s)
pd.Ye = 8.8e9;                                 % Young's modulus (Pa)
pd.B = pd.Ye*(1e3)^3/(12*(1-0.33^2));          % bending stiffness (Pa m^3)
if oo.relaxation_term == 0
    pd.alpha_b = 1.0/(5*pd.td);                % relaxation rate (s^-1)
    pd.kappa_b = 0;                            % relaxation coeff 
elseif oo.relaxation_term == 1
    pd.alpha_b = 0;                            % relaxation rate (s^-1)
    pd.kappa_b = 0;                         % relaxation coeff 
end

% alter default parmaeters 
pd.c_e_reg2 = 0.00/1e3/9.81;        % elastic sheet thickness [m/Pa]
pd.N_reg2 = 1e4; % 1e3              % regularisation pressure for elastic sheet thickness 
pd.u_b = 100/pd.ty;                 % sliding speed [m/s]
pd.sigma = 1e-3;                    % englacial void fraction
pd.h_r = 0.1;                       % roughness height [m]
pd.l_r = 10;                        % roughness length [m]
pd.l_c = 1000;                      % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 1e-3;                      % sheet permeability constant
pd.tau_b = 60e3; %60 kPa            % driving stress [Pa]
pd.melt = pd.G/pd.rho_w/pd.L;       % geothermal heat derived basal melt [m/s]
pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  % geothermal heat + frictional heating derived basal melt [m/s]
pd.meltinterior = ((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; % flux of basal melt up to the ~icedivide (200 km) [m2/s]

% non-dimensionalise
ps = struct;
[ps,pp] = nevis_nondimension(pd,ps,oo);   

%% grid and geometry
load('data/morlighem_for_nevis_140km'); % load Morlighem bedmap (previously collated)
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
gg.n1m = gg.n1;                     % label all edge nodes as boundary nodes for pressure

%% label boundary nodes
gg = nevis_label(gg,gg.n1m);
oo.adjust_boundaries = 1;           % enable option of changing conditions

%% plot grid
% nevis_plot_grid(gg); return;  % check to see what grid looks like

%% initialize variables
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);       % default initialisation
pd.k_f = 0.9;                                     % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);   % initial pressure  k_f*phi_0
N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

%% boundary conditions
% aa.phi = aa.phi_a(gg.nbdy)+k_factor*(aa.phi_0(gg.nbdy)-aa.phi_a(gg.nbdy));    % prescribed boundary pressure
% aa.phi_b = aa.phi_0;             % prescribed boundary pressure at overburden
aa.phi_b = max(aa.phi_0,aa.phi_a); % prescribed boundary pressure at overburden or atmospheric (LAS 18 Nov. 2015)

%% moulins 
oo.random_moulins = 0;
% create new moulin locations:
% [pp.ni_m,pp.sum_m] = nevis_moulins_density_joughin([],[],gg,oo,aa,ps);
% load consistent location moulins for all runs:
% load('nevis/nevis_170207a.mat','pp');
%% lakes
 % a single-point lake
 pp.x_l = [0.5*(x(1)+x(end))];                                   % x-coord of lakes
 pp.y_l = [0.5*(y(1)+y(end))];                                   % y-coord of lakes
 pp.V_l = [1e8/(ps.Q0*ps.t)];                                    % volume of lakes         
 pp.t_drainage = [300*pd.td/ps.t];                               % time of lake drainages (assumed to be the middle time of the Gaussian)
 pp.t_duration = [0.25*pd.td/ps.t];                              % duration of lake drainages, 6hr

 [pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo); % calculate lake catchments 

%% parameters for timesteping
 % hourly timesteps, save timesteps, save moulin pressures
 oo.dt = 0.01/24*pd.td/ps.t; % hourly timesteps
 oo.save_timesteps = 1; % save timesteps
 oo.save_pts_all = 1; 
 oo.pts_ni = pp.ni_l;   

%% surface input
% % to include distributed runoff from e.g. RACMO, define function
% % runoff(t,gg) to return runoff (m/s) at time t (s), at each point on the
% % grid (ie runoff(t,gg) should return a vector of size gg.nIJ-by-1), then
% % include as:
% pp.runoff_function = @(t) runoff(ps.t*t,gg)/ps.m; oo.runoff_function = 1;
% load('nevis/runoff_2009_nevis140.mat');       % load data for year of interest (previously collated)
% RACMO distributed input
oo.distributed_input = 1;                     % If set to 1 turns on RACMO distributed input
pp.meltE = @(t) (60/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 
% pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2009_nevis140)./ps.m;  % distributed input (m/sec)

% RACMO moulin input
% pp.input_function = @(t) runoff_moulins(((t*ps.t)/pd.td),runoff_2009_nevis140,pp.sum_m,gg.Dx(1))./ps.m; % RACMO moulin input (m3/sec)

% %% Timesteps and saving model output
% oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = pp.ni_m;

%% save initial parameters
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
% load('nevis/nevis_22221/0365.mat','vv','tt')
% tspan = (1:1:365)*pd.td/ps.t; % 1 year of timesteps
t_span = (0:2000)*(0.2*pd.td/ps.t);
[tt,vv,info] = nevis_timesteps(t_span,vv,aa,pp,gg,oo);

%% expand/update variables
aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge

save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','oo','tt');

%% Simple animate
% nevis_regional_animation