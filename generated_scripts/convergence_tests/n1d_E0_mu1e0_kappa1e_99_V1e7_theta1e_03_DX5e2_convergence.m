
%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS 1-dimensional model for an idealised ice sheet
% The surface runoff is set to 0, and the bedrock is linear with a slope of 0.01

% Author: Hanwen Zhang  
% Date: 2025-09
format compact

%% setup paths and directories
clear,clc
oo.root = './';                                % filename root
oo.code = '../nevis/src';                      % code directory   
oo.results = 'results';                        % path to the results folders
oo.dataset = 'nevis_regional';                 % dataset name     
oo.casename = 'n1d_E0_mu1e0_kappa1e_99_V1e7_theta1e_03_DX5e2_convergence';                        % case name
                                            % casename
oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
oo.dn = [oo.root, 'data/', oo.dataset, '/'];   % path to the data
addpath(oo.code);                              % add path to code
mkdir(oo.rn);                                  % create directory for results 
oo.Tol_F = 1e-7;                               % tolerance for fsolve

%% parameters
% default parameters 
[pd,oo] = nevis_defaults([],oo);  

oo.evaluate_variables = 1;
oo.input_gaussian = 1;
oo.relaxation_term = 1;                         % 0 is alpha hb, 1 is alpha deltap hb
oo.initial_condition = 1;                       % 1 is default condition from 0365.mat, 0 is using steady-state drainage system, winter or summertime
oo.mean_perms = 1;
oo.modified_mean_perms = 0;
oo.display_residual =0;
pd.alpha_b = 0;                                 % relaxation rate (s^-1)
pd.kappa_b = 1e-99;                             % relaxation coeff

% alter default parameters
pd.mu = 1.0;                                  % water viscosity (Pa s)
pd.Ye = 8.8e9;                                  % Young's modulus (Pa)
pd.B = pd.Ye*(1e3)^3/(12*(1-0.33^2));           % bending stiffness (Pa m^3)
pd.c_e_reg2 = 0.00/1e3/9.81;                    % elastic sheet thickness [m/Pa]
pd.N_reg1 = 1e3;                                % effective pressure cut off for regularization of moulin input term [Pa]
pd.N_reg2 = 1e4; % 1e3                          % regularisation pressure for elastic sheet thickness  [Pa]
pd.u_b = 100/pd.ty;                             % sliding speed [m/s]
pd.sigma = 0e-3;                                % englacial void fraction
pd.h_r = 0.1;                                   % roughness height [m]
pd.l_r = 10;                                    % roughness length [m]
pd.l_c = 1000;                                  % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 1e-3;                                  % sheet permeability constant
pd.tau_b = 60e3; %60 kPa                        % driving stress [Pa]
pd.melt = pd.G/pd.rho_w/pd.L;                   % geothermal heat derived basal melt [m/s]
pd.melt = 0*(pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  % geothermal heat + frictional heating derived basal melt [m/s]
pd.meltinterior = ((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; % flux of basal melt up to the ~icedivide (200 km) [m2/s]
pd.E_lapse = 40/1000/pd.td/10^3;
pd.hb_reg1 = 1e-3;                              % regularizing thickness of blister [m]
pd.hb_reg2 = 1e-3;                              % regularizing thickness of blister [m]
pd.B_reg = 0;

% non-dimensionalise
ps = struct;
[ps,pp] = nevis_nondimension(pd,ps,oo); 

%% grid and geometry
L = 1e5;                               % length of the domain [m]
Dx = 500;
Nx = (L/Dx)+1;
x = linspace(0,(L/ps.x),Nx); 
y = linspace(0,(L/ps.x),1);            % 1-d grid of length 50km 
oo.yperiodic = 1;                      % oo.yperiodic = 1 necessary for a 1-d grid
oo.xperiodic = 0;
gg = nevis_grid(x,y,oo); 
b = 5e4/ps.z - 0.001*gg.nx*ps.x/ps.z;   % bed slope angle = 0.1
s = 200/ps.z + b;

%% mask with minimum ice thickness
H = max(s-b,0);
Hmin = 0/ps.z; 
nout = find(H<=Hmin);
gg = nevis_mask(gg,nout); 
gg.n1m = gg.n1;                                   % label all edge nodes as boundary nodes for pressure

%% label boundary nodes
gg = nevis_label(gg,gg.n1m);
oo.adjust_boundaries = 1;                         % enable option of changing conditions

%% initialize variables
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);         % default initialisation
pd.k_f = 0.9;                                    % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

%% boundary conditions
aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric (LAS 18 Nov. 2015)

%% moulins 
oo.keep_all_moulins = 0;
oo.random_moulins = 0;         
[pp.ni_m,pp.sum_m] = nevis_moulins([0.25*L/ps.x],[0/ps.x],gg,oo);     % one moulin at the lake location

%% supraglacial lakes
pp.x_l = [0.5*L/ps.x];                                          % x-coord of lakes
pp.y_l = [0];                                                   % y-coord of lakes
pp.V_l = [10000000.0/(ps.Q0*ps.t)];                                    % volume of lakes         
pp.t_drainage = [5*pd.td/ps.t];                                % time of lake drainages
pp.t_duration = [0.025*pd.td/ps.t];                             % duration of lake drainages
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments 

%% surface input
oo.surface_runoff = 0;                          % If set to 1 turns on surface runoff
                                                % If set to 0, moulin input is prescribed with a function
oo.RACMO_runoff = 0;                            % If set to 1, use RACMO input as below; 
                                                % If set to 0, use prescribed runoff function 
oo.distributed_input = 0;                       % If set to 1 turns on distributed input
                                                % If set to 0 , input is collected into moulins 
pp.input_function = @(t) 0*(1-exp(-t/(30*pd.td/ps.t)))./(ps.m*ps.x^2);   % prescribed moulin input (m3/sec)

%% Timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.pts_ni = pp.ni_l;                             
oo.t_span = [(0.1:0.1:9.4)*pd.td/ps.t (9.5:0.005:10.5)*pd.td/ps.t (10.6:0.1:100)*pd.td/ps.t];
oo.t_span = [(0:0.005:20)*pd.td/ps.t];

%% save initial parameters
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','vv','oo');
[tt,vv,info] = nevis_timesteps(oo.t_span,vv,aa,pp,gg,oo);

%% expand/update variables
aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); 
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); 
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','oo','tt');
        