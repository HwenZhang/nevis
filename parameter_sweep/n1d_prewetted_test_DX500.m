%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS 1-dimensional model for an idealised ice sheet
% No surface runoff, Gaussian influx at early times, then constant volume
% Linear bed slope, constant ice thickness

% Author: Hanwen Zhang  
% Date: 2025-06-16
format compact

%% setup paths and directories
clear,clc
oo.root = './';                                % filename root
oo.code = '../nevis/src';                      % code directory  
oo.analysis = '../nevis/analysis';             % analysis code directory  
oo.results = 'results';                        % path to the results folders
oo.dataset = 'nevis_regional';                 % dataset name     
oo.casename = 'n1d_prewetted_test_DX500'; % casename
oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
oo.dn = [oo.root, 'data/', oo.dataset, '/'];   % path to the data
addpath(oo.code);                              % add path to code
addpath(oo.analysis);                          % add path to analysis
mkdir(oo.rn);                                  % create directory for results

dt = 1e-3;
Nx = 201;
%% parameters
% default parameters 
[pd,oo] = nevis_defaults([],oo);  

oo.evaluate_variables = 1;
oo.input_gaussian = 0;
oo.relaxation_term = 1;                         % 0 is alpha hb, 1 is alpha deltap hb
oo.initial_condition = 1;                       % 1 is default condition from 0365.mat, 0 is using steady-state drainage system, wither summertime or wintertime
oo.mean_perms = 1;
oo.modified_mean_perms = 0;

% leakage term
if oo.relaxation_term == 0                      % 0: exponential decay: -\alpha_0(1+h/hc+S/Sc) h_b         
    pd.alpha_b = 0.0/(100*pd.td);                % relaxation rate (s^-1)
    pd.kappa_b = 0;                             % relaxation coeff 
    pd.m_l=0;
elseif oo.relaxation_term == 1                  % 1: proportional to pressure diff and thickness: -\kappa/\mu(p_b-p_w)h_b
    pd.alpha_b = 0;                             % relaxation rate (s^-1)
    pd.kappa_b = 1e-60;                         % relaxation coeff 
    pd.alpha_dh = 0;                            % alpha for the change in hb, used in the relaxation term [1/(m s)]
    pd.alpha_ds = 0;                            % alpha for the change in S, used in the relaxation term [1/(m s)]
    pd.m_l=1;
elseif oo.relaxation_term == 2                  % 2: channel control, enhanced at channels: -\alpha_0 (\tanh(S/S_c))
    pd.alpha_b = 1.0/(10*pd.td);                % relaxation rate (s^-1)
    pd.S_crit = 0.1;                            % critical cross section (m^2), below which there is no leakage to the drainage system
    pd.m_l=0;
end

% alter default parameters
pd.mu = 5.0e0;                                  % water viscosity (Pa s)
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
pd.kl_s = 0.0;                                  % leakage dependence on S
pd.kl_h = 0.0;                                  % leakage dependence on h
pd.c0 = 1.0;                                    % constant for leakage dependence on S and h, default is 1.0
pd.hb_reg1 = 1e-3;                              % regularizing thickness of blister [m]
pd.hb_reg2 = 1e-10;                              % regularizing thickness of blister [m]

% non-dimensionalise
ps = struct;
[ps,pp] = nevis_nondimension(pd,ps,oo);   

%% grid and geometry
L = 1e5;                               % length of the domain [m]
x = linspace(0,(L/ps.x),Nx); 
y = linspace(0,(L/ps.x),1);            % 1-d grid of length 50km 
oo.yperiodic = 1;                      % oo.yperiodic = 1 necessary for a 1-d grid
oo.xperiodic = 0;
gg = nevis_grid(x,y,oo); 
b = 5e4/ps.z - 0.0*gg.nx*ps.x/ps.z;   % bed slope angle = 0.1
s = 1000/ps.z + b;

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
pd.k_f = 0.9;                                     % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)
% vv.hb = 0*ones(size(vv.hs));                        % initial blister sheet thickness
% prescribed blister sheet thickness
xb = 0.5*L;                                       % x-coordinate of the blister center (m)
yb = 0;                                           % y-coordinate of the blister center (m)
R0 = 5e3;                                         % initial blister radius (m)
V0 = 1e7;                                         % initial blister volume (m^3)
radius = abs(gg.nx*ps.x - xb);                    % distance from the blister center (m)
% vv.hb = 5e-3*ones(size(radius))/ps.hb;
% vv.hb = 3*V0/(pi*R0^2)/ps.hb*(radius<=R0).*(1-radius.^2/R0^2).^2;   % initial blister sheet thickness
% vv.pb = 192*pd.B*V0*ones(size(vv.hb))/(pi*R0^6)/ps.phi;             % initial blister pressure

%% boundary conditions
% aa.phi = aa.phi_a(gg.nbdy)+k_factor*(aa.phi_0(gg.nbdy)-aa.phi_a(gg.nbdy));    % prescribed boundary pressure
% aa.phi_b = aa.phi_0;                            % prescribed boundary pressure at overburden
aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric (LAS 18 Nov. 2015)

%% moulins 
oo.keep_all_moulins = 0;
oo.random_moulins = 1;         
[pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);  % no moulin at the lake location

%% supraglacial lakes
pp.x_l = [0.5*L/ps.x];                                          % x-coord of lakes
pp.y_l = [0];                                                   % y-coord of lakes
pp.V_l = [1e7/(ps.Q0*ps.t)];                                    % volume of lakes         
pp.t_drainage = [0.5*pd.td/ps.t];                              % time of lake drainages (assumed to be the middle time of the Gaussian)
pp.t_duration = [1.0*pd.td/ps.t];                              % duration of lake drainages, 6hr
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments 

%% surface input
% RACMO distributed input
oo.surface_runoff = 0;                          % If set to 1 turns on surface runoff
                                                % If set to 0, moulin input is prescribed with a function
oo.RACMO_runoff = 0;                            % If set to 1, use RACMO input as below; 
                                                % If set to 0, use prescribed runoff function 
oo.distributed_input = 0;                       % If set to 1 turns on distributed input
                                                % If set to 0 , input is collected into moulins 
pp.input_function = @(t) 0./ps.m;               % RACMO moulin input (m3/sec)

%% Timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.pts_ni = pp.ni_l;                             % save lake pressures
% oo.t_span = [(1:1:44)*pd.td/ps.t (45:0.01:55)*pd.td/ps.t (56:1:2000)*pd.td/ps.t]; % 
oo.t_span = (dt:dt:1.0)*pd.td/ps.t;          % time span for simulation (in ps.t)

%% save initial parameters
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','vv','oo');
[tt,vv,info] = nevis_timesteps(oo.t_span,vv,aa,pp,gg,oo);

%% expand/update variables
aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','oo','tt');

%% Simple animate
% nevis_downslope_power_law