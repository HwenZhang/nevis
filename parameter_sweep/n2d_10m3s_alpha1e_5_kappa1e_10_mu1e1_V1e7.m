%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS 1-dimensional model for an idealised ice sheet
% The surface runoff is set up by a prescribed function
% The surface and bed profiles are the same as in Hewitt 2013

% Author: Hanwen Zhang  
% Date: 2025-05
format compact

%% setup paths and directories
clear,clc
oo.root = './';                                % filename root
oo.code = '../nevis/src';                      % code directory  
oo.results = 'results';                        % path to the results folders
oo.dataset = 'nevis_regional';                 % dataset name     
oo.casename = 'n2d_10m3s_alpha1e_5_kappa1e_10_mu1e1_V1e7';   
           
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
oo.input_gaussian = 1;
oo.relaxation_term = 1;                         % 0 is alpha hb, 1 is alpha deltap hb
oo.initial_condition = 1;                       % 1 is default condition from 0365.mat, 0 is using steady-state drainage system, wither summertime or wintertime
oo.mean_perm = 1;

% leakage term
if oo.relaxation_term == 0                      % 0: exponential decay: -\alpha_0(1+h/hc+S/Sc) h_b         
    pd.alpha_b = 1.0/(10*pd.td);                % relaxation rate (s^-1)
    pd.kappa_b = 0;                             % relaxation coeff 
    pd.m_l = 1;
elseif oo.relaxation_term == 1                  % 1: proportional to pressure diff and thickness: -\kappa/\mu(p_b-p_w)h_b
    pd.alpha_b = 1.0/(1e5*pd.td);               % relaxation rate (s^-1)
    pd.kappa_b = 1e-10;                         % relaxation coeff
    pd.m_l = 1;
elseif oo.relaxation_term == 2                  % 2: channel control, enhanced at channels: -\alpha_0 (\tanh(S/S_c))
    pd.alpha_b = 1.0/(10*pd.td);                % relaxation rate (s^-1)
    pd.S_crit = 0.1;                            % critical cross section (m^2), below which there is no leakage to the drainage system
    pd.m_l=0;
elseif oo.relaxation_term == 3                  % pressure diff control: k[hb]/mu*\Delta p       
    pd.m_l = 0;                                 % no dependence on hb
    pd.kappa_b = 1e-9;
    pd.alpha_b = 0;
end

% alter default parmaeters 
runoff_max = 30;                                % prescribed runoff (mm/day)
pd.mu = 1e1;                                    % water viscosity (Pa s)
pd.Ye = 8.8e9;                                  % Young's modulus (Pa)
pd.B = pd.Ye*(1e3)^3/(12*(1-0.33^2));           % bending stiffness (Pa m^3)
pd.E_lapse = 30/1000/pd.td/10^3;

pd.hb_reg1 = 0;                                 % Regularisation parameter for hb, leakage only occurs when hb >> hb_reg1
pd.hb_reg2 = 1e-6;                              % Regularisation parameter for hb, leakage only occurs when hb >> hb_reg2 [abandoned]
pd.N_reg1 = 1e3;                                % Regularisation parameter for N, (N >> Nreg, input to drainage system; N << -Nreg, input to blister))
pd.alpha_dh = 1.1574e-6;                        % alpha for the change in hb, used in the relaxation term [1/(m s)] 1.1574e-5 ~ 1e0
pd.alpha_ds = 1.5595e-5;                        % alpha for the change in hb, used in the relaxation term [1/(m^2 s)]
% pd.alpha_dh = pd.k_s*ps.h^(pd.alpha_s-1)*ps.phi^pd.beta_s/ps.hb/ps.x^(pd.beta_s+1);
% pd.alpha_ds = 1e2*pd.k_c*ps.S^(pd.alpha_c-1)*ps.phi^(pd.beta_c)/ps.hb/ps.x^(pd.beta_c+2);
pd.alpha_dh = 0;               
pd.alpha_ds = 0; 
% non-dimensionalise
ps = struct;
[ps,pp] = nevis_nondimension(pd,ps,oo);   

%% grid and geometry
L = 5e4;                                     % length of the domain [m]
W = 0.2*L;                                   % width of the domain [m]
x = linspace(0,(L/ps.x),101); 
y = linspace(0,(W/ps.x),20);        
oo.yperiodic = 1;                        % oo.yperiodic = 1 necessary for a 1-d grid
oo.xperiodic = 0;
gg = nevis_grid(x,y,oo); 
b = (0/ps.z)*gg.nx;                      % flat bed
s = (1060/ps.z)*(1-ps.x*gg.nx/L).^0.5;   % ice surface topography 

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
pd.k_f = 0.99;                                    % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

%% boundary conditions
% aa.phi = aa.phi_a(gg.nbdy)+k_factor*(aa.phi_0(gg.nbdy)-aa.phi_a(gg.nbdy));    % prescribed boundary pressure
% aa.phi_b = aa.phi_0;                            % prescribed boundary pressure at overburden
aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric (LAS 18 Nov. 2015)

%% moulins 
oo.keep_all_moulins = 0;
oo.random_moulins = 0;         
[pp.ni_m,pp.sum_m] = nevis_moulins([0.25*L/ps.x],[0.5*W/ps.x],gg,oo);     % one moulin at the lake location

%% supraglacial lakes
pp.x_l = [0.5*L/ps.x];                                          % x-coord of lakes
pp.y_l = [0.5*W/ps.x];                                          % y-coord of lakes
pp.V_l = [1e7/(ps.Q0*ps.t)];                                    % volume of lakes         
pp.t_drainage = [3.0*365*pd.td/ps.t];                           % time of lake drainages (assumed to be the middle time of the Gaussian)
pp.t_duration = [0.25*pd.td/ps.t];                              % duration of lake drainages, 6hr
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments 

%% surface input
% RACMO distributed input
oo.surface_runoff = 0;                          % If set to 1 turns on surface runoff
                                                % If set to 0, moulin input is prescribed with a function
oo.RACMO_runoff = 0;                            % If set to 1, use RACMO input as below; 
                                                % If set to 0, use prescribed runoff function 
oo.distributed_input = 0;                       % If set to 1 turns on distributed input
                                                % If set to 0 , input is collected into moulins 

% load([oo.dn, 'runoff_2009_nevis140.mat']);      % load data for year of interest (previously collated)
% pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2009_nevis140)./ps.m;  % RACMO input (m/sec)
pp.meltE = @(t) (runoff_max/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));
pp.input_function = @(t) 10*(1-exp(-t/(300*pd.td/ps.t)))./(ps.m*ps.x^2);     % RACMO moulin input (100 m3/sec)

%% Timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.pts_ni = [pp.ni_l pp.ni_m];                      % save lake pressures
oo.t_span = [(1:1:3*365-6)*pd.td/ps.t (3*365-5:0.1:3*365+10)*pd.td/ps.t (3*365+11:1:4*365)*pd.td/ps.t];              % time span for simulation (in ps.t)
% oo.t_span = [(1:1:4*365)*pd.td/ps.t];              % time span for simulation (in ps.t)

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
% nevis_regional_animation