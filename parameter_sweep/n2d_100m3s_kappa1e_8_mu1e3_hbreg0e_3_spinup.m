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
oo.casename = 'n2d_100m3s_kappa1e_8_mu1e3_hbreg0e_3_spinup';                    % casename
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
oo.initial_condition = 1;                       % 1 is default condition from 0365.mat, 0 is using steady-state drainage system, winter or summertime
oo.mean_perm = 1;
oo.display_residual =0;
% leakage term
if oo.relaxation_term == 0                      % 0: exponential decay: -\alpha_0(1+h/hc+S/Sc) h_b         
    pd.alpha_b = 1.0/(10*pd.td);                % relaxation rate (s^-1)
    pd.kappa_b = 0;                             % relaxation coeff 
    pd.c0 = 1;
    pd.m_l = 1;
elseif oo.relaxation_term == 1                  % 1: proportional to pressure diff and thickness: -\kappa/\mu(p_b-p_w)h_b
    pd.alpha_b = 0;                             % relaxation rate (s^-1)
    pd.kappa_b = 1e-08;                      % relaxation coeff
    pd.c0 = 1;
    pd.m_l = 1;
end

% alter default parmaeters 
runoff_max = 30;                                % prescribed runoff (mm/day)
moulin_input = 100;                  % prescribed moulin input (m^3/s)
pd.mu = 1000.0;                                  % water viscosity (Pa s)
pd.Ye = 8.8e9;                                  % Young's modulus (Pa)
pd.B = pd.Ye*(1e3)^3/(12*(1-0.33^2));           % bending stiffness (Pa m^3)
pd.E_lapse = 30/1000/pd.td/10^3;

pd.hb_reg1 = 0e-3;                              % Regularisation parameter for hb
pd.hb_reg2 = 1e-3;                              % Regularisation parameter for hb
pd.N_reg1 = 1e3;                                % Regularisation parameter for N

% non-dimensionalise
ps = struct;
[ps,pp] = nevis_nondimension(pd,ps,oo);   

%% grid and geometry
L = 5e4;                                     % length of the domain [m]
W = 0.4*L;                                   % width of the domain [m]
x = linspace(0,(L/ps.x),101); 
y = linspace(0,(W/ps.x),40);        
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
pd.k_f = 0.90;                                    % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

%% boundary conditions
aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric

%% moulins 
oo.keep_all_moulins = 0;
oo.random_moulins = 0;         
[pp.ni_m,pp.sum_m] = nevis_moulins([0.25*L/ps.x],[0.5*W/ps.x],gg,oo);     % one moulin at the lake location

%% supraglacial lakes
pp.x_l = [0.5*L/ps.x];                                          % x-coord of lakes
pp.y_l = [0.5*W/ps.x];                                          % y-coord of lakes
pp.V_l = [0e7/(ps.Q0*ps.t)];                                    % volume of lakes         
pp.t_drainage = [0.5*365*pd.td/ps.t];                           % time of lake drainages
pp.t_duration = [0.25*pd.td/ps.t];                              % duration of lake drainages
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments 

%% surface input
oo.surface_runoff = 0;                          
oo.RACMO_runoff = 0;                            
oo.distributed_input = 0;                       
pp.meltE = @(t) (runoff_max/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));
pp.input_function = @(t) moulin_input*(1-exp(-t/(30*pd.td/ps.t)))./(ps.m*ps.x^2);     

%% Timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.pts_ni = [pp.ni_l pp.ni_m];                      
oo.t_span = (1:1:2*365)*pd.td/ps.t;                 

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
