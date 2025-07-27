%% Script to run NEVIS regional model 
% This script is designed to run the NEVIS 1-dimensional model for an idealised ice sheet
% The surface runoff is set up by a prescribed function
% The surface and bed profiles are the same as in Hewitt 2013

% Author: Hanwen Zhang  
% Date: 2025-05
format compact

%% read in the initial condition
casename = 'n2d_0m3s_kappa1e_8_mu1e1_hbreg5e_3_V1e8_drainage'; % drainage system filename
initname = strrep(casename, '_V1e8_drainage', '_spinup'); % initial condition filename

data = load(['./results/' initname '/' initname]);
pd = data.pd;                                % load parameters from the initial condition
ps = data.ps;                                % load state variables from the initial condition
pp = data.pp;                                % load scaled parameters from the initial condition
aa = data.aa;                                % load state variables from the initial condition
oo = data.oo;                                % load options from the initial condition

oo.casename = casename;                      % drainage system filename
oo.initname = initname;                      % initial condition filename, for spinup

oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
oo.dn = [oo.root, 'data/', oo.dataset, '/'];   % path to the data
addpath(oo.code);                              % add path to code
mkdir(oo.rn);                                  % create directory for results    
pp.c0 = 0;

%% grid and geometry
L = 5e4;                                     % length of the domain [m]
W = 0.4*L;                                   % width of the domain [m]
x = linspace(0,(L/ps.x),101); 
y = linspace(0,(W/ps.x),40);        
oo.yperiodic = 1;                            % oo.yperiodic = 1 necessary for a 1-d grid
oo.xperiodic = 0;
gg = nevis_grid(x,y,oo); 
b = (0/ps.z)*gg.nx;                          % flat bed
s = (1060/ps.z)*(1-ps.x*gg.nx/L).^0.5;       % ice surface topography 

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
init_cond = load(['./results/' oo.initname '/' '0730.mat']); % load initial condition
vv = init_cond.vv;                                % load state variables from the initial condition

%% supraglacial lakes
pp.x_l = [0.5*L/ps.x];                                          
pp.y_l = [0.5*W/ps.x];                                          
pp.V_l = [1e8/(ps.Q0*ps.t)];                                      % 不同 V_l
pp.t_drainage = vv.t+ [20*pd.td/ps.t];                            % time of lake drainages
pp.t_duration = [0.25*pd.td/ps.t];                                % duration of lake drainages
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);            % calculate lake catchments 

%% timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.pts_ni = [pp.ni_l pp.ni_m];                                              
oo.t_span = vv.t + [(1:0.2:36)*pd.td/ps.t (37:1:365)*pd.td/ps.t];            

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
