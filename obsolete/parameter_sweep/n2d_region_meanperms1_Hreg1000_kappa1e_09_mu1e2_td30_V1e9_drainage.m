
%% Script to run NEVIS regional model - DRAINAGE PHASE
% This script loads spinup results and runs drainage simulation
% The surface runoff is set up by a prescribed function
% The surface and bed profiles are the same as in Hewitt 2013

% Author: Auto-generated script  
% Date: 2025-09

format compact

%% read in the initial condition
casename = 'n2d_region_meanperms1_Hreg1000_kappa1e_09_mu1e2_td30_V1e9_drainage'; % drainage system filename
initname = 'n2d_region_meanperms1_Hreg1000_kappa1e_09_mu1e2_spinup';

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
init_cond = load(['./results/' oo.initname '/' '0365.mat']); % load initial condition
vv = init_cond.vv;                                % load state variables from the initial condition

%% supraglacial lakes
% pp.x_l = [0.5*(x(1)+x(end))];                                 % x-coord of lakes
% pp.y_l = [0.5*(y(1)+y(end))];                                 % y-coord of lakes
pp.x_l = 40*(1e3/ps.x);                                         % x-coord of lakes
pp.y_l = -10*(1e3/ps.x);                                        % y-coord of lakes
pp.V_l = [1000000000.0/(ps.Q0*ps.t)];                                    % volume of lakes         
pp.t_drainage = [30*pd.td/ps.t];                      % time of lake drainages (assumed to be the middle time of the Gaussian)
pp.t_duration = [0.25*pd.td/ps.t];                              % duration of lake drainages, 6hr
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments  

%% surface runoff
load([oo.dn '/' 'runoff_2010_nevis140.mat']);                    % load data for year of interest (previously collated)
% RACMO distributed input
oo.surface_runoff = 1;                                           % If set to 1 turns on surface runoff input
                                                                 % 0 is prescribe moulin input with a function
oo.RACMO_runoff = 1;                                             % If set to 1 turns on RACMO surface runoff data
                                                                 % 0 is a simple sinusoidal function
oo.distributed_input = 0;                                        % If set to 1 distributes input across the domain
                                                                 % 0 inputs to discrete moulins

pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2010_nevis140)./ps.m;  % distributed input (m/sec)
oo.pts_ni = [pp.ni_m; pp.ni_l];                                   % points to save moulin and lake pressures

%% timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.t_span = (1:1:365)*pd.td/ps.t; % 1 year of timesteps

% Add GPS station points downstream of the moulin every 5km
% pp.ni_gps = nevis_gps_array((pp.x_l:5e3/ps.x:L/ps.x), pp.y_l, gg, oo); % GPS station points
% oo.pts_ni = [pp.ni_l pp.ni_m pp.ni_gps];    

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
