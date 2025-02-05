 % 5 November 2024: sample nevis run in one dimension

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
oo.casename = 'time_series_blister_vk_concentrated_reg_test';
oo.casename = 'test';
addpath(oo.code);       % add path to code
oo.evaluate_variables = 1;

% Physical set-up
oo.constant_k = 0;      % constant permeability
oo.include_radius = 1;  % include Rb in the Jacobian 
%% parameters
[pd,oo] = nevis_defaults([],oo);
% [ put non-default parameters and options here ]
[ps,pp] = nevis_nondimension(pd);

%% grid and geometry
L = 1e5;                % length of the domain [m]
x = linspace(0,(L/ps.x),101); 
y = linspace(0,(L/ps.x),1); % 1-d grid of length 100km 
oo.yperiodic = 1; % oo.yperiodic = 1 necessary for a 1-d grid
oo.xperiodic = 0;
gg = nevis_grid(x,y,oo); 
b = (0/ps.z)*gg.nx; % flat bed
s = (300/ps.z)*((L/ps.x)-gg.nx); % linear surface topography

%% mask grid
gg = nevis_mask(gg,find(s-b<=0)); 
gg.n1m = gg.n1;                  % margin boundary nodes
gg = nevis_label(gg,gg.n1m);     % label pressure boundary nodes

%% plot grid
% nevis_plot_grid(gg);
% return;

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
% overburden = phi_0 - phi_a;
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden into hydraulic potential, 10% to eff p 
vv.hs = (0.2/ps.hs)*ones(gg.nIJ,1);         % 10cm thick water sheet

%% Surface hydrology

% surface runoff
% runoff function; ramp up input over timescale 30 days (dimensionless)
pp.meltE = @(t) 0.2*(10/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 

% moulins
% [pp.ni_m,pp.sum_m] = nevis_moulins(1000/ps.x,0,gg,oo);  % one moulin at x=1km,y=0  

% lakes
% oo.lake_dict = './lake_data';   % code directory
% pp.x_l = load([oo.lake_dict "/lake_data.mat"],"x_l");        % x-coord of lakes
% pp.y_l = load([oo.lake_dict "/lake_data.mat"],"y_l");        % y-coord of lakes
% pp.V_l = load([oo.lake_dict "/lake_data.mat"],"V_l");        % volume of lakes         
% pp.t_drainage = load([oo.lake_dict "/lake_data.mat"],"t_l"); % time of lake drainages (assumed to be the middle time of the Gaussian)
% pp.t_duration = load([oo.lake_dict "/lake_data.mat"],"delta_t_l"); % duration of lake drainages, 6hr

pp.x_l = [0.5*L/ps.x];                              % x-coord of lakes
pp.y_l = [0];                                       % y-coord of lakes
pp.V_l = [1e8/(ps.Qb_0*ps.t)];                      % volume of lakes         
pp.t_drainage = [20.0];                             % time of lake drainages (assumed to be the middle time of the Gaussian)
pp.t_duration = [0.025];                            % duration of lake drainages, 6hr

[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo); % calculate lake catchments  
vv.Rb(pp.ni_l) = ((ps.Qb_0*ps.t)*pp.V_l/ps.K_0).^0.4/ps.R; % initial blister radius;  

%% parameters for timesteping
% hourly timesteps, save timesteps, save moulin pressures
oo.dt = 0.01/24*pd.td/ps.t; % hourly timesteps
oo.save_timesteps = 1; % save timesteps
oo.save_pts_all = 1; 
% oo.pts_ni = [pp.ni_m;pp.ni_l];  
oo.pts_ni = pp.ni_l;   

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
t_span = (0:750)*(0.4*pd.td/ps.t);
t_span = (0:1500)*(0.2*pd.td/ps.t);
[tt,vv] = nevis_timesteps(t_span,vv,aa,pp,gg,oo);     % save at hourly timesteps

%% plot summary
nevis_summary;

%% expand/update variables
aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge


%% plot
% nevis_1d_plot


%% Simple animate
% nevis_1d_animate;