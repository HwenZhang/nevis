% 5 November 2024: sample nevis run in one dimension

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
addpath(oo.code);       % add path to code

%% parameters
[pd,oo] = nevis_defaults([],oo);
% [ put non-default parameters and options here ]
[ps,pp] = nevis_nondimension(pd);

%% grid and geometry
x = linspace(0,(10000/ps.x),101); y = linspace(0,(10000/ps.x),1); % 1-d grid of length 10km 
oo.yperiodic = 1; % oo.yperiodic = 1 necessary for a 1-d grid
oo.xperiodic = 0;
gg = nevis_grid(x,y,oo); 
b = (0/ps.z)*gg.nx; % flat bed
s = (500/ps.z)*((10000/ps.x)-gg.nx); % linear surface topography

%% mask grid
gg = nevis_mask(gg,find(s-b<=0)); 
gg.n1m = gg.n1;                 % margin boundary nodes
gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes

%% plot grid
% nevis_plot_grid(gg);
% return;

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden 
vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);         % 10cm thick sheet

%% moulins
[pp.ni_m,pp.sum_m] = nevis_moulins(1000/ps.x,0,gg,oo);  % one moulin at x=1km,y=0                 

%% surface input
% runoff function; ramp up input over timescale 30 days
pp.meltE = @(t) (10/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 

% runoff function at each node is precribed in nevis_inputs.m
% r = (1-pp.E_amp*cos(2*pi*t/pp.td)).*max(pp.meltE(t)-pp.E_lapse*aa.s,0); 
% r(gg.nout) = 0;

% hourly timesteps, save timesteps, save moulin pressures
oo.dt = 1/24*pd.td/ps.t; %h ourly timesteps
oo.save_timesteps = 1; % save timesteps
oo.save_pts_all = 1; 
oo.pts_ni = pp.ni_m;   

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
[tt,vv] = nevis_timesteps([0:100]*(pd.td/ps.t),vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot summary
nevis_summary;

%% expand/update variables
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge


%% plot
nevis_1d_plot


%% Simple animate
nevis_1d_animate