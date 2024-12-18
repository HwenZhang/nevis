% 30 July 2014: sample nevis run

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
oo.casename = 'time_series_2d_blister_ck_distributed';
addpath(oo.code);       % add path to code

%% parameters
[pd,oo] = nevis_defaults([],oo);
% [ put non-default parameters and options here ]
[ps,pp] = nevis_nondimension(pd);

%% grid and geometry
% x and y coordinates
L = 1e4;                % length of the domain [m]
x = linspace(0,(L/ps.x),21); y = linspace(0,(L/ps.x),21);
% gg: coords of nodes and edges, oo: optianl flags?
gg = nevis_grid(x,y,oo);
% precribe bed topography
b = (0/ps.z)*gg.nx; % flat bed
% topo high at the centre
s = (500/ps.z)*2*max(0.5^2-(gg.nx-0.5).^2-(gg.ny-0.5).^2,0).^(1/2);

%% mask grid
gg = nevis_mask(gg,find(s-b<=0)); 
gg.n1m = gg.n1;                 % margin boundary nodes
gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes

% %% plot grid
nevis_plot_grid(gg);

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden, 10% eff stress 
vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);         % 10cm thick sheet

%% moulins
% oo.random_moulins = 20;
% [pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);      

%% lakes
[pp.ni_l,pp.sum_l] = nevis_lakes(0.5*L/ps.x,0.5*L/ps.x,gg,oo);  % one lake at the centre

%% surface input
% runoff function; ramp up input over timescale 30 days
pp.meltE = @(t) (100/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));                          

% lake input function
t_input_1 = 20.0;
t_0 = pd.td/ps.t;
sigma = 0.1;
% pp.lake_input_function = @(t) (t>=t_input && t<=t_input+0.1).*(pd.Q_0)/pd.Q_0; % Heaviside
pp.lake_input_function = @(t) (pd.Q_0*t_0)/(sqrt(2*pi)*sigma*pd.Q_0)*...
      exp(-0.5/sigma^2*(t-t_input_1)^2).*(t>=t_input_1-5*sigma && t<=t_input_1+5*sigma);

% dimensionless permeability k as a function of dimensionless sheet thickness h
% pp.k_blister = @(h)     11*pd.mu*pd.k_s*(ps.hs*h).^3/pd.k_bed; % dimensional value/scale 
pp.k_blister = @(h)      4.0; % dimensional value/scale 
oo.blister_distributed = 0;

%% parameters for timesteping
% hourly timesteps, save timesteps, save moulin pressures
oo.dt = 0.01/24*pd.td/ps.t;  % hourly timesteps
oo.save_timesteps = 1;       % save timesteps
oo.save_pts_all = 1;         % save moulin pressures
oo.pts_ni = pp.ni_l;  

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
t_span = (0:750)*(0.4*pd.td/ps.t);
[tt,vv] = nevis_timesteps(t_span,vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot discharge
nevis_plot;

%% plot summary
nevis_summary;

%% animate
nevis_animate([oo.root,oo.fn],450:550,1,0);


