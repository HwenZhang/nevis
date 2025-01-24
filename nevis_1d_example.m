 % 5 November 2024: sample nevis run in one dimension

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
oo.casename = 'time_series_blister_ck_concentrated_reg_1';
oo.include_radius = 1;  % dont calculate radius Rb 
addpath(oo.code);       % add path to code

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
vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);         % 10cm thick water sheet

%% moulins
% [pp.ni_m,pp.sum_m] = nevis_moulins(1000/ps.x,0,gg,oo);  % one moulin at x=1km,y=0         

%% lakes
[pp.ni_l,pp.sum_l] = nevis_lakes(0.5*L/ps.x,0,gg,oo);  % one lake at x=1km,y=0

%% surface input
% runoff function; ramp up input over timescale 30 days (dimensionless)
pp.meltE = @(t) 0.2*(10/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 

% lake input function
t_input_1 = 20.0;
t_0 = pd.td/ps.t;
sigma = 0.025;
% pp.lake_input_function = @(t) (t>=t_input && t<=t_input+0.1).*(pd.Q_0)/pd.Q_0; % Heaviside
pp.lake_input_function = @(t) (pd.Q_0*t_0)/(sqrt(2*pi)*sigma*pd.Q_0)*...
                       exp(-0.5/sigma^2*(t-t_input_1)^2).*(t>=t_input_1-5*sigma && t<=t_input_1+5*sigma);

% dimensionless permeability k as a function of dimensionless sheet thickness h
% pp.k_blister = @(h)     14*pd.mu*pd.k_s*(ps.hs*h).^3/pd.k_bed; % dimensional value/scale 
pp.k_blister = @(h)     1.0; % dimensional value/scale 
% % dimensionless fracture toughness K_c as a function of dimensional sheet thickness h
% pp.K_c = @(h)           pd.K_1c/pd.K_1c;

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
[tt,vv] = nevis_timesteps(t_span,vv,aa,pp,gg,oo);     % save at hourly timesteps

%% plot summary
nevis_summary;

%% expand/update variables
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; 
[vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge


%% plot
% nevis_1d_plot


%% Simple animate
% nevis_1d_animate;