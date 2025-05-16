% 30 July 2014: sample nevis run

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
x = linspace(0,(10000/ps.x),101); y = linspace(0,(10000/ps.x),101);
gg = nevis_grid(x,y,oo);
b = (0/ps.z)*gg.nx;
s = (500/ps.z)*2*max(0.5^2-(gg.nx-0.5).^2-(gg.ny-0.5).^2,0).^(1/2);

%% mask grid
gg = nevis_mask(gg,find(s-b<=0)); 
gg.n1m = gg.n1;                 % margin boundary nodes
gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes

% %% plot grid
% nevis_plot_grid(gg);

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden 
vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);         % 10cm thick sheet

%% moulins
oo.random_moulins = 50;
[pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);                    

%% surface input
pp.meltE = @(t) (100/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));                          % runoff function; ramp up input over timescale 30 days
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = pp.ni_m;   % hourly timesteps, save timesteps, save moulin pressures

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
[tt,vv] = nevis_timesteps([0:50]*(pd.td/ps.t),vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot discharge
nevis_plot;

%% plot summary
nevis_summary;

%% animate
nevis_animate([oo.root,oo.fn],1:50,1,0);


