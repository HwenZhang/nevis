% 29 August 2014: sample run for circular or rectangular domain

format compact;
% clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
addpath(oo.code);

%% parameters
% [ put non-default parameters and options here ]
[pd,oo] = nevis_defaults([],oo);
[ps,pp] = nevis_nondimension(pd);

% %% grid and geometry
% gg = nevis_grid(50,50,1,0,1,0,1,oo);
% b = reshape(0*gg.nx,gg.nIJ,1);
% s = max(reshape(max(1-gg.nx.^2-gg.ny.^2,0).^(1/2),gg.nIJ,1),b);
% 
% %% mask grid
% gg = nevis_mask(gg,find(s-b<=0)); 
% gg.n1m = gg.n1;              % margin boundary nodes
% gg = nevis_label(gg,gg.n1m); % Dirichlet boundary nodes

%% grid and geometry
x = linspace(0,1,50); y = linspace(0,1,50);
gg = nevis_grid(x,y,oo);
b = reshape(0*gg.nx,gg.nIJ,1);
s = max(reshape(max(1-gg.nx.^2,0).^(1/2),gg.nIJ,1),b);

%% mask grid
% gg = nevis_mask(gg,find(s-b<0)); 
gg = nevis_label(gg,'l_e_r_n'); oo.adjust_boundaries = 0; % Dirichlet boundary nodes

%% plot grid
nevis_plot_grid(gg);

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);
vv.hs = (0.1/ps.hs)*gg.nx.^0;

%% boundary conditions
aa.phi = aa.phi_a(gg.nbdy);

%% surface input
pp.meltE = @(t) (10/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 0;
% oo.save_pts_all = 1; oo.pts_ni = pp.ni_m; 

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
[tt,vv] = nevis_timesteps([1:10]*pd.td/ps.t,vv,aa,pp,gg,oo); 
