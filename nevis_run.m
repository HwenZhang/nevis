% 30 July 2014: sample nevis run

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
addpath(oo.code);       % add path to code

%% parameters
% [ put non-default parameters and options here ]
pd = struct;
[pd,oo] = nevis_defaults(pd,oo);
[ps,pp] = nevis_nondimension(pd);

pp.c42 = ps.phi/10^5;
pp.c43 = 10^4/ps.phi;

%% grid and geometry
x = linspace(0,(10000/ps.x),10); y = linspace(0,(10000/ps.x),10);
gg = nevis_grid(x,y,oo);
b = (0/ps.z)*gg.nx;
s = (1000/ps.z)*2*max(0.5^2-(gg.nx-0.5).^2-(gg.ny-0.5).^2,0).^(1/2);

%% mask grid
gg = nevis_mask(gg,find(s-b<=0)); 
gg.n1m = gg.n1;                 % margin boundary nodes
gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes

% %% grid and geometry
% gg = nevis_grid(50,50,1,0,1,0,1,oo);
% b = reshape(0*gg.nx,gg.nIJ,1);
% s = max(reshape(max(1-gg.nx.^2-gg.ny.^2,0).^(1/2),gg.nIJ,1),b);
% 
% %% mask grid
% gg = nevis_mask(gg,find(s-b<=0)); 
% gg.n1m = gg.n1;              % margin boundary nodes
% gg = nevis_label(gg,gg.n1m); % label pressure boundary nodes


% %% grid and geometry
% x = linspace(0,1,50); y = linspace(0,1,50);
% gg = nevis_grid(x,y,oo);
% b = reshape(0*gg.nx,gg.nIJ,1);
% s = max(reshape(max(1-gg.nx.^2,0).^(1/2),gg.nIJ,1),b);
% 
% %% mask grid
% gg = nevis_label(gg,'l_e_r_n');

% %% plot grid
% nevis_plot_grid(gg);

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a); 
vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);

%% moulins
oo.random_moulins = 20;
[pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);                    

%% surface input
pp.meltE = @(t) (100/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));   % runoff function; ramp up input over time
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = pp.ni_m; 

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
[tt,vv] = nevis_timesteps([0:100]*(pd.td/ps.t),vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot discharge
nevis_plot;

%% animate
nevis_animate([oo.root,oo.fn],1:100,1,0);
