% 4 Nov 2025 
% attempt to roughly recreate run from Hewitt 2013 (though uses a different sliding law)

clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
addpath(oo.code);       % add path to code

%% parameters
[pd,oo] = nevis_defaults([],oo);
% [ put non-default parameters and options here ]
pd.u_b = 60/pd.ty;
pd.melt = pd.G/pd.rho_w/pd.L;
[ps,pp] = nevis_nondimension(pd);

%% grid and geometry
x = linspace(0,(50000/ps.x),101); y = linspace(0,(10000/ps.x),21);
gg = nevis_grid(x,y,oo);
b = (0/ps.z)*gg.nx;
s = ((2*10^5/pd.rho_i/pd.g*ps.x)^(1/2)/ps.z)*max(50000/ps.x-gg.nx,0).^(1/2);

%% mask grid
% gg = nevis_mask(gg,[]); 
gg = nevis_mask(gg,find(s-b<=0)); 
% gg = nevis_mask(gg,find(((gg.nx-0.5).^2+(gg.ny-0.5).^2).^(1/2)<0.2)); 
% gg = nevis_mask(gg,unique([gg.bdy.nlbdy; gg.bdy.nrbdy; gg.bdy.ntopbdy; gg.bdy.nbotbdy])); % mask out nodes around boundary
gg.n1m = gg.n1;                 % margin boundary nodes
gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes

%% plot grid
% nevis_plot_grid(gg); return;

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden 
vv.hs = (0.5/ps.hs)*ones(gg.nIJ,1);         % 50cm thick sheet

%% add parameters and boundary labels for ice velocity
pd.n_glen = 1;
pd.A_glen = 1/2/((0.01)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make dimensionless viscosity term a given size;
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity

gg = nevis_label_ice(gg); % add boundary labels needed for ice velocity 
figure(1); clf; nevis_plot_grid_ice(gg); 

%% initial ice velocity
N = ones(gg.nIJ,1); 
u = zeros(gg.eIJ,1); v = zeros(gg.fIJ,1); u(gg.eout2) = NaN; v(gg.fout2) = NaN;
[u,v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
vv.u = u; vv.v = v;
oo.include_ice = 1; % include ice velocity solve

% moulins
% oo.random_moulins = 50;
load(fullfile('data','redundant','moulins_120511.mat')); % load previously saved moulin locations
[pp.ni_m,pp.sum_m] = nevis_moulins(gg.xl+(gg.xr-gg.xl)*x_m,gg.yb+(gg.yt-gg.yb)*y_m,gg,oo);                    

%% surface input
pp.meltE = @(t) 0*t;                          % runoff function
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; % oo.save_pts_all = 1; oo.pts_ni = pp.ni_m;   % hourly timesteps, save timesteps, save moulin pressures

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
[tt,vv] = nevis_timesteps_ice([0:100]*(pd.td/ps.t),vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot summary
nevis_summary;

%% animate
nevis_animate([oo.root,oo.fn],1:100,6,0);

return;

%% now re-run for annual cycle

%% surface input
pd.E_lapse = 60/1000/pd.td/1e3; pp.E_lapse = pd.E_lapse/ps.m*ps.z;
pp.meltE = @(t) (70/1000/pd.td/ps.m)*annual_signal(ps.t*t,135*pd.td,244*pd.td,21*pd.td);     % runoff function, using annual signal
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; % oo.save_pts_all = 1; oo.pts_ni = pp.ni_m;   % hourly timesteps, save timesteps, save moulin pressures

%% timestep 
[tt,vv] = nevis_timesteps_ice([0:7:365]*(pd.td/ps.t),vv,aa,pp,gg,oo);     

%% plot discharge
% nevis_plot;

%% plot summary
nevis_summary;

%% animate
nevis_animate([oo.root,oo.fn],1:100,6,0);


