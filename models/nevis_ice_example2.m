% 10 October 2025 experimental nevis and ice velocity run
% circular geometry with random arrangement of moulins

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

%% add parameters and boundary labels for ice velocity
pd.n_glen = 1;
eps = 0.01; pd.A_glen = 1/2/((eps)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make membrane stress terms of dimensionless size eps in momentum equation;
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity

gg = nevis_label_ice(gg); % add boundary labels needed for ice velocity
figure(1); clf; nevis_plot_grid_ice(gg); 

%% initial ice velocity
N = ones(gg.nIJ,1); 
u = zeros(gg.eIJ,1); v = zeros(gg.fIJ,1); u(gg.eout2) = NaN; v(gg.fout2) = NaN;
[u,v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
[tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
vv.u = u; vv.v = v;

figure(2); clf; imagesc(gg.nx(:,1),gg.ny(1,:),reshape(((gg.nmeanx*u).^2+(gg.nmeany*v).^2).^(1/2),gg.nI,gg.nJ)'); colorbar; shg
% tmp = pcolor(gg.ex,gg.ey,reshape(u,gg.eI,gg.eJ)); colorbar; %set(tmp,'linestyle','none'); 
% imagesc(gg.ex(:,1),gg.ey(1,:),reshape(u,gg.eI,gg.eJ)'); colorbar;

% return;

%% moulins
oo.random_moulins = 50;
[pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);                    

%% surface input
pp.meltE = @(t) (100/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));                          % runoff function; ramp up input over timescale 30 days
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = pp.ni_m;   % hourly timesteps, save timesteps, save moulin pressures

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
oo.include_ice = 1; % include ice velocity solve
[tt,vv] = nevis_timesteps_ice([0:50]*(pd.td/ps.t),vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot discharge
nevis_plot;

%% plot summary
nevis_summary;

%% animate
nevis_animate([oo.root,oo.fn],1:50,1,0);


