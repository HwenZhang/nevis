% 22 Oct 2025
% taken from nevis_140921b, for a greenland margin example

format compact;
clear
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
addpath(oo.code);
script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(script_dir);
load(fullfile(repo_root, 'data', 'topo_b.mat'));  % load topography
oo.reversey = 1;

%% parameters
pd.sigma = 1e-4;
pd.E_lapse = 5/1000/(24*60*60)/10^3; 
pd.phi_s = 0;
[pd,oo] = nevis_defaults(pd,oo);
[ps,pp,sc] = nevis_nondimension(pd);

%% grid and geometry
ps.x0 = min(dd.x); ps.y0 = min(dd.y);
gg = nevis_grid((dd.x-ps.x0)/ps.x,(dd.y-ps.y0)/ps.x,oo); 
b = dd.b/ps.z;
s = dd.s/ps.z;
H = dd.H/ps.z;

%% mask with border
% H = max(s-b,0);
% in_inds = inpolygon(gg.nx,gg.ny,(dd.x_out-ps.x0)/ps.x,(dd.y_out-ps.x0)/ps.x); 
% H(~in_inds) = 0;
Hmin = 15/ps.z; 
gg = nevis_mask(gg,find(H<Hmin)); 
gg.n1m = gg.n1(b(gg.n1)<0); % margin nodes identified with bedrock below sea level
gg.n1m = gg.n1;

%% label boundary nodes
gg = nevis_label(gg,gg.n1m); oo.adjust_boundaries = 1;

%% plot grid
% nevis_plot_grid(gg); 

%% initialize
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a); 
vv.hs = 0.1/ps.hs*ones(gg.nIJ,1);

%% moulins
% oo.random_moulins = 20;
% [pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);                    

%% surface input
pp.meltE = @(t) (10/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t)));   % runoff function; ramp up input over time
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; %oo.save_pts_all = 1; oo.pts_ni = pp.ni_m; 

%% add parameters and boundary labels for ice velocity
pd.n_glen = 1;
pd.A_glen = 1/2/((0.1)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make dimensionless viscosity term a given size
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity

% oo.boundary_method = 'vel_tblr'; % prescribe velocity on all boundary edges of the grid
gg = nevis_label_ice(gg,oo); % add boundary labels needed for ice velocity 
figure(1); clf; nevis_plot_grid_ice(gg); 

%% initial ice velocity
pp.C = 0.1; pp.C2 = 0.1;
N = 1*ones(gg.nIJ,1); 
u = zeros(gg.eIJ,1); v = zeros(gg.fIJ,1); u(gg.eout2) = NaN; v(gg.fout2) = NaN;
[u,v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
[tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b,taudx,taudy] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
vv.u = u; vv.v = v;

figure(2); clf; imagesc(gg.nx(:,1),gg.ny(1,:),reshape(((gg.nmeanx*u).^2+(gg.nmeany*v).^2).^(1/2),gg.nI,gg.nJ)'); colorbar; shg
% tmp = pcolor(gg.ex,gg.ey,reshape(u,gg.eI,gg.eJ)); colorbar; %set(tmp,'linestyle','none'); 
% imagesc(gg.ex(:,1),gg.ey(1,:),reshape(u,gg.eI,gg.eJ)'); colorbar;

oo.include_ice = 1;

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','sc','gg','aa','vv','oo');

%% timestep 
[tt,vv,info] = nevis_timesteps_ice([0:7:1*365]*pd.td/ps.t,vv,aa,pp,gg,oo); 

%% plot summary
nevis_summary;

%% animate
nevis_animate([oo.root,oo.fn],1:53,6,0);
