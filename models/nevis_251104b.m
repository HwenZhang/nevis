% 4 Nov 2025 
% attempt to roughly recreate run from Hewitt 2013 (though uses a different sliding law)
% b continues a for another year with annual input

clear
load('nevis_251104a'); % load previous solution and variables
oo.root = '';           % filename root
oo.fn = mfilename;      % filename
oo.code = '../nevis';   % code directory
addpath(oo.code);       % add path to code        

%% now re-run for annual cycle

%% surface input
pd.E_lapse = 60/1000/pd.td/1e3; pp.E_lapse = pd.E_lapse/ps.m*ps.z;
pp.meltE = @(t) (70/1000/pd.td/ps.m)*annual_signal(ps.t*t,135*pd.td,244*pd.td,21*pd.td);     % runoff function, using annual signal
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = pp.ni_m;   % hourly timesteps, save timesteps, save moulin pressures
oo.dt_max = 0.1; 

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
[tt,vv] = nevis_timesteps_ice([0:7:365]*(pd.td/ps.t),vv,aa,pp,gg,oo);     

%% plot discharge
% nevis_plot;

%% plot summary
nevis_summary;

%% animate
nevis_animate([oo.root,oo.fn],1:53,6,0);



