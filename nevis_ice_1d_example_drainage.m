% 10 October 2025 experimental nevis and ice velocity run
% square-root shaped geometry with random arrangement of moulins

clear 
casename = 'n1d_example_drainage_mu1e1_kappa1e_10_eps0_1';
initname = strrep(casename, '_drainage', '_spinup'); % initial condition filename

%% parameters
data = load(['./results/' initname '/' initname]);
pd = data.pd;                                % load parameters from the initial condition
ps = data.ps;                                % load state variables from the initial condition
pp = data.pp;                                % load scaled parameters from the initial condition
aa = data.aa;                                % load state variables from the initial condition
oo = data.oo;                                % load options from the initial condition

oo.casename = casename;                      % drainage system filename
oo.initname = initname;                      % initial condition filename, for spinup

oo.fn = ['/',oo.casename];                     % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
addpath(oo.code);                              % add path to code
mkdir(oo.rn);                                  % create directory for results    
oo.Tol_F = 1e-8;
oo.cavity_coupling = 0;                        % couple to ice velocity

%% grid and geometry, linear bedslope
L = 100000; % length in m
W = 50000; % width in m
x = linspace(0,(L/ps.x),201); 
y = linspace(0,(W/ps.x),1);
oo.xperiodic = 0;
oo.yperiodic = 1;
gg = nevis_grid(x,y,oo);
b = (0/ps.z)*gg.nx - (0.001*ps.x/ps.z)*(gg.nx);  % linear bedslope
% s = (1000/ps.z)*max(1-gg.nx/max(max(gg.nx)),0).^(1/2);
s = b + (2000/ps.z) - (3000/ps.z).*(gg.nx>0.95*L/ps.x);

%% mask grid
% gg = nevis_mask(gg,[]); 
gg = nevis_mask(gg,find(s-b<=0)); 
% gg = nevis_mask(gg,find(((gg.nx-0.5).^2+(gg.ny-0.5).^2).^(1/2)<0.2)); 
% gg = nevis_mask(gg,unique([gg.bdy.nlbdy; gg.bdy.nrbdy; gg.bdy.ntopbdy; gg.bdy.nbotbdy])); % mask out nodes around boundary
gg.n1m = gg.n1;                 % margin boundary nodes
gg = nevis_label(gg,gg.n1m);    % label pressure boundary nodes

%% plot grid
nevis_plot_grid(gg);

%% initialize
init_cond = load(['./results/' initname '/' '7300.mat']); % load initial condition
vv = init_cond.vv;                                           % load state variables from the initial condition

%% add parameters and boundary labels for ice velocity
% pd.n_glen = 1;
% eps = 0.1; 
% pd.A_glen = 1/2/((eps)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make membrane stress terms of dimensionless size eps in momentum equation
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity
gg = nevis_label_ice(gg); % add boundary labels needed for ice velocity
figure(1); clf; 
nevis_plot_grid_ice(gg); 

%% moulins
oo.keep_all_moulins = 0;
oo.random_moulins = 20;
[pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);                    

%% supraglacial lakes
pp.x_l = [0.5*L/ps.x];                                          % x-coord of lakes
pp.y_l = [0.0*W/ps.x];                                          % y-coord of lakes
pp.V_l = [1e7/(ps.Q0*ps.t)];                                    % volume of lakes         
pp.t_drainage = vv.t + 0.2*365*pd.td/ps.t;                      % time of lake drainages
pp.t_duration = [0.025*pd.td/ps.t];                             % duration of lake drainages
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments 
oo.pts_ni = [pp.ni_l; pp.ni_m];  

%% surface input
pp.meltE = @(t) (0*100/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); % runoff function; ramp up input over timescale 30 days
pp.input_function = @(t) 0*(1-exp(-t/(30*pd.td/ps.t)))./(ps.m*ps.x^2);

%% save initial parameters
save([oo.rn,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
oo.dt = 1/24*pd.td/ps.t; 
oo.save_timesteps = 1; 
oo.save_pts_all = 1; 
oo.t_span = vv.t + (0:0.25:365*1.5)*pd.td/ps.t;         
[tt,vv] = nevis_timesteps_ice(oo.t_span,vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot discharge
% nevis_plot;
save([oo.rn,oo.fn],'pp','pd','ps','gg','aa','vv','oo','tt');

%% animate
% nevis_2d_animation_ice_drainage;
% 
% return;
% 
% figure(3); clf; 
% % imagesc(gg.ex(:,1),gg.ey(1,:),reshape(taudx,gg.eI,gg.eJ)'); colorbar;
% imagesc(gg.nx(:,1),gg.ny(1,:),reshape(tauxx,gg.nI,gg.nJ)'); colorbar;
% 
% [tau1,tau2,theta] = nevis_principal_stress(Txx,Tyy,gg.nmeanc*Txy);
% imagesc(gg.nx(:,1),gg.ny(1,:),reshape(tau1,gg.nI,gg.nJ)'); colorbar;



