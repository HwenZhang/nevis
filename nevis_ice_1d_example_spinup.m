% 10 October 2025 experimental nevis and ice velocity run
% square-root shaped geometry with random arrangement of moulins

clear 
oo.root = './';             % filename root
oo.code = '../nevis/src';   % code directory
oo.results = 'results';     % path to the results folders
addpath(oo.code);           % add path to code
oo.casename = 'n1d_example_spinup_mu1e1_kappa1e_10_eps0_1';
oo.fn = ['/',oo.casename];               % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];      % path to the case results
mkdir(oo.rn);

%% parameters
[pd,oo] = nevis_defaults([],oo);
% [ put non-default parameters and options here ]
oo.evaluate_variables = 1;
oo.input_gaussian = 1;
oo.relaxation_term = 1;                         % 0 is alpha hb, 1 is alpha deltap hb
oo.initial_condition = 0;                       
oo.cavity_coupling = 1;                         % couple to ice velocity

moulin_input = 0;                               % prescribed moulin input (m^3/s)
pd.mu = 1e1;                                    % water viscosity (Pa s)
pd.Ye = 8.8e9;                                  % Young's modulus (Pa)
pd.B = pd.Ye*(1e3)^3/(12*(1-0.33^2));           % bending stiffness (Pa m^3)
pd.E_lapse = 30/1000/pd.td/10^3;

pd.k_s = 1e-3;                                  % sheet roughness parameter

pd.hb_reg1 = 5e-3;                              % Regularisation parameter for hb
pd.hb_reg2 = 1e-3;                              % Regularisation parameter for hb
pd.N_reg1 = 1e4;                                % Regularisation parameter for N
pd.deltap_reg = 1e4;                            % Regularisation parameter for deltap
pd.B_reg = 0;                                   % Regularisation parameter for B

pd.G = 0.01;                                    % geothermal heat flux [J/s/m^2]
pd.melt = pd.G/pd.rho_w/pd.L;                   % geothermal heat derived basal melt [m/s]
pd.alpha_b = 0;                                 % relaxation rate (s^-1)
pd.kappa_b = 1e-10;                             % relaxation coeff
pd.c0 = 1;

% pd.C = 1e99;;
% pd.mu = 0.01;

[ps,pp] = nevis_nondimension(pd,[],oo);

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
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);
vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden 
vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);         % 10cm thick sheet
vv.hb = (0.1/ps.hb)*ones(gg.nIJ,1);         % 10cm thick bed

%% add parameters and boundary labels for ice velocity
pd.n_glen = 1;
eps = 0.1; 
pd.A_glen = 1/2/((eps)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make membrane stress terms of dimensionless size eps in momentum equation
[pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity

gg = nevis_label_ice(gg); % add boundary labels needed for ice velocity
figure(1); clf; 
nevis_plot_grid_ice(gg); 

%% initial ice velocity
N = ones(gg.nIJ,1); 
u = zeros(gg.eIJ,1); 
v = zeros(gg.fIJ,1); 

u(gg.eout2) = NaN; 
v(gg.fout2) = NaN;
[u,v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
[tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
[tau1,tau2,theta] = nevis_principal_stress(Txx,Tyy,gg.nmeanc*Txy);

vv.u = u; 
vv.v = v;

% % tmp = pcolor(gg.ex,gg.ey,reshape(u,gg.eI,gg.eJ)); colorbar; %set(tmp,'linestyle','none'); 
% % imagesc(gg.ex(:,1),gg.ey(1,:),reshape(u,gg.eI,gg.eJ)'); colorbar;
% imagesc(gg.nx(:,1),gg.ny(1,:),reshape(((gg.nmeanx*u).^2+(gg.nmeany*v).^2).^(1/2),gg.nI,gg.nJ)'); colorbar; shg
% 
% return;

oo.include_ice = 1; % include ice velocity solve

%% moulins
oo.keep_all_moulins = 0;
oo.random_moulins = 20;
[pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);                    

%% supraglacial lakes
pp.x_l = [0.5*L/ps.x];                                          % x-coord of lakes
pp.y_l = [0.0*W/ps.x];                                          % y-coord of lakes
pp.V_l = [0e7/(ps.Q0*ps.t)];                                    % volume of lakes         
pp.t_drainage = [0*365*pd.td/ps.t];                             % time of lake drainages
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
oo.t_span = (0:1:365*20)*pd.td/ps.t;         
[tt,vv] = nevis_timesteps_ice(oo.t_span,vv,aa,pp,gg,oo);     % save at daily timesteps

%% plot discharge
% nevis_plot;
save([oo.rn,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% animate
% nevis_2d_animation_ice;

%% plot stresses
figure(3); clf; 
% imagesc(gg.ex(:,1),gg.ey(1,:),reshape(taudx,gg.eI,gg.eJ)'); colorbar;
imagesc(gg.nx(:,1),gg.ny(1,:),reshape(tauxx,gg.nI,gg.nJ)'); colorbar;

[tau1,tau2,theta] = nevis_principal_stress(Txx,Tyy,gg.nmeanc*Txy);
imagesc(gg.nx(:,1),gg.ny(1,:),reshape(tau1,gg.nI,gg.nJ)'); colorbar;



