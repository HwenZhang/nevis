 % 5 November 2024: sample nevis run in one dimension

 clc,clear
 oo.root = './';                                % filename root
 oo.code = '../nevis';                          % code directory  
 oo.results = 'results';                        % path to the results folders
 oo.casename = 'blister_2d_alpha0_2_kappa0_mu1e9_E8_8e9_V1e6_H1_0e3';

                                                % casename
 oo.fn = ['/',oo.casename];                     % filename (same as casename)
 oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
 addpath(oo.code);                              % add path to code

 %% parameters
 [pd,oo] = nevis_defaults([],oo);               % the default parameters

 oo.evaluate_variables = 1;
 oo.use_modified_N = 0;
 oo.input_gaussian = 1;

 pd.alpha_b = 1.0/(5*pd.td);                    % relaxation rate (s^-1)
 pd.mu = 1.0e9;                                 % water viscosity (Pa s)
 pd.Ye = 8.8e9;                                 % Young's modulus (Pa)
 pd.B = pd.Ye*(1e3)^3/(12*(1-0.33^2));
 pd.kappa_b = 0;                                % 

 % nondimensionalisation
 [ps,pp] = nevis_nondimension(pd);
 
 %% grid and geometry
% x and y coordinates
L = 4e4;                                 % domain size (m)
x = linspace(0,(L/ps.x),201); 
y = linspace(0,(L/ps.x),201);
oo.yperiodic = 0;                        % oo.yperiodic = 1 necessary for a 1-d grid
oo.xperiodic = 0;
% gg: coords of nodes and edges, oo: optianl flags?
gg = nevis_grid(x,y,oo);
% precribe bed topography
b = (0/ps.z)*gg.nx; % flat bed
% topo high at the centre
% s = (1000/ps.z)*max(1-(gg.nx-0.5*L/ps.x).^2-(gg.ny-0.5*L/ps.x).^2,0).^(1/2);
s = (1000/ps.z)*ones(size(gg.nx));

 %% mask grid
 gg = nevis_mask(gg,find(s-b<=0)); 
 gg.n1m = gg.n1;                          % margin boundary nodes
 gg = nevis_label(gg,gg.n1m);             % label pressure boundary nodes
 
 %% plot grid
 % nevis_plot_grid(gg);
 % return;
 
 %% initialize
 [aa,vv] = nevis_initialize(b,s,gg,pp,oo);
 % overburden = phi_0 - phi_a;
 vv.phi = aa.phi_a+0.9*(aa.phi_0-aa.phi_a);  % 90% overburden into hydraulic potential, 10% to eff p 
 vv.hs = (0.1/ps.hs)*ones(gg.nIJ,1);         % 10cm thick water sheet

 %% Surface hydrology
 % surface runoff
 % runoff function: ramp up input over timescale 30 days (dimensionless)
 pp.meltE = @(t) (60/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 
 
 % seasonal runoff function over a year (dimensionless)
 % pp.meltE = @(t) (pp.r_m+pp.E_lapse*pp.s_m)/ps.m*(0.5*tanh((t-pp.t_spr)/pp.delta_t)-0.5*tanh((t-pp.t_aut)/pp.delta_t)); 

 % a single-point lake
 pp.x_l = [0.5*L/ps.x];                              % x-coord of lakes
 pp.y_l = [0.5*L/ps.x];                              % y-coord of lakes
 pp.V_l = [1e6/(ps.Q0*ps.t)];                        % volume of lakes         
 pp.t_drainage = [30.0];                             % time of lake drainages (assumed to be the middle time of the Gaussian)
 pp.t_duration = [0.025];                             % duration of lake drainages, 6hr

 [pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo); % calculate lake catchments 

%% moulins
 oo.random_moulins = 0;
% [pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);         

 %% parameters for timesteping
 % hourly timesteps, save timesteps, save moulin pressures
 oo.dt = 0.01/24*pd.td/ps.t; % hourly timesteps
 oo.save_timesteps = 1; % save timesteps
 oo.save_pts_all = 1; 
 oo.pts_ni = pp.ni_l;   
 
 %% save initial parameters
 save([oo.rn,oo.fn],'pp','pd','ps','gg','aa','vv','oo');
 
 %% timestep 
 t_span = (0:2000)*(0.2*pd.td/ps.t);
 [tt,vv] = nevis_timesteps(t_span,vv,aa,pp,gg,oo);     % save at hourly timesteps
 
 %% plot summary
%  nevis_summary;
 
 %% expand/update variables
 aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
 oo.evaluate_variables = 1; 
 [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
 vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
 
 save([oo.rn,oo.fn],'pp','pd','ps','gg','aa','oo','tt');

 %% plot
%  nevis_2d_plot
%  nevis_summary
 
 %% Simple animate
 nevis_demonstration_animation