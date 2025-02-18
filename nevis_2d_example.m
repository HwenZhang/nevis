 % 5 November 2024: sample nevis run in one dimension

 clear
 oo.root = '';           % filename root
 oo.fn = mfilename;      % filename
 oo.code = '../nevis';   % code directory
 oo.casename = 'time_series_blister_vk_concentrated_reg_test';
 oo.casename = 'test_2d_distributed_k0_kbed1e_13_rm60';
%  oo.casename = 'test_2d_no_blister';
%  oo.casename = 'test_2d_concentrated_k0_kbed1e_12_rm60';
 addpath(oo.code);       % add path to code
 oo.evaluate_variables = 1;
 
 % Physical set-up
 oo.constant_k = 1.0;      % constant permeability
 oo.sheet_k = 0.0;         % permeability related to cavity sheet
 oo.channel_k = 1-oo.constant_k-oo.sheet_k; % permeability related to channels
 oo.include_radius = 1;  % include Rb in the Jacobian 
 oo.distributed_lake = 1;

 %% parameters
 [pd,oo] = nevis_defaults([],oo);
 % [ put non-default parameters and options here ]
 [ps,pp] = nevis_nondimension(pd);
 
 %% grid and geometry
 % The geometry is taken from Hewitt 2013
 L = 5e4;                                 % length of the domain [m]
 W = 0.2*L;                               % width of the domain [m]
 x = linspace(0,(L/ps.x),101); 
 y = linspace(0,(W/ps.x),21);        
 oo.yperiodic = 1;                        % oo.yperiodic = 1 necessary for a 1-d grid
 oo.xperiodic = 0;
 gg = nevis_grid(x,y,oo); 
 b = (0/ps.z)*gg.nx;                      % flat bed
 s = (1060/ps.z)*(1-ps.x*gg.nx/L).^0.5;   % ice surface topography 
 
 %% mask grid
 gg = nevis_mask(gg,find(s-b<=0)); 
 gg.n1m = gg.n1;                  % margin boundary nodes
 gg = nevis_label(gg,gg.n1m);     % label pressure boundary nodes
 
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
 
 % lakes
 % oo.lake_dict = './lake_data';   % code directory
 % pp.x_l = load([oo.lake_dict "/lake_data.mat"],"x_l");        % x-coord of lakes
 % pp.y_l = load([oo.lake_dict "/lake_data.mat"],"y_l");        % y-coord of lakes
 % pp.V_l = load([oo.lake_dict "/lake_data.mat"],"V_l");        % volume of lakes         
 % pp.t_drainage = load([oo.lake_dict "/lake_data.mat"],"t_l"); % time of lake drainages (assumed to be the middle time of the Gaussian)
 % pp.t_duration = load([oo.lake_dict "/lake_data.mat"],"delta_t_l"); % duration of lake drainages, 6hr
if ~oo.distributed_lake
    % a single-point lake
    pp.x_l = [0.75*L/ps.x];                             % x-coord of lakes
    pp.y_l = [0.5*W/ps.x];                              % y-coord of lakes
    pp.V_l = [1e8/(ps.Qb_0*ps.t)];                      % volume of lakes         
    pp.t_drainage = [30.0];                             % time of lake drainages (assumed to be the middle time of the Gaussian)
    pp.t_duration = [0.025];                            % duration of lake drainages, 6hr
else
    % distribute the lakes into multiple grid points
    vx = gg.nx(:); vy = gg.ny(:);
    index_lake = ((vx-0.75*L/ps.x).^2 + (vy-0.5*W/ps.x).^2).^0.5 < 3.7e3/ps.x;
    pp.x_l = vx(index_lake);                            % x-coord of lakes
    pp.y_l = vy(index_lake);                            % y-coord of lakes
    pp.V_l = 1e8/sum(index_lake)/(ps.Qb_0*ps.t)*ones(sum(index_lake),1);                      % volume of lakes         
    pp.t_drainage = 30.0*ones(sum(index_lake),1);                             % time of lake drainages (assumed to be the middle time of the Gaussian)
    pp.t_duration = 0.025*ones(sum(index_lake),1);                            % duration of lake drainages, 6hr
end
[pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo); % calculate lake catchments  
vv.Rb(pp.ni_l) = ((ps.Qb_0*ps.t)*pp.V_l/ps.K_0).^0.4/ps.R; % initial blister radius;  

%% moulins
 oo.random_moulins = 0;
 pp.x_l = load("moulins.mat","x_m");        % x-coord of moulins
 pp.y_l = load("moulins.mat","y_m");        % y-coord of moulins
%  [pp.ni_m,pp.sum_m] = nevis_moulins([],[],gg,oo);         
[pp.ni_m,pp.sum_m] = nevis_moulins(pp.x_l.x_m,pp.y_l.y_m,gg,oo);  % one moulin at the lake location

 %% parameters for timesteping
 % hourly timesteps, save timesteps, save moulin pressures
 oo.dt = 0.01/24*pd.td/ps.t; % hourly timesteps
 oo.save_timesteps = 1; % save timesteps
 oo.save_pts_all = 1; 
 % oo.pts_ni = [pp.ni_m;pp.ni_l];  
 oo.pts_ni = pp.ni_l;   
 
 %% save initial parameters
 save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');
 
 %% timestep 
 t_span = (0:1000)*(0.4*pd.td/ps.t);
 % t_span = (0:1500)*(0.2*pd.td/ps.t);
 [tt,vv] = nevis_timesteps(t_span,vv,aa,pp,gg,oo);     % save at hourly timesteps
 
 %% plot summary
 if oo.distributed_lake == 1 
    nevis_summary_distributed; 
 else 
    nevis_summary;
 end
 
 %% expand/update variables
 aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
 oo.evaluate_variables = 1; 
 [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
 vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
 
 %% plot
 nevis_2d_plot
 
 %% Simple animate
 nevis_2d_animate