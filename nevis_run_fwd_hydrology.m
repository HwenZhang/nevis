function [N_new, vv_temp] = nevis_run_fwd_hydrology(C_dim, C_hat, vv_prev, mus)
    % Run the forward hydrology model with the current C and return the updated N
    % Inputs:
    %   C_dim:    dimensional slipperiness field [nIJ x 1]
    %   vv_prev:  (optional) state variables from previous outer iteration
    %             If provided, use as initial condition instead of spinup
    % Outputs:
    %   N_new: updated effective pressure field from hydrology model
    %   vv_temp: final state variables (pass to next iteration)
    
    format compact

    %% read in the initial condition
    oo.root = './';                                % filename root
    oo.code = './src';                             % code directory   
    oo.results = 'results';                        % path to the results folders
    oo.dataset = 'nevis_regional';                 % dataset name
    oo.casename = 'n2d_region_ice_inversion_test'; % casename
    oo.fn = ['/',oo.casename];                     % filename (same as casename)
    oo.rn = [oo.root,oo.results,oo.fn];            % path to the case results
    oo.dn = [oo.root, 'data/', oo.dataset, '/'];   % path to the data
    addpath(oo.code);                              % add path to code
    mkdir(oo.rn);                                  % create directory for results 

    %% parameters
    [pd,oo] = nevis_defaults([],oo);  

    oo.evaluate_variables = 1;
    oo.input_gaussian = 1;
    oo.relaxation_term = 1;                         % 0 is alpha hb, 1 is alpha deltap hb
    oo.initial_condition = 1;                       % 1 is default condition from 0365.mat, 0 is using steady-state drainage system, winter or summertime
    oo.mean_perms = 1;
    oo.modified_mean_perms = 0;
    oo.display_residual = 1;
    oo.N_coupling = 1; % turn on effective pressure coupling                                   
    oo.U_coupling = 1; % turn on basal sliding coupling
    oo.boundary_method = 'stress_l_vel_tbl';
    oo.mask_boundary_method = 'stress_free';
    % oo.step_ice = 0.1;
    % oo.max_iter_new = 100;

    pd.alpha_b = 0;                                 % relaxation rate (s^-1)
    pd.kappa_b = 1e-10;                             % relaxation coeff

    % alter default parmaeters 
    pd.mu = 20.0;                                   % water viscosity (Pa s)
    pd.c_e_reg2 = 0.00/1e3/9.81;                    % elastic sheet thickness [m/Pa]
    pd.u_b = 100/pd.ty;                             % sliding speed [m/s]
    pd.sigma = 1e-3;                                % englacial void fraction
    pd.h_r = 0.1;                                   % roughness height [m]
    pd.l_r = 10;                                    % roughness length [m]
    pd.l_c = 1000;                                  % sheet width contributing to conduit melting [m] default = 10 m
    pd.k_s = 1e-3;                                  % sheet permeability constant
    pd.tau_b = 60e3;                                % driving stress [Pa]
    pd.melt = pd.G/pd.rho_w/pd.L;                   % geothermal heat derived basal melt [m/s]
    pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  
                                                    % geothermal heat + frictional heating derived basal melt [m/s]
    pd.meltinterior = 0*((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; 
                                                    % flux of basal melt up to the ~icedivide (200 km) [m2/s]

    pd.hb_reg1 = 5e-3;                              % Regularisation parameter for hb in blister propagation
    pd.hb_reg2 = 1e-3;                              % Regularisation parameter for hb in blister leakage
    pd.N_reg1 = 1e4;                                % Regularisation parameter for N in blister leakage
    pd.deltap_reg = 1e4;                            % Regularisation parameter for deltap
    pd.B_reg = pd.Ye*(1000)^3/(12*(1-0.33)^2);      % Reg parameter for bending stiffness

    % non-dimensionalise
    ps = struct;
    [ps,pp] = nevis_nondimension(pd,ps,oo);
    
    %% grid and geometry
    load([oo.dn '/' 'morlighem_for_nevis_140km']); % load Morlighem bedmap (previously collated)
    dd = morlighem_for_nevis_140km; dd.skip = 6;
    gg = nevis_grid(dd.X_m(1:dd.skip:end,1)/ps.x,dd.Y_m(1,1:dd.skip:end)/ps.x,oo); 
    b = reshape(dd.B_m(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
    s = reshape(dd.S_m(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
    % save a copy of the bed and surface
    gg.bed = b;
    gg.surface = s;

    x = dd.X_m(1:dd.skip:end,1)/ps.x;
    y = dd.Y_m(1,1:dd.skip:end)/ps.x;

    %% mask with minimum ice thickness
    H = max(s-b,0);
    Hmin = 0/ps.z; 
    % y_max = max(max(gg.ny'));
    nout = find(H<=Hmin);
    if isempty(nout)
        x_max = max(gg.nx);
        y_max = max(gg.ny);
        nout = find(abs(gg.nx - x_max) < 1e-10);
    end
    noutb = nout;
    if isempty(noutb)
        x_max = max(max(gg.nx));
        y_max = max(max(gg.ny'));
        noutb = union(find(abs(gg.ny - y_max) < 1e-10), find(abs(gg.nx - x_max) < 1e-10));
    end
    gg = nevis_mask(gg,nout); 
    gg = nevis_mask_blister(gg,noutb);
    gg.n1m = gg.n1;                                   % label all edge nodes as boundary nodes for pressure

    %% label boundary nodes
    gg = nevis_label(gg,gg.n1m);
    gg = nevis_label_blister(gg,gg.n1_blister,oo);    % label blister boundary nodes
    oo.adjust_boundaries = 1;                         % enable option of changing conditions
    % oo.boundary_method = 'vel_tblr';

    %% add parameters and boundary labels for ice velocity
    % pd.n_glen = 1;
    % eps = 0.05; % ratio of membrane stress to driving stress, used to set the ice softness to ensure that membrane stresses are small in this test case, eps rises, A decreases, and viscosity increases
    % pd.A_glen = 1/2/((eps)*pd.rho_i*pd.g*ps.z*ps.x/pd.u_b); % to make membrane stress terms of dimensionless size eps in momentum equation
    % pd.mu_s = mus;
    % [pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity

    [pd,ps,pp,oo] = nevis_update_parameters_ice(pd,ps,pp,oo); % add parameters etc needed to solve for ice velocity
    gg = nevis_label_ice_test(gg, oo); % add boundary labels needed for ice velocity

    if ~isfield(pp,'eps_reg'), pp.eps_reg = 1e-2; end % regularisation on strain rates
    if ~isfield(pp,'Ub_reg'), pp.Ub_reg = 1e-16; end % regularisation on sliding speed (max-based, matches nevis_velocity)
    if ~isfield(pp,'N_slide_reg'), pp.N_slide_reg = 1e-16; end % regularisation on effective pressure (max-based, matches nevis_velocity)
    if ~isfield(pp,'taud_reg'), pp.taud_reg = 1e-16; end % regularisation on basal shear stress [ may not be needed ? ]
    if ~isfield(pp,'C2'), pp.C2 = 0; end % added power-law coefficient in sliding law

    %% load the slipperiness field for the inversion test
    % load(['./data/C_inversion_results.mat'], 'C_hat_dim');

    %% plot grid
    % nevis_plot_grid_ice(gg); return;                    % check to see what grid looks like

    %% initialize variables
    [aa,vv] = nevis_initialize(b,s,gg,pp,oo);         % default initialisation
    % C_dim is dimensional; convert back to nondimensional C
    aa.C = C_dim * (ps.u_b^(1/pp.n_slide) / ps.tau);
    % check if C_hat == aa.C
    
    % if nargin >= 2 && ~isempty(vv_prev)
    %     % --- Continue from previous outer iteration state ---
    %     fprintf('  Using vv_prev from previous outer iteration as initial condition\n');
    %     vv.phi = vv_prev.phi;
    %     vv.hs  = vv_prev.hs;
    %     vv.u   = vv_prev.u;
    %     vv.v   = vv_prev.v;
    %     vv.Sx = vv_prev.Sx;
    %     vv.Sy = vv_prev.Sy;
    %     vv.Ss = vv_prev.Ss;
    %     vv.Sr = vv_prev.Sr;
    %     vv.hb = vv_prev.hb;
    %     vv.pb = vv_prev.pb;
    % else
    %     % --- Fresh start from spinup state ---
    %     fprintf('  Starting from k_f=0.9 initial condition\n');
    %     pd.k_f = 0.9;                                     % percent overburden (k-factor) 
    %     vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
    %     N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
    %     vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)
    % end

    % --- Fresh start from spinup state ---
    % fprintf('  Starting from k_f=0.9 initial condition\n');
    load(['./data/velocity_inverted.mat'], 'vv_hydro');
    % pd.k_f = 0.9;                                     % percent overburden (k-factor) 
    % vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);     % initial pressure  k_f*phi_0
    % N = aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
    % N = vv_hydro.N;
    N = vv_prev.N;
    vv.phi = aa.phi_0 - N;
    vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

    %% initial ice velocity
    % Option 1: use observed velocity
    load([oo.dn '/' 'measures_for_nevis_140km']); 
    dd = measures_for_nevis_140km;
    un = dd.u_obs_dim/(pd.ty*ps.u);
    vn = dd.v_obs_dim/(pd.ty*ps.u);
    % un = gg.nmeanx2*vv_prev.u;
    % vn = gg.nmeany2*vv_prev.v;
    % check if un(gg.ebdy2) agrees with dd.u_obs_dim(gg.ebdy2)/(pd.ty*ps.u)
    un_filled = un(:);
    vn_filled = vn(:);

    % find nans on es2,fs2 and interpolate with the nearest value
    nan_mask_u = isnan(un_filled);
    if any(nan_mask_u)
        valid_idx = find(~nan_mask_u);
        nan_idx = find(nan_mask_u);
        un_filled(nan_idx) = griddata(gg.nx(valid_idx), gg.ny(valid_idx), ...
                                    un_filled(valid_idx), gg.nx(nan_idx), gg.ny(nan_idx), 'nearest');
    end
    nan_mask_v = isnan(vn_filled);
    if any(nan_mask_v)
        valid_idx = find(~nan_mask_v);
        nan_idx = find(nan_mask_v);
        vn_filled(nan_idx) = griddata(gg.nx(valid_idx), gg.ny(valid_idx), ...
                                    vn_filled(valid_idx), gg.nx(nan_idx), gg.ny(nan_idx), 'nearest');
    end
    vv.u = gg.emean2*un_filled;
    vv.v = gg.fmean2*vn_filled;
    vv.u(gg.eout2) = NaN;
    vv.v(gg.fout2) = NaN;
    aa.u_obs = vv.u;
    aa.v_obs = vv.v;

    % vv.u = reshape(gg.emean2*un(:), gg.eI, gg.nJ);
    % vv.v = reshape(gg.fmean2*vn(:), gg.nI, gg.fJ); %% plot initial velocity field
    % nevis_plot_velocity(gg, vv);

    % Option 2: initialise the velocity field using a simplified SSA model
    % N = ones(gg.nIJ,1); 
    u = vv.u;
    v = vv.v;
    [vv.u,vv.v] = nevis_velocity(aa.H,u,v,N,aa,pp,gg,oo);
    % [tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
    % [tau1,tau2,theta] = nevis_principal_stress(Txx,Tyy,gg.nmeanc*Txy);
    aa.u_obs = vv.u;
    aa.v_obs = vv.v;

    %% boundary conditions
    aa.phi_b = max(aa.phi_0,aa.phi_a);                % prescribed boundary pressure at overburden or atmospheric

    %% moulins 
    oo.density_moulins = 1;
    oo.keep_all_moulins = 0;
    oo.random_lakes = 0;
    % load moulin locations
    tmp = load([oo.dn '/' '/nevis_170207a.mat'],'pp');
    pp.ni_m = tmp.pp.ni_m; pp.sum_m = tmp.pp.sum_m; clear tmp;

    %% supraglacial lakes
    lakes = nevis_import_lakes(2022);                     % import lake data
    % multiple lakes from the catalogue
    pp.x_l = lakes.x_m/ps.x;                          % x-coord of lakes
    pp.y_l = lakes.y_m/ps.x;                          % y-coord of lakes
    pp.V_l = 0*lakes.volume_m3/(ps.Q0*ps.t);          % volume of lakes, set to zero for spinup       
    pp.t_drainage = 0.5*pd.td/ps.t*(lakes.drainage_time_start + lakes.drainage_time_end)-17;  % time of lake drainages (assumed to be the middle time of the Gaussian)
    pp.t_duration = 0.5*pd.td/ps.t*(lakes.drainage_duration_days);  % duration of lake drainages
    [pp.ni_l,pp.sum_l] = nevis_lakes(pp.x_l,pp.y_l,gg,oo);          % calculate lake catchments

    %% surface runoff
    load([oo.dn '/' 'runoff_2022_nevis140.mat']);                   % load data for year of interest (previously collated)
    % RACMO distributed input
    oo.surface_runoff = 1;                                          % If set to 1 turns on surface runoff input
                                                                    % 0 is prescribe moulin input with a function
    oo.RACMO_runoff = 0;                                            % If set to 1 turns on RACMO surface runoff data
                                                                    % 0 is a simple sinusoidal function
    oo.distributed_input = 1;                                       % If set to 1 distributes input across the domain
                                                                    % 0 inputs to discrete moulins

    pp.meltE = @(t) (30/1000/pd.td/ps.m)*(1-exp(-t/(20*pd.td/ps.t))); 
    pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2022_nevis140)./ps.m;  % distributed input (m/sec)
    pp.input_function = @(t) 0;

    %% timestep 
    oo.dt = 1/24*pd.td/ps.t; 
    oo.save_timesteps = 1; 
    oo.save_pts_all = 1; 

    % Add GPS station points downstream of the moulin every 5km
    pp.ni_gps = nevis_gps_array([40e3,40e3,20e3,0,-40e3]/ps.x, [-15e3,-5e3,-15e3,-25e3,-30e3]/ps.x, gg, oo); % GPS station points
    oo.pts_ni = [pp.ni_l' pp.ni_m' pp.ni_gps];    
    oo.t_span = [(1:1:365)*pd.td/ps.t];
    % oo.t_span = [(1:1:59)*pd.td/ps.t (59.5:0.001:60.5)*pd.td/ps.t (61:1:120)*pd.td/ps.t];

    %% save initial parameters
    save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','vv','oo');
    [tt,vv,info] = nevis_timesteps(oo.t_span,vv,aa,pp,gg,oo);

    %% expand/update variables
    aa = nevis_inputs(vv.t,aa,vv,pp,gg,oo);
    oo.evaluate_variables = 1; 
    [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); 
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); 
    save([oo.rn, oo.fn],'pp','pd','ps','gg','aa','oo','tt');

    N_new = aa.phi_0 - vv.phi;
    vv_temp = vv;

end