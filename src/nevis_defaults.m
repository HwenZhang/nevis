function [pd,oo] = nevis_defaults(pd,oo)
% Assign default parameters and options
% Inputs 
%   pd struct of pre-assigned dimensional parameters [optional]
%   oo struct of scales to use for non-dimensionalization [optional]
% Outputs
%   pd struct of parameters
%   oo struct of scaled parameters
%
% 20 July 2014 : taken from hydro_nondimension
% 22 Sept 2014 : changed default melt to 0
% 27 Nov 2024  : Include blister

if nargin<1 || isempty(pd), pd = struct; end
if nargin<2, oo = struct; end

%% Default parameters
if ~isfield(pd,'td'), pd.td = 24*60*60; end                        % seconds in day [s/d]
if ~isfield(pd,'ty'), pd.ty = 24*60*60*365; end                    % seconds in year [s/y]
if ~isfield(pd,'rho_w'), pd.rho_w = 1000; end                      % water density [kg/m^3]
if ~isfield(pd,'rho_i'), pd.rho_i = 910; end                       % ice density [kg/m^3]
if ~isfield(pd,'g'), pd.g = 9.81; end                              % gravitational acceleration [m/s^2]
if ~isfield(pd,'L'), pd.L = 3.35*10^5; end                         % latent heat [J/kg]
if ~isfield(pd,'c'), pd.c = 4.2e3; end                             % specific heat capacity [J/kg/K]
if ~isfield(pd,'gamma_cc'), pd.gamma_cc = 0*7.5e-8; end            % melting point pressure gradient [K/Pa]
if ~isfield(pd,'G'), pd.G = 0.063; end                             % geothermal heat flux [J/s/m^2]

% ice rheology
if ~isfield(pd,'n_Glen'), pd.n_Glen = 3; end                       % exponent for ice rheology
if ~isfield(pd,'A'), pd.A = 6.8*10^(-24); end                      % ice rheological parameter [Pa^(-n)/s]
if ~isfield(pd,'K_s'), pd.K_s = 2*pd.A*pd.n_Glen^(-pd.n_Glen); end % sheet ice rheological parameter [Pa^(-n)/s]
if ~isfield(pd,'K_c'), pd.K_c = pd.K_s; end                        % conduit ice rheological parameter [Pa^(-n)/s]

% flux equation
if ~isfield(pd,'alpha_s'), pd.alpha_s = 3; end                     % cross-sectional area exponent for conduits
if ~isfield(pd,'beta_s'), pd.beta_s = 1; end                       % potential gradient exponent for conduits
if ~isfield(pd,'alpha_c'), pd.alpha_c = 5/4; end                   % sheet thickness exponent for sheet 
if ~isfield(pd,'beta_c'), pd.beta_c = 1/2; end                     % potential gradient exponent for sheet
if ~isfield(pd,'alpha_e'), pd.alpha_e = 3; end                     % sheet thickness exponent for elastic sheet 
if ~isfield(pd,'beta_e'), pd.beta_e = 1; end                       % potential gradient exponent for elastic sheet
if ~isfield(pd,'k_c'), pd.k_c = 0.1; end                           % conduit flux coefficient [m/s/Pa^(1/2)]
if ~isfield(pd,'k_s'), pd.k_s = 0.0001; end                        % sheet flux coefficient [1/Pa/s]
if ~isfield(pd,'k_e'), pd.k_e = 0.0001; end                        % elastic sheet flux coefficient [1/Pa/s]

if ~isfield(pd,'Psi_reg'), pd.Psi_reg = 1e-1; end                  % regularizing potential gradient [Pa/m]
if ~isfield(pd,'h_r'), pd.h_r = 0.1; end                           % roughness height for cavity sheet [m]
if ~isfield(pd,'l_r'), pd.l_r = 10; end                            % roughness length for cavity sheet [m]
if ~isfield(pd,'h_rc'), pd.h_rc = 0; end                           % roughness height for cavity-like conduits [m]
if ~isfield(pd,'S_rc'), pd.S_rc = 0; end                           % max conduit area for cavity-like conduits [m]

if ~isfield(pd,'A_m'), pd.A_m = 10; end                            % moulin cross-sectional area [m^2]
if ~isfield(pd,'l_c'), pd.l_c = 10; end                            % roughness scale for conduit melting [m]
if ~isfield(pd,'u_b'), pd.u_b = 60/pd.ty; end                      % sliding speed [m/s]
if ~isfield(pd,'sigma'), pd.sigma = 0; end                         % englacial void fraction 
if ~isfield(pd,'melt'), pd.melt = 0*((pd.G)/pd.rho_w/pd.L); end    % basal melt rate [m/s]
if ~isfield(pd,'phi_s'), pd.phi_s = -inf; end                      % sea level potential [Pa]
if ~isfield(pd,'gamma_e'), pd.gamma_e = 1; end                     % exponent for elastic sheet
if ~isfield(pd,'c_e_power'), pd.c_e_power = 0; end                 % depth scale for power law elastic sheet [m]
if ~isfield(pd,'N_reg'), pd.N_reg = 10^3; end                      % regularization for ice pressure in power law elastic sheet [Pa]
if ~isfield(pd,'c_e_reg2'), pd.c_e_reg2 = 0; end                   % thickness rate of change for regularization 2 of elastic sheet [m/Pa]
if ~isfield(pd,'N_reg2'), pd.N_reg2 = pd.N_reg; end                % effective pressure cut off for regularization 2 of elastic sheet [Pa]

if ~isfield(pd,'c_e_log'), pd.c_e_log = 0; end                     % depth scale for log regularization of elastic sheet [m] [ obsolete ]
if ~isfield(pd,'c_e_reg1'), pd.c_e_reg1 = 0; end                   % depth scale for log regularization 1 of elastic sheet [m] [ obsolete ]
if ~isfield(pd,'N_reg1'), pd.N_reg1 = pd.N_reg; end                % effective pressure cut off for regularization 1 of elastic sheet [Pa] [ obsolete ]
if ~isfield(pd,'c_e_reg3'), pd.c_e_reg3 = 0; end                   % depth scale for log regularization 3 of elastic sheet [m] [ obsolete ]
if ~isfield(pd,'N_reg3'), pd.N_reg3 = pd.N_reg; end                % effective pressure cut off for regularization 3 of elastic sheet [Pa] [ obsolete ]
if ~isfield(pd,'V_m_reg'), pd.V_m_reg = 0; end                     % regularizing volume on top of moulins [m^3] [ obsolete ]
if ~isfield(pd,'p_w_reg'), pd.p_w_reg = 1*pd.rho_w*pd.g; end       % pressure range for regularization on top of moulins [Pa] [ obsolete ]
if ~isfield(pd,'sigma_log'), pd.sigma_log = 0; end                 % regularizing storage [ obsolete ]
if ~isfield(pd,'N_sigma'), pd.N_sigma = pd.N_reg; end              % effective pressure below which regularizing storage hits in [Pa] [ obsolete ]
if ~isfield(pd,'hb_reg1'), pd.hb_reg1 = 1e-5; end                  % regularizing thickness of blister [m] [ obsolete ]
if ~isfield(pd,'hb_reg2'), pd.hb_reg2 = 1e-2; end                  % regularizing thickness of blister [m] [ obsolete ]

if ~isfield(pd,'p_a_reg'), pd.p_a_reg = 1*pd.rho_w*pd.g; end       % pressure range for regularization of ? 
if ~isfield(pd,'E_lapse'), pd.E_lapse = 60/1000/pd.td/10^3; end    % surface melt lapse rate [m/s/m]
if ~isfield(pd,'E_amp'), pd.E_amp = 0; end                         % diurnal input amplitude 
if ~isfield(pd,'t_spr'), pd.t_spr = 135*pd.td; end                 % seasonal runoff function spring time [s]
if ~isfield(pd,'t_aut'), pd.t_aut = 244*pd.td; end                 % seasonal runoff function autumn time [s]
if ~isfield(pd,'delta_t'), pd.delta_t = 21*pd.td; end              % timescale [s]
if ~isfield(pd,'s_m'), pd.s_m = 500; end                           % peak runoff location [m]
if ~isfield(pd,'r_m'), pd.r_m = 25/1000/pd.td; end                 % base precipitation [m/s]

% blister
if ~isfield(pd,'Q0'), pd.Q_0 = 1e3; end
if ~isfield(pd,'mu'), pd.mu = 1.0e+3; end                          % water viscosity (Pa s)
if ~isfield(pd,'Ye'), pd.Ye = 8.8e+9; end                          % Young's modulus (Pa)
if ~isfield(pd,'kappa_b'), pd.kappa_b = 1e-9; end                  % relaxation coeff ()
if ~isfield(pd,'S_crit'), pd.S_crit = 0.1; end                     % critical cross section (m^2), below which there is no leakage to the drainage system
if ~isfield(pd,'c0'), pd.c0 = 0.0; end
if ~isfield(pd,'deltap_reg'), pd.deltap_reg = 5.0e7; end           % regularization parameter pressure difference (Pa)

if ~isfield(pd,'alpha_b'), pd.alpha_b = 1/(10*pd.td); end          % relaxation rate of the blister (s^-1)
if ~isfield(pd,'alpha_dh'), pd.alpha_dh = 0; end                   % d (relaxation rate) / d h (m^-1 s^-1)
if ~isfield(pd,'alpha_ds'), pd.alpha_ds = 0; end                   % d (relaxation rate) / d S (m^-2 s^-1)
% if ~isfield(pd,'alpha_dh'), pd.alpha_dh = 1.1574e-5; end           % d (relaxation rate) / d h (m^-1 s^-1)
% if ~isfield(pd,'alpha_ds'), pd.alpha_ds = 1.5595e-5; end           % d (relaxation rate) / d S (m^-2 s^-1)

%% Default options
% method options
if ~isfield(oo,'no_channels'), oo.no_channels = 0; end
if ~isfield(oo,'no_sheet'), oo.no_sheet = 0; end
if ~isfield(oo,'noXi'), oo.noXi = 1; end
if ~isfield(oo,'combine_sheet'), oo.combine_sheet = 1; end
if ~isfield(oo,'mean_perms'), oo.mean_perms = 1; end
if ~isfield(oo,'modified_mean_perms'), oo.modified_mean_perms = 0; end % whether to use modified mean permeability
if ~isfield(oo,'input_gaussian'), oo.input_gaussian = 0; end
if ~isfield(oo,'relaxation_term'), oo.relaxation_term = 0; end      % 0 is alpha hb, 1 is alpha deltap hb
if ~isfield(oo,'initial_condition'), oo.initial_condition = 0; end  % 0 is using default file 0365.mat, 1 is using steady-state drainage system
if ~isfield(oo,'include_ice'), oo.include_ice = 0; end              % couple to ice flow
if ~isfield(oo,'include_blister'), oo.include_blister = 1; end      % couple to blister
if ~isfield(oo,'include_pressure'), oo.include_pressure = 1; end    % couple to pressure 
if ~isfield(oo,'verb'), oo.verb = 1; end                            % display if adding or removing Dirichlet boundary points
% if ~isfield(oo,'method'), oo.method = 'lrtb_vel'; end % [ obsolete ]
% if ~isfield(oo,'include_moulins'), oo.include_moulins = 1; end  % [ obsolete ]
% if ~isfield(oo,'include_diag'), oo.include_diag = 1; end  % [ obsolete ]

% input options
if ~isfield(oo,'surface_runoff'), oo.surface_runoff = 0; end   % use surface runoff: 1 for surface runoff, 0 for direct moulin input
if ~isfield(oo,'RACMO_runoff'), oo.RACMO_runoff = 0; end       % use pp.runoff function(t) for runoff
if ~isfield(oo,'distributed_input'), oo.distributed_input = 0; end
if ~isfield(oo,'include_lake'), oo.include_lake = 1; end

% saving options
if ~isfield(oo,'save_timesteps'), oo.save_timesteps = 0; end
if ~isfield(oo,'save_pts_all'), oo.save_pts_all = 0; end
if ~isfield(oo,'code'), oo.root = ''; end
if ~isfield(oo,'root'), oo.root = ''; end
if ~isfield(oo,'fn'), oo.fn = 'var'; end

% plotting options
if ~isfield(oo,'reversey'), oo.reversey = 0; end
if ~isfield(oo,'halfcmap'), oo.halfcmap = 1; end

% timestepping options
if ~isfield(oo,'change_timestep'), oo.change_timestep = 1; end
if ~isfield(oo,'adjust_boundaries'), oo.adjust_boundaries = 0; end
if ~isfield(oo,'Tol_F'), oo.Tol_F = 1e-3; end
if ~isfield(oo,'dt_max'), oo.dt_max = 1e6; end
if ~isfield(oo,'dt_min'), oo.dt_min = 1e-8; end
if ~isfield(oo,'check_Fs'), oo.check_Fs = 1; end
if ~isfield(oo,'Tol_Fs'), oo.Tol_Fs = oo.Tol_F*[1 1 .1 .1 .1 .1]; end

end