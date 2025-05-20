function [ps,pp,sc] = nevis_nondimension(pd,ps,oo)
% Nondimensionalise parameters in pd with scales in ps
% Inputs 
%   pd struct of dimensional parameters
%   ps struct of scales to use for non-dimensionalization [optional]
%   oo struct of options [optional]
% Outputs
%   ps struct of scales
%   pp struct of scaled parameters
%   sc scales for plot [optional; for compatability with older versions]
%
% 20 July 2014 : taken from hydro_nondimension
         
if nargin<2 || isempty(ps), ps = struct; end 
if nargin<3, oo = struct; end

%% Default scales
if ~isfield(ps,'u_b'), ps.u_b = pd.u_b; end             % sliding speed (m/s) [L/T]
if ~isfield(ps,'x'), ps.x = 10*10^3; end                % x scale (m) [L]
if ~isfield(ps,'z'), ps.z = 1*10^3; end                 % z scale (m) [L]
if ~isfield(ps,'phi'), ps.phi = pd.rho_i*pd.g*ps.z; end % potential scale (Pa) [M L^-1 T^-2]
% changed 17 April 2016 from pd.l_c so that it can be zero
if ~isfield(ps,'l_c'), ps.l_c = pd.l_r; end             % channel roughness scale (m) [L]
if ~isfield(ps,'m'), ps.m = 10*10^(-3)/pd.td; end       % basal source scale (m/s) [L T^-1]
if ~isfield(ps,'h'), ps.h = 0.1; end                    % sheet thickness scale (m) [L]
if ~isfield(ps,'sigma'), ps.sigma = 1e-5; end           % englacial void (m) [L]

if ~isfield(ps,'Psi'), ps.Psi = ps.phi/ps.x; end        % phi gradient scale (Pa/m) [M L^-2 T^-2]
if ~isfield(ps,'hs'), ps.hs = ps.h; end                 % cavity sheet thickness (m) [L]
if ~isfield(ps,'hv'), ps.hv = ps.h; end                 % storage sheet thickness (m) [L]
if ~isfield(ps,'he'), ps.he = ps.h; end                 % elastic sheet thickness  (m) [L]
if ~isfield(ps,'hb'), ps.hb = 10*ps.h; end              % blister sheet thickness (m) [L]
if ~isfield(ps,'qs'), ps.qs = ps.m*ps.x; end            % sheet flux scale (m^2/s) [L^2 T^-1]
if ~isfield(ps,'qe'), ps.qe = ps.qs; end                % elastic flux scale (m^2/s) [L^2 T^-1]
if ~isfield(ps,'Q'), ps.Q = ps.qs*100; end              % channel flux scale (m^3/s) [L^3 T^-1]
if ~isfield(ps,'S')                                     % channel c.sectional area (m^2) [L^2]
    ps.S = (ps.Q/pd.k_c/ps.Psi^pd.beta_c).^(1/pd.alpha_c); 
end 
if ~isfield(ps,'t'), ps.t = ps.h*ps.x/ps.qs; end        % time scale (s) [T]
if ~isfield(ps,'Xi'), ps.Xi = ps.Q*ps.Psi; end          % energy dissipation scale [M L T^-3]

% blister
% derived variables
% if ~isfield(ps,'hb'), ps.hb = ps.h; end                 % steady-state radius of the blister (m)
if ~isfield(ps,'Q0'), ps.Q0 = pd.Q_0; end                 % blister inflow scale (m^2/s) [L^3 T^-1]   
if ~isfield(ps,'qb'), ps.qb = ps.qs; end                  % blister flow scale (m^2/s) [L^3 T^-1]   
if ~isfield(ps,'alpha_b'), ps.alpha_b = pd.alpha_b; end   % relaxation rate (s^-1)
    
%% Dimensionless parameters [ many of these can be chosen to be 1 by suitable choice of scales ]
pp.c1 = ps.hs/ps.h;
pp.c2 = ps.he/ps.h;
pp.c3 = ps.hv/ps.h;
pp.c4 = ps.qs*ps.t/ps.h/ps.x;
pp.c5 = ps.qe*ps.t/ps.h/ps.x;
pp.c6 = ps.m*ps.t/ps.h;
pp.c7 = ps.m*ps.t/ps.hb;
pp.c8 = ps.S/ps.h/ps.x;
pp.c9 = ps.Q*ps.t/ps.h/ps.x^2;
pp.c10 = ps.m*ps.t/ps.h;
pp.c11 = ps.Xi*ps.t/ps.h/ps.x/pd.rho_w/pd.L;
pp.c12 = 1;
pp.c13 = ps.Xi*ps.t/ps.S/pd.rho_i/pd.L;
pp.c14 = pd.K_c*ps.t*ps.phi^pd.n_Glen;
pp.c15 = ps.hs/ps.h;
pp.c16 = ps.u_b*pd.h_r/pd.l_r*ps.t/ps.h;
pp.c17 = ps.hs/max(pd.h_r,eps);
pp.c18 = pd.K_s*ps.t*ps.phi^pd.n_Glen;
pp.c19 = ps.Q*ps.phi/ps.Xi/ps.x;
pp.c20 = ps.qs*ps.l_c*ps.phi/ps.Xi/ps.x;
pp.c21 = pd.k_s*ps.hs^pd.alpha_s*ps.Psi^pd.beta_s/ps.qs; 
pp.c22 = pd.k_e*ps.he^pd.alpha_e*ps.Psi^pd.beta_e/ps.qe;
pp.c23 = pd.k_c*ps.S^pd.alpha_c*ps.Psi^pd.beta_c/ps.Q;
pp.c24 = pd.c_e_power/ps.he;
pp.c25 = ps.sigma/pd.rho_w/pd.g*ps.phi/ps.hv;
pp.c26 = pd.rho_w/pd.rho_i*ps.m*ps.t/ps.h;
pp.c27 = pd.rho_w/pd.rho_i;
pp.c28 = pd.c_e_log/ps.he; 
pp.c29 = pd.c_e_reg1/ps.he; 
pp.c30 = pd.N_reg1/ps.phi; 
pp.c31 = pd.c_e_reg2*ps.phi/ps.he; 
pp.c32 = pd.N_reg2/ps.phi; 
pp.c33 = pd.V_m_reg/max(pd.A_m,eps)*pd.rho_w*pd.g/ps.phi; 
pp.c34 = pd.p_w_reg/ps.phi; 
pp.c35 = ps.u_b*pd.h_rc*ps.t/ps.S;
pp.c36 = ps.S/max(pd.S_rc,eps);
pp.c37 = pd.c_e_reg3/ps.he; 
pp.c38 = pd.N_reg3/ps.phi; 
pp.c39 = pd.sigma_log/ps.he;
pp.c40 = pd.N_sigma/ps.phi;
pp.c41 = pd.rho_w*pd.gamma_cc*pd.c;

% blister dimensionless parameters
pp.c42 = ps.hb^3*ps.Psi/(12*pd.mu*ps.qb);
pp.c43 = ps.Q0*ps.t/ps.x^2/ps.hb;
pp.c44 = pd.alpha_b*ps.hb*ps.t/ps.hb;
pp.c45 = ps.qb*ps.t/ps.x/ps.hb;
pp.c48 = 1;
pp.c49 = pd.B*ps.hb/ps.x^4/ps.phi;
pp.c50 = pd.alpha_b*ps.hb*ps.t/ps.h;
if oo.relaxation_term == 1
    pp.c51 = pd.kappa_b*ps.hb*ps.t*ps.phi/ps.hb/pd.mu;
    pp.c52 = pd.kappa_b*ps.hb*ps.t*ps.phi/ps.h/pd.mu;
    pp.m_l = 1.0;
elseif oo.relaxation_term == 0
    pp.c51 = pd.alpha_b*ps.hb*ps.t/ps.hb;
    pp.c52 = pd.alpha_b*ps.hb*ps.t/ps.h;
    pp.m_l = 0.0;
elseif oo.relaxation_term == 2
    pp.c51 = pd.alpha_b*ps.hb*ps.t/ps.hb;
    pp.c52 = pd.alpha_b*ps.hb*ps.t/ps.h;
    pp.m_l = 0.0;
else
    error("Invalid relaxation term!")
end
pp.kl_s = 1e0;    % leakage dependence on S
pp.kl_h = 0;      % leakage dependence on h
pp.c0 = 0;        

pp.n_Glen = pd.n_Glen;
pp.alpha_s = pd.alpha_s;
pp.beta_s = pd.beta_s;
pp.alpha_c = pd.alpha_c;
pp.beta_c = pd.beta_c;
pp.alpha_e = pd.alpha_e;
pp.beta_e = pd.beta_e;
pp.gamma_e = pd.gamma_e;
pp.r = pd.rho_w/pd.rho_i;
pp.sigma = pd.sigma/ps.sigma;
pp.A_m = pd.A_m/ps.x^2/ps.sigma;  
pp.phi_s = pd.phi_s/ps.phi;
pp.td = pd.td/ps.t;
pp.ty = pd.ty/ps.t;

pp.Psi_reg = pd.Psi_reg/ps.Psi; 
pp.N_reg = pd.N_reg/ps.phi;
%pp.N_reg1 = pd.N_reg1/ps.phi;
pp.p_w_reg = pd.p_w_reg/ps.phi;
pp.p_a_reg = pd.p_a_reg/ps.phi;
pp.hb_reg = pd.hb_reg/ps.hb;

pp.E_lapse = pd.E_lapse/ps.m*ps.z;
pp.E_amp = pd.E_amp; 
pp.melt = pd.melt/ps.m;
pp.u_b = pd.u_b/ps.u_b;
pp.l_c = pd.l_c/ps.l_c;
pp.t_spr = pd.t_spr/ps.t; 
pp.t_aut = pd.t_aut/ps.t;    
pp.delta_t = pd.delta_t/ps.t; 
pp.s_m = pd.s_m/ps.z;
pp.r_m = pd.r_m/ps.m;

%% Scalings for plots 
sc.x = ps.x/1000;   % km
sc.z = ps.z;        % m
sc.S = ps.S;        % m
sc.h = ps.h*100;    % cm
sc.phi = ps.phi/10^6; % MPa
sc.m = ps.m*10^3*pd.td;   % mm / d
sc.qs = ps.qs;      % m^2 / s
sc.qe = ps.qe;      % m^2 / s
sc.t = ps.t/pd.td;  % d
sc.Q = ps.Q;        % m^3 / s

end

