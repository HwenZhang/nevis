function [pd,ps,pp,oo] = nevis_update_parameters_ice(pd,ps,pp,oo)
% Add to existing parameter structs to enable solving for ice velocity
% Inputs 
%   pd struct of dimensional parameters
%   ps struct of scales to use for non-dimensionalization 
%   pp struct of scaled parameters
%   oo struct of options [optional]
% Outputs
%   pd struct of dimensional parameters
%   ps struct of scales
%   pp struct of scaled parameters
%
% 15 October 2025

if nargin<4, oo = struct; end

%% Default scales and dimensional parameters (if not already in pd,ps)
ps.u = ps.u_b;
ps.tau = pd.rho_i*pd.g*ps.z^2/ps.x;
ps.N = pd.rho_i*pd.g*ps.z;

pd.u_b = 100/pd.ty;
pd.n_glen = 3;                                              % exponent in ice flow law
pd.A_glen = 1e-24;                                          % coefficient in ice flow law
pd.n_slide = 1;                                             % exponent in sliding law tau_b ~ C u_b^(1/n) / ( 1 + (C/mu)^n u_b/N^n )^(1/n)
pd.C = ps.tau/ps.u^(1/pd.n_slide);                     % coefficient in sliding law tau_b ~ C u_b^(1/n) / ( 1 + (C/mu)^n u_b/N^n )^(1/n)
pd.C2 = 0*ps.tau/ps.u^(1/pd.n_slide);                       % second coefficient in sliding law tau_b ~ C2 u_b^(1/n)
pd.mu_s = 0.25;                                             % Coulomb friction coefficient (make infinite for Weertman sliding)
ps.eta = 1/2*pd.A_glen^(-1/pd.n_glen)*(ps.u/ps.x)^(1/pd.n_glen-1);

%% Dimensionless parameters
pp.n_glen = pd.n_glen;             % exponent for ice viscosity
pp.A_glen = 1;                     % coefficient in ice flow law
pp.n_slide = pd.n_slide;           % exponent in sliding law
pp.C = pd.C*ps.u^(1/pd.n_slide)/ps.tau; % coefficient in sliding law 
pp.C2 = pd.C2*ps.u^(1/pd.n_slide)/ps.tau; % second coefficient in sliding law 
pp.mu = pd.mu_s*ps.N/ps.tau;       % dimensionless Coulomb friction coefficient mu*[N]/[tau] 
pp.c60 = 1;                             
pp.c61 = 1;
pp.c62 = ps.eta*ps.u/pd.rho_i/pd.g/ps.z/ps.x; % dimensionless coefficient in front of extensional stress terms, defined as \varepsilon
pp.c63 = pd.G/pd.rho_w/pd.L/ps.m;     % dimensionless coefficient of melt rate from geothermal heating
pp.c64 = ps.tau*ps.u/pd.rho_w/pd.L/ps.m;    % dimensionless coefficient of melt rate from frictional heating
pp.c65 = ps.u/ps.u_b;             % ratio between sliding speed scale and sliding speed scale in cavity opening term
pp.c66 = 1;

%% Options
oo.include_ice = 1;         % option to include ice velocity solve in nevis_timesteps

end

