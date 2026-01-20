function [pd,ps,pp,oo] = nevis_add_parameters_ice(pd,ps,pp,oo)
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
if ~isfield(ps,'u'), ps.u = ps.u_b; end
if ~isfield(ps,'tau'), ps.tau = pd.rho_i*pd.g*ps.z^2/ps.x; end
if ~isfield(ps,'N'), ps.N = pd.rho_i*pd.g*ps.z; end

if ~isfield(pd,'n_glen'), pd.n_glen = pd.n_Glen; end        % exponent in ice flow law
if ~isfield(pd,'A_glen'), pd.A_glen = pd.A; end             % coefficient in ice flow law
if ~isfield(pd,'n_slide'), pd.n_slide = 1; end              % exponent in sliding law tau_b ~ C u_b^(1/n) / ( 1 + (C/mu)^n u_b/N^n )^(1/n)
if ~isfield(pd,'C'), pd.C = ps.tau/ps.u^(1/pd.n_slide); end % coefficient in sliding law tau_b ~ C u_b^(1/n) / ( 1 + (C/mu)^n u_b/N^n )^(1/n)
if ~isfield(pd,'C2'), pd.C2 = 0; end                        % second coefficient in sliding law tau_b ~ C2 u_b^(1/n)
if ~isfield(pd,'mu_s'), pd.mu_s = 0.1; end                  % Coulomb friction coefficient (make infinite for Weertman sliding)
if ~isfield(ps,'eta'), ps.eta = 1/2*pd.A_glen^(-1/pd.n_glen)*(ps.u/ps.x)^(1/pd.n_glen-1); end

%% Dimensionless parameters
pp.n_glen = pd.n_glen;             % exponent for ice viscosity
pp.A_glen = 1;                     % coefficient in ice flow law
pp.n_slide = pd.n_slide;           % exponent in sliding law
pp.C = pd.C*ps.u^(1/pd.n_slide)/ps.tau; % coefficient in sliding law 
pp.C2 = pd.C2*ps.u^(1/pd.n_slide)/ps.tau; % second coefficient in sliding law 
pp.mu = pd.mu_s*ps.N/ps.tau;       % dimensionless Coulomb friction coefficient mu*[N]/[tau] 
pp.c60 = 1;                             
pp.c61 = 1;
pp.c62 = ps.eta*ps.u/pd.rho_i/pd.g/ps.z/ps.x; % dimensionless coefficient in front of extensional stress terms
pp.c63 = pd.G/pd.rho_w/pd.L/ps.m;     % dimensionless coefficient of melt rate from geothermal heating
pp.c64 = ps.tau*ps.u/pd.rho_w/pd.L/ps.m;    % dimensionless coefficient of melt rate from frictional heating
pp.c65 = ps.u/ps.u_b;             % ratio between sliding speed scale and sliding speed scale in cavity opening term
pp.c66 = 1;

%% Options
% oo.include_ice = 1;         % option to include ice velocity solve in nevis_timesteps

end

