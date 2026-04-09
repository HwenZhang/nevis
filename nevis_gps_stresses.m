%% Compute viscous, elastic, and total stress components at lake and GPS nodes
%  Outputs a .mat file with time series of SSA viscous stress tensor
%  components, Kirchhoff-Love basal-fiber elastic bending stresses derived
%  from blister thickness, their summed total stresses, principal stresses,
%  basal shear stress, and driving stresses at lake and GPS nodes.
clc; clear; close all;

%% Settings
casename = 'n2d_regional_racmo_eps1e_02_kappa5e_11_mu5e0_partition5e_01_k01e_01_spinup';
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];
oo.rn = [oo.root,oo.results,oo.fn];
oo.code = '../nevis/src';
filepath = [oo.rn,'/'];
addpath(oo.code);
formatSpec = '%04d';

%% Time range
tmin_yr = 0.0;
tmax_yr = tmin_yr + 1;
tmin = tmin_yr * 365;  % days
tmax = tmax_yr * 365;

tspan_d = (ps.t/pd.td) * oo.t_span;
index = find(tspan_d >= tmin & tspan_d <= tmax);
if isempty(index)
    error('No frames found in the specified time range [%.1f, %.1f] days', tmin, tmax);
end
fprintf('Found %d frames in time range [%.1f, %.1f] days\n', length(index), tmin, tmax);

%% Extraction point info
has_lakes = isfield(pp, 'ni_l') && ~isempty(pp.ni_l);
has_gps = isfield(pp, 'ni_gps') && ~isempty(pp.ni_gps);

if ~has_lakes && ~has_gps
    error('No lake or GPS node indices found in pp.');
end

ni_l = [];
n_lakes = 0;
lake_x_km = [];
lake_y_km = [];
if has_lakes
    ni_l = pp.ni_l(:);                         % node indices of lake locations
    n_lakes = length(ni_l);
    lake_x_km = (ps.x/1e3) * gg.nx(ni_l);     % lake x-coordinates in km
    lake_y_km = (ps.x/1e3) * gg.ny(ni_l);     % lake y-coordinates in km
end

ni_gps = [];
n_gps = 0;
gps_x_km = [];
gps_y_km = [];
if has_gps
    ni_gps = pp.ni_gps(:);                    % node indices of GPS stations
    n_gps = length(ni_gps);
    gps_x_km = (ps.x/1e3) * gg.nx(ni_gps);    % GPS x-coordinates in km
    gps_y_km = (ps.x/1e3) * gg.ny(ni_gps);    % GPS y-coordinates in km
end

%% Scaling factors
stress_membrane = ps.eta * ps.u / ps.x;     % membrane stress scale [Pa]
stress_basal    = pd.rho_i * pd.g * ps.z^2 / ps.x;  % basal/driving stress scale [Pa]
if isfield(pd, 'nu')
    poisson_ice = pd.nu;
else
    poisson_ice = 0.33;
end
stress_elastic_scale = pd.Ye * ps.z * ps.hb / (2 * (1 - poisson_ice^2) * ps.x^2);

%% Pre-allocate output arrays
n_frames = length(index);
t_days   = zeros(1, n_frames);

% Lake fields [Pa]
tauxx_lake = zeros(n_lakes, n_frames);
tauyy_lake = zeros(n_lakes, n_frames);
tauxy_lake = zeros(n_lakes, n_frames);
tauxx_elastic_lake = zeros(n_lakes, n_frames);
tauyy_elastic_lake = zeros(n_lakes, n_frames);
tauxy_elastic_lake = zeros(n_lakes, n_frames);
tauxx_total_lake = zeros(n_lakes, n_frames);
tauyy_total_lake = zeros(n_lakes, n_frames);
tauxy_total_lake = zeros(n_lakes, n_frames);
Txx_lake = zeros(n_lakes, n_frames);
Tyy_lake = zeros(n_lakes, n_frames);
Txy_lake = zeros(n_lakes, n_frames);
sigma1_lake = zeros(n_lakes, n_frames);
sigma2_lake = zeros(n_lakes, n_frames);
sigma1_elastic_lake = zeros(n_lakes, n_frames);
sigma2_elastic_lake = zeros(n_lakes, n_frames);
sigma1_total_lake = zeros(n_lakes, n_frames);
sigma2_total_lake = zeros(n_lakes, n_frames);
t1_lake = zeros(n_lakes, n_frames);
t2_lake = zeros(n_lakes, n_frames);
t1_elastic_lake = zeros(n_lakes, n_frames);
t2_elastic_lake = zeros(n_lakes, n_frames);
t1_total_lake = zeros(n_lakes, n_frames);
t2_total_lake = zeros(n_lakes, n_frames);
tau_b_lake = zeros(n_lakes, n_frames);
taudx_lake = zeros(n_lakes, n_frames);
taudy_lake = zeros(n_lakes, n_frames);
N_lake = zeros(n_lakes, n_frames);
U_lake = zeros(n_lakes, n_frames);

% GPS fields [Pa]
tauxx_gps = zeros(n_gps, n_frames);
tauyy_gps = zeros(n_gps, n_frames);
tauxy_gps = zeros(n_gps, n_frames);
tauxx_elastic_gps = zeros(n_gps, n_frames);
tauyy_elastic_gps = zeros(n_gps, n_frames);
tauxy_elastic_gps = zeros(n_gps, n_frames);
tauxx_total_gps = zeros(n_gps, n_frames);
tauyy_total_gps = zeros(n_gps, n_frames);
tauxy_total_gps = zeros(n_gps, n_frames);
Txx_gps = zeros(n_gps, n_frames);
Tyy_gps = zeros(n_gps, n_frames);
Txy_gps = zeros(n_gps, n_frames);
sigma1_gps = zeros(n_gps, n_frames);
sigma2_gps = zeros(n_gps, n_frames);
sigma1_elastic_gps = zeros(n_gps, n_frames);
sigma2_elastic_gps = zeros(n_gps, n_frames);
sigma1_total_gps = zeros(n_gps, n_frames);
sigma2_total_gps = zeros(n_gps, n_frames);
t1_gps = zeros(n_gps, n_frames);
t2_gps = zeros(n_gps, n_frames);
t1_elastic_gps = zeros(n_gps, n_frames);
t2_elastic_gps = zeros(n_gps, n_frames);
t1_total_gps = zeros(n_gps, n_frames);
t2_total_gps = zeros(n_gps, n_frames);
tau_b_gps = zeros(n_gps, n_frames);
taudx_gps = zeros(n_gps, n_frames);
taudy_gps = zeros(n_gps, n_frames);
N_gps = zeros(n_gps, n_frames);
U_gps = zeros(n_gps, n_frames);

%% Loop over timesteps
for i_idx = 1:n_frames
    i_t = index(i_idx);
    if mod(i_idx, 50) == 0 || i_idx == 1
        fprintf('Processing frame %d / %d ...\n', i_idx, n_frames);
    end

    % Load timestep
    vva = load([filepath num2str(i_t, formatSpec)], 'vv');
    vva = vva.vv;
    aa = nevis_inputs(vva.t, aa, vva, pp, gg, oo);

    % Time in days
    t_days(i_idx) = vva.t * ps.t / pd.td;

    % Effective pressure (dimensional, Pa)
    N_field = ps.phi * (aa.phi_0 - vva.phi);

    % Compute all stress components
    [tauxx, tauyy, tauxy, Txx, Tyy, Txy, tau_b, taudx, taudy, ...
     sigma1, sigma2, t1, t2] = nevis_stresses(aa.H, vva.u, vva.v, ...
                                               aa.phi_0 - vva.phi, aa, pp, gg, oo);

    % Ice speed at nodes [m/yr]
    uxn = gg.nmeanx2(:,gg.es2) * vva.u(gg.es2);
    vyn = gg.nmeany2(:,gg.fs2) * vva.v(gg.fs2);
    Un  = sqrt(uxn.^2 + vyn.^2);

    Txy_n = gg.nmeanc(:,gg.cs2) * Txy(gg.cs2);

    % Driving stress: defined on edges, interpolate to nodes for point extraction
    taudx_n = gg.nmeanx2(:,gg.es2) * taudx(gg.es2);
    taudy_n = gg.nmeany2(:,gg.fs2) * taudy(gg.fs2);
    tauxy_n = gg.nmeanc(:,gg.cs2) * tauxy(gg.cs2);

    % Kirchhoff-Love basal-fiber bending stress from blister thickness.
    hb_xx = gg.nddx * (gg.eddx * vva.hb);
    hb_yy = gg.nddy * (gg.fddy * vva.hb);
    hb_xy = 0.5 * (gg.cddy * (gg.eddx * vva.hb) + gg.cddx * (gg.fddy * vva.hb));
    hb_xy_n = gg.nmeanc(:,gg.cs2) * hb_xy(gg.cs2);

    curvature_scale = ps.hb / ps.x^2;
    kappa_xx = curvature_scale * hb_xx;
    kappa_yy = curvature_scale * hb_yy;
    kappa_xy = curvature_scale * hb_xy_n;

    H_dim = ps.z * max(aa.H, 0);
    bending_prefactor = pd.Ye * H_dim ./ (2 * (1 - poisson_ice^2));
    shear_prefactor = pd.Ye * H_dim ./ (2 * (1 + poisson_ice));

    tauxx_elastic = -bending_prefactor .* (kappa_xx + poisson_ice * kappa_yy);
    tauyy_elastic = -bending_prefactor .* (kappa_yy + poisson_ice * kappa_xx);
    tauxy_elastic = -shear_prefactor .* kappa_xy;

    [sigma1_elastic, sigma2_elastic, t1_elastic] = ...
        nevis_principal_stress(tauxx_elastic, tauyy_elastic, tauxy_elastic);
    t2_elastic = t1_elastic + pi / 2;

    tauxx_total = stress_membrane * tauxx + tauxx_elastic;
    tauyy_total = stress_membrane * tauyy + tauyy_elastic;
    tauxy_total = stress_membrane * tauxy_n + tauxy_elastic;

    [sigma1_total, sigma2_total, t1_total] = ...
        nevis_principal_stress(tauxx_total, tauyy_total, tauxy_total);
    t2_total = t1_total + pi / 2;

    if has_lakes
        tauxx_lake(:, i_idx) = stress_membrane * tauxx(ni_l);
        tauyy_lake(:, i_idx) = stress_membrane * tauyy(ni_l);
        tauxy_lake(:, i_idx) = stress_membrane * tauxy_n(ni_l);
        tauxx_elastic_lake(:, i_idx) = tauxx_elastic(ni_l);
        tauyy_elastic_lake(:, i_idx) = tauyy_elastic(ni_l);
        tauxy_elastic_lake(:, i_idx) = tauxy_elastic(ni_l);
        tauxx_total_lake(:, i_idx) = tauxx_total(ni_l);
        tauyy_total_lake(:, i_idx) = tauyy_total(ni_l);
        tauxy_total_lake(:, i_idx) = tauxy_total(ni_l);
        Txx_lake(:, i_idx) = stress_membrane * Txx(ni_l);
        Tyy_lake(:, i_idx) = stress_membrane * Tyy(ni_l);
        Txy_lake(:, i_idx) = stress_membrane * Txy_n(ni_l);
        sigma1_lake(:, i_idx) = stress_membrane * sigma1(ni_l);
        sigma2_lake(:, i_idx) = stress_membrane * sigma2(ni_l);
        sigma1_elastic_lake(:, i_idx) = sigma1_elastic(ni_l);
        sigma2_elastic_lake(:, i_idx) = sigma2_elastic(ni_l);
        sigma1_total_lake(:, i_idx) = sigma1_total(ni_l);
        sigma2_total_lake(:, i_idx) = sigma2_total(ni_l);
        t1_lake(:, i_idx) = t1(ni_l);
        t2_lake(:, i_idx) = t2(ni_l);
        t1_elastic_lake(:, i_idx) = t1_elastic(ni_l);
        t2_elastic_lake(:, i_idx) = t2_elastic(ni_l);
        t1_total_lake(:, i_idx) = t1_total(ni_l);
        t2_total_lake(:, i_idx) = t2_total(ni_l);
        tau_b_lake(:, i_idx) = stress_basal * tau_b(ni_l);
        taudx_lake(:, i_idx) = stress_basal * taudx_n(ni_l);
        taudy_lake(:, i_idx) = stress_basal * taudy_n(ni_l);
        N_lake(:, i_idx) = N_field(ni_l);
        U_lake(:, i_idx) = ps.u_b * pd.ty * Un(ni_l);
    end

    if has_gps
        tauxx_gps(:, i_idx) = stress_membrane * tauxx(ni_gps);
        tauyy_gps(:, i_idx) = stress_membrane * tauyy(ni_gps);
        tauxy_gps(:, i_idx) = stress_membrane * tauxy_n(ni_gps);
        tauxx_elastic_gps(:, i_idx) = tauxx_elastic(ni_gps);
        tauyy_elastic_gps(:, i_idx) = tauyy_elastic(ni_gps);
        tauxy_elastic_gps(:, i_idx) = tauxy_elastic(ni_gps);
        tauxx_total_gps(:, i_idx) = tauxx_total(ni_gps);
        tauyy_total_gps(:, i_idx) = tauyy_total(ni_gps);
        tauxy_total_gps(:, i_idx) = tauxy_total(ni_gps);
        Txx_gps(:, i_idx) = stress_membrane * Txx(ni_gps);
        Tyy_gps(:, i_idx) = stress_membrane * Tyy(ni_gps);
        Txy_gps(:, i_idx) = stress_membrane * Txy_n(ni_gps);
        sigma1_gps(:, i_idx) = stress_membrane * sigma1(ni_gps);
        sigma2_gps(:, i_idx) = stress_membrane * sigma2(ni_gps);
        sigma1_elastic_gps(:, i_idx) = sigma1_elastic(ni_gps);
        sigma2_elastic_gps(:, i_idx) = sigma2_elastic(ni_gps);
        sigma1_total_gps(:, i_idx) = sigma1_total(ni_gps);
        sigma2_total_gps(:, i_idx) = sigma2_total(ni_gps);
        t1_gps(:, i_idx) = t1(ni_gps);
        t2_gps(:, i_idx) = t2(ni_gps);
        t1_elastic_gps(:, i_idx) = t1_elastic(ni_gps);
        t2_elastic_gps(:, i_idx) = t2_elastic(ni_gps);
        t1_total_gps(:, i_idx) = t1_total(ni_gps);
        t2_total_gps(:, i_idx) = t2_total(ni_gps);
        tau_b_gps(:, i_idx) = stress_basal * tau_b(ni_gps);
        taudx_gps(:, i_idx) = stress_basal * taudx_n(ni_gps);
        taudy_gps(:, i_idx) = stress_basal * taudy_n(ni_gps);
        N_gps(:, i_idx) = N_field(ni_gps);
        U_gps(:, i_idx) = ps.u_b * pd.ty * Un(ni_gps);
    end
end

%% Save results
outfile = ['./data/gps_stresses_' casename '.mat'];
save(outfile, ...
    'casename', 't_days', ...
    'ni_l', 'n_lakes', 'lake_x_km', 'lake_y_km', ...
    'ni_gps', 'n_gps', 'gps_x_km', 'gps_y_km', ...
    'tauxx_lake', 'tauyy_lake', 'tauxy_lake', ...
    'tauxx_elastic_lake', 'tauyy_elastic_lake', 'tauxy_elastic_lake', ...
    'tauxx_total_lake', 'tauyy_total_lake', 'tauxy_total_lake', ...
    'Txx_lake', 'Tyy_lake', 'Txy_lake', ...
    'sigma1_lake', 'sigma2_lake', ...
    'sigma1_elastic_lake', 'sigma2_elastic_lake', ...
    'sigma1_total_lake', 'sigma2_total_lake', ...
    't1_lake', 't2_lake', ...
    't1_elastic_lake', 't2_elastic_lake', ...
    't1_total_lake', 't2_total_lake', ...
    'tau_b_lake', 'taudx_lake', 'taudy_lake', ...
    'N_lake', 'U_lake', ...
    'tauxx_gps', 'tauyy_gps', 'tauxy_gps', ...
    'tauxx_elastic_gps', 'tauyy_elastic_gps', 'tauxy_elastic_gps', ...
    'tauxx_total_gps', 'tauyy_total_gps', 'tauxy_total_gps', ...
    'Txx_gps', 'Tyy_gps', 'Txy_gps', ...
    'sigma1_gps', 'sigma2_gps', ...
    'sigma1_elastic_gps', 'sigma2_elastic_gps', ...
    'sigma1_total_gps', 'sigma2_total_gps', ...
    't1_gps', 't2_gps', ...
    't1_elastic_gps', 't2_elastic_gps', ...
    't1_total_gps', 't2_total_gps', ...
    'tau_b_gps', 'taudx_gps', 'taudy_gps', ...
    'N_gps', 'U_gps', ...
    'stress_membrane', 'stress_basal', 'stress_elastic_scale', 'poisson_ice');
fprintf('Saved GPS stress data to: %s\n', outfile);
fprintf('  %d lake nodes x %d frames\n', n_lakes, n_frames);
fprintf('  %d GPS stations x %d frames\n', n_gps, n_frames);
fprintf('  Time range: %.1f to %.1f days\n', t_days(1), t_days(end));
