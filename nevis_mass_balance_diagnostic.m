%% Compare E computed on nin (interior only) vs ns (interior+boundary) and Q_out
% This script investigates whether the sawtooth oscillation in Q_out is
% caused by the mismatch: E is summed over ns (= nin + nbdy), but Q_out
% only captures flux through the boundary — so E on boundary nodes is
% double-counted or unaccounted for in the water balance.
%
% Usage: set casename below, then run.

clear; close all;
addpath('./src');

%% ==================== User Settings ====================
casename = 'n2d_region_melt_test_Hreg1000_kappa1e_10_ks1e_03_mu1e1';
suffix = '_spinup';  % '_spinup' or '_drainage'
fullname = [casename suffix];
%% ==================== Load Data ====================
load(['./results/' fullname '/' fullname], 'pp','pd','ps','gg','aa','oo','tt');
oo.fn = ['/' fullname];
oo.rn = [oo.root, oo.results, oo.fn];
oo.code = './src';

%% ==================== Extract time series from tt ====================
% Bug fix: use length(tt), not length(oo.t_span). tt has one entry per
% timestep; saved .mat files are only at t_span intervals — different lengths.
nt = length(tt);
t_days = (ps.t / (24*3600)) * [tt.t];   % dimensional time (days)

% Existing non-dimensional quantities from tt
E_ns_nd   = [tt.E];      % sum over ns, non-dim [ps.m*ps.x^2]
Q_out_nd  = [tt.Q_out];  % total outflow from residual, non-dim [ps.Q]
Q_outQ_nd = [tt.Q_outQ]; % channel outflow, non-dim [ps.Q]
Q_outq_nd = [tt.Q_outq]; % sheet outflow, non-dim [ps.Q]
Q_outb_nd = [tt.Q_outb]; % blister outflow, non-dim [ps.Q]
m_ns_nd   = [tt.m];      % basal melt sum over ns, non-dim [ps.m*ps.x^2]

%% ==================== Recompute E on nin and nbdy separately ====================
% Bug fix: no need to load saved .mat files to get aa.E — it only depends on
% time t. Call nevis_inputs with each tt(i).t directly on all timesteps.
% (Using original aa as base; vv argument unused inside nevis_inputs.)
E_nin_nd  = zeros(1, nt);   % E on interior nodes (nin), non-dim [ps.m*ps.x^2]
E_nbdy_nd = zeros(1, nt);   % E on boundary nodes (nbdy), non-dim [ps.m*ps.x^2]
m_nin_nd  = zeros(1, nt);   % basal melt on nin, non-dim [ps.m*ps.x^2]
m_nbdy_nd = zeros(1, nt);   % basal melt on nbdy, non-dim [ps.m*ps.x^2]

for i = 1:nt
    % Reconstruct aa.E at this timestep's time — no file loading needed
    aa_i = nevis_inputs(tt(i).t, aa, struct, pp, gg, oo);

    % Sum E over interior nodes only
    E_nin_nd(i)  = sum(aa_i.E(gg.nin) .* gg.Dx(gg.nin) .* gg.Dy(gg.nin));
    E_nbdy_nd(i) = sum(aa_i.E(gg.nbdy) .* gg.Dx(gg.nbdy) .* gg.Dy(gg.nbdy));

    % Basal melt split (aa.m is static, but compute per-step for consistency)
    m_nin_nd(i)  = sum(aa_i.m(gg.nin) .* gg.Dx(gg.nin) .* gg.Dy(gg.nin));
    m_nbdy_nd(i) = sum(aa_i.m(gg.nbdy) .* gg.Dx(gg.nbdy) .* gg.Dy(gg.nbdy));
end

% Timestep size from time differences (days)
dt_arr_days = [diff(t_days), NaN];

%% ==================== Dimensional conversion ====================
% Dimensional scales:
%   E  [non-dim] * ps.m * ps.x^2  = m^3/s  (volumetric input rate)
%   Q  [non-dim] * ps.Q            = m^3/s  (volumetric flow rate)
%   dt [non-dim] * ps.t/(24*3600) = days
dim_mQ = ps.m * ps.x^2;           % [m^3/s] scale for E and m
dim_Q  = ps.Q;                     % [m^3/s] scale for channel/sheet flux

% Convert all quantities to m^3/s
E_ns_dim   = E_ns_nd   * dim_mQ;   % total surface input (ns) [m^3/s]
E_nin_dim  = E_nin_nd  * dim_mQ;   % surface input, interior only [m^3/s]
E_nbdy_dim = E_nbdy_nd * dim_mQ;   % surface input, boundary only [m^3/s]
m_ns_dim   = m_ns_nd   * dim_mQ;   % basal melt (ns) [m^3/s]
m_nin_dim  = m_nin_nd  * dim_mQ;   % basal melt (nin) [m^3/s]

Q_out_dim  = Q_out_nd  * dim_Q;    % total outflow from residual [m^3/s]
Q_outQ_dim = Q_outQ_nd * dim_Q;    % channel outflow [m^3/s]
Q_outq_dim = Q_outq_nd * dim_Q;    % sheet outflow [m^3/s]

% Source = E + m, both in m^3/s
source_ns_dim  = E_ns_dim  + m_ns_dim;   % source over ns  [m^3/s]
source_nin_dim = E_nin_dim + m_nin_dim;  % source over nin [m^3/s]

% Water balance imbalance (source + outflow; outflow is negative)
imbalance_ns_dim  = source_ns_dim  + Q_out_dim;
imbalance_nin_dim = source_nin_dim + Q_out_dim;

%% ==================== Plot Comparison ====================
figure('Position', [100 100 1200 900]);

% --- Panel 1: E on ns vs nin vs nbdy (m^3/s) ---
subplot(3,1,1);
plot(t_days, E_ns_dim,   'b-',  'DisplayName', 'E (ns = nin+nbdy)'); hold on;
plot(t_days, E_nin_dim,  'r--', 'DisplayName', 'E (nin only)');
plot(t_days, E_nbdy_dim, 'k:',  'DisplayName', 'E (nbdy only)');
plot(t_days, Q_outq_dim + Q_outQ_dim, 'k-',  'DisplayName', 'Q_outq+Q_outQ');
ylim([0 100])
% Bug fix: Q_outq+Q_outQ is in different units from E — plot separately
% (previously mixed ps.Q and ps.m*ps.x^2 on same y-axis)
ylabel('Surface input E (m^3 s^{-1})');
legend('Location', 'best');
title('Surface input E: interior vs boundary nodes (dimensional)');
grid on;

% --- Panel 2: Water balance — source vs outflow (m^3/s) ---
subplot(3,1,2);
plot(t_days, source_ns_dim,  '--', 'LineWidth', 1.5, 'DisplayName', 'runoff+melting (ns)'); hold on;
plot(t_days, source_nin_dim, 'r-', 'LineWidth', 1.5, 'DisplayName', 'runoff+melting (nin)');
% plot(t_days, -Q_out_dim,     'k-',  'LineWidth', 1.5, 'DisplayName', '-Q_{out} (residual)');
plot(t_days, (Q_outQ_dim + Q_outq_dim), 'm-.', 'LineWidth', 1.5, 'DisplayName', 'sheet and channel outflux');
ylim([0 100])
ylabel('Flow rate (m^3 s^{-1})');
legend('Location', 'best');
title('Source terms vs outflow (dimensional)');
grid on;

% --- Panel 3: Water balance imbalance (m^3/s) ---
subplot(3,1,3);
plot(t_days, imbalance_ns_dim,  'b-',  'DisplayName', 'Imbalance with E(ns)'); hold on;
plot(t_days, imbalance_nin_dim, 'r--', 'DisplayName', 'Imbalance with E(nin)');
yline(0, 'k--');
ylabel('Imbalance (m^3 s^{-1})');
xlabel('Time (days)');
legend('Location', 'best');
title('Water balance imbalance: source + Q_{out}  (should be 0 at steady state)');
grid on;

sgtitle(['Water balance: ' strrep(fullname,'_','\_')], 'FontSize', 14);

%% ==================== Print summary (dimensional, m^3/s) ====================
fprintf('\n===== Summary (m^3/s) =====\n');
fprintf('Mean E(ns):            %.4e m^3/s\n', mean(E_ns_dim));
fprintf('Mean E(nin):           %.4e m^3/s\n', mean(E_nin_dim));
fprintf('Mean E(nbdy):          %.4e m^3/s  (%.2f%% of total)\n', ...
        mean(E_nbdy_dim), 100*mean(E_nbdy_dim)/mean(E_ns_dim));
fprintf('Mean m(ns):            %.4e m^3/s\n', mean(m_ns_dim));
fprintf('Mean m(nin):           %.4e m^3/s\n', mean(m_nin_dim));
fprintf('Mean source(ns):       %.4e m^3/s\n', mean(source_ns_dim));
fprintf('Mean source(nin):      %.4e m^3/s\n', mean(source_nin_dim));
fprintf('Mean -Q_out:           %.4e m^3/s\n', mean(-Q_out_dim));
fprintf('Mean -(Q_outQ+Q_outq): %.4e m^3/s\n', mean(-(Q_outQ_dim+Q_outq_dim)));
fprintf('Mean imbalance (ns):   %.4e m^3/s\n', mean(imbalance_ns_dim));
fprintf('Mean imbalance (nin):  %.4e m^3/s\n', mean(imbalance_nin_dim));
fprintf('Std  Q_out:            %.4e m^3/s\n', std(Q_out_dim));
fprintf('Std  (Q_outQ+Q_outq):  %.4e m^3/s\n', std(Q_outQ_dim+Q_outq_dim));

%% ==================== Additional: dt oscillation check ====================
figure('Position', [100 100 1000 400]);
subplot(1,2,1);
plot(t_days(1:end-1), dt_arr_days(1:end-1), '.-');
ylabel('dt (days)'); xlabel('Time (days)');
title('Timestep size evolution (days)');
grid on;

subplot(1,2,2);
plot(t_days, Q_out_dim,  'b.-', 'DisplayName', 'Q_{out} residual (m^3/s)'); hold on;
plot(t_days, Q_outQ_dim + Q_outq_dim, 'r.-', 'DisplayName', 'Q_{outQ}+Q_{outq} flux (m^3/s)');
ylabel('Outflow (m^3 s^{-1})'); xlabel('Time (days)');
legend('Location', 'best');
title('Q_{out}: residual vs flux-based (m^3/s)');
grid on;

sgtitle('Timestep and outflow diagnostics');

%% ==================== Spatial distribution of E at each saved frame ====================
% Grid coordinates in km
xx = (ps.x/1e3) * gg.nx;
yy = (ps.x/1e3) * gg.ny;
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

path_data = [oo.rn '/'];
formatSpec = '%04d';
n_save = length(oo.t_span);

figure('Position', [100 100 800 600]);
for i = 1:n_save
    fn_i = [path_data num2str(i, formatSpec) '.mat'];
    if ~isfile(fn_i), continue; end
    tmp = load(fn_i, 'vv'); vvi = tmp.vv;

    % Reconstruct E at this time
    aa_i = nevis_inputs(vvi.t, aa, struct, pp, gg, oo);
    E_field = aa_i.E * dim_mQ;  % dimensional E [m^3/s/cell] on each node

    % Reshape to 2D grid
    zE = reshape(E_field, gg.nI, gg.nJ);

    % Plot
    pE = pcolor(xx, yy, zE);
    set(pE, 'linestyle', 'none');
    colorbar;
    colormap(parula);
    xlabel('x (km)'); ylabel('y (km)');
    title(sprintf('E (m^3 s^{-1} per cell),  t = %.1f d', vvi.t * ps.t/(24*3600)));
    axis equal; axis tight;
    drawnow;
    pause(0.3);
end