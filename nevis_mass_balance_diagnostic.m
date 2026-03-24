%% Mass-balance diagnostics for NEVIS outputs
% This script checks approximate mass closure in dimensional units (m^3/s)
% using saved time series tt from one simulation case.
%
% Notes:
% 1) Main-system storage uses hs + he + hv + hm + S.
% 2) Qb_dec is treated as leakage from blister to main system.
% 3) Residuals should be interpreted with the above approximations in mind.

clc; clear; close all;

%% User settings
casename = 'n2d_regional_eps1e_02_kappa1e_10_mu2e1_partition2e_01_spinup_test';
steady_tail_fraction = 0.20;   % use last 20% as steady-state window
r2_map_pair = 'last';          % saved-frame pair for R2 map: 'last' or 1..(nFrames-1)

%% Load result file
case_dir = ['./results/' casename '/'];
result_file = [case_dir casename];

if exist(result_file, 'file')
    % use as-is
elseif exist([result_file '.mat'], 'file')
    result_file = [result_file '.mat'];
else
    cands = dir([case_dir '*.mat']);
    if isempty(cands)
        error('No .mat result file found in %s', case_dir);
    end
    % Prefer file matching case name; otherwise use the only/top candidate.
    names = {cands.name};
    k = find(strcmp(names, [casename '.mat']), 1);
    if isempty(k)
        k = 1;
    end
    result_file = [case_dir cands(k).name];
end

load(result_file);  % expects tt, ps, pd, gg

if ~exist('tt', 'var') || isempty(tt)
    error('No tt struct in %s', result_file);
end
if ~exist('ps', 'var') || ~exist('pd', 'var')
    error('ps/pd not found in %s', result_file);
end

%% Time and dimensionalized fluxes
% t_tt is in model nondimensional time; convert to days and seconds.
t_day = (ps.t / pd.td) * [tt.t];
t_sec = t_day * 86400;

% Dimensional source/sink terms [m^3/s]
M_in = (ps.m * ps.x^2) * [tt.m];                  % basal melt input to main system
E_in = (ps.m * ps.x^2) * [tt.E];                  % runoff input (split internally by Reg_N)
Qb_in_raw = [tt.Qb_in];                           % may be nLake-by-nTime
if isvector(Qb_in_raw)
    Qb_in_raw = reshape(Qb_in_raw, 1, []);
elseif size(Qb_in_raw, 1) == numel(tt) && size(Qb_in_raw, 2) ~= numel(tt)
    % Defensive: if concatenation came out nTime-by-nLake, transpose it.
    Qb_in_raw = Qb_in_raw.';
end

% Aggregate external lake input across lakes to get a 1-by-nTime series.
if numel(pd.Q_0) == 1
    Qb_in = pd.Q_0 * sum(Qb_in_raw, 1, 'omitnan');
else
    q0 = pd.Q_0(:);
    if numel(q0) ~= size(Qb_in_raw, 1)
        error('Size mismatch: pd.Q_0 has %d entries but tt.Qb_in has %d lakes.', numel(q0), size(Qb_in_raw, 1));
    end
    Qb_in = (q0.' * Qb_in_raw);
end
Qb_dec = (ps.h * ps.x^2 / ps.t) * [tt.Qb_dec];    % blister leakage to main system

Qout_main = ps.Q * ([tt.Q_outQ] + [tt.Q_outq]);   % main-system outflow (channel + sheet)
Qout_blister = ps.Q0 * [tt.Q_outb];               % blister outflow (blister-scale flux)

% Some runs store blister terms as NaN when blister is inactive.
if all(isnan(Qout_blister))
    warning('Q_outb is all NaN; treating blister outflow as 0 for diagnostics.');
    Qout_blister = zeros(size(Qout_blister));
else
    Qout_blister(isnan(Qout_blister)) = 0;
end
if all(isnan(Qb_dec))
    warning('Qb_dec is all NaN; treating blister leakage as 0 for diagnostics.');
    Qb_dec = zeros(size(Qb_dec));
else
    Qb_dec(isnan(Qb_dec)) = 0;
end
if all(isnan(Qb_in))
    warning('Qb_in is all NaN; treating blister input as 0 for diagnostics.');
    Qb_in = zeros(size(Qb_in));
else
    Qb_in(isnan(Qb_in)) = 0;
end

%% Dimensional storages
% Main-system storage proxy: cavity + elastic + storage + moulin + channel volume.
hs_scale = ps.h;
if isfield(ps, 'hs'), hs_scale = ps.hs; end
he_scale = ps.h;
if isfield(ps, 'he'), he_scale = ps.he; end
hv_scale = ps.h;
if isfield(ps, 'hv'), hv_scale = ps.hv; end
hm_scale = hv_scale;  % hm uses the same scale factor as hv in R2.

V_hs = (ps.x^2 * hs_scale) * [tt.hs];
V_he = (ps.x^2 * he_scale) * [tt.he];
if isfield(tt, 'hv')
    V_hv = (ps.x^2 * hv_scale) * [tt.hv];
else
    V_hv = zeros(size(V_hs));
    warning('tt.hv not found; using V_hv = 0 in main-storage estimate.');
end
if isfield(tt, 'hm')
    V_hm = (ps.x^2 * hm_scale) * [tt.hm];
else
    V_hm = zeros(size(V_hs));
    warning('tt.hm not found; using V_hm = 0 in main-storage estimate.');
end
V_S  = (ps.x * ps.S) * [tt.S];
V_main = V_hs + V_he + V_hv + V_hm + V_S;

% Blister volume
V_blister = (ps.x^2 * ps.hb) * [tt.Vb];

% Time derivatives [m^3/s]
if numel(t_sec) >= 2
    dV_main_dt = gradient(V_main, t_sec);
    dV_blister_dt = gradient(V_blister, t_sec);
else
    error('Need at least 2 timesteps in tt to compute derivatives.');
end

%% Residuals (positive means apparent excess input)
% Approximate main-system closure:
%   M_in + E_in + Qb_dec - Qout_main - dV_main/dt = 0
res_main = M_in + E_in + Qb_dec - Qout_main - dV_main_dt;

% Approximate blister closure:
%   Qb_in - Qb_dec - Qout_blister - dV_blister/dt = 0
res_blister = Qb_in - Qb_dec - Qout_blister - dV_blister_dt;

% Approximate total closure (Qb_dec cancels):
%   M_in + E_in + Qb_in - Qout_main - Qout_blister - d(V_main+V_blister)/dt = 0
res_total = M_in + E_in + Qb_in - Qout_main - Qout_blister - (dV_main_dt + dV_blister_dt);

%% Exact discrete mass balance over saved frames
% This reconstructs the frame-to-frame storage tendency and the discrete
% source/sink terms used in nevis_backbone, rather than integrating F2 = R2(gg.nin).
exact_t_day = [];
exact_res_main = [];
exact_res_blister = [];
exact_res_total = [];
R2_map_nondim = [];
R2_map_pair_label = '';
R2_map_t_day = NaN;

all_mat = dir([case_dir '*.mat']);
is_frame = false(numel(all_mat), 1);
for i = 1:numel(all_mat)
    is_frame(i) = ~isempty(regexp(all_mat(i).name, '^\d{4}\.mat$', 'once'));
end
frame_files = all_mat(is_frame);

if numel(frame_files) >= 2
    % Ensure source functions are available
    addpath('./src');

    % Sort by frame number
    frame_names = {frame_files.name};
    frame_idx = zeros(numel(frame_names), 1);
    for i = 1:numel(frame_names)
        frame_idx(i) = str2double(frame_names{i}(1:4));
    end
    [~, ord] = sort(frame_idx);
    frame_files = frame_files(ord);

    if ischar(r2_map_pair) || isstring(r2_map_pair)
        if strcmpi(string(r2_map_pair), "last")
            r2_pair_idx = numel(frame_files) - 1;
        else
            error('Unsupported r2_map_pair value: %s', string(r2_map_pair));
        end
    else
        r2_pair_idx = round(r2_map_pair);
    end
    if r2_pair_idx < 1 || r2_pair_idx > numel(frame_files) - 1
        error('r2_map_pair=%d is out of range for %d saved-frame pairs.', r2_pair_idx, numel(frame_files) - 1);
    end

    exact_t_day = nan(numel(frame_files)-1, 1);
    exact_res_main = nan(numel(frame_files)-1, 1);
    exact_res_blister = nan(numel(frame_files)-1, 1);
    exact_res_total = nan(numel(frame_files)-1, 1);

    oo_exact = oo;
    oo_exact.evaluate_jacobian = 0;
    oo_exact.evaluate_variables = 1;
    oo_exact.evaluate_residual = 1;

    oo_state = oo_exact;
    oo_state.evaluate_residual = 0;

    for k = 2:numel(frame_files)
        Sp = load([case_dir frame_files(k-1).name], 'vv');
        Sc = load([case_dir frame_files(k).name], 'vv');
        vv_prev = Sp.vv;
        vv_curr = Sc.vv;

        dt_k = vv_curr.t - vv_prev.t;
        dt_sec_k = ps.t * dt_k;

        aa_prev = nevis_inputs(vv_prev.t, aa, vv_prev, pp, gg, oo_state);
        aa_curr = nevis_inputs(vv_curr.t, aa, vv_curr, pp, gg, oo_exact);

        vv2_prev = nevis_backbone(inf, vv_prev, vv_prev, aa_prev, pp, gg, oo_state);
        [vv2_curr, ~, ~, F2_curr] = nevis_backbone(dt_k, vv_curr, vv_prev, aa_curr, pp, gg, oo_exact);

        area_main = gg.Dx(gg.ns) .* gg.Dy(gg.ns);
        area_blister = gg.Dx(gg.ns_blister) .* gg.Dy(gg.ns_blister);

        V_main_prev_k = ...
            (ps.x^2 * hs_scale) * sum(vv_prev.hs(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x^2 * he_scale) * sum(vv2_prev.he(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x^2 * hv_scale) * sum(vv2_prev.hv(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x^2 * hm_scale) * sum(vv2_prev.hm(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x * ps.S) * ( ...
                sum(vv_prev.Sx(gg.ein) .* (gg.emean(gg.ein,:) * gg.Dx), 'omitnan') + ...
                sum(vv_prev.Sy(gg.fin) .* (gg.fmean(gg.fin,:) * gg.Dy), 'omitnan') + ...
                sum(vv_prev.Ss(gg.cin) .* (gg.cmean(gg.cin,:) * gg.Ds), 'omitnan') + ...
                sum(vv_prev.Sr(gg.cin) .* (gg.cmean(gg.cin,:) * gg.Dr), 'omitnan'));

        V_main_curr_k = ...
            (ps.x^2 * hs_scale) * sum(vv_curr.hs(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x^2 * he_scale) * sum(vv2_curr.he(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x^2 * hv_scale) * sum(vv2_curr.hv(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x^2 * hm_scale) * sum(vv2_curr.hm(gg.ns) .* area_main, 'omitnan') + ...
            (ps.x * ps.S) * ( ...
                sum(vv_curr.Sx(gg.ein) .* (gg.emean(gg.ein,:) * gg.Dx), 'omitnan') + ...
                sum(vv_curr.Sy(gg.fin) .* (gg.fmean(gg.fin,:) * gg.Dy), 'omitnan') + ...
                sum(vv_curr.Ss(gg.cin) .* (gg.cmean(gg.cin,:) * gg.Ds), 'omitnan') + ...
                sum(vv_curr.Sr(gg.cin) .* (gg.cmean(gg.cin,:) * gg.Dr), 'omitnan'));

        V_blister_prev_k = (ps.x^2 * ps.hb) * sum(vv_prev.hb(gg.ns_blister) .* area_blister, 'omitnan');
        V_blister_curr_k = (ps.x^2 * ps.hb) * sum(vv_curr.hb(gg.ns_blister) .* area_blister, 'omitnan');

        dV_main_exact_dt = (V_main_curr_k - V_main_prev_k) / dt_sec_k;
        dV_blister_exact_dt = (V_blister_curr_k - V_blister_prev_k) / dt_sec_k;

        regN_curr = 0.5 - 0.5 * tanh(vv2_curr.N / pp.N_reg1);

        M_in_exact = (ps.m * ps.x^2) * sum(aa_curr.m(gg.ns) .* area_main, 'omitnan');
        E_main_exact = (ps.m * ps.x^2) * ...
            sum((1 - regN_curr(gg.ns)) .* aa_curr.E(gg.ns) .* area_main, 'omitnan');
        E_blister_exact = (ps.m * ps.x^2) * ...
            sum(regN_curr(gg.ns_blister) .* aa_curr.E(gg.ns_blister) .* area_blister, 'omitnan');

        Qb_in_exact = ps.Q0 * sum(vv2_curr.Qb_in, 'omitnan');
        Qb_dec_exact = (ps.h * ps.x^2 / ps.t) * vv2_curr.Qb_dec;
        Qout_main_exact = ps.Q * (vv2_curr.Q_outQ + vv2_curr.Q_outq);
        Qout_blister_exact = ps.Q0 * vv2_curr.Q_outb;

        exact_res_main(k-1) = M_in_exact + E_main_exact + Qb_dec_exact - Qout_main_exact - dV_main_exact_dt;
        exact_res_blister(k-1) = Qb_in_exact + E_blister_exact - Qb_dec_exact - Qout_blister_exact - dV_blister_exact_dt;
        exact_res_total(k-1) = M_in_exact + E_main_exact + E_blister_exact + Qb_in_exact ...
            - Qout_main_exact - Qout_blister_exact - (dV_main_exact_dt + dV_blister_exact_dt);
        exact_t_day(k-1) = (ps.t / pd.td) * vv_curr.t;

        if (k - 1) == r2_pair_idx
            R2_map_nondim = nan(gg.nIJ, 1);
            R2_map_nondim(gg.nin) = F2_curr;
            R2_map_nondim(gg.nbdy) = vv2_curr.R_bdy;
            R2_map_pair_label = sprintf('%s -> %s', frame_files(k-1).name, frame_files(k).name);
            R2_map_t_day = exact_t_day(k-1);
        end
    end
else
    warning('Not enough saved frame files for exact residual check in %s', case_dir);
end

%% Steady-tail diagnostics
n = numel(t_day);
i0 = max(1, floor((1 - steady_tail_fraction) * n));
idx_tail = i0:n;

mean_tail = @(x) mean(x(idx_tail), 'omitnan');
std_tail = @(x) std(x(idx_tail), 0, 'omitnan');

fprintf('\n=== Mass-balance diagnostics: %s ===\n', casename);
fprintf('Steady window: last %.0f%% (%d to %d of %d samples)\n', 100 * steady_tail_fraction, i0, n, n);

fprintf('\n[Main system, approx]\n');
fprintf('  mean(M_in + E_in + Qb_dec) = %.4e m^3/s\n', mean_tail(M_in + E_in + Qb_dec));
fprintf('  mean(Qout_main)            = %.4e m^3/s\n', mean_tail(Qout_main));
fprintf('  mean(dV_main/dt)           = %.4e m^3/s\n', mean_tail(dV_main_dt));
fprintf('  mean residual              = %.4e m^3/s\n', mean_tail(res_main));
fprintf('  std residual               = %.4e m^3/s\n', std_tail(res_main));

fprintf('\n[Blister, approx]\n');
fprintf('  mean(Qb_in)                = %.4e m^3/s\n', mean_tail(Qb_in));
fprintf('  mean(Qb_dec + Qout_blister)= %.4e m^3/s\n', mean_tail(Qb_dec + Qout_blister));
fprintf('  mean(dV_blister/dt)        = %.4e m^3/s\n', mean_tail(dV_blister_dt));
fprintf('  mean residual              = %.4e m^3/s\n', mean_tail(res_blister));
fprintf('  std residual               = %.4e m^3/s\n', std_tail(res_blister));

fprintf('\n[Total, approx]\n');
fprintf('  mean(M_in + E_in + Qb_in)  = %.4e m^3/s\n', mean_tail(M_in + E_in + Qb_in));
fprintf('  mean(Qout_main + Qout_blister) = %.4e m^3/s\n', mean_tail(Qout_main + Qout_blister));
fprintf('  mean(dV_total/dt)          = %.4e m^3/s\n', mean_tail(dV_main_dt + dV_blister_dt));
fprintf('  mean residual              = %.4e m^3/s\n', mean_tail(res_total));
fprintf('  std residual               = %.4e m^3/s\n', std_tail(res_total));

if ~isempty(exact_res_main)
    ns = numel(exact_t_day);
    i0s = max(1, floor((1 - steady_tail_fraction) * ns));
    idx_tail_s = i0s:ns;
    fprintf('\n[Main system, exact discrete]\n');
    fprintf('  samples                    = %d\n', ns);
    fprintf('  mean residual              = %.4e m^3/s\n', mean(exact_res_main(idx_tail_s), 'omitnan'));
    fprintf('  std residual               = %.4e m^3/s\n', std(exact_res_main(idx_tail_s), 0, 'omitnan'));

    fprintf('\n[Blister, exact discrete]\n');
    fprintf('  samples                    = %d\n', ns);
    fprintf('  mean residual              = %.4e m^3/s\n', mean(exact_res_blister(idx_tail_s), 'omitnan'));
    fprintf('  std residual               = %.4e m^3/s\n', std(exact_res_blister(idx_tail_s), 0, 'omitnan'));

    fprintf('\n[Total, exact discrete]\n');
    fprintf('  samples                    = %d\n', ns);
    fprintf('  mean residual              = %.4e m^3/s\n', mean(exact_res_total(idx_tail_s), 'omitnan'));
    fprintf('  std residual               = %.4e m^3/s\n', std(exact_res_total(idx_tail_s), 0, 'omitnan'));
end

if ~isempty(R2_map_nondim)
    R2_node_flux = (ps.h * ps.x^2 / ps.t) * R2_map_nondim .* gg.Dx .* gg.Dy;
    r2_nin_l1 = sum(abs(R2_node_flux(gg.nin)), 'omitnan');
    r2_nbdy_l1 = sum(abs(R2_node_flux(gg.nbdy)), 'omitnan');
    r2_l1_total = r2_nin_l1 + r2_nbdy_l1;
    fprintf('\n[R2 map pair: %s]\n', R2_map_pair_label);
    fprintf('  time                       = %.4f day\n', R2_map_t_day);
    fprintf('  L1(|R2| area), nin         = %.4e m^3/s\n', r2_nin_l1);
    fprintf('  L1(|R2| area), nbdy        = %.4e m^3/s\n', r2_nbdy_l1);
    if r2_l1_total > 0
        fprintf('  nin share                  = %.2f %%\n', 100 * r2_nin_l1 / r2_l1_total);
        fprintf('  nbdy share                 = %.2f %%\n', 100 * r2_nbdy_l1 / r2_l1_total);
    end
end

%% Plot diagnostics
figure('Color', 'w', 'Position', [100 100 1400 900]);

subplot(3,2,1);
plot(t_day, M_in + E_in + Qb_dec, 'k-', 'LineWidth', 1.5); hold on;
plot(t_day, Qout_main, 'b-', 'LineWidth', 1.5);
plot(t_day, dV_main_dt, 'r--', 'LineWidth', 1.2);
xlabel('t [day]'); ylabel('m^3/s'); title('Main terms');
legend('inflow proxy', 'Q_{out,main}', 'dV_{main}/dt', 'Location', 'best');
grid on;

subplot(3,2,2);
plot(t_day, res_main, 'm-', 'LineWidth', 1.5);
hold on;
if ~isempty(exact_res_main)
    plot(exact_t_day, exact_res_main, 'k--', 'LineWidth', 1.2);
    legend('approx residual', 'exact discrete residual', 'Location', 'best');
end
xlabel('t [day]'); ylabel('m^3/s'); title('Main residual');
grid on;

subplot(3,2,3);
plot(t_day, Qb_in, 'k-', 'LineWidth', 1.5); hold on;
plot(t_day, Qb_dec + Qout_blister, 'b-', 'LineWidth', 1.5);
plot(t_day, dV_blister_dt, 'r--', 'LineWidth', 1.2);
xlabel('t [day]'); ylabel('m^3/s'); title('Blister terms');
legend('Q_{b,in}', 'Q_{b,dec}+Q_{out,b}', 'dV_b/dt', 'Location', 'best');
grid on;

subplot(3,2,4);
plot(t_day, res_blister, 'm-', 'LineWidth', 1.5);
hold on;
if ~isempty(exact_res_main)
    plot(exact_t_day, exact_res_blister, 'k--', 'LineWidth', 1.2);
    legend('approx residual', 'exact discrete residual', 'Location', 'best');
end
xlabel('t [day]'); ylabel('m^3/s'); title('Blister residual');
grid on;

subplot(3,2,5);
plot(t_day, M_in + E_in + Qb_in, 'k-', 'LineWidth', 1.5); hold on;
plot(t_day, Qout_main + Qout_blister, 'b-', 'LineWidth', 1.5);
plot(t_day, dV_main_dt + dV_blister_dt, 'r--', 'LineWidth', 1.2);
xlabel('t [day]'); ylabel('m^3/s'); title('Total terms');
legend('total inflow proxy', 'total outflow', 'dV_{total}/dt', 'Location', 'best');
grid on;

subplot(3,2,6);
plot(t_day, res_total, 'm-', 'LineWidth', 1.5);
hold on;
if ~isempty(exact_res_main)
    plot(exact_t_day, exact_res_total, 'k--', 'LineWidth', 1.2);
    legend('approx residual', 'exact discrete residual', 'Location', 'best');
end
xlabel('t [day]'); ylabel('m^3/s'); title('Total residual');
grid on;

sgtitle(['Mass-balance diagnostics: ' casename], 'Interpreter', 'none');

if ~isempty(R2_map_nondim)
    xx = (ps.x / 1e3) * gg.nx;
    yy = (ps.x / 1e3) * gg.ny;

    R2_density_dim = reshape((ps.h / ps.t) * R2_map_nondim, gg.nI, gg.nJ);
    R2_node_flux_dim = reshape((ps.h * ps.x^2 / ps.t) * R2_map_nondim .* gg.Dx .* gg.Dy, gg.nI, gg.nJ);
    R2_node_flux_nin = nan(gg.nIJ, 1);
    R2_node_flux_nbdy = nan(gg.nIJ, 1);
    R2_node_flux_nin(gg.nin) = (ps.h * ps.x^2 / ps.t) * R2_map_nondim(gg.nin) .* gg.Dx(gg.nin) .* gg.Dy(gg.nin);
    R2_node_flux_nbdy(gg.nbdy) = (ps.h * ps.x^2 / ps.t) * R2_map_nondim(gg.nbdy) .* gg.Dx(gg.nbdy) .* gg.Dy(gg.nbdy);
    R2_node_flux_nin_dim = reshape(R2_node_flux_nin, gg.nI, gg.nJ);
    R2_node_flux_nbdy_dim = reshape(R2_node_flux_nbdy, gg.nI, gg.nJ);

    figure('Color', 'w', 'Position', [200 120 1400 900]);
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    pcolor(xx, yy, R2_density_dim);
    shading flat;
    axis equal tight;
    xlabel('x [km]');
    ylabel('y [km]');
    title({'R2 field', sprintf('%s at t = %.3f day', R2_map_pair_label, R2_map_t_day)}, 'Interpreter', 'none');
    cb = colorbar;
    cb.Label.String = 'R2 scaled by ps.h/ps.t [m s^{-1}]';

    nexttile;
    pcolor(xx, yy, R2_node_flux_dim);
    shading flat;
    axis equal tight;
    xlabel('x [km]');
    ylabel('y [km]');
    title('Nodal contribution: R2 * area', 'Interpreter', 'none');
    cb = colorbar;
    cb.Label.String = 'R2 * area * ps.h * ps.x^2 / ps.t [m^3 s^{-1}]';

    nexttile;
    pcolor(xx, yy, R2_node_flux_nin_dim);
    shading flat;
    axis equal tight;
    xlabel('x [km]');
    ylabel('y [km]');
    title('Interior only: nin contribution', 'Interpreter', 'none');
    cb = colorbar;
    cb.Label.String = 'nin-only R2 * area [m^3 s^{-1}]';

    nexttile;
    pcolor(xx, yy, R2_node_flux_nbdy_dim);
    shading flat;
    axis equal tight;
    xlabel('x [km]');
    ylabel('y [km]');
    title('Boundary only: nbdy contribution', 'Interpreter', 'none');
    cb = colorbar;
    cb.Label.String = 'nbdy-only R2 * area [m^3 s^{-1}]';
end
