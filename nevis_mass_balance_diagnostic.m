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
V_hs = (ps.x^2 * ps.h) * [tt.hs];
V_he = (ps.x^2 * ps.h) * [tt.he];
if isfield(tt, 'hv')
    V_hv = (ps.x^2 * ps.h) * [tt.hv];
else
    V_hv = zeros(size(V_hs));
    warning('tt.hv not found; using V_hv = 0 in main-storage estimate.');
end
if isfield(tt, 'hm')
    V_hm = (ps.x^2 * ps.h) * [tt.hm];
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

%% Strict main residual from R2 over saved frames
% This uses the same discrete mass-conservation residual solved in backbone:
% F2 = R2(gg.nin), then integrates it over area and dimensionalizes.
strict_t_day = [];
strict_res_main = [];
strict_res_blister = [];
strict_res_total = [];

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

    strict_t_day = nan(numel(frame_files)-1, 1);
    strict_res_main = nan(numel(frame_files)-1, 1);
    strict_res_blister = nan(numel(frame_files)-1, 1);
    strict_res_total = nan(numel(frame_files)-1, 1);

    oo_strict = oo;
    oo_strict.evaluate_jacobian = 0;
    oo_strict.evaluate_variables = 0;
    oo_strict.evaluate_residual = 1;

    for k = 2:numel(frame_files)
        Sp = load([case_dir frame_files(k-1).name], 'vv');
        Sc = load([case_dir frame_files(k).name], 'vv');
        vv_prev = Sp.vv;
        vv_curr = Sc.vv;

        dt_k = vv_curr.t - vv_prev.t;
        aa_k = nevis_inputs(vv_curr.t, aa, vv_curr, pp, gg, oo_strict);

        [~, ~, ~, F2, ~, ~, ~, ~, ~, F8] = nevis_backbone(dt_k, vv_curr, vv_prev, aa_k, pp, gg, oo_strict);

        % R2 has scale hs/t; multiply by area and dimensional hs/t scale.
        strict_res_main(k-1) = (ps.h * ps.x^2 / ps.t) * ...
            sum(F2 .* gg.Dx(gg.nin) .* gg.Dy(gg.nin), 'omitnan');
        strict_res_blister(k-1) = (ps.hb * ps.x^2 / ps.t) * ...
            sum(F8 .* gg.Dx(gg.nin_blister) .* gg.Dy(gg.nin_blister), 'omitnan');
        strict_res_total(k-1) = strict_res_main(k-1) + strict_res_blister(k-1);
        strict_t_day(k-1) = (ps.t / pd.td) * vv_curr.t;
    end
else
    warning('Not enough saved frame files for strict residual check in %s', case_dir);
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

if ~isempty(strict_res_main)
    ns = numel(strict_t_day);
    i0s = max(1, floor((1 - steady_tail_fraction) * ns));
    idx_tail_s = i0s:ns;
    fprintf('\n[Main system, strict from R2]\n');
    fprintf('  samples                    = %d\n', ns);
    fprintf('  mean residual              = %.4e m^3/s\n', mean(strict_res_main(idx_tail_s), 'omitnan'));
    fprintf('  std residual               = %.4e m^3/s\n', std(strict_res_main(idx_tail_s), 0, 'omitnan'));

    fprintf('\n[Blister, strict from R8]\n');
    fprintf('  samples                    = %d\n', ns);
    fprintf('  mean residual              = %.4e m^3/s\n', mean(strict_res_blister(idx_tail_s), 'omitnan'));
    fprintf('  std residual               = %.4e m^3/s\n', std(strict_res_blister(idx_tail_s), 0, 'omitnan'));

    fprintf('\n[Total, strict from R2+R8]\n');
    fprintf('  samples                    = %d\n', ns);
    fprintf('  mean residual              = %.4e m^3/s\n', mean(strict_res_total(idx_tail_s), 'omitnan'));
    fprintf('  std residual               = %.4e m^3/s\n', std(strict_res_total(idx_tail_s), 0, 'omitnan'));
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
if ~isempty(strict_res_main)
    plot(strict_t_day, strict_res_main, 'k--', 'LineWidth', 1.2);
    legend('approx residual', 'strict R2 residual', 'Location', 'best');
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
if ~isempty(strict_res_main)
    plot(strict_t_day, strict_res_blister, 'k--', 'LineWidth', 1.2);
    legend('approx residual', 'strict R8 residual', 'Location', 'best');
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
if ~isempty(strict_res_main)
    plot(strict_t_day, strict_res_total, 'k--', 'LineWidth', 1.2);
    legend('approx residual', 'strict R2+R8 residual', 'Location', 'best');
end
xlabel('t [day]'); ylabel('m^3/s'); title('Total residual');
grid on;

sgtitle(['Mass-balance diagnostics: ' casename], 'Interpreter', 'none');
