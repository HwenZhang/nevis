%% Plot blister radius power-law for multiple cases (cleaned)
% Define list of case names
casenames = { ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX62_5', ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX125', ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX250', ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX500', ...
};
legends = {'DX=62.5', 'DX=125', 'DX=250', 'DX=500'};
% Fixed threshold for growth phase (days)
inject_time = 200;

% Prepare figure
figure('Position',[100,100,800,600]);
hold on;

% Loop over each case and plot
for i = 1:numel(casenames)
    casename = casenames{i};
    dataFile = fullfile('./results', casename, casename);
    if exist([dataFile, '.mat'], 'file')
        load(dataFile, 'ps', 'tt');
    else
        warning('File %s not found.', dataFile);
        continue;
    end
    % Compute dimensional time (days) and blister radius (m)
    t   = (ps.t/(24*3600)) * [tt.t];
    R_b = ps.x * [tt.Rb];
    % Select growth phase
    idx = t > inject_time;
    t0  = t(idx) - inject_time;
    R0  = R_b(idx);
    % Plot on log-log axes
    loglog(t0, R0, 'LineWidth', 2, 'DisplayName', legends{i});
end

% analytical solution
Q0 = 1e7/(pp.t_duration*ps.t)/1e5;
Rb_a = (8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((t0)*pd.td).^(0.4);
plot(t0, Rb_a, 'k--', 'LineWidth', 2, 'DisplayName', 'Analytical solution');

% Finalize log-log plot
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time since injection (days)');
ylabel('Blister radius (m)');
ylim([1e2, 1e4]);
title('Blister radius vs time for multiple cases');
grid on;
set(gca, 'FontSize', 14);
legend('show', 'Location', 'NorthWest');

% Save figure
% saveas(gcf, './blister_power_law_multiple.png');