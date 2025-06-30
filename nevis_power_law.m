%% Plot blister radius power-law for multiple cases (cleaned)
% Define list of case names
casenames = { ...
    'n1d_Vconst_theta0_kh0_ks0_mu5e0_V1e7_DX62_5_test', ...
    'n1d_Vconst_theta0_kh0_ks0_mu5e0_V1e7_DX125_test', ...
    'n1d_Vconst_theta0_kh0_ks0_mu5e0_V1e7_DX250_test', ...
    'n1d_Vconst_theta0_kh0_ks0_mu5e0_V1e7_DX500_test'
};

casenames = { ...
    'n1d_test_DX31_25', ...
    'n1d_test_DX62_5', ...
    'n1d_test_DX125', ...
    'n1d_test_DX250', ...
    'n1d_test_DX500', ...
};

casenames = { ...
    'test_DX31_25', ...
    'test_DX62_5', ...
    'test_DX125', ...
    'test_DX250', ...
    'test_DX500', ...
};

legends = {'DX=31.25','DX=62.5', 'DX=125', 'DX=250', 'DX=500'};
% legends = {'DX=62.5', 'DX=125', 'DX=250', 'DX=500'};
% Fixed threshold for growth phase (days)
inject_time = 0;

% Prepare figure
figure('Position',[100,100,800,600]);
hold on;

% Loop over each case and plot
for i = 1:numel(casenames)
    casename = casenames{i};
    dataFile = fullfile('./results', casename, casename);
    if exist([dataFile, '.mat'], 'file')
        load(dataFile, 'ps', 'tt', 'pp');
    else
        warning('File %s not found.', dataFile);
        continue;
    end
    % Compute dimensional time (days) and blister radius (m)
    t   = (ps.t/(24*3600)) * [tt.t];
    R_b = ps.x * [tt.Rb];
    h_b = ps.hb * [tt.hb_max];
    % Select growth phase
    idx = t >= inject_time;
    t0  = t(idx) - inject_time;
    h0  = h_b(idx);
    R0 = R_b(idx);
    % R0 = R0 + (1e-1 - R0(1));
    % Plot on log-log axes
    plot(t0, R0, 'LineWidth', 2, 'DisplayName', legends{i});
    disp(num2str(length(t0)));
end

% analytical solution
% Q0 = 1e8/(pp.t_duration*ps.t)/1e5;
% Rb_a = (8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((t0)*pd.td).^(0.4);
% plot(t0, Rb_a, 'k--', 'LineWidth', 2, 'DisplayName', 'Analytical solution');

% Finalize log-log plot
% set(gca, 'XScale', 'linear', 'YScale', 'linear');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time since injection (days)');
ylabel('Blister radius (m)');
% ylim([1e3, 1e4]);
xlim([0, 1]);
title('Blister radius vs time for multiple cases');
grid on;
set(gca, 'FontSize', 14);
legend('show', 'Location', 'NorthWest');

% Save figure
% saveas(gcf, './blister_power_law_multiple.png');