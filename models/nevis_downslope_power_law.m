%% Plot blister radius power-law for multiple cases
% Define list of case names
casenames = { 'n1d_Vconst_theta0_1_kh0_ks0_mu5e0_V1e7_DX31_25', ...
    'n1d_Vconst_theta0_1_kh0_ks0_mu5e0_V1e7_DX62_5', ...
    'n1d_Vconst_theta0_1_kh0_ks0_mu5e0_V1e7_DX125', ...
    'n1d_Vconst_theta0_1_kh0_ks0_mu5e0_V1e7_DX250', ...
    % 添加更多 case 名称
};

% casenames = { 'n1d_Vconst_theta0_01_kh0_ks0_mu5e0_V1e7_DX31_25', ...
%     'n1d_Vconst_theta0_01_kh0_ks0_mu5e0_V1e7_DX62_5', ...
%     'n1d_Vconst_theta0_01_kh0_ks0_mu5e0_V1e7_DX125', ...
%     'n1d_Vconst_theta0_01_kh0_ks0_mu5e0_V1e7_DX250'
% };

legends = {'DX31.25', 'DX62.5', 'DX125', 'DX250'};  % Legends for the cases
% Plot setup
figure('Position',[100,100,400,400]);
hold on;
grid on;
% Fixed injection time (days)
inject_time = 51;

% Loop over each case
for i = 1:numel(casenames)
    casename = casenames{i};
    % Load only necessary data
    dataFile = fullfile('./results', casename, casename);
    load(dataFile, 'ps', 'pd', 'tt');
    
    % Compute dimensional time (days) and blister radius (m)
    t   = (ps.t/(24*60*60)) * [tt.t];
    R_b = ps.x * [tt.Rb];
    
    % Filter to growth phase
    growth = t >= inject_time;
    t0 = t(growth) - inject_time;
    R0 = R_b(growth);
    
    % Plot data curve
    loglog(t0, R0, '-', 'LineWidth', 2, 'DisplayName', legends{i});
end

% Add analytical reference curve (optional)
% Uncomment and adjust if needed:
xf = (0.75)^(2/3) * (pd.rho_w*pd.g*0.1/pd.mu)^(1/3)*(1e7/1e5).^(2/3)*((t - inject_time)*pd.td).^(1/3);
loglog(t - inject_time, xf, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');

% Finalize plot
xlabel('Time (days)');
ylabel('Blister front (m)');
% xlim([0, 500]);   % Set x-axis limits
xticks([1e-2, 1e-1, 1, 10, 100, 1e3, 1e4]); % Set x-axis ticks
title('Blister front vs time with \theta = 0.1');
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
set(gca, 'Box', 'on');
legend('show', 'Location', 'NorthWest');

% Save the figure
% saveas(gcf, './blister_power_law_multiple.png');