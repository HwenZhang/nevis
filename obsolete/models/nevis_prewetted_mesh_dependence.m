% —— Publication-quality plot of blister radius & thickness —— 
% Enable LaTeX interpreter globally
set(groot, 'DefaultTextInterpreter','latex',                    ...
          'DefaultAxesTickLabelInterpreter','latex',           ...
          'DefaultLegendInterpreter','latex',                  ...
          'DefaultAxesFontSize',14,                            ...
          'DefaultLineLineWidth',1.5);

% Case names and legends
casenames = { ...
    'n1d_prewetted_test_DX3_90625', ...
    'n1d_prewetted_test_DX7_8125', ...
    'n1d_prewetted_test_DX15_625', ...
    'n1d_prewetted_test_DX31_25', ...
    'n1d_prewetted_test_DX62_5', ...
    'n1d_prewetted_test_DX125', ...
    'n1d_prewetted_test_DX250', ...
    'n1d_prewetted_test_DX500' ...
};
legends = {'DX=3.90625', 'DX=7.8125', 'DX=15.625', 'DX=31.25', ...
           'DX=62.5', 'DX=125', 'DX=250', 'DX=500'};

%% Prepare figure for radius and thickness plots
fig = figure('Units','normalized','Position',[0.3,0.3,0.4,0.4]); % Normal size, centered
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Preallocate difference metrics
diffR = nan(numel(casenames),1);
diffh = nan(numel(casenames),1);
deltax = nan(numel(casenames),1);

for i = 1:numel(casenames)
    cn = casenames{i};
    dataFile = fullfile('results',cn,cn);
    if ~isfile([dataFile,'.mat'])
        warning('File %s.mat not found, skipping.', dataFile);
        continue;
    end
    S = load(dataFile,'ps','tt','pd');
    % dimensional time (days), radius R_b (m), thickness h_b (m)
    t   = (S.ps.t/(24*3600)) * [S.tt.t];
    R_b = S.ps.x * [S.tt.Rb];
    h_b = S.ps.hb * [S.tt.hb_max];
    
    if i==1
        tref   = t;
        Rbref  = R_b;
        hbref  = h_b;
    end
    
    diffR(i)  = norm(R_b - Rbref, 'fro');
    diffh(i)  = norm(h_b - hbref, 'fro');
    deltax(i) = 7.8125*2^(i-2);
    
    %— Plot radius —
    ax1 = nexttile(1);
    loglog(ax1, t, R_b, 'DisplayName',legends{i});
    hold(ax1,'on');
    xlabel(ax1,'Time since injection (days)');
    ylabel(ax1,'Blister radius (m)');
    title(ax1,'Blister radius vs. time');
    grid(ax1,'on');
    xlim(ax1,[1e-3,1]);    % adjust as needed
    
    %— Plot thickness —
    ax2 = nexttile(2);
    loglog(ax2, t, h_b, 'DisplayName',legends{i});
    hold(ax2,'on');
    xlabel(ax2,'Time since injection (days)');
    ylabel(ax2,'Blister thickness (m)');
    title(ax2,'Blister thickness vs. time');
    grid(ax2,'on');
    xlim(ax2,[1e-3,1]);    % adjust as needed
end

% Analytical solution (on same axes)
Q0 = 1e7/(pd.td)/1e5;
Rb_a = 2.31*(S.pd.B^2*S.pd.hb_reg1*Q0^5/(49*(12*S.pd.mu)^2)).^(1/17) * (t*S.pd.td).^(7/17);
Hb_a = 15*Q0*t*S.pd.td./(16*Rb_a);

loglog(ax1, t, Rb_a, 'k--','DisplayName','Analytical');
loglog(ax2, t, Hb_a, 'k--','DisplayName','Analytical');

% Legends
legend(ax1,'Location','northwest');
legend(ax2,'Location','northwest');

% Save combined figure
exportgraphics(gcf, 'blister_prewetted_time_series.png', 'Resolution', 300);

%% figure for convergence analysis

fig = figure('Units','normalized','Position',[0.3,0.3,0.4,0.4]); % Normal size, centered
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(1);
% Plot the Frobenius norm of radius difference vs mesh size
plot(deltax, diffh, '-o', 'LineWidth', 2, 'DisplayName', '$||R_b - R_{b,ref}||~(m)$');
hold on;
% plot a delta x^(-2) line for reference
x = linspace(min(deltax), max(deltax), 100);
y = 0.0001 * x.^(2);
plot(x, y, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference: $\Delta x^{2}$');
xlabel('$\Delta x (m)$');
ylabel('$||R_b - R_{b,ref}|| (m)$');
title('convergence of radius');
set(ax1, 'XScale', 'log', 'YScale', 'log');
grid on;
set(ax1, 'FontSize', 14);


ax2 = nexttile(2);
% Plot the Frobenius norm of thickness difference vs mesh size
plot(deltax, diffh, '-o', 'LineWidth', 2, 'DisplayName', '$||h_b - h_{b,ref}||~(m)$');
hold on;
% plot a delta x^(-2) line for reference
y = 0.0001 * x.^(2);
plot(x, y, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference: $\Delta x^{2}$');
xlabel('$\Delta x (m)$');
ylabel('$||h_b - h_{b,ref}|| (m)$');
title('convergence of thickness');
set(ax2, 'XScale', 'log', 'YScale', 'log');
grid on;        
set(ax2, 'FontSize', 14);

legend(ax1, 'Location', 'NorthWest');
legend(ax2, 'Location', 'NorthWest');

% Save figure
exportgraphics(gcf, 'blister_prewetted_mesh_dependence.png', 'Resolution', 300);