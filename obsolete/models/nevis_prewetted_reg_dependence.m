% —— Publication-quality plot of blister radius & thickness —— 
% Enable LaTeX interpreter globally
set(groot, 'DefaultTextInterpreter','latex',                    ...
          'DefaultAxesTickLabelInterpreter','latex',           ...
          'DefaultLegendInterpreter','latex',                  ...
          'DefaultAxesFontSize',14,                            ...
          'DefaultLineLineWidth',1.5);

% Case names and legends
casenames = { ...
    'n1d_prewetted_h01e_6', ...
    'n1d_prewetted_h01e_5', ...
    'n1d_prewetted_h01e_4', ...
    'n1d_prewetted_h01e_3', ...
    'n1d_prewetted_h01e_2', ...
    'n1d_prewetted_h01e_1', ...
};
legends = {'h0=10^{-6}', 'h0=10^{-5}', 'h0=10^{-4}', 'h0=10^{-3}', 'h0=10^{-2}', 'h0=10^{-1}'};

%% Prepare figure for radius and thickness plots
fig = figure('Units','normalized','Position',[0.3,0.3,0.4,0.4]); % Normal size, centered
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Preallocate difference metrics
hend = nan(numel(casenames),1);
hreg = nan(numel(casenames),1);

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

    hreg(i) = 1e-7*10^(i);  % Regularizing thickness for this case
    hend(i) = h_b(end-600);
    Rend(i) = R_b(end-600);
    disp(num2str(t(end-600)));
    
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
% exportgraphics(gcf, 'blister_prewetted_time_series.png', 'Resolution', 300);

figure();
loglog(hreg, hend, 'o-');
hold on;
loglog(hreg, 1e-1*hreg.^(-1/17), 'k--');