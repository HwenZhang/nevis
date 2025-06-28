%% Combined plot: two sets of cases in a 1x2 subplot layout
% List of case groups and legends
group1_cases = { ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX62_5', ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX125', ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX250', ...
    'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX500'};
group1_leg = {'DX=62.5 m','DX=125 m','DX=250 m','DX=500 m'};

group2_cases = {'n1d_Vconst_theta0_1_kh0_ks0_mu5e0_V1e7_DX62_5', ...
    'n1d_Vconst_theta0_1_kh0_ks0_mu5e0_V1e7_DX125', ...
    'n1d_Vconst_theta0_1_kh0_ks0_mu5e0_V1e7_DX250', ...
    'n1d_Vconst_theta0_01_kh0_ks0_mu5e0_V1e7_DX500'};
group2_leg = {'DX=62.5 m','DX=125 m','DX=250 m','DX=500 m'};

% Injection times for each group
group1_inj = 200;   % days
group2_inj = 51;    % days

% Create figure with two subplots
figure('Position',[100,100,1000,400]);

%% --- Group 1 subplot ---
ax1 = subplot(1,2,1);
hold on; grid on;

for i = 1:numel(group1_cases)
    casename = group1_cases{i};
    file = fullfile('./results', casename, casename);
    if exist([file '.mat'],'file')
        load(file,'ps','tt');
    else
        warning('Missing %s',casename);
        continue;
    end
    t   = (ps.t/(24*3600))*[tt.t];
    R_b = ps.x * [tt.Rb];
    idx = t > group1_inj;
    t0  = t(idx) - group1_inj;
    R0  = R_b(idx);
    loglog(t0,R0,'LineWidth',2,'DisplayName',group1_leg{i});
end

% analytical solution
load(['./results/' 'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX500' '/' 'n1d_0mm_kh0_ks0_mu5e0_V1e8_DX500'])
Q0 = (ps.Q0*ps.t).*pp.V_l./(pp.t_duration*ps.t)/1e5;
% Rb_a = (8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((t0)*pd.td).^(4/10);
% plot(t0, Rb_a, 'k--', 'LineWidth', 2, 'DisplayName', '$(8.42  B H Q^3 / \mu)^{1/10} t^{2/5}$');

Rb_b = 0.5*(8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((t0)*pd.td).^(4/9);
plot(t0, Rb_b, 'k--', 'LineWidth', 2, 'DisplayName', '$(8.42  B H Q^3 / \mu)^{1/10} t^{2/5}$');
Rb_c = 1.5*(8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((t0)*pd.td).^(4/11);
plot(t0, Rb_c, 'k--', 'LineWidth', 2, 'DisplayName', '$(8.42  B H Q^3 / \mu)^{1/10} t^{2/5}$');

set(ax1, 'XScale','log','YScale','log', ...
         'FontSize',14, 'FontName','Times New Roman', ...
         'LineWidth',1.5, 'TickLength',[0.02 0.025]);
xlabel(ax1,'Time since injection (days)','Interpreter','latex');
ylabel(ax1,'Blister radius (m)','Interpreter','latex');
title(ax1,['$Q_{in}=$ ' num2str(Q0,'%.2e') ' $m^2/s$, with varying $\Delta x$'],'Interpreter','latex','FontSize',16);
% % legend('show','Location','southoutside','Orientation','horizontal','NumColumns',2);
l1 = legend(ax1,'show','Location','northwest','NumColumns',1,'Box','on','FontSize',12);
l1.Interpreter = 'latex';


%% --- Group 2 subplot ---
ax2 = subplot(1,2,2);
hold on; grid on;
for i = 1:numel(group2_cases)
    casename = group2_cases{i};
    file = fullfile('./results', casename, casename);
    if exist([file '.mat'],'file')
        load(file,'ps','pd','tt');
    else
        warning('Missing %s',casename);
        continue;
    end
    t   = (ps.t/(24*3600))*[tt.t];
    R_b = ps.x * [tt.Rb];
    idx = t >= group2_inj;
    t0  = t(idx) - group2_inj;
    R0  = R_b(idx);
    loglog(t0,R0,'LineWidth',2,'DisplayName',group2_leg{i});
end
% Analytical curve for group2
xf = (0.75)^(2/3) * (pd.rho_w*pd.g*0.1/pd.mu)^(1/3) * (1e7/1e5)^(2/3) * ((t - group2_inj)*pd.td).^(1/3);
loglog(t - group2_inj,xf,'k--','LineWidth',2,'DisplayName','$(3/4)^{2/3}(\rho_w g \theta / \mu)^{1/3} V^{2/3} t^{1/3}$');

set(ax2, 'XScale','log','YScale','log', ...
         'FontSize',14, 'FontName','Times New Roman', ...
         'LineWidth',1.5, 'TickLength',[0.02 0.025]);
xlabel(ax2,'Time since injection (days)','Interpreter','latex');
ylabel(ax2,'Blister front (m)','Interpreter','latex');
title(ax2,'$V=100~m^2$, $\theta=0.1$, with varying $\Delta x$','Interpreter','latex','FontSize',16);
legend(ax2,'show','Location','northwest','Interpreter','latex','NumColumns',1,'Box','on','FontSize',12);

axis tight;