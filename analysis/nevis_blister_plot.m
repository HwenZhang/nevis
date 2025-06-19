%% Import necessary libraries
casename = 'n1d_0mm_cg0_00_a0_kh0_ks0_mu5e0_c1_V1e8';
% casename = 'test_2009_140km_mu2e1_kappa0';
load(['./results/' casename '/' casename])
% oo.fn = ['/',casename];                         % filename (same as casename)
% oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
% oo.code = '../nevis/src';                       % code directory  
path = [oo.rn,'/'];
addpath(oo.code);                               % add path to code
  
dt = oo.dt;
tmin = 0*pd.td/ps.t;
tmax = 1000*pd.td/ps.t;
tmin_d = tmin*ps.t/pd.td; 
tmax_d = tmax*ps.t/pd.td;                        % time range for the plot
t_init = 175; t_end = 400;                      % time range for animation

%% colormap
n = 256; % number of colors
cmap = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1); 
        ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];

%% read in the screenshot at the intial timestep
formatSpec = '%04d';
tframe1 = 800.0*pd.td/ps.t; % time frame for the screenshot
nframe1 = round(tframe1/0.05); % frame number
tframe2 = 1200.0*pd.td/ps.t; % time frame for the screenshot
nframe2 = round(tframe2/0.05); % frame number
tframe3 = 1600.0*pd.td/ps.t; % time frame for the screenshot
nframe3 = round(tframe3/0.05); % frame number

vv1 = load([path num2str(nframe1,formatSpec)], 'vv');
vv1 = vv1.vv;
vv2 = load([path num2str(nframe2,formatSpec)], 'vv');
vv2 = vv2.vv;
vv3 = load([path num2str(nframe3,formatSpec)], 'vv');
vv3 = vv3.vv;

xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

Q0 = (ps.Q0*ps.t).*pp.V_l./(pp.t_duration*ps.t)/1e5;
Rb_a1 = (8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((tframe1*ps.t/pd.td-200.0)*pd.td).^(0.4);
Rb_a2 = (8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((tframe2*ps.t/pd.td-200.0)*pd.td).^(0.4);
Rb_a3 = (8.42*pd.B*1000*Q0^3/pd.mu)^(1/10)*((tframe3*ps.t/pd.td-200.0)*pd.td).^(0.4);

hb1_a = (abs(xx*1e3-5e4)<Rb_a1).*15/16*(Q0*(tframe1*ps.t/pd.td-200)*pd.td)/Rb_a1.*(1-((xx*1e3-5e4)/Rb_a1).^2).^2;
hb2_a = (abs(xx*1e3-5e4)<Rb_a2).*15/16*(Q0*(tframe2*ps.t/pd.td-200)*pd.td)/Rb_a2.*(1-((xx*1e3-5e4)/Rb_a2).^2).^2;
hb3_a = (abs(xx*1e3-5e4)<Rb_a3).*15/16*(Q0*(tframe3*ps.t/pd.td-200)*pd.td)/Rb_a3.*(1-((xx*1e3-5e4)/Rb_a3).^2).^2;
pb1_a = 45*pd.B*(Q0*(tframe1*ps.t/pd.td-200)*pd.td)/(2*Rb_a1^5)/1e6; % pressure in MPa
pb2_a = 45*pd.B*(Q0*(tframe2*ps.t/pd.td-200)*pd.td)/(2*Rb_a2^5)/1e6; % pressure in MPa
pb3_a = 45*pd.B*(Q0*(tframe3*ps.t/pd.td-200)*pd.td)/(2*Rb_a3^5)/1e6; % pressure in MPa

%% plot the blister profile
figure;
subplot(2,1,1);
plot(xx, vv1.hb, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
hold on
plot(xx, vv2.hb, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(xx, vv3.hb, 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
plot(xx, hb1_a, '--','LineWidth', 2, 'Color', [0 0.4470 0.7410]);
plot(xx, hb2_a, '--','LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
plot(xx, hb3_a, '--','LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
xlabel('x (km)');
ylabel('blister sheet thickness (m)');
title('Blister Sheet Thickness');
xlim([40 60]);
ylim([0 0.2]);
grid on;
set(gca, 'FontSize', 14);
legend({'t = 800 d', 't = 1400 d', 't = 2000 d', 'analytical'}, 'Location', 'NorthEast');

subplot(2,1,2);
plot(xx, ps.phi/1e6*vv1.pb-ps.phi/1e6, 'LineWidth', 2,'Color', [0 0.4470 0.7410]);
hold on
plot(xx, ps.phi/1e6*vv2.pb-ps.phi/1e6, 'LineWidth', 2,'Color', [0.8500 0.3250 0.0980]);
plot(xx, ps.phi/1e6*vv3.pb-ps.phi/1e6, 'LineWidth', 2,'Color', [0.9290 0.6940 0.1250]);

% yline(pb1_a, 'Color', [0 0.4470 0.7410], 'LineStyle', '--','LineWidth', 2);
% yline(pb2_a, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '--','LineWidth', 2);
% yline(pb3_a, 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--','LineWidth', 2);

xlabel('x (km)');
ylabel('blister sheet pressure (MPa)');
title('Blister Sheet Pressure');
xlim([40 60]);
ylim([-0.02 0.02]);
grid on;
set(gca, 'FontSize', 14);
% Save the figure
saveas(gcf, ['./results/' casename '/blister_profile.png']);