figure();clf;

casename = 'test';
% casename = 'n1d_example_spinup_mu1e1_kappa1e_10_eps0_1';
load(['./results/' casename '/' casename])
tdrain = 365*20+20;
path = ['./results/' casename '/'];
formatSpec = '%04d';

vva = load([path '0002'], 'vv');
vva = vva.vv;
uxn = gg.nmeanx(:,gg.es2)*vva.u(gg.es2); % x-component of ice velocity at nodes
vyn = gg.nmeany(:,gg.fs2)*vva.v(gg.fs2); % y-component of ice velocity at nodes
Un = ((uxn).^2+(vyn).^2).^(1/2); % ice speed at nodes
Un(gg.nout) = NaN;  
vux = ps.u_b*pd.ty*reshape(uxn,gg.nI,gg.nJ); % x-component of ice velocity at nodes in m/yr
vuy = ps.u_b*pd.ty*reshape(vyn,gg.nI,gg.nJ); % y-component of ice velocity at nodes in m/yr

uil = [tt.pts_u];
u_l = ps.u_b*pd.ty*uil(1,:);

xx = gg.nx(:,1)*ps.x/10^3; % x grid in km
xx_m = gg.nx(:,1)*ps.x; % x grid in km

% analytical solution
etac = 1/(2*pd.A_glen);
vux_analytical = pd.ty*pd.rho_i*pd.g*1e3/(12*etac)*xx_m.*(1-xx_m/1e5).^(1/2);

%% ice velocity profile along the flowline
subplot(3,2,1)
centerline_j = round(gg.nJ/2);
vux_pre = vux(:,centerline_j);
vuy_pre = vuy(:,centerline_j);
vel_centerline = plot(xx(:,1),vux_pre, 'r-', 'LineWidth', 2); 
hold on;
vely_centerline = plot(xx(:,1),vuy_pre, 'b-', 'LineWidth', 2); 
vel_analytical = plot(xx(:,1),vux_analytical, 'k--', 'LineWidth', 2);

txt = text(0.02, 0.90, ['Time = ' num2str(vva.t*ps.t/pd.td,'%.1f') ' days'], ...
    'Units', 'normalized', 'FontSize', 10, 'Color', 'k');
xlabel('x (km)','FontSize',10);
ylabel('Ice velocity u (m/yr)','FontSize',10);
title('Ice velocity along centerline','FontSize',10);
set(gca,'FontSize',10);
grid on;
xlim([0 100])
% ylim([0 500]);

subplot(3,2,2)
centerline_j = round(gg.nJ/2);
vux_pre = vux(:,centerline_j);
delta_vel_centerline = plot(xx(:,1),vux_pre, 'r-', 'LineWidth', 2); 
hold on;
xlabel('x (km)','FontSize',10);
ylabel('\Delta u (m/yr)','FontSize',10);
title('Ice velocity anomalies along centerline','FontSize',10);
set(gca,'FontSize',10);
grid on;
xlim([0 100])
ylim([-15 200]);

%% 
tmin = 0*365;
tmax = tmin + 2.5*365;
index = find(oo.t_span>=tmin*pd.td/ps.t & oo.t_span<=tmax*pd.td/ps.t);
t = ps.t/pd.td*[tt.t]-tdrain;
U = ps.u_b*pd.ty*[tt.U];

subplot(3,2,3)
yyaxis left
    pblister = plot(xx(gg.ns),ps.hb*vva.hb(gg.ns),'b-','LineWidth',1.5); 
    ylim([-0.005 0.75]);
    ylabel('h_b [m]')
    title('blister sheet and pb','FontSize',10);
yyaxis right
    % pphi = plot(xx(gg.ns),ps.hs*vva.hs(gg.ns),'r-','LineWidth',1.5);
    zN = (ps.phi/10^6)*reshape((aa.phi_0-vva.phi),gg.nI,gg.nJ);
    pphi = plot(xx, zN, 'r-', 'LineWidth', 1.5);
    hold on
    xlabel('x (km)')
    ylabel('N (kPa)')
    ylim([-0.5 5]);
    grid on
    grid minor

subplot(3,2,4)
yyaxis left
hb_pre = ps.hb*vva.hb(gg.ns);
zN_pre = (ps.phi/10^6)*reshape((aa.phi_0-vva.phi),gg.nI,gg.nJ);
    delta_pblister = plot(xx(gg.ns),ps.hb*vva.hb(gg.ns),'b-','LineWidth',1.5); 
    ylim([-0.01 0.25]);
    ylabel('\Delta h_b [m]')
    title('blister sheet and pb','FontSize',10);
yyaxis right
    % pphi = plot(xx(gg.ns),ps.hs*vva.hs(gg.ns),'r-','LineWidth',1.5);
    zN = (ps.phi/10^6)*reshape((aa.phi_0-vva.phi),gg.nI,gg.nJ);
    delta_pphi = plot(xx, zN, 'r-', 'LineWidth', 1.5);
    hold on
    xlabel('x (km)')
    ylabel('\Delta N (kPa)')
    ylim([-0.5 1.5]);
    grid on
    grid minor

subplot(3,2,5)    
phs = plot(xx(gg.ns),ps.hs*vva.hs(gg.ns),'b-','LineWidth',1.5); 
hold on
ylim([0.0 0.15]);
ylabel('h_s [m]')
title('cavity sheet','FontSize',10);
grid on
grid minor

subplot(3,2,6)
p1 = plot(t, u_l, LineStyle="-", LineWidth=1.0);
hold on
time = plot(t(index(1)), u_l(index(1)), 'ro', 'MarkerFaceColor', 'r');
ylabel('Ice velocity u (m/yr)','FontSize',10);
xlabel('Time (d)', 'FontSize', 10)
title('Averaged ice velocity ','FontSize',10);
set(gca,'FontSize',10);
grid on;
ylim([min(u_l)-1 max(u_l)+1]);
xlim([-10 400]);

v = VideoWriter(['./presentation/animation/' casename '_timeseries'],'MPEG-4');
v.FrameRate = 1;
open(v);

for nframe = index(1:1:end)
    vva = load([path num2str(nframe,formatSpec)], 'vv');
    vva = vva.vv;
    uxn = gg.nmeanx(:,gg.es2)*vva.u(gg.es2); % x-component of ice velocity at nodes
    vyn = gg.nmeany(:,gg.fs2)*vva.v(gg.fs2); % y-component of ice velocity at nodes
    Un = ((uxn).^2+(vyn).^2).^(1/2); % ice speed at nodes
    Un(gg.nout) = NaN;  
    vux = ps.u_b*pd.ty*reshape(uxn,gg.nI,gg.nJ); % x-component of ice velocity at nodes in m/yr
    vuy = ps.u_b*pd.ty*reshape(vyn,gg.nI,gg.nJ); % y-component of ice velocity at nodes in m/yr
    
    vel_centerline.YData = vux(:,centerline_j);
    vely_centerline.YData = vuy(:,centerline_j);
    delta_vel_centerline.YData = vux(:,centerline_j) - vux_pre;

    time.XData = t(nframe);
    time.YData = u_l(nframe);

    % Update the plots for the next frame
    pblister.YData = ps.hb * vva.hb(gg.ns);
    delta_pblister.YData = ps.hb * vva.hb(gg.ns) - hb_pre;
    phs.YData = ps.hs * vva.hs(gg.ns);

    vN = (ps.phi/10^6)*reshape((aa.phi_0-vva.phi),gg.nI,gg.nJ);
    pphi.YData = vN;
    delta_pphi.YData = vN -zN_pre;

    txt.String = ['Time = ' num2str(t(nframe),'%.1f') ' days'];

    refreshdata
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);