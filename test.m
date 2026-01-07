figure();clf;

casename = 'n1d_example_drainage_mu1e1_kappa1e_10_eps0_1';
% casename = 'n1d_example_spinup_mu1e1_kappa1e_10_eps0_1';
load(['./results/' casename '/' casename])
tdrain = 365*20+0.2*365;
path = ['./results/' casename '/'];
formatSpec = '%04d';

vva = load([path '0001'], 'vv');
vva = vva.vv;
uxn = gg.nmeanx(:,gg.es2)*vva.u(gg.es2); % x-component of ice velocity at nodes
vyn = gg.nmeany(:,gg.fs2)*vva.v(gg.fs2); % y-component of ice velocity at nodes
Un = ((uxn).^2+(vyn).^2).^(1/2); % ice speed at nodes
Un(gg.nout) = NaN;  
vux = ps.u_b*pd.ty*reshape(uxn,gg.nI,gg.nJ); % x-component of ice velocity at nodes in m/yr
vuy = ps.u_b*pd.ty*reshape(vyn,gg.nI,gg.nJ); % y-component of ice velocity at nodes in m/yr

uil = [tt.pts_u];
u_l = ps.u_b*pd.ty*uil(1,:);

subplot(3,1,1)
centerline_j = round(gg.nJ/2);
xx = gg.nx(:,1)*ps.x/10^3; % x grid in km
vel_centerline = plot(xx(:,1),vux(:,centerline_j), 'r-', 'LineWidth', 2); 
hold on;
txt = text(0.02, 0.95, ['Time = ' num2str(vva.t*ps.t/pd.td,'%.1f') ' days'], 'Units', 'normalized', 'FontSize', 14, 'Color', 'k');
xlabel('x (km)','FontSize',10);
ylabel('Ice velocity u (m/yr)','FontSize',10);
title('Ice velocity profile along centerline','FontSize',10);
set(gca,'FontSize',10);
grid on;
xlim([0 100])
ylim([0 100]);

tmin = 20*365 + 0.18*365;
tmax = tmin+0.25*365;
index = find(oo.t_span>=tmin*pd.td/ps.t & oo.t_span<=tmax*pd.td/ps.t);
t = ps.t/pd.td*[tt.t]-tdrain;
U = ps.u_b*pd.ty*[tt.U];

subplot(3,1,2)
yyaxis left
    pblister = plot(xx(gg.ns),ps.hb*vva.hb(gg.ns),'b-','LineWidth',1.5); 
    ylim([-0.02 0.05]);
    ylabel('h_b [m]')
    title('blister sheet and pb','FontSize',14);
yyaxis right
    pphi = plot(xx(gg.ns),ps.hs*vva.hs(gg.ns),'r-','LineWidth',1.5);
    hold on
    xlabel('x (km)')
    ylabel('hs (m)')
    % ylim([0 10]);
    grid on
    grid minor

subplot(3,1,3)
p1 = plot(t, u_l, LineStyle="-", LineWidth=1.0);
hold on
time = plot(t(index(1)), u_l(index(1)), 'ro', 'MarkerFaceColor', 'r');
ylabel('Ice velocity u (m/yr)','FontSize',10);
xlabel('Time (d)', 'FontSize', 10)
title('Averaged ice velocity ','FontSize',10);
set(gca,'FontSize',10);
grid on;
ylim([20 45]);
xlim([-10 50]);

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

    disp(ps.u_b*pd.ty*mean(((gg.nmeanx(:,gg.es2)*vva.u(gg.es2)).^2+0*(gg.nmeany(:,gg.fs2)*vva.v(gg.fs2)).^2).^(1/2)));
    vux = ps.u_b*pd.ty*reshape(uxn,gg.nI,gg.nJ); % x-component of ice velocity at nodes in m/yr
    vuy = ps.u_b*pd.ty*reshape(vyn,gg.nI,gg.nJ); % y-component of ice velocity at nodes in m/yr

    vel_centerline.YData = (vux(:,centerline_j).^2 + vuy(:,centerline_j).^2).^(1/2);

    time.XData = t(nframe);
    time.YData = u_l(nframe);

    % Update the plots for the next frame
    pblister.YData = ps.hb * vva.hb(gg.ns);
    pphi.YData = ps.hs*vva.hs(gg.ns);
    
    txt.String = ['Time = ' num2str(t(nframe),'%.1f') ' days'];

    refreshdata
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);