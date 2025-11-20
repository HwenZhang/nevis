figure();clf;

% casename = 'n2d_moulin1e1_eps1e_02_kappa1e_10_mu1e1_V1e7_drainage';
load(['./results/' casename '/' casename])
tdrain = 365*10+0.2*365;

subplot(2,1,1)
centerline_j = round(gg.nJ/2);
vel_centerline = plot(xx(:,1),vux(:,centerline_j), 'r-', 'LineWidth', 2); 
hold on;
xlabel('x (km)','FontSize',10);
ylabel('Ice velocity u (m/yr)','FontSize',10);
title('Ice velocity profile along centerline','FontSize',10);
set(gca,'FontSize',10);
grid on;
ylim([0 50]);

tmin = 10*365 + 0.19*365;
tmax = tmin+1.0*365;
index = find(oo.t_span>=tmin*pd.td/ps.t & oo.t_span<=tmax*pd.td/ps.t);

subplot(2,1,2)
t = ps.t/pd.td*[tt.t]-tdrain;
U = ps.u_b*pd.ty*[tt.U];

p1 = plot(t, U, LineStyle="-", LineWidth=1.0);
hold on
time = plot(t(index(1)), U(index(1)), 'ro', 'MarkerFaceColor', 'r');
ylabel('Ice velocity u (m/yr)','FontSize',10);
xlabel('Time (d)', 'FontSize', 10)
title('Averaged ice velocity ','FontSize',10);
set(gca,'FontSize',10);
grid on;
ylim([5 20]);

v = VideoWriter(['./presentation/animation/' casename '_timeseries'],'MPEG-4');
v.FrameRate = 1;
open(v);
for nframe = index(1:8:end)
    vva = load([path num2str(nframe,formatSpec)], 'vv');
    vva = vva.vv;
    uxn = gg.nmeanx(:,gg.es2)*vva.u(gg.es2); % x-component of ice velocity at nodes
    vyn = gg.nmeany(:,gg.fs2)*vva.v(gg.fs2); % y-component of ice velocity at nodes
    Un = ((uxn).^2+(vyn).^2).^(1/2); % ice speed at nodes
    Un(gg.nout) = NaN;  

    vux = ps.u_b*pd.ty*reshape(uxn,gg.nI,gg.nJ); % x-component of ice velocity at nodes in m/yr
    vuy = ps.u_b*pd.ty*reshape(vyn,gg.nI,gg.nJ); % y-component of ice velocity at nodes in m/yr

    vel_centerline.YData = vux(:,centerline_j);

    time.XData = t(nframe);
    time.YData = U(nframe);

    refreshdata
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);