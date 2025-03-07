% function nevis_plot_summary(tt,vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]
% plot summary of time series in struct tt
%
% 21 August 2014
%% read in the results
casename = oo.casename;
oo.fn = ['/',oo.casename];                      % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
path = [oo.rn,'/'];
load([path casename '.mat'], 'pp','pd','ps','gg','aa','oo','tt');

%% process the variables
t = (ps.t/(24*60*60))*[tt.t];             % dimensional time series (days)
Q_b_in = pd.Q_0*[tt.Qb_in];               % dimensional influx (m^3/s)
Q_b_out = ps.h*ps.x^2/ps.t*[tt.Qb_out];   % dimensional influx (m^3/s)
Q_in = ps.Q*[tt.Q_in];         % dimensional influx (m^3/s)
Q_out = ps.Q*[tt.Q_out];       % dimensional outflux (m^3/s)
m = (ps.m*ps.x^2)*[tt.m];      % dimensional melting rate (m^3/s)
E = (ps.m*ps.x^2)*[tt.E];      % dimensional source terms  (m^3/s)
Q_outQ = ps.Q*[tt.Q_outQ];     % dimensional channel outflux (m^3/s)
Q_outq = ps.Q*[tt.Q_outq];     % dimensional sheet outflux (m^3/s)
h_b = ps.hb*[tt.pts_hb];            %
p_b = ps.phi*[tt.pts_pb];            %

phi = (ps.phi/10^6)*[tt.phi];  % dimensional hydrulic potential (MPa)
N = (ps.phi/10^6)*[tt.N];      % dimensional effective stress (MPa)
hs = ps.x^2*ps.h*[tt.hs];      % integrated hs (m^3)
hs_b = ps.h*[tt.hs_b];         % integrated hs (m^3)
he = ps.x^2*ps.h*[tt.he];      % integrated he (m^3)
p_w = ps.phi/1e6*[tt.pwb];     % dimensional hydraulic potential at the lake (MPa)

Sx_b = ps.S*[tt.Sx_b];
Sy_b = ps.S*[tt.Sy_b];
Ss_b = ps.S*[tt.Ss_b];
Sr_b = ps.S*[tt.Sr_b];
S_ave = (0.25*(Sx_b.^2 + Sy_b.^2 + Ss_b.^2 + Sr_b.^2)).^0.5;

S = ps.x*ps.S*[tt.S];
A = ps.x^2*sum(gg.Dx.*gg.Dy);

if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])    
    pts_phi = (ps.phi/10^6)*[tt.pts_phi];
    pts_hs = ps.hs*[tt.pts_hs];
    pts_N = (ps.phi/10^6)*(aa.phi_0(oo.pts_ni)*[tt.t].^0 - [tt.pts_phi]);
    pts_pw = (ps.phi/10^6)*([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0);
    pts_prat = ([tt.pts_phi]-aa.phi_a(oo.pts_ni)*[tt.t].^0)./...
               (aa.phi_0(oo.pts_ni)*[tt.t].^0-aa.phi_a(oo.pts_ni)*[tt.t].^0);
end

%% make plot
tmin=200; tmax=400; % time range for the plot
f1 = figure(1); clf;
ax(1) = subplot(6,1,1);
    plot(t,Q_b_in,'b-',t,Q_b_out,'r-',LineWidth=1.5);
    hold on;
    plot(t,Q_out,color=[0,0.5,0],LineStyle='-',LineWidth=1.5);
    plot(t,E,color=[0,0,0],LineStyle='--',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Q [ m^3/s ]');
    h=legend('Q_{b,in}','Q_{b,out}','Q_{out}','Q_{in}','NumColumns',2);
    h.Location='southeast';
    text(0.025,0.8,'(a) flux','Units','normalized')

    xlim([tmin tmax])
    set(gca, 'YScale', 'log')
    ylim([1e0 1e3])
    yticks([1e0,1e1,1e2,1e3])
    grid on
    grid minor

ax(2) = subplot(6,1,2); 
    if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])  
    plot(t,pts_N,'-',LineWidth=1.5);
    end
    hold on
    plot(t,N,'-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('N [ MPa ]');
    h = legend('N at the blister','averaged N');
    text(0.025,0.8,'(b) effective pressure','Units','normalized')
    h.Location='southeast';
    xlim([tmin tmax])
    % ylim([0 2])
    grid on    
    grid minor

ax(3) = subplot(6,1,3); 
    yyaxis left
    plot(t,hs./A,'b-',t,he./A,'k-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average h [ m ]');
    ylim([0 0.3])
    text(0.025,0.8,'(c)averaged h and S','Units','normalized')

    yyaxis right
    plot(t,S/A,'r-',LineWidth=1.5);
    ylabel('Channel h [ m ]');
    h=legend('h_{cav}','h_{e}','S','NumColumns',2);
    h.Location='northeast';
    xlim([tmin tmax])
    % ylim([0 0.05])
    grid on
    grid minor

ax(4) = subplot(6,1,4);
    yyaxis left
    plot(t,hs_b,'b-',LineWidth=1.5);    
    % plot(t,hs./A,'b-',t,he./A,'b--',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('h at blister [ m ]');
    % ylim([0.02 0.03])
    text(0.025,0.8,'(d)h and S at the blister','Units','normalized')

    yyaxis right
    plot(t,0.25*(Sx_b+Sy_b+Ss_b+Sr_b),'r-',LineWidth=1.5);
    ylabel('S at blister [ m^2 ]');
    xlim([tmin tmax])
    h=legend('h_{cav}','S','NumColumns',2);
    h.Location='northeast';
    % ylim([0 0.05])
    grid on
    grid minor

ax(5) = subplot(6,1,5); 
    yyaxis left
    plot(t,h_b,'b-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('h [ m ]');
    text(0.025,0.8,'(d) blister thickness and pressure','Units','normalized')

    yyaxis right
    plot(t,p_b/1e6,'r-',LineWidth=1.5);
    hold on

    xlabel('t [ d ]');
    ylabel('p_b [ MPa ]');
    xlim([tmin tmax])
    grid on
    grid minor

ax(6) = subplot(6,1,6); 
    p1 = plot(t,p_b,'k-',LineWidth=1.5); % blister pressure + \rho_i g h
    hold on
    % p2 = plot(t,pts_N,'b-.',LineWidth=1.5);
    p3 = plot(t,p_w,'r-.',LineWidth=1.5); 
    
    xlim([tmin tmax])
    ylim([0 30])
    text(0.025,0.8,'(e) pressure','Units','normalized')
    xlabel('t [ d ]');
    ylabel('N [ MPa ]');
    legend('p_b+\rho_i g h','p_w')
    grid on
    grid minor

set(f1,'Position',[680 235 660 843])
% exportgraphics(gcf,'time_series.jpg','Resolution',600)
% saveas(f1,'time_series.jpeg','jpeg')

img = getframe(gcf);
imwrite(img.cdata, ['./figures/' oo.casename, '.png']);