% function nevis_plot_summary(tt,vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]
% plot summary of time series in struct tt
%
% 21 August 2014

t = (ps.t/(24*60*60))*[tt.t];  % dimensional time series (days)
Q_b_in = pd.Q_0*[tt.Qb_in];        % dimensional influx (m^3/s)
Q_b_out = ps.V/ps.t*[tt.Qb_out];   % dimensional influx (m^3/s)
Q_in = ps.Q*[tt.Q_in];         % dimensional influx (m^3/s)
Q_out = ps.Q*[tt.Q_out];       % dimensional outflux (m^3/s)
m = (ps.m*ps.x^2)*[tt.m];      % dimensional melting rate (m^3/s)
E = (ps.m*ps.x^2)*[tt.E];      % dimensional source terms  (m^3/s)
Q_outQ = ps.Q*[tt.Q_outQ];     % dimensional channel outflux (m^3/s)
Q_outq = ps.Q*[tt.Q_outq];     % dimensional sheet outflux (m^3/s)
V_b = ps.V*[tt.Vb];            %
R_b = ps.R*[tt.Rb];            %

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
% dimensional pressure p_b+\rho_i g h (MPa)
p_b = 0.125*pd.E_e/1e6/pd.alpha_b/(1-pd.nu^2).*V_b./(R_b).^3 + ps.phi/1e6*aa.phi_0(pp.ni_l);

%% make plot
tmin=180; tmax=240; % time range for the plot
f1 = figure(1); clf;
ax(1) = subplot(6,1,1);
    plot(t,Q_b_in,'b-',t,Q_b_out,'r-',t,Q_out,'k-',LineWidth=1.5);
    hold on; 
    xlabel('t [ d ]');
    ylabel('Q [ m^3/s ]');
    legend('Q_{b,in}','Q_{b,out}','Q_{out}','NumColumns',2);
    text(0.025,0.8,'(a) flux','Units','normalized')

    xlim([tmin tmax])
    ylim([0 2000])
    grid on
    grid minor

ax(2) = subplot(6,1,2); 

    plot(t,N,'-','color',0.6*[1 1 1],LineWidth=1.5);
    hold on;
    
    xlabel('t [ d ]');
    ylabel('N [ MPa ]');
    h=legend('Average N');
    text(0.025,0.8,'(b) effective pressure','Units','normalized')
    h.Location='southeast';
    xlim([tmin tmax])
    grid on    
    grid minor

ax(3) = subplot(6,1,3); 
    yyaxis left
    plot(t,hs./A,'b-',t,he./A,'k-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average h [ m ]');
    % ylim([0.02 0.03])
    text(0.025,0.8,'(c)h and S','Units','normalized')

    yyaxis right
    plot(t,1e3*S./A,'r-',LineWidth=1.5);
    ylabel('Average S [ m ]');
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
    ylabel('Average h [ m ]');
    % ylim([0.02 0.03])
    text(0.025,0.8,'(d)h and S at the blister','Units','normalized')

    yyaxis right
    plot(t,Sx_b,'r-',LineWidth=1.5);
    ylabel('Average S [ m ]');
    xlim([tmin tmax])
    h=legend('h_{cav}','S','NumColumns',2);
    h.Location='northeast';
    % ylim([0 0.05])
    grid on
    grid minor

ax(5) = subplot(6,1,5); 
    yyaxis left
    plot(t,V_b,'b-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('V [ m^3 ]');
    text(0.025,0.8,'(d) blister','Units','normalized')

    yyaxis right
    plot(t,R_b,'r-',LineWidth=1.5);
    hold on

    xlabel('t [ d ]');
    ylabel('R [ m ]');
    xlim([tmin tmax])
    % ylim([-200 1000])
    grid on
    grid minor

ax(6) = subplot(6,1,6); 
    p1 = plot(t,p_b,'k-',LineWidth=1.5); % blister pressure + \rho_i g h
    hold on
    % p2 = plot(t,pts_N,'b-.',LineWidth=1.5);
    p3 = plot(t,p_w,'r-.',LineWidth=1.5); 
    
    xlim([tmin tmax])
    ylim([0 50])
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