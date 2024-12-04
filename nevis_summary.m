% function nevis_plot_summary(tt,vv,aa,pp,ps,gg,oo) % [uncomment to make a function ]
% plot summary of time series in struct tt
%
% 21 August 2014

t = (ps.t/(24*60*60))*[tt.t];  % dimensional time series (days)
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
he = ps.x^2*ps.h*[tt.he];      % integrated he (m^3)

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
f1 = figure(1); clf;
ax(1) = subplot(4,1,1);
    plot(t,Q_in,'k-',t,Q_out,'r-',t,m,'b-',t,Q_outQ,'k--',LineWidth=1.5);
    hold on; 
    % plot(t,Q_outQ,'k--', LineWidth=1.5);
    % plot(t,Q_outQ+Q_outq,'k:',LineWidth=1.5)
    xlabel('t [ d ]');
    ylabel('Q [ m^3/s ]');
    legend('Q_{in}','Q_{out}','m','Q_{outQ}','NumColumns',2);
    ylim([0 max(Q_out)]);
    xlim([80 140])
    % ylim([0 2000])
    grid on
ax(2) = subplot(4,1,2); 
    if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])  
    plot(t,pts_N,'-',LineWidth=1.5);
    end
    hold on;
    plot(t,N,'-','color',0.6*[1 1 1],LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average N [ MPa ]');
    xlim([80 140])
    grid on    
ax(3) = subplot(4,1,3); 
    yyaxis left
    plot(t,hs./A,'b-',t,he./A,'b--',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average h [ m ]');
    h=legend('h_{cav}','h_{el}','NumColumns',2);
    h.Location='northeast';
    ylim([2e-3 5e-3])

    yyaxis right
    plot(t,1e3*S./A,'r-',LineWidth=1.5);
    ylabel('Average S [ m ]');
    xlim([80 140])
    grid on
ax(4) = subplot(4,1,4); 
    yyaxis left
    plot(t,V_b,'b-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('V [ m^3 ]');
    % legend('h_{cav}','h_{el}','S')
    grid on

    yyaxis right
    plot(t,R_b,'r-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('R [ m ]');
    xlim([80 140])
    % legend('h_{cav}','h_{el}','S')
    grid on

% ax(5) = subplot(5,1,5); 
%     plot(t,N_b,'k-',LineWidth=1.5);
%     xlabel('t [ d ]');
%     ylabel('N [ MPa ]');
%     % legend('h_{cav}','h_{el}','S')
%     grid on

set(gcf,'Position',[100 100 500 600])
exportgraphics(gcf,'time_series.jpg','Resolution',600)