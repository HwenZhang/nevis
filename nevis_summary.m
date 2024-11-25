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
figure(1); clf;
ax(1) = subplot(3,1,1);
    plot(t,Q_in,'k',t,Q_out,'r',t,m,'b--',t,E,'b',LineWidth=1.5);
    hold on; 
    plot(t,Q_outQ,'r--', LineWidth=1.5);
    plot(t,Q_outQ+Q_outq,'r:',LineWidth=1.5)
    xlabel('t [ d ]');
    ylabel('Q [ m^3/s ]');
    legend('Q_{in}','Q_{out}','m','E');
    ylim([0 max(Q_out)]);

ax(2) = subplot(3,1,2); 
    if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])  
    plot(t,pts_N,'-',LineWidth=1.5);
    end
    hold on;
    plot(t,N,'-','color',0.6*[1 1 1],LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average N [ MPa ]');
    
ax(3) = subplot(3,1,3); 
    plot(t,hs./A,'b-',t,he./A,'k-',t,1e3*S./A,'r-',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average h [ m ]');
    legend('h_{cav}','h_{el}','S')