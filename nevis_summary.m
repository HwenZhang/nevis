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
p_w = ps.phi/1e6*[tt.pwb];     % dimensional hydraulic potential at the lake (MPa)

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
tmin=150; tmax=300; % time range for the plot
f1 = figure(1); clf;
ax(1) = subplot(5,1,1);
    plot(t,Q_in,'k-',t,Q_out,'r-',t,m,'b-',t,Q_outQ,'k--',LineWidth=1.5);
    hold on; 
    % plot(t,Q_outQ,'k--', LineWidth=1.5);
    % plot(t,Q_outQ+Q_outq,'k:',LineWidth=1.5)
    xlabel('t [ d ]');
    ylabel('Q [ m^3/s ]');
    legend('Q_{in}','Q_{out}','m','Q_{outQ}','NumColumns',2);
    ylim([0 max(Q_out)]);
    xlim([tmin tmax])
    % ylim([0 2000])
    grid on
ax(2) = subplot(5,1,2); 
    if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])  
    plot(t,pts_N,'-',LineWidth=1.5);
    end
    hold on;
    plot(t,N,'-','color',0.6*[1 1 1],LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('N [ MPa ]');
    h=legend('N at the blister','Average N');
    h.Location='southeast';
    xlim([tmin tmax])
    grid on    
ax(3) = subplot(5,1,3); 
    yyaxis left
    plot(t,hs./A,'b-',t,he./A,'b--',LineWidth=1.5);
    xlabel('t [ d ]');
    ylabel('Average h [ m ]');
    h=legend('h_{cav}','h_{el}','NumColumns',2);
    h.Location='northeast';
    ylim([0.02 0.03])

    yyaxis right
    plot(t,1e3*S./A,'r-',LineWidth=1.5);
    ylabel('Average S [ m ]');
    xlim([tmin tmax])
    % ylim([0 0.05])
    grid on
ax(4) = subplot(5,1,4); 
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
    xlim([tmin tmax])
    % legend('h_{cav}','h_{el}','S')
    grid on

ax(5) = subplot(5,1,5); 
    p1 = plot(t,p_b,'k-',LineWidth=1.5); % blister pressure + \rho_i g h
    hold on
    p2 = plot(t,pts_N,'b-.',LineWidth=1.5);
    p3 = plot(t,p_w,'r-.',LineWidth=1.5); 
    
    xlim([tmin tmax])
    ylim([0 5])

    xlabel('t [ d ]');
    ylabel('N [ MPa ]');
    legend('p_b+\rho_i g h','N','p_w')
    grid on

set(f1,'Position',[680 235 660 743])
% exportgraphics(gcf,'time_series.jpg','Resolution',600)
% saveas(f1,'time_series.jpeg','jpeg')

img = getframe(gcf);
imwrite(img.cdata, ['./figures/' oo.casename, '.png']);