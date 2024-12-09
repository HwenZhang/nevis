% Script to animate the 1-dimensiona results from nevis
% By Hanwen Zhang on 25/11/2024

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

figure(2);
fs = 12;
path = '/Users/eart0487/nevis/nevis_1d_example/';
formatSpec = '%04d';
vva = load([path num2str(1,formatSpec)], 'vv');
vva = vva.vv;

v = VideoWriter("test.avi");
open(v)

subplot(511);
ani_phi = (ps.phi/10^6)*vva.phi(gg.ns);
ani_phi0 = (ps.phi/10^6)*aa.phi_0(gg.ns);
ani_hs = ps.h*vva.hs(gg.ns);

% hydraulic potential 
    p1 = plot(ps.x/1e3*gg.nx(gg.ns),ani_phi,'-'); 
    set(p1,'LineWidth',1.5,...
    'Color',[0 0 1]); hold on
    p1.YDataSource = 'ani_phi';
    
    p2 = plot(ps.x/1e3*gg.nx(gg.ns),ani_phi0,'-'); 
    set(p2,'LineWidth',1.5,...
    'Color',[1 0 0]);
    p2.YDataSource = 'ani_phi0';

    xlabel('x (km)');
    ylabel('\phi (MPa)');
    xlim([0 10]) 
    ylim([0 10]) 

    % legend
    h=legend([p1,p2],'$\phi$','$\phi_0$');
    h.Interpreter='latex';
    h.FontSize=fs+2;
    h.Location='northeast';
    
    time = ['t=' num2str(vva.t/pd.td,'%.1f'), ' d'];
    ttext = text(0.1,8,time);
    ttext.FontSize=12;
    grid on

% sheet thickness
subplot(512);
    p3 = plot(ps.x/1e3*gg.nx(gg.ns),ani_hs,'-'); 
    set(p3,'LineWidth',1.5,...
    'Color',[0 0 0]); hold on
    p3.YDataSource = 'ani_hs';

    xlabel('x (km)');
    ylabel('h (m)');
    % ylim([0 20]) 

    % legend
    % h=legend(p3,'$h_s$');
    % h.Interpreter='latex';
    % h.FontSize=fs+2;
    % h.Location='northeast';
    grid on

subplot(513);
    t = (ps.t/(24*60*60))*[tt.t];  % dimensional time series (days)
    V_b = ps.V*[tt.Vb];            % dimensional volume (m^3)
    R_b = ps.R*[tt.Rb];            % dimensional radius (m)
    Q_b = ps.V/ps.t*[tt.Qb];       % dimensional flux (m^3/s)
    p_w = ps.phi/1e6*[tt.pwb];     % dimensional hydraulic potential at the lake (MPa)
    S = ps.x*ps.S*[tt.S];          % % dimensional cross sectional area (m^2)
    % dimensional pressure p_b+\rho_i g h (MPa)
    p_b = 0.125*pd.E_e/1e6/pd.alpha_b/(1-pd.nu^2).*V_b./(R_b).^3 + ps.phi/1e6*aa.phi_0(pp.ni_l);
    
    grid on
    yyaxis left
    p4 = plot(t,V_b,'b-'); 
    set(p4,'LineWidth',1.5,...
    'Color',[0 0 1]); hold on
    xlabel('t [ d ]');
    ylabel('Blister volume [ m^3 ]');
    ylim([0 5e8])
    % pointer
    tm = t(1);
    V_m = V_b(1);
    R_m = R_b(1);

    pp1 = plot([tm tm],[0 5e8],'k--',LineWidth=1.5);
    pp1.XDataSource = 'tm';

    yyaxis right
    p5 = plot(t,R_b,'r-'); 
    set(p5,'LineWidth',1.5,...
    'Color',[1 0 0]);
    xlim([0 max(t)])
    ylim([0 5e3])
    xlabel('t [ d ]');
    ylabel('Blister radius [ m ]');
    

subplot(514);
    A = ps.x^2*sum(gg.Dx.*gg.Dy);
    S_m = 1e3*S(1)./A;
    hs = ps.x^2*ps.h*[tt.hs];      % integrated hs (m^3)
    he = ps.x^2*ps.h*[tt.he];      % integrated he (m^3)

    p6 = plot(t,hs./A,'b-','LineWidth',1.5);
    hold on
    p7 = plot(t,he./A,'k-','LineWidth',1.5);
    p8 = plot(t,1e3*S./A,'r-','LineWidth',1.5);

    xlabel('t [ d ]');
    ylabel('Average h [ m ]');
    xlim([0 max(t)])
    h=legend('h_{cav}','h_{el}','S','NumColumns',3);
    h.Location='northwest';
    grid on

    pp3 = plot([tm tm],[0 5],'k--',LineWidth=1.5);
    pp3.XDataSource = 'tm';

subplot(515);
    p9 = plot(t,p_b,'-'); 
    set(p9,'LineWidth',1.5,'Color',[0 0 0]);
    hold on
    p10 = plot(t,p_w,'--'); 
    set(p10,'LineWidth',1.5,'Color',[1 0 0]);
    
    if isfield(tt,'pts_phi') && ~isempty([tt.pts_phi])  
    plot(t,pts_N,'-.',LineWidth=1.5);
    end
    % plot(t,N,'-.','color',0.6*[1 1 1],LineWidth=1.5);

    xlim([0 max(t)])
    ylim([0 5])
    xlabel('t [ d ]');
    ylabel('pressure (MPa)');
    grid on
    
    pp4 = plot([tm tm],[0 1e1],'k--',LineWidth=1.5);
    pp4.XDataSource = 'tm';


set(gcf,'Position',[100 100 500 600])
F(:,1) = getframe; 

%% animate
for i=1:vv.ti_save
    % each xxxx.mat file contains a solution field vv
    vva = load([path num2str(i,formatSpec)], 'vv');
    vva = vva.vv;

    ani_phi = (ps.phi/10^6)*vva.phi(gg.ns);
    ani_phi0 = (ps.phi/10^6)*aa.phi_0(gg.ns);
    ani_hs = ps.h*vva.hs(gg.ns);
    
    % pointer
    [~,it]=min(abs(vva.t*ps.t/(24*60*60)-t));
    tm = [t(it) t(it)];
    V_m = V_b(it);
    R_m = R_b(it);
    S_m = 1e3*S(it)./A;
    Q_m = ps.V/ps.t*Q_b(it);
    set(ttext,{'string'},{['t=' num2str(vva.t*ps.t/(24*60*60),'%.1f'), ' d']})  %notice the column vector of new values

    refreshdata
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame)
    % pause(0.1);
end

close(v)