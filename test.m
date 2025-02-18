% t = (ps.t/(24*60*60))*[tt.t];  % dimensional time series (days)
% Q_b_in = pd.Q_0*[tt.Qb_in];        % dimensional influx (m^3/s)
% Q_b_out = ps.V/ps.t*[tt.Qb_out];   % dimensional influx (m^3/s)
% Q_b_out_test = ps.V/ps.t*[tt.Qb_out_test];   % dimensional influx (m^3/s)
% Q_in = ps.Q*[tt.Q_in];         % dimensional influx (m^3/s)
% Q_out = ps.Q*[tt.Q_out];       % dimensional outflux (m^3/s)
% m = (ps.m*ps.x^2)*[tt.m];      % dimensional melting rate (m^3/s)
% E = (ps.m*ps.x^2)*[tt.E];      % dimensional source terms  (m^3/s)
% Q_outQ = ps.Q*[tt.Q_outQ];     % dimensional channel outflux (m^3/s)
% Q_outq = ps.Q*[tt.Q_outq];     % dimensional sheet outflux (m^3/s)
% V_b = ps.V*[tt.Vb];            %
% R_b = ps.R*[tt.Rb];            %

% phi = (ps.phi/10^6)*[tt.phi];  % dimensional hydrulic potential (MPa)
% N = (ps.phi/10^6)*[tt.N];      % dimensional effective stress (MPa)
% hs = ps.x^2*ps.h*[tt.hs];      % integrated hs (m^3)
% he = ps.x^2*ps.h*[tt.he];      % integrated he (m^3)
% p_w = ps.phi/1e6*[tt.pwb];     % dimensional hydraulic potential at the lake (MPa)

% S = ps.x*ps.S*[tt.S];
% A = ps.x^2*sum(gg.Dx.*gg.Dy);

% t_s = t*pd.td;
% dVdt = diff(V_b)./diff(t_s); % dimensional volume rate

pp.meltE = @(t) 0.2*(10/1000/pd.td/ps.m)*(1-exp(-t/(30*pd.td/ps.t))); 

figure()
vt = linspace(0,300,101)*pd.td/ps.t;
plot(vt,pp.meltE(vt),'b-',LineWidth=1.5);
hold on 