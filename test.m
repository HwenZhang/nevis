figure('Name','Velocity','Position',[100 100 2000 400]);clf
subplot(1,3,1)
zC_hat = reshape(log10(C_hat_dim), gg.nI, gg.nJ);
zC_hat(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, zC_hat); shading flat; colorbar;
caxis([9 12]);
text(0.52, 0.98, '(a) Inverted log_{10}(C)', 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize',12, 'FontWeight','bold', 'Interpreter','tex');
xlabel('x'); ylabel('y'); axis equal tight;

subplot(1,3,2)
Un_obs = sqrt((gg.nmeanx2(:,gg.es2)*u_obs(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_obs(gg.fs2)).^2);
Un_inv = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);
Un_misfit = Un_inv;
Un_misfit(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, ps.u_b*pd.ty*reshape(Un_misfit,gg.nI,gg.nJ)); shading flat; colorbar;
text(0.52, 0.98, '(b) Speed (m/yr)', 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize',12, 'FontWeight','bold');
caxis([0 300]);  % symmetric colorbar
xlabel('x'); ylabel('y'); axis equal tight;

subplot(1,3,3)
Un_obs = sqrt((gg.nmeanx2(:,gg.es2)*u_obs(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_obs(gg.fs2)).^2);
Un_inv = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);
Un_misfit = Un_inv - Un_obs;
Un_misfit(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, ps.u_b*pd.ty*reshape(Un_misfit,gg.nI,gg.nJ)); shading flat; colorbar;
text(0.52, 0.98, '(c) Speed misfit (m/yr)', 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize',12, 'FontWeight','bold');
caxis([-50 50]);  % symmetric colorbar
xlabel('x'); ylabel('y'); axis equal tight;

exportgraphics(gcf, 'vel_inv.png')

%% plot the regularisation and damping terms
c_hat = log(C_hat);
Lc = L * c_hat;                              % gradient of c = log(C)
reg_node = (L' * L) * c_hat;                 % L^T L c  on nodes [nIJ x 1]
damp_node = c_hat - c_prior;                 % c - c_prior on nodes [nIJ x 1]
z_reg  = reshape(reg_node,  gg.nI, gg.nJ);  z_reg(gg.nout)  = NaN;
z_damp = reshape(damp_node, gg.nI, gg.nJ);  z_damp(gg.nout) = NaN;
% combine: alpha * |Lc|^2  +  gamma * |c-c0|^2  per node (for visualisation)
reg_energy  = opts_inv.alpha * reg_node .* c_hat;   % ~ alpha * c^T L^T L c contribution per node
damp_energy = opts_inv.gamma * damp_node.^2;        % gamma * (c-c0)^2 per node
total_penalty = reg_energy + damp_energy;
total_penalty(gg.nout) = NaN;
z = reshape(total_penalty, gg.nI, gg.nJ);

U_nodes_fric = sqrt((gg.nmeanx2(:,gg.es2)*u_inv(gg.es2)).^2 + ...
                    (gg.nmeany2(:,gg.fs2)*v_inv(gg.fs2)).^2);  % [nIJ x 1]
Ub_fric  = max(U_nodes_fric, pp.Ub_reg);      % regularised speed
Nr_fric  = max(N_current,    pp.N_slide_reg);  % regularised N
C_fric   = C_hat;                             % inverted C
mu_fric  = pp.mu * ones(gg.nIJ, 1);           % Coulomb coefficient

% tau_b / Ub  (dimensionless)
taub_over_Ub = Nr_fric .* Ub_fric.^(1/pp.n_slide-1) .* ...
    (mu_fric.^(-pp.n_slide) .* Ub_fric + C_fric.^(-pp.n_slide) .* Nr_fric.^pp.n_slide).^(-1/pp.n_slide);

% tau_b  [dimensional: Pa]
tau_b_dim = taub_over_Ub .* U_nodes_fric * ps.tau;   % [Pa]
tau_b_kPa = tau_b_dim / 1e3;                                   % [kPa]

% Driving stress [kPa] (for comparison)
dsdx_fric = gg.nmeanx2(:,gg.es2) * (gg.eddx(gg.es2,:) * aa.s);
dsdy_fric = gg.nmeany2(:,gg.fs2) * (gg.fddy(gg.fs2,:) * aa.s);
tau_d_kPa = pp.c60 * aa.H .* sqrt(dsdx_fric.^2 + dsdy_fric.^2) * pd.rho_i*pd.g*ps.z^2/ps.x/ 1e3;  % [kPa]

figure('Name','Basal Friction Distribution','Position',[100 100 2000 400]);

% (a) Basal friction tau_b [kPa]
subplot(1,3,1);
z = reshape(tau_b_kPa, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
clim([0 300]);
text(0.52, 0.98, '(a) Basal friction \tau_b [kPa]', 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize',12, 'FontWeight','bold', 'Interpreter','tex');
xlabel('x'); ylabel('y'); 
axis equal tight;

% (b) Driving stress tau_d [kPa]
subplot(1,3,2);
z = reshape(tau_d_kPa, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
clim([0 300]);
text(0.52, 0.98, '(b) Driving stress \tau_d [kPa]', 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize',12, 'FontWeight','bold', 'Interpreter','tex');
xlabel('x'); ylabel('y'); 
axis equal tight;

% (c) Friction fraction tau_b / tau_d
subplot(1,3,3);
fric_frac = tau_b_kPa ./ max(tau_d_kPa, 0.1);  % avoid div by zero
fric_frac(gg.nout) = NaN;
z = reshape(fric_frac, gg.nI, gg.nJ); z(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
clim([0 2]);
text(0.52, 0.98, '(c) Friction fraction \tau_b / \tau_d', 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize',12, 'FontWeight','bold', 'Interpreter','tex');
xlabel('x'); ylabel('y'); axis equal tight;

saveas(gcf, 'tau_inv.png');

%% ============================================================
%  Principal deviatoric stresses and directions
%  ============================================================
% Strain rates on nodes / corners
epsxx = gg.nddx(:,gg.es2)*u_inv(gg.es2);           % du/dx at nodes [nIJ x 1]
epsyy = gg.nddy(:,gg.fs2)*v_inv(gg.fs2);           % dv/dy at nodes [nIJ x 1]
epsxy_c = 0.5*(gg.cddy(:,gg.es2)*u_inv(gg.es2) ... % (du/dy+dv/dx)/2 at corners
             + gg.cddx(:,gg.fs2)*v_inv(gg.fs2));
epsxy = gg.nmeanc(:,gg.cs2)*epsxy_c(gg.cs2);       % averaged to nodes [nIJ x 1]

% Viscosity at nodes (dimensionless)
I2 = pp.eps_reg^2 + epsxx.^2 + epsyy.^2 + epsxx.*epsyy + epsxy.^2;
etabar = pp.A_glen^(-1/pp.n_glen) .* I2.^((1/pp.n_glen-1)/2);

% Deviatoric stress tensor  sigma'_ij = 2 * eta * eps_ij  (depth-integrated: H*sigma')
% Dimensional scale: ps.tau [Pa] (stress scale)
stress_scale = ps.tau / 1e3;  % convert to kPa
txx = 2 * etabar .* epsxx * stress_scale;   % sigma'_xx [kPa]
tyy = 2 * etabar .* epsyy * stress_scale;   % sigma'_yy [kPa]
txy = 2 * etabar .* epsxy * stress_scale;   % sigma'_xy [kPa]

% Principal stresses:  lambda = (txx+tyy)/2 ± sqrt( ((txx-tyy)/2)^2 + txy^2 )
tmean = 0.5*(txx + tyy);
R     = sqrt( (0.5*(txx - tyy)).^2 + txy.^2 );
sigma1 = tmean + R;   % max principal stress [kPa]
sigma2 = tmean - R;   % min principal stress [kPa]

% Principal direction (angle of max principal stress w.r.t. x-axis)
theta = 0.5 * atan2(2*txy, txx - tyy);  % [rad]

% Direction vectors for quiver
dx1 = cos(theta);   dy1 = sin(theta);    % max principal
dx2 = cos(theta+pi/2); dy2 = sin(theta+pi/2); % min principal (perpendicular)

% Mask outside nodes
sigma1(gg.nout) = NaN; sigma2(gg.nout) = NaN;
dx1(gg.nout) = NaN; dy1(gg.nout) = NaN;
dx2(gg.nout) = NaN; dy2(gg.nout) = NaN;

% Reshape to grid
nI = gg.nI; nJ = gg.nJ;
S1 = reshape(sigma1, nI, nJ);
S2 = reshape(sigma2, nI, nJ);

figure('Name','Principal Deviatoric Stresses','Position',[100 100 1600 500]);

% (a) Max principal stress magnitude
subplot(1,3,1);
pcolor(gg.nx, gg.ny, S1); shading flat; colorbar;
title('\sigma''_1 (max principal) [kPa]'); xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, jet);

% (b) Min principal stress magnitude
subplot(1,3,2);
pcolor(gg.nx, gg.ny, S2); shading flat; colorbar;
title('\sigma''_2 (min principal) [kPa]'); xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, jet);

% (c) Principal directions overlaid on effective stress
subplot(1,3,3);
tau_eff = sqrt(0.5*(txx.^2 + tyy.^2) + txy.^2);
tau_eff(gg.nout) = NaN;
pcolor(gg.nx, gg.ny, reshape(tau_eff, nI, nJ)); shading flat; colorbar;
hold on;
% subsample for readability
skip = 4;
ii = 1:skip:nI; jj = 1:skip:nJ;
xq = gg.nx(ii,jj); yq = gg.ny(ii,jj);
% scale arrows by principal magnitude for visibility
s1_q = reshape(sigma1, nI, nJ); s1_q = s1_q(ii,jj);
s2_q = reshape(sigma2, nI, nJ); s2_q = s2_q(ii,jj);
dx1_q = reshape(dx1, nI, nJ); dx1_q = dx1_q(ii,jj);
dy1_q = reshape(dy1, nI, nJ); dy1_q = dy1_q(ii,jj);
dx2_q = reshape(dx2, nI, nJ); dx2_q = dx2_q(ii,jj);
dy2_q = reshape(dy2, nI, nJ); dy2_q = dy2_q(ii,jj);
% Plot as line segments (±direction) — tensile in red, compressive in blue
quiver(xq, yq,  s1_q.*dx1_q,  s1_q.*dy1_q, 0.5, 'r', 'LineWidth', 1);
quiver(xq, yq, -s1_q.*dx1_q, -s1_q.*dy1_q, 0.5, 'r', 'LineWidth', 1, 'ShowArrowHead', 'off');
quiver(xq, yq,  s2_q.*dx2_q,  s2_q.*dy2_q, 0.5, 'b', 'LineWidth', 1);
quiver(xq, yq, -s2_q.*dx2_q, -s2_q.*dy2_q, 0.5, 'b', 'LineWidth', 1, 'ShowArrowHead', 'off');
hold off;
title('\tau_{eff} [kPa] + principal dirs (r=\sigma_1, b=\sigma_2)');
xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, jet);
drawnow;

%% ============================================================
%  Direct estimate of tau_b and C  (no Coulomb cap, n_slide = 1)
%  ============================================================
%  Use SSA force balance with OBSERVED velocities:
%    c60*H*grad(s) + c62*div(H*sigma') = c61*tau_b
%  where sigma' (membrane stress) is computed from observed strain rates
%  and Glen's law viscosity.
%  With no Coulomb cap and n_slide=1:  tau_b = C * Ub, so C = tau_b / Ub.

% --- Clean observed velocities (NaN → 0) ---
u_obs_d = u_obs; u_obs_d(isnan(u_obs_d)) = 0;
v_obs_d = v_obs; v_obs_d(isnan(v_obs_d)) = 0;

% --- Shorthand ---
ns_d = gg.ns2; es_d = gg.es2; fs_d = gg.fs2; cs_d = gg.cs2;
ein_d = gg.ein2; fin_d = gg.fin2; cin_d = gg.cin2;

% --- Renormalize mean operators (same convention as nevis_velocity) ---
temp = gg.emean(:,ns_d)*ones(length(ns_d),1); temp(temp==0)=inf;
emean_d = sparse(1:gg.eIJ,1:gg.eIJ,temp.^(-1),gg.eIJ,gg.eIJ)*gg.emean;
temp = gg.fmean(:,ns_d)*ones(length(ns_d),1); temp(temp==0)=inf;
fmean_d = sparse(1:gg.fIJ,1:gg.fIJ,temp.^(-1),gg.fIJ,gg.fIJ)*gg.fmean;
temp = gg.cmean(:,ns_d)*ones(length(ns_d),1); temp(temp==0)=inf;
cmean_d = sparse(1:gg.cIJ,1:gg.cIJ,temp.^(-1),gg.cIJ,gg.cIJ)*gg.cmean;
temp = gg.nmeanc(:,cs_d)*ones(length(cs_d),1); temp(temp==0)=inf;
nmeanc_d = sparse(1:gg.nIJ,1:gg.nIJ,temp.^(-1),gg.nIJ,gg.nIJ)*gg.nmeanc;

% --- Strain rates from observed velocity ---
epsxx_d = gg.nddx(:,es_d)*u_obs_d(es_d);       % du/dx at nodes
epsyy_d = gg.nddy(:,fs_d)*v_obs_d(fs_d);       % dv/dy at nodes
epsxy_c_d = 0.5*(gg.cddy(:,es_d)*u_obs_d(es_d) ...
               + gg.cddx(:,fs_d)*v_obs_d(fs_d)); % eps_xy at corners
epsxy_n_d = nmeanc_d(:,cs_d)*epsxy_c_d(cs_d);   % averaged to nodes

% --- Viscosity from Glen's law ---
I2_d = pp.eps_reg^2 + epsxx_d.^2 + epsyy_d.^2 + epsxx_d.*epsyy_d + epsxy_n_d.^2;
etabar_d = pp.A_glen^(-1/pp.n_glen) .* I2_d.^((1/pp.n_glen-1)/2);

% --- Membrane stress at nodes (T_xx, T_yy) and corners (T_xy) ---
% T_xx = H*eta*(4*eps_xx + 2*eps_yy),  T_yy = H*eta*(2*eps_xx + 4*eps_yy)
T_xx_d = aa.H .* etabar_d .* (4*epsxx_d + 2*epsyy_d);
T_yy_d = aa.H .* etabar_d .* (2*epsxx_d + 4*epsyy_d);
T_xx_d(gg.nout) = 0;  T_yy_d(gg.nout) = 0;

% T_xy = H_c * eta_c * (du/dy + dv/dx)  at corners
Heta_c_d = cmean_d(cin_d,ns_d) * (aa.H(ns_d) .* etabar_d(ns_d));
T_xy_c_d = zeros(gg.cIJ, 1);
T_xy_c_d(cin_d) = Heta_c_d .* (gg.cddy(cin_d,es_d)*u_obs_d(es_d) ...
                              + gg.cddx(cin_d,fs_d)*v_obs_d(fs_d));

% --- Membrane stress divergence at interior edges ---
% x-edges: d(T_xx)/dx + d(T_xy)/dy
mem_x = pp.c62 * (gg.eddx(ein_d,:)*T_xx_d + gg.eddy(ein_d,:)*T_xy_c_d);
% y-edges: d(T_xy)/dx + d(T_yy)/dy
mem_y = pp.c62 * (gg.fddx(fin_d,:)*T_xy_c_d + gg.fddy(fin_d,:)*T_yy_d);

% --- Driving stress at interior edges ---
driv_x = -pp.c60 * (emean_d(ein_d,ns_d)*aa.H(ns_d)) .* (gg.eddx(ein_d,ns_d)*aa.s(ns_d));
driv_y = -pp.c60 * (fmean_d(fin_d,ns_d)*aa.H(ns_d)) .* (gg.fddy(fin_d,ns_d)*aa.s(ns_d));

% --- Friction traction at interior edges (from SSA force balance) ---
% SSA:  -c61*(tau_b/Ub)*u + c62*membrane + driving = 0
%  =>   (tau_b/Ub)*u = (membrane + driving) / c61
fric_x_total = (mem_x + driv_x) / pp.c61;   % full balance [len(ein) x 1]
fric_y_total = (mem_y + driv_y) / pp.c61;

fric_x_driv  = driv_x / pp.c61;             % driving-only estimate
fric_y_driv  = driv_y / pp.c61;

fric_x_mem   = mem_x / pp.c61;              % membrane contribution
fric_y_mem   = mem_y / pp.c61;

% --- Average friction traction to nodes ---
% Full (with membrane)
fx_full = zeros(gg.eIJ,1); fx_full(ein_d) = fric_x_total;
fy_full = zeros(gg.fIJ,1); fy_full(fin_d) = fric_y_total;
tau_bx_fb = gg.nmeanx2(:,es_d)*fx_full(es_d);
tau_by_fb = gg.nmeany2(:,fs_d)*fy_full(fs_d);
tau_b_fb  = sqrt(tau_bx_fb.^2 + tau_by_fb.^2);  % |tau_b| from force balance

% Driving only (no membrane)
fx_dr = zeros(gg.eIJ,1); fx_dr(ein_d) = fric_x_driv;
fy_dr = zeros(gg.fIJ,1); fy_dr(fin_d) = fric_y_driv;
tau_bx_dr = gg.nmeanx2(:,es_d)*fx_dr(es_d);
tau_by_dr = gg.nmeany2(:,fs_d)*fy_dr(fs_d);
tau_b_dr  = sqrt(tau_bx_dr.^2 + tau_by_dr.^2);   % |tau_d|/c61

% Membrane only
fx_mm = zeros(gg.eIJ,1); fx_mm(ein_d) = fric_x_mem;
fy_mm = zeros(gg.fIJ,1); fy_mm(fin_d) = fric_y_mem;
tau_mx = gg.nmeanx2(:,es_d)*fx_mm(es_d);
tau_my = gg.nmeany2(:,fs_d)*fy_mm(fs_d);

% --- Observed speed on nodes ---
U_obs_direct = sqrt((gg.nmeanx2(:,es_d)*u_obs_d(es_d)).^2 + ...
                     (gg.nmeany2(:,fs_d)*v_obs_d(fs_d)).^2);

% --- Dimensionalize ---
tau_b_fb_kPa  = tau_b_fb * ps.tau / 1e3;        % full balance [kPa]
tau_b_dr_kPa  = tau_b_dr * ps.tau / 1e3;        % driving only [kPa]
mem_effect_kPa = tau_b_fb_kPa - tau_b_dr_kPa;   % membrane adds/removes friction

% --- C estimates (nondimensional, then dimensional) ---
C_fb_nondim = tau_b_fb ./ max(U_obs_direct, 1e-10);   % with membrane
C_fb_dim    = C_fb_nondim * (ps.tau / ps.u_b);
C_dr_nondim = tau_b_dr ./ max(U_obs_direct, 1e-10);   % driving only
C_dr_dim    = C_dr_nondim * (ps.tau / ps.u_b);

% Mask outside
tau_b_fb_kPa(gg.nout) = NaN;  tau_b_dr_kPa(gg.nout) = NaN;
mem_effect_kPa(gg.nout) = NaN;
C_fb_dim(gg.nout) = NaN;  C_dr_dim(gg.nout) = NaN;

% ===== Figure =====
figure('Name','Direct tau_b and C (no Coulomb, n_{slide}=1)','Position',[100 100 1800 700]);

% (a) tau_b from full force balance [kPa]
subplot(2,3,1);
z = reshape(tau_b_fb_kPa, gg.nI, gg.nJ);
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
clim([0 200]);
title('\tau_b^{fb} (driving + membrane) [kPa]'); xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, jet);

% (b) tau_d only [kPa]
subplot(2,3,2);
z = reshape(tau_b_dr_kPa, gg.nI, gg.nJ);
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
clim([0 200]);
title('\tau_d (driving only) [kPa]'); xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, jet);

% (c) Membrane contribution [kPa] (positive = membrane requires extra friction)
subplot(2,3,3);
z = reshape(mem_effect_kPa, gg.nI, gg.nJ);
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
clim([-300 300]);
title('\tau_b^{fb} - \tau_d (membrane effect) [kPa]'); xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, interp1([0 0.5 1],[0 0 1; 1 1 1; 1 0 0],linspace(0,1,256)));

% (d) log10(C_direct) with membrane
subplot(2,3,4);
z = reshape(log10(C_fb_dim), gg.nI, gg.nJ);
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
title('log_{10}(C^{fb}) (with membrane)'); xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, jet);

% (e) log10(C_inv / C_direct)
subplot(2,3,5);
C_inv_dim_node = C_hat_dim; C_inv_dim_node(gg.nout) = NaN;
z = reshape(log10(C_inv_dim_node) - log10(C_fb_dim), gg.nI, gg.nJ);
pcolor(gg.nx, gg.ny, z); shading flat; colorbar;
clim([-2 2]);
title('log_{10}(C^{inv} / C^{fb})'); xlabel('x'); ylabel('y'); axis equal tight;
colormap(gca, interp1([0 0.5 1],[0 0 1; 1 1 1; 1 0 0],linspace(0,1,256)));

% (f) Scatter: C_direct vs C_inverted
subplot(2,3,6);
nin = gg.nin2;
C_fb_plot  = log10(C_fb_dim(nin));
C_inv_plot = log10(C_hat_dim(nin));
U_obs_dim_plot = U_obs_direct(nin) * ps.u_b * pd.ty;  % [m/yr]
scatter(C_fb_plot, C_inv_plot, 8, U_obs_dim_plot, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
lims = [min([C_fb_plot; C_inv_plot]), max([C_fb_plot; C_inv_plot])];
plot(lims, lims, 'k--', 'LineWidth', 1.5);
hold off;
colorbar; colormap(gca, jet);
xlabel('log_{10}(C^{fb})'); ylabel('log_{10}(C^{inv})');
title('C^{fb} vs C^{inv} (color = speed m/yr)');
axis equal; grid on;
drawnow;

fprintf('\n=== Direct Estimate (no Coulomb, n_slide=1, with membrane) ===\n');
fprintf('  tau_b (force balance) [kPa]: mean=%.1f, median=%.1f, max=%.1f\n', ...
    nanmean(tau_b_fb_kPa(nin)), nanmedian(tau_b_fb_kPa(nin)), nanmax(tau_b_fb_kPa(nin)));
fprintf('  tau_d (driving only)  [kPa]: mean=%.1f, median=%.1f, max=%.1f\n', ...
    nanmean(tau_b_dr_kPa(nin)), nanmedian(tau_b_dr_kPa(nin)), nanmax(tau_b_dr_kPa(nin)));
fprintf('  membrane effect       [kPa]: mean=%.1f, median=%.1f, range=[%.1f, %.1f]\n', ...
    nanmean(mem_effect_kPa(nin)), nanmedian(mem_effect_kPa(nin)), ...
    nanmin(mem_effect_kPa(nin)), nanmax(mem_effect_kPa(nin)));
fprintf('  log10(C^{fb}_dim):  mean=%.2f, median=%.2f, range=[%.2f, %.2f]\n', ...
    nanmean(log10(C_fb_dim(nin))), nanmedian(log10(C_fb_dim(nin))), ...
    nanmin(log10(C_fb_dim(nin))), nanmax(log10(C_fb_dim(nin))));
fprintf('  log10(C^{inv}_dim): mean=%.2f, median=%.2f, range=[%.2f, %.2f]\n', ...
    nanmean(log10(C_hat_dim(nin))), nanmedian(log10(C_hat_dim(nin))), ...
    nanmin(log10(C_hat_dim(nin))), nanmax(log10(C_hat_dim(nin))));
fprintf('  Corr(log C^{fb}, log C^{inv}) = %.3f\n', ...
    corr(log10(C_fb_dim(nin)), log10(C_hat_dim(nin)), 'rows', 'complete'));
fprintf('==============================================================\n');