% Export ice velocity and stress fields for specific time slices.
% Called from Python via subprocess. Reads 'export_ice_fields.mat' 
% containing: filenames (cell), casename (string).
% Appends to each .mat file: ice_fields struct with ux, vy, U, 
% sigma1, sigma2, t1, tauxx, tauyy, tauxy, tau_b, stress_scale.

script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
cd(parent_dir);

clear variables -except parent_dir script_dir;
clc;

addpath(genpath('./src'));
python_data = load('export_ice_fields.mat');
casename = python_data.casename;
load(['./results/' casename '/' casename]);

for i = 1:length(python_data.filenames)
    filename = python_data.filenames{i};
    vva = load(filename, 'vv');
    vv = vva.vv;
    aa = nevis_inputs(vv.t, aa, vv, pp, gg, oo);
    pp.deltap_reg = 5.0e7/ps.phi;
    vv2 = nevis_backbone(inf, vv, vv, aa, pp, gg, oo);
    vv2 = nevis_nodedischarge(vv2, aa, pp, gg, oo);
    qnet = ps.qs * (vv2.qs + vv2.qe + vv2.qQ + 0 * vv2.Q);

    % velocity on nodes
    uxn = gg.nmeanx2(:,gg.es2) * vv.u(gg.es2);
    vyn = gg.nmeany2(:,gg.fs2) * vv.v(gg.fs2);
    Un  = sqrt(uxn.^2 + vyn.^2);

    % stresses
    [tauxx,tauyy,tauxy,~,~,~,tau_b,~,~,sigma1,sigma2,t1,~] = ...
        nevis_stresses(aa.H, vv.u, vv.v, aa.phi_0 - vv.phi, aa, pp, gg, oo);

    % elastic bending stress
    if isfield(pd, 'nu')
        poisson_ice = pd.nu;
    else
        poisson_ice = 0.33;
    end
    stress_membrane = ps.eta * ps.u / ps.x;
    tauxy_n = gg.nmeanc(:, gg.cs2) * tauxy(gg.cs2);

    hb_xx = gg.nddx * (gg.eddx * vv.hb);
    hb_yy = gg.nddy * (gg.fddy * vv.hb);
    hb_xy = 0.5 * (gg.cddy * (gg.eddx * vv.hb) + gg.cddx * (gg.fddy * vv.hb));
    hb_xy_n = gg.nmeanc(:, gg.cs2) * hb_xy(gg.cs2);

    curvature_scale = ps.hb / ps.x^2;
    kappa_xx = curvature_scale * hb_xx;
    kappa_yy = curvature_scale * hb_yy;
    kappa_xy = curvature_scale * hb_xy_n;

    H_dim = ps.z * max(aa.H, 0);
    bending_prefactor = pd.Ye * H_dim ./ (2 * (1 - poisson_ice^2));
    shear_prefactor = pd.Ye * H_dim ./ (2 * (1 + poisson_ice));

    tauxx_elastic = -bending_prefactor .* (kappa_xx + poisson_ice * kappa_yy);
    tauyy_elastic = -bending_prefactor .* (kappa_yy + poisson_ice * kappa_xx);
    tauxy_elastic = -shear_prefactor .* kappa_xy;

    [sigma1_elastic, ~, t1_elastic] = ...
        nevis_principal_stress(tauxx_elastic, tauyy_elastic, tauxy_elastic);

    tauxx_total = stress_membrane * tauxx + tauxx_elastic;
    tauyy_total = stress_membrane * tauyy + tauyy_elastic;
    tauxy_total = stress_membrane * tauxy_n + tauxy_elastic;

    [sigma1_total, ~, t1_total] = ...
        nevis_principal_stress(tauxx_total, tauyy_total, tauxy_total);

    % pack into struct
    ice_fields = struct();
    ice_fields.uxn = uxn;
    ice_fields.vyn = vyn;
    ice_fields.Un = Un;
    ice_fields.tauxx = tauxx;
    ice_fields.tauyy = tauyy;
    ice_fields.tauxy = tauxy;
    ice_fields.sigma1 = sigma1;
    ice_fields.sigma2 = sigma2;
    ice_fields.t1 = t1;
    ice_fields.tau_b = tau_b;
    ice_fields.stress_scale = ps.eta * ps.u / ps.x / 1e3;  % to kPa
    ice_fields.vel_scale = ps.u_b * pd.ty;                  % to m/yr
    ice_fields.nout = gg.nout;
    ice_fields.nI = gg.nI;
    ice_fields.nJ = gg.nJ;
    ice_fields.qnet = qnet;
    ice_fields.t_days = vv.t * ps.t / pd.td;
    ice_fields.sigma1_elastic_pa = sigma1_elastic;
    ice_fields.sigma1_total_pa = sigma1_total;
    ice_fields.t1_elastic = t1_elastic;
    ice_fields.t1_total = t1_total;

    save(filename, 'ice_fields', '-append');
    fprintf('Exported ice fields for %s (t=%.1f d)\n', filename, ice_fields.t_days);
end
