function [tauxx,tauyy,tauxy,Txx,Tyy,Txy,tau_b,taudx,taudy,sigma1,sigma2,t1,t2] = nevis_stresses(H,u,v,N,aa,pp,gg,oo)

    % calculate membrane stresses using the shallow shelf approximation for 
    % given ice thickness H, velocity (u,v), and effective pressure N
    % also returns basal shear stress and gravitational driving stresses
    %
    % 15 October 2025 edited from calculate_stresses.m to go together with nevis_velocity_solve.m

    % % fill in any missing required fields and parameters
    % if ~isfield(pp,'C'), pp.C = 1; end % dimensionless sliding coefficient (won't be used if aa contains field C)
    % if ~isfield(pp,'mu'), pp.mu = inf; end % dimensionless Coulomb friction coefficient (won't be used if aa contains field mu)
    % if ~isfield(pp,'eps_reg'), pp.eps_reg = 1e-16; end % regularisation on strain rates
    % if ~isfield(pp,'Ub_reg'), pp.Ub_reg = 1e-16; end % regularisation on sliding speed
    % if ~isfield(pp,'N_slide_reg'), pp.N_slide_reg = 1e-16; end % regularisation on effective pressure in sliding law
    % if ~isfield(pp,'C2'), pp.C2 = 0; end % added power-law coefficient in sliding law
    % if ~isfield(aa,'C'), aa.C = pp.C*ones(gg.nIJ,1); end % dimensionless sliding coefficient 
    % if ~isfield(aa,'mu'), aa.mu = pp.mu*ones(gg.nIJ,1); end % dimensionless Coulomb friction coefficient

    % unpack indices and matrices from gg
    nmeanx = gg.nmeanx2;
    nmeany = gg.nmeany2;
    nmeanc = gg.nmeanc2;
    emean = gg.emean2;
    fmean = gg.fmean2;
    cmean = gg.cmean2;
    ns = gg.ns2;
    es = gg.es2;
    fs = gg.fs2;
    cs = gg.cs2;
    nin = gg.nin2;
    cin = gg.cin2;

    % redfine mean operators to only include active nodes
    temp = nmeanx(:,es)*ones(length(es),1); temp(temp==0) = inf; nmeanx = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeanx;
    temp = nmeany(:,fs)*ones(length(fs),1); temp(temp==0) = inf; nmeany = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeany;
    temp = nmeanc(:,cs)*ones(length(cs),1); temp(temp==0) = inf; nmeanc = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeanc;
    temp = emean(:,ns)*ones(length(ns),1); temp(temp==0) = inf; emean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*emean;
    temp = fmean(:,ns)*ones(length(ns),1); temp(temp==0) = inf; fmean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*fmean;
    temp = cmean(:,ns)*ones(length(ns),1); temp(temp==0) = inf; cmean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*cmean;
    
    % viscosity
    etabar = calculate_viscosity(u,v,pp,gg,oo);

    % membrane stresses, scaled with ps.eta*ps.u/ps.x
    tauxx = 0*ones(gg.nIJ,1); % only used on nbdyx
    tauyy = 0*ones(gg.nIJ,1); % only used on nbdyy
    tauxy = 0*ones(gg.cIJ,1); % only used on cbdy
    tauxx(nin,:) = 2*etabar(nin).*(gg.nddx(nin,:)*u);
    tauyy(nin,:) = 2*etabar(nin).*(gg.nddy(nin,:)*v);
    tauxy(cin,:) = (cmean(cin,ns)*etabar(ns)).*(gg.cddy(cin,:)*u+gg.cddx(cin,:)*v);
    tauxy_nin = nmeanc(:,cin)*tauxy(cin); % tauxy at nin, interpolated from cin

    % depth-integrated membrane stresses, scaled with ps.eta*ps.u/ps.x [= epsilon*pd.rho_i*pd.g*ps.z]
    Txx = H.*(2*tauxx+tauyy);
    Tyy = H.*(tauxx+2*tauyy);
    Txy = (cmean(:,ns)*H(ns)).*tauxy;

    % principal stresses and directions, scaled with ps.eta*ps.u/ps.x
    sigma1 = 0.5*(tauxx+tauyy) + sqrt(0.25*(tauxx-tauyy).^2+tauxy_nin.^2);
    sigma2 = 0.5*(tauxx+tauyy) - sqrt(0.25*(tauxx-tauyy).^2+tauxy_nin.^2);
    t1 = atan(2*tauxy_nin./(tauxx-tauyy))/2; % angle of sigma1 relative to x-axis, in radians   
    t2 = t1 + pi/2; % angle of sigma2 relative to x-axis, in radians

    % sliding speed, scaled with ps.u
    U = ( (nmeanx(:,es)*u(es)).^2+(nmeany(:,fs)*v(fs)).^2 ).^(1/2);

    % basal shear stress, scaled with pd.rho_i*pd.g*ps.z^2/ps.x
    tau_b = slide_fun(U,N,aa.C,pp.mu,pp,gg,oo).*U;

    % driving stress, scaled with pd.rho_i*pd.g*ps.z^2/ps.x
    % s = aa.b+H;
    s = aa.s;
    taudx = 0*ones(gg.eIJ,1); % only used outside domain
    taudy = 0*ones(gg.fIJ,1); % only used outside domain
    taudx(es) = -(emean(es,ns)*H(ns)).*(gg.eddx(es,ns)*s(ns));
    taudy(fs) = -(fmean(fs,ns)*H(ns)).*(gg.fddy(fs,ns)*s(ns));

end

%% viscosity
function eta = calculate_viscosity(u,v,pp,gg,oo)
% calculate viscosity for given depth-averaged velocities u and v
% [ the dimensional factor at the front of eta should be 1/2A^(1/n), after
% scaling with 1/2[A]^(1/n) it becomes A^(-1/n) in terms of dimensionless A
    epsxx = gg.nddx(:,gg.es2)*u(gg.es2);
    epsyy = gg.nddy(:,gg.fs2)*v(gg.fs2);
    epsxy = 1/2*(gg.cddy(:,gg.es2)*u(gg.es2)+gg.cddx(:,gg.fs2)*v(gg.fs2));
    eta = pp.A_glen.^(-1/pp.n_glen).*( pp.eps_reg^2+epsxx.^2+epsyy.^2+epsxx.*epsyy+(gg.nmeanc(:,gg.cs2)*epsxy(gg.cs2)).^2 ).^((1/pp.n_glen-1)/2);
end

%% sliding law
function tau_b_over_Ub = slide_fun(Ub,N,C,mu,pp,gg,oo)
% cavity-based sliding law
% tau_b ~ mu*N for large Ub, tau_b ~ C*Ub^(1/n) for small Ub 
    tau_b_over_Ub = max(N,pp.N_slide_reg).*max(Ub,pp.Ub_reg).^(1/pp.n_slide-1).*(mu.^(-pp.n_slide).*max(Ub,pp.Ub_reg)+C.^(-pp.n_slide).*max(N,pp.N_slide_reg).^pp.n_slide).^(-1/pp.n_slide) + pp.C2*max(Ub,pp.Ub_reg).^(1/pp.n_slide-1);
end