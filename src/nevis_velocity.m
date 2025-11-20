function [u,v] = nevis_velocity(H,u,v,N,aa,pp,gg,oo)
% calculate depth-averaged velocities u,v using the shallow shelf
% approximation for given ice thickness H
% inputs
%   H       ice thickness on nodes [nIJ-by-1]
%   u       boundary values and initial guess for velocity [eIJ-by-1]
%   v       boundary values and initial guess for velocity [fIJ-by-1]
%   aa      prescribed fields and boundary terms [struct]
%   pp      parameters [struct]
%   gg      grid and discretisation operators [struct]
%   oo      options [struct]
% outputs 
%   u       x component of velocity on edges [eIJ-by-1]
%   v       y component of velocity on edges [fIJ-by-1]
%
% 26 July 2022
%   velocity_solve.m taken from iceflow/ice/ice_force_balance_new2.m and stripped down
% 10 Oct 2025 - edited to make nevis_velocity.m, compatibale with nevis

% ITERATION OPTIONS
if ~isfield(oo,'iter_max'), oo.iter_max = 20; end
if ~isfield(oo,'tol_vel'), oo.tol_vel = 1e-4; end
if ~isfield(oo,'tol_eta'), oo.tol_eta = 1e-6; end

% DIAGNOSTIC OPTIONS
if ~isfield(oo,'display_norms'), oo.display_norms = 0; end          % display norms at each iteration
if ~isfield(oo,'verb'), oo.verb = 0; end                            % verbose screen output
% [could add other diagnostic options here]  

% REQUIRED
if ~isfield(pp,'C'), pp.C = 1; end % dimensionless sliding coefficient (won't be used if aa contains field C)
if ~isfield(pp,'mu'), pp.mu = inf; end % dimensionless Coulomb friction coefficient (won't be used if aa contains field mu)
if ~isfield(pp,'A_glen'), pp.A_glen = 1; end % dimensionless Glen's-law coefficient
if ~isfield(pp,'n_glen'), pp.n_glen = pp.n_Glen; end % power-law in ice rheology
if ~isfield(pp,'n_slide'), pp.n_slide = pp.n_glen; end % power-law exponent in sliding law
if ~isfield(pp,'c60'), pp.c60 = 1; end % dimensionless coefficient in front of gravity term in force balance
if ~isfield(pp,'c61'), pp.c61 = 1; end % dimensionless coefficient in front of sliding law in force balance
if ~isfield(pp,'c62'), pp.c62 = 1; end % dimensionless coefficient in front of longitudinal stress in force balance
if ~isfield(pp,'eps_reg'), pp.eps_reg = 1e-16; end % regularisation on strain rates
if ~isfield(pp,'Ub_reg'), pp.Ub_reg = 1e-16; end % regularisation on sliding speed
if ~isfield(pp,'N_slide_reg'), pp.N_slide_reg = 1e-16; end % regularisation on effective pressure in sliding law
if ~isfield(pp,'taud_reg'), pp.taud_reg = 1e-16; end % regularisation on basal shear stress [ may not be needed ? ]
if ~isfield(pp,'C2'), pp.C2 = 0; end % added power-law coefficient in sliding law

% REQRUIED FIELDS
if ~isfield(aa,'C'), aa.C = pp.C*ones(gg.nIJ,1); end % dimensionless sliding coefficient 
if ~isfield(aa,'mu'), aa.mu = pp.mu*ones(gg.nIJ,1); end % dimensionless Coulomb friction coefficient

% unpack indices and matrices from gg
nIJ = gg.nIJ;
eIJ = gg.eIJ;
fIJ = gg.fIJ;
cIJ = gg.cIJ;
ns = gg.ns2;
es = gg.es2;
fs = gg.fs2;
cs = gg.cs2;
nin = gg.nin2;
ein = gg.ein2;
fin = gg.fin2;
cin = gg.cin2;
nbdyx = gg.nbdyx2;
nbdyy = gg.nbdyy2;
ebdy = gg.ebdy2;
fbdy = gg.fbdy2;
cbdy = gg.cbdy2;

nmeanx = gg.nmeanx;
nmeany = gg.nmeany;
nmeanc = gg.nmeanc;
emean = gg.emean;
fmean = gg.fmean;
cmean = gg.cmean;
nddx = gg.nddx;
nddy = gg.nddy;
eddx = gg.eddx;
eddy = gg.eddy;
fddx = gg.fddx;
fddy = gg.fddy;
cddx = gg.cddx;
cddy = gg.cddy;

% redfine mean operators to only include active nodes
temp = nmeanx(:,es)*ones(length(es),1); 
temp(temp==0) = inf;
nmeanx = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeanx;

temp = nmeany(:,fs)*ones(length(fs),1); 
temp(temp==0) = inf; 
nmeany = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeany;

temp = nmeanc(:,cs)*ones(length(cs),1); 
temp(temp==0) = inf; 
nmeanc = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeanc;

temp = emean(:,ns)*ones(length(ns),1); 
temp(temp==0) = inf; 
emean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*emean;

temp = fmean(:,ns)*ones(length(ns),1); 
temp(temp==0) = inf; 
fmean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*fmean;

temp = cmean(:,ns)*ones(length(ns),1); 
temp(temp==0) = inf; 
cmean = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*cmean;

gg.nmeanc = nmeanc; % ? for use in calculate_viscosity 

% required fields
b = aa.b; % bed elevation on nodes [nIJ-by-1]
C = aa.C; % basal friction coefficient on nodes [nIJ-by-1]
mu = aa.mu; % basal friction coefficient on nodes [nIJ-by-1]
Txx = 0*ones(nIJ,1); % membrane stress [only used on nbdyx] [nIJ-by-1]
Tyy = 0*ones(nIJ,1); % membrane stress [only used on nbdyy] [nIJ-by-1]
Txy = 0*ones(cIJ,1); % membrane stress [only used on cbdy] [cIJ-by-1]
s = b+H;      % surface elevation [nIJ-by-1]

% initial guess of velocities ignoring longitudinal stress and for infinite
% effective pressure
if isempty(u)
    taud_x = -pp.c60*(emean(:,ns)*H(ns)).*(eddx(:,:)*s(:));
    taud_y = -pp.c60*(fmean(:,ns)*H(ns)).*(fddy(:,:)*s(:));
    taud = ((gg.nmeanx*taud_x).^2+(gg.nmeany*taud_y).^2).^(1/2);
    u = (gg.emean*(aa.C.^(-1/pp.n_slide).*max(taud,pp.taud_reg).^(1./pp.n_slide-1))).*taud_x;
    v = (gg.fmean*(aa.C.^(-1/pp.n_slide).*max(taud,pp.taud_reg).^(1./pp.n_slide-1))).*taud_y;
end

% speed and viscosity
U = ( (nmeanx(:,es)*u(es)).^2+(nmeany(:,fs)*v(fs)).^2 ).^(1/2); % speed [n_IJ-by-1]
etabar = calculate_viscosity(u,v,pp,gg,oo); % viscosity [nIJ-by-1]

%iterate nonlinear viscosity
iter = 1;
while iter <= oo.iter_max 
    
    % remember current values for convergence check
    U_old = U; 
    etabar_old = etabar;
    
    % setup linear system for velocities
    H_times_etabar_matrix = sparse(nin,nin,H(nin).*etabar(nin),nIJ,nIJ);                   % [nIJ-by-nIJ]
    H_times_etabar_matrix_c = sparse(cin,cin,cmean(cin,ns)*(H(ns).*etabar(ns)),cIJ,cIJ);   % [cIJ-by-cIJ]
 
    Fx_res = -pp.c60*(emean(ein,ns)*H(ns)).*(eddx(ein,ns)*s(ns));
    if ~isempty(ebdy), Fx_res = Fx_res + pp.c62*( eddx(ein,nin)*(4*H_times_etabar_matrix(nin,nin)*nddx(nin,ebdy)) + eddy(ein,cin)*(H_times_etabar_matrix_c(cin,cin)*cddy(cin,ebdy)) )*u(ebdy); end % prescribed boundary velocity
    if ~isempty(fbdy), Fx_res = Fx_res + pp.c62*( eddx(ein,nin)*(2*H_times_etabar_matrix(nin,nin)*nddy(nin,fbdy)) + eddy(ein,cin)*(H_times_etabar_matrix_c(cin,cin)*cddx(cin,fbdy)) )*v(fbdy); end % prescribed boundary velocity
    if ~isempty(nbdyx), Fx_res = Fx_res + pp.c62*( eddx(ein,nbdyx)*( H(nbdyx).*Txx(nbdyx) ) ); end % prescribed boundary stress
    if ~isempty(cbdy), Fx_res = Fx_res + pp.c62*( eddy(ein,cbdy)*( (cmean(cbdy,ns)*H(ns)).*Txy(cbdy) ) ); end % prescribed boundary stress
    Fx_ub = pp.c61*sparse(1:length(ein),1:length(ein),-emean(ein,ns)*( slide_fun(U(ns),N(ns),C(ns),mu(ns),pp,gg,oo) ),length(ein),length(ein)) ...
        + pp.c62*( eddx(ein,nin)*(4*H_times_etabar_matrix(nin,nin)*nddx(nin,ein)) + eddy(ein,cin)*(H_times_etabar_matrix_c(cin,cin)*cddy(cin,ein)) );
    Fx_vb = pp.c62*( eddx(ein,nin)*(2*H_times_etabar_matrix(nin,nin)*nddy(nin,fin)) + eddy(ein,cin)*(H_times_etabar_matrix_c(cin,cin)*cddx(cin,fin)) );
    
    Fy_res = -pp.c60*(fmean(fin,ns)*H(ns)).*(fddy(fin,ns)*s(ns));
    if ~isempty(ebdy), Fy_res = Fy_res + pp.c62*( fddy(fin,nin)*(2*H_times_etabar_matrix(nin,nin)*nddx(nin,ebdy)) + fddx(fin,cin)*(H_times_etabar_matrix_c(cin,cin)*cddy(cin,ebdy)) )*u(ebdy); end % prescribed boundary velocity
    if ~isempty(fbdy), Fy_res = Fy_res + pp.c62*( fddy(fin,nin)*(4*H_times_etabar_matrix(nin,nin)*nddy(nin,fbdy)) + fddx(fin,cin)*(H_times_etabar_matrix_c(cin,cin)*cddx(cin,fbdy)) )*v(fbdy); end % prescribed boundary velocity
    if ~isempty(nbdyy), Fy_res = Fy_res + pp.c62*( fddy(fin,nbdyy)*( H(nbdyy).*Tyy(nbdyy) ) ); end % prescribed boundary stress
    if ~isempty(cbdy), Fy_res = Fy_res + pp.c62*( fddx(fin,cbdy)*( (cmean(cbdy,ns)*H(ns)).*Txy(cbdy) ) ); end % prescribed boundary stress
    Fy_ub = pp.c62*( fddy(fin,nin)*(2*H_times_etabar_matrix(nin,nin)*nddx(nin,ein)) + fddx(fin,cin)*(H_times_etabar_matrix_c(cin,cin)*cddy(cin,ein)) );
    Fy_vb = pp.c61*sparse(1:length(fin),1:length(fin),-fmean(fin,ns)*( slide_fun(U(ns),N(ns),C(ns),mu(ns),pp,gg,oo) ),length(fin),length(fin)) ...
        + pp.c62*( fddy(fin,nin)*(4*H_times_etabar_matrix(nin,nin)*nddy(nin,fin)) + fddx(fin,cin)*(H_times_etabar_matrix_c(cin,cin)*cddx(cin,fin)) );
    
    % solve for velocities 
    rhs = -[Fx_res; Fy_res];
    matrix = [Fx_ub Fx_vb; Fy_ub Fy_vb];
    temp = matrix\rhs;
    u(ein) = temp(1:length(ein)); 
    v(fin) = temp(length(ein)+(1:length(fin)));
    
    % update speed and viscosity
    U = ( (nmeanx(:,es)*u(es)).^2+(nmeany(:,fs)*v(fs)).^2 ).^(1/2); % [n_IJ-by-1]
    etabar = calculate_viscosity(u,v,pp,gg,oo); % [nIJ-by-1]
    
    % check convergence
%      norm_vel = norm((U-U_old)./U,2)/length(U);  % based on relative change of speed
   norm_vel = norm((U(ns)-U_old(ns)),inf)/(norm(U(ns),2)+pp.Ub_reg);  % based on relative change of speed
 %   norm_vel = norm(U-U_old,inf);    % based on absolute change of speed
   norm_eta = norm((etabar(ns)-etabar_old(ns)),2)/norm(etabar(ns),2);  % based on relative change of viscosity
  
    if oo.display_norms
        disp(['norm_vel = ',num2str(norm_vel),' after ',num2str(iter),' iterations']);
        disp(['norm_eta = ',num2str(norm_eta),' after ',num2str(iter),' iterations']);
    end

    % if norm_eta < oo.tol_eta && norm_vel < oo.tol_vel, 
    %     if oo.verb, disp(['Less than tolerance: stopping']); end; 
    %     break,
    % end
    % if norm_eta < oo.tol_eta, 
    %     if oo.verb, disp(['Less than tolerance: stopping']); end; 
    %     break, 
    % end
    if norm_vel < oo.tol_vel 
        if oo.verb, disp(['Less than tolerance: stopping']); end
        break,
    end
    if iter == oo.iter_max, disp(['!!! Max iterations reached: stopping']); break; end
    iter = iter+1; 

end

end

%% visocity
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



