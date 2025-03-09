function [vv2,F,F1,F2,F3,F4,F5,F6,F7,F8,J] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo) 
    %% Inputs
    %   vv      struct containing solution variables
    %   vv0     struct containing current solution variables
    %   aa      prescribed variables and inputs 
    %   pp      parameters
    %   gg      grid and operators
    %   oo      options   
    % Outputs
    %   vv2     struct containing expanded solution variables
    %   F       Residual vector
    %   J       Jacobian matrix
    %
    % IJH 13 August 2014 : largely taken from hydro_timestep_diag
    %   15 Jan 2015 : add aa.lcx etc to channel opening terms; to use as 
    %                 indicator for allowing individual channels or not
    %   26 Oct 2015 : add he to output vv2
    %   18 Set 2016 : edit DXic in jacobian for pressure dependence (c41), edit
    %                 jacobian for sliding opening term in channels (1>c36*Sx)

    %% OPTIONS
    if nargin<7, oo = struct; end
    % don't include sheet equation in residual
    if ~isfield(oo,'no_sheet'), oo.no_sheet = 0; end        
    % don't include channel equations in residual
    if ~isfield(oo,'no_channels'), oo.no_channels = 0; end       
    % use maximum of Xi on nodes at either end rather than mean
    if ~isfield(oo,'maxXi'), oo.maxXi = 0; end                 
    % evaluate Xi rather than taking Xi from aa
    if ~isfield(oo,'includeXi'), oo.includeXi = 1; end   
    % combine cavity and elastic sheet depths hs and he to give one flux qs
    if ~isfield(oo,'combine_sheet'), oo.combine_sheet = 0; end   
    % use arithmetic means for permeability rather than harmonic means
    if ~isfield(oo,'mean_perms'), oo.mean_perms = 0; end                 
     % evaluate expanded solution variables in vv2
    if ~isfield(oo,'evaluate_variables'), oo.evaluate_variables = 0; end   
    % evaluate residual [ otherwise returns empty F ]
    if ~isfield(oo,'evaluate_residual'), oo.evaluate_residual = 1; end 
     % evaluate Jacobian [ otherwise returns emtpy J ]
    if ~isfield(oo,'evaluate_jacobian'), oo.evaluate_jacobian = 1; end     
    % include_blister
    if ~isfield(oo,'include_blister'), oo.include_blister = 1; end 
    % include_radius, if True, radius is included in the Jacobian, if not,
    % radius is calculated outstide the solver
    if ~isfield(oo,'include_pressure'), oo.include_pressure = 1; end 
    % formulation of the permeability
    if ~isfield(oo,'constant_k'), oo.constant_k = 1; end 
    if ~isfield(oo,'sheet_k'), oo.sheet_k = 0; end 
    if ~isfield(oo,'channel_k'), oo.channel_k = 0; end 
    
    %% fill in missing boundary fluxes
    if ~isfield(aa,'phi'), aa.phi = aa.phi_a(gg.nbdy); end

    if ~isfield(aa,'qsx'), aa.qsx = 0*ones(length(gg.ebdy),1); end
    if ~isfield(aa,'qex'), aa.qex = 0*ones(length(gg.ebdy),1); end
    if ~isfield(aa,'Qx'), aa.Qx = 0*ones(length(gg.ebdy),1); end
    
    if ~isfield(aa,'qsy'), aa.qsy = 0*ones(length(gg.fbdy),1); end
    if ~isfield(aa,'qey'), aa.qey = 0*ones(length(gg.fbdy),1); end
    if ~isfield(aa,'Qy'), aa.Qy = 0*ones(length(gg.fbdy),1); end
    
    if ~isfield(aa,'Qs'), aa.Qs = 0*ones(length(gg.cbdy),1); end
    if ~isfield(aa,'Qr'), aa.Qr = 0*ones(length(gg.cbdy),1); end

    if ~isfield(aa,'qbx'), aa.qbx = 0*ones(length(gg.ebdy),1); end   % flow within blister
    if ~isfield(aa,'qby'), aa.qby = 0*ones(length(gg.fbdy),1); end   % flow within blister
    
    %% restrict mean operators 
    % [ to take means over active points must take care to only 
    % numptiply with such points ]

    temp = gg.nmeanx(:,gg.ein)*ones(length(gg.ein),1); temp(temp==0) = inf; 
    % sparse(1:length(temp),1:length(temp),temp.^(-1),...
    % length(temp),length(temp)) (nIxnJ,nIxnJ) diagonal maxtrix
    gg.nmeanxin = sparse(1:length(temp),1:length(temp),temp.^(-1),...
                         length(temp),length(temp))*gg.nmeanx;
    temp = gg.nmeany(:,gg.fin)*ones(length(gg.fin),1); temp(temp==0) = inf; 
    gg.nmeanyin = sparse(1:length(temp),1:length(temp),temp.^(-1),...
                         length(temp),length(temp))*gg.nmeany;
    temp = gg.nmeans(:,gg.cin)*ones(length(gg.cin),1); temp(temp==0) = inf; 
    gg.nmeansin = sparse(1:length(temp),1:length(temp),temp.^(-1),...
                         length(temp),length(temp))*gg.nmeans;
    temp = gg.nmeanr(:,gg.cin)*ones(length(gg.cin),1); temp(temp==0) = inf; 
    gg.nmeanrin = sparse(1:length(temp),1:length(temp),temp.^(-1),...
                         length(temp),length(temp))*gg.nmeanr;

    %% primary variables
    phi = vv.phi;
    hs = vv.hs;
    Sx = vv.Sx;
    Sy = vv.Sy;
    Ss = vv.Ss;
    Sr = vv.Sr;
    hb = vv.hb;  % blister thickness
    pb = vv.pb;  % blister pressure
    
    % channel cross sectional area interpolated at nodes
    Smean = 0.25*(gg.nmeanx*Sx+ gg.nmeany*Sy + gg.nmeans*Ss + gg.nmeanr*Sr);

    %% prescribed fields
    if isfield(vv,'m'), m = vv.m; else m = aa.m; end
    if isfield(vv,'E'), E = vv.E; else E = aa.E; end
    if isfield(vv,'Qb_in'), Qb_in = vv.Qb; else Qb_in = aa.Qb_in; end
    % if isfield(vv,'Ql'), Ql = vv.Ql; else Ql = aa.Ql; end
    if isfield(vv,'Ub'), Ub = vv.Ub; else Ub = aa.Ub; end
    sigma = aa.sigma;

    % boundary potential
    if ~isempty(gg.nbdy)
        phi(gg.nbdy) = aa.phi;
        pb(gg.nbdy) = aa.phi_0(gg.nbdy)-aa.phi_a(gg.nbdy);          % boundary pb = ice overburden
    end
    
    % potential gradients
    Psi_x = gg.eddx*phi;
    Psi_y = gg.fddy*phi;
    Psi_s = gg.cdds*phi;
    Psi_r = gg.cddr*phi;
    Psib_x = gg.eddx*pb;
    Psib_y = gg.fddy*pb;

    % bed potential gradients
    Psi_a_x = gg.eddx*aa.phi_a;
    Psi_a_y = gg.fddy*aa.phi_a;
    Psi_a_s = gg.cdds*aa.phi_a;
    Psi_a_r = gg.cddr*aa.phi_a;
    
    % total potential gradient on nodes [ mean of psi on adjoining edges; 
    % could use meany2 here if reflective boundary conditions ]
    % where nmeany is 1 is for reflective boundary conditions 
    % want to average Psi_y and -Psi_y in Psi
    % nmeany2 = nmeanyin; 
    % if strcmp(opts.method,'r_stress_l_vel_tb_ref'),...
    % nmeany2(nmeany2==1) = 0; end 

    % ||Psi|| magnitude of Psi
    Psi = ( (gg.nmeanxin(:,gg.ein)*Psi_x(gg.ein)).^2+...
            (gg.nmeanyin(:,gg.fin)*Psi_y(gg.fin)).^2 ).^(1/2);

    % moulin storage Am*delta(xm) (dimensionless)
    sigma_m = zeros(gg.nIJ,1);
    if isfield(pp,'ni_m') && ~isempty(pp.ni_m)
        sigma_m(pp.ni_m) = pp.A_m*ones(length(pp.ni_m),1)./...
                           gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
    end

    %effective pressure
    N = aa.phi_0-phi;

    % elastic sheet
    he = he_fun(aa.phi_0-phi,aa.phi_0-aa.phi_a,pp,oo);

    % storage sheet
    hv= hv_fun(phi-aa.phi_a,aa.phi_0-aa.phi_a,pp.c25*sigma,pp,oo);

    % moulin sheet
    hm = hm_fun(phi-aa.phi_a,aa.phi_0-aa.phi_a,pp.c25*sigma_m,pp,oo);
    
    % % upwind directions
    % eup = gg.econnect(:,1).*(Psi_x<=0) + gg.econnect(:,2).*(Psi_x>0);
    % fup = gg.fconnect(:,1).*(Psi_y<=0) + gg.fconnect(:,2).*(Psi_y>0);
    %     
    % % extrapolate he and hv values on boundary nodes from neighbouring nodes [ or perhaps just upstream nodes ? ]     %mean operator of nodes from surrounding nodes [ then adjust to take mean only over interior nodes ]
    % %mean operator over interior nodes
    % nmeann = ( oo.nmeanx(:,ein)*oo.emean(ein,:)+oo.nmeany(:,fin)*oo.fmean(fin,:) ); nmeann = spdiags(zeros(nIJ,1),0,nmeann); 
    % temp = nmeann(:,nin)*ones(length(nin),1); temp(temp==0) = inf; nmeann = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeann;
    % % mean operator over upstream nodes
    % nmeanupn = ( oo.nmeanx(:,ein)*sparse(1:length(ein),eup(ein),ones(length(ein),1),length(ein),nIJ,length(ein))+oo.nmeany(:,fin)*sparse(1:length(fin),fup(fin),ones(length(fin),1),length(fin),nIJ,length(fin)) ); nmeanupn = spdiags(zeros(nIJ,1),0,nmeanupn); 
    % temp = nmeanupn(:,nin)*ones(length(nin),1); temp(temp==0) = inf; nmeanupn = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*nmeanupn;
    % he(nbdy) = nmeann(nbdy,nin)*he(nin);
    % hv(nbdy) = nmeann(nbdy,nin)*hv(nin); 
    % hm(nbdy) = nmeann(nbdy,nin)*hm(nin); 

    % channel fluxes
    Qx = -pp.c23*max(Sx,0).^(pp.alpha_c).*(Psi_x.^2+pp.Psi_reg^2).^((pp.beta_c-1)/2).*Psi_x;
    Qy = -pp.c23*max(Sy,0).^(pp.alpha_c).*(Psi_y.^2+pp.Psi_reg^2).^((pp.beta_c-1)/2).*Psi_y;
    Qs = -pp.c23*max(Ss,0).^(pp.alpha_c).*(Psi_s.^2+pp.Psi_reg^2).^((pp.beta_c-1)/2).*Psi_s;
    Qr = -pp.c23*max(Sr,0).^(pp.alpha_c).*(Psi_r.^2+pp.Psi_reg^2).^((pp.beta_c-1)/2).*Psi_r;

    % sheet permeabilites h^alpha_c [ edited from hydro_timestep_diag ]
    perm_reg = 1e-16;
    if oo.combine_sheet
        % [ include hs and he together as sheet permeability perms and 
        %   have perme set to zero ]    
        hh = hs+he; alf = pp.alpha_s; bet = pp.beta_s;
        if oo.mean_perms % arithmetic mean of  (phi_1*dx1+phi_2*dx2)/(dx1+dx2) [Question here]
            permsx = (gg.Dx(gg.econnect(:,1)).*max(hh(gg.econnect(:,1)),0).^(alf) + ...
                gg.Dx(gg.econnect(:,2)).*max(hh(gg.econnect(:,2)),0).^(alf)).*...
                (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(-1);
            permsy = (gg.Dy(gg.fconnect(:,1)).*max(hh(gg.fconnect(:,1)),0).^(alf) + ...
                gg.Dy(gg.fconnect(:,2)).*max(hh(gg.fconnect(:,2)),0).^(alf)).*...
                (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(-1);
        else % harmonic mean
            permsx = (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(bet).*...
                max(hh(gg.econnect(:,1)),0).^(alf).*max(hh(gg.econnect(:,2)),0).^(alf).*...
                (gg.Dx(gg.econnect(:,1)).*max(hh(gg.econnect(:,2)),0).^(alf/bet) + ...
                gg.Dx(gg.econnect(:,2)).*max(hh(gg.econnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 
            permsy = (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(bet).*...
                max(hh(gg.fconnect(:,1)),0).^(alf).*max(hh(gg.fconnect(:,2)),0).^(alf).*...
                (gg.Dy(gg.fconnect(:,1)).*max(hh(gg.fconnect(:,2)),0).^(alf/bet) + ...
                gg.Dy(gg.fconnect(:,2)).*max(hh(gg.fconnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 
        end
            permex = zeros(gg.eIJ,1);
            permey = zeros(gg.fIJ,1);
    else
        % [ compute sheet permeabilities separately as perms and perme ]
        hh = hs; alf = pp.alpha_s; bet = pp.beta_s;
        if oo.mean_perms
            permsx = (gg.Dx(gg.econnect(:,1)).*max(hh(gg.econnect(:,1)),0).^(alf) +...
                gg.Dx(gg.econnect(:,2)).*max(hh(gg.econnect(:,2)),0).^(alf)).*...
                (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(-1);
            permsy = (gg.Dy(gg.fconnect(:,1)).*max(hh(gg.fconnect(:,1)),0).^(alf) +...
                gg.Dy(gg.fconnect(:,2)).*max(hh(gg.fconnect(:,2)),0).^(alf) ).*...
                (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(-1);
        else
            permsx = (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(bet).*...
                max(hh(gg.econnect(:,1)),0).^(alf).*max(hh(gg.econnect(:,2)),0).^(alf).*...
                ( gg.Dx(gg.econnect(:,1)).*max(hh(gg.econnect(:,2)),0).^(alf/bet) + ...
                gg.Dx(gg.econnect(:,2)).*max(hh(gg.econnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 
            permsy = (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(bet).*...
                max(hh(gg.fconnect(:,1)),0).^(alf).*max(hh(gg.fconnect(:,2)),0).^(alf).*...
                ( gg.Dy(gg.fconnect(:,1)).*max(hh(gg.fconnect(:,2)),0).^(alf/bet) + ...
                gg.Dy(gg.fconnect(:,2)).*max(hh(gg.fconnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 
        end
        hh = he; alf = pp.alpha_e; bet = pp.beta_e;
        if oo.mean_perms
            permex = ( gg.Dx(gg.econnect(:,1)).*max(hh(gg.econnect(:,1)),0).^(alf) + ...
                gg.Dx(gg.econnect(:,2)).*max(hh(gg.econnect(:,2)),0).^(alf) ).*...
                (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(-1);
            permey = ( gg.Dy(gg.fconnect(:,1)).*max(hh(gg.fconnect(:,1)),0).^(alf) + ...
                gg.Dy(gg.fconnect(:,2)).*max(hh(gg.fconnect(:,2)),0).^(alf) ).*...
                (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(-1);
        else
            permex = (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(bet).*...
                max(hh(gg.econnect(:,1)),0).^(alf).*max(hh(gg.econnect(:,2)),0).^(alf).*...
                ( gg.Dx(gg.econnect(:,1)).*max(hh(gg.econnect(:,2)),0).^(alf/bet) + ...
                gg.Dx(gg.econnect(:,2)).*max(hh(gg.econnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 
            permey = (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(bet).*...
                max(hh(gg.fconnect(:,1)),0).^(alf).*max(hh(gg.fconnect(:,2)),0).^(alf).*...
                ( gg.Dy(gg.fconnect(:,1)).*max(hh(gg.fconnect(:,2)),0).^(alf/bet) + ...
                gg.Dy(gg.fconnect(:,2)).*max(hh(gg.fconnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 
        end
    end

    % sheet fluxes [ qsx,qsy cavity sheet flux, qex, qey poroelastic sheet flow ]
    qsx = -pp.c21*permsx.*((gg.emean(:,:)*Psi).^2+pp.Psi_reg^2).^((pp.beta_s-1)/2).*Psi_x;
    qsy = -pp.c21*permsy.*((gg.fmean(:,:)*Psi).^2+pp.Psi_reg^2).^((pp.beta_s-1)/2).*Psi_y;
    qex = -pp.c22*permex.*((gg.emean(:,:)*Psi).^2+pp.Psi_reg^2).^((pp.beta_e-1)/2).*Psi_x;
    qey = -pp.c22*permey.*((gg.fmean(:,:)*Psi).^2+pp.Psi_reg^2).^((pp.beta_e-1)/2).*Psi_y;

    % blister fluxes
    alf = 3.0; bet = 1.0;
    if oo.mean_perms
        permbx = (gg.Dx(gg.econnect(:,1)).*max(hb(gg.econnect(:,1)),0).^alf +...
            gg.Dx(gg.econnect(:,2)).*max(hb(gg.econnect(:,2)),0).^alf).*...
            (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(-1);
        permby = ( gg.Dy(gg.fconnect(:,1)).*max(hb(gg.fconnect(:,1)),0).^(alf) + ...
            gg.Dy(gg.fconnect(:,2)).*max(hb(gg.fconnect(:,2)),0).^(alf) ).*...
            (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(-1);
    else
        permbx = (gg.Dx(gg.econnect(:,1))+gg.Dx(gg.econnect(:,2))).^(bet).*...
                max(hb(gg.econnect(:,1)),0).^(alf).*max(hb(gg.econnect(:,2)),0).^(alf).*...
                ( gg.Dx(gg.econnect(:,1)).*max(hb(gg.econnect(:,2)),0).^(alf/bet) + ...
                gg.Dx(gg.econnect(:,2)).*max(hb(gg.econnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 
        permby = (gg.Dy(gg.fconnect(:,1))+gg.Dy(gg.fconnect(:,2))).^(bet).*...
            max(hb(gg.fconnect(:,1)),0).^(alf).*max(hb(gg.fconnect(:,2)),0).^(alf).*...
            ( gg.Dy(gg.fconnect(:,1)).*max(hb(gg.fconnect(:,2)),0).^(alf/bet) + ...
            gg.Dy(gg.fconnect(:,2)).*max(hb(gg.fconnect(:,1)),0).^(alf/bet) + perm_reg).^(-bet); 

    end
    % permbx, permby are defined on all edges
    qbx = -pp.c42*permbx.*Psib_x; % vector of (gg.eIJ,1)
    qby = -pp.c42*permby.*Psib_y; % vector of (gg.fIJ,1)
    
    % blister to subglacial drainage system term defined on the nodes [ns]
    Qb = pp.c50*hb; % 

    % boundary edge fluxes
    if ~isempty(gg.ebdy)
        qsx(gg.ebdy) = aa.qsx;
        qex(gg.ebdy) = aa.qex;
        qbx(gg.ebdy) = aa.qbx;
        Qx(gg.ebdy) = aa.Qx;
    end

    if ~isempty(gg.fbdy)
        qsy(gg.fbdy) = aa.qsy;
        qey(gg.fbdy) = aa.qey;
        qby(gg.fbdy) = aa.qby;
        Qy(gg.fbdy) = aa.Qy;
    end

    if ~isempty(gg.cbdy)
        Qs(gg.cbdy) = aa.Qs;
        Qr(gg.cbdy) = aa.Qr;
    end
    
    % outside edge fluxes
    qsx(gg.eout) = 0; qex(gg.eout) = 0; qbx(gg.eout) = 0; Qx(gg.eout) = 0; % x-edge
    qsy(gg.fout) = 0; qey(gg.fout) = 0; qby(gg.fout) = 0; Qy(gg.fout) = 0; % y-edge  
    Qs(gg.cout) = 0; Qr(gg.cout) = 0; % corners
    
    % Xi = rho_w * L * M
    % channel dissipation
    Xicx = - pp.c19*Qx.*(Psi_x*(1-pp.c41)+pp.c41*Psi_a_x);
    Xicy = - pp.c19*Qy.*(Psi_y*(1-pp.c41)+pp.c41*Psi_a_y);
    Xics = - pp.c19*Qs.*(Psi_s*(1-pp.c41)+pp.c41*Psi_a_s);
    Xicr = - pp.c19*Qr.*(Psi_r*(1-pp.c41)+pp.c41*Psi_a_r);

    % sheet dissipation
    if oo.includeXi
        Xi = - gg.nmeanxin(:,gg.ein)*(qsx(gg.ein).*Psi_x(gg.ein) + qex(gg.ein).*Psi_x(gg.ein))... 
          - gg.nmeanyin(:,gg.fin)*(qsy(gg.fin).*Psi_y(gg.fin) + qey(gg.fin).*Psi_y(gg.fin));
    else
        Xi = aa.Xi; 
    end
    if oo.maxXi
        Xix = pp.c20*aa.lcx.*max(Xi(gg.econnect(:,1:2)),[],2);
        Xiy = pp.c20*aa.lcy.*max(Xi(gg.fconnect(:,1:2)),[],2);
        Xis = pp.c20*aa.lcs.*max(Xi([gg.fconnect(gg.cconnect(:,1),1)...
                                     gg.fconnect(gg.cconnect(:,2),2)]),[],2);
        Xir = pp.c20*aa.lcr.*max(Xi([gg.fconnect(gg.cconnect(:,1),2)...
                                     gg.fconnect(gg.cconnect(:,2),1)]),[],2);
    else
        Xix = pp.c20*aa.lcx.*(gg.emean(:,:)*Xi);
        Xiy = pp.c20*aa.lcy.*(gg.fmean(:,:)*Xi);
        Xis = pp.c20*aa.lcs.*(gg.cmean(:,:)*Xi);
        Xir = pp.c20*aa.lcr.*(gg.cmean(:,:)*Xi);
    end
    %% Variables
    if oo.evaluate_variables
        vv2 = vv;
        vv2.N = N;
        vv2.he = he;
        vv2.hv = hv;
        vv2.hm = hm;
        vv2.hb = hb;
        vv2.pb = pb;
        % vv2.Sx_n = gg.nmeanx*Sx;
        % vv2.Sy_n = gg.nmeany*Sy;
        % vv2.Ss_n = gg.nmeans*Ss;
        % vv2.Sr_n = gg.nmeanr*Sr;

        % fluxes and potential gradient
        vv2.qsx = qsx;
        vv2.qsy = qsy;
        vv2.qex = qex;
        vv2.qey = qey;
        vv2.qbx = qbx;
        vv2.qby = qby;

        vv2.Qx = Qx;
        vv2.Qy = Qy;
        vv2.Qs = Qs;
        vv2.Qr = Qr;
        vv2.Psi_x = Psi_x;
        vv2.Psi_y = Psi_y;
        vv2.Psi_s = Psi_s;
        vv2.Psi_r = Psi_r;
        vv2.Psi = Psi;
        vv2.Xix = Xix;
        vv2.Xiy = Xiy;
        vv2.Xis = Xis;
        vv2.Xir = Xir;
        vv2.Xi = Xi;
        vv2.Xicx = Xicx;
        vv2.Xicy = Xicy;
        vv2.Xics = Xics;
        vv2.Xicr = Xicr;
        % influx from boundaries
        vv2.Q_in = sum(Qx(gg.ebdy))+sum(Qy(gg.fbdy))+sum(Qs(gg.cbdy))+sum(Qr(gg.cbdy)) ...
                + pp.c4/pp.c9*(sum(qsx(gg.ebdy).*(gg.emean(gg.ebdy,:)*gg.Dy))+...
                               sum(qsy(gg.fbdy).*(gg.fmean(gg.fbdy,:)*gg.Dx))) ...
                + pp.c5/pp.c9*(sum(qex(gg.ebdy).*(gg.emean(gg.ebdy,:)*gg.Dy))+...
                               sum(qey(gg.fbdy).*(gg.fmean(gg.fbdy,:)*gg.Dx)));
    else
        vv2 = struct;
    end

    %% Residual
    if oo.evaluate_residual
        [R1,R2,R3,R4,R5,R6,R7,R8] = residuals();
        
        vv2.R_bdy = R2(gg.nbdy);
        vv2.Qb_in = Qb_in;
        vv2.Qb_dec = sum(Qb.*gg.Dx.*gg.Dy);
        vv2.Q_out = 1/pp.c9*sum( R2(gg.nbdy).*gg.Dx(gg.nbdy).*gg.Dy(gg.nbdy) ); 
        vv2.Qb_out = 1/pp.c43*sum( R8(gg.nbdy).*gg.Dx(gg.nbdy).*gg.Dy(gg.nbdy) ); 
        vv2.Xi = Xi;
        vv2.he = he;
        vv2.Q_in = sum(Qx(gg.ebdy))+sum(Qy(gg.fbdy))+sum(Qs(gg.cbdy))+sum(Qr(gg.cbdy)) ...
                 + pp.c4/pp.c9*(sum(qsx(gg.ebdy).*(gg.emean(gg.ebdy,:)*gg.Dy))+...
                                sum(qsy(gg.fbdy).*(gg.fmean(gg.fbdy,:)*gg.Dx))) ...
                 + pp.c5/pp.c9*(sum(qex(gg.ebdy).*(gg.emean(gg.ebdy,:)*gg.Dy))+...
                                sum(qey(gg.fbdy).*(gg.fmean(gg.fbdy,:)*gg.Dx)));% inflow

        vv2.Q_outQ = -1/pp.c9*sum(pp.c9*((gg.nddx(gg.nbdy,:)*Qx).*gg.Dy(gg.nbdy).^(-1) + ...
                                 (gg.nddy(gg.nbdy,:)*Qy).*gg.Dx(gg.nbdy).^(-1) + ...
                                 (gg.ndds(gg.nbdy,:)*Qs) + (gg.nddr(gg.nbdy,:)*Qr)).*...
                                gg.Dx(gg.nbdy).*gg.Dy(gg.nbdy));  % outflow in channels

        vv2.Q_outq = -1/pp.c9*sum((pp.c4*(gg.nddx(gg.nbdy,:)*qsx + gg.nddy(gg.nbdy,:)*qsy) +...
                                   pp.c5*( gg.nddx(gg.nbdy,:)*qex + gg.nddy(gg.nbdy,:)*qey)).*...
                                gg.Dx(gg.nbdy).*gg.Dy(gg.nbdy));  % outflow in sheet

        vv2.Q_outb = -1/pp.c43*sum((pp.c45*(gg.nddx(gg.nbdy,:)*qbx + gg.nddy(gg.nbdy,:)*qby)).*...
                                gg.Dx(gg.nbdy).*gg.Dy(gg.nbdy));   % outflow in blister
        
        F1 = R1(gg.ns);  % cavity sheet
        F2 = R2(gg.nin); % total mass conservation
        F3 = R3(gg.ein); % channel
        F4 = R4(gg.fin);
        F5 = R5(gg.cin);
        F6 = R6(gg.cin);
        F7 = R7(gg.ns);  % blister thickness
        F8 = R8(gg.nin); % blister mass conservation
        
        if oo.no_channels && oo.no_sheet 
            F =  F2;
        elseif oo.no_channels
            F = [ F1; F2 ];
        elseif oo.no_sheet
            F = [ F2; F3; F4; F5; F6 ];
        else
            if oo.include_blister
                if oo.include_pressure  
                    F = [ F1; F2; F3; F4; F5; F6; F7; F8];
                else
                    F = [ F1; F2; F3; F4; F5; F6; F7];
                end
            else
                F = [ F1; F2; F3; F4; F5; F6];
            end
        end
    else
        F = []; F1 = []; F2 = []; F3 = []; F4 = []; F5 = []; F6 = []; F7 = []; F8 = [];
    end

    %% Jacobian
    if oo.evaluate_jacobian
        J = jacob;
    else
        J = [];
    end

    %% nested functions
    function [R1,R2,R3,R4,R5,R6,R7,R8] = residuals()    
        %% time derivatives [ for use in residuals below ]

        hs_t = pp.c26*m ...
                + pp.c16*Ub.*max(1-pp.c17*hs,0) ...
                - pp.c18*hs.*abs(aa.phi_0-phi).^(pp.n_Glen-1).*(aa.phi_0-phi);  
            % pp.c15*hs_t
        
        h_t = - pp.c4*(gg.nddx(:,:)*qsx + gg.nddy(:,:)*qsy) ... % cavity sheet flux divergence
              - pp.c5*(gg.nddx(:,:)*qex + gg.nddy(:,:)*qey) ... % poroelastic sheet flux divergence
              + pp.c6*m ...                                     % basal melting rate
              + pp.c7*E ...                                     % moulin and lake influx
              + Qb      ...                                     % blister influx
              - pp.c9*((gg.nddx(:,:)*Qx).*gg.Dy.^(-1) + (gg.nddy(:,:)*Qy).*gg.Dx.^(-1)) ... % x,y channel divergence
              - pp.c9*((gg.ndds(:,:)*Qs) + (gg.nddr(:,:)*Qr)) ...                           % s,r channel divergence
              + pp.c11*((gg.nmeanxin(:,gg.ein)*(Xicx(gg.ein)+Xix(gg.ein))).*gg.Dy.^(-1) +...% x,y dissipation
                        (gg.nmeanyin(:,gg.fin)*(Xicy(gg.fin)+Xiy(gg.fin))).*gg.Dx.^(-1)) ...
              + pp.c11*((gg.nmeansin(:,gg.cin)*(Xics(gg.cin)+Xis(gg.cin))).*gg.Ds.*gg.Dx.^(-1).*gg.Dy.^(-1) +... % s,r dissipation
                        (gg.nmeanrin(:,gg.cin)*(Xicr(gg.cin)+Xir(gg.cin))).*gg.Dr.*gg.Dx.^(-1).*gg.Dy.^(-1));    
            % pp.c1*hs_t + pp.c2*hd_t + pp.c3*hv_t + pp.c3*hm_t + pp.c8*channels
    
        Sx_t = - pp.c14*Sx.*abs(gg.emean(:,:)*(aa.phi_0-phi)).^(pp.n_Glen-1).*(gg.emean(:,:)*(aa.phi_0-phi)) ... % ice creep closure
               + pp.c13*(Xicx+Xix) ...  % wall melting opening
               + pp.c35*aa.lcx.*(gg.emean(:,:)*Ub).*max(1-pp.c36*Sx,0); % 
            % pp.c12*Sx_t
    
        Sy_t = - pp.c14*Sy.*abs(gg.fmean(:,:)*(aa.phi_0-phi)).^(pp.n_Glen-1).*(gg.fmean(:,:)*(aa.phi_0-phi)) ...
                + pp.c13*(Xicy+Xiy) ...
                + pp.c35*aa.lcy.*(gg.fmean(:,:)*Ub).*max(1-pp.c36*Sy,0);
            % pp.c12*Sy_t
    
        Ss_t = - pp.c14*Ss.*abs(gg.cmean(:,:)*(aa.phi_0-phi)).^(pp.n_Glen-1).*(gg.cmean(:,:)*(aa.phi_0-phi)) ...
                + pp.c13*(Xics+Xis) ...
                + pp.c35*aa.lcs.*(gg.cmean(:,:)*Ub).*max(1-pp.c36*Ss,0);
            % pp.c12*Ss_t
    
        Sr_t = - pp.c14*Sr.*abs(gg.cmean(:,:)*(aa.phi_0-phi)).^(pp.n_Glen-1).*(gg.cmean(:,:)*(aa.phi_0-phi)) ...
                + pp.c13*(Xicr+Xir) ...
                + pp.c35*aa.lcr.*(gg.cmean(:,:)*Ub).*max(1-pp.c36*Sr,0);
            % pp.c12*Sr_t

        % blister volume
        hb_t = -pp.c45*(gg.nddx(:,:)*qbx + gg.nddy(:,:)*qby) + pp.c43*Qb_in./gg.Dx./gg.Dy -pp.c44*hb;

        % display Vb and Rb in the calculation
        % disp([pp.c42*aa.Qb_in(pp.ni_l); -pp.c43*Vb(pp.ni_l)./(Rb(pp.ni_l)+pp.R_b_reg).^2; -pp.c44*(aa.phi_0(pp.ni_l)-phi(pp.ni_l)).*Rb(pp.ni_l)]);
        % disp([Vb(pp.ni_l) Rb(pp.ni_l)]);

        % blister radius
        % 0 = - pb + (aa.phi_0-aa.phi_a) + pp.c49 * (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*(gg.nddx*gg.eddx + gg.nddy*gg.fddy)*hb;
        % Rb_t = 0;
        
        %% old variables
        hs_old = vv0.hs;
        phi_old = vv0.phi;
        he_old = he_fun(aa.phi_0-phi_old,aa.phi_0-aa.phi_a,pp,oo);
        hv_old = hv_fun(phi_old-aa.phi_a,aa.phi_0-aa.phi_a,pp.c25*sigma,pp,oo);
        hm_old = hm_fun(phi_old-aa.phi_a,aa.phi_0-aa.phi_a,pp.c25*sigma_m,pp,oo);
        Sx_old = vv0.Sx;
        Sy_old = vv0.Sy;
        Ss_old = vv0.Ss;
        Sr_old = vv0.Sr;
        hb_old = vv0.hb;
        
        %% objective residuals    
        R1 = - pp.c15*(hs-hs_old).*dt^(-1) + hs_t; % equation (3), cavity sheet variation
        % equation (10), total conservation of mass
        R2 =  - pp.c1*(hs-hs_old).*dt^(-1) ...
             - pp.c2*(he-he_old).*dt.^(-1) ...
             - pp.c3*(hv-hv_old).*dt.^(-1) ...
             - pp.c3*(hm-hm_old).*dt.^(-1) ...
               - pp.c8*gg.nmeanxin(:,gg.ein)*(Sx(gg.ein)-Sx_old(gg.ein)).*dt^(-1).*gg.Dy.^(-1) ...
               - pp.c8*gg.nmeanyin(:,gg.fin)*(Sy(gg.fin)-Sy_old(gg.fin)).*dt^(-1).*gg.Dx.^(-1) ...
               - pp.c8*gg.nmeansin(:,gg.cin)*(Ss(gg.cin)-Ss_old(gg.cin)).*dt^(-1).*gg.Ds.*gg.Dx.^(-1).*gg.Dy.^(-1) ...
               - pp.c8*gg.nmeanrin(:,gg.cin)*(Sr(gg.cin)-Sr_old(gg.cin)).*dt^(-1).*gg.Dr.*gg.Dx.^(-1).*gg.Dy.^(-1) ...
             + h_t;
        R3 = - pp.c12*(Sx-Sx_old).*dt^(-1) + Sx_t; % equation (7), channel croess sectional area
        R4 = - pp.c12*(Sy-Sy_old).*dt^(-1) + Sy_t;
        R5 = - pp.c12*(Ss-Ss_old).*dt^(-1) + Ss_t;
        R6 = - pp.c12*(Sr-Sr_old).*dt^(-1) + Sr_t;
        R7 = - pb + (aa.phi_0-aa.phi_a) + pp.c49 * (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*(gg.nddx*gg.eddx + gg.nddy*gg.fddy)*hb; % blister pressure
        R8 =        - (hb-hb_old).*dt^(-1) + hb_t; % mass conservation
        
 
    end
    function J = jacob()
        % [ taken from hydro_timestep_diag 13 August 2014, 
        % initial definitions changed to make compatible : should probably update in future ]
        %% options
        opts = oo;
        opts.include_diag = 1;
    
        %% prescribed fields
        lcx = aa.lcx;
        lcy = aa.lcy;
        lcs = aa.lcs;
        lcr = aa.lcr;
        phi_0 = aa.phi_0;
        phi_a = aa.phi_a;
    
        %% extract grid and operators
        nIJ = gg.nIJ;
        eIJ = gg.eIJ;
        fIJ = gg.fIJ;
        nin = gg.nin;
        ein = gg.ein;
        fin = gg.fin;

        nbdy = gg.nbdy;
        ebdy = gg.ebdy;
        fbdy = gg.fbdy;

        ns = gg.ns;
        econnect = gg.econnect;
        fconnect = gg.fconnect;
        cconnect = gg.cconnect;

        nddx = gg.nddx;
        nddy = gg.nddy;
        eddx = gg.eddx;
        fddy = gg.fddy;

        emean = gg.emean;
        fmean = gg.fmean;

        Dx = gg.Dx;
        Dy = gg.Dy;
        Ds = gg.Ds;
        Dr = gg.Dr;
        cIJ = gg.cIJ;
        cin = gg.cin;
        cbdy = gg.cbdy;
        cdds = gg.cdds;
        cddr = gg.cddr;
        ndds = gg.ndds;
        nddr = gg.nddr;
        cmean = gg.cmean;
        nmeanx = gg.nmeanxin;
        nmeany = gg.nmeanyin;
        nmeans = gg.nmeansin;
        nmeanr = gg.nmeanrin;
        nmeany2 = nmeany;
    
        %% extract parameters
        c1 = pp.c1;
        c2 = pp.c2;
        c3 = pp.c3;
        c4 = pp.c4;
        c5 = pp.c5;
        c6 = pp.c6;
        c7 = pp.c7;
        c8 = pp.c8;
        c9 = pp.c9;
        % c10 = pp.c10;
        c11 = pp.c11;
        c12 = pp.c12;
        c13 = pp.c13;
        c14 = pp.c14;
        c15 = pp.c15;
        c16 = pp.c16;
        c17 = pp.c17;
        c18 = pp.c18;
        c19 = pp.c19;
        c20 = pp.c20;
        c21 = pp.c21;
        c22 = pp.c22;
        c23 = pp.c23;
        % c24 = pp.c24;
        c25 = pp.c25; 
        c26 = pp.c26;
        c35 = pp.c35;
        c36 = pp.c36;
        c41 = pp.c41;
        c42 = pp.c42;
        c43 = pp.c43;
        c44 = pp.c44;
        c45 = pp.c45;
        c48 = pp.c48;
        c49 = pp.c49;
        c50 = pp.c50;

        n_Glen = pp.n_Glen;
        alpha_c = pp.alpha_c;
        beta_c = pp.beta_c;
        alpha_s = pp.alpha_s;
        beta_s = pp.beta_s;
        alpha_e = pp.alpha_e;
        beta_e = pp.beta_e;
        Psi_reg = pp.Psi_reg;
        
        %% extract old variables
        hb_old = vv0.hb;
        pb_old = vv0.pb;

        %% adjust size of perms
        permsx = permsx(ein);
        permex = permex(ein);
        permsy = permsy(fin);
        permey = permey(fin);
        permbx = permbx(ein);
        permby = permby(fin);

        %derivative of elastic sheet
        Dhe_phi = sparse(nin,nin, -Dhe_fun(phi_0(nin)-phi(nin), phi_0(nin)-phi_a(nin), pp,opts),nIJ, nIJ, length(nin));
        %derivative of storage sheet
        Dhv_phi = sparse(nin,nin, Dhv_fun(phi(nin)-phi_a(nin),phi_0(nin)-phi_a(nin),c25*sigma(nin),pp,opts), nIJ, nIJ, length(nin));
        %derivative of moulin sheet
        Dhm_phi = sparse(nin,nin, Dhm_fun(phi(nin)-phi_a(nin),phi_0(nin)-phi_a(nin),c25*sigma_m(nin),pp,opts), nIJ, nIJ, length(nin));
        %derivative of Psi 
        DPsi_phi = sparse(1:nIJ,1:nIJ,(nmeanx(:,ein)*Psi_x(ein)).*(Psi.^2+Psi_reg^2).^(-1/2),nIJ,nIJ)...
                   *(nmeanx(:,ein)*eddx(ein,:)) + ...
                   sparse(1:nIJ,1:nIJ,(nmeany2(:,fin)*Psi_y(fin)).*(Psi.^2+Psi_reg^2).^(-1/2),nIJ,nIJ)...
                   *(nmeany2(:,fin)*fddy(fin,:));
        %derivatives of Qx and Qy           
        DQx_phi = sparse(1:length(ein),1:length(ein),...
                  -c23*max(Sx(ein),0).^(alpha_c).*(Psi_x(ein).^2+Psi_reg^2).^((beta_c-1)/2)...
                  .*( beta_c-(beta_c-1)*Psi_reg^2*(Psi_x(ein).^2+Psi_reg^2).^(-1) ), length(ein),length(ein)) * eddx(ein,:);
        DQy_phi = sparse(1:length(fin),1:length(fin),...
                  -c23*max(Sy(fin),0).^(alpha_c).*(Psi_y(fin).^2+Psi_reg^2).^((beta_c-1)/2)...
                  .*( beta_c-(beta_c-1)*Psi_reg^2*(Psi_y(fin).^2+Psi_reg^2).^(-1) ), length(fin),length(fin)) * fddy(fin,:);

        DQx_Sx = sparse(1:length(ein),ein,-c23*max(Sx(ein),0).^(alpha_c-1).*(Psi_x(ein).^2+Psi_reg^2).^((beta_c-1)/2).*Psi_x(ein)*alpha_c ,length(ein),eIJ); 
        DQy_Sy = sparse(1:length(fin),fin,-c23*max(Sy(fin),0).^(alpha_c-1).*(Psi_y(fin).^2+Psi_reg^2).^((beta_c-1)/2).*Psi_y(fin)*alpha_c ,length(fin),fIJ); 
        
        %derivatives of Xi 
        DXicx_phi = - c19*(sparse(1:length(ein),1:length(ein),(1-c41)*Qx(ein),length(ein),length(ein))*eddx(ein,:) + ...
                        sparse(1:length(ein),1:length(ein),Psi_x(ein)*(1-c41)+c41*Psi_a_x(ein),length(ein),length(ein))*(DQx_phi));
        DXicy_phi = - c19*(sparse(1:length(fin),1:length(fin),(1-c41)*Qy(fin),length(fin),length(fin))*fddy(fin,:) + ...
                        sparse(1:length(fin),1:length(fin),Psi_y(fin)*(1-c41)+c41*Psi_a_y(fin),length(fin),length(fin))*(DQy_phi));
        DXicx_Sx = - c19*( sparse(1:length(ein),1:length(ein),(Psi_x(ein)*(1-c41)+Psi_a_x(ein)*c41),length(ein),length(ein))*(DQx_Sx) );
        DXicy_Sy = - c19*( sparse(1:length(fin),1:length(fin),(Psi_y(fin)*(1-c41)+Psi_a_y(fin)*c41),length(fin),length(fin))*(DQy_Sy) );

        if opts.include_diag
            %derivatives of corner Q
            DQs_phi = sparse(1:length(cin),1:length(cin),-c23*max(Ss(cin),0).^(alpha_c).*(Psi_s(cin).^2+Psi_reg^2).^((beta_c-1)/2).*( beta_c-(beta_c-1)*Psi_reg^2*(Psi_s(cin).^2+Psi_reg^2).^(-1) ) ,length(cin),length(cin))*cdds(cin,:);
            DQr_phi = sparse(1:length(cin),1:length(cin),-c23*max(Sr(cin),0).^(alpha_c).*(Psi_r(cin).^2+Psi_reg^2).^((beta_c-1)/2).*( beta_c-(beta_c-1)*Psi_reg^2*(Psi_r(cin).^2+Psi_reg^2).^(-1) ) ,length(cin),length(cin))*cddr(cin,:);
            DQs_Ss = sparse(1:length(cin),cin,-c23*max(Ss(cin),0).^(alpha_c-1).*(Psi_s(cin).^2+Psi_reg^2).^((beta_c-1)/2).*Psi_s(cin)*alpha_c ,length(cin),cIJ); 
            DQr_Sr = sparse(1:length(cin),cin,-c23*max(Sr(cin),0).^(alpha_c-1).*(Psi_r(cin).^2+Psi_reg^2).^((beta_c-1)/2).*Psi_r(cin)*alpha_c ,length(cin),cIJ); 
            
            DXics_phi = - c19*( sparse(1:length(cin),1:length(cin),(1-c41)*Qs(cin),length(cin),length(cin))*cdds(cin,:) + sparse(1:length(cin),1:length(cin),Psi_s(cin)*(1-c41)+Psi_a_s(cin)*c41,length(cin),length(cin))*(DQs_phi) );
            DXicr_phi = - c19*( sparse(1:length(cin),1:length(cin),(1-c41)*Qr(cin),length(cin),length(cin))*cddr(cin,:) + sparse(1:length(cin),1:length(cin),Psi_r(cin)*(1-c41)+Psi_a_r(cin)*c41,length(cin),length(cin))*(DQr_phi) );
            DXics_Ss = - c19*( sparse(1:length(cin),1:length(cin),Psi_s(cin)*(1-c41)+Psi_a_s(cin)*c41,length(cin),length(cin))*(DQs_Ss) );
            DXicr_Sr = - c19*( sparse(1:length(cin),1:length(cin),Psi_r(cin)*(1-c41)+Psi_a_r(cin)*c41,length(cin),length(cin))*(DQr_Sr) );
        else
            DQs_phi = sparse(1,1,0,length(cin),nIJ); 
            DQs_Ss = sparse(1,1,0,length(cin),cIJ); 
            DQr_phi = sparse(1,1,0,length(cin),nIJ); 
            DQr_Sr = sparse(1,1,0,length(cin),cIJ); 
            DXics_phi = sparse(1,1,0,length(cin),nIJ); 
            DXics_Ss = sparse(1,1,0,length(cin),cIJ); 
            DXicr_phi = sparse(1,1,0,length(cin),nIJ); 
            DXicr_Sr = sparse(1,1,0,length(cin),cIJ); 
        end

    % derivatives of permeabilities for sheets
        if opts.combine_sheet
            hh = hs+he;
            if opts.mean_perms
                Dpermsx_hs = sparse(1:length(ein),econnect(ein,1),alpha_s*( Dx(econnect(ein,1)).*max(hh(econnect(ein,1)),0).^(alpha_s-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ) ...
                            +sparse(1:length(ein),econnect(ein,2),alpha_s*( Dx(econnect(ein,2)).*max(hh(econnect(ein,2)),0).^(alpha_s-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ);
                Dpermsy_hs = sparse(1:length(fin),fconnect(fin,1),alpha_s*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,1)),0).^(alpha_s-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ) ...
                            +sparse(1:length(fin),fconnect(fin,2),alpha_s*( Dy(fconnect(fin,2)).*max(hh(fconnect(fin,2)),0).^(alpha_s-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ);
            else
                Dpermsx_hs = sparse(1:length(ein),econnect(ein,1),alpha_s*Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s+alpha_s).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hh(econnect(ein,1)),0).^(alpha_s-1).*( Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ) ...
                            +sparse(1:length(ein),econnect(ein,2),alpha_s*Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s+alpha_s).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hh(econnect(ein,2)),0).^(alpha_s-1).*( Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ); 
                Dpermsy_hs = sparse(1:length(fin),fconnect(fin,1),alpha_s*Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s+alpha_s).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hh(fconnect(fin,1)),0).^(alpha_s-1).*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ) ...
                            +sparse(1:length(fin),fconnect(fin,2),alpha_s*Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s+alpha_s).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hh(fconnect(fin,2)),0).^(alpha_s-1).*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ); 
            end
            Dpermsx_he = Dpermsx_hs;
            Dpermsy_he = Dpermsy_hs;
            Dpermex_he = sparse(1,1,0,length(ein),nIJ); 
            Dpermey_he = sparse(1,1,0,length(fin),nIJ); 
        else
            if opts.mean_perms
                hh = hs;
                Dpermsx_hs = sparse(1:length(ein),econnect(ein,1),alpha_s*( Dx(econnect(ein,1)).*max(hh(econnect(ein,1)),0).^(alpha_s-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ) ...
                            +sparse(1:length(ein),econnect(ein,2),alpha_s*( Dx(econnect(ein,2)).*max(hh(econnect(ein,2)),0).^(alpha_s-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ);
                Dpermsy_hs = sparse(1:length(fin),fconnect(fin,1),alpha_s*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,1)),0).^(alpha_s-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ) ...
                            +sparse(1:length(fin),fconnect(fin,2),alpha_s*( Dy(fconnect(fin,2)).*max(hh(fconnect(fin,2)),0).^(alpha_s-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ);
                hh = he;
                Dpermex_he = sparse(1:length(ein),econnect(ein,1),alpha_s*( Dx(econnect(ein,1)).*max(hh(econnect(ein,1)),0).^(alpha_s-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ) ...
                            +sparse(1:length(ein),econnect(ein,2),alpha_s*( Dx(econnect(ein,2)).*max(hh(econnect(ein,2)),0).^(alpha_s-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ);
                Dpermey_he = sparse(1:length(fin),fconnect(fin,1),alpha_s*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,1)),0).^(alpha_s-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ) ...
                            +sparse(1:length(fin),fconnect(fin,2),alpha_s*( Dy(fconnect(fin,2)).*max(hh(fconnect(fin,2)),0).^(alpha_s-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ);
            else
            % [ mistake ? changed 29 April 2012 ]
                % Dpermsx_hs = sparse(1:length(ein),econnect(ein,1),alpha_s*Dx(econnect(ein,1)).*max(hs(econnect(ein,2)),0).^(alpha_s/beta_s+alpha_s-1).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hs(econnect(ein,1)),0).^(alpha_s).*( Dx(econnect(ein,1)).*max(hs(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hs(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ) ...
                %              +sparse(1:length(ein),econnect(ein,2),alpha_s*Dx(econnect(ein,2)).*max(hs(econnect(ein,1)),0).^(alpha_s/beta_s+alpha_s-1).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hs(econnect(ein,2)),0).^(alpha_s).*( Dx(econnect(ein,1)).*max(hs(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hs(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ); 
                % Dpermsy_hs = sparse(1:length(fin),fconnect(fin,1),alpha_s*Dy(fconnect(fin,1)).*max(hs(fconnect(fin,2)),0).^(alpha_s/beta_s+alpha_s-1).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hs(fconnect(fin,1)),0).^(alpha_s).*( Dy(fconnect(fin,1)).*max(hs(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hs(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ) ...
                %              +sparse(1:length(fin),fconnect(fin,2),alpha_s*Dy(fconnect(fin,2)).*max(hs(fconnect(fin,1)),0).^(alpha_s/beta_s+alpha_s-1).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hs(fconnect(fin,2)),0).^(alpha_s).*( Dy(fconnect(fin,1)).*max(hs(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hs(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ); 
                % Dpermsx_he = sparse(1,1,0,length(ein),nIJ); 
                % Dpermsy_he = sparse(1,1,0,length(fin),nIJ);
                % Dpermex_he = sparse(1:length(ein),econnect(ein,1),alpha_e*Dx(econnect(ein,1)).*max(he(econnect(ein,2)),0).^(alpha_e/beta_e+alpha_e-1).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_e).*max(he(econnect(ein,1)),0).^(alpha_e).*( Dx(econnect(ein,1)).*max(he(econnect(ein,2)),0).^(alpha_e/beta_e) + Dx(econnect(ein,2)).*max(he(econnect(ein,1)),0).^(alpha_e/beta_e) + perm_reg).^(-beta_e-1),length(ein),nIJ) ...
                %              +sparse(1:length(ein),econnect(ein,2),alpha_e*Dx(econnect(ein,2)).*max(he(econnect(ein,1)),0).^(alpha_e/beta_e+alpha_e-1).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_e).*max(he(econnect(ein,2)),0).^(alpha_e).*( Dx(econnect(ein,1)).*max(he(econnect(ein,2)),0).^(alpha_e/beta_e) + Dx(econnect(ein,2)).*max(he(econnect(ein,1)),0).^(alpha_e/beta_e) + perm_reg).^(-beta_e-1),length(ein),nIJ); 
                % Dpermey_he = sparse(1:length(fin),fconnect(fin,1),alpha_e*Dy(fconnect(fin,1)).*max(he(fconnect(fin,2)),0).^(alpha_e/beta_e+alpha_e-1).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_e).*max(he(fconnect(fin,1)),0).^(alpha_e).*( Dy(fconnect(fin,1)).*max(he(fconnect(fin,2)),0).^(alpha_e/beta_e) + Dy(fconnect(fin,2)).*max(he(fconnect(fin,1)),0).^(alpha_e/beta_e) + perm_reg).^(-beta_e-1),length(fin),nIJ) ...
                %              +sparse(1:length(fin),fconnect(fin,2),alpha_e*Dy(fconnect(fin,2)).*max(he(fconnect(fin,1)),0).^(alpha_e/beta_e+alpha_e-1).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_e).*max(he(fconnect(fin,2)),0).^(alpha_e).*( Dy(fconnect(fin,1)).*max(he(fconnect(fin,2)),0).^(alpha_e/beta_e) + Dy(fconnect(fin,2)).*max(he(fconnect(fin,1)),0).^(alpha_e/beta_e) + perm_reg).^(-beta_e-1),length(fin),nIJ); 
                hh = hs;
                Dpermsx_hs = sparse(1:length(ein),econnect(ein,1),alpha_s*Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s+alpha_s).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hh(econnect(ein,1)),0).^(alpha_s-1).*( Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ) ...
                            +sparse(1:length(ein),econnect(ein,2),alpha_s*Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s+alpha_s).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hh(econnect(ein,2)),0).^(alpha_s-1).*( Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ); 
                Dpermsy_hs = sparse(1:length(fin),fconnect(fin,1),alpha_s*Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s+alpha_s).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hh(fconnect(fin,1)),0).^(alpha_s-1).*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ) ...
                            +sparse(1:length(fin),fconnect(fin,2),alpha_s*Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s+alpha_s).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hh(fconnect(fin,2)),0).^(alpha_s-1).*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ); 
                hh = he;
                Dpermex_he = sparse(1:length(ein),econnect(ein,1),alpha_s*Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s+alpha_s).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hh(econnect(ein,1)),0).^(alpha_s-1).*( Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ) ...
                            +sparse(1:length(ein),econnect(ein,2),alpha_s*Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s+alpha_s).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(beta_s).*max(hh(econnect(ein,2)),0).^(alpha_s-1).*( Dx(econnect(ein,1)).*max(hh(econnect(ein,2)),0).^(alpha_s/beta_s) + Dx(econnect(ein,2)).*max(hh(econnect(ein,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(ein),nIJ); 
                Dpermey_he = sparse(1:length(fin),fconnect(fin,1),alpha_s*Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s+alpha_s).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hh(fconnect(fin,1)),0).^(alpha_s-1).*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ) ...
                            +sparse(1:length(fin),fconnect(fin,2),alpha_s*Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s+alpha_s).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(beta_s).*max(hh(fconnect(fin,2)),0).^(alpha_s-1).*( Dy(fconnect(fin,1)).*max(hh(fconnect(fin,2)),0).^(alpha_s/beta_s) + Dy(fconnect(fin,2)).*max(hh(fconnect(fin,1)),0).^(alpha_s/beta_s) + perm_reg).^(-beta_s-1),length(fin),nIJ); 
            end
            Dpermsx_he = sparse(1,1,0,length(ein),nIJ); 
            Dpermsy_he = sparse(1,1,0,length(fin),nIJ); 
        end

        %derivatives of qs
        Dqsx_phi = sparse(1:length(ein),1:length(ein),-c21*permsx.*((emean(ein,:)*Psi).^2+Psi_reg^2).^((beta_s-1)/2),length(ein),length(ein))*eddx(ein,:) ...
                    + sparse(1:length(ein),1:length(ein),-c21*permsx.*((emean(ein,:)*Psi).^2+Psi_reg^2).^((beta_s-3)/2).*(emean(ein,:)*Psi).*Psi_x(ein)*(beta_s-1),length(ein),length(ein))*emean(ein,:)*DPsi_phi ...
                    + sparse(1:length(ein),1:length(ein),-c21*((emean(ein,:)*Psi).^2+Psi_reg^2).^((beta_s-1)/2).*Psi_x(ein) ,length(ein),length(ein))*Dpermsx_he(:,:)*Dhe_phi(:,:);
        Dqsy_phi = sparse(1:length(fin),1:length(fin),-c21*permsy.*((fmean(fin,:)*Psi).^2+Psi_reg^2).^((beta_s-1)/2),length(fin),length(fin))*fddy(fin,:) ...
                    + sparse(1:length(fin),1:length(fin),-c21*permsy.*((fmean(fin,:)*Psi).^2+Psi_reg^2).^((beta_s-3)/2).*(fmean(fin,:)*Psi).*Psi_y(fin)*(beta_s-1),length(fin),length(fin))*fmean(fin,:)*DPsi_phi ...
                    + sparse(1:length(fin),1:length(fin),-c21*((fmean(fin,:)*Psi).^2+Psi_reg^2).^((beta_s-1)/2).*Psi_y(fin) ,length(fin),length(fin))*Dpermsy_he(:,:)*Dhe_phi(:,:);
        Dqsx_hs = sparse(1:length(ein),1:length(ein),-c21*((emean(ein,:)*Psi).^2+Psi_reg^2).^((beta_s-1)/2).*Psi_x(ein),length(ein),length(ein))*Dpermsx_hs(:,:); 
        Dqsy_hs = sparse(1:length(fin),1:length(fin),-c21*((fmean(fin,:)*Psi).^2+Psi_reg^2).^((beta_s-1)/2).*Psi_y(fin) ,length(fin),length(fin))*Dpermsy_hs(:,:); 

        %derivatives of qe
        Dqex_phi = sparse(1:length(ein),1:length(ein),-c22*permex.*((emean(ein,:)*Psi).^2+Psi_reg^2).^((beta_e-1)/2),length(ein),length(ein))*eddx(ein,:) ...
                    + sparse(1:length(ein),1:length(ein),-c22*permex.*((emean(ein,:)*Psi).^2+Psi_reg^2).^((beta_e-3)/2).*(emean(ein,:)*Psi).*Psi_x(ein)*(beta_e-1),length(ein),length(ein))*emean(ein,:)*DPsi_phi ...
                    + sparse(1:length(ein),1:length(ein),-c22*((emean(ein,:)*Psi).^2+Psi_reg^2).^((beta_e-1)/2).*Psi_x(ein) ,length(ein),length(ein))*Dpermex_he(:,:)*Dhe_phi(:,:);
        Dqey_phi = sparse(1:length(fin),1:length(fin),-c22*permey.*((fmean(fin,:)*Psi).^2+Psi_reg^2).^((beta_e-1)/2),length(fin),length(fin))*fddy(fin,:) ...
                    + sparse(1:length(fin),1:length(fin),-c22*permey.*((fmean(fin,:)*Psi).^2+Psi_reg^2).^((beta_e-3)/2).*(fmean(fin,:)*Psi).*Psi_y(fin)*(beta_e-1),length(fin),length(fin))*fmean(fin,:)*DPsi_phi ...
                    + sparse(1:length(fin),1:length(fin),-c22*((fmean(fin,:)*Psi).^2+Psi_reg^2).^((beta_e-1)/2).*Psi_y(fin) ,length(fin),length(fin))*Dpermey_he(:,:)*Dhe_phi(:,:);

        %derivatives of Xi
        % [ using Xi calculated on nodes ]
        DXi_phi = - nmeanx(:,ein)*( sparse(1:length(ein),1:length(ein),qsx(ein)+qex(ein),length(ein),length(ein))*eddx(ein,:) + sparse(1:length(ein),1:length(ein),Psi_x(ein),length(ein),length(ein))*(Dqsx_phi+Dqex_phi) ) ...
                - nmeany(:,fin)*( sparse(1:length(fin),1:length(fin),qsy(fin)+qey(fin),length(fin),length(fin))*fddy(fin,:) + sparse(1:length(fin),1:length(fin),Psi_y(fin),length(fin),length(fin))*(Dqsy_phi+Dqey_phi) );
        DXi_hs = - nmeanx(:,ein)*( sparse(1:length(ein),1:length(ein),Psi_x(ein),length(ein),length(ein))*(Dqsx_hs) ) ...
                - nmeany(:,fin)*( sparse(1:length(fin),1:length(fin),Psi_y(fin),length(fin),length(fin))*(Dqsy_hs) );
        DXix_phi = (1-opts.noXi)*sparse(1:length(ein),1:length(ein),c20*lcx(ein),length(ein),length(ein))*emean(ein,:)*DXi_phi;
        DXix_hs = (1-opts.noXi)*sparse(1:length(ein),1:length(ein),c20*lcx(ein),length(ein),length(ein))*emean(ein,:)*DXi_hs;
        DXiy_phi = (1-opts.noXi)*sparse(1:length(fin),1:length(fin),c20*lcy(fin),length(fin),length(fin))*fmean(fin,:)*DXi_phi;
        DXiy_hs = (1-opts.noXi)*sparse(1:length(fin),1:length(fin),c20*lcy(fin),length(fin),length(fin))*fmean(fin,:)*DXi_hs;

        if opts.include_diag
            %derivatives of corner Xi
            DXis_phi = (1-opts.noXi)*sparse(1:length(cin),1:length(cin),c20*lcs(cin),length(cin),length(cin))*cmean(cin,:)*DXi_phi;
            DXis_hs = (1-opts.noXi)*sparse(1:length(cin),1:length(cin),c20*lcs(cin),length(cin),length(cin))*cmean(cin,:)*DXi_hs;
            DXir_phi = (1-opts.noXi)*sparse(1:length(cin),1:length(cin),c20*lcr(cin),length(cin),length(cin))*cmean(cin,:)*DXi_phi;
            DXir_hs = (1-opts.noXi)*sparse(1:length(cin),1:length(cin),c20*lcr(cin),length(cin),length(cin))*cmean(cin,:)*DXi_hs;
        else
            DXis_phi = sparse(1,1,0,length(cin),nIJ); 
            DXis_hs = sparse(1,1,0,length(cin),nIJ); 
            DXir_phi = sparse(1,1,0,length(cin),nIJ); 
            DXir_hs = sparse(1,1,0,length(cin),nIJ); 
        end

        % derivatives of permeabilities for blisters
        alf = 3.0; bet = 1.0;
        if oo.mean_perms
            Dpermbx_hb = sparse(1:length(ein),econnect(ein,1),alf*( Dx(econnect(ein,1)).*max(hb(econnect(ein,1)),0).^(alf-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ) ...
                       + sparse(1:length(ein),econnect(ein,2),alf*( Dx(econnect(ein,2)).*max(hb(econnect(ein,2)),0).^(alf-1) ).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(-1),length(ein),nIJ);
            Dpermby_hb = sparse(1:length(fin),fconnect(fin,1),alf*( Dy(fconnect(fin,1)).*max(hb(fconnect(fin,1)),0).^(alf-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ) ...
                       + sparse(1:length(fin),fconnect(fin,2),alf*( Dy(fconnect(fin,2)).*max(hb(fconnect(fin,2)),0).^(alf-1) ).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(-1),length(fin),nIJ);
        else
            Dpermbx_hb = sparse(1:length(ein),econnect(ein,1),alf*Dx(econnect(ein,1)).*max(hb(econnect(ein,2)),0).^(alf/bet+alf).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(bet).*max(hb(econnect(ein,1)),0).^(alf-1).*( Dx(econnect(ein,1)).*max(hb(econnect(ein,2)),0).^(alf/bet) + Dx(econnect(ein,2)).*max(hb(econnect(ein,1)),0).^(alf/bet) + perm_reg).^(-bet-1),length(ein),nIJ) ...
                       + sparse(1:length(ein),econnect(ein,2),alf*Dx(econnect(ein,2)).*max(hb(econnect(ein,1)),0).^(alf/bet+alf).*(Dx(econnect(ein,1))+Dx(econnect(ein,2))).^(bet).*max(hb(econnect(ein,2)),0).^(alf-1).*( Dx(econnect(ein,1)).*max(hb(econnect(ein,2)),0).^(alf/bet) + Dx(econnect(ein,2)).*max(hb(econnect(ein,1)),0).^(alf/bet) + perm_reg).^(-bet-1),length(ein),nIJ); 
            Dpermby_hb = sparse(1:length(fin),fconnect(fin,1),alf*Dy(fconnect(fin,1)).*max(hb(fconnect(fin,2)),0).^(alf/bet+alf).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(bet).*max(hb(fconnect(fin,1)),0).^(alf-1).*( Dy(fconnect(fin,1)).*max(hb(fconnect(fin,2)),0).^(alf/bet) + Dy(fconnect(fin,2)).*max(hb(fconnect(fin,1)),0).^(alf/bet) + perm_reg).^(-bet-1),length(fin),nIJ) ...
                       + sparse(1:length(fin),fconnect(fin,2),alf*Dy(fconnect(fin,2)).*max(hb(fconnect(fin,1)),0).^(alf/bet+alf).*(Dy(fconnect(fin,1))+Dy(fconnect(fin,2))).^(bet).*max(hb(fconnect(fin,2)),0).^(alf-1).*( Dy(fconnect(fin,1)).*max(hb(fconnect(fin,2)),0).^(alf/bet) + Dy(fconnect(fin,2)).*max(hb(fconnect(fin,1)),0).^(alf/bet) + perm_reg).^(-bet-1),length(fin),nIJ); 
        end
        DPsibx_pb = eddx(ein,:);
        DPsiby_pb = fddy(fin,:);
        
        % derivatives of blister fluxes
        % qbx = - pp.c42*permbx.*Psib_x(ein) = - pp.c42*permbx.*eddx(ein,:)*pb;                                 %
        % qby = - pp.c42*permby.*Psib_y(fin) = - pp.c42*permby.*fddy(fin,:)*pb;                                 % 

        Dqbx_pb = -c42*permbx.*DPsibx_pb;                                                                       % ein * nIJ
        Dqby_pb = -c42*permby.*DPsiby_pb;                                                                       % fin * nIJ
        Dqbx_hb = sparse(1:length(ein),1:length(ein),-c42*Psib_x(ein),length(ein),length(ein))*Dpermbx_hb(:,:); % ein * nIJ
        Dqby_hb = sparse(1:length(fin),1:length(fin),-c42*Psib_y(fin),length(fin),length(fin))*Dpermby_hb(:,:); % fin * nIJ
        
        %% Derivatives of objective function
        %DF1 
        DF1_hs = sparse(1:length(ns),1:length(ns), -c15*ones(length(ns),1).*dt.^(-1) ...
                        - (ones(length(ns),1)>c17.*hs(ns)).*c16.*c17.*Ub(ns) ... %%
                        - c18.*abs(phi_0(ns)-phi(ns)).^(n_Glen-1).*(phi_0(ns)-phi(ns)) ,length(ns),length(ns)); 
        temp = sparse(nin,1:length(nin), c18*hs(nin).*abs(phi_0(nin)-phi(nin)).^(n_Glen-1)*(n_Glen) ,nIJ,length(nin)); % ns eqns, nin variables 
        DF1_phi = temp(ns,:); 

        %DF2
        temp = sparse(1:length(nin), nin, -c1*ones(length(nin),1).*dt.^(-1) ,length(nin), nIJ); % num of eqn < num of variables

        DF2_hs = + temp(:,ns) ...
                - c4*( nddx(nin,ein)*Dqsx_hs(:,ns) + nddy(nin,fin)*Dqsy_hs(:,ns) ) ...
                + sparse(1:length(nin),1:length(nin),Dy(nin).^(-1),length(nin),length(nin))*( c11*nmeanx(nin,ein)*DXix_hs(:,ns) ) ...
                + sparse(1:length(nin),1:length(nin),Dx(nin).^(-1),length(nin),length(nin))*( c11*nmeany(nin,fin)*DXiy_hs(:,ns) ) ...
                + sparse(1:length(nin),1:length(nin),Ds(nin).*Dx(nin).^(-1).*Dy(nin).^(-1),length(nin),length(nin))*( c11*nmeans(nin,cin)*DXis_hs(:,ns) ) ...
                + sparse(1:length(nin),1:length(nin),Dr(nin).*Dx(nin).^(-1).*Dy(nin).^(-1),length(nin),length(nin))*( c11*nmeanr(nin,cin)*DXir_hs(:,ns) );

        DF2_phi = - c2*Dhe_phi(nin,nin).*dt.^(-1) ...
                - c3*Dhv_phi(nin,nin).*dt.^(-1) ...
                - c3*Dhm_phi(nin,nin).*dt.^(-1) ...
                - c4*( nddx(nin,ein)*Dqsx_phi(:,nin) + nddy(nin,fin)*Dqsy_phi(:,nin) ) ...
                - c5*( nddx(nin,ein)*Dqex_phi(:,nin) + nddy(nin,fin)*Dqey_phi(:,nin) ) ...
                + sparse(1:length(nin),1:length(nin),Dy(nin).^(-1),length(nin),length(nin))*( -c9*nddx(nin,ein)*DQx_phi(:,nin) + c11*nmeanx(nin,ein)*(DXicx_phi(:,nin)+DXix_phi(:,nin)) )...
                + sparse(1:length(nin),1:length(nin),Dx(nin).^(-1),length(nin),length(nin))*( -c9*nddy(nin,fin)*DQy_phi(:,nin) + c11*nmeany(nin,fin)*(DXicy_phi(:,nin)+DXiy_phi(:,nin)) )...
                - c9*ndds(nin,cin)*DQs_phi(:,nin) ...
                - c9*nddr(nin,cin)*DQr_phi(:,nin) ...
                + sparse(1:length(nin),1:length(nin),Ds(nin).*Dx(nin).^(-1).*Dy(nin).^(-1),length(nin),length(nin))*( + c11*nmeans(nin,cin)*(DXics_phi(:,nin)+DXis_phi(:,nin)) ) ...
                + sparse(1:length(nin),1:length(nin),Dr(nin).*Dx(nin).^(-1).*Dy(nin).^(-1),length(nin),length(nin))*( + c11*nmeanr(nin,cin)*(DXicr_phi(:,nin)+DXir_phi(:,nin)) );
                % + sparse(1:length(nin),1:length(nin),-c44*k_b(hs(nin),Smean(nin),pp,oo).*Vb_reg(Vb(nin),pp,oo).*Rb(nin).*Dx(nin).^(-1).*Dy(nin).^(-1),length(nin),length(nin)); % derivative to phi

        DF2_pb = sparse(1:length(nin), 1:length(nin), zeros(length(nin),1), length(nin), length(nin));
        temp   = sparse(1:length(nin), nin, c50*ones(length(nin),1), length(nin), nIJ);
        DF2_hb = temp(:,ns);

        DF2_Sx = sparse(1:length(nin),1:length(nin), Dy(nin).^(-1),length(nin),length(nin))*( -c8*dt^(-1).*nmeanx(nin,ein) ...
                    - c9*nddx(nin,ein)*DQx_Sx(:,ein) + c11*nmeanx(nin,ein)*DXicx_Sx(:,ein) );  

        DF2_Sy = sparse(1:length(nin),1:length(nin),Dx(nin).^(-1),length(nin),length(nin))*( -c8*dt^(-1).*nmeany(nin,fin) ...
                    - c9*nddy(nin,fin)*DQy_Sy(:,fin) + c11*nmeany(nin,fin)*DXicy_Sy(:,fin) ); 
        DF2_Ss = - c9*ndds(nin,cin)*DQs_Ss(:,cin) ...
                    + sparse(1:length(nin),1:length(nin),Ds(nin).*Dx(nin).^(-1).*Dy(nin).^(-1),length(nin),length(nin))*( -c8*dt^(-1).*nmeans(nin,cin) ...
                    + c11*nmeans(nin,cin)*DXics_Ss(:,cin) );
        DF2_Sr = - c9*nddr(nin,cin)*DQr_Sr(:,cin) ...
                    + sparse(1:length(nin),1:length(nin),Dr(nin).*Dx(nin).^(-1).*Dy(nin).^(-1),length(nin),length(nin))*( -c8*dt^(-1).*nmeanr(nin,cin) ...
                    + c11*nmeanr(nin,cin)*DXicr_Sr(:,cin) );
        %DF3            
        DF3_Sx = sparse(1:length(ein),1:length(ein),-c12*ones(length(ein),1).*dt.^(-1) ...
            - c14.*abs(emean(ein,:)*(phi_0-phi)).^(n_Glen-1).*(emean(ein,:)*(phi_0-phi)) ...
            - c36*c35*aa.lcx(ein).*(emean(ein,:)*Ub).*(1>c36*Sx(ein)) ,length(ein),length(ein)) ...
            + c13*DXicx_Sx(:,ein);

        DF3_phi = sparse(1:length(ein),1:length(ein), c14.*Sx(ein).*abs(emean(ein,:)*(phi_0-phi)).^(n_Glen-1)*n_Glen ,length(ein),length(ein))*emean(ein,nin) ...
                + c13*(DXicx_phi(:,nin)+DXix_phi(:,nin));

        DF3_hs = c13*DXix_hs(:,ns);

        %DF4
        DF4_Sy = sparse(1:length(fin),1:length(fin),-c12*ones(length(fin),1).*dt.^(-1) ...
            - c14.*abs(fmean(fin,:)*(phi_0-phi)).^(n_Glen-1).*(fmean(fin,:)*(phi_0-phi)) ...
            - c36*c35*aa.lcy(fin).*(fmean(fin,:)*Ub).*(1>c36*Sy(fin)) ,length(fin),length(fin)) ...
            + c13*DXicy_Sy(:,fin);

        DF4_phi = sparse(1:length(fin),1:length(fin), c14.*Sy(fin).*abs(fmean(fin,:)*(phi_0-phi)).^(n_Glen-1)*n_Glen ,length(fin),length(fin))*fmean(fin,nin) ...
                + c13*(DXicy_phi(:,nin)+DXiy_phi(:,nin));

        DF4_hs = c13*DXiy_hs(:,ns);

        %DF5   
        DF5_Ss = sparse(1:length(cin),1:length(cin),-c12*ones(length(cin),1).*dt.^(-1) ...
            - c14.*abs(cmean(cin,:)*(phi_0-phi)).^(n_Glen-1).*(cmean(cin,:)*(phi_0-phi)) ...
            - c36*c35*aa.lcs(cin).*(cmean(cin,:)*Ub).*(1>c36*Ss(cin)) ,length(cin),length(cin)) ...
            + c13*DXics_Ss(:,cin);

        DF5_phi = sparse(1:length(cin),1:length(cin), c14.*Ss(cin).*abs(cmean(cin,:)*(phi_0-phi)).^(n_Glen-1)*n_Glen ,length(cin),length(cin))*cmean(cin,nin) ...
                + c13*(DXics_phi(:,nin)+DXis_phi(:,nin));

        DF5_hs = c13*DXis_hs(:,ns);

        %DF6           
        DF6_Sr = sparse(1:length(cin),1:length(cin),-c12*ones(length(cin),1).*dt.^(-1) ...
            - c14.*abs(cmean(cin,:)*(phi_0-phi)).^(n_Glen-1).*(cmean(cin,:)*(phi_0-phi)) ...
            - c36*c35*aa.lcr(cin).*(cmean(cin,:)*Ub).*(1>c36*Sr(cin)) ,length(cin),length(cin)) ...
            + c13*DXicr_Sr(:,cin);

        DF6_phi = sparse(1:length(cin),1:length(cin), c14.*Sr(cin).*abs(cmean(cin,:)*(phi_0-phi)).^(n_Glen-1)*n_Glen ,length(cin),length(cin))*cmean(cin,nin) ...
                + c13*(DXicr_phi(:,nin)+DXir_phi(:,nin));

        DF6_hs = c13*DXir_hs(:,ns);

        % DF7
        biharmonic = (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*(gg.nddx*gg.eddx + gg.nddy*gg.fddy);
        DF7_hb = c49 * biharmonic(ns,ns);
        temp = sparse(nin, 1:length(nin), -ones(length(nin),1), nIJ, length(nin)); % ns eqns, nin variables 
        DF7_pb = temp(ns,:);

        % sparse(nin, 1:length(nin), -ones(length(nin),1), nIJ, length(nin)); % ns eqns, nin variables 
        % DF8
        DF8_pb = - c45*(gg.nddx(nin,ein)*Dqbx_pb(:,nin) + gg.nddy(nin,fin)*Dqby_pb(:,nin));

        temp = sparse(1:length(nin), nin, - ones(length(nin),1).*dt.^(-1) - c44*ones(length(nin),1), length(nin), nIJ);
        DF8_hb = temp(:,ns) - c45*(gg.nddx(nin,ein)*Dqbx_hb(:,ns) + gg.nddy(nin,fin)*Dqby_hb(:,ns));

        %% construct Jacobian matrix
        % % [ original way of calculating ]
        % if opts.no_channels && opts.no_sheet  
        %     J =   DF2_phi ;
        % elseif opts.no_channels  
        %     J = [  DF1_hs DF1_phi; ...
        %            DF2_hs DF2_phi ];
        % elseif opts.no_sheet
        %     if opts.include_diag
        %     J = [  DF2_phi DF2_Sx DF2_Sy DF2_Ss DF2_Sr; ...
        %            DF3_phi DF3_Sx sparse(length(ein),length(fin)) sparse(length(ein),length(cin)) sparse(length(ein),length(cin)); ...
        %            DF4_phi sparse(length(fin),length(ein)) DF4_Sy sparse(length(fin),length(cin)) sparse(length(fin),length(cin)); ...
        %            DF5_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) DF5_Ss sparse(length(cin),length(cin)); ...
        %            DF6_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) sparse(length(cin),length(cin)) DF6_Sr ];
        %     else
        %     J = [  DF2_phi DF2_Sx DF2_Sy; ...
        %            DF3_phi DF3_Sx sparse(length(ein),length(fin)); ...
        %            DF4_phi sparse(length(fin),length(ein)) DF4_Sy ];
        %     end
        % else  
        %     if opts.include_diag
        %     J = [  DF1_hs DF1_phi sparse(length(ns),length(ein)) sparse(length(ns),length(fin)) sparse(length(ns),length(cin)) sparse(length(ns),length(cin)); ...
        %            DF2_hs DF2_phi DF2_Sx DF2_Sy DF2_Ss DF2_Sr; ...
        %            DF3_hs DF3_phi DF3_Sx sparse(length(ein),length(fin)) sparse(length(ein),length(cin)) sparse(length(ein),length(cin)); ...
        %            DF4_hs DF4_phi sparse(length(fin),length(ein)) DF4_Sy sparse(length(fin),length(cin)) sparse(length(fin),length(cin)); ...
        %            DF5_hs DF5_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) DF5_Ss sparse(length(cin),length(cin)); ...
        %            DF6_hs DF6_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) sparse(length(cin),length(cin)) DF6_Sr ];
        %     else
        %     J = [  DF1_hs DF1_phi sparse(length(ns),length(ein)) sparse(length(ns),length(fin)); ...
        %            DF2_hs DF2_phi DF2_Sx DF2_Sy; ...
        %            DF3_hs DF3_phi DF3_Sx sparse(length(ein),length(fin)); ...
        %            DF4_hs DF4_phi sparse(length(fin),length(ein)) DF4_Sy ];
        %     end
        % end
    
    % J = [   DF1_hs DF1_phi sparse(length(ns),length(ein)) sparse(length(ns),length(fin)) sparse(length(ns),length(cin)) sparse(length(ns),length(cin)) sparse(length(ns),length(nin)) sparse(length(ns),length(nin)); ...
    %         DF2_hs DF2_phi DF2_Sx DF2_Sy DF2_Ss DF2_Sr sparse(length(nin),length(nin)) sparse(length(nin),length(nin)); ...
    %         DF3_hs DF3_phi DF3_Sx sparse(length(ein),length(fin)) sparse(length(ein),length(cin)) sparse(length(ein),length(cin)) sparse(length(ein),length(nin)) sparse(length(ein),length(nin)); ...
    %         DF4_hs DF4_phi sparse(length(fin),length(ein)) DF4_Sy sparse(length(fin),length(cin)) sparse(length(fin),length(cin)) sparse(length(fin),length(nin)) sparse(length(fin),length(nin)); ...
    %         DF5_hs DF5_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) DF5_Ss sparse(length(cin),length(cin)) sparse(length(cin),length(nin)) sparse(length(cin),length(nin)); ...
    %         DF6_hs DF6_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) sparse(length(cin),length(cin)) DF6_Sr sparse(length(cin),length(nin)) sparse(length(cin),length(nin)); ...
    %         sparse(length(nin),length(ns)) DF7_phi sparse(length(nin),length(ein)) sparse(length(nin),length(fin)) sparse(length(nin),length(cin)) sparse(length(nin),length(cin)) DF7_Vb DF7_Rb; ...
    %         sparse(length(nin),length(ns)) sparse(length(nin),length(nin)) sparse(length(nin),length(ein)) sparse(length(nin),length(fin)) sparse(length(nin),length(cin)) sparse(length(nin),length(cin)) DF8_Vb DF8_Rb ];
    
        %[ alternative that looks nicer : not sure if it will be quicker or slower ]
        J = [   DF1_hs DF1_phi sparse(length(ns),length(ein)) sparse(length(ns),length(fin)) sparse(length(ns),length(cin)) sparse(length(ns),length(cin)) sparse(length(ns),length(ns)) sparse(length(ns),length(nin)); ...
                DF2_hs DF2_phi DF2_Sx DF2_Sy DF2_Ss DF2_Sr DF2_hb DF2_pb; ...
                DF3_hs DF3_phi DF3_Sx sparse(length(ein),length(fin)) sparse(length(ein),length(cin)) sparse(length(ein),length(cin)) sparse(length(ein),length(ns)) sparse(length(ein),length(nin)); ...
                DF4_hs DF4_phi sparse(length(fin),length(ein)) DF4_Sy sparse(length(fin),length(cin)) sparse(length(fin),length(cin)) sparse(length(fin),length(ns)) sparse(length(fin),length(nin)); ...
                DF5_hs DF5_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) DF5_Ss sparse(length(cin),length(cin)) sparse(length(cin),length(ns)) sparse(length(cin),length(nin)); ...
                DF6_hs DF6_phi sparse(length(cin),length(ein)) sparse(length(cin),length(fin)) sparse(length(cin),length(cin)) DF6_Sr sparse(length(cin),length(ns)) sparse(length(cin),length(nin)); ...
                sparse(length(ns),length(ns)) sparse(length(ns),length(nin)) sparse(length(ns),length(ein)) sparse(length(ns),length(fin)) sparse(length(ns),length(cin)) sparse(length(ns),length(cin)) DF7_hb DF7_pb; ...
                sparse(length(nin),length(ns)) sparse(length(nin),length(nin)) sparse(length(nin),length(ein)) sparse(length(nin),length(fin)) sparse(length(nin),length(cin)) sparse(length(nin),length(cin)) DF8_hb DF8_pb ];
        ii = [];
        if ~opts.no_sheet, ii = [ ii 1:length(ns) ]; end % include sheet equation
        ii = [ ii length(ns)+(1:length(nin)) ]; % include the mass conservation equation
        if ~opts.no_channels, ii = [ ii length(ns)+length(nin)+(1:length(ein)) length(ns)+length(nin)+length(ein)+(1:length(fin)) ]; end % include horizontal and vertical channels
        if opts.include_diag, ii = [ ii length(ns)+length(nin)+length(ein)+length(fin)+(1:length(cin)) length(ns)+length(nin)+length(ein)+length(fin)+length(cin)+(1:length(cin)) ]; end % include diagnonal channels
        % include blister equations for Vb and Rb
        if opts.include_blister, ii = [ii length(ns)+length(nin)+length(ein)+length(fin)+length(cin)+length(cin)+(1:length(ns))]; end
        if opts.include_pressure, ii = [ii length(ns)+length(nin)+length(ein)+length(fin)+length(cin)+length(cin)+length(ns)+(1:length(nin))]; end

        J = J(ii,ii);
    end
end
            
%% auxilliary functions [ taken from hydro_timestep_diag 13 August 2014 ]
function out = hm_fun(p_w,p_i,sigma_m,pp,~)
% moulin storage as function of pressure p_w
    sigma_m_reg = pp.c33;
    p_w_reg = pp.c34; 
    r = pp.c27;
    out = sigma_m.*(p_w + sigma_m_reg*exp(-(r*p_i-p_w)/p_w_reg));
    if any(p_w>r*p_i) && sigma_m_reg>0, disp('   Using overflow regulation'); end
end
function out = Dhm_fun(p_w,p_i,sigma_m,pp,~)
% derivative of moulin storage as function of pressure p_w
    sigma_m_reg = pp.c33;
    p_w_reg = pp.c34; 
    r = pp.c27;
    out = sigma_m.*(p_w.^0 + sigma_m_reg*exp(-(r*p_i-p_w)/p_w_reg)/p_w_reg);
end
function out = hv_fun(p_w,p_i,sigma,pp,~)
% storage as function of pressure p_w
    sigma_log = pp.c39;
    N0 = pp.c40;
    N_reg1 = pp.c30;
    
    N = p_i-p_w;
    Sigma_log = - sigma_log*log(min((N+N_reg1)./N0,1));
    ii_neg = N<0; 
    Sigma_log(ii_neg) = -sigma_log*log(min((N_reg1)./N0,1)) -sigma_log*N(ii_neg)/N_reg1;  %has continuous derivative

    out = sigma.*p_w + Sigma_log;
end
function out = Dhv_fun(p_w,p_i,sigma,pp,~)
% derivative of storage as function of pressure p_w
    sigma_log = pp.c39;
    N0 = pp.c40;
    N_reg1 = pp.c30;
    
    N = p_i-p_w;
    DSigma_log = -sigma_log*(N+N_reg1).^(-1).*(N+N_reg1<=N0);
    ii_neg = N<0; DSigma_log(ii_neg) = -sigma_log/N_reg1; 
    
    out = sigma.*p_w.^0 + DSigma_log;
end
function out = he_fun(N,p_i,pp,~)
% elastic sheet depth as function of effective pressure N
    gamma_e = pp.gamma_e;
    c_e_power = pp.c24;
    N_reg = pp.N_reg;
    c_e_log = pp.c28;
    c_e_reg1 = pp.c29;
    N_reg1 = pp.c30;
    c_e_reg2 = pp.c31;
    N_reg2 = pp.c32;
    c_e_reg3 = pp.c37;
    N_reg3 = pp.c38;
    
    he_power = c_e_power*(max(p_i-N,0).*(p_i.^2+N_reg^2).^(-1/2)).^(gamma_e);
    he_log = -c_e_log*log(min((N.^2+N_reg^2).^(1/2).*(p_i.^2+N_reg^2).^(-1/2),1));  %[ effective pressure must be positive ]

    he_reg = - c_e_reg1*log(min((N+N_reg1).*(p_i+N_reg).^(-1),1));
    ii_neg = N<0; 
    he_reg(ii_neg) = -c_e_reg1*log(min((N_reg1).*(p_i(ii_neg)+N_reg).^(-1),1)) -c_e_reg1*N(ii_neg)/N_reg1;  %has continuous derivative

    he_reg2 = -c_e_reg2*min(N,0) + 1/2*N_reg2*c_e_reg2*min(max(1-N/N_reg2,0),1).^2;

    he_reg3 = -c_e_reg3*log(min((N+N_reg1)./N_reg3,1));
    ii_neg = N<0; 
    he_reg3(ii_neg) = -c_e_reg3*log(min((N_reg1)./N_reg3,1)) -c_e_reg3*N(ii_neg)/N_reg1;  %has continuous derivative

    out =  he_power + he_log + he_reg + he_reg2 + he_reg3;
end
function out = Dhe_fun(N,p_i,pp,~)
% derivative of elastic sheet depth as function of effective pressure N
    gamma_e = pp.gamma_e;
    c_e_power = pp.c24;
    N_reg = pp.N_reg;
    c_e_log = pp.c28;
    c_e_reg = pp.c29;
    N_reg1 = pp.c30;
    c_e_reg2 = pp.c31;
    N_reg2 = pp.c32;
    c_e_reg3 = pp.c37;
    N_reg3 = pp.c38;

    Dhe_power = -gamma_e*c_e_power*(max(p_i-N,0).*(p_i.^2+N_reg^2).^(-1/2)).^(gamma_e-1).*(p_i.^2+N_reg^2).^(-1/2).*(N<=p_i);
    Dhe_log = -c_e_log.*N.*(N.^2+N_reg^2).^(-1).*(N<=p_i);  %[ effective pressure must be positive ]

    Dhe_reg = -c_e_reg.*(N+N_reg1).^(-1).*(N+N_reg1<=p_i+N_reg);
    ii_neg = N<0; 
    Dhe_reg(ii_neg) = -c_e_reg/N_reg1;

    Dhe_reg2 = -c_e_reg2*min(max(1-N/N_reg2,0),1);

    Dhe_reg3 = -c_e_reg3*(N+N_reg1).^(-1).*(N+N_reg1<=N_reg3);
    ii_neg = N<0; 
    Dhe_reg3(ii_neg) = -c_e_reg3/N_reg1; 

    out = Dhe_power + Dhe_log + Dhe_reg + Dhe_reg2 + Dhe_reg3;
end