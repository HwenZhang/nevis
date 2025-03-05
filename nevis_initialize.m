function [aa,vv] = nevis_initialize(b,s,gg,pp,oo)
    % Assign default prescribed fields and boundary conditions and initial
    % conditions [ probably need to change initial hs and phi ]
    % Inputs
    %   b bed elevation on nodes
    %   s surface elevation on nodes
    %   gg grid structure
    %   pp parameters
    %   oo options [optional]
    % Outputs
    %   aa structure containing prescribed fields and boundary conditions
    %   vv structure containing initial variables
    %
    % 21 August 2014
    if nargin<5, oo = struct; end
    if ~isfield(pp,'u_b'), pp.u_b = 1; end
    if ~isfield(pp,'l_c'), pp.l_c = 1; end
    if ~isfield(pp,'phi_s'), pp.phi_s = -inf; end
    
    %% reshape inputs
    b = reshape(b,gg.nIJ,1);             % bed profile
    s = reshape(s,gg.nIJ,1);             % surface profile
    
    %% hydraulic potential
    H = max(s-b,0);                      % ice thickness
    phi_a = pp.r*b;                      % atmospheric potential
    phi_0 = phi_a+(s-b);                 % hydraulic potential
    
    %% channel roughness sizes
    lcx = pp.l_c*ones(gg.eIJ,1);         % effective channel width
    lcy = pp.l_c*ones(gg.fIJ,1);
    lcs = pp.l_c*ones(gg.cIJ,1);
    lcr = pp.l_c*ones(gg.cIJ,1);
    
    %% storage capacity [ value taken from pp ]
    sigma = pp.sigma*ones(gg.nIJ,1);
    sigma(gg.nout) = NaN;
    
    %% basal velocity [ value taken from pp ]
    Ub = pp.u_b*ones(gg.nIJ,1); 
    Ub(gg.nout) = 0;
    
    %% basal melting rate [ value taken from pp ]
    m = pp.melt*ones(gg.nIJ,1); 
    m(gg.nout) = 0;
    
    %% input to the moulin [supposed to be the outflow from blister instead]
    E = 0*ones(gg.nIJ,1); 
    E(gg.nout) = 0;

    %% input to and out from the blister
    Qb_in = 0*ones(gg.nIJ,1);
    
    %% variables
    phi = phi_a+0*(phi_0-phi_a); % potential
    hs = 0*ones(gg.nIJ,1);       % cavity thickness
    Sx = 0*ones(gg.eIJ,1);       % cross sectional area of x channels
    Sy = 0*ones(gg.fIJ,1);       % cross sectional area of y channels
    Ss = 0*ones(gg.cIJ,1);       % cross sectional area of diagonal s channels
    Sr = 0*ones(gg.cIJ,1);       % cross sectional area of diagonal r channels
    hb = 0*ones(gg.nIJ,1);       % blister thickness
    pb = 0*ones(gg.nIJ,1);       % blister pressure
    
    %% boundary conditions
    aa.phi = max(phi_a(gg.nbdy),pp.phi_s); % boundary values
    aa.phi_b = max(phi_a,pp.phi_s);        % bed potential
    
    %% put fields and variables in structs
    % precribed field with bound conditions
    aa.phi_0 = phi_0; 
    aa.phi_a = phi_a; 
    aa.lcx = lcx; 
    aa.lcy = lcy; 
    aa.lcs = lcs; 
    aa.lcr = lcr; 
    aa.b = b; 
    aa.s = s; 
    aa.H = H; 
    aa.sigma = sigma; 
    aa.m = m; 
    aa.Ub = Ub; 
    aa.E = E;
    aa.Qb_in = Qb_in;
    
    % solution fields
    vv.phi = phi; 
    vv.hs = hs; 
    vv.Sx = Sx; 
    vv.Sy = Sy; 
    vv.Ss = Ss; 
    vv.Sr = Sr;
    vv.hb = hb;
    vv.pb = pb;

end