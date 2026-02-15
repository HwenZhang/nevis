function gg = nevis_label_ice(gg,oo)
    % update grid with boundary labels for solving for ice velocity
    % 
    % 10 Oct 2025 - edited from ssa/mask_and_label.m and nevis_label, to create grid structure
    %   with boundary node/edge labels for ice velocity solve, which may be
    %   different from those needed by the hydrology model

    if nargin<2, oo = struct; end
    if ~isfield(oo,'boundary_method'), oo.boundary_method = 'default'; end
    % mask_boundary_method: how to treat the irregular mask boundary
    %   'velocity' : prescribe velocity (e.g. u=0) on mask boundary [default]
    %   'stress_free' : prescribe stress-free on mask boundary
    if ~isfield(oo,'mask_boundary_method'), oo.mask_boundary_method = 'velocity'; end

    %% grid point labels [ labelled first along x axis then along y axis, so location i,j becomes i+(j-1)*I ]
    ns = (1:gg.nIJ)';
    es = (1:gg.eIJ)';
    fs = (1:gg.fIJ)';
    cs = (1:gg.cIJ)';

    %% labels of boundary nodes, edges, corners (from nevis_mask)
    nout = gg.nout;
    n1 = gg.n1; % nodes connected by an edge to an outside node : boundary nodes
    e0 = gg.e0; % edges connected to all outside nodes
    e1 = gg.e1; % edges connected to one outside node : boundary edges
    e2 = gg.e2; % edges connected to two edge nodes : along-boundary edges
    f0 = gg.f0; % edges connected to all outside nodes
    f1 = gg.f1; % edges connected to one outside node : boundary edges
    f2 = gg.f2; % edges connected to two edge nodes : along-boundary edges
    c0 = gg.c0; % corners connected to all outside nodes
    c1 = gg.c1; % corners connected to one outside node 
    c2 = gg.c2; % corners connected to two outside nodes
    c3 = gg.c3; % corners connected to three outside nodes

    nx = gg.nx; ny = gg.ny;
    ex = gg.ex; ey = gg.ey;
    fx = gg.fx; fy = gg.fy;
    cx = gg.cx; cy = gg.cy;

    temp = zeros(gg.eIJ,1); 
    temp(e1) = 1; 
    n1x = find(gg.nmeanx*temp>0); 
    n1x = setdiff(n1x,nout); % nodes connected by an e edge to an outside node : x boundary nodes 
    temp = zeros(gg.fIJ,1); 
    temp(f1) = 1; 
    n1y = find(gg.nmeany*temp>0); 
    n1y = setdiff(n1y,nout); % nodes connected by an f edge to an outside node : y boundary nodes 

    %% default boundary labels for whole grid
    nlbdy = gg.bdy.nlbdy; nrbdy = gg.bdy.nrbdy; ntopbdy = gg.bdy.ntopbdy; nbotbdy = gg.bdy.nbotbdy;
    elbdy = gg.bdy.elbdy; erbdy = gg.bdy.erbdy; etopbdy = gg.bdy.etopbdy; ebotbdy = gg.bdy.ebotbdy;
    flbdy = gg.bdy.flbdy; frbdy = gg.bdy.frbdy; ftopbdy = gg.bdy.ftopbdy; fbotbdy = gg.bdy.fbotbdy;
    clbdy = gg.bdy.clbdy; crbdy = gg.bdy.crbdy; ctopbdy = gg.bdy.ctopbdy; cbotbdy = gg.bdy.cbotbdy;
    
    if isfield(oo,'boundary_method') && strcmp(oo.boundary_method,'per_tblr')   
        % [ appropriate for periodic top and bottom boundaries, periodic left and right ]
        ebdy = unique([]);       
        fbdy = unique([]);         
        nbdyx = unique([]); 
        nbdyy = unique([]); 
        cbdy = unique([]); 
        eout = unique([]);
        fout = unique([]);
        cout = unique([]);

    elseif strcmp(oo.boundary_method,'vel_tblr')
        % [ appropriate for prescribing velocity on all boundaries ]
        ebdy = unique([etopbdy; ebotbdy; elbdy; erbdy]);       
        fbdy = unique([ftopbdy; fbotbdy; flbdy; frbdy]);         
        nbdyx = unique([]); 
        nbdyy = unique([]); 
        cbdy = unique([]); 
        eout = unique([]);
        fout = unique([]);
        cout = unique([ctopbdy; cbotbdy; clbdy; crbdy]);

    elseif strcmp(oo.boundary_method,'stress_tblr')
        % [ appropriate for prescribing normal and tangential stress on all boundaries]
        ebdy = unique([]);
        fbdy = unique([]);
        nbdyx = unique([nlbdy; nrbdy]); 
        nbdyy = unique([nbotbdy; ntopbdy]); 
        cbdy = unique([cbotbdy; ctopbdy; clbdy; crbdy]);
        eout = unique([ebotbdy; etopbdy; elbdy; erbdy]);
        fout = unique([fbotbdy; ftopbdy; flbdy; frbdy]);
        cout = unique([]);

    elseif strcmp(oo.boundary_method,'vel_r_stress_tbl')
        % [ appropriate for prescribing velocity on right, stress-free on others ]
        ebdy = unique([erbdy]);
        fbdy = unique([frbdy]);
        nbdyx = unique([nlbdy]); 
        nbdyy = unique([ntopbdy; nbotbdy]); 
        cbdy = unique([clbdy]);
        eout = unique([elbdy]);
        fout = unique([flbdy]);
        cout = unique([crbdy]);

    elseif strcmp(oo.boundary_method, 'stress_l_vel_tbl')
        % [ appropriate for prescribing zero stress on left, velocity on top, bottom, and right ] 
        ebdy = unique([etopbdy; ebotbdy; erbdy]); 
        fbdy = unique([ftopbdy; fbotbdy; frbdy]); 
        nbdyx = unique([nlbdy]); 
        nbdyy = unique([]); 
        cbdy = unique([clbdy]); 
        eout = unique([elbdy]); 
        fout = unique([flbdy]); 
        cout = unique([crbdy]);

    else
        % [ appropriate for prescribing velocity on left boundary, reflective or periodic on top and bottom, and normal and tangential stress on right]
        ebdy = unique([elbdy]);
        fbdy = unique([flbdy]);
        nbdyx = unique([nrbdy]); 
        nbdyy = unique([]); 
        cbdy = unique([crbdy]);
        eout = unique([erbdy]);
        fout = unique([]);
        cout = unique([clbdy]);

    end

    %% modification to boundary labels for masked (irregular) region
    if strcmp(oo.mask_boundary_method, 'stress_free')
        % ============================================================
        % Stress-free on the irregular mask boundary
        % ============================================================
        % The mask boundary nodes (n1x, n1y) are where the domain meets
        % the masked-out region. For stress-free conditions, we prescribe
        % sigma.n = 0 at these nodes instead of prescribing velocity.
        %
        % n1x: nodes adjacent to a masked x-edge (e1) -> stress in x
        % n1y: nodes adjacent to a masked y-edge (f1) -> stress in y
        %
        % The boundary edges (e1, f1) are placed in eout/fout (not ebdy/fbdy)
        % because we do NOT prescribe velocity there; instead the stress
        % equations at nbdyx/nbdyy enforce the stress-free condition.
        
        nbdyx = unique([nbdyx; n1x]);   % stress-free in x on mask boundary
        nbdyy = unique([nbdyy; n1y]);   % stress-free in y on mask boundary
        
        % mask boundary edges -> out (not velocity-prescribed)
        eout = unique([eout; e0; e1; e2]);
        fout = unique([fout; f0; f1; f2]);
        
        % along-boundary edges: prescribe velocity = 0 (tangential)
        % ebdy = unique([ebdy; e2]); 
        % fbdy = unique([fbdy; f2]); 
        
        cout = unique([cout; c0; c3; c2]);
        
        disp('nevis_label_ice: Stress-free BCs on irregular mask boundary');
        
    else
        % ============================================================
        % Velocity (Dirichlet) on the irregular mask boundary [default]
        % ============================================================
        % Prescribe velocity (typically u=0) on all mask boundary edges.
        % This is the original behaviour.
        
        ebdy = unique([ebdy; e1; e2]); 
        fbdy = unique([fbdy; f1; f2]); 
        
        eout = unique([eout; e0]);
        fout = unique([fout; f0]);
        cout = unique([cout; c0; c3; c2]);
        
        disp('nevis_label_ice: Velocity BCs on irregular mask boundary');
    end

    %% restrict labels to active area (i.e. remove points now assigned as out)
    ns = setdiff(ns,nout);   
    es = setdiff(es,eout);   
    fs = setdiff(fs,fout);   
    cs = setdiff(cs,cout);  

    nbdyx = setdiff(nbdyx,nout); 
    nbdyy = setdiff(nbdyy,nout); 
    nbdy = unique([nbdyx; nbdyy]);   
    ebdy = setdiff(ebdy,eout);   
    fbdy = setdiff(fbdy,fout);   
    cbdy = setdiff(cbdy,cout); 

    nin = setdiff(ns,nbdy);   
    ein = setdiff(es,ebdy);  
    fin = setdiff(fs,fbdy);  
    cin = setdiff(cs,cbdy); 
    
    %% diagnostic plot
    figure; clf;
    hold on;
    % external
    plot(nx(nout),ny(nout),'.','color',.9*[1 1 1],'markersize',8);
    plot(ex(eout),ey(eout),'+','color',.9*[1 1 1]);
    plot(fx(fout),fy(fout),'+','color',.9*[1 1 1]);
    plot(cx(cout),cy(cout),'x','color',.9*[1 1 1]);
    % boundary: velocity (Dirichlet)
    plot(ex(ebdy),ey(ebdy),'r+','markersize',6);
    plot(fx(fbdy),fy(fbdy),'b+','markersize',6);
    % boundary: stress (Neumann)
    plot(nx(nbdyx),ny(nbdyx),'bd','markersize',6);
    plot(nx(nbdyy),ny(nbdyy),'gd','markersize',6);
    % interior
    plot(nx(nin),ny(nin),'k.','markersize',4);
    plot(ex(ein),ey(ein),'k+','markersize',3);
    plot(fx(fin),fy(fin),'k+','markersize',3);
    legend('nout','eout','fout','cout','ebdy (vel)','fbdy (vel)',...
           'nbdyx (stress)','nbdyy (stress)','nin','ein','fin',...
           'Location','bestoutside');
    title('Ice velocity boundary labels');
    axis equal; box on;
    hold off;

    % [ possibly want to redefine mean operators here ? ]

    %% assign labels to gg (ice labels with suffix 2 may be slightly different to hydrology labels)
    gg.ns2 = ns;
    gg.es2 = es;
    gg.fs2 = fs;
    gg.cs2 = cs;
    gg.nbdyx2 = nbdyx;
    gg.nbdyy2 = nbdyy;
    gg.ebdy2 = ebdy;
    gg.cbdy2 = cbdy;
    gg.fbdy2 = fbdy;
    gg.nin2 = nin;
    gg.ein2 = ein;
    gg.fin2 = fin;
    gg.cin2 = cin;
    gg.nout2 = nout;
    gg.eout2 = eout;
    gg.fout2 = fout;
    gg.cout2 = cout;
    gg.n1x2 = n1x;
    gg.n1y2 = n1y;

    % redfine mean operators to only include active nodes
    temp = gg.nmeanx(:,es)*ones(length(es),1); 
    temp(temp==0) = inf;
    gg.nmeanx2 = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeanx;

    temp = gg.nmeany(:,fs)*ones(length(fs),1); 
    temp(temp==0) = inf; 
    gg.nmeany2 = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeany;
    temp = gg.nmeanc(:,cs)*ones(length(cs),1); 
    temp(temp==0) = inf; 
    gg.nmeanc2 = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeanc;

    temp = gg.emean(:,ns)*ones(length(ns),1); 
    temp(temp==0) = inf; 
    gg.emean2 = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.emean;
    temp = gg.fmean(:,ns)*ones(length(ns),1); 
    temp(temp==0) = inf; 
    gg.fmean2 = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.fmean;

    temp = gg.cmean(:,ns)*ones(length(ns),1); 
    temp(temp==0) = inf; 
    gg.cmean2 = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.cmean;

end
