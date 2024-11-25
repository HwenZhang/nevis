function [vv] = nevis_nodeoutflow(vv,aa,pp,gg,oo)
% calculate net outflow from nodes, sclaed with ps.x^2*ps.h/ps.t
%   ignores contribution due to changes in storage
%   this is really indended for looking at flux into boundary nodes gg.nbdy
%       -vv.divQ(gg.nbdy) gives flow into boundary nodes from channels
%       -vv.divqs(gg.nbdy) gives flow into boundary nodes from sheet
%
% 1 November 2015 edited from nodedischarge, mostly taken from nevis_backbone, R2

    if ~isfield(vv,'Qx'), 
        if isfield(vv,'nbdy'), gg = nevis_label(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % re-mask
        aa = nevis_inputs(vv.t,aa,pp,gg,oo);
        oo.evaluate_variables = 1; [vv] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
%         oo.evaluate_variables = 1; [vv0] = nevis_backbone(inf,vv0,vv0,aa,pp,gg,oo); % expand solution variables
    end

    %% restrict mean operators [ to take means over active points ; must take care to only numptiply with such points ]
    temp = gg.nmeanx(:,gg.ein)*ones(length(gg.ein),1); temp(temp==0) = inf; nmeanxin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeanx;
    temp = gg.nmeany(:,gg.fin)*ones(length(gg.fin),1); temp(temp==0) = inf; nmeanyin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeany;
    temp = gg.nmeans(:,gg.cin)*ones(length(gg.cin),1); temp(temp==0) = inf; nmeansin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeans;
    temp = gg.nmeanr(:,gg.cin)*ones(length(gg.cin),1); temp(temp==0) = inf; nmeanrin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeanr;

    qsx = vv.qsx;
    qsy = vv.qsy;
    qex = vv.qex;
    qey = vv.qey;
    Qx = vv.Qx;
    Qy = vv.Qy;
    Qs = vv.Qs;
    Qr = vv.Qr;
    Xicx = vv.Xicx;
    Xicy = vv.Xicy;
    Xics = vv.Xics;
    Xicr = vv.Xicr;
    Xix = vv.Xix;
    Xiy = vv.Xiy;
    Xis = vv.Xis;
    Xir = vv.Xir;
%     hs = vv.hs;
%     he = vv.he;
%     hv = vv.hv;
%     hm = vv.hm;
%     Sx = vv.Sx;
%     Sy = vv.Sy;
%     Ss = vv.Ss;
%     Sr = vv.Sr;
%     
%     dt = vv.t-vv0.t;  if dt==0, dt = inf; end % timestep
%     hs_old = vv0.hs;
%     he_old = vv0.he;
%     hv_old = vv0.hv;
%     hm_old = vv0.hm;
%     Sx_old = vv0.Sx;
%     Sy_old = vv0.Sy;
%     Ss_old = vv0.Ss;
%     Sr_old = vv0.Sr;
    
    m = aa.m;
    E = aa.E;
    area = gg.Dx.*gg.Dy; % node area
    
    vv.divqs = pp.c4*( gg.nddx(:,:)*qsx + gg.nddy(:,:)*qsy ).*area;             % area*h_t due to sheet flow
    vv.divqe = pp.c5*( gg.nddx(:,:)*qex + gg.nddy(:,:)*qey ).*area;             % area*h_t due to elastic sheet flow   
    vv.source = - ( pp.c6*( m ) + pp.c7*( E ) ).*area;                          % area*h_t due to source
    vv.divQ = pp.c9*( (gg.nddx(:,:)*Qx).*gg.Dy.^(-1) + (gg.nddy(:,:)*Qy).*gg.Dx.^(-1) ... % area*h_t due to channel flow
             + (gg.ndds(:,:)*Qs) + (gg.nddr(:,:)*Qr) ).*area ;
    vv.diss = - pp.c11*( (nmeanxin(:,gg.ein)*(Xicx(gg.ein)+Xix(gg.ein))).*gg.Dy.^(-1) + (nmeanyin(:,gg.fin)*(Xicy(gg.fin)+Xiy(gg.fin))).*gg.Dx.^(-1) ... % area*h_t due to dissipation
                + (nmeansin(:,gg.cin)*(Xics(gg.cin)+Xis(gg.cin))).*gg.Ds.*gg.Dx.^(-1).*gg.Dy.^(-1) + (nmeanrin(:,gg.cin)*(Xicr(gg.cin)+Xir(gg.cin))).*gg.Dr.*gg.Dx.^(-1).*gg.Dy.^(-1) ).*area;
%     vv.storage = ( - pp.c1*(hs-hs_old).*dt^(-1) ...
%          - pp.c2*(he-he_old).*dt.^(-1) ...
%          - pp.c3*(hv-hv_old).*dt.^(-1) ...
%          - pp.c3*(hm-hm_old).*dt.^(-1) ...
%            - pp.c8*nmeanxin(:,gg.ein)*(Sx(gg.ein)-Sx_old(gg.ein)).*dt^(-1).*gg.Dy.^(-1) ...
%            - pp.c8*nmeanyin(:,gg.fin)*(Sy(gg.fin)-Sy_old(gg.fin)).*dt^(-1).*gg.Dx.^(-1) ...
%            - pp.c8*nmeansin(:,gg.cin)*(Ss(gg.cin)-Ss_old(gg.cin)).*dt^(-1).*gg.Ds.*gg.Dx.^(-1).*gg.Dy.^(-1) ...
%            - pp.c8*nmeanrin(:,gg.cin)*(Sr(gg.cin)-Sr_old(gg.cin)).*dt^(-1).*gg.Dr.*gg.Dx.^(-1).*gg.Dy.^(-1) ).*area;        
            
    vv.outflow = vv.divqs+vv.divqe+vv.divQ+vv.diss+vv.source; %+vv.storage; 

end