function [vv] = nevis_nodedischarge(vv,aa,pp,gg,oo)
% calculate mean fluxes on nodes from values on edges
%  [ previously run nevis_backbone with oo.evalulate_variables = 1 to
%  calculate fluxes ]
%
% 19 April 2014 : taken from hydro_nodedischarge, corrected sign of Qr for y component
% 2 November 2015 : edited to provide means of directed discharge on nodes

    % restrict mean operators [ to take means over active points ; must take care to only numptiply with such points ]
    temp = gg.nmeanx(:,gg.ein)*ones(length(gg.ein),1); temp(temp==0) = inf; nmeanxin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeanx;
    temp = gg.nmeany(:,gg.fin)*ones(length(gg.fin),1); temp(temp==0) = inf; nmeanyin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeany;
    temp = gg.nmeans(:,gg.cin)*ones(length(gg.cin),1); temp(temp==0) = inf; nmeansin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeans;
    temp = gg.nmeanr(:,gg.cin)*ones(length(gg.cin),1); temp(temp==0) = inf; nmeanrin = sparse(1:length(temp),1:length(temp),temp.^(-1),length(temp),length(temp))*gg.nmeanr;

    ein = gg.ein;
    fin = gg.fin;
    cin = gg.cin;
    qsx = vv.qsx;
    qsy = vv.qsy;
    qex = vv.qex;
    qey = vv.qey;
    qbx = vv.qbx;
    qby = vv.qby;
    Qx = vv.Qx;
    Qy = vv.Qy;
    Qs = vv.Qs;
    Qr = vv.Qr;
    Sx = vv.Sx;
    Sy = vv.Sy;
    Ss = vv.Ss;
    Sr = vv.Sr;
    
    % average fluxes on nodes
    % scaled with ps.qs
    vv.qsxn = nmeanxin(:,ein)*qsx(ein); 
    vv.qsyn = nmeanyin(:,fin)*qsy(fin);
    % scaled with ps.qs
    vv.qexn = pp.c5/pp.c4*nmeanxin(:,ein)*qex(ein); 
    vv.qeyn = pp.c5/pp.c4*nmeanyin(:,fin)*qey(fin);
    % scaled with ps.qs
    vv.qbxn = nmeanxin(:,ein)*qbx(ein); 
    vv.qbyn = nmeanyin(:,fin)*qby(fin);

    % channel discharge on nodes converted to 2d flux, scaled with ps.qs
    % [ divide by effective perpendicular area, eg (nmeansin(:,cin)*Qs(cin))./(2*Dx.*Dy./Ds) etc, then x component is qs.*(Dx./Ds) ]
    % [ note dds goes forward in x and y, ddr goes forward in x, backward in y ]
    vv.qQxn = pp.c9/pp.c4*( (nmeanxin(:,ein)*Qx(ein)).*gg.Dy.^(-1) ...
        + (nmeansin(:,cin)*Qs(cin)).*gg.Dy.^(-1)/2 ...
        + (nmeanrin(:,cin)*Qr(cin)).*gg.Dy.^(-1)/2 );
    vv.qQyn = pp.c9/pp.c4*( (nmeanyin(:,fin)*Qy(fin)).*gg.Dx.^(-1) ...
        + (nmeansin(:,cin)*Qs(cin)).*gg.Dx.^(-1)/2 ...
        - (nmeanrin(:,cin)*Qr(cin)).*gg.Dx.^(-1)/2 );
    % channel discharge on nodes, scaled with ps.Q
    % [ as above but multiply by Dy/Dx ]
    vv.Qxn = (nmeanxin(:,ein)*Qx(ein)) ...
        + (nmeansin(:,cin)*Qs(cin))/2 ...
        + (nmeanrin(:,cin)*Qr(cin))/2; ...
    vv.Qyn = (nmeanyin(:,fin)*Qy(fin)) ...
        + (nmeansin(:,cin)*Qs(cin))/2 ...
        - (nmeanrin(:,cin)*Qr(cin))/2;
    
    % absolute values, scaled witn ps.qs
    vv.qs = (vv.qsxn.^2+vv.qsyn.^2).^(1/2);
    vv.qe = (vv.qexn.^2+vv.qeyn.^2).^(1/2);
    vv.qb = (vv.qbxn.^2+vv.qbyn.^2).^(1/2);
    vv.qQ = (vv.qQxn.^2+vv.qQyn.^2).^(1/2);
    vv.Q = (vv.Qxn.^2+vv.Qyn.^2).^(1/2);
    
    % volume of channels at nodes, scaled with ps.h*ps.x^2
    vv.VS = pp.c8*( (nmeanxin(:,ein)*Sx(ein)).*gg.Dx + (nmeanyin(:,fin)*Sy(fin)).*gg.Dy ...
        + (nmeansin(:,cin)*Ss(cin)).*gg.Ds + (nmeanrin(:,cin)*Sr(cin)).*gg.Dr );

end
