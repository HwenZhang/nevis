function nevis_plot_grid(gg,ni)
% plot grid defined in struct gg
% if ni included, highlighs nodes labelled ni
%
% IJH 31 July 2014

if nargin<2, ni = []; end

nx = gg.nx;
ny = gg.ny;
ex = gg.ex;
ey = gg.ey;
fx = gg.fx;
fy = gg.fy;
cx = gg.cx;
cy = gg.cy;

nout = gg.nout;
eout = gg.eout;
fout = gg.fout;
cout = gg.cout;
nbdy = gg.nbdy;
ebdy = gg.ebdy;
fbdy = gg.fbdy;
cbdy = gg.cbdy;
nin = gg.nin;
ein = gg.ein;
fin = gg.fin;
cin = gg.cin;

nbdy_blister = gg.nbdy_blister;
nin_blister = gg.nin_blister;
ebdy_blister = gg.ebdy_blister;
fbdy_blister = gg.fbdy_blister;
cbdy_blister = gg.cbdy_blister;
ein_blister = gg.ein_blister;
fin_blister = gg.fin_blister;
cin_blister = gg.cin_blister;

figure;
    % external
    plot(nx(nout),ny(nout),'.','color',.9*[1 1 1],'markersize',8); hold on;
    plot(ex(eout),ey(eout),'+','color',.9*[1 1 1]);
    plot(fx(fout),fy(fout),'+','color',.9*[1 1 1]);
    plot(cx(cout),cy(cout),'x','color',.9*[1 1 1]);
    
    % bound and internal
    plot(nx(nbdy),ny(nbdy),'r.','markersize',8); 
    plot(nx(nin),ny(nin),'k.','markersize',8);
    plot(ex(ebdy),ey(ebdy),'r+');
    plot(ex(ein),ey(ein),'k+');
    plot(fx(fbdy),fy(fbdy),'g+');
    plot(fx(fin),fy(fin),'k+');
    plot(cx(cbdy),cy(cbdy),'rx');
    plot(cx(cin),cy(cin),'kx');
    
    % plot(nx(nbdy_blister),ny(nbdy_blister),'r.','markersize',8); 
    % plot(nx(nin_blister),ny(nin_blister),'k.','markersize',8);
    % plot(ex(ebdy_blister),ey(ebdy_blister),'r+');
    % plot(ex(ein_blister),ey(ein_blister),'k+');
    % plot(fx(fbdy_blister),fy(fbdy_blister),'g+');
    % plot(fx(fin_blister),fy(fin_blister),'k+');
    % plot(cx(cbdy_blister),cy(cbdy_blister),'rx');
    % plot(cx(cin_blister),cy(cin_blister),'kx');

    % plot(nx(gg.nbdyx2),ny(gg.nbdyx2),'g.','markersize',8); 
    % plot(nx(gg.nbdyy2),ny(gg.nbdyy2),'g.','markersize',8); 

    % get the diff between ns2 and nin2
    % ndiff = setdiff(gg.nin2, gg.nin);
    % plot(nx(ndiff),ny(ndiff),'mo','markersize',8);
    % plot(nx(gg.nin2),ny(gg.nin2),'rx','markersize',8);
    % plot(ex(gg.ebdy2),ey(gg.ebdy2),'r+');
    % plot(ex(gg.ein2),ey(gg.ein2),'k+');
    % plot(fx(gg.fbdy2),fy(gg.fbdy2),'g+');
    % plot(fx(gg.fin2),fy(gg.fin2),'k+');
    % plot(cx(gg.cbdy2),cy(gg.cbdy2),'rx');
    % plot(cx(gg.cin2),cy(gg.cin2),'kx');

    invec = ismember(gg.nin, gg.nin2);
    missing = gg.nin(~invec);
    plot(nx(missing),ny(missing),'bo','markersize',8);

    plot(nx(ni),ny(ni),'go');
    legend('nout','eout','fout','cout','nbdy','nin','ebdy','ein','fbdy','fin','cbdy','cin','ni');
    grid on
   
end