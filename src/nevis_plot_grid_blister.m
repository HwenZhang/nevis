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
    % plot(nx(nbdy),ny(nbdy),'r.','markersize',8); 
    % plot(nx(nin),ny(nin),'k.','markersize',8);
    % plot(ex(ebdy),ey(ebdy),'r+');
    % plot(ex(ein),ey(ein),'k+');
    % plot(fx(fbdy),fy(fbdy),'g+');
    % plot(fx(fin),fy(fin),'k+');
    % plot(cx(cbdy),cy(cbdy),'rx');
    % plot(cx(cin),cy(cin),'kx');
    % 
    plot(nx(nbdy_blister),ny(nbdy_blister),'r.','markersize',8); 
    plot(nx(nin_blister),ny(nin_blister),'k.','markersize',8);
    plot(ex(ebdy_blister),ey(ebdy_blister),'r+');
    plot(ex(ein_blister),ey(ein_blister),'k+');
    plot(fx(fbdy_blister),fy(fbdy_blister),'g+');
    plot(fx(fin_blister),fy(fin_blister),'k+');
    plot(cx(cbdy_blister),cy(cbdy_blister),'rx');
    plot(cx(cin_blister),cy(cin_blister),'kx');
    % % 
    plot(nx(ni),ny(ni),'go');
    grid on
   
end