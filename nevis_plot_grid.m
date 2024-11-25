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
    
    plot(nx(ni),ny(ni),'go');
    grid on
   
end