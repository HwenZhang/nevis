function nevis_plot_grid_ice(gg,ni)
% plot grid defined in struct gg
% if ni included, highlighs nodes labelled ni
%
% IJH 31 July 2014
% 14 Oct 2025 - ice version, edited from nevis_plot_grid

if nargin<2, ni = []; end

nx = gg.nx;
ny = gg.ny;
ex = gg.ex;
ey = gg.ey;
fx = gg.fx;
fy = gg.fy;
cx = gg.cx;
cy = gg.cy;

nout = gg.nout2;
eout = gg.eout2;
fout = gg.fout2;
cout = gg.cout2;
nbdyx = gg.nbdyx2;
nbdyy = gg.nbdyy2;
nbdy = unique([nbdyx; nbdyy]);
ebdy = gg.ebdy2;
fbdy = gg.fbdy2;
cbdy = gg.cbdy2;
nin = gg.nin2;
ein = gg.ein2;
fin = gg.fin2;
cin = gg.cin2;

figure;
    plot(nx(nout),ny(nout),'.','color',.9*[1 1 1],'markersize',8); hold on;
    plot(ex(eout),ey(eout),'x','color',.9*[1 1 1]);
    plot(fx(fout),fy(fout),'x','color',.9*[1 1 1]);
    plot(cx(cout),cy(cout),'+','color',.9*[1 1 1]);
    
    plot(nx(nbdy),ny(nbdy),'r.','markersize',8); 
    plot(nx(nin),ny(nin),'k.','markersize',8);

    plot(ex(ebdy),ey(ebdy),'r+');
    plot(fx(fbdy),fy(fbdy),'b+');
    plot(cx(cbdy),cy(cbdy),'gx');

    
    plot(nx(ni),ny(ni),'go');
   
end