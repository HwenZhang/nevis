function nevis_imagesc(z,gg,ii)
% reshape vector z and plot as scaled image
%   z field to plot
%   gg struct containing grid
%   ii node type to plot on [optional: otherwise guessed from size of z]
%
% IJH 14 August 2014 : taken from hydro_imagesc

    if nargin<3,
        if length(z)==gg.nIJ, ii = 'n'; 
        elseif length(z)==gg.eIJ, ii = 'e';
        elseif length(z)==gg.fIJ, ii = 'f';
        elseif length(z)==gg.cIJ, ii = 'c';
        end
    end
    
    if strcmp(ii,'n'), z = reshape(z,gg.nI,gg.nJ); x = gg.nx; y = gg.ny; z(gg.nout) = NaN;
    elseif strcmp(ii,'e'), z = reshape(z,gg.eI,gg.eJ); x = gg.ex; y = gg.ey; z(gg.eout) = NaN;
    elseif strcmp(ii,'f'), z = reshape(z,gg.fI,gg.fJ); x = gg.fx; y = gg.fy; z(gg.fout) = NaN;
    elseif strcmp(ii,'c'), z = reshape(z,gg.cI,gg.cJ); x = gg.cx; y = gg.cy; z(gg.cout) = NaN;
    end
    
    figure;
    imagesc(x(:,1),y(1,:),z'); 
    colorbar; shg; 
    xlabel('x');
    ylabel('y');
    set(gca,'YDir','normal');
end