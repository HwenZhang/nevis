function reduce_colorbar_size(cc,xfac,yfac,xoff,yoff)
% reduce size of colorbar cc by factor xfac in x direction and yfac in y direction

    % offset factors
    if nargin<5
        xoff = []; yoff = [];
    end
    if isempty(xoff)
        xoff = 0.5*(1-xfac);
    end
    if isempty(yoff)
        yoff = 0.5*(1-yfac);
    end
    
    % keep current axes position
    set(gca,'Position',get(gca,'Position'));
    
    % change position of colorbar
    pos = get(cc,'Position');
    pos(1) = pos(1)+xoff*pos(3);
    pos(2) = pos(2)+yoff*pos(4);
    pos(3) = xfac*pos(3);
    pos(4) = yfac*pos(4);
    set(cc,'Position',pos);
    
end