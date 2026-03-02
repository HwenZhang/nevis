% nevis_plot_velocity.m
% function to plot the ice velocity field for maintenance during timesteps
function nevis_plot_velocity(gg,vv,num)
    u = vv.u;
    v = vv.v;
    hf = findobj('Type','figure','Number',num);
    if isempty(hf)
        hf = figure(num);
        hf.Visible = 'off';
    else
        set(0, 'CurrentFigure', hf);  % make current without raising
    end
    clf(hf);
    ax = axes(hf);
    un = 100*reshape(gg.nmeanx2(:,gg.es2)*u(gg.es2),gg.nI,gg.nJ);
    vn = 100*reshape(gg.nmeany2(:,gg.fs2)*v(gg.fs2),gg.nI,gg.nJ);
    Un = sqrt(un.^2 + vn.^2);
    Un(gg.nout) = NaN;

    hold(ax, 'on');
    pcolor(ax, gg.nx,gg.ny,Un); shading(ax, 'flat'); colorbar(ax);
    clim(ax, [0 300]);
    colormap(ax, jet);
    plot(ax, gg.ex(gg.ebdy2),gg.ey(gg.ebdy2),'r+','markersize',6);
    plot(ax, gg.fx(gg.fbdy2),gg.fy(gg.fbdy2),'b+','markersize',6);
    plot(ax, gg.nx(gg.nout),gg.ny(gg.nout),'k.','markersize',6);
    quiver(ax, gg.nx(1:4:end,1:4:end),gg.ny(1:4:end,1:4:end),un(1:4:end,1:4:end),vn(1:4:end,1:4:end), 'k');
    xlabel(ax, 'x (km)'); ylabel(ax, 'y (km)');
    if num == 1
        title(ax, 'Coupled velocity');
    else
        title(ax, 'SSA velocity');
    end
    axis(ax, 'equal'); axis(ax, 'tight');
    hf.Visible = 'on';
    drawnow
end 