% nevis_plot_velocity.m
% function to plot the ice velocity field for maintenance during timesteps
function nevis_plot_velocity(gg,u,v)
    figure(1); clf;
    un = 100*reshape(gg.nmeanx2(:,gg.es2)*u(gg.es2),gg.nI,gg.nJ);
    vn = 100*reshape(gg.nmeany2(:,gg.fs2)*v(gg.fs2),gg.nI,gg.nJ);
    Un = sqrt(un.^2 + vn.^2);
    Un(gg.nout) = NaN;

    hold on
    pcolor(gg.nx,gg.ny,Un); shading flat; colorbar; hold on
    clim([0 300]);
    plot(gg.ex(gg.ebdy2),gg.ey(gg.ebdy2),'r+','markersize',6);
    plot(gg.fx(gg.fbdy2),gg.fy(gg.fbdy2),'b+','markersize',6);
    plot(gg.nx(gg.nout),gg.ny(gg.nout),'k.','markersize',6);
    quiver(gg.nx(1:4:end,1:4:end),gg.ny(1:4:end,1:4:end),un(1:4:end,1:4:end),vn(1:4:end,1:4:end), 'k');
    xlabel('x (km)'); ylabel('y (km)'); title('Initial ice velocity'); axis equal; axis tight;
    drawnow
end 