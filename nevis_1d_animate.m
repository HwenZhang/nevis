% Script to animate the 1-dimensiona results from nevis
% By Hanwen Zhang on 25/11/2024

figure(4);
path = '/Users/eart0487/nevis/nevis_1d_example/';
formatSpec = '%04d';
vva = load([path num2str(1,formatSpec)], 'vv');
vva = vva.vv;

subplot(211);
ani_phi = vva.phi(gg.ns);
ani_phi0 = aa.phi_0(gg.ns);
ani_hs = vva.hs(gg.ns);

% hydraulic potential 
    p1 = plot(gg.nx(gg.ns),ani_phi,'-'); 
    set(p1,'LineWidth',1.5,...
    'Color',[0 0 1]); hold on
    p1.YDataSource = 'ani_phi';
    
    p2 = plot(gg.nx(gg.ns),ani_phi0,'-'); 
    set(p2,'LineWidth',1.5,...
    'Color',[1 0 0]);
    p2.YDataSource = 'ani_phi';

    xlabel('Distance');
    ylabel('Hydraulic potential');
    ylim([0 1]) 

    % legend
    h=legend([p1,p2],'$\phi$','$\phi_0$');
    h.Interpreter='latex';
    h.FontSize=fs+2;
    h.Location='northeast';
    
    time = ['t=' num2str(vva.t,'%.1f')];
    ttext = text(0.1,0.8,time);
    ttext.FontSize=12;
    grid on

% sheet thickness
subplot(212);
    p3 = plot(gg.nx(gg.ns),ani_hs,'-'); 
    set(p3,'LineWidth',1.5,...
    'Color',[0 0 1]); hold on
    p3.YDataSource = 'ani_hs';

    xlabel('Distance');
    ylabel('Sheet thickness');
    ylim([0 3]) 

    % legend
    h=legend(p3,'$h_s$');
    h.Interpreter='latex';
    h.FontSize=fs+2;
    h.Location='northeast';
    grid on

for i=2:vv.ti_save
    % each xxxx.mat file contains a solution field vv
    vva = load([path num2str(i,formatSpec)], 'vv');
    vva = vva.vv;

    ani_phi = vva.phi(gg.ns);
    ani_phi0 = aa.phi_0(gg.ns);
    ani_hs = vva.hs(gg.ns);
    
    set(ttext,{'string'},{['t=' num2str(vva.t,'%.1f')]})  %notice the column vector of new values

    refreshdata
    drawnow
    % pause(0.01);
end