% Script to plot the 1-dimensiona results from nevis
% By Hanwen Zhang on 25/11/2024

fs = 12; % fontsize
% final solutions
figure(2); clf; 
    subplot(211);
    % hydraulic potential 
        p1 = plot(gg.nx(gg.ns),vv.phi(gg.ns),'-'); 
        set(p1,'LineWidth',1.5,...
        'Color',[0 0 1]); hold on
        
        p2 = plot(gg.nx(gg.ns),aa.phi_0(gg.ns),'-'); 
        set(p2,'LineWidth',1.5,...
        'Color',[1 0 0]);
    
        xlabel('Distance');
        ylabel('Hydraulic potential');
    
        % legend
        h=legend([p1,p2],'$\phi$','$\phi_0$');
        h.Interpreter='latex';
        h.FontSize=fs+2;
        h.Location='northeast';

    subplot(212);
    % hydraulic potential 
        p1 = plot(gg.nx(gg.ns),vv.hs(gg.ns),'-'); 
        set(p1,'LineWidth',1.5,...
        'Color',[0 0 1]); hold on
        
        xlabel('Distance');
        ylabel('Hydraulic potential');
    
        % legend
        h=legend(p1,'$h_s$');
        h.Interpreter='latex';
        h.FontSize=fs+2;
        h.Location='northwest';

figure(3); clf;   
    subplot(211);
    % channel/sheet discharge 
    p1 = plot(gg.nx(gg.ns),vv2.qQ(gg.ns),'-'); 
    set(p1,'LineWidth',1.5,...
    'Color',[0 0 1]); hold on
    
    p2 = plot(gg.nx(gg.ns),vv2.qs(gg.ns)+vv2.qe(gg.ns),'-'); 
    set(p2,'LineWidth',1.5,...
    'Color',[1 0 0]);

    xlabel('Distance');
    ylabel('Hydraulic potential');

    % legend
    h=legend([p1,p2],'$Q_{S}$','$q_{s}+q_{e}$');
    h.Interpreter='latex';
    h.FontSize=fs+2;
    h.Location='southeast';

    % plot(gg.nx(gg.ns),vv2.qQ(gg.ns),'-',gg.nx(gg.ns),vv2.qs(gg.ns)+vv2.qe(gg.ns),'-');  
    xlabel('Distance');
    ylabel('Discharge');