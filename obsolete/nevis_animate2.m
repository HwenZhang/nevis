function nevis_animate2(fn,tis)
% mesh plot of potential
% [ not much use ? ]
% 26 August 2014: taken from nevis_animate

%% load initial timestep
if exist([fn,'/0000']), load([fn,'/0000']); else load(fn); end
load([fn,'/0001']);
if isfield(oo,'reversey'), reversey = oo.reversey; else reversey = 1; end

%% set up figure
figure(1); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 25 25]);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontSize',10);
axes_size = [0.1 0.1 0.8 0.8];

%% plot frames
for i_t = 1:length(tis)
disp(['Frame ',num2str(i_t),' / ',num2str(length(tis)),' ...']);

%% load timestep
load([fn,'/',int2four(tis(i_t))]);

%% extract new variables
if isfield(vv,'nbdy'), gg = nevis_mask2(gg,vv.nbdy); aa.phi = aa.phi_a(gg.nbdy); end % remask
aa = nevis_inputs(vv.t,aa,pp,gg,oo);
oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0; [vv2] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge
nevis_unpack(aa,vv2);
clear tt vv vv2

%% plot
clf; ax = axes('position',axes_size);

    xx = sc.x*gg.nx; yy = sc.x*gg.ny; 
    zz = sc.phi*reshape(phi,gg.nI,gg.nJ); zz(gg.nout) = NaN;
    mesh(xx,yy,zz);
    zlim([25 28]);

    text(0.05,0.95,['t = ',num2str(round(sc.t*t*100)/100)],'units','normalized');
    set(gcf,'Color','w');

    xlabel('x [ km ]');
    ylabel('y [ km ]');
    zlabel('\phi [ MPa ]');

    view(-145,20);
    shg; pause(0.05);
    
end

end