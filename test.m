%% set up a figure
casename = oo.casename;
casename = 'test_2d';
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                      % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
path = [oo.rn,'/'];

%% read in the screenshot at the intial timestep
t_init = 740; t_end = 800; %time range for animation
formatSpec = '%04d';
tframe = 30.0;
nframe = floor(tframe/0.04)+1;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ + vv2.Q);

xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;  
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% get the cross section
figure();
figtiled = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
hb_centre = ps.hb*vva.hb(gg.ny==0.5);
hb_xxxx = (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*(gg.nddx*gg.eddx + gg.nddy*gg.fddy)*vva.hb;
hb_xxxx_centre = ps.hb/ps.x^4*hb_xxxx(gg.ny==0.5);

ax = nexttile(1);
p11 = plot(ax,xx(gg.ny==0.5),hb_centre,LineStyle='-',LineWidth=1.5);
p11.YDataSource = 'hb_centre';
ylabel('h_b (m)')
xlabel('x (km)')
xlim([30 50])
ylim([0 20])

time = ['t=' num2str(vva.t*ps.t/pd.td,'%.1f') ' d'];
ttext = text(32,3.5,time);
ttext.FontSize=15;

ax = nexttile(2);
p22 = plot(ax,xx(gg.ny==0.5),hb_xxxx_centre,LineStyle='-',LineWidth=1.5);
p22.YDataSource = 'hb_xxxx_centre';
ylabel('h_{b,xxxx} (m^-3)')
xlabel('x (km)')
xlim([30 50])
ylim([-1e-10 1e-10])

%% make video
v = VideoWriter(['./results/videos/' 'test'],'MPEG-4');
v.FrameRate = 1;
open(v)
for i_t = t_init:t_end
    disp(['Frame ',num2str(i_t-t_init),' / ',num2str(t_end-t_init),' ...']);
    %% load timestep
    vva = load([path num2str(i_t,formatSpec)], 'vv');
    vva = vva.vv;
    aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
    [vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo); % expand solution variables
    vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo); % calculate node discharge

    hb_centre = ps.hb*vva.hb(gg.ny==0.5);
    hb_xxxx = (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*(gg.nddx*gg.eddx + gg.nddy*gg.fddy)*vva.hb;
    hb_xxxx_centre = ps.hb/ps.x^4*hb_xxxx(gg.ny==0.5);
    
    set(ttext,{'string'},{['t=' num2str(vva.t*ps.t/(24*60*60),'%.1f'), ' d']})  %notice the column vector of new values

    refreshdata
    drawnow
    % pause(0.2)
    frame = getframe(gcf);
    writeVideo(v,frame)
end

close(v)
