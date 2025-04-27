%% set up a figure
% casename = oo.casename;
casename = 'blister_2d_alpha0_2_kappa0_mu1e_3_E8_8e9_V1e7_H1_0e3';
load(['./results/' casename '/' casename])
oo.fn = ['/',casename];                      % filename (same as casename)
oo.rn = [oo.root,oo.results,oo.fn];             % path to the case results
path = [oo.rn,'/'];

%% read in the screenshot at the intial timestep
t_init = 1480; t_end = 1600; %time range for animation
formatSpec = '%04d';
tframe = 30.0;
nframe = floor(tframe/0.02)+1;
vva = load([path num2str(nframe,formatSpec)], 'vv');
vva = vva.vv;
aa = nevis_inputs(vva.t,aa,vva,pp,gg,oo);
[vv2] = nevis_backbone(inf,vva,vva,aa,pp,gg,oo);     % expand solution variables
vv2 = nevis_nodedischarge(vv2,aa,pp,gg,oo);          % calculate node discharge
qnet = ps.qs*(vv2.qs + vv2.qe + vv2.qQ + vv2.Q);
slice = 1.0;
xx = (ps.x/10^3)*gg.nx; % x grid in km
yy = (ps.x/10^3)*gg.ny;   
xx(gg.nout) = NaN;
yy(gg.nout) = NaN;

%% asymptotic solution
% blister radius (only works for a single blister)
nonzeroIdx = find(vva.hb > 1e-3*vva.hb(pp.ni_l));
[~, localIdx] = min(vva.hb(nonzeroIdx));
minidx = nonzeroIdx(localIdx);
if isempty(minidx)
    Rb = 0;
else
    Rb = ps.x*((gg.nx(pp.ni_l)-gg.nx(minidx)).^2 + (gg.ny(pp.ni_l)-gg.ny(minidx)).^2).^(1/2);
end
disp(['Radius of the blister is ' num2str(Rb) '.']);

% asymptotic solution
Vb = ps.x^2*ps.hb*sum(vva.hb.*gg.Dx.*gg.Dy);
hb_centre_analytical = 3*Vb/(pi*Rb^2)*(1-((xx(gg.ny==slice)*1e3-1e4)/Rb).^2).^2 .*(abs(xx(gg.ny==slice)*1e3-1e4)<Rb);
pb_centre_analytical = 192/pi*pd.B*Vb/Rb^6*ones(size(hb_centre_analytical))/1e6;

hb_centre = ps.hb*vva.hb(gg.ny==slice);
hb_xxxx = (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*(gg.nddx*gg.eddx + gg.nddy*gg.fddy)*vva.hb;
hb_xxxx_centre = ps.hb/ps.x^4*hb_xxxx(gg.ny==slice);
pb_centre = pd.B*hb_xxxx_centre/1e6;

%% get the cross section
figure();

figtiled = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax = nexttile(1);
p11 = plot(ax,xx(gg.ny==slice),hb_centre,LineStyle='-',LineWidth=1.5,Marker='o');
p11.YDataSource = 'hb_centre';

hold on
p12 = plot(ax,xx(gg.ny==slice),hb_centre_analytical,LineStyle='--',LineWidth=1.5);
p12.YDataSource = 'hb_centre_analytical';
legend('numerical','analytical','Location','northeast',FontSize=15)

ylabel('h_b (m)')
xlabel('x (km)')
xlim([5 15])
ylim([0 5])

time = ['t=' num2str(vva.t*ps.t/pd.td,'%.1f') ' d'];
ttext = text(0.1,0.8,time,'Units','normalized');
ttext.FontSize=15;

ax = nexttile(2);
p21 = plot(ax,xx(gg.ny==slice),pb_centre,LineStyle='-',LineWidth=1.5);
p21.YDataSource = 'pb_centre';
hold on
p21 = plot(ax,xx(gg.ny==slice),pb_centre_analytical,LineStyle='--',LineWidth=1.5);
p21.YDataSource = 'pb_centre_analytical';
legend('numerical','analytical','Location','northeast',FontSize=15)

ylabel('p_{b} (Pa)')
xlabel('x (km)')
xlim([5 15])
ylim([-10 10])

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

    % blister radius (only works for a single blister)
    nonzeroIdx = find(vva.hb > 1e-3*vva.hb(pp.ni_l));
    [~, localIdx] = min(vva.hb(nonzeroIdx));
    minidx = nonzeroIdx(localIdx);
    if isempty(minidx)
        Rb = 0;
    else
        Rb = ps.x*((gg.nx(pp.ni_l)-gg.nx(minidx)).^2 + (gg.ny(pp.ni_l)-gg.ny(minidx)).^2).^(1/2);
    end
    disp(['Radius of the blister is ' num2str(Rb) '.']);

    % asymptotic solution
    Vb = ps.x^2*ps.hb*sum(vva.hb.*gg.Dx.*gg.Dy);
    hb_centre_analytical = 3*Vb/(pi*Rb^2)*(1-((xx(gg.ny==slice)*1e3-1e4)/Rb).^2).^2 .*(abs(xx(gg.ny==slice)*1e3-1e4)<Rb);
    pb_centre_analytical = 192/pi*pd.B*Vb/Rb^6*ones(size(hb_centre_analytical))/1e6;

    hb_centre = ps.hb*vva.hb(gg.ny==slice);
    hb_xxxx = (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*(gg.nddx*gg.eddx + gg.nddy*gg.fddy)*vva.hb;
    hb_xxxx_centre = ps.hb/ps.x^4*hb_xxxx(gg.ny==slice);
    pb_centre = pd.B*hb_xxxx_centre/1e6;
    
    set(ttext,{'string'},{['t=' num2str(vva.t*ps.t/(24*60*60),'%.1f'), ' d']})  %notice the column vector of new values

    refreshdata
    drawnow
    % pause(0.2)
    frame = getframe(gcf);
    writeVideo(v,frame)
end

close(v)
