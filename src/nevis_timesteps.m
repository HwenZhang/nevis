function [tt,vv,info] = nevis_timesteps(t_span,vv,aa,pp,gg,oo)
% timestep nevis variables over time t_span
% Inputs:
%   t_span  start and end times and times for saving output
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   tt      struct containing time series of average quantities
%   vv      struct containing final solution variables
%   info    information about last computation [optional]
%
% IJH 14 August 2014 : largely taken from hydro_timesteps_diag
%       1 May 2015 : add possibility for time averaging between timesteps
%       26 October 2015 : added calculation of elastic sheet volume to tt
%       17 Nov 2015 : added possibility of prescribing potential aa.phi_b

if ~isfield(pp,'phi_s'), pp.phi_s = -inf; end

% ================= TIMESTEPPING OPTIONS =================
% update ice velocity
if ~isfield(oo,'include_ice'), oo.include_ice = 0; end

% ----------------- adaptive timestepping -----------------
% change timestep based on previous iterations/convergence
if ~isfield(oo,'change_timestep'), oo.change_timestep = 1; end
% allowable timestep
if ~isfield(oo,'dt_max'), oo.dt_max = 1e0; end  % maximum        
if ~isfield(oo,'dt_min'), oo.dt_min = 1e-6; end % mimimum    

% factor to change timestep by
if ~isfield(oo,'dt_factor'), oo.dt_factor = 2; end       
% increase timestep for this number of iterations or less
if ~isfield(oo,'small_iter'), oo.small_iter = 2; end            
% decrease timestep for this number of iterations or more 
if ~isfield(oo,'large_iter'), oo.large_iter = 10; end           

% verbose screen output
if ~isfield(oo,'verb'), oo.verb = 0; end                        
% adjust boundary nodes during timestepping
if ~isfield(oo,'adjust_boundaries'), oo.adjust_boundaries = 1; end  
% frequency of timesteps at which to adjust boundaries
if ~isfield(oo,'dti_boundaries'), oo.dti_boundaries = 4; end 

% ================= INPUT OPTIONS ==================
% use pp.input function(t) for moulin inputs
if ~isfield(oo,'input_function'), oo.input_function = 0; end    
% use pp.runoff function(t) for runoff
if ~isfield(oo,'runoff_function'), oo.runoff_function = 0; end  
% use distributed input for runoff even in presence of moulins
if ~isfield(oo,'distributed_input'), oo.distributed_input = 0; end  
% ================= SAVING OPTIONS =================
% save solution at timesteps in t_span to directory
if ~isfield(oo,'save_timesteps'), oo.save_timesteps = 1; end  
% make directory if it doesn't already exist
% if oo.save_timesteps && ~exist([oo.root,oo.fn],'dir')
%     mkdir([oo.root,oo.fn]); 
% end 
% save output at particular nodes pts_ni
if ~isfield(oo,'save_pts_all'), oo.save_pts_all = 0; end  
% nodes at which to save output
if ~isfield(oo,'pts_ni'), oo.pts_ni = []; end       
% save average pressure between timesteps           
if ~isfield(oo,'save_phi_av'), oo.save_phi_av = 0; end          

% ================= UPDATE BOUNDARIES =================
if oo.adjust_boundaries && isfield(vv,'nbdy')         
    % redefine boundary labels
    gg = nevis_label(gg,vv.nbdy,oo);          
    % fill in missin phi_b for older versions
    if ~isfield(aa,'phi_b'), aa.phi_b = max(aa.phi_a,pp.phi_s); end 
    % boundary conditions
    aa.phi = aa.phi_b(gg.nbdy);                    
end

% ================= TIMESTEPPING =================
t_start = t_span(1);
t_stop = t_span(end);
ti_save = 1;
t_save = t_span(ti_save);
ti = 1;                   
ti_boundaries = 1;

dt1 = oo.dt;    % suggested timestep
t = t_start;    % initial time
vv.t = t;       % initial time

vv.dt = dt1;    % initial timestep
tt = struct;
% if save average potential
if oo.save_phi_av         
    phi_av = vv.phi; % initial average pressure
    phi_min = vv.phi;
    phi_max = vv.phi;
end
% initial expanded variables vv2
[vv2,~,~,~,~,~,~,~,~,~,~] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); 

% time loop
while t<t_stop+oo.dt_min
    %% ice sheet dynamics
    %     if oo.include_ice
    %         % [ this option not operational yet ]
    %         N = max(0,phi_0-phi);  %[ ignore negative effective pressures ]
    %         
    %         [ubar,vbar,eta] = ice_force_balance_new2(b,H,A,N,ubar,vbar,eta,tauxx_bdy,tauyy_bdy,tauxy_bdy,pp,gg,oo,opts); 
    %         [ub,vb,us,vs,qx,qy,taubx,tauby,tauxxbar,tauyybar,tauxybar,tauxz,tauyz,ui,vi,tau] = ice_velocities_new2(b,H,A,N,ubar,vbar,eta,tauxx_bdy,tauyy_bdy,tauxy_bdy,pp,gg,oo,opts);
    %         us(gg.bdy.elbdy,:) = us(gg.bdy.elbdy+1,:); ui(gg.bdy.elbdy,:) = ui(gg.bdy.elbdy+1,:);  % [ xslope may not be well defined on elbdy : extrapolate values from next points ]
    %        %[~,~,wf,~] = ice_velocities_all(b,H,ub,vb,ui,vi,pp,gg,oo,opts); ws = wf(:,1); % vertical surface velocity
    % 
    %        if oo.ice_update,
    %            % update sliding speed [ assumes pd.u was taken as ps.u so that non-dimensional Ub = 1 ]
    %            Ub = ( (oo.nmeanx(:,gg.ein)*ub(gg.ein)).^2 + (oo.nmeany(:,gg.fin)*vb(gg.fin)).^2 ).^(1/2);
    %        end
    %        if oo.melt_update
    %            % update melt rate 
    %            m = ((pd.G+ ps.tau0*ps.u*( oo.nmeanx(:,gg.ein)*(taubx(gg.ein).*ub(gg.ein))+oo.nmeany(:,gg.fin)*(tauby(gg.fin).*vb(gg.fin)) ))/pd.rho_w/pd.L)/ps.m; 
    %        end
    %     end
    %     end

    %% update input
    aa = nevis_inputs(t,aa,vv,pp,gg,oo);

    %% update blister input and output
    % current timestep is dt, current solution is vv, precribed field aa
    % [aa,vv] = nevis_blister_distributed(aa,vv,pp,gg,oo);

    %% time series of average quantities
    tt(ti).ti = ti;
    tt(ti).t = t;
    tt(ti).m = sum(aa.m(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));  % basal melt, scaled with ps.m*ps.x^2
    % supraglacial input, scaled with ps.m*ps.x^2
    tt(ti).E = sum(aa.E(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));  % moulin input, scaled with ps.m*ps.x^2

    tt(ti).Qb_in = aa.Qb_in(pp.ni_l).*gg.Dx(pp.ni_l).*gg.Dy(pp.ni_l);                         % inflow to the blister
    tt(ti).Qb_dec = vv2.Qb_dec;                               % relaxation term of the blister
               
    tt(ti).pwb = vv.phi(pp.ni_l);                             % hydrulic potential at the lake
    tt(ti).hs_b = vv.hs(pp.ni_l);                             % sheet thickness at the lake

    % channel cross sectional area at the nodes
    temp = gg.nmeanx*vv.Sx;
    tt(ti).Sx_b = temp(pp.ni_l);                              % average channel cross section in x-direction
    temp = gg.nmeany*vv.Sy;
    tt(ti).Sy_b = temp(pp.ni_l);
    temp = gg.nmeans*vv.Ss;
    tt(ti).Ss_b = temp(pp.ni_l);
    temp = gg.nmeanr*vv.Sr;
    tt(ti).Sr_b = temp(pp.ni_l);

    tt(ti).Q_in = vv2.Q_in;                                   % inflow, scaled with ps.Q
    tt(ti).Q_out = vv2.Q_out;                                 % outflow, scaled with ps.Q
    tt(ti).Qb_out = vv2.Qb_out;                               % outflow from the blister scaled with ps.Q, from residual

    tt(ti).Q_outQ = vv2.Q_outQ;                               % channel outflow, scaled with ps.Q
    tt(ti).Q_outq = vv2.Q_outq;                               % sheet outflow, scaled with ps.Q
    tt(ti).Q_outb = vv2.Q_outb;                               % blister outflow, scaled with ps.Q, from calculation of the flux term

    tt(ti).phi = mean(vv.phi(gg.ns));                         % mean potential, scaled with ps.phi
    tt(ti).N = mean(aa.phi_0(gg.ns)-vv.phi(gg.ns));           % mean effective pressure, scaled with ps.phi
    tt(ti).hs = sum(vv.hs(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));% total cavity sheet volume, scaled with ps.h*ps.x^2
    tt(ti).Vb = sum(vv.hb(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));% total blister volume, scaled with ps.hb*ps.x^2

    disp(['Volume of the blister is ' num2str(tt(ti).Vb) '.']);
    disp(['The most negative effective pressure is ' num2str(min(aa.phi_0(gg.ns)-vv.phi(gg.ns))) '.']);

    % blister radius (only works for a single blister)
    % [~, maxIdx] = max(vv.hb);
    maxIdx = pp.ni_l(1); % the location of the lake/hydrofracture
    disp(['Max thickness of the blister is ' num2str(max(vv.hb)) '.']);
    disp(['Min thickness of the blister is ' num2str(min(vv.hb)) '.']);
    nonzeroIdx = find(vv.hb > 1e-4); % find non-zero thicknesses
    [~, localIdx] = max((gg.nx(nonzeroIdx)-gg.nx(maxIdx)).^2+(gg.ny(nonzeroIdx)-gg.ny(maxIdx)).^2);
    minidx = nonzeroIdx(localIdx);
    if isempty(minidx)
        tt(ti).Rb = 0;
    else
        tt(ti).Rb = ((gg.nx(maxIdx)-gg.nx(minidx)).^2 + (gg.ny(maxIdx)-gg.ny(minidx)).^2).^(1/2);
    end
    disp(['Radius of the blister is ' num2str(tt(ti).Rb) '.']);
   tt(ti).hb_max = vv.hb(maxIdx); % maximum blister thickness

    % total channel volume, scaled with ps.S*ps.x
    tt(ti).S = sum(vv.Sx(gg.ein).*(gg.emean(gg.ein,:)*gg.Dx))+...
               sum(vv.Sy(gg.fin).*(gg.fmean(gg.fin,:)*gg.Dy))+...
               sum(vv.Ss(gg.cin).*(gg.cmean(gg.cin,:)*gg.Ds))+...
               sum(vv.Sr(gg.cin).*(gg.cmean(gg.cin,:)*gg.Dr)); 
    % total elastic sheet volume, scaled with ps.h*ps.x^2
    tt(ti).he = sum(vv2.he(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns)); 
    if oo.save_pts_all
         tt(ti).pts_phi = vv.phi(oo.pts_ni);
         tt(ti).pts_hs = vv.hs(oo.pts_ni);
         tt(ti).pts_he = vv2.he(oo.pts_ni);
         tt(ti).pts_hb = vv.hb(oo.pts_ni);
         tt(ti).pts_pb = vv.pb(oo.pts_ni);
         tt(ti).pts_hc = vv2.hc(oo.pts_ni);
         if oo.include_ice
            tt(ti).pts_us = vv.us(oo.pts_ni);
         end
    end
    
    %% saving    
    if t >= t_save
        if oo.save_phi_av
            vv.phi_av = phi_av;
            vv.phi_min = phi_min;
            vv.phi_max = phi_max;
            phi_av = vv.phi; % reset average pressure
            phi_min = vv.phi;
            phi_max = vv.phi;
        end
        vv.ti = ti;
        vv.ti_save = ti_save;

        % disp(['nevis_timesteps: Saving ',fn,' ...']);
        % save(fn,'vv','tt','-append');
        disp('nevis_timesteps: Done');
        % save full solution
        if oo.save_timesteps
            fnt = [oo.rn,'/',int2four(ti_save)];
            disp(['nevis_timesteps: Saving ',fnt,' ...']);
            save(fnt,'vv'); 
            disp('nevis_timesteps: Done');
        end
        % [ add other saving options here ]
        if t>t_stop-oo.dt_min, break; end
        ti_save = ti_save+1;
        t_save = t_span(ti_save);  
    elseif t > t_stop-oo.dt_min
        break; 
    end

    %% timestep hydrology
    tic;
    accept = 0;
    decreased = 0;
    increased = 0;
    while ~accept
        dt = dt1; % suggested timestep
        % adjust timestep if this would take too close to t_save
        if t + dt > t_save-oo.dt_min, dt = t_save-t; end 
        disp(['nevis_timesteps: dt = ',num2str(dt),' ...']);
        [vv1,vv2,info] = nevis_timestep(dt,vv,aa,pp,gg,oo);  % time iteration
        
        %% check success and adjust size of timestep [ taken from hydro_timestep_diag ]
        if ~info.failed
            accept = 1; 
            t = t + dt;
            % vv1.hb = max(vv1.hb,0); % set negative hb to zero
            vv = vv1;
            vv.dt = dt;
            vv.t = t;
            comp_time = toc;
            disp(['nevis_timesteps: Done [ ',num2str(comp_time),' s, ',num2str(info.iter_new-1),' iterations ]']);
            
            % add code
            % if max(vv.hb) > 1e-2
            %     figure(2);                    % 或者任何你喜欢的 figure 编号
            %     subplot(2,1,1);
            %     plot(gg.nx, vv.hb, 'o-','LineWidth',1); 
            %     xlabel('x');
            %     ylabel('h_b');
            %     xlim([4.5 5.5]);
            %     title(sprintf('t = %.3g', t));
            %     grid on;

            %     subplot(2,1,2);
            %     kappa = (gg.nddx*gg.eddx + gg.nddy*gg.fddy)*vv.hb;
            %     % dRdt = 1.15*pd.B*1000/(144*pd.mu) * kappa.^3;

            %     plot(gg.nx, kappa, 'o-','LineWidth',1);
            %     xlabel('x');
            %     xlim([4.5 5.5]);
            %     ylabel('dR/dt (m/s)');
            %     grid on;

            %     drawnow;                       % 强制立即刷新
            %     pause(1.0);                    % 暂停 1 秒以便查看图形
            % end
            % end here

        end
        if oo.change_timestep && info.iter_new-1 <= oo.small_iter && ~info.failed && dt1 < oo.dt_max && ~decreased
            % too few iterations
            dt1 = min(oo.dt_factor*dt1,oo.dt_max); increased = 1;
            if oo.verb, disp(['nevis_timesteps: Increase suggested timestep to ',num2str(dt1)]); end
            continue;
        elseif oo.change_timestep && ( info.iter_new-1 >= oo.large_iter || info.failed ) && dt1 > oo.dt_min && ~increased
            % too many iterations
            dt1 = max(dt1/oo.dt_factor,oo.dt_min); decreased = 1;
            if oo.verb, disp(['nevis_timesteps: Decrease suggested timestep to ',num2str(dt1)]); end
            continue;
        end
        % If no convergence
        if ~accept
            comp_time = toc;
            disp(['nevis_timesteps: Failed [ ',num2str(comp_time),' s, ',num2str(info.iter_new-1),' iterations, dt = ',num2str(dt),' ]']);
            return;
        end
    end
    
    %% time averaging of pressure
    if oo.save_phi_av
        phi_av = (phi_av*(t-t_span(ti_save-1)-dt)+vv.phi*dt )/(t-t_span(ti_save-1));
        phi_min = min(phi_min,vv.phi);
        phi_max = max(phi_max,vv.phi);
    end
    
    %% check and adjust boundary nodes
    if oo.adjust_boundaries && ti>=ti_boundaries
        ti_boundaries = ti+oo.dti_boundaries;
        ni1 = gg.nbdy(vv2.R_bdy<0); % Dirichlet nodes with inflow
        if ~isfield(gg,'n1m'), gg.n1m = gg.n1; end % boundary nodes adjacent to margin
        % ni2 = gg.n1m(vv.phi(gg.n1m)-aa.phi_b(gg.n1m)>pp.p_a_reg); % boundary nodes with too high pressure
        ni2 = gg.n1m((vv.phi(gg.n1m)-vv.pb(gg.n1m)-aa.phi_a(gg.n1m)+aa.phi_0(gg.n1m))>pp.p_a_reg); % boundary nodes with too high pressure
        % ni2 = gg.n1m((vv.phi(gg.n1m)-vv.pb(gg.n1m)-aa.phi_a(gg.n1m)+vv.phi(gg.n1m))>pp.p_a_reg); % boundary nodes with too high pressure
        if ~isempty(ni1) || ~isempty(ni2)
            if ~isempty(ni1)
                if oo.verb, disp('nevis_timesteps: Removing Dirichlet indices ...'); disp(ni1); end
            end
            if ~isempty(ni2)
                if oo.verb, disp('nevis_timesteps: Adding Dirichlet indices ...'); disp(ni2); end
            end
            vv.nbdy = union(setdiff(gg.nbdy,ni1),ni2);
            gg = nevis_label(gg,vv.nbdy,oo);               % redefine boundary labels
            aa.phi = aa.phi_b(gg.nbdy);                    % boundary conditions
            % nevis_plot_grid(gg,gg.nbdy);  % for maintenance
        end
    end 
    
    disp(['nevis_timesteps: t = ',num2str(t),' [ / ',num2str(t_stop),' ]']);
    ti = ti + 1; 
end

end