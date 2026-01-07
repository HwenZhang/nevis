function [tt,vv,info] = nevis_timesteps_ice(t_span,vv,aa,pp,gg,oo)
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
%       17 November 2015 : added possibility of prescribing potential aa.phi_b
%       10 October 2025 : added option to solve for ice velocity using
%           nevis_velocity.m

if ~isfield(pp,'phi_s'), pp.phi_s = -inf; end

% TIMESTEPPING OPTIONS
if ~isfield(oo,'include_ice'), oo.include_ice = 0; end          % update ice velocity
if ~isfield(oo,'cavity_coupling'), oo.cavity_coupling = 0; end  % update sliding speed in cavity opening term
if ~isfield(oo,'melt_coupling'), oo.melt_coupling = 0; end      % update sliding speed in basal melt term
if ~isfield(oo,'change_timestep'), oo.change_timestep = 1; end  % change timestep based on previous iterations/convergence
if ~isfield(oo,'dt_max'), oo.dt_max = 1e0; end                  % maximum allowable timestep
if ~isfield(oo,'dt_min'), oo.dt_min = 1e-6; end                 % mimimum allowable timestep
if ~isfield(oo,'dt_factor'), oo.dt_factor = 2; end              % factor to change timestep by
if ~isfield(oo,'small_iter'), oo.small_iter = 2; end            % increase timestep for this number of iterations or less
if ~isfield(oo,'large_iter'), oo.large_iter = 20; end           % decrease timestep for this number of iterations or more 
if ~isfield(oo,'verb'), oo.verb = 0; end                        % verbose screen output
if ~isfield(oo,'adjust_boundaries'), oo.adjust_boundaries = 1; end  % adjust boundary nodes during timestepping
if ~isfield(oo,'dti_boundaries'), oo.dti_boundaries = 4; end    % frequency of timesteps at which to adjust boundaries
% INPUT OPTIONS
if ~isfield(oo,'input_function'), oo.input_function = 0; end    % use pp.input function(t) for moulin inputs
if ~isfield(oo,'runoff_function'), oo.runoff_function = 0; end  % use pp.runoff function(t) for runoff
if ~isfield(oo,'distributed_input'), oo.distributed_input = 0; end  % use distributed input for runoff even in presence of moulins
% SAVING OPTIONS
if ~isfield(oo,'save_timesteps'), oo.save_timesteps = 0; end    % save solution at timesteps in t_span to directory
% if oo.save_timesteps && ~exist([oo.root,oo.fn],'dir'), mkdir([oo.root,oo.fn]); end % make directory if it doesn't already exist
if ~isfield(oo,'save_pts_all'), oo.save_pts_all = 0; end        % save output at particular nodes pts_ni
if ~isfield(oo,'pts_ni'), oo.pts_ni = []; end                   % nodes at which to save output
if ~isfield(oo,'save_phi_av'), oo.save_phi_av = 0; end          % save average pressure between timesteps

% UPDATE BOUNDARIES
if oo.adjust_boundaries && isfield(vv,'nbdy')         
    gg = nevis_label(gg,vv.nbdy,oo);               % redefine boundary labels
    if ~isfield(aa,'phi_b'), aa.phi_b = max(aa.phi_a,pp.phi_s); end % fill in missin phi_b for older versions
    aa.phi = aa.phi_b(gg.nbdy);                    % boundary conditions
end

% TIMESTEPPING
t_start = t_span(1);
t_stop = t_span(end);
ti_save = 1;
t_save = t_span(ti_save); 
ti = 1;
ti_boundaries = 1;
dt1 = oo.dt;    % suggested timestep
t = t_start;
vv.dt = dt1;
vv.t = t;
tt = struct;
if oo.save_phi_av
    phi_av = vv.phi; % initial average pressure
    phi_min = vv.phi;
    phi_max = vv.phi;
end
[vv2,~,~,~,~,~,~,~,~] = nevis_backbone(inf,vv,vv,aa,pp,gg,oo); % initial expanded variables vv2
while t<t_stop+oo.dt_min

    vv.u(gg.ebdy2) = 50/31536000/(1.9026e-06); % set boundary ice velocities
    % vv.v(gg.fbdy2) = 10/31536000/(1.9026e-06);
    %% solve for ice velocity
    if oo.include_ice
        N = max(0,aa.phi_0-vv.phi); 
        [u,v] = nevis_velocity(aa.H,vv.u,vv.v,N,aa,pp,gg,oo);
        vv.u = u; 
        vv.v = v; 
        vv.U = ((gg.nmeanx(:,gg.es2)*vv.u(gg.es2)).^2+0*(gg.nmeany(:,gg.fs2)*vv.v(gg.fs2)).^2).^(1/2);
        if oo.cavity_coupling
            % update sliding speed for cavity opening term (including Ub
            % in vv will cause this to be used in place of aa.Ub)
            vv.Ub = pp.c65*vv.U;  
        end
        if oo.melt_coupling
            % update frictional heating term (including m in vv will cause
            % this to be used ni place of aa.m)
            [~,~,~,~,~,~,tau_b] = nevis_stresses(aa.H,u,v,N,aa,pp,gg,oo);
            vv.m = pp.c63 + pp.c64*( tau_b.*vv.U );
        end
    
    end

    %% update input
    aa = nevis_inputs(t,aa,vv,pp,gg,oo);

    %% time series of average quantities
    tt(ti).ti = ti;
    tt(ti).t = t;
    tt(ti).m = sum(aa.m(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));  % basal melt, scaled with ps.m*ps.x^2
    tt(ti).E = sum(aa.E(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));  % supraglacial input, scaled with ps.m*ps.x^2
    tt(ti).Q_in = vv2.Q_in;                                % inflow, scaled with ps.Q
    tt(ti).Q_out = vv2.Q_out;                              % outflow, scaled with ps.Q
    tt(ti).Q_outQ = vv2.Q_outQ;                            % channel outflow, scaled with ps.Q
    tt(ti).Q_outq = vv2.Q_outq;                            % sheet outflow, scaled with ps.Q
    tt(ti).phi = mean(vv.phi(gg.ns));                      % mean potential, scaled with ps.phi
    tt(ti).N = mean(aa.phi_0(gg.ns)-vv.phi(gg.ns));        % mean effective pressure, scaled with ps.phi
    tt(ti).hs = sum(vv.hs(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns)); % total cavity sheet volume, scaled with ps.h*ps.x^2
    tt(ti).S = sum(vv.Sx(gg.ein).*(gg.emean(gg.ein,:)*gg.Dx))+...
               sum(vv.Sy(gg.fin).*(gg.fmean(gg.fin,:)*gg.Dy))+...
               sum(vv.Ss(gg.cin).*(gg.cmean(gg.cin,:)*gg.Ds))+...
               sum(vv.Sr(gg.cin).*(gg.cmean(gg.cin,:)*gg.Dr)); % total channel volume, scaled with ps.S*ps.x
    tt(ti).he = sum(vv2.he(gg.ns).*gg.Dx(gg.ns).*gg.Dy(gg.ns));% total elastic sheet volume, scaled with ps.h*ps.x^2
    if oo.include_ice
         tt(ti).U = mean(vv.U);                                % mean ice velocity, scaled with ps.u
         disp([[t*10,60*tt(ti).U]]);
    end
    if oo.save_pts_all
         tt(ti).pts_phi = vv.phi(oo.pts_ni);
         tt(ti).pts_hs = vv.hs(oo.pts_ni);
         tt(ti).pts_he = vv2.he(oo.pts_ni);
         if oo.include_ice
            tt(ti).pts_u = vv.U(oo.pts_ni);
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
        % fn = [oo.rn,oo.fn];
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
    elseif t > t_stop-oo.dt_min, break; 
    end

    %% timestep hydrology
    tic;
    accept = 0;
    decreased = 0;
    increased = 0;
    while ~accept
        dt = dt1; %suggested timestep
        if t + dt > t_save-oo.dt_min, dt = t_save-t; end % adjust timestep if this would take too close to t_save
        disp(['nevis_timesteps: dt = ',num2str(dt),' ...']);
        [vv1,vv2,info] = nevis_timestep(dt,vv,aa,pp,gg,oo);

        %% check success and adjust size of timestep [ taken from hydro_timestep_diag ]
        if ~info.failed 
            accept = 1; 
            t = t + dt;
            vv = vv1;
            vv.dt = dt;
            vv.t = t;
            comp_time = toc;
            disp(['nevis_timesteps: Done [ ',num2str(comp_time),' s, ',num2str(info.iter_new-1),' iterations ]']);
        end
        if oo.change_timestep && info.iter_new-1 <= oo.small_iter && ~info.failed && dt1 < oo.dt_max && ~decreased
            dt1 = min(oo.dt_factor*dt1,oo.dt_max); increased = 1;
            if oo.verb, disp(['nevis_timesteps: Increase suggested timestep to ',num2str(dt1)]); end
            continue;
        elseif oo.change_timestep && ( info.iter_new-1 >= oo.large_iter || info.failed ) && dt1 > oo.dt_min && ~increased
            dt1 = max(dt1/oo.dt_factor,oo.dt_min); decreased = 1;
            if oo.verb, disp(['nevis_timesteps: Decrease suggested timestep to ',num2str(dt1)]); end
            continue;
        end
        if ~accept
            comp_time = toc;
            disp(['nevis_timesteps: Failed [ ',num2str(comp_time),' s, ',num2str(info.iter_new-1),' iterations, dt = ',num2str(dt),' ]']);
            return;
        end

    end
    
    %% time averaging of pressure
    if oo.save_phi_av
        phi_av = ( phi_av*(t-t_span(ti_save-1)-dt)+vv.phi*dt )/(t-t_span(ti_save-1));
        phi_min = min(phi_min,vv.phi);
        phi_max = max(phi_max,vv.phi);
    end
    
    %% check and adjust boundary nodes
    if oo.adjust_boundaries && ti>=ti_boundaries
    ti_boundaries = ti+oo.dti_boundaries;
    ni1 = gg.nbdy(vv2.R_bdy<0); % Dirichlet nodes with inflow
    if ~isfield(gg,'n1m'), gg.n1m = gg.n1; end % boundary nodes adjacent to margin
    ni2 = gg.n1m(vv.phi(gg.n1m)>pp.p_a_reg); % boundary nodes with too high pressure
    ni2 = union(ni2, gg.n1m((vv.pb(gg.n1m)+aa.phi_a(gg.n1m)-aa.phi_0(gg.n1m))>pp.p_a_reg));
    if ~isempty(ni1) || ~isempty(ni2)
        if ~isempty(ni1) 
            if oo.verb, disp('nevis_timesteps: Removing Dirichlet indices ...');
            disp(ni1); end
        end
        if ~isempty(ni2)
            if oo.verb, disp('nevis_timesteps: Adding Dirichlet indices ...');
            disp(ni2); end
        end
        vv.nbdy = union(setdiff(gg.nbdy,ni1),ni2);
        gg = nevis_label(gg,vv.nbdy,oo);               % redefine boundary labels
        aa.phi = aa.phi_b(gg.nbdy);                    % boundary conditions
%        nevis_plot_grid(gg,gg.nbdy);  % for maintenance
    end
    end 
    
    disp(['nevis_timesteps: t = ',num2str(t),' [ / ',num2str(t_stop),' ]']);
    ti = ti + 1; 
    
end
end
