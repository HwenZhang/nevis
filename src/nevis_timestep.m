function [vv1,vv2,info] = nevis_timestep(dt,vv,aa,pp,gg,oo)
% advance nevis variables a timestep dt using Newton iteration
% Inputs:
%   dt      timestep
%   vv      struct containing current solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv1     struct containing solution variables
%   vv2     struct containing expanded solution variables
%   info    information about computation
%
% IJH 13 August 2014 : largely taken from hydro_timestep_diag
% HZ 06 June 2024 : adding ice sheet dynamics
num = 10;
% ITERATION OPTIONS
if ~isfield(oo,'Tol_F'), oo.Tol_F = 1e-8; end                       % tolerance on Newton iteration
if ~isfield(oo,'check_Fs'), oo.check_Fs = 0; end                    % check tolerances seprately for residuals 1-6
if ~isfield(oo,'Tol_Fs'), oo.Tol_Fs = oo.Tol_F*ones(1,num); end       % tolerances on Newoton iteration [ for when using check_Fs ]
if length(oo.Tol_Fs)<num
    oo.Tol_Fs = [oo.Tol_Fs oo.Tol_Fs(end)*ones(1,num-length(oo.Tol_Fs))];
end 
if ~isfield(oo,'max_iter_new'), oo.max_iter_new = 50; end           % maximum number of Newton iterations
if ~isfield(oo,'step_new'), oo.step_new = 1; end                    % step size for Newton iteration
if ~isfield(oo,'max2_iter_new'), oo.max2_iter_new = 20; end         % maximum number of Newton iterations before step2 is used
if ~isfield(oo,'step2_new'), oo.step2_new = 0.5*oo.step_new; end    % step2 size for Newton iteration
if ~isfield(oo,'fac2_new'), oo.fac2_new = 1; end                    % factor for norm reduction above which step2 is used
if ~isfield(oo,'step_ice'), oo.step_ice = 1.0; end                  % relaxation factor for ice velocity (R9,R10) Newton step

% ITERATIVE SOLVER OPTIONS (use GMRES+ILU instead of UMFPACK direct solve to avoid memory leaks)
if ~isfield(oo,'iterative_solver'), oo.iterative_solver = 0; end     % 0 = direct, 1 = GMRES with ILU preconditioner
if ~isfield(oo,'direct_solver'), oo.direct_solver = 'lu'; end        % 'lu' = UMFPACK via explicit lu() (factors freed on exit)
                                                                      % 'qr' = SuiteSparseQR via qr() (avoids UMFPACK entirely)
if ~isfield(oo,'gmres_restart'), oo.gmres_restart = 50; end          % GMRES restart parameter
if ~isfield(oo,'gmres_tol'), oo.gmres_tol = 1e-6; end               % GMRES relative tolerance (on true residual)
if ~isfield(oo,'gmres_safety'), oo.gmres_safety = 10; end            % safety factor: internal tol = gmres_tol/gmres_safety
                                                                      %   compensates for left-preconditioning discrepancy
                                                                      %   (GMRES converges on preconditioned residual,
                                                                      %    but relres reports true residual which is larger)
if ~isfield(oo,'gmres_tol_accept'), oo.gmres_tol_accept = 1e-5; end  % accept GMRES result if relres below this (inexact Newton)
if ~isfield(oo,'gmres_maxit'), oo.gmres_maxit = 300; end            % GMRES maximum iterations
if ~isfield(oo,'ilu_type'), oo.ilu_type = 'ilutp'; end              % ILU type: 'nofill', 'crout', 'ilutp'
if ~isfield(oo,'ilu_droptol'), oo.ilu_droptol = 1e-3; end           % ILU drop tolerance

% DIAGNOSTIC OPTIONS
if ~isfield(oo,'plot_residual'), oo.plot_residual = 0; end          % plot residuals at each iteration
if ~isfield(oo,'display_residual'), oo.display_residual = 0; end    % display maximum of residuals at each iteration
if ~isfield(oo,'verb'), oo.verb = 0; end                            % verbose screen output
% [could add other diagnostic options here]  

% ITERATION
vv0 = vv;
vv2 = struct();  % initialise so output is always assigned (even on early abort)
%% struct for solution info
info = struct;
info.failed = 0;            % failure indicator
info.residual_time = 0;     % time spent building residual
info.jacob_time = 0;        % time spent building Jacobian
info.start_time = tic;      % time starter

%% iteration options
Tol_F = oo.Tol_F;                    
check_Fs = oo.check_Fs;              
Tols_F = oo.Tol_Fs;                   
max_iter_new = oo.max_iter_new;     
norm_Fs = NaN*ones(max_iter_new+1,1);
norms_Fs = NaN*ones(max_iter_new+1,num);
for iter_new = 1:max_iter_new+1

    %% safeguard: check solution variables BEFORE any backbone call
    % (backbone internals can overflow/crash if inputs are extreme)
    vars_ok = isfinite_bounded(vv.phi) && isfinite_bounded(vv.hs) && ...
              isfinite_bounded(vv.Sx) && isfinite_bounded(vv.Sy) && ...
              isfinite_bounded(vv.Ss) && isfinite_bounded(vv.Sr);
    if oo.include_blister
        vars_ok = vars_ok && isfinite_bounded(vv.hb) && isfinite_bounded(vv.pb);
    end
    if oo.include_ice
        vars_ok = vars_ok && isfinite_bounded(vv.u) && isfinite_bounded(vv.v);
    end
    % if ~vars_ok
    %     fprintf('**Aborting (iter %d): solution variables contain non-finite or extreme values\n', iter_new);
    %     fprintf('  max|phi|=%.2e, max|hs|=%.2e, max|hb|=%.2e, max|pb|=%.2e\n', ...
    %         max(abs(vv.phi)), max(abs(vv.hs)), max(abs(vv.hb)), max(abs(vv.pb)));
    %     info.failed = 1; break;
    % end

    %% evaluate residual
    tstart = tic;
    oo.evaluate_variables = 0; 
    oo.evaluate_residual = 1; 
    oo.evaluate_jacobian = 0; 
    [vv2,F,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,~] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo);
    info.residual_time = info.residual_time + toc(tstart);

    %% update Xi if not including Xi in iterations
    if iter_new==1 && oo.noXi, aa.Xi = vv2.Xi; oo.includeXi = 0; end % set oo.includeXi = 0 so Xi not evaluated at next iteration

    %% check convergence
    norm_F = norm(F,inf);
    norm_Fs(iter_new) = norm_F;
    norm_F1 = norm(F1,inf); 
    % norm_F2 = norm(F2,inf); 
    norm_F2 = norm(F2,2)/length(F2);
    norm_F3 = norm(F3,inf); 
    norm_F4 = norm(F4,inf);
    norm_F5 = norm(F5,inf);
    norm_F6 = norm(F6,inf);
    norm_F7 = norm(F7,inf);
    norm_F8 = norm(F8,inf);
    if oo.include_ice
        norm_F9 = norm(F9,inf);
        norm_F10 = norm(F10,inf);
    else
        norm_F9 = 0;
        norm_F10 = 0;
    end

    norms_F = [norm_F1 norm_F2 norm_F3 norm_F4 norm_F5 norm_F6 norm_F7 norm_F8 norm_F9 norm_F10];

    norms_Fs(iter_new,:) = norms_F;
    if oo.plot_residual
        figure(1); clf;
        plot(F1,'k.'); hold on;
        plot(F2,'b.');
        plot(F3,'r.');
        plot(F4,'g.');
        plot(F5,'m.');
        plot(F6,'y.');
        plot(F7,'b-');
        plot(F8,'r-');
        plot(F9,'g-');
        plot(F10,'m-');
        shg;
    end

    % if oo.plot_residual
    %     figure(2);clf
    %     scatter(gg.nx(gg.ns_blister),gg.ny(gg.ns_blister),20,abs(F7),'filled'); colorbar; title('Residual F7 on sheet nodes');shg
    % end
    % figure(2);clf
    % scatter(gg.nx(gg.ns_blister),gg.ny(gg.ns_blister),20,abs(F7),'filled'); colorbar; title('Residual F7 on sheet nodes');shg
    % 
    % figure(3);clf
    % scatter(gg.nx(gg.nin_blister),gg.ny(gg.nin_blister),20,abs(F8),'filled'); colorbar; title('Residual F8 on sheet nodes');shg

    if oo.display_residual
        [m1,i1] = max(abs(F1)); 
        [m2,i2] = max(abs(F2)); 
        [m3,i3] = max(abs(F3)); 
        [m4,i4] = max(abs(F4)); 
        [m5,i5] = max(abs(F5)); 
        [m6,i6] = max(abs(F6));
        [m7,i7] = max(abs(F7));
        [m8,i8] = max(abs(F8));
        if oo.include_ice
            [m9,i9] = max(abs(F9));
            [m10,i10] = max(abs(F10));
        else
            m9 = 0; i9 = 0;
            m10 = 0; i10 = 0;
        end
        disp([ iter_new m1 m2 m3 m4 m5 m6 m7 m8 m9 m10]);
        disp([ 0 i1 i2 i3 i4 i5 i6 i7 i8 i9 i10]); 
    end
    if oo.no_channels && oo.no_sheet, iFs = 2; 
    elseif oo.no_channels, iFs = [1 2]; 
    elseif oo.no_sheet, iFs = [2 3 4 5 6];  
    else
        iFs = [1 2 3 4 5 6];
        if oo.include_blister
            if oo.include_pressure
                iFs = [1 2 3 4 5 6 7 8];  
                if oo.include_ice
                    iFs = [1 2 3 4 5 6 7 8 9 10];
                end
            else
                iFs = [1 2 3 4 5 6 7]; 
                if oo.include_ice
                    iFs = [1 2 3 4 5 6 7 9 10];
                end 
            end
        end        
    end
    if check_Fs && all(norms_F(iFs)<=Tols_F(iFs)), 
        if oo.verb, disp(['  norms_F less than tolerance after ',num2str(iter_new-1),' Newton iterations']); end; 
        vv1 = vv;
        info.failed = 0; break, 
    elseif norm_F <= Tol_F && iter_new > 1, 
        if oo.verb, disp(['  norm_F less than tolerance after ',num2str(iter_new-1),' Newton iterations']); end; 
        vv1 = vv;
        info.failed = 0; break, 
    end
    if iter_new > max_iter_new, 
        disp(['**norm_F = ',num2str(norm_F),' exceeds tolerance after ',num2str(iter_new-1),' Newton iterations']); 
        info.failed = 1; break, 
    end

    %% check for divergence — abort early instead of waiting for NaN/crash
    % 1) abort if F contains NaN/Inf
    if any(~isfinite(F))
        fprintf('**Aborting (iter %d): F contains non-finite values (NaN: %d, Inf: %d)\n', ...
            iter_new, sum(isnan(F)), sum(isinf(F)));
        info.failed = 1; break;
    end
    % 2) abort if norm_F is astronomically large (prevents overflow in Jacobian)
    if norm_F > 1e50
        fprintf('**Aborting (iter %d): norm_F = %.2e is dangerously large\n', iter_new, norm_F);
        info.failed = 1; break;
    end
    % 3) abort if norm_F has grown for 3+ consecutive iterations
    if iter_new >= 4
        n_consec_increase = 0;
        for kk = iter_new:-1:2
            if norm_Fs(kk) > norm_Fs(kk-1)
                n_consec_increase = n_consec_increase + 1;
            else
                break;
            end
        end
        if n_consec_increase >= 3
            fprintf('**Aborting: Newton diverging — norm_F increased %d consecutive iterations (%.2e -> %.2e)\n', ...
                n_consec_increase, norm_Fs(iter_new-n_consec_increase), norm_F);
            info.failed = 1; break;
        end
    end
    % 4) abort if norm_F has grown by more than 1e8 relative to initial
    if iter_new >= 2 && norm_F > 1e8 * norm_Fs(1)
        fprintf('**Aborting: Newton diverging — norm_F = %.2e exceeds 1e8 × initial (%.2e)\n', ...
            norm_F, norm_Fs(1));
        info.failed = 1; break;
    end

    %% calculate Jacobian         
    tstart = tic;
    oo.evaluate_variables = 0; oo.evaluate_residual = 0; oo.evaluate_jacobian = 1; 
    [~,~,~,~,~,~,~,~,~,~,~,~,J] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo);
    info.jacob_time = info.jacob_time + toc(tstart);

    %% Newton step [ add line search here ? ]
    if oo.no_channels && oo.no_sheet
        X = vv.phi(gg.nin);
    elseif oo.no_channels
        X = [vv.hs(gg.ns); vv.phi(gg.nin)];
    elseif oo.no_sheet
        X = [vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin)];
    else 
        % include both channel and sheet
        % check if include volume and radius
        if oo.include_blister
            if oo.include_pressure
                X = [vv.hs(gg.ns); vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin); vv.hb(gg.ns_blister); vv.pb(gg.nin_blister)];
                if oo.include_ice
                    X = [vv.hs(gg.ns); vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin); vv.hb(gg.ns_blister); vv.pb(gg.nin_blister); vv.u(gg.ein2); vv.v(gg.fin2)];
                end
            else
                X = [vv.hs(gg.ns); vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin); vv.hb(gg.ns_blister)];
                if oo.include_ice
                    X = [vv.hs(gg.ns); vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin); vv.hb(gg.ns_blister); vv.u(gg.ein2); vv.v(gg.fin2)];
                end
            end
        else
            X = [vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin)];
        end
    end

    % if condest(J) >=1e20, disp(' Aborting Newton step : J is nearly singular'); info.failed = 1; break; end

    % Check J for NaN/Inf/extreme values before linear solve
    % (extreme finite values in J cause UMFPACK to crash with SIGTRAP)
    [~,~,Jv] = find(J);
    max_J = max(abs(Jv));
    if any(~isfinite(Jv))
        fprintf('  Aborting (iter %d): J contains NaN/Inf (max|J|=%.2e)\n', iter_new, max_J);
        info.failed = 1; break;
    end
    if max_J > 1e100
        fprintf('  Aborting (iter %d): J has extreme values (max|J|=%.2e) — would crash UMFPACK\n', iter_new, max_J);
        info.failed = 1; break;
    end
    if oo.iterative_solver
        %% iterative solve: ILU-preconditioned GMRES (avoids UMFPACK memory leak)
        try
            iluOpts = struct('type', oo.ilu_type, 'droptol', oo.ilu_droptol);
            [M1, M2] = ilu(J, iluOpts);
        catch ME_ilu
            % ILU can fail on near-singular matrices; fall back to explicit LU
            if oo.verb, fprintf('  ILU failed (iter %d): %s — falling back to direct (%s)\n', iter_new, ME_ilu.message, oo.direct_solver); end
            [dX, lu_ok] = solve_direct(J, F, oo.direct_solver);
            if ~lu_ok, fprintf('**Aborting (iter %d): direct solve failed\n', iter_new); info.failed = 1; break; end
            M1 = [];  % skip GMRES below
        end
        if ~isempty(M1)
            warnState = warning('off', 'MATLAB:gmres:tooSmallTolerance');
            [dX, gmres_flag, gmres_relres, gmres_iter] = gmres(J, -F, oo.gmres_restart, oo.gmres_tol/oo.gmres_safety, oo.gmres_maxit, M1, M2);
            warning(warnState);
            if gmres_flag == 0 || gmres_relres <= oo.gmres_tol
                % converged, or true relres meets target — use dX as-is
                if oo.verb && gmres_flag ~= 0, fprintf('  GMRES flag=%d, relres=%.2e <= gmres_tol=%.1e — accepting\n', gmres_flag, gmres_relres, oo.gmres_tol); end
            elseif gmres_relres <= oo.gmres_tol_accept
                % not fully converged but good enough for inexact Newton
                if oo.verb, fprintf('  GMRES flag=%d, relres=%.2e — accepting (tol_accept=%.1e)\n', gmres_flag, gmres_relres, oo.gmres_tol_accept); end
            else
                % poor approximation — fall back to direct solve
                if oo.verb, fprintf('  GMRES flag=%d, relres=%.2e — falling back to direct (%s)\n', gmres_flag, gmres_relres, oo.direct_solver); end
                [dX, lu_ok] = solve_direct(J, F, oo.direct_solver);
                if ~lu_ok, fprintf('**Aborting (iter %d): direct solve failed\n', iter_new); info.failed = 1; break; end
            end
        end
    else
        %% direct solve (avoids UMFPACK internal caching / memory leak)
        [dX, lu_ok] = solve_direct(J, F, oo.direct_solver);
        if ~lu_ok, fprintf('**Aborting (iter %d): direct solve failed\n', iter_new); info.failed = 1; break; end
    end

    if any(isnan(dX)), disp(' Aborting Newton step : NaN returned'); info.failed = 1; break; end;
    step = oo.step_new;
    if iter_new>2 && norm_Fs(iter_new)/norm_Fs(iter_new-1)>oo.fac2_new, step = oo.step2_new; end % reduce size of step if previous iterations are slow
    if iter_new>oo.max2_iter_new, step = oo.step2_new; end % reduce size of step if more than max2_iter_new iterations
    % apply per-component relaxation: step_ice for u,v (R9,R10)
    step_vec = step * ones(size(dX));
    % adaptive step_ice: ramp up from step_ice to 1 as iterations progress
    if oo.include_ice
        n_ice = length(gg.ein2) + length(gg.fin2);
        if iter_new <= 3
            % first few iterations: use the (small) prescribed step_ice
            step_ice_iter = oo.step_ice;
        elseif iter_new <= 6
            % intermediate: blend toward 1
            step_ice_iter = min(1, oo.step_ice + (1 - oo.step_ice) * (iter_new - 3) / 3);
        else
            % later iterations: full step
            step_ice_iter = 1;
        end
        step_vec(end-n_ice+1:end) = step_ice_iter;
    end
    X = X + step_vec.*dX;
    
    vv.phi(gg.nbdy) = aa.phi; % boundary conditions
    
    temp1 = 0;
    if ~oo.no_sheet
        temp2 = length(gg.ns); 
        vv.hs(gg.ns) = X(temp1+(1:temp2)); 
        temp1=temp1+temp2;
    end
    
    temp2 = length(gg.nin); 
    vv.phi(gg.nin) = X(temp1+(1:temp2)); 
    temp1=temp1+temp2;

    if ~oo.no_channels
        temp2 = length(gg.ein); vv.Sx(gg.ein) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        temp2 = length(gg.fin); vv.Sy(gg.fin) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        temp2 = length(gg.cin); vv.Ss(gg.cin) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        temp2 = length(gg.cin); vv.Sr(gg.cin) = X(temp1+(1:temp2)); temp1=temp1+temp2;
    end
    % update the blister volume and radius
    if oo.include_blister
        temp2 = length(gg.ns_blister); vv.hb(gg.ns_blister) = X(temp1+(1:temp2)); temp1=temp1+temp2; % vv.hb=max(vv.hb,0);
        if oo.include_pressure
            temp2 = length(gg.nin_blister); vv.pb(gg.nin_blister) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        end
    end
    % update ice sheet velocity
    if oo.include_ice
        temp2 = length(gg.ein2); vv.u(gg.ein2) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        temp2 = length(gg.fin2); vv.v(gg.fin2) = X(temp1+(1:temp2));
    end
end

%% outputs
info.iter_new = iter_new;
info.norm_Fs = norm_Fs;
info.norms_Fs = norms_Fs;
info.stop_time = toc(info.start_time);

% check whether SSA residuals (F9,F10) are above tolerance
info.ice_large = false;
if info.failed && oo.include_ice
    last_norms = norms_Fs(min(iter_new, size(norms_Fs,1)), :);
    info.ice_large = any(last_norms([9 10]) > 10*Tols_F([9 10]));
end

if info.failed
    vv1 = vv0;
    % ensure vv2 is populated even on failure (caller may inspect it)
    if isempty(fieldnames(vv2))
        try
            oo.evaluate_variables = 1; oo.evaluate_residual = 0; oo.evaluate_jacobian = 0;
            [vv2,~,~,~,~,~,~,~,~,~,~,~,~] = nevis_backbone(dt,vv0,vv0,aa,pp,gg,oo);
        catch
            % leave vv2 as empty struct if backbone also fails
        end
    end
end
    
end

function ok = isfinite_bounded(x)
% Check that vector x has no NaN, no Inf, and no extreme values (>1e100)
% that would cause overflow in quadratic/power-law expressions inside backbone
    ok = all(isfinite(x)) && (isempty(x) || max(abs(x)) < 1e50);
end

function [dX, ok] = solve_direct(J, F, method)
% Direct solve with explicit factorization; factors are local variables
% and freed on function return (avoids UMFPACK internal caching leaks).
%   method = 'lu' : sparse LU via UMFPACK  (fast, standard)
%   method = 'qr' : sparse QR via SPQR     (no UMFPACK, more fill-in)
    ok = true;
    try
        if strcmp(method, 'qr')
            [Q, R, E] = qr(J);
            dX = -E * (R \ (Q' * F));
        else  % 'lu' (default)
            [L, U, P, Q] = lu(J);
            dX = -Q * (U \ (L \ (P * F)));
        end
    catch ME
        fprintf('  solve_direct (%s): %s\n', method, ME.message);
        dX = NaN(size(F));
        ok = false;
    end
end
