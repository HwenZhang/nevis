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
num = 8;
% ITERATION OPTIONS
if ~isfield(oo,'Tol_F'), oo.Tol_F = 1e-8; end                       % tolerance on Newton iteration
if ~isfield(oo,'check_Fs'), oo.check_Fs = 0; end                    % check tolerances seprately for residuals 1-6
if ~isfield(oo,'Tol_Fs'), oo.Tol_Fs = oo.Tol_F*ones(1,num); end       % tolerances on Newoton iteration [ for when using check_Fs ]
if length(oo.Tol_Fs)<num
    oo.Tol_Fs = [oo.Tol_Fs oo.Tol_Fs(end)*ones(1,num-length(oo.Tol_Fs))];
end 
if ~isfield(oo,'max_iter_new'), oo.max_iter_new = 20; end            % maximum number of Newton iterations
if ~isfield(oo,'step_new'), oo.step_new = 1; end                    % step size for Newton iteration
if ~isfield(oo,'max2_iter_new'), oo.max2_iter_new = 20; end         % maximum number of Newton iterations before step2 is used
if ~isfield(oo,'step2_new'), oo.step2_new = 0.5*oo.step_new; end    % step2 size for Newton iteration
if ~isfield(oo,'fac2_new'), oo.fac2_new = 1; end                    % factor for norm reduction above which step2 is used

% DIAGNOSTIC OPTIONS
if ~isfield(oo,'plot_residual'), oo.plot_residual = 0; end          % plot residuals at each iteration
if ~isfield(oo,'display_residual'), oo.display_residual = 0; end    % display maximum of residuals at each iteration
if ~isfield(oo,'verb'), oo.verb = 0; end                            % verbose screen output
% [could add other diagnostic options here]  

% ITERATION
vv0 = vv;
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

    %% evaluate residual
    tstart = tic;
    oo.evaluate_variables = 0; 
    oo.evaluate_residual = 1; 
    oo.evaluate_jacobian = 0; 
    [vv2,F,F1,F2,F3,F4,F5,F6,F7,F8,~] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo);
    info.residual_time = info.residual_time + toc(tstart);

    %% update Xi if not including Xi in iterations
    if iter_new==1 && oo.noXi, aa.Xi = vv2.Xi; oo.includeXi = 0; end % set oo.includeXi = 0 so Xi not evaluated at next iteration

    %% check convergence
    norm_F = norm(F,inf);
    norm_Fs(iter_new) = norm_F;
    norm_F1 = norm(F1,inf); 
    norm_F2 = norm(F2,inf); norm_F2 = norm(F2,2)/length(F2);
    norm_F3 = norm(F3,inf); 
    norm_F4 = norm(F4,inf);
    norm_F5 = norm(F5,inf);
    norm_F6 = norm(F6,inf);
    norm_F7 = norm(F7,inf);
    norm_F8 = norm(F8,inf);
    norms_F = [norm_F1 norm_F2 norm_F3 norm_F4 norm_F5 norm_F6 norm_F7 norm_F8];
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
        shg;
    end
    if oo.display_residual
        [m1,i1] = max(abs(F1)); 
        [m2,i2] = max(abs(F2)); 
        [m3,i3] = max(abs(F3)); 
        [m4,i4] = max(abs(F4)); 
        [m5,i5] = max(abs(F5)); 
        [m6,i6] = max(abs(F6));
        [m7,i7] = max(abs(F7));
        [m8,i8] = max(abs(F8));
        disp([ iter_new m1 m2 m3 m4 m5 m6 m7 m8]);
        disp([ 0 i1 i2 i3 i4 i5 i6 i7 i8]); 
    end
    if oo.no_channels && oo.no_sheet, iFs = 2; 
    elseif oo.no_channels, iFs = [1 2]; 
    elseif oo.no_sheet, iFs = [2 3 4 5 6];  
    else
        iFs = [1 2 3 4 5 6];
        if oo.include_blister
            if oo.include_pressure
                iFs = [1 2 3 4 5 6 7 8];  
            else
                iFs = [1 2 3 4 5 6 7];  
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

    %% calculate Jacobian         
    tstart = tic;
    oo.evaluate_variables = 0; oo.evaluate_residual = 0; oo.evaluate_jacobian = 1; 
    [~,~,~,~,~,~,~,~,~,~,J] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo);
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
                X = [vv.hs(gg.ns); vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin); vv.hb(gg.ns); vv.pb(gg.nin)];
            else
                X = [vv.hs(gg.ns); vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin); vv.hb(gg.ns)];
            end
        else
            X = [vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin)];
        end
    end

    % if condest(J) >=1e20, disp(' Aborting Newton step : J is nearly singular'); info.failed = 1; break; end
    dX = -J\F;
    if any(isnan(dX)), disp(' Aborting Newton step : NaN returned'); info.failed = 1; break; end;
    step = oo.step_new;
    if iter_new>2 && norm_Fs(iter_new)/norm_Fs(iter_new-1)>oo.fac2_new, step = oo.step2_new; end % reduce size of step if previous iterations are slow
    if iter_new>oo.max2_iter_new, step = oo.step2_new; end % reduce size of step if more than max2_iter_new iterations
    X = X + step*dX;
    
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
        temp2 = length(gg.ns); vv.hb(gg.ns) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        if oo.include_pressure
            temp2 = length(gg.nin); vv.pb(gg.nin) = X(temp1+(1:temp2)); % temp1=temp1+temp2;
        end
    end
end

%% outputs
info.iter_new = iter_new;
info.norm_Fs = norm_Fs;
info.norms_Fs = norms_Fs;
info.stop_time = toc(info.start_time);

if info.failed
    vv1 = vv0;
end
    
end
