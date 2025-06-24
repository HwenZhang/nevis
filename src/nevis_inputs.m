function aa = nevis_inputs(t,aa,vv,pp,gg,oo)
% aa = nevis_inputs(t,aa,pp,gg,oo)
% update inputs in aa for time t
%
% 21 August 2014
% 11 Dec 2014 : note that if runoff_function includes interp2, it may not 
% work on earlier matlab versions if interpolating NaN values
% 12 November 2015 : edited to use oo.distributed_input option whether
%    oo.input_function is used or not (rather than only if it's not) by
%    reordering

if nargin<5, oo = struct; end
% use pp.input function(t) for moulin inputs
if ~isfield(oo,'surface_runoff'), oo.surface_runoff = 1; end 
if ~isfield(oo,'RACMO_runoff'), oo.RACMO_runoff = 0; end
if ~isfield(oo,'distributed_input'), oo.distributed_input = 0; end

% if include blister and radius
if ~isfield(oo,'include_lake'), oo.include_lake = 1; end

%% runoff function at the surface & input directly to the subglacial hydrological system
if oo.surface_runoff
    % whether to use RACMO surface runoff data
    if oo.RACMO_runoff 
        r = pp.runoff_function(t); 
        r(gg.nout) = 0;
    else
        r = (1-pp.E_amp*cos(2*pi*t/pp.td)).*max(pp.meltE(t)-pp.E_lapse*aa.s,0); 
        r(gg.nout) = 0;      
    end
    % whether to distribute input to the whole domain or descrete input to moulins
    if isfield(pp,'ni_m') && ~isempty(pp.ni_m) && ~oo.distributed_input
        aa.E = 0*ones(gg.nIJ,1);
        % net input to each moulin
        aa.E(pp.ni_m) = (pp.sum_m*(r.*gg.Dx.*gg.Dy))./gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
    else
        aa.E = r; % the surface runoff is distributed across the domain
    end
else
    % Alternatively, prescribe moulin input with a function 
    aa.E(pp.ni_m) = (pp.input_function(t))./gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
end

%% Prescribe blister input
if oo.include_lake && isfield(pp,'ni_l') && ~isempty(pp.ni_l)
    % use the following function to precribe the input from each lake to the local blister
    if oo.input_gaussian
        pp.lake_input_function = @(t) (pp.V_l)./(sqrt(2*pi)*pp.t_duration)...
                        .*exp(-0.5./pp.t_duration.^2.*(t-pp.t_drainage).^2);
                        % .*(t>=pp.t_drainage-5*pp.t_duration).*(t<=pp.t_drainage+5*pp.t_duration);
    else
        pp.lake_input_function = @(t) (pp.V_l)./pp.t_duration...
                        .*(t>=pp.t_drainage-0.5*pp.t_duration).*(t<=pp.t_drainage+0.5*pp.t_duration);
    end
    aa.Qb_in(pp.ni_l) = pp.lake_input_function(t)./gg.Dx(pp.ni_l)./gg.Dy(pp.ni_l);      % add input from lake at point pp.ni_l
end

%% boundary input
if isfield(pp,'Q_in')
    aa.Qx = 0*ones(length(gg.ebdy),1);
    aa.Qx(pp.ei_in) = pp.Q_in(t);
end

end