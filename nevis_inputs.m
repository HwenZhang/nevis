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
if ~isfield(oo,'input_function'), oo.input_function = 0; end 
% use pp.runoff function(t) for runoff
if ~isfield(oo,'runoff_function'), oo.runoff_function = 0; end
% if include blister and radius
if ~isfield(oo,'include_blister'), oo.include_blister = 1; end             % 
if ~isfield(oo,'blister_distributed'), oo.blister_distributed = 1; end     % whether to distribued the blister 
if ~isfield(oo,'include_radius'), oo.include_radius = 1; end               % 

%% runoff function at the surface & input directly to the subglacial hydrological system
% set up the runoff function at each node
if oo.runoff_function
    r = pp.runoff_function(t); 
    r(gg.nout) = 0;
else
    r = (1-pp.E_amp*cos(2*pi*t/pp.td)).*max(pp.meltE(t)-pp.E_lapse*aa.s,0); 
    r(gg.nout) = 0;
end

% sum surface input into moulins and lakes
if isfield(pp,'ni_m') && ~isempty(pp.ni_m) && ~oo.distributed_input
    aa.E = 0*ones(gg.nIJ,1);
    % net input to each moulin
    aa.E(pp.ni_m) = (pp.sum_m*(r.*gg.Dx.*gg.Dy))./gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
else 
    aa.E = r; % the surface runoff is distributed across the domain
end

% Alternatively, prescribe moulin input with a function
if oo.input_function
    aa.E(pp.ni_m) = (pp.input_function(t))./gg.Dx(pp.ni_m)./gg.Dy(pp.ni_m);
else

end

%% Prescribe blister input
% aa.Q_lake is the input to the blister model
Q_lake = 0*ones(gg.nIJ,1);                            % set the lake input to 0 everywhere
if oo.include_blister && isfield(pp,'ni_l') && ~isempty(pp.ni_l)
    % use the following function to precribe the input from each lake to the local blister
    pp.lake_input_function = @(t) (pp.V_l)./(sqrt(2*pi)*pp.t_duration)...
                    .*exp(-0.5./pp.t_duration.^2.*(t-pp.t_drainage).^2)...
                    .*(t>=pp.t_drainage-5*pp.t_duration).*(t<=pp.t_drainage+5*pp.t_duration);

    Q_lake(pp.ni_l) = pp.lake_input_function(t);      % add input from lake at point pp.ni_l
end

%% Update blister radius if not calculated in the solver
if ~oo.include_radius
    vv.Rb = max(vv.Rb, (1/pp.c45*vv.Vb).^0.4); % radius should never decrease unless the volume drops to 0
    disp(['The blister radius is ' num2str(max(vv.Rb)), ' R0.']);
    disp(['The blister volume is ' num2str(max(vv.Vb)), ' V0.']);
end

%% distribute Q_{out} to surrounding grid points
% if oo.blister_distributed  
%     for i_l=pp.ni_l
%         r_b = pp.c47*((gg.nx-gg.nx(i_l)).^2+(gg.ny-gg.ny(i_l)).^2).^(0.5);        % distance from grid points to the lake
%         r_in = r_b <= vv.Rb(i_l);                                                     % if the grid point is in the blister
%         n_r = sum(r_in(:));
%         r_in(i_l) = r_in(i_l) - n_r;
%         Q_lake = Q_lake + 1/n_r*r_in*Q_lake(i_l);
%         disp(['There are ' num2str(n_r) ' grid points within the blister area.']);
%     end
% end
aa.Qb_in = Q_lake;
if aa.Qb_in(pp.ni_l)>0
    disp("Inflow to the blister.");
end
%% boundary input
% boundary input
if isfield(pp,'Q_in')
    aa.Qx = 0*ones(length(gg.ebdy),1);
    aa.Qx(pp.ei_in) = pp.Q_in(t);
end

end