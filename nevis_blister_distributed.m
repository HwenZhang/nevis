function [aa,vv] = nevis_blister_distributed(aa,vv,pp,gg,oo)
    % 5/12/2024 HZ
    % This is a module containing blister physics
    % Input: local pressure of the subglacial drainage system (p_w)
    %        input from the lake (aa.Q_lake), 
    % Update:
    %        the volume (vv.Vb) and corresponding radius (vv.Rb)
    %        the flow rate from the blister to the surrounding, distributed across grid points (aa.Q_b)
    % Return:
    %        aa,vv
    
    % include the input from the lake
    if ~isfield(oo,'include_blister'), oo.include_blister = 1; end    % 
    if ~isfield(oo,'blister_distributed'), oo.blister_distributed = 0; end     % whether to distribued the blister 

    %% input from the lake to the blister
    % a.Q_lake is the input to the blister model
    aa.Q_lake = 0*ones(gg.nIJ,1);                                     % set the lake input to 0 everywhere
    if oo.include_blister && isfield(pp,'ni_l') && ~isempty(pp.ni_l)
        aa.Q_lake(pp.ni_l) = pp.lake_input_function(vv.t);            % add input from lake at point pp.ni_l
    end
    
    %% update volume according to mass conservation
    % inflow rate from lake
    Q_in = pp.c42*aa.Q_lake; 
    % outflow rate by Darcy's law
    Q_out = pp.c43*pp.k_blister(vv.hs).*vv.Vb./(vv.Rb).^2 + ...
            pp.c44*pp.k_blister(vv.hs).*vv.Rb.*(aa.phi_0-vv.phi);

    % Q_out = min(pp.c43*pp.k_blister(vv.hs).*vv.Vb./(vv.Rb).^2,...
    %             pp.c48*pp.k_blister(vv.hs).*vv.Rb) + ...
    %         pp.c44*pp.k_blister(vv.hs).*vv.Rb.*(aa.phi_0-vv.phi);

    % the outflow rate should be smaller than the remaining volume/dt,
    % otherwise get negative volume after extraction
    Q_out = min(Q_out, vv.Vb/vv.dt);
    
    % rate of volume variation
    V_t = Q_in - Q_out;
    
    % update the volume, forward euler
    vv.Vb = vv.Vb + vv.dt.*V_t;
    
    % display the volume of the blister
    disp(['The blister volume is ' num2str(max(vv.Vb)), ' V0.']);
    disp(['The ratio $Q_{out}/Q_{in}$ ' num2str(Q_out(pp.ni_l)/Q_in(pp.ni_l))]);
    
    %% update radius according to toughness/exponential decay
    vv.Rb = (vv.Vb>0).*max(vv.Rb, 1/pp.c45*vv.Vb.^0.4); % radius should never decrease unless the volume drops to 0

    % alternatively, consider the plastic effect
    % sol = fzero(@(x) pp.c45*x^3/(x+pp.Ri)^0.5-vv.Vb(pp.ni_l), 1/pp.c45*vv.Vb(pp.ni_l)^0.4);
    % vv.Rb = (vv.Vb>0).*max(vv.Rb, sol); % radius should never decrease unless the volume drops to 0 !!!!!!!!!!

    % vv.Rb = vv.Rb + (vv.Rb==0).*pp.Ri;                  % keep a minimum radius
    % forward Euler
    % vv.Rb = (vv.Vb>0).*1/pp.c45.*(pp.c45*vv.Rb.^(2.5) + (pp.c45.*vv.Rb.^(2.5)<=vv.Vb).*(V_t>0).*vv.dt.*V_t).^(0.4);
    %% distribute Q_{out} to surrounding grid points
    if oo.blister_distributed
        r_b = pp.c47*((gg.nx-gg.nx(pp.ni_l)).^2+(gg.ny-gg.ny(pp.ni_l)).^2).^(0.5);        % distance from grid points to the lake
        r_in = r_b <= vv.Rb(pp.ni_l);                                                     % if the grid point is in the blister
        n_r = sum(r_in(:));
        % n_r = 1;
        Q_out = 1/n_r*r_in*Q_out(pp.ni_l);
        disp(['There are ' num2str(n_r) ' grid points within the blister area.']);
    end
    
    %% set inflow from the blister to the drainage network
    aa.Qb = Q_out./gg.Dx./gg.Dy;                          % dimensionally, aa.Q_b ~ L T^(-1) ~ [E]

end