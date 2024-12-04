function [aa,vv] = nevis_blister_draft(aa,vv,pp,gg,oo)
% 24/11/2024 HZ
% This is a module containing blister physics
% Input: local pressure of the subglacial drainage system (p_w)
%        input from the lake (aa.Q_lake), 
% Update:
%        the volume (vv.Vb) and corresponding radius (vv.Rb)
%        the flow rate from the blister to the surrounding (aa.Q_b)
% Return:
%        aa,vv

% include the input from the lake
if ~isfield(oo,'lake_input'), oo.lake_input = 1; end     % This is not in function yet!

%% input from the lake to the blister
% a.Q_lake is the input to the blister model
if oo.lake_input && isfield(pp,'ni_l') && ~isempty(pp.ni_l)
    aa.Q_lake = 0*ones(gg.nIJ,1);                        % set the lake input to 0 everywhere
    aa.Q_lake(pp.ni_l) = pp.lake_input_function(vv.t);   % add input from lake at point pp.ni_l
end

%% update volume according to mass conservation
% inflow rate from lake
Q_in = pp.c42*aa.Q_lake; 
% outflow rate by Darcy's law
Q_out = pp.c43*pp.k_blister(vv.hs).*vv.Vb./(vv.Rb).^2 + ...
        pp.c44*pp.k_blister(vv.hs).*vv.Rb.*(aa.phi_0-vv.phi);
% Q_out = Q_in;
% Q_out = zeros(gg.nIJ,1);

% the outflow rate should be smaller than the remaining volume/dt,
% otherwise get negative volume after extraction
Q_out = min(Q_out, vv.Vb/vv.dt);

% rate of volume variation
V_t = Q_in - Q_out;

% update the volume, forward euler
vv.Vb = vv.Vb + vv.dt.*V_t;

% the volume should be non-negative [optional]
% vv.Vb(vv.Vb<0)=0;

% display the volume of the blister
disp(['The blister volume is ' num2str(max(vv.Vb)), ' V0.']);
disp(['The ratio $Q_{out}/Q_{in}$ ' num2str(Q_out(pp.ni_l)/Q_in(pp.ni_l))]);

%% update radius according to toughness/exponential decay
vv.Rb = (vv.Vb>0).*max(vv.Rb, 1/pp.c45*vv.Vb.^0.4); % radius should never decrease unless the volume drops to 0
% forward Euler
% vv.Rb = (vv.Vb>0).*1/pp.c45.*(pp.c45*vv.Rb.^(2.5) + (pp.c45.*vv.Rb.^(2.5)<=vv.Vb).*(V_t>0).*vv.dt.*V_t).^(0.4);
aa.Qb = Q_out./gg.Dx(pp.ni_l)./gg.Dy(pp.ni_l);
% dimensionally, aa.Q_b ~ L T^(-1) ~ [E]

end