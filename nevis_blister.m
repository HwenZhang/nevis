function [aa,vv] = nevis_blister(aa,vv,pp,gg,oo)
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
if ~isfield(oo,'lake_input'), oo.lake_input = 1; end

%% input from the lake to the blister
% a.Q_lake is the input to the blister model
if oo.lake_input && isfield(pp,'ni_l') && ~isempty(pp.ni_l)
    aa.Q_lake = 0*ones(gg.nIJ,1);
    % prescribe input to the lake
    aa.Q_lake(pp.ni_l) = pp.lake_input_function(vv.t);
end

%% update volume according to mass conservation
% Vb = zeros(size(vv.Vb));
% forward euler
Q_in = pp.c42*aa.Q_lake;

% Darcy's law
Q_out = pp.c43.*vv.Vb./vv.Rb; % + pp.c44*vv.Rb.*(aa.phi_0-vv.phi);

% Simplified exponential decay (pp.c42*aa.Q_lake <= vv.Vb/10.0).*
% Q_out = (vv.Vb./5.0);
% Q_out = zeros(gg.nIJ,1);
% the outflux should be smaller than the remaining volume /dt
Q_out = min(Q_out, vv.Vb/vv.dt);

% rate of volume variation
V_t = Q_in - Q_out;

% update the volume
vv.Vb = vv.Vb + vv.dt.*V_t;
% the volume should be non-negative
vv.Vb(vv.Vb<0)=0;

% display the volume of the blister
disp(['The blister volume is ' num2str(max(vv.Vb))]);

%% update radius according to toughness/exponential decay
vv.Rb = (vv.Vb>0).*max(vv.Rb,1/pp.c45*vv.Vb.^0.4);
% forward Euler
% vv.Rb = (vv.Vb>0).*1/pp.c45.*(pp.c45*vv.Rb.^(2.5) + (pp.c45.*vv.Rb.^(2.5)<=vv.Vb).*(V_t>0).*vv.dt.*V_t).^(0.4);
aa.Qb = Q_out./gg.Dx(pp.ni_l)./gg.Dy(pp.ni_l);
% dimensionally, aa.Q_b ~ L T^(-1) ~ [E]

end