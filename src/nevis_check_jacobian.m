function [J,J2] = nevis_check_jacobian(dt,vv,vv0,aa,pp,gg,oo)
% 18 Sept 2016
% check jacobian with numerically calculated entries


% calculate Xi
oo.evaluate_variables = 1;
[vv2,~,~,~,~,~,~,~,~] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo); 
aa.Xi = vv2.Xi; oo.includeXi = 0; 

% jacobian
oo.evaluate_variables = 0; oo.evaluate_jacobian = 1;
[~,F,~,~,~,~,~,~,J] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo); 

% construct vector of variables
if oo.no_channels && oo.no_sheet
    X = vv.phi(gg.nin);
elseif oo.no_channels
    X = [vv.hs(gg.ns); vv.phi(gg.nin)];
elseif oo.no_sheet
    X = [vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin)];
else
    X = [vv.hs(gg.ns); vv.phi(gg.nin); vv.Sx(gg.ein); vv.Sy(gg.fin); vv.Ss(gg.cin); vv.Sr(gg.cin)];
end
% numerically alculated Jacobian
J2 = sparse(size(J,1),size(J,2));
v0 = vv;
X0 = X;
step = 1e-8;
% run through each variable perturbing
for i = 1:length(X),
    % reset vv and X
    X = X0;
    vv = v0;
    % perturb
    X(i) = X(i)+step;
    
    % construct perturbed vv
    temp1 = 0;
    if ~oo.no_sheet
        temp2 = length(gg.ns); vv.hs(gg.ns) = X(temp1+(1:temp2)); temp1=temp1+temp2;
    end
    temp2 = length(gg.nin); vv.phi(gg.nin) = X(temp1+(1:temp2)); temp1=temp1+temp2;
    if ~oo.no_channels
        temp2 = length(gg.ein); vv.Sx(gg.ein) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        temp2 = length(gg.fin); vv.Sy(gg.fin) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        temp2 = length(gg.cin); vv.Ss(gg.cin) = X(temp1+(1:temp2)); temp1=temp1+temp2;
        temp2 = length(gg.cin); vv.Sr(gg.cin) = X(temp1+(1:temp2));
    end
    
    % evaluate pertubed F
    oo.evaluate_jacobian = 0;
    [~,Fstep,~,~,~,~,~,~,~] = nevis_backbone(dt,vv,vv0,aa,pp,gg,oo); 
    J2(:,i) = (Fstep-F)/step;
end

figure(1); clf;
    spy(J);
figure(2); clf;
    spy(J2);
figure(3); clf;
    spy(abs(J-J2)>1e-4);

end
