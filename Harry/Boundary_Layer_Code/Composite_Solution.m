function [x,h,xi] = Composite_Solution(t,Ll,G)

Lr = Ll + 3/2^(2/3) * G^(1/3) * t^(1/3);
x_mid = Lr - 10 * (3 / 2 / G / (Lr-Ll))^(1/2);

if x_mid <= Ll
    x = Ll + (Lr-Ll) * linspace(0,1,200);
else
    x = [Ll + (x_mid-Ll)*linspace(0,1,100).^3, x_mid + (Lr-x_mid) * flip(1-linspace(0,1,100).^3)];
end

load('Sloping_Boundary_Layer.mat','X','H')

Xq = (Lr-x) * (2*G*(Lr-Ll)/3)^(1/2);
h = 3 / (2*(Lr-Ll)) * ( ((x-Ll)/(Lr-Ll)).^(1/2) + interp1(X,H,Xq) - 1);

xi = 2*x/Lr-1;

% plot(x,h)
% hold on
% % scatter(x,h)
% title(num2str(trapz(x,h)))

end