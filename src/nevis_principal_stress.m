function [sigma1,sigma2,theta] = nevis_principal_stress(sigmaxx,sigmayy,sigmaxy)
% calculate principal stress components and angle for two-dimensional
% tensor with diagonal components sigmaxx,sigmayy and off-diagonal sigmaxy
%
% IJH 22 October 2025

    sigma1 = (sigmaxx+sigmayy)/2+(sigmaxy.^2+((sigmaxx-sigmayy)/2).^2).^(1/2);
    sigma2 = (sigmaxx+sigmayy)/2-(sigmaxy.^2+((sigmaxx-sigmayy)/2).^2).^(1/2);
    theta = 1/2*atan2(sigmaxy,(sigmaxx-sigmayy)/2);

end