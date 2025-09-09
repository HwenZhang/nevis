function K = Kern(eta,xi)

K = 4/pi * cot(eta) .* ( (tan(xi)^2 - tan(eta).^2) .* log(abs((tan(xi)+tan(eta))./(tan(xi)-tan(eta)))) - 2 * tan(xi) .* tan(eta) );

K((eta == xi)) = -8/pi * tan(xi);
K((eta == pi/2)) = -16/pi * tan(xi);
K((eta == 0)) = 0;

% K = K/2; % For comparison with Lister

end