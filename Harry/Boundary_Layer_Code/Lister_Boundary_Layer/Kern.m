function K = Kern(eta,y)

K = 2/pi * cot(eta) .* ( (tan(y)^2 - tan(eta).^2) .* log( abs( (tan(y)+tan(eta))./(tan(y)-tan(eta)) ) ) - 2 * tan(y) .* tan(eta) );

K((eta == y)) = -4/pi * tan(y);
K((eta == pi/2)) = -8/pi * tan(y);
K((eta == 0)) = 0;

end