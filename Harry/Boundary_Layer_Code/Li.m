function Li = Li(i,eta,etaq)

N = length(eta);
Li = zeros(1,length(etaq));
if i == 1
    R_ind = logical(etaq<eta(2));
    Li(R_ind) = 1 - etaq(R_ind)/eta(2);
elseif i == N
    L_ind = logical(etaq>eta(N-1));
    Li(L_ind) = 1 - (eta(N)-etaq(L_ind))/(eta(N)-eta(N-1));
else
    L_ind = logical((eta(i-1)<etaq) .* (etaq<eta(i)));
    R_ind = logical((eta(i)<etaq) .* (etaq<eta(i+1)));
    Li(L_ind) = 1 - (eta(i)-etaq(L_ind))/(eta(i)-eta(i-1));
    Li(R_ind) = 1 - (etaq(R_ind)-eta(i))/(eta(i+1)-eta(i));
end

end