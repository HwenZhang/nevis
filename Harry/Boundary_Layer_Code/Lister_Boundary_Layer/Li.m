function Li = Li(i,y,eta)

N = length(y);
Li = zeros(1,length(eta));
if i == 1
    R_ind = logical(eta<=y(2));
    Li(R_ind) = 1 - eta(R_ind)/y(2);
elseif i == N
    L_ind = logical(eta>y(N-1));
    Li(L_ind) = 1 - (y(N)-eta(L_ind))/(y(N)-y(N-1));
else
    L_ind = logical((y(i-1)<eta) .* (eta<=y(i)));
    R_ind = logical((y(i)<eta) .* (eta<y(i+1)));
    Li(L_ind) = 1 - (y(i)-eta(L_ind))/(y(i)-y(i-1));
    Li(R_ind) = 1 - (eta(R_ind)-y(i))/(y(i+1)-y(i));
end

end