function Sloping_Boundary_Layer(N,K)

y = linspace(0,pi/2,N);
% Get A intergal coefficients
A = NaN(N,N);
for i = 2:N-1
    for j = 1:N
        A(i,j) = integral(@(eta) Kern(eta,y(i)) .* Li(j,y,eta), 0 ,pi/2, 'Waypoints',y,'AbsTol',1e-8,'RelTol',1e-8);
    end
end

% Get F integral coefficients
F = NaN(1,N);
for i = 1:N
    F(i) = integral(@(eta) Li(i,y,eta), 0, pi/2, 'Waypoints',y,'AbsTol',1e-8,'RelTol',1e-8);
end

IC = load("Initial_Guess.mat",'y','f');
f0 = interp1(IC.y,IC.f,y,'spline','extrap');


options = optimset('MaxFunEvals',1e5,'MaxIter',1e5);
f_new = fsolve(@(f) Constraints(f,y,A,F,K),f0,options);

% plot(y,f_new)
% xlim([0,pi/2])

h = [0,(A(2:N-1,:)*f_new')',1];
plot(tan(y).^2,h)
xlim([0,5])

end

function Constraints = Constraints(f,y,A,F,K)

    N = length(y);

    % Elasticity/Mass constraints
    A_constraints = f(2:N-1) ./ (tan(y(2:N-1)).^2 .* (1+tan(y(2:end-1)).^2)) + 1 - ((A(2:N-1,:) * f')').^(-2);

    % Fracture constraint
    F_constraint = sqrt(32)/pi*F*f' + K;

    % Origin constraint
    G_constraint = f(1) - 1/(2*K^2);


    % Collect constraints together
    Constraints = [G_constraint,A_constraints,F_constraint];

end

