function Sloping_Boundary_Layer(N)

y = linspace(0,pi/2,N);
% Get A intergal coefficients
A = NaN(N,N);
for i = 2:N-1
    for j = 1:N
        A(i,j) = integral(@(eta) Kern(eta,y(i)) .* Li(j,y,eta) ./ eta.^(2/3), 0 ,pi/2, 'Waypoints',y,'AbsTol',1e-8,'RelTol',1e-8);
    end
end

% Get F integral coefficients
F = NaN(1,N);
for i = 1:N
    F(i) = integral(@(eta) Li(i,y,eta) ./ eta.^(2/3), 0, pi/2, 'Waypoints',y,'AbsTol',1e-8,'RelTol',1e-8);
end

% Load initial guess and solve for solution
IC = load("Initial_Guess_Zero_Toughness.mat",'y','f');
f0 = interp1(IC.y,IC.f,y,'spline','extrap');

options = optimset('MaxFunEvals',1e5,'MaxIter',1e5);
f_new = fsolve(@(f) Constraints(f,y,A,F),f0,options);

% plot(y,f_new)
% xlim([0,pi/2])

% Extract solutions for H and P
H = [0,(A(2:N-1,:)*f_new')',1];
for i = 1:N
    P(i) = -integral(@(eta) 2 * cot(eta) .* interp1(y,f_new,eta) ./ eta.^(2/3),y(i),pi/2,'Waypoints',y(i:N));
end
f = f_new;
X = tan(y).^2;

% Plot results
plot(X,H)
hold on
plot(X,P)
xlim([0,5])

% Save results
save('Sloping_Boundary_Layer','y','f','H','X','P');

end

function Constraints = Constraints(f,y,A,F)

    N = length(y);

    % Elasticity/Mass constraints
    A_constraints = f(2:N-1) ./ (y(2:N-1).^(2/3) .* tan(y(2:N-1)).^2 .* (1+tan(y(2:end-1)).^2)) + 1 - ((A(2:N-1,:) * f')').^(-2);

    % Fracture constraint
    F_constraint = F*f';

    % Origin constraint
    G_constraint = f(1) - 3^(-5/3);

    % Collect constraints together
    Constraints = [G_constraint,A_constraints,F_constraint];

end

