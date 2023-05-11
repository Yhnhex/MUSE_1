% Define the parameters
L = 1.0;  % length of the rod
k = 1.0;  % thermal conductivity
Q = 1.0;  % heat source
T0 = 0.0; % boundary condition at x=0
TL = 1.0; % boundary condition at x=L


% Define the grid
N = 100; % number of grid points
x = linspace(0, L, N);

% Define the initial condition
T_init = [T0; 0];

% Solve the ODE system
[T,~] = ode45(@f, x, T_init);

% Plot the results
plot(x, T(:,1), 'LineWidth', 2);
xlabel('x');
ylabel('T');
title('Temperature Distribution');
grid on;


% Define the function to solve the ODE system
function dydx = f(x, y)
    k = 1.0;  % thermal conductivity
    Q = 1.0;
    % y(1) = T(x), y(2) = T'(x)
    dydx = zeros(2,1);
    dydx(1) = y(2);
    dydx(2) = -Q/k;
end
