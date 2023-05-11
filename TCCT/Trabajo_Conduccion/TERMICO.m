
rho_CFRP = 1650;
rho_HC = 40 ;
CP_CFRP = 930;
CP_HC = 890;
kx_CFRP = 130;
ky_CFRP = 1;
kz_CFRP = 130;
kx_HC = 0.7;
ky_HC = 1.5;
kz_HC = 0.7;
X = 500e-3;
Y_CFRP = 2.5e-3;
Y_HC = 15e-3;
Z = 300e-3;
T0 = 273 ;
q = 50;

R_CFRP = X/(ky_CFRP*Z*Y_CFRP);
R_HC = X/(ky_HC*Z*Y_HC);
R_eq1 = 1/(2/R_CFRP +1/R_HC);
k_eq1 =  X/(R_eq1*(Z*Y_HC + 2*Z*Y_CFRP));

x = linspace(0, X, 100);
[Temperature, dT] = ode45(@f, x, [T0 T0] );

plot(x, T(:,1), 'LineWidth', 2);
xlabel('x');
ylabel('T');
title('Temperature Distribution');
grid on;

function dudx = f(x, u)
    % y(1) = T(x), u(2) = T'(x)
    dudx = zeros(2,1);
    dudx(1) = u(2);
    dudx(2) = -q/k_eq1;
end
