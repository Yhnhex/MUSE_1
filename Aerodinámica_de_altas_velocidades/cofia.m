clc
clear
close all

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.25) % LineWidth
set(groot,'defaultAxesFontSize',15) % Fontsize

syms alpha phi x
RL = 2.7; L = 10;
g = 1.4; M = 10;
Cpmax = 2* (((((g+1)*M)^2)/(4*g*M^2 - 2*(g-1)))^(g/(g-1))*((1-g+2*g*M^2)/(g+1))-1)/(g*M^2);

r = RL.*(x./L).^(1/3);
dr = (RL*x^(-2/3))/(3*L^(1/3));
num = sind(alpha)*cosd(phi)-dr*cosd(alpha);
denom = 1 + dr^2;

Cp = 2*(num^2/denom);

alphas = [5 15 20 25 30];

figure(1)
hold on
grid on
for i = 1:length(alphas)
    Lx = double(solve(subs(subs(Cp,phi,0),alpha,alphas(i)) == 0,x));
    fplot(subs(subs(Cp,phi,0),alpha,alphas(i)),[0 Lx], LineWidth = 1.5)
end
xlim([0 6])
legend('$\alpha = 5^{o}$','$\alpha = 15^{o}$','$\alpha = 20^{o}$','$\alpha = 25^{o}$','$\alpha = 30^{o}$')
hold off

alphas = linspace(0,45,45);
CD = zeros(1,length(alphas));
CD_m = zeros(1,length(alphas));
CL = zeros(1,length(alphas));
CL_m = zeros(1,length(alphas));
CMy = zeros(1,length(alphas));
CMy_m = zeros(1,length(alphas));

for i = 1:length(alphas)
    alpha = deg2rad(alphas(i));

    % Cds
    Cd = @(x,y) -(2/(pi*RL^2)).*((sin(alpha).*cos(y) - (RL.*x.^(-2./3))./(3.*L.^(1./3)).*cos(alpha)).^3).*(RL.*(x./L).^(1./3))./(1+ ((RL.*x.^(-2./3))./(3.*L.^(1./3))).^(2)).^(3./2);
    Cd_m = @(x,y) -(Cpmax/(pi*RL^2)).*((sin(alpha).*cos(y) - (RL.*x.^(-2./3))./(3.*L.^(1./3)).*cos(alpha)).^3).*(RL.*(x./L).^(1./3))./(1+ ((RL.*x.^(-2./3))./(3.*L.^(1./3))).^(2)).^(3./2);
    CD(i) = integral2(Cd,0,L,0,2*pi);
    CD_m(i) = integral2(Cd_m,0,L,0,2*pi);

    % CLs
    Cl = @(x,y) -(2/(pi*RL^2)).*((sin(alpha).*cos(y) - (RL.*x.^(-2./3))./(3.*L.^(1./3)).*cos(alpha)).^2).*(RL.*(x./L).^(1./3)).*((RL.*x.^(-2./3))./(3.*L.^(1./3)).*sin(alpha) + cos(alpha).*cos(y))./(1+ ((RL.*x.^(-2./3))./(3.*L.^(1./3))).^(2)).^(3./2);
    CL(i) = integral2(Cl,0,L,0,2*pi);
    Cl_m = @(x,y) -(Cpmax/(pi*RL^2)).*((sin(alpha).*cos(y) - (RL.*x.^(-2./3))./(3.*L.^(1./3)).*cos(alpha)).^2).*(RL.*(x./L).^(1./3)).*((RL.*x.^(-2./3))./(3.*L.^(1./3)).*sin(alpha) + cos(alpha).*cos(y))./(1+ ((RL.*x.^(-2./3))./(3.*L.^(1./3))).^(2)).^(3./2);
    CL_m(i) = integral2(Cl_m,0,L,0,2*pi);
    
    % Cmy
    Cm = @(x,y) (2/(pi*L*RL^2)).*((sin(alpha).*cos(y) - (RL.*x.^(-2./3))./(3.*L.^(1./3)).*cos(alpha)).^(2)).*(x.*cos(y)+ (RL.*(x./L).^(1./3)).*((RL.*x.^(-2./3))./(3.*L.^(1./3)))*cos(alpha)).*(RL.*(x./L).^(1./3))./(1+ ((RL.*x.^(-2./3))./(3.*L.^(1./3))).^(2)).^(3./2);
    CMy(i) = integral2(Cm,0,L,0,2*pi);
    Cm_m = @(x,y) (Cpmax/(pi*L*RL^2)).*((sin(alpha).*cos(y) - (RL.*x.^(-2./3))./(3.*L.^(1./3)).*cos(alpha)).^(2)).*(x.*cos(y)+ (RL.*(x./L).^(1./3)).*((RL.*x.^(-2./3))./(3.*L.^(1./3)))*cos(alpha)).*(RL.*(x./L).^(1./3))./(1+ ((RL.*x.^(-2./3))./(3.*L.^(1./3))).^(2)).^(3./2);
    CMy_m(i) = integral2(Cm_m,0,L,0,2*pi);
    i = i
end

figure(2)
hold on
grid on
plot(alphas,CD,'k')
plot(alphas,CD_m,'r')

legend('Newton','Newton Modified',Location='southeast')
xlabel('\''Angulo de ataque $\alpha$ [deg]')
ylabel('$C_{D}$')
title('$C_{D}$ para $M=10$')
hold off

figure(3)
hold on
grid on
plot(alphas,CL,'k')
plot(alphas,CL_m,'r')

legend('Newton','Newton Modified',Location='southeast')
xlabel('\''Angulo de ataque $\alpha$ [deg]')
ylabel('$C_{L}$')
title('$C_{L}$ para $M=10$')

hold off

figure(4)
hold on
grid on
plot(alphas,CMy,'k')
plot(alphas,CMy_m,'r',LineStyle='--')

legend('Newton','Newton Modified',Location='southeast')
xlabel('\''Angulo de ataque $\alpha$ [deg]')
ylabel('$C_{My}$')
title('$C_{My}$ para $M=10$')

hold off