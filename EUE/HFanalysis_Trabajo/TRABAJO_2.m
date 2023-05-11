clc
clear 
close all
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.25) % LineWidth
set(groot,'defaultAxesFontSize',25) % Fontsize
%% DATA %%------------------------------------------------------------------------------------------------------------------------------------------------------
p = struct();
a = struct();

p.A = 1.25;           % [m^2]
p.P = 4.52;           % [m]
p.L1 = 1;             % [m]
p.L2 = 1.25;          % [m]
p.t = 0.005;          % [m]
p.rho = 2700;         % [kg/m^3]
p.E = 70E9;           % [Pa]
p.nu = 0.33;          % [-]
a.h = 0.05;           % [m]
a.rho = 1.23;         % [kg/m^3]
a.c0 = 343;           % [m/s]
p.ILF = 0.015;        % [%]
a.CLF = 0.01;           % [%]
p.m = p.rho * p.t * p.A;

% Frequencies and bandwidths
fm = [16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000];
bws = [3.7, 4.6, 5.8, 7.3, 9.2, 11.5, 14.6, 18.3, 22.9, 29, 37, 46, 58, 73, 92, 115, 146, 183, 231, 291, 365, 461, 579, 730, 919, 1156, 1456, 1834, 2307];
w = 2*pi*fm;
bws_w = 2*pi*bws;

% Power associated:
P = zeros(29, 3);
P(:, 1) = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 20 35 50 80 100 150 100 100 100 100]';
P(:, 2) = [4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 4.35 8.7 15.2 21.74 36.96 39.13 45.65 43.47 43.47 43.47 43.47]';
P(:, 3) = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 20 35 50 80 100 150 100 100 100 100]';


%% INMEDIATE CALCULUS %%------------------------------------------------
a.V = p.A*a.h;                                               % [m^3]
a.A = 2 * p.A + 2 * p.L1 * a.h +  2 * p.L2 * a.h;            % [m^2]
a.L = 4 * (p.L1 + p.L2 + a.h) ;                              % [m]
p.D = p.E * p.t ^ 3 / (12 * (1 - p.nu ^ 2));                 % [N·m] (stiffness Al) 
p.n = p.A / (4 * pi) * sqrt(p.rho * p.t / p.D);              % [-] (modal density)
p.fc = (a.c0 ^ 2 / (2 * pi)) * sqrt(p.rho * p.t / p.D);      % [Hz]
p.f11 = a.c0 / (4 * p.fc) * (1 / p.L1 ^ 2 + 1 / p.L2 ^ 2);   % [Hz]

% Variable initialization
a.n = zeros(length(fm), 1); 
a.N = zeros(length(fm), 1);
N_p = zeros(length(fm), 1);
eta_pa = zeros(length(fm), 1);
eta_ap = zeros(length(fm), 1);


E=zeros(16, 5);

%% RESOLUTION %%-----------------------------------------------------------
for i = 1:length(fm)
    a.n(i) = a.V / (pi * a.c0) * (w(i) / a.c0)^2 + a.A / (4 * a.c0) * (w(i) / a.c0) + a.L / (8 * a.c0);
    a.N(i) = a.n(i) * bws_w(i);
    p.N(i) = p.n * bws_w(i);
    eta_pa(i) = p.A * a.rho * a.c0 * sigma(w(i), p, a) / (p.rho * p.A * p.t * w(i));
    eta_ap(i) = p.n / a.n(i) * eta_pa(i);

end
f5 =fm(fm>=315);
% eta_ap = eta_ap(14:end);
% eta_pa = eta_pa(14:end);
% P2 = P(14:end, :);
f_plot = linspace(16, 10000, 9984);
w_plot = 2 * pi * f_plot;
for i = 1:length(fm)
    
    L = [(p.ILF + eta_pa(i)), -eta_ap(i), 0, 0, 0
     -eta_pa(i), (a.CLF + 2*eta_ap(i)), -eta_pa(i), 0, 0
     0, -eta_ap(i), (p.ILF + 2*eta_pa(i)), -eta_ap(i), 0
     0, 0, -eta_pa(i), (a.CLF + 2*eta_ap(i)), -eta_pa(i)
     0, 0, 0, -eta_ap(i), (p.ILF + eta_pa(i))
     ];

    P_syst = [P(i, 1), 0, P(i,2), 0, P(i, 3)]' ./ w(i);
    E(i, :) = L \ P_syst;
end
v_RMS = sqrt(E(:, [1,3,5])./p.m);
P_RMS = sqrt(E(:, [2,4]) * a.rho * a.c0 ^ 2 ./ a.V);
for i = 1:length(f_plot)

    a.n(i) = a.V / (pi * a.c0) * (w_plot(i) / a.c0)^2 + a.A / (4 * a.c0) * (w_plot(i) / a.c0) + a.L / (8 * a.c0);
    eta_pa(i) = p.A * a.rho * a.c0 * sigma(w_plot(i), p, a) / (p.rho * p.A * p.t * w_plot(i));
    eta_ap(i) = p.n / a.n(i) * eta_pa(i);


end








%% FIGURES %%--------------------------------------------------

figure(1)
loglog(fm, (p.N)', Color='k', LineWidth=1.5 );
hold on
loglog( fm, a.N,Color='r', LineWidth=1.5)
yline(5)
xline(315)
grid on
% loglog(fm, );
legend('$$N_p$$' , '$$N_a$$', Location='best')
xlabel('$$f$$ [Hz]')
ylabel('Number of modes')

figure(2)
hold on
grid on
plot(f_plot, log10(eta_ap), LineWidth=1.5, Color='r');
plot(f_plot, log10(eta_pa), LineWidth=1.5, LineStyle="--", Color='k');
legend('$$\eta_{ap}$$' , '$$\eta_{pa}$$', Location='best')
xlabel('$$f$$ [Hz]')
ylabel('CLF')


figure(3)
hold on
grid on
plot(f5, E(14:end, 1), Color='r', LineStyle= "-", LineWidth=1.5);
plot(f5, E(14:end, 2), Color='#7E2F8E', LineStyle= "-", LineWidth=1.5);
plot(f5, E(14:end, 3), Color='b', LineStyle= "-", LineWidth=1.5);
plot(f5, E(14:end, 4), Color='#D95319', LineStyle= "-.", LineWidth=1.5);
plot(f5, E(14:end, 5), Color='k', LineStyle= "--", LineWidth=1.5);

legend('Plate 1' , 'Air layer 1', 'Plate 2' , 'Air layer 2', 'Plate 3', Location='best')
xlabel('$$f$$ [Hz]')
ylabel('$$E$$ [J]')

figure(4)
hold on
grid on
plot(f5, v_RMS(14:end, 1), Color='r', LineStyle= "-", LineWidth=1.5);
plot(f5, v_RMS(14:end, 2), Color='b', LineStyle= "-", LineWidth=1.5);
plot(f5, v_RMS(14:end, 3), Color='k', LineStyle= "--", LineWidth=1.5);

legend('Plate 1' , 'Plate 2' , 'Plate 3', Location='best')
xlabel('$$f$$ [Hz]')
ylabel('$$v_{RMS}$$ [m/s]', Interpreter='latex')

figure(5)
hold on
grid on
plot(f5, P_RMS(14:end, 1), Color='#7E2F8E', LineStyle= "-", LineWidth=1.5);
plot(f5, P_RMS(14:end, 2), Color='#D95319', LineStyle= "--", LineWidth=1.5);

legend( 'Air layer 1', 'Air layer 2', Location='best')
xlabel('$$f$$ [Hz]')
ylabel('$$p_{RMS}$$ [Pa]')






%% FUNCTIONS %%-----------------------------------------------------
function s = sigma(w, p, a)
f = w/(2*pi);
    if p.f11 <= p.fc/2
        Lambda = sqrt(f/p.fc);
        if f<p.f11
            s = 4*p.P^2*(f/a.c0)^2;
        elseif p.f11<f && f <p.fc/2
            s = a.c0*p.P/(p.A*p.fc)*(((1-Lambda^2)*log((1+Lambda)/(1-Lambda)) + 2*Lambda) / (4*pi^2*(1-Lambda^2)^1.5)) + 2* ((2*a.c0)/(p.fc*pi^2))^2*(1-2*Lambda^2)/(p.A*Lambda*sqrt(1-Lambda^2));
        elseif p.fc/2<f && f<p.fc
            s = a.c0*p.P/(p.A*p.fc)*(((1-Lambda^2)*log((1+Lambda)/(1-Lambda)) + 2*Lambda) / (4*pi^2*(1-Lambda^2)^1.5));
        elseif f>p.fc
            s = 1/(sqrt(1-p.fc/f));
        end
    elseif p.f11>p.fc/2
        if f<p.fc
            s = 4*p.A^2*(f/a.c0)^2;
        elseif f>p.fc
            s = 1/(sqrt(1-p.fc/f));
        end
    end
end



