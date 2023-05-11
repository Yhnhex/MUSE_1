tic
clc
clear
close all

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.25) % LineWidth
set(groot,'defaultAxesFontSize',10) % Fontsize


% Constantes
G = 1367;
eta = 29.5;
kbol = 1.380649e-23;
q = 1.6e-19;

Resis = [15.8 16.7 17.2];  % Resistencias


% Definición de los parámetros del sistema
omega = 0.02;               % [rad/s]
Periodo = 2*pi/omega;       % [s]

T0_illuminated = 80;        % [°C]
T0_not_illuminated = -20;   % [°C]

V_oc_0 = 2461.36e-3;        % [V]
I_sc_0 = 506.34e-3;         % [A]
V_mp_0 = 2153.117e-3;       % [V]
I_mp_0 = 495.16e-3;         % [A]

% Definir el ángulo de inclinación de la radiación solar en función del tiempo

t = 0:0.1:Periodo;          % Tiempo de simulación [s]
theta = 111.5 + rad2deg(omega*t);     % [°]

% Temperaturas en función de Theta
Temperatura = zeros(length(theta),1);
for i = 1:length(theta)
    Temperatura(i) = LeyTemperatura(theta(i));
end 

% Gradientes de temperatura
V_oc = (V_oc_0 - (6.2e-3)*Temperatura)*7;       % [V]
I_sc = (I_sc_0 + (0.36e-3)*Temperatura)*2;      % [A]
V_mp = (V_mp_0 - (6.7e-3)*Temperatura)*7;       % [V]
I_mp = (I_mp_0 + (0.24e-3)*Temperatura)*2;      % [A]

betas = I_mp./I_sc;
alphas = V_mp./V_oc;

for j = 1:length(Resis)
    R = Resis(j);
    % Karalkar y Haeefa
    
    C = (1 - betas - alphas)./(2*betas - 1);
    
    IV_KF = zeros(length(alphas),2);
    
    for i = 1:length(alphas)
        alpha = alphas(i);
        beta = betas(i);
        eme = 1 + 1/C(i) + Wfunction(-(log(alpha)*alpha^(-1/C(i))/C(i)))/log(alpha);
        gamma = (2*beta - 1)/((eme-1)*alpha^(eme));
    
        f_KF = @(I) 1 - (1 - gamma)*(I*R)/V_oc(i) - gamma*((I*R)/V_oc(i))^eme - I/I_sc(i);
        I = fzero(f_KF,1);
        IV_KF(i,1) = I;
    end
    IV_KF_T(:,1) = IV_KF(:,1);
    IV_KF_T(:,2) = IV_KF(:,1)*R;
    
    IV_KF(cosd(theta) < cosd(80),1) = 0;
    IV_KF(cosd(theta) < cosd(80),2) = 0;
    IV_KF(:,1) = IV_KF(:,1).*sqrt(cosd(theta'));
    IV_KF(:,2) = IV_KF(:,1)*R;
    
    P_KF = IV_KF(:,1).*IV_KF(:,2);
    P_KF(cosd(theta) < 0) = 0;
    
    
    % Modelo de Das
    
    IV_DAS = zeros(length(alphas),2);
    
    for i = 1:length(alphas)
        alpha = alphas(i);
        beta = betas(i);
        k = Wfunction(beta*log(alpha))/log(alpha);
        h =  (1/beta - 1/k - 1)/alpha;
        f_DAS = @(I) (1 - (I*R/V_oc(i))^k)/(1 + h*(I*R/V_oc(i))) - I/I_sc(i);
        I = fzero(f_DAS,1);
        IV_DAS(i,1) = I;
    end
    
    IV_DAS_T(:,1) = IV_DAS(:,1);
    IV_DAS_T(:,2) = IV_DAS_T(:,1)*R;
    
    IV_DAS(cosd(theta) < cosd(80),1) = 0;
    IV_DAS(cosd(theta) < cosd(80),2) = 0;
    IV_DAS(:,1) = IV_DAS(:,1).*sqrt(cosd(theta'));
    IV_DAS(:,2) = IV_DAS(:,1)*R;
    
    P_DAS = IV_DAS(:,1).*IV_DAS(:,2);
    P_DAS(cosd(theta) < 0) = 0;
    
    % Modelo de Pindado-Cubas
    
    IV_PC_1 = zeros(length(alphas),2);
    
    for i = 1:length(alphas)
        alpha = alphas(i);
        beta = betas(i);
        eta = (1-alpha)/(1-beta)/beta;
    
        f_PC_1 = @(I) 1 - (1.-beta)*((I*R/V_oc(i))/alpha)^(beta/(1-beta)) - I/I_sc(i);
        I1 = real(fzero(f_PC_1,1));
        IV_PC_1(i,1) = I1;
        IV_PC_1(i,2) = I1*R;
    end
    
    % La segunda función se debe aplicar cuando v>alpha, esto es cuando la
    % siguiente operación sea positiva.
    
    signo = IV_PC_1(i,2)./V_oc - alphas;
    
    alpha2 = alphas(signo>0);
    beta2 = betas(signo>0);
    
    IV_PC_2 = zeros(length(alpha2),2);
    
    V_oc_2 = V_oc(signo>0);
    I_sc_2 = I_sc(signo>0);
    
    for i = 1:length(alpha2)
        alpha = alpha2(i);
        beta = beta2(i);
        eta = (1-alpha)/(1-beta)/beta;
    
        f_PC_2 = @(I) beta*(alpha/(I*R/V_oc_2(i)))*(1 - (((I*R/V_oc_2(i)) - alpha)/(1 - alpha))^eta) - I/I_sc_2(i);
        I2 = real(fzero(f_PC_2,1));
        IV_PC_2(i,1) = I2;
    end 
    
    IV_PC = IV_PC_1; 
    
    % Sustituyo en los valores que toca con la segunda función
    IV_PC((signo>0),:) = IV_PC_2(:,:);
    
    IV_PC_T(:,1) = IV_PC(:,1);
    IV_PC_T(:,2) = IV_PC_T(:,1)*R;
    
    IV_PC(cosd(theta) < cosd(80),1) = 0;
    IV_PC(cosd(theta) < cosd(80),2) = 0;
    IV_PC(:,1) = IV_PC(:,1).*sqrt(cosd(theta'));
    IV_PC(:,2) = IV_PC(:,1)*R;
    
    P_PC = IV_PC(:,1).*IV_PC(:,2);
    P_PC(cosd(theta) < 0) = 0;
    
    % 1D-2R
    
    IV_1D2R = zeros(length(theta),2);
    
    for i = 1:length(theta)
        
        V_T = kbol*(Temperatura(i) + 273 )/q;
    
        % CALCULATE APPROX R_S
        a_prox=1;
        fun_R_s = @(R_ss) (7*a_prox*V_T*V_mp(i)*(2*I_mp(i) - I_sc(i)))/((V_mp(i)*I_sc(i)+V_oc(i)*(I_mp(i)-I_sc(i)))*(V_mp(i)-I_mp(i)*R_ss)-7*a_prox*V_T*(V_mp(i)*I_sc(i)-V_oc(i)*I_mp(i)))-exp((V_mp(i)+I_mp(i)*R_ss-V_oc(i))/(7*a_prox*V_T));
        R_s = fzero(fun_R_s, 1.3535);
        
        % CALCULATE R_SHUNT
        
        R_sh= ((V_mp(i)-I_mp(i)*R_s)*(V_mp(i)-R_s*(I_sc(i)-I_mp(i))-7*a_prox*V_T))/((V_mp(i)-I_mp(i)*R_s)*(I_sc(i)-I_mp(i))-7*a_prox*V_T*I_mp(i));
        % CALCULATE REAL a
        R_sh0= R_s+R_sh;
        a =((V_mp(i)-I_mp(i)*R_s)*(V_mp(i)+(I_mp(i)-I_sc(i))*R_sh0))/((V_mp(i)-I_mp(i)*R_sh0)*V_T*7);
        
        
        % CALCULATE REAL R_S
        A = (V_mp(i)+(I_mp(i)-I_sc(i))*R_sh0)*log((V_mp(i)+(I_mp(i)-I_sc(i))*R_sh0)/(V_oc(i)-I_sc(i)*R_sh0));
        B = V_mp(i) - R_sh0*I_mp(i);
        R_s = (A-B)/(A+B)*V_mp(i)/I_mp(i) + B/(A+B)*V_oc(i)/I_mp(i);
        
        % CALCULATE I_pv
        I_pv = (R_sh+R_s)/R_sh*I_sc(i);
        I_o = ((R_sh+R_s)*I_sc(i)-V_oc(i))/(R_sh*exp(V_oc(i)/(7*a*V_T)));
        
        % define function you want to solve for: fun = 0
        fun   = @(I) I_pv - I_o * (exp((I*R + I*R_s)/(7*a*V_T)) - 1) - (I*R + I*R_s)/R_sh - I;
        % get solution for all V
        I_sol = fzero(fun, I_pv);
        IV_1D2R(i,1) = I_sol;
    end
    
    IV_1D2R_T(:,1) = IV_1D2R(:,1);
    IV_1D2R_T(:,2) = IV_1D2R_T(:,1)*R;
    
    IV_1D2R(cosd(theta) < cosd(80),1) = 0;
    IV_1D2R(cosd(theta) < cosd(80),2) = 0;
    IV_1D2R(:,1) = IV_1D2R(:,1).*sqrt(cosd(theta'));
    IV_1D2R(:,2) = IV_1D2R(:,1)*R;
    
    P_1D2R = IV_1D2R(:,1).*IV_1D2R(:,2);
    P_1D2R(cosd(theta) < 0) = 0;

    figure(1)
    subplot(1,2,1)
    hold on
    grid on
    plot(theta,IV_KF(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off
    subplot(1,2,2)
    hold on
    grid on
    plot(theta,P_KF(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

    figure(2)
    subplot(1,2,1)
    hold on
    grid on
    plot(theta,IV_DAS(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off
    subplot(1,2,2)
    hold on
    grid on
    plot(theta,P_DAS(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

    figure(3)
    subplot(1,2,1)
    hold on
    grid on
    plot(theta,IV_PC(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off
    subplot(1,2,2)
    hold on
    grid on
    plot(theta,P_PC(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

    figure(4)
    subplot(1,2,1)
    hold on
    grid on
    plot(theta,IV_1D2R(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off
    subplot(1,2,2)
    hold on
    grid on
    plot(theta,P_1D2R(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

    figure(5)
    hold on
    grid on
    plot(theta,IV_KF_T(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

    figure(6)
    hold on
    grid on
    plot(theta,IV_DAS_T(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

    figure(7)
    hold on
    grid on
    plot(theta,IV_PC_T(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

    figure(8)
    hold on
    grid on
    plot(theta,IV_1D2R_T(:,1),LineWidth=2)
    if j == 3
        legend('$R = 15.8 \Omega$','$R = 16.7 \Omega$', '$R = 17.2 \Omega$')
    end
    hold off

end 
toc







function Temperatura = LeyTemperatura(theta)
T0_illuminated = 80;                % [°C]
T0_not_illuminated = -20;           % [°C]
delta_T = T0_illuminated - T0_not_illuminated;
phi = theta + 90;

Tmax = T0_not_illuminated + delta_T*(1-exp(-180/45));
Tmin = Tmax - (Tmax-T0_not_illuminated)*(1-exp(-(180)/30));

while phi > 360
    phi = phi - 360;
end 
if (phi>= 0) && (phi <= 180)
    Temperatura = Tmin + (T0_illuminated - Tmin)*(1-exp(-phi/45));
elseif (phi > 180) && (phi < 360)
    Temperatura = Tmax - (Tmax-T0_not_illuminated)*(1-exp(-(phi-180)/30));
end 
end

function W = Wfunction(x)
     if (x >= 2e-16) && (x <= 2e-1)

         W = x*exp(0.71116*x^2 - 0.98639*x);

     elseif (x > 0.2) && (x <= 1.2)
         
         W = -1.6579 + 0.1396*(2.9179e5 - (x - 22.8345)^4)^0.25;

     elseif (x > 1.2) && (x <= 10)

         W = -1.2216 + (3.4724e-2)*(1.7091e8 - (x - 114.146)^4)^0.25 ;
    
     elseif (x >= -1e-2) && (x <= -5e-13)

         W = 9.7117e-5*log(-x)^3 + (6.8669e-3)*log(-x)^2 + 1.2*log(-x) -1.1102;

     elseif (x > -5e-13) && (x <= -1e-40)

        W = (1.6705e-6)*log(-x)^3 + (4.4514e-4)*log(-x)^2 + 1.0511*log(-x) - 2.3364;

     elseif (x >= -0.36785) && (x <= -0.27)

         W = -1 -sqrt(42.949*x^2 + 37.694*x + 8.0542);

     elseif (x > -0.27) && (x <= -0.0732) 

        W = 0.14279*log(-x)^3 + 1.04416*log(-x)^2 + 3.92*log(-x) + 1.65795;

     end
end
