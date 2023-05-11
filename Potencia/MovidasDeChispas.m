clc
clear
close all

colors = ["#000000","#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];

% Parámetros de cada uno

RTC_France= readtable('IV_curves.xlsx','Sheet','RTC France','Range', 'B1:B6');
TNJ = readtable('IV_curves.xlsx','Sheet','TNJ','Range', 'B1:B6');
ZTJ= readtable('IV_curves.xlsx','Sheet','ZTJ','Range', 'B1:B6');
P3G30C= readtable('IV_curves.xlsx','Sheet','3G30C','Range', 'B1:B6');
PWP201= readtable('IV_curves.xlsx','Sheet','PWP201','Range', 'B1:B6');
KC200GT2 = readtable('IV_curves.xlsx','Sheet','KC200GT2','Range', 'B1:B6');
SPVSX5 = readtable('IV_curves.xlsx','Sheet','SPVSX5','Range', 'B1:B6');
PSC = readtable('IV_curves.xlsx','Sheet','PSC','Range', 'B1:B6');


% Datos experimentales de cada uno

C_RTC_France = rmmissing(readtable('IV_curves.xlsx','Sheet','RTC France'),1);
C_TNJ = rmmissing(readtable('IV_curves.xlsx','Sheet','TNJ'),1);
C_ZTJ = rmmissing(readtable('IV_curves.xlsx','Sheet','ZTJ'),1);
C_3G30C = rmmissing(readtable('IV_curves.xlsx','Sheet','3G30C'),1);
C_PWP201 = rmmissing(readtable('IV_curves.xlsx','Sheet','PWP201'),1);
C_KC200GT2 = rmmissing(readtable('IV_curves.xlsx','Sheet','KC200GT2'),1);
C_SPVSX5 = rmmissing(readtable('IV_curves.xlsx','Sheet','SPVSX5'),1);
C_PSC = rmmissing(readtable('IV_curves.xlsx','Sheet','PSC'),1);

figure('Name','Datos experimentales')

hold on
grid on

plot(C_RTC_France{:,1},C_RTC_France{:,2},LineWidth=2)
plot(C_TNJ{:,1},C_TNJ{:,2},LineWidth=2)
plot(C_ZTJ{:,1},C_ZTJ{:,2},LineWidth=2)
plot(C_3G30C{:,1},C_3G30C{:,2},LineWidth=2)
plot(C_PWP201{:,1},C_PWP201{:,2},LineWidth=2)
plot(C_KC200GT2{:,1},C_KC200GT2{:,2},LineWidth=2)
plot(C_SPVSX5{:,1},C_SPVSX5{:,2},LineWidth=2)
plot(C_PSC{:,1},C_PSC{:,2},LineWidth=2)


xlabel('Tension [V]','Interpreter','latex')
ylabel('Intensidad experimental [A]','Interpreter','latex','Rotation',90)

legend('RTC France','TNJ','ZTJ','3G30C','PWP201','KC200GT2','SPVSX5','PSC')

hold off

% i = I/Isc , v = V/Voc

alphas = [RTC_France{6,1} TNJ{6,1} ZTJ{6,1} P3G30C{6,1} PWP201{6,1} KC200GT2{6,1} SPVSX5{6,1} PSC{6,1}];
betas = [RTC_France{5,1} TNJ{5,1} ZTJ{5,1} P3G30C{5,1} PWP201{5,1} KC200GT2{5,1} SPVSX5{5,1} PSC{5,1}];
Iscs = [RTC_France{1,1} TNJ{1,1} ZTJ{1,1} P3G30C{1,1} PWP201{1,1} KC200GT2{1,1} SPVSX5{1,1} PSC{1,1}];
Vocs = [RTC_France{4,1} TNJ{4,1} ZTJ{4,1} P3G30C{4,1} PWP201{4,1} KC200GT2{4,1} SPVSX5{4,1} PSC{4,1}];

% Karalkar y Haeefa

gammas = zeros(1,8); emes = zeros(1,8);

syms  v

figure('Name','Karalkar y Haeefa')
hold on
grid on


for i = 1:size(alphas,2)
    
    alpha = alphas(i);
    beta = betas(i);
    
    C = (1 - beta - alpha)/(2*beta - 1);
    
    emes(i) = 1 + 1/C + Wfunction(-(log(alpha)*alpha^(-1/C))/C)/log(alpha);
    gammas(i) = (2*beta - 1)/((emes(i)-1)*alpha^(emes(i)));
    
    i_ad_KF = 1 - (1 - gammas(i))*v - gammas(i)*v^emes(i);

    fplot(i_ad_KF,[0 1], LineWidth=1.2, Color = colors(i))

end
xlabel('$$ v = V/V_{oc} $$ ','Interpreter','latex', FontSize=14)
ylabel('$$ i = I/I_{sc} $$ ','Interpreter','latex', FontSize=14)
title('Karalkar y Haeefa')
xlim([0 1.1])
ylim([0 1.1])
legend('RTC France','TNJ','ZTJ','3G30C','PWP201','KC200GT2','SPVSX5','PSC',Location='southwest')
hold off

% Modelo de Das

hs = zeros(1,8); ks = zeros(1,8);

figure('Name','DAS')
hold on
grid on

for i = 1:size(alphas,2)
    
    alpha = alphas(i);
    beta = betas(i);
    
    k = Wfunction(beta*log(alpha))/log(alpha);
    h =  (1/beta - 1/k - 1)/alpha;

    hs(i) = h; ks(i) = k;

    i_ad_DAS = (1 - v^k)/(1 + h*v);
    
    fplot(i_ad_DAS,[0 1], LineWidth=1.2, Color = colors(i))    

end

xlabel('$$ v = V/V_{oc} $$ ','Interpreter','latex', FontSize=14)
ylabel('$$ i = I/I_{sc} $$ ','Interpreter','latex', FontSize=14)
title('DAS')
xlim([0 1.1])
ylim([0 1.1])
legend('RTC France','TNJ','ZTJ','3G30C','PWP201','KC200GT2','SPVSX5','PSC',Location='southwest')
hold off



% Modelo de Pindado-Cubas

etas = zeros(1,8);

figure('Name','P-C')
hold on
grid on

etas = zeros(1,8);
for i = 1:size(alphas,2)
    
    alpha = alphas(i);
    beta = betas(i);
    
    eta = (1-alpha)/(1-beta)/beta;
    etas(i) = eta;

   
    i_ad_PC_1 = 1 - (1-beta)*(v/alpha)^(beta/(1-beta));

    i_ad_PC_2 = beta*(alpha/v)*(1 - ((v - alpha)/(1 - alpha))^eta);
    
    
    fig(i) = fplot(i_ad_PC_1,[0 alpha], LineWidth=1.2, Color = colors(i));
    fplot(i_ad_PC_2,[alpha 1], LineWidth=1.2, Color = colors(i))

end

xlabel('$$ v = V/V_{oc} $$ ','Interpreter','latex', FontSize=14)
ylabel('$$ i = I/I_{sc} $$ ','Interpreter','latex', FontSize=14)
title('Pindado - Cubas')
xlim([0 1.1])
ylim([0 1.1])
legend(fig([1 2 3 4 5 6 7 8]),'RTC France','TNJ','ZTJ','3G30C','PWP201','KC200GT2','SPVSX5','PSC',Location='southwest')


% Comparativa entre métodos

names = {'RTC France' 'TNJ' 'ZTJ' '3G30C' 'PWP201' 'KC200GT2' 'SPVSX5' 'PSC'};

for i = 1:size(alphas,2)   
    figure('Name',names{i})
    hold on
    grid on
    alpha = alphas(i);
    beta = betas(i);
    eta = etas(i);
    k = ks(i);
    h = hs(i);
    
    i_ad_KF = 1 - (1 - gammas(i))*v - gammas(i)*v^emes(i);
    i_ad_DAS = (1 - v^k)/(1 + h*v);
    i_ad_PC_1 = 1 - (1-beta)*(v/alpha)^(beta/(1-beta));
    i_ad_PC_2 = beta*(alpha/v)*(1 - ((v - alpha)/(1 - alpha))^eta);

    fig(1) = fplot(i_ad_KF,[0 1], LineWidth=1.2);
    fig(2) = fplot(i_ad_DAS,[0 1], LineWidth=1.2);   
    fig(4) = fplot(i_ad_PC_1,[0 alpha], LineWidth=1.2, Color = colors(5));
  
    fig(5) = fplot(i_ad_PC_2,[alpha 1], LineWidth=1.2, Color = colors(5));
    
    if i == 1
       fig(3) = plot(C_RTC_France{:,1}/Vocs(i),C_RTC_France{:,2}/Iscs(i),'--',LineWidth=2);
    elseif i == 2
       fig(3) = plot(C_TNJ{:,1}/Vocs(i),C_TNJ{:,2}/Iscs(i),'--',LineWidth=2);
    elseif i == 3    
       fig(3) = plot(C_ZTJ{:,1}/Vocs(i),C_ZTJ{:,2}/Iscs(i),'--',LineWidth=2);
    elseif i == 4
       fig(3) = plot(C_3G30C{:,1}/Vocs(i),C_3G30C{:,2}/Iscs(i),'--',LineWidth=2);
    elseif i == 5
       fig(3) = plot(C_PWP201{:,1}/Vocs(i),C_PWP201{:,2}/Iscs(i),'--',LineWidth=2);
    elseif i == 6    
       fig(3) = plot(C_KC200GT2{:,1}/Vocs(i),C_KC200GT2{:,2}/Iscs(i),'--',LineWidth=2);
    elseif i == 7    
       fig(3) = plot(C_SPVSX5{:,1}/Vocs(i),C_SPVSX5{:,2}/Iscs(i),'--',LineWidth=2);
    else
       fig(3) = plot(C_PSC{:,1}/Vocs(i),C_PSC{:,2}/Iscs(i),'--',LineWidth=2);
    end 


    legend(fig([1 2 4 3]),'K-F','DAS','P-C','Experimental',Location='southwest')
    xlabel('$$ v = V/V_{oc} $$ ','Interpreter','latex', FontSize=14)
    ylabel('$$ i = I/I_{sc} $$ ','Interpreter','latex', FontSize=14)
    title("Comparación de métodos para " + names{i})
    xlim([0 1.1])
    ylim([0 1.1])
    hold off
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
