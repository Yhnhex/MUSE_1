clear all
tic
%% Definicion de variables

% GEOMETRIA
% Longitud panel [m]
Lx = 0.5;
Lz = 0.3;
% Espesor capas [m]
h_CFRP = 0.0025;
h_AlH = 0.015;
h_tot = 2*h_CFRP+h_AlH;
% Areas [m2]
A_panel = Lx*Lz;
A_heater = 0.1*0.1;
A_latx = Lz*h_tot;
A_latz = Lx*h_tot;

% MATERIALES
% Densidad [kg/m3]
rho_CFRP = 1650;
rho_AlH = 40;
% Capacidad termica [J/kg.K]
C_CFRP = 930;
C_AlH = 890;
% Conductividad termica [W/m.K]
K_CFRP = [130, 1, 130];
K_AlH = [0.7, 1.5, 0.7];
% Emisividad
eps_CFRP = 0.8;

% Heaters
Q = 50;

% Cosas
sigma = 5.6708*10^-8; % [W/m2.K4]
Tsat = 0; % [ºC]

%% Sol. analitica
% Disipacion homogenea
K_eq = (2*K_CFRP(1)*h_CFRP + (K_AlH(1)*h_AlH))/h_tot;
rho_eq = (2*rho_CFRP*h_CFRP + rho_AlH*h_AlH)/h_tot;
C_eq = (2*C_CFRP*h_CFRP + C_AlH*h_AlH)/h_tot;

syms x y z

T_anal1 = 2*Q/(K_eq*A_latx)*x*(Lx-x);

% Disipacion en heaters

M = [0    1  0     0   0     0   0     0   0     0;
     0.05 0 -0.05 -1   0     0   0     0   0     0;
     1    0 -1     0   0     0   0     0   0     0;
     0    0  0.15  1  -0.15 -1   0     0   0     0;
     0    0  1     0  -1     0   0     0   0     0;
     0    0  0     0   0.35  1  -0.35 -1   0     0;
     0    0  0     0   1     0  -1     0   0     0;
     0    0  0     0   0     0   0.45  1  -0.45 -1;
     0    0  0     0   0     0   1     0  -1     0;
     0    0  0     0   0     0   0     0   0.5   1];

q_coefs = [0 -0.05^2/2 -0.05 0.15^2/2 0.15 -0.35^2/2 -0.35 0.45^2/2 0.45 0]'*Q/(K_eq*A_latx);
coefs = M\q_coefs;

T_anal2 = piecewise(x<0.05,           coefs(1)*x ,...
                    x>=0.05 & x<0.15, coefs(3)*x + coefs(4) - x^2*Q/(2*K_eq*A_latx), ...
                    x>=0.15 & x<0.35, coefs(5)*x + coefs(6), ...
                    x>=0.35 & x<0.45, coefs(7)*x + coefs(8) - x^2*Q/(2*K_eq*A_latx), ...
                    x>=0.45 & x<=0.5, coefs(9)*x + coefs(10))*10;


%% Sol. numerica

nx = 500;
dx = Lx/nx;
titer = 150000;
dt = 0.1;

mu = dt/(rho_eq*C_eq);

a = 0;
cont = 0;
time_on = 0;
bucle = false;
conv_mean = false;
conv_RMS = false;

xplot = linspace(0,Lx,nx);
tplot = linspace(1,titer*dt,titer);

T1 = zeros(nx,1);
T2 = zeros(nx,1);
T3 = zeros(nx,1);
T3opt = zeros(nx,1);
T5 = zeros(nx,1);

T5mid = zeros(length(tplot),1);
diff_mean = zeros(length(tplot),1);
diff_RMS = zeros(length(tplot),1);

for t = 1:titer
    
    Aux = T3;

    for i = 2:nx-1

        if (i>=nx/10)&&(i<=3*nx/10) || (i>=7*nx/10)&&(i<=9*nx/10)
            Qnod = Q/(A_latx*(nx/5));
        else 
            Qnod = 0;
        end

        Qradopt = eps_CFRP*sigma*((T3(i)+273.15)^4 - (Tsat+273.15)^4);
%         Qrad = eps_CFRP*sigma*((T3(nx/2)+273.15)^4 - (Tsat+273.15)^4);
        Qrad = 4*eps_CFRP*sigma*(Tsat+273.15)^3*(T3(nx/2) - Tsat);
        Qrad5 = eps_CFRP*sigma*((T5(nx/2)+273.15)^4 - (Tsat+273.15)^4);

        T1(i) = (T1(i+1) + T1(i-1))/2 + dx^2*Q/(K_eq*A_latx*(nx-1));
        T2(i) = (T2(i+1) + T2(i-1))/2 + dx^2*Qnod/K_eq;
        
        T3(i) = (T3(i+1) + T3(i-1))/2 + dx^2*(Qnod/K_eq - Qrad);
        T3opt(i) = (T3opt(i+1) + T3opt(i-1))/2 + dx^2*(Qnod/K_eq - Qradopt);

        if T5(nx/2)<10/nx
            heater_on = true;
        elseif T5(nx/2)>12/nx
            heater_on = false;
        end
        
        if heater_on == true
            T5(i) = (T5(i+1) + T5(i-1))/2 + dx^2*(Qnod/K_eq - Qrad);
        else
            T5(i) = (T5(i+1) + T5(i-1))/2 - dx^2*(Qrad);
        end

        T1(1) = 0; T1(end) = 0;
        T2(1) = 0; T2(end) = 0;
        T3(1) = 0; T3(end) = 0;
    end
    
    if (heater_on == false) && (cont == 0)
        bucle = true;
    end
    
    if bucle == true
        cont = cont + 1;
        if heater_on == true
            time_on = time_on + 1;
        end
    end

    diff_mean(t) = mean(T3 - Aux)*100;
    if (diff_mean(t)/diff_mean(1) < 0.1) && (conv_mean == false)
        t_conv_mean = t;
        conv_mean = true;
    end

    diff_RMS(t) = rms(T3 - Aux)*100;
    if (diff_RMS(t)/diff_RMS(1) < 0.1) && (conv_RMS == false)
        t_conv_RMS = t;
        conv_RMS = true;
    end

    T5mid(t) = T5(nx/2);
    
    if bucle == true
        if (T5mid(t) < T5mid(t-1)) && (T5mid(t-1) > T5mid(t-2))
            a = a + 1;
            t_loop(a) = t;
            if a>1
                ciclo{a} = T5mid(t_loop(a-1):t_loop(a));
            end
        end
    end
end

T1 = T1*2*nx;
T2 = T2*nx;
T3 = T3*nx;
T3opt = T3opt*nx;
T3diff = (T3opt - T3)./T3opt*100;
T5mid = T5mid*nx;

T_anal1_ev = double(subs(T_anal1, x, xplot));
T_anal2_ev = double(subs(T_anal2, x, xplot));

Err1 = (T_anal1_ev - T1')./T_anal1_ev*100 ;
Err2 = (T_anal2_ev - T2')./T_anal2_ev*100 ;

time_on_med = time_on/cont; 
Pmed = Q*time_on_med;





% nx = 51; % número de nodos en el eje x
% ny = 17; % número de nodos en el eje y
% nz = 31; % número de nodos en el eje z
% dx = Lx/(nx-1); % tamaño del paso en x
% dy = h_tot/(ny-1); % tamaño del paso en y
% dz = Lz/(nz-1); % tamaño del paso en z
% dt = 0.1; % tamaño del paso de tiempo
% 
% % Inicializar la temperatura en la malla
% T = zeros(nx, ny, nz);
% T(1, :, :) = 0; % pared inferior a temperatura constante
% T(end, :, :) = 0; % pared superior a temperatura constante
% 
% % Iterar la solución en el tiempo
% for t = 1:100
%     % Calcular la temperatura en el interior de la malla
%     for i = 2:nx-1
%         for j = 1:ny
%             if j*dy <= h_CFRP 
%                alpha_x = K_CFRP(1);
%                alpha_y = K_CFRP(2);
%                alpha_z = K_CFRP(3);
%             elseif (j*dy >= h_CFRP) && (j*dy <= h_CFRP+h_AlH)
%                 alpha_x = K_AlH(1);
%                 alpha_y = K_AlH(2);
%                 alpha_z = K_AlH(3);
%             else
%                 alpha_x = K_CFRP(1);
%                 alpha_y = K_CFRP(2);
%                 alpha_z = K_CFRP(3);
%             end
%             for k = 1:nz
%                 if j == 1
%                     sumy = (T(i, j+1, k) - T(i, j, k))/dy;
%                 elseif j == ny
%                     sumy = (T(i, j, k) - T(i, j-1, k))/dy;
%                 else
%                     sumy = (T(i, j+1, k) - 2*T(i, j, k) + T(i, j-1, k))/(dy^2);
%                 end 
% 
%                 if k == 1
%                     sumz = (T(i, j, k+1) - T(i, j, k))/dy;
%                 elseif k == nz
%                     sumz = (T(i, j, k) - T(i, j, k-1))/dy;
%                 else
%                     sumz = (T(i, j, k+1) - 2*T(i, j, k) + T(i, j, k-1))/(dy^2);
%                 end 
% 
%                 if (i>=5)&&(i<=15)&&(j==ny)&&(k>=5)&&(k<=15)
%                     Qnod = Q/100;
%                 elseif (i>=35)&&(i<=45)&&(j==ny)&&(k>=5)&&(k<=15)
%                     Qnod = Q/100;
%                 else
%                     Qnod = 0;
%                 end
% 
%                 T(i, j, k) = T(i, j, k) + ...
%                     alpha_x*dt*(T(i+1, j, k) - 2*T(i, j, k) + T(i-1, j, k))/(dx^2) + ...
%                     alpha_y*dt*sumy + alpha_z*dt*sumz + Qnod;
%             end
%         end
%     end
%     
%     % Aplicar condiciones de frontera
%     T(1, :, :) = 0; % pared inferior a temperatura constante
%     T(end, :, :) = 0; % pared superior a temperatura constante
%     
%     % Visualizar la solución en un corte transversal
% 
% end
% 
% x = linspace(0, 1, nx);
% y = linspace(0, 1, ny);
% z = linspace(0, 1, nz);
% [X, Y, Z] = meshgrid(x, y, z);
% 
% % Definir las coordenadas de la sección que se desea visualizar
% xi = [0 0.5];
% yi = [0 0.02];
% zi = [0 0.3];
% 
% % Visualizar la sección de la matriz en un mapa de calor
% figure
% slice(X, Y, Z, T, xi, yi, zi)
% colormap(jet)
% colorbar
% xlabel('x')
% ylabel('y')
% zlabel('z')


%% Figuras
close all

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.5) % LineWidth
set(groot,'defaultAxesFontSize',15) % Fontsize


figure(1)
fplot(T_anal1,[0,Lx])
hold on
plot(xplot,T1,LineStyle="--")
title('Perfil de temperatura para disipaci\''on homog\''enea')
xlabel('Longitud del panel [m]')
ylabel('Temperatura [$^\circ$ C]')
legend('Anal\''itica','Numerica')
set(gcf,'Position',[100 100 1000 500])
grid on

figure(2)
plot(xplot,Err1)
title('Error entre soluci\''on anal\''itica y num\''erica para disipaci\''on homog\''enea')
xlabel('Longitud del panel [m]')
ylabel('Error [\%]')
set(gcf,'Position',[100 100 1000 500])
grid on

figure(3)
fplot(T_anal2,[0,Lx])
hold on
plot(xplot,T2,LineStyle="--")
title('Perfil de temperatura para disipaci\''on en heaters')
xlabel('Longitud del panel [m]')
ylabel('Temperatura [$^\circ$ C]')
legend('Anal\''itica','Num\''erica')
set(gcf,'Position',[100 100 1000 500])
grid on

figure(4)
plot(xplot,Err2)
title('Diferencia entre soluci\''on anal\''itica y num\''erica para disipaci\''on en heaters')
xlabel('Longitud del panel [m]')
ylabel('Error [\%]')
set(gcf,'Position',[100 100 1000 500])
grid on

figure(5)
plot(xplot,T3)
hold on
plot(xplot,T3opt,LineStyle="--")
title('Perfil de temperatura con radiaci\''on')
xlabel('Longitud del panel [m]')
ylabel('Temperatura [$^\circ$ C]')
legend('Radiaci\''on linealizada','Radiaci\''on sin linealizar')
set(gcf,'Position',[100 100 1000 500])
grid on

figure(6)
plot(xplot,T3diff)
title('Error cometido por linealizar la radiaci\''on')
xlabel('Longitud del panel [m]')
ylabel('Error [\%]')
set(gcf,'Position',[100 100 1000 500])
grid on

figure(8)
plot(tplot, diff_mean)
hold on
plot(tplot, diff_RMS)
title('Convergencia a estado estacionario')
xlabel('Tiempo [s]')
ylabel('Diferencia [\%]')
legend('Media aritm\''etica','RMS')
set(gcf,'Position',[100 100 1000 500])
grid on

figure(9)
plot(tplot,T5mid)
title('Evoluci\''on de la temperatura del sensor')
xlabel('Tiempo [s]')
ylabel('Temperatura [$^\circ$ C]')
set(gcf,'Position',[100 100 1000 500])
grid on

toc