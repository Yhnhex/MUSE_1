clc
clear
close all

%% TRABAJO DE POTENSIA %% 

Rt = 6378;         % Radio de la Tierra
mu = 3.987*10^5;
J2 = 1.0827*10^(-3);
Rendimiento = 0.3;
A = 0.226*0.34*3;
Al = 0.226*0.226;

n = 3 ; % NÚMERO DE ÓRBITAS QUE QUIERES EVALUAR

f = 0.9;    % Factor de empaquetamiento, CAMBIAR A MANO COMO MACACOS

% Órbita

hs = (512);   % Altitud
a = Rt + hs;         % SMA

alpha = linspace(0,360,1000);
alpha = repmat(alpha,1,n);
alpha_grafica = linspace(0,360*n,1000*n);
betas = [10 20 30 45];

alpha_eclipse = asind(Rt./(a));
T = 2*pi*sqrt(a.^(3)/mu);    % PERIODOS  

% Primer caso: -X apuntando siempre a Nadir y sólo panel frontal

P_O1_SP_tensor = zeros(length(betas),length(alpha),length(hs)); %Creo tensor: cada piso representa una altura, cada fila de las matrices una beta y cada columna un angulo alpha en la orbita
for i = 1:length(hs)
    figure('Name',append('P vs alpha O1SP h = ', num2str(hs(i)) , ' km'))
    hold on
    grid on
    for j = 1:length(betas)
        P_O1_SP = 1367*Rendimiento*cosd(alpha)*f*A*cosd(betas(j));
        P_O1_SP(P_O1_SP<0)=0;
        P_O1_SP_tensor(j,:,i) = P_O1_SP;
        plot(alpha_grafica,P_O1_SP,LineWidth=2)
    end 
    legend('10','20','30','45')
    hold off 
end 

% Segundo caso: siempre orientado hacia el Sol y sólo panel frontal

theta_eclipse = zeros(length(hs),length(betas)); % Cada fila es una altura, cada columna una beta
for i = 1:length(hs)
    for j = 1:length(betas)
        if betas(j) < alpha_eclipse(i)
            theta_eclipse(i,j) = 2*acosd(cosd(alpha_eclipse(i))/cosd(betas(j)));
        end 
    end 
end 
t_eclipse = theta_eclipse.*T./360;
alpha_inicio_eclipse = 180 - theta_eclipse/2;

P_O2_SP_tensor = zeros(length(betas),length(alpha),length(hs));
P_O2_SP = zeros(1,length(alpha));
for i = 1:length(hs)
    figure('Name',append('P vs alpha O2SP h = ', num2str(hs(i)) , ' km'))
    hold on
    grid on
    for j = 1:length(betas)
        P_O2_SP(:) = 1367*Rendimiento*f*A;
        for k = 1:length(alpha)
            if (alpha(k) > alpha_inicio_eclipse(i,j)) && (alpha(k) < (alpha_inicio_eclipse(i,j) + theta_eclipse(i,j)))
            P_O2_SP(k) = 0;
            end
        end 
        P_O2_SP_tensor(j,:,i) = P_O2_SP;
        plot(alpha_grafica,P_O2_SP,LineWidth=2)
    end
    legend('10','20','30','45')
    hold off
end 

for j = 1:length(betas)
    figure('Name',append('P vs alpha O2SP Beta = ', num2str(betas(j))))
    hold on
    grid on

    for i = 1:length(hs)
        plot(alpha_grafica,P_O2_SP_tensor(j,:,i),LineWidth=2)
    end 
    legend('h = 500 km','h = 550 km','h = 600 km','h = 650 km')
    hold off
end

% Tercer caso: -X orientado a nadir y panel lateral

for i = 1:length(hs)
    figure('Name',append('P vs alpha O1PL h = ', num2str(hs(i)) , ' km'))
    hold on
    grid on
    for j = 1:length(betas)
        P_O1_PL = 1367*Rendimiento*cosd(alpha + 90)*f*Al*cosd(betas(j));
        for k = 1:length(alpha)
            if (alpha(k) > 0) && (alpha(k) < (alpha_inicio_eclipse(i,j) + theta_eclipse(i,j)))
            P_O1_PL(k) = 0;
            end
        end 
        P_O1_PL_tot = P_O1_PL + P_O1_SP_tensor(j,:,i);
        plot(alpha_grafica,P_O1_PL_tot,LineWidth=2)
    end
    legend('10','20','30','45')
    hold off
end 


%% Movidas de eclipses y potencias medias, creo cosas grandes para que los surf tengan mazo puntos 

alpha = linspace(0,360,1000); % COMO ALTERNATIVA SE PUEDEN DIVIDIR LOS TIEMPOS DE ECLIPSE Y POTENCIAS MEDIAS ENTRE N

hs_grande = linspace(500,650,50);
betas_grande = linspace(10,45,50);

a_grande = Rt + hs_grande;         % SMA
T_grande = 2*pi*sqrt(a_grande.^(3)/mu);    
alpha_eclipse_grande = asind(Rt./(a_grande));

theta_eclipse_grande = zeros(length(hs_grande),length(betas_grande));
alpha_inicio_eclipse_grande = 180 - theta_eclipse_grande/2;

for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)
        
        theta_eclipse_grande(i,j) = 2*acosd(cosd(alpha_eclipse_grande(i))/cosd(betas_grande(j)));
    
    end 
end 
t_eclipse_grande = theta_eclipse_grande.*T_grande./360;

figure('Name','PEPINO')
surf(hs_grande,betas_grande,t_eclipse_grande)
shading interp

p_media_O1_SP = zeros(length(hs_grande),length(betas_grande));
for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)
        P_O1_SP_grande = 1367*Rendimiento*cosd(alpha)*f*A*cosd(betas_grande(j));
        P_O1_SP_grande(P_O1_SP_grande<0)=0;
        p_media_O1_SP(i,j)= mean(P_O1_SP_grande);
    end 
end 

figure('Name','Potencia media dependiendo de Beta y Altitud para caso O1SP')
surf(betas_grande,hs_grande,p_media_O1_SP)
shading interp

P_O2_SP_grande = zeros(1,length(alpha));

p_media_O2_SP = zeros(length(hs_grande),length(betas_grande));
for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)
        P_O2_SP_grande(:) = 1367*Rendimiento*f*A;
        for k = 1:length(alpha)
            if (alpha(k) > alpha_inicio_eclipse_grande(i,j)) && (alpha(k) < (alpha_inicio_eclipse_grande(i,j) + theta_eclipse_grande(i,j)))
            P_O2_SP_grande(k) = 0;
            end
        end 
        p_media_O2_SP(i,j) = mean(P_O2_SP_grande);
    end
end 

figure('Name','Potencia media dependiendo de Beta y Altitud para caso O2SP')
subplot(1,3,1)
surf(hs_grande,betas_grande,p_media_O2_SP)
view(0,90)
box;
col = colorbar;
shading interp
subplot(1,3,2)
hold on
surf(hs_grande,betas_grande,p_media_O2_SP)
shading interp
set(gca,'Position',[0.05 0.1 0.4 0.8]) % Ajusta el tamaño y posición del primer "subplot"
set(gca,'Position',[0.55 0.1 0.4 0.8]) % Ajusta el tamaño y posición del segundo "subplot"
subplot(1,3,3)
contour(hs_grande,betas_grande,p_media_O2_SP,20,'ShowText',true,"LabelFormat","%0.1f W",LineWidth = 2)
