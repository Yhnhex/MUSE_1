clc
clear
close all

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.25) % LineWidth
set(groot,'defaultAxesFontSize',20) % Fontsize

format long 

%% Constantes 

Rt = 6378;         % Radio de la Tierra [km]
mu = 3.987*10^5; % Para operarse con kilómetros
J2 = 1.0827*10^(-3);
Rendimiento = 0.3;
A = 0.226*0.34*3; % Área de los tres paneles forntales
Al = 0.226*0.226; % Área del panel lateral

n = 2 ;                      % NÚMERO DE ÓRBITAS QUE QUIERES EVALUAR

f = 0.5;                     % Factor de empaquetamiento, CAMBIAR A MANO COMO MACACOS
%% 

%%%%% Parámetros de la órbita %%%%%

hs = (500:50:650);   % Altitud [km]
a = Rt + hs;         % SMA [km]

alpha = linspace(0,360,1000);                % Ángulo del satélite en la órbita
alpha = repmat(alpha,1,n);                   % Reescribo para el nº de vueltas que queramos
alpha_grafica = linspace(0,360*n,1000*n);    % Para las gráficas    
betas = [10 20 30 45];                       % Distintas betas que nos piden simular
phis = linspace(0,330,12);                  % Giros con respecto del eje X del satélite

alpha_eclipse = asind(Rt./(a));              % Ángulo desde la horizontal hasta la intersección del radiovector desde el centro de la órbita con el inicio del cilindro de su sombra
T = 2*pi*sqrt(a.^(3)/mu);                    % Periodos de cada órbita, dependen sólo de la altura (no como el eclipse, por lo que el % de tiempo en sombra será distinto según h y beta)  

%%%%% Primer caso: -X apuntando siempre a Nadir y sólo panel frontal %%%%%

P_O1_SP_tensor = zeros(length(betas),length(alpha),length(hs)); % Creo tensor: cada piso (i) representa una altura, cada fila (j) de las matrices una beta y cada columna un angulo alpha (k) en la orbita
for i = 1:length(hs)      
    figure('Name',append('P vs alpha O1SP h = ', num2str(hs(i)) , ' km'))
    title ('Potencia instantanea con eje $$-X$$ apuntando a Nadir y panel frontal')
    hold on
    
    grid on
    for j = 1:length(betas)
        P_O1_SP = 1367*Rendimiento*cosd(alpha)*f*A*cosd(betas(j));              % El coseno del ángulo es cos(alpha)*cos(beta), lo podemos sacar así o con la función
        P_O1_SP(P_O1_SP<0)=0;
        P_O1_SP_tensor(j,:,i) = P_O1_SP;
        plot(alpha_grafica,P_O1_SP,LineWidth=2)
    end 
    
    xlabel('$\alpha$ [$^{o}$]')
    ylabel('$P$ [W]')
    legend('$\beta = 10^{o}$','$\beta = 20^{o}$','$\beta = 30^{o}$','$\beta = 45^{o}$')
    hold off 
end 

figure('Name','Comprobación de chimichangas')

hold on
sgtitle ('Incidencia con eje $$-X$$ apuntando a Nadir')
grid on
% joputa esta se pone (incidencia y GMAt (por lo s cojones de ines)
for i = 1:length(betas)
    cosTheta = zeros(1,length(alpha));
    Bufarras = zeros(1,length(alpha));
    for j = 1:length(alpha)
        cosTheta(j) = NormalFrontRotation(alpha(j),betas(i));
        Bufarras(j) = acosd(cosTheta(j));
    end 
    %cosTheta(cosTheta<0)=0;
    subplot(1,2,1)
    hold on
    grid on
    plot(alpha_grafica,cosTheta,LineWidth=2)
    xlabel('$$\alpha$$ [$^{o}$]','Interpreter','latex')
        ylabel('cos($$\theta$$)', 'Interpreter','latex')
    xlim([0 360])

    hold off

    subplot(1,2,2)
    hold on
    grid on
    plot(alpha_grafica,Bufarras,LineWidth=2)
    xlabel('$$\alpha$$ [$^{o}$]','Interpreter','latex')
    ylabel('$$\theta$$ [$^{o}$]', 'Interpreter','latex')
    xlim([0 360])
    hold off
end 
legend('$\beta = 10^{o}$','$\beta = 20^{o}$','$\beta = 30^{o}$','$\beta = 45^{o}$')
hold off

%%%%% Segundo caso: siempre orientado hacia el Sol y sólo panel frontal %%%%%

theta_eclipse = zeros(length(hs),length(betas));            % Matriz con todo el ángulo que va a ocupar el satélite. Cada fila es una altura, cada columna una beta
for i = 1:length(hs)
    for j = 1:length(betas)
        if betas(j) < alpha_eclipse(i)
            theta_eclipse(i,j) = 2*acosd(cosd(alpha_eclipse(i))/cosd(betas(j)));
        end 
    end 
end 
t_eclipse = theta_eclipse.*T./360;                          % Tiempo de eclipse
t_eclipse_relativo = t_eclipse./T;                          % Porcentaje del tiempo de la órbita en la que se tendrá un eclipse
alpha_inicio_eclipse = 180 - theta_eclipse/2;               % Ángulo de la órbita en la que empieza el eclipse

P_O2_SP_tensor = zeros(length(betas),length(alpha),length(hs));
P_O2_SP = zeros(1,length(alpha));

for i = 1:length(hs)
    % joputa otro par de ellas tonto
    figure('Name',append('P vs alpha O2SP h = ', num2str(hs(i)) , ' km'))
        title ('Potencia instantanea con eje $$-X$$ apuntando al Sol y panel frontal', 'Interpreter','latex')

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
    xlabel('$\alpha$ [$^{o}$]')
    ylabel('$P$ [W]')
    legend('$\beta = 10^{o}$','$\beta = 20^{o}$','$\beta = 30^{o}$','$\beta = 45^{o}$')
    hold off
end 

for j = 1:length(betas)
    % joputa otro par de ellas tonto
    figure('Name',append('P vs alpha O2SP Beta = ', num2str(betas(j))))
    title ('Potencia instantanea con eje $$-X$$ apuntando al Sol y panel frontal')

    hold on
    grid on

    for i = 1:length(hs)
        plot(alpha_grafica,P_O2_SP_tensor(j,:,i),LineWidth=2)
    end 
    xlabel('$\alpha$ [$^{o}$]')
    ylabel('$P$ [W]')
    legend('$h = 500$ km','$h = 550$ km','$h = 600$ km','$h = 650$ km')
    hold off
end

%%%%% Tercer caso: -X orientado a nadir y panel lateral %%%%%

P_media_O1_lateral = zeros(length(hs),length(betas));
for i = 1:length(hs)

    % joputa otro par de ellas tonto
    figure('Name',append('P vs alpha O1PL h = ', num2str(hs(i)) , ' km'))
    title (append('Potencia instantanea con eje $$-X$$ apuntando a Nadir y paneles frontal y lateral para $h =$ ', num2str(hs(i)) , ' km'))

    hold on
    grid on
    for j = 1:length(betas)
        P_O1_PL = 1367*Rendimiento*cosd(alpha + 90)*f*Al*cosd(betas(j));
        for k = 1:length(alpha)
            if (alpha(k) > 0) && (alpha(k) < (alpha_inicio_eclipse(i,j) + theta_eclipse(i,j)))
            P_O1_PL(k) = 0;
            end
        end
        P_media_O1_lateral(i, j) = mean(P_O1_PL);
        P_O1_PL_tot = P_O1_PL + P_O1_SP_tensor(j,:,i);
        plot(alpha_grafica,P_O1_PL_tot,LineWidth=2)
    end
    xlabel('$\alpha$ [$^{o}$]')
    ylabel('$P$ [W]')
    legend('$\beta = 10^{o}$','$\beta = 20^{o}$','$\beta = 30^{o}$','$\beta = 45^{o}$')
    hold off
end 

%%%%% Cuarto caso: -X orientado a nadir, panel lateral girado varios angulitos %%%%%           

for i = 1:length(hs)
    for j = 1:length(betas)

        beta = betas(j);

        for m = 1:length(phis)
        
        % Potencia generada por el panel lateral cuando está girado en torno al eje z
        P_O3_lateral = 1367*Rendimiento*Al*f*NormalLatRotation(alpha, beta, phis(m));
        P_O3_lateral(P_O3_lateral<0) = 0;

        for k = 1:length(alpha) 
            if (alpha(k) > alpha_inicio_eclipse(i,j)) && (alpha(k) < (alpha_inicio_eclipse(i,j) + theta_eclipse(i,j)))
            P_O3_lateral(k) = 0;
            end
        end

        P_O3_lateral_tensor(j, :, i, m) = P_O3_lateral; % Mismo criterio de definición que para los otros tensores
        P_media_O3_lateral(j, i, m) = mean(P_O3_lateral);
        end
    end 
end

for j = 1:length(betas)
    % joputa DE ESTAS LA   DE MAYOR Y MENOR ALTURA pero mirando a que
    % altura las metes
    figure('Name','Potencia(alpha) para un beta y una altura dadas')
    title( ['Potencia instantanea para $h$ = ', num2str( hs(i) ), ' km; $\beta = $', num2str( betas(j) ),'$^{o}$ y varios  $\varphi$'] )
    hold on
    grid on

    for m = 1:5
        plot(alpha_grafica,P_O3_lateral_tensor(j,:,i,m))
    end 
    xlabel('$\alpha$ [$^{o}$]')
    ylabel('$P$ [W]')
    legend('$\varphi = 0^{o}$', '$\varphi = 30^{o}$','$\varphi = 60^{o}$','$\varphi = 90^{o}$','$\varphi = 120^{o}$')
    hold off
end

% plot(alpha_grafica,P_O3_lateral_tensor(j,:,i,end-2)) Para phi = 270º el panel
% lateral no ve nada de luz

for j = 1:length(betas)
    % de estas otro par de ellas joputa 
    figure('Name','Potencia(alpha) para un beta y una altura dadas con el PANEL FRONTAL')
    title( ['Potencia instantanea con eje $$-X$$ apuntando a Nadir y panel frontal para $h$ = ', num2str( hs(i) ), ' km; $\beta = $', num2str( betas(j) ),'$^{o}$'] )
    hold on
    grid on

    for m = 1:length(phis(1:4))
        plot(alpha_grafica,P_O1_SP_tensor(j,:,i) + P_O3_lateral_tensor(j,:,i,m))
    end 
    legend('$\varphi = 0^{o}$', '$\varphi = 30^{o}$','$\varphi = 60^{o}$','$\varphi = 90^{o}$')
    xlabel('$\alpha$ [$^{o}$]')
    ylabel('$P$ [W]')
    hold off
end

%% Movidas de eclipses y potencias medias, creo cosas grandes para que los surf tengan mazo puntos 

alpha = linspace(0,360,1000);               % COMO ALTERNATIVA SE PUEDEN DIVIDIR LOS TIEMPOS DE ECLIPSE Y POTENCIAS MEDIAS ENTRE N

hs_grande = linspace(500,650,50);
betas_grande = linspace(0,45,50);
phis_grande = linspace(0,360,50);
a_grande = Rt + hs_grande;                  % SMA
T_grande = 2*pi*sqrt(a_grande.^(3)/mu);    
alpha_eclipse_grande = asind(Rt./(a_grande));

theta_eclipse_grande = zeros(length(hs_grande),length(betas_grande));
alpha_inicio_eclipse_grande = 180 - theta_eclipse_grande/2;

for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)
        if betas_grande(j) < alpha_eclipse_grande(i)
        theta_eclipse_grande(i,j) = 2*acosd(cosd(alpha_eclipse_grande(i))/cosd(betas_grande(j)));
        else
        theta_eclipse_grande(i,j) = 0;
        end 
    end 
end 

%%% TIEMPO DE ECLIPSE EN FUNCIÓN DEL BETA Y DE LA ALTURA %%%
t_eclipse_grande = theta_eclipse_grande.*T_grande./360;
t_eclipse_relativo = t_eclipse_grande./T_grande*100;
figure('Name','Tiempos de eclipse')
hold on
sgtitle('Tiempos de eclipse')

subplot(1,2,1)
hold on
grid on
surf(betas_grande,hs_grande,t_eclipse_grande)
%surf(betas_grande,hs_grande,t_eclipse_relativo )
view(45,45)
xlabel('$\beta$ [$^{o}$]')
ylabel('$h$ [km]')
zlabel('$t$ [s]')


%% shading interp
hold off
subplot(1,2,2)
hold on
grid on
[C,h] = contour(betas_grande,hs_grande,t_eclipse_grande,15,'ShowText',true,"LabelFormat","%0.1f s",LineWidth = 3)
%contour(betas_grande,hs_grande,t_eclipse_relativo,15,'ShowText',true,"LabelFormat","%0.1f ",LineWidth = 2)
xlabel('$\beta$ [$^o$]')
ylabel('$h$ [km]')
clabel(C,h,'FontSize',15,'Color','black','labelspacing', 400)
xlabel('$\beta$ [$^{o}$]')
ylabel('$h$ [km]')
zlabel('$t$ [s]')
hold off
hold off

%%% POTENCIA MEDIA DE LAS ÓRBITAS PARA PANEL FRONTAL Y -X SIEMPRE APUNTANDO A NADIR %%% 

% joputa de etstas se ponen una o dos y fin OJO, cAMBIAR LO DE RELATIVO A
% SEGUNDOS Y COMPARAR CON RESULTADOS GMAT
p_media_O1_SP = zeros(length(hs_grande),length(betas_grande));
for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)
        P_O1_SP_grande = 1367*Rendimiento*cosd(alpha)*f*A*cosd(betas_grande(j));
        P_O1_SP_grande(P_O1_SP_grande<0)=0;
        p_media_O1_SP(i,j)= mean(P_O1_SP_grande);
    end 
end 
% joputa esta tb, es la de panel mirando a nadir
figure('Name','Potencia media dependiendo de Beta y Altitud para caso O1SP')
sgtitle('Potencia media con eje $$-X$$ apuntando a Nadir y panel frontal')
subplot(1,2,1)
hold on
grid on
surf(betas_grande,hs_grande,p_media_O1_SP)
view(45,45)
% shading interp
xlabel('$\beta$ [$^o$]')
ylabel('$h$ [km]')
zlabel('$P$ [W]')
hold off
subplot(1,2,2)
hold on
grid on
contour(betas_grande,hs_grande,p_media_O1_SP,15,'ShowText',true,"LabelFormat","%0.1f W",LineWidth = 3)
xlabel('$\beta$ [$^o$]')
ylabel('$h$ [km]')
zlabel('$P$ [W]')
hold off


%%% POTENCIA MEDIA DE LAS ÓRBITAS PARA PANEL FRONTAL Y ORIENTACIÓN PERFECA %%%
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
% joputa esta si pero en complementasion con los t de eclipses de este caso
figure('Name','Potencia media dependiendo de Beta y Altitud para caso O2SP')
hold on
sgtitle('Potencia media con eje $$-X$$ apuntando al Sol y panel frontal')

% subplot(1,3,1)
% surf(betas_grande,hs_grande,p_media_O2_SP)
% view(0,90)
% box;
% col = colorbar;
% % shading interp
subplot(1,2,1)
hold on
grid on
view(45,45)
surf(betas_grande,hs_grande,p_media_O2_SP)
xlabel('$\beta$ [$^o$]')
ylabel('$h$ [km]')
zlabel('$P$ [W]')
% shading interp
% set(gca,'Position',[0.05 0.1 0.4 0.8]) % Ajusta el tamaño y posición del primer "subplot"
% set(gca,'Position',[0.55 0.1 0.4 0.8]) % Ajusta el tamaño y posición del segundo "subplot"
subplot(1,2,2)
contour(hs_grande,betas_grande,p_media_O2_SP,20,'ShowText',true,"LabelFormat","%0.1f W",LineWidth = 3)
ylabel('$\beta$ [$^o$]')
xlabel('$h$ [km]')
zlabel('$P$ [W]')
hold off
%%% POTENCIA MEDIA DE LAS ÓRBITAS PARA 2 PANELES Y -X SIEMPRE APUNTANDO A NADIR %%%

P_O3_lateral_tensor_grande = zeros (length(betas_grande),length(alpha),length(hs),length(phis_grande));
P_O3_frontal_tensor_grande = zeros (length(betas_grande),length(alpha),length(hs),length(phis_grande));
P_media_O3_lateral_grande = zeros(length(betas_grande),length(hs), length(phis_grande));
P_media_O3_grande = zeros(length(betas_grande),length(hs), length(phis_grande)); 


for i = 1:length(hs)         %% BUCLE PARA COMPARAR , NO HACE FALTA hs_grande %%
    for j = 1:length(betas_grande)

        beta = betas_grande(j);

        for m = 1:length(phis_grande)
        
        % Potencia generada por el panel lateral cuando está girado en torno al eje z
        P_O3_lateral_grande = 1367*Rendimiento*Al*f*NormalLatRotation(alpha, beta, phis_grande(m));
        P_O3_frontal_grande = 1367*Rendimiento*A*f*NormalFrontRotation(alpha, beta);
        P_O3_lateral_grande(P_O3_lateral_grande<0) = 0;
        P_O3_frontal_grande(P_O3_frontal_grande<0) = 0;
 
        for k = 1:length(alpha) 
            if (alpha(k) > alpha_inicio_eclipse_grande(i,j)) && (alpha(k) < (alpha_inicio_eclipse_grande(i,j) + theta_eclipse_grande(i,j)))
            P_O3_lateral_grande(k) = 0;
            P_O3_frontal_grande(k) = 0;
            end
        end

        P_O3_lateral_tensor_grande(j, :, i, m) = P_O3_lateral_grande; % Mismo criterio de definición que para los otros tensores
        P_O3_frontal_tensor_grande(j, :, i, m) = P_O3_frontal_grande;

        P_media_O3_lateral_grande(j, i, m) = mean(P_O3_lateral_grande);
        P_media_O3_grande(j, i, m) = mean(P_O3_lateral_grande) + mean(P_O3_frontal_grande);
        end
    end 
end

Resta = zeros(length(betas_grande),length(phis_grande));

for i = 1:length(hs)
    for m = 1:length(phis_grande)
        Resta(:,m) = P_media_O3_lateral_grande(:,i,m) - P_media_O3_lateral_grande(:,i,1);
    end 
% ETS TB joputa PERO SIN LA RESTA HAY QUE EDITARLA MAMAWEBO
    figure('Name','COMPARASION WAPARDA')
    hold on

    grid on
    surf(phis_grande,betas_grande, Resta) 
    view(45,45)
    % shading interp
    title( ['Diferencia entre la potencia instantanea y media para $h$ = ', num2str( hs(i) ), ' km'] )
    ylabel('$\beta$ [$^o$]')
    xlabel('$\varphi$ [$^o$]')
    zlabel('$\Delta P$ [W]')
    xlim([0 360])
    ylim([0 45])
    hold off
   %% ESTA SI joputa
    figure ('Name','COMPARASION (no tan waparda)')
    hold on
    sgtitle( ['Diferencia entre la potencia instantanea y media para $h$ = ', num2str( hs(i) ), ' km'] )
    grid on  

    subplot(1,2,1)
    surf(phis_grande,betas_grande, Resta)
    view(0,90)
    box;
    col = colorbar;
    % shading interp
    ylabel('$\beta$ [$^o$]')
    xlabel('$\varphi$ [$^o$]')
    xlim([0 360])
    ylim([0 45])
    hold off

    subplot(1,2,2)
    hold on
    grid on
    [C,h] = contour(phis_grande,betas_grande, Resta,10,'ShowText',true,"LabelFormat","%0.1f W",LineWidth = 3);
    clabel(C,h,'FontSize',11,'Color','black')
    
    ylabel('$\beta$ [$^o$]')
    xlabel('$\varphi$ [$^o$]')

    hold off
    hold off
end 
%%
P_O3_lateral_tensor_grande_x = zeros (length(betas_grande),length(alpha),length(hs_grande),length(phis_grande));
P_O3_frontal_tensor_grande_x = zeros (length(betas_grande),length(alpha),length(hs_grande),length(phis_grande));
P_media_O3_lateral_grande_x = zeros(length(betas_grande),length(hs_grande), length(phis_grande));
P_media_O3_grande_x = zeros(length(betas_grande),length(hs_grande), length(phis_grande)); 


for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)

        beta = betas_grande(j);

        for m = 1:length(phis_grande)
        
        % Potencia generada por el panel lateral cuando está girado en torno al eje z
        P_O3_lateral_grande_x = 1367*Rendimiento*Al*f*NormalLatRotation(alpha, beta, phis_grande(m));
        P_O3_frontal_grande_x = 1367*Rendimiento*A*f*NormalFrontRotation(alpha, beta);
        P_O3_lateral_grande_x(P_O3_lateral_grande_x<0) = 0;
        P_O3_frontal_grande_x(P_O3_frontal_grande_x<0) = 0;
 
        for k = 1:length(alpha) 
            if (alpha(k) > alpha_inicio_eclipse_grande(i,j)) && (alpha(k) < (alpha_inicio_eclipse_grande(i,j) + theta_eclipse_grande(i,j)))
            P_O3_lateral_grande_x(k) = 0;
            P_O3_frontal_grande_x(k) = 0;
            end
        end

        P_O3_lateral_tensor_grande_x(j, :, i, m) = P_O3_lateral_grande_x; % Mismo criterio de definición que para los otros tensores
        P_O3_frontal_tensor_grande_x(j, :, i, m) = P_O3_frontal_grande_x;

        P_media_O3_lateral_grande_x(j, i, m) = mean(P_O3_lateral_grande_x);
        P_media_O3_grande_x(j, i, m) = mean(P_O3_lateral_grande_x) + mean(P_O3_frontal_grande_x);
        end
    end 
end


%%%%%%%%%%%%%%%%%%ESTA joputa
%%%%%%%%%%%%%%%%%%GRAFICA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%SIIIIUUUUUUUUUU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Potencia media TOTAL en función de beta y de la altitud')

sgtitle('Potencia media total con paneles frontal y lateral')
subplot(1,2,1)
hold on
grid on
surf(hs_grande,betas_grande,P_media_O3_grande_x(:,:,1))
ylabel('$\beta$ [$^o$]')
xlabel('$h$ [km]')
zlabel('$P$ [W]')
view(45,45)
hold off

subplot(1,2,2)
hold on
grid on
contour(hs_grande,betas_grande,P_media_O3_grande_x(:,:,1),15,'ShowText',true,"LabelFormat","%0.1f W",LineWidth = 3)
ylabel('$\beta$ [$^o$]')
xlabel('$h$ [km]')
zlabel('$P$ [W]')
hold off


%% 

for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)

        beta = betas_grande(j);

        for m = 1:length(phis_grande)
        
        % Potencia generada por el panel lateral cuando está girado en torno al eje z
        P_O3_lateral_grande = 1367*Rendimiento*Al*f*NormalLatRotation(alpha, beta, phis_grande(m));
        P_O3_lateral_grande(P_O3_lateral_grande<0) = 0;

        for k = 1:length(alpha) 
            if (alpha(k) > alpha_inicio_eclipse_grande(i,j)) && (alpha(k) < (alpha_inicio_eclipse_grande(i,j) + theta_eclipse_grande(i,j)))
            P_O3_lateral_grande(k) = 0;
            end
        end

        P_O3_lateral_tensor_grande(j, :, i, m) = P_O3_lateral_grande; % Mismo criterio de definición que para los otros tensores
        P_media_O3_lateral_grande(j, i, m) = mean(P_O3_lateral_grande);
        end
    end 
end


phiphis = zeros(length(hs_grande),length(betas_grande));
figure('Name','Ssssuuu')
hold on

grid on
for i = 1:length(hs_grande)
    for j = 1:length(betas_grande)
    [x,xindx] = max(P_media_O3_lateral_grande(j,i,:));
    phiphis(i,j) = phis_grande(xindx);
    end 
end 
surf(betas_grande,hs_grande,phiphis)
xlabel('$\beta$ [$^o$]')
ylabel('$h$ [km]')
zlabel('$\phi$ [$^o$]')
view(45,45)
%% shading interp
hold off

%%%%%%%%%%%%%%% Matrices de rotación %%%%%%%%%%%%%%% 



function cosThetaFront = NormalFrontRotation(alpha,beta)
cosThetaFront =  cosd(alpha)*cosd(beta); 
%                 cos(0.0175*alpha)*sin(0.0175*beta);
%                 -sin(0.0175*alpha)];
end

function cosThetaLat = NormalLatRotation(alpha,beta,phi)
cosThetaLat = (sind(phi)*sind(beta) + cosd(phi)*sind(alpha)*cosd(beta));
              %cos(0.0175*phi)*sin(0.0175*alpha)*sin(0.0175*beta) - sin(0.0175*phi)*cos(0.0175*beta);
              %                                      cos(0.0175*alpha)*cos(0.0175*phi)];  
end 
