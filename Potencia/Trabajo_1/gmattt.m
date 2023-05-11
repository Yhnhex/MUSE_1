clear all

close all
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.25) % LineWidth
set(groot,'defaultAxesFontSize',20) % Fontsize

% c = ['#3498DB'; "#EFB810"; '#C0392B'; '#28B463' ;'#8E44AD'];

%newcolours = ['#000080'; "#EFB810"; '#C0392B'; '#28B463'; '#3498DB'; '#8E44AD'; '#FFA07A'; '#D748DF'; '#59D13D'; '#020202' ];
%colororder(newcolours)
h500 = readmatrix('GMAT\h_500\Incidencias.txt','NumHeaderLines',1);
h550 = readmatrix('GMAT\h_550\Incidencias.txt','NumHeaderLines',1);
h600 = readmatrix('GMAT\h_600\Incidencias.txt','NumHeaderLines',1);
h650 = readmatrix('GMAT\h_650\Incidencias.txt','NumHeaderLines',1);
hs(:,:,1) = h500;
hs(:,:,2) = h550;
hs(:,:,3) = h600;
hs(:,:,4) = h650;
betas=[0 10 20 30 45];
hj=[500 550 600 650];
for H = 1:size(hs, 3)
    alpha = hs(:, 1, H);

    figure(H)
    hold on
    for i = 2:6
        subplot(1,2,1)
        hold on
        grid on
        xlim([0, 360])
        plot(alpha, cosd(hs(:, i, H)))
        legend('$$\beta = 0^{o}$$', '$$\beta = 10^{o}$$', '$$\beta = 20^{o}$$', '$$\beta = 30^{o}$$', '$$\beta = 45^{o}$$', 'Interpreter','latex')
        xlabel('$$\alpha$$ [$^{o}$]','Interpreter','latex')
        ylabel('cos($$\theta$$)', 'Interpreter','latex')
        hold off

        subplot(1,2,2)
        hold on
        grid on
        xlim([0, 360])
        plot(alpha, (hs(:, i, H)))
        legend('$$\beta = 0^{o}$$', '$$\beta = 10^{o}$$', '$$\beta = 20^{o}$$', '$$\beta = 30^{o}$$', '$$\beta = 45^{o}$$', 'Interpreter','latex')
        xlabel('$$\alpha$$ [$^{o}$]','Interpreter','latex')
        ylabel('$$\theta$$ [$^{o}$]', 'Interpreter','latex')
        hold off
            sgtitle(append('Incidencia con eje $-X$ apuntando a Nadir con GMAT para $$h=$$ ', num2str(hj(H)),' km'))

        hold off
    end 
    
end 

t_eclipses=[2139.359 2131.545 2115.087 2118.841; 2128.083 2112.155 2107.321 2105.594; 2092.543 2081.899 2072.306 2063.621; 2023.792 2008.843 1994.982 1982.052; 1830.289 1802.465 1775.7 1749.796]; 


figure('Name','Tiempos de eclipse')
sgtitle('Tiempos de eclipse en funcion de $\beta$ y la altitud $h$')

subplot(1,2,1)
hold on
grid on
%surf(betas_grande,hs_grande,t_eclipse_grande)
surf(betas, hj ,t_eclipses' )
view(45,45)
xlabel('$\beta$ [$^{o}$]')
ylabel('$h$ [km]')
zlabel('$t$ [s]')
%shading interp
hold off
subplot(1,2,2)
hold on
grid on
%contour(betas_grande,hs_grande,t_eclipse_grande,15,'ShowText',true,"LabelFormat","%0.1f s",LineWidth = 2)

[C,hj] = contour(betas, hj,t_eclipses',15,'ShowText',true,"LabelFormat","%0.1f s",LineWidth = 3)
%contour(betas_grande,hs_grande,t_eclipse_relativo,15,'ShowText',true,"LabelFormat","%0.1f ",LineWidth = 2)
xlabel('$\beta$ [$^o$]')
ylabel('$h$ [km]')
clabel(C,hj,'FontSize',15,'Color','black','labelspacing', 400)

hold off
