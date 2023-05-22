clc
clear
close all

% Para personalizar las grÃ¡ficas
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.25) % LineWidth
set(groot,'defaultAxesFontSize',15) % Fontsize


A=xlsread("xpxT.xlsx");
c = ['#3498DB'; "#EFB810"; '#C0392B'; '#28B463' ;'#8E44AD'];
b = [1, 5, 9, 13, 17, 21];
xp = A(1:end-3,b);
p = A(1:end-3, b+1);
xT = A(1:end-3,b+2);
T = A(1:end-3, b+3);


figure(1)
hold on 
grid on
scatter(xT(:,1), T(:, 1) , "filled");
scatter(xT(:,2), T(:, 2), "filled");
%scatter(xT(:,3), T(:, 3), "filled");
scatter(xT(:,4), T(:, 4), "filled");
scatter(xT(:,5), T(:, 5), "filled");
scatter(xT(:,6), T(:, 6), "filled");
xlabel('$$x$$ [m]')
ylabel('$$T$$ [$$^{\circ}$$C]')
legend("$$\alpha = 0^{\circ}$$", "$$\alpha = 10^{\circ}$$", "$$\alpha = 30^{\circ}$$", "$$\alpha = 40^{\circ}$$", "$$\alpha = 50^{\circ}$$")
hold off

figure(2)
hold on 
grid on
scatter(xp(:,1), p(:, 1), "filled");
scatter(xp(:,2), p(:, 2), "filled");
scatter(xp(:,4), p(:, 4), "filled");
scatter(xp(:,5), p(:, 5), "filled");
scatter(xp(:,6), p(:, 6), "filled");
xlabel('$$x$$ [m]')
ylabel('$$p$$ [Pa]')
legend("$$\alpha = 0^{\circ}$$", "$$\alpha = 10^{\circ}$$", "$$\alpha = 30^{\circ}$$", "$$\alpha = 40^{\circ}$$", "$$\alpha = 50^{\circ}$$")


hold off
