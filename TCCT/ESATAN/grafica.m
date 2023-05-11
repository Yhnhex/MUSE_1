set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1.25) % LineWidth
set(groot,'defaultAxesFontSize',20) % Fontsize

M=csvread("Attribute_Chart_1.csv");
figure(1)
hold on
grid on
plot(M(:,1), M(:,2), Color='k', LineStyle='-', Marker='+')
plot(M(:,1), M(:,4), Color='r', LineStyle='-', Marker='o')
plot(M(:,1), M(:,6), Color='b', LineStyle='-', Marker='*')
legend('$T_{PL1}$', '$T_{PL2}$','$T_{PL3}$', 'Interpreter','latex')
xlabel('$t$ [s]', 'Interpreter','latex')
ylabel('$$T^a$$ [$^{\circ} $C]', 'Interpreter','latex')
xlim([0, 6000])
hold off
