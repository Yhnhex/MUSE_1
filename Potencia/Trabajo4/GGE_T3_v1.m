clear all
close all

%% Configuracion de las graficas
linewidth = 3;
fontsize = 20;
line = ["-", "--", "-.", ":", "-", "--", "-.", ":"];

%% Options ajuste curvas
opts = statset('Display','off','TolFun',1e-10);

%% Load Data
load('Data_Est');

%% Ajuste Curvas
%     b(1) = E0
%     b(2) = E1
%     b(3) = Rtd
%     b(4) = E20
%     b(5) = E30
%     b(6) = E21
%     b(7) = E22
%     b(8) = E31

% descarga 5A zona 1: lineal
X1(:,1) = D_5(100:1500, 10); % Energy discharged
X1(:,2) = D_5(100:1500, 2);  % Current
Y1 = D_5(100:1500, 3);       % Voltage

b10 = [20 -1e-6 Rd_5];
b1 = AjusteDescarga1(X1,Y1,b10, opts);

V_est1 = b1(1) + b1(2)*D_5(2:2126, 10) - D_5(2:2126, 2)*b1(3);
graphdisplay1(1, D_5(2:2126, 9), D_5(2:2126, 3), D_5(2:2126, 9), V_est1, ... 
    'V vs Energia descargada', 'Fi [Ws]', 'V [V]', ["V real", "V estimado"], linewidth, fontsize)

% descarga 5A
X2(:,1) = D_5(2:2126, 10); % Energy discharged
X2(:,2) = D_5(2:2126, 2);  % Current
Y2 = D_5(2:2126, 3);       % Voltage

b20 = [24 -1e-6 Rd_5 -1e-10 1e-5];
b2 = AjusteDescarga2(X2,Y2,b20, opts);

V_est2 = b2(1) + b2(2).*X2(:,1) + b2(4)*exp(b2(5).*X2(:,1)) - X2(:,2).*b2(3);
graphdisplay1(2, D_5(2:2126, 9), D_5(2:2126, 3), D_5(2:2126, 9), V_est2, ... 
    'V vs Energia descargada', 'Fi [Ws]', 'V [V]', ["V real", "V estimado"], linewidth, fontsize)


% descarga 5A
X3(:,1) = D_5(3:2126, 10); % Energy discharged
X3(:,2) = D_5(3:2126, 2);  % Current
Y3 = D_5(3:2126, 3);       % Voltage

b30 = [24 -1e-6 Rd_5 -1e-10 1e-5 -1e-13 -1e-12 -1e-7];
% b30 = [24 -1e-6 Rd_5 1e-5 1e-5 1e-5 1e-5 1e-5];
b3 = AjusteDescarga3(X3,Y3,b30, opts);

V_est3 = b3(1) + b3(2).*X3(:,1) + (b3(4)+b3(6)*X3(:,2)+b3(7)*X3(:,2).^2).*exp( (b3(5)+b3(8)*X3(:,2)).*X3(:,1) ) - X3(:,2).*b3(3);
graphdisplay1(3, D_5(3:2126, 9), Y3, D_5(3:2126, 9), V_est3, ... 
    'V vs Energia descargada', 'Fi [Ws]', 'V [V]', ["V bat", "V estimado"], linewidth, fontsize)






%% Functions
%Ajuste modelo estatico bateria descarga ecuacion 1
function [b_] = AjusteDescarga1(X,Y, b0, opts)
    eq1 = @(b, X) b(1) + b(2)*X(:,1) - X(:,2)*b(3);
    
    mdl = fitnlm(X, Y, eq1, b0, 'Options', opts);
    b_ = table2array(mdl.Coefficients(:,1));
end

%Ajuste modelo estatico bateria descarga ecuacion 2
function [b_] = AjusteDescarga2(X,Y,b0, opts)
    eq1 = @(b, X) b(1) + b(2).*X(:,1) + b(4)*exp(b(5).*X(:,1)) - X(:,2).*b(3);
    
    mdl = fitnlm(X, Y, eq1, b0, 'Options', opts);
    b_ = table2array(mdl.Coefficients(:,1));
end

%Ajuste modelo estatico bateria descarga ecuacion 3
function [b_] = AjusteDescarga3(X,Y, b0, opts)
    eq1 = @(b, X) b(1) + b(2).*X(:,1) + (b(4)+b(6)*X(:,2)+b(7)*X(:,2).^2).*exp( (b(5)+b(8)*X(:,2)).*X(:,1) ) - X(:,2).*b(3);
    
    mdl = fitnlm(X, Y, eq1, b0, 'Options', opts);
    b_ = table2array(mdl.Coefficients(:,1));
end









%% Graphs Functions
function graphdisplay0(f, x, y, tit, xlab, ylab, linewidth, fontsize)

    figure(f);
    hold on;
    grid on;
    plot(x, y, 'LineWidth', linewidth)
    hold off;
    xlim([min(x) max(x)])
    ylabel(ylab);
    xlabel(xlab);
    title (tit)
    set(gca, 'FontSize', fontsize);
end

function graphdisplay1(f, x1, y1, x2, y2, tit, xlab, ylab, leg, linewidth, fontsize)

    figure(f);
    hold on;
    grid on;
    plot(x1, y1, 'LineWidth', linewidth)
    plot(x2, y2, 'LineWidth', linewidth)
    hold off;
    
    ylabel(ylab);
    xlabel(xlab);
    legend(leg);
    title (tit)
    set(gca, 'FontSize', fontsize);
end


