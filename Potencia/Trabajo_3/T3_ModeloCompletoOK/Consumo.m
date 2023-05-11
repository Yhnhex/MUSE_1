%% Cargar los datos de consumo

load('Consumo.mat');

% Potencia consumida en el rail BUS.
Wbus = Cons(:, 5) + Cons(:, 11) + Cons(:, 13) + Cons(:, 15) + Cons(:, 17);
Wbusmod = Cons(:, 5) + Cons(:, 11) + Cons(:, 15) + Cons(:, 17);
Wbusnoexp = Cons(:, 5);

% Intensidad del panel solar.
Isp = Cons(:, 3);

% Vector de tiempo
T = Cons(:, 1);
dt = T(2)-T(1);

ZZ = zeros(size(Isp));


%% Datos del modelo de bateria
% Iter 1
% C1 = 93572.2724606597; % [F]
% C2 = 200.030155060056; % [F]
% R1 = 0.0018962621778753; % [R]
% R2 = 6.69645517801497e-05; % [R]
% Rc = 0.115872760572405; % [R]
% Rd = 0.103875281061481; % [R]

% Iter 2
% C1 = 2890.75147425686;
% C2 = 897.404505670572;
% R1 = 0.0182008285048376;
% R2 = 0.000547107095742223;
% Rc = 0.118792125684771;
% Rd = 0.104926383711228;

% Iter 3 (parece que la buena)
% C1 = 2432.00409623172;
% C2 = 440.18125030792;
% R1 = 0.0119470318917375;
% R2 = 0.00018588426519567;
% Rc = 0.229061190046119;
% Rd = 0.103770719943725;

% Iter 4
C1 = 2037.07956638245;
C2 = 282.027037851405;
R1 = 0.00794947264184679;
R2 = 6.21618016293642e-05;
Rc = 0.258351921779644;
Rd = 0.103662535963592;

V0 = 24;
PHI0 = 374.5; % [Wh]

% Datos curvas V-Phi modelo estatico Excel
% Ed0 = 23.53;
% Ed10 = -0.0129;
% Ec0 = 23.963;
% Ec10 = -0.0114;

% Datos curvas V-Phi modelo estatico ajuste bueno
Ed0 = 23.8656693114163;
Ed10 = -0.014039298384087;
Ec0 = 23.7100000000000;
Ec10 = 0.0119000000000000;

% Valores TFG Porras
% Ed0 = 24.45;
% Ed10 = 3.25e-6; % Cambiado signo por ser descarga
% Ec0 = 24.35;
% Ec10 = -2.982e-6;


%% Datos de los conversores DC-DC.

E15max = 0.91;
E5max  = 0.84;
E3max  = 0.78;



%% Calculo de los conversores de 5 V y 3.3 V.

% Potencia consumida.
W5 = Cons(:, 8) + Cons(:, 18);
W3 = Cons(:, 9);

% Intensidad de salida.
I5o = W5 ./  5 ;
I3o = W3 ./ 3.3;

% Eficiencia.
%E5 = E5max .* (1 - exp(-I5o ./ 2));
E5 = eff5( I5o );
%E3 = E3max .* (1 - exp(-I3o ./ 2));
E3 = eff3( I3o );

% Intensidad de entrada.
I5i = ( 5  / 15) .* I5o ./ E5;
I3i = (3.3 / 15) .* I3o ./ E3;

I5i( isnan(I5i) ) = 0;
I3i( isnan(I3i) ) = 0;



%% Calculo de los conversores de +15 V y -15 V.
 
% Potencia consumida.
W15n = Cons(:, 7);
W15p = Cons(:, 6);

% Intensidad de salida.
I15po = W15p ./ 15;
I15no = W15n ./ 15;

% Eficiencia en los railes +15 V y -15 V.
%E15p = E15max .* (1 - exp(-I15po ./ 2));
E15p = eff15( I15po );
%E15n = E15max .* (1 - exp(-I15no ./ 2));
E15n = eff15( I15no );

% La intensidad de entrada no se puede calcular porque V1 = VBUS varia

%% Test_potencia cte
P_test = ones(length(T),1)*5;



%% Funciones de eficiencia de los conversores

function E = eff5(x)
    E = 0.7707 + (0.1116 .* log(x));
end

function E = eff3(x)
    E = 0.7058 + (0.1128 .* log(x));
end


function E = eff15(x)
    E = 0.8549 + (0.1592 .* log(x));
end
