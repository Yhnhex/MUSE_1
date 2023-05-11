clear all
close all

%% Configuracion de las graficas
linewidth = 3;
fontsize = 20;

%% Load Data
% load('Data_Est');
load('Data_Din');

%% Variables
Data0 = 1;
Dataf = length(Bat_din);

% Data
Time = Bat_din(Data0:Dataf,1);
V = Bat_din(Data0:Dataf,2);
I = -Bat_din(Data0:Dataf,3);

Rd = 0.1036;
Rc = 0.1699;
Ed = V + Rd*I;
Ec = V - Rc*I;
Ec_Ed = Ec-Ed;


% Variables a calcular
% u0 = [9e-2, 2e-2, 1000, 200];
u0 = [0.05, 0.05, 1000, 200];

R1 = u0(1);
R2 = u0(2);
C1 = u0(3);
C2 = u0(4);








