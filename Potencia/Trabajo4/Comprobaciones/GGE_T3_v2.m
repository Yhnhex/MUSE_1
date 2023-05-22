clear
close all

%% Configuracion de las graficas
linewidth = 3;
fontsize = 20;

%% Load Data
% load('Data_Est');
% load('Data_Din');

M=xlsread("datos_bateria_6S6P.xlsx", 9, 'B3:D29584');
V = M(3:end, 1);
I = -M(3:end, 2);
t = M(3:end, 3);

R1 = 0.0506063212140136;
R2 = 0.0101258648525245;
R3 = 0.001;

C1 = 200.213974865193;
C2 = 1008.78035302138;
C3 = 100;

Rc = 0.180690162482483;
Rd = 0.0607407915012055;
Ec = V - Rc*I;
Ed = V + Rd*I;
Ec_Ed = Ec - Ed;

% %% Variables
% Data0 = 1;
% Dataf = length(Bat_din);

% % Data
% Time = Bat_din(Data0:Dataf,1);
% V = Bat_din(Data0:Dataf,2);
% I = -Bat_din(Data0:Dataf,3);

% Rd = 0.17552;
% Rc = 0.18028;
% Ed = V + Rd*I;
% Ec = V - Rc*I;
% Ec_Ed = Ec-Ed;


% Variables a calcular
% u0 = [9e-2, 2e-2, 1000, 200];
u0 = [0.05, 0.05, 1000, 200];

R1 = u0(1);
R2 = u0(2);
C1 = u0(3);
C2 = u0(4);








