 clear all
close all

load('Perfiles_Consumo.mat')

%% Valores iniciales

Phi0 = 374.5; 
V0 = 24;

C1 = 2173.48;
C2 = 281.904;
R1 = 0.016495;
R2 = 0.00014539;
Rc = 0.26;
Rd = 0.0171277;

Ed0 = 24.129;
Ed10 = -0.012749;
Ed3 = -7.9149E-11;
Ed4 = 0.072135;

Ec0 = 20.142;
Ec10 = 0.012506;
Ec3 = -9.86E-14;
Ec4 = 0.0779;

%% Cosas sin cable

time = perfildeconsumo(:,1);
dt = time(2) - time(1);

Isp = perfildeconsumo(:,2);
P_m15 = perfildeconsumo(:,5);
P_33 = perfildeconsumo(:,7);

P_Trans_Bus = perfildeconsumo(:,8:11);
P_Trans_5V = perfildeconsumo(:,12);

P_BUS = perfildeconsumo(1,3) + sum(P_Trans_Bus,2);
P_5 = perfildeconsumo(1,6) + sum(P_Trans_5V,2);

I_33 = P_33./3.3;
I_5 = P_5./5;

eff_33 = 0.78*(1 - exp(f_33(I_33)));
eff_5 = 0.84*(1 - exp(f_5(I_5)));

P_33_in = P_33./eff_33;
P_5_in = P_5./eff_5;

P_15 = perfildeconsumo(:,4) + P_33_in + P_5_in;

I_15 = P_15/15;
I_m15 = P_m15/15;

eff_15 = 0.91*( 1 - exp(f_15(I_15)));
eff_m15 = 0.91*( 1 - exp(f_15(I_m15)));

P_15_in = P_15./eff_15;
P_m15_in = P_m15./eff_m15;

P_BUS_tot = P_BUS + P_15_in + P_m15_in;



%% Cosas con cable


L = 1; % [m]
rho = 	1.72e-8; % [Ohm.m] 
A = 0.324e-6; % [m2]
alpha = 3.8e-3; 
dT = 20;

R0_cab = rho*L/A;
R1_cab = R0_cab*(1+alpha*dT);


% Temperatura = 0° C
I_33_cab0 = (3.33 - sqrt(3.33^2 - 4*R0_cab.*P_33))/(2*R0_cab);
I_5_cab0 = (5 - sqrt(5^2 - 4*R0_cab.*P_5))/(2*R0_cab);

eff_33_cab0 = 0.78*(1 - exp(f_33(I_33_cab0)));
eff_5_cab0 = 0.84*(1 - exp(f_5(I_5_cab0)));

P_33_in_cab0 = (P_33 + R0_cab.*I_33_cab0.^2)./eff_33_cab0;
P_5_in_cab0 = (P_5 + R0_cab.*I_5_cab0.^2)./eff_5_cab0;

P_15_cab0 = perfildeconsumo(:,4) + P_33_in_cab0 + P_5_in_cab0;

I_15_cab0 = (15 - sqrt(15^2 - 4*R0_cab.*P_15))/(2*R0_cab);
I_m15_cab0 = (15 - sqrt(15^2 - 4*R0_cab.*P_m15))/(2*R0_cab);

eff_15_cab0 = 0.91*( 1 - exp(f_15(I_15_cab0)));
eff_m15_cab0 = 0.91*( 1 - exp(f_15(I_m15_cab0)));

P_15_in_cab0 = (P_15 + R0_cab.*I_15_cab0.^2)./eff_15_cab0;
P_m15_in_cab0 = (P_m15 + R0_cab.*I_m15_cab0.^2)./eff_m15_cab0;

P_BUS_tot_cab0 = P_BUS + P_15_in_cab0 + P_m15_in_cab0;


% Temperatura = 20° C
I_33_cab1 = (3.33 - sqrt(3.33^2 - 4*R1_cab.*P_33))/(2*R1_cab);
I_5_cab1 = (5 - sqrt(5^2 - 4*R1_cab.*P_5))/(2*R1_cab);

eff_33_cab1 = 0.78*(1 - exp(f_33(I_33_cab1)));
eff_5_cab1 = 0.84*(1 - exp(f_5(I_5_cab1)));

P_33_in_cab1 = (P_33 + R1_cab.*I_33_cab1.^2)./eff_33_cab1;
P_5_in_cab1 = (P_5 + R1_cab.*I_5_cab1.^2)./eff_5_cab1;

P_15_cab1 = perfildeconsumo(:,4) + P_33_in_cab1 + P_5_in_cab1;

I_15_cab1 = (15 - sqrt(15^2 - 4*R1_cab.*P_15))/(2*R1_cab);
I_m15_cab1 = (15 - sqrt(15^2 - 4*R1_cab.*P_m15))/(2*R1_cab);

eff_15_cab1 = 0.91*( 1 - exp(f_15(I_15_cab1)));
eff_m15_cab1 = 0.91*( 1 - exp(f_15(I_m15_cab1)));

P_15_in_cab1 = (P_15 + R1_cab.*I_15_cab1.^2)./eff_15_cab1;
P_m15_in_cab1 = (P_m15 + R1_cab.*I_m15_cab1.^2)./eff_m15_cab1;

P_BUS_tot_cab1 = P_BUS + P_15_in_cab1 + P_m15_in_cab1;

%% diff_0

diff_0 = P_BUS_tot_cab0./P_BUS_tot;
diff_1 = P_BUS_tot_cab1./P_BUS_tot;

figure(1)
plot(diff_0)
figure(2)
plot(diff_1)




%% Caso

Case = 4;

if Case ==1 
    P_BUS_sym = P_BUS_tot;
elseif Case == 2
    P_BUS_sym = P_BUS_tot_cab0;
elseif Case == 3
    P_BUS_sym = P_BUS_tot_cab1;
elseif Case == 4
    P_BUS_sym = zeros(size(P_BUS_tot_cab1));
end

%% Cosas
% 
% syms y
% 
% cosa_chula = y*5/(0.84*(1-exp(f_5(y))));
% 
% fplot(cosa_chula)


% V_33R0 = (3.3 + sqrt(3.3^2 - 4*R0*P_33))/2;
% I_33_R0 = P_33/V_33;
% 
% V_5R0 = (5 + sqrt(5^2 - 4*R0*P_5))/2;
% I_5_R0 = P_5/V_5;
% 
% V_15R0 = (15 + sqrt(15^2 - 4*R0*P_15))/2;
% I_15_R0 = P_15/V_15;
% 
% V_m15R0 = (-15 - sqrt(15^2 - 4*R0*P_m15))/2;
% I_m15_R0 = P_m15/V_m15;
% 
% V_33R1 = (3.3 + sqrt(3.3^2 - 4*R1*P_33))/2;
% I_33_R1 = P_33/V_33;
% 
% V_5R1 = (5 + sqrt(5^2 - 4*R1*P_5))/2;
% I_5_R1 = P_5/V_5;
% 
% V_15R1 = (15 + sqrt(15^2 - 4*R1*P_15))/2;
% I_15_R1 = P_15/V_15;
% 
% V_m15R1 = (-15 - sqrt(15^2 - 4*R1*P_m15))/2;
% I_m15_R1 = P_m15/V_m15;

%% Funciones

function f_33 = f_33(I)
    f_33 = -1.9824.*I.^3 + 4.2363.*I.^2 - 5.073.*I + 0.0212;
end

function f_5 = f_5(I)
    f_5 = -2.235.*I.^3 + 4.9524*I.^2 - 5.9881.*I -0.0321;
end

function f_15 = f_15(I)
    f_15 = 1.4962*I.^3 - 3.7408*I.^2 - 1.3329*I - 0.5245;
end