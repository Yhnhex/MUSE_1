clear all
M=xlsread("datos_bateria_6S6P.xlsx", 9, 'B3:D29584');
V = M(3:end, 1);
I = -M(3:end, 2);
t = M(3:end, 3);


R_1 = 0.01;
R_2 = 0.01;
C_1 = 1000;
C_2 = 1000;


R_c = 0.18028;
R_d = 0.17552;
E_c = V - R_c*I;
E_d = V + R_d*I;
Ec_Ed = E_c - E_d;


plot(t, V)