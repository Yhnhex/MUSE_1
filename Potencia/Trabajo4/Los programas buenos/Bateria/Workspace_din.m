M=xlsread("medidas_bateria_dinamico.xlsx");
t = M(3:end, 1);
V = M(3:end, 2);
I = -M(3:end, 3);


R_1 = 0.01;
R_2 = 0.01;
C_1 = 1000;
C_2 = 1000;


R_c = 0.14415;
R_d = 0.11872;
E_c = V - R_c*I;
E_d = V + R_d*I;
Ec_Ed = E_c - E_d;