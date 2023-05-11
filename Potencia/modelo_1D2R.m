clear all
V_oc =2.7;
I_sc = 0.5202;
V_mp = 2.411;
I_mp = 0.5044;

V_T = 0.026;
% CALCULATE APPROX R_S
   a_prox=1;
   fun_R_s = @(R_ss) (a_prox*V_T*V_mp*(2*I_mp - I_sc))/((V_mp*I_sc+V_oc*(I_mp-I_sc))*(V_mp-I_mp*R_ss)-a_prox*V_T*(V_mp*I_sc-V_oc*I_mp))-exp((V_mp+I_mp*R_ss-V_oc)/a_prox*V_T);
   R_s = fsolve(fun_R_s, 1.3535);

 % CALCULATE R_SHUNT

   R_sh= ((V_mp-I_mp*R_s)*(V_mp-R_s*(I_sc-I_mp)-a_prox*V_T))/((V_mp-I_mp*R_s)*(I_sc-I_mp)-a_prox*V_T*I_mp);
 % CALCULATE REAL a esto realmente ns que paque si no hase na
   R_sh0= R_s+R_sh;
   a =((V_mp-I_mp*R_s)*(V_mp+(I_mp-I_sc)*R_sh0))/((V_mp-I_mp*R_sh0)*V_T);

 % CALCULATE REAL R_S
   A = (V_mp+(I_mp-I_sc)*R_sh0)*log((V_mp+(I_mp-I_sc)*R_sh0)/(V_oc-I_sc*R_sh0));
   B = V_mp - R_sh0*I_mp;
   R_s = (A-B)/(A+B)*V_mp/I_mp + B/(A+B)*V_oc/I_mp;
 
 % CALCULATE I_pv
   I_pv = (R_sh+R_s)/R_sh*I_sc;
   I_o = ((R_sh+R_s)*I_sc-V_oc)/(R_sh*exp(V_oc/(a*V_T)));
 
  
  V     = 0:0.1:V_oc;
 % define function you want to solve for: fun = 0
 fun   = @(I) I_pv - I_o * exp((V + I*R_s)/(a*V_T)) - (V + I*R_s)/R_sh - I;
 % get solution for all V
 I_sol = fsolve(fun, I_pv*ones(size(V)));
plot(V,I_sol)
