C1 = 2173.48;
C2 = 281.904;
R1 = 0.016495;
R2 = 0.00014539;
Rc = 0.26;
Rd = 0.0171277;




P_cam = 5;
P_ADCS = 0.3;
P_OBC = 0.4;
P_TXR = 15;

mode = 1; %1=all on;  2=TXR off, cam on; 3=TXR on, cam off; 4=TXR off, cam off 

if mode == 1
    P_all = P_TXR + P_OBC + P_ADCS + P_cam;
elseif mode == 2
    P_all = P_OBC + P_ADCS + P_cam;
elseif mode == 3
    P_all = P_TXR + P_OBC + P_ADCS;
elseif mode == 4
    P_all = P_OBC + P_ADCS;
end

I_out = P_all./5;
eff_dcdc = 0.9072*(1 - exp(f_5(I_out)));
P_in = P_all ./eff_dcdc;






function f_5 = f_5(I)
    f_5 = -0.4022.*I.^3 + 1.6076*I.^2 - 3.4231.*I -0.0005;
end
