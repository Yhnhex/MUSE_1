from numpy.linalg import inv
from numpy import amax, amin, array, dot, newaxis, logical_and, concatenate, zeros, linspace, ones, flip
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')


matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}

matplotlib.rc('font', **font)
# Input data
rho_CFRP = 1650    # [kg/m^3]
rho_HC = 40        # [kg/m^3] 
C_CFRP = 930       # [J/kg k]
C_HC = 890         # [J/kg k]
kx_CFRP = 130      # [W/m K]
ky_CFRP = 1        # [W/m K]
kz_CFRP = 130      # [W/m K]
kx_HC = 0.7        # [W/m K]
ky_HC = 1.5        # [W/m K] 
kz_HC = 0.7        # [W/m K]
X = 500e-3         # [m]
Y_CFRP = 2.5e-3    # [m]
Y_HC = 15e-3       # [m]
Z = 300e-3         # [m]
T0 = 273           # [K]
eps_CFRP = 0.8
Hp = 50            # [W]
sigma = 5.6708E-8  # [W/m^2.K^4]
t_heat=0

R_CFRP = X/(kx_CFRP*Z*Y_CFRP)
R_HC = X/(kx_HC*Z*Y_HC)
R_eq1 = 1/(2/R_CFRP +1/R_HC)
k_eq =  X/(R_eq1*(Z*Y_HC + 2*Z*Y_CFRP))


Lx = 500 * 1E-3          # Longitud en el eje X [m]
Lz = 300 * 1E-3          # Longitud en el eje Z [m]
H = 2*Y_CFRP + Y_HC      # Longitud en el eje Y [m]

To = 273                 # [K]
rho_eq = (2*rho_CFRP*Y_CFRP + rho_HC*Y_HC)/H
C_eq = (2*C_CFRP*Y_CFRP + C_HC*Y_HC)/H



x =  linspace(0, Lx, 150)

# %% Apartado 1 y 2 analíticos
print(k_eq)

q1 = 2 * Hp / (Lx*Lz)   # [W/m^2]
q4 = q1 / (k_eq*H)
T1_anal = -q4/2 * x**2 + q4 * Lx/2 * x + To # Distribución de temperatura analítica


q_ap2 = 2 * Hp / (200*1E-3*Lz) # Se considera el área total sobre la que se introduce flujo de energía [W/m^2]
q_tot_2 =  q_ap2 / (k_eq*H)

# Analíticamente~~~~~~~~~~~~~~~~~~~~~~~~~~


T2_anal_a = To +  q_tot_2 * 0.1 * x[x<=0.05]

T2_anal_b = - q_tot_2 * (x[ logical_and(x>0.05, x<=0.15)] - 0.05)**2 / 2 +  q_tot_2 * 0.1 * (x[ logical_and(x>0.05, x<=0.15)] - 0.05) + (To +  q_tot_2 * 0.1 * 0.05)

T2_anal_c = (- q_tot_2 * (0.1)**2 / 2 +  q_tot_2 * 0.1 * (0.1) + (To +  q_tot_2 * 0.1 * 0.05) ) *  ones( len(x[ logical_and(x>0.15, x<=0.35)]) ) # Tramo constante

T2_anal_d =  flip(T2_anal_b)

T2_anal_e = To -  q_tot_2 * 0.1 * ( x[x>0.45] - 0.45 - 0.05)

T2_anal =  concatenate( [T2_anal_a, T2_anal_b, T2_anal_c, T2_anal_d, T2_anal_e] ) # Vector con la solución analítica completa

# Numéricamente~~~~~~~~~~~~~~~~~~~~~~~~~~
# %%

nx = 150
dx = Lx/nx
t=1000
nt = 10000
dt = t/nt
xplot = linspace(0, Lx, nx)

T1 =  zeros(nx)
T2 =  zeros(nx)
T3 =  zeros(nx)
T3_nolin =  zeros(nx)
T5 = zeros(nx)
T4 =  zeros([nx, nt])
T5_mat =  zeros([nx, nt])
T42 =  zeros([nx, nt])
Qrad4 = zeros(nt)
Qrad5 = zeros(nt)
Qheat5 = zeros(nt)

diff = zeros(nt)
Qrad4[0] = eps_CFRP * sigma * ((T4[0, 0] + 273.15)**4 - (To)**4)
Qrad5[0] = eps_CFRP * sigma * ((T5_mat[0, 0] + 273.15)**4 - (To)**4)


t_heat=0

for j in range(nt - 1):
    for i in range(1, nx - 1):
        if ((i >= nx/10) and (i <= 3*nx/10) or (i >= 7*nx/10) and (i <= 9*nx/10)):
            Qheat = Hp / (Lz*H * (nx/5))
        else:
            Qheat = 0

        

        Qrad_nolin = eps_CFRP * sigma * ((T3[i] + 273.15)**4 - (To)**4)
        Qrad = 4 * eps_CFRP*sigma * (To)**3 * (T3[nx//2] + 273.15 - To)
        Qrad5[j] = eps_CFRP * sigma * ((T5[i] + 273.15)**4 - (To)**4)
        Qrad4[j] = eps_CFRP * sigma * ((T4[i, j+1] + 273.15)**4 - (To)**4)

        T1[i] = ((T1[i+1] + T1[i-1]) / 2 + dx**2 * Hp / (k_eq * Lz*H * (nx-1))) 
        T2[i] = ((T2[i+1] + T2[i-1]) / 2 + dx**2 * Qheat / k_eq) 

        T3[i] = ((T3[i+1] + T3[i-1]) / 2 + dx**2 * (Qheat/k_eq - Qrad/k_eq)) 
        T3_nolin[i] = ((T3_nolin[i+1] + T3_nolin[i-1]) / 2 + dx**2 * (Qheat/k_eq - Qrad_nolin/k_eq)) 
       
        if T5[nx//2] < 10/nx:
            heater_on = True
        elif T5[nx//2] > 12/nx:
            heater_on = False
    
        if heater_on:
            T5[i] = (T5[i+1] + T5[i-1])/2 + dx**2*(Qheat - Qrad5[j])/k_eq
            t_heat = t_heat + dt
        else:
            T5[i] = (T5[i+1] + T5[i-1])/2 - dx**2*(Qrad)/k_eq
        

    diff[j] = T4[nx//2, j]-T3[nx//2]
    T5_mat[:, j]=T5
    T42[:, j]=T3


        

T1 = T1 * 2 * nx
T2 = T2 * nx
T3 = T3 * nx
T3_nolin = T3_nolin * nx
T42 = T42 * nx

T5_mat = T5_mat * nx
T3diff = (T3_nolin - T3) / T3_nolin

# %%
from numpy import amin
t_heat_av = t_heat/t/100
P_med = t_heat_av*2*Hp

T1_max = amax(T1)
T1_s = T1[nx//2]
T1_max_anal = amax(T1_anal-273)
T1_s_anal = T1_anal[nx//2]-273
T2_max = amax(T2)
T2_s = T2[nx//2]
T2_max_anal = amax(T2_anal-273)
T2_s_anal = T2_anal[nx//2]-273
T3_max = amax(T3)
T3_s = T3[nx//2]
T3_max_nl = amax(T3_nolin)
T3_s_nl = T3_nolin[nx//2]
T42_max = amax(T42)
T42_s = T42[nx//2, -2]
T5_max = amax(T5_mat[nx//2, :])










# %% Representación gráfica
# APARTADO 1
fig, ax = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax.grid(); ax.set_title('')
ax.set_xlabel(r'$x$ [m]'); ax.set_ylabel(r'$T$ [$^\circ$C]')
ax.plot(xplot, T1 , c='r', ls='--', label=r'$T(x)$ diferencias finitas')
ax.plot(x, T1_anal - 273, c='black', label=r'$T(x)$ analitica')
ax.legend(loc = "best")
ax.set_xlim([0, 0.5]); ax.set_ylim([0, 35])
plt.show(); t_heat_av = 0.4892; P_med=t_heat_av*100

# APARTADO 2
fig2, ax2 = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax2.grid(); ax2.set_title('')
ax2.set_xlabel(r'$x$ [m]'); ax2.set_ylabel(r'$T$ [$^\circ$C]')
ax2.plot(xplot, T2 , c='r', ls='--', label=r'$T(x)$ diferencias finitas')
ax2.plot(x, T2_anal - 273, c='black', label=r'$T(x)$ analitica')
ax2.legend(loc = "best")
ax2.set_xlim([0, 0.5]); ax2.set_ylim([0, 30])


# APARTADO 3
fig3, ax3 = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax3.grid(); ax3.set_title('')
ax3.set_xlabel(r'$x$ [m]'); ax3.set_ylabel(r'$T$ [$^\circ$C]')
ax3.plot(xplot, T3_nolin, c='black', label=r'$T(x)$ no lineal')
ax3.plot(xplot, T3 , c='r', ls='--', label=r'$T(x)$ linealizada')
ax3.legend(loc = "best")
ax3.set_xlim([0, 0.5]); ax3.set_ylim([0, 30])

fig32, ax32 = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax32.grid(); ax32.set_title('')
ax32.set_xlabel(r'$x$ [m]'); ax32.set_ylabel(r'$\frac{T_3-T_{3,nl}}{T_3}$ ')
ax32.plot(xplot, T3diff , c='k', label=r'T(x)^{a}')

# ax.plot(xplot, T2 , c='b', marker='s', ls='--', label=r'T(x)^{a}')
# ax.plot(xplot, T3 , c='r', marker=',', ls='--', label=r'T(x)^{a}')
# ax.plot(xplot, T4[:,-1] , c='r', marker=',', ls='--', label=r'T(x)^{a}')
# # ax.plot(xplot, T5[:,-1] , c='g', marker='o', ls='--', label=r'T(x)^{a}')

# ax.plot(x, T2_anal - 273, c='blue', label=r'T(x)^{a}')

# APARTADO 4
fig4, ax4 = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax4.grid(); ax4.set_title('')
ax4.set_xlabel(r'$t$ [s]'); ax4.set_ylabel(r'$T$ [$^\circ$C]')
ax4.plot(linspace(0, t, nt-1), T42[nx//2, 0:nt-1], c='black', label=r'$T(t)$ transitorio')
ax4.axhline(T3[nx//2], c='r', ls='--', label=r'$T(t)$ estacionario')
ax4.legend(loc = "best")
ax.set_xlim([0, 1000]); ax.set_ylim([0, 30])


fig42, ax42 = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax42.grid(); ax42.set_title('')
ax42.set_xlabel(r'$t$ [s]'); ax42.set_ylabel(r'$T$ [$^\circ$C]')
ax42.plot(xplot, T42[:, nt-3], c='black', label=r'$T(x, t=1000s)$ Ap.4')
ax42.plot(xplot, T3, c='r', ls='--', label=r'$T(x)$ Ap.3')
ax42.legend(loc = "best")
ax42.set_xlim([0, 0.5]); ax42.set_ylim([0, 30])
plt.show()


fig43, ax43 = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax43.grid(); ax43.set_title('')
ax43.set_xlabel(r'$t$ [s]'); ax43.set_ylabel(r'$\frac{T_4-T_{3}}{T_3}$')
ax43.plot(linspace(0, t, nt-1), abs((T42[nx//2, 0:nt-1]-T3[nx//2])/T3[nx//2]), c='black', label=r'$T(t)$')
plt.show()


fig5, ax5 = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
ax5.grid(); ax5.set_title('')
ax5.set_xlabel(r'$t$ [s]'); ax5.set_ylabel(r'$T$ [$^\circ$C]')
ax5.axhline(12, c='r', ls='--', label=r'$T_{off}$ ')
ax5.axhline(10, c='b', ls='--', label=r'$T_{on}$ ')
ax5.plot(linspace(0, t, nt-1), T5_mat[nx//2, 0:nt-1], c='black', label=r'$T(t)$')
ax5.legend(loc = "best")
ax5.set_xlim([0, 1000]); ax.set_ylim([0, 20])

# %%
