# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 11:07:08 2023

@author: Rafael Rivero
"""



class Material():
    def __init__(self, h, rho, c, kx, ky, kz):
        self.rho = rho
        self.c = c          # Capacidad térmica específica [ J/(kg·k)]
        self.h = h          # Espesor [m^3]
        self.kx = kx
        self.ky = ky
        self.ky = kz

# %% Inputs

Hp = 50                  # Heater Power [W]

CFRP = Material(rho=1650, h=2.5*1E-3, c=930, kx=130, ky=1, kz=130)
Honey = Material(rho=40, h=15*1E-3, c=890, kx=0.7, ky=1.5, kz=0.7)

Lx = 500 * 1E-3          # Longitud en el eje X [m]
Lz = 300 * 1E-3          # Longitud en el eje Z [m]
H = 2*CFRP.h + Honey.h   # Longitud en el eje Y [m]

To = 273                 # [K]

k_eq = ( 2 * CFRP.h * CFRP.kx + Honey.kx * Honey.h ) / ( H ) # Conducitividad equivalente unidimensional (x) [W/(m·K)]


N = 20 # Número de nodos para la resolución numérica por diferencias finitas

x =  linspace(0, Lx, 20)
x_n =  linspace(0, Lx, N); dx = x_n[1]-x_n[0]

# %% Apartado 1

q1 = 4 * Hp / (Lx*Lz)   # [W/m^2]
q4 = q1 / (k_eq*H)


# Analíticamente~~~~~~~~~~~~~~~~~~~~~~~~~~

T1_anal = -q4/2 * x**2 + q4 * Lx/2 * x + To # Distribución de temperatura analítica

# Numéricamente~~~~~~~~~~~~~~~~~~~~~~~~~~

# % % %  [------------------] [----]   [----]   _R :: Restringido, UR == Conocido, FR == Carga externa desconocida
# % % %  |   KRR  |  KRL    | | UR ]   [ FR ]   _L :: Libre, UL == Desconocido, FL == Carga externa conocida
# % % %  | -------|-------- | |----] = [----]   KRL = transpose(KLR)
# % % %  |   KLR  |  KLL    | | UL ]   [ FL ]   UL = (FL - KLR * UR) / KLL
# % % %  [------------------] [----]   [----]   FR = KRR * UR + KRL * UL

M =  zeros([N, N])

for node, dist in enumerate(x_n):
    if node == 0:
       M[0,0] = -2; M[0,1] = 1 # Invent a medias
    elif  node == N-1:
       M[N-1,N-1] = -2; M[N-1,N-2] = 1 # Invent a medias
    else:
        M[node, node-1] = 1; M[node, node] = -2; M[node, node+1] = 1

M_lib = M[1:-1, 1:-1]
M_res = array( [[M[0,0], M[0,-1]], [M[0,-1], M[-1,-1]]])
M_rl = array([ M[0,1:-1], M[-1, 1:-1] ])  ; M_lr = M_rl.transpose()
T_bound = array( [[To], [To]] )

Q1 = array( [ (-q4 * dx**2 ) for i in range(N-2)] )[ newaxis]; Q1 = Q1.transpose() # Carga térmica distribuida en todos los nodos, excepto en los extremos

T_num1 = dot( inv(M_lib) , ( Q1 - dot(M_lr,  T_bound) ) )
T1 =  concatenate( [ T_bound[0], T_num1[:,0],  T_bound[-1]] )


# %% Apartado 2

q_ap2 = 4 * Hp / (200*1E-3*Lz) # Se considera el área total sobre la que se introduce flujo de energía [W/m^2]
q_tot_2 =  q_ap2 / (k_eq*H)

# Analíticamente~~~~~~~~~~~~~~~~~~~~~~~~~~


T2_anal_a = To +  q_tot_2 * 0.1 * x[x<=0.05]

T2_anal_b = - q_tot_2 * (x[ logical_and(x>0.05, x<=0.15)] - 0.05)**2 / 2 +  q_tot_2 * 0.1 * (x[ logical_and(x>0.05, x<=0.15)] - 0.05) + (To +  q_tot_2 * 0.1 * 0.05)

T2_anal_c = (- q_tot_2 * (0.1)**2 / 2 +  q_tot_2 * 0.1 * (0.1) + (To +  q_tot_2 * 0.1 * 0.05) ) *  ones( len(x[ logical_and(x>0.15, x<=0.35)]) ) # Tramo constante

T2_anal_d =  flip(T2_anal_b)

T2_anal_e = To -  q_tot_2 * 0.1 * ( x[x>0.45] - 0.45 - 0.05)

T2_anal =  concatenate( [T2_anal_a, T2_anal_b, T2_anal_c, T2_anal_d, T2_anal_e] ) # Vector con la solución analítica completa


# Numéricamente~~~~~~~~~~~~~~~~~~~~~~~~~~
Q2 =  zeros([N-2, 1])

for i, n_value in enumerate(x_n[1:-1]):
    if n_value in x_n[ logical_and(x_n>=0.05, x_n<=0.15)]:
        Q2[i] = (- q_tot_2 * dx**2 )
    elif n_value in x_n[ logical_and(x_n>=0.35, x_n<=0.45)]:
        Q2[i] = (- q_tot_2 * dx**2 )



T_num2 = dot( inv(M_lib) , ( Q2 - dot(M_lr,  T_bound) ) )
T2 =  concatenate( [ T_bound[0], T_num2[:,0],  T_bound[-1]] )

#Tsensor = T2()
# %% Apartado 3
eps_CFRP = 0.8
sigma = 5.6708E-8 # [W/m2.K4]

T3_num=zeros([N-2, 1])
T3opt_num=zeros([N-2, 1])
Qrad=zeros([N-2, 1])
Qradopt=zeros([N-2, 1])


for i, n_value in enumerate(x_n[1:-1]):
    if n_value in x_n[ logical_and(x_n>=0.05, x_n<=0.15)]:
        Q2[i] = (- q_tot_2 * dx**2 )
    elif n_value in x_n[ logical_and(x_n>=0.35, x_n<=0.45)]:
        Q2[i] = (- q_tot_2 * dx**2 )

    Qradopt[i] = eps_CFRP*sigma*((T3_num[i]+273.15)**4 - (To)**4)*dx**2
    Qrad[i] = eps_CFRP*sigma*((T3_num[9]+273.15)**4 - (To)**4)*dx**2
    #Qrad5 = eps_CFRP*sigma*((T5(N/2)+273.15)**4 - (To)**4)*dx**2


    T3_num[i] = (T3_num[i+1] + T3_num[i+1])/2 + (Q2[i] - Qrad[i])
    T3opt_num[i] = (T3opt_num[i+1] + T3opt_num[i+1])/2 + (Q2[i] - Qradopt[i])

T3 =  concatenate( [ T_bound[0], T3_num[:,0],  T_bound[-1]] )
T3opt =  concatenate( [ T_bound[0], T3opt_num[:,0],  T_bound[-1]] )









# %% Representación gráfica

fig, ax = plt.subplots(1, 1, figsize=(7,4), constrained_layout=True) #, gridspec_kw={'width_ratios': [2, 1]})
#----------------------------------------------
# ax = axes[0]
ax.grid(); ax.set_title('Apartados 1 y 2')
ax.set_xlabel(r'$x$ [m]'); ax.set_ylabel(r'$T$ [$^\circ$C]')
ax.plot(x_n, T1 - 273, c='k', marker='^', ls='--', label=r'T(x)^{a}')
ax.plot(x_n, T2 - 273, c='b', marker='s', ls='--', label=r'T(x)^{a}')
ax.plot(x_n, T3 - 273, c='r', marker=',', ls='--', label=r'T(x)^{a}')
ax.plot(x, T1_anal - 273, c='black', label=r'T(x)^{a}')
ax.plot(x, T2_anal - 273, c='blue', label=r'T(x)^{a}')

# ax.legend(loc=0, fancybox=False ,edgecolor="black", ncol = 1, fontsize=11)
#--------------------------------------------
# ax = axes[1]
# ax.grid()
# ax.set_xlabel(r'Modes', fontsize=15)
# ax.set_ylabel(r'Energy of each mode', fontsize=15)
# ax.bar(np.arange(5)+1, energiaModos[0:5], width=0.5, bottom=None, align='center', color='black')



# ax.plot(x[1:], modosPropios[:,0], c='black', marker="^", label=r'1st mode')
# ax.plot(x[1:], modosPropios[:,1], c='blue', marker="v", label=r'2nd mode')
# ax.plot(x[1:], modosPropios[:,2], c='red', marker="s", label=r'3rd mode')
# ax.plot(x[1:], modosPropios[:,3], c='green', marker="o", label=r'4th mode')
# ax.plot(x[1:], modosPropios[:,4], c='gray', marker="d", label=r'5th mode')
