
from numpy import linspace, zeros, zeros_like, array, vstack, atleast_1d, tile, logical_and
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_bvp

N = 500
rho_CFRP = 1650
rho_HC = 40 
CP_CFRP = 930
CP_HC = 890
kx_CFRP = 130
ky_CFRP = 1
kz_CFRP = 130
kx_HC = 0.7
ky_HC = 1.5
kz_HC = 0.7
X = 500e-3
Y_CFRP = 2.5e-3
Y_HC = 15e-3
Z = 300e-3
T0 = 273 
q = 50

R_CFRP = X/(kx_CFRP*Z*Y_CFRP)
R_HC = X/(kx_HC*Z*Y_HC)
R_eq1 = 1/(2/R_CFRP +1/R_HC)
k_eq1 =  X/(R_eq1*(Z*Y_HC + 2*Z*Y_CFRP))


print(k_eq1)


def Fourier_1D(x, T, k=k_eq1, Q = 4*q/(X*Z)):
    dTdx = zeros_like(T)
    dTdx[0] = T[1]
    dTdx[1] = -Q/(k*(2*Y_CFRP+Y_HC))
    return dTdx
x = linspace(0, X, N)




y0 = zeros((2,x.size))

# Define the boundary conditions
def bc(Ta, Tb):
    return array([Ta[0]-T0, Tb[0]-T0])

# Define the initial guess for the solution
T_guess = zeros((2, x.size))
T_guess[0,:] = linspace(T0, T0+1, x.size)
# Solve the ODE system
sol = solve_bvp(Fourier_1D, bc, x, T_guess)

# Plot the solution
plt.plot(sol.x, sol.y[0,:]-273)
plt.xlabel('x')
plt.ylabel('Temperature')
plt.show()



# # Section 2

for i
def Fourier_1D_ap2(x, T, k=k_eq1, Q=4*q/(X*Z)):
    # u[0] = T(x), u[1] = T'(x)
    dTdx = zeros_like(T)
    dTdx[0] = T[1]
    if logical_and(x >= 0, x <= 50e-3):
        dTdx[1] = 0
    elif logical_and(x > 50e-3, x <= 150e-3):
        dTdx[1] = - Q / (k * (2 * Y_CFRP + Y_HC))
    elif logical_and(x > 150e-3, x <= 350e-3):
        dTdx[1] = 0
    elif logical_and(x > 350e-3, x <= 450e-3):
        dTdx[1] = - Q / (k * (2 * Y_CFRP + Y_HC))
    elif logical_and(x > 450e-3, x <= 500e-3):
        dTdx[1] = 0

    return dTdx


# Solve the ODE system
sol = solve_bvp(Fourier_1D_ap2, bc, x, T_guess)

# Plot the solution
plt.plot(sol.x, sol.y[0,:]-273)
plt.xlabel('x')
plt.ylabel('Temperature')
plt.show()
