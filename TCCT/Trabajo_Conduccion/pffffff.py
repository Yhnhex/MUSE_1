import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Define the parameters and boundary conditions
L = 1.0    # length of rod
k = 1.0    # thermal conductivity
q = 1.0    # heat source/sink
T0 = 0.0   # left boundary condition
T1 = 1.0   # right boundary condition

# Define the ODE system
def heat_eq(x, T):
    dTdx = np.zeros_like(T)
    dTdx[0] = T[1]
    dTdx[1] = -q/k
    return dTdx

# Define the boundary conditions
def bc(Ta, Tb):
    return np.array([Ta[0]-T0, Tb[0]-T1])

# Define the initial guess for the solution
x = np.linspace(0, L, 100)
T_guess = np.zeros((2, x.size))
T_guess[0,:] = np.linspace(T0, T1, x.size)

# Solve the ODE system
sol = solve_bvp(heat_eq, bc, x, T_guess)

# Plot the solution
plt.plot(sol.x, sol.y[0,:])
plt.xlabel('x')
plt.ylabel('Temperature')
plt.show()
