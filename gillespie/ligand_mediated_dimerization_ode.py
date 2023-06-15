from scipy.integrate import odeint
import numpy as np

# Parameter values
# Rates:
r1 = 0.01
r2 = 0.1

# Initial conditions
# x1 and x2 are the initial displacements; y1 and y2 are the initial velocities
x1 = 100.0
x2 = 0.0

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 10.0
numpoints = 250

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
t_new = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
p = [r1, r2]
w0 = [x1, x2]


def vectorfield(w, t, p):
    """
    Defines the differential equations for the coupled spring-mass system.

    Arguments:
        w :  vector of the state variables:
                  w = [x1,y1,x2,y2]
        t :  time
        p :  vector of the parameters:
                  p = [m1,m2,k1,k2,L1,L2,b1,b2]

    a1 = r1 * x1  # 0  1  0
    a2 = r2 * x2  # 0  0  1
    a3 = r3 * x2  # 0 -1  0
    a4 = r4 * x3  # 0  0 -1
    """
    x1, x2 = w
    r1, r2 = p

    # Create f = (x1',x1',x3'):
    f = [r1 * x1 * x2 - r2 * x1 * (x1 - 1),
         - r1 * x1 * x2 + r2 * x1 * (x1 - 1)]
    return f

"""
# Call the ODE solver.
wsol = odeint(vectorfield, w0, t_new, args=(p,),
              atol=abserr, rtol=relerr)

def give_sol():
    return [wsol, t_new]


def f_x(t):
    return 
"""


def solve_CME():

    k_a = 0.01
    k_d = 0.1

    mol_tot = 100
    time_step = 0.1
    start_time = 0
    end_time = 10

    time_vals = np.arange(start_time,
                          end_time + time_step,
                          time_step)

    x_mat = np.zeroes(end_time / time_step + 1, mol_tot + 1)
    y_mat = x_mat.copy()

    initial_x = np.zeros(mol_tot + 1)
    initial_x[mol_tot] = 1
    x_mat[0, :] = initial_x

    t_num, x_num = np.shape(x_mat)

    for t in range(1, t_num):
        x = np.arange(0, mol_tot + 1)
        y = np.arange(0, mol_tot + 1)
        dxdt = k_a * (x[:-1] - 1) * (y[1:] + 1) * x_mat[t, :-1]
        x_mat[t, :] = x_mat[t - 1, :] + time_step * dxdt