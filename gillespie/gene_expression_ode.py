from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameter values
# Rates:
r1 = 0.5  # 0.01
r2 = 0.5  # 0.1
r3 = 0.3  # 0.0001
r4 = 0.4  # 0.0001

# Initial conditions
# x1 and x2 are the initial displacements; y1 and y2 are the initial velocities
x1 = 100.0
x2 = 1.0
x3 = 0.0

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 100.0
numpoints = 250

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
p = [r1, r2, r3, r4]
w0 = [x1, x2, x3]


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
    x1, x2, x3 = w
    r1, r2, r3, r4 = p

    # Create f = (x1',x1',x3'):
    f = [0,
         r1 * x1 - r3 * x2,
         r2 * x2 - r4 * x3]
    return f


# Call the ODE solver.
wsol = odeint(vectorfield, w0, t, args=(p,),
              atol=abserr, rtol=relerr)

fig, ax = plt.subplots()
ax.plot(t, wsol, label=['DNA',
                        'mRNA',
                        'Protein'])
ax.legend(loc='best', shadow=True)
plt.show()
