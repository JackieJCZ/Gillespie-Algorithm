import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie
from scipy.interpolate import interp1d
from scipy.stats import entropy

plt.rcParams.update({'font.size': 20})

# Ligand-Mediated Dimerization
LMD_stoch_subst = np.array([[1, 1],
                            [2, 0]])
LMD_stoch_prods = np.array([[2, 0],
                            [1, 1]])
LMD_rates = np.array([0.01, 0.1])
LMD_X0 = np.array([100, 0])
LMD_t_max = 10

a = LMD_rates[1] + LMD_rates[0] * (LMD_X0[0] + LMD_X0[1])
b = LMD_rates[0] + LMD_rates[1]
c = (LMD_rates[0] * LMD_X0[1] / LMD_X0[0]
     + LMD_rates[1] * (1 / LMD_X0[0] - 1))

stoptime = 10.0
numpoints = 250
t_new = np.array([stoptime * float(i) / (numpoints - 1)
                 for i in range(numpoints - 10)])
iters = 10


def interpolate(t_vals, x_vals):
    inter_func = interp1d(t_vals, x_vals)
    return inter_func(t_new)


def f_x(t):
    return a / (b + c * np.exp(-a * t))


def f_y(t):
    return LMD_X0[0] + LMD_X0[1] - f_x(t)


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
    x = f_x(t)
    y = f_y(t)

    # Create f = (x1',x1',x3'):
    f = [r1 * (x - 1) * (y + 1) * x1
         - (r1 * y + r2 * (x - 1)) * x * x1
         + r2 * x]
    return f


def call_main():
    fig, ax = plt.subplots(1, 2,
                           sharex=True,
                           sharey=True,
                           figsize=(12, 5))

    all_sols_x0 = np.zeros((iters, len(t_new)))
    all_sols_x1 = np.zeros((iters, len(t_new)))

    ax[0].set_xlabel('time / s')
    ax[0].set_ylabel('Molecule number')
    ax[1].set_xlabel('time / s')
    fig.suptitle('Trajectories over time')

    for i in range(iters):
        elapsed_time, X, t = gillespie(LMD_X0,
                                       LMD_rates,
                                       LMD_stoch_subst,
                                       LMD_stoch_prods,
                                       LMD_t_max)

        for j in range(len(LMD_X0)):
            ax[0].step(t, X[:, j],
                       linewidth=0.5)

            ax[1].step(t, X[:, j],
                       linewidth=0.5)

            all_sols_x0[i, :] = interpolate(t,
                                            X[:, 0])
            all_sols_x1[i, :] = interpolate(t,
                                            X[:, 1])

    all_sols_x0 = np.mean(all_sols_x0, axis=0)
    all_sols_x1 = np.mean(all_sols_x1, axis=0)

    ax[0].plot(t_new, all_sols_x0,
               '--',
               c='black',
               alpha=0.75,
               linewidth=1.5,
               label='Dimer mean')
    ax[0].plot(t_new, all_sols_x1,
               '--',
               c='blue',
               alpha=0.75,
               linewidth=1.5,
               label='Ligand mean')

    ax[1].plot(t_new, all_sols_x0,
               c='black',
               alpha=0.75,
               linewidth=1.5,
               label='Dimer mean')
    ax[1].plot(t_new, all_sols_x1,
               c='blue',
               alpha=0.75,
               linewidth=1.5,
               label='Ligand mean')

    y = f_x(t_new)

    ax[0].plot(t_new, y, label='Dimer',
               c='black',
               linewidth=1.5)
    ax[0].plot(t_new, (LMD_X0[0] + LMD_X0[1] - y),
               label='Ligand',
               c='blue',
               linewidth=1.5)

    ax[1].axhline(a / b,
                  linestyle='--',
                  c='black',
                  alpha=0.5,
                  linewidth=1.5)
    ax[1].axhline((LMD_X0[0] + LMD_X0[1] - a / b),
                  linestyle='--',
                  c='blue',
                  alpha=0.5,
                  linewidth=1.5)

    ax[0].locator_params(axis='x', nbins=6)
    ax[0].legend(loc='best', shadow=False)
    ax[1].legend(loc='best', shadow=False)
    fig.subplots_adjust(wspace=0)
    plt.show()

call_main()
