import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Gene Expression:

# TOGGLE SWITCH:

GE_stoch_subst = np.array([[1, 0, 0],
                           [0, 1, 0],
                           [0, 1, 0],
                           [0, 0, 1]])
GE_stoch_prods = np.array([[1, 1, 0],
                           [0, 1, 1],
                           [0, 0, 0],
                           [0, 0, 0]])
# GE_rates = np.array([0.01, 0.1, 0.0001, 0.0001])
GE_rates = np.array([0.5, 0.5, 0.3, 0.4])
GE_X0 = np.array([100, 1, 0])
GE_t_max = 100

elapsed_time, X, t = gillespie(GE_X0, GE_rates,
                               GE_stoch_subst,
                               GE_stoch_prods,
                               GE_t_max)

selected = None


def call_base_gillespie():

    print(
        f"Elapsed time: {elapsed_time}ms"
    )

    fig, ax = plt.subplots()
    if selected is None:
        for i in range(len(GE_X0)):
            ax.step(t, X[:, i], label='X' + str(i))
    else:
        for i in selected:
            ax.step(t, X[:, i], label='X' + str(i))
    ax.legend(loc='best', shadow=True)
    plt.show()


# call_base_gillespie()


def run_multiple(iters):
    vals = []
    for i in range(iters):
        e, X, t = gillespie(GE_X0, GE_rates,
                            GE_stoch_subst,
                            GE_stoch_prods,
                            GE_t_max)
        vals.append(X[-1, :])
    print(np.mean(vals, axis=0))

beta = GE_rates[0] * GE_X0[0]
gamma = GE_rates[2]
# run_multiple(100)

from scipy.special import factorial
import numpy as np

def Pm(n):
    vals = [np.exp(-beta / gamma)]
    for i in range(1, n + 1):
        vals.append(vals[-1] * (beta / gamma) / i)
    return vals

n = 400
x_vals = [i for i in range(n + 1)]

plt.plot(x_vals, Pm(n))