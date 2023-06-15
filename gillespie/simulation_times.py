import numpy as np
from gillespie import gillespie
from matplotlib import pyplot as plt

# Degradation of Biomolecule:
DB_stoch_subst = np.array([[1]])
DB_stoch_prods = np.array([[0]])
DB_rates = np.array([0.5])
DB_X0 = np.array([100])
DB_t_max = 11

# Base gillespie

time_vals = [1,10,100,1000,10000]
x_vals = [1,10,100,1000,10000]

def new_gill(X0, t):
    e_v = []
    for i in range(10):
        e, X, tp = gillespie(np.array([X0]), DB_rates,
                            DB_stoch_subst,
                            DB_stoch_prods,
                            t)
        e_v.append(e)
    return np.mean(e_v)

x_val, t_val = np.meshgrid(x_vals, time_vals)
x_val_ans = []
for x in x_vals:
    # t_vals = np.power(10, np.linspace(0, 4, 10))
    # x_plot = np.vectorize(new_gill)(x, t_vals)
    # plt.plot(t_vals, x_plot, label=f'X0={x}')
    # plt.legend()

    x_val_ans.append(new_gill(x, 10))

plt.plot(x_vals, x_val_ans)