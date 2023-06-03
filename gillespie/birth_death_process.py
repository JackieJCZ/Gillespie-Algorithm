import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Birth - Death Process:
BD_stoch_subst = np.array([[1, 0],
                           [0, 1]])
BD_stoch_prods = np.array([[0, 1],
                           [1, 0]])
BD_rates = np.array([0.5, 0.5])
BD_X0 = np.array([100, 0])
BD_t_max = 5


def stationary_point(X, num=20):
    return np.average(X[-num:],
                      axis=0)


elapsed_time, X, t = gillespie(BD_X0,
                               BD_rates,
                               BD_stoch_subst,
                               BD_stoch_prods,
                               BD_t_max)
selected = [1]


# # Standard plot of distributions over time


def call_base_gillespie():
    print(
        f"Elapsed time: {elapsed_time}ms"
    )

    fig, ax = plt.subplots()
    if selected is None:
        for i in range(len(BD_X0)):
            ax.plot(t, X[:, i], '-*', label='X' + str(i))
    else:
        for i in selected:
            ax.plot(t, X[:, i], '-*', label='X' + str(i))
    ax.legend(loc='best', shadow=True)
    plt.show()


call_base_gillespie()


# # Stationary distribution of the base algorithm


def call_stationary():
    end_vals = []
    for i in range(10000):
        elapsed_time, X, t = gillespie(BD_X0,
                                       BD_rates,
                                       BD_stoch_subst,
                                       BD_stoch_prods,
                                       BD_t_max)
        end_vals.append(X[-1][1])
    plt.hist(end_vals,
             range=(0, np.max(end_vals)),
             bins=int(np.max(end_vals)))


call_stationary()


# # Imshow plot of stationary points for different k values


def find_stationary(ks, kd):
    e, X, t = gillespie(BD_X0,
                        [ks, kd],
                        BD_stoch_subst,
                        BD_stoch_prods,
                        BD_t_max)
    return stationary_point(X)[1]


def call_imshow():
    ks, kd = np.mgrid[0:1.1:0.01, 0:1.1:0.01]

    D = np.vectorize(find_stationary)(ks, kd)

    plt.imshow(D,
               origin='lower',
               extent=(0, 1, 0, 1))
    plt.colorbar()
    plt.show()


# call_imshow()
