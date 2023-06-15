import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial
from gillespie import gillespie

# Birth - Death Process:
BD_stoch_subst = np.array([[0],
                           [1]])
BD_stoch_prods = np.array([[1],
                           [0]])
BD_rates = np.array([0.2, 0.001])
BD_X0 = np.array([100])
BD_t_max = 1000

λ = BD_rates[0] / BD_rates[1]


def poisson(x):
    return np.power(λ, x) * np.exp(-λ) / factorial(x)


def stationary_point(X, num=20):
    return np.average(X[-num:],
                      axis=0)


# # Standard plot of distributions over time

selected = None


def call_base_gillespie():

    elapsed_time, X, t = gillespie(BD_X0,
                                   BD_rates,
                                   BD_stoch_subst,
                                   BD_stoch_prods,
                                   BD_t_max)

    print(
        f"Elapsed time: {elapsed_time}ms"
    )

    fig, ax = plt.subplots()
    if selected is None:
        for i in range(len(BD_X0)):
            ax.step(t, X[:, i], '-*', label='X' + str(i))
    else:
        for i in selected:
            ax.step(t, X[:, i], '-*', label='X' + str(i))
    ax.legend(loc='best', shadow=True)
    plt.show()


call_base_gillespie()


# # Stationary distribution of the base algorithm


def call_stationary():
    iters = 10000
    end_vals = []
    for i in range(iters):
        elapsed_time, X, t = gillespie(BD_X0,
                                       BD_rates,
                                       BD_stoch_subst,
                                       BD_stoch_prods,
                                       BD_t_max)
        end_vals.append(X[-1][0])
    max_val = int(np.max(end_vals))
    p_points = np.linspace(0, max_val, 100)
    fig, ax = plt.subplots()
    ax.hist(end_vals,
            range=(0, max_val),
            bins=max_val + 1)
    ax.plot(p_points, poisson(p_points) * iters)


# call_stationary()


# # Imshow plot of stationary points for different k values


def find_stationary(ks, kd):
    e, X, t = gillespie(BD_X0,
                        [ks, kd],
                        BD_stoch_subst,
                        BD_stoch_prods,
                        BD_t_max)
    return stationary_point(X)[0]


def call_imshow():
    step = 0.01
    ks, kd = np.mgrid[0.1:1 + step:step, 0.1:1 + step:step]

    D = np.vectorize(find_stationary)(ks, kd)
    _max = np.max(D)

    fig, axes = plt.subplots(2, 1)
    ax1, ax2 = axes
    im = ax1.imshow(D,
                    origin='lower',
                    extent=(0, 1, 0, 1),
                    vmin=0,
                    vmax=_max)

    D_theory = ks / kd

    ax2.imshow(D_theory,
               origin='lower',
               extent=(0, 1, 0, 1),
               vmin=0,
               vmax=_max)

    fig.colorbar(im, ax=axes.ravel().tolist())


# call_imshow()
