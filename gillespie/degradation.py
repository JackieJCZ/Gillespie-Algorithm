import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie
from scipy.special import binom
from scipy.stats import entropy
from scipy.interpolate import interp1d

# Degradation of Biomolecule:
DB_stoch_subst = np.array([[1]])
DB_stoch_prods = np.array([[0]])
DB_rates = np.array([0.5])
DB_X0 = np.array([100])
DB_t_max = 11

# Base gillespie


def call_base_gillespie():
    elapsed_time, X, t = gillespie(DB_X0, DB_rates,
                                   DB_stoch_subst,
                                   DB_stoch_prods,
                                   DB_t_max)
    selected = None

    print(
        f"Elapsed time: {elapsed_time}ms"
    )

    fig, ax = plt.subplots()
    if selected is None:
        for i in range(len(DB_X0)):
            ax.plot(t, X[:, i], '-*', label='X' + str(i))
    else:
        for i in selected:
            ax.plot(t, X[:, i], '-*', label='X' + str(i))
    ax.legend(loc='best', shadow=True)
    plt.show()

# call_base_gillespie()

# Multiple gillespies to compare rate constant.


def call_multi_k_gillespie():
    rates = np.linspace(0.1, 1, 10)
    fig, ax = plt.subplots()
    for rate in rates:
        e, X, t = gillespie(DB_X0,
                            [rate],
                            DB_stoch_subst,
                            DB_stoch_prods,
                            30)
        ax.plot(t, X[:, 0], label=f'k = {rate :.1f}')
    ax.legend(loc='best', shadow=True)
    plt.show()


# call_multi_k_gillespie()

# Degradation time vs k.


def g(n, n_0, k, t):
    return (binom(n_0, n)
            * np.exp(-k * n * t)
            * np.power((1 - np.exp(-k * t)), (n_0 - n)))


def call_runtime_against_k():
    iters = 10000
    fig, ax = plt.subplots()
    end_vals = []
    for i in range(iters):
        e, X, t = gillespie(DB_X0,
                            DB_rates,
                            DB_stoch_subst,
                            DB_stoch_prods,
                            5)
        end_vals.append(X[-1, 0])
    max_val = int(np.max(end_vals))
    p_points = np.linspace(-0.5, max_val + 0.5, 100)
    bins = np.array([i for i in range(max_val + 2)]) - 0.5
    ax.set_xlabel("Molecule number")
    ax.set_ylabel("Probability")
    ax.set_title("Distribution of Molecules")
    pdf = ax.hist(end_vals,
                  density=True,
                  range=(0, max_val + 1),
                  bins=bins,
                  ec=(0, 0, 0, 0.5),
                  fc=(0, 0, 0, 0))
    mu_1 = DB_X0[0] * np.exp(- DB_rates[0] * 5)
    mu_2 = sum([a * b for a, b in enumerate(pdf[0])])
    ax.hist(end_vals,
            density=True,
            range=(0, max_val + 1),
            bins=bins,
            histtype='step',
            linewidth=1.5,
            ec='black',
            label='Realised\n distribution')
    ax.plot(p_points,
            g(p_points,
              DB_X0[0],
              DB_rates[0],
              5),
            linewidth=1.5,
            c='red',
            label='Theoretical\n distribution'
            )
    print(mu_1, mu_2)
    ax.axvline(x=mu_1,
               linewidth=1,
               c='red')
    ax.axvline(x=mu_2,
               linewidth=1.5,
               linestyle='--',
               c='black')
    ax.legend(loc='upper left',
              bbox_to_anchor=(0.53, 1),
              fontsize="16")
    ax.locator_params(axis='x', nbins=int(max_val/4))
    fig.subplots_adjust(hspace=0.2)
    pdf_compare = np.array([i for i in range(max_val + 1)])
    pdf_compare = g(pdf_compare,
                    DB_X0[0],
                    DB_rates[0],
                    5)
    print(entropy(pdf[0], pdf_compare))


call_runtime_against_k()


def f(t, n_0, k):
    return n_0 * np.exp(-k * t)


def call_mean_trajectory():

    stoptime = 10.0
    numpoints = 250
    t_new = [stoptime * float(i) / (numpoints - 1)
             for i in range(numpoints)]
    iters = 10

    def interpolate(t_vals, x_vals):
        t_vals = np.append(t_vals, 100)
        x_vals = np.append(x_vals, 0)
        inter_func = interp1d(t_vals, x_vals)
        return inter_func(np.array(t_new))

    all_sols = np.zeros((iters, len(t_new)))

    fig, ax = plt.subplots()

    plt.rcParams.update({'font.size': 20})
    ax.set_xlabel('time / s')
    ax.set_ylabel('Molecules of A')
    ax.set_title('Trajectories of A over time')
    for i in range(iters):
        e, X, t = gillespie(DB_X0,
                            DB_rates,
                            DB_stoch_subst,
                            DB_stoch_prods,
                            DB_t_max)
        ax.step(t, X[:, 0],
                linewidth=0.5)
        all_sols[i, :] = interpolate(t,
                                     X[:, 0])

    all_sols = np.mean(all_sols, axis=0)

    ax.plot(t_new, all_sols,
            '--',
            c='black',
            linewidth=1,
            label='Simulation\n mean')

    times = np.linspace(0, DB_t_max, 100)
    ax.plot(times,
            f(times, DB_X0[0], DB_rates[0]),
            linewidth=1.5,
            c="black",
            alpha=0.75,
            label='Theoretical\n trajectory')
    ax.locator_params(axis='x', nbins=6)
    ax.legend(loc='best', shadow=False)
    fig.subplots_adjust(hspace=0.2)
    plt.show()


call_mean_trajectory()
