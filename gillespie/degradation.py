import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Degradation of Biomolecule:
DB_stoch_subst = np.array([[1]])
DB_stoch_prods = np.array([[0]])
DB_rates = np.array([0.5])
DB_X0 = np.array([1000])
DB_t_max = 10

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


call_multi_k_gillespie()

# Degradation time vs k.


def call_runtime_against_k():
    iters = 10
    rates = np.linspace(0.1, 1, 100)
    fig, ax = plt.subplots()
    times = []
    for rate in rates:
        k_time = []
        for i in range(iters):
            e, X, t = gillespie(DB_X0,
                                [rate],
                                DB_stoch_subst,
                                DB_stoch_prods,
                                100)
            k_time.append(t[-1])
        times.append(np.average(k_time))
    ax.plot(rates, times)
    plt.show()


# call_runtime_against_k()
