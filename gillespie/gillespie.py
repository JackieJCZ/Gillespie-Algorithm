# import modules and libraries
import numpy as np
from scipy.special import binom, factorial
import matplotlib.pyplot as plt
from timeit import default_timer as timer

# ****************************
# Input
# ****************************
# Degradation of Biomolecule:
DB_stoch_subst = np.array([[2, 0, 0, 0]])
DB_stoch_prods = np.array([[0, 1, 3, 2]])
DB_rates = np.array([0.1])
DB_X0 = np.array([1000, 0, 0, 0])
DB_t_max = 1
# ****************************
# Birth - Death Process:
BD_stoch_subst = np.array([[1, 0],
                           [0, 1],
                           [0, 1]])
BD_stoch_prods = np.array([[0, 1],
                           [1, 0],
                           [0, 0]])
BD_rates = np.array([0.5, 0.5, 0.01])
BD_X0 = np.array([1000, 20])
BD_t_max = 10
# ****************************
# Protein Complex Formation:
PCF_stoch_subst = np.array([[]])
PCF_stoch_prods = np.array([[]])
PCF_rates = np.array([])
PCF_X0 = np.array([])
PCF_t_max = 10000
# ****************************
# Ligand-Mediated Dimerization
LMD_stoch_subst = np.array([[]])
LMD_stoch_prods = np.array([[]])
LMD_rates = np.array([])
LMD_X0 = np.array([])
LMD_t_max = 10000
# ****************************
# Gene Expression:
GE_stoch_subst = np.array([[]])
GE_stoch_prods = np.array([[]])
GE_rates = np.array([])
GE_X0 = np.array([])
GE_t_max = 10000
# ****************************

# IGNORE #
def new_bin(i, j):
    if np.any(i <= 0):
        return 0
    else:
        return binom(i, j)


def propensity_calc_2(X, M, N):
    a = np.zeros(M)
    return a
# IGNORE #


def propensity_calc(X, M, N, stoch_subst, rates):
    a = np.zeros(M)
    for i in range(M):
        h_i = 1
        for j in range(N):
            if stoch_subst[i, j] == 0:
                continue
            else:
                if X[j] <= 0:
                    h_i = 0
                    break
                else:
                    h_i = (h_i * binom(X[j], stoch_subst[i, j])
                           * factorial(stoch_subst[i, j]))
        a[i] = h_i * rates[i]
    a0 = sum(a)
    return a, a0


def find_τ_r(a, a0):
    r1, r2 = np.random.random(2)
    τ = (1 / a0) * np.log(1 / r1)  # Pick a τ
    µ = 0                          # Pick a reaction
    N = r2 * a0 - a[µ]
    while N > 0:
        µ += 1
        N -= a[µ]
    return τ, µ


def gillespie(X0, rates, stoch_subst, stoch_prods, t_max, r_max=1000000):

    start = timer()

    stoch = stoch_prods - stoch_subst
    M, N = np.shape(stoch)

    X = np.zeros((r_max, N))
    X[0] = X0
    t = [0]    # The vector of times
    r_num = 0  # The number of reactions computed

    while t[-1] <= t_max:
        current_species = X[r_num]
        a, a0 = propensity_calc(current_species,
                                M,
                                N,
                                stoch_subst,
                                rates)
        τ, r = find_τ_r(a, a0)

        t += [t[-1] + τ]
        X[r_num + 1] = current_species + np.transpose(stoch[r])
        r_num += 1

    X = X[:r_num]
    t = np.array(t)[:r_num]

    print(
        f"Elapsed time: {timer() - start}ms"
        )

    fig, ax = plt.subplots()
    for i in range(N):
        ax.plot(t, X[:, i], '-*', label='X' + str(i))
    ax.legend(loc='best', shadow=True)
    plt.show()


gillespie(DB_X0, DB_rates, DB_stoch_subst, DB_stoch_prods, DB_t_max)
gillespie(BD_X0, BD_rates, BD_stoch_subst, BD_stoch_prods, BD_t_max)
# gillespie(PCF_X0, PCF_rates, PCF_stoch_subst, PCF_stoch_prods, PCF_t_max)
# gillespie(LMD_X0, LMD_rates, LMD_stoch_subst, LMD_stoch_prods, LMD_t_max)
# gillespie(GE_X0, GE_rates, GE_stoch_subst, GE_stoch_prods, GE_t_max)
