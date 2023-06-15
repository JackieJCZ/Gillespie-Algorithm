# import modules and libraries
import numpy as np
from scipy.special import binom, factorial
from timeit import default_timer as timer


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


def gillespie(X0, rates, stoch_subst, stoch_prods,
              t_max, r_max=1000000):

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

    return timer() - start, X, t
