import numpy as np
from scipy.special import factorial, binom
import math


def choose_t_r(a0, a):
    """Choose the next time and reaction happening randomly."""
    r1 = np.random.random()
    r2 = np.random.random()
    T = (1/a0)*math.log(1/r1)

    mu = 0
    N = r2*a0 - a[mu]
    while (N > 0):
        mu += 1
        N -= a[mu]
    next_r = mu
    return T, next_r


def gillespie(init, rates, stoch_subst, stoch_prods, tmax, nrmax):
    """The main calculation of the Gillespie alogrithm."""
    stoch = stoch_prods + stoch_subst
    num_rxn, num_spec = np.shape(stoch)
    current_species = init.copy()
    current_t = 0
    t_count = 0
    largenum = 1000000
    store_t = np.zeros(largenum)
    store_mols = np.zeros((largenum, num_spec))
    store_r = np.zeros(largenum)
    store_t[t_count] = current_t
    store_mols[t_count, :] = current_species
    # The main loop over the time given
    while (current_t < tmax):
        a = np.ones(num_rxn)
        for i in range(num_rxn):
            hi = 1
            for j in range(len(init)):
                if (stoch_subst[i, j]):
                    if (current_species[j]):
                        hi = hi*binom(current_species[j],
                                      np.absolute(
                                      stoch_subst[i, j]))*factorial(
                                      np.absolute(stoch_subst[i, j]))
                    else:
                        hi = 0
            a[i] = hi*rates[i]
        a0 = sum(a)
        T, next_r = choose_t_r(a0, a)
        current_t += T
        current_species += np.transpose(stoch[next_r, :])
        t_count += 1
        store_t[t_count] = current_t
        store_mols[t_count, :] = current_species
        store_r[t_count] = next_r
    # Store the values accordingly to return
    store_t = store_t[:t_count]
    store_mols = store_mols[:t_count, :]
    store_r = store_r[:t_count]
    return store_t, store_mols, store_r
