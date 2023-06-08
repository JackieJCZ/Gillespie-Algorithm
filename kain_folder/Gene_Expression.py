import numpy as np
from scipy.special import factorial, binom
import matplotlib.pyplot as plt
import math
import scipy


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


def my_gillespie(init, rates, stoch_subst, stoch_prods, tmax, nrmax):
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


genesub = np.array([[-1, 0, 0], [0, -1, 0], [0, -1, 0], [0, 0, -1]])
genepro = np.array([[1, 1, 0], [0, 1, 1], [0, 0, 0], [0, 0, 0]])
rates = np.array([0.5, 0.5, 0.3, 0.4])
init = np.array([100, 1, 0])
tmax = 30
nrmax = 100
store_t, store_mols, store_r = my_gillespie(init, rates, genesub,
                                            genepro, tmax, nrmax)
plt.step(store_t, store_mols, label=['DNA', 'RNA', 'protein'])
plt.xlabel('time / s')
plt.ylabel('molecule number')
plt.legend(title='Populations:', fontsize=12)
plt.title('Stochastic occurance of one reaction', fontsize=22)
plt.rcParams.update({'font.size': 15})
plt.show()
# Iterate 10 times for average
iter = 10
store_DNA = np.zeros(iter)
store_RNA = np.zeros(iter)
store_protein = np.zeros(iter)
fig, ax = plt.subplots(1, 2,
                       figsize=(10, 5),
                       sharex=True,
                       sharey=True)
for i in range(0, iter):
    t_temp, mol_temp, r_temp = my_gillespie(init, rates,
                                            genesub, genepro,
                                            tmax, nrmax)
    ax[0].step(t_temp, mol_temp[:, 1], '-', linewidth=0.5)
    ax[0].set_xlabel('time / s')
    ax[0].set_ylabel('molecule number')
    ax[1].step(t_temp, mol_temp[:, 2], '-', linewidth=0.5)
    ax[1].set_xlabel('time / s')


def F(t, x):
    """Rate equation of the gene expression."""
    ret = np.zeros(len(init))
    for i in range(0, len(genepro)):
        ret -= genesub[i].dot(x)*rates[i]*(genepro[i]+genesub[i])
    return ret


# Solve the rate equation and plot with the simulated results
results = scipy.integrate.solve_ivp(F, (0, tmax), init)
ax[0].plot(results.t, results.y[1], '-', c='black', label='mRNA',
           linewidth=2.0)
ax[0].legend(loc='lower right')
ax[1].plot(results.t, results.y[2], '-', c='black', label='protein',
           linewidth=2.0)
ax[1].legend(loc='lower right')
fig.subplots_adjust(hspace=0.5,
                    wspace=0)
fig.suptitle('Comparison of simulated molecule number with expected',
             fontsize=22)
plt.rcParams.update({'font.size': 15})
plt.show()
