import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
from GA import gillespie
from Theoretical_Func import p_mn


# Initial condtions
genesub = np.array([[-1, 0, 0], [0, -1, 0], [0, -1, 0], [0, 0, -1]])
genepro = np.array([[1, 1, 0], [0, 1, 1], [0, 0, 0], [0, 0, 0]])
rates = np.array([0.5, 0.5, 0.3, 0.4])
init = np.array([100, 1, 0])
tmax = 30
nrmax = 100
# Iterate for a given times for comparison with expected values
iter = 10000
stored_val = np.zeros((iter, 2))

fig, ax = plt.subplots(1, 2,
                       figsize=(10, 5),
                       sharex=True,
                       sharey=True)

for i in range(0, iter):
    t_temp, mol_temp, r_temp = gillespie(init, rates,
                                         genesub, genepro,
                                         tmax, nrmax)
    stored_val[i, ] = [mol_temp[-1, 1], mol_temp[-1, 2]]
    ax[0].step(t_temp, mol_temp[:, 1], '-', linewidth=0.5)
    ax[0].set_xlabel('time / s')
    ax[0].set_ylabel('molecule number')
    ax[1].step(t_temp, mol_temp[:, 2], '-', linewidth=0.5)
    ax[1].set_xlabel('time / s')


def RRE(t, x):
    """Rate equation of the gene expression."""
    ret = np.zeros(len(init))
    for i in range(0, len(genepro)):
        ret -= genesub[i].dot(x)*rates[i]*(genepro[i]+genesub[i])
    return ret


# Solve the rate equation and plot with the simulated results
results = scipy.integrate.solve_ivp(RRE, (0, tmax), init)
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

# Solve the CME for m and n
start_n, end_n, start_m, end_m = 155, 255, 135, 205
m_vals = np.arange(start_m, end_m + 1, 1)
n_vals = np.arange(start_n, end_n + 1, 1)
m, n = np.meshgrid(m_vals, n_vals)
p_mn = p_mn(rates, start_n, end_n, start_m, end_m)
p_n, p_m = scipy.stats.contingency.margins(p_mn)


def getkw(mi, ma):
    '''Get keyword arguments for plots with given extreme values.'''
    bins = np.array([i for i in range(mi, ma + 5, 5)]) - 2.5
    return dict(density=True, linewidth=1.5, bins=bins, ec='black')


# Plot the marginal distribution of mRNA and compare with simulation
min_m = int(min(m_vals))
max_m = int(max(m_vals))
kwargs = getkw(min_m, max_m)
pdf_m = plt.hist(stored_val[:, 0], fc='white', **kwargs)
plt.hist(stored_val[:, 0], label='Realised\n distribution',
         histtype='step', **kwargs)
mu_m = sum([a * b for a, b in enumerate(pdf_m[0])])
plt.axvline(x=mu_m, linewidth=1, c='black')
plt.plot(m_vals, p_m[0], label='Theoretical\n distribution',
         c='red', linewidth=1.5)
plt.xlabel('molecule number', fontsize=20)
plt.ylabel('probability', fontsize=20)
plt.xlim(min_m-5, max_m+5)
plt.legend(fontsize=12)
plt.title('Steady State Distribution of Simulated\nand Theoretical mRNA Level',
          fontsize=20)
plt.rcParams.update({'font.size': 20})
plt.show()

# Plot the marginal distribution of protein and compare with simulation
min_n = int(min(n_vals))
max_n = int(max(n_vals))
kwargs = getkw(min_n, max_n)
pdf_n = plt.hist(stored_val[:, 1], fc='white', **kwargs)
plt.hist(stored_val[:, 1], label='Realised\n distribution',
         histtype='step', **kwargs)
mu_n = sum([a * b for a, b in enumerate(pdf_n[0])])
plt.axvline(x=mu_n, linewidth=1, c='black')
plt.plot(n_vals, p_n.T[0], label='Theoretical\n distribution',
         c='red', linewidth=1.5)
plt.xlabel('molecule number', fontsize=20)
plt.ylabel('probability', fontsize=20)
plt.xlim(min_n-5, max_n+5)
plt.legend(fontsize=12)
plt.title('Steady State Distribution of Simulated\n'
          'and Theoretical Protein Level', fontsize=20)
plt.rcParams.update({'font.size': 20})
plt.show()

# Create the joint distribution plot for mRNA and protein using seaborn
bins = [[x - 0.5 for x in range(start_m, end_m + 2)],
        [x - 0.5 for x in range(start_n, end_n + 2)]]
p = sns.jointplot(x=stored_val[:, 0], y=stored_val[:, 1], kind="hist",
                  joint_kws=dict(bins=30))
plt.xlabel('mRNA number')
plt.ylabel('protein number')
p.fig.suptitle('Realised Joint Probability Distribution\n'
               'of mRNA and Protein')
p.fig.subplots_adjust(top=0.85)
plt.show()
# The expected plot of distributions using imshow
plt.imshow(p_mn,
           extent=[start_m, end_m, start_n, end_n], origin='lower')
plt.xlabel('mRNA number')
plt.ylabel('protein number')
plt.title('Theoretical Joint Probability Distribution\nof mRNA and Protein',
          fontsize=20)
plt.colorbar()
