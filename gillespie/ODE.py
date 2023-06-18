from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt
from gillespie import gillespie
from scipy.stats import entropy

plt.rcParams.update({'font.size': 20})

# Ligand-Mediated Dimerization
LMD_stoch_subst = np.array([[1, 1],
                            [2, 0]])
LMD_stoch_prods = np.array([[2, 0],
                            [1, 1]])
LMD_rates = np.array([0.01, 0.1])
LMD_X0 = np.array([100, 0])
LMD_t_max = 1

n = 100  # total molecules
states = np.array([[n - i, i] for i in range(n + 1)])


def cme_func(x, t, k_a, k_d):
    func = np.zeros(n + 2)
    for i in range(n + 1):
        func[i] = (k_a * (n - 1 - i) * (i + 1) * x[i + 1]
                   - (k_a * i + k_d * (n - 1 - i)) * (n - i) * x[i]
                   + k_d * (n - i) * (n + 1 - i) * x[i - 1])
    return func


k_a = 0.01
k_d = 0.1
ts = np.linspace(0, LMD_t_max, 101)
y0 = np.zeros(n + 2)
y0[0] = 1

res = odeint(cme_func, y0, ts, args=(k_a, k_d))

fig, ax = plt.subplots()
ax.plot(np.array([i for i in range(n + 1)]),
        res[-1, -2::-1],
        'ro',
        ms=3,
        mec='r',
        label='Theoretical\n distribution')
ax.vlines(np.array([i for i in range(n + 1)]),
          0,
          res[-1, -2::-1],
          colors='r',
          linestyles='-',
          lw=1,
          alpha=0.5)

iters = 100
end_vals = []
for i in range(iters):
    e, X, t = gillespie(LMD_X0,
                        LMD_rates,
                        LMD_stoch_subst,
                        LMD_stoch_prods,
                        LMD_t_max)
    end_vals.append(X[-1, 0])
max_val = int(np.max(end_vals))
p_points = np.linspace(-0.5, max_val + 0.5, 100)
bins = np.array([i for i in range(max_val + 2)]) - 0.5
ax.set_xlabel("Dimer number")
ax.set_ylabel("Probability")
ax.set_title("Distribution of Dimers")
pdf = ax.hist(end_vals,
              density=True,
              range=(0, max_val + 1),
              bins=bins,
              ec=(0, 0, 0, 0.5),
              fc=(0, 0, 0, 0))
ax.hist(end_vals,
        density=True,
        range=(0, max_val + 1),
        bins=bins,
        histtype='step',
        linewidth=1.5,
        ec='black',
        label='Realised\n distribution')
ax.set_xlim([-5, 30])
ax.legend(loc='upper right',
          bbox_to_anchor=(0.47, 0.99),
          fontsize=16)

results = res[-1, -2::-1]
new_pdf = pdf[0]
new_pdf += (new_pdf == 0) * 0.000001
new_pdf = list(pdf[0]) + [0.000001 for i in range(len(results)
                                                  - len(pdf[0]))]
results += (results < 0.000001) * 0.000001
print(entropy(new_pdf,
              results))


solve_mat = np.zeros((n, n))
for i in range(1, n - 1):
    solve_mat[i, i + 1] = k_a * (n - 1 - i) * (i + 1)
    solve_mat[i, i] = - (k_a * i + k_d * (n - 1 - i)) * (n - i)
    solve_mat[i, i - 1] = k_d * (n - i) * (n + 1 - i)
solve_mat[0, 0] = k_d * (n - 1) * n
solve_mat[0, 1] = k_a * (n - 1)
solve_mat[n - 1, :] = np.ones(n)

res = np.zeros(n)
res[-1] = 1

stat_vals = np.linalg.solve(solve_mat, res)
stat_vals = list(stat_vals) + [0]

fig2, ax2 = plt.subplots()
ax2.plot(np.array([i for i in range(n + 1)]),
         stat_vals[::-1],
         'ro',
         ms=3,
         mec='r',
         label='Theoretical\n distribution')
ax2.vlines(np.array([i for i in range(n + 1)]),
           0,
           stat_vals[::-1],
           colors='r',
           linestyles='-',
           lw=1,
           alpha=0.5)
ax2.set_xlim([-5, 30])

end_vals_2 = []
for i in range(iters):
    e, X, t = gillespie(LMD_X0,
                        LMD_rates,
                        LMD_stoch_subst,
                        LMD_stoch_prods,
                        10)
    end_vals_2.append(X[-1, 0])
max_val_2 = int(np.max(end_vals_2))
p_points_2 = np.linspace(-0.5, max_val_2 + 0.5, 100)
bins_2 = np.array([i for i in range(max_val_2 + 2)]) - 0.5
ax2.set_xlabel("Dimer number")
ax2.set_ylabel("Probability")
ax2.set_title("Distribution of Dimers")
pdf2 = ax2.hist(end_vals_2,
                density=True,
                range=(0, max_val_2 + 1),
                bins=bins_2,
                ec=(0, 0, 0, 0.5),
                fc=(0, 0, 0, 0))
ax2.hist(end_vals_2,
         density=True,
         range=(0, max_val_2 + 1),
         bins=bins_2,
         histtype='step',
         linewidth=1.5,
         ec='black',
         label='Realised\n distribution')
ax2.set_xlim([-5, 30])
ax2.legend(loc='upper left',
           bbox_to_anchor=(0.47, 0.99),
           fontsize=18)

results_2 = np.array(stat_vals[::-1])
new_pdf_2 = pdf2[0]
new_pdf_2 += (new_pdf_2 == 0) * 0.000001
new_pdf_2 = list(pdf2[0]) + [0.000001 for i in range(len(results_2)
                                                     - len(pdf2[0]))]
results_2 += (results_2 < 0.000001) * 0.000001
print(entropy(new_pdf_2,
              results_2))