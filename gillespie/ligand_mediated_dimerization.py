import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Ligand-Mediated Dimerization
LMD_stoch_subst = np.array([[1, 1, 0, 0, 0],
                            [0, 0, 1, 0, 0],
                            [0, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 1]])
LMD_stoch_prods = np.array([[0, 0, 1, 0, 0],
                            [1, 1, 0, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 1],
                            [0, 0, 0, 1, 0]])
LMD_rates = np.array([0.01, 0.1, 0.001, 10, 10, 1])
LMD_X0 = np.array([10, 100, 0, 0, 0])
LMD_t_max = 100


elapsed_time, X, t = gillespie(LMD_X0,
                               LMD_rates,
                               LMD_stoch_subst,
                               LMD_stoch_prods,
                               LMD_t_max)
selected = [1, 2, 3]

print(
    f"Elapsed time: {elapsed_time}ms"
)

fig, ax = plt.subplots()
if selected is None:
    for i in range(len(LMD_X0)):
        ax.plot(t, X[:, i], '-*', label='X' + str(i))
else:
    for i in selected:
        ax.plot(t, X[:, i], '-*', label='X' + str(i))
ax.legend(loc='best', shadow=True)
plt.show()
