import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Gene Expression:

# TOGGLE SWITCH:

GE_stoch_subst = np.array([[1, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, 0, 1, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, 0, 1, 0, 0],
                           [0, 0, 0, 0, 1, 0],
                           [0, 0, 0, 0, 0, 1]])
GE_stoch_prods = np.array([[1, 0, 1, 0, 0, 0],
                           [0, 1, 0, 1, 0, 0],
                           [0, 0, 1, 0, 1, 0],
                           [0, 0, 0, 1, 0, 1],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0]])
GE_rates = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
GE_X0 = np.array([100, 100, 1000, 10, 100, 100])
GE_t_max = 1000

elapsed_time, X, t = gillespie(GE_X0, GE_rates,
                               GE_stoch_subst,
                               GE_stoch_prods,
                               GE_t_max)

selected = [2, 3]

print(
    f"Elapsed time: {elapsed_time}ms"
)

fig, ax = plt.subplots()
if selected is None:
    for i in range(len(GE_X0)):
        ax.plot(t, X[:, i], label='X' + str(i))
else:
    for i in selected:
        ax.plot(t, X[:, i], label='X' + str(i))
ax.legend(loc='best', shadow=True)
plt.show()
