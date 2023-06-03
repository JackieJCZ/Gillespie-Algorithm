import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Gene Expression:
GE_stoch_subst = np.array([[]])
GE_stoch_prods = np.array([[]])
GE_rates = np.array([])
GE_X0 = np.array([])
GE_t_max = 10000

elapsed_time, X, t = gillespie(GE_X0, GE_rates,
                               GE_stoch_subst,
                               GE_stoch_prods,
                               GE_t_max)

selected = None

print(
    f"Elapsed time: {elapsed_time}ms"
)

fig, ax = plt.subplots()
if selected is None:
    for i in range(len(GE_X0)):
        ax.plot(t, X[:, i], '-*', label='X' + str(i))
else:
    for i in selected:
        ax.plot(t, X[:, i], '-*', label='X' + str(i))
ax.legend(loc='best', shadow=True)
plt.show()
