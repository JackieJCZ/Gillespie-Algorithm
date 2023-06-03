import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Degradation of Biomolecule:
DB_stoch_subst = np.array([[1]])
DB_stoch_prods = np.array([[0]])
DB_rates = np.array([0.5])
DB_X0 = np.array([1000])
DB_t_max = 10


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
