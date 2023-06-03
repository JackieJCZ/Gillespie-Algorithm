import numpy as np
import matplotlib.pyplot as plt
from gillespie import gillespie

# Protein Complex Formation:
PCF_stoch_subst = np.array([[]])
PCF_stoch_prods = np.array([[]])
PCF_rates = np.array([])
PCF_X0 = np.array([])
PCF_t_max = 10000


elapsed_time, X, t = gillespie(PCF_X0, PCF_rates,
                               PCF_stoch_subst,
                               PCF_stoch_prods,
                               PCF_t_max)
selected = None

print(
    f"Elapsed time: {elapsed_time}ms"
)

fig, ax = plt.subplots()
if selected is None:
    for i in range(len(PCF_X0)):
        ax.plot(t, X[:, i], '-*', label='X' + str(i))
else:
    for i in selected:
        ax.plot(t, X[:, i], '-*', label='X' + str(i))
ax.legend(loc='best', shadow=True)
plt.show()
