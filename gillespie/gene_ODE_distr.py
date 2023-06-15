from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from scipy.special import factorial
from scipy.stats import poisson

m = np.linspace(1, 200, 200)
n = np.linspace(1, 200, 200)
m, n = np.meshgrid(m, n)

k_2 = 0.5
k_4 = 0.4

p_points = poisson.pmf(n, k_2 * m / k_4)


def p_n(m, n):
    return (1 / factorial(n)
            * np.power(k_2 * m / k_4, n)
            * np.exp(- k_2 * m / k_4))


fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(m, n, p_points,
                       cmap=plt.cm.cividis)
fig.colorbar(surf, shrink=0.5, aspect=8)
plt.show()


plt.imshow(p_points)
plt.colorbar()
