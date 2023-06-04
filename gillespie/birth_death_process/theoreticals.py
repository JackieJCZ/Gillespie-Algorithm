import numpy as np
import matplotlib.pyplot as plt

k1, k2 = np.mgrid[0.1:1.01:0.01, 0.1:1.01:0.01]

D = k1 / k2

plt.imshow(D,
           origin='lower',
           extent=(0, 1, 0, 1))
plt.colorbar()
plt.show()
