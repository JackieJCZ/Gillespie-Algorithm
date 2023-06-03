import numpy as np
import matplotlib.pyplot as plt

k1, k2 = np.mgrid[0.1:1.1:0.01, 0.1:1.1:0.01]

D = (1 - k1 / k2) * k1 / k2

plt.imshow(D, origin='lower')
