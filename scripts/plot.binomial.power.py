from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.binom
import statsmodels.stats
read_size = range(0,50000,1)#read size
effect_size = np.arange(0.5,0.99,0.01) # transpose
pwoer =  = np.cos(x ** 2 + y ** 2)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(x, y, z,cmap='viridis', edgecolor='none')
ax.set_title('Surface plot')
plt.show()