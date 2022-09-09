import numpy as np 
import matplotlib.pyplot as plt 

xm, am = np.loadtxt("test_max_relative_error.txt", unpack=True, delimiter=',')

plt.plot(xm, am, 'bo')

plt.xscale('log');plt.yscale('log')
plt.legend()
plt.show()