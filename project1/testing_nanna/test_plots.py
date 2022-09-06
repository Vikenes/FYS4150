import numpy as np
import matplotlib.pyplot as plt


n_steps = [10, 100, 1000]
path = "folder/"

# Problem 7

fig, ax = plt.subplots()

for n in n_steps:
    file = f'num_sol_{n}steps.txt'
    x, v = np.loadtxt(path + file, unpack=True)
    ax.plot(x, v, ls='--', label=r'$n_{steps}=%i$'%n)

u = 1 - (1-np.exp(-10))*x - np.exp(-10*x)
ax.plot(x, u, lw=1.2, c='k', label=r'$u(x)$', alpha=0.6)

ax.legend()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$v(x)$')
plt.show()

# Problem 8

