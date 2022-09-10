from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 

xm, am = np.loadtxt("test_max_relative_error.txt", unpack=True, delimiter=',')

n6, gt6 = np.loadtxt("general_thomas_timed6.txt", unpack=True, delimiter=",")
n6, st6 = np.loadtxt("special_thomas_timed6.txt", unpack=True, delimiter=",")

n7, gt7 = np.loadtxt("general_thomas_timed7.txt", unpack=True, delimiter=",")
n7, st7 = np.loadtxt("special_thomas_timed7.txt", unpack=True, delimiter=",")

# plt.plot(xm, am, 'bo')

# plt.plot(n_s, gt, 'ro', label='gt')
# plt.plot(n_s, st, 'bo', label='st')

plt.plot(n6, gt6, 'rx', label='gt6')
plt.plot(n6, st6, 'bx', label='st6')


plt.plot(n7, gt7, 'gx', label='gt7')
plt.plot(n7, st7, 'kx', label='st7')


plt.xscale('log')
# plt.yscale('log')
plt.legend()
plt.show()