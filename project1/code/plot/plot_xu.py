import numpy as np 
import matplotlib.pyplot as plt 
import os 

here = os.path.abspath(".")
data_path = here + "/../../output/data/"
plot_path = here + "/../../output/plots/"
# print(here)
# exit()

x, v = np.loadtxt(data_path + "x_u.txt", unpack=True)

file = plot_path + "u_x.pdf"

fig, ax = plt.subplots()
ax.plot(x,v)
fig.savefig(file)